"""
ORF Finder Module

Run Prodigal to predict ORFs within IS elements and identify noncoding regions.
"""

import os
import subprocess
import tempfile
from typing import Dict, List, Optional, Tuple
from pathlib import Path


def run_prodigal(fasta_path: str, output_prefix: str = None,
                 meta_mode: bool = True, quiet: bool = True) -> Tuple[str, str]:
    """
    Run Prodigal on a FASTA file.

    Args:
        fasta_path: Path to input FASTA file
        output_prefix: Prefix for output files (default: same as input)
        meta_mode: Use metagenomic mode (-p meta) for short sequences
        quiet: Suppress Prodigal output

    Returns:
        Tuple of (gff_path, protein_path)
    """
    if output_prefix is None:
        output_prefix = fasta_path.rsplit('.', 1)[0]

    gff_path = f"{output_prefix}.gff"
    protein_path = f"{output_prefix}.faa"
    genes_path = f"{output_prefix}.fna"

    cmd = [
        "prodigal",
        "-i", fasta_path,
        "-o", gff_path,
        "-f", "gff",           # GFF output format
        "-a", protein_path,    # Protein translations
        "-d", genes_path,      # Gene nucleotide sequences
    ]

    if meta_mode:
        cmd.extend(["-p", "meta"])  # Metagenomic mode for short sequences

    if quiet:
        cmd.extend(["-q"])

    # Run prodigal (assumes module is loaded)
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"Prodigal failed: {result.stderr}")

    return gff_path, protein_path


def parse_gff(gff_path: str) -> List[Dict]:
    """
    Parse Prodigal GFF output.

    Args:
        gff_path: Path to GFF file

    Returns:
        List of gene dictionaries with start, end, strand, attributes
    """
    genes = []

    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            if not line.strip():
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            seqid, source, feature, start, end, score, strand, phase, attributes = parts

            if feature != 'CDS':
                continue

            # Parse attributes
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key] = value

            gene = {
                'seqid': seqid,
                'start': int(start),
                'end': int(end),
                'strand': strand,
                'score': float(score) if score != '.' else None,
                'phase': int(phase) if phase != '.' else 0,
                'gene_id': attr_dict.get('ID', ''),
                'partial': attr_dict.get('partial', '00'),
                'start_type': attr_dict.get('start_type', ''),
                'rbs_motif': attr_dict.get('rbs_motif', ''),
                'rbs_spacer': attr_dict.get('rbs_spacer', ''),
                'gc_cont': float(attr_dict.get('gc_cont', 0)),
                'conf': float(attr_dict.get('conf', 0)),
            }
            genes.append(gene)

    return genes


def parse_protein_fasta(fasta_path: str) -> Dict[str, Dict]:
    """
    Parse Prodigal protein FASTA output.

    Args:
        fasta_path: Path to protein FASTA file

    Returns:
        Dictionary mapping gene_id to protein info (sequence, header details)
    """
    proteins = {}
    current_id = None
    current_seq = []
    current_info = {}

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous
                if current_id:
                    proteins[current_id] = {
                        'sequence': ''.join(current_seq),
                        'length': len(''.join(current_seq)),
                        **current_info
                    }

                # Parse header: >seqid_genenum # start # end # strand # info
                parts = line[1:].split(' # ')
                current_id = parts[0]
                current_info = {
                    'start': int(parts[1]) if len(parts) > 1 else 0,
                    'end': int(parts[2]) if len(parts) > 2 else 0,
                    'strand': int(parts[3]) if len(parts) > 3 else 1,
                    'info': parts[4] if len(parts) > 4 else ''
                }
                current_seq = []
            else:
                current_seq.append(line)

        # Don't forget last entry
        if current_id:
            proteins[current_id] = {
                'sequence': ''.join(current_seq),
                'length': len(''.join(current_seq)),
                **current_info
            }

    return proteins


def calculate_noncoding_regions(is_length: int, genes: List[Dict]) -> List[Dict]:
    """
    Calculate noncoding regions given IS element length and gene positions.

    Args:
        is_length: Total length of IS element
        genes: List of gene dictionaries with 'start' and 'end'

    Returns:
        List of noncoding region dictionaries
    """
    if not genes:
        # Entire sequence is noncoding
        return [{
            'start': 1,
            'end': is_length,
            'length': is_length,
            'type': 'intergenic'
        }]

    # Sort genes by start position
    sorted_genes = sorted(genes, key=lambda x: x['start'])

    noncoding = []

    # Region before first gene
    first_gene_start = sorted_genes[0]['start']
    if first_gene_start > 1:
        noncoding.append({
            'start': 1,
            'end': first_gene_start - 1,
            'length': first_gene_start - 1,
            'type': '5_prime_utr'
        })

    # Intergenic regions between genes
    for i in range(len(sorted_genes) - 1):
        current_end = sorted_genes[i]['end']
        next_start = sorted_genes[i + 1]['start']

        if next_start > current_end + 1:
            noncoding.append({
                'start': current_end + 1,
                'end': next_start - 1,
                'length': next_start - current_end - 1,
                'type': 'intergenic'
            })

    # Region after last gene
    last_gene_end = sorted_genes[-1]['end']
    if last_gene_end < is_length:
        noncoding.append({
            'start': last_gene_end + 1,
            'end': is_length,
            'length': is_length - last_gene_end,
            'type': '3_prime_utr'
        })

    return noncoding


def extract_subsequence(sequence: str, start: int, end: int) -> str:
    """Extract subsequence using 1-based coordinates."""
    return sequence[start - 1:end]


def annotate_is_element(is_id: str, is_sequence: str,
                        temp_dir: str = None,
                        keep_files: bool = False) -> Dict:
    """
    Annotate a single IS element with ORF predictions and noncoding regions.

    Args:
        is_id: Identifier for this IS element
        is_sequence: DNA sequence of the IS element
        temp_dir: Directory for temporary files (default: system temp)
        keep_files: Keep intermediate files

    Returns:
        Dictionary with ORF and noncoding region annotations
    """
    if temp_dir is None:
        temp_dir = tempfile.gettempdir()

    # Create temporary FASTA file
    fasta_path = os.path.join(temp_dir, f"{is_id}.fna")
    output_prefix = os.path.join(temp_dir, f"{is_id}_prodigal")

    with open(fasta_path, 'w') as f:
        f.write(f">{is_id}\n")
        # Write sequence in 60-char lines
        for i in range(0, len(is_sequence), 60):
            f.write(is_sequence[i:i+60] + "\n")

    try:
        # Run Prodigal
        gff_path, protein_path = run_prodigal(fasta_path, output_prefix)

        # Parse results
        genes = parse_gff(gff_path)
        proteins = parse_protein_fasta(protein_path)

        # Filter genes for this IS element
        is_genes = [g for g in genes if g['seqid'] == is_id]

        # Calculate noncoding regions
        noncoding = calculate_noncoding_regions(len(is_sequence), is_genes)

        # Build ORF annotations with protein sequences
        orfs = []
        for gene in is_genes:
            gene_key = f"{is_id}_{len(orfs) + 1}"

            # Find matching protein
            protein_seq = ''
            for prot_id, prot_info in proteins.items():
                if prot_id.startswith(is_id) and prot_info['start'] == gene['start']:
                    protein_seq = prot_info['sequence']
                    break

            orf = {
                'orf_id': gene_key,
                'start': gene['start'],
                'end': gene['end'],
                'strand': gene['strand'],
                'length_nt': gene['end'] - gene['start'] + 1,
                'length_aa': len(protein_seq),
                'dna_sequence': extract_subsequence(is_sequence, gene['start'], gene['end']),
                'protein_sequence': protein_seq,
                'partial': gene['partial'],
                'start_type': gene['start_type'],
                'rbs_motif': gene['rbs_motif'],
                'gc_content': gene['gc_cont'],
                'confidence': gene['conf']
            }
            orfs.append(orf)

        # Add sequences to noncoding regions
        for nc in noncoding:
            nc['sequence'] = extract_subsequence(is_sequence, nc['start'], nc['end'])

        result = {
            'is_id': is_id,
            'is_length': len(is_sequence),
            'num_orfs': len(orfs),
            'num_noncoding_regions': len(noncoding),
            'total_coding_bp': sum(orf['length_nt'] for orf in orfs),
            'total_noncoding_bp': sum(nc['length'] for nc in noncoding),
            'coding_fraction': sum(orf['length_nt'] for orf in orfs) / len(is_sequence) if is_sequence else 0,
            'orfs': orfs,
            'noncoding_regions': noncoding
        }

        return result

    finally:
        # Cleanup
        if not keep_files:
            for ext in ['.fna', '_prodigal.gff', '_prodigal.faa', '_prodigal.fna']:
                path = os.path.join(temp_dir, f"{is_id}{ext}")
                if os.path.exists(path):
                    os.remove(path)


def annotate_is_elements(is_records: List[Dict], temp_dir: str = None) -> List[Dict]:
    """
    Annotate multiple IS elements with ORF predictions.

    Args:
        is_records: List of IS element dictionaries (from ISExtractor)
        temp_dir: Directory for temporary files

    Returns:
        List of annotated IS element dictionaries
    """
    results = []

    for record in is_records:
        is_id = record['is_id']
        is_sequence = record['is_element']['sequence']

        if not is_sequence:
            continue

        try:
            annotation = annotate_is_element(is_id, is_sequence, temp_dir)

            # Merge with original record
            annotated = {
                **record,
                'orf_annotation': annotation
            }
            results.append(annotated)

        except Exception as e:
            print(f"Warning: Failed to annotate {is_id}: {e}")
            results.append({
                **record,
                'orf_annotation': {'error': str(e)}
            })

    return results


if __name__ == "__main__":
    import json
    import sys

    # Example usage
    if len(sys.argv) < 2:
        print("Usage: python orf_finder.py <is_sequence>")
        print("       python orf_finder.py --test")
        sys.exit(1)

    if sys.argv[1] == "--test":
        # Test with a sample IS element
        test_seq = "ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA"
        result = annotate_is_element("test_is", test_seq * 10)
        print(json.dumps(result, indent=2))
    else:
        result = annotate_is_element("input_is", sys.argv[1])
        print(json.dumps(result, indent=2))
