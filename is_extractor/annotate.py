"""
IS Element Annotation Script

Annotate IS elements with:
1. Protein sequences (transposases) - with start, end, strand
2. Noncoding regions - with start, end, sequence

Usage:
    python annotate.py <sample_dir> [output.json]

Requires:
    module load bio/prodigal/2.6.3-gcc-11.4.0
"""

import os
import json
import subprocess
import tempfile
from typing import Dict, List, Optional
from pathlib import Path


def parse_fasta(fasta_path: str) -> Dict[str, str]:
    """Parse FASTA file into {header: sequence} dict."""
    sequences = {}
    current_header = None
    current_seq = []

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    sequences[current_header] = ''.join(current_seq)
                current_header = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            sequences[current_header] = ''.join(current_seq)

    return sequences


def run_prodigal(fasta_path: str, output_prefix: str) -> tuple:
    """Run Prodigal and return paths to output files."""
    gff_path = f"{output_prefix}.gff"
    protein_path = f"{output_prefix}.faa"
    genes_path = f"{output_prefix}.fna"

    cmd = [
        "prodigal",
        "-i", fasta_path,
        "-o", gff_path,
        "-f", "gff",
        "-a", protein_path,
        "-d", genes_path,
        "-p", "meta",  # metagenomic mode for short sequences
        "-q"           # quiet
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Prodigal failed: {result.stderr}")

    return gff_path, protein_path, genes_path


def parse_prodigal_output(gff_path: str, protein_path: str, is_id: str) -> List[Dict]:
    """
    Parse Prodigal GFF and protein FASTA to extract ORF info.

    Returns list of ORF dictionaries with:
        - orf_id, start, end, strand, protein_sequence
    """
    # Parse GFF for coordinates
    orfs = []
    with open(gff_path, 'r') as f:
        orf_num = 0
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'CDS':
                continue

            orf_num += 1
            orfs.append({
                'orf_id': f"{is_id}_orf{orf_num}",
                'start': int(parts[3]),
                'end': int(parts[4]),
                'strand': parts[6],
                'protein_sequence': ''  # Will fill from protein file
            })

    # Parse protein FASTA
    if os.path.exists(protein_path):
        proteins = parse_fasta(protein_path)
        # Match proteins to ORFs by position
        for orf in orfs:
            for prot_id, prot_seq in proteins.items():
                # Prodigal header format: >seqid_num # start # end # strand # info
                if f"# {orf['start']} # {orf['end']} #" in prot_id or prot_id.startswith(is_id):
                    # Try to match by parsing header
                    pass
            # Simpler: assign in order (Prodigal outputs in same order)

        # Just assign in order since Prodigal outputs match GFF order
        prot_list = list(proteins.values())
        for i, orf in enumerate(orfs):
            if i < len(prot_list):
                orf['protein_sequence'] = prot_list[i]

    return orfs


def calculate_noncoding(is_length: int, orfs: List[Dict]) -> List[Dict]:
    """
    Calculate noncoding regions from IS length and ORF positions.

    Returns list of noncoding region dictionaries with:
        - region_id, start, end, length, type
    """
    if not orfs:
        return [{
            'region_id': 'nc_1',
            'start': 1,
            'end': is_length,
            'length': is_length,
            'type': 'full_noncoding'
        }]

    # Sort ORFs by start position
    sorted_orfs = sorted(orfs, key=lambda x: x['start'])

    noncoding = []
    region_num = 0

    # 5' region before first ORF
    first_start = sorted_orfs[0]['start']
    if first_start > 1:
        region_num += 1
        noncoding.append({
            'region_id': f'nc_{region_num}',
            'start': 1,
            'end': first_start - 1,
            'length': first_start - 1,
            'type': '5_prime'
        })

    # Intergenic regions between ORFs
    for i in range(len(sorted_orfs) - 1):
        current_end = sorted_orfs[i]['end']
        next_start = sorted_orfs[i + 1]['start']

        if next_start > current_end + 1:
            region_num += 1
            noncoding.append({
                'region_id': f'nc_{region_num}',
                'start': current_end + 1,
                'end': next_start - 1,
                'length': next_start - current_end - 1,
                'type': 'intergenic'
            })

    # 3' region after last ORF
    last_end = sorted_orfs[-1]['end']
    if last_end < is_length:
        region_num += 1
        noncoding.append({
            'region_id': f'nc_{region_num}',
            'start': last_end + 1,
            'end': is_length,
            'length': is_length - last_end,
            'type': '3_prime'
        })

    return noncoding


def annotate_single_is(is_id: str, is_sequence: str,
                       temp_dir: str = None) -> Dict:
    """
    Annotate a single IS element.

    Args:
        is_id: Identifier for this IS element
        is_sequence: DNA sequence
        temp_dir: Directory for temp files

    Returns:
        Dictionary with 'proteins' and 'noncoding_regions'
    """
    if temp_dir is None:
        temp_dir = tempfile.gettempdir()

    is_length = len(is_sequence)

    # Write temporary FASTA
    fasta_path = os.path.join(temp_dir, f"{is_id}.fna")
    output_prefix = os.path.join(temp_dir, f"{is_id}_prodigal")

    with open(fasta_path, 'w') as f:
        f.write(f">{is_id}\n")
        for i in range(0, len(is_sequence), 60):
            f.write(is_sequence[i:i+60] + "\n")

    try:
        # Run Prodigal
        gff_path, protein_path, genes_path = run_prodigal(fasta_path, output_prefix)

        # Parse ORFs
        orfs = parse_prodigal_output(gff_path, protein_path, is_id)

        # Calculate noncoding regions
        noncoding = calculate_noncoding(is_length, orfs)

        # Add sequences to noncoding regions
        for nc in noncoding:
            nc['sequence'] = is_sequence[nc['start']-1:nc['end']]

        # Build result
        result = {
            'is_id': is_id,
            'is_length': is_length,
            'proteins': [
                {
                    'protein_id': orf['orf_id'],
                    'start': orf['start'],
                    'end': orf['end'],
                    'strand': orf['strand'],
                    'length_aa': len(orf['protein_sequence']),
                    'protein_sequence': orf['protein_sequence']
                }
                for orf in orfs
            ],
            'noncoding_regions': [
                {
                    'region_id': nc['region_id'],
                    'start': nc['start'],
                    'end': nc['end'],
                    'length': nc['length'],
                    'type': nc['type'],
                    'sequence': nc['sequence']
                }
                for nc in noncoding
            ]
        }

        return result

    finally:
        # Cleanup temp files
        for ext in ['.fna', '_prodigal.gff', '_prodigal.faa', '_prodigal.fna']:
            path = os.path.join(temp_dir, f"{is_id}{ext}")
            if os.path.exists(path):
                os.remove(path)


def annotate_is_records(is_records: List[Dict], temp_dir: str = None) -> List[Dict]:
    """
    Annotate multiple IS element records.

    Args:
        is_records: List of IS element dicts from extract_is_elements()
        temp_dir: Directory for temp files

    Returns:
        List of annotated records with 'annotation' field added
    """
    results = []

    for record in is_records:
        is_id = record['is_id']
        is_seq = record['is_element']['sequence']

        if not is_seq:
            continue

        try:
            annotation = annotate_single_is(is_id, is_seq, temp_dir)

            # Merge annotation into record
            annotated_record = {
                **record,
                'annotation': annotation
            }
            results.append(annotated_record)

        except Exception as e:
            print(f"Warning: Failed to annotate {is_id}: {e}")
            results.append({
                **record,
                'annotation': {'error': str(e)}
            })

    return results


def main():
    """Main function for command-line usage."""
    import sys

    if len(sys.argv) < 2:
        print(__doc__)
        print("\nExample:")
        print("  module load bio/prodigal/2.6.3-gcc-11.4.0")
        print("  python annotate.py /path/to/sample_dir output.json")
        sys.exit(1)

    sample_dir = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None

    # Import extractor
    from extractor import extract_is_elements

    print(f"Extracting IS elements from: {sample_dir}")
    is_records = extract_is_elements(sample_dir)
    print(f"Found {len(is_records)} IS elements")

    print("Running Prodigal annotation...")
    annotated = annotate_is_records(is_records)
    print(f"Annotated {len(annotated)} IS elements")

    # Output
    if output_file:
        with open(output_file, 'w') as f:
            json.dump(annotated, f, indent=2)
        print(f"Results saved to: {output_file}")
    else:
        # Print summary
        print("\n" + "=" * 70)
        print("ANNOTATION SUMMARY")
        print("=" * 70)

        for record in annotated[:5]:
            ann = record.get('annotation', {})
            if 'error' in ann:
                print(f"\n{record['is_id']}: ERROR - {ann['error']}")
                continue

            print(f"\n{record['is_id']} ({ann['is_length']} bp):")
            print(f"  Proteins: {len(ann['proteins'])}")
            for prot in ann['proteins']:
                print(f"    - {prot['protein_id']}: {prot['start']}-{prot['end']} ({prot['strand']}) {prot['length_aa']} aa")
            print(f"  Noncoding regions: {len(ann['noncoding_regions'])}")
            for nc in ann['noncoding_regions']:
                print(f"    - {nc['region_id']}: {nc['start']}-{nc['end']} ({nc['length']} bp) [{nc['type']}]")

        if len(annotated) > 5:
            print(f"\n... and {len(annotated) - 5} more")


if __name__ == "__main__":
    main()
