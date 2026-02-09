"""
IS Element Extractor

Extract IS element DNA sequences and flanking regions from MGEfinder output.
Flanking regions are extracted from the assembly contig where the IS sits,
ensuring upstream_flank → IS_element → downstream_flank are contiguous.
"""

import os
import sys
import csv
from typing import Dict, List, Optional
from pathlib import Path

# Some IS sequences are very long, increase CSV field size limit
csv.field_size_limit(sys.maxsize)


def parse_fasta(fasta_path: str) -> Dict[str, str]:
    """Parse a FASTA file into a dictionary of {header: sequence}."""
    sequences = {}
    current_header = None
    current_seq = []

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header is not None:
                    sequences[current_header] = ''.join(current_seq)
                # Extract just the sequence ID (first word after >)
                current_header = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

        # Don't forget the last sequence
        if current_header is not None:
            sequences[current_header] = ''.join(current_seq)

    return sequences


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    comp = str.maketrans('ACGTacgt', 'TGCAtgca')
    return seq.translate(comp)[::-1]


class ISExtractor:
    """Extract IS elements and flanking regions from MGEfinder output."""

    def __init__(self, sample_dir: str, flank_size: int = 80):
        """
        Initialize extractor for a sample directory.

        Args:
            sample_dir: Path to MGEfinder sample directory
                       (containing 00.genome, 00.assembly, 03.results, etc.)
            flank_size: Size of flanking regions to extract (default 80bp)
        """
        self.sample_dir = Path(sample_dir)
        self.flank_size = flank_size
        self.sample_id = self.sample_dir.name

        # Define file paths
        self.genome_path = self.sample_dir / "00.genome" / "genome.fna"
        self.assembly_path = self.sample_dir / "00.assembly" / f"{self.sample_id}.fna"
        self.clusterseq_path = self.sample_dir / "03.results" / "genome" / "01.clusterseq.genome.tsv"
        self.genotype_path = self.sample_dir / "03.results" / "genome" / "02.genotype.genome.tsv"
        self.fasta_path = self.sample_dir / "03.results" / "genome" / "04.makefasta.genome.repr_seqs.fna"

        # Cache for loaded data
        self._genome_seqs = None
        self._assembly_seqs = None
        self._is_sequences = None

    def _load_genome(self) -> Dict[str, str]:
        """Load reference genome sequences."""
        if self._genome_seqs is None:
            if not self.genome_path.exists():
                raise FileNotFoundError(f"Genome file not found: {self.genome_path}")
            self._genome_seqs = parse_fasta(str(self.genome_path))
        return self._genome_seqs

    def _load_assembly(self) -> Dict[str, str]:
        """Load assembly contig sequences."""
        if self._assembly_seqs is None:
            if not self.assembly_path.exists():
                self._assembly_seqs = {}
            else:
                self._assembly_seqs = parse_fasta(str(self.assembly_path))
        return self._assembly_seqs

    def _load_is_sequences(self) -> Dict[str, str]:
        """Load IS element sequences from FASTA."""
        if self._is_sequences is None:
            if self.fasta_path.exists():
                self._is_sequences = parse_fasta(str(self.fasta_path))
            else:
                self._is_sequences = {}
        return self._is_sequences

    def _parse_clusterseq(self) -> List[dict]:
        """Parse clusterseq TSV file."""
        records = []
        if not self.clusterseq_path.exists():
            return records

        with open(self.clusterseq_path, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                records.append(row)
        return records

    def _parse_genotype(self) -> Dict[str, dict]:
        """Parse genotype TSV file, keyed by seqid."""
        records = {}
        if not self.genotype_path.exists():
            return records

        with open(self.genotype_path, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                seqid = row.get('seqid', '')
                if seqid:
                    records[seqid] = row
        return records

    def check_results_exist(self) -> bool:
        """Check if MGEfinder results exist for this sample."""
        return self.clusterseq_path.exists()

    def has_mge_detected(self) -> bool:
        """Check if any MGEs were detected (not just empty results)."""
        no_mge_file = self.sample_dir / "03.results" / "genome" / "NO_MGE_DETECTED.txt"
        if no_mge_file.exists():
            return False
        return self.check_results_exist()

    def _detect_strand(self, is_sequence: str, contig_seq: str,
                       start: int, end: int) -> str:
        """
        Detect whether the IS element is on the forward or reverse strand
        of the assembly contig.

        Args:
            is_sequence: The inferred IS element sequence
            contig_seq: The full assembly contig sequence
            start: Start position (0-based)
            end: End position (0-based, exclusive)

        Returns:
            '+' for forward strand, '-' for reverse strand, '?' if unknown
        """
        if not is_sequence or not contig_seq or start >= end:
            return '?'

        sub = contig_seq[start:end].upper()
        is_upper = is_sequence.upper()

        if sub == is_upper:
            return '+'
        if reverse_complement(sub) == is_upper:
            return '-'
        return '?'

    def _extract_flanking_from_assembly(self, contig_seq: str, start: int,
                                         end: int, strand: str) -> dict:
        """
        Extract flanking regions from the assembly contig.

        The returned upstream/downstream are always in the IS element's
        reading frame:
          upstream_flank → IS_element → downstream_flank

        For forward strand (+):
          upstream  = contig[start-flank : start]
          downstream = contig[end : end+flank]

        For reverse strand (-):
          upstream  = revcomp(contig[end : end+flank])
          downstream = revcomp(contig[start-flank : start])

        Args:
            contig_seq: Full assembly contig sequence
            start: IS start on contig (0-based)
            end: IS end on contig (0-based, exclusive)
            strand: '+' or '-'

        Returns:
            Dict with 'upstream' and 'downstream' flanking info
        """
        contig_len = len(contig_seq)
        flank = self.flank_size

        # Extract raw flanking from contig coordinates
        left_start = max(0, start - flank)
        left_seq = contig_seq[left_start:start]
        right_seq = contig_seq[end:min(contig_len, end + flank)]

        if strand == '-':
            # IS is on reverse strand: flip and revcomp
            upstream_seq = reverse_complement(right_seq)
            upstream_contig_start = end + 1                    # 1-based
            upstream_contig_end = min(contig_len, end + flank)

            downstream_seq = reverse_complement(left_seq)
            downstream_contig_start = left_start + 1           # 1-based
            downstream_contig_end = start
        else:
            # Forward strand (or unknown): use as-is
            upstream_seq = left_seq
            upstream_contig_start = left_start + 1             # 1-based
            upstream_contig_end = start

            downstream_seq = right_seq
            downstream_contig_start = end + 1                  # 1-based
            downstream_contig_end = min(contig_len, end + flank)

        return {
            'upstream': {
                'sequence': upstream_seq,
                'length': len(upstream_seq),
                'contig_start': upstream_contig_start,
                'contig_end': upstream_contig_end,
            },
            'downstream': {
                'sequence': downstream_seq,
                'length': len(downstream_seq),
                'contig_start': downstream_contig_start,
                'contig_end': downstream_contig_end,
            }
        }

    def extract(self) -> List[Dict]:
        """
        Extract all IS elements with their sequences and flanking regions.

        Flanking regions are extracted from the assembly contig where the IS
        element sits, in the IS element's reading frame. This ensures:
            upstream_flank → IS_element → downstream_flank
        are contiguous and on the same contig.

        Returns:
            List of dictionaries, one per IS element.
        """
        if not self.check_results_exist():
            return []

        if not self.has_mge_detected():
            return []

        # Load data
        assembly_seqs = self._load_assembly()
        clusterseq_records = self._parse_clusterseq()
        genotype_records = self._parse_genotype()

        results = []

        for record in clusterseq_records:
            seqid = record.get('seqid', '')
            cluster = record.get('cluster', '')
            group = record.get('group', '')

            # Get IS element sequence from inferred_seq column
            is_sequence = record.get('inferred_seq', '')

            # Parse location from 'loc' column (e.g., "SAMN30697274.contig00004:127020-249243")
            loc = record.get('loc', '')
            is_contig = ''
            is_start = 0
            is_end = 0

            if ':' in loc and '-' in loc:
                contig_part, coords = loc.rsplit(':', 1)
                is_contig = contig_part
                start_end = coords.split('-')
                if len(start_end) == 2:
                    is_start = int(start_end[0])
                    is_end = int(start_end[1])

            # Get insertion site info from genotype
            genotype = genotype_records.get(seqid, {})
            ref_contig = genotype.get('contig', '')
            pos_5p = int(genotype.get('pos_5p', 0)) if genotype.get('pos_5p') else 0
            pos_3p = int(genotype.get('pos_3p', 0)) if genotype.get('pos_3p') else 0

            # Detect strand and extract flanking from assembly contig
            strand = '?'
            flanking_upstream = {
                'sequence': '', 'length': 0,
                'contig_start': 0, 'contig_end': 0,
            }
            flanking_downstream = {
                'sequence': '', 'length': 0,
                'contig_start': 0, 'contig_end': 0,
            }

            if is_contig and is_contig in assembly_seqs and is_start < is_end:
                contig_seq = assembly_seqs[is_contig]
                strand = self._detect_strand(is_sequence, contig_seq, is_start, is_end)
                flanking = self._extract_flanking_from_assembly(
                    contig_seq, is_start, is_end, strand)
                flanking_upstream = flanking['upstream']
                flanking_downstream = flanking['downstream']

            # Build IS element dict
            is_element = {
                'sequence': is_sequence,
                'length': len(is_sequence),
                'contig': is_contig,
                'start': is_start,
                'end': is_end,
                'strand': strand,
                'method': record.get('method', ''),
                'source_file': str(self.clusterseq_path),
            }

            # Build final record
            is_record = {
                'is_id': f"{self.sample_id}_{seqid}",
                'sample': self.sample_id,
                'seqid': seqid,
                'cluster': cluster,
                'group': group,
                'pair_id': record.get('pair_id', ''),
                'confidence': genotype.get('conf', ''),
                'is_element': is_element,
                'flanking_upstream': {
                    **flanking_upstream,
                    'source_file': str(self.assembly_path),
                    'contig': is_contig,
                },
                'flanking_downstream': {
                    **flanking_downstream,
                    'source_file': str(self.assembly_path),
                    'contig': is_contig,
                },
                'insertion_site': {
                    'contig': ref_contig,
                    'pos_5p': pos_5p,
                    'pos_3p': pos_3p,
                    'source_file': str(self.genotype_path)
                }
            }

            results.append(is_record)

        return results


def extract_is_elements(sample_dir: str, flank_size: int = 80) -> List[Dict]:
    """
    Convenience function to extract IS elements from a sample directory.

    Args:
        sample_dir: Path to MGEfinder sample directory
        flank_size: Size of flanking regions (default 80bp)

    Returns:
        List of IS element dictionaries
    """
    extractor = ISExtractor(sample_dir, flank_size)
    return extractor.extract()


def extract_batch(samples_dir: str, sample_ids: Optional[List[str]] = None,
                  flank_size: int = 80) -> Dict[str, List[Dict]]:
    """
    Extract IS elements from multiple samples.

    Args:
        samples_dir: Base directory containing sample subdirectories
        sample_ids: List of sample IDs to process (None = all)
        flank_size: Size of flanking regions (default 80bp)

    Returns:
        Dictionary mapping sample_id -> list of IS element records
    """
    samples_path = Path(samples_dir)
    results = {}

    if sample_ids is None:
        # Find all sample directories with results
        sample_ids = [
            d.name for d in samples_path.iterdir()
            if d.is_dir() and (d / "03.results").exists()
        ]

    for sample_id in sample_ids:
        sample_dir = samples_path / sample_id
        if sample_dir.exists():
            extractor = ISExtractor(str(sample_dir), flank_size)
            if extractor.has_mge_detected():
                results[sample_id] = extractor.extract()

    return results


if __name__ == "__main__":
    # Example usage
    import json
    import sys

    if len(sys.argv) < 2:
        print("Usage: python extractor.py <sample_dir> [flank_size]")
        sys.exit(1)

    sample_dir = sys.argv[1]
    flank_size = int(sys.argv[2]) if len(sys.argv) > 2 else 80

    results = extract_is_elements(sample_dir, flank_size)
    print(json.dumps(results, indent=2))
