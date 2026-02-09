"""
IS Element Extractor

Extract IS element DNA sequences and flanking regions from MGEfinder output.
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


def extract_flanking(genome_seq: str, pos_5p: int, pos_3p: int,
                     flank_size: int = 80) -> Dict[str, dict]:
    """
    Extract flanking regions from genome sequence.

    Args:
        genome_seq: Full genome sequence string
        pos_5p: 5' position (1-based)
        pos_3p: 3' position (1-based)
        flank_size: Size of flanking region to extract (default 80bp)

    Returns:
        Dictionary with upstream and downstream flanking info
    """
    genome_len = len(genome_seq)

    # Determine insertion site boundaries (handle both orientations)
    left_pos = min(pos_5p, pos_3p)
    right_pos = max(pos_5p, pos_3p)

    # Calculate flanking region coordinates (convert to 0-based)
    upstream_start = max(0, left_pos - 1 - flank_size)
    upstream_end = left_pos - 1
    downstream_start = right_pos  # already 0-based after -1+1
    downstream_end = min(genome_len, right_pos + flank_size)

    return {
        'upstream': {
            'sequence': genome_seq[upstream_start:upstream_end],
            'start': upstream_start + 1,  # back to 1-based
            'end': upstream_end,
            'length': upstream_end - upstream_start
        },
        'downstream': {
            'sequence': genome_seq[downstream_start:downstream_end],
            'start': downstream_start + 1,  # back to 1-based
            'end': downstream_end,
            'length': downstream_end - downstream_start
        }
    }


class ISExtractor:
    """Extract IS elements and flanking regions from MGEfinder output."""

    def __init__(self, sample_dir: str, flank_size: int = 80):
        """
        Initialize extractor for a sample directory.

        Args:
            sample_dir: Path to MGEfinder sample directory
                       (containing 00.genome, 03.results, etc.)
            flank_size: Size of flanking regions to extract (default 80bp)
        """
        self.sample_dir = Path(sample_dir)
        self.flank_size = flank_size
        self.sample_id = self.sample_dir.name

        # Define file paths
        self.genome_path = self.sample_dir / "00.genome" / "genome.fna"
        self.clusterseq_path = self.sample_dir / "03.results" / "genome" / "01.clusterseq.genome.tsv"
        self.genotype_path = self.sample_dir / "03.results" / "genome" / "02.genotype.genome.tsv"
        self.fasta_path = self.sample_dir / "03.results" / "genome" / "04.makefasta.genome.repr_seqs.fna"

        # Cache for loaded data
        self._genome_seqs = None
        self._is_sequences = None

    def _load_genome(self) -> Dict[str, str]:
        """Load reference genome sequences."""
        if self._genome_seqs is None:
            if not self.genome_path.exists():
                raise FileNotFoundError(f"Genome file not found: {self.genome_path}")
            self._genome_seqs = parse_fasta(str(self.genome_path))
        return self._genome_seqs

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

    def extract(self) -> List[Dict]:
        """
        Extract all IS elements with their sequences and flanking regions.

        Returns:
            List of dictionaries, one per IS element, containing:
            - is_id: Unique identifier for this IS element
            - sample: Sample ID
            - seqid: Sequence ID from MGEfinder
            - cluster: Cluster ID
            - group: Group ID
            - is_element: Dict with sequence and location info
            - flanking_upstream: Dict with upstream flanking sequence and location
            - flanking_downstream: Dict with downstream flanking sequence and location
        """
        if not self.check_results_exist():
            return []

        if not self.has_mge_detected():
            return []

        # Load data
        genome_seqs = self._load_genome()
        is_fasta_seqs = self._load_is_sequences()
        clusterseq_records = self._parse_clusterseq()
        genotype_records = self._parse_genotype()

        results = []

        for record in clusterseq_records:
            seqid = record.get('seqid', '')
            cluster = record.get('cluster', '')
            group = record.get('group', '')

            # Get IS element sequence
            # Try from clusterseq first (inferred_seq column)
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

            # Get insertion site info from genotype for flanking regions
            genotype = genotype_records.get(seqid, {})
            ref_contig = genotype.get('contig', '')
            pos_5p = int(genotype.get('pos_5p', 0)) if genotype.get('pos_5p') else 0
            pos_3p = int(genotype.get('pos_3p', 0)) if genotype.get('pos_3p') else 0

            # Build IS element dict
            is_element = {
                'sequence': is_sequence,
                'length': len(is_sequence),
                'source_file': str(self.clusterseq_path),
                'contig': is_contig,
                'start': is_start,
                'end': is_end,
                'method': record.get('method', '')
            }

            # Extract flanking regions from reference genome
            flanking_upstream = {
                'sequence': '',
                'start': 0,
                'end': 0,
                'length': 0,
                'source_file': '',
                'contig': ''
            }
            flanking_downstream = {
                'sequence': '',
                'start': 0,
                'end': 0,
                'length': 0,
                'source_file': '',
                'contig': ''
            }

            if ref_contig and ref_contig in genome_seqs and pos_5p > 0 and pos_3p > 0:
                genome_seq = genome_seqs[ref_contig]
                flanking = extract_flanking(genome_seq, pos_5p, pos_3p, self.flank_size)

                flanking_upstream = {
                    'sequence': flanking['upstream']['sequence'],
                    'start': flanking['upstream']['start'],
                    'end': flanking['upstream']['end'],
                    'length': flanking['upstream']['length'],
                    'source_file': str(self.genome_path),
                    'contig': ref_contig
                }
                flanking_downstream = {
                    'sequence': flanking['downstream']['sequence'],
                    'start': flanking['downstream']['start'],
                    'end': flanking['downstream']['end'],
                    'length': flanking['downstream']['length'],
                    'source_file': str(self.genome_path),
                    'contig': ref_contig
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
                'flanking_upstream': flanking_upstream,
                'flanking_downstream': flanking_downstream,
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
