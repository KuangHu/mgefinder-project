#!/usr/bin/env python3
"""
Extract IS element sequences to create a reference for Circle-Map analysis.
Filters by length and removes duplicates.

Usage:
    python extract_is_reference.py --sample-dir samples/SAMN00002238 --output is_reference.fna
"""

import argparse
import os
import hashlib


def extract_is_reference(sample_dir, output_fasta, min_length=500, max_length=5000):
    """
    Extract IS sequences from MGEfinder output for use as Circle-Map reference.

    Args:
        sample_dir: Path to sample directory
        output_fasta: Output FASTA file path
        min_length: Minimum IS length to include
        max_length: Maximum IS length to include

    Returns:
        Number of sequences extracted
    """
    sample_id = os.path.basename(sample_dir.rstrip('/'))

    # Input: MGEfinder IS sequences
    is_fasta = f"{sample_dir}/03.results/genome/04.makefasta.genome.all_seqs.fna"

    if not os.path.exists(is_fasta):
        print(f"No IS sequences found for {sample_id}")
        return 0

    # Check if file is empty or contains NO_MGE marker
    no_mge_file = f"{sample_dir}/03.results/genome/NO_MGE_DETECTED.txt"
    if os.path.exists(no_mge_file):
        print(f"No MGEs detected for {sample_id}")
        return 0

    # Parse FASTA manually (avoid BioPython dependency)
    sequences = []
    seen_hashes = set()
    current_id = None
    current_seq = []

    with open(is_fasta, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence
                if current_id and current_seq:
                    seq = ''.join(current_seq)
                    seq_len = len(seq)
                    seq_hash = hashlib.md5(seq.encode()).hexdigest()

                    if min_length <= seq_len <= max_length and seq_hash not in seen_hashes:
                        sequences.append({
                            'id': f"{sample_id}_{current_id}",
                            'seq': seq,
                            'length': seq_len
                        })
                        seen_hashes.add(seq_hash)

                # Start new sequence
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

        # Don't forget last sequence
        if current_id and current_seq:
            seq = ''.join(current_seq)
            seq_len = len(seq)
            seq_hash = hashlib.md5(seq.encode()).hexdigest()

            if min_length <= seq_len <= max_length and seq_hash not in seen_hashes:
                sequences.append({
                    'id': f"{sample_id}_{current_id}",
                    'seq': seq,
                    'length': seq_len
                })

    # Write output
    if sequences:
        os.makedirs(os.path.dirname(output_fasta) or '.', exist_ok=True)
        with open(output_fasta, 'w') as f:
            for seq_data in sequences:
                f.write(f">{seq_data['id']} length={seq_data['length']}\n")
                # Write sequence in 80-char lines
                seq = seq_data['seq']
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + '\n')

        print(f"Extracted {len(sequences)} unique IS sequences for {sample_id}")
        print(f"  Length range: {min(s['length'] for s in sequences)} - {max(s['length'] for s in sequences)} bp")
    else:
        print(f"No IS sequences passed filters for {sample_id}")

    return len(sequences)


def main():
    parser = argparse.ArgumentParser(
        description='Extract IS element sequences for Circle-Map analysis')
    parser.add_argument('--sample-dir', required=True,
                        help='Path to sample directory')
    parser.add_argument('--output', required=True,
                        help='Output FASTA file path')
    parser.add_argument('--min-length', type=int, default=500,
                        help='Minimum IS length (default: 500)')
    parser.add_argument('--max-length', type=int, default=5000,
                        help='Maximum IS length (default: 5000)')

    args = parser.parse_args()

    count = extract_is_reference(
        args.sample_dir,
        args.output,
        args.min_length,
        args.max_length
    )

    exit(0 if count > 0 else 1)


if __name__ == '__main__':
    main()
