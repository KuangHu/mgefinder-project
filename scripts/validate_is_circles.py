#!/usr/bin/env python3
"""
Validate Circle-Map results against IS element lengths.
True circularized IS: circle coordinates ≈ full IS length

Usage:
    python validate_is_circles.py \
        --circles is_circles.bed \
        --reference is_reference.fna \
        --output circularized_is.tsv
"""

import argparse
import csv
import os


def parse_fasta_lengths(fasta_file):
    """Parse FASTA file and return dict of {seq_id: length}."""
    lengths = {}
    current_id = None
    current_len = 0

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    lengths[current_id] = current_len
                current_id = line[1:].split()[0]
                current_len = 0
            else:
                current_len += len(line)

        if current_id:
            lengths[current_id] = current_len

    return lengths


def validate_circles(circle_bed, is_reference, output_tsv,
                     length_tolerance=0.1, min_split_reads=2,
                     start_tolerance=0.05, end_tolerance=0.05):
    """
    Filter circles that match IS element lengths.

    Criteria for TRUE circularized IS:
    1. Circle length within ±tolerance of IS reference length
    2. Circle starts near position 0 (head)
    3. Circle ends near IS length (tail)
    4. At least min_split_reads supporting junction

    Args:
        circle_bed: Circle-Map BED output file
        is_reference: IS reference FASTA file
        output_tsv: Output TSV file for validated circles
        length_tolerance: Allowed deviation from IS length (0.1 = 10%)
        min_split_reads: Minimum split reads required
        start_tolerance: Max start position as fraction of IS length
        end_tolerance: Min end position as fraction of IS length (from end)

    Returns:
        List of validated circle dictionaries
    """

    # Load IS reference lengths
    is_lengths = parse_fasta_lengths(is_reference)
    print(f"Loaded {len(is_lengths)} IS sequences from reference")

    if not os.path.exists(circle_bed) or os.path.getsize(circle_bed) == 0:
        print("No circles detected by Circle-Map")
        return []

    # Parse circles and validate
    all_circles = []
    valid_circles = []

    with open(circle_bed, 'r') as f:
        for line in f:
            if not line.strip():
                continue

            parts = line.strip().split('\t')
            is_id = parts[0]
            circle_start = int(parts[1])
            circle_end = int(parts[2])
            split_reads = int(parts[3]) if len(parts) > 3 else 0
            discordant = int(parts[4]) if len(parts) > 4 else 0
            score = float(parts[5]) if len(parts) > 5 else 0

            circle_length = circle_end - circle_start

            all_circles.append({
                'is_id': is_id,
                'circle_start': circle_start,
                'circle_end': circle_end,
                'circle_length': circle_length,
                'split_reads': split_reads,
                'discordant_reads': discordant,
                'score': score
            })

            if is_id not in is_lengths:
                continue

            is_length = is_lengths[is_id]
            length_ratio = circle_length / is_length if is_length > 0 else 0

            # Validation criteria
            is_full_length = (1 - length_tolerance) <= length_ratio <= (1 + length_tolerance)
            starts_at_head = circle_start <= is_length * start_tolerance
            ends_at_tail = circle_end >= is_length * (1 - end_tolerance)
            has_support = split_reads >= min_split_reads

            # Determine validation status
            if is_full_length and starts_at_head and ends_at_tail and has_support:
                status = 'CONFIRMED_CIRCULAR_IS'
            elif is_full_length and has_support:
                status = 'LIKELY_CIRCULAR_IS'
            elif has_support:
                status = 'PARTIAL_CIRCLE'
            else:
                status = 'LOW_CONFIDENCE'
                continue  # Skip low confidence

            valid_circles.append({
                'is_id': is_id,
                'is_length': is_length,
                'circle_start': circle_start,
                'circle_end': circle_end,
                'circle_length': circle_length,
                'length_ratio': round(length_ratio, 3),
                'split_reads': split_reads,
                'discordant_reads': discordant,
                'score': score,
                'starts_at_head': starts_at_head,
                'ends_at_tail': ends_at_tail,
                'validation': status
            })

    print(f"Total circles from Circle-Map: {len(all_circles)}")
    print(f"Validated circles: {len(valid_circles)}")

    # Count by status
    status_counts = {}
    for vc in valid_circles:
        status = vc['validation']
        status_counts[status] = status_counts.get(status, 0) + 1

    for status, count in sorted(status_counts.items()):
        print(f"  {status}: {count}")

    # Write output
    if valid_circles:
        fieldnames = ['is_id', 'is_length', 'circle_start', 'circle_end',
                      'circle_length', 'length_ratio', 'split_reads',
                      'discordant_reads', 'score', 'validation']

        with open(output_tsv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t',
                                    extrasaction='ignore')
            writer.writeheader()
            writer.writerows(valid_circles)

        print(f"\nResults written to: {output_tsv}")
    else:
        # Write empty file with header
        with open(output_tsv, 'w') as f:
            f.write("is_id\tis_length\tcircle_start\tcircle_end\tcircle_length\t"
                    "length_ratio\tsplit_reads\tdiscordant_reads\tscore\tvalidation\n")
        print(f"\nNo validated circles found. Empty file written to: {output_tsv}")

    return valid_circles


def main():
    parser = argparse.ArgumentParser(
        description='Validate Circle-Map results for circularized IS elements')
    parser.add_argument('--circles', required=True,
                        help='Circle-Map BED output file')
    parser.add_argument('--reference', required=True,
                        help='IS reference FASTA file')
    parser.add_argument('--output', required=True,
                        help='Output TSV file')
    parser.add_argument('--tolerance', type=float, default=0.1,
                        help='Length tolerance (default: 0.1 = 10%%)')
    parser.add_argument('--min-reads', type=int, default=2,
                        help='Minimum split reads (default: 2)')

    args = parser.parse_args()

    valid = validate_circles(
        args.circles,
        args.reference,
        args.output,
        length_tolerance=args.tolerance,
        min_split_reads=args.min_reads
    )

    # Exit with count of confirmed circles
    confirmed = sum(1 for v in valid if v['validation'] == 'CONFIRMED_CIRCULAR_IS')
    print(f"\n=== CONFIRMED CIRCULAR IS ELEMENTS: {confirmed} ===")


if __name__ == '__main__':
    main()
