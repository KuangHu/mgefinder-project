#!/usr/bin/env python3
"""
Analyze ENA sample metadata file.

Reports:
- Total sample count
- Field completeness (% non-empty, excluding placeholder values)
- Coverage of master_list samples
- Common placeholder values found
"""

import csv
import sys
from pathlib import Path
from collections import Counter

csv.field_size_limit(sys.maxsize)

# Values to treat as "missing"
MISSING_VALUES = {
    '', 'missing', 'not applicable', 'not collected', 'not provided',
    'unknown', 'n/a', 'na', 'none', '-', 'null', 'unspecified',
    'not determined', 'not available'
}


def analyze_metadata(
    metadata_file: str,
    master_list_file: str = None,
    sample_size: int = None
):
    """
    Analyze metadata file completeness and coverage.

    Args:
        metadata_file: Path to metadata TSV
        master_list_file: Optional path to master list for coverage check
        sample_size: If set, only analyze first N rows (for quick preview)
    """
    print(f"Analyzing: {metadata_file}")
    print("=" * 60)

    # Load master list sample IDs if provided
    master_samples = None
    if master_list_file and Path(master_list_file).exists():
        print(f"Loading master list from {master_list_file}...")
        master_samples = set()
        with open(master_list_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                master_samples.add(row['Sample_ID'])
        print(f"  Master list contains {len(master_samples)} samples")

    # Initialize counters
    field_counts = {}  # {field: {'total': n, 'filled': n, 'values': Counter}}
    total_rows = 0
    matched_master = 0

    print(f"\nScanning metadata file...")

    with open(metadata_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        headers = reader.fieldnames

        # Initialize counters for each field
        for h in headers:
            field_counts[h] = {
                'total': 0,
                'filled': 0,
                'placeholders': Counter(),
                'sample_values': []
            }

        for row in reader:
            total_rows += 1

            # Check if in master list
            if master_samples and row.get('sample_accession') in master_samples:
                matched_master += 1

            # Count field completeness
            for field, value in row.items():
                if field not in field_counts:
                    continue

                field_counts[field]['total'] += 1
                value_lower = value.strip().lower() if value else ''

                if value_lower in MISSING_VALUES:
                    if value_lower:  # Track non-empty placeholders
                        field_counts[field]['placeholders'][value.strip()] += 1
                else:
                    field_counts[field]['filled'] += 1
                    # Keep some sample values
                    if len(field_counts[field]['sample_values']) < 5:
                        field_counts[field]['sample_values'].append(value.strip()[:50])

            if sample_size and total_rows >= sample_size:
                print(f"  (Stopped at {sample_size} rows for preview)")
                break

            if total_rows % 1000000 == 0:
                print(f"  Processed {total_rows:,} rows...")

    # Print results
    print(f"\n{'=' * 60}")
    print(f"SUMMARY")
    print(f"{'=' * 60}")
    print(f"Total samples: {total_rows:,}")

    if master_samples:
        coverage = (matched_master / len(master_samples)) * 100
        print(f"Master list coverage: {matched_master:,} / {len(master_samples):,} ({coverage:.1f}%)")

    print(f"\n{'=' * 60}")
    print(f"FIELD COMPLETENESS (excluding placeholder values)")
    print(f"{'=' * 60}")
    print(f"{'Field':<25} {'Filled':>12} {'Total':>12} {'%':>8}")
    print("-" * 60)

    for field in headers:
        fc = field_counts[field]
        pct = (fc['filled'] / fc['total'] * 100) if fc['total'] > 0 else 0
        print(f"{field:<25} {fc['filled']:>12,} {fc['total']:>12,} {pct:>7.1f}%")

    print(f"\n{'=' * 60}")
    print(f"SAMPLE VALUES PER FIELD")
    print(f"{'=' * 60}")
    for field in headers:
        samples = field_counts[field]['sample_values']
        if samples:
            print(f"\n{field}:")
            for s in samples[:3]:
                print(f"  - {s}")

    print(f"\n{'=' * 60}")
    print(f"COMMON PLACEHOLDER VALUES")
    print(f"{'=' * 60}")
    all_placeholders = Counter()
    for field in headers:
        all_placeholders.update(field_counts[field]['placeholders'])

    for val, count in all_placeholders.most_common(10):
        print(f"  '{val}': {count:,}")

    return {
        'total_rows': total_rows,
        'field_counts': field_counts,
        'matched_master': matched_master
    }


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Analyze ENA metadata file")
    parser.add_argument("--metadata", required=True, help="Metadata TSV file")
    parser.add_argument("--master", help="Master list file for coverage check")
    parser.add_argument("--preview", type=int, help="Only analyze first N rows")
    args = parser.parse_args()

    analyze_metadata(args.metadata, args.master, args.preview)
