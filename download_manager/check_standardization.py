#!/usr/bin/env python3
"""
Check standardization of metadata fields.
Shows unique value counts and patterns for each field.
"""

import csv
import sys
import re
from collections import Counter

csv.field_size_limit(sys.maxsize)

MISSING_VALUES = {
    '', 'missing', 'not applicable', 'not collected', 'not provided',
    'unknown', 'n/a', 'na', 'none', '-', 'null', 'unspecified',
    'not determined', 'not available'
}

def check_field_standardization(metadata_file: str, fields_to_check: list, sample_limit: int = None):
    """Check how standardized each field is."""

    field_values = {f: Counter() for f in fields_to_check}
    total_rows = 0

    print(f"Scanning metadata file...")

    with open(metadata_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')

        for row in reader:
            total_rows += 1

            for field in fields_to_check:
                val = row.get(field, '').strip()
                val_lower = val.lower()

                if val_lower not in MISSING_VALUES:
                    field_values[field][val] += 1

            if sample_limit and total_rows >= sample_limit:
                break

            if total_rows % 1000000 == 0:
                print(f"  Processed {total_rows:,} rows...")

    print(f"\nTotal rows scanned: {total_rows:,}\n")

    # Analyze each field
    for field in fields_to_check:
        values = field_values[field]
        unique_count = len(values)
        total_filled = sum(values.values())

        print("=" * 70)
        print(f"FIELD: {field}")
        print("=" * 70)
        print(f"Unique values: {unique_count:,}")
        print(f"Filled entries: {total_filled:,}")

        # Top values
        print(f"\nTop 15 values:")
        for val, count in values.most_common(15):
            pct = count / total_filled * 100 if total_filled > 0 else 0
            display_val = val[:60] + "..." if len(val) > 60 else val
            print(f"  {count:>8,} ({pct:>5.1f}%)  {display_val}")

        # Check for patterns
        print(f"\nPattern analysis:")

        # Numeric check
        numeric_count = sum(c for v, c in values.items() if re.match(r'^-?\d+\.?\d*$', v))
        if numeric_count > 0:
            print(f"  - Pure numeric values: {numeric_count:,} ({numeric_count/total_filled*100:.1f}%)")

        # Values with units
        with_units = sum(c for v, c in values.items() if re.search(r'\d+\s*[a-zA-Z°%]+', v))
        if with_units > 0:
            print(f"  - Values with units: {with_units:,} ({with_units/total_filled*100:.1f}%)")

        # Range values (e.g., "10-20" or "10 to 20")
        ranges = sum(c for v, c in values.items() if re.search(r'\d+\s*(-|to)\s*\d+', v.lower()))
        if ranges > 0:
            print(f"  - Range values: {ranges:,} ({ranges/total_filled*100:.1f}%)")

        # Free text (long strings)
        free_text = sum(c for v, c in values.items() if len(v) > 50)
        if free_text > 0:
            print(f"  - Long text (>50 chars): {free_text:,} ({free_text/total_filled*100:.1f}%)")

        # Standardization score
        top_10_coverage = sum(c for _, c in values.most_common(10)) / total_filled * 100 if total_filled > 0 else 0
        top_100_coverage = sum(c for _, c in values.most_common(100)) / total_filled * 100 if total_filled > 0 else 0

        print(f"\nStandardization indicators:")
        print(f"  - Top 10 values cover: {top_10_coverage:.1f}%")
        print(f"  - Top 100 values cover: {top_100_coverage:.1f}%")

        if unique_count < 100:
            print(f"  => HIGHLY STANDARDIZED (only {unique_count} unique values)")
        elif top_100_coverage > 90:
            print(f"  => MOSTLY STANDARDIZED (top 100 cover {top_100_coverage:.1f}%)")
        elif top_100_coverage > 50:
            print(f"  => PARTIALLY STANDARDIZED")
        else:
            print(f"  => FREE TEXT / NOT STANDARDIZED")

        print()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--fields", nargs="+", required=True)
    parser.add_argument("--limit", type=int, help="Row limit for quick check")
    args = parser.parse_args()

    check_field_standardization(args.metadata, args.fields, args.limit)
