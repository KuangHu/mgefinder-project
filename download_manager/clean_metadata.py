#!/usr/bin/env python3
"""
Clean ENA metadata file - keep only useful, standardized columns.

Keeps: sample_accession, scientific_name, host, temperature, ph, altitude
Drops: country, collection_date, salinity, depth, environment_*, isolation_source
"""

import csv
import sys
from pathlib import Path

csv.field_size_limit(sys.maxsize)

KEEP_COLUMNS = [
    'sample_accession',
    'scientific_name',
    'host',
    'temperature',
    'ph',
    'altitude'
]

MISSING_VALUES = {
    '', 'missing', 'not applicable', 'not collected', 'not provided',
    'unknown', 'n/a', 'na', 'none', '-', 'null', 'unspecified',
    'not determined', 'not available', 'missing: third party data'
}


def clean_metadata(input_file: str, output_file: str):
    """
    Clean metadata file, keeping only useful columns.

    Args:
        input_file: Path to full metadata TSV
        output_file: Path for cleaned output TSV
    """
    print(f"Input:  {input_file}")
    print(f"Output: {output_file}")
    print(f"Keeping columns: {KEEP_COLUMNS}")
    print("-" * 60)

    total_rows = 0

    with open(input_file, 'r') as fin, open(output_file, 'w', newline='') as fout:
        reader = csv.DictReader(fin, delimiter='\t')
        writer = csv.DictWriter(fout, fieldnames=KEEP_COLUMNS, delimiter='\t')
        writer.writeheader()

        for row in reader:
            total_rows += 1

            # Extract and clean values
            clean_row = {}
            for col in KEEP_COLUMNS:
                val = row.get(col, '').strip()
                # Normalize missing values to empty string
                if val.lower() in MISSING_VALUES:
                    val = ''
                clean_row[col] = val

            writer.writerow(clean_row)

            if total_rows % 1000000 == 0:
                print(f"  Processed {total_rows:,} rows...")

    print("-" * 60)
    print(f"Done! Wrote {total_rows:,} rows to {output_file}")

    # Show file sizes
    input_size = Path(input_file).stat().st_size / (1024 * 1024)
    output_size = Path(output_file).stat().st_size / (1024 * 1024)
    reduction = (1 - output_size / input_size) * 100

    print(f"Size: {input_size:.1f} MB -> {output_size:.1f} MB ({reduction:.1f}% reduction)")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Clean ENA metadata file")
    parser.add_argument("--input", required=True, help="Input metadata TSV")
    parser.add_argument("--output", required=True, help="Output cleaned TSV")
    args = parser.parse_args()

    clean_metadata(args.input, args.output)
