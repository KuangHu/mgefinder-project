"""
Master List Builder - Merge matched samples with metadata from ENA files.

Creates a unified TSV containing Sample_ID, Scientific_Name, Ref_FTP, and Read_FTP
for all samples that have both paired-end reads and assembly data.
"""

import csv
import sys
from pathlib import Path
from typing import Optional

# Increase CSV field size limit for long URLs
csv.field_size_limit(sys.maxsize)


def build_master_list(
    matched_samples_file: str,
    reads_file: str,
    assembly_file: str,
    output_file: str = "master_list.tsv"
) -> int:
    """
    Build a master list merging matched samples with their metadata.

    Args:
        matched_samples_file: Path to file containing matched sample IDs
        reads_file: Path to ENA reads TSV file
        assembly_file: Path to ENA assembly TSV file
        output_file: Path for output master list TSV

    Returns:
        Number of complete pairs written to output
    """
    # Load matched sample IDs
    print("1. Loading matched sample IDs...")
    valid_samples = set()
    with open(matched_samples_file, 'r') as f:
        for line in f:
            valid_samples.add(line.strip())
    print(f"   Loaded {len(valid_samples)} valid samples.")

    # Store merged data: {sample_id: {'ref_url': '', 'read_url': '', 'name': ''}}
    data_map = {}

    # Process assembly file (extract reference URL)
    print("2. Parsing Assembly Metadata...")
    with open(assembly_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            sid = row['sample_accession']
            if sid in valid_samples:
                # Only take first occurrence to avoid duplicates
                if sid not in data_map:
                    data_map[sid] = {
                        'name': row['scientific_name'],
                        'ref_url': row['submitted_ftp'],
                        'read_url': None
                    }

    # Process reads file (extract paired-end fastq URLs)
    print("3. Parsing Reads Metadata...")
    skipped_non_paired = 0
    with open(reads_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            sid = row['sample_accession']
            if sid in data_map and data_map[sid]['read_url'] is None:
                ftp = row['fastq_ftp']
                if ftp and ftp.strip():
                    urls = ftp.split(';')
                    # Only accept paired-end reads (must have both _1.fastq.gz and _2.fastq.gz)
                    if len(urls) >= 2 and '_1.fastq.gz' in ftp and '_2.fastq.gz' in ftp:
                        r1 = [u for u in urls if '_1.fastq.gz' in u]
                        r2 = [u for u in urls if '_2.fastq.gz' in u]
                        if r1 and r2:
                            data_map[sid]['read_url'] = f"{r1[0]};{r2[0]}"
                        else:
                            skipped_non_paired += 1
                    else:
                        skipped_non_paired += 1
    print(f"   Skipped {skipped_non_paired} non-paired-end samples.")

    # Write output
    print(f"4. Writing to {output_file}...")
    count = 0
    with open(output_file, 'w') as f:
        f.write("Sample_ID\tScientific_Name\tRef_FTP\tRead_FTP\n")
        for sid, info in data_map.items():
            if info['ref_url'] and info['read_url']:
                f.write(f"{sid}\t{info['name']}\t{info['ref_url']}\t{info['read_url']}\n")
                count += 1

    print(f"Done! Successfully merged {count} complete pairs into {output_file}")
    return count


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Build master list from ENA metadata")
    parser.add_argument("--matched", required=True, help="Matched samples file")
    parser.add_argument("--reads", required=True, help="ENA reads TSV file")
    parser.add_argument("--assembly", required=True, help="ENA assembly TSV file")
    parser.add_argument("--output", default="master_list.tsv", help="Output file")
    args = parser.parse_args()

    build_master_list(args.matched, args.reads, args.assembly, args.output)
