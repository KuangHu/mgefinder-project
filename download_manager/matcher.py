"""
Sample Matcher - Find samples that have both reads and assembly data.

This module identifies "golden pairs" - samples that exist in both the reads
and assembly datasets from ENA, making them suitable for MGEfinder analysis.
"""

import subprocess
from pathlib import Path


def match_samples(
    reads_file: str,
    assembly_file: str,
    output_file: str = "matched_samples.txt",
    sample_col: int = 2
) -> int:
    """
    Find samples that have both reads and assembly data.

    Args:
        reads_file: Path to ENA reads TSV file
        assembly_file: Path to ENA assembly TSV file
        output_file: Path for matched samples output
        sample_col: Column number containing sample IDs (1-indexed, default=2)

    Returns:
        Number of matched samples found
    """
    reads_path = Path(reads_file)
    assembly_path = Path(assembly_file)
    output_path = Path(output_file)

    if not reads_path.exists():
        raise FileNotFoundError(f"Reads file not found: {reads_file}")
    if not assembly_path.exists():
        raise FileNotFoundError(f"Assembly file not found: {assembly_file}")

    # Extract and sort sample IDs from reads file
    print("1. Extracting sample IDs from reads file...")
    sorted_reads = output_path.parent / "sorted_reads_samples.txt"
    cmd_reads = f"cut -f{sample_col} '{reads_file}' | tail -n +2 | sort > '{sorted_reads}'"
    subprocess.run(cmd_reads, shell=True, check=True)

    # Extract and sort sample IDs from assembly file
    print("2. Extracting sample IDs from assembly file...")
    sorted_assembly = output_path.parent / "sorted_assembly_samples.txt"
    cmd_assembly = f"cut -f{sample_col} '{assembly_file}' | tail -n +2 | sort > '{sorted_assembly}'"
    subprocess.run(cmd_assembly, shell=True, check=True)

    # Find intersection
    print("3. Finding intersection (samples with both reads and assembly)...")
    cmd_match = f"comm -12 '{sorted_reads}' '{sorted_assembly}' > '{output_file}'"
    subprocess.run(cmd_match, shell=True, check=True)

    # Count results
    with open(output_file, 'r') as f:
        count = sum(1 for _ in f)

    print(f"Found {count} matched samples.")
    return count


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Find samples with both reads and assembly")
    parser.add_argument("--reads", required=True, help="ENA reads TSV file")
    parser.add_argument("--assembly", required=True, help="ENA assembly TSV file")
    parser.add_argument("--output", default="matched_samples.txt", help="Output file")
    args = parser.parse_args()

    match_samples(args.reads, args.assembly, args.output)
