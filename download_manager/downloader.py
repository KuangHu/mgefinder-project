"""
Sample Downloader - Download reference genomes and reads from ENA.

Creates MGEfinder-compatible directory structure:
    SampleID/
        00.reference/
            genome.fna
        00.reads/
            SampleID_R1.fastq.gz
            SampleID_R2.fastq.gz
"""

import csv
import subprocess
import sys
from pathlib import Path
from typing import Optional, Callable

csv.field_size_limit(sys.maxsize)


def download_samples(
    target_file: str,
    output_dir: str = ".",
    skip_existing: bool = True,
    progress_callback: Optional[Callable[[str, int, int], None]] = None
) -> dict:
    """
    Download samples from ENA to MGEfinder directory structure.

    Args:
        target_file: Path to target samples TSV
        output_dir: Base directory for downloads
        skip_existing: Skip samples that already have complete data
        progress_callback: Optional callback(sample_id, current, total)

    Returns:
        Dict with counts: {'success': n, 'failed': n, 'skipped': n}
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    stats = {'success': 0, 'failed': 0, 'skipped': 0}

    # Load samples
    samples = []
    with open(target_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            samples.append(row)

    total = len(samples)
    print(f"Starting download for {total} samples...")
    print("-" * 50)

    for idx, row in enumerate(samples, 1):
        sample_id = row['Sample_ID']
        name = row['Scientific_Name']
        ref_url = row['Ref_FTP']
        read_urls = row['Read_FTP']

        if progress_callback:
            progress_callback(sample_id, idx, total)

        print(f"[{idx}/{total}] Processing: {sample_id} ({name})")

        base_dir = output_path / sample_id
        ref_dir = base_dir / "00.reference"
        reads_dir = base_dir / "00.reads"

        # Check if already complete
        if skip_existing:
            genome_file = ref_dir / "genome.fna"
            r1_file = reads_dir / f"{sample_id}_R1.fastq.gz"
            r2_file = reads_dir / f"{sample_id}_R2.fastq.gz"

            if genome_file.exists() and r1_file.exists() and r2_file.exists():
                print(f"   Skipping (already exists)")
                stats['skipped'] += 1
                continue

        # Create directories
        ref_dir.mkdir(parents=True, exist_ok=True)
        reads_dir.mkdir(parents=True, exist_ok=True)

        try:
            # Download reference genome
            if not ref_url.startswith(('ftp://', 'http://')):
                ref_url = f"ftp://{ref_url}"

            print(f"   -> [Ref] Downloading genome...")
            genome_gz = ref_dir / "genome.fna.gz"
            subprocess.run(
                ["wget", "-q", "--show-progress", "-O", str(genome_gz), ref_url],
                check=True
            )

            # Decompress reference
            subprocess.run(["gunzip", "-f", str(genome_gz)], check=True)

            # Download reads
            r1_url, r2_url = read_urls.split(';')

            if not r1_url.startswith(('ftp://', 'http://')):
                r1_url = f"ftp://{r1_url}"
            if not r2_url.startswith(('ftp://', 'http://')):
                r2_url = f"ftp://{r2_url}"

            print(f"   -> [R1] Downloading Read 1...")
            r1_file = reads_dir / f"{sample_id}_R1.fastq.gz"
            subprocess.run(
                ["wget", "-q", "--show-progress", "-O", str(r1_file), r1_url],
                check=True
            )

            print(f"   -> [R2] Downloading Read 2...")
            r2_file = reads_dir / f"{sample_id}_R2.fastq.gz"
            subprocess.run(
                ["wget", "-q", "--show-progress", "-O", str(r2_file), r2_url],
                check=True
            )

            print(f"   Done: {sample_id}")
            stats['success'] += 1

        except subprocess.CalledProcessError as e:
            print(f"   FAILED: {e}")
            stats['failed'] += 1

        print("-" * 50)

    print(f"\nDownload complete!")
    print(f"  Success: {stats['success']}")
    print(f"  Failed:  {stats['failed']}")
    print(f"  Skipped: {stats['skipped']}")

    return stats


def download_single_sample(
    sample_id: str,
    ref_url: str,
    read_urls: str,
    output_dir: str = "."
) -> bool:
    """
    Download a single sample.

    Args:
        sample_id: Sample accession ID
        ref_url: FTP URL for reference genome
        read_urls: Semicolon-separated R1;R2 FTP URLs
        output_dir: Base directory for download

    Returns:
        True if successful, False otherwise
    """
    output_path = Path(output_dir)
    base_dir = output_path / sample_id
    ref_dir = base_dir / "00.reference"
    reads_dir = base_dir / "00.reads"

    ref_dir.mkdir(parents=True, exist_ok=True)
    reads_dir.mkdir(parents=True, exist_ok=True)

    try:
        # Normalize URLs
        if not ref_url.startswith(('ftp://', 'http://')):
            ref_url = f"ftp://{ref_url}"

        # Download and decompress reference
        genome_gz = ref_dir / "genome.fna.gz"
        subprocess.run(["wget", "-q", "-O", str(genome_gz), ref_url], check=True)
        subprocess.run(["gunzip", "-f", str(genome_gz)], check=True)

        # Download reads
        r1_url, r2_url = read_urls.split(';')
        if not r1_url.startswith(('ftp://', 'http://')):
            r1_url = f"ftp://{r1_url}"
        if not r2_url.startswith(('ftp://', 'http://')):
            r2_url = f"ftp://{r2_url}"

        subprocess.run(
            ["wget", "-q", "-O", str(reads_dir / f"{sample_id}_R1.fastq.gz"), r1_url],
            check=True
        )
        subprocess.run(
            ["wget", "-q", "-O", str(reads_dir / f"{sample_id}_R2.fastq.gz"), r2_url],
            check=True
        )

        return True

    except subprocess.CalledProcessError:
        return False


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Download samples from ENA")
    parser.add_argument("--input", required=True, help="Target samples TSV file")
    parser.add_argument("--output", default=".", help="Output directory")
    parser.add_argument("--no-skip", action="store_true", help="Don't skip existing")
    args = parser.parse_args()

    download_samples(args.input, args.output, skip_existing=not args.no_skip)
