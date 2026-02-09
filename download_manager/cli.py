#!/usr/bin/env python3
"""
Command-line interface for the Download Manager.

Usage:
    python -m download_manager.cli match --reads FILE --assembly FILE
    python -m download_manager.cli build --matched FILE --reads FILE --assembly FILE
    python -m download_manager.cli subset --input FILE --size 10000 --batch-name batch3
    python -m download_manager.cli download --input FILE --output DIR
    python -m download_manager.cli register --input FILE --batch-name batch1
    python -m download_manager.cli status
    python -m download_manager.cli pipeline  # Run full workflow
"""

import argparse
import sys
from pathlib import Path

from .matcher import match_samples
from .master_list import build_master_list
from .subset import select_subset, load_exclude_file
from .downloader import download_samples
from .registry import register_from_file, show_status
from .config import (
    DEFAULT_READS_FILE,
    DEFAULT_ASSEMBLY_FILE,
    DEFAULT_MASTER_LIST,
    DEFAULT_TARGET_FILE,
    DEFAULT_REGISTRY_FILE,
)


def cmd_match(args):
    """Run sample matching step."""
    match_samples(
        reads_file=args.reads,
        assembly_file=args.assembly,
        output_file=args.output
    )


def cmd_build(args):
    """Run master list building step."""
    build_master_list(
        matched_samples_file=args.matched,
        reads_file=args.reads,
        assembly_file=args.assembly,
        output_file=args.output
    )


def cmd_subset(args):
    """Run subset selection step."""
    exclude = load_exclude_file(args.exclude) if args.exclude else None
    registry = str(args.registry) if args.registry else None

    select_subset(
        master_list_file=args.input,
        output_file=args.output,
        target_size=args.size,
        seed=args.seed,
        exclude_ids=exclude,
        registry_file=registry,
        batch_name=args.batch_name,
        no_register=args.no_register,
    )


def cmd_download(args):
    """Run download step."""
    download_samples(
        target_file=args.input,
        output_dir=args.output,
        skip_existing=not args.no_skip
    )


def cmd_register(args):
    """Register samples from a file into the registry."""
    registry = str(args.registry)
    n = register_from_file(
        registry_file=registry,
        input_file=args.input,
        batch_name=args.batch_name,
        date_added=args.date,
    )
    print(f"Registered {n:,} new samples as '{args.batch_name}'")
    print()
    show_status(registry)


def cmd_status(args):
    """Show registry status."""
    show_status(str(args.registry))


def cmd_pipeline(args):
    """Run full pipeline with default paths."""
    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    matched_file = work_dir / "matched_samples.txt"
    master_file = work_dir / "master_list.tsv"
    target_file = work_dir / f"target_{args.size}.tsv"
    registry = str(args.registry)

    print("=" * 60)
    print("STEP 1: Matching samples")
    print("=" * 60)
    match_samples(
        reads_file=str(DEFAULT_READS_FILE),
        assembly_file=str(DEFAULT_ASSEMBLY_FILE),
        output_file=str(matched_file)
    )

    print("\n" + "=" * 60)
    print("STEP 2: Building master list")
    print("=" * 60)
    build_master_list(
        matched_samples_file=str(matched_file),
        reads_file=str(DEFAULT_READS_FILE),
        assembly_file=str(DEFAULT_ASSEMBLY_FILE),
        output_file=str(master_file)
    )

    print("\n" + "=" * 60)
    print(f"STEP 3: Selecting {args.size} samples")
    print("=" * 60)
    select_subset(
        master_list_file=str(master_file),
        output_file=str(target_file),
        target_size=args.size,
        seed=args.seed,
        registry_file=registry,
        batch_name=args.batch_name,
        no_register=args.no_register,
    )

    if not args.no_download:
        print("\n" + "=" * 60)
        print("STEP 4: Downloading samples")
        print("=" * 60)
        download_samples(
            target_file=str(target_file),
            output_dir=str(work_dir / "samples"),
            skip_existing=True
        )

    print("\n" + "=" * 60)
    print("Pipeline complete!")
    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(
        description="ENA Download Manager for MGEfinder",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # Match command
    p_match = subparsers.add_parser('match', help='Find samples with both reads and assembly')
    p_match.add_argument('--reads', required=True, help='ENA reads TSV file')
    p_match.add_argument('--assembly', required=True, help='ENA assembly TSV file')
    p_match.add_argument('--output', default='matched_samples.txt', help='Output file')
    p_match.set_defaults(func=cmd_match)

    # Build command
    p_build = subparsers.add_parser('build', help='Build master list from metadata')
    p_build.add_argument('--matched', required=True, help='Matched samples file')
    p_build.add_argument('--reads', required=True, help='ENA reads TSV file')
    p_build.add_argument('--assembly', required=True, help='ENA assembly TSV file')
    p_build.add_argument('--output', default='master_list.tsv', help='Output file')
    p_build.set_defaults(func=cmd_build)

    # Subset command
    p_subset = subparsers.add_parser('subset',
        help='Select diverse subset of samples (auto-excludes previous batches)')
    p_subset.add_argument('--input', required=True, help='Master list TSV file')
    p_subset.add_argument('--output', default='target_samples.tsv', help='Output file')
    p_subset.add_argument('--size', type=int, default=10000, help='Target sample count')
    p_subset.add_argument('--seed', type=int, help='Random seed')
    p_subset.add_argument('--batch-name', help='Batch label (e.g., "batch3"). '
                          'Required to register picks.')
    p_subset.add_argument('--registry', type=Path, default=DEFAULT_REGISTRY_FILE,
                          help=f'Sample registry file (default: {DEFAULT_REGISTRY_FILE})')
    p_subset.add_argument('--exclude', help='Additional exclude file (one ID per line)')
    p_subset.add_argument('--no-register', action='store_true',
                          help="Don't register new picks in the registry")
    p_subset.set_defaults(func=cmd_subset)

    # Download command
    p_download = subparsers.add_parser('download', help='Download samples from ENA')
    p_download.add_argument('--input', required=True, help='Target samples TSV file')
    p_download.add_argument('--output', default='.', help='Output directory')
    p_download.add_argument('--no-skip', action='store_true', help="Don't skip existing")
    p_download.set_defaults(func=cmd_download)

    # Register command
    p_register = subparsers.add_parser('register',
        help='Register samples from a file into the registry')
    p_register.add_argument('--input', required=True,
                           help='Sample list: TSV with Sample_ID column, or '
                                'plain text with one ID per line')
    p_register.add_argument('--batch-name', required=True,
                           help='Batch label (e.g., "batch1")')
    p_register.add_argument('--registry', type=Path, default=DEFAULT_REGISTRY_FILE,
                           help=f'Registry file (default: {DEFAULT_REGISTRY_FILE})')
    p_register.add_argument('--date', default=None,
                           help='Date to record (default: today)')
    p_register.set_defaults(func=cmd_register)

    # Status command
    p_status = subparsers.add_parser('status',
        help='Show sample registry summary')
    p_status.add_argument('--registry', type=Path, default=DEFAULT_REGISTRY_FILE,
                         help=f'Registry file (default: {DEFAULT_REGISTRY_FILE})')
    p_status.set_defaults(func=cmd_status)

    # Pipeline command
    p_pipeline = subparsers.add_parser('pipeline', help='Run full workflow')
    p_pipeline.add_argument('--work-dir', default='.', help='Working directory')
    p_pipeline.add_argument('--size', type=int, default=10000, help='Target sample count')
    p_pipeline.add_argument('--seed', type=int, help='Random seed')
    p_pipeline.add_argument('--batch-name', help='Batch label for registry')
    p_pipeline.add_argument('--registry', type=Path, default=DEFAULT_REGISTRY_FILE,
                           help=f'Registry file (default: {DEFAULT_REGISTRY_FILE})')
    p_pipeline.add_argument('--no-register', action='store_true',
                           help="Don't register picks")
    p_pipeline.add_argument('--no-download', action='store_true',
                           help='Skip download step')
    p_pipeline.set_defaults(func=cmd_pipeline)

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    args.func(args)


if __name__ == "__main__":
    main()
