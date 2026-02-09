#!/usr/bin/env python3
"""
Verify IS extraction results.

Loads an is_elements.json file and verifies that:
  1. upstream_flank + IS_sequence + downstream_flank is contiguous on the
     assembly contig (for assembly-mapped IS elements).
  2. ORF and noncoding regions tile the IS sequence completely (no gaps,
     no overlaps beyond gene overlap).
  3. Each ORF's dna_sequence matches the IS subsequence at its coordinates.
  4. Each noncoding region's sequence matches the IS subsequence.

Usage:
    python verify_is_extraction.py <sample_dir>
    python verify_is_extraction.py <is_elements.json> --sample-dir <sample_dir>
    python verify_is_extraction.py --batch <samples_dir> [--max-samples N]
    python verify_is_extraction.py --live <sample_dir>  (run extractor + verify, no JSON needed)
"""

import argparse
import json
import os
import sys
from pathlib import Path

# Add project root to path so is_extractor can be imported
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from is_extractor.extractor import parse_fasta, reverse_complement


def load_assembly(sample_dir: str) -> dict:
    """Load assembly FASTA for a sample."""
    sample_dir = Path(sample_dir)
    sample_id = sample_dir.name
    assembly_path = sample_dir / "00.assembly" / f"{sample_id}.fna"
    if assembly_path.exists():
        return parse_fasta(str(assembly_path))
    return {}


def verify_contiguity(record: dict, assembly: dict) -> dict:
    """
    Verify that upstream_flank + IS_sequence + downstream_flank is contiguous
    on the assembly contig.

    Returns dict with 'pass', 'skip', 'reason' keys.
    """
    contig_name = record['is_element'].get('contig', '')
    strand = record['is_element'].get('strand', '?')
    start = record['is_element'].get('start', 0)
    end = record['is_element'].get('end', 0)

    # Skip if no assembly coordinates
    if not contig_name or strand == '?' or start >= end:
        return {'pass': None, 'skip': True,
                'reason': f'no assembly coords (contig={contig_name}, strand={strand})'}

    if contig_name not in assembly:
        return {'pass': None, 'skip': True,
                'reason': f'contig {contig_name} not in assembly'}

    contig_seq = assembly[contig_name]
    contig_len = len(contig_seq)

    up_seq = record['flanking_upstream'].get('sequence', '')
    is_seq = record['is_element'].get('sequence', '')
    down_seq = record['flanking_downstream'].get('sequence', '')

    if not is_seq:
        return {'pass': None, 'skip': True, 'reason': 'empty IS sequence'}

    # Build expected contig region.
    # For forward strand (+):
    #   contig order: upstream_flank | IS | downstream_flank
    #   region = contig[start - len(up) : end + len(down)]
    #
    # For reverse strand (-):
    #   contig order: revcomp(downstream) | revcomp(IS) | revcomp(upstream)
    #   The upstream (in IS frame) comes from contig RIGHT of IS,
    #   the downstream (in IS frame) comes from contig LEFT of IS.
    #   region = contig[start - len(down) : end + len(up)]
    if strand == '+':
        region_start = max(0, start - len(up_seq))
        region_end = min(contig_len, end + len(down_seq))
    else:
        region_start = max(0, start - len(down_seq))
        region_end = min(contig_len, end + len(up_seq))

    contig_region = contig_seq[region_start:region_end]

    # Concatenate extracted sequences (in IS reading frame)
    concat = up_seq + is_seq + down_seq

    if strand == '+':
        match = (concat == contig_region)
    else:
        match = (reverse_complement(concat) == contig_region)

    if match:
        return {'pass': True, 'skip': False, 'reason': 'contiguous on contig'}
    else:
        return {
            'pass': False, 'skip': False,
            'reason': (f'MISMATCH: concat len={len(concat)}, '
                       f'contig region len={len(contig_region)}, strand={strand}')
        }


def verify_orf_tiling(record: dict) -> dict:
    """
    Verify that ORFs + noncoding regions tile the IS sequence completely.

    Checks:
      - All positions from 1..is_length are covered by ORFs or noncoding regions
      - Each ORF dna_sequence matches IS subsequence
      - Each noncoding region sequence matches IS subsequence

    Returns dict with 'pass', 'skip', 'errors' keys.
    """
    annotation = record.get('orf_annotation')
    if not annotation or 'error' in annotation:
        return {'pass': None, 'skip': True,
                'reason': 'no orf_annotation or annotation error'}

    is_seq = record['is_element'].get('sequence', '')
    is_length = len(is_seq)
    orfs = annotation.get('orfs', [])
    noncoding = annotation.get('noncoding_regions', [])
    errors = []

    # 1. Check each ORF's dna_sequence matches IS subsequence
    for orf in orfs:
        orf_id = orf.get('orf_id', '?')
        orf_start = orf.get('start', 0)  # 1-based
        orf_end = orf.get('end', 0)      # 1-based
        orf_dna = orf.get('dna_sequence', '')

        if not orf_dna:
            errors.append(f'{orf_id}: empty dna_sequence')
            continue

        expected = is_seq[orf_start - 1:orf_end]
        if orf_dna != expected:
            errors.append(
                f'{orf_id}: dna_sequence mismatch at [{orf_start}:{orf_end}] '
                f'(got len={len(orf_dna)}, expected len={len(expected)})'
            )

    # 2. Check each noncoding region's sequence matches IS subsequence
    for nc in noncoding:
        nc_start = nc.get('start', 0)
        nc_end = nc.get('end', 0)
        nc_seq = nc.get('sequence', '')
        nc_type = nc.get('type', '?')

        if not nc_seq:
            errors.append(f'noncoding {nc_type} [{nc_start}:{nc_end}]: empty sequence')
            continue

        expected = is_seq[nc_start - 1:nc_end]
        if nc_seq != expected:
            errors.append(
                f'noncoding {nc_type} [{nc_start}:{nc_end}]: sequence mismatch '
                f'(got len={len(nc_seq)}, expected len={len(expected)})'
            )

    # 3. Check that ORFs + noncoding cover every position in [1..is_length]
    #    Collect all intervals
    intervals = []
    for orf in orfs:
        intervals.append((orf.get('start', 0), orf.get('end', 0), 'orf'))
    for nc in noncoding:
        intervals.append((nc.get('start', 0), nc.get('end', 0), 'noncoding'))

    if intervals:
        intervals.sort()
        covered = set()
        for s, e, _ in intervals:
            for pos in range(s, e + 1):
                covered.add(pos)

        expected_positions = set(range(1, is_length + 1))
        missing = expected_positions - covered
        if missing:
            # Report ranges for readability
            missing_sorted = sorted(missing)
            ranges = []
            range_start = missing_sorted[0]
            range_end = missing_sorted[0]
            for pos in missing_sorted[1:]:
                if pos == range_end + 1:
                    range_end = pos
                else:
                    ranges.append(f'{range_start}-{range_end}')
                    range_start = pos
                    range_end = pos
            ranges.append(f'{range_start}-{range_end}')
            errors.append(f'uncovered positions: {", ".join(ranges[:5])}'
                          f'{" ..." if len(ranges) > 5 else ""} '
                          f'({len(missing)} bp total)')

    passed = len(errors) == 0
    return {'pass': passed, 'skip': False, 'errors': errors}


def verify_sample(sample_dir: str, json_path: str = None,
                  verbose: bool = False) -> dict:
    """
    Run all verifications for a single sample.

    Args:
        sample_dir: Path to sample directory
        json_path: Path to JSON file (default: {sample_dir}/is_extraction/is_elements.json)
        verbose: Print per-element details

    Returns:
        Summary dict with counts
    """
    sample_dir = Path(sample_dir)
    sample_id = sample_dir.name

    if json_path is None:
        json_path = sample_dir / "is_extraction" / "is_elements.json"
    else:
        json_path = Path(json_path)

    if not json_path.exists():
        return {'sample': sample_id, 'error': f'JSON not found: {json_path}'}

    with open(json_path) as f:
        records = json.load(f)

    if not records:
        return {'sample': sample_id, 'total': 0, 'note': 'empty JSON'}

    # Load assembly
    assembly = load_assembly(str(sample_dir))

    contiguity_pass = 0
    contiguity_fail = 0
    contiguity_skip = 0
    tiling_pass = 0
    tiling_fail = 0
    tiling_skip = 0
    tiling_errors = []

    for record in records:
        is_id = record.get('is_id', '?')

        # Test 1: contiguity
        c_result = verify_contiguity(record, assembly)
        if c_result['skip']:
            contiguity_skip += 1
            if verbose:
                print(f'  {is_id} contiguity: SKIP ({c_result["reason"]})')
        elif c_result['pass']:
            contiguity_pass += 1
            if verbose:
                print(f'  {is_id} contiguity: PASS')
        else:
            contiguity_fail += 1
            if verbose:
                print(f'  {is_id} contiguity: FAIL ({c_result["reason"]})')

        # Test 2: ORF/noncoding tiling
        t_result = verify_orf_tiling(record)
        if t_result['skip']:
            tiling_skip += 1
            if verbose:
                print(f'  {is_id} tiling: SKIP ({t_result.get("reason", "")})')
        elif t_result['pass']:
            tiling_pass += 1
            if verbose:
                print(f'  {is_id} tiling: PASS')
        else:
            tiling_fail += 1
            if verbose:
                for err in t_result['errors']:
                    print(f'  {is_id} tiling: FAIL - {err}')
            tiling_errors.extend(
                [(is_id, e) for e in t_result.get('errors', [])])

    summary = {
        'sample': sample_id,
        'total': len(records),
        'contiguity': {
            'pass': contiguity_pass,
            'fail': contiguity_fail,
            'skip': contiguity_skip,
        },
        'tiling': {
            'pass': tiling_pass,
            'fail': tiling_fail,
            'skip': tiling_skip,
        },
    }
    if tiling_errors:
        summary['tiling_errors'] = tiling_errors

    return summary


def print_summary(summary: dict):
    """Print a sample verification summary."""
    sample = summary.get('sample', '?')
    total = summary.get('total', 0)

    if 'error' in summary:
        print(f'{sample}: ERROR - {summary["error"]}')
        return

    if total == 0:
        print(f'{sample}: no IS elements')
        return

    c = summary['contiguity']
    t = summary['tiling']

    c_status = 'OK' if c['fail'] == 0 else 'FAIL'
    t_status = 'OK' if t['fail'] == 0 else 'FAIL'

    print(f'{sample}: {total} IS elements | '
          f'contiguity: {c["pass"]} pass, {c["fail"]} fail, {c["skip"]} skip [{c_status}] | '
          f'tiling: {t["pass"]} pass, {t["fail"]} fail, {t["skip"]} skip [{t_status}]')

    if 'tiling_errors' in summary:
        for is_id, err in summary['tiling_errors'][:5]:
            print(f'  {is_id}: {err}')
        if len(summary['tiling_errors']) > 5:
            print(f'  ... and {len(summary["tiling_errors"]) - 5} more errors')


def main():
    parser = argparse.ArgumentParser(
        description='Verify IS extraction results against assembly contigs.')
    parser.add_argument('input', nargs='?',
                        help='Sample directory or is_elements.json file')
    parser.add_argument('--sample-dir', metavar='DIR',
                        help='Sample directory (required when input is a JSON file)')
    parser.add_argument('--batch', metavar='SAMPLES_DIR',
                        help='Run on all samples in a directory')
    parser.add_argument('--live', metavar='SAMPLE_DIR',
                        help='Run extractor on a sample and verify (no JSON needed)')
    parser.add_argument('--max-samples', type=int, default=0,
                        help='Maximum number of samples to check (0=all)')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Print per-element details')

    args = parser.parse_args()

    if not args.input and not args.batch and not args.live:
        parser.print_help()
        sys.exit(1)

    if args.live:
        # Live mode: run extractor then verify against assembly
        from is_extractor.extractor import extract_is_elements
        sample_dir = Path(args.live)
        sample_id = sample_dir.name
        assembly = load_assembly(str(sample_dir))

        print(f'Running extractor on {sample_id}...')
        records = extract_is_elements(str(sample_dir))
        print(f'Extracted {len(records)} IS elements')
        print()

        c_pass = c_fail = c_skip = 0
        for record in records:
            is_id = record.get('is_id', '?')
            result = verify_contiguity(record, assembly)
            if result['skip']:
                c_skip += 1
                if args.verbose:
                    print(f'  {is_id} contiguity: SKIP ({result["reason"]})')
            elif result['pass']:
                c_pass += 1
                if args.verbose:
                    print(f'  {is_id} contiguity: PASS (strand={record["is_element"]["strand"]})')
            else:
                c_fail += 1
                print(f'  {is_id} contiguity: FAIL ({result["reason"]})')

        status = 'OK' if c_fail == 0 else 'FAIL'
        print(f'\n{sample_id}: contiguity {c_pass} pass, {c_fail} fail, {c_skip} skip [{status}]')
        sys.exit(1 if c_fail > 0 else 0)

    elif args.batch:
        # Batch mode: check all samples with existing JSON files
        samples_dir = Path(args.batch)
        samples_checked = 0
        total_pass = 0
        total_fail = 0
        total_skip = 0
        failed_samples = []

        for sample_path in sorted(samples_dir.iterdir()):
            json_path = sample_path / "is_extraction" / "is_elements.json"
            if not json_path.exists():
                continue

            summary = verify_sample(str(sample_path), verbose=args.verbose)
            print_summary(summary)

            if 'contiguity' in summary:
                total_pass += summary['contiguity']['pass']
                total_fail += summary['contiguity']['fail']
                total_skip += summary['contiguity']['skip']
                if summary['contiguity']['fail'] > 0:
                    failed_samples.append(summary['sample'])

            samples_checked += 1
            if args.max_samples and samples_checked >= args.max_samples:
                break

        print()
        print(f'=== BATCH SUMMARY ===')
        print(f'Samples checked: {samples_checked}')
        print(f'Contiguity totals: {total_pass} pass, {total_fail} fail, {total_skip} skip')
        if failed_samples:
            print(f'Failed samples: {", ".join(failed_samples[:10])}'
                  f'{"..." if len(failed_samples) > 10 else ""}')
        else:
            print(f'All assembly-mapped IS elements passed contiguity check.')

    else:
        # Single sample or JSON file
        input_path = Path(args.input)

        if input_path.suffix == '.json':
            # JSON file — need sample_dir for assembly
            if args.sample_dir:
                sample_dir = args.sample_dir
            else:
                # Try to infer: {sample_dir}/is_extraction/is_elements.json
                sample_dir = str(input_path.parent.parent)
            summary = verify_sample(sample_dir, str(input_path),
                                    verbose=args.verbose)
        else:
            # Sample directory
            summary = verify_sample(str(input_path), verbose=args.verbose)

        print_summary(summary)

        # Exit code: 0 if all pass, 1 if any fail
        if 'contiguity' in summary:
            if summary['contiguity']['fail'] > 0:
                sys.exit(1)
            if summary['tiling']['fail'] > 0:
                sys.exit(1)


if __name__ == '__main__':
    main()
