#!/usr/bin/env python3
"""
Run full IS element extraction + annotation on samples.
Uses multiprocessing to parallelize across CPUs.

Processes samples from a list or auto-discovers incomplete ones:
  Step 1: Extract IS elements (sequence, flanking, insertion site)
  Step 2: Annotate with Prodigal (ORFs, protein sequences, noncoding regions)

Output:
  - Per-sample JSON: samples/{SAMPLE}/is_extraction/is_elements.json
  - Combined summary TSV: {output-dir}/is_elements_summary.tsv
"""

import argparse
import os
import sys
import json
import csv
import time
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# Some IS sequences are very long, increase CSV field size limit
csv.field_size_limit(sys.maxsize)

# Add project root to path
sys.path.insert(0, '/global/home/users/kh36969/mgefinder_project')


def parse_args():
    parser = argparse.ArgumentParser(
        description='Extract and annotate IS elements from MGEfinder samples.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  # Auto-discover incomplete samples
  python %(prog)s --samples-dir /path/to/samples

  # Use explicit sample list
  python %(prog)s --samples-dir /path/to/samples --sample-list samples_with_is.txt

  # Custom output dir and worker count
  python %(prog)s --samples-dir /path/to/samples --output-dir /path/to/results --workers 16
""")
    parser.add_argument('--samples-dir', required=True,
                        help='Path to directory containing sample subdirectories')
    parser.add_argument('--sample-list', default=None,
                        help='Path to file with one sample ID per line. '
                             'If omitted, auto-discovers samples that have IS FASTA '
                             'but no is_extraction/is_elements.json')
    parser.add_argument('--output-dir', default=None,
                        help='Directory for summary TSV output '
                             '(default: {samples-dir}/../is_extraction_results)')
    parser.add_argument('--workers', type=int, default=30,
                        help='Number of parallel workers (default: 30)')
    parser.add_argument('--flank-size', type=int, default=80,
                        help='Flanking sequence size in bp (default: 80)')
    return parser.parse_args()


# These globals are set in main() and read by worker processes via fork
_SAMPLES_DIR = None
_FLANK_SIZE = None


def process_sample(sample_id):
    """
    Process a single sample: extract IS elements + annotate with Prodigal.
    Each worker runs in its own process with its own temp directory.

    Returns a dict with status and results for the main process to collect.
    """
    # Must set CSV limit in each worker process
    import csv as _csv
    _csv.field_size_limit(sys.maxsize)

    from is_extractor import extract_is_elements, annotate_is_elements

    sample_dir = _SAMPLES_DIR / sample_id
    temp_dir = tempfile.mkdtemp(prefix=f'is_{sample_id}_')

    try:
        # Step 1: Extract IS elements
        results = extract_is_elements(str(sample_dir), flank_size=_FLANK_SIZE)

        if not results:
            return {'sample_id': sample_id, 'status': 'empty', 'records': []}

        # Step 2: Annotate with Prodigal
        annotated = annotate_is_elements(results, temp_dir=temp_dir)

        # Save per-sample JSON
        sample_out_dir = sample_dir / 'is_extraction'
        sample_out_dir.mkdir(parents=True, exist_ok=True)
        json_path = sample_out_dir / 'is_elements.json'
        with open(json_path, 'w') as jf:
            json.dump(annotated, jf, indent=2)

        # Build summary rows for TSV
        summary_rows = []
        for rec in annotated:
            orf_ann = rec.get('orf_annotation', {})
            has_annotation = 'error' not in orf_ann and 'orfs' in orf_ann

            summary_rows.append({
                'is_id': rec['is_id'],
                'sample': rec['sample'],
                'seqid': rec['seqid'],
                'cluster': rec['cluster'],
                'group': rec['group'],
                'pair_id': rec['pair_id'],
                'confidence': rec['confidence'],
                'is_length': rec['is_element']['length'],
                'is_method': rec['is_element']['method'],
                'is_contig': rec['is_element']['contig'],
                'is_start': rec['is_element']['start'],
                'is_end': rec['is_element']['end'],
                'insertion_contig': rec['insertion_site']['contig'],
                'pos_5p': rec['insertion_site']['pos_5p'],
                'pos_3p': rec['insertion_site']['pos_3p'],
                'upstream_length': rec['flanking_upstream']['length'],
                'downstream_length': rec['flanking_downstream']['length'],
                'num_orfs': orf_ann.get('num_orfs', 0) if has_annotation else '',
                'num_noncoding_regions': orf_ann.get('num_noncoding_regions', 0) if has_annotation else '',
                'total_coding_bp': orf_ann.get('total_coding_bp', 0) if has_annotation else '',
                'total_noncoding_bp': orf_ann.get('total_noncoding_bp', 0) if has_annotation else '',
                'coding_fraction': f"{orf_ann.get('coding_fraction', 0):.3f}" if has_annotation else '',
            })

        n_orfs = sum(
            rec.get('orf_annotation', {}).get('num_orfs', 0)
            for rec in annotated
            if 'error' not in rec.get('orf_annotation', {})
        )

        return {
            'sample_id': sample_id,
            'status': 'ok',
            'n_is': len(annotated),
            'n_orfs': n_orfs,
            'records': summary_rows,
        }

    except Exception as e:
        return {'sample_id': sample_id, 'status': 'error', 'error': str(e), 'records': []}

    finally:
        # Cleanup temp dir
        try:
            for f in os.listdir(temp_dir):
                os.remove(os.path.join(temp_dir, f))
            os.rmdir(temp_dir)
        except OSError:
            pass


def discover_samples(samples_dir):
    """Auto-discover samples that have IS FASTA but no is_elements.json yet."""
    is_fasta_name = '04.makefasta.genome.all_seqs.fna'
    samples = []
    for entry in sorted(os.listdir(samples_dir)):
        sample_dir = samples_dir / entry
        if not sample_dir.is_dir():
            continue
        is_fasta = sample_dir / '03.results' / 'genome' / is_fasta_name
        no_mge = sample_dir / '03.results' / 'genome' / 'NO_MGE_DETECTED.txt'
        json_out = sample_dir / 'is_extraction' / 'is_elements.json'
        if is_fasta.exists() and is_fasta.stat().st_size > 0 \
                and not no_mge.exists() and not json_out.exists():
            samples.append(entry)
    return samples


def main():
    global _SAMPLES_DIR, _FLANK_SIZE

    args = parse_args()
    start_time = time.time()

    _SAMPLES_DIR = Path(args.samples_dir).resolve()
    _FLANK_SIZE = args.flank_size

    output_dir = Path(args.output_dir) if args.output_dir \
        else _SAMPLES_DIR.parent / 'is_extraction_results'
    output_dir.mkdir(parents=True, exist_ok=True)

    # Build sample list
    if args.sample_list:
        with open(args.sample_list) as f:
            sample_ids = [line.strip() for line in f if line.strip()]
    else:
        print("No --sample-list provided, auto-discovering samples...")
        sample_ids = discover_samples(_SAMPLES_DIR)

    # Skip already-completed samples
    todo = []
    skipped = 0
    for sid in sample_ids:
        json_path = _SAMPLES_DIR / sid / 'is_extraction' / 'is_elements.json'
        if json_path.exists():
            skipped += 1
        else:
            todo.append(sid)

    print(f"Total samples: {len(sample_ids)}")
    print(f"Already completed: {skipped}")
    print(f"To process: {len(todo)}")
    print(f"Workers: {args.workers}")
    print(f"Flank size: {_FLANK_SIZE} bp")
    print(f"Output directory: {output_dir}")
    print(f"Steps: extraction + Prodigal annotation")
    print(flush=True)

    if not todo:
        print("Nothing to do!")
        return

    # Counters
    total_is = 0
    total_orfs = 0
    processed = 0
    failed = 0
    failed_samples = []
    empty = 0

    # Summary TSV — append to existing file so reruns don't lose previous results
    summary_path = output_dir / 'is_elements_summary.tsv'
    summary_fields = [
        'is_id', 'sample', 'seqid', 'cluster', 'group', 'pair_id',
        'confidence', 'is_length', 'is_method', 'is_contig',
        'is_start', 'is_end',
        'insertion_contig', 'pos_5p', 'pos_3p',
        'upstream_length', 'downstream_length',
        'num_orfs', 'num_noncoding_regions',
        'total_coding_bp', 'total_noncoding_bp', 'coding_fraction'
    ]

    write_header = not summary_path.exists() or summary_path.stat().st_size == 0
    with open(summary_path, 'a', newline='') as summary_f:
        writer = csv.DictWriter(summary_f, fieldnames=summary_fields, delimiter='\t')
        if write_header:
            writer.writeheader()

        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = {executor.submit(process_sample, sid): sid for sid in todo}

            for future in as_completed(futures):
                result = future.result()
                sid = result['sample_id']

                if result['status'] == 'ok':
                    processed += 1
                    total_is += result['n_is']
                    total_orfs += result['n_orfs']
                    for row in result['records']:
                        writer.writerow(row)
                    summary_f.flush()
                elif result['status'] == 'empty':
                    empty += 1
                else:
                    failed += 1
                    failed_samples.append(sid)
                    print(f"  ERROR [{sid}]: {result.get('error', 'unknown')}", file=sys.stderr)

                done = processed + empty + failed
                if done % 100 == 0:
                    elapsed = time.time() - start_time
                    rate = done / elapsed if elapsed > 0 else 0
                    eta = (len(todo) - done) / rate if rate > 0 else 0
                    print(f"  Progress: {done}/{len(todo)} "
                          f"({processed} ok, {empty} empty, {failed} failed) "
                          f"| {total_is} IS, {total_orfs} ORFs "
                          f"| {elapsed:.0f}s elapsed, ~{eta:.0f}s ETA",
                          flush=True)

    elapsed = time.time() - start_time

    print(f"\n{'='*60}")
    print(f"DONE in {elapsed:.1f}s ({elapsed/60:.1f} min)")
    print(f"{'='*60}")
    print(f"  Samples with IS extracted: {processed}")
    print(f"  Samples empty (no IS): {empty}")
    print(f"  Samples failed: {failed}")
    print(f"  Total IS elements: {total_is}")
    print(f"  Total ORFs predicted: {total_orfs}")
    print(f"  Summary TSV: {summary_path}")
    print(f"  Per-sample JSONs: samples/{{SAMPLE}}/is_extraction/is_elements.json")

    if failed_samples:
        failed_path = output_dir / 'failed_samples.txt'
        with open(failed_path, 'w') as f:
            for s in failed_samples:
                f.write(s + '\n')
        print(f"  Failed sample list: {failed_path}")


if __name__ == '__main__':
    main()
