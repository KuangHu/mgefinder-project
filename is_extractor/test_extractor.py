#!/usr/bin/env python
"""
Test script for IS Extractor module.

Validates that the dictionary location info can be used to extract
the same sequences from the original source files.
"""

import sys
from pathlib import Path
from extractor import ISExtractor, parse_fasta


def extract_sequence_from_file(source_file: str, contig: str,
                                start: int, end: int) -> str:
    """
    Extract sequence directly from source file using coordinates.

    Args:
        source_file: Path to FASTA file
        contig: Contig/sequence name
        start: Start position (1-based)
        end: End position (1-based, inclusive)

    Returns:
        Extracted sequence string
    """
    sequences = parse_fasta(source_file)

    if contig not in sequences:
        return f"ERROR: Contig '{contig}' not found in {source_file}"

    seq = sequences[contig]
    # Convert to 0-based indexing
    extracted = seq[start-1:end]
    return extracted


def test_flanking_sequences(sample_dir: str, verbose: bool = True):
    """
    Test that flanking sequences in dictionary match source files.

    Args:
        sample_dir: Path to sample directory
        verbose: Print detailed output
    """
    print(f"Testing sample: {sample_dir}")
    print("=" * 70)

    extractor = ISExtractor(sample_dir)

    if not extractor.has_mge_detected():
        print("No MGE detected in this sample.")
        return True

    results = extractor.extract()
    print(f"Found {len(results)} IS elements\n")

    all_passed = True
    tested = 0
    passed = 0

    for i, record in enumerate(results):
        is_id = record['is_id']

        # Test upstream flanking
        upstream = record['flanking_upstream']
        if upstream['sequence'] and upstream['source_file']:
            tested += 1
            extracted_up = extract_sequence_from_file(
                upstream['source_file'],
                upstream['contig'],
                upstream['start'],
                upstream['end']
            )

            up_match = (extracted_up == upstream['sequence'])
            if up_match:
                passed += 1
            else:
                all_passed = False

            if verbose or not up_match:
                print(f"[{is_id}] Upstream flanking:")
                print(f"  Location: {upstream['contig']}:{upstream['start']}-{upstream['end']}")
                print(f"  Dict seq:      {upstream['sequence'][:50]}...")
                print(f"  Extracted seq: {extracted_up[:50]}...")
                print(f"  Match: {'PASS' if up_match else 'FAIL'}")
                print()

        # Test downstream flanking
        downstream = record['flanking_downstream']
        if downstream['sequence'] and downstream['source_file']:
            tested += 1
            extracted_down = extract_sequence_from_file(
                downstream['source_file'],
                downstream['contig'],
                downstream['start'],
                downstream['end']
            )

            down_match = (extracted_down == downstream['sequence'])
            if down_match:
                passed += 1
            else:
                all_passed = False

            if verbose or not down_match:
                print(f"[{is_id}] Downstream flanking:")
                print(f"  Location: {downstream['contig']}:{downstream['start']}-{downstream['end']}")
                print(f"  Dict seq:      {downstream['sequence'][:50]}...")
                print(f"  Extracted seq: {extracted_down[:50]}...")
                print(f"  Match: {'PASS' if down_match else 'FAIL'}")
                print()

        # Only show first 3 in verbose mode
        if verbose and i >= 2:
            print(f"... (showing first 3 of {len(results)} elements)")
            print()
            break

    # Summary
    print("=" * 70)
    print(f"SUMMARY: {passed}/{tested} flanking sequence tests passed")

    if all_passed:
        print("STATUS: ALL TESTS PASSED")
    else:
        print("STATUS: SOME TESTS FAILED")

    return all_passed


def test_is_element_sequence(sample_dir: str, verbose: bool = True):
    """
    Test that IS element sequences exist in the clusterseq file.

    Note: Same IS element can appear at multiple insertion sites, so we
    verify that the sequence EXISTS in the source file (not exact row match).

    Args:
        sample_dir: Path to sample directory
        verbose: Print detailed output
    """
    print(f"\nTesting IS element sequences: {sample_dir}")
    print("=" * 70)

    extractor = ISExtractor(sample_dir)

    if not extractor.has_mge_detected():
        print("No MGE detected in this sample.")
        return True

    results = extractor.extract()

    # Read all sequences from clusterseq file
    import csv
    clusterseq_path = extractor.clusterseq_path
    all_sequences = set()

    with open(clusterseq_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            seq = row.get('inferred_seq', '')
            if seq:
                all_sequences.add(seq)

    print(f"Loaded {len(all_sequences)} unique sequences from clusterseq file\n")

    all_passed = True
    tested = 0
    passed = 0

    for i, record in enumerate(results):
        is_id = record['is_id']
        is_elem = record['is_element']

        if is_elem['sequence']:
            tested += 1

            # Check if sequence exists in source file
            match = is_elem['sequence'] in all_sequences

            if match:
                passed += 1
            else:
                all_passed = False

            if verbose or not match:
                print(f"[{is_id}] IS element sequence:")
                print(f"  Length: {is_elem['length']} bp")
                print(f"  Location: {is_elem['contig']}:{is_elem['start']}-{is_elem['end']}")
                print(f"  Dict seq: {is_elem['sequence'][:50]}...")
                print(f"  Found in source: {'PASS' if match else 'FAIL'}")
                print()

        if verbose and i >= 2:
            print(f"... (showing first 3 of {len(results)} elements)")
            print()
            break

    print("=" * 70)
    print(f"SUMMARY: {passed}/{tested} IS element sequence tests passed")

    if all_passed:
        print("STATUS: ALL TESTS PASSED")
    else:
        print("STATUS: SOME TESTS FAILED")

    return all_passed


def run_all_tests(sample_dir: str, verbose: bool = True):
    """Run all validation tests."""
    print("\n" + "=" * 70)
    print("IS EXTRACTOR VALIDATION TESTS")
    print("=" * 70 + "\n")

    test1 = test_flanking_sequences(sample_dir, verbose)
    test2 = test_is_element_sequence(sample_dir, verbose)

    print("\n" + "=" * 70)
    print("FINAL RESULT")
    print("=" * 70)

    if test1 and test2:
        print("ALL TESTS PASSED - Dictionary info correctly maps to source files")
        return True
    else:
        print("SOME TESTS FAILED - Check output above for details")
        return False


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python test_extractor.py <sample_dir> [--quiet]")
        print()
        print("Example:")
        print("  python test_extractor.py /path/to/samples/SAMN18872762")
        sys.exit(1)

    sample_dir = sys.argv[1]
    verbose = "--quiet" not in sys.argv

    success = run_all_tests(sample_dir, verbose)
    sys.exit(0 if success else 1)
