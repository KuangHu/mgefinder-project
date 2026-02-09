"""
Sample Registry - Track which biosamples have been picked across batches.

Maintains a central TSV file so that future subset selections automatically
exclude previously picked samples. Each entry records the sample ID, which
batch it was assigned to, its species, and the date it was added.

Registry format (TSV):
    Sample_ID    Batch    Scientific_Name    Date_Added
"""

import csv
import os
import sys
from datetime import date
from pathlib import Path
from typing import Optional

csv.field_size_limit(sys.maxsize)

REGISTRY_FIELDS = ['Sample_ID', 'Batch', 'Scientific_Name', 'Date_Added']


def load_registry(registry_file: str) -> dict:
    """
    Load the sample registry and return a dict of {Sample_ID: row_dict}.

    Returns empty dict if the file doesn't exist yet.
    """
    path = Path(registry_file)
    if not path.exists():
        return {}

    registry = {}
    with open(path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            registry[row['Sample_ID']] = row
    return registry


def load_registry_ids(registry_file: str) -> set:
    """Load just the sample IDs from the registry (fast, low memory)."""
    path = Path(registry_file)
    if not path.exists():
        return set()

    ids = set()
    with open(path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            ids.add(row['Sample_ID'])
    return ids


def append_to_registry(
    registry_file: str,
    samples: list,
    batch_name: str,
    date_added: Optional[str] = None,
) -> int:
    """
    Append new sample entries to the registry.

    Args:
        registry_file: Path to registry TSV
        samples: List of dicts, each with at least 'Sample_ID' and optionally
                 'Scientific_Name'. Can also be a list of (sample_id, species) tuples.
        batch_name: Label for this batch (e.g., "batch1", "batch2")
        date_added: Date string (default: today's date)

    Returns:
        Number of samples appended
    """
    path = Path(registry_file)
    if date_added is None:
        date_added = date.today().isoformat()

    # Load existing IDs to avoid duplicates
    existing_ids = load_registry_ids(registry_file)

    write_header = not path.exists() or path.stat().st_size == 0
    count = 0

    with open(path, 'a', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=REGISTRY_FIELDS, delimiter='\t')
        if write_header:
            writer.writeheader()

        for sample in samples:
            if isinstance(sample, dict):
                sid = sample.get('Sample_ID', '')
                species = sample.get('Scientific_Name', '')
            elif isinstance(sample, (list, tuple)) and len(sample) >= 2:
                sid, species = sample[0], sample[1]
            else:
                sid = str(sample)
                species = ''

            if sid in existing_ids:
                continue

            writer.writerow({
                'Sample_ID': sid,
                'Batch': batch_name,
                'Scientific_Name': species,
                'Date_Added': date_added,
            })
            existing_ids.add(sid)
            count += 1

    return count


def register_from_file(
    registry_file: str,
    input_file: str,
    batch_name: str,
    date_added: Optional[str] = None,
) -> int:
    """
    Register samples from an input file (target TSV or plain text ID list).

    Supports two formats:
      - TSV with header containing 'Sample_ID' column (and optionally 'Scientific_Name')
      - Plain text with one sample ID per line

    Returns:
        Number of new samples registered
    """
    path = Path(input_file)
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    samples = []

    # Detect format by reading first line
    with open(path, 'r') as f:
        first_line = f.readline().strip()

    if '\t' in first_line and 'Sample_ID' in first_line:
        # TSV with header
        with open(path, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                samples.append({
                    'Sample_ID': row['Sample_ID'],
                    'Scientific_Name': row.get('Scientific_Name', ''),
                })
    else:
        # Plain text (one ID per line)
        with open(path, 'r') as f:
            for line in f:
                sid = line.strip()
                if sid:
                    samples.append({'Sample_ID': sid, 'Scientific_Name': ''})

    return append_to_registry(registry_file, samples, batch_name, date_added)


def show_status(registry_file: str):
    """Print a summary of the registry contents."""
    path = Path(registry_file)
    if not path.exists():
        print(f"Registry not found: {registry_file}")
        print("Run 'register' to create it.")
        return

    registry = load_registry(registry_file)
    total = len(registry)

    if total == 0:
        print("Registry is empty.")
        return

    # Count per batch
    batch_counts = {}
    for row in registry.values():
        batch = row.get('Batch', 'unknown')
        batch_counts[batch] = batch_counts.get(batch, 0) + 1

    # Count unique species
    species = set()
    for row in registry.values():
        name = row.get('Scientific_Name', '')
        if name:
            species.add(name)

    print(f"Registry: {registry_file}")
    print(f"Total samples: {total:,}")
    print(f"Unique species: {len(species):,}")
    print(f"Batches: {len(batch_counts)}")
    print()
    for batch in sorted(batch_counts.keys()):
        print(f"  {batch}: {batch_counts[batch]:,} samples")
