"""
Subset Selector - Select a diverse subset of samples for analysis.

Uses a diversity-first strategy:
1. First, pick one representative per species (maximize diversity)
2. Then, randomly fill remaining slots from the pool

Integrates with the sample registry to automatically exclude previously
picked samples and register new picks.
"""

import csv
import random
import sys
from pathlib import Path
from typing import Optional, Set

from .registry import load_registry_ids, append_to_registry

# Increase CSV field size limit
csv.field_size_limit(sys.maxsize)


def load_exclude_file(exclude_file: str) -> set:
    """Load sample IDs from a plain-text exclude file (one ID per line)."""
    ids = set()
    with open(exclude_file, 'r') as f:
        for line in f:
            sid = line.strip()
            if sid:
                ids.add(sid)
    return ids


def select_subset(
    master_list_file: str,
    output_file: str = "target_samples.tsv",
    target_size: int = 10000,
    seed: Optional[int] = None,
    exclude_ids: Optional[Set[str]] = None,
    registry_file: Optional[str] = None,
    batch_name: Optional[str] = None,
    no_register: bool = False,
) -> int:
    """
    Select a diverse subset of samples from the master list.

    Strategy:
    1. Exclude previously picked samples (from registry + exclude_ids)
    2. Pick 1 representative per unique species (Round 1 - diversity)
    3. If under target, randomly fill from remaining samples
    4. If over target (unlikely), randomly subsample
    5. Optionally register new picks in the registry

    Args:
        master_list_file: Path to master list TSV
        output_file: Path for output subset TSV
        target_size: Target number of samples to select
        seed: Random seed for reproducibility (optional)
        exclude_ids: Set of sample IDs to exclude (optional)
        registry_file: Path to sample registry TSV (optional).
                       If provided, auto-excludes registered samples.
        batch_name: Label for this batch in the registry (e.g., "batch3").
                    Required if registry_file is set and no_register is False.
        no_register: If True, don't append new picks to the registry.

    Returns:
        Number of samples selected
    """
    if seed is not None:
        random.seed(seed)

    # Build combined exclusion set
    all_excluded = set()
    if exclude_ids:
        all_excluded.update(exclude_ids)

    if registry_file:
        registry_ids = load_registry_ids(registry_file)
        all_excluded.update(registry_ids)
        print(f"Registry: {len(registry_ids):,} samples already picked")

    if all_excluded:
        print(f"Total exclusions: {len(all_excluded):,} sample IDs")

    print("1. Reading master list and grouping by species...")
    species_map = {}  # { "Escherichia coli": [row1, row2...], ... }
    skipped_excluded = 0

    with open(master_list_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        headers = reader.fieldnames

        for row in reader:
            sid = row['Sample_ID']
            if sid in all_excluded:
                skipped_excluded += 1
                continue

            name = row['Scientific_Name']
            if name not in species_map:
                species_map[name] = []
            species_map[name].append(row)

    total_available = sum(len(rows) for rows in species_map.values())
    unique_species_count = len(species_map)
    print(f"   Excluded {skipped_excluded:,} previously picked samples")
    print(f"   Available: {total_available:,} samples across {unique_species_count:,} species")

    if total_available == 0:
        print("ERROR: No samples available after exclusions!")
        return 0

    if total_available < target_size:
        print(f"   Warning: Only {total_available:,} available, less than target {target_size:,}")

    selected_samples = []

    # Round 1: Pick 1 representative per species
    print("2. Picking 1 representative per species...")
    remaining_pool = []

    for name, rows in species_map.items():
        chosen = random.choice(rows)
        selected_samples.append(chosen)

        rows.remove(chosen)
        remaining_pool.extend(rows)

    current_count = len(selected_samples)
    print(f"   Round 1 selected {current_count:,} samples (high diversity).")

    # Round 2: Fill to target if needed
    if current_count < target_size:
        needed = target_size - current_count
        print(f"3. Need {needed:,} more to reach {target_size:,}. Randomly filling...")

        if needed > len(remaining_pool):
            print("   Warning: Not enough samples total! Taking everything.")
            selected_samples.extend(remaining_pool)
        else:
            random_fill = random.sample(remaining_pool, needed)
            selected_samples.extend(random_fill)

    # Handle case where we have more species than target (unlikely)
    elif current_count > target_size:
        print(f"3. Too many unique species! Subsampling to {target_size:,}...")
        selected_samples = random.sample(selected_samples, target_size)

    # Write output
    print(f"4. Writing {len(selected_samples):,} samples to {output_file}...")
    with open(output_file, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=headers, delimiter='\t')
        writer.writeheader()
        writer.writerows(selected_samples)

    # Register in registry
    if registry_file and not no_register:
        if not batch_name:
            print("   Warning: --batch-name not set, skipping registry update")
        else:
            registry_samples = [
                {'Sample_ID': row['Sample_ID'],
                 'Scientific_Name': row['Scientific_Name']}
                for row in selected_samples
            ]
            n_registered = append_to_registry(
                registry_file, registry_samples, batch_name)
            print(f"5. Registered {n_registered:,} new samples in registry "
                  f"as '{batch_name}'")

    print(f"Done! Selected {len(selected_samples):,} samples.")
    return len(selected_samples)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Select diverse subset of samples")
    parser.add_argument("--input", required=True, help="Master list TSV file")
    parser.add_argument("--output", default="target_samples.tsv", help="Output file")
    parser.add_argument("--size", type=int, default=10000, help="Target sample count")
    parser.add_argument("--seed", type=int, help="Random seed for reproducibility")
    parser.add_argument("--exclude", help="File of sample IDs to exclude (one per line)")
    parser.add_argument("--registry", help="Path to sample registry TSV")
    parser.add_argument("--batch-name", help="Batch label for registry")
    parser.add_argument("--no-register", action="store_true",
                        help="Don't update the registry with new picks")
    args = parser.parse_args()

    exclude = load_exclude_file(args.exclude) if args.exclude else None

    select_subset(
        args.input, args.output, args.size, args.seed,
        exclude_ids=exclude,
        registry_file=args.registry,
        batch_name=args.batch_name,
        no_register=args.no_register,
    )
