"""
Download Manager Module for MGEfinder Project

This module handles downloading bacterial genome data from ENA (European Nucleotide Archive).

Workflow:
1. match_samples() - Find samples with both reads and assembly data
2. build_master_list() - Merge metadata into a unified master list
3. select_subset() - Select a diverse subset for analysis (auto-excludes previous batches)
4. download_samples() - Download the actual files

Sample registry tracks all picked samples across batches to prevent duplicates.
"""

from .matcher import match_samples
from .master_list import build_master_list
from .subset import select_subset
from .downloader import download_samples
from .registry import load_registry, append_to_registry, register_from_file

__all__ = [
    'match_samples', 'build_master_list', 'select_subset', 'download_samples',
    'load_registry', 'append_to_registry', 'register_from_file',
]
