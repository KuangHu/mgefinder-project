"""
Configuration for Download Manager.

Default paths for ENA data files on the cluster.
"""

from pathlib import Path

# Base path for ENA database files
ENA_BASE_PATH = Path("/global/home/groups/pc_rubinlab/databases/kuang/ena")

# Project home directory
PROJECT_DIR = Path("/global/home/users/kh36969/mgefinder_project")

# Default file paths
DEFAULT_READS_FILE = ENA_BASE_PATH / "ena_all_bacteria_reads.tsv"
DEFAULT_ASSEMBLY_FILE = ENA_BASE_PATH / "ena_all_bacteria_assembly.tsv"
DEFAULT_MASTER_LIST = ENA_BASE_PATH / "master_list.tsv"
DEFAULT_TARGET_FILE = ENA_BASE_PATH / "target_10k_uniform.tsv"
DEFAULT_REGISTRY_FILE = PROJECT_DIR / "sample_registry.tsv"

# MGEfinder directory structure
REFERENCE_DIR = "00.reference"
READS_DIR = "00.reads"
GENOME_FILENAME = "genome.fna"
