#!/bin/bash
#SBATCH --job-name=is_extract
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --exclusive

# ============================================================================
# Full IS Element Extraction + Prodigal Annotation
# Uses 30 parallel workers on a full 32-CPU node
#
# Usage:
#   bash run_is_extraction.sh /path/to/mgefinder_batch
#
# Submits a single-node SLURM job that runs the Python extraction script
# with multiprocessing (not an array job — IS extraction is lightweight).
# ============================================================================

set -euo pipefail

BASE_DIR="${1:-}"

if [[ -z "$BASE_DIR" ]]; then
    echo "Usage: bash $0 /path/to/mgefinder_batch"
    echo ""
    echo "Submits a SLURM job to extract and annotate IS elements."
    echo "The Python script auto-discovers incomplete samples unless"
    echo "a sample list is found at {BASE_DIR}/samples_with_is.txt."
    echo ""
    echo "Example:"
    echo "  bash $0 /global/scratch/users/kh36969/mgefinder_batch1"
    exit 1
fi

BASE_DIR="$(realpath "$BASE_DIR")"
SAMPLES_DIR="${BASE_DIR}/samples"
SCRIPTS_DIR="${BASE_DIR}/scripts"
LOG_DIR="${BASE_DIR}/jobs"

mkdir -p "$LOG_DIR"

# If running under SLURM already, just run the Python script directly
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    echo "Starting IS extraction + annotation at $(date)"
    echo "Host: $(hostname)"
    echo "CPUs available: $(nproc)"
    echo "Base directory: ${BASE_DIR}"

    # Load Prodigal for ORF prediction
    module load bio/prodigal/2.6.3-gcc-11.4.0
    echo "Prodigal version: $(prodigal -v 2>&1)"

    python3 "${SCRIPTS_DIR}/run_is_extraction.py" \
        --samples-dir "$SAMPLES_DIR"

    echo "Finished at $(date)"
    exit 0
fi

# Launcher mode: submit as SLURM job
echo "Submitting IS extraction job..."
echo "  Base directory: ${BASE_DIR}"
echo "  Samples directory: ${SAMPLES_DIR}"

SCRIPT_PATH="$(realpath "$0")"
JOB_ID=$(sbatch \
    --output="${LOG_DIR}/is_extract_%j.out" \
    --error="${LOG_DIR}/is_extract_%j.err" \
    --parsable \
    "$SCRIPT_PATH" "$BASE_DIR")

echo "  Job ID: ${JOB_ID}"
echo ""
echo "Monitor with:"
echo "  squeue -j ${JOB_ID}"
echo "  tail -f ${LOG_DIR}/is_extract_${JOB_ID}.out"
