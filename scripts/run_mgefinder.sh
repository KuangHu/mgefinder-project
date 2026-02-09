#!/bin/bash -l
#SBATCH --job-name=mgefinder
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

#===============================================================================
# run_mgefinder.sh
#
# Run MGEfinder de novo workflow on prepared samples.
# Supports SLURM array jobs for parallel processing.
#
# Usage:
#   # Launcher mode - discover samples and submit jobs
#   bash run_mgefinder.sh /path/to/samples [--dry-run] [--max-jobs N]
#
#   # Direct run on single sample
#   bash run_mgefinder.sh /path/to/samples --sample SAMEA2273967
#===============================================================================

set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# MODE DETECTION: launcher vs. SLURM worker
# ─────────────────────────────────────────────────────────────────────────────

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    # =========================================================================
    # WORKER MODE: process a single sample
    # =========================================================================

    BASE_DIR="$1"
    BATCH_NUM="${2:-1}"
    SAMPLE_LIST="${BASE_DIR}/scripts/.mgefinder_queue_batch${BATCH_NUM}.txt"
    CONDA_ENV="mgefinder"
    THREADS="${SLURM_CPUS_PER_TASK:-8}"
    MEMORY_MB="${SLURM_MEM_PER_NODE:-32000}"

    SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
    SAMPLE_DIR="${BASE_DIR}/${SAMPLE}"

    echo "============================================"
    echo "MGEfinder De Novo Workflow"
    echo "============================================"
    echo "Date: $(date)"
    echo "Batch: ${BATCH_NUM}, Job ID: ${SLURM_JOB_ID}, Task: ${SLURM_ARRAY_TASK_ID}"
    echo "Sample: ${SAMPLE}"
    echo "Directory: ${SAMPLE_DIR}"
    echo "Threads: ${THREADS}"
    echo "Memory: ${MEMORY_MB} MB"
    echo "============================================"

    if [[ ! -d "$SAMPLE_DIR" ]]; then
        echo "ERROR: Sample directory does not exist: $SAMPLE_DIR"
        exit 1
    fi

    # --- Validate required inputs ---
    GENOME_FILE="${SAMPLE_DIR}/00.genome/genome.fna"
    ASSEMBLY_FILE="${SAMPLE_DIR}/00.assembly/${SAMPLE}.fna"
    BAM_FILE="${SAMPLE_DIR}/00.bam/${SAMPLE}.genome.bam"

    echo "Checking required files..."
    for f in "$GENOME_FILE" "$ASSEMBLY_FILE" "$BAM_FILE"; do
        if [[ ! -s "$f" ]]; then
            echo "ERROR: Missing or empty: $f"
            exit 1
        fi
        echo "  OK: $(basename "$f") ($(ls -lh "$f" | awk '{print $5}'))"
    done
    echo ""

    # --- Activate conda ---
    if [[ -f "$(conda info --base 2>/dev/null)/etc/profile.d/conda.sh" ]]; then
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate "$CONDA_ENV"
    else
        export PATH="$HOME/.conda/envs/${CONDA_ENV}/bin:$PATH"
    fi

    echo "Environment: $CONDA_ENV"
    echo "MGEfinder: $(mgefinder --version 2>&1 | grep -v CHECKING | grep -v Current | grep -v Expected | grep -v '###' | head -1 || echo 'installed')"
    echo "============================================"

    # --- Run MGEfinder ---
    echo "Starting MGEfinder de novo workflow..."
    echo ""

    cd "$SAMPLE_DIR"

    mgefinder workflow denovo \
        --cores "$THREADS" \
        --memory "$MEMORY_MB" \
        --keep-going \
        --rerun-incomplete \
        --sensitive \
        . \
        2>&1 | tee mgefinder.log

    EXIT_CODE=${PIPESTATUS[0]}

    echo ""
    echo "============================================"
    echo "Results"
    echo "============================================"

    if [[ $EXIT_CODE -eq 0 ]]; then
        echo "MGEfinder completed successfully!"

        # List output files
        echo ""
        echo "Output files:"
        find . -maxdepth 2 -name "*.fna" -o -name "*.tsv" -o -name "*.csv" 2>/dev/null | \
            grep -v "00\." | head -20 | while read f; do
            echo "  $f ($(ls -lh "$f" 2>/dev/null | awk '{print $5}'))"
        done

        echo ""
        echo "SUCCESS: ${SAMPLE}"
    else
        # Check if this is a "no termini found" case (not a real error)
        if grep -q "No termini found in the input file" mgefinder.log 2>/dev/null; then
            echo ""
            echo "NOTE: No mobile genetic elements (MGEs) detected in this sample."
            echo "This is a valid biological result, not an error."
            echo ""

            # Create empty placeholder results so downstream processing works
            mkdir -p 03.results/genome

            # Create empty results files with headers only
            echo -e "cluster_id\tseq_id\tmethod\tloc" > 03.results/genome/01.clusterseq.genome.tsv
            echo -e "sample\tcluster_id\tpair_id\tmethod\tcontig\tpos_5p\tpos_3p\tgenotype" > 03.results/genome/02.genotype.genome.tsv
            echo -e "cluster_id\tnum_unique_sites\tnum_samples\tsamples" > 03.results/genome/03.summarize.genome.clusters.tsv
            echo -e "group_id\tcluster_id\tmethod" > 03.results/genome/03.summarize.genome.groups.tsv
            touch 03.results/genome/04.makefasta.genome.all_seqs.fna
            touch 03.results/genome/04.makefasta.genome.repr_seqs.fna

            # Create a marker file indicating no MGEs found
            echo "NO_MGE_DETECTED" > 03.results/genome/NO_MGE_DETECTED.txt
            echo "Sample: ${SAMPLE}" >> 03.results/genome/NO_MGE_DETECTED.txt
            echo "Date: $(date)" >> 03.results/genome/NO_MGE_DETECTED.txt
            echo "Reason: No termini/insertion sites found in the input data" >> 03.results/genome/NO_MGE_DETECTED.txt

            echo "Created empty placeholder result files."
            echo ""
            echo "SUCCESS (NO_MGE): ${SAMPLE}"
        else
            echo "MGEfinder failed with exit code: $EXIT_CODE"
            echo "Check log: ${SAMPLE_DIR}/mgefinder.log"
            echo ""
            echo "FAILED: ${SAMPLE}"
            exit 1
        fi
    fi

    echo "Completed: $(date)"
    echo "============================================"
    exit 0
fi

# ─────────────────────────────────────────────────────────────────────────────
# Check for direct single-sample mode
# ─────────────────────────────────────────────────────────────────────────────

for arg in "$@"; do
    if [[ "$arg" == "--sample" ]]; then
        # Direct single sample mode
        BASE_DIR=""
        SAMPLE=""
        THREADS=8
        MEMORY_MB=32000

        while [[ $# -gt 0 ]]; do
            case "$1" in
                --sample) SAMPLE="$2"; shift 2 ;;
                --threads) THREADS="$2"; shift 2 ;;
                --memory) MEMORY_MB="$2"; shift 2 ;;
                *) BASE_DIR="$1"; shift ;;
            esac
        done

        if [[ -z "$BASE_DIR" ]] || [[ -z "$SAMPLE" ]]; then
            echo "Usage: bash $0 /path/to/samples --sample SAMPLE_ID [--threads N] [--memory MB]"
            exit 1
        fi

        SAMPLE_DIR="${BASE_DIR}/${SAMPLE}"
        CONDA_ENV="mgefinder"

        echo "============================================"
        echo "MGEfinder De Novo Workflow (Direct Mode)"
        echo "============================================"
        echo "Sample: ${SAMPLE}"
        echo "Directory: ${SAMPLE_DIR}"
        echo "Threads: ${THREADS}"
        echo "Memory: ${MEMORY_MB} MB"
        echo "============================================"

        # Activate conda
        if [[ -f "$(conda info --base 2>/dev/null)/etc/profile.d/conda.sh" ]]; then
            source "$(conda info --base)/etc/profile.d/conda.sh"
            conda activate "$CONDA_ENV"
        else
            export PATH="$HOME/.conda/envs/${CONDA_ENV}/bin:$PATH"
        fi

        cd "$SAMPLE_DIR"

        mgefinder workflow denovo \
            --cores "$THREADS" \
            --memory "$MEMORY_MB" \
            --keep-going \
            --rerun-incomplete \
            --sensitive \
            . \
            2>&1 | tee mgefinder.log

        EXIT_CODE=${PIPESTATUS[0]}

        if [[ $EXIT_CODE -eq 0 ]]; then
            echo ""
            echo "SUCCESS: ${SAMPLE}"
            exit 0
        else
            # Check if this is a "no termini found" case (not a real error)
            if grep -q "No termini found in the input file" mgefinder.log 2>/dev/null; then
                echo ""
                echo "NOTE: No mobile genetic elements (MGEs) detected in this sample."
                echo "This is a valid biological result, not an error."
                echo ""

                # Create empty placeholder results so downstream processing works
                mkdir -p 03.results/genome

                # Create empty results files with headers only
                echo -e "cluster_id\tseq_id\tmethod\tloc" > 03.results/genome/01.clusterseq.genome.tsv
                echo -e "sample\tcluster_id\tpair_id\tmethod\tcontig\tpos_5p\tpos_3p\tgenotype" > 03.results/genome/02.genotype.genome.tsv
                echo -e "cluster_id\tnum_unique_sites\tnum_samples\tsamples" > 03.results/genome/03.summarize.genome.clusters.tsv
                echo -e "group_id\tcluster_id\tmethod" > 03.results/genome/03.summarize.genome.groups.tsv
                touch 03.results/genome/04.makefasta.genome.all_seqs.fna
                touch 03.results/genome/04.makefasta.genome.repr_seqs.fna

                # Create a marker file indicating no MGEs found
                echo "NO_MGE_DETECTED" > 03.results/genome/NO_MGE_DETECTED.txt
                echo "Sample: ${SAMPLE}" >> 03.results/genome/NO_MGE_DETECTED.txt
                echo "Date: $(date)" >> 03.results/genome/NO_MGE_DETECTED.txt
                echo "Reason: No termini/insertion sites found in the input data" >> 03.results/genome/NO_MGE_DETECTED.txt

                echo "Created empty placeholder result files."
                echo ""
                echo "SUCCESS (NO_MGE): ${SAMPLE}"
                exit 0
            else
                echo "MGEfinder failed with exit code: $EXIT_CODE"
                echo "Check log: ${SAMPLE_DIR}/mgefinder.log"
                echo ""
                echo "FAILED: ${SAMPLE}"
                exit 1
            fi
        fi
    fi
done

# =========================================================================
# LAUNCHER MODE: discover samples and submit SLURM array job(s)
# =========================================================================

BASE_DIR=""
DRY_RUN=false
MAX_JOBS=10
FORCE=false
BATCH_SIZE=1000

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY_RUN=true; shift ;;
        --max-jobs) MAX_JOBS="$2"; shift 2 ;;
        --force) FORCE=true; shift ;;
        --batch-size) BATCH_SIZE="$2"; shift 2 ;;
        *) BASE_DIR="$1"; shift ;;
    esac
done

if [[ -z "$BASE_DIR" ]]; then
    echo "Usage: bash $0 /path/to/samples [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  --dry-run       Show what would be submitted, don't submit"
    echo "  --max-jobs N    Max concurrent SLURM array tasks (default: 10)"
    echo "  --force         Re-run MGEfinder even if outputs exist"
    echo "  --batch-size N  Max array size per job (default: 1000)"
    echo ""
    echo "Single sample mode:"
    echo "  bash $0 /path/to/samples --sample SAMPLE_ID [--threads N] [--memory MB]"
    exit 1
fi

BASE_DIR="$(realpath "$BASE_DIR")"
QUEUE_FILE="${BASE_DIR}/scripts/.mgefinder_queue.txt"
mkdir -p "${BASE_DIR}/scripts"

echo "============================================"
echo "MGEfinder Sample Discovery"
echo "============================================"
echo "Base directory: ${BASE_DIR}"
echo "Scanning for prepared samples..."
echo ""

# --- Discover samples ready for MGEfinder ---
> "$QUEUE_FILE"
TOTAL=0
READY=0
NEED_PROCESSING=0
ALREADY_DONE=0

for dir in "${BASE_DIR}"/SAM*/; do
    [[ ! -d "$dir" ]] && continue
    SAMPLE=$(basename "$dir")
    TOTAL=$((TOTAL + 1))

    # Check if sample is prepared (has genome, assembly, bam)
    GENOME="${dir}00.genome/genome.fna"
    ASSEMBLY="${dir}00.assembly/${SAMPLE}.fna"
    BAM="${dir}00.bam/${SAMPLE}.genome.bam"

    if [[ ! -s "$GENOME" ]] || [[ ! -s "$ASSEMBLY" ]] || [[ ! -s "$BAM" ]]; then
        continue
    fi

    READY=$((READY + 1))

    # Check if MGEfinder already ran (look for output files)
    if [[ "$FORCE" == "false" ]]; then
        # Check for MGEfinder output directories/files
        if [[ -d "${dir}01.find" ]] || [[ -f "${dir}mgefinder.log" ]]; then
            # Check if it completed successfully
            if [[ -f "${dir}mgefinder.log" ]] && grep -q "Complete" "${dir}mgefinder.log" 2>/dev/null; then
                ALREADY_DONE=$((ALREADY_DONE + 1))
                continue
            fi
        fi
    fi

    echo "$SAMPLE" >> "$QUEUE_FILE"
    NEED_PROCESSING=$((NEED_PROCESSING + 1))
done

echo "Results:"
echo "  Total sample dirs:    ${TOTAL}"
echo "  Prepared (ready):     ${READY}"
echo "  Already processed:    ${ALREADY_DONE}"
echo "  Need MGEfinder:       ${NEED_PROCESSING}"
echo "============================================"

if [[ "$NEED_PROCESSING" -eq 0 ]]; then
    echo "Nothing to do."
    exit 0
fi

echo ""
echo "Samples to process:"
cat "$QUEUE_FILE" | sed 's/^/  /'
echo ""

# --- Calculate batches ---
NUM_BATCHES=$(( (NEED_PROCESSING + BATCH_SIZE - 1) / BATCH_SIZE ))

echo "Creating ${NUM_BATCHES} batch queue file(s)..."
for ((batch=1; batch<=NUM_BATCHES; batch++)); do
    start=$(( (batch - 1) * BATCH_SIZE + 1 ))
    end=$(( batch * BATCH_SIZE ))
    [[ $end -gt $NEED_PROCESSING ]] && end=$NEED_PROCESSING

    batch_file="${BASE_DIR}/scripts/.mgefinder_queue_batch${batch}.txt"
    sed -n "${start},${end}p" "$QUEUE_FILE" > "$batch_file"

    batch_size=$(wc -l < "$batch_file")
    echo "  Batch ${batch}: ${batch_size} samples"
done
echo ""

if [[ "$DRY_RUN" == "true" ]]; then
    echo "[DRY RUN] Would submit:"
    for ((batch=1; batch<=NUM_BATCHES; batch++)); do
        batch_file="${BASE_DIR}/scripts/.mgefinder_queue_batch${batch}.txt"
        batch_size=$(wc -l < "$batch_file")
        echo "  sbatch --array=1-${batch_size}%${MAX_JOBS} $0 ${BASE_DIR} ${batch}"
    done
    exit 0
fi

# --- Submit SLURM jobs ---
SCRIPT_PATH="$(realpath "$0")"
SUBMITTED_JOBS=()

echo "Submitting ${NUM_BATCHES} batch(es)..."
for ((batch=1; batch<=NUM_BATCHES; batch++)); do
    batch_file="${BASE_DIR}/scripts/.mgefinder_queue_batch${batch}.txt"
    batch_size=$(wc -l < "$batch_file")

    JOB_ID=$(sbatch \
        --array="1-${batch_size}%${MAX_JOBS}" \
        --parsable \
        "$SCRIPT_PATH" "$BASE_DIR" "$batch")

    SUBMITTED_JOBS+=("$JOB_ID")
    echo "  Batch ${batch}: Job ${JOB_ID} (${batch_size} tasks, max ${MAX_JOBS} concurrent)"
done

echo ""
echo "============================================"
echo "Submitted ${#SUBMITTED_JOBS[@]} job(s): ${SUBMITTED_JOBS[*]}"
echo ""
echo "Monitor:"
echo "  squeue -u $USER"
echo "  sacct -j ${SUBMITTED_JOBS[0]} --format=JobID,State,Elapsed,MaxRSS"
echo "============================================"
