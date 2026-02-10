#!/bin/bash -l
#SBATCH --job-name=is_circle
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --exclusive

#===============================================================================
# detect_circularized_is.sh
#
# Detect circularized IS elements using Circle-Map.
# Each SLURM array task processes a chunk of 32 samples on a full node,
# running them in parallel (1 CPU per sample).
#
# Usage:
#   # Launcher mode - discover samples and submit jobs
#   bash detect_circularized_is.sh /path/to/mgefinder_batch [--dry-run] [--max-jobs N]
#
#   # Force reprocessing of all eligible samples
#   bash detect_circularized_is.sh /path/to/mgefinder_batch --force
#===============================================================================

set -euo pipefail

SAMPLES_PER_NODE=32

# ─────────────────────────────────────────────────────────────────────────────
# SINGLE SAMPLE PROCESSING FUNCTION
# ─────────────────────────────────────────────────────────────────────────────
process_sample() {
    local SAMPLE="$1"
    local SAMPLES_DIR="$2"
    local SCRIPTS_DIR="$3"
    local THREADS=1

    local SAMPLE_DIR="$SAMPLES_DIR/$SAMPLE"
    local WORK_DIR="$SAMPLE_DIR/is_circle_analysis"

    echo "[${SAMPLE}] Starting at $(date '+%H:%M:%S')"

    # Create work directory
    mkdir -p "$WORK_DIR"
    cd "$WORK_DIR"

    # ── Check Prerequisites ──
    local IS_FASTA="$SAMPLE_DIR/03.results/genome/04.makefasta.genome.all_seqs.fna"
    local NO_MGE_FILE="$SAMPLE_DIR/03.results/genome/NO_MGE_DETECTED.txt"

    if [[ -f "$NO_MGE_FILE" ]]; then
        echo "[${SAMPLE}] No MGEs detected, skipping"
        echo "SKIPPED:NO_MGE" > status.txt
        return 0
    fi

    if [[ ! -f "$IS_FASTA" ]] || [[ ! -s "$IS_FASTA" ]]; then
        echo "[${SAMPLE}] No IS sequences found, skipping"
        echo "SKIPPED:NO_IS_FILE" > status.txt
        return 0
    fi

    local READS_R1="$SAMPLE_DIR/00.reads/${SAMPLE}_R1.fastq.gz"
    local READS_R2="$SAMPLE_DIR/00.reads/${SAMPLE}_R2.fastq.gz"

    if [[ ! -f "$READS_R1" ]] || [[ ! -f "$READS_R2" ]]; then
        echo "[${SAMPLE}] FASTQ files not found, skipping"
        echo "SKIPPED:NO_FASTQ" > status.txt
        return 0
    fi

    # ── Step 1: Extract IS Reference ──
    python "$SCRIPTS_DIR/extract_is_reference.py" \
        --sample-dir "$SAMPLE_DIR" \
        --output is_reference.fna \
        --min-length 500 \
        --max-length 5000 2>/dev/null

    if [[ ! -s is_reference.fna ]]; then
        echo "[${SAMPLE}] No valid IS sequences (filtered by length)"
        echo "SKIPPED:NO_VALID_IS" > status.txt
        return 0
    fi

    # ── Step 2: Align Reads to IS Reference ──
    bwa index is_reference.fna 2>/dev/null

    bwa mem -t $THREADS is_reference.fna "$READS_R1" "$READS_R2" 2>bwa.log | \
        samtools sort -@ $THREADS -o is_aligned.bam -

    samtools index is_aligned.bam

    local MAPPED_READS
    MAPPED_READS=$(samtools view -c -F 4 is_aligned.bam)

    if [[ $MAPPED_READS -lt 100 ]]; then
        echo "[${SAMPLE}] Too few reads mapped ($MAPPED_READS), skipping"
        echo "SKIPPED:LOW_MAPPING" > status.txt
        rm -f *.bam *.bam.bai
        return 0
    fi

    # ── Step 3: Run Circle-Map ──
    samtools sort -n -@ $THREADS -o is_aligned.qname.bam is_aligned.bam

    Circle-Map ReadExtractor \
        -i is_aligned.qname.bam \
        -o circular_candidates.bam 2>readextractor.log

    local CANDIDATES
    CANDIDATES=$(samtools view -c circular_candidates.bam 2>/dev/null || echo "0")

    if [[ "$CANDIDATES" -lt 10 ]]; then
        echo "[${SAMPLE}] No circular candidates ($CANDIDATES reads)"
        echo "COMPLETED:0:NO_CANDIDATES" > status.txt
        echo -e "is_id\tis_length\tcircle_start\tcircle_end\tcircle_length\tlength_ratio\tsplit_reads\tdiscordant_reads\tscore\tvalidation" > circularized_is.tsv
        rm -f *.bam *.bam.bai
        return 0
    fi

    samtools sort -@ $THREADS -o circular_candidates.sorted.bam circular_candidates.bam
    samtools index circular_candidates.sorted.bam

    Circle-Map Realign -t $THREADS \
        -i circular_candidates.sorted.bam \
        -qbam is_aligned.qname.bam \
        -sbam is_aligned.bam \
        -fasta is_reference.fna \
        -o is_circles.bed 2>realign.log

    if [[ ! -s is_circles.bed ]]; then
        echo "[${SAMPLE}] No circles detected"
        echo "COMPLETED:0:NO_CIRCLES" > status.txt
        echo -e "is_id\tis_length\tcircle_start\tcircle_end\tcircle_length\tlength_ratio\tsplit_reads\tdiscordant_reads\tscore\tvalidation" > circularized_is.tsv
        rm -f *.bam *.bam.bai
        return 0
    fi

    # ── Step 4: Validate Results ──
    python "$SCRIPTS_DIR/validate_is_circles.py" \
        --circles is_circles.bed \
        --reference is_reference.fna \
        --output circularized_is.tsv \
        --tolerance 0.1 \
        --min-reads 2 2>/dev/null

    # ── Summary ──
    if [[ -s circularized_is.tsv ]]; then
        local CONFIRMED LIKELY TOTAL
        CONFIRMED=$(grep -c "CONFIRMED_CIRCULAR_IS" circularized_is.tsv || echo "0")
        LIKELY=$(grep -c "LIKELY_CIRCULAR_IS" circularized_is.tsv || echo "0")
        TOTAL=$((CONFIRMED + LIKELY))
        echo "[${SAMPLE}] Done: ${CONFIRMED} confirmed, ${LIKELY} likely circular IS"
        echo "COMPLETED:$TOTAL:CONFIRMED=$CONFIRMED,LIKELY=$LIKELY" > status.txt
    else
        echo "[${SAMPLE}] Done: no circularized IS"
        echo "COMPLETED:0:NO_CIRCLES" > status.txt
    fi

    # ── Cleanup ──
    rm -f *.bam *.bam.bai is_reference.fna.* 2>/dev/null
    return 0
}

export -f process_sample

# ─────────────────────────────────────────────────────────────────────────────
# MODE DETECTION: launcher vs. SLURM worker
# ─────────────────────────────────────────────────────────────────────────────

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    # =========================================================================
    # WORKER MODE: process a chunk of samples on a full node
    # =========================================================================

    BASE_DIR="$1"
    BATCH_NUM="${2:-1}"
    SAMPLES_DIR="${BASE_DIR}/samples"
    SCRIPTS_DIR="${BASE_DIR}/scripts"
    SAMPLE_LIST="${SCRIPTS_DIR}/.is_circle_queue_batch${BATCH_NUM}.txt"

    echo "=============================================="
    echo "IS Circle Detection — Node Worker"
    echo "=============================================="
    echo "Batch: ${BATCH_NUM}, Task: ${SLURM_ARRAY_TASK_ID}"
    echo "Host: $(hostname), CPUs: $(nproc)"
    echo "Start: $(date)"
    echo "=============================================="

    # Activate conda environment
    set +u
    source $(conda info --base)/etc/profile.d/conda.sh
    conda activate mgefinder
    set -u

    # Calculate sample range for this task
    CHUNK_START=$(( (SLURM_ARRAY_TASK_ID - 1) * SAMPLES_PER_NODE + 1 ))
    CHUNK_END=$(( SLURM_ARRAY_TASK_ID * SAMPLES_PER_NODE ))
    TOTAL_IN_BATCH=$(wc -l < "$SAMPLE_LIST")
    [[ $CHUNK_END -gt $TOTAL_IN_BATCH ]] && CHUNK_END=$TOTAL_IN_BATCH

    echo "Processing samples ${CHUNK_START}-${CHUNK_END} of ${TOTAL_IN_BATCH}"
    echo ""

    # Extract chunk of samples
    CHUNK_SAMPLES=$(sed -n "${CHUNK_START},${CHUNK_END}p" "$SAMPLE_LIST")
    NUM_SAMPLES=$(echo "$CHUNK_SAMPLES" | wc -l)

    # Process all samples in parallel (1 CPU each, up to SAMPLES_PER_NODE concurrent)
    RUNNING=0
    DONE=0
    FAILED=0

    for SAMPLE in $CHUNK_SAMPLES; do
        (
            process_sample "$SAMPLE" "$SAMPLES_DIR" "$SCRIPTS_DIR"
        ) &

        RUNNING=$((RUNNING + 1))

        # Wait if we've hit the parallelism limit
        if [[ $RUNNING -ge $SAMPLES_PER_NODE ]]; then
            wait -n 2>/dev/null || true
            RUNNING=$((RUNNING - 1))
        fi
    done

    # Wait for all remaining jobs
    wait

    echo ""
    echo "=============================================="
    echo "Node complete: processed $NUM_SAMPLES samples"
    echo "End: $(date)"
    echo "=============================================="
    exit 0
fi

# =========================================================================
# LAUNCHER MODE: discover samples and submit SLURM array job(s)
# =========================================================================

# --- Parse arguments ---
BASE_DIR=""
DRY_RUN=false
MAX_JOBS=50
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
    echo "Usage: bash $0 /path/to/mgefinder_batch [OPTIONS]"
    echo ""
    echo "Discovers samples ready for IS circle detection and submits SLURM array jobs."
    echo "Each array task uses a full node (32 CPUs) to process ${SAMPLES_PER_NODE} samples in parallel."
    echo ""
    echo "Options:"
    echo "  --dry-run       Show what would be submitted, don't submit"
    echo "  --max-jobs N    Max concurrent SLURM array tasks (default: 50)"
    echo "  --force         Re-process samples even if circularized_is.tsv exists"
    echo "  --batch-size N  Max samples per queue file (default: 1000)"
    echo ""
    echo "Examples:"
    echo "  bash $0 /global/scratch/users/kh36969/mgefinder_batch1 --dry-run"
    echo "  bash $0 /global/scratch/users/kh36969/mgefinder_batch2 --max-jobs 30"
    exit 1
fi

BASE_DIR="$(realpath "$BASE_DIR")"
SAMPLES_DIR="${BASE_DIR}/samples"
SCRIPTS_DIR="${BASE_DIR}/scripts"
QUEUE_FILE="${SCRIPTS_DIR}/.is_circle_queue.txt"
LOG_DIR="${BASE_DIR}/jobs/is_circle_logs"

mkdir -p "$SCRIPTS_DIR" "$LOG_DIR"

echo "============================================"
echo "IS Circle Detection — Sample Discovery"
echo "============================================"
echo "Base directory: ${BASE_DIR}"
echo "Samples per node: ${SAMPLES_PER_NODE}"
echo "Scanning for eligible samples..."
echo ""

# --- Discover samples ---
> "$QUEUE_FILE"
TOTAL=0
HAS_IS=0
HAS_FASTQ=0
ALREADY_DONE=0
NEED_PROCESSING=0

for dir in "${SAMPLES_DIR}"/SAM*/; do
    [[ ! -d "$dir" ]] && continue
    SAMPLE=$(basename "$dir")
    TOTAL=$((TOTAL + 1))

    IS_FASTA="${dir}03.results/genome/04.makefasta.genome.all_seqs.fna"
    NO_MGE="${dir}03.results/genome/NO_MGE_DETECTED.txt"

    if [[ -f "$NO_MGE" ]]; then
        continue
    fi
    if [[ ! -f "$IS_FASTA" ]] || [[ ! -s "$IS_FASTA" ]]; then
        continue
    fi
    HAS_IS=$((HAS_IS + 1))

    if [[ ! -f "${dir}00.reads/${SAMPLE}_R1.fastq.gz" ]] || \
       [[ ! -f "${dir}00.reads/${SAMPLE}_R2.fastq.gz" ]]; then
        continue
    fi
    HAS_FASTQ=$((HAS_FASTQ + 1))

    if [[ "$FORCE" == "false" ]]; then
        if [[ -f "${dir}is_circle_analysis/circularized_is.tsv" ]]; then
            ALREADY_DONE=$((ALREADY_DONE + 1))
            continue
        fi
    fi

    echo "$SAMPLE" >> "$QUEUE_FILE"
    NEED_PROCESSING=$((NEED_PROCESSING + 1))
done

echo "Results:"
echo "  Total sample dirs:      ${TOTAL}"
echo "  With IS sequences:      ${HAS_IS}"
echo "  With IS + FASTQ:        ${HAS_FASTQ}"
echo "  Already completed:      ${ALREADY_DONE}"
echo "  Need processing:        ${NEED_PROCESSING}"
echo "============================================"

if [[ "$NEED_PROCESSING" -eq 0 ]]; then
    echo "Nothing to do. All eligible samples already have circularized_is.tsv."
    exit 0
fi

echo ""
echo "Queue file: ${QUEUE_FILE}"
echo "First 10 samples in queue:"
head -10 "$QUEUE_FILE" | sed 's/^/  /'
[[ "$NEED_PROCESSING" -gt 10 ]] && echo "  ... and $((NEED_PROCESSING - 10)) more"
echo ""

# --- Split into batch files and calculate array sizes ---
NUM_BATCHES=$(( (NEED_PROCESSING + BATCH_SIZE - 1) / BATCH_SIZE ))

echo "Creating ${NUM_BATCHES} batch queue file(s)..."
for ((batch=1; batch<=NUM_BATCHES; batch++)); do
    start=$(( (batch - 1) * BATCH_SIZE + 1 ))
    end=$(( batch * BATCH_SIZE ))
    [[ $end -gt $NEED_PROCESSING ]] && end=$NEED_PROCESSING

    batch_file="${SCRIPTS_DIR}/.is_circle_queue_batch${batch}.txt"
    sed -n "${start},${end}p" "$QUEUE_FILE" > "$batch_file"

    batch_count=$(wc -l < "$batch_file")
    # Each array task processes SAMPLES_PER_NODE samples
    array_size=$(( (batch_count + SAMPLES_PER_NODE - 1) / SAMPLES_PER_NODE ))
    echo "  Batch ${batch}: ${batch_count} samples → ${array_size} array tasks (${SAMPLES_PER_NODE} samples/node)"
done
echo ""

if [[ "$DRY_RUN" == "true" ]]; then
    echo "[DRY RUN] Would submit ${NUM_BATCHES} batch(es):"
    for ((batch=1; batch<=NUM_BATCHES; batch++)); do
        batch_file="${SCRIPTS_DIR}/.is_circle_queue_batch${batch}.txt"
        batch_count=$(wc -l < "$batch_file")
        array_size=$(( (batch_count + SAMPLES_PER_NODE - 1) / SAMPLES_PER_NODE ))
        echo "  Batch ${batch}: sbatch --array=1-${array_size}%${MAX_JOBS} (${batch_count} samples, ${array_size} nodes)"
    done
    echo ""
    echo "Full sample list saved to: ${QUEUE_FILE}"
    exit 0
fi

# --- Submit SLURM array job(s) ---
SCRIPT_PATH="$(realpath "$0")"
SUBMITTED_JOBS=()

echo "Submitting ${NUM_BATCHES} batch(es)..."
echo ""

for ((batch=1; batch<=NUM_BATCHES; batch++)); do
    batch_file="${SCRIPTS_DIR}/.is_circle_queue_batch${batch}.txt"
    batch_count=$(wc -l < "$batch_file")
    array_size=$(( (batch_count + SAMPLES_PER_NODE - 1) / SAMPLES_PER_NODE ))

    JOB_ID=$(sbatch \
        --array="1-${array_size}%${MAX_JOBS}" \
        --output="${LOG_DIR}/is_circle_%A_%a.out" \
        --error="${LOG_DIR}/is_circle_%A_%a.err" \
        --parsable \
        "$SCRIPT_PATH" "$BASE_DIR" "$batch")

    SUBMITTED_JOBS+=("$JOB_ID")

    echo "Batch ${batch}/${NUM_BATCHES}:"
    echo "  Job ID: ${JOB_ID}"
    echo "  Array: 1-${array_size} (${batch_count} samples, ${SAMPLES_PER_NODE}/node, max ${MAX_JOBS} concurrent nodes)"
    echo ""

    sleep 0.5
done

echo "============================================"
echo "Summary"
echo "============================================"
echo "Total samples to process: ${NEED_PROCESSING}"
echo "Batches submitted: ${NUM_BATCHES}"
echo "Job IDs: ${SUBMITTED_JOBS[*]}"
echo ""
echo "Monitor with:"
echo "  squeue -u $USER"
echo "  squeue -j ${SUBMITTED_JOBS[0]}"
echo ""
echo "Check progress:"
echo "  ls -1 ${SAMPLES_DIR}/SAM*/is_circle_analysis/circularized_is.tsv 2>/dev/null | wc -l"
echo "============================================"
