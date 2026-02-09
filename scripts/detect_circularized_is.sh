#!/bin/bash -l
#SBATCH --job-name=is_circle
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

#===============================================================================
# detect_circularized_is.sh
#
# Detect circularized IS elements using Circle-Map.
# Supports SLURM array jobs with automatic sample discovery and batching.
#
# Usage:
#   # Launcher mode - discover samples and submit jobs
#   bash detect_circularized_is.sh /path/to/mgefinder_batch [--dry-run] [--max-jobs N]
#
#   # Force reprocessing of all eligible samples
#   bash detect_circularized_is.sh /path/to/mgefinder_batch --force
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
    SAMPLES_DIR="${BASE_DIR}/samples"
    SCRIPTS_DIR="${BASE_DIR}/scripts"
    SAMPLE_LIST="${SCRIPTS_DIR}/.is_circle_queue_batch${BATCH_NUM}.txt"
    THREADS="${SLURM_CPUS_PER_TASK:-8}"

    SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

    if [[ -z "$SAMPLE" ]]; then
        echo "ERROR: No sample found at line $SLURM_ARRAY_TASK_ID"
        exit 1
    fi

    SAMPLE_DIR="$SAMPLES_DIR/$SAMPLE"
    WORK_DIR="$SAMPLE_DIR/is_circle_analysis"

    echo "=============================================="
    echo "Detecting Circularized IS Elements"
    echo "=============================================="
    echo "Sample: $SAMPLE"
    echo "Batch: ${BATCH_NUM}, Job ID: ${SLURM_JOB_ID}, Task: ${SLURM_ARRAY_TASK_ID}"
    echo "Start: $(date)"
    echo "=============================================="

    # Activate conda environment (temporarily disable -u for conda scripts with unbound variables)
    set +u
    source $(conda info --base)/etc/profile.d/conda.sh
    conda activate mgefinder
    set -u

    # Create work directory
    mkdir -p "$WORK_DIR"
    cd "$WORK_DIR"

    # ============================================
    # Check Prerequisites
    # ============================================

    # Check for IS sequences
    IS_FASTA="$SAMPLE_DIR/03.results/genome/04.makefasta.genome.all_seqs.fna"
    NO_MGE_FILE="$SAMPLE_DIR/03.results/genome/NO_MGE_DETECTED.txt"

    if [[ -f "$NO_MGE_FILE" ]]; then
        echo "No MGEs detected for this sample"
        echo "SKIPPED:NO_MGE" > status.txt
        exit 0
    fi

    if [[ ! -f "$IS_FASTA" ]] || [[ ! -s "$IS_FASTA" ]]; then
        echo "No IS sequences found for $SAMPLE"
        echo "SKIPPED:NO_IS_FILE" > status.txt
        exit 0
    fi

    # Check for FASTQ files
    READS_R1="$SAMPLE_DIR/00.reads/${SAMPLE}_R1.fastq.gz"
    READS_R2="$SAMPLE_DIR/00.reads/${SAMPLE}_R2.fastq.gz"

    if [[ ! -f "$READS_R1" ]] || [[ ! -f "$READS_R2" ]]; then
        echo "FASTQ files not found (may have been cleaned up)"
        echo "SKIPPED:NO_FASTQ" > status.txt
        exit 0
    fi

    # ============================================
    # Step 1: Extract IS Reference
    # ============================================
    echo ""
    echo "Step 1: Extracting IS reference sequences..."

    python "$SCRIPTS_DIR/extract_is_reference.py" \
        --sample-dir "$SAMPLE_DIR" \
        --output is_reference.fna \
        --min-length 500 \
        --max-length 5000

    if [[ ! -s is_reference.fna ]]; then
        echo "No valid IS sequences extracted (all filtered by length)"
        echo "SKIPPED:NO_VALID_IS" > status.txt
        exit 0
    fi

    IS_COUNT=$(grep -c "^>" is_reference.fna)
    echo "Extracted $IS_COUNT IS sequences"

    # ============================================
    # Step 2: Align Reads to IS Reference
    # ============================================
    echo ""
    echo "Step 2: Aligning reads to IS-only reference..."

    # Index reference
    bwa index is_reference.fna 2>/dev/null

    # Align
    bwa mem -t $THREADS is_reference.fna "$READS_R1" "$READS_R2" 2>bwa.log | \
        samtools sort -@ $THREADS -o is_aligned.bam -

    samtools index is_aligned.bam

    # Check alignment rate
    TOTAL_READS=$(samtools view -c is_aligned.bam)
    MAPPED_READS=$(samtools view -c -F 4 is_aligned.bam)
    echo "Total reads: $TOTAL_READS"
    echo "Mapped to IS: $MAPPED_READS"

    if [[ $MAPPED_READS -lt 100 ]]; then
        echo "Too few reads mapped to IS reference"
        echo "SKIPPED:LOW_MAPPING" > status.txt
        rm -f *.bam *.bam.bai
        exit 0
    fi

    # ============================================
    # Step 3: Run Circle-Map
    # ============================================
    echo ""
    echo "Step 3: Running Circle-Map on IS reference..."

    # Sort by read name for Circle-Map
    samtools sort -n -@ $THREADS -o is_aligned.qname.bam is_aligned.bam

    # Extract circular candidate reads
    echo "  Extracting circular candidate reads..."
    Circle-Map ReadExtractor \
        -i is_aligned.qname.bam \
        -o circular_candidates.bam 2>readextractor.log

    # Check if any candidates found
    CANDIDATES=$(samtools view -c circular_candidates.bam 2>/dev/null || echo "0")
    echo "  Circular candidate reads: $CANDIDATES"

    if [[ "$CANDIDATES" -lt 10 ]]; then
        echo "  No circular candidates found"
        echo "COMPLETED:0:NO_CANDIDATES" > status.txt
        # Still create empty output file
        echo -e "is_id\tis_length\tcircle_start\tcircle_end\tcircle_length\tlength_ratio\tsplit_reads\tdiscordant_reads\tscore\tvalidation" > circularized_is.tsv
        rm -f *.bam *.bam.bai
        exit 0
    fi

    # Sort and index candidates
    samtools sort -@ $THREADS -o circular_candidates.sorted.bam circular_candidates.bam
    samtools index circular_candidates.sorted.bam

    # Run Circle-Map Realign
    echo "  Running Circle-Map Realign..."
    Circle-Map Realign -t $THREADS \
        -i circular_candidates.sorted.bam \
        -qbam is_aligned.qname.bam \
        -sbam is_aligned.bam \
        -fasta is_reference.fna \
        -o is_circles.bed 2>realign.log

    # Check Circle-Map output
    if [[ ! -s is_circles.bed ]]; then
        echo "  No circles detected by Circle-Map"
        echo "COMPLETED:0:NO_CIRCLES" > status.txt
        echo -e "is_id\tis_length\tcircle_start\tcircle_end\tcircle_length\tlength_ratio\tsplit_reads\tdiscordant_reads\tscore\tvalidation" > circularized_is.tsv
        rm -f *.bam *.bam.bai
        exit 0
    fi

    CIRCLE_COUNT=$(wc -l < is_circles.bed)
    echo "  Raw circles detected: $CIRCLE_COUNT"

    # ============================================
    # Step 4: Validate Results
    # ============================================
    echo ""
    echo "Step 4: Validating circles against IS lengths..."

    python "$SCRIPTS_DIR/validate_is_circles.py" \
        --circles is_circles.bed \
        --reference is_reference.fna \
        --output circularized_is.tsv \
        --tolerance 0.1 \
        --min-reads 2

    # ============================================
    # Summary
    # ============================================
    echo ""
    echo "=============================================="

    if [[ -s circularized_is.tsv ]]; then
        CONFIRMED=$(grep -c "CONFIRMED_CIRCULAR_IS" circularized_is.tsv || echo "0")
        LIKELY=$(grep -c "LIKELY_CIRCULAR_IS" circularized_is.tsv || echo "0")
        TOTAL=$((CONFIRMED + LIKELY))

        echo "RESULTS:"
        echo "  Confirmed circular IS: $CONFIRMED"
        echo "  Likely circular IS: $LIKELY"
        echo "COMPLETED:$TOTAL:CONFIRMED=$CONFIRMED,LIKELY=$LIKELY" > status.txt
    else
        echo "No circularized IS elements detected"
        echo "COMPLETED:0" > status.txt
    fi

    # ============================================
    # Cleanup intermediate files
    # ============================================
    echo ""
    echo "Cleaning up intermediate files..."
    rm -f *.bam *.bam.bai
    rm -f is_reference.fna.*  # BWA index files

    # Keep:
    # - is_reference.fna (for verification)
    # - is_circles.bed (raw Circle-Map output)
    # - circularized_is.tsv (validated results)
    # - status.txt
    # - *.log files

    echo "=============================================="
    echo "Complete: $SAMPLE"
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
    echo ""
    echo "Options:"
    echo "  --dry-run       Show what would be submitted, don't submit"
    echo "  --max-jobs N    Max concurrent SLURM array tasks (default: 50)"
    echo "  --force         Re-process samples even if circularized_is.tsv exists"
    echo "  --batch-size N  Max array size per job (default: 1000)"
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

    # Must have IS FASTA (non-empty, no NO_MGE marker)
    IS_FASTA="${dir}03.results/genome/04.makefasta.genome.all_seqs.fna"
    NO_MGE="${dir}03.results/genome/NO_MGE_DETECTED.txt"

    if [[ -f "$NO_MGE" ]]; then
        continue
    fi
    if [[ ! -f "$IS_FASTA" ]] || [[ ! -s "$IS_FASTA" ]]; then
        continue
    fi
    HAS_IS=$((HAS_IS + 1))

    # Must have FASTQ files
    if [[ ! -f "${dir}00.reads/${SAMPLE}_R1.fastq.gz" ]] || \
       [[ ! -f "${dir}00.reads/${SAMPLE}_R2.fastq.gz" ]]; then
        continue
    fi
    HAS_FASTQ=$((HAS_FASTQ + 1))

    # Check if already completed
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

# --- Calculate number of batches needed ---
NUM_BATCHES=$(( (NEED_PROCESSING + BATCH_SIZE - 1) / BATCH_SIZE ))

# --- Split queue file into batch files ---
echo "Creating ${NUM_BATCHES} batch queue file(s)..."
for ((batch=1; batch<=NUM_BATCHES; batch++)); do
    start=$(( (batch - 1) * BATCH_SIZE + 1 ))
    end=$(( batch * BATCH_SIZE ))
    [[ $end -gt $NEED_PROCESSING ]] && end=$NEED_PROCESSING

    batch_file="${SCRIPTS_DIR}/.is_circle_queue_batch${batch}.txt"
    sed -n "${start},${end}p" "$QUEUE_FILE" > "$batch_file"

    batch_size=$(wc -l < "$batch_file")
    echo "  Batch ${batch}: ${batch_size} samples (lines ${start}-${end})"
done
echo ""

if [[ "$DRY_RUN" == "true" ]]; then
    echo "[DRY RUN] Would submit ${NUM_BATCHES} batch(es):"
    for ((batch=1; batch<=NUM_BATCHES; batch++)); do
        batch_file="${SCRIPTS_DIR}/.is_circle_queue_batch${batch}.txt"
        batch_size=$(wc -l < "$batch_file")
        echo "  Batch ${batch}: sbatch --array=1-${batch_size}%${MAX_JOBS} $0 ${BASE_DIR} ${batch}"
        echo "            (${batch_size} tasks)"
    done
    echo ""
    echo "Full sample list saved to: ${QUEUE_FILE}"
    exit 0
fi

# --- Submit SLURM array job(s) in batches ---
SCRIPT_PATH="$(realpath "$0")"
SUBMITTED_JOBS=()

echo "Submitting ${NUM_BATCHES} batch(es)..."
echo ""

for ((batch=1; batch<=NUM_BATCHES; batch++)); do
    batch_file="${SCRIPTS_DIR}/.is_circle_queue_batch${batch}.txt"
    batch_size=$(wc -l < "$batch_file")

    JOB_ID=$(sbatch \
        --array="1-${batch_size}%${MAX_JOBS}" \
        --output="${LOG_DIR}/is_circle_%A_%a.out" \
        --error="${LOG_DIR}/is_circle_%A_%a.err" \
        --parsable \
        "$SCRIPT_PATH" "$BASE_DIR" "$batch")

    SUBMITTED_JOBS+=("$JOB_ID")

    echo "Batch ${batch}/${NUM_BATCHES}:"
    echo "  Job ID: ${JOB_ID}"
    echo "  Array: 1-${batch_size} (${batch_size} tasks, max ${MAX_JOBS} concurrent)"
    echo "  Queue: .is_circle_queue_batch${batch}.txt"
    echo ""

    # Small delay to avoid overwhelming scheduler
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
echo "  sacct -j ${SUBMITTED_JOBS[0]} --format=JobID,State,Elapsed"
echo ""
echo "Check progress:"
echo "  ls -1 ${SAMPLES_DIR}/SAM*/is_circle_analysis/circularized_is.tsv 2>/dev/null | wc -l"
echo "============================================"
