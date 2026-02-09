#!/bin/bash -l
#SBATCH --job-name=prepare_sample
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32      # <--- 修改这里：填满节点 (lr6 主要是 32 核)
#SBATCH --mem=0                 # <--- 修改这里：使用全部内存 (避免 32G 不够用)

#===============================================================================
# prepare_all.sh (Fixed batched version)
#
# Auto-discovers sample directories and submits SLURM array jobs.
# Automatically splits large arrays into batches to respect cluster limits.
#===============================================================================

set -euo pipefail

# ───────────────────────────────────────────────────────────────────────────
# MODE DETECTION: launcher vs. SLURM worker
# ───────────────────────────────────────────────────────────────────────────

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    # =====================================================================
    # WORKER MODE: process a single sample
    # =====================================================================
    
    BASE_DIR="$1"
    BATCH_NUM="${2:-1}"  # Batch number passed as argument
    SAMPLE_LIST="${BASE_DIR}/scripts/.prepare_queue_batch${BATCH_NUM}.txt"
    CONDA_ENV="mgefinder"
    THREADS="${SLURM_CPUS_PER_TASK:-8}"
    
    SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
    SAMPLE_DIR="${BASE_DIR}/${SAMPLE}"
    
    echo "============================================"
    echo "Prepare Sample for MGEfinder"
    echo "============================================"
    echo "Date: $(date)"
    echo "Batch: ${BATCH_NUM}, Job ID: ${SLURM_JOB_ID}, Task: ${SLURM_ARRAY_TASK_ID}"
    echo "Sample: ${SAMPLE}"
    echo "Directory: ${SAMPLE_DIR}"
    echo "Threads: ${THREADS}"
    echo "============================================"
    
    if [[ ! -d "$SAMPLE_DIR" ]]; then
        echo "ERROR: Sample directory does not exist: $SAMPLE_DIR"
        exit 1
    fi
    
    # --- Define paths ---
    READS_DIR="${SAMPLE_DIR}/00.reads"
    REF_DIR="${SAMPLE_DIR}/00.reference"
    GENOME_DIR="${SAMPLE_DIR}/00.genome"
    BAM_DIR="${SAMPLE_DIR}/00.bam"
    ASSEMBLY_DIR="${SAMPLE_DIR}/00.assembly"
    
    R1_FILE="${READS_DIR}/${SAMPLE}_R1.fastq.gz"
    R2_FILE="${READS_DIR}/${SAMPLE}_R2.fastq.gz"
    REF_FILE="${REF_DIR}/genome.fna"
    REF_GZ="${REF_DIR}/genome.fna.gz"
    
    # --- Decompress reference if needed ---
    if [[ ! -f "$REF_FILE" ]] && [[ -f "$REF_GZ" ]]; then
        echo "Decompressing reference genome..."
        gunzip -k "$REF_GZ"
    fi
    
    # --- Validate inputs ---
    for f in "$R1_FILE" "$R2_FILE" "$REF_FILE"; do
        if [[ ! -s "$f" ]]; then
            echo "ERROR: Missing or empty: $f"
            exit 1
        fi
    done
    
    echo "Input files:"
    echo "  R1: $(ls -lh "$R1_FILE" | awk '{print $5}')"
    echo "  R2: $(ls -lh "$R2_FILE" | awk '{print $5}')"
    echo "  Reference: $(ls -lh "$REF_FILE" | awk '{print $5}')"
    echo "============================================"
    
    # --- Activate conda ---
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV"
    
    echo "Conda: $CONDA_ENV"
    echo "  bwa: $(bwa 2>&1 | grep Version || echo 'installed')"
    echo "  samtools: $(samtools --version | head -1)"
    echo "  skesa: $(skesa --version 2>&1 || echo 'unknown')"
    echo "============================================"
    
    GENOME_REF="${GENOME_DIR}/genome.fna"
    BAM_FILE="${BAM_DIR}/${SAMPLE}.genome.bam"
    ASSEMBLY_FILE="${ASSEMBLY_DIR}/${SAMPLE}.fna"
    
    # --- Step 1: 00.genome ---
    echo "[Step 1/5] Setting up 00.genome..."
    mkdir -p "$GENOME_DIR"
    if [[ ! -f "$GENOME_REF" ]]; then
        cp "$REF_FILE" "$GENOME_REF"
        echo "  Copied reference"
    else
        echo "  Already exists, skipping"
    fi
    
    # --- Step 2: BWA index ---
    echo "[Step 2/5] Building BWA index..."
    if [[ ! -f "${GENOME_REF}.bwt" ]]; then
        bwa index "$GENOME_REF" 2>&1
        echo "  Index created"
    else
        echo "  Already exists, skipping"
    fi
    
    # --- Step 3: SKESA assembly ---
    echo "[Step 3/5] Assembling with SKESA..."
    mkdir -p "$ASSEMBLY_DIR"
    if [[ -s "$ASSEMBLY_FILE" ]]; then
        echo "  Already exists ($(grep -c '^>' "$ASSEMBLY_FILE") contigs), skipping"
    else
        MEMORY_GB="${SKESA_MEMORY_GB:-30}"
        skesa \
            --reads "$R1_FILE","$R2_FILE" \
            --cores "$THREADS" \
            --memory "$MEMORY_GB" \
            --contigs_out "$ASSEMBLY_FILE"
        
        if [[ $? -ne 0 ]] || [[ ! -s "$ASSEMBLY_FILE" ]]; then
            echo "ERROR: SKESA assembly failed"
            exit 1
        fi
        
        N_CONTIGS=$(grep -c '^>' "$ASSEMBLY_FILE")
        TOTAL_BP=$(grep -v '^>' "$ASSEMBLY_FILE" | tr -d '\n' | wc -c)
        echo "  Done: ${N_CONTIGS} contigs, ${TOTAL_BP} bp"
    fi
    
    # --- Step 4: BWA MEM alignment ---
    echo "[Step 4/5] Aligning with BWA MEM..."
    mkdir -p "$BAM_DIR"
    if [[ -f "$BAM_FILE" ]] && [[ -f "${BAM_FILE}.bai" ]]; then
        echo "  Already exists, skipping"
    else
        bwa mem \
            -t "$THREADS" \
            "$GENOME_REF" \
            "$R1_FILE" \
            "$R2_FILE" \
            2> "${BAM_DIR}/${SAMPLE}.bwa.log" \
        | samtools sort \
            -@ "$THREADS" \
            -o "$BAM_FILE" \
            -

        if [[ $? -ne 0 ]]; then
            echo "ERROR: Alignment failed"
            cat "${BAM_DIR}/${SAMPLE}.bwa.log"
            exit 1
        fi

        echo "  Done"
        cat "${BAM_DIR}/${SAMPLE}.bwa.log"
    fi
    
    # --- Step 5: Index BAM ---
    echo "[Step 5/5] Indexing BAM..."
    if [[ ! -f "${BAM_FILE}.bai" ]]; then
        samtools index -@ "$THREADS" "$BAM_FILE"
        echo "  Done"
    else
        echo "  Already exists, skipping"
    fi
    
    # --- Verify ---
    echo "============================================"
    echo "Verification"
    echo "============================================"
    ALL_OK=true
    for f in "$GENOME_REF" "${GENOME_REF}.bwt" "$ASSEMBLY_FILE" "$BAM_FILE" "${BAM_FILE}.bai"; do
        if [[ -s "$f" ]]; then
            echo "  OK $(basename "$f") ($(ls -lh "$f" | awk '{print $5}'))"
        else
            echo "  FAIL $(basename "$f")"
            ALL_OK=false
        fi
    done
    echo ""
    
    if $ALL_OK; then
        echo "READY: ${SAMPLE}"
    else
        echo "FAILED: ${SAMPLE}"
        exit 1
    fi
    
    echo "Completed: $(date)"
    echo "============================================"
    exit 0
fi

# =========================================================================
# LAUNCHER MODE: discover samples and submit SLURM array job(s)
# =========================================================================

# --- Parse arguments ---
BASE_DIR=""
DRY_RUN=false
MAX_JOBS=20
FORCE=false
BATCH_SIZE=1000  # Max array size per job

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
    echo "Usage: bash $0 /path/to/mgefinder [--dry-run] [--max-jobs N] [--force] [--batch-size N]"
    echo ""
    echo "Options:"
    echo "  --dry-run       Show what would be submitted, don't submit"
    echo "  --max-jobs N    Max concurrent SLURM array tasks (default: 20)"
    echo "  --force         Re-process samples even if outputs already exist"
    echo "  --batch-size N  Max array size per job (default: 1000)"
    exit 1
fi

BASE_DIR="$(realpath "$BASE_DIR")"
QUEUE_FILE="${BASE_DIR}/scripts/.prepare_queue.txt"
LOG_DIR="${BASE_DIR}/scripts/logs"
mkdir -p "$LOG_DIR"

echo "============================================"
echo "MGEfinder Sample Discovery"
echo "============================================"
echo "Base directory: ${BASE_DIR}"
echo "Scanning for samples..."
echo ""

# --- Discover samples ---
TMPDIR_SCAN="${BASE_DIR}/scripts/.scan_tmp"
mkdir -p "$TMPDIR_SCAN"

echo "  Finding samples with reads..."
find "$BASE_DIR" -maxdepth 3 -name "*_R1.fastq.gz" -path "*/SAM*/00.reads/*" -size +0c \
    | sed 's|.*/\(SAM[^/]*\)/00.reads/.*|\1|' \
    | sort -u > "${TMPDIR_SCAN}/has_r1.txt"

TOTAL=$(wc -l < "${TMPDIR_SCAN}/has_r1.txt")
echo "  Found ${TOTAL} samples with R1 reads"

echo "  Checking R2 and reference files..."
> "$QUEUE_FILE"
MISSING_INPUT=0
READY=0
NEED_PROCESSING=0

while IFS= read -r SAMPLE; do
    dir="${BASE_DIR}/${SAMPLE}/"
    
    # Quick checks for R2 and reference
    if [[ ! -s "${dir}00.reads/${SAMPLE}_R2.fastq.gz" ]]; then
        MISSING_INPUT=$((MISSING_INPUT + 1))
        continue
    fi
    
    if [[ ! -s "${dir}00.reference/genome.fna" ]] && [[ ! -s "${dir}00.reference/genome.fna.gz" ]]; then
        MISSING_INPUT=$((MISSING_INPUT + 1))
        continue
    fi
    
    # Check if already fully prepared
    if [[ "$FORCE" == "false" ]]; then
        if [[ -s "${dir}00.bam/${SAMPLE}.genome.bam" ]] && \
           [[ -s "${dir}00.bam/${SAMPLE}.genome.bam.bai" ]] && \
           [[ -s "${dir}00.assembly/${SAMPLE}.fna" ]] && \
           [[ -s "${dir}00.genome/genome.fna" ]]; then
            READY=$((READY + 1))
            continue
        fi
    fi
    
    echo "$SAMPLE" >> "$QUEUE_FILE"
    NEED_PROCESSING=$((NEED_PROCESSING + 1))
done < "${TMPDIR_SCAN}/has_r1.txt"

rm -rf "$TMPDIR_SCAN"

echo "Results:"
echo "  Total sample dirs:    ${TOTAL}"
echo "  Missing inputs:       ${MISSING_INPUT}"
echo "  Already prepared:     ${READY}"
echo "  Need processing:      ${NEED_PROCESSING}"
echo "============================================"

if [[ "$NEED_PROCESSING" -eq 0 ]]; then
    echo "Nothing to do. All samples with valid inputs are already prepared."
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
echo "Creating ${NUM_BATCHES} batch queue files..."
for ((batch=1; batch<=NUM_BATCHES; batch++)); do
    start=$(( (batch - 1) * BATCH_SIZE + 1 ))
    end=$(( batch * BATCH_SIZE ))
    [[ $end -gt $NEED_PROCESSING ]] && end=$NEED_PROCESSING
    
    batch_file="${BASE_DIR}/scripts/.prepare_queue_batch${batch}.txt"
    sed -n "${start},${end}p" "$QUEUE_FILE" > "$batch_file"
    
    batch_size=$(wc -l < "$batch_file")
    echo "  Batch ${batch}: ${batch_size} samples (lines ${start}-${end})"
done
echo ""

if [[ "$DRY_RUN" == "true" ]]; then
    echo "[DRY RUN] Would submit ${NUM_BATCHES} batch(es):"
    for ((batch=1; batch<=NUM_BATCHES; batch++)); do
        batch_file="${BASE_DIR}/scripts/.prepare_queue_batch${batch}.txt"
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
    batch_file="${BASE_DIR}/scripts/.prepare_queue_batch${batch}.txt"
    batch_size=$(wc -l < "$batch_file")
    
    JOB_ID=$(sbatch \
        --array="1-${batch_size}%${MAX_JOBS}" \
        --parsable \
        "$SCRIPT_PATH" "$BASE_DIR" "$batch")
    
    SUBMITTED_JOBS+=("$JOB_ID")
    
    echo "Batch ${batch}/${NUM_BATCHES}:"
    echo "  Job ID: ${JOB_ID}"
    echo "  Array: 1-${batch_size} (${batch_size} tasks, max ${MAX_JOBS} concurrent)"
    echo "  Queue: .prepare_queue_batch${batch}.txt"
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
echo "  squeue -u $USER                           # all your jobs"
echo "  squeue -j ${SUBMITTED_JOBS[0]}            # first batch"
echo "  tail -f ${LOG_DIR}/prepare_sample_${SUBMITTED_JOBS[0]}_1.out  # first task log"
echo "  sacct -j ${SUBMITTED_JOBS[0]} --format=JobID,State,Elapsed    # first batch summary"
echo ""
echo "Check progress:"
echo "  ls -1 ${BASE_DIR}/*/00.bam/*.bam 2>/dev/null | wc -l  # completed BAM files"
echo "============================================"
