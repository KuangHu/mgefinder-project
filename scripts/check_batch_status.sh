#!/bin/bash
# =============================================================================
# check_batch_status.sh - Check MGEfinder batch download/processing status
#
# Usage:
#   bash check_batch_status.sh /path/to/batch/folder
#   bash check_batch_status.sh /path/to/batch/folder --quick    # Skip slow counts
#   bash check_batch_status.sh /path/to/batch/folder --verbose  # Show more details
#
# Example:
#   bash check_batch_status.sh /global/scratch/users/kh36969/mgefinder_batch2
# =============================================================================

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Parse arguments
BATCH_DIR="${1:-.}"
QUICK_MODE=false
VERBOSE=false

for arg in "$@"; do
    case $arg in
        --quick) QUICK_MODE=true ;;
        --verbose) VERBOSE=true ;;
    esac
done

# Validate directory
if [[ ! -d "$BATCH_DIR" ]]; then
    echo -e "${RED}Error: Directory not found: $BATCH_DIR${NC}"
    exit 1
fi

BATCH_DIR=$(cd "$BATCH_DIR" && pwd)
SAMPLES_DIR="$BATCH_DIR/samples"
JOBS_DIR="$BATCH_DIR/jobs"

echo "============================================================"
echo -e "${BLUE}MGEfinder Batch Status Report${NC}"
echo "============================================================"
echo "Batch directory: $BATCH_DIR"
echo "Report time: $(date)"
echo ""

# -----------------------------------------------------------------------------
# 1. Sample Folder Counts
# -----------------------------------------------------------------------------
echo -e "${YELLOW}=== Sample Downloads ===${NC}"

if [[ -d "$SAMPLES_DIR" ]]; then
    # Count sample folders (SAM* pattern)
    TOTAL_FOLDERS=$(ls -d "$SAMPLES_DIR"/SAM* 2>/dev/null | wc -l)
    echo "Sample folders: $TOTAL_FOLDERS"

    # Count files using find (faster than loops)
    echo ""
    echo "File counts:"
    REF_COUNT=$(find "$SAMPLES_DIR" -path '*/00.reference/genome.fna' -type f 2>/dev/null | wc -l)
    ASM_COUNT=$(find "$SAMPLES_DIR" -path '*/00.assembly/*.fna' -type f 2>/dev/null | wc -l)
    R1_COUNT=$(find "$SAMPLES_DIR" -name '*_R1.fastq.gz' -path '*/00.reads/*' -type f 2>/dev/null | wc -l)
    R2_COUNT=$(find "$SAMPLES_DIR" -name '*_R2.fastq.gz' -path '*/00.reads/*' -type f 2>/dev/null | wc -l)

    printf "  %-20s %6d\n" "References:" "$REF_COUNT"
    printf "  %-20s %6d\n" "Assemblies:" "$ASM_COUNT"
    printf "  %-20s %6d\n" "R1 reads:" "$R1_COUNT"
    printf "  %-20s %6d\n" "R2 reads:" "$R2_COUNT"

    # Count fully complete samples (unless quick mode)
    if [[ "$QUICK_MODE" == false ]]; then
        echo ""
        echo -n "Counting fully complete samples (ref+asm+R1+R2)... "
        COMPLETE=0
        for d in "$SAMPLES_DIR"/SAM*/; do
            [[ ! -d "$d" ]] && continue
            s=$(basename "$d")
            if [[ -f "$d/00.reference/genome.fna" ]] && \
               [[ -f "$d/00.assembly/${s}.fna" ]] && \
               [[ -f "$d/00.reads/${s}_R1.fastq.gz" ]] && \
               [[ -f "$d/00.reads/${s}_R2.fastq.gz" ]]; then
                ((COMPLETE++))
            fi
        done
        echo -e "${GREEN}$COMPLETE${NC}"

        if [[ $TOTAL_FOLDERS -gt 0 ]]; then
            PCT=$(echo "scale=1; $COMPLETE * 100 / $TOTAL_FOLDERS" | bc)
            echo "  Completion rate: ${PCT}%"
        fi
    else
        echo ""
        echo "(Use without --quick to count fully complete samples)"
    fi

    # Count BAM files (prepare step output)
    BAM_COUNT=$(find "$SAMPLES_DIR" -name '*.bam' -path '*/00.bam/*' -type f 2>/dev/null | wc -l)
    if [[ $BAM_COUNT -gt 0 ]]; then
        echo ""
        echo -e "${YELLOW}=== Prepare Step ===${NC}"
        printf "  %-20s %6d\n" "BAM files:" "$BAM_COUNT"
    fi

    # Count results (mgefinder step output)
    RESULTS_COUNT=$(find "$SAMPLES_DIR" -name '01.clusterseq.genome.tsv' -path '*/03.results/*' -type f 2>/dev/null | wc -l)
    NO_MGE_COUNT=$(find "$SAMPLES_DIR" -name 'NO_MGE_DETECTED.txt' -path '*/03.results/*' -type f 2>/dev/null | wc -l)
    if [[ $RESULTS_COUNT -gt 0 ]] || [[ $NO_MGE_COUNT -gt 0 ]]; then
        echo ""
        echo -e "${YELLOW}=== MGEfinder Results ===${NC}"
        printf "  %-20s %6d\n" "With MGEs found:" "$RESULTS_COUNT"
        printf "  %-20s %6d\n" "No MGEs detected:" "$NO_MGE_COUNT"
        printf "  %-20s %6d\n" "Total processed:" "$((RESULTS_COUNT + NO_MGE_COUNT))"
    fi
else
    echo -e "${RED}Samples directory not found: $SAMPLES_DIR${NC}"
fi

# -----------------------------------------------------------------------------
# 2. Job Status (from SLURM logs)
# -----------------------------------------------------------------------------
echo ""
echo -e "${YELLOW}=== Job Logs Analysis ===${NC}"

if [[ -d "$JOBS_DIR" ]]; then
    # Find all log directories
    LOG_DIRS=$(find "$JOBS_DIR" -type d -name 'logs' 2>/dev/null)

    if [[ -n "$LOG_DIRS" ]]; then
        TOTAL_LOGS=0
        COMPLETE_LOGS=0
        ERROR_LOGS=0

        for log_dir in $LOG_DIRS; do
            if [[ -d "$log_dir" ]]; then
                # Count logs with "Complete:" marker
                complete=$(grep -l "Complete:" "$log_dir"/*.out 2>/dev/null | wc -l)
                errors=$(grep -l "ERROR:" "$log_dir"/*.out 2>/dev/null | wc -l)
                total=$(ls "$log_dir"/*.out 2>/dev/null | wc -l)

                TOTAL_LOGS=$((TOTAL_LOGS + total))
                COMPLETE_LOGS=$((COMPLETE_LOGS + complete))
                ERROR_LOGS=$((ERROR_LOGS + errors))
            fi
        done

        printf "  %-25s %6d\n" "Total job logs:" "$TOTAL_LOGS"
        printf "  %-25s %6d\n" "Completed successfully:" "$COMPLETE_LOGS"
        printf "  %-25s %6d\n" "With errors:" "$ERROR_LOGS"

        # Show error summary if verbose
        if [[ "$VERBOSE" == true ]] && [[ $ERROR_LOGS -gt 0 ]]; then
            echo ""
            echo "  Top error types:"
            for log_dir in $LOG_DIRS; do
                grep -h "ERROR:" "$log_dir"/*.out 2>/dev/null
            done | sort | uniq -c | sort -rn | head -10 | while read count msg; do
                printf "    %5d  %s\n" "$count" "$msg"
            done
        fi
    else
        echo "  No log directories found"
    fi
else
    echo "  Jobs directory not found: $JOBS_DIR"
fi

# -----------------------------------------------------------------------------
# 3. SLURM Job Status (if squeue available)
# -----------------------------------------------------------------------------
echo ""
echo -e "${YELLOW}=== Current SLURM Jobs ===${NC}"

if command -v squeue &> /dev/null; then
    USER=$(whoami)
    RUNNING=$(squeue -u "$USER" --noheader 2>/dev/null | wc -l)

    if [[ $RUNNING -gt 0 ]]; then
        echo "  Running/pending jobs: $RUNNING"
        echo ""
        echo "  Job breakdown:"
        squeue -u "$USER" --format="%.30j" --noheader 2>/dev/null | \
            sed 's/_[0-9]*$//' | sort | uniq -c | sort -rn | head -10 | \
            while read count name; do
                printf "    %5d  %s\n" "$count" "$name"
            done
    else
        echo "  No running/pending jobs"
    fi
else
    echo "  (squeue not available)"
fi

# -----------------------------------------------------------------------------
# 4. Disk Usage
# -----------------------------------------------------------------------------
echo ""
echo -e "${YELLOW}=== Disk Usage ===${NC}"

if [[ -d "$SAMPLES_DIR" ]]; then
    USAGE=$(du -sh "$SAMPLES_DIR" 2>/dev/null | cut -f1)
    echo "  Samples directory: $USAGE"
fi

if [[ -d "$BATCH_DIR/reference_genomes" ]]; then
    REF_USAGE=$(du -sh "$BATCH_DIR/reference_genomes" 2>/dev/null | cut -f1)
    echo "  Reference genomes: $REF_USAGE"
fi

# -----------------------------------------------------------------------------
# 5. Summary
# -----------------------------------------------------------------------------
echo ""
echo "============================================================"
echo -e "${BLUE}Summary${NC}"
echo "============================================================"

if [[ -d "$SAMPLES_DIR" ]] && [[ "$QUICK_MODE" == false ]]; then
    if [[ $COMPLETE -gt 0 ]]; then
        echo -e "Ready for prepare step: ${GREEN}$COMPLETE samples${NC}"
    fi
    if [[ $BAM_COUNT -gt 0 ]]; then
        echo -e "Ready for mgefinder:    ${GREEN}$BAM_COUNT samples${NC}"
    fi
    if [[ $RESULTS_COUNT -gt 0 ]] || [[ $NO_MGE_COUNT -gt 0 ]]; then
        echo -e "Fully processed:        ${GREEN}$((RESULTS_COUNT + NO_MGE_COUNT)) samples${NC}"
    fi
fi

echo ""
echo "Report complete."
