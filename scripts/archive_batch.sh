#!/bin/bash
#===============================================================================
# archive_batch.sh
#
# Archive a batch directory into a tar.gz, excluding large data files
# (FASTQ reads, BAM alignments, BAM indexes, BWA indexes).
#
# Usage:
#   bash archive_batch.sh /path/to/mgefinder_batchN [/path/to/output.tar.gz]
#
# Output defaults to: {BASE_DIR}/../mgefinder_batchN_results.tar.gz
#===============================================================================

set -euo pipefail

BASE_DIR="${1:-}"

if [[ -z "$BASE_DIR" ]]; then
    echo "Usage: bash $0 /path/to/mgefinder_batchN [/path/to/output.tar.gz]"
    echo ""
    echo "Archives the batch directory excluding FASTQ, BAM, and BWA index files."
    echo ""
    echo "Excluded file types:"
    echo "  *.fastq.gz    (raw reads, ~1-3 GB per sample)"
    echo "  *.bam         (alignments)"
    echo "  *.bam.bai     (BAM indexes)"
    echo "  *.bwt *.pac *.ann *.amb *.sa  (BWA indexes, regenerable)"
    echo ""
    echo "Example:"
    echo "  bash $0 /global/scratch/users/kh36969/mgefinder_batch1"
    echo "  bash $0 /global/scratch/users/kh36969/mgefinder_batch2 /path/to/batch2_results.tar.gz"
    exit 1
fi

BASE_DIR="$(realpath "$BASE_DIR")"
BATCH_NAME="$(basename "$BASE_DIR")"
OUTPUT="${2:-$(dirname "$BASE_DIR")/${BATCH_NAME}_results.tar.gz}"

# Verify input exists
if [[ ! -d "$BASE_DIR" ]]; then
    echo "ERROR: Directory not found: $BASE_DIR"
    exit 1
fi

# Check output location is writable
OUTPUT_DIR="$(dirname "$OUTPUT")"
if [[ ! -d "$OUTPUT_DIR" ]]; then
    echo "ERROR: Output directory does not exist: $OUTPUT_DIR"
    exit 1
fi

echo "============================================"
echo "Archive Batch: ${BATCH_NAME}"
echo "============================================"
echo "Source:  ${BASE_DIR}"
echo "Output:  ${OUTPUT}"
echo ""

# Show what we're excluding
echo "Scanning directory..."
TOTAL_SIZE=$(du -sh "$BASE_DIR" 2>/dev/null | cut -f1)
echo "  Total size (with data files): ${TOTAL_SIZE}"

# Count excluded files
N_FASTQ=$(find "$BASE_DIR" -name "*.fastq.gz" 2>/dev/null | wc -l)
N_BAM=$(find "$BASE_DIR" -name "*.bam" -o -name "*.bam.bai" 2>/dev/null | wc -l)
N_BWA=$(find "$BASE_DIR" -name "*.bwt" -o -name "*.pac" -o -name "*.ann" -o -name "*.amb" -o -name "*.sa" 2>/dev/null | wc -l)

echo "  FASTQ files to skip: ${N_FASTQ}"
echo "  BAM/BAI files to skip: ${N_BAM}"
echo "  BWA index files to skip: ${N_BWA}"
echo ""

# Estimate archive size (everything except excluded files)
echo "Estimating archive size (this may take a moment)..."
ARCHIVE_SIZE=$(du -sh --exclude="*.fastq.gz" --exclude="*.bam" --exclude="*.bam.bai" \
    --exclude="*.bwt" --exclude="*.pac" --exclude="*.ann" --exclude="*.amb" --exclude="*.sa" \
    "$BASE_DIR" 2>/dev/null | cut -f1)
echo "  Estimated archive content: ${ARCHIVE_SIZE}"
echo ""

echo "Creating archive..."
echo "  (This may take a while for large batches)"
echo ""

NCPU=$(nproc 2>/dev/null || echo 1)
if command -v pigz &>/dev/null && [[ "$NCPU" -gt 1 ]]; then
    echo "  Using pigz with ${NCPU} CPUs for parallel compression"
    tar cf - \
        -C "$(dirname "$BASE_DIR")" \
        --exclude="*.fastq.gz" \
        --exclude="*.bam" \
        --exclude="*.bam.bai" \
        --exclude="*.bwt" \
        --exclude="*.pac" \
        --exclude="*.ann" \
        --exclude="*.amb" \
        --exclude="*.sa" \
        "$BATCH_NAME" | pigz -p "$NCPU" > "$OUTPUT"
else
    echo "  Using single-threaded gzip"
    tar czf "$OUTPUT" \
        -C "$(dirname "$BASE_DIR")" \
        --exclude="*.fastq.gz" \
        --exclude="*.bam" \
        --exclude="*.bam.bai" \
        --exclude="*.bwt" \
        --exclude="*.pac" \
        --exclude="*.ann" \
        --exclude="*.amb" \
        --exclude="*.sa" \
        "$BATCH_NAME"
fi

ARCHIVE_SIZE_ACTUAL=$(du -sh "$OUTPUT" | cut -f1)

echo "============================================"
echo "Done!"
echo "============================================"
echo "Archive: ${OUTPUT}"
echo "Size:    ${ARCHIVE_SIZE_ACTUAL}"
echo ""
echo "To extract:"
echo "  tar xzf ${OUTPUT} -C /path/to/destination"
echo ""
echo "To list contents:"
echo "  tar tzf ${OUTPUT} | head -50"
echo "============================================"
