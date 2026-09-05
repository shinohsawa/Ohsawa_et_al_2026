#!/bin/bash
# ==============================================================================
# SRA Data Download and Formatting Pipeline (with Prefetch)
# Project: Ohsawa et al. (2026)
# ==============================================================================

set -e

# ====================================================================
# 0. Setup & Paths
# ====================================================================
# Set paths relative to the project root (assuming execution from the script directory)
METADATA="../metadata.csv"
OUTDIR="../rawdata"
TMPDIR="../rawdata/tmp_download" # Temporary directory for raw downloads

# Create necessary directories
mkdir -p "$OUTDIR"
mkdir -p "$TMPDIR"

echo "Starting robust download pipeline..."

# ====================================================================
# Pipeline Execution Loop
# ====================================================================
# Filter "RNA-Seq" entries from the metadata and process them in a loop
grep "RNA-Seq" "$METADATA" | while IFS=, read -r Run SampleName LibraryName AssayType
do
    # Skip if files already exist (allows resuming interrupted downloads)
    if [ -f "$OUTDIR/${LibraryName}_R1.fq.gz" ] && [ -f "$OUTDIR/${LibraryName}_R2.fq.gz" ]; then
        echo "Skipping $LibraryName (Already downloaded)"
        continue
    fi

    echo "=========================================================="
    echo " Processing: $LibraryName (SRA: $Run)"
    echo "=========================================================="

    # 1. Pre-fetch the SRA file for stability (Downloads .sra to a local cache/folder)
    echo " [1/3] Running prefetch to secure data locally..."
    prefetch "$Run" --max-size 50G

    # 2. Extract paired-end FASTQ files from the prefetched data
    echo " [2/3] Extracting FASTQ files with fasterq-dump..."
    fasterq-dump --split-files --threads 8 -O "$TMPDIR" "$Run"

    # Locate the extracted files (e.g., SRR123456_1.fastq)
    R1_FILE=$(ls "$TMPDIR"/*_1.fastq 2>/dev/null | head -n 1)
    R2_FILE=$(ls "$TMPDIR"/*_2.fastq 2>/dev/null | head -n 1)

    # 3. Compress (gzip) and save to the output directory using the LibraryName format
    echo " [3/3] Compressing and renaming..."
    if [ -n "$R1_FILE" ]; then
        gzip -c "$R1_FILE" > "$OUTDIR/${LibraryName}_R1.fq.gz"
    fi

    if [ -n "$R2_FILE" ]; then
        gzip -c "$R2_FILE" > "$OUTDIR/${LibraryName}_R2.fq.gz"
    fi

    # 4. Clean up temporary FASTQ files and the prefetch cache directory for this sample
    rm -f "$TMPDIR"/*.fastq
    rm -rf "$Run" 2>/dev/null || true

    echo " -> Successfully saved as ${LibraryName}_R1/R2.fq.gz"
done

# Remove the temporary directory
rmdir "$TMPDIR" 2>/dev/null || true

echo "=========================================================="
echo " All RNA-Seq datasets have been downloaded and formatted! "
echo "=========================================================="