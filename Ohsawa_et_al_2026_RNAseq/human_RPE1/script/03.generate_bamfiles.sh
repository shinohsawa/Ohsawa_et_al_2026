#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Directory and File Settings
REF_GENOME="../reference/STAR_index"

FASTP_DIR="../fastp"
BAM_DIR="../mapped_results"
THREADS=12

# Create output directory
mkdir -p "$BAM_DIR"

echo "=========================================================="
echo " Starting STAR Alignment for Human (TE analysis ready) "
echo "=========================================================="

# 🌟 Loop through _1.clean.fastq.gz files cleaned by Fastp
for R1_FILE in "$FASTP_DIR"/*_1.clean.fastq.gz; do
    
    # Skip if no files found
    [ -e "$R1_FILE" ] || continue

    # Extract sample name
    BASENAME=$(basename "$R1_FILE")
    SAMPLE="${BASENAME%_1.clean.fastq.gz}"
    R2_FILE="${FASTP_DIR}/${SAMPLE}_2.clean.fastq.gz"

    echo "----------------------------------------------------------"
    echo "Processing sample: $SAMPLE"

    if [[ ! -f "$R2_FILE" ]]; then
        echo "Error: Cleaned Read 2 file not found for $SAMPLE!"
        continue
    fi

    # Define output prefix for STAR
    STAR_PREFIX="${BAM_DIR}/${SAMPLE}_"
    SORTED_BAM="${STAR_PREFIX}Aligned.sortedByCoord.out.bam"

    # Skip if BAM already exists
    if [ -f "$SORTED_BAM" ]; then
        echo "  Already mapped $SAMPLE, skipping..."
        continue
    fi

    # Set up temporary directory for stable operation
    STAR_TMP="/tmp/star_${SAMPLE}"
    rm -rf "$STAR_TMP"

    # ==========================================================
    # 🌟 Run STAR mapping (Human & TEtranscripts-specific settings)
    # ==========================================================
    echo "  Running STAR..."
    STAR --runThreadN "$THREADS" \
         --genomeDir "$REF_GENOME" \
         --outTmpDir "$STAR_TMP" \
         --readFilesIn "$R1_FILE" "$R2_FILE" \
         --readFilesCommand gunzip -c \
         --outFilterMultimapNmax 100 \
         --winAnchorMultimapNmax 100 \
         --outFilterMismatchNmax 3 \
         --outSAMstrandField intronMotif \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix "$STAR_PREFIX" \
         --limitBAMsortRAM 10000000000

    # Create index for the sorted BAM file
    echo "  Indexing BAM file..."
    samtools index -@ "$THREADS" "$SORTED_BAM"

    # Clean up temporary STAR files to save space
    rm -rf "$STAR_TMP"
    rm -f "${STAR_PREFIX}Log.progress.out" "${STAR_PREFIX}SJ.out.tab"

    echo "Finished mapping: $SAMPLE"
done

echo "=========================================================="
echo " Alignment completed. Generating Bamtools Stats... "
echo "=========================================================="

STATS_FILE="${BAM_DIR}/bamtools_stats_summary.txt"
> "$STATS_FILE" # Clear the file if it exists

for BAM_FILE in "$BAM_DIR"/*Aligned.sortedByCoord.out.bam; do
    [ -e "$BAM_FILE" ] || continue
    
    echo "Analyzing $BAM_FILE ..."
    echo "========================================" >> "$STATS_FILE"
    echo "Sample: $(basename "$BAM_FILE")" >> "$STATS_FILE"
    
    # Assuming bamtools is installed
    if command -v bamtools &> /dev/null; then
        bamtools stats -in "$BAM_FILE" >> "$STATS_FILE"
    else
        echo "bamtools not found, using samtools flagstat instead" >> "$STATS_FILE"
        samtools flagstat "$BAM_FILE" >> "$STATS_FILE"
    fi
    
    echo "--- Mapping reads per chromosome ---" >> "$STATS_FILE"
    samtools idxstats "$BAM_FILE" >> "$STATS_FILE"
    echo -e "\n" >> "$STATS_FILE"
done

echo "=========================================================="
echo " All processing complete! Stats saved to $STATS_FILE "
echo "=========================================================="