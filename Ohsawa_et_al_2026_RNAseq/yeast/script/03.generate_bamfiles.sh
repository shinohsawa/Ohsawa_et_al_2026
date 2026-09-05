#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Directory and File Settings
REF_GENOME="../reference/star_index"
REF_GTF="../reference/cerevisiae_albicans_NLS_mcherry.gtf"
FASTP_DIR="../fastp"
BAM_DIR="../bamfiles"
THREADS=12

# Create output directory
mkdir -p "$BAM_DIR"

echo "=========================================================="
echo " Starting STAR Alignment for all samples "
echo "=========================================================="

# Loop through all cleaned Read 1 files automatically
for R1_FILE in "$FASTP_DIR"/*_R1.clean.fastq.gz; do
    
    # Skip if no files found
    [ -e "$R1_FILE" ] || continue

    # Extract sample name
    BASENAME=$(basename "$R1_FILE")
    SAMPLE="${BASENAME%_R1.clean.fastq.gz}"
    R2_FILE="${FASTP_DIR}/${SAMPLE}_R2.clean.fastq.gz"

    echo "----------------------------------------------------------"
    echo "Processing sample: $SAMPLE"

    if [[ ! -f "$R2_FILE" ]]; then
        echo "Error: Cleaned Read 2 file not found for $SAMPLE!"
        continue
    fi

    # Define output prefix for STAR
    STAR_PREFIX="${BAM_DIR}/${SAMPLE}_"
    SORTED_BAM="${STAR_PREFIX}Aligned.sortedByCoord.out.bam"

    # Skip if BAM already exists (useful for resuming interrupted jobs)
    if [ -f "$SORTED_BAM" ]; then
        echo "  Already mapped $SAMPLE, skipping..."
        continue
    fi

    # CRITICAL FIX: Use the system's local /tmp directory which supports FIFO pipes
    STAR_TMP="/tmp/star_${SAMPLE}"
    rm -rf "$STAR_TMP"

    # Run STAR Mapping (Default multimap settings: max 10)
    echo "  Running STAR..."
    STAR --runThreadN "$THREADS" \
         --twopassMode Basic \
         --genomeDir "$REF_GENOME" \
         --sjdbGTFfile "$REF_GTF" \
         --outTmpDir "$STAR_TMP" \
         --readFilesIn "$R1_FILE" "$R2_FILE" \
         --readFilesCommand gunzip -c \
         --sjdbOverhang 149 \
         --outFilterMismatchNmax 3 \
         --chimSegmentMin 2 \
         --alignIntronMax 299999 \
         --outSAMstrandField intronMotif \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix "$STAR_PREFIX"

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

for BAM_FILE in "$BAM_DIR"/*.bam; do
    [ -e "$BAM_FILE" ] || continue
    
    echo "Analyzing $BAM_FILE ..."
    echo "========================================" >> "$STATS_FILE"
    echo "Sample: $(basename "$BAM_FILE")" >> "$STATS_FILE"
    
    bamtools stats -in "$BAM_FILE" >> "$STATS_FILE"
    
    echo "--- Mapping reads per chromosome ---" >> "$STATS_FILE"
    samtools idxstats "$BAM_FILE" >> "$STATS_FILE"
    echo -e "\n" >> "$STATS_FILE"
done

echo "=========================================================="
echo " All processing complete! Stats saved to $STATS_FILE "
echo "=========================================================="