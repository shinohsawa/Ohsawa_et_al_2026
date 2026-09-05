#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# ==========================================
# Directory and File Settings 
# ==========================================
BASE_DIR=".."
CLEAN_DIR="${BASE_DIR}/clean_data"
BAM_DIR="${BASE_DIR}/mapped_results"
INDEX_DIR="${BASE_DIR}/STAR_index"

# 🌟 Flies have a small genome, so we use the server's local /tmp where FIFOs can be created!
TMP_DIR="/tmp/sosawa_star_shi2020"
mkdir -p "${BAM_DIR}" "${TMP_DIR}"

# Reference genome (Drosophila melanogaster BDGP6.32)
REF_FASTA="${BASE_DIR}/reference/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa"
REF_GTF="${BASE_DIR}/reference/Drosophila_melanogaster.BDGP6.32.109.gtf"

# Number of threads to use
THREADS=8

# ==========================================================
# Step 1: Genome Index Generation 
# ==========================================================
if [ ! -f "${INDEX_DIR}/SA" ]; then
    echo "=========================================================="
    echo " 1. Generating STAR genome index for Drosophila..."
    echo "=========================================================="
    mkdir -p "${INDEX_DIR}"
    
    STAR --runThreadN "${THREADS}" \
         --runMode genomeGenerate \
         --genomeDir "${INDEX_DIR}" \
         --genomeFastaFiles "${REF_FASTA}" \
         --sjdbGTFfile "${REF_GTF}" \
         --sjdbOverhang 49 \
         --genomeSAindexNbases 12
         
    echo "Genome index generated successfully!"
else
    echo "=========================================================="
    echo " 1. Genome index already exists. Skipping..."
    echo "=========================================================="
fi

# ==========================================================
# Step 2: STAR Alignment (Single-end & TE-specific settings)
# ==========================================================
echo "=========================================================="
echo " 2. Starting STAR Alignment for Drosophila (TE analysis ready)"
echo "=========================================================="

# Automatically target LC1, LC2, S4, S8 (Single-end) in the clean_data folder
for CLEAN_FILE in "${CLEAN_DIR}"/*.clean.fastq.gz; do
    
    # Skip if no files found
    [ -e "${CLEAN_FILE}" ] || continue

    # Extract sample name (e.g., LC1_1.clean.fastq.gz -> LC1_1)
    BASENAME=$(basename "${CLEAN_FILE}")
    SAMPLE="${BASENAME%.clean.fastq.gz}"

    echo "----------------------------------------------------------"
    echo "Processing sample: ${SAMPLE}"

    # Define output prefix and sorted BAM path
    STAR_PREFIX="${BAM_DIR}/${SAMPLE}_"
    SORTED_BAM="${STAR_PREFIX}Aligned.sortedByCoord.out.bam"

    # Skip if BAM already exists and is not empty
    if [[ -s "${SORTED_BAM}" ]]; then
        echo "  Already mapped ${SAMPLE}, skipping..."
        continue
    fi

    # Set up temporary directory for stable operation (specify a large scratch space!)
    STAR_TMP="${TMP_DIR}/star_${SAMPLE}"
    rm -rf "${STAR_TMP}"

    # ==========================================================
    # 🌟 Run STAR mapping (TEtranscripts-specific parameters!)
    # ==========================================================
    echo "  Running STAR..."
    STAR --runThreadN "${THREADS}" \
         --genomeDir "${INDEX_DIR}" \
         --outTmpDir "${STAR_TMP}" \
         --readFilesIn "${CLEAN_FILE}" \
         --readFilesCommand zcat \
         --outFilterMultimapNmax 100 \
         --winAnchorMultimapNmax 100 \
         --outFilterMismatchNmax 3 \
         --outSAMstrandField intronMotif \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix "${STAR_PREFIX}" \
         --limitBAMsortRAM 10000000000

    # 🌟 Automatically generate an index (.bai) immediately after mapping!
    echo "  Indexing BAM file..."
    samtools index -@ "${THREADS}" "${SORTED_BAM}"

    # Clean up temporary STAR files to save space
    rm -rf "${STAR_TMP}"
    rm -f "${STAR_PREFIX}Log.progress.out" "${STAR_PREFIX}SJ.out.tab"

    echo "Finished mapping: ${SAMPLE}"
done

# ==========================================================
# Step 3: Mapping Summary Stats (Automatic summary aggregation!)
# ==========================================================
echo "=========================================================="
echo " 3. Alignment completed. Generating Mapping Stats Summary..."
echo "=========================================================="

STATS_FILE="${BAM_DIR}/mapping_stats_summary.txt"
> "${STATS_FILE}" # Clear the file if it exists

for BAM_FILE in "${BAM_DIR}"/*Aligned.sortedByCoord.out.bam; do
    [ -e "${BAM_FILE}" ] || continue
    
    echo "Analyzing $(basename "${BAM_FILE}") ..."
    echo "========================================" >> "${STATS_FILE}"
    echo "Sample: $(basename "${BAM_FILE}")" >> "${STATS_FILE}"
    
    # Output basic mapping rates using samtools flagstat
    samtools flagstat "${BAM_FILE}" >> "${STATS_FILE}"
    
    echo "--- Mapping reads per chromosome ---" >> "${STATS_FILE}"
    samtools idxstats "${BAM_FILE}" >> "${STATS_FILE}"
    echo -e "\n" >> "${STATS_FILE}"
done

echo "=========================================================="
echo " 🎉 All processing complete! Stats saved to ${STATS_FILE}"
echo "=========================================================="