#!/bin/bash
set -e

# ==========================================================
# 04.tecount.sh : selected Drosophila samples (Single-End)
# ==========================================================

# Directory settings
BAM_DIR="../mapped_results"
OUT_DIR="../te_counts"
mkdir -p "${OUT_DIR}"

# 🌟 HPC safety measure: Just like STAR, use local /tmp where FIFO/temporary file creation is 100% guaranteed!
export TMPDIR="/tmp/sosawa_tecount_shi2020"
mkdir -p "${TMPDIR}"

# Drosophila melanogaster Reference Files
GTF="../reference/Drosophila_melanogaster.BDGP6.32.109.gtf"
TE_GTF="../reference/BDGP6_rmsk_TE.gtf"

THREADS=8

echo "=========================================================="
echo " Starting TEcount for selected Drosophila samples "
echo "=========================================================="

# Accurately target successfully mapped BAMs in mapped_results
for bam in "${BAM_DIR}"/*_Aligned.sortedByCoord.out.bam; do
    
    
    FILENAME=$(basename "${bam}")
    SAMPLE="${FILENAME%_Aligned.sortedByCoord.out.bam}"
    
    # 🌟 Inherited from the Human version: Automatically skip already completed samples!
    if [ -f "${OUT_DIR}/${SAMPLE}.cntTable" ]; then
        echo "✅ ${SAMPLE} is already finished. Skipping..."
        echo "----------------------------------------"
        continue
    fi

    echo "▶️ Processing ${SAMPLE}..."
    
    # 🌟 Fast pre-sorting with 8-core samtools power before passing to TEcount!
    TEMP_NAME_BAM="${TMPDIR}/${SAMPLE}_namesorted.bam"
    echo "   1/2: Running fast name-sorting with samtools (8 threads)..."
    # Specify memory limit (-m 2G) and temporary prefix (-T)
    samtools sort -n -@ ${THREADS} -m 2G -T "${TMPDIR}/${SAMPLE}_sort_tmp" "${bam}" -o "${TEMP_NAME_BAM}"

    echo "   2/2: Running TEcount..."
    # 🌟 Since it's already pre-sorted, remove the dangerous --sortByPos flag!
    TEcount -b "${TEMP_NAME_BAM}" \
            --GTF "${GTF}" \
            --TE "${TE_GTF}" \
            --format BAM \
            --mode multi \
            --project "${OUT_DIR}/${SAMPLE}"

    # Immediately clean up used temporary files to protect disk space
    rm -f "${TEMP_NAME_BAM}"

    echo "🎉 Finished ${SAMPLE}"
    echo "----------------------------------------"
done

echo "=========================================================="
echo " 🚀 All TE count processing completed successfully! "
echo "=========================================================="