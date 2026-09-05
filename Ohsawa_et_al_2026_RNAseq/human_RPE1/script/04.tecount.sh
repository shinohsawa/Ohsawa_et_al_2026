#!/bin/bash
set -e

# Directory settings
BAM_DIR="../mapped_results"
OUT_DIR="../te_counts"
mkdir -p "${OUT_DIR}"

# 🌟 Specify temporary directory
export TMPDIR="/mnt/biol_bc_neurohr_scratch_1/sosawa/tmp_tecount"
mkdir -p "${TMPDIR}"

# Reference files
GTF="../reference/gencode.v44.primary_assembly.annotation.gtf"
TE_GTF="../reference/GRCh38_GENCODE_rmsk_TE.gtf"

THREADS=8

for bam in "${BAM_DIR}"/*_Aligned.sortedByCoord.out.bam; do
    sample=$(basename "${bam}" _Aligned.sortedByCoord.out.bam)
    
    # Automatically skip already completed samples (e.g., day2_r1)!
    if [ -f "${OUT_DIR}/${sample}.cntTable" ]; then
        echo "✅ ${sample} is already finished. Skipping..."
        echo "----------------------------------------"
        continue
    fi

    echo "▶️ Processing ${sample}..."
    
    # 🌟 Pre-sort with samtools before passing to TEcount!
    TEMP_NAME_BAM="${TMPDIR}/${sample}_namesorted.bam"
    echo "   1/2: Running fast name-sorting with samtools (8 threads)..."
    samtools sort -n -@ ${THREADS} -m 2G -T "${TMPDIR}/${sample}_sort_tmp" "${bam}" -o "${TEMP_NAME_BAM}"

    echo "   2/2: Running TEcount..."
    TEcount -b "${TEMP_NAME_BAM}" \
            --GTF "${GTF}" \
            --TE "${TE_GTF}" \
            --format BAM \
            --mode multi \
            --project "${OUT_DIR}/${sample}"

    # Immediately clean up used temporary files
    rm -f "${TEMP_NAME_BAM}"

    echo "🎉 Finished ${sample}"
    echo "----------------------------------------"
done

echo "=========================================="
echo " 🚀 All TE count processing completed successfully!"
echo "=========================================="