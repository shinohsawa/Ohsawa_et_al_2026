#!/bin/bash
# Stop script on error
set -e

# ==============================================================================
# 0. Configuration
# ==============================================================================
THREADS=12
ENV_NAME="chipseq_env"

# Directory settings
REF_DIR="../reference"
CLEAN_DIR="../clean_data"
BAM_DIR="../sorted_bams"
BW_CPM_DIR="../bigwigs_cpm"
BW_COMP_DIR="../bigwigs_log2FC"
QUANT_DIR="../spikein_quantification" # For quantification data output

# Reference settings (Combined genome)
INDEX_PREFIX="${REF_DIR}/combined_Sc_Cg_mCherry_index"
BED_FILE="${REF_DIR}/Sc_and_IGS2.bed"

mkdir -p "$BAM_DIR" "$BW_CPM_DIR" "$BW_COMP_DIR" "$QUANT_DIR"

echo "=========================================================="
echo " 🌟 STEP 0: Conda Environment Check & Setup "
echo "=========================================================="
CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"

if conda env list | grep -q "^${ENV_NAME} "; then
    echo "✅ Conda environment '${ENV_NAME}' is ready."
else
    echo "❌ Error: Environment '${ENV_NAME}' not found."
    exit 1
fi
conda activate "$ENV_NAME"

echo "=========================================================="
echo " 🌟 STEP 1: Bowtie2 Index Check (Build Skipped) "
echo "=========================================================="
if [ -f "${INDEX_PREFIX}.1.bt2" ]; then
    echo "✅ Bowtie2 index found: ${INDEX_PREFIX}"
else
    echo "❌ Error: Bowtie2 index not found at ${INDEX_PREFIX}"
    exit 1
fi

echo "=========================================================="
echo " 🌟 STEP 2: Mapping & Sorting (Sc + Cg + mCherry) "
echo "=========================================================="
for R1_FILE in "$CLEAN_DIR"/*_R1.clean.fastq.gz; do
    [ -e "$R1_FILE" ] || continue
    
    BASENAME=$(basename "$R1_FILE")
    SAMPLE="${BASENAME%_R1.clean.fastq.gz}"
    R2_FILE="${CLEAN_DIR}/${SAMPLE}_R2.clean.fastq.gz"
    FINAL_BAM="${BAM_DIR}/${SAMPLE}_sorted.bam"

    if [ -f "${FINAL_BAM}.bai" ]; then
        echo "✅ Already mapped and sorted: $SAMPLE"
        continue
    fi

    echo "▶️ Mapping & Sorting sample: $SAMPLE"
    bowtie2 -p "$THREADS" -x "$INDEX_PREFIX" -1 "$R1_FILE" -2 "$R2_FILE" 2> "${BAM_DIR}/${SAMPLE}_bowtie2.log" | \
    samtools sort -@ "$THREADS" -o "$FINAL_BAM" -
    
    samtools index "$FINAL_BAM"
    echo "🎉 Finished Mapping: $SAMPLE"
done

echo "=========================================================="
echo " 🌟 STEP 3: Spike-in Quantification (MT excluded) "
echo "=========================================================="
QUANT_FILE="${QUANT_DIR}/scale_factors.txt"

# Create header
echo -e "Sample\tIP_Sc\tIP_Cg\tIP_Plasmid(2089-4630)\tInput_Sc\tInput_Cg\tInput_Plasmid(2089-4630)\tOccupancy_Ratio" > "$QUANT_FILE"

# Regex for exact matching of Roman numerals I to XVI (Ensembl nuclear chromosome names)
SC_NUCLEAR_REGEX="^(I|II|III|IV|V|VI|VII|VIII|IX|X|XI|XII|XIII|XIV|XV|XVI)$"

for IP_BAM in "$BAM_DIR"/IP_*_sorted.bam; do
    [ -e "$IP_BAM" ] || continue
    
    BASENAME=$(basename "$IP_BAM")
    CORE_NAME=${BASENAME#IP_}
    CORE_NAME=${CORE_NAME%_sorted.bam}
    
    INPUT_BAM="${BAM_DIR}/input_${CORE_NAME}_sorted.bam"
    
    if [ ! -f "$INPUT_BAM" ]; then
        echo "⚠️ Warning: Input not found for $CORE_NAME. Skipping."
        continue
    fi
    
    # Cg and plasmid aggregation
    ip_cg=$(samtools idxstats "$IP_BAM" | awk -F '\t' '$1 ~ /Chr/ {sum += $3} END {print sum+0}')
    in_cg=$(samtools idxstats "$INPUT_BAM" | awk -F '\t' '$1 ~ /Chr/ {sum += $3} END {print sum+0}')
    
    ip_plasmid=$(samtools view -c "$IP_BAM" "pGN470withoutScevisiaeorigin:2089-4630")
    in_plasmid=$(samtools view -c "$INPUT_BAM" "pGN470withoutScevisiaeorigin:2089-4630")
    
    # Sc aggregation (Exclude mitochondrial MT, count only nuclear genome I-XVI)
    ip_sc=$(samtools idxstats "$IP_BAM" | awk -F '\t' -v regex="$SC_NUCLEAR_REGEX" '$1 ~ regex {sum += $3} END {print sum+0}')
    in_sc=$(samtools idxstats "$INPUT_BAM" | awk -F '\t' -v regex="$SC_NUCLEAR_REGEX" '$1 ~ regex {sum += $3} END {print sum+0}')
    
    # Calculate Ratio
    if [ "$ip_cg" -eq 0 ] || [ "$in_cg" -eq 0 ]; then
        ratio="NA"
        echo "⚠️ $CORE_NAME: C. glabrata read count is 0."
    else
        ratio=$(awk -v ip_sc="$ip_sc" -v ip_cg="$ip_cg" -v ip_plasmid="$ip_plasmid" \
                    -v in_sc="$in_sc" -v in_cg="$in_cg" -v in_plasmid="$in_plasmid" \
                    'BEGIN { printf "%.5f", ((ip_sc + ip_plasmid) / ip_cg) / ((in_sc + in_plasmid) / in_cg) }')
    fi
    
    # Write results
    echo -e "${CORE_NAME}\t${ip_sc}\t${ip_cg}\t${ip_plasmid}\t${in_sc}\t${in_cg}\t${in_plasmid}\t${ratio}" >> "$QUANT_FILE"
    echo "✅ Calculation complete: $CORE_NAME (Ratio = $ratio)"
done

echo "✅ Saved quantification data to: $QUANT_FILE"

echo "=========================================================="
echo " 🌟 STEP 4: CPM BigWig Generation (Sc Relative Only) "
echo "=========================================================="
for BAM in "$BAM_DIR"/*_sorted.bam; do
    [ -e "$BAM" ] || continue
    BASENAME=$(basename "$BAM")
    SAMPLE="${BASENAME%_sorted.bam}"
    BW_OUT="${BW_CPM_DIR}/${SAMPLE}_CPM.bw"

    if [ -f "$BW_OUT" ]; then
        echo "✅ Already created CPM BigWig: $SAMPLE"
        continue
    fi

    echo "▶️ Generating CPM BigWig: $SAMPLE (Ignoring Cg and mCherry for CPM scale)"
    # Automatically get chromosomes to ignore for CPM calculation (Cg and mCherry)
    IGNORE_CHRS=$(samtools idxstats "$BAM" | awk '$1 ~ /Chr|pGN470/ {print $1}')
    
    bamCoverage -b "$BAM" -o "$BW_OUT" \
                --normalizeUsing CPM \
                --ignoreForNormalization $IGNORE_CHRS \
                --centerReads --binSize 10 -p "$THREADS"
done

echo "=========================================================="
echo " 🌟 STEP 5: IP vs Input (log2FC) BigWig Generation "
echo "=========================================================="
for IP_BAM in "$BAM_DIR"/IP_*_sorted.bam; do
    [ -e "$IP_BAM" ] || continue

    BASENAME=$(basename "$IP_BAM")
    SAMPLE_CORE="${BASENAME#IP_}"
    SAMPLE_CORE="${SAMPLE_CORE%_sorted.bam}"
    
    INPUT_BAM="${BAM_DIR}/input_${SAMPLE_CORE}_sorted.bam"
    COMP_BW="${BW_COMP_DIR}/${SAMPLE_CORE}_log2FC.bw"

    if [[ ! -f "$INPUT_BAM" ]]; then
        echo "⚠️ Warning: Input BAM not found for $SAMPLE_CORE! Skipping."
        continue
    fi

    if [ -f "$COMP_BW" ]; then
        echo "✅ Already calculated log2FC: $SAMPLE_CORE"
        continue
    fi

    echo "▶️ Calculating log2FC (IP/Input): $SAMPLE_CORE"
    IGNORE_CHRS=$(samtools idxstats "$IP_BAM" | awk '$1 ~ /Chr|pGN470/ {print $1}')

    bamCompare -b1 "$IP_BAM" -b2 "$INPUT_BAM" -o "$COMP_BW" \
               --scaleFactorsMethod readCount \
               --ignoreForNormalization $IGNORE_CHRS \
               --operation log2 --binSize 10 -p "$THREADS"
done

echo "=========================================================="
echo " 🌟 STEP 6: Summary Matrix & Tab File Generation"
echo "=========================================================="
SUMMARY_DIR="../summary_from_log2FC_bigwig"
mkdir -p "$SUMMARY_DIR"

SUMMARY_NPZ="${SUMMARY_DIR}/ChIP_Pol2_summary.npz"
SUMMARY_TAB="${SUMMARY_DIR}/ChIP_Pol2_gene_scores.tab"

if [ -f "$SUMMARY_TAB" ] && [ -f "$SUMMARY_NPZ" ]; then
    echo "✅ Already generated summary files in ${SUMMARY_DIR}"
else
    shopt -s nullglob
    bw_files=("$BW_COMP_DIR"/*_log2FC.bw)
    shopt -u nullglob

    if [ ${#bw_files[@]} -eq 0 ]; then
        echo "⚠️ Warning: No log2FC BigWig files found. Skipping summary."
    else
        echo "▶️ Running multiBigwigSummary..."
        multiBigwigSummary BED-file \
            -b "${bw_files[@]}" \
            -o "$SUMMARY_NPZ" \
            --outRawCounts "$SUMMARY_TAB" \
            --BED "$BED_FILE" \
            -p "$THREADS"
            
        echo "🎉 Finished generating summary files!"
    fi
fi

echo "=========================================================="
echo " 🚀 SUCCESS: All processes completed! "
echo "=========================================================="