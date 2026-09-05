#!/bin/bash
set -e

# Create output directory and define paths
OUT_DIR="../TSS_Metagene_Occupancy"
mkdir -p "$OUT_DIR"

BW_DIR="../bigwigs_log2FC"
BED="../reference/Sc_gene_IGS2.bed"
THREADS=12

echo "=========================================================="
echo " 🌟 Metagene Matrix Generation for R (reps.gz)"
echo "=========================================================="

# 1. Check if BED file exists
if [ ! -f "$BED" ]; then
    echo "❌ Error: BED file not found at $BED"
    exit 1
fi

# ==========================
# 1. 3225 (EtOH: 3, EST: 3)
# ==========================
echo "▶️ Processing Strain: 3225 (3 EtOH reps, 3 EST reps)..."
computeMatrix scale-regions -R "$BED" \
  -S "$BW_DIR"/3225_EtOH_1_log2FC.bw "$BW_DIR"/3225_EtOH_2_log2FC.bw "$BW_DIR"/3225_EtOH_3_log2FC.bw \
     "$BW_DIR"/3225_EST_1_log2FC.bw  "$BW_DIR"/3225_EST_2_log2FC.bw  "$BW_DIR"/3225_EST_3_log2FC.bw \
  -b 500 -a 500 --regionBodyLength 1000 --binSize 10 \
  -o "$OUT_DIR/matrix_3225_reps.gz" -p "$THREADS"

# ==========================
# 2. 3229 (EtOH: 3, EST: 3)
# ==========================
echo "▶️ Processing Strain: 3229 (3 EtOH reps, 3 EST reps)..."
computeMatrix scale-regions -R "$BED" \
  -S "$BW_DIR"/3229_EtOH_1_log2FC.bw "$BW_DIR"/3229_EtOH_2_log2FC.bw "$BW_DIR"/3229_EtOH_3_log2FC.bw \
     "$BW_DIR"/3229_EST_1_log2FC.bw  "$BW_DIR"/3229_EST_2_log2FC.bw  "$BW_DIR"/3229_EST_3_log2FC.bw \
  -b 500 -a 500 --regionBodyLength 1000 --binSize 10 \
  -o "$OUT_DIR/matrix_3229_reps.gz" -p "$THREADS"

# ==========================
# 3. 3506 (★EtOH: 2, EST: 3) *Excluding 3506_EtOH_2!
# ==========================
echo "▶️ Processing Strain: 3506 (2 EtOH reps, 3 EST reps)..."
computeMatrix scale-regions -R "$BED" \
  -S "$BW_DIR"/3506_EtOH_1_log2FC.bw "$BW_DIR"/3506_EtOH_3_log2FC.bw \
     "$BW_DIR"/3506_EST_1_log2FC.bw  "$BW_DIR"/3506_EST_2_log2FC.bw  "$BW_DIR"/3506_EST_3_log2FC.bw \
  -b 500 -a 500 --regionBodyLength 1000 --binSize 10 \
  -o "$OUT_DIR/matrix_3506_reps.gz" -p "$THREADS"

echo "=========================================================="
echo " 🎉 SUCCESS: All matrix_..._reps.gz files generated! "
echo "=========================================================="