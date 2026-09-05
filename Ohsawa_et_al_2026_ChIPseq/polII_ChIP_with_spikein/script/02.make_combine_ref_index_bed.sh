#!/bin/bash
# Exit immediately if a command exits with a non-zero status
set -e

# ==============================================================================
# 0. Configuration (Check directories and filenames according to your environment)
# ==============================================================================
# Target directory (Use "." if running the script inside the Combined_Genome folder)
REF_DIR="../reference" 
THREADS=12

# --- Input files ---
SC_FASTA="${REF_DIR}/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa"
CG_FASTA="${REF_DIR}/C_glabrata_CBS138_current_chromosomes.fasta"
MCHERRY_FASTA="${REF_DIR}/pGN470withoutScevisiaeorigin.fa"

GTF_FILE="${REF_DIR}/Saccharomyces_cerevisiae.R64-1-1.113.gtf"
IGS2_FILE="${REF_DIR}/IGS2.bed"

# --- Output files ---
COMBINED_FASTA="${REF_DIR}/combined_Sc_Cg_mCherry.fasta"
INDEX_PREFIX="${REF_DIR}/combined_Sc_Cg_mCherry_index"

TEMP_BED="${REF_DIR}/temp_gtf_genes.bed"
# Rename the output BED file to exclude mCherry from the name
FINAL_BED="${REF_DIR}/Sc_and_IGS2.bed" 

echo "=========================================================="
echo " 🌟 STEP 1: Combining FASTA files and creating index"
echo "=========================================================="

echo "▶️ Combining 3 genomes (Sc, Cg, mCherry) into a single FASTA..."
cat "$SC_FASTA" "$CG_FASTA" "$MCHERRY_FASTA" > "$COMBINED_FASTA"
echo "✅ Combination complete: $COMBINED_FASTA"

echo "▶️ Building Bowtie2 index (Threads: $THREADS) ..."
bowtie2-build --threads "$THREADS" "$COMBINED_FASTA" "$INDEX_PREFIX"
echo "✅ Index build complete: $INDEX_PREFIX"


echo "=========================================================="
echo " 🌟 STEP 2: Creating BED file for aggregation (Sc + IGS2)"
echo "=========================================================="

echo "▶️ Extracting gene regions from GTF into a simple 0-based BED format..."
awk 'BEGIN{OFS="\t"} $3=="gene" {
    # Extract gene_id from the 9th column
    match($0, /gene_id "([^"]+)"/, a);
    gene = a[1] ? a[1] : "gene_" NR;
    # Subtract 1 from $4 to make it 0-based
    print $1, $4-1, $5, gene, "0", $7
}' "$GTF_FILE" > "$TEMP_BED"

echo "▶️ Merging IGS2 (replaced with XII) with GTF data and sorting..."
# 1. TEMP_BED (S. cerevisiae)
# 2. IGS2.bed (Replace NC_001144 with XII)
# Combine these two and sort by chromosome name (1st column) and start position (2nd column)
(
    cat "$TEMP_BED"
    sed 's/ref|NC_001144|/XII/' "$IGS2_FILE"
) | sort -k1,1 -k2,2n > "$FINAL_BED"

echo "🧹 Deleting temporary files..."
rm "$TEMP_BED"

echo "=========================================================="
echo " 🎉 SUCCESS: Reference and BED file preparation completed successfully!"
echo "=========================================================="
echo "- Combined FASTA: $COMBINED_FASTA"
echo "- Bowtie2 Index:  $INDEX_PREFIX"
echo "- Final BED:      $FINAL_BED"
echo ""
echo "--- Verifying special elements added to BED ---"
grep -E "IGS2" "$FINAL_BED" || true