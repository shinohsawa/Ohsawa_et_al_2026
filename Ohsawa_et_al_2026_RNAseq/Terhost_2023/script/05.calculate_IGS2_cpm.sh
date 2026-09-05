#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

echo "============================================================"
echo " Calculating endogenous Sc specific CPM for IGS2 region "
echo "============================================================"

# Launch Rscript to execute the code embedded up to the EOF tag
/usr/bin/Rscript - << 'EOF'
library(Rsubread)

# 1. Set paths
setwd("/mnt/biol_bc_neurohr_scratch_1/sosawa/NGS/Ohsawa_et_al_2026_RNAseq/Terhorst_2023/")

cat("1. Running featureCounts for SAF region (IGS2)...\n")
saf_file <- "/mnt/biol_bc_neurohr_scratch_1/sosawa/NGS/Ohsawa_et_al_2026_RNAseq/yeast/reference/IGS2.saf"

if (!file.exists(saf_file)) {
    stop("Error: SAF file not found at ", saf_file)
}

bam_folder <- "./bamfiles/"
bam_files <- list.files(bam_folder, pattern = "\\.bam$", full.names = TRUE)

# Run featureCounts (Single-End & Multi-Mapping Rescue for IGS2)
fc_results <- featureCounts(
    files = bam_files,
    annot.ext = saf_file,
    isPairedEnd = FALSE,            # 🌟 Set to FALSE for single-end reads
    countMultiMappingReads = TRUE,  # 🌟 Essential for repeat regions like IGS2
    allowMultiOverlap = TRUE,
    fraction = TRUE,                # 🌟 Fractionally assign multi-mapping reads
    isGTFAnnotationFile = FALSE,
    nthreads = 12,
    tmpDir = "/tmp"
)

# 2. Clean column names (remove long "_Aligned.sortedByCoord.out.bam" suffix)
raw_IGS2 <- fc_results$counts
colnames(raw_IGS2) <- gsub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", basename(colnames(raw_IGS2)))

cat("2. Calculating endogenous mRNA specific CPM...\n")
# Load the whole gene count data created earlier
sc_counts <- read.csv("./gene_counts/raw_counts.csv", row.names=1)
colnames(sc_counts) <- gsub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", colnames(sc_counts))

# 🌟 CRITICAL: Completely exclude noise genes before calculating the library size (denominator)!
filter_pattern <- "^RDN|^YLR154|^YLR162|^YLR157|^YLR159|^YLR156|^ETS|^ITS|^snR|^SCR1|^LSR1"
sc_counts_filtered <- sc_counts[!grepl(filter_pattern, rownames(sc_counts)), ]

# Use the total read counts of pure mRNA only as the library size (denominator)
sc_lib_sizes <- colSums(sc_counts_filtered)
sc_lib_sizes <- sc_lib_sizes[colnames(raw_IGS2)] # Ensure the sample order matches

# Compute CPM (Counts Per Million)
cpm_matrix <- t(t(raw_IGS2) / sc_lib_sizes) * 1e6
cpm_df <- as.data.frame(cpm_matrix)

# 3. Data formatting and saving
# Add a "region" column and move it to the leftmost position
cpm_df$region <- rownames(cpm_df)
cpm_df <- cpm_df[, c("region", setdiff(colnames(cpm_df), "region"))]

output_file <- "./gene_counts/IGS2_cpm.csv"
write.csv(cpm_df, output_file, row.names=FALSE)

cat("\n=== SUCCESS: IGS2 CPM file saved to", output_file, "===\n")
EOF

echo "All processes completed!"