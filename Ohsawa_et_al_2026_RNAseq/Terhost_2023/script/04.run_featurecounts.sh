Bash
#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

echo "=========================================================="
echo " Starting featureCounts via Rsubread (SE + Multimap Rescue) "
echo "=========================================================="

Rscript - << 'EOF'

# Load necessary library
library(Rsubread)

# 1. Update the working directory path to the "Terhorst_2023" project
#setwd("/mnt/biol_bc_neurohr_scratch_1/sosawa/NGS/Ohsawa_et_al_2026_RNAseq/Terhorst_2023/")

# Create output directory
dir.create("./gene_counts", showWarnings = FALSE)

# Define file paths
bam_folder <- "./bamfiles/"
bam_files <- list.files(bam_folder, pattern = "\\.bam$", full.names = TRUE)

# 2. Update the GTF file path to the yeast reference genome
gtf <- "./reference/Saccharomyces_cerevisiae.R64-1-1.113.gtf"

# Run featureCounts (Single-end configuration with multi-mapping rescue settings!)
fc_results <- featureCounts(
    files = bam_files,
    annot.ext = gtf,
    isPairedEnd = FALSE,           # 🌟 Set to FALSE for single-end reads
    countMultiMappingReads = TRUE, # 🌟 Rescue multi-mapping reads!
    fraction = TRUE,               # 🌟 Fractionally assign reads mapped to multiple locations to avoid bias!
    allowMultiOverlap = TRUE,      # 🌟 Allow counting of reads that overlap multiple genes in the dense yeast genome!
    isGTFAnnotationFile = TRUE,
    nthreads = 12,
    tmpDir = "/tmp"
)

# Save counts to CSV
write.csv(fc_results$counts, file = "./gene_counts/raw_counts.csv", row.names = TRUE)

# Save summary statistics to CSV
write.csv(fc_results$stat, file = "./gene_counts/raw_counts_summary.csv", row.names = TRUE)

EOF

echo "=========================================================="
echo " featureCounts completed successfully! "
echo " Results are saved in ./gene_counts/ "
echo "=========================================================="