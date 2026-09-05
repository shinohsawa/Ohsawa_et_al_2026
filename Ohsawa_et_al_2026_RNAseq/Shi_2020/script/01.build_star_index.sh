#!/bin/bash
# ==============================================================================
# Build STAR Genome Index for Drosophila melanogaster
# ==============================================================================

# Exit immediately if a command exits with a non-zero status
set -e

# ==========================================
# Directory and File Settings
# ==========================================
BASE_DIR=".."
REF_DIR="${BASE_DIR}/reference"
INDEX_DIR="${BASE_DIR}/STAR_index"

# Reference Genome and Annotation
REF_FASTA="${REF_DIR}/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa"
REF_GTF="${REF_DIR}/Drosophila_melanogaster.BDGP6.32.109.gtf"

# Number of threads to use
THREADS=8

# Overhang parameter: Ideally "ReadLength - 1"
# Assuming 50bp reads based on the provided alignment script
OVERHANG=49

# Parameter calculation for small/medium genomes: min(14, log2(GenomeLength)/2 - 1)
# Drosophila genome is ~140MB. log2(140M)/2 - 1 ≈ 12.5 -> We use 12.
SA_INDEX_NBASES=12

# ==========================================
# Step 1: Download Reference Files
# ==========================================
mkdir -p "${REF_DIR}"
cd "${REF_DIR}"

echo "=========================================================="
echo " 1. Checking and Downloading Reference Files"
echo "=========================================================="

# Download and unzip FASTA
if [ ! -f "Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa" ]; then
    echo "▶️ Downloading Drosophila FASTA from Ensembl..."
    curl -O ftp://ftp.ensembl.org/pub/release-109/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa.gz
    echo "▶️ Unzipping FASTA..."
    gunzip Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa.gz
else
    echo "✅ FASTA file already exists."
fi

# Download and unzip GTF
if [ ! -f "Drosophila_melanogaster.BDGP6.32.109.gtf" ]; then
    echo "▶️ Downloading Drosophila GTF from Ensembl..."
    curl -O ftp://ftp.ensembl.org/pub/release-109/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.109.gtf.gz
    echo "▶️ Unzipping GTF..."
    gunzip Drosophila_melanogaster.BDGP6.32.109.gtf.gz
else
    echo "✅ GTF file already exists."
fi

# Return to script directory
cd - > /dev/null


# ==========================================
# Step 2: Build STAR Index
# ==========================================
mkdir -p "${INDEX_DIR}"

echo "=========================================================="
echo " 2. Starting STAR Genome Index Generation for Drosophila"
echo "=========================================================="
echo "FASTA:    ${REF_FASTA}"
echo "GTF:      ${REF_GTF}"
echo "Output:   ${INDEX_DIR}"
echo "Overhang: ${OVERHANG}"
echo "----------------------------------------------------------"

# Run STAR in genome generation mode
STAR --runThreadN "${THREADS}" \
     --runMode genomeGenerate \
     --genomeDir "${INDEX_DIR}" \
     --genomeFastaFiles "${REF_FASTA}" \
     --sjdbGTFfile "${REF_GTF}" \
     --sjdbOverhang "${OVERHANG}" \
     --genomeSAindexNbases "${SA_INDEX_NBASES}"

echo "=========================================================="
echo " 🎉 STAR Genome Index successfully built! "
echo "=========================================================="