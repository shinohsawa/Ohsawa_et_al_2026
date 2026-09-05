#!/bin/bash
# ==============================================================================
# Download Reference and Build STAR Genome Index for Human (GRCh38 / GENCODE v44)
# Project: Ohsawa et al. (2026)
# ==============================================================================

# Exit immediately if a command exits with a non-zero status
set -e

# ==========================================
# Directory and File Settings
# ==========================================
REF_DIR="../reference"
FASTA="${REF_DIR}/GRCh38.primary_assembly.genome.fa"
GTF="${REF_DIR}/gencode.v44.primary_assembly.annotation.gtf"
OUT_INDEX="${REF_DIR}/STAR_index"

# Create directories if they do not exist
mkdir -p "${REF_DIR}"
mkdir -p "${OUT_INDEX}"

# Number of threads to use
THREADS=12

# ==========================================================
# Step 1: Download Human Reference Genome and Annotation
# ==========================================================
FASTA_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz"

echo "=========================================================="
echo " 1. Downloading GENCODE v44 References "
echo "=========================================================="

# Download and unzip FASTA if it doesn't exist
if [ ! -f "${FASTA}" ]; then
    echo "Downloading FASTA file from GENCODE..."
    wget -c -O "${FASTA}.gz" "${FASTA_URL}"
    echo "Unzipping FASTA file..."
    gunzip "${FASTA}.gz"
else
    echo "FASTA file already exists. Skipping download."
fi

# Download and unzip GTF if it doesn't exist
if [ ! -f "${GTF}" ]; then
    echo "Downloading GTF file from GENCODE..."
    wget -c -O "${GTF}.gz" "${GTF_URL}"
    echo "Unzipping GTF file..."
    gunzip "${GTF}.gz"
else
    echo "GTF file already exists. Skipping download."
fi

# ==========================================================
# Step 2: Generate STAR Genome Index
# ==========================================================
echo "=========================================================="
echo " 2. Starting STAR Index Generation (RAM-saving mode) "
echo "=========================================================="
echo "FASTA: ${FASTA}"
echo "GTF:   ${GTF}"
echo "OUT:   ${OUT_INDEX}"
echo "----------------------------------------------------------"

# Clean up any previous failed index attempts
rm -rf "${OUT_INDEX}"/*

# Build STAR index (using --genomeSAsparseD 2 to halve RAM requirements)
STAR --runMode genomeGenerate \
     --runThreadN "${THREADS}" \
     --genomeDir "${OUT_INDEX}" \
     --genomeFastaFiles "${FASTA}" \
     --sjdbGTFfile "${GTF}" \
     --sjdbOverhang 149 \
     --genomeSAsparseD 2 \
     --limitGenomeGenerateRAM 25000000000

echo "=========================================================="
echo " Human STAR Index successfully created in: ${OUT_INDEX} "
echo "=========================================================="