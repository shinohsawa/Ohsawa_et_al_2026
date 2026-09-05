#!/bin/bash
# ==============================================================================
# Build STAR Genome Index for Saccharomyces cerevisiae
# Optimized for Single-End 40bp reads
# ==============================================================================

set -e

# Directory and File Settings
REF_DIR="../reference"
INDEX_DIR="${REF_DIR}/star_index"

# Note: Please ensure this matches the actual name of your FASTA file
FASTA="${REF_DIR}/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa"
GTF="${REF_DIR}/Saccharomyces_cerevisiae.R64-1-1.113.gtf"

THREADS=12

# Overhang parameter: Ideally "ReadLength - 1"
# Since the reads are 40bp, this is set to 39. 
# (If your reads are actually 50bp or 150bp, change this to 49 or 149 respectively)
OVERHANG=39

# Parameter calculation for small genomes: min(14, log2(GenomeLength)/2 - 1)
# S. cerevisiae genome is ~12MB. log2(12M)/2 - 1 ≈ 10.7 -> We use 10.
SA_INDEX_NBASES=10

# Create output directory for the index
mkdir -p "$INDEX_DIR"

echo "=========================================================="
echo " Starting STAR Genome Index Generation "
echo "=========================================================="
echo "FASTA:    $FASTA"
echo "GTF:      $GTF"
echo "Output:   $INDEX_DIR"
echo "Overhang: $OVERHANG (optimized for 40bp reads)"
echo "----------------------------------------------------------"

# Run STAR in genome generation mode
STAR --runThreadN "$THREADS" \
     --runMode genomeGenerate \
     --genomeDir "$INDEX_DIR" \
     --genomeFastaFiles "$FASTA" \
     --sjdbGTFfile "$GTF" \
     --sjdbOverhang "$OVERHANG" \
     --genomeSAindexNbases "$SA_INDEX_NBASES"

echo "=========================================================="
echo " STAR Genome Index successfully built! "
echo "=========================================================="