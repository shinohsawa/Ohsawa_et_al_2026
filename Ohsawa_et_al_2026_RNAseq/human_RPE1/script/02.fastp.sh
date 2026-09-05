#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Directory settings
RAW_DIR="../rawdata"
FASTP_DIR="../fastp"

# Create output directory for Fastp
mkdir -p "$FASTP_DIR"

# Number of threads to use
THREADS=12

echo "=========================================================="
echo " Starting Fastp processing for all samples in $RAW_DIR "
echo "=========================================================="

# 🌟 Loop through the organized _1.fq.gz files
for R1_FILE in "$RAW_DIR"/*_1.fq.gz; do
    
    # Skip if file does not exist
    [ -e "$R1_FILE" ] || continue

    # Extract sample name by removing the suffix
    BASENAME=$(basename "$R1_FILE")
    SAMPLE="${BASENAME%_1.fq.gz}"

    # Define the corresponding Read 2 file path
    R2_FILE="${RAW_DIR}/${SAMPLE}_2.fq.gz"

    # Define output file paths
    OUTPUT_R1="${FASTP_DIR}/${SAMPLE}_1.clean.fastq.gz"
    OUTPUT_R2="${FASTP_DIR}/${SAMPLE}_2.clean.fastq.gz"
    REPORT_HTML="${FASTP_DIR}/${SAMPLE}_fastp.html"
    REPORT_JSON="${FASTP_DIR}/${SAMPLE}_fastp.json"

    echo "----------------------------------------------------------"
    echo "Processing sample: $SAMPLE"

    # Check if the corresponding Read 2 file exists
    if [[ ! -f "$R2_FILE" ]]; then
        echo "  Error: Read 2 file ($R2_FILE) not found for sample $SAMPLE!"
        continue
    fi

    # ★ Skip functionality ★
    # Skip if the clean R1 and R2 files already exist and are not empty
    if [[ -s "$OUTPUT_R1" && -s "$OUTPUT_R2" ]]; then
        echo "  Already processed $SAMPLE, skipping..."
        continue
    fi

    # Run Fastp
    echo "  Running Fastp..."
    fastp -i "$R1_FILE" -I "$R2_FILE" \
          -o "$OUTPUT_R1" -O "$OUTPUT_R2" \
          --thread "$THREADS" -q 30 \
          --detect_adapter_for_pe \
          --length_required 20 \
          --html "$REPORT_HTML" \
          --json "$REPORT_JSON" \
          2> /dev/null

    echo "Finished processing: $SAMPLE"
done

echo "=========================================================="
echo " All samples processed successfully! "
echo " Cleaned files and reports are saved in: $FASTP_DIR "
echo "=========================================================="