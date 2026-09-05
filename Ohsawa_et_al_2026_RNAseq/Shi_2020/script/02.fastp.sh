#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Directory settings 
RAW_DIR="../raw_data"
FASTP_DIR="../clean_data"
REPORT_DIR="${FASTP_DIR}/reports"

# Create output directories for Fastp and reports
mkdir -p "$FASTP_DIR" "$REPORT_DIR"

# Number of threads to use (Set to 8 threads for HPC)
THREADS=8

echo "=========================================================="
echo " Starting Fastp processing for selected samples in $RAW_DIR "
echo " Target: LC1, LC2, S4, S8 (Single-End mode)"
echo "=========================================================="

# 🌟  Loop restricted to files starting with LC1, LC2, S4, or S8
for RAW_FILE in "$RAW_DIR"/LC[12]_*.fastq.gz "$RAW_DIR"/S[48]_*.fastq.gz; do
    
    # Skip if file does not exist (Safety mechanism)
    [ -e "$RAW_FILE" ] || continue

    # Extract sample name by removing the suffix (e.g., LC1_1.fastq.gz -> LC1_1)
    BASENAME=$(basename "$RAW_FILE")
    SAMPLE="${BASENAME%.fastq.gz}"

    # Define output file paths (Only one file since it's single-end!)
    OUTPUT_FILE="${FASTP_DIR}/${SAMPLE}.clean.fastq.gz"
    REPORT_HTML="${REPORT_DIR}/${SAMPLE}_fastp.html"
    REPORT_JSON="${REPORT_DIR}/${SAMPLE}_fastp.json"

    echo "----------------------------------------------------------"
    echo "Processing sample: $SAMPLE"

    # ★ Skip functionality (Safety guard perfectly matching the yeast/human version!) ★
    # Skip if the clean file already exists and is not empty
    if [[ -s "$OUTPUT_FILE" ]]; then
        echo "  Already processed $SAMPLE, skipping..."
        continue
    fi

    # Run Fastp (Single-end mode: specifying only -i and -o)
    echo "  Running Fastp..."
    fastp -i "$RAW_FILE" \
          -o "$OUTPUT_FILE" \
          --thread "$THREADS" -q 30 \
          --length_required 20 \
          --html "$REPORT_HTML" \
          --json "$REPORT_JSON" \
          2> /dev/null

    echo "Finished processing: $SAMPLE"
done

echo "=========================================================="
echo " Selected samples processed successfully! "
echo " Cleaned files and reports are saved in: $FASTP_DIR "
echo "=========================================================="