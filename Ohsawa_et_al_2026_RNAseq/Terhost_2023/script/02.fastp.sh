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
echo " Starting Fastp processing for SE samples in $RAW_DIR "
echo "=========================================================="


for FASTQ_FILE in "$RAW_DIR"/*.fastq.gz; do
    
    # Skip if file does not exist
    [ -e "$FASTQ_FILE" ] || continue

    # Extract sample name by removing the suffix
    BASENAME=$(basename "$FASTQ_FILE")
    SAMPLE="${BASENAME%.fastq.gz}"

    # Define output file paths 
    OUTPUT_CLEAN="${FASTP_DIR}/${SAMPLE}.clean.fastq.gz"
    REPORT_HTML="${FASTP_DIR}/${SAMPLE}_fastp.html"
    REPORT_JSON="${FASTP_DIR}/${SAMPLE}_fastp.json"

    echo "----------------------------------------------------------"
    echo "Processing sample: $SAMPLE"

   
    if [[ -s "$OUTPUT_CLEAN" ]]; then
        echo "  Already processed $SAMPLE, skipping..."
        continue
    fi

    # Run Fastp (Single-End モード)
    echo "  Running Fastp..."
    fastp -i "$FASTQ_FILE" \
          -o "$OUTPUT_CLEAN" \
          --thread "$THREADS" -q 30 \
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