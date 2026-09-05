#!/bin/bash
set -e

# ==========================================================
# 00.data_download.sh: SRA Data Download, Rename, and Parallel Compression Pipeline
# ==========================================================

# Directory settings by yourself
BASE_DIR=".."
RAW_DIR="${BASE_DIR}/raw_data"
META_DIR="${BASE_DIR}/metadata"
TMP_DIR="/mnt/biol_bc_neurohr_scratch_1/sosawa/tmp_sra"

mkdir -p "${RAW_DIR}" "${META_DIR}" "${TMP_DIR}"

echo "=========================================================="
echo " 1. Fetching metadata (SRP239059)"
echo "=========================================================="
cd "${META_DIR}"

# Fetch RunInfo and create a list of SRR numbers
esearch -db sra -query "SRP239059" | efetch -format runinfo > runinfo.csv
cut -d ',' -f 1 runinfo.csv | grep "^SRR" > srr_list.txt

echo "Number of samples fetched: $(wc -l < srr_list.txt)"
echo "----------------------------------------------------------"

echo "=========================================================="
echo " 2. Downloading SRA data and converting to FASTQ"
echo "=========================================================="
cd "${RAW_DIR}"

# 1. Download SRA files
prefetch --option-file "${META_DIR}/srr_list.txt" --max-size 100G

# 2. FASTQ conversion (🌟 Use -t to specify a large scratch space for temporary files to prevent crashes!)
cat "${META_DIR}/srr_list.txt" | xargs -I {} fasterq-dump {} --split-files -e 8 -t "${TMP_DIR}"

echo "FASTQ conversion completed! Cleaning up the original SRA directories..."
rm -rf ./SRR*/
echo "----------------------------------------------------------"

echo "=========================================================="
echo " 3. Automatic renaming to sample names"
echo "=========================================================="
# 🌟 Generate and execute a Python script on the fly
cat << 'EOF' > rename_fastq.py
import csv
import glob
import os
import re

meta_dir = "../metadata"
sra_result_file = os.path.join(meta_dir, "sra_result.csv")
runinfo_file = os.path.join(meta_dir, "runinfo.csv")

# 1. Create a mapping table from SRX (Experiment ID) -> Sample Name using sra_result.csv
srx_to_name = {}
if os.path.exists(sra_result_file):
    with open(sra_result_file, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            srx = row.get('Experiment Accession')
            title = row.get('Experiment Title', '')
            match = re.search(r':\s*([^;]+);', title)
            if srx and match:
                srx_to_name[srx] = match.group(1).strip()
else:
    print(f"⚠️ Warning: {sra_result_file} not found. Renaming might be skipped.")

# 2. Read SRR (filename) -> SRX mapping from runinfo.csv to determine the final name
srr_to_name = {}
if os.path.exists(runinfo_file):
    with open(runinfo_file, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            srr = row.get('Run')
            srx = row.get('Experiment')
            if srr and srx and srx in srx_to_name:
                srr_to_name[srr] = srx_to_name[srx]

# 3. Batch rename the actual files (Supports _1.fastq / _2.fastq!)
count = 0
for fastq in glob.glob('SRR*.fastq'):
    # Separate SRR123456 and _1.fastq from SRR123456_1.fastq
    match = re.match(r'(SRR\d+)(.*\.fastq)', fastq)
    if match:
        srr, suffix = match.groups()
        if srr in srr_to_name:
            base_name = srr_to_name[srr]
            new_name = f"{base_name}{suffix}"
            
            if os.path.exists(new_name):
                new_name = f"{base_name}_{srr}{suffix}"
                
            print(f"✅ Renamed successfully: [{fastq}] -> {new_name}")
            os.rename(fastq, new_name)
            count += 1
        else:
            print(f"⚠️ Warning: Corresponding name for {srr} not found.")

print(f"\nCompleted: Successfully renamed {count} files in total!")
EOF

python3 rename_fastq.py
rm -f rename_fastq.py
echo "----------------------------------------------------------"

echo "=========================================================="
echo " 4. Parallel compression of FASTQ files (Super fast with 8 cores!)"
echo "=========================================================="
# 🌟 Compress 8 files simultaneously using -P 8!
ls *.fastq | xargs -I {} -P 8 gzip {}

echo "=========================================================="
echo " 🎉 All processes completed successfully! Data is located in raw_data/ "
echo "=========================================================="