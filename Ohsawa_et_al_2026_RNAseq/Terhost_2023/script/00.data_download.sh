#!/bin/bash

# Exit if any command fails or unset variables are used
set -eu

# Target directory layout tuning
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR/.."

RAW_DIR="./rawdata"
mkdir -p "$RAW_DIR"

THREADS=8
SAMPLE_SHEET="./reference/sample_sheet.xlsx"

if [[ ! -f "$SAMPLE_SHEET" ]]; then
    echo "ERROR: Sample sheet not found at $SAMPLE_SHEET"
    exit 1
fi

echo "=========================================================="
echo " Starting SRA Download and Auto-Rename Pipeline "
echo " Reading metadata from $SAMPLE_SHEET "
echo "=========================================================="

# Create temporary robust python parser for xlsx extraction
cat << 'EOF' > ./reference/tmp_xlsx_parser.py
import sys

def run():
    file_path = sys.argv[1]
    # Primary engine: try high-fidelity pandas layer
    try:
        import pandas as pd
        df = pd.read_excel(file_path)
        for _, row in df.iterrows():
            srr = str(row.iloc[0]).strip()
            new_name = str(row.iloc[4]).strip()
            if srr.startswith("SRR") and new_name and "sample" not in srr:
                print(f"{srr}\t{new_name}")
        return
    except Exception:
        pass

    # Secondary engine: zero-dependency native zip/xml parse fallback
    try:
        import zipfile
        import xml.etree.ElementTree as ET
        with zipfile.ZipFile(file_path) as z:
            strings = []
            if "xl/sharedStrings.xml" in z.namelist():
                s_root = ET.fromstring(z.read("xl/sharedStrings.xml"))
                for si in s_root.findall("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}si"):
                    t = si.find("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}t")
                    strings.append(t.text if t is not None else "")
            
            w_root = ET.fromstring(z.read("xl/worksheets/sheet1.xml"))
            sheet_data = w_root.find("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}sheetData")
            for row in sheet_data.findall("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}row"):
                cells = { "".join([char for char in c.get("r") if not char.isdigit()]): c for c in row.findall("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}c") }
                if "A" in cells and "E" in cells:
                    cA, cE = cells["A"], cells["E"]
                    vA = cA.find("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}v")
                    vE = cE.find("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}v")
                    srr = vA.text if vA is not None else ""
                    new_name = vE.text if vE is not None else ""
                    if cA.get("t") == "s" and srr.isdigit(): srr = strings[int(srr)]
                    if cE.get("t") == "s" and new_name.isdigit(): new_name = strings[int(new_name)]
                    srr, new_name = srr.strip(), new_name.strip()
                    if srr.startswith("SRR") and new_name and "sample" not in srr:
                        print(f"{srr}\t{new_name}")
    except Exception:
        sys.exit(1)

if __name__ == "__main__":
    run()
EOF

# Extract pairs and download loops
python3 ./reference/tmp_xlsx_parser.py "$SAMPLE_SHEET" | while read -r line; do
    if [[ -z "$line" ]]; then continue; fi
    
    srr=$(echo "$line" | cut -f1)
    new_name=$(echo "$line" | cut -f2)

    OUTPUT_FASTQ="${RAW_DIR}/${new_name}.fastq"
    OUTPUT_GZ="${RAW_DIR}/${new_name}.fastq.gz"

    if [[ -f "$OUTPUT_GZ" ]]; then
        echo "Already downloaded: $new_name ($srr), skipping..."
        continue
    fi

    echo "----------------------------------------------------------"
    echo "Downloading $srr -> $new_name ..."

    prefetch "$srr" --max-size 50G
    fasterq-dump "$srr" -O "$RAW_DIR" -e "$THREADS"

    if [[ -f "${RAW_DIR}/${srr}.fastq" ]]; then
        mv "${RAW_DIR}/${srr}.fastq" "$OUTPUT_FASTQ"
        echo "Compressing $OUTPUT_FASTQ to .gz ..."
        gzip "$OUTPUT_FASTQ"
    else
        echo "ERROR: Failed to find downloaded fastq for $srr"
        rm -f ./reference/tmp_xlsx_parser.py
        exit 1
    fi

    sleep 2
    rm -rf "$srr" || true
    echo "Finished: $new_name"
done

# Cleanup temporary parser script
rm -f ./reference/tmp_xlsx_parser.py

echo "=========================================================="
echo " All downloads and renaming completed successfully! "
echo " Saved files can be found in: $RAW_DIR "
echo "=========================================================="
