#!/bin/bash

# エラー発生時にスクリプトを停止する安全装置
set -e

# ディレクトリ設定（scriptフォルダから見た相対パス）
RAW_DIR="../rawdata"
FASTP_DIR="../clean_data" # マッピング前の綺麗なデータ用フォルダ

# 出力先フォルダの作成
mkdir -p "$FASTP_DIR"

# スレッド数（HPC環境に合わせて8〜12で設定）
THREADS=12

echo "=========================================================="
echo " 🌟 Starting Fastp processing for all samples in $RAW_DIR "
echo "=========================================================="

# 今回のFGCZのファイル形式「_R1_001.fastq.gz」に合わせて検索
for R1_FILE in "$RAW_DIR"/*_R1_001.fastq.gz; do
    
    # ファイルが存在しない場合はスキップ
    [ -e "$R1_FILE" ] || continue

    # ファイル名から「_R1_001.fastq.gz」を削ってサンプル名（ベース名）を抽出
    BASENAME=$(basename "$R1_FILE")
    SAMPLE="${BASENAME%_R1_001.fastq.gz}"

    # R2のファイルパスを自動生成
    R2_FILE="${RAW_DIR}/${SAMPLE}_R2_001.fastq.gz"

    # 出力ファイルの定義
    OUTPUT_R1="${FASTP_DIR}/${SAMPLE}_R1.clean.fastq.gz"
    OUTPUT_R2="${FASTP_DIR}/${SAMPLE}_R2.clean.fastq.gz"
    REPORT_HTML="${FASTP_DIR}/${SAMPLE}_fastp.html"
    REPORT_JSON="${FASTP_DIR}/${SAMPLE}_fastp.json"

    echo "----------------------------------------------------------"
    echo "▶️ Processing sample: $SAMPLE"

    # ペアとなるR2ファイルが存在するかチェック
    if [[ ! -f "$R2_FILE" ]]; then
        echo "⚠️ Error: Read 2 file ($R2_FILE) not found for sample $SAMPLE! Skipping..."
        continue
    fi

    # ★安全スキップ機能★ 既にクリーンデータが存在し、かつ空ファイル(0バイト)でない場合はスキップ
    if [[ -s "$OUTPUT_R1" && -s "$OUTPUT_R2" ]]; then
        echo "✅ Already processed $SAMPLE, skipping..."
        continue
    fi

    # Fastpの実行
    echo "   Running Fastp..."
    fastp -i "$R1_FILE" -I "$R2_FILE" \
          -o "$OUTPUT_R1" -O "$OUTPUT_R2" \
          --thread "$THREADS" -q 30 \
          --detect_adapter_for_pe \
          --length_required 20 \
          --html "$REPORT_HTML" \
          --json "$REPORT_JSON" \
          2> /dev/null

    echo "🎉 Finished processing: $SAMPLE"
done

echo "=========================================================="
echo " 🚀 All samples processed successfully! "
echo " Cleaned files and reports are saved in: $FASTP_DIR "
echo "=========================================================="