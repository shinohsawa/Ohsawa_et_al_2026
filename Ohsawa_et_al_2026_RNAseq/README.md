# RNA-seq Bioinformatics Pipeline for Ohsawa et al., 2026

This directory contains the custom RNA-seq bioinformatics pipelines used in **"[Insert your manuscript title here]"** (Ohsawa et al., 2026, submitted to *Nature*).

The pipelines cover RNA-seq analysis across multiple datasets and model systems (*Saccharomyces cerevisiae*, *Drosophila melanogaster*, and Human RPE1 cells). The analyses focus on sub-telomeric gene expression, transposable element (TE) quantification, and the calculation of absolute mRNA concentrations normalized by cell volume.

## 📊 Data Availability
- **Yeast (Own data) & Human RPE1 RNA-seq Data:** Deposited in NCBI GEO under accession number `[GSEXXXXXX]`.
- **Drosophila RNA-seq Data:** Publicly available dataset from Shi et al., 2020 (SRA: `SRP239059`).
- **Yeast cdc28-13 RNA-seq Data:** Publicly available dataset from Terhost et al., 2023 (SRA: `SRR22926215–SRR22926293`).

## 💻 System Requirements & Dependencies
- **Command-line tools:** `fastp`, `STAR`, `samtools`, `featureCounts` (Subread), `TEcount` (TEtranscripts).
- **Python / R dependencies:** `Python 3`, `DESeq2`, `ggplot2`, `dplyr`, `rstatix`, `biomaRt`, `openxlsx`, `ggpubr`.

---

## 🚀 Pipeline Overview & Execution Order

The analyses are divided into four independent sub-directories. All scripts are configured to run from within their respective `script/` directories.

### 1. Yeast RNA-seq Pipeline (`yeast/script/`)
Core pipeline for mapping yeast reads (including *C. albicans* spike-in), running featureCounts, and calculating absolute mRNA concentrations incorporating cell volume.

- `00.data_download.sh` - Downloads/prepares raw sequencing data.
- `02.fastp.sh` - Quality control and adapter trimming.
- `03.generate_bamfiles.sh` - Maps reads to the combined yeast genome using STAR.
- `04.run_featurecounts.sh` - Read counting using `featureCounts`.
- `05.split_counts_by_genome.py` - Python script to separate counts between *S. cerevisiae* and spike-in *C. albicans*.
- `06.calculate_IGS2_cpm.sh` - Calculates endogenous mRNA-specific CPM for the IGS2 region.
- `07.plot.R` / `07.plot_original.R` - Core DESeq2 analysis, PCA, and differential expression visualizations.
- `08.IGS2_plot_for_D.R` - Specific visualizations for IGS2 expression.
- `09.global_analysis.R` - Calculates absolute mRNA concentration (`[mRNA] = Counts(K) / fL`) using cell volume data.

### 2. Yeast cdc28-13 Pipeline (`Terhost_2023/script/`)
Re-analysis of public dataset (Terhost et al., 2023) to validate conserved transcription responses.

- `00.data_download.sh` - Downloads raw SRA data (`SRR22926215–SRR22926293`).
- `01.build_star_index.sh` - Builds STAR index for the reference genome.
- `02.fastp.sh` - Quality control.
- `03.generate_bamfiles.sh` - STAR alignment.
- `04.run_featurecounts.sh` - Gene quantification.
- `05.calculate_IGS2_cpm.sh` - IGS2 CPM calculation for comparison.
- `06.plot.R` - DESeq2 analysis and plotting.

### 3. Drosophila RNA-seq Pipeline (`Shi_2020/script/`)
Analyzes *D. melanogaster* RNA-seq data (Shi et al., 2020) with a specific focus on Transposable Elements (TEs).

- `00.data_download.sh` - Downloads SRA data (`SRP239059`) and automatically renames files using metadata.
- `01.build_star_index.sh` - Builds STAR index for the *Drosophila* genome (BDGP6).
- `02.fastp.sh` - Quality control (Single-end mode).
- `03.generate_bamfiles.sh` - Maps reads using TEtranscripts-optimized STAR parameters.
- `04.tecount.sh` - Pre-sorts BAMs by name and quantifies TE/gene expression using `TEcount`.
- `05.deseq2_analysis.R` - Differential expression analysis, PCA, and TE class-specific volcano/boxplots.

### 4. Human RPE1 RNA-seq Pipeline (`human_RPE1/script/`)
Analyzes RNA-seq data from Human RPE1 cells to observe conserved sub-telomeric and TE responses.

- `00.data_download.sh` - Data preparation.
- `01.build_star_index_human.sh` - Builds STAR index for the Human genome (hg38/GENCODE v44).
- `02.fastp.sh` - Read cleaning (Paired-end).
- `03.generate_bamfiles.sh` - STAR alignment with TE-specific configurations.
- `04.tecount.sh` - TE and gene quantification.
- `05.deseq2_analysis.R` - Runs DESeq2, extracts sub-telomeric coordinates dynamically via `biomaRt`, and plots expression comparisons.

---
*For any questions regarding the pipeline, please contact [Shin Ohsawa/sosawa@ethz.ch] or open an issue.*