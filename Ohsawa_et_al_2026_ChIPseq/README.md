# ChIP-seq Bioinformatics Pipeline for Ohsawa et al., 2026

This directory contains the custom ChIP-seq bioinformatics pipeline used to analyze RNA Polymerase II (Pol II) occupancy in **"[Insert your manuscript title here]"** (Ohsawa et al., 2026, submitted to *Nature*).

To achieve highly accurate quantitative comparisons across conditions, this pipeline utilizes a multi-species spike-in normalization strategy mapping reads to a combined genome (*S. cerevisiae* + *C. glabrata* + mCherry plasmid).

## 📊 Data Availability
- **Pol II ChIP-seq Raw Data:** Deposited in NCBI GEO under accession number `[GSEXXXXXX]`.

## 💻 System Requirements & Dependencies
- **Command-line tools:** `fastp`, `Bowtie2`, `samtools`, `deepTools` (bamCoverage, bamCompare, multiBigwigSummary, computeMatrix).
- **R packages:** `ggplot2`, `dplyr`, `tidyr`, `rstatix`, `rtracklayer`, `GenomicRanges`, `ggrepel`.

---

## 🚀 Pipeline Overview & Execution Order

The scripts are configured to be executed sequentially from within the `script/` directory.

### Yeast ChIP-seq Pipeline (Pol II Occupancy)
- `00.data_download.sh` - Downloads and prepares the raw ChIP-seq sequencing data.
- `01.fastp.sh` - Read cleaning and quality control for paired-end ChIP-seq data.
- `02.make_combine_ref_index_bed.sh` - Combines reference genomes (Sc, Cg, mCherry), builds Bowtie2 indices, and extracts gene regions into a customized BED file (including IGS2 modifications).
- `03.map_and_bw.sh` - Bowtie2 mapping, exact spike-in quantification (excluding mitochondrial sequences), and generation of log2FC(IP/Input) BigWig files via `deepTools`.
- `04.general_plot.R` - Aggregates ChIP scores, merges with corresponding RNA-seq data, and generates paired boxplots for Sub-telomeres, Ty1/Ty2 retrotransposons, and iESR clusters.
- `05.make_metagene_matrix.sh` - Generates metagene matrices for TSS/TES regions using `computeMatrix`.
- `06.metagene_plot.R` - Plots normalized Pol II metagene profiles, centering medians by BED regions.
- `07.IGS2_track.R` - Synthesizes average waveforms from replicate BigWig files and plots customized dual-overlay genome tracks for Chromosome XII (including flat-topping for distinct visualization).

---
*For any questions regarding the pipeline, please contact [Shin Ohsawa/ sosawa@ethz.ch] or open an issue.*