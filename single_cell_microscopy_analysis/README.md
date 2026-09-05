# Single-Cell Microscopy Analysis and Biophysical N/C Ratio Modeling

This module constitutes the core image analysis and biophysical modeling pipeline for **Ohsawa *et al.* (2026)**. It processes single-cell morphological and fluorescence data extracted via segmentation pipelines, converts relative fluorescence intensities into absolute molecule counts ($n$) using Coomassie Brilliant Blue (CBB) calibration, and optimizes the theoretical colloidal capacity ($N_t$) governing nuclear-to-cytoplasmic (N/C) volume scaling. Finally, it calculates the physical nuclear envelope tension ($\sigma$) based on volumetric constraints.

---

## 📁 Directory Structure

    single_cell_microscopy_analysis/
    ├── README.md               # This document
    ├── script/
    │   └── single_cell_microscopy_analysis_for_paper.R  # Main analysis script
    ├── rawdata/                # Raw input datasets (segmentation outputs, CBB, vacuole data from "Ohsawa et al., 2026" experiment)
    └── results/                # Output destination for generated plots and shared RData variables

---

## 🛠 System Requirements & R Dependencies

The script was developed and tested using **R (v4.3+)**. To run the analysis, the following R packages are required:

    install.packages(c("dplyr", "tidyr", "readr", "ggplot2", "stringr", "scales", "writexl", "readxl"))

---

## 🚀 How to Run the Analysis

To ensure complete reproducibility and avoid file path errors, please follow these steps:

1. **Open the Project:** Download or clone the repository and open the **`.Rproj`** file located at the root directory (`Ohsawa_et_al_2026/.Rproj`) using RStudio.
2. **Execute the Script:** Open the R script and run it entirely.
3. **Relative Paths:** The script strictly uses relative paths (`./rawdata/` and `./results/`) assuming execution within the RStudio Project environment. **Do not modify working directory paths (`setwd`)** unless running outside of RStudio.

---

## 📊 Input Data Description (`rawdata/`)

This directory must contain the following empirical datasets generated from our experimental pipelines:

### 1. Microscopy Segmentation Outputs (Replicates 1-3)
* **`rep1_251015/`, `rep2_251016/`, `rep3_251017/`**
  * Contains single-cell geometric features (area, volume) and fluorescence intensity metrics extracted via image segmentation tools (Cell-ACDC).
  * Key files: `AllPos_acdc_output_cell.csv`, `AllPos_acdc_output_nuclear.csv`, and `position_basename.csv` (for genotype and induction status mapping).

### 2. Absolute Calibration & Cellular Parameters
* **`250917_SO_qunatification_of_CBB.xlsx`**
  * **Description:** CBB staining quantification data used to compute the universal conversion factor, mapping arbitrary fluorescence units to absolute 3xFP molecule counts per cell.
* **`summary_df_vac.rda`**
  * **Description:** Pre-calculated vacuolar volume fractions ($\phi_v$) under uninduced and induced states for each strain.

---

## 📈 Expected Output (`results/`)

Running the script will populate the `results/` directory with the following outputs. 
*Note: The generated `.RData` files are essential inputs for the downstream `colloid_particle_model/` pipeline.*

### Shared Variables (Passed to downstream pipelines)
| File Name | Description |
| --- | --- |
| **`n_n.RData`** | The calculated absolute baseline copy number of 3xFP molecules in a single uninduced cell. |
| **`best_Nt_from_NES.RData`** | The optimized intracellular colloidal capacity ($N_t$) fitted from the 3xFP-NES scaling data. |
| **`SourceData_SingleCell_Microscopy.xlsx`** | The complete compiled single-cell dataset with geometry, fluorescence, absolute copy numbers ($n$), and membrane tension ($\sigma$) values. |

### Generated PDF Plots (Main & Supplementary Figures)
| File Name | Description / Corresponding Model |
| --- | --- |
| **`3xFP_NLS_model_paper.pdf`** | Best-fit theoretical scaling curve (LOESS optimization) for the 3xFP-NLS system over single-cell data points. |
| **`3xFP_NES_model_for_paper.pdf`** | Best-fit theoretical scaling curve (LOESS optimization) for the 3xFP-NES system over single-cell data points. |
| **`Expected_Nt_paper.pdf`** | Boxplot comparing the distribution of expected $N_t$ calculated independently from individual NLS and NES cells. |
| **`3xFP_NLS_membrane_tension_scatter_plot_for_paper.pdf`** | Scatter plot showing calculated nuclear envelope tension ($\sigma$, mN/m) mapped against absolute 3xFP-NLS abundance ($n$). |
| **`3xFP_NLS_membrane_tension_box_plot_for_paper.pdf`** | Boxplot comparing structural nuclear envelope tension ($\sigma$) between uninduced and induced 3xFP-NLS states. |

---

## 🔬 Author & Citation

Please refer to the main manuscript **Ohsawa *et al.* (2026)** for detailed mathematical derivations, image acquisition settings, and physical scaling assumptions.