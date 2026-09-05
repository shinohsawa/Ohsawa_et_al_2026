# Estimation of Intracellular Colloidal Particle Abundance and Biophysical Modeling of N/C Ratio Scaling

This module is part of the analytical pipeline for **Ohsawa *et al.* (2026)**. It estimates the absolute abundance of intracellular freely diffusing colloidal particles ($N_t$) and their average molecular weight by integrating proteome-wide stoichiometry datasets. Furthermore, it performs biophysical and theoretical simulations of nuclear-to-cytoplasmic (N/C) volume ratio scaling under varying vacuolar volume fractions ($\phi_v$) and nuclear membrane tension ($\sigma$) constraints.

---

## 📁 Directory Structure

    colloid_particle_model/
    ├── README.md               # This document
    ├── script/
    │   └── colloidal_proteome_and_nc_scaling_model.R  # Main simulation script
    ├── data/                   # Input proteome datasets and baseline variables
    └── results/                # Output destination for generated PDF plots

---

## 🛠 System Requirements & R Dependencies

The script was developed and tested using **R (v4.3+)**. To run the analysis, the following R packages are required:

    install.packages(c("dplyr", "tidyr", "readr", "ggplot2", "stringr", "scales", "writexl"))

---

## 🚀 How to Run the Analysis

To ensure complete reproducibility and avoid file path errors, please follow these steps:

1. **Open the Project:** Download or clone the repository and open the **`.Rproj`** file located at the root directory (`Ohsawa_et_al_2026/.Rproj`) using RStudio.
2. **Execute the Script:** Open `script/colloidal_proteome_and_nc_scaling_model.R` and run the entire script.
3. **Relative Paths:** The script strictly uses relative paths (`../data/` and `../results/`) assuming execution within the RStudio Project environment. **Do not modify working directory paths (`setwd`)** unless running outside of RStudio.

---

## 📊 Input Data Description (`data/`)

This directory contains the input datasets required to estimate osmotically active colloidal particles and model N/C ratio scaling. The origins, citations, and filtering criteria for each file are detailed below:

### 1. External Database Files (Proteome & Stoichiometry)

* **`paxDB_4932-WHOLE_ORGANISM-integrated_clean.txt`**
  * **Origin:** Whole-organism protein abundance dataset for *Saccharomyces cerevisiae* downloaded from **PaxDb** (v4.1).
  * **Reference:** [https://doi.org/10.1093/nar/gkaf1066](https://doi.org/10.1093/nar/gkaf1066) / [https://pax-db.org/](https://pax-db.org/)
  * **Usage:** Used alongside protein mass data to calculate absolute copy numbers assuming a total cellular protein mass of **5 pg** (BNID 106225; von der Haar *et al.*, 2002).

* **`20230308_complexPortal_559292.tsv`**
  * **Origin:** Curated multiprotein complex stoichiometry data from the **Complex Portal** database (Release: March 2023, Taxonomy ID: 559292).
  * **Reference:** [https://doi.org/10.1093/nar/gkae1085](https://doi.org/10.1093/nar/gkae1085)
  * **Usage:** Used to correct for protein complex assembly. Complex abundances are calculated as the median abundance of their subunits (scaled by stoichiometry), and individual subunits are subsequently removed from the monomer pool.

* **`yeast_proteome_559292.tsv` & `protein_masses.tsv`**
  * **Origin:** Reference proteome metadata, molecular weights, and Cellular Component annotations from Gene Ontology (**GO_C**) extracted from the **UniProt Knowledgebase (UniProtKB)**.
  * **Reference:** [https://doi.org/10.1093/nar/gkae1010](https://doi.org/10.1093/nar/gkae1010)
  * **Usage:** Provides molecular masses and GO_C annotations used to isolate soluble, freely diffusing particles by excluding membrane-associated ("membrane") and vacuolar ("vacu") components, as well as applying a size exclusion threshold ($\le$ 40 kDa).

* **`yeast_subunit.tsv`**
  * **Origin:** Curated from "Subunit structure" annotations in **UniProtKB** along with specific literature-based stoichiometry adjustments.
  * **Usage:** Adjusts effective particle counts for homo-oligomers (homodimers, trimers, and tetramers) not described in the Complex Portal, with specific stoichiometric corrections applied to nucleosomes and cytoskeletal proteins.

---

## 📈 Expected Output (`results/`)

Running the script will automatically generate and save **9 PDF plots** into the `results/` directory, corresponding to main and supplementary figures in the manuscript:

| File Name | Description / Corresponding Model |
| --- | --- |
| **`colloid_particle_assumption.pdf`** | Stepwise estimation of total intracellular colloidal particle count ($N_t$) across progressive stoichiometric and GO_C filtering steps. |
| **`particle_size.pdf`** | Abundance-weighted histogram distribution of proteome-wide particle masses compared against the 3xFP reporter size. |
| **`3xFP_NLS_model.pdf`** | Theoretical N/C ratio scaling curves for 3xFP-NLS across various hypothetical colloidal capacities ($N_t$). |
| **`3xFP_NES_model.pdf`** | Theoretical N/C ratio scaling curves for 3xFP-NES across various hypothetical colloidal capacities ($N_t$). |
| **`3xFP_NLS_model_vac.pdf`** | NLS scaling simulation incorporating fixed vacuolar volume fractions ($\phi_v = 0 \text{ to } 0.4$). |
| **`3xFP_NES_model_vac.pdf`** | NES scaling simulation incorporating fixed vacuolar volume fractions ($\phi_v = 0 \text{ to } 0.4$). |
| **`3xFP_NLS_model_vac_change.pdf`** | Dynamic NLS scaling model accounting for vacuolar volume expansion upon estradiol induction. |
| **`3xFP_NES_model_vac_change.pdf`** | Dynamic NES scaling model accounting for vacuolar volume expansion upon estradiol induction. |
| **`3xFP_NLS_model_with_nuc_tension.pdf`** | Biophysical scaling model predicting N/C ratio behavior under increasing nuclear envelope tension ($\sigma$). |

---

## 🔬 Author & Citation

Please refer to the main manuscript **Ohsawa *et al.* (2026)** for detailed mathematical derivations, biophysical assumptions, and experimental methods.