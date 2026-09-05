library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(stringr)
library(scales)

# ====================================================================
# 0. Setup and Paths
# ====================================================================
# IMPORTANT FOR REPRODUCIBILITY (Method 2 + Method 3):
# Please open the RStudio project file (.Rproj) located in the root directory 
# before running this script, or set your working directory to the 'script/' folder.
# setwd("/Volumes/Shin/paper_making/NLS-mCherry/github/Ohsawa_et_al_2026/colloid_particle_model/script") # Commented out for GitHub submission
data_path <- "../data/"
figure_path <- "../results/"

# List of chromatin structure IDs to be excluded
chromatin_structure_ids <- c(
  "CPX-1610", "CPX-1614", "CPX-2566", "CPX-1611", "CPX-1612", "CPX-1613", # Nucleosomes
  "P53551", # Histone H1
  "P36012"  # CSE4 (Centromeric Histone H3)
)

# ====================================================================
# 1. Data Loading & ID Matching
# ====================================================================
protein_n <- read_tsv(paste0(data_path, "paxDB_4932-WHOLE_ORGANISM-integrated_clean.txt"))
proteome <- read_tsv(paste0(data_path, "yeast_proteome_559292.tsv")) %>% 
  mutate(string_external_id = str_remove(xref_string, ";$"))
proteome_mass <- read_tsv(paste0(data_path, "protein_masses.tsv"))

protein_n <- protein_n %>%
  mutate(string_suffix = str_remove(string_external_id, "^4932\\."))

proteome_long <- proteome %>%
  filter(!is.na(gene_names)) %>%
  separate_rows(gene_names, sep = " ") %>%
  distinct(gene_names, .keep_all = TRUE)

lookup_table <- protein_n %>%
  filter(!string_external_id %in% proteome$string_external_id) %>%
  select(string_suffix) %>%
  distinct() %>%
  inner_join(proteome_long, by = c("string_suffix" = "gene_names"))

matched <- protein_n %>%
  filter(string_external_id %in% proteome$string_external_id) %>%
  left_join(proteome %>% select(string_external_id, accession, protein_name, gene_names, go_c),
            by = "string_external_id")

rescued <- protein_n %>%
  filter(!string_external_id %in% proteome$string_external_id) %>%
  left_join(lookup_table, by = "string_suffix")

# ====================================================================
# 2. Absolute Abundance Calculation
# ====================================================================
M <- 5E-12 # grams, or 5 pg, as per BioNumbers

final_protein_n <- bind_rows(matched, rescued) %>%
  left_join(proteome_mass %>% select(accession, mass), by = "accession") %>%
  mutate(
    MW_average = sum(abundance * 1e-6 * mass, na.rm = TRUE),
    N_total = (M * 6.02e23) / MW_average,
    abundance_n = abundance * 1e-6 * N_total
  ) %>%
  select(accession, protein_name, abundance, go_c, mass, abundance_n, gene_names = string_suffix)

tot_all_proteins <- sum(final_protein_n$abundance_n, na.rm = TRUE)

# ====================================================================
# 3. Complexes Assembly Correction
# ====================================================================
complexes <- read_tsv(paste0(data_path, "20230308_complexPortal_559292.tsv"))[, c(1, 2, 19, 8)]
colnames(complexes) <- c("complex_id", "complex_name", "uniprot_ids_stoichiometry", "go_c")

complexes_long <- complexes %>% 
  filter(!grepl("vacu", go_c, ignore.case = TRUE)) %>% 
  separate_rows(uniprot_ids_stoichiometry, sep = "\\|") %>% 
  mutate(copies = str_extract(uniprot_ids_stoichiometry, "\\((\\d+)\\)") %>% 
           str_remove_all("[()]") %>%
           as.numeric(),
         id = str_remove(uniprot_ids_stoichiometry, "\\(.*\\)")) %>% 
  mutate(copies = ifelse(copies == 0, 1, copies)) %>% 
  filter(!is.na(uniprot_ids_stoichiometry)) %>% 
  inner_join(final_protein_n %>% select(id = accession, abundance_n), by = "id") %>% 
  inner_join(proteome_mass %>% select(id = accession, mass), by = "id") %>% 
  group_by(complex_id, complex_name, go_c) %>%
  reframe(complex_mass = sum(mass * copies, na.rm = TRUE),
          abundance_component = abundance_n / copies, 
          abundance_n = median(abundance_component, na.rm = TRUE)) %>% 
  mutate(comp = ifelse((grepl("cytoplasm", go_c, ignore.case = TRUE) | grepl("cytosol", go_c, ignore.case = TRUE)) & grepl("nucleus", go_c, ignore.case = TRUE), "shared",  
                       ifelse(grepl("nucleus", go_c, ignore.case = TRUE), "nucleus", "cytoplasm")))

# 💥 BUG FIX: Exception handling to prevent Actin (P60010) from being absorbed into chromatin complexes and depleted
proteins_in_complexes <- complexes %>% 
  separate_rows(uniprot_ids_stoichiometry, sep = "\\|") %>% 
  mutate(id = str_remove(uniprot_ids_stoichiometry, "\\(.*\\)")) %>%
  distinct(id) %>%
  filter(id != "P60010")

# ====================================================================
# 4. Homo-oligomers Correction
# ====================================================================
proteome_subunit <- read_tsv(paste0(data_path, "yeast_subunit.tsv")) %>% 
  mutate(norm_factor = ifelse(grepl("dimer", cc_subunit, ignore.case = TRUE), 2,
                              ifelse(grepl("trimer", cc_subunit, ignore.case = TRUE), 3,
                                     ifelse(grepl("tetramer", cc_subunit, ignore.case = TRUE), 4, 1))))

proteins_n_no_complex <- final_protein_n %>% 
  filter(!accession %in% proteins_in_complexes$id) %>% 
  left_join(proteome_subunit %>% select(accession, norm_factor), by = "accession") %>%
  mutate(norm_factor = replace_na(norm_factor, 1),
         abundance_n = abundance_n / norm_factor,
         mass = mass * norm_factor) %>%
  mutate(comp = ifelse((grepl("cytoplasm", go_c, ignore.case = TRUE) | grepl("cytosol", go_c, ignore.case = TRUE)) & grepl("nucleus", go_c, ignore.case = TRUE), "shared",  
                       ifelse(grepl("nucleus", go_c, ignore.case = TRUE), "nucleus", "cytoplasm"))) %>% 
  select(accession, mass, abundance_n, comp, norm_factor, gene_names, protein_name, go_c) %>% 
  filter(!grepl("ribosomal subunit", protein_name, ignore.case = TRUE))

# ====================================================================
# 5. Recombination, Filtering & Empirical Adjustments
# ====================================================================
combined_with_complexes <- proteins_n_no_complex %>% 
  bind_rows(complexes_long %>% 
              distinct(accession = complex_id, mass = complex_mass, abundance_n, comp, protein_name = complex_name, go_c) %>% 
              mutate(norm_factor = 1, gene_names = NA)) %>%
  
  # --- GUIDO'S EMPIRICAL CORRECTIONS (Based on Yeast Data) ---
  # 1. Remove Histones/Nucleosomes entirely (<1% free fraction)
  filter(!accession %in% chromatin_structure_ids) %>% 
  
  # 2. Scale Actin & Tubulin based on empirical soluble fractions
  mutate(
    abundance_n = case_when(
      accession == "P60010" ~ abundance_n * 0.33,
      accession %in% c("CPX-1424", "CPX-1425", "P05217", "P05218", "P02557") ~ abundance_n * 0.10,
      TRUE ~ abundance_n
    )
  )
# -----------------------------------------------------------

val_monomer <- tot_all_proteins
val_monomer_and_complex <- sum(combined_with_complexes$abundance_n, na.rm = TRUE)

val_without_small <- combined_with_complexes %>% 
  filter(mass > 40000) %>% 
  summarise(n = sum(abundance_n, na.rm = TRUE)) %>% 
  pull(n)

val_without_membrane <- combined_with_complexes %>% 
  filter(mass > 40000) %>% 
  filter(!grepl("membrane", go_c, ignore.case = TRUE)) %>% 
  summarise(n = sum(abundance_n, na.rm = TRUE)) %>% 
  pull(n)

# ====================================================================
# 6. Visualization
# ====================================================================
colloidal_df <- data.frame(
  condition = c("monomer", "monomer_and_complex", "without_small_protein", "without_membrane_protein"),
  colloidal_particle = c(val_monomer, val_monomer_and_complex, val_without_small, val_without_membrane)
) %>%
  mutate(particle_rounded = round(colloidal_particle, -5),
         condition = factor(condition, levels = c("monomer", "monomer_and_complex", "without_small_protein", "without_membrane_protein")))

ggplot(colloidal_df) +
  geom_bar(aes(x = condition, y = particle_rounded, fill = condition), stat = "identity", color = "black") +
  geom_text(aes(x = condition, y = particle_rounded / 2, label = paste0(round(particle_rounded / 1e6), " M")), color = "black") +
  scale_y_continuous(labels = label_number(scale = 1e-6, suffix = "M")) +
  scale_fill_manual(values = c("#5ABAB6", "#FFE7FF", "#DDADDA", "#E71DEE")) +
  theme_classic()  

ggsave(paste0(figure_path, "colloid_particle_assumption.pdf"), width = 7, height = 4)

# ====================================================================
# 7. Average Particle Size Calculation & Visualization
# ====================================================================
mw_p_n <- 8.8e4

combined_with_complexes %>% 
  mutate(mass = case_when(protein_name == "60S cytosolic large ribosomal subunit" ~ 2.2e6, TRUE ~ mass)) %>%
  filter(protein_name != "40S cytosolic small ribosomal subunit") %>% 
  mutate(overall_mean = weighted.mean(mass, w = abundance_n, na.rm = TRUE), mw_p_n = mw_p_n) %>%
  ggplot(aes(x = mass, weight = abundance_n)) +
  geom_histogram(bins = 50, alpha = 0.5, position = "identity", fill = "lightblue") +
  geom_segment(aes(x = overall_mean, xend = overall_mean, y = 0, yend = Inf), color = "black", linetype = "dashed", size = 0.5) +
  geom_text(aes(x = overall_mean+4e4, y = Inf, label = paste0("average particle size: ",round(overall_mean/1000), " kDa")), angle = 90, vjust = -0.5, hjust = 1, color = "black", size = 4, check_overlap = TRUE) +
  geom_segment(aes(x = mw_p_n, xend = mw_p_n, y = 0, yend = Inf), color = "magenta", linetype = "dashed", size = 0.5) +
  geom_text(aes(x = mw_p_n, y = Inf, label = paste0("3xFP size: ", round(mw_p_n/1000), " kDa")), angle = 90, vjust = -0.5, hjust = 1, color = "magenta", size = 4, check_overlap = TRUE) +
  scale_x_log10() +
  theme_classic() +
  labs(x = "Particle mass (Da)", y = "Weighted count")

ggsave(paste0(figure_path, "particle_size.pdf"), width = 6, height = 4)

# ====================================================================
# 8. Time-series Proteomics Analysis (GN2019 Replicates)
# ====================================================================
GN2019_proteome_rep <- read.csv(paste0(data_path, "1-s2.0-S0092867419300510-mmc3_clean_rep.csv"))
rep_data <- GN2019_proteome_rep
colnames(rep_data) <- str_remove(colnames(rep_data), "^X")           
colnames(rep_data) <- str_replace(colnames(rep_data), "_rep\\.?0?", "_rep") 

rep_long <- rep_data %>%
  pivot_longer(cols = -string_suffix, names_to = c("time", "rep"), names_pattern = "(.*)h_(.*)") %>%
  mutate(time = paste0(time, "h"))

rep_1h_base <- rep_long %>% filter(time == "1h") %>% select(string_suffix, rep, val_1h = value)
rep_fc <- rep_long %>%
  left_join(rep_1h_base, by = c("string_suffix", "rep")) %>%
  mutate(fc = ifelse(is.na(value / val_1h) | is.infinite(value / val_1h), 1, value / val_1h))

final_protein_time <- final_protein_n %>%
  inner_join(rep_fc, by = c("gene_names" = "string_suffix")) %>%
  mutate(abundance_base = abundance_n, abundance_n = abundance_base * fc) %>%
  select(accession, protein_name, go_c, mass, gene_names, time, rep, abundance_n)

complexes_time <- complexes %>% 
  filter(!grepl("vacu", go_c, ignore.case = TRUE)) %>% 
  separate_rows(uniprot_ids_stoichiometry, sep = "\\|") %>% 
  mutate(copies = as.numeric(str_remove_all(str_extract(uniprot_ids_stoichiometry, "\\((\\d+)\\)"), "[()]")),
         id = str_remove(uniprot_ids_stoichiometry, "\\(.*\\)")) %>% 
  mutate(copies = ifelse(is.na(copies) | copies == 0, 1, copies)) %>% 
  filter(!is.na(uniprot_ids_stoichiometry)) %>% 
  inner_join(final_protein_time %>% select(id = accession, time, rep, abundance_n, mass), by = "id", relationship = "many-to-many") %>% 
  group_by(complex_id, complex_name, go_c, time, rep) %>%
  summarise(complex_mass = sum(mass * copies, na.rm = TRUE), abundance_n = median(abundance_n / copies, na.rm = TRUE), .groups = "drop") %>% 
  mutate(comp = ifelse((grepl("cytoplasm", go_c, ignore.case = TRUE) | grepl("cytosol", go_c, ignore.case = TRUE)) & grepl("nucleus", go_c, ignore.case = TRUE), "shared", 
                       ifelse(grepl("nucleus", go_c, ignore.case = TRUE), "nucleus", 
                              ifelse(grepl("cytoplasm", go_c, ignore.case = TRUE) | grepl("cytosol", go_c, ignore.case = TRUE), "cytoplasm", "other"))))

proteins_n_no_complex_time <- final_protein_time %>% 
  filter(!accession %in% proteins_in_complexes$id) %>% 
  left_join(proteome_subunit %>% select(accession, norm_factor), by = "accession") %>% 
  mutate(norm_factor = replace_na(norm_factor, 1),
         abundance_n = abundance_n / norm_factor,
         mass = mass * norm_factor) %>%
  mutate(comp = ifelse((grepl("cytoplasm", go_c, ignore.case = TRUE) | grepl("cytosol", go_c, ignore.case = TRUE)) & grepl("nucleus", go_c, ignore.case = TRUE), "shared", 
                       ifelse(grepl("nucleus", go_c, ignore.case = TRUE), "nucleus", 
                              ifelse(grepl("cytoplasm", go_c, ignore.case = TRUE) | grepl("cytosol", go_c, ignore.case = TRUE), "cytoplasm", "other")))) %>% 
  filter(!grepl("ribosomal subunit", protein_name, ignore.case = TRUE))

table_quantities_time_rep <- proteins_n_no_complex_time %>% 
  bind_rows(complexes_time %>% select(accession = complex_id, mass = complex_mass, abundance_n, comp, protein_name = complex_name, go_c, time, rep) %>% mutate(norm_factor = 1, gene_names = NA)) %>% 
  
  # --- GUIDO'S EMPIRICAL CORRECTIONS (Time-series) ---
  # 1. Remove Histones/Nucleosomes entirely
  filter(!accession %in% chromatin_structure_ids) %>% 
  
  # 2. Scale Actin & Tubulin based on empirical soluble fractions
  mutate(
    abundance_n = case_when(
      accession == "P60010" ~ abundance_n * 0.33,
      accession %in% c("CPX-1424", "CPX-1425", "P05217", "P05218", "P02557") ~ abundance_n * 0.10,
      TRUE ~ abundance_n
    )
  ) %>%
  # ---------------------------------------------------
filter(
  mass > 40000, 
  !grepl("membrane", go_c, ignore.case = TRUE)
) %>% 
  group_by(time, rep, comp) %>%
  reframe(tot_abundance_n = sum(abundance_n, na.rm = TRUE))


# Statistical Analysis (T-test)
ratio_data_rep <- table_quantities_time_rep %>% 
  group_by(time, rep) %>% 
  summarise( 
    nuc = sum(tot_abundance_n[comp == "nucleus"], na.rm = TRUE), 
    cyto = sum(tot_abundance_n[comp == "cytoplasm"], na.rm = TRUE), 
    .groups = "drop" 
  ) %>% 
  mutate(nc_ratio = 100 * nuc / (nuc + cyto))

ratio_data_rep$time <- factor(ratio_data_rep$time, levels = c("1h", "3h", "5h", "7h"))
base_1h_ratios <- ratio_data_rep$nc_ratio[ratio_data_rep$time == "1h"]

t_test_results <- ratio_data_rep %>%
  filter(time != "1h") %>%
  group_by(time) %>%
  summarise(mean_ratio = mean(nc_ratio), p_value = t.test(nc_ratio, base_1h_ratios)$p.value, .groups = "drop") %>%
  mutate(significance = symnum(p_value, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")))

# Plot Stacked Bar with Stats
plot_stack_summary <- table_quantities_time_rep %>%
  group_by(time, rep) %>%
  mutate(vol_ratio = tot_abundance_n / sum(tot_abundance_n)) %>%
  group_by(time, comp) %>%
  summarise(mean_vol_ratio = mean(vol_ratio), .groups = "drop") %>%
  mutate(comp = factor(comp, levels = c("other", "shared", "nucleus", "cytoplasm")))

label_data <- ratio_data_rep %>%
  group_by(time) %>%
  summarise(mean_ratio = mean(nc_ratio), .groups = "drop") %>%
  left_join(t_test_results %>% select(time, p_value), by = "time") %>%
  mutate(
    p_text = ifelse(is.na(p_value), "", paste0("p = ", signif(p_value, 2))),
    label_text = ifelse(p_text == "", paste0(round(mean_ratio, 1), "%"), paste0(round(mean_ratio, 1), "%\n", p_text))
  )

comp_colors <- c("cytoplasm" = "#5ABAB6", "nucleus" = "#E71DEE", "shared" = "#DDADDA", "other" = "#FFE7FF")

p_stacked_stat <- ggplot(plot_stack_summary, aes(x = time, y = mean_vol_ratio, fill = comp)) +
  geom_bar(stat = "identity", position = "fill", color = "black", width = 0.7) +
  scale_fill_manual(values = comp_colors) +
  geom_text(data = label_data, aes(x = time, y = 1.05, label = label_text), inherit.aes = FALSE, color = "black", fontface = "bold", size = 4.5, lineheight = 0.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), breaks = seq(0, 1, by = 0.25)) +
  labs(y = "Particle Ratio (%)", fill = "Compartment") +
  theme_classic() 
p_stacked_stat
ggsave(paste0(figure_path, "time_series_stacked_with_stats.pdf"), plot = p_stacked_stat, width = 6.5, height = 5)

# ==============================================================================
# 9. Model Plotting (7 PDFs)
# ==============================================================================
# Load external and previous script results required for modeling
nc_ratios_df <- read.csv("../../single_cell_microscopy_analysis/results/nc_ratio_baselines.csv")

# Load variables generated in single_cell_microscopy_analysis_for_paper.R
load("../../single_cell_microscopy_analysis/results/best_Nt_from_NES.RData")
load("../../single_cell_microscopy_analysis/results/n_n.RData")

NCratio_mean_sc_3xFP_NLS    <- nc_ratios_df$value[nc_ratios_df$variable == "NCratio_mean_sc_3xFP_NLS"]
NCratio_mean_sc_3xFP_NES    <- nc_ratios_df$value[nc_ratios_df$variable == "NCratio_mean_sc_3xFP_NES"]
NCratio_mean_sc_3xFP_NLS_ui <- nc_ratios_df$value[nc_ratios_df$variable == "NCratio_mean_sc_3xFP_NLS_ui"]
NCratio_mean_sc_3xFP_NES_ui <- nc_ratios_df$value[nc_ratios_df$variable == "NCratio_mean_sc_3xFP_NES_ui"]

# ------------------------------------------------------------------------------
# 9.1. 3xFP_NLS_model.pdf
# ------------------------------------------------------------------------------
N_t <- sort(colloidal_df$particle_rounded)
n_nls <- seq(0, 5E7, by = 1E5)
r0_NLS <- NCratio_mean_sc_3xFP_NLS_ui
R_v <- 0

df5_model1 <- expand.grid(R_v = R_v, N_t = N_t, n_nls = n_nls) %>% 
  mutate(exp_obs_NC = r0_NLS + n_nls * ((1 - R_v - r0_NLS) / N_t),
         N_t_1e6 = factor(N_t / 1e+06))

ggplot() +
  geom_hline(yintercept = NCratio_mean_sc_3xFP_NLS, linetype = "dashed") +
  geom_line(data = df5_model1, aes(x = n_nls, y = exp_obs_NC, color = N_t_1e6), size = 0.8) +  
  scale_color_manual(values = c("#E71DEE" , "#DDADDA", "#FFE7FF","#5ABAB6")) +
  annotate("text", x = 1.1e5, y = 0.38, label = "italic(phi[v*','*0]) ==~italic(phi[v*','*i]) == 0", parse = TRUE, size = 5, hjust = 0, color = "black") +
  scale_x_log10(limits = c(1e+05, 1e+07)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  labs(x = expression(italic(n) * ", number of 3xFP-NLS"), y = expression(italic(R)[pred] * ", N/C ratio")) +
  theme_classic() 

ggsave(paste0(figure_path, "3xFP_NLS_model.pdf"), width = 4, height = 4)

# ------------------------------------------------------------------------------
# 9.2. 3xFP_NES_model.pdf
# ------------------------------------------------------------------------------
n_nes <- seq(0, 5E7, by = 1E5)
r0_NES <- NCratio_mean_sc_3xFP_NES_ui
R_v <- 0

df6_model2 <- expand.grid(R_v = R_v, N_t = N_t, n_nes = n_nes) %>% 
  mutate(exp_obs_NC = r0_NES * (1 - (n_nes / N_t)),
         N_t_1e6 = factor(N_t / 1e+06))

ggplot() +
  geom_hline(yintercept = NCratio_mean_sc_3xFP_NES, linetype = "dashed") +
  geom_line(data = df6_model2, aes(x = n_nes, y = exp_obs_NC, color = N_t_1e6), size = 0.8) +  
  scale_color_manual(values = c("#E71DEE" , "#DDADDA", "#FFE7FF","#5ABAB6")) +
  annotate("text", x = 1.1e5, y = 0.07, label = "italic(phi[v*','*0]) == ~ italic(phi[v*','*i]) == 0", parse = TRUE, size = 5, hjust = 0, color = "black") +
  scale_x_log10(limits = c(1e+05, 1e+07)) +
  scale_y_continuous(limits = c(0, 0.075)) +
  labs(x = expression(italic(n) * ", number of 3xFP-NES"), y = expression(italic(R)[pred] * ", N/C ratio")) +
  theme_classic() 

ggsave(paste0(figure_path, "3xFP_NES_model.pdf"), width = 4, height = 4)

# ------------------------------------------------------------------------------
# 9.3. 3xFP_NLS_model_vac.pdf
# ------------------------------------------------------------------------------
R_v_seq <- seq(0, 0.4, by = 0.2)

df5_model3 <- expand.grid(R_v = R_v_seq, N_t = N_t, n_nls = n_nls) %>% 
  filter(N_t > 1.4e+07 & N_t < 2.1e+07) %>%
  mutate(exp_obs_NC = r0_NLS + n_nls * ((1 - R_v - r0_NLS) / N_t),
         N_t_1e6 = N_t / 1e+06)

ggplot() +
  geom_hline(yintercept = NCratio_mean_sc_3xFP_NLS, linetype = "dashed") +
  geom_line(data = df5_model3, aes(x = n_nls, y = exp_obs_NC, color = as.factor(R_v)), size = 0.8) +  
  scale_color_manual(values = c("#E71DEE","#FF5778","#FFCD4B"), name = "R_v") +
  annotate("text", x = 1.1e5, y = 0.38, label = paste0("italic(N)[t] == '", round(mean(as.numeric(as.character(df5_model3$N_t_1e6)))), " M'"), parse = TRUE, size = 5, hjust = 0, color = "black") +
  scale_x_log10(limits = c(1e+05, 1e+07)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  labs(x = expression(italic(n) * ", number of 3xFP-NLS"), y = expression(italic(R)[pred] * ", N/C ratio")) +
  theme_classic() 

ggsave(paste0(figure_path, "3xFP_NLS_model_vac.pdf"), width = 4, height = 4)

# ------------------------------------------------------------------------------
# 9.4. 3xFP_NES_model_vac.pdf
# ------------------------------------------------------------------------------
df6_model4 <- expand.grid(R_v = R_v_seq, N_t = N_t, n_nes = n_nes) %>% 
  filter(N_t > 1.4e+07 & N_t < 2.1e+07) %>%
  mutate(exp_obs_NC = r0_NES * (1 - (n_nes / N_t)),
         N_t_1e6 = N_t / 1e+06,
         R_v = factor(R_v))

ggplot() +
  geom_hline(yintercept = NCratio_mean_sc_3xFP_NES, linetype = "dashed") +
  geom_line(data = df6_model4, aes(x = n_nes, y = exp_obs_NC, color = R_v, linewidth = R_v)) +  
  scale_color_manual(values = c("#E71DEE","#FF5778","#FFCD4B")) +
  scale_linewidth_manual(values = c(1.5, 1, 0.5)) +
  annotate("text", x = 1.1e5, y = 0.07, label = paste0("italic(N)[t] == '", round(mean(as.numeric(as.character(df6_model4$N_t_1e6)))), " M'"), parse = TRUE, size = 5, hjust = 0, color = "black") +
  scale_x_log10(limits = c(1e+05, 1e+07)) +
  scale_y_continuous(limits = c(0, 0.075)) +
  labs(x = expression(italic(n) * ", number of 3xFP-NES"), y = expression(italic(R)[pred] * ", N/C ratio")) +
  theme_classic() 

ggsave(paste0(figure_path, "3xFP_NES_model_vac.pdf"), width = 4, height = 4)

# ------------------------------------------------------------------------------
# 9.5. 3xFP_NLS_model_vac_change.pdf
# ------------------------------------------------------------------------------
R_v2_baseline <- 0.1
n_start_vac <- 3e+04
n_end_vac <- 1e+07
N_t_val <- 1.4e+07
target_Rv1_list <- c(0, 0.1, 0.2)
n_nls_vec <- 10^seq(log10(n_start_vac), log10(n_end_vac), length.out = 100)

df5_model5 <- expand.grid(target_Rv1 = target_Rv1_list, n_nls = n_nls_vec, N_t = N_t_val) %>%
  mutate(
    N_t_1e6 = N_t / 1e6,
    R_v1 = R_v2_baseline + (target_Rv1 - R_v2_baseline) * (n_nls - n_start_vac) / (n_end_vac - n_start_vac),
    exp_obs_NC = ((n_nls - n_start_vac) / N_t) * (1 - r0_NLS / (1 - R_v2_baseline)) * (1 - R_v1) + r0_NLS * (1 - R_v1) / (1 - R_v2_baseline)
  )

ggplot() +
  geom_hline(yintercept = NCratio_mean_sc_3xFP_NLS, linetype = "dashed", color = "gray50") +
  geom_line(data = df5_model5, aes(x = n_nls, y = exp_obs_NC, color = as.factor(target_Rv1)), size = 0.8) +  
  scale_color_manual(values = c("#E71DEE","#FF5778","#FFCD4B"), name = bquote(atop(italic(phi)[v*','*i], (italic(n)[i] == .(n_end_vac / 1e6) ~ "M")))) +
  annotate("text", x = 1.1e5, y = 0.35, label = paste0("atop(italic(N)[t] == '", round(N_t_val / 1e6, 1), " M', ", "italic(phi)[v*','*0] == '", R_v2_baseline, "' ~ (italic(n)[0] == '", n_start_vac / 1e6, " M'))"), parse = TRUE, size = 4, hjust = 0, color = "black") +
  scale_x_log10(limits = c(1e+05, 1e+07)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  labs(x = expression(italic(n) * ", number of 3xFP-NLS"), y = expression(italic(R)[pred] * ", N/C ratio")) +
  theme_classic()

ggsave(paste0(figure_path, "3xFP_NLS_model_vac_change.pdf"), width = 4, height = 4)

# ------------------------------------------------------------------------------
# 9.6. 3xFP_NES_model_vac_change.pdf
# ------------------------------------------------------------------------------
n_nes_vec <- 10^seq(log10(n_start_vac), log10(n_end_vac), length.out = 100)

df6_model6 <- expand.grid(target_Rv1 = target_Rv1_list, n_nes = n_nes_vec, N_t = N_t_val) %>%
  mutate(
    N_t_1e6 = N_t / 1e6,
    R_v1 = R_v2_baseline + (target_Rv1 - R_v2_baseline) * (n_nes - n_start_vac) / (n_end_vac - n_start_vac),
    exp_obs_NC = (1 - (n_nes - n_start_vac) / N_t) * (1 - R_v1) * (r0_NES / (1 - R_v2_baseline))
  )

ggplot() +
  geom_hline(yintercept = NCratio_mean_sc_3xFP_NES, linetype = "dashed", color = "gray50") +
  geom_line(data = df6_model6, aes(x = n_nes, y = exp_obs_NC, color = as.factor(target_Rv1)), size = 0.8) +  
  scale_color_manual(values = c("#E71DEE","#FF5778","#FFCD4B"), name = bquote(atop(italic(phi)[v*','*i], (italic(n)[i] == .(n_end_vac / 1e6) ~ "M")))) +
  annotate("text", x = 1.1e5, y = 0.07, label = paste0("atop(italic(N)[t] == '", round(N_t_val / 1e6, 1), " M', ", "italic(phi)[v*','*0] == '", R_v2_baseline, "' ~ (italic(n)[0] == '", n_start_vac / 1e6, " M'))"), parse = TRUE, size = 4, hjust = 0, color = "black") +
  scale_x_log10(limits = c(1e+05, 1e+07)) +
  scale_y_continuous(limits = c(0, 0.075)) +
  labs(x = expression(italic(n) * ", number of 3xFP-NES"), y = expression(italic(R)[pred] * ", N/C ratio")) +
  theme_classic()

ggsave(paste0(figure_path, "3xFP_NES_model_vac_change.pdf"), width = 4, height = 4)

# ------------------------------------------------------------------------------
# 9.7. 3xFP_NLS_model_with_nuc_tension.pdf
# ------------------------------------------------------------------------------
kT <- 4.11e-21
N_t_tension <- best_Nt_from_NES 
n_c_ratio_seq <- seq(r0_NLS, 0.4, by = 0.0001)
sigma_val_mN_m <- c(0, 0.1, 0.5, 1)
cell_vol_m3 <- 50 * 1e-18

df5_model7 <- expand.grid(N_t = N_t_tension, n_c_ratio = n_c_ratio_seq, sigma_val_N_m = sigma_val_mN_m / 1e3) %>% 
  mutate(
    nuc_vol_m3 = cell_vol_m3 * n_c_ratio,
    sigma_val_mN_m = sigma_val_N_m * 1e3,
    n_fp = ((sigma_val_N_m * (1 - n_c_ratio)) / ((kT / (2 * (nuc_vol_m3^(2/3)))) * ((3 / (4 * pi))^(1/3))) + N_t_tension * (n_c_ratio - r0_NLS)) / (1 - r0_NLS)
  )

ggplot() +
  geom_hline(yintercept = NCratio_mean_sc_3xFP_NLS, linetype = "dashed") +
  geom_vline(xintercept = n_n, linetype = "dashed") +
  geom_line(data = df5_model7, aes(x = n_fp, y = n_c_ratio, color = as.factor(sigma_val_mN_m)), size = 0.8) +
  scale_color_manual(values = c("#E71DEE", "#DDADDA", "#FFE7FF", "#5ABAB6"), name = "sigma") +
  scale_x_log10(limits = c(1e+05, 1e+07)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  labs(x = expression(italic(n) * ", number of 3xFP-NLS"), y = expression(italic(R)[pred] * ", N/C ratio")) +
  annotate("text", x = 1.1e5, y = 0.35, label = paste0("atop(italic(N)[t] == '", round(N_t_tension / 1e6, 1), " M', italic(phi)[v*','*0] == italic(phi)[v*','*i] ~ '=' ~ 0)"), parse = TRUE, size = 4, hjust = 0, color = "black") +
  theme_classic() 

ggsave(paste0(figure_path, "3xFP_NLS_model_with_nuc_tension.pdf"), width = 4, height = 4)

# ==============================================================================
print("All 10 plots have been successfully outputted to the results folder!")
# ==============================================================================