library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rstatix)
library(ggpubr)
library(readxl)
library(ggrepel)

# setwd("/Volumes/biol_bc_neurohr_scratch_1/sosawa/NGS/Ohsawa_et_al_2026_ChIPseq/polII_ChIP_with_spikein/results/")

# ==========================================
# pol II occupancy from all reads
# ==========================================

# 1. Load data
load("../../../Ohsawa_et_al_2026_RNAseq/yeast/yeast_DESeq2_results/A_df_final.RData")
head(df_final)

# ==============================================================
# 1. Load and format ChIP-seq score file (.tab)
# ==============================================================
chip_tab_path <- "../summary_from_log2FC_bigwig/ChIP_Pol2_gene_scores.tab"

# Read data (headers contain "'" so we handle them)
chip_data <- read.delim(chip_tab_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# Remove extra quotes from column names
colnames(chip_data) <- gsub("'", "", colnames(chip_data))

# Ensure the first column is named "chr" (sometimes it is "#chr", etc.)
colnames(chip_data)[1] <- "chr"

# ==============================================================
# 2. Load BED file containing gene names
# ==============================================================
bed_path <- "../reference/Sc_and_IGS2.bed" 
bed_data <- read.table(bed_path, header = FALSE, stringsAsFactors = FALSE, sep = "\t")

# Set BED column names (4th column is the extracted gene_name)
colnames(bed_data) <- c("chr", "start", "end", "gene", "score", "strand")

# ==============================================================
# 3. Merge gene names into ChIP data using coordinates as keys!
# ==============================================================
# Join using coordinates from BED file and chip_data
chip_with_genes <- inner_join(bed_data %>% dplyr::select(chr, start, end, gene), 
                              chip_data, 
                              by = c("chr", "start", "end"))

# --- Function: Create PCA plot with the exact same design as RNA-seq ---
create_rnaseq_style_pca <- function(data, cols_to_use, plot_title) {
  # 1. Create PCA matrix (remove NAs, transpose so rows=samples, cols=genes)
  pca_mat <- data %>%
    dplyr::select(all_of(cols_to_use)) %>%
    drop_na() %>%
    t()
  
  # 2. Run Principal Component Analysis (prcomp)
  pca_res <- prcomp(pca_mat, center = TRUE, scale. = FALSE)
  
  # 3. Calculate Proportion of Variance for each principal component (round to integer: equivalent to RNA-seq's round(100 * attr...))
  percentVar <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)))
  
  # 4. Create dataframe for ggplot and fix factor levels
  pcaData <- as.data.frame(pca_res$x) %>%
    mutate(Sample_Name = rownames(.)) %>%
    mutate(
      strain    = sapply(strsplit(Sample_Name, "_"), `[`, 1),
      treatment = sapply(strsplit(Sample_Name, "_"), `[`, 2),
      replicate = sapply(strsplit(Sample_Name, "_"), `[`, 3)
    ) %>%
    mutate(
      # Fix order to be exactly the same as RNA-seq!
      strain    = factor(strain, levels = c("3225", "3229", "3506")),
      treatment = factor(treatment, levels = c("EtOH", "EST"))
    )
  
  # 5. Create plot with design identical to the RNA-seq script!
  p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = strain, shape = treatment)) +
    geom_point(size = 5) +
    # * Uncomment the line below to add small replicate numbers
    # geom_text_repel(aes(label = replicate), size = 3.5, show.legend = FALSE) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    scale_color_manual(values = c("3225" = "#3853A3", "3229" = "gray", "3506" = "orange")) +
    scale_shape_manual(values = c("EtOH" = 16, "EST" = 17)) +
    ggtitle(plot_title) +
    theme_test() +  # ★ Same theme as RNA-seq!
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      legend.title = element_blank(), # Remove legend title for a cleaner look
      legend.position = "right",
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 12)
    )
  
  return(p)
}

# --- 1. Get all sample columns (including outlier 3506_EtOH_2) ---
all_sample_cols <- names(chip_with_genes)[5:ncol(chip_with_genes)]

# --- 2. Get clean sample columns excluding the outlier ---
clean_sample_cols <- all_sample_cols[all_sample_cols != "3506_EtOH_2_log2FC.bw"]

# --- 3. Create PCA and save as PDF (6x6 inch size, same as RNA-seq) ---
cat("▶️ Generating PCA Plot (Before Outlier Removal)...\n")
p_pca_all <- create_rnaseq_style_pca(chip_with_genes, all_sample_cols, "ChIP-seq PCA (All Samples)")
ggsave("PCA_with_all_sample.pdf", plot = p_pca_all, width = 6, height = 6)

cat("▶️ Generating PCA Plot (After Outlier Removal)...\n")
p_pca_clean <- create_rnaseq_style_pca(chip_with_genes, clean_sample_cols, "ChIP-seq PCA (without 3506_EtOH_2)")
ggsave("PCA_without_3xFP-NESmock_rep2.pdf", plot = p_pca_clean, width = 6, height = 6)

# ==============================================================
# ★ Set exclusion rule for the 8 subsequent plots ★
# ==============================================================
# Finalize sample_cols as clean_sample_cols here!
sample_cols <- clean_sample_cols

cat("✅ Outlier 3506_EtOH_2 removed. Ready for downstream plotting!\n")


chip_long <- chip_with_genes %>%
  rowwise() %>%
  # Calculate the mean of all samples per gene as ChIP_baseMean
  mutate(ChIP_baseMean = mean(c_across(all_of(sample_cols)), na.rm = TRUE)) %>%
  ungroup() %>%
  # Remove unnecessary coordinate columns and convert to long format (Tidy Data)
  dplyr::select(gene, ChIP_baseMean, all_of(sample_cols)) %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "Sample",
    values_to = "ChIP_Value"
  ) %>%
  # Extract Strain, Treatment (EST/EtOH), and Replicate (Rep) from sample names
  mutate(
    Strain    = sapply(strsplit(Sample, "_"), `[`, 1),
    Treatment = sapply(strsplit(Sample, "_"), `[`, 2),
    Replicate = sapply(strsplit(Sample, "_"), `[`, 3)
  )

# ==============================================================
# 2. Run "unpaired t-test" in bulk by gene and strain
# ==============================================================
chip_stats <- chip_long %>%
  group_by(gene, ChIP_baseMean, Strain) %>%
  # ▼ Changed here! Set paired = FALSE for independent sample comparison ▼
  t_test(ChIP_Value ~ Treatment, paired = FALSE, alternative = "two.sided", var.equal = FALSE) %>%
  ungroup()

# ==============================================================
# 3. Calculate exact [Mean EST - Mean EtOH] for each strain,
#    and normalize (shift median to 0) across BED regions simultaneously
# ==============================================================
chip_means <- chip_long %>%
  group_by(gene, Strain, Treatment) %>%
  # Means are calculated correctly even with uneven sample numbers (e.g., 3 EST, 2 EtOH)
  summarise(mean_val = mean(ChIP_Value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Treatment, values_from = mean_val) %>%
  
  # 1. First, calculate pure subtraction (EST - EtOH)
  mutate(Raw_ChIP_log2FC = EST - EtOH) %>%
  
  # 2. ★ Group by Strain here and shift the median to 0 ★
  group_by(Strain) %>%
  mutate(
    ChIP_log2FC = Raw_ChIP_log2FC - median(Raw_ChIP_log2FC, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  
  # Extract only necessary columns (Raw_ChIP_log2FC is no longer needed)
  dplyr::select(gene, Strain, ChIP_log2FC)

# ==============================================================
# 4. Merge statistical results and format Contrast names to match RNA-seq
# ==============================================================
chip_formatted <- chip_stats %>%
  # Merge accurate log2FC
  left_join(chip_means, by = c("gene", "Strain")) %>%
  # Calculate multiple testing correction (padj) for each strain
  group_by(Strain) %>%
  mutate(ChIP_padj = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  # Perfectly match Contrast names from RNA-seq data (df_final)
  mutate(
    Contrast = case_when(
      Strain == "3225" ~ "A3225_EST_vs_EtOH",
      Strain == "3229" ~ "A3229_EST_vs_EtOH",
      Strain == "3506" ~ "A3506_EST_vs_EtOH",
      TRUE ~ NA_character_
    )
  ) %>%
  # Select and rename necessary columns cleanly
  dplyr::select(
    gene, 
    Contrast, 
    ChIP_baseMean, 
    ChIP_log2FC, 
    ChIP_pvalue = p, 
    ChIP_padj
  )

# ==============================================================
# 5. Merge with RNA-seq data (df_final) by gene & Contrast!
# ==============================================================
df_combined_omits <- left_join(df_final, chip_formatted, by = c("gene", "Contrast"))

# Verify
head(df_combined_omits %>% 
       dplyr::select(gene, Contrast, significant, log2FoldChange, padj, ChIP_baseMean, ChIP_log2FC, ChIP_padj))

## retrotransposon
GO0032197 <- read_xlsx("../../../Ohsawa_et_al_2026_RNAseq/yeast/reference/GO0032197_retrotransposition.xlsx")
unique_GO <- GO0032197 %>%
  distinct(gene, .keep_all = TRUE)

df_combined_omits <- left_join(df_combined_omits, unique_GO, by = "gene")
df_combined_omits$GO0032197[is.na(df_combined_omits$GO0032197)] <- "NO"
df_combined_omits$GO0032197 <- factor(df_combined_omits$GO0032197, levels = c("NO", "Ty1", "Ty2", "Ty3", "Ty4", "Ty5", "other"))
df_combined_omits <- df_combined_omits %>% arrange(GO0032197)

# ==============================================================
# 6. Verify final results
# ==============================================================
# Display top rows to check if RNA-seq and ChIP-seq values align correctly per Contrast
head(df_combined_omits %>% 
       dplyr::select(gene, Contrast, significant, log2FoldChange, padj, ChIP_baseMean, ChIP_log2FC, ChIP_padj))

# scatter plot for correlation between RNA-seq vs ChIP-seq
ggplot(df_combined_omits, aes(x=log2FoldChange,y=ChIP_log2FC))+
  geom_smooth(method = "lm", color = "lightblue")+
  geom_point(alpha = 0.2, aes(color=significant))+
  scale_color_manual(values = c("blue","gray","red"))+
  facet_wrap(~Contrast)+
  theme_classic()
ggsave(paste0("scatterplot_ChiPseq_vs_RNA_seq.pdf"), width = 8, height = 4)

### MA plot
label_data_TPM <- df_combined_omits %>%
  dplyr::count(Contrast,significant) %>%
  filter(significant %in% c("up", "down")) %>%
  mutate(
    label = paste0(significant, ": ", n),
    x = max(log10(df_combined_omits$meanTPM_EtOH + 1), na.rm = TRUE) - 2,
    y = ifelse(significant == "up",
               max(df_combined_omits$log2FoldChange, na.rm = TRUE) - 1,
               min(df_combined_omits$log2FoldChange, na.rm = TRUE) + 1),
    color = ifelse(significant == "up", "red", "blue")
  )

# 1. Define background data (adjust range/colors here if needed)
bg_data <- data.frame(
  xmin = c(0.5, 1.5, 2.5),
  xmax = c(1.5, 2.5, 3.5),
  fill_color = c("lightblue", "lightgray", "lightpink")
)

ggplot(df_combined_omits, aes(x = log10(meanTPM_EtOH + 1), y = ChIP_log2FC)) + 
  geom_rect(
    data = bg_data,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill_color),
    alpha = 0.3,          # Transparency (0.3~0.5 recommended so points remain visible)
    inherit.aes = FALSE   # Do not inherit x, y from the main data
  ) +
  # Use the literal color names from the dataframe for fill
  scale_fill_identity() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  
  # --- 1. Scatter plot of all data ---
  geom_point(alpha = 0.5, size = 1, aes(color = significant)) +
  geom_smooth(method = "lm", color = "lightblue")+
  # Color settings for scatter plot points
  scale_color_manual(values = c("blue", "gray", "red")) +
  
  theme_classic() +
  facet_wrap(vars(Contrast)) +
  labs(x = "log10(meanTPM_EtOH  + 1)", y = "log2FC[estradiol/EtOH]")

ggsave(paste0("MA_plot_TPM.pdf"), width = 8, height = 4)

# 1. Perform statistical test first (rstatix)
df_binned <- df_combined_omits %>%
  # Calculate log10(TPM + 1)
  mutate(log10_TPM = log10(meanTPM_EtOH + 1)) %>%
  
  # Filter by specified range (0.5 < x < 4.0)
  filter(log10_TPM > 0.5 & log10_TPM < 4.0) %>%
  
  # Create bins (intervals) in increments of 0.5
  mutate(TPM_bin = cut(
    log10_TPM, 
    breaks = seq(0.5, 3.5, by = 1), # From 0.5 to 4.0 in 0.5 increments
    include.lowest = TRUE # Include lowest boundary
  )) %>%
  
  # Remove NAs that failed binning (just in case)
  filter(!is.na(TPM_bin))

# Group by bin (TPM_bin) and perform t-test against A3229 baseline
stat.test <- df_binned %>%
  group_by(TPM_bin) %>%
  t_test(ChIP_log2FC ~ Contrast, ref.group = "A3229_EST_vs_EtOH") %>%
  adjust_pvalue(method = "BH") %>% # Apply correction if needed (e.g., "BH")
  add_significance("p") %>%
  # Calculate positions on the plot (dodge matches boxplot widths)
  add_xy_position(fun = "max", x = "TPM_bin", dodge = 0.8) %>%
  mutate(label = paste0("p = ", signif(p.adj, 2)))

# Verify
print(stat.test)

# 2. Create Plot
ggplot(df_binned, aes(x = TPM_bin, y = ChIP_log2FC, fill = Contrast)) +
  # === Background bands (using annotate) ===
  # X-axis is a Factor, so coordinates are 1, 2, 3... from the left
  
  # 1st bin (0.5 < TPM < 1.5) -> x coords 0.5 ~ 1.5
  annotate("rect", xmin = 0.55, xmax = 1.45, ymin = -Inf, ymax = Inf, 
           fill = "lightblue", alpha = 0.3) +
  
  # 2nd bin (1.5 < TPM < 2.5) -> x coords 1.5 ~ 2.5
  annotate("rect", xmin = 1.55, xmax = 2.45, ymin = -Inf, ymax = Inf, 
           fill = "lightgray", alpha = 0.3) +
  
  # 3rd bin (2.5 < TPM < 3.5) -> x coords 2.5 ~ 3.5
  annotate("rect", xmin = 2.55, xmax = 3.45, ymin = -Inf, ymax = Inf, 
           fill = "lightpink", alpha = 0.3) +
  # Boxplot (match width and position_dodge for a clean look)
  geom_boxplot(outlier.size = 1, outlier.color = "gray", outlier.alpha = 0.7, alpha = 0.7, width = 0.7, position = position_dodge(0.8)) +
  
  # Y-axis adjustment (if necessary)
  coord_cartesian(ylim = c(-3, 4)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  scale_fill_manual(values = c("#3853A3", "gray", "orange")) +
  theme_classic() +
  
  labs(
    x = "Basal Expression Level (log10[TPM_EtOH + 1])",
    y = "log2FC[estradiol/EtOH]"
  ) +
  
  # 3. Display pre-calculated P-values (instead of stat_compare_means)
  stat_pvalue_manual(
    stat.test, 
    label = "label", 
    tip.length = 0.01,
    hide.ns = TRUE  # Hide non-significant (ns) values
  )
ggsave(paste0("Boxplot_TPMbinned.pdf"), width = 6, height = 4)

### t.test between Contrasts
# ==============================================================
# Subtelomere Boxplot (Paired Test Version)
# ==============================================================
df_combined_omits_subtel <- subset(df_combined_omits, df_combined_omits$subtel == "subtel")

# --- 1. Extract common genes to properly create pairs ---
df_subtel_paired <- df_combined_omits_subtel %>%
  dplyr::select(gene, Contrast, ChIP_log2FC) %>%
  tidyr::pivot_wider(names_from = Contrast, values_from = ChIP_log2FC) %>%
  na.omit() %>%
  tidyr::pivot_longer(cols = -gene, names_to = "Contrast", values_to = "ChIP_log2FC") %>%
  left_join(df_combined_omits_subtel %>% distinct(gene, Contrast, .keep_all = TRUE) %>% dplyr::select(-ChIP_log2FC), by = c("gene", "Contrast")) %>%
  mutate(Contrast = factor(Contrast, levels = c("A3225_EST_vs_EtOH", "A3229_EST_vs_EtOH", "A3506_EST_vs_EtOH"))) %>%
  arrange(Contrast, gene) # ★ Super important for proper pairing

# --- 2. Statistical Test ---
df_combined_omits_allsubtel.stat.test <- df_subtel_paired %>%
  rstatix::t_test(ChIP_log2FC ~ Contrast, paired = TRUE, ref.group = "A3229_EST_vs_EtOH") %>%
  adjust_pvalue(method = "none") %>%
  add_significance() %>%
  add_xy_position(fun = "max", x = "Contrast", dodge = 0.9) %>%
  mutate(label = paste0("p = ", signif(p.adj, 2)))

# --- 3. Create Plot ---
ggplot(df_subtel_paired, aes(x = Contrast, y = ChIP_log2FC)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  geom_boxplot(alpha = 0.7, aes(fill = Contrast), outlier.color = "gray", outlier.alpha = 0.7) + 
  scale_fill_manual(values=c("#3853A3", "gray", "orange")) +
  theme_classic() +
  stat_pvalue_manual(
    data = df_combined_omits_allsubtel.stat.test,
    label = "label",
    tip.length = 0.01,
    bracket.size = 0.5
  ) +
  labs(y = "log2FC[induced/mock]", x = NULL)
ggsave("Boxplot_subtel.pdf", width = 4, height = 4.5)


# ==============================================================
# Retrotransposon (Ty1/Ty2) Boxplot (Paired Test Version)
# ==============================================================
df_combined_omits_GO0032197 <- subset(df_combined_omits, df_combined_omits$GO0032197 != "NO" & df_combined_omits$GO0032197 != "other")
df_combined_omits_GO0032197_ty1_ty2 <- subset(df_combined_omits, df_combined_omits$GO0032197 %in% c("Ty1", "Ty2"))

# --- 1. Extract common genes to properly create pairs ---
df_ty_paired <- df_combined_omits_GO0032197_ty1_ty2 %>%
  dplyr::select(gene, GO0032197, Contrast, ChIP_log2FC) %>%
  tidyr::pivot_wider(names_from = Contrast, values_from = ChIP_log2FC) %>%
  na.omit() %>%
  tidyr::pivot_longer(cols = -c(gene, GO0032197), names_to = "Contrast", values_to = "ChIP_log2FC") %>%
  left_join(df_combined_omits_GO0032197_ty1_ty2 %>% distinct(gene, Contrast, .keep_all = TRUE) %>% dplyr::select(-ChIP_log2FC), by = c("gene", "GO0032197", "Contrast")) %>%
  mutate(GO0032197 = factor(GO0032197, levels = c("Ty1", "Ty2"))) %>%
  mutate(Contrast = factor(Contrast, levels = c("A3225_EST_vs_EtOH", "A3229_EST_vs_EtOH", "A3506_EST_vs_EtOH"))) %>%
  arrange(GO0032197, Contrast, gene)

# --- 2. Statistical Test (added paired = TRUE) ---
df_combined_omits_GO0032197_ty1_ty2.stat.test <- df_ty_paired %>%
  group_by(GO0032197) %>%
  rstatix::t_test(ChIP_log2FC ~ Contrast, paired = TRUE, ref.group = "A3229_EST_vs_EtOH") %>%
  adjust_pvalue(method = "none") %>%
  add_significance() %>%
  add_xy_position(fun = "max", x = "Contrast", dodge = 0.9) %>%
  mutate(label = paste0("p = ", signif(p.adj, 2)))

# --- 3. Create Plot ---
ggplot(df_ty_paired, aes(x = Contrast, y = ChIP_log2FC, fill = Contrast)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  geom_boxplot( outlier.color = "gray", outlier.alpha = 0.7, alpha = 0.7, width = 0.7, position = position_dodge(0.8)) +
  stat_pvalue_manual(data = df_combined_omits_GO0032197_ty1_ty2.stat.test, label = "label", tip.length = 0.01) +
  scale_fill_manual(values=c("#3853A3","gray", "orange")) +
  facet_wrap(~GO0032197) +
  theme_classic() +
  labs(y = "log2FC[estradiol/EtOH]")
ggsave("Boxplot_retrotranspon.pdf", width = 5, height = 4)


# ==============================================================
# iESR Boxplot (Paired Test Version)
# ==============================================================
iESR <- read_xlsx("../../../Ohsawa_et_al_2026_RNAseq/yeast/reference/iESR_cluster_ESR_clusters_UPDATED_2017_from_SGD.xlsx")
unique_GO <- iESR %>% distinct(gene, .keep_all = TRUE)

df_combined_omits <- left_join(df_combined_omits, unique_GO, by = "gene")
df_combined_omits$ESR[is.na(df_combined_omits$ESR)] <- "NO"
df_combined_omits_iESR <- subset(df_combined_omits, df_combined_omits$ESR != "NO")

# --- 1. Extract common genes to properly create pairs ---
df_iesr_paired <- df_combined_omits_iESR %>%
  dplyr::select(gene, ESR, Contrast, ChIP_log2FC) %>%
  tidyr::pivot_wider(names_from = Contrast, values_from = ChIP_log2FC) %>%
  na.omit() %>%
  tidyr::pivot_longer(cols = -c(gene, ESR), names_to = "Contrast", values_to = "ChIP_log2FC") %>%
  left_join(df_combined_omits_iESR %>% distinct(gene, Contrast, .keep_all = TRUE) %>% dplyr::select(-ChIP_log2FC), by = c("gene", "ESR", "Contrast")) %>%
  mutate(Contrast = factor(Contrast, levels = c("A3225_EST_vs_EtOH", "A3229_EST_vs_EtOH", "A3506_EST_vs_EtOH"))) %>%
  arrange(ESR, Contrast, gene)

# --- 2. Statistical Test (added paired = TRUE) ---
df_combined_omits_iESR.stat.test <- df_iesr_paired %>%
  group_by(ESR) %>%
  rstatix::t_test(ChIP_log2FC ~ Contrast, paired = TRUE, ref.group = "A3229_EST_vs_EtOH") %>%
  adjust_pvalue(method = "none") %>%
  add_significance() %>%
  add_xy_position(fun = "max", x = "Contrast", dodge = 0.9) %>%
  mutate(label = paste0("p = ", signif(p.adj, 2)))

# --- 3. Create Plot ---
ggplot(df_iesr_paired, aes(x = Contrast, y = ChIP_log2FC)) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(alpha = 0.7, aes(fill = Contrast), outlier.color = "gray", outlier.alpha = 0.7) +
  stat_pvalue_manual(data = df_combined_omits_iESR.stat.test, label = "label", tip.length = 0.01) +
  scale_fill_manual(values=c("#3853A3", "lightgray", "orange")) +
  theme_classic() +
  facet_wrap(~ESR) +
  labs(y = "log2FC[estradiol/EtOH]")
ggsave("Boxplot_iESR.pdf", width = 5, height = 4)


# ==============================================================
# 1. Normalize by shifting median to 0 for each sample using all gene values
# ==============================================================
target_genes <- c("IGS2-1", "IGS2-2")

# Get score columns (exclude 3506_EtOH_2)
all_score_cols <- names(chip_with_genes)[grepl("log2FC\\.bw$", names(chip_with_genes))]
score_cols <- all_score_cols[all_score_cols != "3506_EtOH_2_log2FC.bw"]

# ★ Key point! Subtract the overall median for each column (sample) ★
chip_normalized <- chip_with_genes %>%
  mutate(across(all_of(score_cols), ~ . - median(., na.rm = TRUE)))

# ==============================================================
# 2. Extract IGS2-1 and IGS2-2 data and convert to long (Tidy) format
# ==============================================================
df_igs2_norm_long <- chip_normalized %>%
  # Extract IGS2 from normalized data
  filter(gene %in% target_genes) %>%
  dplyr::select(gene, all_of(score_cols)) %>%
  pivot_longer(
    cols = all_of(score_cols),
    names_to = "Sample",
    values_to = "ChIP_Value"
  ) %>%
  mutate(
    Strain    = sapply(strsplit(Sample, "_"), `[`, 1),
    Treatment = sapply(strsplit(Sample, "_"), `[`, 2),
    Replicate = sapply(strsplit(Sample, "_"), `[`, 3)
  ) %>%
  mutate(Treatment = factor(Treatment, levels = c("EtOH", "EST")))

# ==============================================================
# 3. T-test and P-value coordinate calculation on normalized data
# ==============================================================
stat.test.norm <- df_igs2_norm_long %>%
  group_by(gene, Strain) %>%
  t_test(ChIP_Value ~ Treatment, paired = FALSE, var.equal = FALSE) %>%
  add_significance() %>%
  add_xy_position(x = "Strain", dodge = 0.8) %>%
  mutate(label = paste0("p = ", signif(p, 2)))

print(stat.test.norm)

# ==============================================================
# 4. Create boxplot (normalized data)
# ==============================================================
p_igs2_norm <- ggplot(df_igs2_norm_long, aes(x = Strain, y = ChIP_Value, fill = Treatment)) +
  
  # Boxplot
  geom_boxplot(
    position = position_dodge(0.8), 
    width = 0.6, 
    alpha = 0.7, 
    outlier.shape = NA, 
    color = "black"
  ) +
  
  # Dots
  geom_point(
    aes(group = Treatment), 
    position = position_jitterdodge(jitter.width = 0, dodge.width = 0.8), 
    show.legend = FALSE, 
    alpha = 0.7, 
    size = 2,
    color = "gray"
  ) +
  
  # P-values
  stat_pvalue_manual(
    stat.test.norm,
    label = "label",
    tip.length = 0.01,
    bracket.nudge.y = 0.05
  ) +
  
  # Design settings
  scale_fill_manual(values = c("EtOH" = "gray", "EST" = "firebrick")) +
  facet_wrap(~gene) + 
  theme_classic() +
  labs(
    x = "Strain",
    y = "Normalized ChIP-seq signal"
  ) +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold", size = 12)
  )

print(p_igs2_norm)

# Save as PDF
ggsave("Boxplot_IGS2.pdf", plot = p_igs2_norm, width = 5, height = 4)

#========================================
# absolute pol II occupancy with spikein
#=======================================

# 1. Load data
df <- read.delim("../spikein_quantification/scale_factors.txt", header = TRUE)

# 2. Split Sample column to extract Strain, Treatment, Replicate
df <- df %>%
  mutate(
    Strain = sapply(strsplit(Sample, "_"), `[`, 1),
    Treatment = sapply(strsplit(Sample, "_"), `[`, 2),
    Replicate = sapply(strsplit(Sample, "_"), `[`, 3) # Get Replicate numbers (1,2,3)
  ) %>%
  mutate(
    Strain = factor(Strain, levels = c("3225", "3229", "3506"), labels = c("3xFP-NLS", "empty", "3xFP-NES"))
  )

# 3. Convert Long data to Wide, calculate EST / EtOH ratio per replicate
df_ratio <- df %>%
  select(Strain, Treatment, Replicate, Occupancy_Ratio) %>%
  pivot_wider(names_from = Treatment, values_from = Occupancy_Ratio) %>%
  mutate(Relative_Ratio = EST / EtOH)

# 4. One-sample t-test (test if Ratio is significantly different from 1.0 (mu = 1))
stat_test_one_sample <- df_ratio %>%
  group_by(Strain) %>%
  t_test(Relative_Ratio ~ 1, mu = 1.0) %>%
  add_significance() %>%
  # Shift y.position up slightly to make asterisks visible
  add_xy_position(x = "Strain", fun = "max") %>%
  mutate(y.position = max(df_ratio$Relative_Ratio)*0.98,
         p.value = paste0("p = ",signif(p, 2)))

print(stat_test_one_sample)

# 5. Point plot (3 Ratio points + Mean/SE + Baseline 1.0)
p <- ggplot(df_ratio, aes(x = Strain, y = Relative_Ratio)) +
  
  # Draw theoretical baseline (no change = 1.0) with a dashed line
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "gray40", size = 0.8) +
  
  # Draw error bars and mean points (Mean ± SE)
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1, size = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 5, aes(colour = Strain)) +
  
  # Draw individual replicate data points (3 points) with slight transparency
  # geom_jitter(width = 0.1, size = 3, color = "gray40", alpha = 0.7) +
  scale_color_manual(values = c("#3853a3", "#bebebe", "#ffa500"))+
  # Y-axis settings (start from 0)
  scale_y_continuous(limits = c(0.7, max(df_ratio$Relative_Ratio) ), expand = expansion(mult = c(0, 0))) +
  
  # Add 1-sample t-test p-values (asterisks)
  stat_pvalue_manual(stat_test_one_sample, label = "p.value", tip.length = 0) +
  
  # Theme and label settings
  theme_classic() +
  labs(
    x = "Strain",
    y = "Pol II occupancy ratio (induced/mock)"
  ) 
print(p)

# Save as PDF
ggsave("./paired_relative_ratio.pdf", plot = p, width = 4, height = 4)