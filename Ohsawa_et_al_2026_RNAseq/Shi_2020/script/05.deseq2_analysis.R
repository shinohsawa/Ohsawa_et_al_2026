library(DESeq2)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)
library(purrr)
library(scales)
library(readr)
library(tidyr)
library(openxlsx)
library(rstatix)
library(ggpubr)

# ==============================================================================
# 0. Directory Setup (Optimized for HPC environments and specified output paths)
# ==============================================================================
# Automatically set paths assuming the script is run from the script/ directory
# setwd("/Volumes/biol_bc_neurohr_scratch_1/sosawa/NGS/Ohsawa_et_al_2026_RNAseq/Shi_2020/script/") # Commented out for reproducibility

# Define output folders (based on DESeq2_results/ in the same directory level as script/)
base_out_dir <- "../DESeq2_results"
comp_out_dir <- file.path(base_out_dir, "comparison")
plot_out_dir <- file.path(base_out_dir, "plot")

# Automatically create folders if they don't exist (Error prevention)
dir.create(base_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(comp_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_out_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. Automatic Merging of Count Table (.cntTable) Files
# ==============================================================================
cat("▶️ 1/6: Merging raw count tables from te_counts_output...\n")
cnt_files <- list.files("../te_counts", pattern = "\\.cntTable$", full.names = TRUE)
if (length(cnt_files) == 0) {
  stop("Error: No .cntTable files found in ../te_counts_output/")
}

merged_counts <- cnt_files %>%
  map(function(f) {
    sample_name <- tools::file_path_sans_ext(basename(f))
    
    # 🌟 Fix Point: Add quote = "" to prevent automatic removal of double quotes!
    read_tsv(f, col_names = c("raw_ID", "count"), skip = 1, show_col_types = FALSE, quote = "") %>%
      mutate(
        Category = ifelse(str_detect(raw_ID, '"'), "gene", "TE"),
        ID = str_replace_all(raw_ID, '"', ''),
        Sample = sample_name
      ) %>%
      dplyr::select(ID, Category, Sample, count)
  }) %>%
  bind_rows() %>%
  pivot_wider(names_from = Sample, values_from = count, values_fill = 0)

# ==============================================================================
# 2. Sample Selection (small vs large) and Metadata Construction
# ==============================================================================
cat("▶️ 2/6: Filtering samples for small (S4, S8) vs large (LC1, LC2)...\n")
counts_selected <- merged_counts %>%
  dplyr::select(ID, Category, starts_with(c("S4", "S8", "LC1", "LC2")))

sample_names <- colnames(counts_selected)[3:ncol(counts_selected)]
condition <- ifelse(str_detect(sample_names, "^S"), "small", "large")
replicate <- str_extract(sample_names, "[0-9]+$")

sample_info <- data.frame(
  row.names = sample_names,
  condition = factor(condition, levels = c("small", "large")),
  replicate = replicate
)

# Low expression filter (Count > 1 in at least 3 samples)
counts_selected_filtered <- counts_selected[rowSums(counts_selected[, 3:ncol(counts_selected)] > 1) >= 3, ]

# ==============================================================================
# 3. Running DESeq2 Analysis
# ==============================================================================
cat("▶️ 3/6: Running DESeq2 differential expression analysis...\n")
counts_matrix <- counts_selected_filtered %>%
  column_to_rownames("ID") %>%
  dplyr::select(-Category) %>%
  as.matrix()

dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = sample_info,
  design = ~ condition
)
dds <- DESeq(dds)

category_info <- counts_selected %>%
  dplyr::select(ID, Category) %>%
  distinct()

res_df <- results(dds, contrast = c("condition", "large", "small")) %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  left_join(category_info, by = "ID") %>%
  mutate(
    significant = case_when(
      log2FoldChange >= 1 & padj <= 0.05 ~ "up",
      log2FoldChange <= -1 & padj <= 0.05 ~ "down",
      TRUE ~ "ns"
    )
  ) %>%
  filter(!is.na(padj)) %>%
  separate(col = ID, into = c("TE_name", "Superfamily", "Class"), sep = ":", remove = FALSE, fill = "right")

# 📊 OUTPUT 1: PCA.pdf
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = group)) +
  geom_point(size = 5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(values = c("red", "lightblue")) +
  theme_test()
ggsave(file.path(plot_out_dir, "PCA.pdf"), width = 6, height = 6)

# ==============================================================================
# 4. TPM Calculation & Saving to RData
# ==============================================================================
cat("▶️ 4/6: Calculating TPM values and saving RData...\n")
ensembl <- biomaRt::useMart(
  biomart = "ensembl", 
  dataset = "dmelanogaster_gene_ensembl", 
  host = "https://jul2022.archive.ensembl.org"
)

target_genes <- rownames(counts_matrix)
geneName <- biomaRt::getBM(
  attributes = c("ensembl_gene_id", "start_position", "end_position"), 
  filters = 'ensembl_gene_id', values = target_genes, mart = ensembl
) %>% filter(!duplicated(ensembl_gene_id)) %>% column_to_rownames("ensembl_gene_id")

gene_lengths_kb <- abs(geneName$start_position - geneName$end_position)[match(target_genes, rownames(geneName))] / 1000

rpk <- counts_matrix / gene_lengths_kb
rpk_clean <- na.omit(as.data.frame(rpk))
scale_factors <- colSums(rpk_clean) / 1e6
tpm_matrix <- t(t(rpk_clean) / scale_factors)

# 🌟 As specified, output to the same level as the comparison folder (directly under DESeq2_results/)!
save(tpm_matrix, file = file.path(base_out_dir, "tpm_matrix.RData"))

df_tpm_calc <- as.data.frame(tpm_matrix) %>%
  rownames_to_column("ID") %>%
  mutate(mean_small = rowMeans(dplyr::select(., starts_with("S"))))

df_merged <- res_df %>%
  left_join(dplyr::select(df_tpm_calc, ID, mean_small), by = "ID")

# 📊 OUTPUT 2: MA_plot_TPM_significant.pdf
label_data_TPM <- df_merged %>%
  dplyr::count(significant) %>%
  filter(significant %in% c("up", "down")) %>%
  mutate(
    label = paste0(significant, ": ", n),
    x = max(log10(df_merged$mean_small + 1), na.rm = TRUE) - 1,
    y = ifelse(significant == "up", max(df_merged$log2FoldChange, na.rm = TRUE) - 1, min(df_merged$log2FoldChange, na.rm = TRUE) + 1),
    color = ifelse(significant == "up", "red", "blue")
  )

bg_data <- data.frame(
  xmin = c(0.5, 1.5, 2.5), xmax = c(1.5, 2.5, 3.5),
  fill_color = c("lightblue", "lightgray", "lightpink")
)

ggplot(df_merged, aes(x = log10(mean_small + 1), y = log2FoldChange)) + 
  geom_rect(data = bg_data, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill_color), alpha = 0.3, inherit.aes = FALSE) +
  scale_fill_identity() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  geom_point(alpha = 0.5, size = 1, aes(color = significant)) +
  geom_smooth(method = lm , color = "lightblue")+
  scale_color_manual(values = c("blue", "gray", "red")) +
  geom_text(data = label_data_TPM, aes(x = x, y = y, label = label), color = label_data_TPM$color, inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.5) +
  theme_classic() +
  labs(x = "log10(meanTPM_small + 1)", y = "log2FC[large/small]")
ggsave(file.path(plot_out_dir, "MA_plot_TPM_significant.pdf"), width = 6, height = 4)

# 📊 OUTPUT 3: Boxplot_TPMbinned.pdf
df_binned <- df_merged %>%
  mutate(log10_TPM = log10(mean_small + 1)) %>%
  filter(log10_TPM > 0.5 & log10_TPM < 4.0) %>%
  mutate(TPM_bin = cut(log10_TPM, breaks = seq(0.5, 3.5, by = 1), include.lowest = TRUE)) %>%
  filter(!is.na(TPM_bin))

stat.test_one <- df_binned %>%
  group_by(TPM_bin) %>%
  t_test(log2FoldChange ~ 1, mu = 0, detailed = TRUE) %>%
  adjust_pvalue(method = "none") %>%
  add_significance() %>%
  add_xy_position(fun = "max", x = "TPM_bin", dodge = 0.75) %>%
  mutate(label = paste0("p = ", signif(p.adj, 2)))

ggplot(df_binned, aes(x = TPM_bin, y = log2FoldChange)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  annotate("rect", xmin = 0.55, xmax = 1.45, ymin = -Inf, ymax = Inf, fill = "lightblue", alpha = 0.3) +
  annotate("rect", xmin = 1.55, xmax = 2.45, ymin = -Inf, ymax = Inf, fill = "lightgray", alpha = 0.3) +
  annotate("rect", xmin = 2.55, xmax = 3.45, ymin = -Inf, ymax = Inf, fill = "lightpink", alpha = 0.3) +
  geom_boxplot(outlier.size = 1, outlier.color = "gray", outlier.alpha = 0.7, alpha = 0.7, width = 0.5, aes(fill = Category )) +
  theme_classic() +
  scale_fill_manual(values = c("#fd7eb7"))+
  stat_pvalue_manual(stat.test_one, label = "label", tip.length = NA, remove.bracket = TRUE, hide.ns = TRUE) +
  labs(x = "TPM bin (small)", y = "log2FC[large/small]")
ggsave(file.path(plot_out_dir, "Boxplot_TPMbinned.pdf"), width = 4, height = 4)

# ==============================================================================
# 5. Sub-telomeric Region Analysis & Heatmap
# ==============================================================================
cat("▶️ 5/6: Annotating sub-telomeric regions and plotting heatmaps...\n")
unique_genes <- res_df %>% filter(Category == "gene") %>% pull(ID) %>% unique()

geneName_coords <- biomaRt::getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position"), 
  filters = 'ensembl_gene_id', values = unique_genes, mart = ensembl
) %>% filter(!duplicated(ensembl_gene_id))
colnames(geneName_coords) <- c("ID", "gene_name", "chr_name", "start_pos", "end_pos")

target_chrs <- c("2L", "2R", "3L", "3R", "4", "X", "Y")
dm6_info <- data.frame(
  chr_name = target_chrs,
  size     = c(23513712, 25286936, 28110227, 32079331, 1348131, 23542271, 3667352),
  tel_side = c("start", "end", "start", "end", "start", "start", "start")
)

df_coords <- res_df %>%
  left_join(geneName_coords, by = "ID") %>%
  filter(chr_name %in% target_chrs) %>%
  mutate(chr_name = factor(chr_name, levels = target_chrs), region_type = "body")

subtel_window <- 80000
for (i in 1:nrow(dm6_info)) {
  chr <- dm6_info$chr_name[i]
  chr_len <- dm6_info$size[i]
  side <- dm6_info$tel_side[i]
  idx_chr <- which(df_coords$chr_name == chr & !is.na(df_coords$start_pos))
  if (length(idx_chr) > 0) {
    idx_start <- idx_chr[df_coords$start_pos[idx_chr] <= subtel_window]
    df_coords$region_type[idx_start] <- ifelse(side == "start", "subtel", "centr")
    idx_end <- idx_chr[df_coords$start_pos[idx_chr] >= (chr_len - subtel_window)]
    df_coords$region_type[idx_end] <- ifelse(side == "start", "centr", "subtel")
  }
}

# 📊 OUTPUT 4: all_contrast_subtel_boxplot.pdf
df_all_subtel <- df_coords %>% filter(region_type == "subtel")
df_all_allsubtel.stat.test_one <- df_all_subtel %>%
  t_test(log2FoldChange ~ 1, mu = 0, detailed = TRUE) %>%
  add_xy_position(fun = "max", x = "Category") %>%
  mutate(label = paste0("p = ", signif(p, 2)))

p_label <- df_all_allsubtel.stat.test_one$label[1]

# 2. Calculate a position that is "2" higher than the maximum value of the boxplot (top of the whisker)
y_pos <- max(df_all_subtel$log2FoldChange, na.rm = TRUE) + 2

# 3. Draw the graph
ggplot(df_all_subtel, aes(x = "large_vs_small", y = log2FoldChange)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  geom_boxplot(alpha = 0.7, fill = "#fd7eb7", outlier.color = "gray") +
  # 🌟 Instead of using stat_pvalue_manual, use annotate to force drawing at "X=1 (center)"!
  annotate("text", x = 1, y = y_pos, label = p_label, size = 5) +
  
  # scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + # Ensure 15% margin at the top
  theme_classic() +
  labs(x = "Contrast", y = "log2FC[large/small] in Subtelomere")
ggsave(file.path(plot_out_dir, "all_contrast_subtel_boxplot.pdf"), height = 4, width = 1.6)

# 📊 OUTPUT 5: subtel_chr_heatmap.pdf
df_subtel_chr_stat <- df_all_subtel %>%
  group_by(chr_name) %>% filter(n() > 1) %>% t_test(log2FoldChange ~ 1, mu = 0) %>%
  adjust_pvalue(method = "none") %>% add_significance()

df_subtel_all_stat <- df_all_subtel %>%
  t_test(log2FoldChange ~ 1, mu = 0) %>%
  mutate(chr_name = "all", p.adj.signif = case_when(p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns"))

df_heatmap_data <- bind_rows(df_subtel_chr_stat, df_subtel_all_stat) %>%
  left_join(bind_rows(df_all_subtel %>% group_by(chr_name) %>% summarise(mean_val = mean(log2FoldChange, na.rm = TRUE)),
                      data.frame(chr_name = "all", mean_val = mean(df_all_subtel$log2FoldChange, na.rm = TRUE))), by = "chr_name") %>%
  mutate(chr_name = factor(chr_name, levels = c("all", target_chrs)),
         tile_label = case_when(chr_name %in% c("X", "Y") ~ "no gene", p.adj.signif == "ns" ~ "", TRUE ~ p.adj.signif)) %>%
  complete(chr_name)

ggplot(df_heatmap_data, aes(x = "large_vs_small", y = chr_name, fill = mean_val)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, na.value = "grey90", oob = scales::squish, name = "Mean log2FC") +
  scale_x_discrete(position = "top") + scale_y_discrete(limits = rev) + 
  geom_text(aes(label = ifelse(tile_label == "no gene", "", tile_label)), color = "black", size = 6) +
  geom_text(aes(label = ifelse(tile_label == "no gene", "no gene", "")), color = "black", size = 3.5) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), panel.grid = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file.path(plot_out_dir, "subtel_chr_heatmap.pdf"), device = "pdf", height = 4, width = 2.5)

# ==============================================================================
# 6. Transposable Elements (TE) Classification Analysis
# ==============================================================================
cat("▶️ 6/6: Analyzing Transposable Elements (TE) and generating volcano/boxplots...\n")
df_long_TE_subset <- res_df %>% filter(Category == "TE")
sorted_levels <- names(sort(table(df_long_TE_subset$Class), decreasing = TRUE))
df_long_TE_subset$Class <- factor(df_long_TE_subset$Class, levels = sorted_levels)
df_long_TE_binned <- df_long_TE_subset %>% group_by(Class) %>% filter(n() >= 2) %>% ungroup()

# 📊 OUTPUT 6: Boxplot_TE_class.pdf
ggplot(df_long_TE_binned, aes(x = Class, y = log2FoldChange)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  geom_boxplot(aes(fill = Class), outlier.colour = "gray") +
  scale_fill_brewer(palette = "Pastel1") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "TE Class", y = "log2FC[large/small]")
ggsave(file.path(plot_out_dir, "Boxplot_TE_class.pdf"), width = 4, height = 4)

# 📊 OUTPUT 7: volcanoplot_TE_class_color.pdf
df_long_TE_ns <- df_long_TE_subset %>% filter(significant == "ns")
df_long_TE_sign <- df_long_TE_subset %>% filter(significant != "ns")

ggplot() + 
  geom_hline(yintercept = -log10(0.05), color = "gray", linetype = "dashed") +
  geom_vline(xintercept = -1, color = "gray", linetype = "dashed") +
  geom_vline(xintercept = 1, color = "gray", linetype = "dashed") +
  geom_point(data = df_long_TE_ns, aes(y = -log10(padj), x = log2FoldChange), alpha = 0.5, size = 1, color = "gray") +
  geom_point(data = df_long_TE_sign, aes(y = -log10(padj), x = log2FoldChange, color = Class), size = 2) +
  scale_color_brewer(palette = "Pastel1") + theme_classic() +
  labs(x = "log2FoldChange", y = "-log10(padj)")
ggsave(file.path(plot_out_dir, "volcanoplot_TE_class_color.pdf"), width = 5, height = 4)

# ==============================================================================
# Save Merged Excel Data
# ==============================================================================
write.xlsx(res_df, file = file.path(comp_out_dir, "All_Contrasts_Merged.xlsx"))
cat("\n🎉 [Success] All 7 cleaned figures, tables, and RData saved to DESeq2_results/ successfully!\n")