library(DESeq2)
library(openxlsx)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)
library(purrr)
library(ggrepel)
library(ggpubr)
library(rstatix)
library(readr)
library(tidyr)
library(rtracklayer)
library(GenomicRanges)

# ==============================================================================
# 0. Directory Setup (Folder structure perfectly matched with the Drosophila analysis)
# ==============================================================================
# Set the script directory as the working directory (change if necessary)
# setwd("/Volumes/biol_bc_neurohr_scratch_1/sosawa/NGS/Ohsawa_et_al_2026_RNAseq/human_RPE1/script/")

base_out_dir <- "../DESeq2_results"
comp_out_dir <- file.path(base_out_dir, "comparison")
plot_out_dir <- file.path(base_out_dir, "plot")

# Automatically create folders
dir.create(base_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(comp_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_out_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. Python script alternative: Automatic merging and classification of count tables (.cntTable)
# ==============================================================================
cat("▶️ 1/6: Merging raw count tables and identifying Gene/TE...\n")
cnt_files <- list.files("../te_counts", pattern = "\\.cntTable$", full.names = TRUE)

merged_counts <- cnt_files %>%
  map(function(f) {
    sample_name <- tools::file_path_sans_ext(basename(f))
    # Remove double quotes like in Python, classify as "gene" if starting with ENSG, otherwise "TE"
    read_tsv(f, col_names = c("raw_ID", "count"), skip = 1, show_col_types = FALSE, quote = "") %>%
      mutate(
        ID = str_replace_all(raw_ID, '"', ''),
        category = ifelse(str_starts(ID, "ENSG"), "gene", "TE"),
        Sample = sample_name
      ) %>%
      dplyr::select(ID, category, Sample, count)
  }) %>%
  bind_rows() %>%
  pivot_wider(names_from = Sample, values_from = count, values_fill = 0)

# ==============================================================================
# 2. Data filtering and metadata construction
# ==============================================================================
cat("▶️ 2/6: Filtering data (keeping counts >= 1 in at least 4 samples)...\n")

# Extract count columns (from 3rd column onwards) and filter
counts_only <- merged_counts %>% dplyr::select(-ID, -category)
keep <- rowSums(counts_only >= 1) >= 4
filtered_data <- merged_counts[keep, ]

count_matrix <- filtered_data %>% dplyr::select(-ID, -category) %>% as.matrix()
rownames(count_matrix) <- filtered_data$ID

# Metadata (day2 vs day7)
sample_names <- colnames(count_matrix)
condition <- ifelse(grepl("day2", sample_names), "day2", "day7")
colData <- data.frame(row.names = sample_names, condition = factor(condition, levels = c("day2", "day7")))

# ==============================================================================
# 3. DESeq2 Analysis & 📊 OUTPUT 1: PCA Plot
# ==============================================================================
cat("▶️ 3/6: Running DESeq2 and generating PCA...\n")
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = colData, design = ~ condition)
dds <- DESeq(dds)

# PCA Plot
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(show.legend = FALSE) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_classic() +
  scale_color_manual(values = c("day2" = "gray", "day7" = "coral2"))
ggsave(file.path(plot_out_dir, "PCA_plot_standard_norm.pdf"), plot = p_pca, width = 5, height = 4)

# ==============================================================================
# 4. Results extraction, TPM calculation & RData output
# ==============================================================================
cat("▶️ 4/6: Extracting results and calculating TPM...\n")
res <- results(dds, contrast = c("condition", "day7", "day2"))

df_long <- as.data.frame(res) %>%
  rownames_to_column("ID") %>%
  left_join(dplyr::select(filtered_data, ID, category), by = "ID") %>%
  filter(!is.na(padj)) %>%
  separate(col = ID, into = c("TE_name", "Superfamily", "Class"), sep = ":", remove = FALSE, fill = "right") %>%
  mutate(significant = case_when(
    log2FoldChange >= 1 & padj <= 0.05 ~ "up",
    log2FoldChange <= -1 & padj <= 0.05 ~ "down",
    TRUE ~ "ns"
  ))

# Load GTF files (adjust paths to your actual environment)
gtf_human <- import("../reference/gencode.v44.primary_assembly.annotation.gtf")
gtf_te    <- import("../reference/GRCh38_GENCODE_rmsk_TE.gtf")

# Calculate exon lengths
exons_human <- gtf_human[gtf_human$type == "exon"]
reduced_human <- IRanges::reduce(splitAsList(exons_human, exons_human$gene_id))
human_lengths <- sum(width(reduced_human))

exons_te <- if("exon" %in% gtf_te$type) gtf_te[gtf_te$type == "exon"] else gtf_te
reduced_te <- IRanges::reduce(splitAsList(exons_te, exons_te$gene_id))
te_lengths <- sum(width(reduced_te))

all_lengths_bp <- c(human_lengths, te_lengths)
gene_lengths_kb <- all_lengths_bp[rownames(count_matrix)] / 1000

# Calculate TPM matrix
rpk <- count_matrix / gene_lengths_kb
rpk[is.na(rpk)] <- 0
scale_factors <- colSums(rpk) / 1e6
tpm_matrix <- t(t(rpk) / scale_factors)

# Save TPM matrix
save(tpm_matrix, file = file.path(base_out_dir, "tpm_matrix.RData"))

# Create df_merged
df_tpm_calc <- as.data.frame(tpm_matrix) %>%
  rownames_to_column("ID") %>%
  mutate(meanTPM_day2 = rowMeans(dplyr::select(., contains("day2"))))

df_merged <- left_join(df_long, dplyr::select(df_tpm_calc, ID, meanTPM_day2), by = "ID")

# 📊 OUTPUT 2: MA Plot (TPM)
label_data_TPM <- df_merged %>% dplyr::count(significant) %>% filter(significant %in% c("up", "down")) %>%
  mutate(label = paste0(significant, ": ", n),
         x = max(log10(df_merged$meanTPM_day2 + 1), na.rm = TRUE) - 2,
         y = ifelse(significant == "up", max(df_merged$log2FoldChange, na.rm = TRUE) - 1, min(df_merged$log2FoldChange, na.rm = TRUE) + 1),
         color = ifelse(significant == "up", "red", "blue"))

# 🌟 Fix Point 1: Add the 0~0.5 interval to bg_data and create a new "alpha_val" column
bg_data <- data.frame(
  xmin       = c(0, 0.5, 1.5, 2.5), 
  xmax       = c(0.5, 1.5, 2.5, 3.5), 
  fill_color = c("lightblue", "lightblue", "lightgray", "lightpink"),
  alpha_val  = c(0.7, 0.3, 0.3, 0.3)  # Specify 0.7 for 0~0.5, and 0.3 for the others!
)

p_ma <- ggplot(df_merged, aes(x = log10(meanTPM_day2 + 1), y = log2FoldChange)) + 
  
  # 🌟 Fix Point 2: Put alpha = alpha_val inside aes(), and remove the fixed alpha = 0.3
  geom_rect(data = bg_data, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill_color, alpha = alpha_val), inherit.aes = FALSE) +
  
  # 🌟 Fix Point 3: Add scale_alpha_identity() to let the numeric values be recognized directly as transparency
  scale_fill_identity() + 
  scale_alpha_identity() + 
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  geom_point(alpha = 0.5, size = 1, aes(color = significant)) + 
  scale_color_manual(values = c("blue", "gray", "red")) +
  geom_smooth(method = "lm",  color = "lightblue") +
  geom_text(data = label_data_TPM, aes(x = x, y = y, label = label), color = label_data_TPM$color, inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.5) +
  theme_classic()

p_ma 
ggsave(file.path(plot_out_dir, "MA_plot_TPM_significant.pdf"), plot = p_ma, width = 6, height = 4)
# 📊 OUTPUT 3: Boxplot TPM Binned
df_binned <- df_merged %>% mutate(log10_TPM = log10(meanTPM_day2 + 1)) %>% filter(log10_TPM >= 0 & log10_TPM <= 3.5) %>%
  mutate(TPM_bin = cut(log10_TPM, breaks = c(0, 0.5, 1.5, 2.5, 3.5), include.lowest = TRUE)) %>% filter(!is.na(TPM_bin))

stat.test_one <- df_binned %>% group_by(TPM_bin) %>% t_test(log2FoldChange ~ 1, mu = 0, detailed = TRUE) %>%
  adjust_pvalue(method = "none") %>% add_significance() %>% add_xy_position(fun = "max", x = "TPM_bin", dodge = 0.75) %>%
  mutate(label = paste0("p = ", signif(p.adj, 2)))

p_box <- ggplot(df_binned, aes(x = TPM_bin, y = log2FoldChange)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  annotate("rect", xmin = 0.55, xmax = 1.45, ymin = -Inf, ymax = Inf, fill = "lightblue", alpha = 0.7) +
  annotate("rect", xmin = 1.55, xmax = 2.45, ymin = -Inf, ymax = Inf, fill = "lightblue", alpha = 0.3) +
  annotate("rect", xmin = 2.55, xmax = 3.45, ymin = -Inf, ymax = Inf, fill = "lightgray", alpha = 0.3) +
  annotate("rect", xmin = 3.55, xmax = 4.45, ymin = -Inf, ymax = Inf, fill = "lightpink", alpha = 0.3) +
  geom_boxplot(fill = "#9dc3e6", outlier.size = 1, outlier.color = "gray", outlier.alpha = 0.7, alpha = 0.7, width = 0.7) +
  theme_classic() + stat_pvalue_manual(stat.test_one, label = "label", tip.length = NA, remove.bracket = TRUE, hide.ns = TRUE)
p_box
ggsave(file.path(plot_out_dir, "Boxplot_TPMbinned.pdf"), plot = p_box, width = 4, height = 4)

# ==============================================================================
# 5. Sub-telomere Analysis & Heatmap
# ==============================================================================
cat("▶️ 5/6: Annotating sub-telomeric regions...\n")
hg38_lengths <- read.table("../reference/STAR_index/chrNameLength.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(hg38_lengths) <- c("seqnames", "chr_length")
main_chrs <- paste0("chr", c(1:22, "X", "Y"))
hg38_lengths <- hg38_lengths %>% filter(seqnames %in% main_chrs)

df_genes <- as.data.frame(gtf_human[gtf_human$type == "gene"]) %>% mutate(base_id = sub("\\..*", "", gene_id))
df_genes_subtel <- df_genes %>% inner_join(hg38_lengths, by = "seqnames") %>%
  mutate(is_subtel = (start <= 500000) | (end >= (chr_length - 500000))) %>%
  dplyr::select(base_id, is_subtel) %>% distinct(base_id, .keep_all = TRUE)

centromeres_hg38 <- data.frame(seqnames = main_chrs,
                               centromere_pos = c(122026459, 92326171, 90504854, 49660117, 46405641, 58830166, 58054331, 43838887, 47367679, 39258202, 51644205, 34856694, 16000000, 16000000, 17000000, 35335801, 22263006, 15460898, 24681782, 26369569, 11288129, 13728592, 58100000, 10300000))

df_chr_arm <- df_genes %>% dplyr::select(base_id, seqnames, start) %>% distinct(base_id, .keep_all = TRUE) %>%
  inner_join(centromeres_hg38, by = "seqnames") %>% mutate(chr_no = as.character(seqnames), arm = ifelse(start < centromere_pos, "left", "right")) %>% dplyr::select(base_id, chr_no, arm)

df_merged <- df_merged %>% mutate(base_id = sub("\\..*", "", ID)) %>% left_join(df_genes_subtel, by = "base_id") %>%
  mutate(subtel = ifelse(is_subtel == TRUE, "subtel", "no")) %>% left_join(df_chr_arm, by = "base_id") %>%
  mutate(chr_no = ifelse(category == "TE" | is.na(chr_no), "TE", chr_no), arm = ifelse(category == "TE" | is.na(arm), "TE", arm)) %>% dplyr::select(-base_id, -is_subtel)

# 📊 OUTPUT 4: subtel_boxplot.pdf
df_subtel <- df_merged %>% filter(subtel == "subtel") %>% mutate(Contrast = "Day7_vs_Day2")
df_subtel.stat.test_one <- df_subtel %>% t_test(log2FoldChange ~ 1, mu = 0, detailed = TRUE) %>%
  add_xy_position(fun = "max", x = "Contrast") %>% mutate(label = paste0("p = ", signif(p, 2)))

y_pos <- max(df_subtel$log2FoldChange, na.rm = TRUE) + 0.5
p_subtel <- ggplot(df_subtel, aes(x = Contrast, y = log2FoldChange)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  geom_boxplot(alpha = 0.7, fill = "#9dc3e6", outlier.color = "gray") +
  annotate("text", x = 1, y = y_pos, label = df_subtel.stat.test_one$label[1], size = 5) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + theme_classic()
p_subtel
ggsave(file.path(plot_out_dir, "subtel_boxplot.pdf"), plot = p_subtel, height = 4, width = 1.6)

# ==============================================================================
# 📊 OUTPUT 5: subtel_chr_heatmap.pdf
# ==============================================================================

# ① Test per chromosome (n > 1)
df_subtel_chr_stat <- df_subtel %>% 
  group_by(Contrast, chr_name = chr_no) %>% 
  filter(n() > 1) %>% 
  t_test(log2FoldChange ~ 1, mu = 0, detailed = TRUE) %>% 
  adjust_pvalue(method = "none") %>% 
  add_significance()

# ② Test for all chromosomes combined (all)
df_subtel_all_stat <- df_subtel %>% 
  group_by(Contrast) %>% 
  t_test(log2FoldChange ~ 1, mu = 0, detailed = TRUE) %>% 
  adjust_pvalue(method = "none") %>% 
  add_significance() %>% 
  mutate(chr_name = "all")

# ③ Combine test results
df_combined_stat <- bind_rows(df_subtel_chr_stat, df_subtel_all_stat)

# ④ Mean values per chromosome
df_means_chr <- df_subtel %>% 
  group_by(Contrast, chr_name = chr_no) %>% 
  summarise(mean_log2FC = mean(log2FoldChange, na.rm = TRUE), .groups = "drop")

# ⑤ Overall mean values
df_means_all <- df_subtel %>% 
  group_by(Contrast) %>% 
  summarise(mean_log2FC = mean(log2FoldChange, na.rm = TRUE), .groups = "drop") %>% 
  mutate(chr_name = "all")

# ⑥ Combine mean values
df_means_combined <- bind_rows(df_means_chr, df_means_all)

# ⑦ Data integration and formatting (★ Apply the perfect logic from the previous analysis)
chr_levels <- c("all", paste0("chr", 1:22), "chrX", "chrY")

df_heatmap_data <- df_combined_stat %>%
  left_join(df_means_combined, by = c("Contrast", "chr_name")) %>%
  # Sort and complete based on factor levels
  mutate(chr_name = factor(chr_name, levels = chr_levels)) %>%
  complete(chr_name, Contrast) %>%
  mutate(
    # Convert P-values to asterisks (convert NA to ns)
    p.adj.signif = case_when(
      is.na(p) ~ "ns",
      p < 0.0001 ~ "****",
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    # Determine tile labels
    tile_label = case_when(
      is.na(mean_log2FC) ~ "no gene",
      TRUE ~ ""
    ),
    # Determine 3-row layout (works perfectly since it is calculated after chr14 is generated by complete)
    plot_row = factor(case_when(
      chr_name %in% c("all", paste0("chr", 1:8)) ~ "Row 1",
      chr_name %in% paste0("chr", 9:16)          ~ "Row 2",
      TRUE                                       ~ "Row 3"
    ), levels = c("Row 1", "Row 2", "Row 3"))
  )

# ⑧ Draw the 3-row layout heatmap
p_heatmap <- ggplot(df_heatmap_data, aes(x = chr_name, y = Contrast, fill = mean_log2FC)) + 
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, na.value = "grey90", oob = scales::squish, name = "Mean log2FC\n(Subtel)") +
  scale_x_discrete(position = "top") + 
  scale_y_discrete(limits = rev) +
  geom_text(aes(label = ifelse(p.adj.signif == "ns", "", p.adj.signif)), color = "black", size = 5, vjust = 0.7) +
  geom_text(aes(label = tile_label), color = "black", size = 3) + 
  facet_wrap(~ plot_row, ncol = 1, scales = "free_x") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(hjust = 0.5), 
    panel.grid = element_blank(), 
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(), 
    strip.background = element_blank(), 
    strip.text = element_blank()
  )

print(p_heatmap)
ggsave(file.path(plot_out_dir, "subtel_chr_heatmap.pdf"), plot = p_heatmap, device = "pdf", width = 7.5, height = 3)
# ==============================================================================
# 📊 OUTPUT 5: subtel_chr_heatmap.pdf (Subtel vs Non-Subtel Two-sample t-test)
# ==============================================================================

# ① Test per chromosome (Subtel vs Non-Subtel)
df_chr_stat <- df_merged %>%
  # Check if there are at least 2 genes for both "subtel" and "no" per chromosome
  group_by(Contrast, chr_name = chr_no, subtel) %>%
  mutate(n_genes = n()) %>%
  ungroup() %>%
  group_by(Contrast, chr_name) %>%
  filter(min(n_genes) > 1 & length(unique(subtel)) == 2) %>% 
  t_test(log2FoldChange ~ subtel, detailed = TRUE) %>% 
  adjust_pvalue(method = "none") %>% 
  add_significance()

# ② Test for all chromosomes combined (all) (Subtel vs Non-Subtel)
df_all_stat <- df_merged %>%
  mutate(chr_name = "all") %>%
  group_by(Contrast, chr_name) %>%
  t_test(log2FoldChange ~ subtel, detailed = TRUE) %>%
  adjust_pvalue(method = "none") %>%
  add_significance()

# ③ Combine test results
df_combined_stat <- bind_rows(df_chr_stat, df_all_stat)

# ④ Mean values per chromosome (use the mean log2FC of sub-telomeric regions for heatmap coloring)
df_means_chr <- df_merged %>%
  filter(subtel == "subtel") %>%
  group_by(Contrast, chr_name = chr_no) %>%
  summarise(mean_log2FC = mean(log2FoldChange, na.rm = TRUE), .groups = "drop")

# ⑤ Overall mean values
df_means_all <- df_merged %>%
  filter(subtel == "subtel") %>%
  group_by(Contrast) %>%
  summarise(mean_log2FC = mean(log2FoldChange, na.rm = TRUE), .groups = "drop") %>%
  mutate(chr_name = "all")

# ⑥ Combine mean values
df_means_combined <- bind_rows(df_means_chr, df_means_all)

# ⑦ Data integration and formatting
chr_levels <- c("all", paste0("chr", 1:22), "chrX", "chrY")

df_heatmap_data <- df_combined_stat %>%
  left_join(df_means_combined, by = c("Contrast", "chr_name")) %>%
  mutate(chr_name = factor(chr_name, levels = chr_levels)) %>%
  complete(chr_name, Contrast) %>%
  mutate(
    # Convert P-values to asterisks (convert NA to ns)
    p.adj.signif = case_when(
      is.na(p) ~ "ns",
      p < 0.0001 ~ "****",
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    # Determine tile labels
    tile_label = case_when(
      is.na(mean_log2FC) ~ "no gene",
      TRUE ~ ""
    ),
    # Determine 3-row layout
    plot_row = factor(case_when(
      chr_name %in% c("all", paste0("chr", 1:8)) ~ "Row 1",
      chr_name %in% paste0("chr", 9:16)          ~ "Row 2",
      TRUE                                       ~ "Row 3"
    ), levels = c("Row 1", "Row 2", "Row 3"))
  )

# ⑧ Draw the 3-row layout heatmap
p_heatmap <- ggplot(df_heatmap_data, aes(x = chr_name, y = Contrast, fill = mean_log2FC)) + 
  geom_tile(color = "white") +
  # Change legend title to indicate it's a "Subtel vs Non-Subtel" comparison
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, na.value = "grey90", oob = scales::squish, name = "Mean log2FC\n(Subtel)") +
  scale_x_discrete(position = "top") + 
  scale_y_discrete(limits = rev) +
  geom_text(aes(label = ifelse(p.adj.signif == "ns", "", p.adj.signif)), color = "black", size = 5, vjust = 0.7) +
  geom_text(aes(label = tile_label), color = "black", size = 3) + 
  facet_wrap(~ plot_row, ncol = 1, scales = "free_x") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(hjust = 0.5), 
    panel.grid = element_blank(), 
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(), 
    strip.background = element_blank(), 
    strip.text = element_blank()
  )

print(p_heatmap)
ggsave(file.path(plot_out_dir, "subtel_chr_heatmap_twosample.pdf"), plot = p_heatmap, device = "pdf", width = 7.5, height = 3)
# ==============================================================================
# 6. Transposable Elements (TE) Classification Analysis
# ==============================================================================
cat("▶️ 6/6: Analyzing Transposable Elements (TE) and generating volcano/boxplots...\n")
df_long_TE <- df_merged %>% filter(category == "TE", Class %in% c("LTR", "DNA", "LINE", "SINE", "Unknown", "Satellite")) %>%
  mutate(Class = factor(Class, levels = c("LTR", "LINE", "DNA", "Satellite", "Unknown", "SINE")))

# 📊 OUTPUT 6: volcanoplot_TE_class_color.pdf
p_volcano_te <- ggplot() + 
  geom_hline(yintercept = -log10(0.05), color = "black", linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), color = "black", linetype = "dashed") +
  geom_point(data = filter(df_long_TE, significant == "ns"), aes(y = -log10(padj), x = log2FoldChange), alpha = 0.5, size = 1, color = "gray") +
  geom_point(data = filter(df_long_TE, significant != "ns"), aes(y = -log10(padj), x = log2FoldChange, color = Class), size = 2) +
  scale_color_brewer(palette = "Pastel1") + theme_classic()
ggsave(file.path(plot_out_dir, "volcanoplot_TE_class_color.pdf"), plot = p_volcano_te, width = 5, height = 4)


# 📊 OUTPUT 7: Boxplot_TE_class.pdf
p_box_te <- df_long_TE %>% group_by(Class) %>% filter(n() >= 2) %>% ungroup() %>%
  ggplot(aes(x = Class, y = log2FoldChange)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  geom_boxplot(aes(fill = Class), outlier.colour = "gray") +
  scale_fill_brewer(palette = "Pastel1") + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(plot_out_dir, "Boxplot_TE_class.pdf"), plot = p_box_te, width = 4.5, height = 4)

# ==============================================================================
# Save Merged Excel Data
# ==============================================================================
write.xlsx(df_merged, file = file.path(comp_out_dir, "All_Contrasts_Merged.xlsx"))
cat("\n🎉 [Success] All 7 human-specific figures and RData saved to DESeq2_results/ successfully!\n")