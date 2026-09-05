# ==============================================================================
# Yeast RNA-seq: Terhorst 2023 Master Pipeline
# (rRNA Filter -> DESeq2 -> TPM -> PCA -> DE Export -> Final 6 Plots)
# ==============================================================================

# Load necessary libraries once
library(DESeq2)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(rstatix)
library(biomaRt)
library(readxl)

# Change directory according to your environment
# setwd("/Volumes/biol_bc_neurohr_scratch_1/sosawa/NGS/Ohsawa_et_al_2026_RNAseq/Terhorst_2023/") 

# Create necessary directories at once
dir.create("./DESeq2_results", showWarnings = FALSE)
dir.create("./DESeq2_results/comparison", showWarnings = FALSE)
dir.create("./DESeq2_results/plots", recursive = TRUE, showWarnings = FALSE)

cat("=" %+% strrep("=", 79) %+% "\n")
cat("Yeast RNA-seq Fully Integrated Pipeline (DESeq2 + TPM + Plots)\n")
cat(strrep("=", 80) %+% "\n\n")

# ==============================================================================
# PART 1: Load and prepare data (Noise Filter)
# ==============================================================================
cat("\n[PART 1] Loading and filtering data...\n")

cdc28 <- read.csv("./gene_counts/raw_counts.csv", header = TRUE, row.names = 1, check.names = FALSE)
colnames(cdc28) <- gsub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", colnames(cdc28))

# Apply 6-layer noise filter
filter_pattern <- "^RDN|^YLR154|^YLR162|^YLR157|^YLR159|^YLR156|^ETS|^ITS|^snR|^SCR1|^LSR1"
noise_rows <- grepl(filter_pattern, rownames(cdc28))
cdc28 <- cdc28[!noise_rows, ]
cat("- Filtered clean mRNA genes:", nrow(cdc28), "genes\n")

unique_genes <- unique(row.names(cdc28))
cdc28 <- round(cdc28)

# Create metadata
parse_sample_info <- function(sample_names) {
  strain <- ifelse(grepl("WT", sample_names), "WT", "cdc28")
  treatment <- ifelse(grepl("2h", sample_names), "2h", "8h")
  replicate <- gsub(".*_(rep[1-3])$", "\\1", sample_names)
  data.frame(sample = sample_names, strain = strain, treatment = treatment, replicate = replicate, stringsAsFactors = FALSE)
}
sample_info <- parse_sample_info(colnames(cdc28))

# Low expression gene filter
cdc28_filtered <- cdc28[rowSums(cdc28 > 1) >= 3, ]
cat("- Filtered to", nrow(cdc28_filtered), "expressed genes\n")

# ==============================================================================
# PART 2: BiomaRt Fetch (Run only once!)
# ==============================================================================
cat("\n[PART 2] Fetching Gene Info via BiomaRt (Length, Chr, Position)...\n")

ensembl <- biomaRt::useMart(biomart ="fungi_mart", dataset="scerevisiae_eg_gene", host="https://fungi.ensembl.org")
geneName <- biomaRt::getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "description", "chromosome_name", "start_position", "end_position", "strand", "transcript_length"), 
  filters = 'ensembl_gene_id', values = unique_genes, mart = ensembl
)
geneName <- geneName[!duplicated(geneName$ensembl_gene_id), ]
rownames(geneName) <- geneName$ensembl_gene_id
colnames(geneName) <- c("gene", "gene_name", "description", "chr_name", "start_pos", "end_pos", "strand_orient", "transcript_length")

# ==============================================================================
# PART 3: DESeq2 Analysis and TPM Calculation
# ==============================================================================
cat("\n[PART 3] Running DESeq2 & Calculating TPM...\n")

col_data <- data.frame(
  strain = factor(sample_info$strain, levels = c("WT", "cdc28")),
  treatment = factor(sample_info$treatment, levels = c("2h", "8h")),
  row.names = sample_info$sample
)

dds <- DESeqDataSetFromMatrix(countData = cdc28_filtered, colData = col_data, design = ~ strain + treatment)
dds <- DESeq(dds)
normalized_counts <- counts(dds, normalized = TRUE)

# Calculate TPM
gene_lengths_bp <- geneName$transcript_length[match(rownames(cdc28_filtered), geneName$gene)]
valid_idx <- !is.na(gene_lengths_bp)
gene_lengths_kb <- gene_lengths_bp / 1000

rpk <- cdc28_filtered
rpk[valid_idx, ] <- cdc28_filtered[valid_idx, ] / gene_lengths_kb[valid_idx]
rpk[!valid_idx, ] <- NA 

scale_factors <- colSums(rpk, na.rm = TRUE) / 1e6
tpm_matrix <- t(t(rpk) / scale_factors)

save(tpm_matrix, file = "./DESeq2_results/tpm_matrix.RData")

# ==============================================================================
# PART 4: PCA Plot
# ==============================================================================
cat("\n[PART 4] Generating PCA Plot...\n")

vsd_subset <- vst(dds, blind = FALSE)
pcaData_subset <- plotPCA(vsd_subset, intgroup = c("strain", "treatment"), returnData = TRUE)
percentVar_subset <- round(100 * attr(pcaData_subset, "percentVar"))

p_pca <- ggplot(pcaData_subset, aes(PC1, PC2, color = strain, shape = treatment)) +
  geom_point(size = 5) +
  xlab(paste0("PC1: ", percentVar_subset[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_subset[2], "% variance")) +
  scale_color_manual(values=c("blue","orange")) +
  ggtitle("Samples: PCA Plot (rRNA Filtered Normalized)") +
  theme_test() 
ggsave("./DESeq2_results/plots/pca.pdf", plot = p_pca, width = 6, height = 6)

# ==============================================================================
# PART 5: DE Analysis (8h vs 2h) & Collect Results in Memory
# ==============================================================================
cat("\n[PART 5] Differential expression analysis (8h vs 2h)...\n")

get_group_means <- function(norm_counts, group1_samps, group2_samps) {
  list(mean1 = rowMeans(norm_counts[, group1_samps, drop = FALSE]), 
       mean2 = rowMeans(norm_counts[, group2_samps, drop = FALSE]))
}

df_list <- list() # Export to Excel while generating a list in memory for plotting

for (strain in unique(sample_info$strain)) {
  strain_cols <- sample_info$sample[sample_info$strain == strain]
  dds_strain <- dds[, strain_cols]
  dds_strain$strain <- droplevels(dds_strain$strain)
  design(dds_strain) <- ~ treatment
  dds_strain <- DESeq(dds_strain, quiet = TRUE)
  res <- results(dds_strain, contrast = c("treatment", "8h", "2h"))
  
  t8h_samps <- sample_info$sample[sample_info$strain == strain & sample_info$treatment == "8h"]
  t2h_samps <- sample_info$sample[sample_info$strain == strain & sample_info$treatment == "2h"]
  
  means <- get_group_means(counts(dds_strain, normalized = TRUE), t2h_samps, t8h_samps)
  
  # Create results dataframe
  df <- as.data.frame(res)
  df$gene <- rownames(df)
  df$mean_group1 <- means$mean1[rownames(df)]
  df$mean_group2 <- means$mean2[rownames(df)]
  df <- df[, c("gene", "baseMean", "mean_group1", "mean_group2", "log2FoldChange", "lfcSE", "pvalue", "padj")]
  df <- df[order(df$padj), ]
  
  # Export to Excel
  fname <- paste0(strain, "_8h_vs_2h")
  write.xlsx(df, file = paste0("./DESeq2_results/comparison/", fname, ".xlsx"))
  
  # Save to list
  df$Contrast <- fname
  df_list[[fname]] <- df
}

# Combine lists to create df_all
df_all <- bind_rows(df_list)
df_all$Contrast <- factor(df_all$Contrast, levels = c("WT_8h_vs_2h", "cdc28_8h_vs_2h"))

# ==============================================================================
# PART 6: Build Plotting Data (TPM append & BiomaRt merge)
# ==============================================================================
cat("\n[PART 6] Preparing Data for Final Plots...\n")

# Significance labels
df_all$significant <- ifelse(df_all$log2FoldChange >= 1 & df_all$padj <= 0.05, "up",
                             ifelse(df_all$log2FoldChange <= -1 & df_all$padj <= 0.05, "down", "ns"))
df_all$significant[is.na(df_all$significant)] <- "ns"

# Calculate average TPM_2h and append directly to df_all (Super smart implementation!)
wt_2h_samps <- sample_info$sample[sample_info$strain == "WT" & sample_info$treatment == "2h"]
cdc28_2h_samps <- sample_info$sample[sample_info$strain == "cdc28" & sample_info$treatment == "2h"]

meanTPM_WT_2h <- rowMeans(tpm_matrix[, wt_2h_samps, drop=FALSE], na.rm=TRUE)
meanTPM_cdc28_2h <- rowMeans(tpm_matrix[, cdc28_2h_samps, drop=FALSE], na.rm=TRUE)

df_all <- df_all %>%
  mutate(meanTPM_2h = case_when(
    grepl("WT", Contrast) ~ meanTPM_WT_2h[gene],
    grepl("cdc28", Contrast) ~ meanTPM_cdc28_2h[gene],
    TRUE ~ NA_real_
  ))

# Merge with BiomaRt info (fetched in PART 2)
df_all <- left_join(df_all, geneName, by="gene")

# Factorize chromosomes & calculate Subtelomeres
df_all$chr_name <- factor(df_all$chr_name, levels = c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "Mito"))

genome_size <- read.csv("../yeast/reference/sacCer3.chrom.sizes.csv")
df_all$subtel <- "no"
for (i in seq_len(nrow(genome_size))) {
  chr <- genome_size$chr_name[i]
  chr_len <- genome_size$size[i]
  subtel_idx <- which(df_all$chr_name == chr & (df_all$start_pos <= 30000 | df_all$start_pos >= (chr_len - 30000)))
  df_all$subtel[subtel_idx] <- "subtel"
}

df_all_noMito <- subset(df_all, chr_name != "Mito")

# ==============================================================================
# PART 7: Generating Final 6 Plots
# ==============================================================================
cat("\n[PART 7] Generating Final Summary Plots...\n")

# [Plot 1] MA Plot (TPM)
label_data_TPM <- df_all %>% count(Contrast, significant) %>% filter(significant %in% c("up", "down")) %>%
  mutate(label = paste0(significant, ": ", n),
         x = max(log10(df_all$meanTPM_2h + 1), na.rm = TRUE) - 2,
         y = ifelse(significant == "up", max(df_all$log2FoldChange, na.rm = TRUE) - 1, min(df_all$log2FoldChange, na.rm = TRUE) + 1),
         color = ifelse(significant == "up", "red", "blue"))

bg_data <- data.frame(xmin = c(0.5, 1.5, 2.5), xmax = c(1.5, 2.5, 3.5), fill_color = c("lightblue", "lightgray", "lightpink"))

p1 <- ggplot(df_all, aes(x = log10(meanTPM_2h + 1), y = log2FoldChange)) + 
  geom_rect(data = bg_data, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill_color), alpha = 0.3, inherit.aes = FALSE) +
  scale_fill_identity() + geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_point(alpha = 0.5, size = 1, aes(color = significant)) + scale_color_manual(values = c("blue", "gray", "red")) +
  geom_smooth(method = lm, color ="lightblue")+
  geom_text(data = label_data_TPM, aes(x = x, y = y, label = label), color = label_data_TPM$color, inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.5) +
  coord_cartesian(xlim = c(0, 4)) + theme_classic() + facet_wrap(vars(Contrast)) + labs(x = "log10(meanTPM_2h + 1)", y = "log2FC[8h/2h]")

ggsave("./DESeq2_results/plots/MA_plot_TPM_significant.pdf", plot = p1, width = 7, height = 4)


# [Plot 2] TPM Binned Boxplot
df_binned <- df_all %>% mutate(log10_TPM = log10(meanTPM_2h + 1)) %>% filter(log10_TPM > 0.5 & log10_TPM < 4.0) %>%
  mutate(TPM_bin = cut(log10_TPM, breaks = seq(0.5, 3.5, by = 1), include.lowest = TRUE)) %>% filter(!is.na(TPM_bin))

stat_test_binned <- df_binned %>% group_by(TPM_bin) %>% t_test(log2FoldChange ~ Contrast) %>%
  adjust_pvalue(method = "BH") %>% add_significance("p") %>% add_xy_position(fun = "max", x = "TPM_bin", dodge = 0.8) %>%
  mutate(label = paste0("p = ", signif(p.adj, 2)))

p2 <- ggplot(df_binned, aes(x = TPM_bin, y = log2FoldChange, fill = Contrast)) +
  annotate("rect", xmin = 0.55, xmax = 1.45, ymin = -Inf, ymax = Inf, fill = "lightblue", alpha = 0.3) +
  annotate("rect", xmin = 1.55, xmax = 2.45, ymin = -Inf, ymax = Inf, fill = "lightgray", alpha = 0.3) +
  annotate("rect", xmin = 2.55, xmax = 3.45, ymin = -Inf, ymax = Inf, fill = "lightpink", alpha = 0.3) +
  geom_boxplot(outlier.size = 1, outlier.color = "gray", outlier.alpha = 0.7, alpha = 0.7, width = 0.7, position = position_dodge(0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) + scale_fill_manual(values = c("gray", "#91cf50")) +
  theme_classic() + labs(x = "basal expression level (log10[TPM_2h + 1])", y = "log2FC[8h/2h]") +
  stat_pvalue_manual(stat_test_binned, label = "label", tip.length = 0.01, hide.ns = TRUE)

ggsave("./DESeq2_results/plots/boxplot_TPMbinned.pdf", plot = p2, width = 6, height = 4)


# [Plot 3] Subtelomere Boxplot
df_subtel_only <- subset(df_all_noMito, subtel == "subtel")

stat_test_subtel <- df_subtel_only %>% t_test(log2FoldChange ~ Contrast) %>% adjust_pvalue(method = "none") %>%
  add_significance() %>% add_xy_position(fun = "max", x = "Contrast", dodge = 0.9) %>% mutate(label = paste0("p = ", signif(p.adj, 2)))

p3 <- ggplot(df_subtel_only, aes(x = Contrast, y = log2FoldChange)) + geom_boxplot(alpha = 0.7, aes(fill = Contrast)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("gray", "#91cf50")) + stat_pvalue_manual(data = stat_test_subtel, label = "label", tip.length = 0.01) +
  theme_classic() + labs(y = "log2FC[8h/2h]")

ggsave("./DESeq2_results/plots/boxplot_subtel.pdf", plot = p3, device = "pdf", height = 4, width = 4)


# [Plot 4] Subtelomere Chromosome Heatmap
stat_test_chr <- df_all_noMito %>% group_by(Contrast, chr_name) %>% t_test(log2FoldChange ~ subtel) %>% adjust_pvalue(method = "none") %>% add_significance()
stat_test_all <- df_all_noMito %>% group_by(Contrast) %>% t_test(log2FoldChange ~ subtel) %>% adjust_pvalue(method = "none") %>% add_significance() %>% mutate(chr_name = "All")
df_combined_stat <- bind_rows(stat_test_chr, stat_test_all)

df_means_chr <- df_all_noMito %>% group_by(Contrast, chr_name, subtel) %>% summarise(mean_val = mean(log2FoldChange, na.rm = TRUE), .groups = "drop") %>% pivot_wider(names_from = subtel, values_from = mean_val, names_prefix = "mean_log2FC_")
df_means_all <- df_all_noMito %>% group_by(Contrast, subtel) %>% summarise(mean_val = mean(log2FoldChange, na.rm = TRUE), .groups = "drop") %>% pivot_wider(names_from = subtel, values_from = mean_val, names_prefix = "mean_log2FC_") %>% mutate(chr_name = "All")
df_combined_stat <- df_combined_stat %>% left_join(bind_rows(df_means_chr, df_means_all), by = c("Contrast", "chr_name"))

# ★ Change 1: Swapped x and y specifications
p4 <- ggplot(df_combined_stat, aes(x = Contrast, y = chr_name, fill = mean_log2FC_subtel)) +
  geom_tile(color = "white") + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, oob = scales::squish, name = "Mean log2FC\n(Subtel)") +
  
  # Display X-axis (Contrast) on top, and Y-axis (chr_name) from top to bottom using rev
  scale_x_discrete(position = "top") + 
  scale_y_discrete(limits = rev) +
  
  geom_text(aes(label = ifelse(p.adj.signif == "ns", "", p.adj.signif)), color = "black", size = 6) +
  theme_minimal() + 
  theme(
    # If Contrast text is too long and overlaps, adjust here (e.g., angle = 45, hjust = 0)
    axis.text.x = element_text(angle = 0, hjust = 0.5), 
    panel.grid = element_blank(), 
    axis.title.x = element_blank(), 
    axis.title.y = element_blank()
  )
p4
# ★ Change 2: Changed save size from horizontal (2x12) to vertical (e.g., 10x4)
ggsave("./DESeq2_results/plots/heatmap_subtel_chr.pdf", plot = p4, device = "pdf", height = 4, width = 4)

# [Plot 5] Retrotransposon (Ty1/Ty2) Boxplot
GO0032197 <- read_xlsx("../yeast/reference/GO0032197_retrotransposition.xlsx") %>% distinct(gene, .keep_all = TRUE)
df_all <- left_join(df_all, GO0032197, by = "gene")
df_all$GO0032197[is.na(df_all$GO0032197)] <- "NO"

df_ty1_ty2 <- subset(df_all, GO0032197 %in% c("Ty1", "Ty2"))
df_ty1_ty2$GO0032197 <- factor(df_ty1_ty2$GO0032197, levels = c("Ty1", "Ty2"))

stat_test_ty <- df_ty1_ty2 %>% group_by(GO0032197) %>% t_test(log2FoldChange ~ Contrast) %>%
  adjust_pvalue(method = "none") %>% add_significance() %>% add_xy_position(fun = "max", x = "Contrast", dodge = 0.9) %>% mutate(label = paste0("p = ", signif(p.adj, 2)))

p5 <- ggplot(df_ty1_ty2, aes(x = Contrast, y = log2FoldChange)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) + geom_boxplot(alpha = 0.7, aes(fill = Contrast), outlier.colour = "gray") +
  stat_pvalue_manual(data = stat_test_ty, label = "label", tip.length = 0.01) + scale_fill_manual(values = c("gray", "#91cf50")) +
  facet_wrap(~GO0032197) + theme_classic() +
  #coord_cartesian(ylim = c(-1, 2)) + 
  labs(y = "log2FC[8h/2h]")

ggsave("./DESeq2_results/plots/boxplot_retrotranspon_ty1_ty2.pdf", plot = p5, width = 4.5, height = 4)


# [Plot 6] iESR Boxplot
iESR <- read_xlsx("../yeast/reference/iESR_cluster_ESR_clusters_UPDATED_2017_from_SGD.xlsx") %>% distinct(gene, .keep_all = TRUE)
df_all <- left_join(df_all, iESR, by = "gene")
df_all$ESR[is.na(df_all$ESR)] <- "NO"
df_iESR <- subset(df_all, ESR != "NO")

stat_test_iesr <- df_iESR %>% group_by(ESR) %>% t_test(log2FoldChange ~ Contrast) %>%
  adjust_pvalue(method = "none") %>% add_significance() %>% add_xy_position(fun = "max", x = "Contrast", dodge = 0.9) %>% mutate(label = paste0("p = ", signif(p.adj, 2)))

p6 <- ggplot(df_iESR, aes(x = Contrast, y = log2FoldChange)) + geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(alpha = 0.7, aes(fill = Contrast)) + stat_pvalue_manual(data = stat_test_iesr, label = "label", tip.length = 0.01) +
  scale_fill_manual(values = c("gray", "#91cf50")) + theme_classic() + facet_wrap(~ESR) + labs(y = "log2FC[8h/2h]")

ggsave("./DESeq2_results/plots/boxplot_iESR.pdf", plot = p6, width = 4, height = 4)

cat("\n=== EVERYTHING COMPLETED AWESOMELY IN ONE SHOT! ===\n")