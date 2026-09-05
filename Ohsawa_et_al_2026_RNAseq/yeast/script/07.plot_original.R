library(DESeq2)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(rstatix)
library(biomaRt)
library(stringr)
library(clusterProfiler)
library(org.Sc.sgd.db)
library(enrichplot)

# ==============================================================================
# 1. Initial Setup and Directory Creation
# ==============================================================================
setwd("/Volumes/biol_bc_neurohr_scratch_1/sosawa/NGS/Ohsawa_et_al_2026_RNAseq/yeast/") 

dir.create("./yeast_DESeq2_results", showWarnings = FALSE)
dir.create("./yeast_DESeq2_results/plots/PCA", recursive = TRUE, showWarnings = FALSE)
dir.create("./yeast_DESeq2_results/plots/GO", recursive = TRUE, showWarnings = FALSE)
dir.create("./yeast_DESeq2_results/DE_results", recursive = TRUE, showWarnings = FALSE)
dir.create("./yeast_DESeq2_results/plots/Summary_Plots", recursive = TRUE, showWarnings = FALSE)

cat("=================================================================\n")
cat(" Yeast RNA-seq Full Pipeline (Global Factor Ordering v9)\n")
cat("=================================================================\n\n")

# ------------------------------------------------------------------------------
# 🎨 User Configuration Area
# ------------------------------------------------------------------------------
group_configs <- list(
  "A" = list(
    pca_colors   = c("blue", "gray", "orange", "darkgreen"),
    ma_colors    = c("blue", "gray", "red"),
    box_colors   = c("blue", "gray", "orange"),
    esr_colors   = c("#3853A3", "lightgray", "orange"),
    ref_contrast = "A3229_EST_vs_EtOH",
    ma_width = 8, ma_height = 4, box_tpm_width = 6, box_tpm_height = 4, heat_width = 12, heat_height = 3, box_sub_width = 4, box_sub_height = 4, box_ty_width = 5, box_ty_height = 4, box_esr_width = 5, box_esr_height = 4
  ),
  "B" = list(
    pca_colors   = c("blue", "lightblue"),
    ma_colors    = c("blue", "gray", "red"),
    box_colors   = c("blue", "lightblue"),
    esr_colors   = c("blue", "lightblue"),
    ref_contrast = "B1684_EST_vs_EtOH",
    ma_width = 6, ma_height = 4, box_tpm_width = 5, box_tpm_height = 4, heat_width = 12, heat_height = 2.5, box_sub_width = 3, box_sub_height = 4, box_ty_width = 4, box_ty_height = 4, box_esr_width = 4, box_esr_height = 4
  ),
  "C" = list(
    pca_colors   = c("blue", "gray", "red"),
    ma_colors    = c("blue", "gray", "red"),
    box_colors   = c("blue", "gray", "red"),
    esr_colors   = c("#3853A3", "lightgray", "red"),
    ref_contrast = "C3229_EST_vs_EtOH",
    ma_width = 8, ma_height = 4, box_tpm_width = 6, box_tpm_height = 4, heat_width = 12, heat_height = 3, box_sub_width = 4, box_sub_height = 4, box_ty_width = 5, box_ty_height = 4, box_esr_width = 5, box_esr_height = 4
  ),
  "D" = list(
    pca_colors   = c("lightgray", "darkgray", "#FFD166", "black", "#FF9F1C"),
    ma_colors    = c("blue", "gray", "red"),
    box_colors   = c("lightgray", "#FFD166", "darkgray", "#FF9F1C"),
    esr_colors   = c("gray", "#FFD166", "gray", "#FF9F1C"),
    ref_contrast = "t4h_EtOH_vs_t2h", 
    ma_width = 6, ma_height = 4, box_tpm_width = 8, box_tpm_height = 4, heat_width = 12, heat_height = 4, box_sub_width = 4, box_sub_height = 4, box_ty_width = 7, box_ty_height = 4, box_esr_width = 5.5, box_esr_height = 4
  )
)

# ==============================================================================
# 2. Global Data Loading (BiomaRt & IGS2 CPM)
# ==============================================================================
cat("[GLOBAL] Loading master counts and fetching Gene Info via BiomaRt...\n")
scerevisiae <- read.csv("./gene_counts/raw_counts_Scerevisiae.csv", header = TRUE, row.names = 1)
raw_counts_global <- round(scerevisiae)
colnames(raw_counts_global) <- sub("_Aligned.*", "", colnames(raw_counts_global))

unique_genes <- rownames(raw_counts_global)
ensembl <- biomaRt::useMart(biomart ="fungi_mart", dataset="scerevisiae_eg_gene", host="https://fungi.ensembl.org")
geneName <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description","chromosome_name", "start_position", "end_position","strand", "transcript_length"), filters = 'ensembl_gene_id', values = unique_genes, mart = ensembl)
geneName <- geneName[!duplicated(geneName$ensembl_gene_id),]
colnames(geneName) <- c("gene","gene_name","description","chr_name", "start_pos","end_pos","strand_orient", "transcript_length")

igs2_file <- "./gene_counts/IGS2_cpm.csv" 
if (!file.exists(igs2_file) && file.exists("./gene_counts/IGS2_subtel_sorted_cpm.csv")) igs2_file <- "./gene_counts/IGS2_subtel_sorted_cpm.csv"

if (file.exists(igs2_file)) {
  cat("[GLOBAL] Loading IGS2 CPM data...\n")
  df_igs2_long <- read.csv(igs2_file) %>% filter(region %in% c("IGS2-1", "IGS2-2")) %>% pivot_longer(cols = -region, names_to = "sample", values_to = "cpm") %>% mutate(sample = sub("_Aligned.*", "", sample))
} else {
  df_igs2_long <- NULL
  cat("[WARNING] IGS2 CPM file not found. Skipping IGS2 plots.\n")
}

# ==============================================================================
# 🚀 Master Loop for Groups (A, B, C, D)
# ==============================================================================
for (grp in names(group_configs)) {
  cat(paste0("\n\n=================================================================\n"))
  cat(paste0(" >>> Starting Pipeline for Group: ", grp, " <<<\n"))
  cat(paste0("=================================================================\n"))
  
  cfg <- group_configs[[grp]]
  
  target_samples <- grep(paste0("^", grp), colnames(raw_counts_global), value = TRUE)
  if (grp == "B") { target_samples <- grep("B448", target_samples, value = TRUE, invert = TRUE) }
  if (grp == "D") { target_samples <- grep("D37_4h_EST_2", target_samples, value = TRUE, invert = TRUE) }
  
  if (length(target_samples) == 0) { cat(paste0("No samples found for Group ", grp, "! Skipping...\n")); next }
  
  counts_target <- raw_counts_global[, target_samples, drop = FALSE]
  counts_target <- counts_target[rowSums(counts_target > 1) >= 1, ]
  
  if (grp == "D") {
    col_data <- data.frame(sample = target_samples, group = factor(sub("_[0-9]+$", "", sub("^D37_", "t", target_samples)), levels = c("t2h_EtOH", "t4h_EtOH", "t4h_EST", "t6h_EtOH", "t6h_EST")), replicate = factor(sub(".*_([0-9]+)$", "\\1", target_samples)), row.names = target_samples)
    col_data$group <- relevel(col_data$group, ref = "t2h_EtOH")
  } else {
    col_data <- data.frame(sample = target_samples, strain = factor(sub("_(EST|EtOH)_[0-9]+$", "", target_samples)), treatment = factor(sub(".*_(EST|EtOH)_[0-9]+$", "\\1", target_samples), levels = c("EtOH", "EST")), replicate = factor(sub(".*_(EST|EtOH)_([0-9]+)$", "\\2", target_samples)), row.names = target_samples)
    
    # 🚨 ここで PCA 用の順番も固定！ (B3225 を B1684 の左に設定)
    if (grp == "A") { col_data$strain <- factor(col_data$strain, levels = c("A3225", "A3229", "A3506")) }
    if (grp == "B") { col_data$strain <- factor(col_data$strain, levels = c("B3225", "B1684")) }
    if (grp == "C") { col_data$strain <- factor(col_data$strain, levels = c("C3225", "C3229", "C448")) }
  }
  
  gene_lengths_bp <- geneName$transcript_length[match(rownames(counts_target), geneName$gene)]
  gene_lengths_kb <- gene_lengths_bp / 1000
  rpk <- counts_target / gene_lengths_kb
  scale_factors <- colSums(rpk, na.rm = TRUE) / 1e6
  tpm_matrix <- t(t(rpk) / scale_factors)
  save(tpm_matrix, file = paste0("./yeast_DESeq2_results/", grp, "_tpm_matrix.RData"))
  
  # ==============================================================================
  # PART 3: DESeq2 Analysis & PCA Plotting
  # ==============================================================================
  cat(paste0("[", grp, "] Running DESeq2 & Generating PCA...\n"))
  if (grp == "D") {
    dds <- DESeqDataSetFromMatrix(countData = counts_target, colData = col_data, design = ~ group)
    dds <- DESeq(dds, quiet = TRUE); normalized_counts <- counts(dds, normalized = TRUE); vsd <- vst(dds, blind = FALSE)
    pcaData <- plotPCA(vsd, intgroup = c("group"), returnData = TRUE); percentVar <- round(100 * attr(pcaData, "percentVar"))
    pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = group)) + geom_point(size = 5) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + scale_color_manual(values = cfg$pca_colors) + theme_test()
  } else {
    dds <- DESeqDataSetFromMatrix(countData = counts_target, colData = col_data, design = ~ strain + treatment)
    dds <- DESeq(dds, quiet = TRUE); normalized_counts <- counts(dds, normalized = TRUE); vsd <- vst(dds, blind = FALSE)
    pcaData <- plotPCA(vsd, intgroup = c("strain", "treatment"), returnData = TRUE); percentVar <- round(100 * attr(pcaData, "percentVar"))
    pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = strain, shape = treatment)) + geom_point(size = 5) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + scale_color_manual(values = cfg$pca_colors) + theme_test()
  }
  ggsave(paste0("./yeast_DESeq2_results/plots/PCA/", grp, "_PCA_plot.pdf"), pca_plot, width = 6, height = 6)
  
  # ==============================================================================
  # PART 4: Differential Expression (DE) Analysis Output
  # ==============================================================================
  cat(paste0("[", grp, "] Executing Differential Expression Contrasts...\n"))
  if (grp == "D") {
    for (g in c("t4h_EtOH", "t4h_EST", "t6h_EtOH", "t6h_EST")) {
      res <- results(dds, contrast = c("group", g, "t2h_EtOH")); df <- as.data.frame(res); df$gene <- rownames(df)
      base_cols <- c("gene", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")
      df <- cbind(df[, base_cols], as.data.frame(normalized_counts[, rownames(subset(col_data, group == "t2h_EtOH")), drop=FALSE]), as.data.frame(normalized_counts[, rownames(subset(col_data, group == g)), drop=FALSE]))
      write.xlsx(df[order(df$padj), ], file = paste0("./yeast_DESeq2_results/DE_results/", grp, "_", g, "_vs_t2h.xlsx"))
    }
  } else {
    for (strain in unique(col_data$strain)) {
      strain_samples <- col_data[col_data$strain == strain, ]
      if (!"EST" %in% strain_samples$treatment || !"EtOH" %in% strain_samples$treatment) next
      dds_strain <- dds[, rownames(strain_samples)]; dds_strain$strain <- droplevels(dds_strain$strain); design(dds_strain) <- ~ treatment; dds_strain <- DESeq(dds_strain, quiet = TRUE)
      res <- results(dds_strain, contrast = c("treatment", "EST", "EtOH")); norm_use <- counts(dds_strain, normalized = TRUE)
      df <- as.data.frame(res); df$gene <- rownames(df); base_cols <- c("gene", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")
      df <- cbind(df[, base_cols], as.data.frame(norm_use[, rownames(subset(strain_samples, treatment == "EtOH")), drop=FALSE]), as.data.frame(norm_use[, rownames(subset(strain_samples, treatment == "EST")), drop=FALSE]))
      write.xlsx(df[order(df$padj), ], file = paste0("./yeast_DESeq2_results/DE_results/", grp, "_", strain, "_EST_vs_EtOH.xlsx"))
    }
  }
  
  # ==============================================================================
  # PART 5: Combined Dataframe Integration
  # ==============================================================================
  cat(paste0("[", grp, "] Building integrated master dataset dynamically...\n"))
  df_list <- list()
  if (grp == "D") {
    for (g in c("t4h_EtOH", "t4h_EST", "t6h_EtOH", "t6h_EST")) { file_path <- paste0("./yeast_DESeq2_results/DE_results/", grp, "_", g, "_vs_t2h.xlsx"); if(file.exists(file_path)) df_list[[g]] <- read.xlsx(file_path) %>% dplyr::select(1:6) %>% mutate(Contrast = paste0(g, "_vs_t2h")) }
  } else {
    for (strain in unique(col_data$strain)) { file_path <- paste0("./yeast_DESeq2_results/DE_results/", grp, "_", strain, "_EST_vs_EtOH.xlsx"); if(file.exists(file_path)) df_list[[strain]] <- read.xlsx(file_path) %>% dplyr::select(1:6) %>% mutate(Contrast = paste0(strain, "_EST_vs_EtOH")) }
  }
  df_nonspikein <- bind_rows(df_list)
  
  # 🚨 ここで全プロット用の順番も固定！ (B3225 を B1684 の左に設定)
  if (grp == "A") { df_nonspikein$Contrast <- factor(df_nonspikein$Contrast, levels = c("A3225_EST_vs_EtOH", "A3229_EST_vs_EtOH", "A3506_EST_vs_EtOH")) }
  if (grp == "B") { df_nonspikein$Contrast <- factor(df_nonspikein$Contrast, levels = c("B3225_EST_vs_EtOH", "B1684_EST_vs_EtOH")) }
  if (grp == "C") { df_nonspikein$Contrast <- factor(df_nonspikein$Contrast, levels = c("C3225_EST_vs_EtOH", "C3229_EST_vs_EtOH", "C448_EST_vs_EtOH")) }
  if (grp == "D") { df_nonspikein$Contrast <- factor(df_nonspikein$Contrast, levels = c("t4h_EtOH_vs_t2h", "t4h_EST_vs_t2h", "t6h_EtOH_vs_t2h", "t6h_EST_vs_t2h")) }
  
  df_nonspikein$significant <- ifelse(df_nonspikein$log2FoldChange >= 1 & df_nonspikein$padj <= 0.05, "up", ifelse(df_nonspikein$log2FoldChange <= -1 & df_nonspikein$padj <= 0.05, "down", "ns"))
  df_nonspikein$significant[is.na(df_nonspikein$significant)] <- "ns"
  df_nonspikein <- left_join(df_nonspikein, geneName, by="gene")
  df_nonspikein$chr_name <- factor(df_nonspikein$chr_name, levels = c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "Mito"))
  
  genome_size <- read.csv("./reference/sacCer3.chrom.sizes.csv"); df_nonspikein$subtel <- "no"
  for (i in seq_len(nrow(genome_size))) { chr <- genome_size$chr_name[i]; chr_len <- genome_size$size[i]; subtel_idx <- which(df_nonspikein$chr_name == chr & (df_nonspikein$start_pos <= 30000 | df_nonspikein$start_pos >= (chr_len - 30000))); df_nonspikein$subtel[subtel_idx] <- "subtel" }
  
  df_tpm_calc <- as.data.frame(tpm_matrix); df_tpm_calc$gene <- rownames(df_tpm_calc)
  if (grp == "D") {
    df_tpm_calc$mean_t2_EtOH <- rowMeans(df_tpm_calc[, grep("D37_2h_EtOH", colnames(df_tpm_calc)), drop=FALSE], na.rm=TRUE)
    df_final <- df_nonspikein %>% left_join(df_tpm_calc %>% dplyr::select(gene, mean_t2_EtOH), by = "gene") %>% mutate(meanTPM_EtOH = mean_t2_EtOH) %>% dplyr::select(-mean_t2_EtOH)
  } else {
    for(strain in unique(col_data$strain)){ col_pattern <- paste0(strain, "_EtOH"); matched_cols <- grep(col_pattern, colnames(df_tpm_calc), value = TRUE); if(length(matched_cols) > 0){ df_tpm_calc[[paste0("mean_", strain, "_EtOH")]] <- rowMeans(df_tpm_calc[, matched_cols, drop=FALSE], na.rm=TRUE) } }
    tpm_long <- df_tpm_calc %>% dplyr::select(gene, starts_with("mean_")) %>% pivot_longer(cols = starts_with("mean_"), names_to = "tmp_col", values_to = "meanTPM_EtOH") %>% mutate(strain_match = gsub("mean_|_EtOH", "", tmp_col))
    df_nonspikein <- df_nonspikein %>% mutate(strain_match = sub("_(EST_vs_EtOH|EST_vs_EtOH)", "", Contrast))
    df_final <- left_join(df_nonspikein, tpm_long, by = c("gene", "strain_match")) %>% dplyr::select(-strain_match)
  }
  
  
  # ==============================================================================
  # ★ 追加: 各グループ (特に Group A) の df_final を RData として保存 ★
  # ==============================================================================
  # 1. 各グループごとの名前で保存 (例: ./yeast_DESeq2_results/A_df_final.RData)
  save_path <- paste0("./yeast_DESeq2_results/", grp, "_df_final.RData")
  save(df_final, file = save_path)
  cat(paste0("[", grp, "] Successfully saved master dataframe to: ", save_path, "\n"))
  df_final <- df_final %>% filter(!is.na(log2FoldChange))
  # ==============================================================================
  # PART 6: Publication-Ready Summary Plots Generation
  # ==============================================================================
  cat(paste0("[", grp, "] Plotting Final 6 Biological Figures...\n"))
  plot_dir <- "./yeast_DESeq2_results/plots/Summary_Plots/"
  bg_data <- data.frame(xmin = c(0.5, 1.5, 2.5), xmax = c(1.5, 2.5, 3.5), fill_color = c("lightblue", "lightgray", "lightpink"))
  
  # 1. MA plot 
  label_data_TPM <- df_final %>% dplyr::count(Contrast, significant) %>% filter(significant %in% c("up", "down")) %>% mutate(label = paste0(significant, ": ", n), x = max(log10(df_final$meanTPM_EtOH + 1), na.rm=T) - 2, y = ifelse(significant == "up", max(df_final$log2FoldChange, na.rm=T) - 1, min(df_final$log2FoldChange, na.rm=T) + 1), color = ifelse(significant == "up", "red", "blue"))
  if (grp == "D") {
    df_4h <- df_final %>% filter(grepl("t4h", Contrast)); label_4h <- label_data_TPM %>% filter(grepl("t4h", Contrast))
    p1_4h <- ggplot(df_4h, aes(x = log10(meanTPM_EtOH + 1), y = log2FoldChange)) + geom_rect(data = bg_data, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill_color), alpha = 0.3, inherit.aes = FALSE) + scale_fill_identity() + geom_hline(yintercept = 0, linetype = "dashed") + geom_point(alpha = 0.5, size = 1, aes(color = significant)) + geom_smooth(method = lm, color = "lightblue")+ scale_color_manual(values = c("blue", "gray", "red")) + geom_text(data = label_4h, aes(x = x, y = y, label = label), color = label_4h$color, inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.5) + theme_classic() + facet_wrap(vars(Contrast), nrow = 1) + labs(x = "log10(meanTPM_2h + 1)", y = "log2FC[estradiol/EtOH]")
    ggsave(paste0(plot_dir, grp, "_MA_plot_TPM_significant_4h.pdf"), plot = p1_4h, width = cfg$ma_width, height = cfg$ma_height)
    df_6h <- df_final %>% filter(grepl("t6h", Contrast)); label_6h <- label_data_TPM %>% filter(grepl("t6h", Contrast))
    p1_6h <- ggplot(df_6h, aes(x = log10(meanTPM_EtOH + 1), y = log2FoldChange)) + geom_rect(data = bg_data, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill_color), alpha = 0.3, inherit.aes = FALSE) + scale_fill_identity() + geom_hline(yintercept = 0, linetype = "dashed") + geom_point(alpha = 0.5, size = 1, aes(color = significant)) + geom_smooth(method = lm, color = "lightblue")+ scale_color_manual(values = c("blue", "gray", "red")) + geom_text(data = label_6h, aes(x = x, y = y, label = label), color = label_6h$color, inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.5) + theme_classic() + facet_wrap(vars(Contrast), nrow = 1) + labs(x = "log10(meanTPM_2h + 1)", y = "log2FC[estradiol/EtOH]")
    ggsave(paste0(plot_dir, grp, "_MA_plot_TPM_significant_6h.pdf"), plot = p1_6h, width = cfg$ma_width, height = cfg$ma_height)
  } else {
    p1 <- ggplot(df_final, aes(x = log10(meanTPM_EtOH + 1), y = log2FoldChange)) + geom_rect(data = bg_data, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill_color), alpha = 0.3, inherit.aes = FALSE) + scale_fill_identity() + geom_hline(yintercept = 0, linetype = "dashed") +geom_point(alpha = 0.5, size = 1, aes(color = significant)) +geom_smooth(method = lm, color = "lightblue")+  scale_color_manual(values = c("blue", "gray", "red")) + geom_text(data = label_data_TPM, aes(x = x, y = y, label = label), color = label_data_TPM$color, inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.5) + theme_classic() + facet_wrap(vars(Contrast), nrow = 1) + labs(x = "log10(meanTPM_EtOH + 1)", y = "log2FC[estradiol/EtOH]")
    ggsave(paste0(plot_dir, grp, "_MA_plot_TPM_significant.pdf"), plot = p1, width = cfg$ma_width, height = cfg$ma_height)
  }
  
  # 2. Boxplot TPM binned 
  df_binned <- df_final %>% mutate(log10_TPM = log10(meanTPM_EtOH + 1)) %>% filter(log10_TPM > 0.5 & log10_TPM < 4.0) %>% mutate(TPM_bin = cut(log10_TPM, breaks = seq(0.5, 3.5, by = 1), include.lowest = TRUE)) %>% filter(!is.na(TPM_bin))
  df_binned_safe <- df_binned %>% group_by(TPM_bin, Contrast) %>% filter(n() >= 2) %>% ungroup()
  if (grp == "D") { d_comps <- list(c("t4h_EtOH_vs_t2h", "t4h_EST_vs_t2h"), c("t6h_EtOH_vs_t2h", "t6h_EST_vs_t2h")); stat_test_binned <- df_binned_safe %>% group_by(TPM_bin) %>% t_test(log2FoldChange ~ Contrast, comparisons = d_comps) %>% adjust_pvalue(method = "BH") %>% add_significance("p") %>% add_xy_position(x = "TPM_bin", dodge = 0.8) %>% mutate(label = paste0("p = ", signif(p.adj, 2)))
  } else { stat_test_binned <- df_binned_safe %>% droplevels() %>% group_by(TPM_bin) %>% filter(cfg$ref_contrast %in% Contrast & n_distinct(Contrast) >= 2) %>% t_test(log2FoldChange ~ Contrast, ref.group = cfg$ref_contrast) %>% adjust_pvalue(method = "BH") %>% add_significance("p") %>% add_xy_position(fun = "max", x = "TPM_bin", dodge = 0.8) %>% mutate(label = paste0("p = ", signif(p.adj, 2))) }
  p2 <- ggplot(df_binned, aes(x = TPM_bin, y = log2FoldChange, fill = Contrast)) + annotate("rect", xmin = 0.55, xmax = 1.45, ymin = -Inf, ymax = Inf, fill = "lightblue", alpha = 0.3) + annotate("rect", xmin = 1.55, xmax = 2.45, ymin = -Inf, ymax = Inf, fill = "lightgray", alpha = 0.3) + annotate("rect", xmin = 2.55, xmax = 3.45, ymin = -Inf, ymax = Inf, fill = "lightpink", alpha = 0.3) + geom_boxplot(outlier.size = 1, outlier.color = "gray", outlier.alpha = 0.7, alpha = 0.7, width = 0.7, position = position_dodge(0.8)) + geom_hline(yintercept = 0, linetype = "dashed") + scale_fill_manual(values = cfg$box_colors) + theme_classic() + labs(x = "Basal Expression Level (log10[TPM_EtOH + 1])", y = "log2FC[estradiol/EtOH]") + stat_pvalue_manual(stat_test_binned, label = "label", tip.length = 0.01, hide.ns = TRUE)
  ggsave(paste0(plot_dir, grp, "_boxplot_TPMbinned.pdf"), plot = p2, width = cfg$box_tpm_width, height = cfg$box_tpm_height)
  
  # 3. Subtelomere Heatmap
  df_nonspikein_noMito <- subset(df_final, chr_name != "Mito")
  df_safe_chr <- df_nonspikein_noMito %>% group_by(Contrast, chr_name, subtel) %>% filter(n() >= 2) %>% ungroup() %>% group_by(Contrast, chr_name) %>% filter(n_distinct(subtel) == 2) %>% ungroup()
  df_stat_chr <- df_safe_chr %>% droplevels() %>% group_by(Contrast, chr_name) %>% t_test(log2FoldChange ~ subtel) %>% adjust_pvalue(method = "none") %>% add_significance()
  df_safe_all <- df_nonspikein_noMito %>% group_by(Contrast, subtel) %>% filter(n() >= 2) %>% ungroup() %>% group_by(Contrast) %>% filter(n_distinct(subtel) == 2) %>% ungroup()
  df_stat_all <- df_safe_all %>% droplevels() %>% group_by(Contrast) %>% t_test(log2FoldChange ~ subtel) %>% adjust_pvalue(method = "none") %>% add_significance() %>% mutate(chr_name = "All")
  df_combined_stat <- bind_rows(df_stat_chr, df_stat_all); df_means_chr <- df_nonspikein_noMito %>% group_by(Contrast, chr_name, subtel) %>% summarise(mean_val = mean(log2FoldChange, na.rm = TRUE), .groups = "drop") %>% pivot_wider(names_from = subtel, values_from = mean_val, names_prefix = "mean_log2FC_"); df_means_all <- df_nonspikein_noMito %>% group_by(Contrast, subtel) %>% summarise(mean_val = mean(log2FoldChange, na.rm = TRUE), .groups = "drop") %>% pivot_wider(names_from = subtel, values_from = mean_val, names_prefix = "mean_log2FC_") %>% mutate(chr_name = "All")
  df_combined_stat <- left_join(df_combined_stat, bind_rows(df_means_chr, df_means_all), by = c("Contrast", "chr_name"))
  p3 <- ggplot(df_combined_stat, aes(x = chr_name, y = Contrast, fill = mean_log2FC_subtel)) + geom_tile(color = "white") + scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, oob = scales::squish, name = "Mean log2FC\n(Subtel)") + scale_x_discrete(position = "top") + scale_y_discrete(limits = rev) + geom_text(aes(label = ifelse(p.adj.signif == "ns", "", p.adj.signif)), color = "black", size = 6) + theme_minimal() + theme(axis.text.x = element_text(angle = 0, hjust = 0.5), panel.grid = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
  ggsave(paste0(plot_dir, grp, "_heatmap_subtel_chr.pdf"), plot = p3, device = "pdf", width = cfg$heat_width, height = cfg$heat_height)
  
  # 4. All Subtelomere Boxplot
  df_subtel_only <- subset(df_nonspikein_noMito, subtel == "subtel")
  if (grp == "D") {
    df_sub_4h <- df_subtel_only %>% filter(grepl("t4h", Contrast)) %>% droplevels(); stat_sub_4h <- df_sub_4h %>% group_by(Contrast) %>% filter(n() >= 2) %>% ungroup() %>% t_test(log2FoldChange ~ Contrast, ref.group = "t4h_EtOH_vs_t2h") %>% adjust_pvalue(method = "none") %>% add_significance() %>% add_xy_position(fun = "max", x = "Contrast", dodge = 0.9) %>% mutate(label = paste0("p = ", signif(p.adj, 2)))
    p4_4h <- ggplot(df_sub_4h, aes(x = Contrast, y = log2FoldChange)) + geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) + geom_boxplot(alpha = 0.7, aes(fill = Contrast), outlier.color = "gray", outlier.alpha = 0.7) + scale_fill_manual(values=cfg$box_colors[1:2]) + stat_pvalue_manual(data = stat_sub_4h, label = "label", tip.length = 0.01, bracket.size = 0.5) + theme_classic() + labs(y = "log2FC[estradiol/EtOH]", x = NULL); ggsave(paste0(plot_dir, grp, "_boxplot_subtel_4h.pdf"), plot = p4_4h, width = cfg$box_sub_width, height = cfg$box_sub_height)
    df_sub_6h <- df_subtel_only %>% filter(grepl("t6h", Contrast)) %>% droplevels(); stat_sub_6h <- df_sub_6h %>% group_by(Contrast) %>% filter(n() >= 2) %>% ungroup() %>% t_test(log2FoldChange ~ Contrast, ref.group = "t6h_EtOH_vs_t2h") %>% adjust_pvalue(method = "none") %>% add_significance() %>% add_xy_position(fun = "max", x = "Contrast", dodge = 0.9) %>% mutate(label = paste0("p = ", signif(p.adj, 2)))
    p4_6h <- ggplot(df_sub_6h, aes(x = Contrast, y = log2FoldChange)) + geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) + geom_boxplot(alpha = 0.7, aes(fill = Contrast), outlier.color = "gray", outlier.alpha = 0.7) + scale_fill_manual(values=cfg$box_colors[3:4]) + stat_pvalue_manual(data = stat_sub_6h, label = "label", tip.length = 0.01, bracket.size = 0.5) + theme_classic() + labs(y = "log2FC[estradiol/EtOH]", x = NULL); ggsave(paste0(plot_dir, grp, "_boxplot_subtel_6h.pdf"), plot = p4_6h, width = cfg$box_sub_width, height = cfg$box_sub_height)
  } else {
    df_sub_safe <- df_subtel_only %>% group_by(Contrast) %>% filter(n() >= 2) %>% ungroup(); stat_all_subtel <- df_sub_safe %>% droplevels() %>% filter(cfg$ref_contrast %in% Contrast & n_distinct(Contrast) >= 2) %>% t_test(log2FoldChange ~ Contrast, paired = TRUE, ref.group = cfg$ref_contrast) %>% adjust_pvalue(method = "none") %>% add_significance() %>% add_xy_position(fun = "max", x = "Contrast", dodge = 0.9) %>% mutate(label = paste0("p = ", signif(p.adj, 2)))
    p4 <- ggplot(df_subtel_only, aes(x = Contrast, y = log2FoldChange)) + geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) + geom_boxplot(alpha = 0.7, aes(fill = Contrast), outlier.color = "gray", outlier.alpha = 0.7) + scale_fill_manual(values=cfg$box_colors) + stat_pvalue_manual(data = stat_all_subtel, label = "label", tip.length = 0.01, bracket.size = 0.5) + theme_classic() + labs(y = "log2FC[estradiol/EtOH]", x = NULL); ggsave(paste0(plot_dir, grp, "_boxplot_subtel.pdf"), plot = p4, width = cfg$box_sub_width, height = cfg$box_sub_height)
  }
  
  # 5. Ty1/Ty2 Boxplot
  GO0032197 <- read.xlsx("./reference/GO0032197_retrotransposition.xlsx") %>% distinct(gene, .keep_all = TRUE); df_ty <- left_join(df_final, GO0032197, by = "gene"); df_ty$GO0032197[is.na(df_ty$GO0032197)] <- "NO"; df_ty1_ty2 <- subset(df_ty, GO0032197 %in% c("Ty1", "Ty2")); df_ty_safe <- df_ty1_ty2 %>% group_by(GO0032197, Contrast) %>% filter(n() >= 2) %>% ungroup()
  if (grp == "D") { d_comps <- list(c("t4h_EtOH_vs_t2h", "t4h_EST_vs_t2h"), c("t6h_EtOH_vs_t2h", "t6h_EST_vs_t2h")); stat_ty <- df_ty_safe %>% group_by(GO0032197) %>% t_test(log2FoldChange ~ Contrast, comparisons = d_comps) %>% adjust_pvalue(method = "none") %>% add_significance() %>% add_xy_position(x = "Contrast", dodge = 0.9) %>% mutate(label = paste0("p = ", signif(p.adj, 2)))
  } else { stat_ty <- df_ty_safe %>% droplevels() %>% group_by(GO0032197) %>% filter(cfg$ref_contrast %in% Contrast & n_distinct(Contrast) >= 2) %>% t_test(log2FoldChange ~ Contrast, paired = TRUE, ref.group = cfg$ref_contrast) %>% adjust_pvalue(method = "none") %>% add_significance() %>% add_xy_position(fun = "max", x = "Contrast", dodge = 0.9) %>% mutate(label = paste0("p = ", signif(p.adj, 2))) }
  p5 <- ggplot(df_ty1_ty2, aes(x = Contrast, y = log2FoldChange, fill = Contrast)) + geom_hline(yintercept = 0, linetype = "dashed") + geom_boxplot(outlier.size = 1, outlier.alpha = 0.7, alpha = 0.7, width = 0.7, position = position_dodge(0.8)) + stat_pvalue_manual(stat_ty, label = "label", tip.length = 0.01) + scale_fill_manual(values=cfg$box_colors) + facet_wrap(~GO0032197) + theme_classic() + labs(y = "log2FC[estradiol/EtOH]"); ggsave(paste0(plot_dir, grp, "_boxplot_retrotranspon_ty1_ty2.pdf"), plot = p5, width = cfg$box_ty_width, height = cfg$box_ty_height)
  
  # 6. iESR Boxplot
  iESR <- read.xlsx("./reference/iESR_cluster_ESR_clusters_UPDATED_2017_from_SGD.xlsx") %>% distinct(gene, .keep_all = TRUE); df_iesr <- left_join(df_final, iESR, by = "gene"); df_iesr$ESR[is.na(df_iesr$ESR)] <- "NO"; df_iesr_filtered <- subset(df_iesr, ESR != "NO")
  
  if (grp == "D") {
    # --- 4hの分離 ---
    df_iesr_4h <- df_iesr_filtered %>% filter(grepl("t4h", Contrast)) %>% droplevels(); stat_iesr_4h <- df_iesr_4h %>% group_by(ESR, Contrast) %>% filter(n() >= 2) %>% ungroup() %>% group_by(ESR) %>% filter("t4h_EtOH_vs_t2h" %in% Contrast & n_distinct(Contrast) >= 2) %>% t_test(log2FoldChange ~ Contrast, ref.group = "t4h_EtOH_vs_t2h") %>% adjust_pvalue(method = "none") %>% add_significance() %>% add_xy_position(fun = "max", x = "Contrast", dodge = 0.9) %>% mutate(label = paste0("p = ", signif(p.adj, 2)))
    p6_4h <- ggplot(df_iesr_4h, aes(x = Contrast, y = log2FoldChange)) + geom_hline(yintercept = 0, linetype = "dashed") + geom_boxplot(alpha = 0.7, aes(fill = Contrast)) + stat_pvalue_manual(stat_iesr_4h, label = "label", tip.length = 0.01) + scale_fill_manual(values=cfg$esr_colors[1:2]) + theme_classic() + facet_wrap(~ESR) + labs(y = "log2FC[estradiol/EtOH]"); ggsave(paste0(plot_dir, grp, "_boxplot_iESR_4h.pdf"), plot = p6_4h, width = cfg$box_esr_width, height = cfg$box_esr_height)
    
    # --- 6hの分離 ---
    df_iesr_6h <- df_iesr_filtered %>% filter(grepl("t6h", Contrast)) %>% droplevels(); stat_iesr_6h <- df_iesr_6h %>% group_by(ESR, Contrast) %>% filter(n() >= 2) %>% ungroup() %>% group_by(ESR) %>% filter("t6h_EtOH_vs_t2h" %in% Contrast & n_distinct(Contrast) >= 2) %>% t_test(log2FoldChange ~ Contrast, ref.group = "t6h_EtOH_vs_t2h") %>% adjust_pvalue(method = "none") %>% add_significance() %>% add_xy_position(fun = "max", x = "Contrast", dodge = 0.9) %>% mutate(label = paste0("p = ", signif(p.adj, 2)))
    p6_6h <- ggplot(df_iesr_6h, aes(x = Contrast, y = log2FoldChange)) + geom_hline(yintercept = 0, linetype = "dashed") + geom_boxplot(alpha = 0.7, aes(fill = Contrast)) + stat_pvalue_manual(stat_iesr_6h, label = "label", tip.length = 0.01) + scale_fill_manual(values=cfg$esr_colors[3:4]) + theme_classic() + facet_wrap(~ESR) + labs(y = "log2FC[estradiol/EtOH]"); ggsave(paste0(plot_dir, grp, "_boxplot_iESR_6h.pdf"), plot = p6_6h, width = cfg$box_esr_width, height = cfg$box_esr_height)
    
  } else {
    # --- D以外のグループ (通常処理) ---
    df_iesr_safe <- df_iesr_filtered %>% group_by(ESR, Contrast) %>% filter(n() >= 2) %>% ungroup(); stat_iesr <- df_iesr_safe %>% droplevels() %>% group_by(ESR) %>% filter(cfg$ref_contrast %in% Contrast & n_distinct(Contrast) >= 2) %>% t_test(log2FoldChange ~ Contrast, paired = TRUE, ref.group = cfg$ref_contrast) %>% adjust_pvalue(method = "none") %>% add_significance() %>% add_xy_position(fun = "max", x = "Contrast", dodge = 0.9) %>% mutate(label = paste0("p = ", signif(p.adj, 2)))
    p6 <- ggplot(df_iesr_filtered, aes(x = Contrast, y = log2FoldChange)) + geom_hline(yintercept = 0, linetype = "dashed") + geom_boxplot(alpha = 0.7, aes(fill = Contrast)) + stat_pvalue_manual(stat_iesr, label = "label", tip.length = 0.01) + scale_fill_manual(values=cfg$esr_colors) + theme_classic() + facet_wrap(~ESR) + labs(y = "log2FC[estradiol/EtOH]"); ggsave(paste0(plot_dir, grp, "_boxplot_iESR.pdf"), plot = p6, width = cfg$box_esr_width, height = cfg$box_esr_height)
  }
  # ==============================================================================
  # PART 7: GO Enrichment Analysis (Group A only)
  # ==============================================================================
  if (grp == "A") {
    cat(paste0("[", grp, "] Running GO Enrichment Analysis...\n"))
    go_plot_dir <- "./yeast_DESeq2_results/plots/GO/"
    for (current_contrast in unique(df_final$Contrast)) {
      df_contrast <- subset(df_final, Contrast == current_contrast)
      for (ont_type in c("BP", "MF", "CC")) {
        up_genes <- df_contrast$gene[df_contrast$significant == "up"]
        if (length(up_genes) > 0) { ego_up <- enrichGO(gene = up_genes, universe = df_contrast$gene, OrgDb = org.Sc.sgd.db, keyType = "ORF", ont = ont_type, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = FALSE); if (!is.null(ego_up) && nrow(ego_up) > 0) { p_go_up <- dotplot(ego_up, showCategory = 15, title = paste0("GO (Up) - ", ont_type, ": ", current_contrast)); ggsave(paste0(go_plot_dir, grp, "_GO_Up_", current_contrast, "_", ont_type, ".pdf"), plot = p_go_up, width = 8, height = 6); write.xlsx(as.data.frame(ego_up), file = paste0("./yeast_DESeq2_results/DE_results/", grp, "_GO_Up_", current_contrast, "_", ont_type, ".xlsx")) } }
        down_genes <- df_contrast$gene[df_contrast$significant == "down"]
        if (length(down_genes) > 0) { ego_down <- enrichGO(gene = down_genes, universe = df_contrast$gene, OrgDb = org.Sc.sgd.db, keyType = "ORF", ont = ont_type, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = FALSE); if (!is.null(ego_down) && nrow(ego_down) > 0) { p_go_down <- dotplot(ego_down, showCategory = 15, title = paste0("GO (Down) - ", ont_type, ": ", current_contrast)); ggsave(paste0(go_plot_dir, grp, "_GO_Down_", current_contrast, "_", ont_type, ".pdf"), plot = p_go_down, width = 8, height = 6); write.xlsx(as.data.frame(ego_down), file = paste0("./yeast_DESeq2_results/DE_results/", grp, "_GO_Down_", current_contrast, "_", ont_type, ".xlsx")) } }
      }
    }
  } else { cat(paste0("[", grp, "] Skipping GO Enrichment Analysis.\n")) }
  
  # ==============================================================================
  # PART 8: 📊 IGS2 CPM Barplots
  # ==============================================================================
  if (!is.null(df_igs2_long)) {
    cat(paste0("[", grp, "] Plotting IGS2 Barplots...\n"))
    df_grp_igs2 <- df_igs2_long %>% filter(sample %in% rownames(col_data)) %>% left_join(col_data, by = "sample")
    
    if (nrow(df_grp_igs2) > 0) {
      if (grp == "A") { df_grp_igs2$strain <- factor(df_grp_igs2$strain, levels = c("A3225", "A3229", "A3506")) }
      # 🚨 順番固定 (B3225 を B1684 の左に配置)
      if (grp == "B") { df_grp_igs2$strain <- factor(df_grp_igs2$strain, levels = c("B3225", "B1684")) }
      if (grp == "C") { df_grp_igs2$strain <- factor(df_grp_igs2$strain, levels = c("C3225", "C3229", "C448")) }
      
      if (grp == "D") {
        df_grp_igs2 <- df_grp_igs2 %>% mutate(time = sub("_.*", "", group), treatment = sub(".*_", "", group)) %>% filter(time != "t2h")
        df_grp_igs2$treatment <- factor(df_grp_igs2$treatment, levels = c("EtOH", "EST"))
        stat_igs2 <- df_grp_igs2 %>% group_by(region, time) %>% filter(n_distinct(treatment) >= 2) %>% t_test(cpm ~ treatment) %>% adjust_pvalue(method = "none") %>% add_significance() %>% add_xy_position(fun = "max", x = "time", dodge = 0.8) %>% mutate(label = paste0("p = ", signif(p, 2)))
        x_var <- "time"
      } else {
        stat_igs2 <- df_grp_igs2 %>% group_by(region, strain) %>% filter(n_distinct(treatment) >= 2) %>% t_test(cpm ~ treatment) %>% adjust_pvalue(method = "none") %>% add_significance() %>% add_xy_position(fun = "max", x = "strain", dodge = 0.8) %>% mutate(label = paste0("p = ", signif(p, 2)))
        x_var <- "strain"
      }
      
      p_igs2 <- ggplot(df_grp_igs2, aes(x = .data[[x_var]], y = cpm, fill = treatment)) +
        stat_summary(fun = mean, geom = "bar", position = position_dodge(0.8), width = 0.7, color = "black") +
        stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(0.8), width = 0.2) +
        geom_point(aes(group = treatment), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), show.legend = FALSE, alpha = 0.7, size = 2, color = "black") +
        scale_fill_manual(values = c("EtOH" = "gray", "EST" = "firebrick")) +
        facet_wrap(~region, scales = "free_x", space = "free_x") +
        theme_classic() +
        labs(y = "Normalized Counts (CPM)", x = ifelse(grp == "D", "Time", "Strain")) +
        theme(strip.background = element_rect(fill = "white", color = "black"), strip.text = element_text(face = "bold", size = 12))
      
      if(nrow(stat_igs2) > 0) {
        y_nudge <- max(df_grp_igs2$cpm, na.rm=TRUE) * 0.05
        p_igs2 <- p_igs2 + stat_pvalue_manual(stat_igs2, label = "label", tip.length = 0.01, bracket.nudge.y = y_nudge)
      }
      ggsave(paste0(plot_dir, grp, "_IGS2_Barplot_Stats.pdf"), plot = p_igs2, width = 5, height = 4)
    }
  }
}
cat("\n=== ALL PIPELINES COMPLETED AWESOMELY AND SUCCESSFULLY! ===\n")