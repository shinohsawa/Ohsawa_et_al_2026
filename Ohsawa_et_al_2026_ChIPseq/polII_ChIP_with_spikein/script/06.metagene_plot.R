library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
# setwd("/Volumes/biol_bc_neurohr_scratch_1/sosawa/NGS/Ohsawa_et_al_2026_ChIPseq/polII_ChIP/results/")

# ==========================================
# 1. Function to automatically read data per strain, [Normalize by BED regions], and average replicates
# ==========================================
# Added arguments n_etoh (number of EtOH replicates) and n_est (number of EST replicates)
process_strain <- function(file_path, strain_name, n_etoh, n_est) {
  
  # Read gz file directly (output from computeMatrix)
  mat <- fread(file_path, skip=1)
  
  # Number of bins per sample is fixed at 200
  bins <- 200
  
  # =========================================================
  # ★ BED region-based normalization (Median-centering per Replicate)
  # For each replicate (200 bins), calculate the median signal of all genes,
  # and subtract it from all values to align the baseline to zero.
  # =========================================================
  
  # 1. Normalize EtOH replicates
  for (r in 1:n_etoh) {
    cols_rep <- (7 + (r - 1) * bins) : (6 + r * bins)
    rep_median <- median(as.matrix(mat[, ..cols_rep]), na.rm = TRUE)
    for (col_name in names(mat)[cols_rep]) {
      set(mat, j = col_name, value = mat[[col_name]] - rep_median)
    }
  }
  
  # 2. Normalize EST replicates
  for (r in 1:n_est) {
    cols_rep <- (7 + (n_etoh * bins) + (r - 1) * bins) : (6 + (n_etoh * bins) + r * bins)
    rep_median <- median(as.matrix(mat[, ..cols_rep]), na.rm = TRUE)
    for (col_name in names(mat)[cols_rep]) {
      set(mat, j = col_name, value = mat[[col_name]] - rep_median)
    }
  }
  # =========================================================
  
  # Re-acquire indices for all EtOH columns and all EST columns
  idx_etoh <- 7:(6 + (n_etoh * bins))
  idx_est  <- (7 + (n_etoh * bins)) : (6 + (n_etoh * bins) + (n_est * bins))
  
  # Sort in descending order by the "mean of all replicates" for the normalized EtOH condition
  mat$EtOH_mean <- rowMeans(mat[, ..idx_etoh], na.rm = TRUE)
  df_genome <- mat
  
  df_genome <- df_genome[order(-EtOH_mean)]
  
  # --- Internal Function: Extract specific % and convert to long format ---
  get_profile <- function(start_pct, end_pct, label_name) {
    n <- nrow(df_genome)
    start_row <- max(1, floor(n * (start_pct / 100)))
    end_row <- min(n, ceiling(n * (end_pct / 100)))
    
    subset_df <- df_genome[start_row:end_row, ]
    col_means_all <- colMeans(subset_df[, 7:ncol(subset_df)], na.rm = TRUE)
    
    # ★ Average the corresponding bins across replicates ★
    etoh_bin_means <- numeric(bins)
    est_bin_means  <- numeric(bins)
    
    for (b in 1:bins) {
      rep_idx_etoh <- sapply(1:n_etoh, function(r) (r - 1) * bins + b)
      etoh_bin_means[b] <- mean(col_means_all[rep_idx_etoh], na.rm = TRUE)
      
      rep_idx_est <- sapply(1:n_est, function(r) (n_etoh * bins) + (r - 1) * bins + b)
      est_bin_means[b] <- mean(col_means_all[rep_idx_est], na.rm = TRUE)
    }
    
    data.frame(
      Strain = strain_name,
      Group = label_name,
      Bin = 1:bins,
      `EtOH` = etoh_bin_means,
      `EST`  = est_bin_means,
      check.names = FALSE
    ) %>%
      pivot_longer(cols = c(`EtOH`, `EST`),
                   names_to = "Treatment", values_to = "Occupancy")
  }
  
  bind_rows(
    get_profile(0, 10, "1 Top 10%"),
    get_profile(10, 50, "2. Mid 10-50%"),
    get_profile(50, 100, "3. Bottom 50%")
  )
}

# ==========================================
# 2. Read data for 3 strains at once and combine into a single dataframe
# ==========================================
# ★ Specify Spike-in normalized _reps.gz here (pass the number of replicates to the function)
plot_data <- bind_rows(
  process_strain("../TSS_Metagene_Occupancy/matrix_3225_reps.gz", "A3225", n_etoh = 3, n_est = 3),
  process_strain("../TSS_Metagene_Occupancy/matrix_3229_reps.gz", "A3229",    n_etoh = 3, n_est = 3),
  process_strain("../TSS_Metagene_Occupancy/matrix_3506_reps.gz", "A3506", n_etoh = 2, n_est = 3) # EtOH has 2 replicates!
)

# Fix the order (Factor levels) for Treatment and Strain
plot_data$Treatment <- factor(plot_data$Treatment, levels=c("EtOH", "EST"))
plot_data$Strain <- factor(plot_data$Strain, levels=c("A3225", "A3229", "A3506"))

# ==========================================
# 3. Set X-axis labels (-0.5kb, start, stop, +0.5kb)
# ==========================================
x_breaks <- c(1, 50, 150, 200)
x_labels <- c("-0.5 kb", "start", "stop", "+0.5 kb")

# ==========================================
# 4. Draw horizontal x vertical panel plots (facet_grid) using ggplot
# ==========================================
p <- ggplot(plot_data, aes(x = Bin, y = Occupancy, color = Strain, alpha = Treatment)) +
  
  # Fill the background from -500bp (Bin 1) to start (Bin 50) with light gray
  annotate("rect", xmin = 1, xmax = 50, ymin = -Inf, ymax = Inf, fill = "gray80", alpha = 0.4) +
  
  # Draw dashed lines at start (Bin 50) and stop (Bin 150)
  geom_vline(xintercept = c(50, 150), linetype = "dashed", color = "gray40", size = 0.5) +
  
  # Main waveform data
  geom_line(size = 1.0) +
  
  # Set colors and transparency
  scale_color_manual(values = c("#3853a3", "#bebebe", "#ffa500")) +
  scale_alpha_manual(values = c(0.5, 1)) +
  
  scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = c(0, 0)) +
  
  # ★ Remove limits = c(0, NA) since normalization produces negative values
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  
  # Arrange Group vertically and Strain horizontally
  facet_grid(Group ~ Strain, scales = "free_y") + 
  
  # Overall design theme
  theme_classic() +
  labs(
    x = "",
    y = "Normalized Pol II log2FC\n(Median-centered per BED regions)"
  ) +
  theme(
    legend.title = element_blank(),
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

print(p)

# Save as PDF
ggsave("./metagene_plot_percentiles.pdf", plot = p, width = 7.5, height = 4)