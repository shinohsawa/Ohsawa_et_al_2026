library(dplyr)
library(tidyr)
library(rstatix)
library(ggpubr)
library(ggplot2)
library(stringr)

# setwd("/Volumes/biol_bc_neurohr_scratch_1/sosawa/NGS/Ohsawa_et_al_2026_RNAseq/yeast/") # Commented out for GitHub submission

# ==============================================================================
# 1. Data Loading and Formatting
# ==============================================================================

# Load data
df <- read.csv("gene_counts/IGS2_cpm.csv")
df_long <- df %>%
  dplyr::select(region, starts_with("D37_")) %>%
  pivot_longer(
    cols = starts_with("D37_"),
    names_to = "sample_name",
    values_to = "cpm"
  )

# Remove unnecessary strings from sample names and extract metadata
df_clean <- df_long %>%
  mutate(
    # Remove the "_Aligned.sortedByCoord.out" suffix
    sample_clean = str_replace(sample_name, "_Aligned\\.sortedByCoord\\.out$", "")
  ) %>%
  # Split the sample name (e.g., "D37_4h_EST_1") into 4 columns
  separate(sample_clean, into = c("Strain", "Time", "Condition", "Replicate"), sep = "_") %>%
  # Create a new Treatment column for the X-axis with the specified names
  mutate(
    Treatment = case_when(
      Time == "2h" & Condition == "EtOH" ~ "2h",
      Time == "4h" & Condition == "EtOH" ~ "4h",
      Time == "4h" & Condition == "EST"  ~ "4hEST",
      Time == "6h" & Condition == "EtOH" ~ "6h",
      Time == "6h" & Condition == "EST"  ~ "6hEST",
      TRUE ~ "Other"
    )
  )

# Convert to factor and specify the display order
df_clean$Treatment <- factor(df_clean$Treatment, levels = c("2h", "4h", "4hEST", "6h", "6hEST"))

# Verification
print("Data breakdown for analysis:")
print(table(df_clean$region, df_clean$Treatment))

# ==============================================================================
# 2. T-test Execution
# ==============================================================================
# Perform t-test by Treatment for each region (using 2h as the reference group)
stat.test <- df_clean %>%
  group_by(region) %>%
  t_test(cpm ~ Treatment, ref.group = "2h") %>%
  adjust_pvalue(method = "none") %>%
  add_significance() %>%
  add_xy_position(x = "Treatment", dodge = 0.8)

# Create P-value labels (e.g., p = 0.05)
stat.test <- stat.test %>%
  mutate(label = paste0("p = ", signif(p, 2)))

# ==============================================================================
# 3. Barplot Generation (ggplot2)
# ==============================================================================
p <- ggplot(df_clean, aes(x = Treatment, y = cpm, fill = Treatment)) +
  # Bar plot for mean values
  stat_summary(fun = mean, geom = "bar", position = position_dodge(0.8), width = 0.7, color = "black", alpha = 0.7) +
  # Error bars (Mean +/- Standard Error)
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(0.8), width = 0.2) +
  # Individual data points
  geom_point(position = position_dodge(0.8), show.legend = FALSE, alpha = 0.6) +
  
  # Add statistical analysis results (P-values)
  stat_pvalue_manual(
    stat.test,
    label = "label",
    tip.length = 0.01,
    step.increase = 0.08 # Stagger the brackets to prevent overlapping P-values
  ) +
  
  # Specified color palette
  scale_fill_manual(values = c("gray", "greenyellow", "green4", "hotpink", "maroon4")) +
  
  # Facet by region
  facet_wrap(~region, scales = "free_x", space = "free_x") +
  
  theme_classic() +
  labs(
    x = "Time & Treatment",
    y = "Normalized Counts (CPM)",
    fill = "Condition"
  ) +
  
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "none" # Hide legend for a cleaner look since colors map directly to the X-axis
  )

# Print and save the plot
print(p)
ggsave("./yeast_DESeq2_results/plots/Summary_Plots/D_IGS2_Barplot_Stats_with_2h.pdf", plot = p, width = 5.5, height = 6)