library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
# setwd("/Volumes/biol_bc_neurohr_scratch_1/sosawa/NGS/Ohsawa_et_al_2026_ChIPseq/polII_ChIP/results//")

# ==============================================================
# 1. Set the target genomic region to extract
# ==============================================================
target_region <- GRanges(
  seqnames = "XII",                  # Chromosome name
  ranges = IRanges(start = 440000, end = 480001) # Target coordinates
)

# ==============================================================
# 2. Load BigWig files of all individual samples and synthesize a "Pooled (Average)" waveform in R!
# ==============================================================
# Directory containing BigWig files (adjust according to your environment)
bw_dir <- "../bigwigs_cpm/"

strains <- c("3225", "3229", "3506")
treatments <- c("EtOH", "EST")
mol_types <- c("IP", "input") # * Note the lowercase if created as input_3225_... in the terminal!

df_list <- list()

for (st in strains) {
  for (tr in treatments) {
    for (mt in mol_types) {
      
      # Search for corresponding individual replicate files
      # e.g., IP_3225_EtOH_1_CPM.bw, input_3229_EST_2_CPM.bw, etc.
      pattern <- paste0("^", mt, "_", st, "_", tr, "_.*_CPM\\.bw$")
      target_files <- list.files(path = bw_dir, pattern = pattern, full.names = TRUE)
      
      # ⚠️ Exclude the outlier 3506_EtOH_2 from the search results!
      target_files <- target_files[!grepl("3506_EtOH_2", target_files)]
      
      if (length(target_files) > 0) {
        
        # Load all corresponding replicates into a list
        rep_data_list <- lapply(target_files, function(file) {
          bw_data <- import(file, which = target_region)
          as.data.frame(bw_data) %>% dplyr::select(seqnames, start, end, score)
        })
        
        # Vertically bind data from all replicates and calculate the mean for each coordinate (start)!
        # This results in a "perfect average waveform", equivalent to creating a pooled.bw
        df_pooled <- bind_rows(rep_data_list) %>%
          group_by(seqnames, start, end) %>%
          summarise(score = mean(score, na.rm = TRUE), .groups = "drop") %>%
          mutate(
            Strain = st, 
            Treatment = tr, 
            MoleculeType = ifelse(mt == "input", "Input", "IP") # Capitalize Input for display
          )
        
        # Store in the list
        df_list[[paste0(st, "_", tr, "_", mt)]] <- df_pooled
        
      } else {
        warning(paste("File not found:", st, tr, mt))
      }
    }
  }
}

# Combine into a single massive dataframe
df_all_raw <- bind_rows(df_list)

# ==============================================================
# 3. Data processing and "Forced Y-axis upper limit cut" (Flat-topping)
# ==============================================================
df_plot_final <- df_all_raw %>%
  mutate(
    StrainIdx = factor(Strain, levels = c("3225", "3229", "3506")),
    TypeIdx = factor(MoleculeType, levels = c("IP", "Input")),
    TreatmentIdx = factor(Treatment, levels = c("EtOH", "EST")),
    
    # ★ Point 1: Cut spikes (abnormal peaks) exceeding the specified limit here ★
    # This prevents the entire graph from being squashed by local peaks
    score = case_when(
      MoleculeType == "IP" & score > 500 ~ 500,
      MoleculeType == "Input" & score > 1000 ~ 1000,
      # Align values below zero (if any) to zero
      score < 0 ~ 0,
      TRUE ~ score
    )
  )

# ==============================================================
# 4. Create dummy data (Blank Data) for Y-axis scale expansion
# ==============================================================
# Create "invisible points" to forcefully expand the graph scale to 500 or 1000
df_blank_max <- expand.grid(
  StrainIdx = factor(c("3225", "3229", "3506"), levels = c("3225", "3229", "3506")),
  TypeIdx = factor(c("IP", "Input"), levels = c("IP", "Input")),
  TreatmentIdx = factor("EtOH", levels = c("EtOH", "EST")) # Arbitrary factor
) %>%
  mutate(
    start = 450000, # Arbitrary coordinate within the X-axis range
    score = ifelse(TypeIdx == "IP", 500, 1000) # Plot points at 500 for IP and 1000 for Input
  )

# Create Y=0 dummy data and combine
df_blank_0 <- df_blank_max %>% mutate(score = 0)
df_blank <- bind_rows(df_blank_max, df_blank_0)

# ==============================================================
# 5. Draw with ggplot (make Y-axes independent using facet_wrap)
# ==============================================================
p_dual_overlay <- ggplot(mapping = aes(x = start, y = score)) +
  
  # 1. [Background] Draw EST data with filled areas (specify 3 colors)
  geom_area(
    data = filter(df_plot_final, TreatmentIdx == "EST"),
    aes(fill = StrainIdx), 
    alpha = 0.8, 
    position = "identity"
  ) +
  
  # 2. [Foreground] Draw EtOH data with black lines
  geom_line(
    data = filter(df_plot_final, TreatmentIdx == "EtOH"),
    color = "black", 
    size = 0.6
  ) +
  
  # Dummy data for Y-axis scale expansion (IP:0-500, Input:0-1000)
  geom_blank(data = df_blank, aes(x = start, y = score)) +
  
  # Zero baseline
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 0.4) +
  
  # ★ Specified layout: 1 column vertically, independent Y-axes, labels on the right ★
  facet_wrap(~ TypeIdx + StrainIdx, ncol = 1, scales = "free_y", strip.position = "right") + 
  
  # Color settings
  scale_fill_manual(values = c("3225" = "#3853a3", "3229" = "#bebebe", "3506" = "#ffa500")) +
  
  # Design adjustments
  theme_classic() +
  labs(
    x = "Genomic Position (Chromosome XII)",
    y = "Pol II ChIP-seq Signal (Averaged CPM)",
    fill = "Strain (EST)"
  ) +
  
  # ★ Specified theme: White background and vertical text for labels ★
  theme(
    strip.background = element_rect(fill = "white", color = "gray20"),
    # angle = -90 makes the text easier to read from top to bottom
    strip.text.y = element_text(face = "bold", size = 11, angle = -90), 
    panel.spacing = unit(1, "lines")
    # legend.position = "bottom" is commented out as specified (displays on the right by default)
  ) +
  scale_x_continuous(labels = scales::comma, expand = c(0, 0))+
  
  
  # Thin out the displayed ticks to about 3 (e.g., 0, 500, 1000) to prevent Y-axis overlap (1000, 750, 500...)
  scale_y_continuous(breaks = scales::pretty_breaks(n =1))
print(p_dual_overlay)

# Save as PDF
ggsave("IGS2_track.pdf", plot = p_dual_overlay, width = 6, height = 4)