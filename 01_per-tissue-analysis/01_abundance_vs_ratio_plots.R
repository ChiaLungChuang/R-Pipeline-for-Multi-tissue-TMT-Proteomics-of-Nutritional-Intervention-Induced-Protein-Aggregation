# ============================================================
# Project: Multi-tissue Proteomics — Dietary Intervention
# Script:  01_abundance_vs_ratio_plots.R
# Purpose: Per-tissue scatter plots of protein Abundance fold
#          change (Treatment vs. Control) vs. Insolubility
#          Ratio fold change (Treatment vs. Control). Point
#          size encodes -log10(p-value) of ratio significance.
#          Generates labeled (extreme values) and R-squared
#          variants. Also generates a combined multi-tissue
#          overlay plot.
# Input:   Per-tissue Excel files with two sheets each:
#            (Ratio_INSOL_SOL) Treatment vs Control
#            (SUM_INSOL_SOL) Treatment vs Control
# Output:  Treatment_<Tissue>_extreme.pdf
#          Treatment_<Tissue>_extreme_R.pdf
#          Treatment_<Tissue>_combined.csv
#          Treatment_<Tissue>_processed_all.csv
#          Treatment_multi_tissue_comparison.pdf
# Author:  Chia-Lung Chuang
# Updated: 2026-03
# NOTE:    Scripts only — raw data not included. Underlying
#          dataset is tied to a manuscript in preparation.
# ============================================================

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpmisc)
library(readxl)
library(scales)

# NOTE: Set your working directory to the folder containing the input files
# setwd("path/to/your/data")

# ============================================================
# Tissue Color Palette
# ============================================================

tissue_colors <- list(
  tissue_a = "darkgreen",
  tissue_b = "blue4",
  tissue_c = "purple4",
  tissue_d = "darkgoldenrod",
  tissue_e = "darkred"
)

# ============================================================
# Core Function: Per-Tissue Abundance vs. Ratio Plot
# ============================================================

# Reads the Ratio (insolubility ratio) and Sum (total abundance)
# sheets from a per-tissue Excel file, merges them, and generates:
#   1. Scatter plot with extreme-value gene labels
#   2. Same plot with Pearson R² annotation
# Significance threshold: p <= 0.05 for insolubility ratio.
# Point size: -log10(p-value), minimum 0.5 for non-significant.
# Labeling: significant proteins with |ratio FC| > 1 or
#   |abundance FC| > defined thresholds.

create_tissue_plots <- function(file_path, tissue_name, tissue_color,
                                label_ratio_up   =  1.0,
                                label_ratio_down = -1.0,
                                label_sum_up     =  0.4,
                                label_sum_down   = -1.0) {
  
  # --- Load and merge sheets ---
  df_ratio <- read_excel(file_path,
                         sheet = "(Ratio_INSOL_SOL) Treatment vs Control") %>%
    select(AC, GN, DE, `Log2Fold(Control/Treatment)`, `p-value`) %>%
    rename(Log2FC_ratio  = `Log2Fold(Control/Treatment)`,
           p_value_ratio = `p-value`) %>%
    mutate(AC_GN = paste(AC, GN, sep = "_"))
  
  df_sum <- read_excel(file_path,
                       sheet = "(SUM_INSOL_SOL) Treatment vs Control") %>%
    select(AC, GN, DE, `Log2Fold(Control/Treatment)`, `p-value`) %>%
    rename(Log2FC_sum  = `Log2Fold(Control/Treatment)`,
           p_value_sum = `p-value`) %>%
    mutate(AC_GN = paste(AC, GN, sep = "_")) %>%
    select(-AC, -GN, -DE)
  
  final_df <- full_join(df_ratio, df_sum, by = "AC_GN")
  write.csv(final_df,
            paste0("Treatment_", toupper(tissue_name), "_combined.csv"),
            row.names = FALSE)
  
  # --- Prepare plot data ---
  plot_data <- final_df %>%
    mutate(
      p_value_ratio     = as.numeric(p_value_ratio),
      p_value_sum       = as.numeric(p_value_sum),
      significant_ratio = p_value_ratio <= 0.05,
      significant_sum   = p_value_sum   <= 0.05,
      size_ratio        = ifelse(p_value_ratio <= 0.05, -log10(p_value_ratio), 0.5),
      Gene_Name         = GN,
      point_alpha       = ifelse(significant_ratio, 0.4, 0.05),
      Label_extreme     = case_when(
        is.na(GN) | GN == "NA"  ~ NA_character_,
        significant_ratio &
          (Log2FC_ratio >  label_ratio_up  |
             Log2FC_ratio <  label_ratio_down |
             Log2FC_sum   >  label_sum_up    |
             Log2FC_sum   <  label_sum_down)  ~ GN,
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!grepl("co|CON_", AC_GN))
  
  write.csv(plot_data,
            paste0("Treatment_", toupper(tissue_name), "_processed_all.csv"),
            row.names = FALSE)
  
  # --- Base plot ---
  base_plot <- ggplot(plot_data,
                      aes(x = Log2FC_sum, y = Log2FC_ratio)) +
    geom_point(aes(size = size_ratio, alpha = point_alpha),
               color = tissue_color) +
    scale_alpha_identity() +
    scale_size_continuous(
      range  = c(0.4, 9),
      name   = "P-Value Size (-log10)",
      breaks = c(1.3, 2, 3, 4, 5),
      limits = c(0.5, 8),
      labels = c("1.30 (p<0.05)", "2", "3", "4", "5")
    ) +
    scale_x_continuous(labels = label_number(accuracy = 0.1)) +
    labs(
      title   = paste0(str_to_title(tissue_name),
                       ": Abundance FC vs. Insolubility FC (Treatment vs. Control)"),
      x       = "Log2 Abundance fold change (Treatment vs. Control)",
      y       = "Log2 Insolubility Ratio fold change (Treatment vs. Control)",
      caption = paste0("Point size: -log10(p-value) of insolubility ratio.\n",
                       "Small points (size=0.5) indicate p > 0.05.")
    ) +
    theme_minimal() +
    theme(
      axis.text   = element_text(size = 14, face = "bold"),
      axis.title  = element_text(size = 14, face = "bold"),
      legend.text  = element_text(size = 12, face = "italic"),
      legend.title = element_text(size = 14, face = "bold"),
      plot.title   = element_text(size = 16, face = "bold")
    )
  
  # --- Labeled variant ---
  p_extreme <- base_plot +
    geom_text_repel(
      data         = function(x) subset(x, !is.na(Label_extreme)),
      aes(label    = Label_extreme),
      size          = 5, fontface = "bold",
      segment.color = "grey50",
      max.overlaps  = Inf,
      nudge_x = 0.15, nudge_y = 0.15,
      box.padding = 0.5, point.padding = 0.5,
      min.segment.length = 0.1,
      na.rm = TRUE
    )
  ggsave(paste0("Treatment_", toupper(tissue_name), "_extreme.pdf"),
         p_extreme, width = 7200, height = 6400, units = "px", dpi = 600)
  
  # --- R-squared variant ---
  p_extreme_R <- p_extreme +
    stat_poly_eq(
      aes(label = paste(after_stat(rr.label))),
      rr.digits = 3, parse = TRUE,
      color = "black", show.legend = FALSE,
      label.y = 0.05, label.x = 0.9
    )
  ggsave(paste0("Treatment_", toupper(tissue_name), "_extreme_R.pdf"),
         p_extreme_R, width = 7200, height = 6400, units = "px", dpi = 600)
  
  return(list(
    plot_data    = plot_data,
    base_plot    = base_plot,
    extreme_plot = p_extreme,
    r2_plot      = p_extreme_R
  ))
}

# ============================================================
# Run Per-Tissue Analysis
# ============================================================

tissue_a_results <- create_tissue_plots(
  "Treatment_TissueA.xlsx", "tissue_a", tissue_colors$tissue_a)

tissue_b_results <- create_tissue_plots(
  "Treatment_TissueB.xlsx", "tissue_b", tissue_colors$tissue_b)

tissue_c_results <- create_tissue_plots(
  "Treatment_TissueC.xlsx", "tissue_c", tissue_colors$tissue_c)

tissue_d_results <- create_tissue_plots(
  "Treatment_TissueD.xlsx", "tissue_d", tissue_colors$tissue_d)

tissue_e_results <- create_tissue_plots(
  "Treatment_TissueE.xlsx", "tissue_e", tissue_colors$tissue_e)

# ============================================================
# Multi-Tissue Overlay Plot
# ============================================================

# Combines all five tissues into a single scatter plot.
# Non-significant proteins are rendered gray; significant
# proteins are colored by tissue. A single R² value is
# computed across all proteins in all tissues.

create_multi_tissue_plot <- function(tissue_results_list,
                                     tissue_names,
                                     tissue_colors) {
  tissue_labels <- setNames(
    paste("Tissue", LETTERS[seq_along(tissue_names)]),
    tissue_names
  )
  
  combined <- bind_rows(lapply(seq_along(tissue_results_list), function(i) {
    tissue_results_list[[i]]$plot_data %>%
      mutate(
        Tissue      = tissue_names[i],
        point_alpha = ifelse(significant_ratio, 0.2, 0.01)
      )
  }))
  
  r2_all <- round(summary(lm(Log2FC_ratio ~ Log2FC_sum, data = combined))$r.squared, 3)
  
  ggplot(combined, aes(x = Log2FC_sum, y = Log2FC_ratio)) +
    geom_point(aes(size = size_ratio, alpha = point_alpha, color = Tissue)) +
    annotate("text", x = -2.5, y = 2.5,
             label = paste0("R² = ", r2_all),
             hjust = 0, size = 5, fontface = "bold") +
    scale_alpha_identity() +
    scale_size_continuous(
      range  = c(0.4, 5),
      name   = "P-Value Size (-log10)",
      breaks = c(1.3, 2, 3, 4, 5),
      limits = c(0.5, 8),
      labels = c("1.30 (p<0.05)", "2", "3", "4", "5")
    ) +
    scale_color_manual(
      values = unlist(tissue_colors[tissue_names]),
      labels = tissue_labels,
      name   = "Tissue"
    ) +
    labs(
      title   = "All Tissues: Abundance FC vs. Insolubility FC (Treatment vs. Control)",
      x       = "Log2 Abundance fold change (Treatment vs. Control)",
      y       = "Log2 Insolubility Ratio fold change (Treatment vs. Control)",
      caption = paste0("Point size: -log10(p-value). ",
                       "Small points: p > 0.05.")
    ) +
    theme_minimal() +
    theme(
      axis.text    = element_text(size = 12, face = "bold"),
      axis.title   = element_text(size = 18, face = "bold"),
      legend.text  = element_text(size = 14),
      legend.title = element_text(size = 16, face = "bold"),
      plot.title   = element_text(size = 24, face = "bold")
    )
}

multi_tissue_plot <- create_multi_tissue_plot(
  list(tissue_a_results, tissue_b_results, tissue_c_results,
       tissue_d_results, tissue_e_results),
  c("tissue_a", "tissue_b", "tissue_c", "tissue_d", "tissue_e"),
  tissue_colors
)

ggsave("Treatment_multi_tissue_comparison.pdf",
       multi_tissue_plot,
       width = 10000, height = 8000, units = "px", dpi = 600)