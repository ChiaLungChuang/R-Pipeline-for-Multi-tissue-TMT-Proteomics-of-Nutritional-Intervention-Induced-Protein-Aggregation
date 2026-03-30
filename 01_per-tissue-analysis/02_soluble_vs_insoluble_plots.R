# ============================================================
# Project: Multi-tissue Proteomics — Dietary Intervention
# Script:  02_soluble_vs_insoluble_plots.R
# Purpose: Per-tissue and combined scatter plots of Soluble
#          fraction FC vs. Insoluble fraction FC (Treatment vs.
#          Control). Point size encodes -log10(p-value) of
#          insolubility ratio. Generates labeled plots with
#          Pearson R² annotations for each tissue and a
#          combined multi-tissue overlay.
# Input:   Per-tissue Excel files with three sheets each:
#            Soluble (Control vs Treatment)
#            Insol (Control vs Treatment)
#            (Ratio_INSOL_SOL) Treatment vs Control
# Output:  <tissue>_soluble_vs_insoluble.pdf
#          <tissue>_processed_data.csv
#          all_tissues_soluble_vs_insoluble.pdf
#          all_tissues_combined_clean_data.csv
# Author:  Chia-Lung Chuang
# Updated: 2026-03
# NOTE:    Scripts only — raw data not included. Underlying
#          dataset is tied to a manuscript in preparation.
# ============================================================

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(readxl)

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
# Shared Plot Theme
# ============================================================

cr_theme <- theme_minimal() +
  theme(
    axis.text    = element_text(size = 14),
    axis.title   = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    plot.title   = element_text(size = 16, face = "bold")
  )

# ============================================================
# Helper: Detect and Rename Log2Fold Column
# ============================================================

# Handles variable column naming across tissues by detecting
# any column starting with "Log2Fold" and renaming it.

rename_log2fold <- function(df, new_name) {
  col <- grep("Log2Fold", names(df), value = TRUE)[1]
  if (is.na(col)) stop(paste("Log2Fold column not found. Available:", paste(names(df), collapse = ", ")))
  df %>% rename(!!new_name := all_of(col))
}

rename_pvalue <- function(df, new_name) {
  col <- grep("p-value", names(df), value = TRUE, ignore.case = TRUE)[1]
  if (is.na(col)) stop(paste("p-value column not found."))
  df %>% rename(!!new_name := all_of(col))
}

# ============================================================
# Core Function: Per-Tissue Soluble vs. Insoluble Plot
# ============================================================

# For each tissue, reads soluble, insoluble, and ratio sheets,
# merges them on protein accession + gene name, and plots:
#   Soluble FC (x) vs. Insoluble FC (y)
# Significant proteins (p_value_ratio < 0.05) colored by tissue;
# non-significant proteins shown in gray.
# R² computed separately for all proteins and significant-only.

create_tissue_plots <- function(file_path, tissue_name, tissue_color,
                                label_sol_up    =  0.6,
                                label_sol_down  = -1.0,
                                label_insol_up  =  0.4,
                                label_insol_down = -1.0) {
  
  cat("Processing:", tissue_name, "\n")
  
  # --- Load sheets ---
  sol_sheet   <- read_excel(file_path, sheet = "Soluble (Control vs Treatment)")
  insol_sheet <- read_excel(file_path, sheet = "Insol(Control vs Treatment)")
  ratio_sheet <- read_excel(file_path, sheet = "(Ratio_INSOL_SOL) Treatment vs Control")
  
  # Rename key columns robustly
  sol_sheet   <- rename_log2fold(sol_sheet,   "Log2FC_soluble")
  sol_sheet   <- rename_pvalue(sol_sheet,      "p_value_soluble")
  insol_sheet <- rename_log2fold(insol_sheet, "Log2FC_insoluble")
  insol_sheet <- rename_pvalue(insol_sheet,    "p_value_insoluble")
  ratio_sheet <- rename_log2fold(ratio_sheet, "Log2FC_ratio")
  ratio_sheet <- rename_pvalue(ratio_sheet,    "p_value_ratio")
  
  # --- Merge ---
  sol_data <- sol_sheet %>%
    select(`Protein Accession`, GN, Description,
           Log2FC_soluble, p_value_soluble)
  
  insol_data <- insol_sheet %>%
    select(`Protein Accession`, GN,
           Log2FC_insoluble, p_value_insoluble)
  
  ratio_data <- ratio_sheet %>%
    select(AC, GN, Log2FC_ratio, p_value_ratio) %>%
    rename(`Protein Accession` = AC)
  
  combined <- sol_data %>%
    inner_join(insol_data, by = c("Protein Accession", "GN")) %>%
    left_join(ratio_data,  by = c("Protein Accession", "GN")) %>%
    mutate(
      AC_GN            = paste(`Protein Accession`, GN, sep = "_"),
      p_value_soluble  = as.numeric(p_value_soluble),
      p_value_insoluble = as.numeric(p_value_insoluble),
      p_value_ratio    = as.numeric(p_value_ratio),
      significant_ratio    = p_value_ratio   <= 0.05,
      significant_soluble  = p_value_soluble <= 0.05,
      significant_insoluble = p_value_insoluble <= 0.05,
      size_ratio   = ifelse(p_value_ratio < 0.05, -log10(p_value_ratio), 0.5),
      point_alpha  = ifelse(significant_ratio, 0.4, 0.1),
      Label = case_when(
        is.na(GN) | GN == "NA" ~ NA_character_,
        significant_ratio &
          (Log2FC_soluble  >  label_sol_up   |
             Log2FC_soluble  <  label_sol_down  |
             Log2FC_insoluble > label_insol_up  |
             Log2FC_insoluble < label_insol_down) ~ GN,
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(Log2FC_soluble), !is.na(Log2FC_insoluble)) %>%
    filter(!grepl("CON_|co", `Protein Accession`, ignore.case = TRUE))
  
  # --- R² values ---
  r2_all <- round(summary(lm(Log2FC_insoluble ~ Log2FC_soluble,
                             data = combined))$r.squared, 3)
  sig_only <- combined %>% filter(significant_ratio)
  r2_sig   <- if (nrow(sig_only) > 2) {
    round(summary(lm(Log2FC_insoluble ~ Log2FC_soluble,
                     data = sig_only))$r.squared, 3)
  } else NA
  
  x_min <- min(combined$Log2FC_soluble,  na.rm = TRUE)
  y_max <- max(combined$Log2FC_insoluble, na.rm = TRUE)
  
  # --- Plot ---
  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
    geom_point(data = combined %>% filter(!significant_ratio),
               aes(x = Log2FC_soluble, y = Log2FC_insoluble,
                   size = size_ratio),
               color = "gray70", alpha = 0.1) +
    geom_point(data = sig_only,
               aes(x = Log2FC_soluble, y = Log2FC_insoluble,
                   size = size_ratio),
               color = as.character(tissue_color), alpha = 0.4) +
    geom_text_repel(
      data = function(x) subset(combined, !is.na(Label)),
      aes(x = Log2FC_soluble, y = Log2FC_insoluble, label = Label),
      size = 4, fontface = "bold",
      segment.color = "grey50",
      max.overlaps = 30, box.padding = 0.5, point.padding = 0.5,
      min.segment.length = 0.1, force = 3, na.rm = TRUE
    ) +
    annotate("text",
             x = x_min + 0.2, y = y_max - 0.2,
             label = paste0("All proteins: R² = ", r2_all),
             hjust = 0, color = "black", size = 4) +
    annotate("text",
             x = x_min + 0.2, y = y_max - 0.5,
             label = paste0("Significant proteins: R² = ",
                            ifelse(is.na(r2_sig), "N/A", r2_sig)),
             hjust = 0, color = tissue_color, size = 4) +
    scale_size_continuous(
      range  = c(0.5, 8),
      name   = "P-Value (-log10)",
      breaks = c(1.3, 2, 3, 4, 5),
      limits = c(0.5, 8),
      labels = c("1.30 (p<0.05)", "2", "3", "4", "5")
    ) +
    labs(
      title   = paste0(str_to_title(tissue_name),
                       ": Soluble vs. Insoluble FC (Treatment vs. Control)"),
      x       = "Log2 Soluble fold change (Treatment vs. Control)",
      y       = "Log2 Insoluble fold change (Treatment vs. Control)",
      caption = paste0("Point size: -log10(p-value) of insolubility ratio.\n",
                       "Gray: p >= 0.05. Colored: p < 0.05.")
    ) +
    cr_theme
  
  ggsave(paste0(tissue_name, "_soluble_vs_insoluble.pdf"),
         p, width = 12, height = 10, dpi = 600)
  
  # --- Save data ---
  write.csv(combined,
            paste0(tissue_name, "_processed_data.csv"), row.names = FALSE)
  
  clean_data <- combined %>%
    select(`Protein Accession`, GN, Description,
           Log2FC_soluble, p_value_soluble,
           Log2FC_insoluble, p_value_insoluble,
           AC_GN, Log2FC_ratio, p_value_ratio, size_ratio)
  write.csv(clean_data,
            paste0(tissue_name, "_clean_data.csv"), row.names = FALSE)
  
  cat("  Done:", tissue_name, "\n\n")
  
  return(list(
    plot_data    = combined,
    clean_data   = clean_data,
    plot         = p
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
# Combined Multi-Tissue Plot
# ============================================================

# Stacks all five tissues into one dataframe, renders
# non-significant proteins gray and significant proteins
# colored by tissue.

create_combined_plot <- function(tissue_results_list, tissue_colors) {
  tissue_names <- names(tissue_results_list)
  
  # Align column sets across tissues before binding
  all_cols <- Reduce(union, lapply(tissue_results_list,
                                   function(x) names(x$plot_data)))
  aligned <- lapply(seq_along(tissue_results_list), function(i) {
    d <- tissue_results_list[[i]]$plot_data
    missing <- setdiff(all_cols, names(d))
    for (col in missing) d[[col]] <- NA
    d[, all_cols] %>% mutate(tissue = tissue_names[i])
  })
  
  combined_all <- bind_rows(aligned)
  
  non_sig <- combined_all %>% filter(!significant_ratio)
  sig_pts  <- combined_all %>% filter( significant_ratio)
  
  r2_all <- round(summary(lm(Log2FC_insoluble ~ Log2FC_soluble,
                             data = combined_all))$r.squared, 3)
  r2_sig <- round(summary(lm(Log2FC_insoluble ~ Log2FC_soluble,
                             data = sig_pts))$r.squared, 3)
  
  x_min <- min(combined_all$Log2FC_soluble,  na.rm = TRUE)
  y_max <- max(combined_all$Log2FC_insoluble, na.rm = TRUE)
  
  p_combined <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
    geom_point(data = non_sig,
               aes(x = Log2FC_soluble, y = Log2FC_insoluble,
                   size = size_ratio),
               color = "gray70", alpha = 0.1) +
    geom_point(data = sig_pts,
               aes(x = Log2FC_soluble, y = Log2FC_insoluble,
                   color = tissue, size = size_ratio),
               alpha = 0.3) +
    annotate("text",
             x = x_min + 0.2, y = y_max - 0.2,
             label = paste0("All proteins: R² = ", r2_all),
             hjust = 0, color = "black", size = 4) +
    annotate("text",
             x = x_min + 0.2, y = y_max - 0.5,
             label = paste0("Significant proteins: R² = ", r2_sig),
             hjust = 0, color = "darkblue", size = 4) +
    scale_size_continuous(
      range  = c(0.5, 8),
      name   = "P-Value (-log10)",
      breaks = c(1.3, 2, 3, 4, 5),
      limits = c(0.5, 8),
      labels = c("1.30 (p<0.05)", "2", "3", "4", "5")
    ) +
    scale_color_manual(
      values = unlist(tissue_colors[tissue_names]),
      labels = paste("Tissue", LETTERS[seq_along(tissue_names)]),
      name   = "Tissue"
    ) +
    labs(
      title   = "All Tissues: Soluble vs. Insoluble FC (Treatment vs. Control)",
      x       = "Log2 Soluble fold change (Treatment vs. Control)",
      y       = "Log2 Insoluble fold change (Treatment vs. Control)",
      caption = "Gray: p >= 0.05. Colored: significant (p < 0.05) by tissue."
    ) +
    cr_theme
  
  write.csv(
    combined_all %>%
      select(`Protein Accession`, GN, Description,
             Log2FC_soluble, p_value_soluble,
             Log2FC_insoluble, p_value_insoluble,
             AC_GN, Log2FC_ratio, p_value_ratio,
             size_ratio, tissue),
    "all_tissues_combined_clean_data.csv",
    row.names = FALSE
  )
  
  ggsave("all_tissues_soluble_vs_insoluble.pdf",
         p_combined, width = 12, height = 10, dpi = 600)
  
  return(p_combined)
}

all_results <- list(
  tissue_a = tissue_a_results,
  tissue_b = tissue_b_results,
  tissue_c = tissue_c_results,
  tissue_d = tissue_d_results,
  tissue_e = tissue_e_results
)

combined_plot <- create_combined_plot(all_results, tissue_colors)