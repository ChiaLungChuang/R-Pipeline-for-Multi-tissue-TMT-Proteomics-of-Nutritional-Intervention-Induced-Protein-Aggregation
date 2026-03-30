# ============================================================
# Project: Multi-tissue Proteomics — Dietary Intervention
# Script:  07_fraction_comparison_plots.R
# Purpose: Per-tissue and combined scatter plots comparing
#          protein fold changes across three axes:
#            (A) Soluble FC vs. Insolubility Ratio FC
#            (B) Soluble FC vs. Insoluble FC (direct fractions)
#          Generates labeled (top 20) and clean variants with
#          Pearson R² annotations. Also exports per-tissue
#          and combined processed CSVs in a standardized format
#          for use by downstream scripts (biophysical
#          correlations, network visualization).
# Input:   <tissue>_clean_data.csv
#            (output of Script 02, soluble vs. insoluble)
# Output:  Treatment_<Tissue>_soluble_vs_ratio.pdf
#          Treatment_<Tissue>_soluble_vs_insoluble.pdf
#          all_tissues_soluble_vs_ratio.pdf
#          all_tissues_soluble_vs_insoluble.pdf
#          Treatment_<Tissue>_processed_all_Treatment.csv
#          Treatment_all_tissues_processed.csv
# Author:  Chia-Lung Chuang
# Updated: 2026-03
# NOTE:    Scripts only — raw data not included. Underlying
#          dataset is tied to a manuscript in preparation.
# ============================================================

library(tidyverse)
library(ggrepel)

# NOTE: Set your working directory to the folder containing *_clean_data.csv files
# setwd("path/to/your/data")

# ============================================================
# Tissue Color Palette
# ============================================================

tissue_colors <- c(
  tissue_a = "darkgreen",
  tissue_b = "blue4",
  tissue_c = "purple4",
  tissue_d = "darkgoldenrod",
  tissue_e = "darkred"
)

tissue_labels <- c(
  tissue_a = "Tissue A", tissue_b = "Tissue B",
  tissue_c = "Tissue C", tissue_d = "Tissue D",
  tissue_e = "Tissue E"
)

# ============================================================
# 1. Load Clean Data
# ============================================================

# Reads *_clean_data.csv produced by Script 02.
# Standardizes column names and adds a Tissue column.

load_clean_data <- function(tissue_name) {
  f <- paste0(tissue_name, "_clean_data.csv")
  if (!file.exists(f)) { warning("Not found: ", f); return(NULL) }
  d <- read_csv(f, show_col_types = FALSE)
  if (!"AC" %in% names(d) && "Protein.Accession" %in% names(d))
    d <- rename(d, AC = `Protein.Accession`)
  if (!"Gene_Name" %in% names(d) && "GN" %in% names(d))
    d <- rename(d, Gene_Name = GN)
  d %>% mutate(Tissue = tissue_name)
}

tissue_names <- c("tissue_a", "tissue_b", "tissue_c", "tissue_d", "tissue_e")
all_data <- Filter(Negate(is.null), lapply(setNames(tissue_names, tissue_names),
                                           load_clean_data))

# ============================================================
# 2. Shared Plotting Helpers
# ============================================================

r2_annotate <- function(data, x_col, y_col, color = "black", vjust = 0) {
  # Compute and add R² as annotate() call on a ggplot object
  if (nrow(data) < 3) return(NULL)
  r2  <- round(summary(lm(as.formula(paste(y_col, "~", x_col)),
                          data = data))$r.squared, 3)
  x_pos <- quantile(data[[x_col]], 0.75, na.rm = TRUE)
  y_pos <- max(data[[y_col]],   na.rm = TRUE) * (0.95 - 0.08 * vjust)
  annotate("text", x = x_pos, y = y_pos,
           label = paste0("R² = ", r2), hjust = 0,
           size = 4, color = color)
}

save_pdf <- function(p, f, w = 12, h = 10)
  ggsave(f, p, width = w, height = h, dpi = 600)

cr_theme <- theme_minimal() +
  theme(axis.text   = element_text(size = 14, face = "bold"),
        axis.title  = element_text(size = 14, face = "bold"),
        plot.title  = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"))

# ============================================================
# 3. Soluble FC vs. Insolubility Ratio FC
# ============================================================

create_soluble_vs_ratio <- function(d, tissue_name, tissue_color) {
  sig    <- d %>% filter(p_value_ratio <= 0.05) %>%
    mutate(size_pt = -log10(p_value_ratio))
  nonsig <- d %>% filter(p_value_ratio >  0.05) %>%
    mutate(size_pt = 0.5)
  
  # Top 20 by combined fold-change magnitude
  top20 <- sig %>%
    mutate(rank_score = abs(Log2FC_ratio) + abs(Log2FC_soluble)) %>%
    slice_max(rank_score, n = 20, with_ties = FALSE) %>%
    mutate(Label = Gene_Name)
  
  r2_all <- r2_annotate(d,   "Log2FC_soluble", "Log2FC_ratio", "black",    0)
  r2_sig <- r2_annotate(sig, "Log2FC_soluble", "Log2FC_ratio", tissue_color, 1)
  
  base <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
    geom_point(data = nonsig,
               aes(x = Log2FC_soluble, y = Log2FC_ratio, size = size_pt),
               color = "gray70", alpha = 0.1) +
    geom_point(data = sig,
               aes(x = Log2FC_soluble, y = Log2FC_ratio, size = size_pt),
               color = tissue_color, alpha = 0.4) +
    scale_size_continuous(range = c(0.5, 8), name = "P-Value (-log10)",
                          breaks = c(1.3, 2, 3, 4, 5), limits = c(0.5, 8),
                          labels = c("1.30 (p<0.05)", "2", "3", "4", "5")) +
    labs(
      title   = paste0(tissue_labels[tissue_name],
                       ": Soluble FC vs. Insolubility Ratio FC (Treatment vs. Control)"),
      x       = "Log2 Soluble fold change (Treatment vs. Control)",
      y       = "Log2 Insolubility Ratio FC (Treatment vs. Control)",
      caption = paste0("n significant = ", nrow(sig))
    ) + cr_theme
  
  if (!is.null(r2_all)) base <- base + r2_all
  if (!is.null(r2_sig)) base <- base + r2_sig
  
  # Labeled variant
  p_lab <- base + geom_text_repel(
    data = top20,
    aes(x = Log2FC_soluble, y = Log2FC_ratio, label = Label),
    size = 3, max.overlaps = 30, box.padding = 0.5, seed = 42
  )
  
  save_pdf(base,  paste0("Treatment_", toupper(tissue_name), "_soluble_vs_ratio.pdf"))
  save_pdf(p_lab, paste0("Treatment_", toupper(tissue_name), "_soluble_vs_ratio_labeled.pdf"))
  
  list(plot = base, sig_data = sig %>% mutate(Tissue = tissue_name),
       all_data = d %>% mutate(Tissue = tissue_name))
}

# ============================================================
# 4. Soluble FC vs. Insoluble FC (Direct Fractions)
# ============================================================

create_soluble_vs_insoluble <- function(d, tissue_name, tissue_color) {
  sig    <- d %>%
    filter(p_value_soluble <= 0.05 | p_value_insoluble <= 0.05) %>%
    mutate(size_pt = -log10(pmin(p_value_soluble, p_value_insoluble, na.rm = TRUE)))
  nonsig <- d %>%
    filter(p_value_soluble > 0.05, p_value_insoluble > 0.05) %>%
    mutate(size_pt = 0.5)
  
  top20 <- sig %>%
    mutate(rank_score = abs(Log2FC_insoluble) + abs(Log2FC_soluble)) %>%
    slice_max(rank_score, n = 20, with_ties = FALSE) %>%
    mutate(Label = Gene_Name)
  
  r2_all <- r2_annotate(d,   "Log2FC_soluble", "Log2FC_insoluble", "black",     0)
  r2_sig <- r2_annotate(sig, "Log2FC_soluble", "Log2FC_insoluble", tissue_color, 1)
  
  base <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
    geom_point(data = nonsig,
               aes(x = Log2FC_soluble, y = Log2FC_insoluble, size = size_pt),
               color = "gray70", alpha = 0.1) +
    geom_point(data = sig,
               aes(x = Log2FC_soluble, y = Log2FC_insoluble, size = size_pt),
               color = tissue_color, alpha = 0.4) +
    scale_size_continuous(range = c(0.5, 8), name = "P-Value (-log10)",
                          breaks = c(1.3, 2, 3, 4, 5), limits = c(0.5, 8),
                          labels = c("1.30 (p<0.05)", "2", "3", "4", "5")) +
    labs(
      title   = paste0(tissue_labels[tissue_name],
                       ": Soluble FC vs. Insoluble FC (Treatment vs. Control)"),
      x       = "Log2 Soluble fold change (Treatment vs. Control)",
      y       = "Log2 Insoluble fold change (Treatment vs. Control)",
      caption = paste0("n significant = ", nrow(sig))
    ) + cr_theme
  
  if (!is.null(r2_all)) base <- base + r2_all
  if (!is.null(r2_sig)) base <- base + r2_sig
  
  p_lab <- base + geom_text_repel(
    data = top20,
    aes(x = Log2FC_soluble, y = Log2FC_insoluble, label = Label),
    size = 3, max.overlaps = 30, box.padding = 0.5, seed = 42
  )
  
  save_pdf(base,  paste0("Treatment_", toupper(tissue_name), "_soluble_vs_insoluble.pdf"))
  save_pdf(p_lab, paste0("Treatment_", toupper(tissue_name), "_soluble_vs_insoluble_labeled.pdf"))
  
  list(plot = base, sig_data = sig %>% mutate(Tissue = tissue_name),
       all_data = d %>% mutate(Tissue = tissue_name))
}

# ============================================================
# 5. Run Per-Tissue Analysis and Export Processed CSVs
# ============================================================

ratio_results <- list()
insol_results <- list()

for (tn in names(all_data)) {
  d     <- all_data[[tn]]
  color <- tissue_colors[tn]
  
  ratio_results[[tn]] <- create_soluble_vs_ratio(d, tn, color)
  insol_results[[tn]] <- create_soluble_vs_insoluble(d, tn, color)
  
  # Export standardized processed CSV for use by Scripts 04/06
  export_d <- d %>%
    mutate(
      size_ratio = ifelse(p_value_ratio <= 0.05, -log10(p_value_ratio), 0.5)
    )
  write.csv(export_d,
            paste0("Treatment_", toupper(tn), "_processed_all_Treatment.csv"),
            row.names = FALSE)
  cat("Exported processed CSV:", tn, "\n")
}

# ============================================================
# 6. Combined Multi-Tissue Plots
# ============================================================

build_combined_plot <- function(results_list, x_col, y_col,
                                x_label, y_label, out_file) {
  sig_all <- bind_rows(lapply(results_list, `[[`, "sig_data"))
  all_all <- bind_rows(lapply(results_list, `[[`, "all_data"))
  
  nonsig_all <- all_all %>%
    anti_join(sig_all %>% select(AC, Tissue), by = c("AC", "Tissue")) %>%
    mutate(size_pt = 0.5)
  sig_all <- sig_all %>%
    mutate(size_pt = -log10(pmin(p_value_ratio, 0.05, na.rm = TRUE)))
  
  r2_all <- if (nrow(all_all) > 2) round(
    summary(lm(as.formula(paste(y_col, "~", x_col)),
               data = all_all))$r.squared, 3) else NA
  r2_sig <- if (nrow(sig_all) > 2) round(
    summary(lm(as.formula(paste(y_col, "~", x_col)),
               data = sig_all))$r.squared, 3) else NA
  
  x_min <- min(all_all[[x_col]], na.rm = TRUE)
  y_max <- max(all_all[[y_col]], na.rm = TRUE)
  
  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
    geom_point(data = nonsig_all,
               aes(x = .data[[x_col]], y = .data[[y_col]], size = size_pt),
               color = "gray70", alpha = 0.1) +
    geom_point(data = sig_all,
               aes(x = .data[[x_col]], y = .data[[y_col]],
                   color = Tissue, size = size_pt),
               alpha = 0.4) +
    scale_color_manual(values = tissue_colors, labels = tissue_labels,
                       name = "Tissue") +
    scale_size_continuous(range = c(0.5, 8), name = "P-Value (-log10)",
                          breaks = c(1.3, 2, 3, 4, 5), limits = c(0.5, 8)) +
    annotate("text", x = x_min + 0.2, y = y_max - 0.2,
             label = paste0("All: R² = ", r2_all),
             hjust = 0, size = 4, color = "black") +
    annotate("text", x = x_min + 0.2, y = y_max - 0.5,
             label = paste0("Significant: R² = ", r2_sig),
             hjust = 0, size = 4, color = "darkblue") +
    labs(title   = paste0("All Tissues: ", x_label, " vs. ", y_label),
         x = x_label, y = y_label) +
    cr_theme
  
  save_pdf(p, out_file, w = 14, h = 12)
  cat("Saved:", out_file, "\n")
  p
}

build_combined_plot(
  ratio_results,
  "Log2FC_soluble", "Log2FC_ratio",
  "Log2 Soluble FC (Treatment vs. Control)",
  "Log2 Insolubility Ratio FC (Treatment vs. Control)",
  "all_tissues_soluble_vs_ratio.pdf"
)

build_combined_plot(
  insol_results,
  "Log2FC_soluble", "Log2FC_insoluble",
  "Log2 Soluble FC (Treatment vs. Control)",
  "Log2 Insoluble FC (Treatment vs. Control)",
  "all_tissues_soluble_vs_insoluble.pdf"
)

# ============================================================
# 7. Export Combined Processed CSV
# ============================================================

bind_rows(lapply(names(all_data), function(tn) {
  all_data[[tn]] %>%
    mutate(Tissue = tn,
           size_ratio = ifelse(p_value_ratio <= 0.05,
                               -log10(p_value_ratio), 0.5))
})) %>%
  write.csv("Treatment_all_tissues_processed.csv", row.names = FALSE)

cat("Saved: Treatment_all_tissues_processed.csv\n")

# Summary table
summary_tbl <- bind_rows(lapply(names(ratio_results), function(tn) {
  d <- ratio_results[[tn]]
  data.frame(
    Tissue              = tn,
    Total_Proteins      = nrow(d$all_data),
    Significant_Ratio   = nrow(d$sig_data),
    Significant_Soluble = sum(all_data[[tn]]$p_value_soluble <= 0.05, na.rm = TRUE),
    Significant_Insol   = sum(all_data[[tn]]$p_value_insoluble <= 0.05, na.rm = TRUE)
  )
}))
write.csv(summary_tbl, "Treatment_fraction_comparison_summary.csv", row.names = FALSE)
cat("Summary:\n"); print(summary_tbl)