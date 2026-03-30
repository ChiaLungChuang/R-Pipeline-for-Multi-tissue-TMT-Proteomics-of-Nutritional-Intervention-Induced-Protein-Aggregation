# ============================================================
# Project: Multi-tissue Proteomics — Dietary Intervention
# Script:  06_biophysical_correlations.R
# Purpose: Correlates protein insolubility fold change
#          (Treatment vs. Control) against two biophysical
#          properties: protein half-life (days) and thermal
#          stability (melting temperature, °C). Both properties
#          are joined to the per-tissue processed data via
#          protein accession ID matching. Generates per-tissue
#          scatter plots (labeled and R²-annotated variants)
#          and a multi-tissue overlay for half-life.
# Input:   Treatment_<Tissue>_processed_all_Treatment.csv
#            (output of Script 01)
#          protein_half-life.xlsx  — columns: Proteins (ID),
#            per-tissue half-life columns
#          meltTemperature.xlsx   — columns: Protein_ID,
#            gene_name, meltPoint
# Output:  Treatment_<Tissue>_halflife_*.pdf/csv
#          Treatment_multi_tissue_halflife.pdf
#          Treatment_halflife_correlations.csv
#          Treatment_<Tissue>_melttemp_*.pdf/csv
# Author:  Chia-Lung Chuang
# Updated: 2026-03
# NOTE:    Scripts only — raw data not included. Underlying
#          dataset is tied to a manuscript in preparation.
#          Half-life data available for Tissues A, B, C only.
#          Melting temperature data available for Tissue B only.
# ============================================================

library(tidyverse)
library(ggrepel)
library(readxl)
library(ggpmisc)

# NOTE: Set your working directory to the folder containing all input files
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

# ============================================================
# 1. Load Per-Tissue Processed Data
# ============================================================

load_tissue <- function(tissue_name) {
  f <- paste0("Treatment_", toupper(tissue_name), "_processed_all_Treatment.csv")
  if (!file.exists(f)) { warning("Not found: ", f); return(NULL) }
  d <- read_csv(f, show_col_types = FALSE)
  if (!"AC" %in% names(d) && "Protein.Accession" %in% names(d))
    d <- rename(d, AC = `Protein.Accession`)
  if (!"Gene_Name" %in% names(d) && "GN" %in% names(d))
    d <- rename(d, Gene_Name = GN)
  # Extract UniProt accession from pipe-delimited format (e.g. sp|P12345|GENE_MOUSE)
  d %>% mutate(ID = str_extract(AC, "(?<=\\|)[A-Z0-9]+(\\-\\d+)?(?=\\|)"))
}

tissue_data <- list(
  tissue_a = load_tissue("tissuea"),
  tissue_b = load_tissue("tissueb"),
  tissue_c = load_tissue("tissuec"),
  tissue_d = load_tissue("tissued"),
  tissue_e = load_tissue("tissuee")
)
tissue_data <- Filter(Negate(is.null), tissue_data)

# ============================================================
# 2. Shared Plot Builder
# ============================================================

# Core scatter plot: biophysical property (x) vs. insolubility FC (y).
# Significant points (p_ratio ≤ 0.05) sized by -log10(p), colored by
# tissue; non-significant points shown as small gray dots.

build_biophys_plot <- function(plot_data, x_col, tissue_name, tissue_color,
                               x_label, x_lim = NULL,
                               label_col = NULL,
                               title_suffix = "") {
  sig     <- plot_data %>% filter(p_value_ratio <= 0.05) %>%
    mutate(size_pt = -log10(p_value_ratio), pt_alpha = 0.7)
  nonsig  <- plot_data %>% filter(p_value_ratio >  0.05) %>%
    mutate(size_pt = 0.5,                   pt_alpha = 0.3)
  
  r2_all <- if (nrow(plot_data) > 2)
    round(summary(lm(as.formula(paste("Log2FC_ratio ~", x_col)),
                     data = plot_data))$r.squared, 3) else NA
  
  p <- ggplot() +
    geom_point(data = nonsig,
               aes(x = .data[[x_col]], y = Log2FC_ratio, size = size_pt),
               color = "gray70", alpha = 0.15) +
    geom_point(data = sig,
               aes(x = .data[[x_col]], y = Log2FC_ratio, size = size_pt),
               color = tissue_color, alpha = 0.7) +
    scale_size_continuous(
      range  = c(0.4, 9), name = "P-Value (-log10)",
      breaks = c(1.3, 2, 3, 4, 5), limits = c(0.5, 8),
      labels = c("1.30 (p<0.05)", "2", "3", "4", "5")
    ) +
    scale_x_continuous(labels = scales::label_number(accuracy = 0.1),
                       limits = x_lim) +
    labs(
      title   = paste0(str_to_title(tissue_name), ": ", x_label,
                       " vs. Insolubility FC", title_suffix),
      x       = x_label,
      y       = "Log2 Insolubility Ratio FC (Treatment vs. Control)",
      caption = paste0("Point size = -log10(p-value). Small points: p > 0.05.\n",
                       "n significant = ", nrow(sig))
    ) +
    theme_minimal() +
    theme(axis.text   = element_text(size = 14, face = "bold"),
          axis.title  = element_text(size = 14, face = "bold"),
          plot.title  = element_text(size = 16, face = "bold"),
          legend.text = element_text(size = 12))
  
  if (!is.na(r2_all)) {
    x_pos <- if (!is.null(x_lim)) x_lim[2] * 0.7 else
      quantile(plot_data[[x_col]], 0.85, na.rm = TRUE)
    y_pos <- max(plot_data$Log2FC_ratio, na.rm = TRUE) * 0.9
    p <- p + annotate("text", x = x_pos, y = y_pos,
                      label = paste0("R² = ", r2_all),
                      hjust = 0, size = 4.5, fontface = "bold")
  }
  
  if (!is.null(label_col) && label_col %in% names(plot_data)) {
    label_pts <- plot_data %>% filter(!is.na(.data[[label_col]]))
    if (nrow(label_pts) > 0)
      p <- p + geom_text_repel(
        data = label_pts,
        aes(x = .data[[x_col]], y = Log2FC_ratio,
            label = .data[[label_col]]),
        size = 3, segment.color = "grey50",
        max.overlaps = Inf, box.padding = 0.5,
        point.padding = 0.5, na.rm = TRUE
      )
  }
  p
}

save_pdf <- function(p, file, w = 7200, h = 6400)
  ggsave(file, p, width = w, height = h, units = "px", dpi = 600)

# ============================================================
# 3. Half-Life Analysis
# ============================================================

# Half-life data available for Tissues A, B, C (update column
# names below to match the actual column headers in your file).
#
# Mapping: tissue_name → column name in protein_half-life.xlsx
halflife_col_map <- c(
  tissue_a = "protein.half.life..tissue_a.",
  tissue_b = "protein.half.life..tissue_b.",
  tissue_c = "protein.half.life..tissue_c."
)

run_halflife_analysis <- function() {
  hl_raw <- read_excel("protein_half-life.xlsx")
  hl_raw <- hl_raw %>%
    mutate(ID = str_extract(Proteins, "^[A-Z0-9]+(?=[-_]|$)"))
  
  results <- list()
  
  for (tissue_name in names(halflife_col_map)) {
    if (!tissue_name %in% names(tissue_data)) next
    hl_col <- halflife_col_map[tissue_name]
    if (!hl_col %in% names(hl_raw)) {
      cat("Half-life column not found:", hl_col, "— skipping\n"); next
    }
    
    d <- tissue_data[[tissue_name]] %>%
      left_join(hl_raw %>% select(ID, Protein_Half_Life = all_of(hl_col)),
                by = "ID") %>%
      filter(!is.na(Protein_Half_Life), !is.na(Log2FC_ratio)) %>%
      mutate(
        Label_hl = case_when(
          p_value_ratio <= 0.05 & Protein_Half_Life > 5  ~ Gene_Name,
          p_value_ratio <= 0.05 & Protein_Half_Life < 0  ~ Gene_Name,
          TRUE ~ NA_character_
        )
      )
    
    write.csv(d, paste0("Treatment_", toupper(tissue_name),
                        "_halflife.csv"), row.names = FALSE)
    cat(tissue_name, "— half-life matched:", nrow(d), "proteins\n")
    
    # No-label + R²
    p_base <- build_biophys_plot(d, "Protein_Half_Life", tissue_name,
                                 tissue_colors[tissue_name],
                                 "Protein Half-Life (days)",
                                 title_suffix = " (Treatment vs. Control)")
    p_r2   <- p_base + stat_poly_eq(
      aes(x = Protein_Half_Life, y = Log2FC_ratio,
          label = paste(after_stat(rr.label))),
      formula = y ~ x, rr.digits = 3, parse = TRUE,
      color = "black", label.y = 0.1, label.x = 0.8
    )
    save_pdf(p_base, paste0("Treatment_", toupper(tissue_name), "_halflife_nolabel.pdf"))
    save_pdf(p_r2,   paste0("Treatment_", toupper(tissue_name), "_halflife_nolabel_R.pdf"))
    
    # Labeled + R²
    p_lab  <- build_biophys_plot(d, "Protein_Half_Life", tissue_name,
                                 tissue_colors[tissue_name],
                                 "Protein Half-Life (days)",
                                 label_col = "Label_hl",
                                 title_suffix = " (Treatment vs. Control)")
    p_labR <- p_lab + stat_poly_eq(
      aes(x = Protein_Half_Life, y = Log2FC_ratio,
          label = paste(after_stat(rr.label))),
      formula = y ~ x, rr.digits = 3, parse = TRUE,
      color = "black", label.y = 0.1, label.x = 0.8
    )
    save_pdf(p_lab,  paste0("Treatment_", toupper(tissue_name), "_halflife_labeled.pdf"))
    save_pdf(p_labR, paste0("Treatment_", toupper(tissue_name), "_halflife_labeled_R.pdf"))
    
    results[[tissue_name]] <- d
  }
  
  # Multi-tissue overlay
  if (length(results) >= 2) {
    combined <- bind_rows(lapply(names(results), function(tn)
      results[[tn]] %>% mutate(Tissue = tn))) %>%
      mutate(size_pt  = ifelse(p_value_ratio <= 0.05,
                               -log10(p_value_ratio), 0.5),
             pt_alpha = ifelse(p_value_ratio <= 0.05, 0.5, 0.15))
    
    r2_combined <- round(summary(lm(Log2FC_ratio ~ Protein_Half_Life,
                                    data = combined))$r.squared, 3)
    
    p_multi <- ggplot(combined,
                      aes(x = Protein_Half_Life, y = Log2FC_ratio,
                          color = Tissue, size = size_pt, alpha = pt_alpha)) +
      geom_point() +
      scale_alpha_identity() +
      scale_color_manual(values = tissue_colors,
                         labels = c(tissue_a = "Tissue A", tissue_b = "Tissue B",
                                    tissue_c = "Tissue C")) +
      scale_size_continuous(range = c(0.4, 7), name = "P-Value (-log10)",
                            breaks = c(1.3, 2, 3, 4, 5), limits = c(0.5, 8)) +
      annotate("text",
               x = quantile(combined$Protein_Half_Life, 0.8, na.rm = TRUE),
               y = max(combined$Log2FC_ratio, na.rm = TRUE) * 0.9,
               label = paste0("R² = ", r2_combined),
               size = 5, fontface = "bold") +
      labs(title = "Multi-tissue: Protein Half-Life vs. Insolubility FC",
           x = "Protein Half-Life (days)",
           y = "Log2 Insolubility Ratio FC (Treatment vs. Control)") +
      theme_minimal() +
      theme(axis.text  = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 14, face = "bold"),
            plot.title = element_text(size = 16, face = "bold"))
    
    save_pdf(p_multi, "Treatment_multi_tissue_halflife.pdf", w = 8000, h = 6400)
    
    # Correlation summary
    cor_summary <- bind_rows(lapply(names(results), function(tn) {
      d <- results[[tn]]
      ct <- cor.test(d$Protein_Half_Life, d$Log2FC_ratio)
      data.frame(Tissue = tn, r = round(ct$estimate, 4),
                 p_value = ct$p.value, n = nrow(d))
    }))
    write.csv(cor_summary, "Treatment_halflife_correlations.csv", row.names = FALSE)
    cat("\nHalf-life correlations:\n"); print(cor_summary)
  }
}

run_halflife_analysis()

# ============================================================
# 4. Melting Temperature Analysis
# ============================================================

# Melting temperature data currently available for Tissue B only.
# Add additional tissues to tissue_melt_tissues if data is available.

tissue_melt_tissues <- c("tissue_b")

run_melttemp_analysis <- function() {
  mt_raw <- read_excel("meltTemperature.xlsx") %>%
    mutate(ID = str_extract(Protein_ID, "^[A-Z0-9]+(?=_)"))
  
  for (tissue_name in tissue_melt_tissues) {
    if (!tissue_name %in% names(tissue_data)) next
    
    d <- tissue_data[[tissue_name]] %>%
      left_join(mt_raw %>% select(ID, Melt_Temperature = meltPoint),
                by = "ID") %>%
      filter(!is.na(Melt_Temperature), !is.na(Log2FC_ratio)) %>%
      mutate(
        Label_mt = case_when(
          p_value_ratio <= 0.05 & Melt_Temperature > 60 ~ Gene_Name,
          p_value_ratio <= 0.05 & Melt_Temperature < 42 ~ Gene_Name,
          TRUE ~ NA_character_
        )
      )
    
    write.csv(d, paste0("Treatment_", toupper(tissue_name),
                        "_melttemp.csv"), row.names = FALSE)
    cat(tissue_name, "— melting temp matched:", nrow(d), "proteins\n")
    
    # No-label
    p_base <- build_biophys_plot(d, "Melt_Temperature", tissue_name,
                                 tissue_colors[tissue_name],
                                 "Protein Melting Temperature (°C)",
                                 title_suffix = " (Treatment vs. Control)")
    p_r2   <- p_base + stat_poly_eq(
      aes(x = Melt_Temperature, y = Log2FC_ratio,
          label = paste(after_stat(rr.label))),
      formula = y ~ x, rr.digits = 3, parse = TRUE,
      color = "black", label.y = 0.1, label.x = 0.8
    )
    save_pdf(p_base, paste0("Treatment_", toupper(tissue_name), "_melttemp_nolabel.pdf"))
    save_pdf(p_r2,   paste0("Treatment_", toupper(tissue_name), "_melttemp_nolabel_R.pdf"))
    
    # Labeled
    p_lab  <- build_biophys_plot(d, "Melt_Temperature", tissue_name,
                                 tissue_colors[tissue_name],
                                 "Protein Melting Temperature (°C)",
                                 label_col = "Label_mt",
                                 title_suffix = " (Treatment vs. Control)")
    p_labR <- p_lab + stat_poly_eq(
      aes(x = Melt_Temperature, y = Log2FC_ratio,
          label = paste(after_stat(rr.label))),
      formula = y ~ x, rr.digits = 3, parse = TRUE,
      color = "black", label.y = 0.1, label.x = 0.8
    )
    save_pdf(p_lab,  paste0("Treatment_", toupper(tissue_name), "_melttemp_labeled.pdf"))
    save_pdf(p_labR, paste0("Treatment_", toupper(tissue_name), "_melttemp_labeled_R.pdf"))
  }
}

run_melttemp_analysis()