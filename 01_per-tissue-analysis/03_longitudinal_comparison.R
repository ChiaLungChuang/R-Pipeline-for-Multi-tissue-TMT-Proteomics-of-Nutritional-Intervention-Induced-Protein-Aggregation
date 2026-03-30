# ============================================================
# Project: Multi-tissue Proteomics — Dietary Intervention
# Script:  03_longitudinal_comparison.R
# Purpose: Cross-study integration comparing protein
#          insolubility changes from a dietary intervention
#          (Treatment vs. Control) against longitudinal
#          aging changes (Timepoint2 vs. Timepoint1) across
#          matched tissues. Classifies proteins as concordant
#          (same direction) or discordant (reversed) between
#          the two studies. Generates scatter plots with
#          directional color-coding and a concordance summary.
# Input:   longitudinal_<Tissue>.csv — processed outputs from
#            the companion longitudinal aging project
#          Treatment_<Tissue>_processed_all_Treatment.csv —
#            processed outputs from Script 01
# Output:  <Tissue>_vs_<Tissue>_insolubility_Treatment.pdf
#          <Tissue>_vs_<Tissue>_insolubility_data_Treatment.csv
#          Treatment_Longitudinal_Concordance_Summary.csv
#          shared_significant_proteins_Treatment.csv
# Author:  Chia-Lung Chuang
# Updated: 2026-03
# NOTE:    Scripts only — raw data not included. Underlying
#          datasets are tied to manuscripts in preparation.
#          Longitudinal aging CSVs are outputs of a separate
#          companion analysis pipeline.
# ============================================================

library(tidyverse)
library(ggrepel)
library(ggpmisc)

# NOTE: Set your working directory to the folder containing all input CSVs
# setwd("path/to/your/data")

# ============================================================
# Study Design Mappings
# ============================================================

# Maps each Treatment tissue to its longitudinal counterpart(s).
# Some Treatment tissues correspond to multiple sub-regions
# in the longitudinal dataset.

tissue_mappings <- list(
  "Tissue_A" = c("Tissue_A1", "Tissue_A2", "Tissue_A3"),
  "Tissue_B" = c("Tissue_B"),
  "Tissue_C" = c("Tissue_C"),
  "Tissue_D" = c("Tissue_D"),
  "Tissue_E" = c("Tissue_E1", "Tissue_E2")
)

# File mappings: longitudinal aging data (Timepoint2 vs Timepoint1)
longitudinal_file_mappings <- c(
  "Tissue_A1" = "longitudinal_TissueA1.csv",
  "Tissue_A2" = "longitudinal_TissueA2.csv",
  "Tissue_A3" = "longitudinal_TissueA3.csv",
  "Tissue_B"  = "longitudinal_TissueB.csv",
  "Tissue_C"  = "longitudinal_TissueC.csv",
  "Tissue_D"  = "longitudinal_TissueD.csv",
  "Tissue_E1" = "longitudinal_TissueE1.csv",
  "Tissue_E2" = "longitudinal_TissueE2.csv"
)

# File mappings: dietary intervention processed data
treatment_file_mappings <- c(
  "Tissue_A" = "Treatment_TISSUEA_processed_all_Treatment.csv",
  "Tissue_B" = "Treatment_TISSUEB_processed_all_Treatment.csv",
  "Tissue_C" = "Treatment_TISSUEC_processed_all_Treatment.csv",
  "Tissue_D" = "Treatment_TISSUED_processed_all_Treatment.csv",
  "Tissue_E" = "Treatment_TISSUEE_processed_all_Treatment.csv"
)

# ============================================================
# Data Loaders
# ============================================================

# Load longitudinal aging data and extract insolubility FC
# columns for Timepoint2 vs. Timepoint1 comparison.
load_longitudinal_data <- function(tissue_label, data_dir = ".") {
  file_path <- file.path(data_dir, longitudinal_file_mappings[tissue_label])
  if (!file.exists(file_path)) {
    warning("File not found: ", file_path)
    return(data.frame())
  }
  
  d <- read.csv(file_path)
  
  # Detect Timepoint2 vs. Timepoint1 ratio and abundance columns
  ratio_col   <- grep("Ratio_FC.*timepoint2.*timepoint1",
                      colnames(d), value = TRUE, ignore.case = TRUE)[1]
  abund_col   <- grep("Abundance_FC.*timepoint2.*timepoint1",
                      colnames(d), value = TRUE, ignore.case = TRUE)[1]
  ratio_p_col <- grep("ra\\.p\\.value.*timepoint2.*timepoint1",
                      colnames(d), value = TRUE, ignore.case = TRUE)[1]
  abund_p_col <- grep("ab\\.p\\.value.*timepoint2.*timepoint1|total\\.p\\.value.*timepoint2.*timepoint1",
                      colnames(d), value = TRUE, ignore.case = TRUE)[1]
  
  if (is.na(ratio_col) || is.na(ratio_p_col)) {
    warning("Required columns not found for ", tissue_label)
    return(data.frame())
  }
  
  if (is.na(abund_p_col)) abund_p_col <- ratio_p_col
  
  d %>%
    mutate(
      Ratio_FC_Longitudinal     = .data[[ratio_col]],
      Abundance_FC_Longitudinal = if (!is.na(abund_col)) .data[[abund_col]] else NA_real_,
      p_value_Ratio_Longitudinal    = .data[[ratio_p_col]],
      p_value_Abundance_Longitudinal = .data[[abund_p_col]],
      Tissue_Label = tissue_label
    ) %>%
    select(GN_ID, Gene_Name,
           Ratio_FC_Longitudinal, Abundance_FC_Longitudinal,
           p_value_Ratio_Longitudinal, p_value_Abundance_Longitudinal,
           Tissue_Label)
}

# Load dietary intervention processed data (from Script 01)
load_treatment_data <- function(tissue_label, data_dir = ".") {
  file_path <- file.path(data_dir, treatment_file_mappings[tissue_label])
  if (!file.exists(file_path)) {
    warning("File not found: ", file_path)
    return(data.frame())
  }
  
  d <- read.csv(file_path)
  
  d %>%
    mutate(GN_ID = paste(GN, AC, sep = "_")) %>%
    select(GN_ID, Gene_Name = GN,
           Ratio_FC_Treatment    = Log2FC_ratio,
           p_value_Ratio_Treatment = p_value_ratio,
           Tissue_Label = tissue_label)
}

# ============================================================
# Comparison Plot Function
# ============================================================

# Plots longitudinal insolubility FC (x) vs. treatment
# insolubility FC (y) for proteins significant in both
# datasets. Color encodes concordance pattern.

create_comparison_plot <- function(longitudinal_data, treatment_data,
                                   p_threshold = 0.05,
                                   output_file = NULL,
                                   longitudinal_label = "Longitudinal",
                                   treatment_label = "Treatment") {
  merged <- longitudinal_data %>%
    inner_join(treatment_data, by = "GN_ID") %>%
    mutate(Gene_Name = coalesce(Gene_Name.x, Gene_Name.y)) %>%
    select(-Gene_Name.x, -Gene_Name.y)
  
  sig <- merged %>%
    filter(p_value_Ratio_Longitudinal <= p_threshold,
           p_value_Ratio_Treatment    <= p_threshold) %>%
    mutate(
      consistency = case_when(
        Ratio_FC_Longitudinal > 0 & Ratio_FC_Treatment < 0 ~ "Longitudinal Up, Treatment Down",
        Ratio_FC_Longitudinal < 0 & Ratio_FC_Treatment > 0 ~ "Longitudinal Down, Treatment Up",
        Ratio_FC_Longitudinal > 0 & Ratio_FC_Treatment > 0 ~ "Both Increased",
        Ratio_FC_Longitudinal < 0 & Ratio_FC_Treatment < 0 ~ "Both Decreased",
        TRUE ~ "Other"
      )
    )
  
  if (nrow(sig) == 0) {
    warning("No proteins significant in both datasets for: ",
            longitudinal_label, " vs. ", treatment_label)
    return(NULL)
  }
  
  # Label top 10 by combined fold-change magnitude per category
  label_data <- sig %>%
    filter(!is.na(Gene_Name), consistency != "Other") %>%
    mutate(magnitude = abs(Ratio_FC_Longitudinal) + abs(Ratio_FC_Treatment)) %>%
    group_by(consistency) %>%
    slice_max(order_by = magnitude, n = 10, with_ties = FALSE) %>%
    ungroup()
  
  cor_result <- cor.test(sig$Ratio_FC_Longitudinal, sig$Ratio_FC_Treatment)
  
  p <- ggplot(sig, aes(x = Ratio_FC_Longitudinal, y = Ratio_FC_Treatment)) +
    geom_point(aes(color    = consistency,
                   size     = -log10(pmin(p_value_Ratio_Longitudinal,
                                          p_value_Ratio_Treatment))),
               alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    stat_poly_eq(
      formula = y ~ x,
      aes(label = paste(after_stat(rr.label))),
      parse = TRUE, color = "black",
      label.y = 1, label.x = 0.5, size = 4
    ) +
    geom_text_repel(
      data = label_data,
      aes(label = Gene_Name),
      size = 3, max.overlaps = Inf, box.padding = 0.4, seed = 42
    ) +
    scale_color_manual(values = c(
      "Both Increased"                 = "#4DAF4A",
      "Both Decreased"                 = "#377EB8",
      "Longitudinal Up, Treatment Down" = "#E41A1C",
      "Longitudinal Down, Treatment Up" = "#FF7F00",
      "Other"                          = "gray50"
    )) +
    scale_size_continuous(range = c(1, 6)) +
    labs(
      title   = paste0("Protein Insolubility: ",
                       longitudinal_label, " vs. ", treatment_label),
      x       = "Log2 Insolubility FC (Timepoint2 / Timepoint1)",
      y       = "Log2 Insolubility FC (Treatment / Control)",
      color   = "Change Pattern",
      size    = "-log10(p-value)",
      caption = paste0("n = ", nrow(sig),
                       " proteins significant (p < ", p_threshold, ") in both studies")
    ) +
    theme_minimal() +
    theme(
      plot.title   = element_text(face = "bold", size = 14),
      axis.title   = element_text(size = 12),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
  
  if (!is.null(output_file))
    ggsave(output_file, p, width = 12, height = 8)
  
  list(
    plot        = p,
    statistics  = sig %>% count(consistency) %>% rename(n = n),
    correlation = cor_result,
    data        = sig
  )
}

# ============================================================
# Run All Tissue Comparisons
# ============================================================

run_all_comparisons <- function(data_dir = ".",
                                output_dir = "longitudinal_vs_treatment_results") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  all_results   <- list()
  summary_stats <- data.frame()
  
  get_count <- function(stats_df, category) {
    if (is.null(stats_df) || !all(c("consistency", "n") %in% colnames(stats_df)))
      return(0)
    idx <- which(stats_df$consistency == category)
    if (length(idx)) stats_df$n[idx] else 0
  }
  
  for (treatment_tissue in names(tissue_mappings)) {
    longitudinal_tissues <- tissue_mappings[[treatment_tissue]]
    
    treatment_data <- load_treatment_data(treatment_tissue, data_dir)
    if (nrow(treatment_data) == 0) next
    
    for (long_tissue in longitudinal_tissues) {
      cat("Comparing:", long_tissue, "vs", treatment_tissue, "\n")
      
      long_data <- load_longitudinal_data(long_tissue, data_dir)
      if (nrow(long_data) == 0) next
      
      comp_name   <- paste(long_tissue, treatment_tissue, sep = "_vs_")
      output_file <- file.path(output_dir,
                               paste0(comp_name, "_insolubility_Treatment.pdf"))
      
      result <- create_comparison_plot(
        long_data, treatment_data,
        output_file          = output_file,
        longitudinal_label   = long_tissue,
        treatment_label      = treatment_tissue
      )
      
      if (!is.null(result$data)) {
        write.csv(result$data,
                  file.path(output_dir, paste0(comp_name, "_data.csv")),
                  row.names = FALSE)
      }
      
      all_results[[comp_name]] <- result
      
      summary_stats <- bind_rows(summary_stats, data.frame(
        Comparison      = comp_name,
        Total_Proteins  = if (!is.null(result$statistics)) sum(result$statistics$n) else NA,
        Both_Increased  = get_count(result$statistics, "Both Increased"),
        Both_Decreased  = get_count(result$statistics, "Both Decreased"),
        Long_Up_Trt_Down = get_count(result$statistics, "Longitudinal Up, Treatment Down"),
        Long_Down_Trt_Up = get_count(result$statistics, "Longitudinal Down, Treatment Up"),
        Other           = get_count(result$statistics, "Other")
      ))
    }
  }
  
  write.csv(summary_stats,
            file.path(output_dir, "Treatment_Longitudinal_Summary.csv"),
            row.names = FALSE)
  
  list(results = all_results, summary = summary_stats)
}

# ============================================================
# Concordance Analysis
# ============================================================

# For proteins significant in BOTH studies, calculates the
# % that change in the same direction (concordant) vs.
# opposite direction (discordant).

calculate_concordance <- function(longitudinal_data, treatment_data,
                                  p_threshold = 0.05) {
  merged <- longitudinal_data %>%
    inner_join(treatment_data %>%
                 select(Gene_Name, Ratio_FC_Treatment, p_value_Ratio_Treatment),
               by = "Gene_Name", relationship = "many-to-many")
  
  both_sig <- merged %>%
    filter(p_value_Ratio_Longitudinal < p_threshold,
           p_value_Ratio_Treatment    < p_threshold)
  
  if (nrow(both_sig) == 0) return(NULL)
  
  cat_data <- both_sig %>%
    mutate(
      Pattern = case_when(
        Ratio_FC_Longitudinal > 0 & Ratio_FC_Treatment > 0 ~ "Concordant_BothUp",
        Ratio_FC_Longitudinal < 0 & Ratio_FC_Treatment < 0 ~ "Concordant_BothDown",
        Ratio_FC_Longitudinal > 0 & Ratio_FC_Treatment < 0 ~ "Discordant_LongUp_TrtDown",
        Ratio_FC_Longitudinal < 0 & Ratio_FC_Treatment > 0 ~ "Discordant_LongDown_TrtUp",
        TRUE ~ "Other"
      ),
      Concordance = case_when(
        grepl("Concordant", Pattern)  ~ "Concordant",
        grepl("Discordant", Pattern)  ~ "Discordant",
        TRUE ~ "Other"
      )
    )
  
  n_total       <- nrow(cat_data)
  pattern_counts <- table(cat_data$Pattern)
  
  list(
    n_both_sig     = n_total,
    n_concordant   = sum(cat_data$Concordance == "Concordant"),
    n_discordant   = sum(cat_data$Concordance == "Discordant"),
    pct_concordant = 100 * sum(cat_data$Concordance == "Concordant") / n_total,
    pct_discordant = 100 * sum(cat_data$Concordance == "Discordant") / n_total,
    pattern_counts = pattern_counts
  )
}

run_concordance_analysis <- function(data_dir = ".",
                                     output_dir = "longitudinal_vs_treatment_results",
                                     p_threshold = 0.05) {
  concordance_summary <- data.frame()
  
  for (treatment_tissue in names(tissue_mappings)) {
    long_tissues   <- tissue_mappings[[treatment_tissue]]
    treatment_data <- load_treatment_data(treatment_tissue, data_dir)
    if (nrow(treatment_data) == 0) next
    
    for (long_tissue in long_tissues) {
      cat("Concordance:", long_tissue, "vs", treatment_tissue, "\n")
      
      long_data <- load_longitudinal_data(long_tissue, data_dir)
      if (nrow(long_data) == 0) next
      
      conc <- calculate_concordance(long_data, treatment_data, p_threshold)
      
      if (is.null(conc)) {
        cat("  No proteins significant in both datasets\n\n")
        conc <- list(n_both_sig = 0, n_concordant = 0, n_discordant = 0,
                     pct_concordant = NA, pct_discordant = NA,
                     pattern_counts = integer(0))
      } else {
        cat(sprintf("  Both significant: %d | Concordant: %.1f%% | Discordant: %.1f%%\n\n",
                    conc$n_both_sig, conc$pct_concordant, conc$pct_discordant))
      }
      
      row <- data.frame(
        Treatment_Tissue      = treatment_tissue,
        Longitudinal_Tissue   = long_tissue,
        Comparison            = paste(long_tissue, "vs", treatment_tissue),
        N_Both_Significant    = conc$n_both_sig,
        N_Concordant          = conc$n_concordant,
        Pct_Concordant        = conc$pct_concordant,
        N_Discordant          = conc$n_discordant,
        Pct_Discordant        = conc$pct_discordant,
        stringsAsFactors = FALSE
      )
      
      for (patt in names(conc$pattern_counts)) {
        row[[paste0("N_", patt)]]   <- as.numeric(conc$pattern_counts[patt])
        row[[paste0("Pct_", patt)]] <- 100 * as.numeric(conc$pattern_counts[patt]) /
          max(conc$n_both_sig, 1)
      }
      
      concordance_summary <- bind_rows(concordance_summary, row)
    }
  }
  
  write.csv(concordance_summary,
            file.path(output_dir, "Treatment_Longitudinal_Concordance_Summary.csv"),
            row.names = FALSE)
  
  cat("\n=== CONCORDANCE SUMMARY ===\n")
  print(concordance_summary %>%
          select(Comparison, N_Both_Significant, Pct_Concordant, Pct_Discordant),
        row.names = FALSE)
  
  concordance_summary
}

# ============================================================
# Execute
# ============================================================

comparison_results <- run_all_comparisons(
  data_dir   = ".",
  output_dir = "longitudinal_vs_treatment_results"
)

concordance_results <- run_concordance_analysis(
  data_dir   = ".",
  output_dir = "longitudinal_vs_treatment_results"
)