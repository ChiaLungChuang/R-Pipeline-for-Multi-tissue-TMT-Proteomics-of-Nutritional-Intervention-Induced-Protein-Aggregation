# ============================================================
# Project: Multi-tissue Proteomics — Dietary Intervention
# Script:  05_upset_plots.R
# Purpose: UpSet plots showing intersection of significant
#          proteins across tissues, for both insolubility
#          ratio and protein abundance. Generates separate
#          plots for all proteins, proteins increased under
#          treatment, and proteins decreased under treatment.
#          Also exports bidirectional protein summaries and
#          a comprehensive cross-tissue intersection report.
# Input:   Treatment_<Tissue>_processed_all_Treatment.csv
#            (output of Script 01)
# Output:  UpSet PDFs (insolubility/abundance × direction),
#          bidirectional protein CSVs,
#          protein_analysis_summary.xlsx
# Author:  Chia-Lung Chuang
# Updated: 2026-03
# NOTE:    Scripts only — raw data not included. Underlying
#          dataset is tied to a manuscript in preparation.
# Dependencies: BiocManager not required.
#   install.packages(c("ComplexUpset","writexl","ggVennDiagram"))
# ============================================================

library(tidyverse)
library(ComplexUpset)
library(writexl)

# NOTE: Set your working directory to the folder containing the input CSVs
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
# 1. Load Data
# ============================================================

load_tissue_data <- function(tissue_name) {
  f <- paste0("Treatment_", toupper(tissue_name), "_processed_all_Treatment.csv")
  if (!file.exists(f)) { warning("File not found: ", f); return(NULL) }
  d <- read_csv(f, show_col_types = FALSE)
  if (!"AC" %in% names(d) && "Protein.Accession" %in% names(d))
    d <- rename(d, AC = `Protein.Accession`)
  if (!"Gene_Name" %in% names(d) && "GN" %in% names(d))
    d <- rename(d, Gene_Name = GN)
  d
}

# Combine all tissues into a tidy long-format table
build_protein_table <- function(fc_col, pval_col, p_threshold = 0.05) {
  tissue_names <- c("tissuea","tissueb","tissuec","tissued","tissuee")
  tissue_labels <- c(tissue_a = "Tissue A", tissue_b = "Tissue B",
                     tissue_c = "Tissue C", tissue_d = "Tissue D",
                     tissue_e = "Tissue E")
  
  bind_rows(lapply(tissue_names, function(tn) {
    d <- load_tissue_data(tn)
    if (is.null(d)) return(NULL)
    label <- tissue_labels[paste0(sub("tissue", "tissue_", tn))]
    d %>%
      filter(.data[[pval_col]] <= p_threshold) %>%
      select(AC, Gene_Name, value = all_of(fc_col)) %>%
      mutate(Tissue    = label,
             Direction = ifelse(value > 0, "Increased", "Decreased"))
  }))
}

insol_data <- build_protein_table("Log2FC_ratio", "p_value_ratio")
abund_data <- build_protein_table("Log2FC_sum",   "p_value_sum")

# ============================================================
# 2. UpSet Plot Functions
# ============================================================

# Converts long-format protein table to binary presence matrix
make_binary_matrix <- function(data) {
  data %>%
    distinct(AC, Tissue) %>%
    mutate(present = TRUE) %>%
    pivot_wider(names_from = Tissue, values_from = present,
                values_fill = FALSE)
}

# Creates colored dot queries per tissue
make_queries <- function(tissues) {
  lapply(tissues, function(t)
    upset_query(set = t, color = tissue_colors[t],
                fill  = tissue_colors[t],
                only_components = "intersections_matrix"))
}

# Core UpSet plot builder
create_upset <- function(data, file, title, bar_color = "#2c3e50",
                         width = 16, height = 10) {
  binary <- make_binary_matrix(data)
  tissues <- setdiff(names(binary), "AC")
  if (length(tissues) == 0 || nrow(binary) < 2) {
    cat("Insufficient data for UpSet plot:", title, "\n"); return(invisible(NULL))
  }
  
  queries <- make_queries(intersect(tissues, names(tissue_colors)))
  n_total <- n_distinct(binary$AC)
  
  p <- upset(
    binary, tissues,
    name          = "Tissue",
    width_ratio   = 0.15,
    min_size      = 1,
    queries       = queries,
    matrix        = intersection_matrix(
      geom = geom_point(size = 5),
      segment = geom_segment(linewidth = 0.8)
    ),
    base_annotations = list(
      "Intersection size" = intersection_size(
        counts = TRUE,
        mapping = aes(fill = bar_color)
      ) + scale_fill_identity() + theme_minimal()
    ),
    themes       = upset_default_themes(text = element_text(size = 13)),
    stripes      = "white",
    sort_sets    = "descending",
    sort_intersections = "descending"
  ) + ggtitle(title, subtitle = paste(n_total, "proteins"))
  
  pdf(file, width = width, height = height)
  print(p); dev.off()
  cat("Saved:", file, "\n")
}

# ============================================================
# 3. Generate All UpSet Plots
# ============================================================

# --- Insolubility ---
create_upset(insol_data,
             "Treatment_insolubility_all_upset.pdf",
             "Insolubility: All Significant Proteins",
             bar_color = "#2c3e50")

create_upset(filter(insol_data, Direction == "Increased"),
             "Treatment_insolubility_increased_upset.pdf",
             "Insolubility: Proteins Increased Under Treatment",
             bar_color = "#e74c3c")

create_upset(filter(insol_data, Direction == "Decreased"),
             "Treatment_insolubility_decreased_upset.pdf",
             "Insolubility: Proteins Decreased Under Treatment",
             bar_color = "#3498db")

# --- Abundance ---
create_upset(abund_data,
             "Treatment_abundance_all_upset.pdf",
             "Abundance: All Significant Proteins",
             bar_color = "#2c3e50")

create_upset(filter(abund_data, Direction == "Increased"),
             "Treatment_abundance_increased_upset.pdf",
             "Abundance: Proteins Increased Under Treatment",
             bar_color = "#e74c3c")

create_upset(filter(abund_data, Direction == "Decreased"),
             "Treatment_abundance_decreased_upset.pdf",
             "Abundance: Proteins Decreased Under Treatment",
             bar_color = "#3498db")

# ============================================================
# 4. Bidirectional Protein Analysis
# ============================================================

# Proteins that are increased in some tissues but decreased in
# others — candidates for tissue-specific regulatory effects.

find_bidirectional <- function(data, prefix) {
  summary <- data %>%
    group_by(AC, Gene_Name) %>%
    summarize(
      n_tissues      = n_distinct(Tissue),
      tissues        = paste(sort(unique(Tissue)), collapse = ", "),
      increased_in   = paste(sort(Tissue[Direction == "Increased"]), collapse = ", "),
      decreased_in   = paste(sort(Tissue[Direction == "Decreased"]), collapse = ", "),
      is_bidirectional = n_distinct(Direction) > 1,
      .groups = "drop"
    ) %>%
    filter(is_bidirectional)
  
  if (nrow(summary) > 0) {
    write.csv(summary,
              paste0(prefix, "_bidirectional_summary.csv"), row.names = FALSE)
    details <- data %>%
      filter(AC %in% summary$AC) %>%
      arrange(AC, Tissue)
    write.csv(details,
              paste0(prefix, "_bidirectional_detail.csv"), row.names = FALSE)
    cat("Bidirectional proteins:", nrow(summary), "—", prefix, "\n")
  } else {
    cat("No bidirectional proteins found —", prefix, "\n")
  }
  summary
}

insol_bidir <- find_bidirectional(insol_data, "Treatment_insolubility")
abund_bidir <- find_bidirectional(abund_data, "Treatment_abundance")

# ============================================================
# 5. Cross-Tissue Intersection Statistics
# ============================================================

# For each possible tissue combination (1–5 tissues), counts
# how many proteins appear in that exact intersection.

count_intersections <- function(data) {
  tissues <- unique(data$Tissue)
  results <- data.frame(Tissues = character(), Count = integer(),
                        stringsAsFactors = FALSE)
  
  for (k in seq_along(tissues)) {
    combos <- combn(tissues, k, simplify = FALSE)
    for (combo in combos) {
      sets <- lapply(combo, function(t)
        data %>% filter(Tissue == t) %>% pull(AC) %>% unique())
      intersection <- Reduce(intersect, sets)
      results <- rbind(results, data.frame(
        Tissues = paste(combo, collapse = " & "),
        Count   = length(intersection)
      ))
    }
  }
  arrange(results, desc(Count))
}

# ============================================================
# 6. Comprehensive Summary Report
# ============================================================

tissue_summary <- function(data) {
  data %>%
    group_by(Tissue) %>%
    summarize(
      Total    = n_distinct(AC),
      Increased = n_distinct(AC[Direction == "Increased"]),
      Decreased = n_distinct(AC[Direction == "Decreased"]),
      Pct_Increased = round(Increased / Total * 100, 1),
      Pct_Decreased = round(Decreased / Total * 100, 1),
      .groups = "drop"
    )
}

excel_sheets <- list(
  "Insolubility by Tissue"     = tissue_summary(insol_data),
  "Abundance by Tissue"        = tissue_summary(abund_data),
  "Insol Intersections"        = count_intersections(insol_data),
  "Abund Intersections"        = count_intersections(abund_data),
  "Insol Increased Intersect"  = count_intersections(filter(insol_data, Direction == "Increased")),
  "Insol Decreased Intersect"  = count_intersections(filter(insol_data, Direction == "Decreased")),
  "Abund Increased Intersect"  = count_intersections(filter(abund_data, Direction == "Increased")),
  "Abund Decreased Intersect"  = count_intersections(filter(abund_data, Direction == "Decreased")),
  "Insol Bidirectional"        = if (nrow(insol_bidir) > 0) insol_bidir else data.frame(Note = "None found"),
  "Abund Bidirectional"        = if (nrow(abund_bidir) > 0) abund_bidir else data.frame(Note = "None found")
)

write_xlsx(excel_sheets, "protein_analysis_summary.xlsx")
cat("Saved: protein_analysis_summary.xlsx\n")