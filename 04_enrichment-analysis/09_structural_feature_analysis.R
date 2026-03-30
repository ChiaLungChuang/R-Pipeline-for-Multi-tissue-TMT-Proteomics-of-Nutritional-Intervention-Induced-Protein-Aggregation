# ============================================================
# Project: Multi-tissue Proteomics — Dietary Intervention
# Script:  09_structural_feature_analysis.R
# Purpose: Tests whether 8 biophysical/structural protein
#          features (length, charge, secondary structure,
#          disorder, phase separation, amyloid prediction)
#          differ between treatment-responsive and unchanged
#          proteins. Two complementary approaches:
#            (A) Mann-Whitney enrichment test: treatment-
#                affected vs. unchanged proteins (background)
#            (B) Pearson correlation: feature vs. Log2FC,
#                all proteins and in subgroups
#          Results exported as heatmaps (effect size, Pearson
#          r, row-Z-score, significant-protein enrichment,
#          cross-feature correlation) and an Excel workbook
#          with protein-level detail tables per tissue.
# Input:   <tissue>_clean_data.csv
#            (output of Script 02)
#          reference_structural_features.csv — columns:
#            name (UniProt ID), len, chg, ssa, ssb, ssc,
#            dis, pip, plc (from Molzahn et al. PNAS)
# Output:  structural_features_enrichment.xlsx
#          structural_features_correlation.xlsx
#          Fig3D_correlation_heatmap.pdf
#          Fig3E_zscore_clustered.pdf
#          Fig3F_significant_enrichment.pdf
#          Fig3G_crosscorr_<Tissue>.pdf
#          Enrichment_effect_size.pdf
# Author:  Chia-Lung Chuang
# Updated: 2026-03
# NOTE:    Scripts only — raw data not included. Underlying
#          dataset is tied to a manuscript in preparation.
#          CR_Protective / CR_Promotes labels in this script
#          map to Decreased / Increased in treatment group;
#          rename per your biological interpretation.
# ============================================================

library(tidyverse)
library(pheatmap)
library(writexl)
library(scales)

# NOTE: Set your working directory to the folder containing all input files
# setwd("path/to/your/data")

# ============================================================
# Tissue and Feature Configuration
# ============================================================

tissue_colors <- c(
  tissue_a = "darkgreen", tissue_b = "blue4",
  tissue_c = "purple4",  tissue_d = "darkgoldenrod",
  tissue_e = "darkred"
)

tissue_labels <- c(
  tissue_a = "Tissue A", tissue_b = "Tissue B",
  tissue_c = "Tissue C", tissue_d = "Tissue D",
  tissue_e = "Tissue E"
)

# Short column names from the reference file → display labels
feature_map <- c(
  len = "Protein Length",
  chg = "% Charged AA",
  ssa = "% Alpha Helix",
  ssb = "% Beta Sheet",
  ssc = "% Coiled Coil",
  dis = "% Disordered",
  pip = "Phase Separation Score",
  plc = "Amyloid Prediction"
)

features   <- names(feature_map)
feat_names <- unname(feature_map)

# ============================================================
# 1. Load Reference Data (Structural Features)
# ============================================================

# Expected columns: name (UniProt ID), len, chg, ssa, ssb,
#   ssc, dis, pip, plc
# Adjust file name/path to match your copy of the reference.

load_reference <- function(file = "reference_structural_features.csv") {
  d <- read.csv(file, skip = 1, stringsAsFactors = FALSE)
  for (col in features)
    d[[col]] <- as.numeric(as.character(d[[col]]))
  
  # Ensure plc is 0/1 binary
  if ("plc" %in% names(d)) {
    vals <- tolower(as.character(d$plc))
    d$plc <- case_when(
      vals %in% c("1","yes","y","true","positive","amyloid") ~ 1,
      vals %in% c("0","no","n","false","negative","none")    ~ 0,
      TRUE ~ suppressWarnings(as.numeric(vals))
    )
  }
  d
}

ref_data <- load_reference()
cat("Reference loaded:", nrow(ref_data), "proteins\n")

# ============================================================
# 2. Load Per-Tissue Data
# ============================================================

load_tissue <- function(tissue_name) {
  f <- paste0(tissue_name, "_clean_data.csv")
  if (!file.exists(f)) { warning("Not found: ", f); return(NULL) }
  d <- read_csv(f, show_col_types = FALSE)
  if (!"Gene_Name" %in% names(d) && "GN" %in% names(d)) d <- rename(d, Gene_Name = GN)
  
  # Extract UniProt ID from pipe-delimited Protein Accession
  if ("Protein.Accession" %in% names(d) || "Protein Accession" %in% names(d)) {
    ac_col <- if ("Protein.Accession" %in% names(d)) "Protein.Accession" else "Protein Accession"
    d$name <- sapply(strsplit(d[[ac_col]], "\\|"),
                     function(x) if (length(x) >= 2) x[2] else x[1])
  }
  
  # Classify (Treatment vs. Control):
  #   Protective = treatment decreases insolubility (negative Log2FC_ratio)
  #   Promotes   = treatment increases insolubility (positive Log2FC_ratio)
  d <- d %>% mutate(
    CR_Effect = case_when(
      p_value_ratio < 0.05 & Log2FC_ratio < 0 ~ "CR_Protective",
      p_value_ratio < 0.05 & Log2FC_ratio > 0 ~ "CR_Promotes",
      TRUE ~ "No_Change"
    )
  )
  d$Tissue <- tissue_name
  d
}

tissue_names <- c("tissue_a","tissue_b","tissue_c","tissue_d","tissue_e")
tissue_data <- Filter(Negate(is.null),
                      setNames(lapply(tissue_names, load_tissue), tissue_names))

# Merged list (tissue data + structural features)
merged_data <- lapply(names(tissue_data), function(tn) {
  merge(tissue_data[[tn]], ref_data, by = "name", all.x = FALSE)
}) %>% setNames(names(tissue_data))

# ============================================================
# 3. Helper: Safe Statistical Tests
# ============================================================

safe_mwu <- function(target, bg, min_n = 3) {
  t <- na.omit(target); b <- na.omit(bg)
  if (length(t) < min_n || length(b) < min_n)
    return(list(success = FALSE))
  if (sd(t, na.rm = TRUE) == 0 && sd(b, na.rm = TRUE) == 0)
    return(list(success = FALSE))
  r <- tryCatch(wilcox.test(t, b, exact = FALSE), error = function(e) NULL)
  if (is.null(r)) return(list(success = FALSE))
  list(success     = TRUE,
       p.value     = r$p.value,
       effect_size = median(t) - median(b),
       n_target    = length(t),
       n_bg        = length(b))
}

safe_cor <- function(x, y, min_n = 10) {
  idx <- !is.na(x) & !is.na(y) & is.finite(x) & is.finite(y)
  if (sum(idx) < min_n || sd(x[idx]) == 0) return(list(success = FALSE))
  r <- tryCatch(cor.test(x[idx], y[idx], method = "pearson"),
                error = function(e) NULL)
  if (is.null(r) || is.na(r$p.value)) return(list(success = FALSE))
  list(success  = TRUE,
       r        = as.numeric(r$estimate),
       p.value  = r$p.value,
       n        = sum(idx))
}

# ============================================================
# 4. Enrichment Analysis (Mann-Whitney U)
# ============================================================

# Compares CR_Protective + CR_Promotes (target) vs. No_Change (background)

run_enrichment <- function() {
  all_rows <- list()
  
  for (tn in names(merged_data)) {
    d <- merged_data[[tn]]
    target <- d %>% filter(CR_Effect %in% c("CR_Protective","CR_Promotes"))
    bg     <- d %>% filter(CR_Effect == "No_Change")
    
    for (feat in features) {
      feat_type <- if (feat == "plc") "binary" else "continuous"
      
      if (feat_type == "binary") {
        # Fisher's exact for binary amyloid prediction
        t1 <- sum(target[[feat]] == 1, na.rm = TRUE)
        t0 <- sum(target[[feat]] == 0, na.rm = TRUE)
        b1 <- sum(bg[[feat]] == 1, na.rm = TRUE)
        b0 <- sum(bg[[feat]] == 0, na.rm = TRUE)
        if (any(c(t1+t0, b1+b0) < 3)) next
        res <- tryCatch(
          fisher.test(matrix(c(t1,t0,b1,b0), nrow = 2)),
          error = function(e) NULL)
        if (is.null(res)) next
        all_rows[[length(all_rows)+1]] <- data.frame(
          Tissue = tn, Feature = feature_map[feat],
          Test = "Fisher", N_Target = t1+t0, N_Background = b1+b0,
          Effect_Size = t1/(t1+t0) - b1/(b1+b0),
          P_Value = res$p.value)
      } else {
        res <- safe_mwu(target[[feat]], bg[[feat]])
        if (!res$success) next
        all_rows[[length(all_rows)+1]] <- data.frame(
          Tissue = tn, Feature = feature_map[feat],
          Test = "MWU", N_Target = res$n_target, N_Background = res$n_bg,
          Effect_Size = res$effect_size, P_Value = res$p.value)
      }
    }
  }
  
  bind_rows(all_rows) %>%
    group_by(Tissue) %>%
    mutate(P_Adjusted_BH = p.adjust(P_Value, "BH")) %>%
    ungroup() %>%
    mutate(Sig_Stars = case_when(
      P_Value < 0.001 ~ "***",
      P_Value < 0.01  ~ "**",
      P_Value < 0.05  ~ "*",
      TRUE            ~ ""))
}

enrichment_results <- run_enrichment()

# ============================================================
# 5. Correlation Analysis (Pearson r)
# ============================================================

# Three sub-analyses per tissue × feature:
#   All_Proteins, Significant_Only, CR_Protective_Only, CR_Promotes_Only

run_correlations <- function() {
  all_rows <- list()
  
  for (tn in names(merged_data)) {
    d <- merged_data[[tn]]
    subsets <- list(
      All_Proteins       = d,
      Significant_Only   = d %>% filter(CR_Effect != "No_Change"),
      CR_Protective_Only = d %>% filter(CR_Effect == "CR_Protective"),
      CR_Promotes_Only   = d %>% filter(CR_Effect == "CR_Promotes")
    )
    
    for (analysis in names(subsets)) {
      sub <- subsets[[analysis]]
      if (nrow(sub) < 10) next
      for (feat in features) {
        res <- safe_cor(sub[[feat]], sub$Log2FC_ratio)
        if (!res$success) next
        all_rows[[length(all_rows)+1]] <- data.frame(
          Tissue = tn, Analysis = analysis,
          Feature = feature_map[feat],
          r = res$r, P_Value = res$p.value, N = res$n)
      }
    }
  }
  
  bind_rows(all_rows) %>%
    group_by(Tissue, Analysis) %>%
    mutate(P_Adjusted_BH = p.adjust(P_Value, "BH")) %>%
    ungroup() %>%
    mutate(Sig_Stars = case_when(
      P_Value < 0.001 ~ "***",
      P_Value < 0.01  ~ "**",
      P_Value < 0.05  ~ "*",
      TRUE            ~ ""))
}

correlation_results <- run_correlations()

# ============================================================
# 6. Build Heatmap Matrices
# ============================================================

to_matrix <- function(df, val_col, tissue_col = "Tissue",
                      feature_col = "Feature") {
  df %>%
    select(all_of(c(tissue_col, feature_col, val_col))) %>%
    pivot_wider(names_from = all_of(tissue_col),
                values_from = all_of(val_col)) %>%
    column_to_rownames(feature_col) %>%
    as.matrix()
}

sig_label_matrix <- function(df, val_col, sig_col = "Sig_Stars",
                             tissue_col = "Tissue",
                             feature_col = "Feature") {
  vals <- to_matrix(df, val_col, tissue_col, feature_col)
  sigs <- to_matrix(df, sig_col, tissue_col, feature_col)
  matrix(sprintf("%.2f\n%s", vals, sigs),
         nrow = nrow(vals), ncol = ncol(vals),
         dimnames = dimnames(vals))
}

# Correlation heatmap (All_Proteins analysis)
cor_all <- correlation_results %>% filter(Analysis == "All_Proteins")
cor_mat  <- to_matrix(cor_all, "r")
cor_labs <- sig_label_matrix(cor_all, "r")

# Z-score matrix (row-standardised)
zscore_mat <- t(scale(t(cor_mat)))
zscore_mat[is.nan(zscore_mat)] <- 0

# Enrichment effect size matrix
enrich_mat  <- to_matrix(enrichment_results, "Effect_Size")
enrich_labs <- sig_label_matrix(enrichment_results, "Effect_Size")

# ============================================================
# 7. Visualizations
# ============================================================

# Rename matrix columns to display labels
rename_cols <- function(mat) {
  colnames(mat) <- tissue_labels[colnames(mat)]
  mat
}

cor_mat_d  <- rename_cols(cor_mat)
zscore_d   <- rename_cols(zscore_mat)
enrich_d   <- rename_cols(enrich_mat)
cor_labs_d <- { m <- cor_labs; colnames(m) <- tissue_labels[colnames(m)]; m }

heat_theme_args <- list(
  cluster_rows = TRUE, cluster_cols = FALSE,
  border_color = "grey80", fontsize = 11,
  fontsize_row = 9, fontsize_col = 10, angle_col = 45
)

# --- Fig 3D: Pearson r heatmap ---
pdf("Fig3D_correlation_heatmap.pdf", width = 10, height = 7)
do.call(pheatmap, c(list(
  mat             = cor_mat_d,
  color           = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
  breaks          = seq(-0.5, 0.5, length.out = 101),
  display_numbers = cor_labs_d, number_format = "%s",
  main            = "Pearson r: Structural Feature vs. Insolubility FC"
), heat_theme_args))
dev.off()

# --- Fig 3E: Z-score heatmap (clustered) ---
valid_rows <- apply(zscore_d, 1, function(x) !any(is.na(x)))
if (sum(valid_rows) >= 2) {
  zd_clean  <- zscore_d[valid_rows, ]
  pdf("Fig3E_zscore_clustered.pdf", width = 10, height = 8)
  do.call(pheatmap, c(list(
    mat             = zd_clean,
    color           = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
    breaks          = seq(-3, 3, length.out = 101),
    main            = "Z-score of Feature-FC Correlations (row-standardised)"
  ), heat_theme_args))
  dev.off()
}

# --- Enrichment effect size heatmap ---
pdf("Enrichment_effect_size.pdf", width = 10, height = 7)
do.call(pheatmap, c(list(
  mat             = enrich_d,
  color           = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
  display_numbers = enrich_labs, number_format = "%s",
  main            = "Enrichment: Effect Size (Treatment-affected vs. Unchanged)"
), heat_theme_args))
dev.off()

# --- Fig 3F: Correlation of significance with features ---
# Binary outcome: is protein significant (p_ratio < 0.05)?
sig_cor_rows <- list()
for (tn in names(merged_data)) {
  d <- merged_data[[tn]] %>%
    mutate(Is_Sig = as.numeric(p_value_ratio < 0.05))
  for (feat in features) {
    res <- safe_cor(d[[feat]], d$Is_Sig)
    if (!res$success) next
    sig_cor_rows[[length(sig_cor_rows)+1]] <- data.frame(
      Tissue = tn, Feature = feature_map[feat],
      r = res$r, P_Value = res$p.value, N = res$n)
  }
}
sig_cor <- bind_rows(sig_cor_rows) %>%
  mutate(Sig_Stars = case_when(P_Value < 0.001~"***", P_Value < 0.01~"**",
                               P_Value < 0.05~"*", TRUE~""))

if (nrow(sig_cor) > 0) {
  sig_mat  <- rename_cols(to_matrix(sig_cor, "r"))
  sig_labs <- { m <- sig_label_matrix(sig_cor, "r"); colnames(m) <- tissue_labels[colnames(m)]; m }
  pdf("Fig3F_significant_enrichment.pdf", width = 10, height = 7)
  do.call(pheatmap, c(list(
    mat             = sig_mat,
    color           = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
    display_numbers = sig_labs, number_format = "%s",
    main            = "Fig 3F: Feature Association with Significant Proteins"
  ), heat_theme_args))
  dev.off()
}

# --- Fig 3G: Cross-feature correlation per tissue ---
for (tn in names(merged_data)) {
  cr_d <- merged_data[[tn]] %>%
    filter(CR_Effect %in% c("CR_Protective","CR_Promotes")) %>%
    select(all_of(features)) %>%
    select(where(~ sum(!is.na(.)) > 2 && var(., na.rm = TRUE) > 0))
  
  if (ncol(cr_d) < 2 || nrow(cr_d) < 10) next
  cm <- cor(cr_d, use = "pairwise.complete.obs")
  rownames(cm) <- colnames(cm) <- feature_map[colnames(cm)]
  
  pdf(paste0("Fig3G_crosscorr_", tn, ".pdf"), width = 9, height = 9)
  pheatmap(cm,
           color          = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
           breaks         = seq(-1, 1, length.out = 101),
           display_numbers = TRUE, number_format = "%.2f",
           main           = paste0("Cross-feature Correlations — ",
                                   tissue_labels[tn]),
           border_color   = "grey80", fontsize = 10,
           cluster_rows   = TRUE, cluster_cols = TRUE)
  dev.off()
}

# ============================================================
# 8. Export Excel Workbooks
# ============================================================

# Protein-level detail tables per tissue
prot_sheets <- list()
for (tn in names(merged_data)) {
  d <- merged_data[[tn]] %>%
    filter(CR_Effect %in% c("CR_Protective","CR_Promotes")) %>%
    select(any_of(c("name","Gene_Name","CR_Effect","Log2FC_ratio",
                    "p_value_ratio","Log2FC_soluble","p_value_soluble",
                    "Log2FC_insoluble","p_value_insoluble")),
           all_of(features)) %>%
    rename(any_of(c(Protein_Length = "len", Percent_Charged_AA = "chg",
                    Percent_Alpha_Helix = "ssa", Percent_Beta_Sheet = "ssb",
                    Percent_Coiled_Coil = "ssc", Percent_Disordered = "dis",
                    Phase_Separation_Score = "pip", Amyloid_Prediction = "plc")))
  if (nrow(d) == 0) next
  prot_sheets[[paste0(tn, "_all")]] <- d
  prot_sheets[[paste0(tn, "_protective")]] <- d %>% filter(CR_Effect == "CR_Protective")
  prot_sheets[[paste0(tn, "_promotes")]]   <- d %>% filter(CR_Effect == "CR_Promotes")
}

write_xlsx(c(
  list(Enrichment_Results = enrichment_results),
  prot_sheets
), "structural_features_enrichment.xlsx")

write_xlsx(list(
  All_Correlations   = correlation_results,
  Significant_Only   = filter(correlation_results, P_Value < 0.05),
  Sig_Feature_Corr   = sig_cor
), "structural_features_correlation.xlsx")

cat("\nOutputs saved:\n")
cat("  structural_features_enrichment.xlsx\n")
cat("  structural_features_correlation.xlsx\n")
cat("  Fig3D_correlation_heatmap.pdf\n")
cat("  Fig3E_zscore_clustered.pdf\n")
cat("  Fig3F_significant_enrichment.pdf\n")
cat("  Fig3G_crosscorr_<tissue>.pdf\n")
cat("  Enrichment_effect_size.pdf\n")