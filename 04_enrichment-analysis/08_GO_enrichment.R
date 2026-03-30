# ============================================================
# Project: Multi-tissue Proteomics — Dietary Intervention
# Script:  08_GO_enrichment.R
# Purpose: GO enrichment analysis (ORA via clusterProfiler)
#          for protein abundance and insolubility ratio changes
#          across 5 tissues. Runs per-tissue enrichment for
#          increased and decreased proteins (Treatment vs.
#          Control), exports dot plots (GeneRatio and fold
#          enrichment), barplots, GO network maps, and
#          Excel workbooks. Also generates cross-tissue
#          comparison heatmaps and a unified summary report.
# Input:   Treatment_<Tissue>_processed_all_Treatment.csv
#            (output of Script 07)
# Output:  GO_results/<Tissue>/  — per-tissue CSVs + PDFs
#          GO_results/cross_tissue_*_heatmap.pdf
#          GO_results/<data_type>_all_tissues_summary.txt
# Author:  Chia-Lung Chuang
# Updated: 2026-03
# NOTE:    Scripts only — raw data not included. Underlying
#          dataset is tied to a manuscript in preparation.
#          Organism: mouse (org.Mm.eg.db). Change OrgDb and
#          keyType if analyzing a different species.
#          Analysis uses p-value-only filtering (padj < 0.05,
#          no fold-change threshold applied).
# ============================================================

library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(pheatmap)
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

tissue_labels <- c(
  tissue_a = "Tissue A", tissue_b = "Tissue B",
  tissue_c = "Tissue C", tissue_d = "Tissue D",
  tissue_e = "Tissue E"
)

# ============================================================
# 1. Load Per-Tissue Processed Data
# ============================================================

load_tissue <- function(tissue_name) {
  f <- paste0("Treatment_", toupper(tissue_name), "_processed_all_Treatment.csv")
  if (!file.exists(f)) { warning("Not found: ", f); return(NULL) }
  d <- read_csv(f, show_col_types = FALSE)
  if (!"Gene_Name" %in% names(d) && "GN" %in% names(d)) d <- rename(d, Gene_Name = GN)
  d
}

tissue_names <- c("tissue_a","tissue_b","tissue_c","tissue_d","tissue_e")
tissue_data <- Filter(Negate(is.null),
                      setNames(lapply(tissue_names, load_tissue), tissue_names))

# ============================================================
# 2. Significant Gene Lists
# ============================================================

# Returns list: $increased, $decreased, $background gene symbols
# Direction is Treatment vs. Control:
#   increased = positive Log2FC (increased under treatment)
#   decreased = negative Log2FC (decreased under treatment)

get_gene_lists <- function(d, fc_col, pval_col, p_threshold = 0.05) {
  list(
    increased  = d %>% filter(.data[[pval_col]] <= p_threshold,
                              .data[[fc_col]] > 0) %>% pull(Gene_Name) %>% unique(),
    decreased  = d %>% filter(.data[[pval_col]] <= p_threshold,
                              .data[[fc_col]] < 0) %>% pull(Gene_Name) %>% unique(),
    background = d %>% pull(Gene_Name) %>% unique()
  )
}

# ============================================================
# 3. Core GO Enrichment Function
# ============================================================

# Runs enrichGO, saves dotplot (GeneRatio + fold enrichment),
# barplot, and optional network plot to a PDF.

run_go <- function(genes, background, ont, title, pdf_prefix) {
  
  if (length(genes) < 5) {
    message("  Skipping (", length(genes), " genes)")
    return(NULL)
  }
  
  res <- tryCatch(
    enrichGO(gene          = genes,
             universe      = background,
             OrgDb         = org.Mm.eg.db,
             keyType       = "SYMBOL",
             ont           = ont,
             pAdjustMethod = "BH",
             pvalueCutoff  = 0.05,
             qvalueCutoff  = 0.2,
             readable      = TRUE),
    error = function(e) { message("  enrichGO error: ", e$message); NULL }
  )
  
  if (is.null(res) || nrow(res@result) == 0) {
    message("  No significant terms for ", title, " — ", ont)
    return(NULL)
  }
  
  # Save CSV
  write.csv(res@result, paste0(pdf_prefix, "_", ont, ".csv"), row.names = FALSE)
  
  # Helper: parse "n/N" ratio string → numeric
  parse_ratio <- function(r) {
    v <- strsplit(r, "/")[[1]]
    as.numeric(v[1]) / as.numeric(v[2])
  }
  
  # Prepare plot data (top 20 by padj)
  pd <- res@result %>%
    arrange(p.adjust) %>%
    slice_head(n = 20) %>%
    mutate(
      GeneRatio_num    = sapply(GeneRatio, parse_ratio),
      BgRatio_num      = sapply(BgRatio,   parse_ratio),
      fold_enrichment  = GeneRatio_num / BgRatio_num,
      log10_padj       = -log10(p.adjust)
    )
  
  dot_theme <- theme_minimal() +
    theme(plot.title  = element_text(size = 13, face = "bold"),
          axis.text.y = element_text(size = 9, hjust = 1),
          plot.margin = margin(20, 20, 20, 20))
  
  make_dot <- function(x_col, x_lab) {
    ggplot(pd, aes(x = .data[[x_col]],
                   y = reorder(Description, .data[[x_col]]),
                   size = Count, color = log10_padj)) +
      geom_point() +
      scale_color_viridis_c(name = "-log10(adj.p)", option = "plasma",
                            begin = 0.2, end = 0.9) +
      scale_size(name = "Gene Count", range = c(3, 10)) +
      labs(title = paste(title, "-", ont), x = x_lab, y = "") +
      dot_theme
  }
  
  bar_p <- ggplot(pd, aes(x = reorder(Description, log10_padj),
                          y = log10_padj, fill = GeneRatio_num)) +
    geom_col() + coord_flip() +
    scale_fill_viridis_c(name = "Gene Ratio", option = "viridis") +
    geom_text(aes(label = Count), hjust = -0.2, size = 3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "red", alpha = 0.7) +
    labs(title = paste(title, "-", ont, "(by significance)"),
         x = "", y = "-log10(adj.p)") +
    dot_theme
  
  pdf(paste0(pdf_prefix, "_", ont, "_plots.pdf"), width = 12, height = 10)
  print(make_dot("GeneRatio_num",   "Gene Ratio"))
  print(make_dot("fold_enrichment", "Fold Enrichment"))
  print(bar_p)
  
  # Network map
  tryCatch({
    net <- enrichplot::pairwise_termsim(res)
    print(enrichplot::emapplot(net,
                               showCategory = min(30, nrow(res@result)),
                               color = "p.adjust", layout = "kk"))
  }, error = function(e) message("  Network plot skipped: ", e$message))
  
  dev.off()
  res
}

# ============================================================
# 4. Per-Tissue Analysis
# ============================================================

run_tissue_analysis <- function(tissue_name, d, data_type) {
  
  fc_col   <- if (data_type == "abundance") "Log2FC_sum"   else "Log2FC_ratio"
  pval_col <- if (data_type == "abundance") "p_value_sum"  else "p_value_ratio"
  
  genes    <- get_gene_lists(d, fc_col, pval_col)
  results  <- list()
  dir.create(paste0("GO_results/", tissue_name), recursive = TRUE, showWarnings = FALSE)
  
  for (direction in c("increased","decreased")) {
    gene_list <- genes[[direction]]
    if (length(gene_list) < 5) next
    results[[direction]] <- list()
    
    for (ont in c("BP","MF","CC")) {
      prefix <- paste0("GO_results/", tissue_name, "/", data_type, "_", direction)
      title  <- paste(tissue_labels[tissue_name], direction, data_type)
      results[[direction]][[ont]] <- run_go(gene_list, genes$background,
                                            ont, title, prefix)
    }
  }
  
  # Excel summary
  sheets <- list()
  for (dir in names(results))
    for (ont in names(results[[dir]]))
      if (!is.null(results[[dir]][[ont]]))
        sheets[[paste(tissue_name, dir, ont, sep = "_")]] <-
    results[[dir]][[ont]]@result
  
  if (length(sheets) > 0)
    write_xlsx(sheets, paste0("GO_results/", tissue_name, "/",
                              data_type, "_GO_results.xlsx"))
  
  list(sig_proteins = genes, results = results)
}

# ============================================================
# 5. Cross-Tissue Comparison Heatmap
# ============================================================

# For each direction × ontology, builds a -log10(padj) matrix
# of the top 20 terms per tissue, keeps terms in ≥2 tissues,
# and saves a pheatmap.

cross_tissue_heatmap <- function(all_results, data_type, direction, ont,
                                 top_n = 20) {
  
  term_data <- list()
  
  for (tissue in names(all_results)) {
    res <- all_results[[tissue]]$results[[direction]][[ont]]
    if (is.null(res)) next
    go_d <- res@result %>% arrange(p.adjust) %>% head(top_n)
    for (i in seq_len(nrow(go_d))) {
      term <- go_d$Description[i]
      if (!term %in% names(term_data)) term_data[[term]] <- list()
      term_data[[term]][[tissue]] <- -log10(go_d$p.adjust[i])
    }
  }
  
  tissues <- names(all_results)
  unique_terms <- names(term_data)
  if (length(unique_terms) == 0) return(invisible(NULL))
  
  mat <- matrix(0, nrow = length(unique_terms), ncol = length(tissues),
                dimnames = list(unique_terms, tissues))
  for (term in unique_terms)
    for (tissue in tissues)
      if (!is.null(term_data[[term]][[tissue]]))
        mat[term, tissue] <- term_data[[term]][[tissue]]
  
  mat <- mat[rowSums(mat > 0) >= 2, , drop = FALSE]
  if (nrow(mat) == 0) return(invisible(NULL))
  
  mat <- mat[order(rowSums(mat), decreasing = TRUE), , drop = FALSE]
  if (nrow(mat) > 30) mat <- mat[1:30, ]
  
  # Clean column names for display
  colnames(mat) <- tissue_labels[colnames(mat)]
  
  direction_label <- if (direction == "increased") "Increased" else "Decreased"
  f <- paste0("GO_results/cross_tissue_", data_type, "_",
              direction, "_", ont, "_heatmap.pdf")
  
  pdf(f, width = 12, height = max(8, nrow(mat) * 0.3))
  pheatmap(mat,
           color          = colorRampPalette(c("gray90","darkred"))(100),
           cluster_rows   = FALSE, cluster_cols = FALSE,
           display_numbers = TRUE, number_format = "%.1f",
           number_color   = "black", fontsize_number = 7,
           fontsize_row   = 9, fontsize_col   = 10,
           main = paste0("GO ", ont, " (", direction_label, ") — ",
                         toupper(data_type)),
           angle_col = 45)
  dev.off()
  cat("Saved:", f, "\n")
}

# ============================================================
# 6. Unified Summary Text
# ============================================================

write_unified_summary <- function(all_results, data_type) {
  f <- paste0("GO_results/", data_type, "_all_tissues_summary.txt")
  sink(f)
  cat("GO ENRICHMENT SUMMARY —", toupper(data_type), "\n")
  cat(strrep("=", 60), "\n\n")
  
  cat(sprintf("%-12s %8s %8s %12s\n", "Tissue","Increased","Decreased","Background"))
  for (tn in names(all_results)) {
    g <- all_results[[tn]]$sig_proteins
    cat(sprintf("%-12s %8d %8d %12d\n", tn,
                length(g$increased), length(g$decreased), length(g$background)))
  }
  cat("\n")
  
  for (tn in names(all_results)) {
    cat("\n", toupper(tissue_labels[tn]), "\n", strrep("-", 40), "\n\n")
    for (direction in c("increased","decreased")) {
      res_bp <- all_results[[tn]]$results[[direction]][["BP"]]
      if (is.null(res_bp)) next
      cat("  Top", direction, "BP terms:\n")
      top5 <- head(res_bp@result, 5)
      for (i in seq_len(nrow(top5))) {
        pf <- ifelse(top5$p.adjust[i] < 0.001, "p<0.001",
                     ifelse(top5$p.adjust[i] < 0.01,  "p<0.01", "p<0.05"))
        cat(sprintf("    %d. %s (%s, n=%d)\n", i,
                    top5$Description[i], pf, top5$Count[i]))
      }
    }
  }
  sink()
  cat("Summary saved:", f, "\n")
}

# ============================================================
# 7. Run All Analyses
# ============================================================

dir.create("GO_results", showWarnings = FALSE)

for (data_type in c("abundance","insolubility")) {
  cat("\n====", toupper(data_type), "====\n")
  
  all_results <- list()
  for (tn in names(tissue_data)) {
    cat("\n--", tn, "--\n")
    all_results[[tn]] <- run_tissue_analysis(tn, tissue_data[[tn]], data_type)
  }
  
  # Cross-tissue heatmaps
  for (direction in c("increased","decreased"))
    for (ont in c("BP","MF","CC"))
      cross_tissue_heatmap(all_results, data_type, direction, ont)
  
  write_unified_summary(all_results, data_type)
  
  saveRDS(all_results, paste0("GO_results/", data_type, "_results.rds"))
}