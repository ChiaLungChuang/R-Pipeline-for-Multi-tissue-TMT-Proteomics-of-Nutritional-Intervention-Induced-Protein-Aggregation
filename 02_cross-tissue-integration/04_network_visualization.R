# ============================================================
# Project: Multi-tissue Proteomics — Dietary Intervention
# Script:  04_network_visualization.R
# Purpose: Bipartite protein-tissue network visualization for
#          5 tissues. Constructs networks for (1) insolubility
#          ratio and (2) protein abundance, with gene nodes
#          colored by fold change gradient or tissue. Also
#          overlays mRNA expression (TPM) on the insolubility
#          network to highlight proteins with discordant
#          protein vs. transcriptional changes.
# Input:   Treatment_<Tissue>_processed_all_Treatment.csv
#            (output of Script 01)
#          mRNA_levels.xlsx — average TPM per gene per tissue
# Output:  Network PDFs (labeled and clean variants),
#          node/edge CSVs, inter-organ gene lists,
#          protein direction summary CSVs
# Author:  Chia-Lung Chuang
# Updated: 2026-03
# NOTE:    Scripts only — raw data not included. Underlying
#          dataset is tied to a manuscript in preparation.
# ============================================================

library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(ggrepel)
library(readxl)
library(scales)
library(ggnewscale)

# NOTE: Set your working directory to the folder containing the input CSVs
# setwd("path/to/your/data")

# ============================================================
# Tissue Color Palette
# ============================================================

tissue_colors <- c(
  tissue_a    = "darkgreen",
  tissue_b    = "blue4",
  tissue_c    = "purple4",
  tissue_d    = "darkgoldenrod",
  tissue_e    = "darkred",
  inter_organ = "grey50"
)

tissue_labels <- c(
  tissue_a = "Tissue A",
  tissue_b = "Tissue B",
  tissue_c = "Tissue C",
  tissue_d = "Tissue D",
  tissue_e = "Tissue E"
)

# ============================================================
# 1. Load Per-Tissue Processed Data
# ============================================================

load_tissue_data <- function(tissue_name) {
  f <- paste0("Treatment_", toupper(tissue_name), "_processed_all_Treatment.csv")
  if (!file.exists(f)) { warning("File not found: ", f); return(NULL) }
  read_csv(f, show_col_types = FALSE)
}

tissue_data <- list(
  tissue_a = load_tissue_data("tissuea"),
  tissue_b = load_tissue_data("tissueb"),
  tissue_c = load_tissue_data("tissuec"),
  tissue_d = load_tissue_data("tissued"),
  tissue_e = load_tissue_data("tissuee")
)
tissue_data <- Filter(Negate(is.null), tissue_data)
tissue_names <- names(tissue_data)

# ============================================================
# 2. Load mRNA Data
# ============================================================

# mRNA_levels.xlsx contains one sheet per tissue sub-region
# with columns: Symbol (gene name), average TPMs
# Sheet names should be mapped to the anonymized tissue labels
# below — update the mapping to match your actual sheet names.

load_mrna_data <- function(file = "mRNA_levels.xlsx") {
  sheet_to_tissue <- c(
    tissue_a1 = "tissue_a",
    tissue_a2 = "tissue_a",
    tissue_a3 = "tissue_a",
    tissue_b  = "tissue_b",
    tissue_c  = "tissue_c",
    tissue_d  = "tissue_d",
    tissue_e1 = "tissue_e",
    tissue_e2 = "tissue_e"
  )
  bind_rows(lapply(excel_sheets(file), function(s) {
    df <- read_excel(file, sheet = s)
    df %>%
      rename(Gene_Name = Symbol, average_TPMs = `average TPMs`) %>%
      mutate(Tissue = sheet_to_tissue[tolower(s)])
  })) %>%
    filter(!is.na(Tissue))
}

mrna_data <- tryCatch(load_mrna_data(), error = function(e) {
  warning("mRNA data not loaded: ", e$message); NULL
})

# ============================================================
# 3. Core Network Builder
# ============================================================

# Creates a bipartite tidygraph: gene nodes ↔ tissue nodes.
# Gene nodes carry avg_Log2FC (mean across tissues present in).
# Inter-organ genes = significant in ≥2 tissues.

build_network <- function(tissue_data_list, tissue_names,
                          fc_col   = "Log2FC_ratio",
                          pval_col = "p_value_ratio",
                          p_threshold = 0.05) {
  
  all_proteins <- bind_rows(lapply(seq_along(tissue_data_list), function(i) {
    d <- tissue_data_list[[i]]
    tissue <- tissue_names[i]
    
    # Standardize column names across script 01 variants
    if (!"AC" %in% names(d) && "Protein.Accession" %in% names(d))
      d <- rename(d, AC = `Protein.Accession`)
    if (!"Gene_Name" %in% names(d) && "GN" %in% names(d))
      d <- rename(d, Gene_Name = GN)
    
    d %>%
      filter(.data[[pval_col]] <= p_threshold) %>%
      select(AC, Gene_Name, value = all_of(fc_col)) %>%
      mutate(Tissue = tissue,
             Direction = ifelse(value > 0, "Increased", "Decreased"))
  }))
  
  gene_meta <- all_proteins %>%
    group_by(AC, Gene_Name) %>%
    summarize(
      n_tissues   = n_distinct(Tissue),
      tissues     = paste(sort(unique(Tissue)), collapse = ", "),
      gene_type   = ifelse(n_tissues > 1, "inter_tissue", "intra_tissue"),
      Tissue      = ifelse(n_tissues > 1, "inter_organ", first(Tissue)),
      Direction   = ifelse(n_distinct(Direction) == 1, first(Direction), "Mixed"),
      avg_Log2FC  = mean(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  edges    <- all_proteins %>% select(from = AC, to = Tissue) %>% distinct()
  vertices <- bind_rows(
    gene_meta %>% rename(name = AC) %>%
      mutate(type = "gene",
             node_type = ifelse(gene_type == "inter_tissue", "Inter-tissue", "Intra-tissue")),
    distinct(all_proteins, name = Tissue) %>%
      mutate(Tissue = name, type = "tissue", gene_type = NA,
             Gene_Name = NA, Direction = NA, n_tissues = NA,
             node_type = "Tissue", avg_Log2FC = 0)
  )
  
  graph <- tbl_graph(nodes = vertices, edges = edges, directed = FALSE) %>%
    activate(nodes) %>%
    mutate(degree = centrality_degree(),
           betweenness = centrality_betweenness())
  
  list(graph = graph, vertices = vertices, edges = edges,
       all_proteins = all_proteins, gene_meta = gene_meta)
}

# ============================================================
# 4. Network Plot Function
# ============================================================

# Unified plotter for both insolubility and abundance networks.
# use_gradient = TRUE: gene nodes colored by avg_Log2FC gradient
# use_gradient = FALSE: gene nodes colored by tissue
# alpha_genes: transparency for gene nodes (0–1)

plot_network <- function(net, layout_type = "kk",
                         with_labels  = TRUE,
                         use_gradient = TRUE,
                         alpha_genes  = 0.7,
                         title = "Protein-Tissue Network") {
  set.seed(42)
  
  p <- ggraph(net$graph, layout = layout_type) +
    geom_edge_link(alpha = 0.1) +
    
    # Tissue nodes
    geom_node_point(
      data = function(x) filter(x, type == "tissue"),
      aes(color = Tissue, size = node_type)
    ) +
    scale_color_manual(values = tissue_colors, name = "Tissue",
                       labels = tissue_labels) +
    scale_size_manual(
      values = c("Tissue" = 10, "Inter-tissue" = 3, "Intra-tissue" = 3),
      name = "Node Type"
    )
  
  if (use_gradient) {
    p <- p +
      ggnewscale::new_scale_color() +
      geom_node_point(
        data = function(x) filter(x, type == "gene"),
        aes(color = avg_Log2FC, size = node_type),
        alpha = alpha_genes
      ) +
      scale_color_gradientn(
        colors = c("blue4","royalblue","slateblue","mediumpurple",
                   "orchid","indianred","red4"),
        values = scales::rescale(c(-2,-1,-0.5, 0, 0.5, 1, 2)),
        limits = c(-2, 2), oob = squish,
        name   = "Log2FC\n(Treatment/Control)",
        guide  = guide_colorbar(title.position = "top",
                                barwidth = 1, barheight = 8)
      )
  } else {
    p <- p +
      ggnewscale::new_scale_color() +
      geom_node_point(
        data = function(x) filter(x, type == "gene"),
        aes(color = Tissue, size = node_type),
        alpha = alpha_genes
      ) +
      scale_color_manual(values = tissue_colors, guide = "none")
  }
  
  p <- p +
    theme_void() +
    labs(title = title) +
    theme(
      legend.position  = "right",
      legend.box       = "vertical",
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 20),
      plot.margin      = margin(0.5, 0.5, 0.5, 0.5, "cm"),
      legend.title     = element_text(size = 14),
      legend.text      = element_text(size = 12)
    )
  
  if (with_labels) {
    p <- p +
      geom_node_text(
        data  = function(x) filter(x, node_type == "Inter-tissue"),
        aes(label = Gene_Name),
        repel = TRUE, size = 3, color = "black",
        max.overlaps = Inf, box.padding = 0.5
      )
  }
  
  p
}

# ============================================================
# 5. mRNA Overlay Network
# ============================================================

# Extends the basic network to color gene nodes by log10(TPM).
# Genes with TPM = 0 (undetectable at mRNA level) are shown
# in black — these are candidates for post-transcriptional
# regulation by the dietary intervention.

plot_mrna_network <- function(net, mrna_data,
                              layout_type  = "kk",
                              with_labels  = TRUE,
                              alpha_genes  = 0.7,
                              title = "Insolubility Network with mRNA Levels") {
  # Join mRNA TPM to gene nodes
  mrna_avg <- mrna_data %>%
    group_by(Gene_Name) %>%
    summarize(average_TPMs = max(average_TPMs, na.rm = TRUE), .groups = "drop")
  
  graph_with_mrna <- net$graph %>%
    activate(nodes) %>%
    left_join(mrna_avg, by = "Gene_Name") %>%
    mutate(
      average_TPMs = replace_na(average_TPMs, 0),
      mrna_label   = ifelse(type == "gene" & average_TPMs < 1, Gene_Name, NA_character_)
    )
  
  set.seed(42)
  p <- ggraph(graph_with_mrna, layout = layout_type) +
    geom_edge_link(alpha = 0.1) +
    
    # Tissue nodes
    geom_node_point(
      data = function(x) filter(x, type == "tissue"),
      aes(color = Tissue, size = node_type)
    ) +
    scale_color_manual(values = tissue_colors, name = "Tissue",
                       labels = tissue_labels) +
    scale_size_manual(
      values = c("Tissue" = 10, "Inter-tissue" = 3, "Intra-tissue" = 3),
      name = "Node Type"
    ) +
    
    # Gene nodes: black if TPM = 0, gradient if TPM > 0
    ggnewscale::new_scale_color() +
    geom_node_point(
      data  = function(x) filter(x, type == "gene", average_TPMs == 0),
      aes(size = node_type), color = "black", alpha = alpha_genes
    ) +
    ggnewscale::new_scale_color() +
    geom_node_point(
      data  = function(x) filter(x, type == "gene", average_TPMs > 0),
      aes(color = log10(average_TPMs), size = node_type),
      alpha = alpha_genes
    ) +
    scale_color_gradientn(
      colors   = c("blue","grey","red","orange","green"),
      name     = "log10(mRNA TPM)",
      na.value = "grey50",
      guide    = guide_colorbar(title.position = "top",
                                barwidth = 1, barheight = 8)
    ) +
    theme_void() +
    labs(
      title    = title,
      subtitle = "Large nodes = tissues. Gene nodes: black = undetectable mRNA, colored = detectable."
    ) +
    theme(
      legend.position = "right",
      legend.box      = "vertical",
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 20),
      plot.subtitle   = element_text(hjust = 0.5, size = 14),
      plot.margin     = margin(0.5, 0.5, 0.5, 0.5, "cm"),
      legend.title    = element_text(size = 14),
      legend.text     = element_text(size = 12)
    )
  
  if (with_labels) {
    p <- p +
      geom_node_text(
        data  = function(x) filter(x, node_type == "Inter-tissue"),
        aes(label = Gene_Name),
        repel = TRUE, size = 3, color = "black",
        max.overlaps = Inf, box.padding = 0.5
      ) +
      geom_node_text(
        data  = function(x) filter(x, !is.na(mrna_label)),
        aes(label = mrna_label),
        repel = TRUE, size = 3, color = "black",
        max.overlaps = Inf, box.padding = 0.5
      )
  }
  
  p
}

# ============================================================
# 6. Export Helpers
# ============================================================

export_network <- function(net, prefix) {
  net$graph %>% activate(nodes) %>% as_tibble() %>%
    write.csv(paste0(prefix, "_nodes.csv"), row.names = FALSE)
  net$graph %>% activate(edges) %>% as_tibble() %>%
    write.csv(paste0(prefix, "_edges.csv"), row.names = FALSE)
  inter <- net$gene_meta %>%
    filter(gene_type == "inter_tissue") %>% arrange(desc(n_tissues))
  write.csv(inter, paste0(prefix, "_inter_tissue_genes.csv"), row.names = FALSE)
  cat("Exported:", prefix, "—", nrow(inter), "inter-tissue genes\n")
  inter
}

analyze_directions <- function(net, prefix) {
  d <- net$all_proteins %>%
    group_by(AC, Gene_Name) %>%
    summarize(
      N_tissues       = n_distinct(Tissue),
      Tissues         = paste(sort(unique(Tissue)), collapse = ", "),
      Direction_detail = paste(Tissue, Direction, sep = ": ", collapse = "; "),
      Overall         = ifelse(n_distinct(Direction) == 1, first(Direction), "Mixed"),
      .groups = "drop"
    )
  write.csv(d, paste0(prefix, "_direction_summary.csv"), row.names = FALSE)
  write.csv(filter(d, Overall == "Increased"),
            paste0(prefix, "_increased.csv"), row.names = FALSE)
  write.csv(filter(d, Overall == "Decreased"),
            paste0(prefix, "_decreased.csv"), row.names = FALSE)
  write.csv(filter(d, Overall == "Mixed"),
            paste0(prefix, "_mixed.csv"), row.names = FALSE)
  cat("Direction summary saved:", prefix, "\n")
}

# ============================================================
# 7. Build and Save All Networks
# ============================================================

save_pdf <- function(p, file, w = 14400, h = 12800)
  ggsave(file, p, width = w, height = h, units = "px", dpi = 600)

# --- Insolubility Ratio Network ---
cat("\nBuilding insolubility ratio network...\n")
insol_net <- build_network(tissue_data, tissue_names,
                           fc_col = "Log2FC_ratio", pval_col = "p_value_ratio")

save_pdf(plot_network(insol_net, with_labels = TRUE,  use_gradient = FALSE,
                      title = "Insolubility Ratio Network — Treatment vs. Control"),
         "Treatment_insolubility_network_labeled.pdf")

save_pdf(plot_network(insol_net, with_labels = FALSE, use_gradient = TRUE,
                      title = "Insolubility Ratio Network — Treatment vs. Control"),
         "Treatment_insolubility_network_gradient.pdf")

export_network(insol_net, "Treatment_insolubility_network")
analyze_directions(insol_net, "Treatment_insolubility")

# --- Abundance Network ---
cat("\nBuilding abundance network...\n")
abund_net <- build_network(tissue_data, tissue_names,
                           fc_col = "Log2FC_sum", pval_col = "p_value_sum")

save_pdf(plot_network(abund_net, with_labels = TRUE,  use_gradient = FALSE,
                      title = "Abundance Network — Treatment vs. Control"),
         "Treatment_abundance_network_labeled.pdf")

save_pdf(plot_network(abund_net, with_labels = FALSE, use_gradient = TRUE,
                      title = "Abundance Network — Treatment vs. Control"),
         "Treatment_abundance_network_gradient.pdf")

export_network(abund_net, "Treatment_abundance_network")
analyze_directions(abund_net, "Treatment_abundance")

# --- mRNA Overlay Network ---
if (!is.null(mrna_data)) {
  cat("\nBuilding mRNA overlay network...\n")
  
  save_pdf(plot_mrna_network(insol_net, mrna_data, with_labels = TRUE,
                             title = "Insolubility Network with mRNA Levels"),
           "Treatment_insolubility_network_mRNA_labeled.pdf")
  
  save_pdf(plot_mrna_network(insol_net, mrna_data, with_labels = FALSE,
                             title = "Insolubility Network with mRNA Levels"),
           "Treatment_insolubility_network_mRNA_clean.pdf")
  
  # Export: inter-tissue genes with low/undetectable mRNA
  mrna_avg <- mrna_data %>%
    group_by(Gene_Name) %>%
    summarize(average_TPMs = max(average_TPMs, na.rm = TRUE), .groups = "drop")
  
  insol_net$gene_meta %>%
    filter(gene_type == "inter_tissue") %>%
    left_join(mrna_avg, by = "Gene_Name") %>%
    mutate(average_TPMs = replace_na(average_TPMs, 0)) %>%
    filter(average_TPMs < 1) %>%
    arrange(average_TPMs) %>%
    write.csv("Treatment_inter_tissue_genes_low_mRNA.csv", row.names = FALSE)
  
  cat("mRNA overlay networks saved.\n")
} else {
  cat("mRNA data not available — skipping overlay networks.\n")
}
