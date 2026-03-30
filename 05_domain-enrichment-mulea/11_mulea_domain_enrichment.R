# ============================================================
# Project: Multi-tissue Proteomics — Dietary Intervention
# Script:  11_mulea_domain_enrichment.R
# Purpose: Protein domain enrichment analysis using mulea ORA
#          (overrepresentation analysis with empirical FDR)
#          against two ontologies:
#            (A) InterPro domains — GMT file from InterPro
#            (B) GO Molecular Function — via muleaData
#          Runs three comparisons per tissue:
#            All_CR_Effects, Treatment_Protective_Only,
#            Treatment_Promotes_Only
#          Results exported to Excel workbooks.
#          Note on mulea: standard hypergeometric p-values
#          can be unreliable for small gene sets; mulea's
#          empirical FDR (eFDR) via permutation is used
#          instead. Known issue: eFDR = 0 for very strong
#          enrichments — these are replaced with 1e-10 for
#          downstream log-transformation.
# Input:   Treatment_<Tissue>_processed_all_Treatment.csv
#            (output of Script 07)
#          mouse_interpro_domains.gmt.txt — InterPro GMT
#          domain_annotations_cache.csv  — optional cache
#            with domain_id, domain_name, domain_type cols
# Output:  mulea_InterPro_results.xlsx
#          mulea_GO_MF_results.xlsx
# Author:  Chia-Lung Chuang
# Updated: 2026-03
# NOTE:    Scripts only — raw data not included. Underlying
#          dataset is tied to a manuscript in preparation.
#          Organism: mouse. Change OrgDb/muleaData resource
#          for different species.
#          Treatment_Protective = decreased insolubility
#          Treatment_Promotes   = increased insolubility
# ============================================================

library(tidyverse)
library(mulea)
library(writexl)

# NOTE: Set your working directory to the folder containing all input files
# setwd("path/to/your/data")

# ============================================================
# Tissue Configuration
# ============================================================

tissue_labels <- c(
  tissue_a = "Tissue A", tissue_b = "Tissue B",
  tissue_c = "Tissue C", tissue_d = "Tissue D",
  tissue_e = "Tissue E"
)

comparisons <- list(
  All_Treatment_Effects       = c("Treatment_Protective","Treatment_Promotes"),
  Treatment_Protective_Only   = c("Treatment_Protective"),
  Treatment_Promotes_Only     = c("Treatment_Promotes")
)

# ============================================================
# 1. Load Per-Tissue Data
# ============================================================

load_tissue <- function(tissue_name) {
  f <- paste0("Treatment_", toupper(tissue_name),
              "_processed_all_Treatment.csv")
  if (!file.exists(f)) { warning("Not found: ", f); return(NULL) }
  d <- read_csv(f, show_col_types = FALSE)
  if (!"Gene_Name" %in% names(d) && "GN" %in% names(d)) d <- rename(d, Gene_Name = GN)
  
  # Extract UniProt ID from pipe-delimited Protein Accession
  ac_col <- intersect(c("Protein.Accession","Protein Accession"), names(d))[1]
  if (!is.na(ac_col))
    d$UniProt_ID <- sapply(strsplit(d[[ac_col]], "\\|"),
                           function(x) if (length(x) >= 2) x[2] else x[1])
  
  # Direction classification (Treatment vs. Control):
  #   Protective = treatment decreases insolubility (negative Log2FC_ratio)
  #   Promotes   = treatment increases insolubility (positive Log2FC_ratio)
  d <- d %>% mutate(
    Treatment_Effect = case_when(
      p_value_ratio < 0.05 & Log2FC_ratio < 0 ~ "Treatment_Protective",
      p_value_ratio < 0.05 & Log2FC_ratio > 0 ~ "Treatment_Promotes",
      TRUE ~ "No_Change"
    )
  )
  d$Tissue <- tissue_name
  d
}

tissue_names <- c("tissue_a","tissue_b","tissue_c","tissue_d","tissue_e")
tissue_data <- Filter(Negate(is.null),
                      setNames(lapply(tissue_names, load_tissue), tissue_names))

# ============================================================
# 2. Core mulea ORA Runner
# ============================================================

# Runs ORA for all tissues × comparisons against a given GMT.
# Returns a list with two elements:
#   $all_results  — combined data.frame of all tests
#   $excel_sheets — named list ready for write_xlsx()

run_mulea_ora <- function(gmt, gene_map = NULL,
                          n_perm = 10000, min_genes = 5) {
  
  all_rows  <- list()
  xl_sheets <- list()
  
  for (tn in names(tissue_data)) {
    d          <- tissue_data[[tn]]
    background <- unique(na.omit(d$UniProt_ID))
    cat("  Tissue:", tn, "| background:", length(background), "\n")
    
    for (comp_name in names(comparisons)) {
      cats  <- comparisons[[comp_name]]
      genes <- unique(na.omit(d$UniProt_ID[d$Treatment_Effect %in% cats]))
      if (length(genes) < min_genes) {
        cat("   ", comp_name, ": skipped (<", min_genes, " proteins)\n"); next
      }
      
      res <- tryCatch({
        set.seed(42)
        ora_model <- ora(gmt                         = gmt,
                         element_names               = genes,
                         background_element_names    = background,
                         p_value_adjustment_method   = "eFDR",
                         number_of_permutations      = n_perm,
                         nthreads                    = 1)
        run_test(ora_model)
      }, error = function(e) { message("   Error: ", e$message); NULL })
      
      if (is.null(res) || nrow(res) == 0) next
      
      # Annotate overlapping proteins with gene names (if map provided)
      if (!is.null(gene_map) && "list_of_values" %in% names(gmt)) {
        gmt_lookup <- setNames(gmt$list_of_values, gmt$ontology_id)
        res$Overlapping_Proteins <- sapply(res$ontology_id, function(did) {
          ids    <- gmt_lookup[[did]]
          overlap <- intersect(ids, genes)
          labels  <- sapply(overlap, function(id) {
            gn <- gene_map[id]
            if (is.na(gn)) id else paste0(gn, " (", id, ")")
          })
          paste(labels, collapse = "; ")
        })
      }
      
      res <- res %>%
        # eFDR = 0 arises for very strong enrichments; replace for log-safety
        mutate(eFDR    = ifelse(eFDR    == 0, 1e-10, eFDR),
               p_value = ifelse(p_value == 0, 1e-10, p_value),
               Tissue     = tn,
               Comparison = comp_name,
               N_Test     = length(genes),
               N_Background = length(background)) %>%
        filter(eFDR < 0.05) %>%
        arrange(eFDR)
      
      n_sig <- nrow(res)
      cat("   ", comp_name, ":", n_sig, "significant domains\n")
      if (n_sig == 0) next
      
      all_rows[[length(all_rows)+1]] <- res
      xl_name <- substr(paste0(substr(tn,1,8), "_", comp_name), 1, 31)
      xl_sheets[[xl_name]] <- res
    }
  }
  
  list(all_results  = if (length(all_rows) > 0) bind_rows(all_rows)
       else data.frame(),
       excel_sheets = xl_sheets)
}

# ============================================================
# 3. InterPro Domain Analysis
# ============================================================

cat("\n=== InterPro Domain Enrichment ===\n")

# Load GMT
gmt_raw <- read_tsv("mouse_interpro_domains.gmt.txt",
                    show_col_types = FALSE)
gmt_list <- lapply(strsplit(gmt_raw$list_of_values, ", "), as.character)
names(gmt_list) <- gmt_raw$ontology_id
gmt_interpro    <- list_to_gmt(gmt_list)
gmt_interpro$ontology_name <- gmt_raw$ontology_name[
  match(gmt_interpro$ontology_id, gmt_raw$ontology_id)]
gmt_interpro <- as.data.frame(gmt_interpro)

# Optional domain metadata (type: Family, Domain, Repeat, …)
domain_meta <- NULL
if (file.exists("domain_annotations_cache.csv")) {
  domain_meta <- read_csv("domain_annotations_cache.csv",
                          show_col_types = FALSE) %>%
    distinct(domain_id, .keep_all = TRUE) %>%
    rename(ontology_id = domain_id)
  cat("Domain metadata cache loaded:", nrow(domain_meta), "entries\n")
}

# Build gene → name map from all tissues
gene_map_interpro <- unlist(lapply(tissue_data, function(d)
  setNames(d$Gene_Name, d$UniProt_ID)), use.names = FALSE)
gene_map_interpro <- setNames(
  gene_map_interpro,
  unlist(lapply(tissue_data, `[[`, "UniProt_ID"), use.names = FALSE))
gene_map_interpro <- gene_map_interpro[!duplicated(names(gene_map_interpro))]

interpro_out <- run_mulea_ora(gmt_interpro, gene_map = gene_map_interpro)

# Join domain type metadata
if (!is.null(domain_meta) && nrow(interpro_out$all_results) > 0) {
  interpro_out$all_results <- interpro_out$all_results %>%
    left_join(domain_meta %>% select(ontology_id, domain_type),
              by = "ontology_id") %>%
    mutate(domain_type = replace_na(domain_type, "Unknown"))
  
  interpro_out$excel_sheets <- lapply(interpro_out$excel_sheets, function(s) {
    s %>% left_join(domain_meta %>% select(ontology_id, domain_type),
                    by = "ontology_id") %>%
      mutate(domain_type = replace_na(domain_type, "Unknown"))
  })
}

readme_ip <- data.frame(
  Column = c("Tissue","Comparison","ontology_id","ontology_name","domain_type",
             "eFDR","p_value","nr_common_with_tested_elements",
             "Overlapping_Proteins","N_Test","N_Background"),
  Description = c(
    "Anonymised tissue name (tissue_a … tissue_e)",
    "All_Treatment_Effects | Treatment_Protective_Only | Treatment_Promotes_Only",
    "InterPro accession (e.g. IPR012345)",
    "Domain name from InterPro",
    "Family | Domain | Repeat | Active_site | etc.",
    "Empirical FDR (permutation-based). Values of 0 replaced with 1e-10.",
    "Raw hypergeometric p-value (zeros also replaced with 1e-10).",
    "Number of treatment-affected proteins containing this domain",
    "Gene names and UniProt IDs of overlapping proteins",
    "Total proteins in the test set for this comparison",
    "Total proteins in the background"
  ))

write_xlsx(c(list(README = readme_ip,
                  All_Significant = interpro_out$all_results),
             interpro_out$excel_sheets),
           "mulea_InterPro_results.xlsx")
cat("Saved: mulea_InterPro_results.xlsx\n")

# ============================================================
# 4. GO Molecular Function Analysis
# ============================================================

cat("\n=== GO Molecular Function Enrichment ===\n")

# Requires muleaData and ExperimentHub
gmt_gomf <- tryCatch({
  library(muleaData); library(ExperimentHub)
  eh      <- ExperimentHub()
  hub     <- query(eh, "muleaData")
  matches <- hub[grep("GO.*Mus_musculus.*UniprotID", hub$title, TRUE)]
  mf_hits <- matches[grep("MF|Molecular", matches$title, TRUE)]
  if (length(mf_hits) == 0) stop("No GO MF resource found")
  raw <- as.data.frame(eh[[names(mf_hits)[1]]])
  names(raw) <- sub("ontologyId",   "ontology_id",   names(raw))
  names(raw) <- sub("ontologyName", "ontology_name", names(raw))
  names(raw) <- sub("listOfValues", "list_of_values", names(raw))
  raw <- raw[!is.na(raw$ontology_id) & !is.na(raw$list_of_values), ]
  cat("GO MF resource loaded:", nrow(raw), "terms\n")
  raw
}, error = function(e) {
  message("GO MF unavailable (", e$message, ") — skipping GO analysis")
  NULL
})

if (!is.null(gmt_gomf)) {
  gomf_out <- run_mulea_ora(gmt_gomf)
  
  readme_go <- data.frame(
    Column      = c("Tissue","Comparison","ontology_id","ontology_name",
                    "eFDR","p_value","nr_common_with_tested_elements"),
    Description = c(
      "Anonymised tissue name",
      "All_Treatment_Effects | Treatment_Protective_Only | Treatment_Promotes_Only",
      "GO accession (e.g. GO:0003674)",
      "GO term name",
      "Empirical FDR (zeros replaced with 1e-10)",
      "Raw p-value (zeros replaced with 1e-10)",
      "Proteins with this GO term in the test set"
    ))
  
  write_xlsx(c(list(README = readme_go,
                    All_Significant = gomf_out$all_results),
               gomf_out$excel_sheets),
             "mulea_GO_MF_results.xlsx")
  cat("Saved: mulea_GO_MF_results.xlsx\n")
}

cat("\nDone. Outputs:\n")
cat("  mulea_InterPro_results.xlsx\n")
if (!is.null(gmt_gomf)) cat("  mulea_GO_MF_results.xlsx\n")