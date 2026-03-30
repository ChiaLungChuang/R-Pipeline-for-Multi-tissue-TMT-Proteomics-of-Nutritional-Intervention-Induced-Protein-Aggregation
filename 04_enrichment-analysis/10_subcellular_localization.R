# ============================================================
# Project: Multi-tissue Proteomics — Dietary Intervention
# Script:  10_subcellular_localization.R
# Purpose: Annotates treatment-responsive proteins with GO
#          Cellular Component terms (via org.Mm.eg.db), then
#          summarizes and visualizes how protein subcellular
#          location relates to the direction of insolubility
#          change (treatment-protective vs. promotes
#          aggregation) across 5 tissues.
#          Outputs: per-tissue Excel with GO annotations,
#          chi-square compartment comparison, violin/beeswarm
#          plots, enrichment bar charts, heatmaps, bubble
#          plots, per-tissue volcano plots, and an interactive
#          HTML explorer.
# Input:   Treatment_<Tissue>_processed_all_Treatment.csv
#            (output of Script 07)
# Output:  CR_proteins_subcellular_localization.xlsx
#          subcellular_chi_square_tests.csv
#          subcellular_plots/  — all PDFs and HTML
# Author:  Chia-Lung Chuang
# Updated: 2026-03
# NOTE:    Scripts only — raw data not included. Underlying
#          dataset is tied to a manuscript in preparation.
#          Organism: mouse (org.Mm.eg.db).
# ============================================================

library(tidyverse)
library(org.Mm.eg.db)
library(GO.db)
library(writexl)
library(pheatmap)
library(ggrepel)
library(ggbeeswarm)
library(plotly)
library(htmlwidgets)

# NOTE: Set your working directory to the folder containing the input CSVs
# setwd("path/to/your/data")

# ============================================================
# Tissue Configuration
# ============================================================

tissue_labels <- c(
  tissue_a = "Tissue A", tissue_b = "Tissue B",
  tissue_c = "Tissue C", tissue_d = "Tissue D",
  tissue_e = "Tissue E"
)

effect_colors <- c(
  "CR Protective"          = "#2166AC",
  "CR Promotes Aggregation" = "#B2182B"
)

# ============================================================
# 1. Load Per-Tissue Data
# ============================================================

load_tissue <- function(tissue_name) {
  f <- paste0("Treatment_", toupper(tissue_name), "_processed_all_Treatment.csv")
  if (!file.exists(f)) { warning("Not found: ", f); return(NULL) }
  d <- read_csv(f, show_col_types = FALSE)
  if (!"Gene_Name" %in% names(d) && "GN" %in% names(d)) d <- rename(d, Gene_Name = GN)
  
  # Treatment/Control direction (Treatment vs. Control format):
  #   Positive Log2FC_ratio = more insoluble under treatment → promotes aggregation
  #   Negative Log2FC_ratio = less insoluble under treatment → protective
  d %>% mutate(
    Tissue    = tissue_name,
    Direction = case_when(
      p_value_ratio < 0.05 & Log2FC_ratio > 0 ~ "Increased_in_Treatment",
      p_value_ratio < 0.05 & Log2FC_ratio < 0 ~ "Decreased_in_Treatment",
      TRUE ~ "Not_Significant"
    ),
    CR_Effect = case_when(
      Direction == "Decreased_in_Treatment" ~ "CR Protective",
      Direction == "Increased_in_Treatment" ~ "CR Promotes Aggregation",
      TRUE ~ "Not_Significant"
    )
  )
}

tissue_names <- c("tissue_a","tissue_b","tissue_c","tissue_d","tissue_e")
tissue_data <- Filter(Negate(is.null),
                      setNames(lapply(tissue_names, load_tissue), tissue_names))

# ============================================================
# 2. GO Cellular Component Annotation
# ============================================================

# Fetch CC terms for all unique gene symbols in one batch.

get_go_cc <- function(gene_symbols) {
  genes <- unique(gene_symbols[!is.na(gene_symbols) & gene_symbols != ""])
  tryCatch({
    entrez <- mapIds(org.Mm.eg.db, keys = genes, column = "ENTREZID",
                     keytype = "SYMBOL", multiVals = "first")
    go_raw <- select(org.Mm.eg.db,
                     keys   = entrez[!is.na(entrez)],
                     columns = c("SYMBOL","GO","ONTOLOGY"),
                     keytype = "ENTREZID") %>%
      filter(ONTOLOGY == "CC")
    go_terms <- select(GO.db, keys = unique(go_raw$GO),
                       columns = c("GOID","TERM"), keytype = "GOID")
    go_raw %>%
      left_join(go_terms, by = c("GO" = "GOID")) %>%
      group_by(SYMBOL) %>%
      summarize(
        GO_CC_Terms = paste(unique(TERM[!is.na(TERM)]), collapse = "; "),
        GO_CC_IDs   = paste(unique(GO[!is.na(GO)]),    collapse = "; "),
        .groups = "drop"
      )
  }, error = function(e) {
    warning("GO annotation failed: ", e$message)
    data.frame(SYMBOL = character(), GO_CC_Terms = character(), GO_CC_IDs = character())
  })
}

all_genes    <- unique(unlist(lapply(tissue_data, `[[`, "Gene_Name")))
go_cc        <- get_go_cc(all_genes)

# ============================================================
# 3. Binary Compartment Columns
# ============================================================

compartment_patterns <- list(
  Nucleus              = c("nucleus","nuclear","chromosome","nucleolus",
                           "nucleoplasm","chromatin"),
  Cytoplasm            = c("cytoplasm","cytosol","cytoplasmic"),
  Mitochondrion        = "mitochondri",
  Endoplasmic_Reticulum = c("endoplasmic reticulum","\\ber\\b"),
  Golgi                = "golgi",
  Plasma_Membrane      = c("plasma membrane","cell membrane","cell surface"),
  Ribosome             = c("ribosom","ribosomal"),
  Lysosome             = c("lysosome","lysosomal"),
  Cytoskeleton         = c("cytoskeleton","actin","tubulin","microtubule","filament"),
  Extracellular        = c("extracellular","secreted","extracellular space",
                           "extracellular region")
)
compartments <- names(compartment_patterns)

add_compartment_cols <- function(df) {
  for (comp in compartments) {
    patterns <- compartment_patterns[[comp]]
    df[[comp]] <- as.integer(
      !is.na(df$GO_CC_Terms) &
        Reduce(`|`, lapply(patterns, grepl, tolower(df$GO_CC_Terms)))
    )
  }
  df
}

# Annotate significant proteins per tissue
annotated <- lapply(names(tissue_data), function(tn) {
  sig <- tissue_data[[tn]] %>% filter(Direction != "Not_Significant")
  sig %>%
    left_join(go_cc, by = c("Gene_Name" = "SYMBOL")) %>%
    add_compartment_cols()
}) %>% setNames(names(tissue_data))

# ============================================================
# 4. Save Per-Tissue Excel with Annotations
# ============================================================

sheets <- list()
for (tn in names(annotated)) {
  d <- annotated[[tn]]
  sheets[[paste0(tn, "_All")]]         <- d
  sheets[[paste0(tn, "_Decreased")]]   <- d %>% filter(Direction == "Decreased_in_Treatment") %>% arrange(Log2FC_ratio)
  sheets[[paste0(tn, "_Increased")]]   <- d %>% filter(Direction == "Increased_in_Treatment") %>% arrange(desc(Log2FC_ratio))
}
write_xlsx(sheets, "CR_proteins_subcellular_localization.xlsx")
cat("Saved: CR_proteins_subcellular_localization.xlsx\n")

# ============================================================
# 5. Chi-Square: Compartment Comparison Between Directions
# ============================================================

chi_rows <- list()
for (tn in names(annotated)) {
  d     <- annotated[[tn]]
  dec_d <- d %>% filter(Direction == "Decreased_in_Treatment")
  inc_d <- d %>% filter(Direction == "Increased_in_Treatment")
  if (nrow(dec_d) < 5 || nrow(inc_d) < 5) next
  for (comp in compartments) {
    dec_y <- sum(dec_d[[comp]], na.rm = TRUE)
    inc_y <- sum(inc_d[[comp]], na.rm = TRUE)
    ct <- matrix(c(dec_y, nrow(dec_d)-dec_y, inc_y, nrow(inc_d)-inc_y),
                 nrow = 2)
    if (any(ct < 5)) next
    res <- tryCatch(chisq.test(ct), error = function(e) NULL)
    if (is.null(res)) next
    chi_rows[[length(chi_rows)+1]] <- data.frame(
      Tissue = tn, Compartment = comp,
      Decreased_N = dec_y,  Decreased_Pct = round(dec_y/nrow(dec_d)*100, 1),
      Increased_N = inc_y,  Increased_Pct = round(inc_y/nrow(inc_d)*100, 1),
      Chi_sq = round(res$statistic, 3),
      P_Value = round(res$p.value, 4),
      Significant = res$p.value < 0.05)
  }
}
chi_results <- bind_rows(chi_rows)
write.csv(chi_results, "subcellular_chi_square_tests.csv", row.names = FALSE)
cat("Saved: subcellular_chi_square_tests.csv\n")

# ============================================================
# 6. Prepare Long-Format Data for Plotting
# ============================================================

dir.create("subcellular_plots", showWarnings = FALSE)

all_prot <- bind_rows(annotated)

violin_data <- all_prot %>%
  pivot_longer(cols = all_of(compartments),
               names_to = "Compartment", values_to = "Present") %>%
  filter(Present == 1, CR_Effect != "Not_Significant") %>%
  mutate(
    Compartment  = gsub("_"," ", Compartment),
    Log2FC_abs   = abs(Log2FC_ratio),
    Tissue_Label = tissue_labels[Tissue]
  )

enrichment_data <- violin_data %>%
  group_by(Tissue_Label, CR_Effect, Compartment) %>%
  summarize(Count = n(), .groups = "drop") %>%
  group_by(Tissue_Label, CR_Effect) %>%
  mutate(Total = sum(Count), Percentage = 100 * Count / Total) %>%
  ungroup()

# Complete grid to expose zero-count combinations
complete_grid <- expand.grid(
  Tissue_Label = unname(tissue_labels),
  CR_Effect    = c("CR Protective","CR Promotes Aggregation"),
  Compartment  = gsub("_"," ", compartments),
  stringsAsFactors = FALSE
)
enrichment_complete <- complete_grid %>%
  left_join(enrichment_data,
            by = c("Tissue_Label","CR_Effect","Compartment")) %>%
  mutate(across(c(Count,Percentage,Total), ~replace_na(., 0)))

# ============================================================
# 7. Visualization Functions
# ============================================================

save_pdf <- function(p, name, w = 14, h = 10)
  ggsave(file.path("subcellular_plots", name), p,
         width = w, height = h, dpi = 300)

cr_theme <- theme_minimal(base_size = 10) +
  theme(axis.text.x   = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        strip.text    = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# --- 1. Violin ---
p_violin <- ggplot(violin_data,
                   aes(x = CR_Effect, y = Log2FC_abs, fill = CR_Effect)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  facet_wrap(~Compartment, ncol = 5, scales = "free_y") +
  scale_fill_manual(values = effect_colors) +
  labs(title = "Insolubility |FC| by Compartment",
       y = "|Log2 Insolubility Ratio FC|", x = "", fill = "") +
  cr_theme
save_pdf(p_violin, "1_violin.pdf", w = 16, h = 10)

# --- 2. Beeswarm (tissue rows) ---
p_bee <- ggplot(violin_data,
                aes(x = Compartment, y = Log2FC_abs, color = CR_Effect)) +
  geom_quasirandom(dodge.width = 0.8, alpha = 0.5, size = 0.8) +
  stat_summary(fun = median, geom = "crossbar", aes(group = CR_Effect),
               position = position_dodge(0.8), width = 0.4,
               linewidth = 0.5, color = "black") +
  facet_wrap(~Tissue_Label, ncol = 1) +
  scale_color_manual(values = effect_colors) +
  labs(title = "Individual Proteins by Compartment",
       y = "|Log2 Insolubility Ratio FC|", x = "", color = "") +
  cr_theme
save_pdf(p_bee, "2_beeswarm.pdf", w = 14, h = 14)

# --- 3. Enrichment bars ---
p_enrich <- ggplot(enrichment_complete,
                   aes(x = Compartment, y = Percentage, fill = CR_Effect)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = ifelse(Count > 0, Count, "")),
            position = position_dodge(0.7), vjust = -0.4, size = 2) +
  facet_wrap(~Tissue_Label, ncol = 3) +
  scale_fill_manual(values = effect_colors) +
  labs(title = "Compartment Enrichment by Tissue",
       x = "Compartment", y = "% of Proteins", fill = "") +
  cr_theme
save_pdf(p_enrich, "3_enrichment.pdf", w = 16, h = 10)

# --- 4. Heatmap ---
hm_mat <- enrichment_complete %>%
  mutate(Group = paste(Tissue_Label, CR_Effect, sep = " — ")) %>%
  select(Compartment, Group, Percentage) %>%
  pivot_wider(names_from = Group, values_from = Percentage, values_fill = 0) %>%
  column_to_rownames("Compartment") %>%
  as.matrix()

pdf(file.path("subcellular_plots","4_heatmap.pdf"), width = 14, height = 8)
pheatmap(hm_mat,
         color = colorRampPalette(c("white","#FDB462","#E31A1C"))(100),
         cluster_rows = TRUE, cluster_cols = FALSE,
         main  = "Subcellular Localization — % of Proteins (color) with counts",
         fontsize = 8, angle_col = 45, border_color = "grey80")
dev.off()

# --- 5. Bubble ---
p_bubble <- ggplot(enrichment_complete,
                   aes(x = Tissue_Label, y = Compartment,
                       size = Count, color = CR_Effect)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(range = c(1,15), name = "Count") +
  scale_color_manual(values = effect_colors, name = "") +
  labs(title = "Subcellular Distribution", x = "", y = "") +
  cr_theme
save_pdf(p_bubble, "5_bubble.pdf", w = 12, h = 8)

# --- 6. Compartment-focused bars ---
p_comp <- ggplot(enrichment_complete,
                 aes(x = Tissue_Label, y = Percentage, fill = CR_Effect)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = ifelse(Count > 0, Count, "")),
            position = position_dodge(0.7), vjust = -0.4, size = 2) +
  facet_wrap(~Compartment, ncol = 5, scales = "free_y") +
  scale_fill_manual(values = effect_colors) +
  labs(title = "Tissue Comparison Within Each Compartment",
       x = "", y = "% of Proteins", fill = "") +
  cr_theme
save_pdf(p_comp, "6_compartment_focused.pdf", w = 16, h = 10)

# --- 7. Per-tissue Volcano plots ---
comp_priority <- c("Mitochondrion","Nucleus","Endoplasmic Reticulum","Golgi",
                   "Lysosome","Plasma Membrane","Ribosome","Cytoskeleton",
                   "Extracellular","Cytoplasm","Other")

volcano_prot <- all_prot %>%
  pivot_longer(all_of(compartments), names_to = "Comp_raw", values_to = "Present") %>%
  mutate(Compartment = gsub("_"," ", Comp_raw)) %>%
  group_by(across(-c(Comp_raw, Present, Compartment))) %>%
  summarize(
    Compartment = first(Compartment[Present == 1 &
                                      Compartment %in% comp_priority]) %||% "Other",
    .groups = "drop"
  )

comp_colors <- setNames(
  c("#E41A1C","#377EB8","#4DAF4A","#E6AB02","#984EA3",
    "#F781BF","#A65628","#008080","#FF7F00","#666666","#D9D9D9"),
  c("Mitochondrion","Nucleus","Endoplasmic Reticulum","Golgi","Lysosome",
    "Plasma Membrane","Ribosome","Cytoskeleton","Extracellular","Cytoplasm","Other"))

for (tn in names(tissue_data)) {
  vd <- volcano_prot %>% filter(Tissue == tn, CR_Effect != "Not_Significant")
  if (nrow(vd) == 0) next
  top10 <- vd %>% arrange(p_value_ratio) %>% head(10)
  
  pv <- ggplot(vd, aes(x = Log2FC_ratio, y = -log10(p_value_ratio),
                       color = Compartment)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_text_repel(data = top10, aes(label = Gene_Name), color = "black",
                    size = 3, min.segment.length = 0,
                    box.padding = 0.5, segment.color = "grey50") +
    scale_color_manual(values = comp_colors, na.value = "grey80") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey60") +
    labs(title   = paste("Volcano —", tissue_labels[tn]),
         subtitle = "Colored by subcellular compartment",
         x = "Log2 Insolubility Ratio FC (Treatment vs. Control)",
         y = "-log10(p-value)") +
    theme_classic(base_size = 11)
  
  ggsave(file.path("subcellular_plots", paste0("7_volcano_", tn, ".pdf")),
         pv, width = 10, height = 7, dpi = 300)
}

# --- 8. Interactive HTML explorer ---
int_data <- violin_data %>%
  mutate(Tooltip = paste0(
    "<b>Gene:</b> ", Gene_Name, "<br>",
    "<b>Tissue:</b> ", Tissue_Label, "<br>",
    "<b>Log2FC:</b> ", round(Log2FC_ratio, 3), "<br>",
    "<b>p-value:</b> ", formatC(p_value_ratio, format="e", digits=2), "<br>",
    "<b>Compartment:</b> ", Compartment
  ))

p_int <- ggplot(int_data,
                aes(x = Compartment, y = Log2FC_ratio,
                    color = CR_Effect, text = Tooltip)) +
  geom_jitter(width = 0.3, alpha = 0.6, size = 1.2) +
  facet_wrap(~Tissue_Label) +
  scale_color_manual(values = effect_colors) +
  theme_light(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Interactive Protein Explorer",
       y = "Log2 Insolubility Ratio FC", color = "")

saveWidget(ggplotly(p_int, tooltip = "text"),
           file.path("subcellular_plots","8_interactive_explorer.html"),
           selfcontained = TRUE)

cat("\nAll outputs saved to subcellular_plots/\n")
cat("  1_violin.pdf\n  2_beeswarm.pdf\n  3_enrichment.pdf\n")
cat("  4_heatmap.pdf\n  5_bubble.pdf\n  6_compartment_focused.pdf\n")
cat("  7_volcano_<tissue>.pdf (one per tissue)\n")
cat("  8_interactive_explorer.html\n")