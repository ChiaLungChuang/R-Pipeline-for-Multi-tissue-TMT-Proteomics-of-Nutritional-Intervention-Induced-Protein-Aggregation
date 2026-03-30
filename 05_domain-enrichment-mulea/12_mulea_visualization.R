# ============================================================
# Project: Multi-tissue Proteomics — Dietary Intervention
# Script:  12_mulea_visualization.R
# Purpose: Visualization and direction consistency analysis
#          for mulea domain enrichment results (Script 11).
#          Cross-tissue plots: overview bar chart, enrichment
#          presence heatmap, -log10(eFDR) significance heatmap,
#          bubble plot, domain type distribution, direction
#          comparison (protective vs. promotes), mirrored bar
#          chart, and tissue-sharing summary.
#          Per-tissue plots: 4-panel comprehensive PDF for
#          each tissue (overview, top domains, direction, type).
#          Direction analysis: identifies domains that switch
#          direction between tissues vs. consistent domains;
#          exports direction matrix and comparison tables.
# Input:   mulea_InterPro_results.xlsx
#            (output of Script 11)
# Output:  mulea_plots/  — all cross-tissue PDFs
#          mulea_tissue_plots/ — per-tissue PDFs
#          mulea_direction_analysis/ — direction CSVs + plots
# Author:  Chia-Lung Chuang
# Updated: 2026-03
# NOTE:    Scripts only — raw data not included. Underlying
#          dataset is tied to a manuscript in preparation.
# ============================================================

library(tidyverse)
library(pheatmap)
library(patchwork)
library(openxlsx)
library(scales)

# NOTE: Set your working directory to the folder containing mulea_InterPro_results.xlsx
# setwd("path/to/your/data")

# ============================================================
# Configuration
# ============================================================

tissue_labels <- c(
  tissue_a = "Tissue A", tissue_b = "Tissue B",
  tissue_c = "Tissue C", tissue_d = "Tissue D",
  tissue_e = "Tissue E"
)
all_tissues <- names(tissue_labels)

direction_colors <- c(
  "All_Treatment_Effects"     = "#984EA3",
  "Treatment_Protective_Only" = "#2166AC",
  "Treatment_Promotes_Only"   = "#B2182B"
)

# ============================================================
# 1. Load Results
# ============================================================

if (!file.exists("mulea_InterPro_results.xlsx"))
  stop("Run Script 11 first to generate mulea_InterPro_results.xlsx")

wb          <- loadWorkbook("mulea_InterPro_results.xlsx")
sheet_names <- setdiff(names(wb), "README")

all_results <- bind_rows(lapply(sheet_names, function(s)
  read.xlsx("mulea_InterPro_results.xlsx", sheet = s)))

cat("Loaded:", nrow(all_results), "significant enrichments\n")
cat("Tissues:", paste(unique(all_results$Tissue), collapse=", "), "\n\n")

# Translate tissue codes to display labels
all_results <- all_results %>%
  mutate(Tissue_Label = ifelse(Tissue %in% names(tissue_labels),
                               tissue_labels[Tissue], Tissue))

# eFDR safety: replace 0 with 1e-10 (should already be fixed, but guard)
all_results <- all_results %>%
  mutate(eFDR    = ifelse(eFDR    == 0, 1e-10, eFDR),
         p_value = ifelse(p_value == 0, 1e-10, p_value))

# ============================================================
# 2. Helpers
# ============================================================

dir.create("mulea_plots",          showWarnings = FALSE)
dir.create("mulea_tissue_plots",   showWarnings = FALSE)
dir.create("mulea_direction_plots",showWarnings = FALSE)

save_pdf <- function(p, path, w = 10, h = 7)
  ggsave(path, p, width = w, height = h, dpi = 300)

cr_theme <- theme_minimal(base_size = 11) +
  theme(axis.text.x   = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        panel.grid.minor = element_blank())

# Complete grid helper: ensures all tissue × domain cells present
complete_grid <- function(df, domain_col, value_col, fill_val = NA) {
  expand.grid(
    Domain = unique(df[[domain_col]]),
    Tissue = all_tissues,
    stringsAsFactors = FALSE
  ) %>%
    left_join(df %>% rename(Domain = all_of(domain_col),
                            Value  = all_of(value_col)),
              by = c("Domain","Tissue")) %>%
    mutate(Value = replace_na(Value, fill_val))
}

# ============================================================
# 3. Cross-Tissue Plots
# ============================================================

# --- 1. Overview: enrichment counts per tissue × comparison ---
overview <- all_results %>%
  count(Tissue_Label, Comparison)

p1 <- ggplot(overview, aes(x = Tissue_Label, y = n, fill = Comparison)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = n), position = position_dodge(0.7),
            vjust = -0.4, size = 3) +
  scale_fill_manual(values = direction_colors,
                    labels = c("All Treatment","Protective Only","Promotes Only")) +
  labs(title = "Enriched InterPro Domains per Tissue",
       subtitle = "eFDR < 0.05", x = "", y = "N Domains", fill = "") +
  cr_theme
save_pdf(p1, "mulea_plots/1_overview_barplot.pdf")

# --- 2. Presence heatmap: top 20 domains × all tissues ---
top20 <- all_results %>%
  filter(Comparison == "All_Treatment_Effects") %>%
  count(ontology_name) %>%
  slice_max(n, n = 20) %>%
  pull(ontology_name)

presence <- expand.grid(ontology_name = top20, Tissue = all_tissues,
                        stringsAsFactors = FALSE) %>%
  left_join(all_results %>%
              filter(Comparison == "All_Treatment_Effects") %>%
              select(ontology_name, Tissue) %>%
              mutate(Present = 1),
            by = c("ontology_name","Tissue")) %>%
  mutate(Present = replace_na(Present, 0),
         Tissue_Label = tissue_labels[Tissue])

p2 <- ggplot(presence,
             aes(x = Tissue_Label, y = reorder(ontology_name, Present, sum),
                 fill = factor(Present))) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_manual(values = c("0"="grey92","1"="#E41A1C"),
                    labels = c("Not enriched","Enriched"), name = "") +
  labs(title = "Top 20 Domains Across Tissues",
       subtitle = "All_Treatment_Effects — grey = not significant",
       x = "", y = "InterPro Domain") +
  cr_theme + theme(panel.grid = element_blank())
save_pdf(p2, "mulea_plots/2_presence_heatmap.pdf", h = 10)

# --- 3. Significance heatmap: -log10(eFDR), capped at 5 ---
top30_names <- all_results %>%
  filter(Comparison == "All_Treatment_Effects") %>%
  group_by(ontology_name) %>%
  summarize(min_eFDR = min(eFDR)) %>%
  slice_min(min_eFDR, n = 30) %>%
  pull(ontology_name)

efdr_wide <- expand.grid(ontology_name = top30_names, Tissue = all_tissues,
                         stringsAsFactors = FALSE) %>%
  left_join(all_results %>%
              filter(Comparison == "All_Treatment_Effects") %>%
              select(ontology_name, Tissue, eFDR),
            by = c("ontology_name","Tissue")) %>%
  mutate(eFDR = replace_na(eFDR, 1),
         log_score = pmin(-log10(eFDR), 5)) %>%
  pivot_wider(names_from = Tissue, values_from = log_score,
              values_fill = 0) %>%
  column_to_rownames("ontology_name")

# Rename columns to display labels
colnames(efdr_wide) <- tissue_labels[colnames(efdr_wide)]

pdf("mulea_plots/3_significance_heatmap.pdf", width = 10, height = 12)
pheatmap(efdr_wide,
         color = colorRampPalette(c("white","yellow","orange","red"))(100),
         cluster_rows = TRUE, cluster_cols = FALSE,
         fontsize = 9, fontsize_row = 8, angle_col = 45,
         main = "Top 30 Domains: -log10(eFDR), capped at 5\nWhite = not significant")
dev.off()

# --- 4. Bubble plot ---
bubble <- expand.grid(ontology_name = top30_names, Tissue = all_tissues,
                      stringsAsFactors = FALSE) %>%
  left_join(all_results %>%
              filter(Comparison == "All_Treatment_Effects") %>%
              select(ontology_name, Tissue, eFDR,
                     Count = nr_common_with_tested_elements),
            by = c("ontology_name","Tissue")) %>%
  mutate(Count      = replace_na(Count, 0),
         log_score  = ifelse(is.na(eFDR), 0, pmin(-log10(eFDR), 5)),
         Tissue_Lab = tissue_labels[Tissue])

p4 <- ggplot(bubble,
             aes(x = Tissue_Lab,
                 y = reorder(ontology_name, log_score, max),
                 size = Count, color = log_score)) +
  geom_point(alpha = 0.75) +
  scale_size_continuous(range = c(1,12), name = "Protein Count") +
  scale_color_gradientn(
    colors = c("grey85","yellow","orange","red"),
    values = rescale(c(0, 0.5, 2, 5)),
    name   = "-log10(eFDR)", limits = c(0,5), oob = squish) +
  labs(title = "Top 30 Domains: Significance Bubble Plot",
       subtitle = "Grey = not significant", x = "", y = "InterPro Domain") +
  cr_theme + theme(axis.text.y = element_text(size = 7))
save_pdf(p4, "mulea_plots/4_bubble_plot.pdf", h = 12)

# --- 5. Domain type distribution ---
if ("domain_type" %in% names(all_results)) {
  type_ct <- all_results %>%
    filter(Comparison == "All_Treatment_Effects") %>%
    count(Tissue_Label, domain_type)
  
  p5 <- ggplot(type_ct, aes(x = Tissue_Label, y = n, fill = domain_type)) +
    geom_col() +
    geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Domain Type Distribution", x = "", y = "N Domains",
         fill = "Type") +
    cr_theme
  save_pdf(p5, "mulea_plots/5_domain_types.pdf")
}

# --- 6. Protective vs. Promotes count comparison ---
dir_ct <- all_results %>%
  filter(Comparison != "All_Treatment_Effects") %>%
  count(Tissue_Label, Comparison)

p6 <- ggplot(dir_ct, aes(x = Tissue_Label, y = n, fill = Comparison)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = n), position = position_dodge(0.7),
            vjust = -0.4, size = 3) +
  scale_fill_manual(values = direction_colors[2:3],
                    labels = c("Protective","Promotes")) +
  labs(title = "Protective vs. Promotes: Enriched Domains",
       x = "", y = "N Domains", fill = "") +
  cr_theme
save_pdf(p6, "mulea_plots/6_protective_vs_promotes.pdf")

# --- 7. Mirrored bar: top 15 per direction (all tissues combined) ---
mirrored <- all_results %>%
  filter(Comparison %in% c("Treatment_Protective_Only","Treatment_Promotes_Only")) %>%
  group_by(Comparison) %>%
  slice_min(eFDR, n = 15, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    LogP        = -log10(eFDR),
    Signed_LogP = ifelse(Comparison == "Treatment_Protective_Only", -LogP, LogP),
    Label       = paste0(ontology_name, " (n=",
                         nr_common_with_tested_elements, ")")
  )

p7 <- ggplot(mirrored, aes(x = reorder(Label, abs(Signed_LogP)),
                           y = Signed_LogP, fill = Comparison)) +
  geom_col(width = 0.8) +
  geom_hline(yintercept = 0, linewidth = 1) +
  coord_flip() +
  scale_fill_manual(values = direction_colors[2:3],
                    labels = c("Protective (Decreased)","Promotes (Increased)")) +
  labs(title = "Top 15 Domains: Protective vs. Promotes (All Tissues)",
       x = "Domain (protein count)", y = "Signed -log10(eFDR)", fill = "") +
  cr_theme
save_pdf(p7, "mulea_plots/7_mirrored_barplot.pdf", h = 8)

# --- 8. Tissue-sharing summary ---
sharing <- all_results %>%
  filter(Comparison != "All_Treatment_Effects") %>%
  group_by(ontology_name) %>%
  summarize(N_Tissues = n_distinct(Tissue), .groups = "drop") %>%
  count(N_Tissues) %>%
  mutate(Pattern = paste0(N_Tissues, " tissue(s)"))

p8 <- ggplot(sharing, aes(x = reorder(Pattern, N_Tissues), y = n)) +
  geom_col(fill = "#984EA3", width = 0.7) +
  geom_text(aes(label = n), vjust = -0.4, size = 4) +
  labs(title = "Domain Tissue-Sharing (Protective + Promotes)",
       x = "N Tissues Enriched", y = "N Domains") +
  cr_theme
save_pdf(p8, "mulea_plots/8_tissue_sharing.pdf", w = 8, h = 6)

# ============================================================
# 4. Per-Tissue Comprehensive Plots
# ============================================================

make_tissue_plots <- function(tn) {
  d      <- all_results %>% filter(Tissue == tn)
  tlab   <- tissue_labels[tn]
  if (nrow(d) == 0) { cat("  No enrichments —", tn, "\n"); return(invisible()) }
  
  # Panel A: count overview
  pA <- ggplot(d %>% count(Comparison),
               aes(x = Comparison, y = n, fill = Comparison)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = n), vjust = -0.4, size = 4) +
    scale_fill_manual(values = direction_colors) +
    labs(title = paste(tlab, "— Enrichment Overview"),
         x = "", y = "N Domains") +
    cr_theme + theme(legend.position = "none")
  
  # Panel B: top 20 by significance
  top20d <- d %>% filter(Comparison == "All_Treatment_Effects") %>%
    slice_min(eFDR, n = 20, with_ties = FALSE)
  pB <- if (nrow(top20d) > 0) {
    ggplot(top20d,
           aes(x = reorder(ontology_name, -log10(eFDR)),
               y = -log10(eFDR))) +
      geom_col(fill = "#984EA3", width = 0.7) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
      geom_text(aes(label = nr_common_with_tested_elements),
                angle = 90, hjust = -0.2, size = 2.5) +
      coord_flip() +
      labs(title = paste(tlab, "— Top Domains"),
           x = "", y = "-log10(eFDR)") +
      cr_theme
  } else ggplot() + annotate("text",x=0,y=0, label="No All_Treatment enrichments") + theme_void()
  
  # Panel C: mirrored direction
  dir_d <- d %>%
    filter(Comparison %in% c("Treatment_Protective_Only","Treatment_Promotes_Only")) %>%
    group_by(Comparison) %>%
    slice_min(eFDR, n = 10, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(SignedLP  = ifelse(Comparison == "Treatment_Protective_Only",
                              -(-log10(eFDR)), -log10(eFDR)),
           Label     = paste0(ontology_name, " (",
                              nr_common_with_tested_elements, ")"))
  pC <- if (nrow(dir_d) > 0) {
    ggplot(dir_d, aes(x = reorder(Label, abs(SignedLP)),
                      y = SignedLP, fill = Comparison)) +
      geom_col(width = 0.8) +
      geom_hline(yintercept = 0, linewidth = 1) +
      coord_flip() +
      scale_fill_manual(values = direction_colors[2:3]) +
      labs(title = paste(tlab, "— Direction Comparison"),
           x = "", y = "Signed -log10(eFDR)", fill = "") +
      cr_theme + theme(axis.text.y = element_text(size = 7))
  } else ggplot() + annotate("text",x=0,y=0, label="No directional enrichments") + theme_void()
  
  # Panel D: domain type
  pD <- if ("domain_type" %in% names(d)) {
    type_d <- d %>% filter(Comparison == "All_Treatment_Effects") %>% count(domain_type)
    ggplot(type_d, aes(x = reorder(domain_type, n), y = n, fill = domain_type)) +
      geom_col(width = 0.7) +
      geom_text(aes(label = n), hjust = -0.2, size = 3.5) +
      coord_flip() +
      scale_fill_brewer(palette = "Set2") +
      labs(title = paste(tlab, "— Domain Types"), x = "", y = "N") +
      cr_theme + theme(legend.position = "none")
  } else ggplot() + theme_void()
  
  combined <- (pA + pD) / (pB | pC)
  save_pdf(combined,
           paste0("mulea_tissue_plots/", tn, "_comprehensive.pdf"),
           w = 16, h = 12)
  cat("  Saved:", tn, "_comprehensive.pdf\n")
}

cat("Creating per-tissue plots...\n")
for (tn in all_tissues) make_tissue_plots(tn)

# ============================================================
# 5. Direction Consistency Analysis
# ============================================================

cat("\nDirection consistency analysis...\n")

dir_data <- all_results %>%
  filter(Comparison %in% c("Treatment_Protective_Only",
                           "Treatment_Promotes_Only")) %>%
  mutate(Direction = ifelse(Comparison == "Treatment_Protective_Only",
                            "Protective","Promotes"))

# Per-domain, check if it appears in both directions across all tissues
dir_summary <- dir_data %>%
  group_by(ontology_id, ontology_name, Direction) %>%
  summarize(N_Tissues = n(), Tissues = paste(sort(unique(Tissue)), collapse=", "),
            .groups = "drop") %>%
  pivot_wider(names_from = Direction,
              values_from = c(N_Tissues, Tissues),
              values_fill = list(N_Tissues = 0, Tissues = ""))

# Switching = enriched in BOTH directions
switching <- dir_summary %>%
  filter(N_Tissues_Protective > 0 & N_Tissues_Promotes > 0) %>%
  arrange(desc(N_Tissues_Protective + N_Tissues_Promotes))

cat("  Switching domains (both directions):", nrow(switching), "\n")
cat("  Consistent protective:", sum(dir_summary$N_Tissues_Protective > 0 &
                                      dir_summary$N_Tissues_Promotes == 0), "\n")
cat("  Consistent promotes:",   sum(dir_summary$N_Tissues_Promotes  > 0 &
                                      dir_summary$N_Tissues_Protective == 0), "\n")

# Direction matrix (one row per domain, one col per tissue)
dir_wide <- dir_data %>%
  select(ontology_id, ontology_name, Tissue, Direction) %>%
  pivot_wider(names_from = Tissue, values_from = Direction, values_fill = "-")

# Save
write.csv(dir_summary, "mulea_direction_plots/direction_summary.csv", row.names = FALSE)
write.csv(dir_wide,    "mulea_direction_plots/direction_matrix.csv",   row.names = FALSE)
if (nrow(switching) > 0)
  write.csv(switching, "mulea_direction_plots/switching_domains.csv",  row.names = FALSE)

# Switching domain heatmap
if (nrow(switching) > 0) {
  switch_plot_d <- dir_data %>%
    filter(ontology_name %in% switching$ontology_name) %>%
    mutate(Tissue_Lab = tissue_labels[Tissue])
  
  p_switch <- ggplot(switch_plot_d,
                     aes(x = Tissue_Lab, y = ontology_name, fill = Direction)) +
    geom_tile(color = "white", linewidth = 0.8) +
    scale_fill_manual(values = c("Protective"="#2166AC","Promotes"="#B2182B")) +
    labs(title = "Switching Domains: Opposite Directions in Different Tissues",
         subtitle = "These domains show tissue-specific functional roles",
         x = "", y = "InterPro Domain", fill = "") +
    cr_theme + theme(panel.grid = element_blank())
  
  save_pdf(p_switch, "mulea_direction_plots/switching_domains_heatmap.pdf",
           w = 10, h = max(6, nrow(switching) * 0.3))
}

cat("\nAll outputs saved:\n")
cat("  mulea_plots/           — 8 cross-tissue plots\n")
cat("  mulea_tissue_plots/    — per-tissue 4-panel PDFs\n")
cat("  mulea_direction_plots/ — direction consistency CSVs and heatmap\n")