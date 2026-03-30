# R Pipeline for Multi-tissue TMT Proteomics of Nutritional Intervention-Induced Protein Aggregation

Analytical pipeline for a multi-tissue proteomics study examining how nutritional intervention affects protein insolubility (aggregation propensity) and protein abundance across five mouse tissues.

> **Note:** This repository contains scripts only. Raw mass spectrometry data and processed protein quantification files are not included, as the underlying manuscript is in preparation. Tissue identities and condition labels have been anonymized pending publication.

---

## Study Overview

This project investigates the relationship between nutritional intervention and protein aggregation across multiple tissues using TMT (tandem mass tag) quantitative proteomics. For each tissue, proteins were fractionated into soluble and insoluble fractions. The **insolubility ratio** (insoluble/soluble) serves as a proxy for protein aggregation propensity. Changes in this ratio — as well as changes in total protein abundance — were compared between treatment and control groups.

**Five tissues analyzed:** Tissue A, Tissue B, Tissue C, Tissue D, Tissue E

**Three quantitative readouts per protein:**
- Soluble fraction fold change
- Insoluble fraction fold change
- Insolubility ratio fold change (primary outcome)

**Statistical approach:** Welch's t-test with Benjamini–Hochberg FDR correction; significance threshold p < 0.05 (unadjusted p-value, consistent with prior work in the field). No fold-change threshold was applied.

---

## Repository Structure

```
CR-proteomics-mouse-aging/
├── 01_per-tissue-analysis/
│   ├── 01_abundance_vs_ratio_plots.R
│   ├── 02_soluble_vs_insoluble_plots.R
│   └── 03_longitudinal_comparison.R
├── 02_cross-tissue-integration/
│   ├── 04_network_visualization.R
│   └── 05_upset_plots.R
├── 03_biophysical-correlations/
│   ├── 06_biophysical_correlations.R
│   └── 07_fraction_comparison_plots.R
├── 04_enrichment-analysis/
│   ├── 08_GO_enrichment.R
│   ├── 09_structural_feature_analysis.R
│   └── 10_subcellular_localization.R
└── 05_domain-enrichment-mulea/
    ├── 11_mulea_domain_enrichment.R
    └── 12_mulea_visualization.R
```

---

## Scripts

### `01_per-tissue-analysis/`

**`01_abundance_vs_ratio_plots.R`** — Per-tissue scatter plots of total protein abundance fold change (x) vs. insolubility ratio fold change (y). Significant proteins are colored by tissue; non-significant proteins shown as gray background. Generates individual tissue plots and a multi-tissue overlay. Also exports standardized processed CSVs used by downstream scripts.

**`02_soluble_vs_insoluble_plots.R`** — Per-tissue and combined scatter plots comparing soluble fraction fold change vs. insoluble fraction fold change directly. Annotates top proteins by combined fold-change magnitude. Useful for distinguishing proteins that change in both fractions vs. only one.

**`03_longitudinal_comparison.R`** — Integrates nutritional intervention data with a parallel longitudinal aging dataset to assess concordance. Identifies proteins where the nutritional intervention recapitulates or opposes age-related changes in insolubility.

### `02_cross-tissue-integration/`

**`04_network_visualization.R`** — Builds bipartite protein–tissue networks using `ggraph`/`tidygraph`. Gene nodes connect to tissue nodes when a protein shows significant insolubility or abundance changes in that tissue. Inter-tissue proteins (significant in ≥2 tissues) are highlighted. Includes a third network variant overlaying mRNA expression levels (TPM) to identify proteins with discordant transcriptional vs. post-translational regulation. Exports node/edge CSVs and inter-tissue gene lists.

**`05_upset_plots.R`** — UpSet plots (via `ComplexUpset`) showing how significant proteins are distributed across tissues for both insolubility and abundance, separated by direction (increased vs. decreased under treatment). Includes bidirectional protein analysis (proteins increased in some tissues, decreased in others) and a comprehensive intersection summary exported to Excel.

### `03_biophysical-correlations/`

**`06_biophysical_correlations.R`** — Correlates insolubility ratio fold change with two biophysical properties: protein half-life (available for three tissues) and thermal stability / melting temperature (available for one tissue). Generates per-tissue scatter plots in four variants (labeled/unlabeled × with/without Pearson R²) and a multi-tissue overlay for half-life. Pearson correlations are reported with 95% confidence intervals.

**`07_fraction_comparison_plots.R`** — Two complementary scatter plot types per tissue: (A) soluble FC vs. insolubility ratio FC, and (B) soluble FC vs. insoluble FC. Annotates top 20 proteins by combined fold-change magnitude. Also exports standardized `Treatment_<Tissue>_processed_all_Treatment.csv` files consumed by Scripts 04, 06, 08–12.

### `04_enrichment-analysis/`

**`08_GO_enrichment.R`** — Gene Ontology overrepresentation analysis (ORA) via `clusterProfiler::enrichGO`. Runs BP, MF, and CC ontologies per tissue for proteins increased and decreased under treatment. Outputs dotplots (GeneRatio and fold enrichment), barplots, GO network maps (via `enrichplot`), per-tissue Excel workbooks, cross-tissue comparison heatmaps (-log10 padj), and a unified text summary. Mouse organism database (`org.Mm.eg.db`); BH FDR correction; padj < 0.05.

**`09_structural_feature_analysis.R`** — Tests whether eight biophysical/structural protein features (sequence length, charged AA%, α-helix%, β-sheet%, coiled-coil%, disordered%, phase separation score, amyloid prediction) differ between treatment-responsive and unchanged proteins. Two approaches: (A) Mann-Whitney U enrichment test comparing treatment-affected vs. unchanged proteins; (B) Pearson correlation of each feature against Log2FC across three protein subsets. Visualized as heatmaps matching a published figure style (correlation r, Z-score, significant-protein enrichment, cross-feature correlation). Reference data from Molzahn et al. *PNAS*.

**`10_subcellular_localization.R`** — Annotates treatment-responsive proteins with GO Cellular Component terms via `org.Mm.eg.db`, classifying each protein into 10 major compartments. Generates violin plots, beeswarm plots, enrichment bar charts, percentage heatmaps, bubble plots, and per-tissue volcano plots colored by compartment. Chi-square tests compare compartment distributions between aggregation-protective and aggregation-promoting proteins. Includes an interactive HTML explorer (plotly).

### `05_domain-enrichment-mulea/`

**`11_mulea_domain_enrichment.R`** — Protein domain enrichment using `mulea` ORA with empirical FDR (permutation-based, 10,000 permutations), which is more reliable than standard hypergeometric p-values for proteomics gene sets. Runs against two ontologies: (A) InterPro domains (GMT file) and (B) GO Molecular Function (via `muleaData`). Three comparisons per tissue: all treatment-affected proteins, treatment-protective only, treatment-promotes only. Known limitation: eFDR = 0 for very strong enrichments is replaced with 1e-10 for downstream log-transformation. Results exported to Excel.

**`12_mulea_visualization.R`** — Comprehensive visualization of mulea results. Cross-tissue: enrichment count overview, domain presence heatmap, -log10(eFDR) significance heatmap, bubble plot, domain type distribution, protective vs. promotes comparison, mirrored bar chart, and tissue-sharing summary. Per-tissue: 4-panel comprehensive PDF per tissue. Direction consistency analysis: identifies "switching domains" (enriched in opposite directions across tissues, indicating tissue-specific functional roles) vs. consistently directional domains; exports direction matrix CSV.

---

## Dependencies

```r
# CRAN
install.packages(c(
  "tidyverse", "ggrepel", "ggpmisc", "pheatmap", "patchwork",
  "scales", "writexl", "openxlsx", "ComplexUpset", "ggbeeswarm",
  "plotly", "htmlwidgets", "mulea", "corrplot"
))

# Bioconductor
BiocManager::install(c(
  "clusterProfiler", "org.Mm.eg.db", "GO.db", "enrichplot",
  "ggraph", "tidygraph", "ggnewscale",
  "muleaData", "ExperimentHub", "fgsea"
))
```

**R version:** Developed under R 4.3+. Package versions not pinned; minor adjustments may be needed for older or newer versions.

---

## Data Requirements

Each script expects input files in the working directory. The primary inputs are:

| File pattern | Generated by | Used by |
|---|---|---|
| `<tissue>_clean_data.csv` | External preprocessing | Scripts 07, 09 |
| `Treatment_<Tissue>_processed_all_Treatment.csv` | Script 07 | Scripts 04, 06, 08–12 |
| `mRNA_levels.xlsx` | External | Script 04 |
| `protein_half-life.xlsx` | External | Script 06 |
| `meltTemperature.xlsx` | External | Script 06 |
| `reference_structural_features.csv` | External (Molzahn et al.) | Script 09 |
| `mouse_interpro_domains.gmt.txt` | External (InterPro) | Script 11 |
| `domain_annotations_cache.csv` | Optional | Script 11 |
| `mulea_InterPro_results.xlsx` | Script 11 | Script 12 |

Scripts are designed to be run sequentially (01 → 12), though each folder is largely self-contained once the processed CSVs from Script 07 are available.

---

## Key Analytical Decisions

- **No fold-change threshold** was applied in any enrichment analysis; significance is determined by p-value alone (p < 0.05). This is consistent across GO, structural feature, subcellular localization, and domain enrichment analyses.
- **Direction convention:** Positive Log2FC = increased under treatment (promotes aggregation); negative Log2FC = decreased under treatment (protective against aggregation).
- **Background for ORA:** The background gene set for all enrichment analyses is the full set of proteins detected in that tissue (not the entire proteome), consistent with standard practice for proteomics ORA.
- **mulea empirical FDR** is used instead of standard BH FDR for domain enrichment because small gene sets in proteomics experiments can produce unreliable hypergeometric p-values.

---

## Citation

Manuscript in preparation. Repository will be updated with full citation upon publication.

---

## Author

**Chia-Lung Chuang**  
Postdoctoral Research Associate  
Department of Developmental Neurobiology  
St. Jude Children's Research Hospital  
Memphis, TN

PI: Dr. Fabio Demontis
