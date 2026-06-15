# Mitochondrial Adaptive Failure in Hepatic Insulin Resistance
### Transcriptomic Analysis of GSE245301 (Mouse Liver RNA-seq)

---

## Overview

This repository contains the complete bioinformatics pipeline and results for a transcriptomic study of mitochondrial dysfunction in diet-induced hepatic insulin resistance in mice.

**Dataset:** [GSE245301](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245301) — Mouse liver RNA-seq across three metabolic states:

| Group | Description |
|-------|-------------|
| **Control** | Lean, chow-fed mice (normal diet) |
| **OB-IS** | Obese, insulin-sensitive (HFD 4 weeks) |
| **OB-IR** | Obese, insulin-resistant (HFD 12 weeks) |

**Core finding:** Progressive suppression of oxidative phosphorylation, fatty acid metabolism, and mitochondrial stress response pathways accompanies the transition from insulin sensitivity to insulin resistance in obese liver.

---

## Repository Structure

```
mito/
├── Scripts/
│   └── GSE245301_pipeline.R          # Complete analysis pipeline (Stages 1–8)
│
├── Data/                              # Processed data outputs
│   ├── GSE245301_sample_table.csv    # Sample metadata and group assignments
│   ├── GSE245301_raw_count_matrix.csv
│   ├── GSE245301_VST_normalised_matrix.csv
│   ├── GSE245301_DEG_CompA_Control_vs_OB_IR.csv
│   ├── GSE245301_DEG_CompB_OB_IS_vs_OB_IR.csv
│   ├── Top_DEGs_OBIS_vs_OBIR.csv
│   ├── Top_DEGs_Annotated.csv
│   └── Top50_DEGs_OBIS_vs_OBIR.csv
│
├── Figures/                           # Main publication figures
│   ├── GSE245301_Figure1_PCA.png     # PCA of VST-normalised counts
│   ├── WGCNA_ModuleTraitHeatmap.png  # Module–trait correlation heatmap
│   └── FLOW.svg                       # Analysis workflow diagram
│
├── Publication_Figures/               # Final publication-quality figures
│   ├── Publication_Volcano_Plot.png  # Volcano plot — OB-IS vs OB-IR
│   ├── Publication_Heatmap.png       # Top 30 DEG heatmap
│   ├── Publication_Heatmap_FIXED.png # Corrected version
│   └── Top30_Heatmap_Genes.csv
│
├── GSEA_backup/                       # Gene Set Enrichment Analysis results
│   ├── GSEA_Figure5_Hallmark_Barplot.png
│   ├── GSEA_Figure6_KEGG_Barplot.png
│   ├── GSEA_Hallmark_Significant.csv
│   ├── GSEA_KEGG_Significant.csv
│   ├── GSEA_Key_Story_Pathways.csv
│   ├── GSEA_Mitochondrial_Stress_Pathways.csv
│   ├── HALLMARK_OXIDATIVE_PHOSPHORYLATION_Enrichment.png
│   ├── HALLMARK_FATTY_ACID_METABOLISM_Enrichment.png
│   ├── HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY_Enrichment.png
│   ├── HALLMARK_MTORC1_SIGNALING_Enrichment.png
│   ├── HALLMARK_TNFA_SIGNALING_VIA_NFKB_Enrichment.png
│   └── HALLMARK_UNFOLDED_PROTEIN_RESPONSE_Enrichment.png
│
├── WGCNA_backup/                      # Weighted Gene Co-expression Network Analysis
│   ├── ModuleTraitHeatmap_FINAL.png
│   ├── Eigengene_Trajectory_FINAL.png
│   ├── Top50_HubGenes_BlueModule.csv
│   ├── ME1_gene_list_FINAL.csv
│   └── ModuleTrait_Summary_FINAL.csv
│
├── Network_Analysis/                  # STRING protein interaction network
│   ├── Hub_Gene_Network_FIXED.png
│   ├── Hub_Genes.csv
│   └── STRING_Interactions.csv
│
├── Enrichment/Functional_Enrichment/  # GO and KEGG enrichment
│   ├── GO_Dotplot.png
│   ├── KEGG_Dotplot.png
│   ├── GO_Biological_Process.csv
│   └── KEGG_Pathways.csv
│
├── Relaxed_Enrichment/                # Relaxed threshold enrichment (padj < 0.1)
│   ├── Relaxed_GO_Dotplot.png
│   ├── Relaxed_KEGG_Dotplot.png
│   ├── Relaxed_GO_Enrichment.csv
│   └── Relaxed_KEGG_Enrichment.csv
│
├── GSE132800_validation/              # Cross-species validation (human FPKM data)
│   ├── signature_validation_plot.png
│   ├── cross_species_heatmap.png
│   ├── cross_species_translational_matrix.csv
│   ├── validation_metrics_summary.csv
│   └── *.R                            # Validation scripts
│
└── Drug_Discovery_NEK7/
    └── 6NPY.cif.gz                    # NEK7 crystal structure (PDB: 6NPY)
```

---

## Analysis Pipeline

The full pipeline is in `Scripts/GSE245301_pipeline.R` and runs end-to-end in 8 stages:

| Stage | Description |
|-------|-------------|
| **1–4** | GEO download → count matrix → DESeq2 differential expression |
| **5** | Ensembl ID → gene symbol annotation (`org.Mm.eg.db`) |
| **6** | STRING hub gene network analysis |
| **7** | GO and KEGG functional enrichment (`clusterProfiler`) |
| **8** | Volcano plot, expression heatmap, GSEA (Hallmark, KEGG, Reactome) + WGCNA |

**Before running:** set `outdir <- "."` at the top of the script (or the path to your working directory containing the `Data/` folder).

---

## Key Results

### Differential Expression (OB-IS vs OB-IR)
- Hundreds of significantly differentially expressed genes (padj < 0.05, |log2FC| > 1)
- Strong downregulation of mitochondrial and metabolic gene programmes in OB-IR

### Top GSEA Hallmark Pathways

| Pathway | NES | padj |
|---------|-----|------|
| OXIDATIVE_PHOSPHORYLATION | −2.80 | 8.8 × 10⁻²⁹ |
| FATTY_ACID_METABOLISM | −2.38 | 6.6 × 10⁻¹² |
| REACTIVE_OXYGEN_SPECIES_PATHWAY | −1.60 | 0.038 |
| MTORC1_SIGNALING | −1.40 | 0.044 |
| TNFA_SIGNALING_VIA_NFKB | +1.51 | 0.026 |
| IL6_JAK_STAT3_SIGNALING | +1.89 | 0.004 |

Negative NES = suppressed in OB-IR; Positive NES = activated in OB-IR.

### WGCNA
- 91 co-expression modules identified across 12,796 genes
- The blue module (ME1) shows the strongest negative correlation with OB-IR and contains hub genes involved in mitochondrial electron transport and fatty acid oxidation

---

## Dependencies

```r
# Core
library(GEOquery)
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)

# Annotation
library(org.Mm.eg.db)
library(AnnotationDbi)

# Enrichment & GSEA
library(clusterProfiler)
library(enrichplot)
library(fgsea)
library(msigdbr)

# Network
library(STRINGdb)
library(igraph)
library(ggraph)

# WGCNA
library(WGCNA)
library(matrixStats)
```

Install all packages:
```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "DESeq2", "org.Mm.eg.db", "AnnotationDbi",
                        "clusterProfiler", "enrichplot", "EnhancedVolcano"))
install.packages(c("WGCNA", "fgsea", "msigdbr", "STRINGdb",
                   "igraph", "ggraph", "matrixStats", "tidyverse",
                   "pheatmap", "RColorBrewer", "ggrepel", "ashr"))
```

---

## Data Availability

Raw data is publicly available on NCBI GEO:
- Primary dataset: [GSE245301](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245301)
- Validation dataset: [GSE132800](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132800)

---

## Citation

If you use this pipeline or results, please cite the original GEO dataset (GSE245301) and this repository.

---

## Author

Analysis pipeline developed for the study of mitochondrial adaptive failure in diet-induced hepatic insulin resistance.
