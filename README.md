# Mitochondrial Adaptive Failure in Hepatic Insulin Resistance
<<<<<<< HEAD
### Transcriptomic Analysis of GSE245301 — Mouse Liver RNA-seq
=======
### Transcriptomic Analysis of GSE245301 (Mouse Liver RNA-seq)
>>>>>>> deedabbafd3bf0a39190f88dde71384fe6862c5c

---

## Overview

<<<<<<< HEAD
A complete R-based bioinformatics pipeline analysing mitochondrial dysfunction in diet-induced hepatic insulin resistance, applied to publicly available mouse liver RNA-seq data.

**Dataset:** [GSE245301](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245301)

| Group | Description |
|-------|-------------|
| **Control** | Lean chow-fed mice |
| **OB-IS** | Obese, insulin-sensitive (HFD 4 weeks) |
| **OB-IR** | Obese, insulin-resistant (HFD 12 weeks) |

**Core finding:** The transition from insulin-sensitive to insulin-resistant obesity drives progressive suppression of oxidative phosphorylation (NES = −2.80, padj = 8.8×10⁻²⁹), fatty acid metabolism, and the ROS pathway — alongside activation of TNF-α/NF-κB and IL-6/JAK-STAT3 inflammatory signalling.
=======
This repository contains the complete bioinformatics pipeline and results for a transcriptomic study of mitochondrial dysfunction in diet-induced hepatic insulin resistance in mice.

**Dataset:** [GSE245301](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245301) — Mouse liver RNA-seq across three metabolic states:

| Group | Description |
|-------|-------------|
| **Control** | Lean, chow-fed mice (normal diet) |
| **OB-IS** | Obese, insulin-sensitive (HFD 4 weeks) |
| **OB-IR** | Obese, insulin-resistant (HFD 12 weeks) |

**Core finding:** Progressive suppression of oxidative phosphorylation, fatty acid metabolism, and mitochondrial stress response pathways accompanies the transition from insulin sensitivity to insulin resistance in obese liver.
>>>>>>> deedabbafd3bf0a39190f88dde71384fe6862c5c

---

## Repository Structure

```
<<<<<<< HEAD
├── 1_Scripts/          R scripts — main pipeline + validation
├── 2_Data/             Raw and processed count matrices, DEG tables
├── 3_Figures/          All plots (numbered in analysis order)
├── 4_Results/          Result tables — DEGs, GSEA, WGCNA, network
├── 5_Validation/       Cross-species validation (GSE132800)
└── README.md
=======
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
>>>>>>> deedabbafd3bf0a39190f88dde71384fe6862c5c
```

---

<<<<<<< HEAD
## Figures

| File | Description |
|------|-------------|
| `Fig0_Analysis_Workflow.svg` | Full pipeline overview |
| `Fig1_PCA.png` | PCA of VST-normalised counts |
| `Fig2_Volcano_Plot.png` | Volcano plot — OB-IS vs OB-IR |
| `Fig3_DEG_Heatmap.png` | Top 30 differentially expressed genes |
| `Fig4_GSEA_Hallmark_Barplot.png` | GSEA Hallmark pathways barplot |
| `Fig5_GSEA_KEGG_Barplot.png` | GSEA KEGG pathways barplot |
| `Fig6_GSEA_OxPhos_Enrichment.png` | Oxidative phosphorylation enrichment plot |
| `Fig7_GSEA_FattyAcid_Enrichment.png` | Fatty acid metabolism enrichment plot |
| `Fig8_GSEA_ROS_Enrichment.png` | Reactive oxygen species pathway enrichment |
| `Fig9_GSEA_TNFa_Enrichment.png` | TNF-α/NF-κB signalling enrichment |
| `Fig10_GSEA_mTORC1_Enrichment.png` | mTORC1 signalling enrichment |
| `Fig11_GSEA_UPR_Enrichment.png` | Unfolded protein response enrichment |
| `Fig12_GO_Dotplot.png` | GO Biological Process dotplot |
| `Fig13_KEGG_Dotplot.png` | KEGG pathway dotplot |
| `Fig14_KEGG_GSEA_Dotplot.png` | KEGG GSEA dotplot |
| `Fig15_KEGG_GSEA_Ridgeplot.png` | KEGG GSEA ridgeplot |
| `Fig16_WGCNA_Eigengene_Trajectory.png` | Module eigengene trajectory across groups |
| `Fig17_WGCNA_ModuleTrait_Heatmap.png` | WGCNA module–trait correlation heatmap |
| `Fig18_Hub_Gene_Network.png` | STRING hub gene protein interaction network |

---

## Pipeline Stages

| Stage | Description |
|-------|-------------|
| 1–4 | GEO download → count matrix → DESeq2 differential expression |
| 5 | Ensembl ID → gene symbol annotation (`org.Mm.eg.db`) |
| 6 | STRING hub gene network analysis |
| 7 | GO and KEGG functional enrichment (`clusterProfiler`) |
| 8 | Volcano plot, heatmap, GSEA (Hallmark/KEGG/Reactome), WGCNA |

**To run:** set `outdir <- "."` at the top of `1_Scripts/GSE245301_pipeline.R`.
=======
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
>>>>>>> deedabbafd3bf0a39190f88dde71384fe6862c5c

---

## Key Results

<<<<<<< HEAD
=======
### Differential Expression (OB-IS vs OB-IR)
- Hundreds of significantly differentially expressed genes (padj < 0.05, |log2FC| > 1)
- Strong downregulation of mitochondrial and metabolic gene programmes in OB-IR

### Top GSEA Hallmark Pathways

>>>>>>> deedabbafd3bf0a39190f88dde71384fe6862c5c
| Pathway | NES | padj |
|---------|-----|------|
| OXIDATIVE_PHOSPHORYLATION | −2.80 | 8.8 × 10⁻²⁹ |
| FATTY_ACID_METABOLISM | −2.38 | 6.6 × 10⁻¹² |
<<<<<<< HEAD
| REACTIVE_OXYGEN_SPECIES | −1.60 | 0.038 |
=======
| REACTIVE_OXYGEN_SPECIES_PATHWAY | −1.60 | 0.038 |
>>>>>>> deedabbafd3bf0a39190f88dde71384fe6862c5c
| MTORC1_SIGNALING | −1.40 | 0.044 |
| TNFA_SIGNALING_VIA_NFKB | +1.51 | 0.026 |
| IL6_JAK_STAT3_SIGNALING | +1.89 | 0.004 |

<<<<<<< HEAD
Negative NES = suppressed in OB-IR. Positive NES = activated in OB-IR.

WGCNA identified 91 co-expression modules across 12,796 genes. The blue module (ME1) shows the strongest negative correlation with insulin resistance and contains mitochondrial electron transport chain hub genes.
=======
Negative NES = suppressed in OB-IR; Positive NES = activated in OB-IR.

### WGCNA
- 91 co-expression modules identified across 12,796 genes
- The blue module (ME1) shows the strongest negative correlation with OB-IR and contains hub genes involved in mitochondrial electron transport and fatty acid oxidation
>>>>>>> deedabbafd3bf0a39190f88dde71384fe6862c5c

---

## Dependencies

```r
<<<<<<< HEAD
BiocManager::install(c("GEOquery", "DESeq2", "org.Mm.eg.db", "AnnotationDbi",
                        "clusterProfiler", "enrichplot", "EnhancedVolcano"))

install.packages(c("WGCNA", "fgsea", "msigdbr", "STRINGdb", "igraph",
                   "ggraph", "matrixStats", "tidyverse", "pheatmap",
                   "RColorBrewer", "ggrepel", "ashr"))
=======
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
>>>>>>> deedabbafd3bf0a39190f88dde71384fe6862c5c
```

---

## Data Availability

<<<<<<< HEAD
- Primary dataset: [GSE245301](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245301)
- Validation dataset: [GSE132800](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132800)
=======
Raw data is publicly available on NCBI GEO:
- Primary dataset: [GSE245301](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245301)
- Validation dataset: [GSE132800](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132800)

---

## Citation

If you use this pipeline or results, please cite the original GEO dataset (GSE245301) and this repository.

---

## Author

Analysis pipeline developed for the study of mitochondrial adaptive failure in diet-induced hepatic insulin resistance.
>>>>>>> deedabbafd3bf0a39190f88dde71384fe6862c5c
