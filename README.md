# Mitochondrial Adaptive Failure in Hepatic Insulin Resistance
### Transcriptomic Analysis of GSE245301 — Mouse Liver RNA-seq

## Overview

A complete R bioinformatics pipeline analysing mitochondrial dysfunction in diet-induced hepatic insulin resistance, applied to publicly available mouse liver RNA-seq data.

**Dataset:** [GSE245301](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245301)

| Group | Description |
|-------|-------------|
| **Control** | Lean chow-fed mice |
| **OB-IS** | Obese, insulin-sensitive (HFD 4 weeks) |
| **OB-IR** | Obese, insulin-resistant (HFD 12 weeks) |

**Core finding:** The transition to insulin resistance drives suppression of oxidative phosphorylation (NES = -2.80), fatty acid metabolism, and ROS pathway, alongside activation of TNF-a and IL-6/JAK-STAT3 inflammatory signalling.

## Repository Structure

- 1_Scripts/ — R scripts, main pipeline and validation
- 2_Data/ — Raw and processed count matrices, DEG tables
- 3_Figures/ — All plots numbered in analysis order
- 4_Results/ — Result tables, DEGs, GSEA, WGCNA, network
- 5_Validation/ — Cross-species validation GSE132800

## Key Results

| Pathway | NES | padj |
|---------|-----|------|
| OXIDATIVE_PHOSPHORYLATION | -2.80 | 8.8e-29 |
| FATTY_ACID_METABOLISM | -2.38 | 6.6e-12 |
| REACTIVE_OXYGEN_SPECIES | -1.60 | 0.038 |
| MTORC1_SIGNALING | -1.40 | 0.044 |
| TNFA_SIGNALING_VIA_NFKB | +1.51 | 0.026 |
| IL6_JAK_STAT3_SIGNALING | +1.89 | 0.004 |

WGCNA identified 91 co-expression modules across 12,796 genes. The blue module (ME1) shows the strongest negative correlation with insulin resistance.

## Dependencies

```r
BiocManager::install(c("GEOquery","DESeq2","org.Mm.eg.db","AnnotationDbi","clusterProfiler","enrichplot","EnhancedVolcano"))
install.packages(c("WGCNA","fgsea","msigdbr","STRINGdb","igraph","ggraph","matrixStats","tidyverse","pheatmap","RColorBrewer","ggrepel","ashr"))
```

## Data

- Primary: [GSE245301](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245301)
- Validation: [GSE132800](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132800)

*This is a learning project — feedback and suggestions are very welcome.*
