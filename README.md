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

**Core finding:** The transition to insulin resistance (OB-IS → OB-IR) drives suppression of oxidative phosphorylation (NES = -2.76) and fatty-acid metabolism, alongside activation of IL-6/JAK-STAT3 inflammatory signalling. WGCNA identifies a 1,449-gene "repair/adaptive" (turquoise) module that is active in OB-IS but collapses in OB-IR, and which contains the mitophagy/antioxidant genes *Map1lc3b*, *Gpx1*, *Gclc*, *Tomm7* and *Ppargc1b*.

## Repository Structure

- 1_Scripts/ — R scripts: `00_download_data.R` (public GEO download), `GSE245301_pipeline.R` (main pipeline) and validation scripts
- 2_Data/ — Raw count matrix, sample table, VST matrix, ranked GSEA list
- 3_Figures/ — Publication figures numbered in analysis order
- 4_Results/ — DEG, GO/KEGG, GSEA, STRING network and WGCNA result tables
- 5_Validation/ — Independent cross-species validation (GSE132800)

The pipeline runs from a scratch working directory and writes intermediate
outputs there; the curated tables and figures committed under `2_Data/`,
`3_Figures/` and `4_Results/` are the final copies.

## Key Results

GSEA Hallmark, OB-IS vs OB-IR (positive NES = higher in OB-IR):

| Pathway | NES | padj |
|---------|-----|------|
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | -2.76 | 1.7e-28 |
| HALLMARK_FATTY_ACID_METABOLISM | -2.33 | 6.5e-12 |
| HALLMARK_ADIPOGENESIS | -1.98 | 4.9e-07 |
| HALLMARK_PEROXISOME | -2.00 | 4.8e-05 |
| HALLMARK_IL6_JAK_STAT3_SIGNALING | +1.87 | 0.0018 |

WGCNA detected 27 co-expression modules across the 12,796 filtered genes.
The turquoise module is the "repair" module (OB-IR correlation r = -0.85,
padj = 0.003; OB-IS r = +0.62); hub genes include *C8g*, *Chchd2*, *Rmdn1*,
*Nudt16l1* and *Acat3*.

## Dependencies

```r
BiocManager::install(c("GEOquery","DESeq2","org.Mm.eg.db","AnnotationDbi","clusterProfiler","enrichplot","EnhancedVolcano"))
install.packages(c("WGCNA","fgsea","msigdbr","STRINGdb","igraph","ggraph","matrixStats","tidyverse","pheatmap","RColorBrewer","ggrepel","ashr"))
```

## Data

- Primary: [GSE245301](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245301)
- Validation: [GSE132800](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132800)

*This is a learning project — feedback and suggestions are very welcome.*


## Plain-English Explanation
See [EXPLAIN_PLAIN_ENGLISH.md](EXPLAIN_PLAIN_ENGLISH.md) for a non-technical summary of what was done, steps, and verification.
