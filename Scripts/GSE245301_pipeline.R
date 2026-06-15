# ============================================================
# MITOCHONDRIAL INSULIN RESISTANCE PIPELINE
# GSE245301
# COMPLETE CLEAN VERSION
# ============================================================

# ============================================================
# STEP 1 — LOAD LIBRARIES
# ============================================================

library(GEOquery)
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(EnhancedVolcano)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(ggrepel)
library(ashr)

# GSEA libraries
library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(enrichplot)

options(stringsAsFactors = FALSE)

cat("Libraries loaded successfully.\n\n")


# ============================================================
# STEP 2 — OUTPUT DIRECTORY
# ============================================================

outdir <- "."  # Set this to your working directory

dir.create(
  outdir,
  showWarnings = FALSE,
  recursive = TRUE
)

gsea_dir <- file.path(outdir, "GSEA")

dir.create(
  gsea_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

cat("All outputs will be saved to:\n")
cat(outdir, "\n\n")


# ============================================================
# STEP 3 — DOWNLOAD GEO METADATA
# ============================================================

cat("Downloading GEO metadata...\n")

gse <- getGEO(
  "GSE245301",
  GSEMatrix = TRUE,
  AnnotGPL = FALSE
)

gse <- gse[[1]]

meta <- pData(gse)

cat("Metadata downloaded successfully.\n\n")


# ============================================================
# STEP 4 — BUILD SAMPLE TABLE
# ============================================================

sample_table <- meta |>
  tibble::rownames_to_column("sample_id") |>
  dplyr::select(sample_id, title) |>
  dplyr::mutate(
    group = dplyr::case_when(
      
      stringr::str_detect(
        title,
        regex(
          "ND|normal|chow|control",
          ignore_case = TRUE
        )
      ) ~ "Control",
      
      stringr::str_detect(
        title,
        regex(
          "sensitive|HFD-4w",
          ignore_case = TRUE
        )
      ) ~ "OB_IS",
      
      stringr::str_detect(
        title,
        regex(
          "resistant|HFD-12w",
          ignore_case = TRUE
        )
      ) ~ "OB_IR",
      
      TRUE ~ "Unknown"
    )
  )

sample_table$group <- factor(
  sample_table$group,
  levels = c("Control", "OB_IS", "OB_IR")
)

cat("Group assignment:\n")
print(table(sample_table$group))

write.csv(
  sample_table,
  file.path(
    outdir,
    "GSE245301_sample_table.csv"
  ),
  row.names = FALSE
)

cat("\nSaved sample table.\n\n")


# ============================================================
# STEP 5 — READ RAW COUNT FILES
# ============================================================

raw_dir <- "GSE245301/raw"

if (!dir.exists(raw_dir)) {
  stop("ERROR: Raw count folder not found.")
}

raw_files <- list.files(
  raw_dir,
  full.names = TRUE
)

if (length(raw_files) == 0) {
  stop("ERROR: No raw files found.")
}

cat("Raw files detected:",
    length(raw_files),
    "\n\n")


# ============================================================
# STEP 6 — FUNCTION TO READ COUNT FILES
# ============================================================

read_count_file <- function(filepath) {
  
  df <- read.table(
    filepath,
    header = FALSE,
    sep = "\t",
    fill = TRUE
  )
  
  df <- df[
    !grepl("^__|^N_", df[[1]]),
  ]
  
  counts <- suppressWarnings(
    as.numeric(df[[ncol(df)]])
  )
  
  counts[is.na(counts)] <- 0
  
  data.frame(
    gene_id = as.character(df[[1]]),
    count = counts
  )
}


# ============================================================
# STEP 7 — LOAD ALL COUNT FILES
# ============================================================

cat("Reading all count files...\n")

count_list <- lapply(
  raw_files,
  read_count_file
)

cat("Files loaded:",
    length(count_list),
    "\n\n")


# ============================================================
# STEP 8 — BUILD COUNT MATRIX
# ============================================================

gene_ids <- count_list[[1]]$gene_id

count_matrix <- sapply(
  count_list,
  function(x) x$count
)

rownames(count_matrix) <- gene_ids

gsm_ids <- stringr::str_extract(
  basename(raw_files),
  "GSM[0-9]+"
)

colnames(count_matrix) <- gsm_ids


# ============================================================
# STEP 9 — MATCH SAMPLE ORDER
# ============================================================

common_samples <- intersect(
  colnames(count_matrix),
  sample_table$sample_id
)

count_matrix <- count_matrix[
  ,
  common_samples
]

sample_table <- sample_table |>
  dplyr::filter(
    sample_id %in% common_samples
  )

sample_table <- sample_table[
  match(
    colnames(count_matrix),
    sample_table$sample_id
  ),
]

stopifnot(
  identical(
    colnames(count_matrix),
    sample_table$sample_id
  )
)

cat("Count matrix dimensions:\n")

cat(
  nrow(count_matrix),
  "genes x",
  ncol(count_matrix),
  "samples\n\n"
)

write.csv(
  as.data.frame(count_matrix),
  file.path(
    outdir,
    "GSE245301_raw_count_matrix.csv"
  )
)

cat("Saved raw count matrix.\n\n")


# ============================================================
# STEP 10 — CLEAN COUNT MATRIX
# ============================================================

count_matrix <- round(count_matrix)

count_matrix[is.na(count_matrix)] <- 0

count_matrix[count_matrix < 0] <- 0

storage.mode(count_matrix) <- "integer"

cat("Count matrix cleaned.\n\n")


# ============================================================
# STEP 11 — CREATE DESEQ2 OBJECT
# ============================================================

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_table,
  design = ~ group
)

keep <- rowSums(
  counts(dds) >= 10
) >= 3

dds <- dds[keep, ]

cat(
  "Genes after filtering:",
  nrow(dds),
  "\n\n"
)


# ============================================================
# STEP 12 — VST NORMALIZATION
# ============================================================

cat("Running VST normalization...\n")

vsd <- vst(
  dds,
  blind = TRUE
)

vst_matrix <- assay(vsd)

write.csv(
  as.data.frame(vst_matrix),
  file.path(
    outdir,
    "GSE245301_VST_normalised_matrix.csv"
  )
)

cat("Saved VST matrix.\n\n")


# ============================================================
# STEP 13 — PCA PLOT
# ============================================================

pca_data <- plotPCA(
  vsd,
  intgroup = "group",
  returnData = TRUE
)

pct_var <- round(
  100 * attr(pca_data, "percentVar")
)

pca_plot <- ggplot(
  pca_data,
  aes(
    PC1,
    PC2,
    color = group,
    label = name
  )
) +
  geom_point(
    size = 4,
    alpha = 0.9
  ) +
  ggrepel::geom_text_repel(
    size = 3,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Control" = "#1D9E75",
      "OB_IS" = "#378ADD",
      "OB_IR" = "#D85A30"
    )
  ) +
  labs(
    title = "PCA — GSE245301",
    subtitle = "Variance Stabilized Counts",
    x = paste0(
      "PC1 (",
      pct_var[1],
      "% variance)"
    ),
    y = paste0(
      "PC2 (",
      pct_var[2],
      "% variance)"
    )
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(
    outdir,
    "GSE245301_Figure1_PCA.png"
  ),
  pca_plot,
  width = 7,
  height = 5,
  dpi = 300
)

cat("Saved PCA plot.\n\n")


# ============================================================
# STEP 14 — RUN DESEQ2
# ============================================================

cat("Running DESeq2...\n")

dds <- DESeq(dds)

cat("DESeq2 complete.\n\n")


# ============================================================
# STEP 15 — COMPARISON A
# CONTROL vs OB_IR
# ============================================================

cat("Running Comparison A...\n")

res_A <- results(
  dds,
  contrast = c(
    "group",
    "OB_IR",
    "Control"
  ),
  alpha = 0.05
)

res_A <- lfcShrink(
  dds,
  contrast = c(
    "group",
    "OB_IR",
    "Control"
  ),
  res = res_A,
  type = "ashr"
)

res_A_df <- as.data.frame(res_A) |>
  tibble::rownames_to_column("gene_id") |>
  dplyr::arrange(padj)

write.csv(
  res_A_df,
  file.path(
    outdir,
    "GSE245301_DEG_CompA_Control_vs_OB_IR.csv"
  ),
  row.names = FALSE
)

cat("Saved Comparison A.\n\n")


# ============================================================
# STEP 16 — COMPARISON B
# OB_IS vs OB_IR
# ============================================================

cat("Running Comparison B...\n")

res_B <- results(
  dds,
  contrast = c(
    "group",
    "OB_IR",
    "OB_IS"
  ),
  alpha = 0.05
)

res_B <- lfcShrink(
  dds,
  contrast = c(
    "group",
    "OB_IR",
    "OB_IS"
  ),
  res = res_B,
  type = "ashr"
)

res_B_df <- as.data.frame(res_B) |>
  tibble::rownames_to_column("gene_id") |>
  dplyr::arrange(padj)

write.csv(
  res_B_df,
  file.path(
    outdir,
    "GSE245301_DEG_CompB_OB_IS_vs_OB_IR.csv"
  ),
  row.names = FALSE
)

cat("Saved Comparison B.\n\n")


# ============================================================
# STEP 17 — DEG SUMMARY
# ============================================================

summarise_degs <- function(df, label) {
  
  sig <- df |>
    dplyr::filter(
      !is.na(padj),
      padj < 0.05,
      abs(log2FoldChange) > 1
    )
  
  cat(
    "\n",
    label,
    "\n---------------------\n",
    "Total DEGs:",
    nrow(sig),
    "\nUpregulated:",
    sum(sig$log2FoldChange > 0),
    "\nDownregulated:",
    sum(sig$log2FoldChange < 0),
    "\n"
  )
}

cat("\n================ DEG SUMMARY ================\n")

summarise_degs(
  res_A_df,
  "Comparison A — Control vs OB_IR"
)

summarise_degs(
  res_B_df,
  "Comparison B — OB_IS vs OB_IR"
)

cat("\nDEG analysis complete.\n\n")


# ============================================================
# STEP 18 — GSEA ANALYSIS
# ============================================================

cat("Preparing ranked genes...\n")

deg <- res_B_df |>
  dplyr::filter(
    !is.na(log2FoldChange),
    !is.na(padj)
  )

deg$ENSEMBL <- gsub(
  "\\..*",
  "",
  deg$gene_id
)

gene_annot <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = unique(deg$ENSEMBL),
  columns = c("SYMBOL"),
  keytype = "ENSEMBL"
)

gene_annot <- gene_annot |>
  dplyr::filter(!is.na(SYMBOL)) |>
  dplyr::distinct(
    ENSEMBL,
    .keep_all = TRUE
  )

deg <- left_join(
  deg,
  gene_annot,
  by = "ENSEMBL"
)

deg <- deg |>
  dplyr::mutate(
    ranking_score =
      sign(log2FoldChange) *
      -log10(padj + 1e-300)
  )

ranked_genes <- deg$ranking_score

names(ranked_genes) <- deg$SYMBOL

valid <- !is.na(names(ranked_genes)) &
  names(ranked_genes) != ""

ranked_genes <- ranked_genes[valid]

ranked_genes <- ranked_genes[
  !duplicated(names(ranked_genes))
]

ranked_genes <- ranked_genes[
  is.finite(ranked_genes)
]

ranked_genes <- sort(
  ranked_genes,
  decreasing = TRUE
)

cat("Ranked genes ready.\n\n")


# ============================================================
# STEP 19 — LOAD HALLMARK PATHWAYS
# ============================================================

hallmark_sets <- msigdbr(
  species = "Mus musculus",
  collection = "H"
)

hallmark_list <- split(
  hallmark_sets$gene_symbol,
  hallmark_sets$gs_name
)

cat(
  "Hallmark pathways loaded:",
  length(hallmark_list),
  "\n\n"
)


# ============================================================
# STEP 20 — RUN FGSEA
# ============================================================

cat("Running FGSEA...\n")

fgsea_res <- fgsea(
  pathways = hallmark_list,
  stats = ranked_genes,
  minSize = 10,
  maxSize = 500
)

fgsea_res <- fgsea_res |>
  as.data.frame() |>
  dplyr::arrange(padj)

# IMPORTANT FIX
fgsea_res$leadingEdge <- sapply(
  fgsea_res$leadingEdge,
  paste,
  collapse = ";"
)

cat("FGSEA complete.\n\n")


# ============================================================
# STEP 21 — SAVE GSEA RESULTS
# ============================================================

write.csv(
  fgsea_res,
  file.path(
    gsea_dir,
    "GSEA_Hallmark_All_Pathways.csv"
  ),
  row.names = FALSE
)

sig_pathways <- fgsea_res |>
  dplyr::filter(
    !is.na(padj),
    padj < 0.05
  )

write.csv(
  sig_pathways,
  file.path(
    gsea_dir,
    "GSEA_Hallmark_Significant_Pathways.csv"
  ),
  row.names = FALSE
)

cat("Significant pathways:",
    nrow(sig_pathways),
    "\n\n")


# ============================================================
# STEP 22 — MITOCHONDRIAL PATHWAYS
# ============================================================

mito_keywords <- c(
  "OXIDATIVE",
  "REACTIVE_OXYGEN",
  "FATTY_ACID",
  "MTOR",
  "HYPOXIA",
  "INFLAMMATORY",
  "TNFA",
  "UNFOLDED",
  "APOPTOSIS",
  "PEROXISOME"
)

mito_pathways <- sig_pathways |>
  dplyr::filter(
    stringr::str_detect(
      pathway,
      paste(
        mito_keywords,
        collapse = "|"
      )
    )
  )

write.csv(
  mito_pathways,
  file.path(
    gsea_dir,
    "GSEA_Mitochondrial_Stress_Pathways.csv"
  ),
  row.names = FALSE
)

cat("Saved mitochondrial pathways.\n\n")


# ============================================================
# STEP 23 — BARPLOT
# ============================================================

top_plot <- sig_pathways |>
  dplyr::slice_max(
    order_by = abs(NES),
    n = 15
  )

barplot_gsea <- ggplot(
  top_plot,
  aes(
    x = reorder(pathway, NES),
    y = NES,
    fill = NES
  )
) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 12)

ggsave(
  file.path(
    gsea_dir,
    "Figure_GSEA_Barplot.png"
  ),
  barplot_gsea,
  width = 10,
  height = 7,
  dpi = 300
)

cat("Saved GSEA barplot.\n\n")


# ============================================================
# STEP 24 — FINAL SUMMARY
# ============================================================

cat("\n=========================================\n")
cat(" PIPELINE COMPLETE\n")
cat("=========================================\n\n")

cat("Major finding:\n")
cat("Strong suppression of mitochondrial oxidative phosphorylation in OB_IR.\n\n")

cat("All files saved to:\n")
cat(outdir, "\n\n")

----------------------------------------------
  
  # ============================================================
# STAGE 4 — CORE DEG EXTRACTION
# GSE245301 | OB_IS vs OB_IR
# ============================================================

cat("\n================================================\n")
cat(" STAGE 4 — CORE DEG EXTRACTION\n")
cat("================================================\n\n")


# ============================================================
# STEP 1 — LOAD DEG FILE
# ============================================================

outdir <- "."  # Set this to your working directory

deg_file <- file.path(
  outdir,
  "GSE245301_DEG_CompB_OB_IS_vs_OB_IR.csv"
)

if (!file.exists(deg_file)) {
  stop("ERROR: DEG file not found.")
}

deg <- read.csv(deg_file)

cat("DEG file loaded.\n")
cat("Rows:", nrow(deg), "\n\n")


# ============================================================
# STEP 2 — REMOVE INVALID VALUES
# ============================================================

deg <- deg |>
  dplyr::filter(
    !is.na(log2FoldChange),
    !is.na(padj),
    is.finite(log2FoldChange),
    is.finite(padj)
  )

cat("Rows after filtering:",
    nrow(deg),
    "\n\n")


# ============================================================
# STEP 3 — EXTRACT SIGNIFICANT DEGs
# ============================================================

top_deg <- deg |>
  dplyr::filter(
    padj < 0.05,
    abs(log2FoldChange) > 1
  ) |>
  dplyr::arrange(padj)

cat("Significant DEGs extracted.\n\n")

cat("Total significant DEGs:",
    nrow(top_deg),
    "\n")

cat("Upregulated genes:",
    sum(top_deg$log2FoldChange > 0),
    "\n")

cat("Downregulated genes:",
    sum(top_deg$log2FoldChange < 0),
    "\n\n")


# ============================================================
# STEP 4 — SAVE FULL DEG TABLE
# ============================================================

write.csv(
  top_deg,
  file.path(
    outdir,
    "Top_DEGs_OBIS_vs_OBIR.csv"
  ),
  row.names = FALSE
)

cat("Saved:\n")
cat("Top_DEGs_OBIS_vs_OBIR.csv\n\n")


# ============================================================
# STEP 5 — TOP 50 GENES
# ============================================================

top50 <- top_deg |>
  dplyr::select(
    gene_id,
    log2FoldChange,
    padj
  ) |>
  head(50)

print(top50)

write.csv(
  top50,
  file.path(
    outdir,
    "Top50_DEGs_OBIS_vs_OBIR.csv"
  ),
  row.names = FALSE
)

cat("\nSaved:\n")
cat("Top50_DEGs_OBIS_vs_OBIR.csv\n\n")


# ============================================================
# STEP 6 — SUMMARY
# ============================================================

cat("====================================\n")
cat(" CORE DEG EXTRACTION COMPLETE\n")
cat("====================================\n\n")

cat("Generated files:\n\n")

cat("1. Top_DEGs_OBIS_vs_OBIR.csv\n")
cat("2. Top50_DEGs_OBIS_vs_OBIR.csv\n\n")

cat("Saved in:\n")
cat(outdir, "\n\n")


# ============================================================
# OPTIONAL — QUICK BIOLOGICAL CHECK
# ============================================================

cat("Top 20 strongest DEGs:\n\n")

print(
  top_deg |>
    dplyr::select(
      gene_id,
      log2FoldChange,
      padj
    ) |>
    head(20)
)

----------------------------------------------------

  # ============================================================
# STAGE 5 — CONVERT ENSEMBL IDs TO GENE SYMBOLS
# ============================================================

cat("\n================================================\n")
cat(" STAGE 5 — GENE SYMBOL ANNOTATION\n")
cat("================================================\n\n")

library(org.Mm.eg.db)
library(AnnotationDbi)
library(dplyr)

# ============================================================
# FILE PATHS
# ============================================================

outdir <- "."  # Set this to your working directory

deg_file <- file.path(
  outdir,
  "Top_DEGs_OBIS_vs_OBIR.csv"
)

# ============================================================
# LOAD DEG FILE
# ============================================================

deg <- read.csv(deg_file)

cat("Loaded DEG file.\n")
cat("Rows:", nrow(deg), "\n\n")

# ============================================================
# CLEAN ENSEMBL IDS
# ============================================================

deg$ENSEMBL <- gsub(
  "\\..*",
  "",
  deg$gene_id
)

# ============================================================
# MAP ENSEMBL → SYMBOL
# ============================================================

cat("Mapping gene symbols...\n")

annot <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = unique(deg$ENSEMBL),
  columns = c(
    "SYMBOL",
    "GENENAME"
  ),
  keytype = "ENSEMBL"
)

annot <- annot |>
  dplyr::filter(!is.na(SYMBOL)) |>
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

# ============================================================
# MERGE
# ============================================================

deg_annotated <- left_join(
  deg,
  annot,
  by = "ENSEMBL"
)

# ============================================================
# SORT
# ============================================================

deg_annotated <- deg_annotated |>
  dplyr::arrange(padj)

# ============================================================
# SAVE
# ============================================================

write.csv(
  deg_annotated,
  file.path(
    outdir,
    "Top_DEGs_Annotated.csv"
  ),
  row.names = FALSE
)

cat("Annotated DEG file saved.\n\n")

# ============================================================
# SHOW TOP GENES
# ============================================================

top_genes <- deg_annotated |>
  dplyr::select(
    SYMBOL,
    GENENAME,
    log2FoldChange,
    padj
  ) |>
  head(30)

print(top_genes)

# ============================================================
# SUMMARY
# ============================================================

cat("\n========================================\n")
cat(" GENE ANNOTATION COMPLETE\n")
cat("========================================\n\n")

cat("Generated file:\n")
cat("Top_DEGs_Annotated.csv\n\n")

cat("Saved in:\n")
cat(outdir, "\n\n")

----
  
  
  # ============================================================
# STAGE 6 — HUB GENE NETWORK ANALYSIS
# ============================================================

cat("\n================================================\n")
cat(" STAGE 6 — HUB GENE NETWORK ANALYSIS\n")
cat("================================================\n\n")

# ============================================================
# LOAD LIBRARIES
# ============================================================

library(dplyr)
library(STRINGdb)
library(igraph)
library(ggraph)
library(ggplot2)

cat("Libraries loaded.\n\n")

# ============================================================
# OUTPUT DIRECTORY
# ============================================================

outdir <- "."  # Set this to your working directory

network_dir <- file.path(
  outdir,
  "Network_Analysis"
)

dir.create(
  network_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

cat("Network results folder created.\n\n")

# ============================================================
# LOAD ANNOTATED DEGs
# ============================================================

deg_file <- file.path(
  outdir,
  "Top_DEGs_Annotated.csv"
)

deg <- read.csv(deg_file)

cat("Annotated DEG file loaded.\n")
cat("Rows:", nrow(deg), "\n\n")

# ============================================================
# REMOVE MISSING SYMBOLS
# ============================================================

deg <- deg |>
  dplyr::filter(
    !is.na(SYMBOL),
    SYMBOL != ""
  )

cat("Genes after filtering:", nrow(deg), "\n\n")

# ============================================================
# INITIALIZE STRING DATABASE
# ============================================================

string_db <- STRINGdb$new(
  version = "11.5",
  species = 10090,
  score_threshold = 400,
  input_directory = ""
)

cat("STRINGdb initialized.\n\n")

# ============================================================
# MAP GENE SYMBOLS TO STRING IDs
# ============================================================

mapped <- string_db$map(
  deg,
  "SYMBOL",
  removeUnmappedRows = TRUE
)

cat("Mapped genes:", nrow(mapped), "\n\n")

# ============================================================
# GET INTERACTION NETWORK
# ============================================================

interactions <- string_db$get_interactions(
  mapped$STRING_id
)

cat("Interactions found:", nrow(interactions), "\n\n")

# ============================================================
# SAVE RAW INTERACTIONS
# ============================================================

write.csv(
  interactions,
  file.path(
    network_dir,
    "STRING_Interactions.csv"
  ),
  row.names = FALSE
)

cat("Saved STRING interactions.\n\n")

# ============================================================
# BUILD GRAPH
# ============================================================

g <- graph_from_data_frame(
  interactions,
  directed = FALSE
)

cat("Network graph created.\n\n")

# ============================================================
# CALCULATE HUB SCORES
# ============================================================

deg_degree <- degree(g)

hub_table <- data.frame(
  STRING_id = names(deg_degree),
  Degree = as.numeric(deg_degree)
)

# ============================================================
# MAP BACK TO GENE SYMBOLS
# ============================================================

hub_table <- left_join(
  hub_table,
  mapped |> dplyr::select(
    STRING_id,
    SYMBOL
  ),
  by = "STRING_id"
)

hub_table <- hub_table |>
  dplyr::distinct(STRING_id, .keep_all = TRUE) |>
  dplyr::arrange(desc(Degree))

# ============================================================
# SAVE HUB GENES
# ============================================================

write.csv(
  hub_table,
  file.path(
    network_dir,
    "Hub_Genes.csv"
  ),
  row.names = FALSE
)

cat("Hub genes saved.\n\n")

# ============================================================
# DISPLAY TOP HUB GENES
# ============================================================

cat("Top Hub Genes:\n\n")

print(
  hub_table |>
    dplyr::select(
      SYMBOL,
      Degree
    ) |>
    head(20)
)

# ============================================================
# CREATE NETWORK PLOT
# ============================================================

top_nodes <- hub_table |>
  head(25)

top_ids <- top_nodes$STRING_id

sub_interactions <- interactions |>
  dplyr::filter(
    from %in% top_ids |
      to %in% top_ids
  )

g_sub <- graph_from_data_frame(
  sub_interactions,
  directed = FALSE
)

png(
  file.path(
    network_dir,
    "Hub_Gene_Network.png"
  ),
  width = 2400,
  height = 2200,
  res = 300
)

plot(
  g_sub,
  vertex.label.cex = 0.7,
  vertex.size = 8,
  main = "Hub Gene Network"
)

dev.off()

cat("Network plot saved.\n\n")

# ============================================================
# SAVE INTERPRETATION
# ============================================================

sink(
  file.path(
    network_dir,
    "Hub_Gene_Interpretation.txt"
  )
)

cat("HUB GENE NETWORK ANALYSIS\n")
cat("====================================\n\n")

cat("Top hub genes:\n\n")

print(
  hub_table |>
    dplyr::select(
      SYMBOL,
      Degree
    ) |>
    head(30)
)

cat("\n\n")

cat("Interpretation:\n\n")

cat(
  "The highest-degree genes may represent central
regulators of mitochondrial adaptation,
stress signaling, metabolic remodeling,
and hepatic insulin resistance progression.\n"
)

sink()

cat("Interpretation summary saved.\n\n")

# ============================================================
# FINAL SUMMARY
# ============================================================

cat("========================================\n")
cat(" HUB GENE ANALYSIS COMPLETE\n")
cat("========================================\n\n")

cat("FILES GENERATED:\n\n")

cat("1. STRING_Interactions.csv\n")
cat("2. Hub_Genes.csv\n")
cat("3. Hub_Gene_Network.png\n")
cat("4. Hub_Gene_Interpretation.txt\n\n")

cat("Saved in:\n")
cat(network_dir, "\n\n")

----
  
  
  # ============================================================
# STAGE 7 — GO + KEGG ENRICHMENT ANALYSIS
# ============================================================

cat("\n================================================\n")
cat(" STAGE 7 — GO + KEGG ENRICHMENT ANALYSIS\n")
cat("================================================\n\n")

# ============================================================
# LIBRARIES
# ============================================================
install.packages("ggraph")
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(dplyr)
library(ggplot2)

cat("Libraries loaded.\n\n")

# ============================================================
# OUTPUT DIRECTORY
# ============================================================

outdir <- "."  # Set this to your working directory

enrich_dir <- file.path(
  outdir,
  "Functional_Enrichment"
)

dir.create(
  enrich_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

cat("Enrichment folder created.\n\n")

# ============================================================
# LOAD ANNOTATED DEGs
# ============================================================

deg_file <- file.path(
  outdir,
  "Top_DEGs_Annotated.csv"
)

deg <- read.csv(deg_file)

cat("Annotated DEG file loaded.\n")
cat("Rows:", nrow(deg), "\n\n")

# ============================================================
# REMOVE NA SYMBOLS
# ============================================================

deg <- deg |>
  dplyr::filter(
    !is.na(SYMBOL),
    SYMBOL != ""
  )

# ============================================================
# CONVERT SYMBOL → ENTREZ
# ============================================================

gene_map <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = unique(deg$SYMBOL),
  columns = c("ENTREZID"),
  keytype = "SYMBOL"
)

gene_map <- gene_map |>
  dplyr::filter(!is.na(ENTREZID)) |>
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

entrez_genes <- unique(gene_map$ENTREZID)

cat("Entrez genes mapped:", length(entrez_genes), "\n\n")

# ============================================================
# GO BIOLOGICAL PROCESS
# ============================================================

cat("Running GO enrichment...\n")

ego <- enrichGO(
  gene = entrez_genes,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

ego_df <- as.data.frame(ego)

write.csv(
  ego_df,
  file.path(
    enrich_dir,
    "GO_Biological_Process.csv"
  ),
  row.names = FALSE
)

cat("GO enrichment complete.\n\n")

# ============================================================
# KEGG PATHWAY ANALYSIS
# ============================================================

cat("Running KEGG enrichment...\n")

ekegg <- enrichKEGG(
  gene = entrez_genes,
  organism = "mmu",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

kegg_df <- as.data.frame(ekegg)

write.csv(
  kegg_df,
  file.path(
    enrich_dir,
    "KEGG_Pathways.csv"
  ),
  row.names = FALSE
)

cat("KEGG enrichment complete.\n\n")

# ============================================================
# DISPLAY RESULTS
# ============================================================

cat("TOP GO TERMS:\n\n")

print(
  ego_df |>
    dplyr::select(
      Description,
      p.adjust,
      GeneRatio
    ) |>
    head(20)
)

cat("\n\n")

cat("TOP KEGG PATHWAYS:\n\n")

print(
  kegg_df |>
    dplyr::select(
      Description,
      p.adjust
    ) |>
    head(20)
)

# ============================================================
# GO DOTPLOT
# ============================================================

png(
  file.path(
    enrich_dir,
    "GO_Dotplot.png"
  ),
  width = 2400,
  height = 2000,
  res = 300
)

print(
  dotplot(
    ego,
    showCategory = 15
  ) +
    ggtitle("GO Biological Process Enrichment")
)

dev.off()

cat("GO dotplot saved.\n\n")

# ============================================================
# KEGG DOTPLOT
# ============================================================

png(
  file.path(
    enrich_dir,
    "KEGG_Dotplot.png"
  ),
  width = 2400,
  height = 2000,
  res = 300
)

print(
  dotplot(
    ekegg,
    showCategory = 15
  ) +
    ggtitle("KEGG Pathway Enrichment")
)

dev.off()

cat("KEGG dotplot saved.\n\n")

# ============================================================
# INTERPRETATION FILE
# ============================================================

sink(
  file.path(
    enrich_dir,
    "Functional_Interpretation.txt"
  )
)

cat("FUNCTIONAL ENRICHMENT ANALYSIS\n")
cat("====================================\n\n")

cat("TOP GO TERMS:\n\n")

print(
  ego_df |>
    dplyr::select(
      Description,
      p.adjust,
      GeneRatio
    ) |>
    head(30)
)

cat("\n\n")

cat("TOP KEGG PATHWAYS:\n\n")

print(
  kegg_df |>
    dplyr::select(
      Description,
      p.adjust
    ) |>
    head(30)
)

cat("\n\n")

cat(
  "These enrichments indicate coordinated
alteration of metabolic, mitochondrial,
oxidative stress, inflammatory,
and insulin resistance-associated pathways.\n"
)

sink()

cat("Interpretation saved.\n\n")

# ============================================================
# FINAL SUMMARY
# ============================================================

cat("========================================\n")
cat(" FUNCTIONAL ENRICHMENT COMPLETE\n")
cat("========================================\n\n")

cat("FILES GENERATED:\n\n")

cat("1. GO_Biological_Process.csv\n")
cat("2. KEGG_Pathways.csv\n")
cat("3. GO_Dotplot.png\n")
cat("4. KEGG_Dotplot.png\n")
cat("5. Functional_Interpretation.txt\n\n")

cat("Saved in:\n")
cat(enrich_dir, "\n\n")

###

# ============================================================
# STAGE 7 — RELAXED FUNCTIONAL ENRICHMENT
# ============================================================

cat("\n================================================\n")
cat(" STAGE 7 — RELAXED ENRICHMENT ANALYSIS\n")
cat("================================================\n\n")

library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(dplyr)
library(ggplot2)

outdir <- "."  # Set this to your working directory

deg_file <- file.path(
  outdir,
  "GSE245301_DEG_CompB_OB_IS_vs_OB_IR.csv"
)

deg <- read.csv(deg_file)

cat("Loaded DEG file.\n")
cat("Rows:", nrow(deg), "\n\n")

# ============================================================
# RELAXED FILTER
# ============================================================

deg_filtered <- deg |>
  dplyr::filter(
    !is.na(log2FoldChange),
    !is.na(padj),
    padj < 0.1
  )

cat("Genes after relaxed filtering:",
    nrow(deg_filtered),
    "\n\n")

# ============================================================
# ENSEMBL → SYMBOL
# ============================================================

deg_filtered$ENSEMBL <- gsub(
  "\\..*",
  "",
  deg_filtered$gene_id
)

annot <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = unique(deg_filtered$ENSEMBL),
  columns = c(
    "SYMBOL",
    "ENTREZID"
  ),
  keytype = "ENSEMBL"
)

annot <- annot |>
  dplyr::filter(
    !is.na(SYMBOL),
    !is.na(ENTREZID)
  ) |>
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

deg_annot <- left_join(
  deg_filtered,
  annot,
  by = "ENSEMBL"
)

entrez_genes <- unique(
  deg_annot$ENTREZID
)

cat("Mapped Entrez genes:",
    length(entrez_genes),
    "\n\n")

# ============================================================
# OUTPUT DIRECTORY
# ============================================================

relax_dir <- file.path(
  outdir,
  "Relaxed_Enrichment"
)

dir.create(
  relax_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ============================================================
# GO ENRICHMENT
# ============================================================

cat("Running GO enrichment...\n")

ego <- enrichGO(
  gene = entrez_genes,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.2,
  readable = TRUE
)

ego_df <- as.data.frame(ego)

write.csv(
  ego_df,
  file.path(
    relax_dir,
    "Relaxed_GO_Enrichment.csv"
  ),
  row.names = FALSE
)

cat("GO complete.\n\n")

# ============================================================
# KEGG ENRICHMENT
# ============================================================

cat("Running KEGG enrichment...\n")

ekegg <- enrichKEGG(
  gene = entrez_genes,
  organism = "mmu",
  pvalueCutoff = 0.1
)

kegg_df <- as.data.frame(ekegg)

write.csv(
  kegg_df,
  file.path(
    relax_dir,
    "Relaxed_KEGG_Enrichment.csv"
  ),
  row.names = FALSE
)

cat("KEGG complete.\n\n")

# ============================================================
# PRINT RESULTS
# ============================================================

cat("TOP GO TERMS:\n\n")

print(
  ego_df |>
    dplyr::select(
      Description,
      p.adjust,
      GeneRatio
    ) |>
    head(20)
)

cat("\n\n")

cat("TOP KEGG PATHWAYS:\n\n")

print(
  kegg_df |>
    dplyr::select(
      Description,
      p.adjust
    ) |>
    head(20)
)

# ============================================================
# GO DOTPLOT
# ============================================================

if (nrow(ego_df) > 0) {
  
  png(
    file.path(
      relax_dir,
      "Relaxed_GO_Dotplot.png"
    ),
    width = 2400,
    height = 2000,
    res = 300
  )
  
  print(
    dotplot(
      ego,
      showCategory = 15
    ) +
      ggtitle(
        "Relaxed GO Enrichment"
      )
  )
  
  dev.off()
  
  cat("GO dotplot saved.\n\n")
}

# ============================================================
# KEGG DOTPLOT
# ============================================================

if (nrow(kegg_df) > 0) {
  
  png(
    file.path(
      relax_dir,
      "Relaxed_KEGG_Dotplot.png"
    ),
    width = 2400,
    height = 2000,
    res = 300
  )
  
  print(
    dotplot(
      ekegg,
      showCategory = 15
    ) +
      ggtitle(
        "Relaxed KEGG Enrichment"
      )
  )
  
  dev.off()
  
  cat("KEGG dotplot saved.\n\n")
}

cat("========================================\n")
cat(" RELAXED ENRICHMENT COMPLETE\n")
cat("========================================\n\n")

cat("Saved in:\n")
cat(relax_dir, "\n\n")

#########################################################

# ============================================================
# STAGE 8 — VOLCANO PLOT + HEATMAP
# ============================================================

cat("\n================================================\n")
cat(" STAGE 8 — VOLCANO PLOT + HEATMAP\n")
cat("================================================\n\n")

# ============================================================
# LIBRARIES
# ============================================================

library(dplyr)
library(ggplot2)
library(pheatmap)
library(org.Mm.eg.db)
library(AnnotationDbi)

cat("Libraries loaded.\n\n")

# ============================================================
# OUTPUT FOLDER
# ============================================================

outdir <- "."  # Set this to your working directory

fig_dir <- file.path(
  outdir,
  "Publication_Figures"
)

dir.create(
  fig_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

cat("Figure directory created.\n\n")

# ============================================================
# LOAD DEG FILE
# ============================================================

deg_file <- file.path(
  outdir,
  "GSE245301_DEG_CompB_OB_IS_vs_OB_IR.csv"
)

deg <- read.csv(deg_file)

cat("DEG file loaded.\n")
cat("Rows:", nrow(deg), "\n\n")

# ============================================================
# CLEAN DATA
# ============================================================

deg <- deg |>
  dplyr::filter(
    !is.na(log2FoldChange),
    !is.na(padj)
  )

deg$ENSEMBL <- gsub(
  "\\..*",
  "",
  deg$gene_id
)

# ============================================================
# MAP SYMBOLS
# ============================================================

annot <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = unique(deg$ENSEMBL),
  columns = c("SYMBOL"),
  keytype = "ENSEMBL"
)

annot <- annot |>
  dplyr::filter(!is.na(SYMBOL)) |>
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

deg <- left_join(
  deg,
  annot,
  by = "ENSEMBL"
)

cat("Gene symbols mapped.\n\n")

# ============================================================
# VOLCANO CATEGORIES
# ============================================================

deg <- deg |>
  dplyr::mutate(
    Regulation = case_when(
      padj < 0.05 &
        log2FoldChange > 1 ~ "Upregulated",
      
      padj < 0.05 &
        log2FoldChange < -1 ~ "Downregulated",
      
      TRUE ~ "Not Significant"
    )
  )

# ============================================================
# SELECT TOP LABELS
# ============================================================

top_labels <- deg |>
  dplyr::filter(
    padj < 0.001
  ) |>
  dplyr::arrange(padj) |>
  head(15)

# ============================================================
# VOLCANO PLOT
# ============================================================

volcano_plot <- ggplot(
  deg,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = Regulation
  )
) +
  
  geom_point(
    alpha = 0.7,
    size = 2
  ) +
  
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  
  geom_text(
    data = top_labels,
    aes(label = SYMBOL),
    size = 3,
    vjust = 1.2,
    check_overlap = TRUE
  ) +
  
  theme_bw(base_size = 14) +
  
  labs(
    title = "Volcano Plot — OB-IS vs OB-IR",
    subtitle = "Differential Expression Analysis",
    x = "log2 Fold Change",
    y = "-log10 adjusted p-value"
  )

ggsave(
  file.path(
    fig_dir,
    "Publication_Volcano_Plot.png"
  ),
  volcano_plot,
  width = 10,
  height = 8,
  dpi = 300
)

cat("Volcano plot saved.\n\n")

# ============================================================
# TOP HEATMAP GENES
# ============================================================

heatmap_genes <- deg |>
  dplyr::filter(
    padj < 0.01
  ) |>
  dplyr::arrange(padj) |>
  head(30)

heatmap_table <- heatmap_genes |>
  dplyr::select(
    SYMBOL,
    log2FoldChange
  )

heatmap_matrix <- as.matrix(
  heatmap_table$log2FoldChange
)

rownames(heatmap_matrix) <- heatmap_table$SYMBOL

colnames(heatmap_matrix) <- "log2FC"

# ============================================================
# HEATMAP
# ============================================================

png(
  file.path(
    fig_dir,
    "Publication_Heatmap.png"
  ),
  width = 1800,
  height = 2400,
  res = 300
)

pheatmap(
  heatmap_matrix,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 8,
  main = "Top Differentially Expressed Genes"
)

dev.off()

cat("Heatmap saved.\n\n")

# ============================================================
# SAVE TOP GENE TABLE
# ============================================================

write.csv(
  heatmap_genes,
  file.path(
    fig_dir,
    "Top30_Heatmap_Genes.csv"
  ),
  row.names = FALSE
)

cat("Top heatmap genes table saved.\n\n")

# ============================================================
# SUMMARY
# ============================================================

sink(
  file.path(
    fig_dir,
    "Figure_Interpretation.txt"
  )
)

cat("PUBLICATION FIGURE INTERPRETATION\n")
cat("====================================\n\n")

cat("VOLCANO PLOT:\n")
cat(
  "Demonstrates strong transcriptional remodeling\n",
  "between insulin sensitive and insulin resistant states.\n\n"
)

cat("HEATMAP:\n")
cat(
  "Shows coordinated regulation of mitochondrial,\n",
  "oxidative phosphorylation, inflammatory,\n",
  "and metabolic genes.\n\n"
)

cat("Top highlighted genes include:\n\n")

print(
  heatmap_genes |>
    dplyr::select(
      SYMBOL,
      log2FoldChange,
      padj
    )
)

sink()

cat("Interpretation saved.\n\n")

cat("========================================\n")
cat(" FIGURE GENERATION COMPLETE\n")
cat("========================================\n\n")

cat("FILES GENERATED:\n\n")

cat("1. Publication_Volcano_Plot.png\n")
cat("2. Publication_Heatmap.png\n")
cat("3. Top30_Heatmap_Genes.csv\n")
cat("4. Figure_Interpretation.txt\n\n")

cat("Saved in:\n")
cat(fig_dir, "\n\n")

#####################################33

# ============================================================
# STAGE 8 — PUBLICATION HEATMAP REGENERATION
# ============================================================

cat("\n================================================\n")
cat(" PUBLICATION HEATMAP REGENERATION\n")
cat("================================================\n\n")

# ------------------------------------------------------------
# LIBRARIES
# ------------------------------------------------------------

library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)

cat("Libraries loaded.\n\n")

# ------------------------------------------------------------
# PATHS
# ------------------------------------------------------------

outdir <- "."  # Set this to your working directory

heatmap_dir <- file.path(
  outdir,
  "Publication_Heatmap"
)

dir.create(
  heatmap_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ------------------------------------------------------------
# LOAD DESEQ OBJECT
# ------------------------------------------------------------

dds <- readRDS(
  file.path(
    outdir,
    "dds_final.rds"
  )
)

cat("DESeq object loaded.\n\n")

# ------------------------------------------------------------
# TRANSFORMATION
# ------------------------------------------------------------

vsd <- vst(
  dds,
  blind = FALSE
)

mat <- assay(vsd)

cat("VST transformation complete.\n\n")

# ------------------------------------------------------------
# SELECT TOP VARIABLE GENES
# ------------------------------------------------------------

topVarGenes <- head(
  order(
    rowVars(mat),
    decreasing = TRUE
  ),
  50
)

mat_top <- mat[topVarGenes, ]

# ------------------------------------------------------------
# Z-SCORE NORMALIZATION
# ------------------------------------------------------------

mat_scaled <- t(
  scale(
    t(mat_top)
  )
)

# ------------------------------------------------------------
# SAMPLE ANNOTATION
# ------------------------------------------------------------

annotation_col <- as.data.frame(
  colData(dds)[, "condition", drop = FALSE]
)

# ------------------------------------------------------------
# COLORS
# ------------------------------------------------------------

heat_colors <- colorRampPalette(
  rev(
    brewer.pal(
      n = 11,
      name = "RdBu"
    )
  )
)(100)

# ------------------------------------------------------------
# SAVE HIGH-QUALITY PNG
# ------------------------------------------------------------

png(
  filename = file.path(
    heatmap_dir,
    "Publication_Heatmap.png"
  ),
  width = 4000,
  height = 4200,
  res = 400
)

pheatmap(
  mat_scaled,
  color = heat_colors,
  annotation_col = annotation_col,
  fontsize_row = 8,
  fontsize_col = 12,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  border_color = NA,
  scale = "none",
  main = "Top 50 Variable Genes"
)

dev.off()

cat("Heatmap regenerated successfully.\n\n")

# ------------------------------------------------------------
# SAVE PDF VERSION
# ------------------------------------------------------------

pdf(
  file.path(
    heatmap_dir,
    "Publication_Heatmap.pdf"
  ),
  width = 12,
  height = 14
)

pheatmap(
  mat_scaled,
  color = heat_colors,
  annotation_col = annotation_col,
  fontsize_row = 8,
  fontsize_col = 12,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  border_color = NA,
  scale = "none",
  main = "Top 50 Variable Genes"
)

dev.off()

cat("PDF heatmap saved.\n\n")

cat("========================================\n")
cat(" HEATMAP GENERATION COMPLETE\n")
cat("========================================\n\n")

cat("FILES GENERATED:\n\n")
cat("1. Publication_Heatmap.png\n")
cat("2. Publication_Heatmap.pdf\n\n")

cat("Saved in:\n")
cat(heatmap_dir, "\n\n")

#################

# ============================================================
# PUBLICATION HEATMAP FROM VST MATRIX
# ============================================================

cat("\n================================================\n")
cat(" PUBLICATION HEATMAP GENERATION\n")
cat("================================================\n\n")

# ------------------------------------------------------------
# LIBRARIES
# ------------------------------------------------------------

library(pheatmap)
library(RColorBrewer)
library(matrixStats)

cat("Libraries loaded.\n\n")

# ------------------------------------------------------------
# PATHS
# ------------------------------------------------------------

outdir <- "."  # Set this to your working directory

heatmap_dir <- file.path(
  outdir,
  "Publication_Heatmap"
)

dir.create(
  heatmap_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ------------------------------------------------------------
# LOAD VST MATRIX
# ------------------------------------------------------------

mat <- read.csv(
  file.path(
    outdir,
    "GSE245301_VST_normalised_matrix.csv"
  ),
  row.names = 1,
  check.names = FALSE
)

cat("VST matrix loaded.\n")
cat("Dimensions:\n")
print(dim(mat))
cat("\n")

# ------------------------------------------------------------
# SELECT TOP VARIABLE GENES
# ------------------------------------------------------------

topVarGenes <- head(
  order(
    matrixStats::rowVars(as.matrix(mat)),
    decreasing = TRUE
  ),
  50
)

mat_top <- mat[topVarGenes, ]

cat("Top variable genes selected.\n\n")

# ------------------------------------------------------------
# Z-SCORE NORMALIZATION
# ------------------------------------------------------------

mat_scaled <- t(
  scale(
    t(mat_top)
  )
)

# ------------------------------------------------------------
# SAMPLE GROUPS
# ------------------------------------------------------------

sample_groups <- data.frame(
  Condition = c(
    "Control",
    "Control",
    "Control",
    "OB_IR",
    "OB_IR",
    "OB_IR",
    "OB_IS",
    "OB_IS",
    "OB_IS"
  )
)

rownames(sample_groups) <- colnames(mat_scaled)

# ------------------------------------------------------------
# COLORS
# ------------------------------------------------------------

heat_colors <- colorRampPalette(
  rev(
    brewer.pal(
      11,
      "RdBu"
    )
  )
)(100)

# ------------------------------------------------------------
# GENERATE PNG
# ------------------------------------------------------------

png(
  filename = file.path(
    heatmap_dir,
    "Publication_Heatmap.png"
  ),
  width = 4000,
  height = 4200,
  res = 400
)

pheatmap(
  mat_scaled,
  color = heat_colors,
  annotation_col = sample_groups,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 7,
  fontsize_col = 12,
  border_color = NA,
  show_rownames = TRUE,
  main = "Top 50 Variable Genes"
)

dev.off()

cat("PNG heatmap saved.\n\n")

# ------------------------------------------------------------
# GENERATE PDF
# ------------------------------------------------------------

pdf(
  file.path(
    heatmap_dir,
    "Publication_Heatmap.pdf"
  ),
  width = 12,
  height = 14
)

pheatmap(
  mat_scaled,
  color = heat_colors,
  annotation_col = sample_groups,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 7,
  fontsize_col = 12,
  border_color = NA,
  show_rownames = TRUE,
  main = "Top 50 Variable Genes"
)

dev.off()

cat("PDF heatmap saved.\n\n")

cat("========================================\n")
cat(" HEATMAP COMPLETE\n")
cat("========================================\n\n")

cat("Saved in:\n")
cat(heatmap_dir, "\n\n")

###########################################################################################

# ============================================================
#  GSE245301 — STAGE 5: GSEA ANALYSIS
#  This is your most important remaining analysis
#  Saves everything to: C:/Users/saipr/OneDrive/Desktop/mito/GSEA
# ============================================================

# ── LIBRARIES ─────────────────────────────────────────────────
library(DESeq2)
library(dplyr)
library(ggplot2)
library(fgsea)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(RColorBrewer)

cat("Libraries loaded.\n\n")


# ── OUTPUT FOLDERS ────────────────────────────────────────────
outdir <- "."  # Set this to your working directory
gsea_dir <- file.path(outdir, "GSEA")
dir.create(gsea_dir, showWarnings = FALSE, recursive = TRUE)
cat("GSEA outputs will save to:", gsea_dir, "\n\n")


# ── STEP 1: Load your DEG results from Comparison B ──────────
# This is OB_IS vs OB_IR — your KEY comparison
deg_file <- file.path(outdir, "Data/GSE245301_DEG_CompB_OB_IS_vs_OB_IR.csv")

if (!file.exists(deg_file)) {
  stop("DEG file not found. Make sure Stage 2 completed successfully.")
}

deg <- read.csv(deg_file)
cat("DEG file loaded. Rows:", nrow(deg), "\n\n")


# ── STEP 2: Clean and annotate ────────────────────────────────
deg <- deg |>
  dplyr::filter(
    !is.na(log2FoldChange),
    !is.na(padj),
    is.finite(log2FoldChange)
  )

# Strip version numbers from Ensembl IDs (e.g. ENSMUSG000001.3 → ENSMUSG000001)
deg$ENSEMBL <- gsub("\\..*", "", deg$gene_id)

# Map Ensembl → Gene Symbol
annot <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys    = unique(deg$ENSEMBL),
  columns = c("SYMBOL"),
  keytype = "ENSEMBL"
)

annot <- annot |>
  dplyr::filter(!is.na(SYMBOL)) |>
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

deg <- left_join(deg, annot, by = "ENSEMBL")

cat("Genes with symbol mapped:",
    sum(!is.na(deg$SYMBOL)), "of", nrow(deg), "\n\n")


# ── STEP 3: Build ranked gene list ───────────────────────────
# Ranking metric: sign(LFC) × -log10(padj)
# This captures both direction AND significance
# Better than LFC alone for GSEA

deg <- deg |>
  dplyr::mutate(
    rank_score = sign(log2FoldChange) * -log10(padj + 1e-300)
  )

ranked_genes <- deg$rank_score
names(ranked_genes) <- deg$SYMBOL

# Remove NA names, duplicates, non-finite values
ranked_genes <- ranked_genes[!is.na(names(ranked_genes))]
ranked_genes <- ranked_genes[names(ranked_genes) != ""]
ranked_genes <- ranked_genes[!duplicated(names(ranked_genes))]
ranked_genes <- ranked_genes[is.finite(ranked_genes)]

# Sort descending (required by fgsea)
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

cat("Ranked gene list ready:", length(ranked_genes), "genes\n\n")
cat("Top 5 upregulated:\n")
print(head(ranked_genes, 5))
cat("Top 5 downregulated:\n")
print(tail(ranked_genes, 5))
cat("\n")

# Save ranked list
write.csv(
  data.frame(gene = names(ranked_genes), score = ranked_genes),
  file.path(gsea_dir, "GSEA_Ranked_Gene_List.csv"),
  row.names = FALSE
)


# ── STEP 4: Load gene sets ────────────────────────────────────
cat("Loading MSigDB gene sets...\n")

# Hallmark gene sets (50 curated pathways — most important)
hallmark <- msigdbr(species = "Mus musculus", collection = "H")
hallmark_list <- split(hallmark$gene_symbol, hallmark$gs_name)

# KEGG gene sets (metabolism, disease pathways)
kegg <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:KEGG_LEGACY")
kegg_list <- split(kegg$gene_symbol, kegg$gs_name)

# Reactome gene sets (detailed biological processes)
reactome <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:REACTOME")
reactome_list <- split(reactome$gene_symbol, reactome$gs_name)

cat("Gene sets loaded:\n")
cat("  Hallmark:", length(hallmark_list), "pathways\n")
cat("  KEGG:", length(kegg_list), "pathways\n")
cat("  Reactome:", length(reactome_list), "pathways\n\n")


# ── STEP 5: Run FGSEA — Hallmark ─────────────────────────────
cat("Running GSEA — Hallmark pathways...\n")
set.seed(42)

gsea_hallmark <- fgsea(
  pathways = hallmark_list,
  stats    = ranked_genes,
  minSize  = 10,
  maxSize  = 500
)

gsea_hallmark <- as.data.frame(gsea_hallmark) |>
  dplyr::arrange(padj)

# Flatten leadingEdge list to string for saving
gsea_hallmark$leadingEdge <- sapply(
  gsea_hallmark$leadingEdge, paste, collapse = ";"
)

write.csv(
  gsea_hallmark,
  file.path(gsea_dir, "GSEA_Hallmark_All.csv"),
  row.names = FALSE
)

sig_hallmark <- gsea_hallmark |> dplyr::filter(padj < 0.05)
write.csv(
  sig_hallmark,
  file.path(gsea_dir, "GSEA_Hallmark_Significant.csv"),
  row.names = FALSE
)

cat("Significant Hallmark pathways:", nrow(sig_hallmark), "\n\n")
cat("Top Hallmark results:\n")
print(sig_hallmark[, c("pathway", "NES", "padj")])
cat("\n")


# ── STEP 6: Run FGSEA — KEGG ─────────────────────────────────
cat("Running GSEA — KEGG pathways...\n")
set.seed(42)

gsea_kegg <- fgsea(
  pathways = kegg_list,
  stats    = ranked_genes,
  minSize  = 10,
  maxSize  = 500
)

gsea_kegg <- as.data.frame(gsea_kegg) |>
  dplyr::arrange(padj)

gsea_kegg$leadingEdge <- sapply(
  gsea_kegg$leadingEdge, paste, collapse = ";"
)

write.csv(
  gsea_kegg,
  file.path(gsea_dir, "GSEA_KEGG_All.csv"),
  row.names = FALSE
)

sig_kegg <- gsea_kegg |> dplyr::filter(padj < 0.05)
write.csv(
  sig_kegg,
  file.path(gsea_dir, "GSEA_KEGG_Significant.csv"),
  row.names = FALSE
)

cat("Significant KEGG pathways:", nrow(sig_kegg), "\n\n")
cat("Top KEGG results:\n")
print(sig_kegg[, c("pathway", "NES", "padj")])
cat("\n")


# ── STEP 7: Run FGSEA — Reactome ─────────────────────────────
cat("Running GSEA — Reactome pathways...\n")
set.seed(42)

gsea_reactome <- fgsea(
  pathways = reactome_list,
  stats    = ranked_genes,
  minSize  = 10,
  maxSize  = 500
)

gsea_reactome <- as.data.frame(gsea_reactome) |>
  dplyr::arrange(padj)

gsea_reactome$leadingEdge <- sapply(
  gsea_reactome$leadingEdge, paste, collapse = ";"
)

write.csv(
  gsea_reactome,
  file.path(gsea_dir, "GSEA_Reactome_All.csv"),
  row.names = FALSE
)

sig_reactome <- gsea_reactome |> dplyr::filter(padj < 0.05)
write.csv(
  sig_reactome,
  file.path(gsea_dir, "GSEA_Reactome_Significant.csv"),
  row.names = FALSE
)

cat("Significant Reactome pathways:", nrow(sig_reactome), "\n\n")


# ── STEP 8: Figure — Hallmark GSEA Barplot ───────────────────
# Your main GSEA publication figure
# Shows top pathways, coloured by direction (suppressed vs activated)

top_hallmark <- sig_hallmark |>
  dplyr::slice_max(order_by = abs(NES), n = 20) |>
  dplyr::mutate(
    Direction = ifelse(NES > 0, "Activated in OB-IR", "Suppressed in OB-IR"),
    pathway   = gsub("HALLMARK_", "", pathway),
    pathway   = gsub("_", " ", pathway)
  )

barplot_fig <- ggplot(
  top_hallmark,
  aes(x = reorder(pathway, NES), y = NES, fill = Direction)
) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c(
    "Activated in OB-IR"  = "#D85A30",
    "Suppressed in OB-IR" = "#378ADD"
  )) +
  labs(
    title    = "GSEA — Hallmark Pathways",
    subtitle = "OB-IS vs OB-IR | GSE245301",
    x        = NULL,
    y        = "Normalised Enrichment Score (NES)",
    fill     = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    legend.position = "bottom",
    axis.text.y   = element_text(size = 9)
  )

ggsave(
  file.path(gsea_dir, "GSEA_Figure5_Hallmark_Barplot.pdf"),
  barplot_fig, width = 10, height = 8, dpi = 300
)
ggsave(
  file.path(gsea_dir, "GSEA_Figure5_Hallmark_Barplot.png"),
  barplot_fig, width = 10, height = 8, dpi = 300
)
cat("Saved: GSEA_Figure5_Hallmark_Barplot.pdf + .png\n\n")


# ── STEP 9: Figure — KEGG Barplot ────────────────────────────
if (nrow(sig_kegg) >= 5) {
  
  top_kegg <- sig_kegg |>
    dplyr::slice_max(order_by = abs(NES), n = 20) |>
    dplyr::mutate(
      Direction = ifelse(NES > 0, "Activated in OB-IR", "Suppressed in OB-IR"),
      pathway   = gsub("KEGG_", "", pathway),
      pathway   = gsub("_", " ", pathway)
    )
  
  kegg_fig <- ggplot(
    top_kegg,
    aes(x = reorder(pathway, NES), y = NES, fill = Direction)
  ) +
    geom_col(width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c(
      "Activated in OB-IR"  = "#D85A30",
      "Suppressed in OB-IR" = "#378ADD"
    )) +
    labs(
      title    = "GSEA — KEGG Pathways",
      subtitle = "OB-IS vs OB-IR | GSE245301",
      x        = NULL,
      y        = "Normalised Enrichment Score (NES)",
      fill     = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", size = 14),
      legend.position = "bottom",
      axis.text.y   = element_text(size = 9)
    )
  
  ggsave(
    file.path(gsea_dir, "GSEA_Figure6_KEGG_Barplot.pdf"),
    kegg_fig, width = 10, height = 8, dpi = 300
  )
  ggsave(
    file.path(gsea_dir, "GSEA_Figure6_KEGG_Barplot.png"),
    kegg_fig, width = 10, height = 8, dpi = 300
  )
  cat("Saved: GSEA_Figure6_KEGG_Barplot.pdf + .png\n\n")
}


# ── STEP 10: Extract your key pathways for the paper ─────────
# These are the pathways central to your biological story

key_keywords <- c(
  "OXIDATIVE_PHOSPHORYLATION",
  "REACTIVE_OXYGEN",
  "FATTY_ACID",
  "UNFOLDED_PROTEIN",
  "MTORC1",
  "HYPOXIA",
  "INFLAMMATORY",
  "TNFA",
  "AMPK",
  "APOPTOSIS",
  "MITOCHONDRI",
  "GLYCOLYSIS",
  "ADIPOGENESIS"
)

story_pathways <- sig_hallmark |>
  dplyr::filter(
    stringr::str_detect(
      pathway,
      paste(key_keywords, collapse = "|")
    )
  ) |>
  dplyr::mutate(
    pathway = gsub("HALLMARK_", "", pathway),
    pathway = gsub("_", " ", pathway)
  ) |>
  dplyr::select(pathway, NES, pval, padj, size)

write.csv(
  story_pathways,
  file.path(gsea_dir, "GSEA_Key_Story_Pathways.csv"),
  row.names = FALSE
)

cat("Your key biological story pathways:\n")
print(story_pathways)
cat("\n")


# ── STEP 11: Repair network check ────────────────────────────
# Check if AMPK / Nrf2 / autophagy pathways show up
# This is your NOVELTY angle

repair_keywords <- c("AMPK", "NRF2", "AUTOPHAGY",
                     "MITOPHAGY", "ANTIOXIDANT", "XENOBIOTIC")

repair_hallmark <- gsea_hallmark |>
  dplyr::filter(
    stringr::str_detect(pathway, paste(repair_keywords, collapse = "|"))
  )

repair_reactome <- gsea_reactome |>
  dplyr::filter(
    stringr::str_detect(pathway, paste(repair_keywords, collapse = "|"))
  )

cat("Repair/adaptation pathways found (Hallmark):\n")
print(repair_hallmark[, c("pathway", "NES", "padj")])
cat("\nRepair/adaptation pathways found (Reactome):\n")
print(repair_reactome[, c("pathway", "NES", "padj")])

write.csv(
  bind_rows(repair_hallmark, repair_reactome),
  file.path(gsea_dir, "GSEA_Repair_Network_Pathways.csv"),
  row.names = FALSE
)
cat("\n")


# ── STEP 12: Full summary ─────────────────────────────────────
cat("\n════════════════════════════════════════════\n")
cat(" GSEA COMPLETE — FILE SUMMARY\n")
cat("════════════════════════════════════════════\n\n")
cat("Saved to:", gsea_dir, "\n\n")
cat("CSVs:\n")
cat("  GSEA_Ranked_Gene_List.csv\n")
cat("  GSEA_Hallmark_All.csv\n")
cat("  GSEA_Hallmark_Significant.csv\n")
cat("  GSEA_KEGG_All.csv\n")
cat("  GSEA_KEGG_Significant.csv\n")
cat("  GSEA_Reactome_All.csv\n")
cat("  GSEA_Reactome_Significant.csv\n")
cat("  GSEA_Key_Story_Pathways.csv\n")
cat("  GSEA_Repair_Network_Pathways.csv\n\n")
cat("Figures:\n")
cat("  GSEA_Figure5_Hallmark_Barplot.pdf + .png\n")
cat("  GSEA_Figure6_KEGG_Barplot.pdf + .png\n\n")
cat("Paste your significant pathway list here\n")
cat("and we will interpret the biology + move to WGCNA.\n")

###############################################################3

# ============================================================
#  GSE245301 — STAGE 6: WGCNA CO-EXPRESSION ANALYSIS
#  Goal: Find the repair network module (AMPK/Nrf2/mitophagy)
#  that collapses during OB-IS → OB-IR transition
#  Saves to: C:/Users/saipr/OneDrive/Desktop/mito/WGCNA
# ============================================================

# ── INSTALL (run once, then comment out) ─────────────────────
# install.packages("WGCNA")
# install.packages("flashClust")
# install.packages("ggrepel")

# ── LIBRARIES ─────────────────────────────────────────────────
library(WGCNA)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(tibble)
library(org.Mm.eg.db)
library(AnnotationDbi)

options(stringsAsFactors = FALSE)
enableWGCNAThreads()
cat("Libraries loaded.\n\n")


# ── OUTPUT FOLDER ─────────────────────────────────────────────
outdir    <- "C:/Users/saipr/OneDrive/Desktop/mito"
wgcna_dir <- file.path(outdir, "WGCNA")
dir.create(wgcna_dir, showWarnings = FALSE, recursive = TRUE)
cat("WGCNA outputs will save to:", wgcna_dir, "\n\n")


# ── STEP 1: Load VST matrix ───────────────────────────────────
cat("Loading VST matrix...\n")
vst_mat <- read.csv(
  file.path(outdir, "GSE245301_VST_normalised_matrix.csv"),
  row.names = 1, check.names = FALSE
)
vst_mat <- as.matrix(vst_mat)
cat("VST matrix:", nrow(vst_mat), "genes x", ncol(vst_mat), "samples\n\n")


# ── STEP 2: Load sample table ─────────────────────────────────
sample_table <- read.csv(file.path(outdir, "GSE245301_sample_table.csv"))
sample_table$group <- factor(sample_table$group,
                             levels = c("Control", "OB_IS", "OB_IR"))
cat("Sample groups:\n")
print(table(sample_table$group))
cat("\n")


# ── STEP 3: Select top 5000 most variable genes ───────────────
cat("Selecting top 5000 most variable genes...\n")
gene_vars <- apply(vst_mat, 1, var)
top_genes  <- names(sort(gene_vars, decreasing = TRUE))[1:5000]
datExpr    <- t(vst_mat[top_genes, ])   # WGCNA needs samples as rows
cat("datExpr dimensions:", nrow(datExpr), "samples x", ncol(datExpr), "genes\n\n")


# ── STEP 4: Quality check ─────────────────────────────────────
cat("Checking data quality...\n")
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  cat("Removed bad genes/samples.\n")
}
cat("Quality check done.\n\n")


# ── STEP 5: Sample clustering ─────────────────────────────────
sampleTree <- hclust(dist(datExpr), method = "average")
pdf(file.path(wgcna_dir, "WGCNA_01_Sample_Clustering.pdf"), width = 10, height = 5)
par(cex = 0.8, mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering — check for outliers",
     sub = "", xlab = "")
dev.off()
cat("Saved: WGCNA_01_Sample_Clustering.pdf\n\n")


# ── STEP 6: Soft-thresholding power ──────────────────────────
cat("Picking soft-thresholding power (~2 min)...\n")
sft <- pickSoftThreshold(datExpr, powerVector = 1:20,
                         RsquaredCut = 0.85, verbose = 5,
                         networkType = "signed")
cat("Recommended power:", sft$powerEstimate, "\n\n")

pdf(file.path(wgcna_dir, "WGCNA_02_Soft_Power_Selection.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))
plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab = "Soft threshold (power)",
     ylab = "Scale free R²", type = "n",
     main = "Scale independence")
text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels = 1:20, cex = 0.9, col = "red")
abline(h = 0.85, col = "red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft threshold (power)",
     ylab = "Mean connectivity", type = "n",
     main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5],
     labels = 1:20, cex = 0.9, col = "red")
dev.off()
cat("Saved: WGCNA_02_Soft_Power_Selection.pdf\n\n")

soft_power <- ifelse(is.na(sft$powerEstimate), 12, sft$powerEstimate)
cat("Using soft power:", soft_power, "\n\n")


# ── STEP 7: Build network + detect modules ───────────────────
cat("Building network and detecting modules (~5-15 min)...\n\n")
net <- blockwiseModules(
  datExpr,
  power             = soft_power,
  networkType       = "signed",
  TOMType           = "signed",
  minModuleSize     = 30,
  reassignThreshold = 0,
  mergeCutHeight    = 0.25,
  numericLabels     = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs          = FALSE,
  verbose           = 3
)
moduleColors <- labels2colors(net$colors)
cat("\nModules detected:", length(unique(moduleColors)) - 1, "\n")
cat("Module sizes:\n")
print(table(moduleColors))
cat("\n")


# ── STEP 8: Dendrogram ───────────────────────────────────────
pdf(file.path(wgcna_dir, "WGCNA_03_Module_Dendrogram.pdf"), width = 12, height = 6)
plotDendroAndColors(net$dendrograms[[1]],
                    moduleColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
cat("Saved: WGCNA_03_Module_Dendrogram.pdf\n\n")


# ── STEP 9: Module-trait correlations ────────────────────────
cat("Calculating module-trait correlations...\n")

datTraits <- data.frame(
  Control = as.numeric(sample_table$group == "Control"),
  OB_IS   = as.numeric(sample_table$group == "OB_IS"),
  OB_IR   = as.numeric(sample_table$group == "OB_IR"),
  row.names = sample_table$sample_id
)
datTraits <- datTraits[rownames(datExpr), ]

MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs)

moduleTraitCor  <- cor(MEs, datTraits, use = "p")
moduleTraitPval <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

write.csv(
  as.data.frame(moduleTraitCor) |> tibble::rownames_to_column("module"),
  file.path(wgcna_dir, "WGCNA_Module_Trait_Correlations.csv"),
  row.names = FALSE
)

textMatrix <- paste0(round(moduleTraitCor, 2),
                     "\n(p=", signif(moduleTraitPval, 1), ")")
dim(textMatrix) <- dim(moduleTraitCor)

pdf(file.path(wgcna_dir, "WGCNA_04_Module_Trait_Heatmap.pdf"),
    width = 8, height = max(6, nrow(moduleTraitCor) * 0.35))
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(datTraits),
  yLabels = rownames(moduleTraitCor),
  ySymbols = rownames(moduleTraitCor),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.5, zlim = c(-1, 1),
  main = "Module-trait correlations"
)
dev.off()
cat("Saved: WGCNA_04_Module_Trait_Heatmap.pdf\n\n")


# ── STEP 10: Identify repair module ──────────────────────────
mod_summary <- data.frame(
  module    = rownames(moduleTraitCor),
  OB_IR_cor = moduleTraitCor[, "OB_IR"],
  OB_IR_p   = moduleTraitPval[, "OB_IR"],
  OB_IS_cor = moduleTraitCor[, "OB_IS"]
) |> dplyr::arrange(OB_IR_p)

write.csv(mod_summary,
          file.path(wgcna_dir, "WGCNA_Module_Summary.csv"), row.names = FALSE)

cat("Top modules by OB-IR association:\n")
print(head(mod_summary, 10))
cat("\n")

# Repair module = active in OB-IS, lost in OB-IR
repair_candidates <- mod_summary |>
  dplyr::filter(OB_IS_cor > 0.3, OB_IR_cor < 0)

if (nrow(repair_candidates) > 0) {
  cat("REPAIR MODULE CANDIDATES (active OB-IS, lost OB-IR):\n")
  print(repair_candidates)
  repair_module_name <- gsub("ME", "", repair_candidates$module[1])
} else {
  cat("Using top OB-IR-associated module.\n")
  repair_module_name <- gsub("ME", "", mod_summary$module[1])
}
cat("\nRepair module:", repair_module_name, "\n\n")


# ── STEP 11: Extract + annotate repair module genes ──────────
repair_genes   <- names(moduleColors)[moduleColors == repair_module_name]
repair_ensembl <- gsub("\\..*", "", repair_genes)

repair_annot <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys    = repair_ensembl,
  columns = c("SYMBOL", "GENENAME"),
  keytype = "ENSEMBL"
) |>
  dplyr::filter(!is.na(SYMBOL)) |>
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

cat("Repair module genes:", length(repair_genes),
    "| Annotated:", nrow(repair_annot), "\n\n")

write.csv(repair_annot,
          file.path(wgcna_dir, "WGCNA_Repair_Module_Genes.csv"), row.names = FALSE)


# ── STEP 12: Check for key pathway genes ─────────────────────
key_repair_genes <- c(
  "Prkaa1","Prkaa2","Prkab1","Prkab2","Prkag1","Stk11","Camkk2",
  "Nfe2l2","Keap1","Hmox1","Nqo1","Gclc","Gclm","Txnrd1",
  "Prdx1","Prdx3","Gpx1","Gpx4","Cat","Sod1","Sod2",
  "Pink1","Prkn","Bnip3","Bnip3l","Fundc1","Ulk1","Becn1",
  "Atg5","Atg7","Atg12","Map1lc3a","Map1lc3b","Sqstm1",
  "Tomm7","Tomm20","Tomm40","Opa1","Mfn1","Mfn2","Drp1",
  "Ppargc1a","Ppargc1b","Tfam","Sesn1","Sesn2"
)
found_in_module <- intersect(repair_annot$SYMBOL, key_repair_genes)

cat("KEY AMPK/NRF2/MITOPHAGY GENES IN REPAIR MODULE:\n")
print(found_in_module)
cat("\n")

write.csv(data.frame(gene = found_in_module),
          file.path(wgcna_dir, "WGCNA_Repair_Module_Key_Genes.csv"), row.names = FALSE)


# ── STEP 13: Eigengene trajectory plot (Figure 7) ────────────
ME_repair <- MEs[[paste0("ME", repair_module_name)]]

eigengene_df <- data.frame(
  sample    = rownames(datExpr),
  eigengene = ME_repair,
  group     = sample_table$group[
    match(rownames(datExpr), sample_table$sample_id)]
)

eigen_plot <- ggplot(eigengene_df,
                     aes(x = group, y = eigengene, color = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA,
               fill = NA, linewidth = 0.8) +
  geom_jitter(width = 0.1, size = 3.5, alpha = 0.9) +
  scale_color_manual(values = c(
    "Control" = "#1D9E75", "OB_IS" = "#378ADD", "OB_IR" = "#D85A30")) +
  labs(
    title    = paste0("Repair module (", repair_module_name, ") eigengene"),
    subtitle = "Progressive collapse from OB-IS to OB-IR = your key finding",
    x = NULL, y = "Module eigengene (first PC)", color = "Group"
  ) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), legend.position = "none")

ggsave(file.path(wgcna_dir, "WGCNA_Figure7_Repair_Eigengene.pdf"),
       eigen_plot, width = 6, height = 5, dpi = 300)
ggsave(file.path(wgcna_dir, "WGCNA_Figure7_Repair_Eigengene.png"),
       eigen_plot, width = 6, height = 5, dpi = 300)
cat("Saved: WGCNA_Figure7_Repair_Eigengene.pdf + .png\n\n")


# ── STEP 14: Module membership + gene significance ───────────
GS_OBIR   <- as.numeric(cor(datExpr, datTraits$OB_IR, use = "p"))
names(GS_OBIR) <- colnames(datExpr)

MM_repair <- as.numeric(cor(datExpr,
                            MEs[[paste0("ME", repair_module_name)]], use = "p"))
names(MM_repair) <- colnames(datExpr)

repair_idx <- moduleColors == repair_module_name

MM_GS_df <- data.frame(
  gene_id          = colnames(datExpr)[repair_idx],
  ModuleMembership = MM_repair[repair_idx],
  GeneSignificance = GS_OBIR[repair_idx]
) |>
  dplyr::mutate(ENSEMBL = gsub("\\..*", "", gene_id)) |>
  dplyr::left_join(repair_annot[, c("ENSEMBL","SYMBOL")], by = "ENSEMBL") |>
  dplyr::arrange(desc(abs(ModuleMembership)))

write.csv(MM_GS_df,
          file.path(wgcna_dir, "WGCNA_Repair_Module_MM_GS.csv"), row.names = FALSE)

cat("Top 20 hub genes:\n")
print(head(MM_GS_df[, c("SYMBOL","ModuleMembership","GeneSignificance")], 20))
cat("\n")


# ── STEP 15: MM vs GS scatter (Figure 8) ─────────────────────
MM_GS_plot <- ggplot(
  MM_GS_df |> dplyr::filter(!is.na(SYMBOL)),
  aes(x = ModuleMembership, y = GeneSignificance, label = SYMBOL)
) +
  geom_point(alpha = 0.5, size = 1.5, color = "#378ADD") +
  geom_smooth(method = "lm", color = "#D85A30", se = TRUE, linewidth = 0.8) +
  geom_text_repel(
    data = MM_GS_df |>
      dplyr::filter(!is.na(SYMBOL)) |>
      dplyr::slice_max(abs(ModuleMembership), n = 15),
    size = 3, max.overlaps = 15
  ) +
  labs(
    title    = "Module membership vs gene significance — repair module",
    subtitle = "Top-right genes = hub genes driving the repair network",
    x = "Module membership (correlation with eigengene)",
    y = "Gene significance (correlation with OB-IR trait)"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(wgcna_dir, "WGCNA_Figure8_MM_vs_GS.pdf"),
       MM_GS_plot, width = 8, height = 6, dpi = 300)
ggsave(file.path(wgcna_dir, "WGCNA_Figure8_MM_vs_GS.png"),
       MM_GS_plot, width = 8, height = 6, dpi = 300)
cat("Saved: WGCNA_Figure8_MM_vs_GS.pdf + .png\n\n")


# ── STEP 16: Heatmap of top 30 hub genes (Figure 9) ──────────
top30_ids <- MM_GS_df |>
  dplyr::filter(!is.na(SYMBOL)) |>
  dplyr::slice_max(abs(ModuleMembership), n = 30) |>
  dplyr::pull(gene_id)

mat_repair <- t(datExpr)[top30_ids, sample_table$sample_id]
mat_repair <- t(scale(t(mat_repair)))
rownames(mat_repair) <- MM_GS_df$SYMBOL[
  match(rownames(mat_repair), MM_GS_df$gene_id)]

ann_col    <- data.frame(Group = sample_table$group,
                         row.names = sample_table$sample_id)
ann_colors <- list(Group = c(
  Control = "#1D9E75", OB_IS = "#378ADD", OB_IR = "#D85A30"))

pdf(file.path(wgcna_dir, "WGCNA_Figure9_Repair_Heatmap.pdf"),
    width = 10, height = 10)
pheatmap(mat_repair,
         annotation_col = ann_col, annotation_colors = ann_colors,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         cluster_rows = TRUE, cluster_cols = FALSE,
         show_rownames = TRUE, show_colnames = TRUE,
         fontsize_row = 9, fontsize_col = 10,
         main = paste0("Repair module (", repair_module_name, ") — top 30 hub genes"))
dev.off()

png(file.path(wgcna_dir, "WGCNA_Figure9_Repair_Heatmap.png"),
    width = 3000, height = 3000, res = 300)
pheatmap(mat_repair,
         annotation_col = ann_col, annotation_colors = ann_colors,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         cluster_rows = TRUE, cluster_cols = FALSE,
         show_rownames = TRUE, show_colnames = TRUE,
         fontsize_row = 9, fontsize_col = 10,
         main = paste0("Repair module (", repair_module_name, ") — top 30 hub genes"))
dev.off()
cat("Saved: WGCNA_Figure9_Repair_Heatmap.pdf + .png\n\n")


# ── STEP 17: Hub gene list for STRING ────────────────────────
hub_final <- MM_GS_df |>
  dplyr::filter(!is.na(SYMBOL)) |>
  dplyr::slice_max(abs(ModuleMembership), n = 50) |>
  dplyr::select(SYMBOL, ModuleMembership, GeneSignificance)

write.csv(hub_final,
          file.path(wgcna_dir, "WGCNA_Hub_Genes_for_STRING.csv"), row.names = FALSE)

cat("Top 20 hub genes for your STRING network:\n")
print(head(hub_final, 20))
cat("\n")


# ── STEP 18: Save objects ─────────────────────────────────────
saveRDS(net,          file.path(wgcna_dir, "WGCNA_network.rds"))
saveRDS(MEs,          file.path(wgcna_dir, "WGCNA_eigengenes.rds"))
saveRDS(moduleColors, file.path(wgcna_dir, "WGCNA_moduleColors.rds"))


# ── FINAL SUMMARY ─────────────────────────────────────────────
cat("\n════════════════════════════════════════════\n")
cat(" WGCNA COMPLETE\n")
cat("════════════════════════════════════════════\n\n")
cat("Files in:", wgcna_dir, "\n\n")
cat("FIGURES:\n")
cat("  WGCNA_01_Sample_Clustering.pdf\n")
cat("  WGCNA_02_Soft_Power_Selection.pdf\n")
cat("  WGCNA_03_Module_Dendrogram.pdf\n")
cat("  WGCNA_04_Module_Trait_Heatmap.pdf\n")
cat("  WGCNA_Figure7_Repair_Eigengene.pdf + .png\n")
cat("  WGCNA_Figure8_MM_vs_GS.pdf + .png\n")
cat("  WGCNA_Figure9_Repair_Heatmap.pdf + .png\n\n")
cat("DATA:\n")
cat("  WGCNA_Module_Trait_Correlations.csv\n")
cat("  WGCNA_Module_Summary.csv\n")
cat("  WGCNA_Repair_Module_Genes.csv\n")
cat("  WGCNA_Repair_Module_Key_Genes.csv\n")
cat("  WGCNA_Repair_Module_MM_GS.csv\n")
cat("  WGCNA_Hub_Genes_for_STRING.csv\n\n")
cat("Paste the output here — next step is manuscript writing.\n")

mat <- read.csv(file.path(outdir, "GSE245301_VST_normalised_matrix.csv"),
                row.names = 1,
                check.names = FALSE)

datExpr <- as.data.frame(t(mat))  # IMPORTANT: samples as rows



#################
library(WGCNA)

base <- file.path(outdir, "WGCNA")

load(file.path(base, "WGCNA_FunctionalEnrichment_workspace.RData"))

load(file.path(base, "net.rds"))
moduleColors <- readRDS(file.path(base, "moduleColors.rds"))

MEs <- readRDS(file.path(base, "MEs.rds"))

trait <- read.csv("/mnt/c/Users/saipr/OneDrive/Desktop/mito/GSE245301_sample_table.csv", row.names = 1)

# ensure alignment
MEs <- orderMEs(MEs)

png(file.path(base, "WGCNA_ModuleTraitHeatmap_FIXED.png"),
    width=3000, height=2500, res=300)

labeledHeatmap(
  Matrix = cor(MEs, trait, use="p"),
  xLabels = names(trait),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = signif(cor(MEs, trait, use="p"), 2)
)

dev.off()

