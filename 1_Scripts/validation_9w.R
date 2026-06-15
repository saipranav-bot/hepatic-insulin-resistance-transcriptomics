# validate_14w.R - 14-Week Cohort External Validation Pipeline
dir.create("GSE132800_validation", showWarnings = FALSE, recursive = TRUE)

required_packages <- c("tidyverse", "broom", "data.table")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos="https://cloud.r-project.org")
  }
}

library(tidyverse)
library(broom)
library(data.table)

dest_file <- "GSE132800_validation/GSE132800_fpkm_14w.txt.gz"
url_14w <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132800/suppl/GSE132800_fpkm_14w.txt.gz"

# 1. Download Block (Targeting 14w)
if(!file.exists(dest_file) || file.info(dest_file)$size == 0) {
  print("Step 1: Downloading 14-week source data matrix from NCBI...")
  tryCatch({
    download.file(url_14w, destfile = dest_file, mode = "wb")
  }, error = function(e) {
    stop("Network Error: Could not connect to NCBI server. Check your internet connection.")
  })
} else {
  print("Step 1: Verified cached 14-week FPKM file exists locally.")
}

# 2. Parser Block
print("Step 2: Loading and cleaning matrix columns...")
fpkm_matrix <- data.table::fread(dest_file, data.table = FALSE)
colnames(fpkm_matrix)[1] <- "GeneSymbol"

# 3. Clean Gene Symbols via Pipe Splitting
fpkm_matrix$GeneSymbol <- sapply(strsplit(as.character(fpkm_matrix[[1]]), "\\|"), `[`, 1)
fpkm_matrix$GeneSymbol_lower <- tolower(fpkm_matrix$GeneSymbol)

user_signature <- c("Nek7", "Tomm7", "Etfa", "Acaa1a", "Ets1", "Cox7a2", "Junb")
target_genes_lower <- tolower(user_signature)

print("Step 3: Filtering matrix for signature genes...")
expr_long <- fpkm_matrix %>%
  filter(GeneSymbol_lower %in% target_genes_lower) %>%
  select(-GeneSymbol_lower)

# Dynamically find the ID column name to drop it safely (handles variations like ID, V1, etc.)
id_col <- colnames(expr_long)[1]
expr_long <- expr_long %>%
  pivot_longer(cols = -c(GeneSymbol, !!sym(id_col)), names_to = "SampleID", values_to = "FPKM") %>%
  select(-!!sym(id_col))

# 4. Adaptive Sample Classification (Targeting 14w naming schema)
expr_long <- expr_long %>%
  mutate(Group = case_when(
    grepl("_ob", SampleID, ignore.case = TRUE) ~ "IR",
    grepl("_wt", SampleID, ignore.case = TRUE) ~ "Control",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Group))

# 5. Statistical Computations
print("Step 4: Running Welch Two-Sample T-Tests for 14-week timeline...")
validation_results <- expr_long %>%
  group_by(GeneSymbol) %>%
  do({
    df <- .
    if(length(unique(df$Group)) == 2 && var(df$FPKM, na.rm = TRUE) > 0) {
      tryCatch({
        broom::tidy(t.test(FPKM ~ Group, data = df))
      }, error = function(e) {
        data.frame(statistic = NA, p.value = NA)
      })
    } else {
      data.frame(statistic = NA, p.value = NA)
    }
  }) %>%
  ungroup()

final_output <- validation_results %>%
  select(GeneSymbol, statistic, p.value) %>%
  mutate(
    p.value = as.numeric(p.value),
    Significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE ~ "NS"
    )
  )

print("=============================================")
print("  VERIFIED 14-WEEK EXTERNAL VALIDATION METRICS ")
print("=============================================")
print(as.data.frame(final_output))
print("=============================================")

write.csv(final_output, "GSE132800_validation/validation_14w_metrics_summary.csv", row.names = FALSE)
print("14-Week pipeline complete! Output saved safely.")
