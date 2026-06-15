# validate_14w.R - Production-Grade 14-Week Cohort External Validation
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

print("Step 1: Reading locally cached 14-week file...")
fpkm_matrix <- data.table::fread(dest_file, data.table = FALSE)

# Store the raw original ID column name to manipulate safely
original_id_col <- colnames(fpkm_matrix)[1]

print("Step 2: Splitting composite IDs to extract clean Gene Symbols...")
# Safely extract "Csprs" from "Csprs|XLOC_000025"
fpkm_matrix$GeneSymbol <- sapply(strsplit(as.character(fpkm_matrix[[original_id_col]]), "\\|"), `[`, 1)
fpkm_matrix$GeneSymbol_lower <- tolower(fpkm_matrix$GeneSymbol)

user_signature <- c("Nek7", "Tomm7", "Etfa", "Acaa1a", "Ets1", "Cox7a2", "Junb")
target_genes_lower <- tolower(user_signature)

print("Step 3: Filtering matrix for your signature genes...")
# Isolate matching rows first
filtered_matrix <- fpkm_matrix %>%
  filter(GeneSymbol_lower %in% target_genes_lower)

# Re-shape to long format while strictly keeping GeneSymbol intact
expr_long <- filtered_matrix %>%
  select(-all_of(c(original_id_col, "GeneSymbol_lower"))) %>% 
  pivot_longer(cols = -GeneSymbol, names_to = "SampleID", values_to = "FPKM")

print(paste("Successfully tracked", length(unique(expr_long$GeneSymbol)), "out of 7 genes in long format."))

# 4. Sample Classification
expr_long <- expr_long %>%
  mutate(Group = case_when(
    grepl("_ob", SampleID, ignore.case = TRUE) ~ "IR",
    grepl("_wt", SampleID, ignore.case = TRUE) ~ "Control",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Group))

# 5. Statistical Computations
print("Step 4: Running Welch Two-Sample T-Tests...")
validation_results <- expr_long %>%
  group_by(GeneSymbol) %>%
  do({
    df <- .
    if(length(unique(df$Group)) == 2) {
      broom::tidy(t.test(FPKM ~ Group, data = df))
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
