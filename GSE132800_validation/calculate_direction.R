library(tidyverse)
library(data.table)

user_signature <- c("Nek7", "Tomm7", "Etfa", "Acaa1a", "Ets1", "Cox7a2", "Junb")
target_genes_lower <- tolower(user_signature)

process_file <- function(filepath, timepoint) {
  mat <- data.table::fread(filepath, data.table = FALSE)
  orig_col <- colnames(mat)[1]
  mat$GeneSymbol <- sapply(strsplit(as.character(mat[[orig_col]]), "\\|"), `[`, 1)
  
  mat %>%
    filter(tolower(GeneSymbol) %in% target_genes_lower) %>%
    select(-all_of(c(orig_col))) %>%
    pivot_longer(cols = -GeneSymbol, names_to = "SampleID", values_to = "FPKM") %>%
    mutate(Group = case_when(
      grepl("_ob", SampleID, ignore.case = TRUE) ~ "IR",
      grepl("_wt", SampleID, ignore.case = TRUE) ~ "Control",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(Group)) %>%
    group_by(GeneSymbol, Group) %>%
    summarise(Mean_FPKM = mean(FPKM, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = Group, values_from = Mean_FPKM) %>%
    mutate(
      Log2FC = log2((IR + 0.01) / (Control + 0.01)),
      Direction = if_else(Log2FC > 0, "UP in IR", "DOWN in IR"),
      Timepoint = timepoint
    )
}

df_9w <- process_file("GSE132800_validation/GSE132800_fpkm_9w.txt.gz", "9 Weeks")
df_14w <- process_file("GSE132800_validation/GSE132800_fpkm_14w.txt.gz", "14 Weeks")

final_summary <- bind_rows(df_9w, df_14w) %>%
  arrange(GeneSymbol, Timepoint)

print("=========================================================")
print("     PROGRESSION DIRECTIONALITY MATRIX (9w vs 14w)       ")
print("=========================================================")
print(as.data.frame(final_summary), row.names = FALSE)
print("=========================================================")
