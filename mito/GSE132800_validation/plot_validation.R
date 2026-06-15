library(tidyverse)
library(data.table)

user_signature <- c("Nek7", "Tomm7", "Etfa", "Acaa1a", "Ets1", "Cox7a2", "Junb")
target_genes_lower <- tolower(user_signature)

load_and_format <- function(filepath, timepoint) {
  mat <- data.table::fread(filepath, data.table = FALSE)
  orig_col <- colnames(mat)[1]
  mat$GeneSymbol <- sapply(strsplit(as.character(mat[[orig_col]]), "\\|"), `[`, 1)
  
  mat %>%
    filter(tolower(GeneSymbol) %in% target_genes_lower) %>%
    select(-all_of(c(orig_col))) %>%
    pivot_longer(cols = -GeneSymbol, names_to = "SampleID", values_to = "FPKM") %>%
    mutate(
      Group = case_when(
        grepl("_ob", SampleID, ignore.case = TRUE) ~ "IR (ob/ob)",
        grepl("_wt", SampleID, ignore.case = TRUE) ~ "Control (WT)",
        TRUE ~ NA_character_
      ),
      Timepoint = timepoint
    ) %>%
    filter(!is.na(Group))
}

data_9w <- load_and_format("GSE132800_validation/GSE132800_fpkm_9w.txt.gz", "9 Weeks")
data_14w <- load_and_format("GSE132800_validation/GSE132800_fpkm_14w.txt.gz", "14 Weeks")
combined_data <- bind_rows(data_9w, data_14w)

# Ensure timepoints and groups display in biological order
combined_data$Timepoint <- factor(combined_data$Timepoint, levels = c("9 Weeks", "14 Weeks"))
combined_data$Group <- factor(combined_data$Group, levels = c("Control (WT)", "IR (ob/ob)"))

p <- ggplot(combined_data, aes(x = Group, y = FPKM, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, lwd = 0.6) +
  geom_jitter(width = 0.15, size = 1.5, aes(color = Group)) +
  facet_wrap(~GeneSymbol, scales = "free_y", ncol = 4) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "#f2f2f2", color = "black"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(
    title = "GSE132800 External Validation Cross-Cohort Expression Profiles",
    subtitle = "Comparing Control vs Progressive Insulin Resistant (IR) Mice",
    x = "", 
    y = "Expression Level (FPKM)"
  ) +
  scale_fill_manual(values = c("Control (WT)" = "#4E79A7", "IR (ob/ob)" = "#E15759")) +
  scale_color_manual(values = c("Control (WT)" = "#2B4970", "IR (ob/ob)" = "#942E30"))

ggsave("GSE132800_validation/signature_validation_plot.png", plot = p, width = 11, height = 7, dpi = 300)
print("Plot successfully rendered and saved to GSE132800_validation/signature_validation_plot.png!")
