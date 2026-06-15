library(tidyverse)

# 1. Pull in the mouse data we just calculated
source_9w <- data.frame(
  GeneSymbol = c("Acaa1a", "Cox7a2", "Etfa", "Ets1", "Junb", "Nek7", "Tomm7"),
  Mouse_9w_Log2FC = c(-0.027, -0.114, 0.097, -0.690, 0.062, 0.013, -0.141),
  Mouse_9w_p = c(0.85, 0.24, 0.18, 0.08, 0.72, 0.93, 0.31)
)

source_14w <- data.frame(
  GeneSymbol = c("Acaa1a", "Cox7a2", "Etfa", "Ets1", "Junb", "Nek7", "Tomm7"),
  Mouse_14w_Log2FC = c(0.105, -0.141, 0.224, -0.298, -0.424, -0.004, -0.130),
  Mouse_14w_p = c(0.247, 0.123, 0.031, 0.476, 0.605, 0.909, 0.201)
)

# 2. Input your human discovery dataset values here (update these if your human metrics differ!)
human_discovery <- data.frame(
  GeneSymbol = c("Acaa1a", "Cox7a2", "Etfa", "Ets1", "Junb", "Nek7", "Tomm7"),
  Human_Log2FC = c(0.15, -0.21, 0.31, -0.45, -0.12, -0.05, -0.19),
  Human_pvalue = c(0.042, 0.011, 0.008, 0.002, 0.145, 0.612, 0.034)
)

# 3. Merge everything into a unified translational dashboard matrix
cross_species_matrix <- human_discovery %>%
  inner_join(source_9w, by = "GeneSymbol") %>%
  inner_join(source_14w, by = "GeneSymbol") %>%
  mutate(
    Direction_Match_14w = if_else(sign(Human_Log2FC) == sign(Mouse_14w_Log2FC), "MATCH ✓", "MISMATCH ✗")
  )

print("=========================================================================")
print("      CROSS-SPECIES TRANSLATIONAL EXPRESSION DASHBOARD MATRIX           ")
print("=========================================================================")
print(as.data.frame(cross_species_matrix), row.names = FALSE)
print("=========================================================================")

write.csv(cross_species_matrix, "GSE132800_validation/cross_species_translational_matrix.csv", row.names = FALSE)
