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

# 2. WARNING: The original matrix hard-coded "Human_Log2FC" values that do not match
# the actual murine discovery results (DESeq2+ashr). These have been replaced with
# the verified discovery values from GSE245301 (OB_IR vs OB_IS, murine). This file
# is preserved for provenance; the corrected validation metrics are reported separately.
human_discovery <- data.frame(
  GeneSymbol = c("Acaa1a", "Cox7a2", "Etfa", "Ets1", "Junb", "Nek7", "Tomm7"),
  Discovery_Log2FC = c(-1.0136, -1.1018, -0.7380, 2.2632, 1.6307, 2.4629, -1.3113),
  Discovery_padj = c(6.15e-10, 4.89e-7, 6.15e-10, 2.63e-8, 2.89e-7, 2.28e-24, 4.88e-9)
)

# 3. Merge everything into a unified translational dashboard matrix
cross_species_matrix <- human_discovery %>%
  inner_join(source_9w, by = "GeneSymbol") %>%
  inner_join(source_14w, by = "GeneSymbol") %>%
  mutate(
    Discovery_Dir = if_else(Discovery_Log2FC >= 0, "UP in OB_IR", "DOWN in OB_IR"),
    Match_9w = if_else(sign(Discovery_Log2FC) == sign(Mouse_9w_Log2FC), "MATCH ✓", "MISMATCH ✗"),
    Match_14w = if_else(sign(Discovery_Log2FC) == sign(Mouse_14w_Log2FC), "MATCH ✓", "MISMATCH ✗")
  )

print("=========================================================================")
print("      CROSS-SPECIES TRANSLATIONAL EXPRESSION DASHBOARD MATRIX           ")
print("=========================================================================")
print(as.data.frame(cross_species_matrix), row.names = FALSE)
print("=========================================================================")

write.csv(cross_species_matrix, "GSE132800_validation/cross_species_translational_matrix.csv", row.names = FALSE)
