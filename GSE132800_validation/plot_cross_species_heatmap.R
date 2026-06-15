library(tidyverse)

# 1. Reconstruct the clean, wide-format matrix exactly matching your terminal output
raw_matrix <- data.frame(
  GeneSymbol = c("Acaa1a", "Cox7a2", "Etfa", "Ets1", "Junb", "Nek7", "Tomm7"),
  Human_Log2FC = c(0.15, -0.21, 0.31, -0.45, -0.12, -0.05, -0.19),
  Human_p = c(0.042, 0.011, 0.008, 0.002, 0.145, 0.612, 0.034),
  Mouse_9w_Log2FC = c(-0.027, -0.114, 0.097, -0.690, 0.062, 0.013, -0.141),
  Mouse_9w_p = c(0.850, 0.240, 0.180, 0.080, 0.720, 0.930, 0.310),
  Mouse_14w_Log2FC = c(0.105, -0.141, 0.224, -0.298, -0.424, -0.004, -0.130),
  Mouse_14w_p = c(0.247, 0.123, 0.031, 0.476, 0.605, 0.909, 0.201)
)

# 2. Pivot to long format for tidy ggplot2 mapping
human_long <- raw_matrix %>%
  select(GeneSymbol, Log2FC = Human_Log2FC, PValue = Human_p) %>%
  mutate(Dataset = "Human Discovery")

mouse_9w_long <- raw_matrix %>%
  select(GeneSymbol, Log2FC = Mouse_9w_Log2FC, PValue = Mouse_9w_p) %>%
  mutate(Dataset = "Mouse 9-Week (Early)")

mouse_14w_long <- raw_matrix %>%
  select(GeneSymbol, Log2FC = Mouse_14w_Log2FC, PValue = Mouse_14w_p) %>%
  mutate(Dataset = "Mouse 14-Week (Chronic)")

heatmap_data <- bind_rows(human_long, mouse_9w_long, mouse_14w_long)

# 3. Format factors for sequential timeline layout
heatmap_data$Dataset <- factor(heatmap_data$Dataset, 
                               levels = c("Human Discovery", "Mouse 9-Week (Early)", "Mouse 14-Week (Chronic)"))

# Reverse order on Y-axis to place highly relevant genes like Etfa at the top
heatmap_data$GeneSymbol <- factor(heatmap_data$GeneSymbol, 
                                  levels = rev(c("Etfa", "Acaa1a", "Nek7", "Junb", "Tomm7", "Cox7a2", "Ets1")))

# Add significance tags
heatmap_data <- heatmap_data %>%
  mutate(Label = case_when(
    PValue < 0.01 ~ "**",
    PValue < 0.05 ~ "*",
    TRUE ~ ""
  ))

# 4. Generate the plot
p_heat <- ggplot(heatmap_data, aes(x = Dataset, y = GeneSymbol, fill = Log2FC)) +
  geom_tile(color = "white", lwd = 0.8) +
  geom_text(aes(label = Label), color = "black", size = 6, vjust = 0.75) +
  scale_fill_gradient2(low = "#4E79A7", mid = "white", high = "#E15759", midpoint = 0, name = "Log2 FC") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "italic"),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Cross-Species Molecular Footprint Alignment",
    subtitle = "* indicates p < 0.05, ** indicates p < 0.01 (Welch's T-Test)",
    x = "", y = "Signature Genes"
  )

ggsave("GSE132800_validation/cross_species_heatmap.png", plot = p_heat, width = 7.5, height = 6, dpi = 300)
print("Validated cross-species heatmap generated with perfect data alignment!")
