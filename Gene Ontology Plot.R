# Load necessary libraries
library("ggplot2")
library("dplyr")
library("ggtext")  # For improved text rendering

# Read the file
go_data <- read.table("High_vs_low.GOseq.enriched", 
                      header = TRUE, 
                      sep = "\t", 
                      stringsAsFactors = FALSE, 
                      fill = TRUE) # Fill missing values

# Assign "Unidentified" to missing terms and ontology, and create the combined label with GO term
go_data <- go_data %>%
  mutate(
    term = ifelse(is.na(term), "Unidentified", term),
    ontology = ifelse(is.na(ontology), "Unidentified", ontology)
  ) %>%
  mutate(term_with_go = paste(term, "(", category, ")", sep = ""))

# Prepare the data (filter for significant terms and take the top 20 by FDR)
go_data_filtered <- go_data %>%
  filter(over_represented_FDR < 0.05) %>%
  mutate(
    `-log10(FDR)` = -log10(over_represented_FDR),
    numDEInCat = round(numDEInCat)  # Round the number of DE genes
  ) %>%
  arrange(desc(`-log10(FDR)`)) %>%
  slice_head(n = 20)

# Enhanced scatter plot
ggplot(go_data_filtered, aes(
  x = `-log10(FDR)`, 
  y = reorder(term_with_go, `-log10(FDR)`), 
  size = numDEInCat, 
  color = ontology
)) +
  geom_point(alpha = 0.8, stroke = 1.5) +  # Added stroke for better point definition
  scale_color_manual(values = c("BP" = "#1f78b4", "CC" = "#33a02c", "MF" = "#e31a1c", "Unidentified" = "grey")) +
  scale_size_continuous(
    range = c(3, 10),
    breaks = c(min(go_data_filtered$numDEInCat), 
               median(go_data_filtered$numDEInCat), 
               max(go_data_filtered$numDEInCat))
  ) +
  theme_minimal(base_size = 24) +
  theme(
    panel.grid.major.y = element_line(color = "#f0f0f0", linetype = "dashed"),
    panel.grid.major.x = element_line(color = "#f0f0f0", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 10)),
    axis.text.y = element_text(size = 24, color = "black"),
    axis.text.x = element_text(size = 24, color = "black"),
    axis.ticks = element_line(color = "black"), # Add tick marks
    axis.ticks.length = unit(0.25, "cm"),       # Adjust tick mark length
    legend.position = "right",
    legend.key.size = unit(1, "cm"),  # Larger legend keys
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 20),
    plot.title = element_text(size = 26, face = "bold", hjust = 0.5, margin = margin(b = 15)),
    axis.line = element_line(linewidth = 0.8, color = "black") # Use linewidth for axis lines
  ) +
  labs(
    title = "Significant Top 20 Enriched GO Terms",
    x = expression(bold(-log[10](FDR))),  # Bold mathematical notation
    y = "GO Term",
    color = "Ontology",
    size = "Number of DEGs"  # Updated legend title
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 5)),  # Consistent legend point size
    size = guide_legend(title.position = "top")  # Move size legend title to top
  )

# Save the plot
ggsave("GO_enrichment_with_unidentified_ontology.png", width = 20, height = 14, dpi = 300, bg = "white")

