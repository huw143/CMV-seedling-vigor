# Load necessary libraries
library("readxl")
library("dplyr")
library("ggplot2")
library("scales")  # For label formatting

# Read the trinotate.xls file
trinotate_data <- read.csv("filtered_Trinotate_results.csv", header = TRUE)

# Check column names to ensure proper selection of `sprot_Top_BLASTX_hit`
# Replace 'sprot_Top_BLASTX_hit' with the exact column name if it differs
if (!"sprot_Top_BLASTX_hit" %in% colnames(trinotate_data)) {
  stop("The 'sprot_Top_BLASTX_hit' column is missing in the provided file.")
}

# Extract species identifier (e.g., "HORVU") from the `sprot_Top_BLASTX_hit` column
trinotate_data <- trinotate_data %>%
  mutate(
    species = ifelse(
      sprot_Top_BLASTX_hit == ".",
      NA, # Handle missing data
      sub(".*?\\^(.*?)_\\^.*", "\\1", sprot_Top_BLASTX_hit) # Extract species identifier
    )
  )




# Summarize occurrences by species
species_summary <- trinotate_data %>%
  filter(!is.na(species)) %>%
  group_by(species) %>%
  summarize(count = n(), .groups = "drop")


# Extract species identifier (e.g., "CATRO", "PROFR", "FAGES") from the `species` column
species_summary <- species_summary %>%
  mutate(
    species_identifier = ifelse(
      is.na(species),
      NA, # Handle missing data
      sub(".*?_(.*?)\\^.*", "\\1", species) # Extract species identifier after "_"
    )
  )




# Summarize occurrences of each species identifier
species_count <- species_summary %>%
  group_by(species_identifier) %>%
  summarize(count = sum(count), .groups = "drop") %>%
  arrange(desc(count)) %>%  # Sort by count in descending order  (1536 species in total)
  slice_head(n = 26)       # Keep only the top 20 species


str(species_count)

# Assuming `species_count` is your tibble
filtered_species_count <- species_count %>%
  filter(!species_identifier %in% c("ECOLI", "ECO57", "ECO24", "HUMAN", "NEUCR", "ECOL6"))



# Create a lookup table for species identifiers, scientific names, and common names
# Create a lookup table for species identifiers, scientific names, and common names
species_info <- data.frame(
  species_identifier = c("ARATH", "PEA", "ORYSJ", "SOYBN", "MEDTR",
                         "TOBAC", "SOLLC", "SOLTU", "ORYSI", "MEDSA",
                         "LOTJA", "PETHY", "NICPL", "VICFA", "SPIOL",
                         "VITVI","RICCO", "CICAR","MAIZE","PHAVU"),
  scientific_name = c(
    "Arabidopsis thaliana","Pisum sativum", "Oryza sativa japonica",
    "Glycine max","Medicago truncatula", "Nicotiana tabacum", 
    "Solanum lycopersicum", "Solanum tuberosum", "Oryza sativa indica",
    "Medicago sativa", "Lotus japonicus", "Petunia hybrida",
    "Nicotiana plumbaginifolia", "Vicia faba", "Spinacia oleracea",
    "Vitis vinifera", "Ricinus communis", "Cicer arietinum",
    "Zea mays","Phaseolus vulgaris"),
  common_name = c(
    "Thale cress", "Pea", 
    "Rice (Japonica)", "Soybean", 
    "Barrelclover", "Tobacco", 
    "Tomato", "Potato", 
    "Rice (Indica)", "Alfalfa",
    "Japanese Lotus", "Petunia",
    "Wild Tobacco", "Faba bean", 
    "Spinach", "Grape",
    "Castor bean", "Chickpea", 
    "Maize", "Kidney bean")
)

# View the updated table
print(species_info)


# Combine scientific and common names into a single column
species_info <- species_info %>%
  mutate(
    scientific_common_name = paste0(scientific_name, " (", common_name, ")")
  )

# Add scientific and common names to the species_count data
filtered_species_count <- filtered_species_count %>%
  left_join(species_info, by = "species_identifier")




p <- ggplot(filtered_species_count, aes(x = reorder(scientific_common_name, count), y = count, fill = scientific_common_name)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.3) +
  scale_fill_viridis_d(option = "viridis") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(color = "#f0f0f0", linetype = "dashed"),
    panel.grid.major.x = element_line(color = "#f0f0f0", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linetype = "solid"),
    axis.text = element_text(face = "plain", color = "black"),
    axis.title = element_text(face = "bold", color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    axis.ticks = element_line(color = "black")
  ) +
  labs(
    x = "Species (Scientific Name - Common Name)",
    y = "Number of Occurrences",
    title = "Distribution of Top 20 Species"
  ) +
  coord_flip() +
  scale_y_continuous(
    labels = label_number(scale_cut = cut_short_scale()),
    breaks = scales::pretty_breaks(n = 10),
    expand = expansion(mult = c(0, 0.1))  # Reduced left (lower) margin
  ) +
  scale_x_discrete(
    expand = expansion(add = c(0.5, 0.5))
  ) +
  geom_text(aes(label = scales::number(count, big.mark = ",", accuracy = 1)), 
            hjust = -0.2,  
            color = "black", 
            fontface = "bold", 
            size = 3.5) +
  expand_limits(y = 0)  # Ensure y-axis starts exactly at 0


p <- p + theme(plot.margin = margin(1, 1, 1, 1, "cm"))
ggsave("species_count_final.png", 
       plot = p,
       units = "in", 
       width = 12, 
       height = 8,
       bg = "white", 
       dpi = 500)




write.csv(filtered_species_count, "filtered_species_count.csv")
