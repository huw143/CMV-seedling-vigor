# Load necessary library
library("readxl")  # For reading Excel files
library("dplyr")   # For data manipulation
library("stringr")             # Load the package

# Step 1: Load the database
trinotate_data <- read.delim("myTrinotate.xls", header = TRUE)  # For TSV

# Step 2: Extract relevant fields from the 'sprot_Top_BLASTX_hit' column
# Assuming 'sprot_Top_BLASTX_hit' is a string formatted as 'IBB_HORVU^...^66.7%ID^E:6.88e-42^...'
trinotate_data <- trinotate_data %>%
  mutate(
    E_value = as.numeric(sub("E:", "", str_extract(sprot_Top_BLASTX_hit, "E:[^\\^]+"))),
    Identity = as.numeric(sub("%ID", "", str_extract(sprot_Top_BLASTX_hit, "[0-9.]+%ID"))),
    Coverage = as.numeric(sub("%", "", str_extract(sprot_Top_BLASTX_hit, "[0-9.]+(?=%)")))
  )

# Step 3: Apply filtering criteria
filtered_data <- trinotate_data %>%
  filter(
    E_value <= 1e-10,  # Stricter threshold
    Identity >= 50,    # Or change to 30 for distant homologs
    Coverage >= 70     # Query or subject coverage ≥ 70%
  )

# Step 4: Save or view the filtered database
# View(filtered_data)  # View the filtered data in RStudio
write.csv(filtered_data, "filtered_Trinotate_results.csv", row.names = FALSE)  # Save as CSV
