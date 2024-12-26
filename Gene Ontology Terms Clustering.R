library("simplifyEnrichment")
library("dplyr")
library("GO.db")
library("ggplot2")

# Read the file
go_data <- read.table("High_vs_low.GOseq.enriched", 
                      header = TRUE, 
                      sep = "\t", 
                      stringsAsFactors = FALSE, 
                      fill = TRUE) # Fill missing values

# Check the first few rows to ensure it loaded correctly
head(go_data)

# Remove rows where ontology is NA
go_data <- go_data[!is.na(go_data$ontology), ]


# Prepare the data (filter for significant terms and take the top 20 by FDR)
go_data_filtered <- go_data %>%
  filter(over_represented_FDR < 0.05) %>%
  mutate(`-log10(FDR)` = -log10(over_represented_FDR)) %>%
  arrange(desc(`-log10(FDR)`)) 

str(go_data_filtered)


# Extract the "category" column based on "BP" in "ontology"
bp_category_vector <- go_data_filtered$category[go_data_filtered$ontology == "BP"]
# Convert the extracted data to a vector
bp_category_vector <- as.vector(bp_category_vector)



go_id = bp_category_vector
head(go_id)
mat = GO_similarity(go_id, ont = 'BP')

png("simplifyGO_BP.png", width = 3200, height = 2400, res = 300)
df_BP <- simplifyGO(mat, column_title="Similarity matrix of 170 Biological Process GO terms")
dev.off()
write.csv(df_BP, "df_BP.csv")

# Interactively visualize the clustering results
# cl = binary_cut(mat)
# export_to_shiny_app(mat, cl)




# Extract the "category" column based on "MF" in "ontology"
MF_category_vector <- go_data_filtered$category[go_data_filtered$ontology == "MF"]
# Convert the extracted data to a vector
MF_category_vector <- as.vector(MF_category_vector)

go_id = MF_category_vector
head(go_id)
mat = GO_similarity(go_id, ont = 'MF')

png("simplifyGO_MF.png", width = 3200, height = 2400, res = 300)
df_MF <- simplifyGO(mat, column_title="Similarity matrix of 101 Molecular Function GO terms")
dev.off()
write.csv(df_MF, "df_MF.csv")





# Extract the "category" column based on "MF" in "ontology"
CC_category_vector <- go_data_filtered$category[go_data_filtered$ontology == "CC"]
# Convert the extracted data to a vector
CC_category_vector <- as.vector(CC_category_vector)

go_id = CC_category_vector
head(go_id)
mat = GO_similarity(go_id, ont = 'CC')

png("simplifyGO_CC.png", width = 3200, height = 2400, res = 300)
df_CC <- simplifyGO(mat, column_title="Similarity matrix of 27 Molecular Function GO terms")
dev.off()
write.csv(df_CC, "df_CC.csv")