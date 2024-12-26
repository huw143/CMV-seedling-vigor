# This analysis used 3322 DEGs

# KEGG Pathway Analysis for Non-Model Organism RNA-seq
# Required Libraries
library("clusterProfiler")
library("dplyr")
library("ggplot2")
library("pathview")
library("readxl")
library("KEGGREST")
library("tidyr")


# Read the Trinotate annotation file
# Assumes the file is a tab-separated file with KEGG annotations
# Modify the column names as per your specific Trinotate output
kegg_annotation <- read.delim("myTrinotate.xls", header = TRUE, stringsAsFactors = FALSE)

# load DEGs
differential_expression_results <- read_excel("differential_expression_results.xlsx")

# check the structure of your dataframes
str(differential_expression_results)
str(kegg_annotation)


# Create a mapping dataframe
mapping_df <- data.frame(
  gene_id = differential_expression_results$ID,  # Use the gene ID column from DE results
  log2FoldChange = differential_expression_results$log2FoldChange,
  kegg_ko = kegg_annotation$EggNM.KEGG_ko[match(differential_expression_results$ID, kegg_annotation$X.gene_id)]  # Adjust column names as needed
)

mapping_df <- mapping_df %>%
  mutate(kegg_ko = strsplit(kegg_ko, ",")) %>%  # Split multiple KOs into lists
  unnest(kegg_ko)  # Create one row per KO



# Remove rows without KEGG annotations
mapping_df <- mapping_df[!is.na(mapping_df$kegg_ko), ]


# Create log2 fold change vector with KEGG IDs
log2_fold_change <- setNames(
  mapping_df$log2FoldChange,
  mapping_df$kegg_ko
)

# Clean up KO IDs (remove 'ko:' prefix, ensure correct format)
names(log2_fold_change) <- gsub("ko:", "", names(log2_fold_change))
names(log2_fold_change) <- gsub("K0*", "K", names(log2_fold_change))

# Filter to ensure only valid K-number IDs
log2_fold_change <- log2_fold_change[grepl("^K\\d+$", names(log2_fold_change))]




# KEGG Pathway Enrichment Analysis:
# Perform KEGG pathway enrichment
kegg_enrichment <- enrichKEGG(
  gene = names(log2_fold_change),  # Use the KEGG IDs
  organism = "ko",  # For non-model organisms
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)



# Visualize enrichment results
# Barplot of enriched pathways
kegg_bar <- barplot(kegg_enrichment, 
                    showCategory = 20, 
                    title = "Top 20 Enriched KEGG Pathways")

# Dotplot of enriched pathways
kegg_dot <- dotplot(kegg_enrichment, 
                    showCategory = 20, 
                    title = "Top 20 Enriched KEGG Pathways")




# Save plots
ggsave("kegg_barplot.png", kegg_bar, width = 10, height = 8)
ggsave("kegg_dotplot.png", kegg_dot, width = 10, height = 8)

# Save enrichment results
write.csv(kegg_enrichment@result, "kegg_enrichment_results.csv")

print(kegg_bar)
print(kegg_dot)




# Optional: Pathway view for top pathways
# You'll need the KEGG pathway ID and gene expression data
# This is an example - adjust according to your specific data
top_pathways <- kegg_enrichment@result$ID[1:10]  # Top 5 pathways

for(pathway in top_pathways) {
  tryCatch({
    pathview(
      gene.data = log2_fold_change,
      pathway.id = pathway,
      species = "ko",
      out.suffix = paste0("pathway_", pathway)
    )
  }, error = function(e) {
    print(paste("Could not generate pathway for", pathway))
  })
}

