# R-scripts-porcine

# =========Install the libraries ==============

BiocManager::install("org.Ss.eg.db")

# Load packages
library(clusterProfiler)
library(org.Ss.eg.db) # If you are using pig
library(enrichplot)
library(GOSemSim)
library(DOSE)
library(ggplot2)

# Your symbol gene list
gene_symbols <- c("gene name") # add the gene names

# Remove duplicates
gene_symbols <- unique(gene_symbols)

# 1. Convert to ENTREZID
gene_df <- bitr(gene_symbols,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Ss.eg.db)

# Check for successful conversions
print(gene_df)


#2. GO (Biological Process) Enrichment
ego <- enrichGO(gene          = gene_df$ENTREZID,
                OrgDb         = org.Ss.eg.db,
                ont           = "BP",
                pAdjustMethod = "fdr",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2,
                readable      = TRUE)

# Check results
head(ego)

# 3. Create a semantic network with GOSemSim
# Generate semantic data from the pig

ssGO <- godata('org.Ss.eg.db', ont = "BP")  # If it gives an error, you can try with ont = "BP"

# Generate semantic similarity matrix
ego <- pairwise_termsim(ego, semData = ssGO)

# 4. Network type map (emapplot)
emapplot(ego, showCategory = 20, layout = "kk")

ggsave("GO_network_strength-legs.png", width = 8, height = 6, dpi = 300)

#Save the results to a png file
ggsave("name.png", width = 8, height = 6, dpi = 300)

# Save the results to a CSV file
write.csv(ego@result, "GO_name_result.csv", row.names = FALSE)


