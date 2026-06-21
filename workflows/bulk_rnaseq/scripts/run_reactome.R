# Reactome Pathway Enrichment for TNBC Manuscript
options(bitmapType='cairo')

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos="http://cran.us.r-project.org")
}
if (!requireNamespace("ReactomePA", quietly = TRUE)) {
  BiocManager::install("ReactomePA", update=FALSE, ask=FALSE)
}
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler", update=FALSE, ask=FALSE)
}
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db", update=FALSE, ask=FALSE)
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos="https://cloud.r-project.org")
}

library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# Set directories
results_dir <- "/home/suleimanhajizadeh/Documents/GitHub/Bioinformatics-analysis/workflows/bulk_rnaseq/results"
deg_file <- file.path(results_dir, "significant_tnbc_genes.csv")

print("1. Loading TNBC significant DEGs...")
degs <- read.csv(deg_file)

print("2. Mapping gene symbols to Entrez IDs...")
entrez_ids <- bitr(degs$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

print("3. Executing Reactome Pathway Enrichment Analysis...")
ereact <- enrichPathway(gene = entrez_ids$ENTREZID, 
                        organism = "human", 
                        pvalueCutoff = 0.05, 
                        readable = TRUE)

if (!is.null(ereact) && nrow(ereact@result) > 0) {
  print("4. Saving Reactome pathway enrichment summary...")
  write.csv(as.data.frame(ereact), file.path(results_dir, "Reactome_pathways_summary.csv"), row.names=FALSE)
  
  print("5. Generating Reactome Pathway Dotplot...")
  dot_plot <- dotplot(ereact, showCategory=15, title="Top 15 Enriched Reactome Pathways") +
              theme(axis.text.y = element_text(size = 12, face="bold"))
  
  ggsave(file.path(results_dir, "Reactome_Pathway_Dotplot.png"), plot=dot_plot, width=11, height=8, dpi=300)
  print("SUCCESS: Reactome enrichment results and figure saved successfully!")
} else {
  print("WARNING: No significant Reactome pathways enriched!")
}
