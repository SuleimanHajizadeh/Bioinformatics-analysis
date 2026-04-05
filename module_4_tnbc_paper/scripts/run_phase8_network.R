# TNBC Phase 8: Systems Biology & Co-expression Network Inference
options(bitmapType='cairo')

if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph", repos="http://cran.us.r-project.org")
if (!requireNamespace("ggraph", quietly = TRUE)) install.packages("ggraph", repos="http://cran.us.r-project.org")
if (!requireNamespace("tidygraph", quietly = TRUE)) install.packages("tidygraph", repos="http://cran.us.r-project.org")

library(DESeq2)
library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)

setwd("/bioinformatics/module_4_tnbc_paper/")

print("1. Loading Transcriptomic Matrices for Network Inference...")
counts_data <- read.csv("tnbc_gene_counts.csv", row.names = 1, check.names = FALSE)
coldata <- read.csv("tnbc_phenotypes.csv", row.names = 1)
coldata <- coldata[colnames(counts_data), , drop = FALSE]
coldata$Treatment <- factor(coldata$Treatment, levels = c("Normal", "Tumor"))

dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = coldata, design = ~ Treatment)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
vsd <- vst(dds, blind=FALSE)
vsd_mat <- assay(vsd)

print("2. Isolating Top Structural Genes for the Network...")
res <- read.csv("tnbc_deg_results.csv", row.names = 1)
res_ord <- res[order(res$padj), ]
sig_genes_top <- rownames(head(res_ord[!is.na(res_ord$padj) & abs(res_ord$log2FoldChange) > 2, ], 50))

print("3. Computing Pearson Mathematical Correlations...")
# Create a correlation matrix among our top 50 DEGs
cor_matrix <- cor(t(vsd_mat[sig_genes_top, ]), method="pearson")

# Keep only extremely strong interactions (R > 0.85 or R < -0.85)
cor_matrix[abs(cor_matrix) < 0.85] <- 0
diag(cor_matrix) <- 0

print("4. Building Topological Spider-Web (iGraph)...")
graph <- graph_from_adjacency_matrix(cor_matrix, mode="undirected", weighted=TRUE)

# We remove isolated "lonely" genes that don't talk to anyone
isolated <- which(degree(graph) == 0)
graph <- delete_vertices(graph, isolated)

# Calculate Network Centrality (Degree) to find the "Commander" Hub Genes
V(graph)$degree <- degree(graph)
V(graph)$label <- V(graph)$name

print("5. Generating Systems Biology Network Graph...")
png("Systems_Biology_Network.png", width=10, height=8, res=300, units="in")

set.seed(123)
# Beautiful Spider Web layout
ggraph(graph, layout = "fr") + 
  geom_edge_link(aes(edge_alpha = abs(weight)), edge_colour = "gray50") + 
  geom_node_point(aes(size = degree, color = degree)) + 
  geom_node_text(aes(label = label), repel = TRUE, size = 5, fontface = "bold") +
  scale_color_gradient(low = "lightblue", high = "darkred") +
  theme_void() +
  labs(title = "Systems Biology: TNBC Hub-Gene Co-expression Network",
       subtitle = "Node size/color explicitly indicates the Master Regulator status (Hub Centrality).") +
  theme(plot.title = element_text(face="bold", size=18, hjust = 0.5),
        plot.subtitle = element_text(size=14, hjust = 0.5),
        legend.position = "none")

dev.off()

print("====== PHASE 8: NETWORK INFERENCE COMPLETE! HUB GENES FOUND! ======")
