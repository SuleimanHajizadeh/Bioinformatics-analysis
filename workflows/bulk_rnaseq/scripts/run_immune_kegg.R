options(bitmapType='cairo')
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

setwd("/bioinformatics/module_4_tnbc_paper/")

print("1. Data və VST Yüklənir...")
counts_data <- read.csv("tnbc_gene_counts.csv", row.names = 1, check.names = FALSE)
coldata <- read.csv("tnbc_phenotypes.csv", row.names = 1)
coldata <- coldata[colnames(counts_data), , drop = FALSE]
coldata$Treatment <- factor(coldata$Treatment, levels = c("Normal", "Tumor"))

dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = coldata, design = ~ Treatment)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
vsd <- vst(dds, blind=FALSE)

print("2. İmmun T-Hüceyrə Profilinin Çıxarılması...")
immune_markers <- c("CD3D", "CD3E", "CD4", "CD8A", "CD8B", "FOXP3", "IFNG", "PRF1", 
                    "GZMB", "PDCD1", "CTLA4", "LAG3", "TIGIT", "CD163", "MS4A1", "CD19")

available_markers <- intersect(immune_markers, rownames(vsd))
immune_mat <- assay(vsd)[available_markers, ]
immune_mat <- immune_mat - rowMeans(immune_mat)

annotation_col <- data.frame(Treatment=coldata$Treatment)
rownames(annotation_col) <- rownames(coldata)
ann_colors <- list(Treatment = c(Normal = "#28a745", Tumor = "#dc3545"))
my_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)

print("3. İmmun Xəritəsinin (Heatmap) Dizaynı...")
png("Immune_Landscape_Heatmap.png", width=8, height=6, res=300, units="in")
pheatmap(immune_mat, 
         annotation_col = annotation_col, 
         annotation_colors = ann_colors,
         color = my_colors,
         show_rownames = TRUE,
         scale="row",
         cluster_cols = TRUE,
         fontsize_row=11,
         main="Immune Landscape: Tumor vs Normal (Hot/Cold)")
dev.off()

print("4. KEGG Pathway Analizi İcra Edilir...")
sig_genes <- read.csv("significant_tnbc_genes.csv")
gene_list <- sig_genes$Gene
entrez_ids <- bitr(gene_list, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

ekegg <- enrichKEGG(gene         = entrez_ids$ENTREZID,
                    organism     = 'hsa',
                    pvalueCutoff = 0.05)

if (!is.null(ekegg) && nrow(ekegg@result) > 0) {
  ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  
  dot_plot <- dotplot(ekegg, showCategory=15, title="Top 15 Enriched KEGG Pathways") +
     theme(axis.text.y = element_text(size = 12, face="bold"))
  
  ggsave("KEGG_Pathway_Dotplot.png", plot=dot_plot, width=10, height=8, dpi=300)
  write.csv(as.data.frame(ekegg), "KEGG_pathways_summary.csv", row.names=FALSE)
  print("KEGG Şəkli Uğurla Yazıldı!")
} else {
  print("XƏBƏRDARLIQ: Heç bir KEGG yolu tapılmadı!")
}

print("====== KEGG VƏ İMMUN PROFİL TAMAMLANDI! MƏQALƏ YAZILA BİLƏR! ======")
