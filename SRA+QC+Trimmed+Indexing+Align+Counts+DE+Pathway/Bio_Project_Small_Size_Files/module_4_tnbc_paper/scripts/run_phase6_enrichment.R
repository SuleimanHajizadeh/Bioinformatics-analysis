options(bitmapType='cairo')
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)

setwd("/bioinformatics/module_4_tnbc_paper/")

print("1. DESeq2 Nəticələri Yüklənir...")
res <- read.csv("tnbc_deg_results.csv", row.names = 1)
res$Gene <- rownames(res)

print("2. Volcano Qrafiki Üçün Qütbləşmə Filtrləri Qoyulur...")
# add a column for significance classification
res$Significance <- "Not Significant"
res$Significance[res$log2FoldChange > 2.0 & res$padj < 0.05] <- "Up-regulated"
res$Significance[res$log2FoldChange < -2.0 & res$padj < 0.05] <- "Down-regulated"

# specific top genes labeling (Requested by user: AGO2, ZFAS1 etc.)
top_genes <- res[order(res$padj), ]
label_genes <- head(top_genes$Gene, 20)
res$Label <- ifelse(res$Gene %in% label_genes | res$Gene %in% c("AGO2", "ZFAS1"), res$Gene, "")

print("3. Estetik Volcano Plot (Magma Partlayışı) Qrafiki Yaradılır...")
volcano_plot <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.6) +
  scale_color_manual(values = c("Down-regulated" = "#2c7bb6", 
                                "Not Significant" = "grey", 
                                "Up-regulated" = "#d7191c")) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(aes(label = Label), 
                  size = 4.5, 
                  box.padding = 0.6, 
                  max.overlaps = 30, 
                  fontface = "bold.italic",
                  color = "black") +
  theme_minimal() +
  labs(title = "Volcano Plot (TNBC vs Normal)",
       x = "Log2 Fold Change",
       y = "-Log10(adj. p-value)") +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.title = element_text(face="bold", size=14))

ggsave("Volcano_plot.png", plot=volcano_plot, width=9, height=8, dpi=300)

print("4. Mütləq Gen Namizədlərinin Seçilməsi...")
sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 2.0)

if(nrow(sig_genes) == 0) {
  stop("CATASTROPHIC ERROR! Səhv-Nəticə Nəzarəti: Heç bir gen kriteryalara girmədi. Model yanlışdır!")
}

print(sprintf("UĞUR! Ümumi %d Təhlükəli Süzülmüş Kəskin Gen Tapıldı!", nrow(sig_genes)))
write.csv(sig_genes, "significant_tnbc_genes.csv", row.names=FALSE)

print("5. Xərçəng Mexanizminin Öyrənilməsi (GO Pathway Enrichment)...")
# Convert official gene symbols to Entrez IDs for clusterProfiler
gene_list <- sig_genes$Gene
entrez_ids <- bitr(gene_list, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

print("   - GO: Bioloji Proses Analizi Başlayır...")
ego <- enrichGO(gene          = entrez_ids$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP", # Biological Process
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2,
                readable      = TRUE)

if (!is.null(ego) && nrow(ego@result) > 0) {
  print("6. Xəstəliyin Mexanik Mexanizm Qrafiki (DotPlot) Yaradılır...")
  dot_plot <- dotplot(ego, showCategory=15, title="Top 15 Enriched GO Pathways") +
     theme(axis.text.y = element_text(size = 12, face="bold"))
  ggsave("GO_Pathway_Dotplot.png", plot=dot_plot, width=11, height=8, dpi=300)
  write.csv(as.data.frame(ego), "GO_pathways_summary.csv", row.names=FALSE)
} else {
  print("XƏBƏRDARLIQ: Heç bir GO yolu tapılmadı və ya qırılmalar lokal səviyyədədir!")
}

print("====== MƏRHƏLƏ 6 BİTİB ÇATDIRILDI! MƏQALƏ TƏMƏLİ QURAŞDIRILIR! ======")
