# TNBC Mərhələ 5: DESeq2 Pipeline
options(bitmapType='cairo')
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

setwd("/bioinformatics/module_4_tnbc_paper/")

print("1. Məlumatlar yüklənir...")
counts_data <- read.csv("tnbc_gene_counts.csv", row.names = 1, check.names = FALSE)
coldata <- read.csv("tnbc_phenotypes.csv", row.names = 1)

# Matris sütunları ilə Xəstə Fenotiplərini güzgü kimi simmetriyaya gətiririk
coldata <- coldata[colnames(counts_data), , drop = FALSE]

# Referans (Baseline) səviyyəsini 'Normal' təyin edirik
coldata$Treatment <- factor(coldata$Treatment, levels = c("Normal", "Tumor"))

print("2. DESeqDataSet qurulur və Ölü genlər (low counts) filtrdən keçirilir...")
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = coldata,
                              design = ~ Treatment)

# Ən azı 10 oxunuşu olmayan yatmış genləri silirik ki, riyazi model təmizlənsin
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

print("3. DESeq2 modeli hesablayır (P-Adj, Log2FoldChange)...")
dds <- DESeq(dds)

# Nəticələri Tumor-Normal qütbləşməsinə görə çıxarırıq
res <- results(dds, contrast=c("Treatment", "Tumor", "Normal"))
resOrdered <- res[order(res$padj),]

print("4. Nəticə Matrisi (tnbc_deg_results.csv) yaddaşa yazılır...")
write.csv(as.data.frame(resOrdered), file="tnbc_deg_results.csv")

print("5. PCA Qrafiki (Principal Component Analysis) şifrəsi partladılır...")
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Treatment"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaPlot <- ggplot(pcaData, aes(PC1, PC2, color=Treatment)) +
  geom_point(size=6, alpha=0.8, stroke=1, color="black", shape=21, aes(fill=Treatment)) +
  theme_minimal() +
  labs(title="Principal Component Analysis (TNBC vs Normal)",
       x=paste0("PC1: ",percentVar[1],"% variance"),
       y=paste0("PC2: ",percentVar[2],"% variance")) +
  scale_fill_manual(values=c("Normal"="#28a745", "Tumor"="#dc3545")) +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=16),
        axis.title = element_text(face="bold", size=14),
        legend.title = element_blank(),
        legend.text = element_text(size=12))

ggsave("PCA_plot.png", plot=pcaPlot, width=8, height=6, dpi=300)

print("6. Heatmap Xəritəsi (Ən Kritik 50 Marker Gen) vizuallaşdırılır...")
top_genes <- rownames(head(resOrdered, 50))
mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat) 

annotation_col <- data.frame(Treatment=coldata$Treatment)
rownames(annotation_col) <- rownames(coldata)
ann_colors <- list(Treatment = c(Normal = "#28a745", Tumor = "#dc3545"))
my_colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

png("Heatmap_plot.png", width=10, height=8, res=300, units="in")
pheatmap(mat, 
         annotation_col = annotation_col, 
         annotation_colors = ann_colors,
         color = my_colors,
         show_rownames = TRUE,
         scale="row",
         fontsize_row=7,
         main="Top 50 Differentially Expressed Genes")
dev.off()

print("====== XƏRÇƏNGİN ŞİFRƏSİ QIRILDI! MƏRHƏLƏ 5 DÜYMƏSİ TAMAMLANDI! ======")
