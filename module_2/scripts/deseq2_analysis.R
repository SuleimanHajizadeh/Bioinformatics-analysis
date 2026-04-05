# 1. DESeq2 və Vizuallaşdırma Paketlərinin Yüklənməsi
library(DESeq2)
library(ggplot2)
library(pheatmap)

# 2. Faylların Oxunması və Təmizlənməsi
counts_file <- "/bioinformatics/module_2/hcc1395_hisat2/final_counts.csv"
meta_file   <- "/bioinformatics/module_2/hcc1395_hisat2/metadata.csv"

# Count matrix yüklənib sətir adları (genes) təyin edilir
counts_data <- read.csv(counts_file, row.names = 1)
# Sütun adları qısa və aydın edək (featureCounts ".bam" sözlərini əlavə edir deyə)
colnames(counts_data) <- c("SRR1553607", "SRR1553606", "SRR1553608")

# Metadata yüklənir
col_data <- read.csv(meta_file, row.names = 1)

# XƏBƏRDARLIQ: DataFrame-dəki sütunlarla Metadata sətirlərinin eyni sıra ilə olduğundan əmin olmaq!
col_data <- col_data[match(colnames(counts_data), rownames(col_data)), , drop = FALSE]

# 3. DESeqDataSet Obyektinin Yaradılması
# "Tumor vs Normal" müqayisəsi etdiyimiz üçün "Condition" dizayn kimi qeyd edilir
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = col_data,
                              design = ~ Condition)

# 4. Filtrləmə (Filtering):
# Çox aşağı ifadələnən və ya sıfır "count"-a sahib genləri xaric edirik (Bütün xəstələrdə cəmi piki < 10 olanlar silinir)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# 5. Normalizasiya və Variasiyanın Sabitləşdirilməsi (VST)
# PCA və Məsafə qrafikləri çəkmək üçün məlumatları median normalization edirik
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

# ==========================================
# 6. Keyfiyyət Nəzarəti (QC) - QRAFİKLƏR
# ==========================================

# A) PCA (Principal Components Analysis) Plot
pca_plot <- plotPCA(vsd, intgroup="Condition") +
  theme_minimal() +
  ggtitle("RNA-Seq Nümunələrinin PCA Qrafiki (Tumor vs Normal)")

ggsave("/bioinformatics/module_2/hcc1395_hisat2/QC_PCA_plot.png", plot=pca_plot, width=7, height=5)

# B) Nümunələr Arası Məsafə (Sample Distance Heatmap)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- NULL

png("/bioinformatics/module_2/hcc1395_hisat2/QC_Distance_Heatmap.png", width=600, height=500)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         main = "Nümunələr Arasındakı Məsafə (Daha tünd rəng = Daha oxşar)")
dev.off()

# ==========================================
# 7. DESeq2 Analizinin Özü (Differential Expression)
# ==========================================
dds <- DESeq(dds)

# Tumorları Normallara qarşı müqayisə edirik
res <- results(dds, contrast=c("Condition", "Tumor", "Normal"))

# Nəticələri padj (p-value adjusted) yəni etibarlılığa görə sıralayırıq
res_ordered <- res[order(res$padj),]

# Nəticəni CSV kimi yaddaşa yazırıq
write.csv(as.data.frame(res_ordered), "/bioinformatics/module_2/hcc1395_hisat2/DESeq2_results.csv")

# Terminal ekranında ən güclü və etibarlı (Top 10) diferensial genləri göstər
cat("\n=== ƏN ÇOX FƏRQLƏNƏN 10 GEN (Tumor vs Normal) ===\n")
print(head(res_ordered, 10))
