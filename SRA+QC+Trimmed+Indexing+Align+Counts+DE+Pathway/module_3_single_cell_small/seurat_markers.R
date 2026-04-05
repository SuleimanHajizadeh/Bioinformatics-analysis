# 1. Lazımi paketlərin yüklənməsi
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

# 2. Əvvəlki mərhələdə klasterləşdirilmiş Seurat obyektini açırıq
pbmc <- readRDS("/bioinformatics/module_3_single_cell/pbmc_clustered.rds")

# Yoxlamaq üçün idents (hansı klasterlərə bölündüyünü) təyin edirik
Idents(pbmc) <- "seurat_clusters"

cat("Markerlərin hesablanması (FindAllMarkers) başlayır... Bu bir qədər vaxt apara bilər.\n")

# 3. İSTİFADƏÇİ TƏLƏBİ: Markerlərin Tapılması
# Yalnız o klasterdə dominant olan (ancaq pozitiv, only.pos) genlər:
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# 4. İSTİFADƏÇİ TƏLƏBİ: Top 10 Markerlərin Seçilməsi və Saxlanması (LogFC dərəcəsinə görə)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10, file = "/bioinformatics/module_3_single_cell/pbmc_top10_markers.csv", row.names = FALSE)
cat("Hər klaster üçün Top 10 marker pbmc_top10_markers.csv olaraq yadda saxlanıldı.\n")

# Hər klasterin vizuallaşdırma üçün Top 3 markerini təyin edirik ki, DotPlot qarışıq olmasın
top3 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
distinct_features <- unique(top3$gene)

# 5. İSTİFADƏÇİ TƏLƏBİ: Vizuallaşdırma (DotPlot)
png("/bioinformatics/module_3_single_cell/pbmc_marker_dotplot.png", width = 1000, height = 500, res = 120)
DotPlot(pbmc, features = distinct_features) + 
  RotatedAxis() + 
  ggtitle("PBMC 3k: Klasterlər Üzrə Ən Güclü Marker Genlərin İfadəsi") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
invisible(dev.off())
cat("DotPlot pbmc_marker_dotplot.png olaraq yadda saxlanıldı.\n\n")

# 6. İSTİFADƏÇİ TƏLƏBİ: MSigDB üçün Hər klasterin Ən Güclü 1 Geninin Ekran Çıxarışı
top1 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC) %>% select(cluster, gene, avg_log2FC, p_val_adj)
cat("=== KLASTERLƏRİN ƏN GÜCLÜ (MSigDB-lik) MARKER GENLƏRİ ===\n")
print(as.data.frame(top1))
