# 1. Paketlərin Yüklənməsi
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))

# 2. Matrisin oxunması və İlkin Obyektin yaradılması
data_dir <- "/bioinformatics/module_3_single_cell/filtered_gene_bc_matrices/hg19/"
pbmc.data <- Read10X(data.dir = data_dir)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# Mühüm: Mitoxondrial faizin hesablanması
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# 3. İSTİFADƏÇİ TƏLƏBİ: Hüceyrə Filtrləməsi (Visual QC-yə əsasən)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# 4. İSTİFADƏÇİ TƏLƏBİ: Normalizasiya və Dəyişkən Genlərin Tapılması
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# 5. İSTİFADƏÇİ TƏLƏBİ: Ölçülərin Endirilməsi və Klasterləşdirmə
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# PCA (Principal Component Analysis)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Klasterlərin tapılması (İlk 10 Dimensiya əsasında)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# UMAP (Uniform Manifold Approximation and Projection)
pbmc <- RunUMAP(pbmc, dims = 1:10)

# 6. İSTİFADƏÇİ TƏLƏBİ: Nəticə və Vizuallaşdırma
png("/bioinformatics/module_3_single_cell/pbmc_umap_islands.png", width = 800, height = 600, res=120)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("PBMC 3k - Hüceyrə Adaları (UMAP)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
invisible(dev.off())

# Seurat obyektinin gələcək annotasiya üçün yadda saxlanılması
saveRDS(pbmc, file = "/bioinformatics/module_3_single_cell/pbmc_clustered.rds")

# Statistik hesabat
cat("\n=== UMAP ANALİZİ TAMAMLANDI ===\n")
cat("Filtrdən sonra qalan Sağlam Hüceyrələr:", ncol(pbmc), "\n")
cat("Aşkar edilən Spesifik Genlər:", nrow(pbmc), "\n")
cat("Yaradılan Hüceyrə Qrupları (Klasterlər):", length(levels(Idents(pbmc))), "\n")
cat("UMAP faylı saxlanıldı: /bioinformatics/module_3_single_cell/pbmc_umap_islands.png\n")
