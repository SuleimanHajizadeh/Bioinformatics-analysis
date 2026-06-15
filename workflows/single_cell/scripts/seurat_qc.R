# 1. Paketləri yükləyirik
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))

# 2. 10x Genomics PBMC 3k Matrisini oxuyuruq
data_dir <- "/bioinformatics/module_3_single_cell/filtered_gene_bc_matrices/hg19/"
pbmc.data <- Read10X(data.dir = data_dir)

# 3. Seurat Obyektini Yaxşı Hüceyrələrlə Yaradırıq
# min.cells = 3 (Gen ən azı 3 hüceyrədə olmalıdır)
# min.features = 200 (Həyatda olan hüceyrə ən azı 200 növ gen ifadə etməlidir)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# 4. Mitoxondrial QC (Keyfiyyət Nəzarəti)
# İnsan genomunda (hg19) mitoxondrial genlər "MT-" ilə başlayır
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# 5. Vizuallaşdırma (Violin Plot)
png("/bioinformatics/module_3_single_cell/pbmc_qc_violin.png", width = 1000, height = 600, res = 120)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
invisible(dev.off())

# 6. Nəticələrin Çıxarılması
cat("Analiz Uğurla Bitdi!\n")
cat("Daxil edilən Hüceyrə sayı (Cells):", ncol(pbmc), "\n")
cat("Daxil edilən Gen sayı (Genes):", nrow(pbmc), "\n")
cat("Violin qrafiki yadda saxlanıldı: /bioinformatics/module_3_single_cell/pbmc_qc_violin.png\n")
