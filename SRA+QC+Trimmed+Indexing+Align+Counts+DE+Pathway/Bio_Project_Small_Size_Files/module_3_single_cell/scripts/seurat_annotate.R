# 1. Lazımi paketləri yükləyirik
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))

# 2. Klasterləşmiş və Markerləri hesablanmış obyekti açırıq
pbmc <- readRDS("/bioinformatics/module_3_single_cell/pbmc_clustered.rds")

# Yoxlama üçün mövcud kimlikləri (0-dan 8-ə) təyin edirik
Idents(pbmc) <- "seurat_clusters"

# 3. İSTİFADƏÇİ TƏLƏBİ: Yeni Bioloji Adların (Cell Types) Təyini
# Qeyd: Bu ardıcıllıq klaster nömrələri olan 0, 1, 2, ..., 8 sırasına uyğun olmalıdır
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B cells", 
                     "CD8 T", "FCGR3A+ Mono", "NK cells", "Dendritic cells", "Platelets")

names(new.cluster.ids) <- levels(pbmc)

# Obyekti yeni adlarla "Rename" (Adlandır) edirik
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# 4. İSTİFADƏÇİ TƏLƏBİ: Yekun Vizuallaşdırma (UMAP)
png("/bioinformatics/module_3_single_cell/pbmc_final_annotation.png", width = 1000, height = 600, res = 120)
# 'label = TRUE' qrafikin üstünə adları yazacaq
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 4.5) + 
  ggtitle("PBMC 3k: Bioloji Klasterlərin İdentifikasiyası (UMAP)") +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)) +
  NoLegend() # Rənglər onsuzda yazılarla eynidir deyə qıraq cədvəlini ləğv edirik
invisible(dev.off())

# 5. İSTİFADƏÇİ TƏLƏBİ: Son Obyektin Saxlanılması
saveRDS(pbmc, file = "/bioinformatics/module_3_single_cell/pbmc_final.rds")

cat("Hüceyrələr uğurla rəsmi MSigDB adları ilə adlandırıldı!\n")
cat("Yekun qrafik yadda saxlanıldı: /bioinformatics/module_3_single_cell/pbmc_final_annotation.png\n")
cat("Tamamlanmış R Obyekti gələcək araşdırmalar üçün yadda saxlanıldı: pbmc_final.rds\n")
