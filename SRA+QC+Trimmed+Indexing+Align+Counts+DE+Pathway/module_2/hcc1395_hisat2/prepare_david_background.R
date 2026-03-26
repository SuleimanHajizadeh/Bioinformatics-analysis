# 1. Lazımi paketləri yükləmək (əgər yoxdursa: BiocManager::install("org.Hs.eg.db"))
library(DESeq2)

# Əvvəlki analizdən nəticələri oxuyuruq
res_file <- "/bioinformatics/module_2/hcc1395_hisat2/DESeq2_results.csv"
res <- read.csv(res_file, row.names=1)

# ==========================================
# 1. BACKGROUND (FON) SİYAHISININ YARADILMASI
# ==========================================
# DESeq2 analizində iştirak edən, yəni "küy" filtrindən (count>=10) keçən BÜTÜN genlər fon hesab olunur.
background_genes <- rownames(res)

write.table(background_genes, 
            file="/bioinformatics/module_2/hcc1395_hisat2/Background_Genes.txt", 
            row.names=FALSE, col.names=FALSE, quote=FALSE)

cat("Fon (Background) gen siyahısı uğurla yaradıldı: Background_Genes.txt\n")

# ==========================================
# 2. İDENTİFİKATORLARIN ÇEVRİLMƏSİ (Gene Symbol -> Entrez ID)
# ==========================================
# DAVID tez-tez Rəsmi Gen Adlarını (məs. TOB2) tam dəstəkləsə də, ideal analiz üçün
# onları rəqəmsal "Entrez ID"-lərə çevirmək tövsiyə olunur.

# Bioconductor paketi vasitəsilə insan gen məlumat bazasını yükləyirik
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  # Serverdə yüklənməyibsə, xəbərdarlıq edirik: (BiocManager::install("org.Hs.eg.db"))
  cat("Qeyd: İdentifikator çevirməsi üçün 'org.Hs.eg.db' paketi tələb olunur.\n")
} else {
  library(org.Hs.eg.db)
  
  # Yalnız əhəmiyyətli (Significant) genlərimizi oxuyaq
  sig_genes <- read.table("/bioinformatics/module_2/hcc1395_hisat2/Significant_Genes.txt", header=FALSE)$V1
  
  # Symbol -> Entrez ID çevrilməsi
  gene_conversions <- mapIds(org.Hs.eg.db,
                             keys = as.character(sig_genes),
                             column = "ENTREZID",
                             keytype = "SYMBOL",
                             multiVals = "first")
  
  # NA (tapılmayan) dəyərləri təmizləyirik
  entrez_ids <- na.omit(unique(gene_conversions))
  
  write.table(entrez_ids, 
              file="/bioinformatics/module_2/hcc1395_hisat2/Significant_Genes_Entrez.txt", 
              row.names=FALSE, col.names=FALSE, quote=FALSE)
  
  cat("Entrez ID çevrilməsi tamamlandı: Significant_Genes_Entrez.txt\n")
}
