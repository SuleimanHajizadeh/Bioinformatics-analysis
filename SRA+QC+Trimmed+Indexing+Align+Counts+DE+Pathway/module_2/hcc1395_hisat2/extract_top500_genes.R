# DESeq2 nəticələrini oxuyaq
res_file <- "/bioinformatics/module_2/hcc1395_hisat2/DESeq2_results.csv"
res <- read.csv(res_file, row.names=1)

# Genləri sadəcə p-dəyərinə (pvalue) görə ən vacibdən başlayaraq sıralayırıq
res_sorted <- res[order(res$pvalue), ]

# Həmçinin bir qədər artma/azalma (log2FoldChange) fərqi olanları seçək
res_filtered <- res_sorted[abs(res_sorted$log2FoldChange) > 0.5, ]

# Siyahıdan ən əhəmiyyətli Top 500 geni götürək (DAVID üçün optimal rəqəm)
top500 <- head(rownames(res_filtered), 500)

write.table(top500, file="/bioinformatics/module_2/hcc1395_hisat2/Significant_Genes.txt", 
            row.names=FALSE, col.names=FALSE, quote=FALSE)

cat("500 Genlik geniş siyahı uğurla Significant_Genes.txt faylına yazıldı!\n")

# Həmçinin Entrez ID versiyasını yeniləyək
suppressMessages(library(org.Hs.eg.db))
gene_conversions <- mapIds(org.Hs.eg.db, keys = as.character(top500), column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
entrez_ids <- na.omit(unique(gene_conversions))

write.table(entrez_ids, file="/bioinformatics/module_2/hcc1395_hisat2/Significant_Genes_Entrez.txt", 
            row.names=FALSE, col.names=FALSE, quote=FALSE)
cat("Entrez ID-lər də yeniləndi.\n")
