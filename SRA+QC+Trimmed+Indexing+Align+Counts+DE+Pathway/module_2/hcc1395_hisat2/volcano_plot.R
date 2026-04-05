# 1. Paketləri yüklə
library(ggplot2)

# 2. DESeq2 Nəticələrini Oxu
res_file <- "/bioinformatics/module_2/hcc1395_hisat2/DESeq2_results.csv"
res <- read.csv(res_file, row.names=1)

# NA olanları 1-ə bərabərlə (xətanın qarşısını almaq üçün)
res$padj[is.na(res$padj)] <- 1

# 3. Volcano Plot Üçün Rəngləndirmə (Significance Thresholds)
# padj < 0.05 və |log2FoldChange| > 1 olanlar "Significant" qəbul edilir
res$Significant <- "NO"
res$Significant[res$padj < 0.05 & res$log2FoldChange > 1] <- "UP"
res$Significant[res$padj < 0.05 & res$log2FoldChange < -1] <- "DOWN"

# 4. Volcano Plot-un Çəkilməsi
p <- ggplot(data=res, aes(x=log2FoldChange, y=-log10(padj), color=Significant)) +
  geom_point(alpha=0.8, size=3) +
  scale_color_manual(values=c("DOWN"="blue", "NO"="grey", "UP"="red")) +
  theme_minimal() +
  labs(title="Volcano Plot (Tumor vs Normal)",
       x="log2 Fold Change",
       y="-log10(Adjusted p-value)") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="black") +
  geom_vline(xintercept=c(-1, 1), linetype="dashed", color="black")

ggsave("/bioinformatics/module_2/hcc1395_hisat2/Volcano_Plot.png", plot=p, width=7, height=6)

# 5. Əhəmiyyətli (Significant) Genləri DAVID/Reactome üçün Təmiz Mətnə Yazmaq
# "NO" olmayan (yəni UP və DOWN) genlərin adlarını (rownames) çıxarırıq
sig_genes <- rownames(res[res$Significant != "NO", ])

# Əgər faylda heç bir əhəmiyyətli gen yoxdursa belə (bu kiçik subset olduğu üçün), siyahını yaradırıq.
if(length(sig_genes) == 0) {
  # Nümunə naminə Top 10 geni yazaq ki, boş qalmasın
  sig_genes <- head(rownames(res[order(res$pvalue), ]), 10)
}

write.table(sig_genes, file="/bioinformatics/module_2/hcc1395_hisat2/Significant_Genes.txt", 
            row.names=FALSE, col.names=FALSE, quote=FALSE)

cat("Volcano Plot və Gen siyahısı uğurla yaradıldı!\n")
