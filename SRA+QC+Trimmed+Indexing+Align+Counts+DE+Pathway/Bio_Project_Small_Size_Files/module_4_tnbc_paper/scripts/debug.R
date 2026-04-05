options(bitmapType='cairo')
library(DESeq2)
setwd("/bioinformatics/module_4_tnbc_paper/")
counts_data <- read.csv("tnbc_gene_counts.csv", row.names = 1, check.names = FALSE)
coldata <- read.csv("tnbc_phenotypes.csv", row.names = 1)
coldata <- coldata[colnames(counts_data), , drop = FALSE]
coldata$Treatment <- factor(coldata$Treatment, levels = c("Normal", "Tumor"))
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = coldata, design = ~ Treatment)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
cat("Starting DESeq...\n")
dds <- DESeq(dds)
cat("DESeq Done!\n")
res <- results(dds, contrast=c("Treatment", "Tumor", "Normal"))
cat("results Done!\n")
write.csv(as.data.frame(res), file="debug_test.csv")
cat("write Done!\n")
