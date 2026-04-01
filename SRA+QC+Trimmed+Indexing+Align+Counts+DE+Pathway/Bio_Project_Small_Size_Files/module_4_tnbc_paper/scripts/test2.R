library(DESeq2)
counts_data <- read.csv("tnbc_gene_counts.csv", row.names = 1, check.names = FALSE)
coldata <- read.csv("tnbc_phenotypes.csv", row.names = 1)
coldata <- coldata[colnames(counts_data), , drop = FALSE]
coldata$Treatment <- factor(coldata$Treatment, levels = c("Normal", "Tumor"))
print("Trying DESeqDataSetFromMatrix")
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = coldata, design = ~ Treatment)
print("Success dds!")
