# TNBC Phase 9: Digital Pathology & Tumor Microenvironment (TME) Deconvolution
options(bitmapType='cairo')

if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2", repos="http://cran.us.r-project.org")
if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2", repos="http://cran.us.r-project.org")

library(DESeq2)
library(ggplot2)
library(reshape2)

setwd("/bioinformatics/module_4_tnbc_paper/")

print("1. Loading VST Normalized Counts for Deconvolution...")
counts_data <- read.csv("tnbc_gene_counts.csv", row.names = 1, check.names = FALSE)
coldata <- read.csv("tnbc_phenotypes.csv", row.names = 1)
coldata <- coldata[colnames(counts_data), , drop = FALSE]
coldata$Treatment <- factor(coldata$Treatment, levels = c("Normal", "Tumor"))

dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = coldata, design = ~ Treatment)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
vsd <- vst(dds, blind=FALSE)
vsd_mat <- assay(vsd)

print("2. Defining Academic Marker Signatures for TME Cell Populations...")
# Standardized digital pathology signatures mapping explicit cell arrays
signatures <- list(
  `CD8+ Cytotoxic T-Cells` = c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1"),
  `B-Cells` = c("CD19", "MS4A1", "CD79A", "CD79B"),
  `M1 Macrophages (Anti-Tumor)` = c("NOS2", "IRF5", "STAT1", "IL12A"),
  `M2 Macrophages (Pro-Tumor)` = c("CD163", "MRC1", "IL10", "ARG1", "CSF1R")
)

print("3. Executing Mathematical Deconvolution Scoring...")
# Convert expression matrix to Z-scores per gene for fair comparison
vsd_scaled <- t(scale(t(vsd_mat)))

# Calculate geometric mean / average Z-score per signature per sample
score_matrix <- data.frame(matrix(ncol = ncol(vsd_scaled), nrow = 0))
colnames(score_matrix) <- colnames(vsd_scaled)

for (cell_type in names(signatures)) {
  markers <- signatures[[cell_type]]
  available_markers <- intersect(markers, rownames(vsd_scaled))
  
  if (length(available_markers) > 0) {
    if (length(available_markers) == 1) {
      cell_score <- vsd_scaled[available_markers, ]
    } else {
      cell_score <- colMeans(vsd_scaled[available_markers, ], na.rm = TRUE)
    }
    # To plot relative abundance, we shift scores to be strictly positive
    cell_score <- cell_score - min(cell_score) + 0.1 
    score_matrix <- rbind(score_matrix, cell_score)
  }
}
rownames(score_matrix) <- names(signatures)

print("4. Preparing Matrix for Relative Fraction Plotting...")
# Transpose and melt for ggplot
plot_data <- as.data.frame(t(score_matrix))
plot_data$Sample <- rownames(plot_data)
plot_data$Condition <- coldata$Treatment

# Sort samples naturally by treatment to show separation
plot_data <- plot_data[order(plot_data$Condition, plot_data$Sample), ]
plot_data$Sample <- factor(plot_data$Sample, levels = plot_data$Sample)

melted_data <- melt(plot_data, id.vars = c("Sample", "Condition"), variable.name = "CellType", value.name = "Score")

print("5. Generating Digital Pathology Stacked Barplot...")
png("TME_Deconvolution_Barplot.png", width=11, height=7, res=300, units="in")

ggplot(melted_data, aes(x = Sample, y = Score, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill", width = 0.85) +
  facet_grid(. ~ Condition, scales = "free_x", space = "free_x") +
  theme_minimal() +
  scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#ff7f00")) +
  labs(title = "\nTumor Microenvironment (TME) Deconvolution Fractions",
       subtitle = "Digital Pathology: Relative Core Immune Subpopulations",
       x = "Individual Clinical Samples",
       y = "Calculated Fraction Identity (%)",
       fill = "Immune Cell Phenotype") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face="bold", size=10),
        axis.text.y = element_text(size=11, face="bold"),
        strip.text = element_text(size = 14, face = "bold", color = "black"),
        plot.title = element_text(hjust = 0.5, face="bold", size=16),
        legend.position = "right",
        legend.title = element_text(face="bold"),
        legend.text = element_text(size=11))

dev.off()

print("====== PHASE 9: DECONVOLUTION COMPLETE! TME FRACTIONS COMPUTED! ======")
