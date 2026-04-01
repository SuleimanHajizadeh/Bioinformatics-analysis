# TNBC Phase 7: Machine Learning Prediction Pipeline
options(bitmapType='cairo')

if (!requireNamespace("randomForest", quietly = TRUE)) install.packages("randomForest", repos="http://cran.us.r-project.org")
if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC", repos="http://cran.us.r-project.org")

library(DESeq2)
library(randomForest)
library(pROC)
library(ggplot2)

setwd("/bioinformatics/module_4_tnbc_paper/")

print("1. Data Loading for Machine Learning...")
counts_data <- read.csv("tnbc_gene_counts.csv", row.names = 1, check.names = FALSE)
coldata <- read.csv("tnbc_phenotypes.csv", row.names = 1)
coldata <- coldata[colnames(counts_data), , drop = FALSE]
coldata$Treatment <- factor(coldata$Treatment, levels = c("Normal", "Tumor"))

# Get VST data to use as ML features
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = coldata, design = ~ Treatment)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
vsd <- vst(dds, blind=FALSE)
vsd_mat <- assay(vsd)

# Load Top DEGs for Feature Selection (We use the top 20 marker genes as our ML inputs)
res <- read.csv("tnbc_deg_results.csv", row.names = 1)
res_ordered <- res[order(res$padj), ]
top_genes <- rownames(head(res_ordered[!is.na(res_ordered$padj), ], 20))

print("2. Preparing Data Matrix for Random Forest...")
ml_data <- as.data.frame(t(vsd_mat[top_genes, ]))
colnames(ml_data) <- make.names(colnames(ml_data))
ml_data$Class <- coldata$Treatment

set.seed(42)  # For reproducibility
print("3. Training Random Forest Classifier...")
# We use the expression of 20 genes to predict if the patient has TNBC or is Normal
rf_model <- randomForest(Class ~ ., data = ml_data, ntree = 500, importance = TRUE)

print("4. Plotting Variable Importance...")
png("ML_Variable_Importance.png", width=8, height=6, res=300, units="in")
# MeanDecreaseGini tells us which genes are the most critical for the AI's decision
varImpPlot(rf_model, main="Random Forest: Top 20 Predictor Genes (TNBC vs Normal)", type=2, col="blue", pch=19, cex=1.2)
dev.off()

print("5. Generating ROC Curve & AUC...")
# We use Out-Of-Bag (OOB) predictions since our dataset is 20 samples (acts like internal cross-validation)
rf_probs <- predict(rf_model, type = "prob")[, "Tumor"]
roc_obj <- roc(ml_data$Class, rf_probs, levels=c("Normal", "Tumor"))

png("ML_ROC_Curve.png", width=6, height=6, res=300, units="in")
plot(roc_obj, main=paste0("Random Forest ROC Curve (AUC = ", round(auc(roc_obj), 3), ")"), 
     col="#d7191c", lwd=4, legacy.axes=TRUE, print.auc=FALSE)
text(0.4, 0.4, paste0("Area Under Curve (AUC) = ", round(auc(roc_obj), 3)), cex=1.3, col="black", font=2)
dev.off()

print("====== PHASE 7: MACHINE LEARNING COMPLETE! AI HAS LEARNED TO DIAGNOSE TNBC! ======")
