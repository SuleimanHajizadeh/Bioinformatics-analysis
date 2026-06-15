# ==============================================================================
# TNBC Phase 7 — NESTED LOOCV + PERMUTATION TEST (Metodoloji Cəhətdən Düzgün)
# ==============================================================================
#
# Dataset: SRP042620 (NCBI SRA)
# Nümunələr: 10 TNBC Tumor + 10 Normal (20 nümunə)
# Gen sayı: ~43,000 (DESeq2-dən keçmiş raw counts)
#
# Bu skript original run_phase7_ml.R-dəki DATA LEAKAGE xətasını aradan qaldırır.
#
# Orijinal xəta:
#   - tnbc_deg_results.csv bütün 20 xəstənin məlumatından hesablanıb (qlobal DESeq2)
#   - Sonra OOB predictions alınıb — lakin gen seçimi artıq CV-dən əvvəl edildiyi
#     üçün AUC = 0.97 süni olaraq şişib (Feature Selection Leakage)
#
# Düzgün yanaşma:
#   - DESeq2 və ya sürətli Welch t-test əsaslı gen seçimi HƏR LOOCV DÖVRÜNÜN
#     DAXİLİNDƏ (train fold-unda) aparılır
#   - Test xəstəsi (left-out sample) gen seçiminə heç vaxt toxunmur
#   - Permutasiya testi modelin statistik əhəmiyyətini sübut edir
#
# Müəllif: Antigravity AI
# Tarix: 2026-05-24
# ==============================================================================

options(bitmapType = 'cairo')
library(DESeq2)
library(randomForest)
library(pROC)
library(ggplot2)

setwd("/bioinformatics/module_4_tnbc_paper/")

cat("\n==============================================================\n")
cat("TNBC Robust ML Pipeline — Nested LOOCV + Permutation Test\n")
cat("==============================================================\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 1. VERİLƏNLƏRİN YÜKLƏNMƏSİ VƏ ÖN EMAL
# ─────────────────────────────────────────────────────────────────────────────
# VST (Variance Stabilizing Transformation): RNA-seq count məlumatlarını
# normal paylanmaya yaxınlaşdırır. Bu, t-test və Random Forest üçün əsasdır.
cat("[1/5] Verilənlər yüklənir və VST transformasiyası tətbiq olunur...\n")

counts_data <- read.csv("data/tnbc_gene_counts.csv", row.names = 1, check.names = FALSE)
coldata     <- read.csv("data/tnbc_phenotypes.csv",  row.names = 1)
coldata     <- coldata[colnames(counts_data), , drop = FALSE]
coldata$Treatment <- factor(coldata$Treatment, levels = c("Normal", "Tumor"))

dds  <- DESeqDataSetFromMatrix(countData = counts_data, colData = coldata, design = ~ Treatment)
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]

# VST — log2 benzər sabit dispersiyalı transformasiya (ML üçün daha uyğundur)
vsd     <- vst(dds, blind = TRUE) # blind=TRUE: qrup məlumatı olmadan normallaşdır
vsd_mat <- assay(vsd) # Genlər sətirdə, nümunələr sütunda: dim = [n_genes x 20]

# Maşın öyrənməsi üçün matrixi transponə edirik: [20 nümunə x n_gen]
X <- t(vsd_mat)  # X: 20 sətir (nümunə), sütunlar (genlər)
y <- coldata$Treatment  # Factor: "Tumor" / "Normal"

cat(sprintf("   → Matris ölçüsü: %d nümunə × %d gen\n\n", nrow(X), ncol(X)))

# ─────────────────────────────────────────────────────────────────────────────
# 2. ULTRA-OPTİMALLAŞDIRILMIŞ WELCH T-TEST (Vektorlaşdırılmış)
# ─────────────────────────────────────────────────────────────────────────────
# Adi apply(X, 2, t.test) ~43,000 gen üçün çox yavaşdır.
# Bu funksiya matris cəbri ilə eyni anda bütün genləri hesablayır.
fast_ttest_pvalues <- function(X_mat, y_fac) {
  g1 <- y_fac == levels(y_fac)[1]
  g2 <- y_fac == levels(y_fac)[2]
  n1 <- sum(g1); n2 <- sum(g2)
  
  m1 <- colMeans(X_mat[g1, ]); m2 <- colMeans(X_mat[g2, ])
  v1 <- colSums(sweep(X_mat[g1, ], 2, m1, "-")^2) / (n1 - 1)
  v2 <- colSums(sweep(X_mat[g2, ], 2, m2, "-")^2) / (n2 - 1)
  
  se     <- sqrt(v1/n1 + v2/n2)
  t_stat <- (m1 - m2) / se
  df     <- (v1/n1 + v2/n2)^2 / ((v1/n1)^2/(n1-1) + (v2/n2)^2/(n2-1))
  
  p_vals        <- 2 * pt(-abs(t_stat), df = df)
  names(p_vals) <- colnames(X_mat)
  return(p_vals)
}

# ─────────────────────────────────────────────────────────────────────────────
# 3. METOD A — ORIJINAL (SƏHV): Feature Selection Leakage İlə
# ─────────────────────────────────────────────────────────────────────────────
# Bu, original run_phase7_ml.R-in simulyasiyasıdır.
# DESeq2 bütün 20 nümunədən hesablanıb, sonra OOB AUC götürülüb.
cat("[2/5] METOD A (Orijinal / Leakage) hesablanır...\n")

# Bütün 20 nümunə üzərindən gen seçimi (Leakage!)
p_vals_global <- fast_ttest_pvalues(X, y)
top50_leakage <- names(sort(p_vals_global))[1:50]

preds_leakage <- numeric(nrow(X))
for (i in 1:nrow(X)) {
  X_train <- X[-i, top50_leakage, drop = FALSE]
  y_train <- y[-i]
  X_test  <- X[ i, top50_leakage, drop = FALSE]
  
  rf  <- randomForest(X_train, y_train, ntree = 500, mtry = 7, nodesize = 3)
  preds_leakage[i] <- predict(rf, X_test, type = "prob")[, "Tumor"]
}
auc_leakage <- auc(roc(y, preds_leakage, quiet = TRUE))
cat(sprintf("   → [LEAKAGE] AUC: %.4f  ← Şişirdilmiş, elmi cəhətdən etibarsız!\n\n", auc_leakage))

# ─────────────────────────────────────────────────────────────────────────────
# 4. METOD B — DÜZGÜN: Nested LOOCV
# ─────────────────────────────────────────────────────────────────────────────
# Statistik prinsip: Hər LOOCV dövrünə (iterasiyaya) baxanda, test nümunəsinin
# heç bir məlumatı gen seçiminə, normallaşdırmaya və ya model öyrənməsinə
# daxil olmamalıdır. Bu, "Optimistic Bias"-dan azaddır.
cat("[3/5] METOD B (Nested LOOCV — Düzgün) hesablanır...\n")

set.seed(42)
n <- nrow(X)
preds_nested <- numeric(n)

for (i in 1:n) {
  X_train_raw <- X[-i, ]
  y_train     <- y[-i]
  X_test_raw  <- X[ i, , drop = FALSE]
  
  # Gen seçimi YALNIZ train fold-u üzərindən (test nümunəsi görünmür)
  # Bu, "Within-fold feature selection" — Nested CV-nin əsas prinsipi
  p_vals_fold      <- fast_ttest_pvalues(X_train_raw, y_train)
  selected_genes   <- names(sort(p_vals_fold))[1:50]
  
  X_train_sel <- X_train_raw[, selected_genes, drop = FALSE]
  X_test_sel  <- X_test_raw[,  selected_genes, drop = FALSE]
  
  # mtry = sqrt(50) ≈ 7: Random Forest üçün standart feature sampling
  rf_nested <- randomForest(X_train_sel, y_train, ntree = 500, mtry = 7, nodesize = 3)
  preds_nested[i] <- predict(rf_nested, X_test_sel, type = "prob")[, "Tumor"]
}

roc_nested <- roc(y, preds_nested, quiet = TRUE)
auc_nested <- auc(roc_nested)
cat(sprintf("   → [NESTED CV] Dürüst AUC: %.4f\n\n", auc_nested))

# ─────────────────────────────────────────────────────────────────────────────
# 5. PERMUTASİYA TESTİ
# ─────────────────────────────────────────────────────────────────────────────
# Null hypothesis (H₀): Modelin müşahidə etdiyi AUC bioloji siqnaldan deyil,
# təsadüfi nümunələnmə variasiyasından irəli gəlir.
# Əgər p < 0.05 → H₀ rədd edilir → Model statistik olaraq əhəmiyyətlidir.
cat("[4/5] Permutasiya Testi (n=50) icra olunur...\n")

set.seed(123)
n_perm <- 50
null_aucs <- numeric(n_perm)

for (p in 1:n_perm) {
  y_perm  <- sample(y)
  preds_p <- numeric(n)
  
  for (i in 1:n) {
    X_train_raw  <- X[-i, ]
    y_train_perm <- y_perm[-i]
    X_test_raw   <- X[ i, , drop = FALSE]
    
    p_vals_perm    <- fast_ttest_pvalues(X_train_raw, y_train_perm)
    sel_perm       <- names(sort(p_vals_perm))[1:50]
    X_tr_sel       <- X_train_raw[, sel_perm, drop = FALSE]
    X_te_sel       <- X_test_raw[,  sel_perm, drop = FALSE]
    
    rf_p        <- randomForest(X_tr_sel, y_train_perm, ntree = 500, mtry = 7, nodesize = 3)
    preds_p[i]  <- predict(rf_p, X_te_sel, type = "prob")[, "Tumor"]
  }
  null_aucs[p] <- auc(roc(y_perm, preds_p, quiet = TRUE))
  cat(sprintf("   Permutasiya %2d/50: AUC = %.4f\n", p, null_aucs[p]))
}

p_value <- sum(null_aucs >= auc_nested) / n_perm
cat(sprintf("\n   → Permutasiya Orta AUC: %.4f\n", mean(null_aucs)))
cat(sprintf("   → Empirik P-dəyəri: %.4f\n\n", p_value))

# ─────────────────────────────────────────────────────────────────────────────
# 6. ROC ƏYRİSİ QRAFIKI (NESTED CV)
# ─────────────────────────────────────────────────────────────────────────────
cat("[5/5] ROC əyrisi qrafiki yaradılır...\n")

roc_data <- data.frame(
  FPR        = 1 - roc_nested$specificities,
  TPR        = roc_nested$sensitivities
)

roc_plot <- ggplot(roc_data, aes(x = FPR, y = TPR)) +
  geom_line(color = "#2e86de", linewidth = 2) +
  geom_abline(linetype = "dashed", color = "gray50") +
  annotate("text", x = 0.5, y = 0.2,
           label = sprintf("Nested CV AUC = %.3f\nPermutasiya P = %.3f", auc_nested, p_value),
           size = 5, fontface = "bold", color = "#d63031") +
  labs(
    title    = "ROC Curve — Nested LOOCV (Dürüst Qiymətləndirmə)",
    subtitle = "TNBC vs Normal — SRP042620 (10+10 nümunə, top 50 gen)",
    x        = "False Positive Rate (1 - Specificity)",
    y        = "True Positive Rate (Sensitivity)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40")
  )

ggsave("results/ML_ROC_Curve_NestedCV.png", plot = roc_plot, width = 7, height = 6, dpi = 300)

# ─────────────────────────────────────────────────────────────────────────────
# 7. YEKUN HESABAT
# ─────────────────────────────────────────────────────────────────────────────
cat("\n====================================================\n")
cat("            YEKUN HESABAT (Final Report)\n")
cat("====================================================\n")
cat(sprintf("  Nümunə sayı:       20 (10 Tumor + 10 Normal)\n"))
cat(sprintf("  Toplam gen sayı:   %d (VST-dən sonra filtr: count >= 10)\n", ncol(X)))
cat(sprintf("  Seçilən gen sayı:  50 (LOOCV daxilindəki Welch t-test)\n\n"))
cat(sprintf("  [LEAKAGE]   Orijinal (Yanlış) AUC:  %.4f  ← Elmi cəhətdən etibarsız\n", auc_leakage))
cat(sprintf("  [NESTED CV] Dürüst AUC:             %.4f  ← Elmi cəhətdən etibarlı\n\n", auc_nested))
cat(sprintf("  Permutasiya P-dəyəri:  %.4f\n", p_value))

if (p_value < 0.05) {
  cat("  ✅ Nəticə: Model statistik olaraq əhəmiyyətlidir (p < 0.05).\n")
  cat("     TNBC-nin gen ifadəsi profili gerçəkdən ayırd edilə bilir.\n")
} else {
  cat("  ⚠️  Nəticə: Model şans səviyyəsindədir (p >= 0.05).\n")
  cat("     Nümunə sayını artırın və ya gen seçim metodunu yenidən nəzərdən keçirin.\n")
}
cat("====================================================\n")
