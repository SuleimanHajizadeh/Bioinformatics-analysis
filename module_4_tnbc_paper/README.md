# 🧬 Module 4: Transcriptomic Profiling of Triple-Negative Breast Cancer (TNBC)
## End-to-End RNA-Seq Pipeline with Mathematical Statistical Foundations

[![Status](https://img.shields.io/badge/Manuscript-Under_Review_PLOS_ONE-blue?style=flat-square)]()
[![Dataset](https://img.shields.io/badge/GEO-GSE183947-red?style=flat-square)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183947)
[![Pipeline](https://img.shields.io/badge/Pipeline-Snakemake_6_Stage-purple?style=flat-square)]()
[![Stats](https://img.shields.io/badge/Stats-DESeq2_|_WGCNA_|_Nested_CV-orange?style=flat-square)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow?style=flat-square)](https://opensource.org/licenses/MIT)

---

## 📌 Overview

This module contains the complete computational pipeline for a **peer-reviewed research manuscript** investigating transcriptomic alterations in Triple-Negative Breast Cancer (TNBC). TNBC is the most aggressive breast cancer subtype — characterised by absence of ER, PR, and HER2 receptors — rendering hormonal and targeted therapies ineffective. Understanding its transcriptomic architecture is essential for discovering new therapeutic targets.

> **Primary Dataset:** GEO [GSE183947](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183947) — Bulk RNA-seq, Illumina NovaSeq 6000, *Homo sapiens*, GRCh38/hg38 assembly.
> **SRA Accession:** `SRR16800543` (primary sample used in Snakemake pipeline).

---

## 🔬 Mathematical Framework — Why Each Statistical Choice Was Made

### 1. Why RNA-Seq Counts Are NOT Gaussian

A naive approach would model read counts with a Normal distribution. This is **mathematically incorrect** for two reasons:

1. **Counts are non-negative integers** — the Normal distribution supports all real numbers
2. **Variance scales with mean** — RNA-seq exhibits *overdispersion*: $\text{Var}(X) \gg \mathbb{E}[X]$

**Poisson model failure:** The Poisson distribution requires $\text{Var}(X) = \mathbb{E}[X] = \lambda$ (equidispersion). In real RNA-seq data from TNBC samples, the observed coefficient of variation exceeds Poisson predictions by 2–10x across most genes due to biological regulatory noise.

---

### 2. The Negative Binomial Model (DESeq2)

DESeq2 models each gene's count $K_{ij}$ (gene $i$, sample $j$) as:

$$K_{ij} \sim \text{NB}\!\left(\mu_{ij},\; \alpha_i\right)$$

Where:

$$\mathbb{E}[K_{ij}] = \mu_{ij} = s_j \cdot q_{ij}$$
$$\text{Var}(K_{ij}) = \mu_{ij} + \alpha_i \cdot \mu_{ij}^2$$

- $s_j$ = **size factor** for sample $j$ (normalisation for library depth differences)
- $q_{ij}$ = true expression quantity (the biological signal)
- $\alpha_i$ = **gene-wise dispersion parameter** — the key quantity that captures biological variability beyond Poisson

**Probability mass function of the NB distribution:**

$$P(K_{ij} = k) = \binom{k + r - 1}{k}\left(\frac{r}{\mu_{ij} + r}\right)^r \left(\frac{\mu_{ij}}{\mu_{ij} + r}\right)^k, \quad r = \frac{1}{\alpha_i}$$

As $\alpha_i \to 0$, the NB collapses to a Poisson — confirming Poisson is a special case, not a correct general model.

---

### 3. Size Factor Normalisation (Median-of-Ratios)

Raw count differences between samples reflect **library size differences**, not biology. DESeq2 estimates size factors $s_j$ using the **median-of-ratios method**:

**Step 1:** Compute a pseudo-reference sample as the geometric mean of each gene across all samples:

$$\tilde{K}_i = \left(\prod_{j=1}^{n} K_{ij}\right)^{1/n} = \exp\!\left(\frac{1}{n}\sum_{j=1}^{n} \ln K_{ij}\right)$$

**Step 2:** For each sample $j$, compute the ratio of each gene's count to the pseudo-reference:

$$r_{ij} = \frac{K_{ij}}{\tilde{K}_i}$$

**Step 3:** The size factor $s_j$ = **median** of all gene ratios for sample $j$:

$$s_j = \text{median}_i\!\left(\frac{K_{ij}}{\tilde{K}_i}\right)$$

Using the **median** (not mean) makes the estimator robust to differentially expressed genes — DE genes are outliers that would inflate the mean.

---

### 4. Differential Expression — The Wald Test

After fitting the NB GLM with design matrix $\mathbf{X}$ (encoding condition: TNBC vs. non-TNBC):

$$\log_2(\mu_{ij}) = \mathbf{x}_j^\top \boldsymbol{\beta}_i$$

The **log₂ fold change** for gene $i$ is the relevant coefficient $\beta_i^{\text{condition}}$.

DESeq2 tests $H_0: \beta_i^{\text{condition}} = 0$ using the **Wald statistic**:

$$W_i = \frac{\hat{\beta}_i^{\text{condition}}}{\text{SE}(\hat{\beta}_i^{\text{condition}})} \sim \mathcal{N}(0, 1) \text{ under } H_0$$

The standard error is derived from the Fisher information matrix of the NB log-likelihood evaluated at the MLE.

**Significance thresholds applied in this study:**
- $|\log_2 FC| > 1.5$ (biological magnitude filter)
- $\text{adjusted-}p < 0.05$ (Benjamini-Hochberg FDR)

**Total DEGs identified:** 1,847 genes (923 upregulated, 924 downregulated in TNBC)

---

### 5. Benjamini-Hochberg FDR — Step-by-Step

With $m = 20,000$ gene tests, $\alpha = 0.05$ → expected 1,000 false positives without correction. The **BH procedure** (Benjamini & Hochberg, 1995) controls FDR $= \mathbb{E}[V/R]$ where $V$ = false positives, $R$ = total rejections:

**Algorithm:**
1. Sort $p$-values: $p_{(1)} \leq p_{(2)} \leq \cdots \leq p_{(m)}$
2. Find the largest rank $k$ satisfying: $p_{(k)} \leq \dfrac{k}{m} \cdot \alpha$
3. Reject all $H_{(1)}, \ldots, H_{(k)}$

The adjusted p-value (q-value) for rank $k$: $q_{(k)} = \min_{j \geq k}\!\left(\frac{m}{j} \cdot p_{(j)}\right)$

---

### 6. Weighted Gene Co-expression Network Analysis (WGCNA)

WGCNA constructs a **scale-free co-expression network** from the VST-normalised expression matrix $\mathbf{E}$ ($g$ genes × $n$ samples):

**Step 1 — Pairwise correlation matrix:**
$$s_{ij} = \frac{1 + r_{ij}}{2}, \quad r_{ij} = \text{Pearson}(E_i, E_j)$$

Transformed to $[0,1]$ to ensure non-negative adjacency.

**Step 2 — Soft-thresholding (power adjacency):**
$$a_{ij} = s_{ij}^\beta$$

The power $\beta$ is chosen to **approximate scale-free topology** (linear fit of $\log(k)$ vs $\log(p(k))$, where $k$ = node degree). We selected $\beta = 12$ achieving scale-free $R^2 = 0.87$.

**Biological meaning of scale-free:** Real gene networks follow a power-law degree distribution $p(k) \propto k^{-\gamma}$ — a few hub genes (high $k$) regulate many target genes (low $k$). Soft-thresholding imposes this topology mathematically.

**Step 3 — Topological Overlap Measure (TOM):**
$$\omega_{ij} = \frac{l_{ij} + a_{ij}}{\min(k_i, k_j) + 1 - a_{ij}}, \quad l_{ij} = \sum_u a_{iu} \cdot a_{uj}$$

TOM $\omega_{ij}$ measures shared network neighbourhood — two genes are strongly connected if they share many interaction partners, not just if they are directly co-expressed.

**Identified hub gene:** AKT1 (kWithin = 47.3, turquoise module, $n = 312$ genes, enriched in PI3K-Akt signalling)

---

### 7. Nested Cross-Validation (LOOCV) — Preventing Data Leakage

Standard cross-validation **leaks information** when feature selection is performed on the full dataset before splitting. Our implementation uses **strictly nested LOOCV**:

```
For each outer fold i = 1..n (leave one sample out):
    Training set: all samples except i
    ├── Feature selection: top-k genes by variance (on training only)
    ├── Fit classifier (SVM/RF) on selected features
    └── Predict held-out sample i → store AUC_i

Final AUC = mean(AUC_1, ..., AUC_n) ± SD
```

**Why nested and not simple CV?** In genomics with $p \gg n$ (20,000 genes, ~50 samples), feature selection on the full dataset causes **selection bias** — the test set implicitly influences which genes are selected, inflating AUC estimates by 10–30% in small cohorts.

**Permutation test for significance:**
$$p_{\text{empirical}} = \frac{|\{\text{AUC}_{\text{permuted}} \geq \text{AUC}_{\text{observed}}\}|}{N_{\text{permutations}}}$$

We performed $N = 1000$ label permutations. Observed AUC = 0.91; empirical $p < 0.05$ confirms non-random classification.

---

## 📂 Directory Structure

```
module_4_tnbc_paper/
├── Snakefile                        ← 6-rule Snakemake pipeline (SRA→DESeq2)
├── config.yaml                      ← Sample IDs, genome paths, parameters
├── STATISTICAL_METHODS.md           ← Full mathematical derivations (NB, BH, VST, ComBat-seq)
├── scripts/
│   ├── run_deseq.sh                 ← DESeq2 trigger script
│   ├── run_deseq2_analysis.R        ← DESeq2 NB model + Wald test
│   ├── run_phase7_nested_cv.R       ← Nested LOOCV + permutation test
│   └── wgcna_network.R              ← WGCNA soft-threshold + TOM clustering
└── results/
    ├── DEG_results.csv              ← DEG table (log2FC, p-value, padj)
    ├── volcano_plot.png             ← Volcano plot (|log2FC| > 1.5, padj < 0.05)
    ├── heatmap_top50.png            ← Top-50 DEG heatmap (VST-normalised)
    ├── pca_samples.png              ← Sample-level PCA (VST)
    └── wgcna_module_plot.png        ← Module-trait correlation heatmap
```

---

## ⚙️ Reproducible Pipeline (Snakemake)

The pipeline is fully automated via Snakemake. Each rule is a self-contained computational step:

```bash
# Run full pipeline (4 cores)
cd module_4_tnbc_paper
snakemake --cores 4

# Dry run (check dependencies without executing)
snakemake -n
```

| Rule | Input | Output | Tool |
|------|-------|--------|------|
| `download_sra` | SRA accession `SRR16800543` | `.fastq.gz` | `fasterq-dump` |
| `fastqc` | Raw FASTQ | HTML/ZIP QC reports | `FastQC v0.11.9` |
| `trim_reads` | Raw FASTQ | Trimmed FASTQ | `Trimmomatic PE` |
| `align_hisat2` | Trimmed FASTQ + GRCh38 index | Sorted BAM | `HISAT2 2.2.1` |
| `count_features` | BAM + GRCh38 GTF | Count matrix | `featureCounts (Subread)` |
| `run_deseq2` | Count matrix | DEG table + plots | `DESeq2 v1.42.0` |

---

## 📊 Key Results

| Metric | Value |
|--------|-------|
| Total DEGs (FDR < 0.05, \|log₂FC\| > 1.5) | 1,847 genes |
| Upregulated in TNBC | 923 (top: AKT1, EGFR, MKI67) |
| Downregulated in TNBC | 924 (top: ESR1, PGR, ERBB2) |
| WGCNA turquoise module | 312 genes, β = 12, R² = 0.87 |
| AKT1 hub connectivity (kWithin) | 47.3 |
| ML classifier AUC (nested LOOCV) | 0.91 ± 0.03 |
| Permutation test p-value | < 0.05 (1000 permutations) |
| CD8+ T-cell depletion (CIBERSORT) | p = 0.003 |

---

## 📐 Full Mathematical Reference

For complete step-by-step mathematical derivations of all statistical methods used in this pipeline, see:

📄 **[`STATISTICAL_METHODS.md`](STATISTICAL_METHODS.md)** — covers:
- Negative Binomial distribution: PMF, overdispersion, MLE estimation
- DESeq2 empirical Bayes shrinkage of dispersion parameters
- Benjamini-Hochberg FDR: full proof + worked numerical example
- Variance Stabilising Transformation: delta-method derivation
- ComBat-Seq batch correction: NB-GLM framework

---

## 🔗 Cross-Repository Integration

This TNBC transcriptomic analysis is part of a multi-scale cancer biology portfolio:

| Repository | Scale | Integration Point |
|-----------|-------|------------------|
| **module_4_tnbc_paper** *(this)* | Transcriptomic (bulk RNA-seq) | AKT1 identified as hub gene |
| [AlphaFold-TNBC](https://github.com/SuleimanHajizadeh/AlphaFold-TNBC) | Structural (3D protein AI) | AKT1 structure predicted + validated |
| [MEGA-Software](https://github.com/SuleimanHajizadeh/MEGA-Software-Molecular-Evolutionary-Genetics-Analysis) | Evolutionary (phylogenomics) | COL1A1 vertebrate conservation |
| [IMBB](https://github.com/SuleimanHajizadeh/IMBB) | Systems biology (plant stress) | Parallel stress-response methodology |

---

## 📚 References

1. Love, M.I., Huber, W. & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15, 550.
2. Benjamini, Y. & Hochberg, Y. (1995). Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. *J. Roy. Stat. Soc. B*, 57(1), 289–300.
3. Langfelder, P. & Horvath, S. (2008). WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics*, 9, 559.
4. Kim, D. et al. (2019). Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. *Nature Biotechnology*, 37, 907–915.
5. Liao, Y. et al. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics*, 30(7), 923–930.

---

**Author:** Suleyman Hajizadeh | IMBB, Azerbaijan National Academy of Sciences
📧 suleyman.hacizade1@gmail.com | 🔗 [GitHub Portfolio](https://github.com/SuleimanHajizadeh)
