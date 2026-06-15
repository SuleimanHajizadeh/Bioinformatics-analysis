# Statistical Methodology Primer: RNA-Seq Differential Expression Analysis

This document provides the mathematical and statistical foundations for the transcriptomics pipeline implemented in this repository. In top-tier research institutions, understanding the underlying mathematical frameworks of software packages (e.g., `DESeq2`, `edgeR`) is critical for ensuring experimental reproducibility and avoiding false discovery artifacts.

---

## 1. The Negative Binomial Model for RNA-Seq Count Data

### Why Normal Distribution Fails
For high-throughput sequencing, the raw data consists of discrete read counts mapped to genomic features (genes/transcripts). Traditional microarray analyses assumed a continuous Normal (Gaussian) distribution:
$$Y_i \sim \mathcal{N}(\mu, \sigma^2)$$

However, read counts are:
1. **Discrete:** Counts are non-negative integers ($k \in \{0, 1, 2, \dots\}$).
2. **Skewed:** The distributions of lowly expressed genes are highly skewed toward zero.
3. **Heteroscedastic:** The variance is not constant across different expression levels.

### The Poisson Assumption and Overdispersion
If reads are sampled independently at random from a large pool, the probability of mapping a read to a specific gene can be modeled by a **Poisson distribution**:
$$P(X = k) = \frac{\lambda^k e^{-\lambda}}{k!}$$

A fundamental property of the Poisson distribution is that the mean is equal to the variance:
$$\text{E}[X] = \text{Var}(X) = \lambda$$

In real biological replicates, this assumption is violated. Due to biological variation (e.g., cell-to-cell differences, environmental shifts, genetic variations), the variance of read counts exceeds the mean, a phenomenon known as **overdispersion**:
$$\text{Var}(X) > \text{E}[X]$$

Using a Poisson model under overdispersion underestimates the variance, leading to artificially low p-values and an inflated False Positive rate.

### The Negative Binomial Solution
To model overdispersion, we represent biological variance by letting the Poisson parameter $\lambda$ itself vary according to a Gamma distribution ($\lambda \sim \text{Gamma}(\alpha, \beta)$). Compounding a Poisson distribution with a Gamma distribution yields the **Negative Binomial (NB) distribution**:

$$X \sim \text{NB}(\mu, \alpha)$$

Under the Negative Binomial distribution, the mean-variance relationship is defined as:
$$\text{Var}(X) = \mu + \alpha \mu^2$$

where:
* $\mu$ is the mean read count.
* $\alpha$ is the **dispersion parameter**, reflecting the biological coefficient of variation ($BCV = \sqrt{\alpha}$).

As the dispersion parameter $\alpha \to 0$, the Negative Binomial distribution collapses into the Poisson distribution ($\text{Var}(X) \to \mu$).

In `DESeq2`, the dispersion $\alpha_g$ for each gene $g$ is estimated using a maximum likelihood framework and shrunk toward a curve representing the global mean-dispersion trend across all genes, utilizing an empirical Bayes approach.

---

## 2. Multiple Testing Correction: Benjamini-Hochberg FDR

### The Multiple Testing Problem
In a typical bulk RNA-seq experiment, statistical tests (e.g., Wald test) are performed simultaneously for each of $m \approx 20,000$ genes. If we select a standard significance threshold of $\alpha = 0.05$:

$$\text{Expected False Positives} = m \times \alpha = 20,000 \times 0.05 = 1,000$$

We would identify roughly 1,000 genes as "differentially expressed" purely by random chance, making downstream biological interpretation impossible.

### Family-Wise Error Rate (FWER) vs. False Discovery Rate (FDR)
* **Bonferroni Correction (FWER):** Limits the probability of making *at least one* false positive across the entire experiment. The adjusted significance threshold is:
  $$\alpha_{\text{Bonf}} = \frac{\alpha}{m} = \frac{0.05}{20,000} = 2.5 \times 10^{-6}$$
  This is extremely conservative and leads to high False Negative rates (missing true biological signals).
* **Benjamini-Hochberg (FDR):** Controls the expected *proportion* of false discoveries among the rejected null hypotheses:
  $$\text{FDR} = \text{E}\left[\frac{V}{R}\right] \quad (\text{where } R \text{ is total discoveries and } V \text{ is false positives})$$

### The Benjamini-Hochberg Algorithm
For a target FDR level $q$ (typically 0.05):
1. Sort the raw p-values of all $m$ tests in ascending order:
   $$P_{(1)} \le P_{(2)} \le \dots \le P_{(m)}$$
2. Find the largest index $k$ such that:
   $$P_{(k)} \le \frac{k}{m} q$$
3. Reject the null hypothesis $H_0$ for all genes with $P_{(i)}$ where $i \le k$.

The adjusted p-value (or q-value) for gene $i$ is calculated as:
$$q_i = \min \left( \min_{j \ge i} \left\{ \frac{m \cdot P_{(j)}}{j} \right\}, 1 \right)$$

This approach controls the FDR at $q$ while maximizing statistical power.

---

## 3. Data Transformation: Variance Stabilizing Transformation (VST) vs. Regularized Log (rlog)

Raw read counts are heteroscedastic: genes with high mean counts have high variance. However, standard machine learning methods (e.g., PCA, clustering, Random Forests) require homoscedastic data (variance independent of the mean).

To solve this, `DESeq2` provides two transformations:
1. **Variance Stabilizing Transformation (VST):** Computes a variance-stabilizing transformation from the fitted mean-variance trend:
   $$y = \int \frac{1}{\sqrt{h(\mu)}} d\mu$$
   where $h(\mu) = \mu + \alpha \mu^2$ is the variance function.
2. **Regularized Logarithmic Transformation (rlog):** Fits a generalized linear model (GLM) for each gene, using a ridge-like penalty to shrink log-fold changes for genes with low counts or high dispersion towards zero.

### Comparison and Choice
| Feature | VST | rlog |
|---------|-----|------|
| **Computation Speed** | Extremely fast ($O(N)$) | Slow for large sample sizes ($O(N^2)$) |
| **Sensitivity to Outliers** | Moderately robust | Highly robust (strong shrinkage) |
| **Use Case** | Medium to large sample sizes ($N > 30$) | Small datasets ($N < 15$) |

**Decision:** In our pipeline, **VST** is selected for data normalization prior to PCA and clustering due to its computational efficiency and optimal stabilization of high-throughput sequencing variance.

---

## 4. Addressing Batch Effects: ComBat-seq and GLM Design

Batch effects are systematic technical variations introduced when samples are processed in different groups, days, or sequencing lanes.

### Designing the GLM formula
When biological replicates are nested in different batches, batch effects can be accounted for directly inside the `DESeq2` generalized linear model (GLM) design formula:
$$\log(\mu_{ij}) = \beta_0 + \beta_{\text{batch}} X_{i,\text{batch}} + \beta_{\text{condition}} X_{i,\text{condition}}$$

Including the batch factor in the design formula controls for the batch-associated variance, enabling unbiased estimation of $\beta_{\text{condition}}$.

### Empirical Bayes Batch Correction (ComBat-seq)
For algorithms that require batch-corrected expression matrices directly (e.g., WGCNA, clustering), adjusting the model design is insufficient. In such cases, we apply **ComBat-seq**, an empirical Bayes approach specifically designed for count matrices. 

ComBat-seq models the counts using a Negative Binomial regression and estimates batch parameters, returning corrected counts:
$$X_{ij}^{\text{adjusted}} = g(X_{ij}, \hat{\gamma}_i, \hat{\delta}_i)$$
where $\gamma_i$ and $\delta_i$ represent the batch-specific location and scale parameters. This approach preserves the integer structure of the count data, allowing downstream programs to process the adjusted matrix without violating count distributional assumptions.
