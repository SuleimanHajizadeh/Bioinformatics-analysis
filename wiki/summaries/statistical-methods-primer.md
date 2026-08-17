---
type: summary
title: "Statistical Methodology Primer: RNA-Seq Differential Expression Analysis"
authors: "Bioinformatics Pipeline Team"
year: 2026
raw_source: "docs/STATISTICAL_METHODS.md"
ingested: 2026-08-17
tags: [statistics, rna-seq, deseq2, normalization, multiple-testing]
---

# Statistical Methodology Primer: RNA-Seq Differential Expression Analysis

## Core Thesis / Executive Summary
High-throughput RNA sequencing produces discrete, overdispersed, and heteroscedastic read count data. Standard Gaussian/Normal models and simple Poisson assumptions fail to model biological variability accurately. The robust mathematical framework uses the **Negative Binomial (NB)** distribution coupled with empirical Bayes dispersion shrinkage, **Benjamini-Hochberg FDR** correction for multiple testing, and **Variance Stabilizing Transformation (VST)** for downstream machine learning and clustering.

## Key Takeaways
- **Overdispersion**: Biological variation results in $\text{Var}(X) > \text{E}[X]$; Poisson underestimation leads to high false positive rates.
- **Negative Binomial Model**: $X \sim \text{NB}(\mu, \alpha)$ with $\text{Var}(X) = \mu + \alpha \mu^2$, where $\alpha$ is the dispersion parameter.
- **Dispersion Shrinkage**: DESeq2 uses an empirical Bayes approach to shrink gene-wise dispersion estimates toward a global trend curve.
- **FDR Control**: Benjamini-Hochberg step-up procedure controls the False Discovery Rate at $q=0.05$ while maintaining higher statistical power than Bonferroni FWER.
- **Data Transformation**: VST ($O(N)$) stabilizes variance across the mean without suffering from the $O(N^2)$ computational overhead of `rlog` on larger cohorts.

## Touched Wiki Pages
- Concept: [[negative-binomial-model]]
- Concept: [[multiple-testing-fdr]]
- Concept: [[variance-stabilizing-transformation]]
- Entity: [[deseq2]]
- Synthesis: [[transcriptomic-analysis-framework]]
