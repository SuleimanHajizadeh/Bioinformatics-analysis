---
type: entity
title: DESeq2
tags: [tools, software, bioconductor, r, differential-expression]
created: 2026-08-17
updated: 2026-08-17
sources:
  - "[[statistical-methods-primer]]"
aliases: [DESeq2 Package, Love et al. 2014]
---

# DESeq2

## Overview
**DESeq2** is an open-source Bioconductor R package designed for differential gene expression analysis of RNA-Seq count data. Developed by Michael Love, Wolfgang Huber, and Simon Anders, it employs shrinkage estimation for dispersions and fold changes to improve stability and interpretability.

## Core Features & Methodology
- **Normalization**: Median-of-ratios method to calculate sample-specific size factors.
- **Dispersion Estimation**: Empirical Bayes shrinkage pulling individual gene dispersions toward a trend curve over the mean.
- **Hypothesis Testing**: Generalized Linear Model (GLM) fitting and Wald / Likelihood Ratio Tests (LRT).
- **Log2 Fold Change Shrinkage**: `apeglm` or normal shrinkage estimators to compress high variance from lowly expressed genes.

## Related Concepts & Entities
- [[negative-binomial-model]]
- [[multiple-testing-fdr]]
- [[variance-stabilizing-transformation]]
- [[triple-negative-breast-cancer]]

## Citations & Sources
- [[statistical-methods-primer]]
