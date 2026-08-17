---
type: concept
title: Variance Stabilizing Transformation (VST)
tags: [normalization, transformation, rna-seq, homoscedasticity]
created: 2026-08-17
updated: 2026-08-17
sources:
  - "[[statistical-methods-primer]]"
aliases: [VST, Homoscedastic Transformation]
---

# Variance Stabilizing Transformation (VST)

## Overview
Raw RNA-Seq count matrices exhibit strong mean-dependent variance (heteroscedasticity). Distance-based statistical algorithms (PCA, hierarchical clustering, machine learning classifiers) assume homoscedastic errors. **Variance Stabilizing Transformation (VST)** maps count data to a log2-like continuous scale where variance is approximately independent of the mean.

## Key Formulation

Given the variance function $h(\mu) = \mu + \alpha \mu^2$, the transformation calculates:
$$y = \int \frac{1}{\sqrt{h(\mu)}} d\mu$$

### VST vs. Regularized Log (rlog)
| Property | VST | rlog |
| :--- | :--- | :--- |
| **Complexity** | $O(N)$ (linear scaling) | $O(N^2)$ (GLM fitting per gene) |
| **Outlier Robustness** | Moderate | High |
| **Recommended Cohort** | $N \ge 30$ samples | $N < 15$ samples |

## Analytical Relevance
- Essential pre-processing step for Unsupervised PCA visualization, WGCNA network generation, and Random Forest feature matrices.

## Related Concepts & Entities
- [[negative-binomial-model]]
- [[deseq2]]
- [[wgcna-coexpression-networks]]

## Citations & Sources
- [[statistical-methods-primer]]
