---
type: concept
title: Nested Cross-Validation and Feature Selection Leakage in Omics
tags: [machine-learning, cross-validation, data-leakage, random-forest, statistics]
created: 2026-08-17
updated: 2026-08-17
sources:
  - "[[tnbc-analytics-suite]]"
aliases: [Nested CV, Feature Selection Leakage, LOOCV Leakage Prevention]
---

# Nested Cross-Validation and Feature Selection Leakage in Omics

## Overview
In high-dimensional bioinformatics datasets where feature count vastly exceeds sample size ($p \gg n$, e.g., 20,000 genes vs. 20 patients), performing feature selection or differential expression filtering on the entire dataset *prior* to cross-validation causes **Data Leakage**. This results in severely inflated, over-optimistic performance metrics (e.g., AUC near 1.0). **Nested Cross-Validation** isolates feature selection strictly within training folds.

## Mechanisms & Mathematical Formulation

### The Data Leakage Trap
$$\text{Naive Pipeline}: \quad \mathcal{D}_{\text{all}} \xrightarrow{\text{Global Feature Selection}} \text{Top } k \text{ Genes} \xrightarrow{\text{Cross-Validation}} \text{Classifier Training \& Evaluation}$$
In this naive setup, the test folds have already influenced which genes were selected, violating independence assumptions.

### Methodologically Correct Nested Workflow
In each fold $i$ of Leave-One-Out Cross-Validation (LOOCV):
1. **Partition**: $\mathcal{D}_{\text{train}}^{(-i)} = \mathcal{D} \setminus \{x_i\}, \quad \mathcal{D}_{\text{test}}^{(i)} = \{x_i\}$.
2. **Internal Feature Selection**: Run two-sample t-test / DESeq2 strictly on $\mathcal{D}_{\text{train}}^{(-i)}$ to identify top $k$ features $\mathcal{S}^{(i)}$.
3. **Train**: Fit model (e.g., Random Forest) solely on $\mathcal{D}_{\text{train}}^{(-i)}[\mathcal{S}^{(i)}]$.
4. **Predict**: Evaluate the unseen sample $x_i[\mathcal{S}^{(i)}]$.
5. **Statistical Validation**: Run label Permutation Tests (e.g., 100 iterations) to verify that observed AUC exceeds random chance ($p < 0.05$).

## Analytical Relevance
- Implemented in `run_phase7_nested_cv.R` to yield an honest validation metric ($\text{AUC} \approx 0.91$) on small sample cohorts.

## Related Concepts & Entities
- [[random-forest-classifier]]
- [[triple-negative-breast-cancer]]
- [[variance-stabilizing-transformation]]

## Citations & Sources
- [[tnbc-analytics-suite]]
