---
type: entity
title: Random Forest Machine Learning Classifier
tags: [tools, machine-learning, algorithms, classification, random-forest]
created: 2026-08-17
updated: 2026-08-17
sources:
  - "[[tnbc-analytics-suite]]"
aliases: [Random Forest, Breiman Random Forest, RF Classifier]
---

# Random Forest Machine Learning Classifier

## Overview
**Random Forest** is an ensemble learning method constructed by fitting multiple decision trees on bootstrap samples of the training data and aggregating their predictions (bagging) while selecting a random subset of features at each candidate split.

## Strengths in Transcriptomics
- **Non-linear Interactions**: Naturally captures high-order gene-gene regulatory interactions.
- **Robustness to Overfitting**: Bagging and random feature subspaces mitigate noise.
- **Variable Importance**: Measures mean decrease in Gini impurity or accuracy to rank candidate biomarker genes.

## Repository Implementation
- Used in `run_phase7_nested_cv.R` embedded within a Leave-One-Out Cross-Validation (LOOCV) loop to classify TNBC vs. Normal samples based on VST-transformed gene expression, achieving $\text{AUC} \approx 0.91$.

## Related Concepts & Entities
- [[nested-cross-validation-ml]]
- [[triple-negative-breast-cancer]]
- [[variance-stabilizing-transformation]]
- [[tnbc-analytics-suite]]
