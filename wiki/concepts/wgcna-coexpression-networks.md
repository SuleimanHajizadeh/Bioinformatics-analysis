---
type: concept
title: Weighted Gene Co-expression Network Analysis (WGCNA)
tags: [wgcna, systems-biology, coexpression, hub-genes, tnbc]
created: 2026-08-17
updated: 2026-08-17
sources:
  - "[[statistical-methods-primer]]"
aliases: [WGCNA, Co-expression Networks, Hub Gene Analysis]
---

# Weighted Gene Co-expression Network Analysis (WGCNA)

## Overview
**Weighted Gene Co-expression Network Analysis (WGCNA)** is a systems biology method for finding clusters (modules) of highly correlated genes, summarizing such clusters using module eigengenes or intramodular hub genes, and relating modules to clinical traits.

## Key Principles & Algorithm Steps

1. **Similarity Matrix**: Pairwise Pearson/Spearman correlation between gene expression profiles:
   $$s_{ij} = |\text{cor}(x_i, x_j)|$$
2. **Adjacency Calculation (Scale-Free Criterion)**: Raising similarity to power $\beta$ to penalize weak correlations:
   $$a_{ij} = s_{ij}^\beta$$
   Power $\beta$ is picked using the scale-free topology criterion ($R^2 \ge 0.85$).
3. **Topological Overlap Matrix (TOM)**: Measures network interconnectedness beyond pairwise relationships:
   $$\text{TOM}_{ij} = \frac{l_{ij} + a_{ij}}{\min(k_i, k_j) + 1 - a_{ij}}$$
4. **Hierarchical Clustering & Dynamic Tree Cut**: Partitions TOM distance into discrete co-expression color modules.
5. **Intramodular Connectivity ($k_{\text{within}}$)**: Identifies central hub genes that drive module function.

## Application in Triple-Negative Breast Cancer
- In this repository's TNBC workflow ($\beta = 12, R^2 = 0.87$), the **turquoise module** (312 genes) showed marked enrichment in PI3K-Akt oncogenic signaling.
- **Top Hub Gene**: *AKT1* ($k_{\text{within}} = 47.3$).

## Related Concepts & Entities
- [[triple-negative-breast-cancer]]
- [[deseq2]]
- [[variance-stabilizing-transformation]]

## Citations & Sources
- [[statistical-methods-primer]]
