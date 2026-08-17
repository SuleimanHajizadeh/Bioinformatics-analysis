---
type: concept
title: Tumor Microenvironment (TME) Immune Deconvolution
tags: [immunology, tme, deconvolution, cibersort, digital-pathology]
created: 2026-08-17
updated: 2026-08-17
sources:
  - "[[tnbc-analytics-suite]]"
aliases: [TME Deconvolution, Digital Pathology, In Silico Cell Fractioning]
---

# Tumor Microenvironment (TME) Immune Deconvolution

## Overview
Bulk tumor transcriptomes reflect a composite mixture of malignant epithelial cells, stromal fibroblasts, endothelial cells, and infiltrating immune cell populations. **In Silico Immune Deconvolution** computationally infers the relative or absolute abundances of immune cell subtypes from bulk RNA-Seq expression matrices using defined marker gene panels.

## Deconvolution Methodology in This Workflow

### Marker Panels
- **CD8+ Cytotoxic T-Cells**: `CD8A`, `CD8B`, `GZMA`, `GZMB`, `PRF1` (direct cytolytic anti-tumor effectors).
- **B-Cells**: `CD19`, `MS4A1` (CD20), `CD79A`, `CD79B` (humoral response).
- **M1 Macrophages (Pro-inflammatory / Anti-Tumor)**: `NOS2`, `IRF5`, `STAT1`, `IL12A`.
- **M2 Macrophages (Immunosuppressive / Pro-Tumor)**: `CD163`, `MRC1` (CD206), `IL10`, `ARG1`, `CSF1R`.

### Mathematical Scoring
1. Scale normalized VST count matrix to Z-scores across all samples:
   $$z_{gj} = \frac{x_{gj} - \bar{x}_g}{\sigma_g}$$
2. Compute composite population scores per cell signature $k$:
   $$\text{Score}_{k, j} = \frac{1}{|S_k|} \sum_{g \in S_k} z_{gj}$$
3. Calibrate scores to non-negative relative abundance profiles across Tumor vs. Normal sample strata.

## Clinical Relevance in Triple-Negative Breast Cancer
- Elevated CD8+ T-cell infiltration is a primary predictive biomarker for response to immune checkpoint inhibitors (anti-PD-1 / anti-PD-L1 therapies).
- High M2 / M1 macrophage ratios correlate with immune evasion, tumor vascularization, and poor clinical prognosis.

## Related Concepts & Entities
- [[triple-negative-breast-cancer]]
- [[deseq2]]
- [[single-cell-clustering-umap]]

## Citations & Sources
- [[tnbc-analytics-suite]]
