---
type: entity
title: Seurat R Package for Single-Cell Genomics
tags: [tools, software, r, single-cell, scrna-seq, bioconductor]
created: 2026-08-17
updated: 2026-08-17
sources:
  - "[[single-cell-seurat-workflow]]"
aliases: [Seurat, Satija Lab Seurat]
---

# Seurat R Package for Single-Cell Genomics

## Overview
**Seurat** is an R package designed for quality control, analysis, and exploration of single-cell RNA-seq data developed by Rahul Satija's laboratory at the New York Genome Center. It provides comprehensive tools for multimodal data integration, clustering, and differential expression analysis.

## Key Capabilities & Data Structures
- **`SeuratObject`**: Multi-layered container storing raw counts, normalized assays, dimensionality reductions (PCA, UMAP, t-SNE), and cell-level metadata.
- **Normalization**: Standard `LogNormalize` and advanced regularized negative binomial models (`SCTransform`).
- **Feature Selection**: Finding highly variable features across cells using variance-stabilizing transformation.
- **Graph Clustering**: Louvain and Leiden community detection algorithms.

## Repository Implementation
- Drives the single-cell transcriptomics pipeline (`workflows/single_cell/scripts/`), taking 10x Genomics PBMC datasets through QC, UMAP clustering, and cell-type annotation.

## Related Concepts & Entities
- [[single-cell-clustering-umap]]
- [[variance-stabilizing-transformation]]
- [[single-cell-seurat-workflow]]
