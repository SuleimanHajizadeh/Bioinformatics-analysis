---
type: summary
title: "Single-Cell Transcriptomics Analysis Pipeline with Seurat"
authors: "Suleyman Hajizadeh"
year: 2026
raw_source: "workflows/single_cell/scripts/"
ingested: 2026-08-17
tags: [single-cell, scrna-seq, seurat, 10x-genomics, umap, clustering]
---

# Single-Cell Transcriptomics Analysis Pipeline with Seurat

## Core Thesis / Executive Summary
An end-to-end R workflow utilizing **Seurat** for 10x Genomics scRNA-Seq data (PBMC 3k benchmark). The pipeline encompasses cell-level quality control (filtering low feature counts and high mitochondrial content), log-normalization, highly variable gene selection, dimensionality reduction (PCA and UMAP), Louvain graph-based clustering, and differential biomarker extraction for cell-type annotation.

## Pipeline Architecture & Parameters
1. **Quality Control (`seurat_qc.R`)**:
   - Initial filtering: `min.cells = 3`, `min.features = 200`.
   - Mitochondrial threshold: `PercentageFeatureSet(pbmc, pattern = "^MT-")`.
   - Stringent cutoff: $200 < \text{nFeature\_RNA} < 2500$ and $\text{percent.mt} < 5\%$.
2. **Clustering & Embedding (`seurat_cluster.R`)**:
   - Normalization: `LogNormalize` (scale factor 10,000).
   - Feature Selection: Top 2,000 variable genes via `vst`.
   - Scaling & PCA: Centered/scaled on all genes, top 10 principal components retained.
   - Graph Clustering: Shared Nearest Neighbor (SNN) graph with `FindNeighbors(dims=1:10)` and Louvain modularity optimization (`FindClusters(resolution=0.5)`).
   - Non-linear Embedding: UMAP projection on top 10 PCs.
3. **Marker Discovery (`seurat_markers.R`) & Annotation (`seurat_annotate.R`)**:
   - Wilcoxon Rank-Sum test for cluster-specific marker genes.
   - Assignment of canonical PBMC identities (CD4+ T-Cells, CD8+ T-Cells, B-Cells, CD14+ Monocytes, NK Cells, FCGR3A+ Monocytes, Dendritic Cells, Platelets).

## Touched Wiki Pages
- Concept: [[single-cell-clustering-umap]]
- Entity: [[seurat]]
- Synthesis: [[multi-omics-analytical-ecosystem]]
