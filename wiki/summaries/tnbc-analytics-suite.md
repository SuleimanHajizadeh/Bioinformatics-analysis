---
type: summary
title: "TNBC Advanced Transcriptomic Analytics Suite (Phases 5-9)"
authors: "Suleyman Hajizadeh"
year: 2026
raw_source: "workflows/bulk_rnaseq/scripts/"
ingested: 2026-08-17
tags: [tnbc, deseq2, clusterprofiler, random-forest, wgcna, cibersort, deconvolution]
---

# TNBC Advanced Transcriptomic Analytics Suite (Phases 5-9)

## Core Thesis / Executive Summary
A multi-stage downstream computational suite in R for analyzing Triple-Negative Breast Cancer (TNBC) vs. Normal cohorts across five analytical dimensions:
1. **Phase 5 (DESeq2)**: Low count filtering ($\ge 10$ reads), Wald testing, PCA clustering, and Top-50 DEG heatmap.
2. **Phase 6 (Enrichment)**: Volcano plot (|log2FC| > 2.0, padj < 0.05) and clusterProfiler GO/KEGG pathway over-representation.
3. **Phase 7 (Machine Learning)**: Nested Leave-One-Out Cross-Validation (LOOCV) Random Forest classifier with feature selection within folds to prevent data leakage, plus permutation testing.
4. **Phase 8 (Co-expression Network)**: Pearson correlation adjacency matrix on top DEGs ($|R| \ge 0.85$), topological filtering, and degree centrality hub gene extraction using `igraph` & `ggraph`.
5. **Phase 9 (Microenvironment Deconvolution)**: Geometric mean Z-score digital pathology scoring for CD8+ Cytotoxic T-Cells, B-Cells, M1 Macrophages, and M2 Macrophages across tumor vs. normal phenotypes.

## Touched Wiki Pages
- Concept: [[nested-cross-validation-ml]]
- Concept: [[immune-deconvolution-tme]]
- Concept: [[wgcna-coexpression-networks]]
- Entity: [[triple-negative-breast-cancer]]
- Entity: [[deseq2]]
- Entity: [[random-forest-classifier]]
- Synthesis: [[transcriptomic-analysis-framework]]
