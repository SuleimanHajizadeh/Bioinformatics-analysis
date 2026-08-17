---
type: summary
title: "HCC1395 Breast Cancer Reference Cell Line RNA-Seq Benchmark"
authors: "Bioinformatics Practice Team"
year: 2026
raw_source: "exploratory_analysis/practice_modules/module_2/scripts/deseq2_analysis.R"
ingested: 2026-08-17
tags: [benchmark, hcc1395, seqc2, deseq2, hisat2, tnbc, cell-line]
---

# HCC1395 Breast Cancer Reference Cell Line RNA-Seq Benchmark

## Core Thesis / Executive Summary
Implementation of standard quality control, HISAT2 alignment, and DESeq2 differential analysis on the gold-standard **HCC1395** (primary TNBC cell line) and **HCC1395BL** (matched normal Epstein-Barr virus-transformed B-lymphoblastoid line) reference benchmark derived from the FDA SEQC2 Consortium.

## Methodological Workflow
- **Data Ingestion**: Raw featureCounts table (`final_counts.csv`) for replicates SRR1553606, SRR1553607, SRR1553608.
- **Filtering**: Removal of low-count genes ($\sum \text{counts} < 10$).
- **Variance Stabilization**: `varianceStabilizingTransformation` (`vsd`, `blind=FALSE`).
- **Quality Control**:
  - Sample Distance Matrix Heatmap (`dist(t(assay(vsd)))`) confirming distinct tumor vs. normal clustering.
  - Principal Component Analysis (PCA) demonstrating clear separation along PC1.
- **Statistical Model**: Wald test with Benjamini-Hochberg adjustment, extracting DEGs with $|\log_2\text{FC}| > 2$ and $q < 0.05$.

## Touched Wiki Pages
- Entity: [[hcc1395-cell-line]]
- Entity: [[deseq2]]
- Concept: [[variance-stabilizing-transformation]]
- Concept: [[negative-binomial-model]]
