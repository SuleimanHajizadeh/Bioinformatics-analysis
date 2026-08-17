---
type: concept
title: Pathway Enrichment Analysis (ORA vs. GSEA)
tags: [pathway, enrichment, gsea, ora, clusterprofiler, gene-ontology, kegg, reactome]
created: 2026-08-17
updated: 2026-08-17
sources:
  - "[[tnbc-analytics-suite]]"
  - "[[statistical-methods-primer]]"
aliases: [Pathway Enrichment, ORA, GSEA, clusterProfiler Analysis]
---

# Pathway Enrichment Analysis (ORA vs. GSEA)

## Overview
High-throughput transcriptomics produces lists of differentially expressed genes (DEGs). **Pathway Enrichment Analysis** statistically evaluates whether predefined biological gene sets (Gene Ontology terms, KEGG pathways, Reactome cascades, MSigDB Hallmark signatures) are significantly over-represented among the detected genes.

## Methodological Approaches

### 1. Over-Representation Analysis (ORA)
- **Input**: Hard-thresholded DEG list (e.g., $|\log_2\text{FC}| > 2, q < 0.05$).
- **Statistical Test**: Hypergeometric distribution / Fisher's Exact Test:
  $$P(X \ge k) = \sum_{i=k}^{\min(n, K)} \frac{\binom{K}{i} \binom{N-K}{n-i}}{\binom{N}{n}}$$
  where $N$ is total genome background, $K$ is total genes in the pathway, $n$ is total DEGs, and $k$ is DEGs in the pathway.
- **Package**: `clusterProfiler::enrichGO`, `enrichKEGG`.

### 2. Gene Set Enrichment Analysis (GSEA)
- **Input**: Non-thresholded, full ranked gene list ordered by differential metric (e.g., $\text{sign}(\log_2\text{FC}) \times -\log_{10}(\text{padj})$).
- **Algorithm**: Calculates an Enrichment Score (ES) using a weighted Kolmogorov-Smirnov running sum statistic.
- **Advantage**: Detects subtle, coordinated pathway-level shifts without arbitrary p-value cutoffs.

## Application in TNBC Workflow
- In `run_phase6_enrichment.R` and `run_reactome.R`, ORA identified hyperactivation of the PI3K-Akt signaling cascade, mitotic cell cycle checkpoint bypass, and extracellular matrix remodeling in [[triple-negative-breast-cancer]].

## Related Concepts & Entities
- [[clusterprofiler]]
- [[pi3k-akt-mtor-signaling-pathway]]
- [[deseq2]]
- [[triple-negative-breast-cancer]]

## Citations & Sources
- [[tnbc-analytics-suite]]
