---
type: entity
title: clusterProfiler Bioconductor Package
tags: [tools, software, bioconductor, r, pathway-enrichment, gsea, gene-ontology]
created: 2026-08-17
updated: 2026-08-17
sources:
  - "[[tnbc-analytics-suite]]"
aliases: [clusterProfiler, Yu et al. 2012]
---

# clusterProfiler Bioconductor Package

## Overview
**clusterProfiler** is an R/Bioconductor package developed by Guangchuang Yu for comparing and visualising functional profiles among gene clusters using Gene Ontology (GO), KEGG, Reactome, and MSigDB gene sets.

## Key Capabilities
- Over-representation analysis (`enrichGO`, `enrichKEGG`).
- Universal Gene Set Enrichment Analysis (`GSEA`).
- Rich visualization suite (`dotplot`, `cnetplot`, `emapplot`, `gseaplot2`).

## Repository Implementation
- Integrated into `run_phase6_enrichment.R` and `run_immune_kegg.R` to decipher dysregulated cancer cascades downstream of [[deseq2]] results.

## Related Concepts & Entities
- [[pathway-enrichment-gsea]]
- [[deseq2]]
- [[triple-negative-breast-cancer]]
