---
type: entity
title: HCC1395 / HCC1395BL Breast Cancer Genomic Benchmark Pair
tags: [cell-line, tnbc, seqc2, benchmark, genomics, rna-seq]
created: 2026-08-17
updated: 2026-08-17
sources:
  - "[[hcc1395-rnaseq-benchmark]]"
aliases: [HCC1395, HCC1395BL, SEQC2 Cell Line]
---

# HCC1395 / HCC1395BL Breast Cancer Genomic Benchmark Pair

## Overview
**HCC1395** is a primary human Triple-Negative Breast Cancer cell line established from a 43-year-old female patient (ductal carcinoma, Stage 1, Grade 3). **HCC1395BL** is the matched normal lymphoblastoid cell line established from peripheral blood lymphocytes of the same patient.

## Reference Standard in Bioinformatics
- Adopted by the FDA-led **SEQC2 (Sequencing Quality Control Phase II)** Consortium as the international reference standard for somatic and germline variant detection, bulk RNA-Seq transcript quantification, and reproducibility benchmarking.

## Repository Implementation
- Analyzed in `exploratory_analysis/practice_modules/module_2/scripts/deseq2_analysis.R` to validate HISAT2 alignment fidelity, DESeq2 size-factor scaling, and distance clustering.

## Related Concepts & Entities
- [[hcc1395-rnaseq-benchmark]]
- [[deseq2]]
- [[triple-negative-breast-cancer]]
