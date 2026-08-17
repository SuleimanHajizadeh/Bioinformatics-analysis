---
type: entity
title: Triple-Negative Breast Cancer (TNBC)
tags: [oncology, breast-cancer, tnbc, clinical-genomics, biomarkers]
created: 2026-08-17
updated: 2026-08-17
sources:
  - "[[statistical-methods-primer]]"
aliases: [TNBC, Basal-like Breast Cancer]
---

# Triple-Negative Breast Cancer (TNBC)

## Overview
**Triple-Negative Breast Cancer (TNBC)** represents approximately 15–20% of all invasive breast carcinomas and is characterized by the absence of immunohistochemical expression of Estrogen Receptor (*ESR1*), Progesterone Receptor (*PGR*), and lack of Human Epidermal Growth Factor Receptor 2 (*ERBB2* / HER2) gene amplification.

## Molecular Characteristics & Expression Signatures
- **Negative Receptors**: *ESR1* (-), *PGR* (-), *ERBB2* (-).
- **Common Upregulated Hubs / Biomarkers**:
  - *MKI67* (Ki-67 proliferation index marker)
  - *EGFR* (Epidermal Growth Factor Receptor)
  - *AKT1* (PI3K/Akt pathway central kinase)
- **Clinical Prognosis**: High proliferative index, aggressive clinical trajectory, elevated risk of early visceral recurrence, and lack of targeted hormonal / anti-HER2 therapies.

## Repository Analytical Focus
- Characterized in dataset **GSE183947** with 1,847 significant DEGs identified using [[deseq2]].
- Co-expression network analysis via [[wgcna-coexpression-networks]] identified the *AKT1*-driven turquoise module.

## Related Concepts & Entities
- [[wgcna-coexpression-networks]]
- [[deseq2]]
- [[transcriptomic-analysis-framework]]

## Citations & Sources
- [[statistical-methods-primer]]
