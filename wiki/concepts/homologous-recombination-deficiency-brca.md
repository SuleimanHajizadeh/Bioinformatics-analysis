---
type: concept
title: Homologous Recombination Deficiency (HRD) and Synthetic Lethality
tags: [oncology, genetics, brca1, brca2, hrd, parp-inhibitors, dna-repair, synthetic-lethality]
created: 2026-08-17
updated: 2026-08-17
sources:
  - "[[tnbc-targeted-and-immune-therapeutics]]"
  - "[[clinical-genomics-pipeline]]"
aliases: [HRD, BRCAness, Synthetic Lethality, PARP Trapping]
---

# Homologous Recombination Deficiency (HRD) and Synthetic Lethality

## Overview
**Homologous Recombination Deficiency (HRD)** occurs when cells lose the ability to accurately repair DNA double-strand breaks (DSBs) through the high-fidelity homologous recombination pathway, predominantly due to germline or somatic mutations in [[brca1-brca2]], *PALB2*, or *RAD51C*. HRD is present in approximately 50% of [[triple-negative-breast-cancer]] tumors.

## Synthetic Lethality Mechanism

```mermaid
graph TD
    A[Endogenous / Chemotherapy DNA Single-Strand Breaks] -->|Repaired by| B[PARP1/2 Base Excision Repair]
    B -->|PARP Inhibitor: Olaparib/Talazoparib| C[Trapped PARP1-DNA Complexes]
    C -->|Replication Fork Collapse| D[DNA Double-Strand Breaks]
    D -->|In Normal Cells: Intact BRCA1/2| E[High-Fidelity HR Repair: Cell Survival]
    D -->|In HRD/TNBC Tumor Cells: Inactive BRCA1/2| F[Error-Prone NHEJ Repair: Chromosomal Catastrophe & Apoptosis]
```

## Genomic Biomarkers of HRD
1. **Loss of Heterozygosity (LOH)**: Chromosomal segments showing single parental allele loss.
2. **Telomeric Allelic Imbalance (TAI)**: Unequal parental contributions extending to telomeres.
3. **Large-Scale State Transitions (LST)**: Chromosomal breaks between adjacent regions $\ge 10\text{ Mb}$.
- Combined composite HRD score $\ge 42$ identifies patients responsive to platinum chemotherapy and PARP inhibitors.

## Related Concepts & Entities
- [[brca1-brca2]]
- [[triple-negative-breast-cancer]]
- [[tnbc-targeted-and-immune-therapeutics]]
- [[clinvar]]

## Citations & Sources
- [[tnbc-targeted-and-immune-therapeutics]]
