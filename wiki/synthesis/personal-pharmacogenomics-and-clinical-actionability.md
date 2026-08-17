---
type: synthesis
title: Personal Pharmacogenomics and Direct-to-Consumer Clinical Actionability
tags: [synthesis, pharmacogenomics, precision-medicine, cpic, clinical-actionability]
created: 2026-08-17
updated: 2026-08-17
sources:
  - "[[clinical-genomics-pipeline]]"
  - "[[statistical-methods-primer]]"
aliases: [Personal PGx Synthesis, DTC Genomics Clinical Interpretation]
---

# Personal Pharmacogenomics and Direct-to-Consumer Clinical Actionability

## Overview
Direct-to-Consumer (DTC) SNP microarrays generate over 600,000 genotyped loci per individual. This synthesis outlines the computational translation of raw consumer arrays into standardized clinical interpretations via [[clinvar]], CPIC guidelines, and specialized pharmacogenomic panels.

```mermaid
graph TD
    A[Raw DTC Microarray: 609K SNPs] --> B[[clinical-genomics-pipeline]]
    B --> C[Behavioral & Lifestyle Traits: [[cyp1a2]], [[comt]], *ACTN3*]
    B --> D[ClinVar Medical Overlap: 21,894 Loci]
    D --> E[ACMG High-Risk Pathogenic Alleles: [[brca1-brca2]]]
    D --> F[Pharmacogenomic ADME Alerts: *SLCO1B1*, *DPYD*, *TPMT*]
    E --> G[Clinical Genetics Referral]
    F --> H[EHR Pharmacogenomics Prescription Guardrails]
```

## Clinical Stratification Table

| Category | Target Gene | Representative rsID | Phenotype / Impact | Actionability Level |
| :--- | :--- | :--- | :--- | :--- |
| **Oncological Risk** | [[brca1-brca2]] | ClinVar Pathogenic variants | Loss of homologous recombination repair; high lifetime breast/ovarian cancer risk | **Tier 1 (High Actionability)** |
| **Pharmacogenomics** | `SLCO1B1` | `rs4149056` (521T>C) | Impaired hepatic statin uptake; elevated risk of rhabdomyolysis | **CPIC Level A (Dosing Adjustment)** |
| **Metabolic Trait** | [[cyp1a2]] | `rs762551` (*1F) | Fast vs. Slow caffeine & theophylline clearance | **Lifestyle / Dietary Advice** |
| **Neurochemical** | [[comt]] | `rs4680` (Val158Met) | "Worrier" (Met/Met) vs. "Warrior" (Val/Val) dopamine metabolism | **Cognitive & Ergonomic Planning** |

## Related Wiki Pages
- [[clinical-genomics-pipeline]]
- [[pharmacogenomics-drug-metabolism]]
- [[variant-pathogenicity-annotation]]
- [[clinvar]]
- [[brca1-brca2]]
- [[cyp1a2]]
- [[comt]]
