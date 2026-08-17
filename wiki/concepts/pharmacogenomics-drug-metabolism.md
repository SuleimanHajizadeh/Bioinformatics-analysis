---
type: concept
title: Clinical Pharmacogenomics (PGx) and Drug Metabolism
tags: [pharmacogenomics, pgx, cpic, metabolism, cyp450, precision-medicine]
created: 2026-08-17
updated: 2026-08-17
sources:
  - "[[clinical-genomics-pipeline]]"
aliases: [Pharmacogenomics, PGx, Drug Metabolism, CPIC Guidelines]
---

# Clinical Pharmacogenomics (PGx) and Drug Metabolism

## Overview
**Pharmacogenomics (PGx)** investigates how individual genomic variations affect drug absorption, distribution, metabolism, and excretion (ADME), dictating therapeutic efficacy and adverse drug reaction (ADR) risk according to Clinical Pharmacogenetics Implementation Consortium (CPIC) guidelines.

## Key Enzymatic Families & Genomic Loci

### 1. Phase I Cytochrome P450 Enzymes (CYP450)
- **[[cyp1a2]] (`rs762551`, *1F allele)**: Alters caffeine and theophylline metabolic clearance; high inductivity in A/A fast metabolizers.
- **`CYP2D6` / `CYP2C19`**: Regulate activation of prodrugs (e.g., Codeine, Tamoxifen) and clearance of cardiovascular and psychiatric medications.

### 2. Transporters and Neurotransmission
- **`SLCO1B1` (`rs4149056`, 521T>C)**: Hepatic uptake transporter; C allele causes statin accumulation and severe statin-induced myopathy.
- **[[comt]] (`rs4680`, Val158Met)**: Degrades synaptic catecholamines (dopamine/norepinephrine); Met/Met (A/A) "Worrier" genotype exhibits 3-4x reduced enzymatic activity.

## Analytical Implementation in Repository
- `analyze_dna.py` directly infers metabolic phenotypes for *CYP1A2*, *COMT*, *SLCO1B1*, and *OPRM1* from personal array chips.

## Related Concepts & Entities
- [[variant-pathogenicity-annotation]]
- [[cyp1a2]]
- [[comt]]
- [[clinvar]]

## Citations & Sources
- [[clinical-genomics-pipeline]]
