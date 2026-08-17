---
type: entity
title: NIH ClinVar Database
tags: [databases, ncbi, clinvar, clinical-genomics, variants, pathogenicity]
created: 2026-08-17
updated: 2026-08-17
sources:
  - "[[clinical-genomics-pipeline]]"
aliases: [ClinVar, NCBI ClinVar]
---

# NIH ClinVar Database

## Overview
**ClinVar** is a freely accessible public archive maintained by the National Center for Biotechnology Information (NCBI) at the National Institutes of Health (NIH). It aggregates relationships among human genetic variations and phenotypes, along with supporting clinical evidence and classification tiers (Pathogenic, Benign, VUS, Drug Response).

## Core Data Schema & Attributes
- **`RS# (dbSNP)`**: Standardized Reference SNP cluster identifier.
- **`ClinicalSignificance`**: ACMG clinical classification (Pathogenic, Likely Pathogenic, Benign, Drug Response, Risk Factor).
- **`PhenotypeList`**: Associated diseases, syndromes, or drug response categories.
- **`ReviewStatus`**: Star rating system indicating assertion criteria and expert panel review consensus.

## Repository Implementation
- The script `clinvar_annotator.py` downloads `variant_summary.txt.gz`, filters for non-zero RSIDs, and performs an inner merge with personal 609K SNP array CSV data to extract medically actionable findings.

## Related Concepts & Entities
- [[variant-pathogenicity-annotation]]
- [[clinical-genomics-pipeline]]
