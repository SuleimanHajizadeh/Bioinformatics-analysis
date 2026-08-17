---
type: concept
title: Clinical Variant Pathogenicity Annotation
tags: [genomics, clinvar, dbsnp, snp-array, pathogenicity, acmg]
created: 2026-08-17
updated: 2026-08-17
sources:
  - "[[clinical-genomics-pipeline]]"
aliases: [ClinVar Annotation, ACMG Classification, Variant Interpretation]
---

# Clinical Variant Pathogenicity Annotation

## Overview
High-throughput personal genomics microarrays interrogate hundreds of thousands of single nucleotide polymorphisms (SNPs). **Variant Pathogenicity Annotation** systematically maps observed genomic alleles to standardized clinical archives (such as NIH ClinVar and dbSNP) to identify disease predispositions, carrier states, and pharmacogenomic drug response profiles according to ACMG/AMP clinical guidelines.

## Methodological Architecture

### 1. dbSNP Identifier Mapping
- Personal genotyping arrays (e.g., Illumina OmniExpress / Global Screening Array used by MyHeritage) report genetic variation indexed by **rsID** (e.g., `rs1815739`).
- Parsing strips non-numeric prefixes to allow high-speed numeric lookups against ClinVar database indices (`RS# (dbSNP)`).

### 2. Clinical Significance Categories (ClinVar)
Variants are stratified by ACMG clinical tiers:
- **Pathogenic / Likely Pathogenic**: Strong clinical evidence linking allele to monogenic disease risk or metabolic failure.
- **Drug Response / Pharmacogenomic**: Variants altering drug metabolism rates (e.g., *CYP1A2* caffeine/theophylline clearance, *SLCO1B1* statin-induced myopathy risk).
- **Risk Factor / Protective**: Polygenic modifiers influencing disease susceptibility (e.g., *FTO* adiposity risk, *CCR5-Delta32* HIV resistance).
- **Benign / Likely Benign**: Polymorphisms with high population frequencies and no deleterious functional effect.
- **Variant of Uncertain Significance (VUS)**: Inconclusive clinical evidence.

## Analytical Pipelines in Repository
- `analyze_dna.py`: Rapid targeted trait dictionary matching for lifestyle, athletic, and metabolic phenotypes.
- `clinvar_annotator.py`: Full 609K SNP array inner join with ClinVar bulk archive `variant_summary.txt.gz`.

## Related Concepts & Entities
- [[clinvar]]
- [[multi-omics-analytical-ecosystem]]

## Citations & Sources
- [[clinical-genomics-pipeline]]
