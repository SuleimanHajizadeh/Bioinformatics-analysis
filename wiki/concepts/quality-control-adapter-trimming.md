---
type: concept
title: Next-Generation Sequencing Quality Control & Adapter Trimming
tags: [ngs, quality-control, fastqc, trimmomatic, preprocessing, phred-score]
created: 2026-08-17
updated: 2026-08-17
sources:
  - "[[bulk-rnaseq-snakemake-pipeline]]"
aliases: [NGS Quality Control, Read Trimming, Phred Scoring, Adapter Removal]
---

# Next-Generation Sequencing Quality Control & Adapter Trimming

## Overview
Raw high-throughput sequencing reads contain technical artifacts, base-calling errors, and synthetic adapter sequences that compromise alignment accuracy. Systematic pre-alignment quality control and trimming ensure high mapping fidelity.

## Key Metrics & Mathematical Definition

### 1. Phred Quality Score ($Q$)
Quantifies the probability $P$ of an incorrect base call:
$$Q = -10 \log_{10}(P) \iff P = 10^{-Q/10}$$
- **$Q = 20$**: $1\%$ error probability ($99\%$ base call accuracy).
- **$Q = 30$**: $0.1\%$ error probability ($99.9\%$ base call accuracy - gold standard).

### 2. Trimming Strategy in Snakemake Workflow
Implemented via `Trimmomatic PE`:
- `ILLUMINACLIP:adapters.fa:2:30:10`: Identifies and cuts standard TruSeq/Illumina adapter sequences with a seed mismatch limit of 2 and palindrome clip threshold of 30.
- `LEADING:3 TRAILING:3`: Cuts low-quality bases at read termini ($Q < 3$).
- `SLIDINGWINDOW:4:15`: Scans read in 4-base windows, clipping when mean quality drops below $Q = 15$.
- `MINLEN:36`: Discards trimmed reads shorter than 36 bp to prevent spurious non-specific genomic alignments.

## Related Concepts & Entities
- [[fastqc-trimmomatic]]
- [[hisat2]]
- [[snakemake]]
- [[bulk-rnaseq-snakemake-pipeline]]

## Citations & Sources
- [[bulk-rnaseq-snakemake-pipeline]]
