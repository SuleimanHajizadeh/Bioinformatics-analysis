---
type: entity
title: FastQC & Trimmomatic Pre-processing Utilities
tags: [tools, quality-control, preprocessing, fastqc, trimmomatic, ngs]
created: 2026-08-17
updated: 2026-08-17
sources:
  - "[[bulk-rnaseq-snakemake-pipeline]]"
aliases: [FastQC, Trimmomatic, Pre-processing Pipeline]
---

# FastQC & Trimmomatic Pre-processing Utilities

## Overview
**FastQC** (Babraham Bioinformatics) provides comprehensive visual quality metrics on raw high-throughput sequencing data (per-base quality, GC content, sequence duplication, adapter contamination). **Trimmomatic** (Bolger et al.) is a flexible paired-end read trimmer designed specifically for Illumina NGS data.

## Repository Implementation
- Configured in `workflows/bulk_rnaseq/Snakefile` (rules `fastqc` and `trim_reads`) and standalone shell scripts `run_qc.sh` and `run_trimmomatic.sh`.

## Related Concepts & Entities
- [[quality-control-adapter-trimming]]
- [[hisat2]]
- [[snakemake]]
