---
type: entity
title: featureCounts Read Summarization Utility
tags: [tools, subread, quantification, rna-seq, counting]
created: 2026-08-17
updated: 2026-08-17
sources:
  - "[[bulk-rnaseq-snakemake-pipeline]]"
aliases: [featureCounts, Subread featureCounts]
---

# featureCounts Read Summarization Utility

## Overview
**featureCounts** is an ultrafast and highly efficient read summarization program developed as part of the Subread software suite (Liao et al.). It counts mapped sequencing reads against genomic features such as genes, exons, promoters, and genomic intervals defined in GTF/GFF format.

## Key Properties
- Orders of magnitude faster than `HTSeq-count` while maintaining low memory footprint.
- Native support for paired-end alignments (`-p`), multi-threading (`-T`), and multi-mapping filtering.

## Repository Implementation
- Rule `count_features` in `workflows/bulk_rnaseq/Snakefile` and standalone script `run_featurecounts.sh` quantify mapped BAM files into `counts.txt` for input to [[deseq2]].

## Related Concepts & Entities
- [[snakemake]]
- [[hisat2]]
- [[deseq2]]
- [[bulk-rnaseq-snakemake-pipeline]]
