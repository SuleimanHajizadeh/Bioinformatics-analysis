---
type: entity
title: HISAT2 Spliced Read Aligner
tags: [tools, alignment, genomics, rna-seq, hisat2]
created: 2026-08-17
updated: 2026-08-17
sources:
  - "[[bulk-rnaseq-snakemake-pipeline]]"
aliases: [HISAT2, Hierarchical Indexing for Spliced Alignment of Transcripts]
---

# HISAT2 Spliced Read Aligner

## Overview
**HISAT2** (Hierarchical Indexing for Spliced Alignment of Transcripts 2) is a fast and sensitive alignment program for mapping next-generation sequencing reads (both DNA and RNA) against a reference genome. It uses a graph-based indexing structure based on the Hierarchical Graph Ferragina-Manzini (HGFM) index.

## Key Capabilities
- Accurately aligns spliced RNA-seq reads across exon-exon junctions.
- Highly memory-efficient index supporting human reference genome GRCh38.

## Repository Implementation
- Integrated into `workflows/bulk_rnaseq/Snakefile` (rule `align_hisat2`) and bash script `run_hisat2.sh`, piping raw alignments directly into `samtools` for coordinate sorting and index generation.

## Related Concepts & Entities
- [[snakemake]]
- [[featurecounts]]
- [[bulk-rnaseq-snakemake-pipeline]]
