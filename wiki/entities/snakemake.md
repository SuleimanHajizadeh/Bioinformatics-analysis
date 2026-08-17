---
type: entity
title: Snakemake Workflow Management System
tags: [tools, workflow-engine, python, reproducibility, pipelines]
created: 2026-08-17
updated: 2026-08-17
sources:
  - "[[bulk-rnaseq-snakemake-pipeline]]"
aliases: [Snakemake, Workflow Engine]
---

# Snakemake Workflow Management System

## Overview
**Snakemake** is a Python-based workflow management engine designed to create reproducible and scalable data analyses. It uses a human-readable domain-specific language (DSL) based on Python to specify rules, inputs, outputs, and shell executions.

## Key Capabilities
- **Directed Acyclic Graph (DAG)**: Automatically calculates dependency trees based on filenames and wildcards.
- **Resilience & Caching**: Reruns only rules whose input files have changed or whose targets are missing.
- **Resource Management**: Dynamic scaling across local multicore CPUs, SLURM clusters, or cloud instances.

## Repository Implementation
- Orchestrates the primary bulk RNA-Seq pipeline (`workflows/bulk_rnaseq/Snakefile`) from raw SRA fetching through Trimmomatic, HISAT2 alignment, featureCounts, and DESeq2 handoff.

## Related Concepts & Entities
- [[hisat2]]
- [[featurecounts]]
- [[deseq2]]
- [[bulk-rnaseq-snakemake-pipeline]]
