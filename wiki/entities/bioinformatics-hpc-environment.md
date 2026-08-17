---
type: entity
title: Bioinformatics High-Performance Computing & Software Environment
tags: [infrastructure, hpc, centos, linux, conda, compilers, environment]
created: 2026-08-17
updated: 2026-08-17
sources:
  - "[[statistical-methods-primer]]"
  - "[[bulk-rnaseq-snakemake-pipeline]]"
aliases: [HPC Environment, CentOS AlmaLinux System, Bioinformatics Server Stack]
---

# Bioinformatics High-Performance Computing & Software Environment

## Overview
The computational infrastructure supporting this bioinformatics repository consists of an Enterprise Linux distribution (**CentOS / AlmaLinux 9.x**) equipped with development toolchains, mathematical acceleration libraries (BLAS/LAPACK, Armadillo, SuperLU, ARPACK), OpenMPI parallel processing, and R/Python conda environments.

## Core Software & Library Stack
- **Compilers & Math Acceleration**: GCC, Clang, GFortran, Armadillo 12.6, SuperLU 6.0, ARPACK 3.8.
- **Workflow & Container Engines**: [[snakemake]], Podman/Docker, Apptainer/Singularity.
- **High-Throughput Sequencing Aligners & Tools**: [[hisat2]], [[featurecounts]], `samtools`, `bedtools`, `Trimmomatic`, `FastQC`, `sra-tools`.
- **R / Bioconductor Environment**: [[deseq2]], [[seurat]], `WGCNA`, `clusterProfiler`, `pheatmap`, `ggplot2`.
- **Python ML & Structural Stack**: `scikit-learn`, `pandas`, `numpy`, [[alphafold2]], [[autodock-vina]], `biopython`, `PyMOL`.

## Related Concepts & Entities
- [[snakemake]]
- [[hisat2]]
- [[featurecounts]]
- [[deseq2]]
- [[alphafold2]]
- [[autodock-vina]]
