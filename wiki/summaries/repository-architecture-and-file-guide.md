---
type: summary
title: "Master Repository Architecture and File Documentation Guide"
authors: "Bioinformatics Pipeline Team"
year: 2026
raw_source: "docs/REPOSITORY_DOCUMENTATION.md"
ingested: 2026-08-17
tags: [repository, architecture, documentation, snakemake, seurat, clinvar, pipelines]
---

# Master Repository Architecture and File Documentation Guide

## Core Thesis / Executive Summary
A complete architectural reference mapping every directory, execution pipeline, script, configuration file, and manuscript in this repository. Bridges bulk transcriptomics (Snakemake, HISAT2, DESeq2, WGCNA, ML), single-cell genomics (Seurat), personal clinical genomics (ClinVar), and AI structural pharmacology (AlphaFold2, AutoDock Vina).

## Key Components Cataloged
1. **Bulk RNA-Seq Pipeline**: `workflows/bulk_rnaseq/` (10-stage workflow from SRA to immune deconvolution).
2. **Single-Cell Workflow**: `workflows/single_cell/` (Seurat QC, clustering, UMAP, marker annotation).
3. **Clinical Genotyping**: `clinical_genomics/myheritage_snp/` (`analyze_dna.py` & `clinvar_annotator.py`).
4. **Research Manuscript**: `manuscripts/manuscript_draft.md` (Multi-scale TNBC AKT1 inhibitor modeling).
5. **LLM Wiki**: `wiki/` (Persistent, compounding Obsidian knowledge base of 50+ interlinked pages).

## Touched Wiki Pages
- Concept: [[negative-binomial-model]]
- Concept: [[single-cell-clustering-umap]]
- Concept: [[variant-pathogenicity-annotation]]
- Entity: [[snakemake]]
- Entity: [[seurat]]
- Entity: [[clinvar]]
- Entity: [[akt1-kinase]]
- Synthesis: [[multi-omics-analytical-ecosystem]]
