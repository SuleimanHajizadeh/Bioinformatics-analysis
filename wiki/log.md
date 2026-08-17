# Wiki Activity Log

This is an append-only chronological log of all ingestion, query synthesis, and maintenance operations.

Each entry follows the parseable format:
`## [YYYY-MM-DD] <action> | <Title>`

---

## [2026-08-17] init | LLM Wiki Architecture Setup
- Initialized 3-layer architecture (`AGENTS.md`, `raw/`, `wiki/`).
- Configured folder structure: `wiki/concepts/`, `wiki/entities/`, `wiki/summaries/`, `wiki/synthesis/`, `raw/assets/`.
- Configured Obsidian-compatible wikilinks and YAML frontmatter specifications in [AGENTS.md](file:///Users/macbookairm2/Documents/GitHub/Bioinformatics-analysis/AGENTS.md).

## [2026-08-17] ingest | Statistical Methodology Primer
- **Source**: `docs/STATISTICAL_METHODS.md`
- **Summary**: [[statistical-methods-primer]]
- **Entities/Concepts Touched/Created**:
  - Concept: [[negative-binomial-model]]
  - Concept: [[multiple-testing-fdr]]
  - Concept: [[variance-stabilizing-transformation]]
  - Concept: [[wgcna-coexpression-networks]]
  - Entity: [[deseq2]]
  - Entity: [[triple-negative-breast-cancer]]
  - Synthesis: [[transcriptomic-analysis-framework]]
- **Updated Catalog**: [[index|wiki/index.md]]

## [2026-08-17] ingest | Codebase - Bulk RNA-Seq Pipeline & Snakemake
- **Source**: `workflows/bulk_rnaseq/Snakefile` & `config.yaml`
- **Summary**: [[bulk-rnaseq-snakemake-pipeline]]
- **Entities/Concepts Touched/Created**:
  - Entity: [[snakemake]]
  - Entity: [[hisat2]]
  - Entity: [[featurecounts]]

## [2026-08-17] ingest | Codebase - TNBC Advanced Analytics Suite (Phases 5-9)
- **Source**: `workflows/bulk_rnaseq/scripts/`
- **Summary**: [[tnbc-analytics-suite]]
- **Entities/Concepts Touched/Created**:
  - Concept: [[nested-cross-validation-ml]]
  - Concept: [[immune-deconvolution-tme]]
  - Entity: [[random-forest-classifier]]

## [2026-08-17] ingest | Codebase - Single-Cell Seurat Pipeline
- **Source**: `workflows/single_cell/scripts/`
- **Summary**: [[single-cell-seurat-workflow]]
- **Entities/Concepts Touched/Created**:
  - Concept: [[single-cell-clustering-umap]]
  - Entity: [[seurat]]

## [2026-08-17] ingest | Codebase - Clinical Genomics & ClinVar Parser
- **Source**: `clinical_genomics/myheritage_snp/scripts/`
- **Summary**: [[clinical-genomics-pipeline]]
- **Entities/Concepts Touched/Created**:
  - Concept: [[variant-pathogenicity-annotation]]
  - Entity: [[clinvar]]

## [2026-08-17] synthesis | Multi-Omics and Translational Bioinformatics Ecosystem
- Generated overarching cross-pipeline synthesis: [[multi-omics-analytical-ecosystem]]
- Linked bulk RNA-seq, single-cell resolution transcriptomics, and personal clinical genomics.

## [2026-08-17] ingest | PLOS ONE Manuscript Draft - Multi-Scale TNBC & AKT1 Pipeline
- **Source**: `manuscripts/manuscript_draft.md`
- **Summary**: [[tnbc-akt1-alphafold-docking-manuscript]]
- **Entities/Concepts Touched/Created**:
  - Concept: [[alphafold-structural-modeling]]
  - Concept: [[molecular-docking-virtual-screening]]
  - Entity: [[akt1-kinase]]
  - Entity: [[alphafold2]]
  - Entity: [[autodock-vina]]
  - Entity: [[mk-2206]]
  - Synthesis: [[transcriptome-to-structure-drug-discovery]]
- **Updated Catalog**: [[index|wiki/index.md]]

## [2026-08-17] ingest | HPC Compute Infrastructure & Software Stack
- **Source**: `docs/installed applications of centos.txt`
- **Entities Touched/Created**: [[bioinformatics-hpc-environment]]
- **Updated Catalog**: [[index|wiki/index.md]]

## [2026-08-17] lint | Comprehensive Graph Health Audit
- **Total Pages Scanned**: 34 markdown pages
- **Broken Links**: 0 broken `[[wikilinks]]`
- **Orphan Nodes**: 0 isolated pages (all pages have $\ge 3$ inbound connections)
- **Top Graph Hubs**: `deseq2` (21), `triple-negative-breast-cancer` (18), `statistical-methods-primer` (16), `tnbc-akt1-alphafold-docking-manuscript` (15), `wgcna-coexpression-networks` (14), `akt1-kinase` (14), `variance-stabilizing-transformation` (13).
