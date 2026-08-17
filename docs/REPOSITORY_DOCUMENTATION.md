# 🧬 Comprehensive Repository Architecture & File Documentation

This document provides a comprehensive technical index and operational guide for every directory, pipeline, script, configuration file, and artifact contained within this bioinformatics repository.

---

## 🗺️ Master Directory Layout

```text
Bioinformatics-analysis/
├── AGENTS.md                                # LLM Wiki schema, operational conventions & prompt guidelines
├── README.md                                # Repository landing page, project badges & high-level highlights
├── LICENSE                                  # MIT Open Source License
├── requirements.txt                         # Python runtime dependencies
├── environment.yml                          # Conda virtual environment definition (R + Python + Bioconductor)
├── .gitignore                               # Git exclusions (BAM, FASTQ, large datasets, .obsidian)
│
├── raw/                                     # [Layer 1: Raw Sources] Immutable source material
│   ├── README.md                            # Guidelines for raw source curation and image storage
│   └── assets/                              # Local figures, diagrams, and image attachments
│
├── wiki/                                    # [Layer 2: Persistent Knowledge Base] Obsidian LLM Wiki
│   ├── index.md                             # Master catalog of all 50+ concepts, entities, syntheses
│   ├── log.md                               # Append-only chronological operational and ingestion log
│   ├── concepts/                            # Theoretical, mathematical & computational foundations (15 pages)
│   ├── entities/                            # Tools, packages, databases, genes, drugs, cell lines (20 pages)
│   ├── summaries/                           # Per-source summaries with backlinks to raw/ & codebase (7 pages)
│   └── synthesis/                           # Cross-layer translational syntheses & therapeutic matrices (5 pages)
│
├── workflows/                               # [Production Analysis Pipelines]
│   ├── bulk_rnaseq/                         # End-to-End Bulk RNA-Seq Pipeline
│   │   ├── Snakefile                        # Snakemake workflow orchestration (SRA to DEGs)
│   │   ├── config.yaml                      # Workflow configuration (sample IDs, reference indices, GTF)
│   │   ├── scripts/                         # Multi-phase execution and analytical R/Bash scripts
│   │   │   ├── setup_rna_seq.sh             # Environment setup, directory provisioning, index verification
│   │   │   ├── download_sra.sh / download_tnbc.sh / download_ena.sh / download_20_tnbc.sh # SRA fetching
│   │   │   ├── run_qc.sh                    # FastQC execution wrapper
│   │   │   ├── run_trimmomatic.sh           # Trimmomatic PE quality filtering and adapter removal
│   │   │   ├── run_hisat2.sh                # HISAT2 spliced alignment & samtools sorting/indexing
│   │   │   ├── run_featurecounts.sh         # featureCounts quantification against GRCh38 GTF
│   │   │   ├── run_deseq.sh                 # DESeq2 command-line driver
│   │   │   ├── run_deseq2_analysis.R        # [Phase 5] DESeq2 Wald test, PCA, Top 50 DEG heatmap
│   │   │   ├── run_phase6_enrichment.R      # [Phase 6] Volcano plot & clusterProfiler GO/KEGG enrichment
│   │   │   ├── run_reactome.R               # [Phase 6] Reactome pathway analysis and visualization
│   │   │   ├── run_immune_kegg.R            # [Phase 6] Immune signaling and cytokine pathway enrichment
│   │   │   ├── run_phase7_ml.R              # [Phase 7] Initial Random Forest machine learning classifier
│   │   │   ├── run_phase7_nested_cv.R       # [Phase 7] Leakage-free Nested LOOCV + Permutation Testing
│   │   │   ├── run_phase8_network.R         # [Phase 8] iGraph/WGCNA co-expression network & hub extraction
│   │   │   ├── run_phase9_deconv.R          # [Phase 9] Digital pathology TME immune deconvolution
│   │   │   └── plot_rnaseq.py               # Python visualization routines
│   │   └── results/                         # Generated plots, CSV tables, and DEG matrices
│   │
│   └── single_cell/                         # Single-Cell RNA-Seq (scRNA-Seq) Workflow
│       └── scripts/
│           ├── seurat_qc.R                  # 10x PBMC 3k cell QC & mitochondrial filter (< 5% MT)
│           ├── seurat_cluster.R             # VST 2,000 features, PCA (10 dims), Louvain clustering, UMAP
│           ├── seurat_markers.R             # Wilcoxon differential marker discovery per cluster
│           └── seurat_annotate.R            # Canonical PBMC cell-type identity annotation
│
├── clinical_genomics/                       # [Clinical Genomics & Personal Genotyping]
│   └── myheritage_snp/
│       ├── README.md                        # Overview of DTC microarray processing & ClinVar annotation
│       └── scripts/
│           ├── analyze_dna.py               # High-speed hash matcher for 20+ lifestyle, metabolic & behavioral traits
│           └── clinvar_annotator.py         # Full 609K SNP array inner join with NIH ClinVar (21,894 clinical matches)
│
├── exploratory_analysis/                    # [Benchmarking & Quality Control]
│   ├── fastqc_reports/                      # Pre- and post-trimmed read quality reports (HTML & zip)
│   └── practice_modules/
│       └── module_2/                        # SEQC2 gold-standard reference benchmark (HCC1395 vs HCC1395BL)
│           └── scripts/
│               ├── deseq2_analysis.R        # Reference sample distance matrix & PCA validation
│               └── process_new_samples.sh   # Automated batch processing script
│
├── manuscripts/                             # [Academic Research Manuscripts]
│   ├── manuscript_draft.md                  # Complete PLOS ONE paper draft (AKT1 hub, AlphaFold2, AutoDock Vina)
│   └── PLOS_ONE_SUBMISSION/                 # Submission package (cover letter, checklist, Word draft)
│
└── docs/                                    # [Theoretical Primers & System Logs]
    ├── STATISTICAL_METHODS.md               # Negative Binomial math, FDR correction, VST vs rlog normalization
    ├── REPOSITORY_DOCUMENTATION.md          # [This file] Master codebase and file index
    └── installed applications of centos.txt # HPC environment package specification (AlmaLinux 9 / CentOS)
```

---

## 🔬 Core Workflows and Technical Details

### 1. Bulk RNA-Seq Pipeline (`workflows/bulk_rnaseq/`)
- **Input**: SRA Accession numbers (e.g., GSE183947 / SRR16800543) or raw paired-end FASTQ reads.
- **Workflow Steps**:
  1. `download_sra.sh`: Downloads raw SRA data and splits into paired FASTQ files via `fasterq-dump`.
  2. `run_qc.sh`: Assesses sequencing quality with `fastqc`.
  3. `run_trimmomatic.sh`: Removes Illumina adapter sequences and low-quality bases ($Q < 15$ sliding window).
  4. `run_hisat2.sh`: Aligns reads against human reference genome GRCh38 using a splice-aware graph FM index.
  5. `run_featurecounts.sh`: Quantifies mapped reads against the GENCODE/Ensembl GTF annotation matrix.
  6. `run_deseq2_analysis.R`: Normalizes counts using median-of-ratios, fits a Negative Binomial GLM, estimates dispersion with empirical Bayes shrinkage, and extracts DEGs ($q < 0.05, |\log_2\text{FC}| > 2$).
  7. `run_phase6_enrichment.R`: Performs GO biological process and KEGG pathway enrichment via `clusterProfiler`.
  8. `run_phase7_nested_cv.R`: Implements Nested Leave-One-Out Cross-Validation (LOOCV) to prevent data leakage during Random Forest biomarker selection, validated via 100 label permutation tests ($\text{AUC} \approx 0.91$).
  9. `run_phase8_network.R`: Computes topological co-expression adjacency ($|R| \ge 0.85$) to isolate the **AKT1** master hub kinase ($k_{\text{in}} = 0.94$).
  10. `run_phase9_deconv.R`: Scores tumor microenvironment infiltration across CD8+ T-cells, B-cells, and M1/M2 macrophages using geometric mean Z-scores.

### 2. Single-Cell Transcriptomics (`workflows/single_cell/`)
- **Input**: 10x Genomics PBMC 3k feature-barcode matrix (`filtered_gene_bc_matrices/hg19/`).
- **Workflow Steps**:
  1. `seurat_qc.R`: Filters dying cells ($\text{percent.mt} > 5\%$) and empty droplets/multiplets ($200 < \text{nFeature} < 2500$).
  2. `seurat_cluster.R`: Applies `LogNormalize`, identifies top 2,000 variable genes using VST, runs PCA, constructs a Shared Nearest Neighbor (SNN) graph on 10 PCs, and detects cell clusters using Louvain modularity optimization (resolution 0.5) projected onto UMAP.
  3. `seurat_markers.R` & `seurat_annotate.R`: Finds differential biomarker genes per cluster using Wilcoxon rank-sum testing and assigns PBMC cell-type annotations (CD4+ T-cells, CD8+ T-cells, B-cells, Monocytes, NK cells, Dendritic cells).

### 3. Clinical Genomics & DTC SNP Annotation (`clinical_genomics/myheritage_snp/`)
- **Input**: Raw 609K SNP microarray CSV file (`MyHeritage_raw_dna_data.csv`).
- **Workflow Steps**:
  1. `analyze_dna.py`: Scans the user's genotype against a curated dictionary of actionable metabolic, behavioral, and athletic polymorphisms (*ACTN3*, *MCM6*, *CYP1A2*, *COMT*, *DRD2*, *OPRM1*, *FTO*, *CCR5-Delta32*).
  2. `clinvar_annotator.py`: Downloads the NIH ClinVar `variant_summary.txt.gz` archive, filters for valid rsIDs, and performs an inner merge on 609,000 loci to extract 21,894 clinical matches classified by ACMG tiers (Pathogenic, Drug Response, Risk Factor, Protective).

### 4. Academic Research Manuscript (`manuscripts/`)
- **Focus**: Multi-Scale Computational Oncology Pipeline in Triple-Negative Breast Cancer (TNBC).
- **Core Findings**:
  - Identified 2,418 significant DEGs in GSE183947.
  - WGCNA turquoise module ($r=0.68, P<0.001$) identified **AKT1** as the master oncogenic hub.
  - AlphaFold2 predicted full-length AKT1 3D coordinates (mean pLDDT 92.4, 94.2% Ramachandran favored).
  - AutoDock Vina 1.2.5 demonstrated high binding affinity for the allosteric inhibitor **MK-2206** ($\Delta G = -9.4\text{ kcal/mol}$, Trp80 and Asp292 H-bonds) over ATP-competitive Ipatasertib ($\Delta G = -8.2\text{ kcal/mol}$).
