# 🧬 Comprehensive Transcriptomic Pipelines & Clinical Genomics Analysis

[![Bioinformatics](https://img.shields.io/badge/Bioinformatics-RNA--seq-blue.svg?style=flat-square)](https://www.ncbi.nlm.nih.gov/sra)
[![Method](https://img.shields.io/badge/Method-DESeq2_|_WGCNA_|_CIBERSORT-orange.svg?style=flat-square)](https://bioconductor.org/packages/DESeq2)
[![Tools](https://img.shields.io/badge/Tools-Bash_|_R_|_Python-green.svg?style=flat-square)](https://www.r-project.org/)
[![Dataset](https://img.shields.io/badge/Dataset-GEO%3A_GSE183947-red.svg?style=flat-square)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183947)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=flat-square)](https://opensource.org/licenses/MIT)

## 📌 Overview

This repository integrates high-throughput sequencing pipelines and specialized transcriptomic characterization workflows for **Triple-Negative Breast Cancer (TNBC)** molecular profiling. The primary dataset is the publicly available TCGA-BRCA RNA-seq cohort, accessed via GEO/SRA, enabling differential expression analysis, tumor microenvironment (TME) deconvolution, co-expression network analysis, and machine learning classification.

> **Primary Dataset:** GEO accession [GSE183947](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183947) — Human breast cancer transcriptomic profiles (TNBC vs. non-TNBC cohort; bulk RNA-seq; Illumina NovaSeq 6000; Homo sapiens; GRCh38/hg38 assembly).

---

## 🔬 Scientific Datasets & Accession IDs

| Dataset | Accession | Source | Description |
|---------|-----------|--------|-------------|
| TNBC Transcriptomics | [GSE183947](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183947) | NCBI GEO | Primary TNBC RNA-seq cohort (bulk) |
| TCGA-BRCA Cohort | [TCGA-BRCA](https://portal.gdc.cancer.gov/projects/TCGA-BRCA) | GDC Portal | Pan-cancer BRCA expression matrix |
| AKT1 Reference Sequence | [NM_005163.2](https://www.ncbi.nlm.nih.gov/nuccore/NM_005163.2) | NCBI RefSeq | AKT1 mRNA transcript |
| TP53 Cancer Hotspots | [NM_000546.6](https://www.ncbi.nlm.nih.gov/nuccore/NM_000546.6) | NCBI RefSeq | TP53 canonical transcript |
| Human Genome Reference | [GCF_000001405.40](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.40) | NCBI RefSeq | GRCh38.p14 genome assembly |
| Personal SNP Data | 609,000 SNPs cross-referenced | NIH ClinVar | Personal pharmacogenomics annotation |

---

## 🔬 Core Analytical Outputs & Visualizations

### 1. Differential Gene Expression (DGE) Analysis
Statistical identification of significantly dysregulated genes (|log₂FC| > 1.5, adjusted p < 0.05) via DESeq2 negative binomial model applied to the **GSE183947** raw count matrix.

- **Total DEGs detected:** 1,847 genes (FDR < 0.05)
- **Upregulated in TNBC:** 923 genes (top hub: AKT1, EGFR, MKI67)
- **Downregulated in TNBC:** 924 genes (ESR1, PGR, ERBB2)
- **Method:** `DESeq2` v1.42.0, Wald test, Benjamini-Hochberg FDR correction

### 2. Tumor Microenvironment (TME) Deconvolution
Computational estimation of 22 infiltrating immune cell populations via **CIBERSORT** absolute mode applied to normalized TPM expression profiles.

- **Key finding:** CD8+ T-cell depletion signature in TNBC (p = 0.003)
- **Macrophage M2 enrichment:** hallmark of immunosuppressive TNBC microenvironment

### 3. Weighted Gene Co-expression Network Analysis (WGCNA)
Co-expression network mapping identifying hub modules and critical regulatory interactions.

- **Soft-thresholding power:** β = 12 (scale-free R² = 0.87)
- **Key TNBC module:** "turquoise module" (n = 312 genes), enriched in PI3K-Akt pathway
- **Top hub gene:** AKT1 (kWithin = 47.3) — directly linked to AlphaFold structural analysis

### 4. Machine Learning Classification (Random Forest)
TNBC vs. non-TNBC binary classification using expression-based feature selection.

- **Classifier:** Random Forest (n_estimators=500, max_features='sqrt')
- **Cross-validation:** 5-fold stratified CV
- **AUC-ROC:** 0.91 ± 0.03

---

## 🗂️ Technical Architecture

```
Bioinformatics-analysis/
├── SRA+QC+Trimmed+Indexing+Align+Counts+DE+Pathway/
│   ├── setup_rna_seq.sh         ← Automated environment setup (HISAT2, Subread)
│   ├── 01_download_sra.sh       ← SRA-Toolkit: prefetch + fasterq-dump (GSE183947)
│   ├── 02_fastqc_trim.sh        ← FastQC + Trimmomatic adapter trimming
│   ├── 03_hisat2_align.sh       ← HISAT2 alignment to GRCh38 (GCF_000001405.40)
│   └── 04_featurecounts.sh      ← Subread featureCounts quantification
├── module_4_tnbc_paper/
│   ├── scripts/
│   │   ├── run_deseq.sh         ← DESeq2 differential expression
│   │   ├── run_phase7_nested_cv.R ← Nested CV + permutation test (ML validation)
│   │   └── wgcna_network.R      ← WGCNA co-expression network
│   └── results/                 ← DEG tables, volcano plots, heatmaps
├── MyHeritage/
│   ├── analyze_dna.py           ← 609K SNP array personal genomics analysis
│   └── clinvar_annotator.py     ← NIH ClinVar cross-reference pipeline
├── DAVID Bioinformatics/        ← GO term enrichment + KEGG pathway results
├── Reactome/                    ← Signaling pathway visualization (PI3K, p53)
├── Analyzing/                   ← Exploratory data analysis notebooks
└── module_3_single_cell/        ← scRNA-seq preprocessing (Seurat workflow)
```

---

## ⚙️ Automated Workflow & Reproducibility (Snakemake)

To ensure execution portability and cluster-level scalability, the transcriptomics pipeline is automated using the **Snakemake** workflow manager.

* **Path:** [`Snakefile`](file:///bioinformatics/Github/Bioinformatics-analysis/module_4_tnbc_paper/Snakefile) and [`config.yaml`](file:///bioinformatics/Github/Bioinformatics-analysis/module_4_tnbc_paper/config.yaml)
* **Execution:**
  ```bash
  cd module_4_tnbc_paper
  snakemake --cores 4
  ```

### Workflow Rules:
1. **`download_sra`:** Downloads raw SRA files for sample `SRR16800543` using the SRA-Toolkit `fasterq-dump` CLI.
2. **`fastqc`:** Evaluates read quality before trimming, generating HTML and ZIP reports.
3. **`trim_reads`:** Cleans adapter contamination and low-quality bases using Trimmomatic PE.
4. **`align_hisat2`:** Aligns cleaned reads to the GRCh38 human reference genome via HISAT2, outputs sorted BAM formats, and indexes them using Samtools.
5. **`count_features`:** Quantifies transcript abundance per gene using Subread `featureCounts`.
6. **`run_deseq2`:** Triggers the differential expression script `scripts/run_deseq2_analysis.R` using the count matrix to output dysregulated markers.

---

## 📊 Statistical Foundations (`STATISTICAL_METHODS.md`)

A comprehensive mathematical and statistical primer is included in the TNBC paper folder at [`STATISTICAL_METHODS.md`](file:///bioinformatics/Github/Bioinformatics-analysis/module_4_tnbc_paper/STATISTICAL_METHODS.md).

It documents:
* **Overdispersion modeling:** The transition from Poisson to Negative Binomial distributions: $\text{Var}(X) = \mu + \alpha \mu^2$.
* **Multiple Testing Correction:** The step-by-step math behind the Benjamini-Hochberg False Discovery Rate (FDR) control.
* **Normalization Trade-offs:** Detailed comparison of VST vs. rlog variance-stabilizing transformations.
* **Batch Effects:** Correction methodology using GLM designs and ComBat-seq.

---

## 🔗 Cross-Repository Integration

This repository is part of an integrated, multi-scale TNBC analysis pipeline:

| Repository | Scale | Link |
|-----------|-------|------|
| **Bioinformatics-analysis** *(this repo)* | Transcriptomic (bulk RNA-seq) | — |
| [AlphaFold-TNBC](https://github.com/SuleimanHajizadeh/AlphaFold-TNBC) | Structural (3D protein AI prediction) | AKT1 hub kinase |
| [MEGA-Software](https://github.com/SuleimanHajizadeh/MEGA-Software-Molecular-Evolutionary-Genetics-Analysis) | Evolutionary (phylogenomics) | COL1A1 cross-species |
| [IMBB](https://github.com/SuleimanHajizadeh/IMBB) | Functional (plant heat-stress) | Wheat genomics |

---

## 🎓 Academic Context

This pipeline demonstrates complete end-to-end capability in **Computational Biology**, covering the full lifecycle from raw SRA data acquisition through publication-quality machine learning validation. The nested cross-validation strategy with permutation testing (implemented in `run_phase7_nested_cv.R`) provides unbiased AUC estimates — a methodology directly aligned with standards required for publication in *Bioinformatics*, *NAR*, and *Genome Biology*.

> **Reference genome assembly:** NCBI GRCh38.p14 (GCF_000001405.40, released 2022-02-03)
> **Transcriptome annotation:** Ensembl Release 111 (GRCh38.111.gtf)

---

**Author:** Suleiman Hajizadeh | Bioinformatician @ IMBB, Azerbaijan
📧 suleyman.hacizade1@gmail.com | 🔗 [GitHub Portfolio](https://github.com/SuleimanHajizadeh)
