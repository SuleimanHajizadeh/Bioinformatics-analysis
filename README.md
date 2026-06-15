# 🧬 Transcriptomic Pipelines, Clinical Genomics & Machine Learning Profile

[![Bioinformatics](https://img.shields.io/badge/Bioinformatics-RNA--seq-blue.svg?style=flat-square)](https://www.ncbi.nlm.nih.gov/sra)
[![Method](https://img.shields.io/badge/Method-DESeq2_|_WGCNA_|_CIBERSORT-orange.svg?style=flat-square)](https://bioconductor.org/packages/DESeq2)
[![Tools](https://img.shields.io/badge/Tools-Bash_|_R_|_Python-green.svg?style=flat-square)](https://www.r-project.org/)
[![Dataset](https://img.shields.io/badge/Dataset-GEO%3A_GSE183947-red.svg?style=flat-square)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183947)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=flat-square)](https://opensource.org/licenses/MIT)

## 📌 Overview

This repository contains production-grade pipelines and analysis modules for **transcriptomic profiling, clinical oncology classification, and clinical genomics annotation**. The primary focus is the molecular characterization of **Triple-Negative Breast Cancer (TNBC)**, using a robust computational framework that spans from raw sequencing reads to peer-reviewed manuscript submission.

---

## 🗂️ Technical Architecture & Repository Layout

Designed in accordance with scientific best practices for computational workflows:

```
Bioinformatics-analysis/
├── requirements.txt            # Python dependencies
├── environment.yml             # Conda environment definition
├── .gitignore                  # Git untracked pattern file
├── LICENSE
├── README.md
│
├── workflows/                  # Reusable analysis pipelines
│   ├── bulk_rnaseq/            # Snakemake + R bulk RNA-Seq pipeline
│   │   ├── Snakefile           # Snakemake workflow rules (SRA to DEGs)
│   │   ├── config.yaml         # Configuration file for workflow rules
│   │   ├── scripts/            # Shell and R scripts (DESeq2, WGCNA, ML classification)
│   │   └── results/            # Volcano plots, co-expression networks, pathways
│   └── single_cell/            # Seurat single-cell transcriptomics workflow
│       ├── scripts/
│       └── results/
│
├── clinical_genomics/          # Variant annotation and personal genomics
│   └── myheritage_snp/         # NIH ClinVar annotation of 609K personal SNP array data
│       ├── data/               # Input CSVs (Suleiman and Huseyn cohorts)
│       ├── scripts/            # analyze_dna.py & clinvar_annotator.py
│       └── results/            # ClinVar match reports and CSV logs
│
├── exploratory_analysis/       # Quality control reports and practice modules
│   ├── fastqc_reports/         # Pre- and post-trimmed Read Quality reports
│   └── practice_modules/       # Exploratory bioinformatics scripts and coursework
│
├── manuscripts/                # Peer-reviewed journal submission materials
│   └── PLOS_ONE_SUBMISSION/    # Manuscript, cover letter, and participant checklists
│
└── docs/                       # Theoretical and statistical primers
    ├── STATISTICAL_METHODS.md  # Poisson/Negative Binomial math, FDR & VST normalization
    └── DNT_Gunu_Bioinformatika.pdf
```

---

## 🔬 Core Analytical Modules

### 1. Differential Gene Expression (DGE) Analysis
Statistical identification of significantly dysregulated genes (|log₂FC| > 1.5, adjusted p < 0.05) via the **DESeq2** negative binomial model.
*   **Dataset:** [GSE183947](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183947) (Illumina NovaSeq 6000 raw counts).
*   **Total DEGs detected:** 1,847 genes (FDR < 0.05)
*   **Upregulated in TNBC:** 923 genes (top hubs: *AKT1*, *EGFR*, *MKI67*)
*   **Downregulated in TNBC:** 924 genes (*ESR1*, *PGR*, *ERBB2*)

### 2. Weighted Gene Co-expression Network Analysis (WGCNA)
System-level mapping of breast cancer subtypes to identify hub genes and regulatory modules.
*   **Soft-thresholding power:** β = 12 (scale-free topology model fit R² = 0.87)
*   **Key TNBC module:** "turquoise module" (n = 312 genes), highly enriched in the PI3K-Akt signaling pathway.
*   **Top hub gene:** *AKT1* (kWithin = 47.3) — directly linked to downstream AlphaFold structural modeling.

### 3. Machine Learning Classification (Random Forest)
Robust TNBC vs. non-TNBC binary classification using transcriptomic expression-based feature selection.
*   **Validation Strategy:** Nested 5-Fold Stratified Cross-Validation with permutation testing.
*   **Performance:** AUC-ROC = 0.91 ± 0.03.

### 4. Clinical Variant Annotation (ClinVar Parser)
A high-throughput Python script mapping personal genomic SNP microarray data (609K SNPs) to the NIH ClinVar database to identify pathogenic and drug-response alleles.

---

## ⚙️ Automated Workflow & Reproducibility (Snakemake)

The bulk RNA-Seq pipeline is automated using the **Snakemake** workflow manager to ensure scalability and reproducibility:

```bash
# Navigate to the bulk RNA-seq workspace
cd workflows/bulk_rnaseq

# Execute the workflow using 4 cores
snakemake --cores 4
```

### Workflow Pipeline Steps:
1. **`download_sra`:** Downloads raw SRA sequence files via SRA-Toolkit `fasterq-dump`.
2. **`fastqc`:** Evaluates read quality before trimming, generating HTML and ZIP reports.
3. **`trim_reads`:** Cleans adapter contamination and low-quality bases using Trimmomatic PE.
4. **`align_hisat2`:** Aligns cleaned reads to the GRCh38 human reference genome via HISAT2.
5. **`count_features`:** Quantifies transcript abundance per gene using Subread `featureCounts`.
6. **`run_deseq2`:** Performs differential expression analysis, outputting dysregulated markers.

---

## 🔗 Cross-Repository Integration

This portfolio is integrated with other specialized scientific repositories:

| Scale | Repository | Focus |
|-------|------------|-------|
| **Transcriptomic (Bulk/scRNA-seq)** | [Bioinformatics-analysis](https://github.com/SuleimanHajizadeh/Bioinformatics-analysis) *(this repo)* | Molecular profiling & ML classification |
| **Structural Biology (3D Protein AI)** | [computational-structural-biology](https://github.com/SuleimanHajizadeh/computational-structural-biology) | AKT1/STN7 kinase modeling & molecular docking |
| **Phylogenomics & Evolution** | [MEGA-Software-Genetics](https://github.com/SuleimanHajizadeh/MEGA-Software-Molecular-Evolutionary-Genetics-Analysis) | COL1A1 multi-species sequence alignment |

---

## 🎓 Academic Context

This repository is optimized for academic peer review. The statistical foundation document [`docs/STATISTICAL_METHODS.md`](./docs/STATISTICAL_METHODS.md) contains the detailed equations and rationales for:
*   Poisson vs. Negative Binomial distribution modeling: $\text{Var}(X) = \mu + \alpha \mu^2$
*   Benjamini-Hochberg False Discovery Rate (FDR) control math
*   VST vs. rlog variance-stabilizing transformation comparisons

The manuscripts folder [`manuscripts/PLOS_ONE_SUBMISSION/`](./manuscripts/PLOS_ONE_SUBMISSION/) contains the revised manuscript draft and supplementary documentation prepared for submission to *PLOS ONE*.

---

**Author:** Suleiman Hajizadeh | Bioinformatician @ IMBB, Azerbaijan  
📧 suleyman.hacizade1@gmail.com | 🔗 [GitHub Portfolio](https://github.com/SuleimanHajizadeh)
