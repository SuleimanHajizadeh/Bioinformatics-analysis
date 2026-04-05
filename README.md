# Bioinformatics-Analysis: Decoding the Molecular Blueprint of Triple-Negative Breast Cancer (TNBC)

[![Bioinformatics](https://img.shields.io/badge/Bioinformatics-Pipeline-blue.svg)](https://github.com/SuleimanHajizadeh/bioinformatics-analysis)
[![Status](https://img.shields.io/badge/Status-Manuscript--Ready-green.svg)](https://github.com/SuleimanHajizadeh/bioinformatics-analysis)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 🔬 Overview

This repository hosts a comprehensive multi-omics research portfolio and advanced bioinformatic pipelines dedicated to decoding the molecular drivers of Triple-Negative Breast Cancer (TNBC). By leveraging **Machine Learning**, **Spatial Deconvolution**, and **Systems Biology**, this project identifies critical transcriptomic signatures that define the tumor microenvironment (TME) and guide precision immunotherapy strategies.

### 🌟 Key Research Highlights (Module 4)

*   **Immune Landscape Priming:** Identification of an "Immune-Hot" phenotype defined by the co-expression of the chemokine **CXCL13** and checkpoint receptors **PDCD1 (PD-1)** and **CTLA4**.
*   **Pathway Hyperactivation:** Discovery of the **PI3K-Akt signaling cascade** as a generalized master regulatory network driving malignant proliferation in TNBC.
*   **TME Deconvolution:** Mathematical isolation of **M2 Macrophage dominance** against the lack of cytotoxic CD8+ T-cell infiltration, effectively stratifying patients for checkpoint blockade therapy.

---

## 📁 Repository Structure

The project is structured into modular bioinformatics pipelines and specialized research modules:

### [Module 4: TNBC Research Project (Final Paper)](./module_4_tnbc_paper/)
The core of this portfolio, containing the manuscript drafts, final figures, and summarized results for submission to high-impact journals (AACR, Elsevier, Springer Nature).
- **Manuscripts:** Specialized versions for *Cancer Research*, *Computers in Biology and Medicine*, and *Scientific Reports*.
- **Figures:** PCA plots, Volcano plots, Heatmaps, and KEGG/GO pathway enrichments.
- **Data:** Cleaned counts and phenotypic metadata (GSE58135/SRP042620).

### [Module 2: Bulk RNA-Seq Pipeline (GSE58135)](./module_2/)
A high-performance command-line pipeline for processing raw RNA-Seq biopsies.
- **Workflow:** FastQC -> Trimmomatic -> HISAT2 (GRCh38) -> featureCounts -> DESeq2.
- **Automation:** Bash drivers for batch processing large sets of SRA accessions.

### [Module 3: Single-Cell RNA-Seq (PBMC Analysis)](./module_3_single_cell/)
Advanced single-cell transcriptomics using the R Seurat framework.
- **Analysis:** Quality Control (QC), t-SNE/UMAP clustering, and differential marker identification for cell type annotation.

---

## 🛠️ Technical Stack & Workflow

### Bioinformatics Environment
The analysis is conducted within a dedicated Conda environment (`rna_seq_env`) on a Linux server (AlmaLinux 9).

| Stage | Tools Used |
| :--- | :--- |
| **Quality Control** | FastQC (v0.12.1), MultiQC |
| **Read Trimming** | Trimmomatic (v0.39) |
| **Mapping/Alignment** | HISAT2 (v2.2.1) using GRCh38 |
| **Quantification** | Subread/featureCounts (v2.0.6) |
| **Differential Expression** | R/DESeq2 (v1.38.3) |
| **Functional Enrichment** | R/clusterProfiler (v4.6.2) |
| **Visualization** | ggplot2, pheatmap, cowplot |

---

## 📊 Key Visualizations

The research includes high-impact visualizations ready for publication:
- **PCA Analysis:** Confirming distinct transcriptomic clusters between normal and malignant biopsies.
- **Volcano Plots:** Identifying 450 robust DEGs (|Log2FC| > 2, FDR < 0.05).
- **Immune Heatmaps:** Decoding the spatial deconvolution of checkpoint signatures.
- **KEGG Dotplots:** Visualizing the hyperactivation of PI3K-Akt and Cell Cycle pathways.

---

## 🚀 Setup & Execution

### 1. Clone the Portfolio
```bash
git clone https://github.com/SuleimanHajizadeh/bioinformatics-analysis.git
cd bioinformatics-analysis
```

### 2. Environment Configuration
The core dependencies can be found in `environment.yml`:
```bash
conda env create -f module_4_tnbc_paper/environment.yml
conda activate rna_seq_env
```

### 3. Run Analysis
To execute the DESeq2 pipeline:
```bash
bash module_4_tnbc_paper/scripts/run_deseq.sh
```

---

## 📧 Contact & Collaboration

**Suleiman Hajizadeh**  
Lead Investigator  
*Western Caspian University, Azerbaijan*  
- **Email:** suleyman.hacizade1@gmail.com  
- **Research Topic:** Precision Oncology and Systems Biology

---
> [!NOTE]  
> All analysis scripts are designed for reproducibility. For individual file access to those >100MB (RDS data), please contact the author directly or use the provided RSYNC/SCP commands.
