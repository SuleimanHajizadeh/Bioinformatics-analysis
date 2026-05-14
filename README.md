# 🧬 Comprehensive Transcriptomic Pipelines & Data Analysis

[![Bioinformatics](https://img.shields.io/badge/Bioinformatics-RNA--seq-blue.svg?style=flat-square)]()
[![Method](https://img.shields.io/badge/Method-DESeq2_|_ML_|_Network_Analysis-orange.svg?style=flat-square)]()
[![Tools](https://img.shields.io/badge/Tools-Bash_|_R_|_Python-green.svg?style=flat-square)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=flat-square)](https://opensource.org/licenses/MIT)

## 📌 Overview

This repository integrates high-throughput sequencing pipelines and specialized transcriptomic characterization workflows. The primary focus is the molecular profiling of tumor microenvironments and clinical genomics datasets, utilizing differential expression analysis, tumor microenvironment (TME) deconvolution, and systems biology network mapping.

---

## 🔬 Core Analytical Outputs & Visualizations

Our transcriptomic analysis pipelines yielded significant data-driven insights. Below are key visual outputs generated from raw sequencing datasets.

### 1. Differential Expression & Immunological Profiling
Statistical identification of significantly dysregulated genes through bulk RNA-seq analysis.

<p align="center">
  <img src="module_4_tnbc_paper/Volcano_plot.png" alt="Volcano Plot" width="600"/>
</p>

### 2. Tumor Microenvironment (TME) Deconvolution
Computational estimation of infiltrating immune cell populations from bulk transcriptomic profiles.

<p align="center">
  <img src="module_4_tnbc_paper/TME_Deconvolution_Barplot.png" alt="TME Deconvolution" width="600"/>
</p>

### 3. Systems Biology Network
Co-expression network mapping highlighting hub genes and critical regulatory interactions in the generated datasets.

<p align="center">
  <img src="module_4_tnbc_paper/Systems_Biology_Network.png" alt="Systems Biology Network" width="600"/>
</p>

---

## 🗂️ Technical Inventory and Architecture

### 1. Advanced Transcriptomic Profiling ([`module_4_tnbc_paper/`](./module_4_tnbc_paper/))
Contains downstream R scripts, results matrices, and high-resolution graphical outputs representing complex molecular profiles.

### 2. End-to-End RNA-Seq Pipeline ([`SRA+QC+.../`](./SRA+QC+Trimmed+Indexing+Align+Counts+DE+Pathway/))
A systematic, shell-driven environment for raw high-throughput data processing.
- **`setup_rna_seq.sh`**: Automated environment and directory configuration.
- Includes pre-processing modules: SRA downloading, FastQC validation, read trimming (Trimmomatic), and alignment (HISAT2).

### 3. Personal Clinical Genomics ([`MyHeritage/`](./MyHeritage/))
Python-based algorithms (`analyze_dna.py`, `clinvar_annotator.py`) designed to cross-reference 609K personal SNPs against the NIH ClinVar database for pharmacogenomic annotation.

### 4. Pathway Topology Mapping ([`Reactome/`](./Reactome/))
Spatial visualization of molecular signaling pathways and functional biological systems.

---

## ⚙️ Pipeline Implementation & Reproducibility

### Environment and Dependencies
The project utilizes a dedicated Conda environment on an **AlmaLinux 9** server.
- **Mapping:** HISAT2 (GRCh38 Indexing)
- **Quantification:** Subread (`featureCounts`)
- **Differential Expression:** R / `DESeq2`
- **Machine Learning:** `scikit-learn` (PCA, Random Forest)

### Execution Trace
```bash
# 1. Initialize repository structure and raw data handling
bash SRA+QC+Trimmed+Indexing+Align+Counts+DE+Pathway/setup_rna_seq.sh

# 2. Run Downstream Analysis
bash module_4_tnbc_paper/scripts/run_deseq.sh
```

---

## 🎓 Academic Context
This repository acts as a comprehensive technical portfolio demonstrating end-to-end capability in **Computational Biology**. It covers the entire lifecycle of high-throughput data analysis: from raw sequence quality control to publication-ready statistical and machine learning insights.

---

**Author:** Suleiman Hajizadeh | Bioinformatician @ IMBB, Azerbaijan
📧 suleyman.hacizade1@gmail.com
