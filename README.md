# Multi-Omics Framework: Transcriptomic Pipelines and Breast Carcinoma Characterization

[![Bioinformatics](https://img.shields.io/badge/Bioinformatics-Pipeline-blue.svg)](https://github.com/SuleimanHajizadeh/bioinformatics-analysis)
[![Status](https://img.shields.io/badge/Status-Manuscript--Ready-green.svg)](https://github.com/SuleimanHajizadeh/bioinformatics-analysis)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository integrates high-throughput sequencing pipelines and specialized transcriptomic characterization workflows. The primary focus is the molecular profiling of Triple-Negative Breast Cancer (TNBC), utilizing differential expression analysis (DESeq2), tumor microenvironment (TME) deconvolution, and pathway topology mapping across bulk and single-cell datasets.

---

## Technical Inventory and Project Architecture

The analytical assets are categorized by their role in the transcriptomic z-chain:

### 1. High-Resolution TNBC Characterization ([Module 4](./module_4_tnbc_paper/))
Contains primary research manuscript drafts and publication-ready metadata.
- **Key Signatures:** Analysis of **CXCL13**, **PDCD1**, and **CTLA4** immune-hot phenotypes.
- **Regulatory Networks:** Characterization of PI3K-Akt signaling hyperactivation in malignant biopsies.
- **Data Source:** Processed results from **GSE142731** and **GSE58135**.

### 2. QC Benchmarking and Technical Validation ([Analyzing/](./Analyzing/))
Technical validation reports for dataset integrity.
- **Samples:** FastQC metrics for **HBR_Rep1** and **HBR_Rep2** (HRCC-ordered ERCC-Mix2 transcripts).
- **Format:** Comprehensive HTML and ZIP archives for sequencing quality assessment (R1/R2 reads).

### 3. Functional Annotation and Target Reporting ([DAVID Bioinformatics/](./DAVID%20Bioinformatics/))
Detailed functional and biological process reports (Dated 2026-03-25) for core targets:
- **Malignant Markers:** `AKT1`, `AURKA`, `BIRC5`.
- **References:** `ACTB` (Actin Beta) housekeeping gene reports in CSV, PDF, and XLSX formats.

### 4. Pathway Topology and Visualization ([Reactome/](./Reactome/))
Spatial mapping of molecular signaling pathways.
- **Assets:** Reacfoam visualizations and pathway screenshots (March 28, 2026).
- **Pathways:** Integrated analysis of signaling cascades via the Reactome Database.

### 5. Automated Pipeline Archives ([SRA+QC+...](./SRA+QC+Trimmed+Indexing+Align+Counts+DE+Pathway/))
A systematic, shell-driven environment for end-to-end RNA-Seq processing.
- **Core Script:** `setup_rna_seq.sh` - Automated environment and directory configuration.
- **Modular Layout:** Archive of `module_1`, `module_2`, and `module_3_single_cell_small` benchmarking projects.

---

## Pipeline Implementation

### Environment and Dependencies
The project utilizes a dedicated **Conda** environment on an **AlmaLinux 9** server.
- **Mapping:** HISAT2 (GRCh38 Indexing).
- **Quantification:** Subread (featureCounts).
- **Normalization/DE:** R/DESeq2.

### Execution Trace
To initialize the research environment and run the core quantification:
```bash
# Initialize repository structure
bash SRA+QC+Trimmed+Indexing+Align+Counts+DE+Pathway/setup_rna_seq.sh

# Run Differential Expression (Module 4)
bash module_4_tnbc_paper/scripts/run_deseq.sh
```

---

## Contact

**Suleiman Hajizadeh**  
Lead Investigator @ Western Caspian University, Azerbaijan  
- **Email:** suleyman.hacizade1@gmail.com  
- **Expertise:** Omics Integration, Precision Oncology, Pathogen Genomics

---
> [!NOTE]
> All analytical results are accompanied by their corresponding source metrics (FastQC/DAVID). For access to high-volume RDS data (>100MB), refer to the provided remote synchronization protocols.
