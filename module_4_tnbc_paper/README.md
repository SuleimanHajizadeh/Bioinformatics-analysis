# 🧬 Transcriptomic Profiling and Machine Learning Analysis of Triple-Negative Breast Cancer (TNBC)

[![Status](https://img.shields.io/badge/Status-Under_Review-blue.svg)]()
[![Method](https://img.shields.io/badge/Method-RNA--Seq-orange.svg)]()
[![Tools](https://img.shields.io/badge/Tools-R_|_Python-green.svg)]()

## 📌 Project Overview
This repository contains the complete computational biology pipeline and analytical scripts used for the research manuscript detailing the transcriptomic profiling of Triple-Negative Breast Cancer (TNBC) subtypes.

The study leverages next-generation sequencing (NGS) data to perform differential expression analysis, identify critical prognostic biomarkers, and construct robust machine learning classification models. The manuscript is currently under review at **PLOS ONE**.

## 📊 Methodology & Pipeline
The computational workflow is entirely reproducible and includes:
1. **Raw Data Processing:** Quality control (FastQC) and alignment (HISAT2) of raw FASTQ reads.
2. **Differential Expression Analysis:** Utilizing `DESeq2` in R to identify statistically significant transcriptomic alterations between TNBC and non-TNBC control tissues.
3. **Dimensionality Reduction:** Principal Component Analysis (PCA) and hierarchical clustering (Heatmaps) to explore transcriptomic variance.
4. **Functional Enrichment:** Gene Ontology (GO) and KEGG pathway analyses to elucidate the biological systems driving TNBC progression.
5. **Machine Learning:** Implementing Random Forest classifiers to evaluate the predictive power and Variable Importance of the top differentially expressed genes.

## 📂 Repository Structure
* `/scripts` - R and Python scripts for data processing and statistical modeling.
* `/data` - Contains processed count matrices and normalized datasets.
* `/results` - Generated analytical outputs.
* `/PLOS_ONE_SUBMISSION` - Formatted manuscript files, tables, and figures prepared for academic publishing.

## 🔬 Key Visualizations Included
The repository contains high-resolution analytical plots:
- **Differential Expression:** `Volcano_plot.png` & `Heatmap_plot.png`
- **Dimensionality Reduction:** `PCA_plot.png`
- **Predictive Modeling:** `ML_ROC_Curve.png` & `ML_Variable_Importance.png`
- **Systems Biology:** `Systems_Biology_Network.png` & `TME_Deconvolution_Barplot.png`

## 🎓 Academic Context
This research project was designed, executed, and documented by Süleyman Hajızadə as a demonstration of independent capability in **Computational Biology, Transcriptomics, and Machine Learning**. It serves as a foundational component for advanced academic applications (e.g., MPhil in Computational Biology).
