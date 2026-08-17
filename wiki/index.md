# Wiki Catalog & Index

Welcome to the **LLM Wiki** — a persistent, compounding knowledge base for transcriptomic analysis, clinical genomics, and computational biology.

---

## 📊 Overview of Wiki Structure

| Section | Description | Item Count |
| :--- | :--- | :--- |
| **[[#Concepts\|Concepts]]** | Theoretical foundations, algorithms, statistical & biophysical principles | 10 |
| **[[#Entities\|Entities]]** | Software tools, packages, databases, genes, drugs, algorithms | 12 |
| **[[#Syntheses\|Syntheses]]** | Cross-domain summaries, multi-source frameworks, theses | 3 |
| **[[#Summaries\|Summaries]]** | Source-specific summaries linking back to `raw/` and codebase | 6 |

---

## 🧠 Concepts

- [[negative-binomial-model]]: Mathematical formulation of overdispersed count data in RNA-Seq.
- [[multiple-testing-fdr]]: Benjamini-Hochberg FDR control vs. Bonferroni FWER for large-scale omics.
- [[variance-stabilizing-transformation]]: Homoscedastic data transformations (VST vs. rlog) for machine learning & PCA.
- [[wgcna-coexpression-networks]]: Weighted gene co-expression network analysis, scale-free topology, and hub gene identification.
- [[nested-cross-validation-ml]]: Preventing feature selection data leakage in high-dimensional omics ($p \gg n$) using nested LOOCV & permutation tests.
- [[single-cell-clustering-umap]]: 10x Genomics scRNA-seq workflow (QC, VST variable genes, PCA, Louvain graph clustering, UMAP).
- [[immune-deconvolution-tme]]: In silico digital pathology deconvolution of tumor microenvironment cell populations.
- [[variant-pathogenicity-annotation]]: Mapping personal genomic SNP arrays against NIH ClinVar for clinical significance and pharmacogenomics.
- [[alphafold-structural-modeling]]: 3D macromolecular prediction, pLDDT local accuracy, PAE inter-domain metrics, and Ramachandran validation.
- [[molecular-docking-virtual-screening]]: Physics of binding affinity ($\Delta G_{\text{bind}}$), allosteric vs. orthosteric inhibition, and virtual screening.

---

## 🏷️ Entities

- [[deseq2]]: Bioconductor R package for differential gene expression using empirical Bayes dispersion shrinkage.
- [[triple-negative-breast-cancer]]: Molecular characteristics, receptor profiles (*ESR1*, *PGR*, *ERBB2*), and oncogenic drivers.
- [[snakemake]]: Python-based workflow management engine for reproducible computational pipelines.
- [[hisat2]]: Fast and sensitive spliced alignment program based on hierarchical graph FM indexing.
- [[featurecounts]]: Ultrafast read summarization utility for genomic feature quantification.
- [[seurat]]: Leading R toolkit for single-cell genomics, QC, multimodal integration, and clustering.
- [[clinvar]]: NCBI public archive of human genetic variants, clinical assertions, and evidence tiers.
- [[random-forest-classifier]]: Ensemble tree classifier with permutation testing for robust transcriptomic biomarker discovery.
- [[akt1-kinase]]: RAC-alpha serine/threonine-protein kinase (UniProt P31749) master oncogenic hub in TNBC.
- [[alphafold2]]: DeepMind neural structural prediction system.
- [[autodock-vina]]: Fast, accurate molecular docking and virtual screening engine.
- [[mk-2206]]: Potent, highly selective allosteric inhibitor of AKT kinases ($\Delta G = -9.4\text{ kcal/mol}$).

---

## 🔬 Syntheses

- [[transcriptomic-analysis-framework]]: Unified end-to-end framework linking raw counts, count modeling, WGCNA, and ML classification.
- [[multi-omics-analytical-ecosystem]]: Cross-layer synthesis connecting bulk RNA-Seq, single-cell resolution transcriptomics, and personal clinical genomics.
- [[transcriptome-to-structure-drug-discovery]]: Multi-scale framework connecting whole-genome RNA-Seq co-expression networks to AlphaFold structural biology and virtual docking.

---

## 📚 Summaries

- [[statistical-methods-primer]]: Summary of mathematical and statistical methodology in RNA-Seq differential expression (`docs/STATISTICAL_METHODS.md`).
- [[bulk-rnaseq-snakemake-pipeline]]: Ingestion of the Snakemake workflow orchestration and rules (`workflows/bulk_rnaseq/Snakefile`).
- [[tnbc-analytics-suite]]: Ingestion of Phases 5–9 downstream R analytics (DESeq2, clusterProfiler, Nested CV Random Forest, iGraph/WGCNA, CIBERSORTx).
- [[single-cell-seurat-workflow]]: Ingestion of the Seurat single-cell transcriptomics pipeline (`workflows/single_cell/scripts/`).
- [[clinical-genomics-pipeline]]: Ingestion of personal SNP array analysis and ClinVar variant annotation scripts (`clinical_genomics/myheritage_snp/scripts/`).
- [[tnbc-akt1-alphafold-docking-manuscript]]: Ingestion of the PLOS ONE manuscript connecting TNBC RNA-Seq, WGCNA, AlphaFold2 AKT1 modeling, and AutoDock Vina inhibitor screening (`manuscripts/manuscript_draft.md`).

---

## 📜 System Log

Track recent changes and ingestion records in [[log|wiki/log.md]].
