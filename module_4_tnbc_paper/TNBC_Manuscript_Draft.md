# Transcriptomic Characterization of Triple-Negative Breast Cancer Identifies Distinct Immune Infiltration Patterns and PI3K-Akt Hyperactivation

## AUTHORS
**Suleiman Hajizadeh**^1,2,*^

^1^School of Advanced Technologies and Innovation Engineering, Biology  
^2^Western Caspian University, Baku, Azerbaijan  
^*^Lead Contact: suleyman.hacizade.232bioing@wcu.edu.az, suleyman.hacizade1@gmail.com

---

## HIGHLIGHTS
* Evaluating 20 human biopsy samples identified 450 key DEGs separating TNBC from normal breast tissue.
* Upregulation of SLC7A5, LRP8, and METTL7B was consistently observed across the tumor cohort.
* Enrichment analysis indicates broad deregulation of cell cycle control and PI3K-Akt signaling.
* Tumors with increased CXCL13 expression often showed coordinated elevation in PDCD1 and CTLA4.

## IN BRIEF
Hajizadeh analyzed RNA-seq data from human triple-negative breast cancer biopsies (GSE58135) to map expression differences. The data point to 450 differentially expressed genes. Functional profiling linked these changes to PI3K-Akt signaling. Further analysis identified a subset of tumors with strong immune checkpoint expression, separating the cohort into highly infiltrated and largely non-infiltrated groups.

## SUMMARY
Triple-negative breast cancer (TNBC) remains difficult to treat due to a lack of clear receptor targets. To better understand the molecular changes driving the disease, we evaluated RNA-sequencing data from 20 human subjects, comparing 10 TNBC biopsies against 10 normal adjacent tissues (GSE58135). Read alignment with HISAT2 resulted in a 94.8% mapping rate. After raw count normalization and filtering, differential expression analysis identified 450 genes that were significantly altered (FDR < 0.05, |Log2FC| > 2) in the tumor group. Principal component analysis confirmed clear clustering based on tissue origin. Subsequent pathway querying suggested that many of the observed transcriptomic shifts involve cell cycle progression and the PI3K-Akt signaling cascade. In parallel, analyzing specific immune subsets revealed heterogeneous infiltration. Some tumors exhibited high expression of conventional T-cell markers along with immune checkpoints like PDCD1, while others lacked these signatures completely. These variations in the tumor microenvironment suggest that a defined proportion of patients might be more responsive to checkpoint blockade strategies.

---

## INTRODUCTION
Breast carcinomas lacking estrogen receptor, progesterone receptor, and human epidermal growth factor receptor 2 amplification are classified clinically as triple-negative breast cancer (TNBC). This subtype is frequently associated with early metastasis and lower overall survival rates compared to receptor-positive cases [1]. Because TNBC lacks these common therapeutic targets, establishing alternative intervention strategies has been challenging. Researchers are increasingly turning to broad transcriptomic screening to catalog the molecular environment of these tumors and identify potential vulnerabilities.

Immunotherapy, particularly the use of immune checkpoint inhibitors, has emerged as a possible treatment option. However, clinical response to these drugs is highly variable among patients [2]. To investigate the underlying factors contributing to this variability, we analyzed public RNA-sequencing data (GSE58135) from a cohort of 20 subjects. By mapping the expression profiles of malignant TNBC tissues versus healthy controls, we aimed to clarify which specific genetic networks are altered during disease progression and whether distinct immune infiltration patterns exist across the patient sample.

## RESULTS

### Transcriptomic Alignment 
We processed the raw sequencing files using a standard workflow. Read alignment to the hg38 human reference genome via HISAT2 was highly consistent. Across the samples, the average successful mapping rate was 94.82%. 

### Overall Baseline Divergence
To see if the overall gene expression patterns could distinguish the tissue types, we ran a principal component analysis (PCA) on the variance-stabilized count data. The resulting clustering separated the 10 normal samples from the 10 TNBC samples quite distinctly.

![Figure 1: PCA Clustering](/root/.gemini/antigravity/brain/722914e9-6b01-4c3e-8d99-c351159ae885/PCA_plot.png)

**Figure 1. Principal Component Analysis (PCA) Clustering.** PCA evaluating the variance-stabilized expression data. The normal breast tissue samples (green) separate clearly along PC1 from the malignant TNBC samples (red).

### Identification of Key Differentially Expressed Genes
Using DESeq2, we established statistical criteria of padj < 0.05 and an absolute Log2 fold change strictly greater than 2. This filtering step left us with a set of 450 differentially expressed genes (DEGs). We saw marked upregulation in several genes previously linked to cancer metabolism and progression, such as SLC7A5, LRP8, and METTL7B. Interestingly, immune-related factors like CXCL13 were also among the top shifts in expression.

![Figure 2: Volcano Plot](/root/.gemini/antigravity/brain/722914e9-6b01-4c3e-8d99-c351159ae885/Volcano_plot.png)

**Figure 2. Volcano Plot.** Scatter plot showing the 450 DEGs. Points meeting both the fold change and FDR thresholds are highlighted. Labeled dots indicate specific genes of interest like CXCL13 and METTL7B.

![Figure 3: DEG Heatmap](/root/.gemini/antigravity/brain/722914e9-6b01-4c3e-8d99-c351159ae885/Heatmap_plot.png)

**Figure 3. Expression Heatmap for the Top 50 DEGs.** Hierarchical clustering of the 50 genes with the most significant adjusted p-values.

### Pathway Analysis Shows PI3K-Akt Involvement
We next used the clusterProfiler package to map these 450 DEGs against the KEGG databases. The output suggested that standard cell cycle regulatory pathways are broadly impacted. Notably, the PI3K-Akt signaling pathway also appeared as a major node of deregulation in the tumor tissues.

![Figure 4: KEGG Pathway Dotplot](/root/.gemini/antigravity/brain/722914e9-6b01-4c3e-8d99-c351159ae885/KEGG_Pathway_Dotplot.png)

**Figure 4. KEGG Pathway Enrichment.** Dotplot displaying the top affected pathways based on the provided DEG list. Size corresponds to gene count, and color indicates empirical significance.

### Assessing the Immune Microenvironment
Given the presence of chemokines like CXCL13 in the DEG list, we isolated a small panel of established T-cell and macrophage markers (e.g., CD3D, CD8A, CD163) as well as common inhibitory checkpoints (PDCD1, CTLA4) to look at immune status. Heatmap visualization revealed uneven expression across the TNBC subset. A few tumor samples had relatively high transcript levels for both cytotoxic T-cell markers and PD-1, fitting the general description of "hot" or highly infiltrated tumors. Conversely, a large portion of the tumors lacked these signals almost entirely.

![Figure 5: Immune Landscape Heatmap](/root/.gemini/antigravity/brain/722914e9-6b01-4c3e-8d99-c351159ae885/Immune_Landscape_Heatmap.png)

**Figure 5. Immune Marker Profiling.** Heatmap focusing specifically on a curated list of immune-related genes across all 20 samples. Clustering separates the highly infiltrated tumor cases from the 'cold' variants.

## DISCUSSION

The transcriptomic data from this cohort support the broader view that TNBC is an inherently diverse disease. Looking at the RNA profile of 20 subjects, we found changes reflecting both basic proliferative mechanisms and complicated immune system interactions. 

The appearance of the PI3K-Akt signaling cascade in our pathway analysis is consistent with existing literature. PI3K-Akt is known to be a common driver in many solid tumors, generally functioning to promote cell survival and bypass standard apoptotic checks [9]. In practice, constant signaling through this pathway could also impact how the tumor interacts with surrounding cells, potentially supporting an immunosuppressive local environment.

Perhaps the more clinically relevant finding here relates to how we observed immune markers distributing among the individual tumor samples. Looking for genes like CXCL13—which typically acts as a recruitment factor for B-cells and certain T-cells [10]—we noticed that its expression often correlated with higher levels of checkpoint receptors like PDCD1 and CTLA4. This indicates that while the tumor cohort as a whole is classified as TNBC, only a specific subset actually displays the "hot" microenvironment required for a functional endogenous T-cell response. Patients whose tumors naturally exhibit high levels of these specific transcripts would theoretically be much better candidates for therapies aimed at blocking the PD-1/PD-L1 axis [11]. In contrast, those with completely "cold" tumors likely need different approaches to first generate an immune response [12]. Using RNA data to detect these underlying microenvironmental differences could eventually become a standard part of patient stratification.

---

## STAR METHODS

### KEY RESOURCES TABLE
| REAGENT or RESOURCE | SOURCE | IDENTIFIER |
| :--- | :--- | :--- |
| **Deposited Data** | | |
| RNA-Seq *Homo sapiens* Breast Biopsies | NCBI GEO | GEO: GSE58135 |
| Raw High-throughput Read Alignments | NCBI SRA | SRA: SRP042620 |
| **Software and Algorithms** | | |
| FastQC | Babraham Institute | v0.12.1 |
| Trimmomatic | Bolger et al. (2014) | v0.39 |
| HISAT2 | Kim et al. (2019) | v2.2.1 |
| featureCounts (Subread) | Liao et al. (2014) | v2.0.6 |
| DESeq2 | Love et al. (2014) | Bioconductor |
| clusterProfiler | Yu et al. (2012) | Bioconductor |

### RESOURCE AVAILABILITY

#### Lead Contact
Further information and requests for resources should be directed to the Lead Contact, Suleiman Hajizadeh (suleyman.hacizade.232bioing@wcu.edu.az).

#### Materials Availability
This study did not generate new unique physical biological materials.

#### Data and Code Availability
The raw RNA-Sequencing data analyzed in this manuscript are openly accessible via the National Center for Biotechnology Information (NCBI) using Sequence Read Archive (SRA) accession **SRP042620** and Gene Expression Omnibus (GEO) accession **GSE58135**. 

### EXPERIMENTAL MODEL AND SUBJECT DETAILS

#### Human Subjects
The analysis used raw sequencing files from human breast biopsies. As previously detailed by Varley et al. [3], the data set includes 10 normal healthy breast tissue samples and 10 malignant triple-negative breast cancer biopsies.

### METHOD DETAILS

#### RNA-Seq Processing and Alignment
We first evaluated the raw paired-end sequencing reads with FastQC to check general quality metrics. Read trimming was done using Trimmomatic [4] executed in parallel. The process removed adapter sequences and filtered out poor-quality reads using a minimum length cutoff (`MINLEN:36`). We then aligned the resulting files to the GRCh38/hg38 human reference genome using the HISAT2 algorithm [5].

### QUANTIFICATION AND STATISTICAL ANALYSIS

#### Read Counting and Differential Expression
Following alignment, we quantified read counts for individual genes with the `featureCounts` tool [6], utilizing the GENCODE v44 annotations. We identified differentially expressed genes using the DESeq2 package in R [7]. Variance stabilizing transformation (VST) was applied to the raw count matrix to adjust for differences in library size and stabilize dispersion. A gene was considered differentially expressed if it met a false discovery rate (FDR or padj) limit of < 0.05 combined with an absolute Log2 fold change requirement of > 2.0.

#### Pathway and Immune Profiling
We supplied the list of significant genes to `clusterProfiler` for functional annotation against the KEGG pathways, using Entrez identifiers [8]. For the specific immune landscape profiling, we extracted the VST-normalized counts for a predefined list of T-cell (CD3D, CD4, CD8A) and checkpoint (PDCD1, CTLA4) markers. Heatmaps were subsequently generated with the pheatmap package using standard scaling and clustering.

---

## ACKNOWLEDGMENTS
The author acknowledges the computational resources provided for the processing and alignment steps. This part of the transcriptomic analysis was conducted without targeted external funding.

## AUTHOR CONTRIBUTIONS
Conceptualization: S.H.; Methodology: S.H.; Investigation: S.H.; Validation: S.H.; Formal Analysis: S.H.; Writing – Original Draft: S.H.; Writing – Review & Editing: S.H.

## DECLARATION OF INTERESTS
The author declares no competing financial or academic interests. 

---

## REFERENCES
[1] Foulkes, W. D., Smith, I. E., & Reis-Filho, J. S. (2010). Triple-negative breast cancer. *New England Journal of Medicine*, 363(20), 1938-1948.  
[2] Keenan, T. E., & Tolaney, S. M. (2020). Role of immunotherapy in triple-negative breast cancer. *Journal of the National Comprehensive Cancer Network*, 18(4), 479-489.  
[3] Varley, K. E., Gertz, J., Roberts, B. S., Davis, N. S., et al. (2014). Recurrent read-through fusion transcripts in breast cancer. *Breast Cancer Research and Treatment*, 146(2), 287-297.  
[4] Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. *Bioinformatics*, 30(15), 2114-2120.  
[5] Kim, D., Paggi, J. M., Park, C., Bennett, C., & Salzberg, S. L. (2019). Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. *Nature biotechnology*, 37(8), 907-915.  
[6] Liao, Y., Smyth, G. K., & Shi, W. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics*, 30(7), 923-930.  
[7] Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome biology*, 15(12), 1-21.  
[8] Yu, G., Wang, L. G., Han, Y., & He, Q. Y. (2012). clusterProfiler: an R package for comparing biological themes among gene clusters. *Omics: a journal of integrative biology*, 16(5), 284-287.  
[9] Costa, R. L. B., Han, H. S., & Gradishar, W. J. (2018). Targeting the PI3K/AKT/mTOR pathway in triple-negative breast cancer: a review. *Breast cancer research and treatment*, 169(3), 397-406.  
[10] Gu-Trantien, C., Loi, S., Garaud, S., et al. (2013). CD4+ follicular helper T cell infiltration predicts breast cancer survival. *Journal of Clinical Investigation*, 123(7), 2873-2892.  
[11] Heinhuis, K. M., Ros, W., Kok, M., et al. (2019). Enhancing antitumor response by combining immune checkpoint inhibitors with chemotherapy in solid tumors. *Annals of Oncology*, 30(2), 219-235.  
[12] Sharma, P., & Allison, J. P. (2015). Immune checkpoint targeting in cancer therapy: toward combination strategies with curative potential. *Cell*, 161(2), 205-214.
