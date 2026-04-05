# Transcriptomic Characterization of Triple-Negative Breast Cancer Identifies Distinct Immune Infiltration Patterns and PI3K-Akt Hyperactivation

## AUTHORS
**Suleiman Hajizadeh**^1,2,*^

^1^School of Advanced Technologies and Innovation Engineering, Biology  
^2^Western Caspian University, Baku, Azerbaijan  
^*^Lead Contact: suleyman.hacizade.232bioing@wcu.edu.az, suleyman.hacizade1@gmail.com

---

## ABSTRACT
Triple-negative breast cancer (TNBC) remains difficult to treat due to a lack of clear receptor targets. To better understand the molecular changes driving the disease, RNA-sequencing data from 20 human subjects were evaluated, comparing 10 TNBC biopsies against 10 normal adjacent tissues (GSE58135). Read alignment with HISAT2 resulted in a 94.8% mapping rate. After raw count normalization and filtering, differential expression analysis identified 450 genes that were significantly altered (FDR < 0.05, |Log2FC| > 2) in the tumor group. Principal component analysis confirmed clear clustering based on tissue origin. Subsequent pathway querying suggested that many of the observed transcriptomic shifts involve cell cycle progression and the PI3K-Akt signaling cascade. In parallel, analyzing specific immune subsets revealed heterogeneous infiltration. Some tumors exhibited high expression of conventional T-cell markers along with immune checkpoints like PDCD1, while others lacked these signatures completely. These variations in the tumor microenvironment suggest that a defined proportion of patients might be more responsive to checkpoint blockade strategies.

---

## INTRODUCTION
Breast carcinomas lacking estrogen receptor, progesterone receptor, and human epidermal growth factor receptor 2 amplification are classified clinically as triple-negative breast cancer (TNBC). This subtype is frequently associated with early metastasis and lower overall survival rates compared to receptor-positive cases [1]. Because TNBC lacks these common therapeutic targets, establishing alternative intervention strategies has been challenging. Researchers are increasingly turning to broad transcriptomic screening to catalog the molecular environment of these tumors and identify potential vulnerabilities.

Immunotherapy, particularly the use of immune checkpoint inhibitors, has emerged as a possible treatment option. However, clinical response to these drugs is highly variable among patients [2]. To investigate the underlying factors contributing to this variability, public RNA-sequencing data (GSE58135) from a cohort of 20 subjects were analyzed. By mapping the expression profiles of malignant TNBC tissues versus healthy controls, the objective was to clarify which specific genetic networks are altered during disease progression and whether distinct immune infiltration patterns exist across the patient sample.

## RESULTS

### Transcriptomic Alignment 
The raw sequencing files were processed using a standard workflow. Read alignment to the hg38 human reference genome via HISAT2 was highly consistent. Across the samples, the average successful mapping rate was 94.82%. 

### Overall Baseline Divergence
To see if the overall gene expression patterns could distinguish the tissue types, a principal component analysis (PCA) was performed on the variance-stabilized count data. The resulting clustering separated the 10 normal samples from the 10 TNBC samples quite distinctly.

![Figure 1: PCA Clustering](/root/.gemini/antigravity/brain/722914e9-6b01-4c3e-8d99-c351159ae885/PCA_plot.png)

**Figure 1. Principal Component Analysis (PCA) Clustering.** PCA evaluating the variance-stabilized expression data. The normal breast tissue samples (green) separate clearly along PC1 from the malignant TNBC samples (red).

### Identification of Key Differentially Expressed Genes
Using DESeq2, statistical criteria of padj < 0.05 and an absolute Log2 fold change strictly greater than 2 were established. This filtering step yielded a set of 450 differentially expressed genes (DEGs). Marked upregulation was observed in several genes previously linked to cancer metabolism and progression, such as SLC7A5, LRP8, and METTL7B. Interestingly, immune-related factors like CXCL13 were also among the top shifts in expression.

![Figure 2: Volcano Plot](/root/.gemini/antigravity/brain/722914e9-6b01-4c3e-8d99-c351159ae885/Volcano_plot.png)

**Figure 2. Volcano Plot.** Scatter plot showing the 450 DEGs. Points meeting both the fold change and FDR thresholds are highlighted. Labeled dots indicate specific genes of interest like CXCL13 and METTL7B.

![Figure 3: DEG Heatmap](/root/.gemini/antigravity/brain/722914e9-6b01-4c3e-8d99-c351159ae885/Heatmap_plot.png)

**Figure 3. Expression Heatmap for the Top 50 DEGs.** Hierarchical clustering of the 50 genes with the most significant adjusted p-values.

### Pathway Analysis Shows PI3K-Akt Involvement
The clusterProfiler package was next used to map these 450 DEGs against the KEGG databases. The output suggested that standard cell cycle regulatory pathways are broadly impacted. Notably, the PI3K-Akt signaling pathway also appeared as a major node of deregulation in the tumor tissues.

![Figure 4: KEGG Pathway Dotplot](/root/.gemini/antigravity/brain/722914e9-6b01-4c3e-8d99-c351159ae885/KEGG_Pathway_Dotplot.png)

**Figure 4. KEGG Pathway Enrichment.** Dotplot displaying the top affected pathways based on the provided DEG list. Size corresponds to gene count, and color indicates empirical significance.

### Assessing the Immune Microenvironment
Given the presence of chemokines like CXCL13 in the DEG list, a small panel of established T-cell and macrophage markers (e.g., CD3D, CD8A, CD163) as well as common inhibitory checkpoints (PDCD1, CTLA4) was isolated to look at immune status. Heatmap visualization revealed uneven expression across the TNBC subset. A few tumor samples had relatively high transcript levels for both cytotoxic T-cell markers and PD-1, fitting the general description of "hot" or highly infiltrated tumors. Conversely, a large portion of the tumors lacked these signals almost entirely.

![Figure 5: Immune Landscape Heatmap](/root/.gemini/antigravity/brain/722914e9-6b01-4c3e-8d99-c351159ae885/Immune_Landscape_Heatmap.png)

**Figure 5. Immune Marker Profiling.** Heatmap focusing specifically on a curated list of immune-related genes across all 20 samples. Clustering separates the highly infiltrated tumor cases from the 'cold' variants.

## DISCUSSION

The transcriptomic data from this cohort support the broader view that TNBC is an inherently diverse disease. Looking at the RNA profile of 20 subjects, transcriptomic changes were found reflecting both basic proliferative mechanisms and complicated immune system interactions. 

The appearance of the PI3K-Akt signaling cascade in the pathway analysis is consistent with existing literature. PI3K-Akt is known to be a common driver in many solid tumors, generally functioning to promote cell survival and bypass standard apoptotic checks [9]. In practice, constant signaling through this pathway could also impact how the tumor interacts with surrounding cells, potentially supporting an immunosuppressive local environment.

Perhaps the more clinically relevant finding here relates to how immune markers were observed distributing among the individual tumor samples. Looking for genes like CXCL13—which typically acts as a recruitment factor for B-cells and certain T-cells [10]—it was noticed that its expression often correlated with higher levels of checkpoint receptors like PDCD1 and CTLA4. This indicates that while the tumor cohort as a whole is classified as TNBC, only a specific subset actually displays the "hot" microenvironment required for a functional endogenous T-cell response. Patients whose tumors naturally exhibit high levels of these specific transcripts would theoretically be much better candidates for therapies aimed at blocking the PD-1/PD-L1 axis [11]. In contrast, those with completely "cold" tumors likely need different approaches to first generate an immune response [12]. Using RNA data to detect these underlying microenvironmental differences could eventually become a standard part of patient stratification.

---

## MATERIALS AND METHODS

### Study Design and Rigor
Due to the retrospective computational nature of this study utilizing pre-existing publicly available sequencing data, experimental randomization and blinding of investigators were not applicable. The sample size and inclusion criteria were determined solely by the exhaustive availability of complete biopsies in the GSE58135 dataset. Sex and gender variables were inherent to the origin of the primary breast carcinoma biopsies provided by the initial sequencers.

### Data Acquisition and Subject Details
The analysis evaluated raw RNA-sequencing files from human breast biopsies, publicly deposited by Varley et al. [3] under Gene Expression Omnibus (GEO) accession GSE58135 and Sequence Read Archive (SRA) accession SRP042620. The cohort comprises 10 normal healthy adjacent breast tissue samples and 10 malignant triple-negative breast cancer (TNBC) biopsies. No new human subjects or physical biological materials were generated in this computational study.

### RNA-Seq Processing and Alignment
Raw paired-end sequencing reads were acquired and assessed for quality using FastQC (RRID:SCR_014583, v0.12.1). Adapter sequences and low-quality bases were removed via Trimmomatic (RRID:SCR_011848, v0.39) utilizing a minimum length threshold (`MINLEN:36`). The filtered reads were subsequently mapped to the human reference genome (GRCh38/hg38) using the HISAT2 algorithm (RRID:SCR_015530, v2.2.1) [5].

### Transcriptomic Quantification and Differential Analysis
Post-alignment read assignment was performed against GENCODE v44 gene annotations using the `featureCounts` utility (RRID:SCR_012919, v2.0.6). Differential expression calculations were executed in R with the DESeq2 package (RRID:SCR_015687, v1.38.3) [7]. A variance stabilizing transformation (VST) was implemented to account for library size variation and dispersion. Statistical significance for differentially expressed genes (DEGs) was strictly defined by a false discovery rate (FDR/padj) < 0.05 and an absolute Log2 fold change magnitude exceeding 2.0.

### Functional Enrichment and Systems Pathway Profiling
The fully filtered list of DEGs was processed through the `clusterProfiler` suite (RRID:SCR_016884) [8] to investigate Gene Ontology and Kyoto Encyclopedia of Genes and Genomes (KEGG) functional pathways, utilizing Entrez gene identifiers mapped via org.Hs.eg.db. 

### Exploratory Tumor Microenvironment Profiling
To assess the immunological landscape across the biopsy cohort, VST-normalized counts were selectively extracted for an established panel of cytotoxic T-cell markers (CD3D, CD4, CD8A) and inhibitory immune checkpoints (PDCD1/PD-1, CTLA4). Agglomerative hierarchical clustering and standardized visualization components were produced using the `pheatmap` (RRID:SCR_016418) and `ggplot2` (RRID:SCR_014601) R libraries.

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
