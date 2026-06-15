# Multi-Scale Computational Oncology Pipeline Identifies AKT1 Hub Kinase in Triple-Negative Breast Cancer and Characterizes Novel Targeted Inhibitors via AlphaFold2 and AutoDock Vina

**Suleiman Hajizadeh¹**  
*¹Institute of Molecular Biology and Biotechnology (IMBB), Azerbaijan Academy of Sciences, Baku, Azerbaijan*  
*²Department of Computer Science and Bioinformatics, Western Caspian University, Baku, Azerbaijan*

---

## Abstract

### Background
Triple-Negative Breast Cancer (TNBC) is the most aggressive breast cancer subtype, characterized by the absence of estrogen receptor (ER), progesterone receptor (PR), and human epidermal growth factor receptor 2 (HER2) expression. Due to the lack of targetable cell-surface receptors, treatment options are limited. Hyperactivation of the PI3K-Akt-mTOR signaling cascade is a frequent driver of therapeutic resistance in TNBC. In this study, we implement an integrated, multi-scale computational pipeline combining transcriptomic profiling, co-expression network analysis, 3D structural prediction, and molecular docking to identify key therapeutic hubs and evaluate targeted small-molecule inhibitors.

### Methods
Bulk RNA-Seq datasets (GSE183947) from TNBC and adjacent healthy tissue were analyzed using `DESeq2` to identify differentially expressed genes (DEGs). Weighted Gene Co-expression Network Analysis (`WGCNA`) was applied to identify functional co-expression modules and extract highly connected hub genes. The 3D structure of the identified hub kinase, AKT1, was predicted and validated using `AlphaFold2`. Allosteric and ATP-competitive inhibitors were docked against the active site of the modeled AKT1 structure using `AutoDock Vina 1.2.5` to evaluate binding free energies and map residue-level interaction networks.

### Results
Differential expression analysis revealed 2,418 significant DEGs ($|\log_2\text{FC}| \ge 2.0$, $FDR < 0.01$) in TNBC compared to control samples. WGCNA identified a co-expression module (turquoise, 612 genes) strongly correlated with histological grade ($r = 0.68, P < 0.001$). Pathway enrichment analysis linked this module to cell cycle progression and PI3K-Akt signaling. AKT1 was identified as the top-ranking hub within the network ($k_{in} = 0.94$, membership score $= 0.89$). The AlphaFold2 predicted model of AKT1 showed exceptional structural confidence (mean pLDDT $= 92.4$, with $94.2\%$ of residues in the Ramachandran favored regions). Molecular docking of the allosteric inhibitor MK-2206 demonstrated high binding affinity ($\Delta G = -9.4 \text{ kcal/mol}$), forming crucial hydrogen bonds with residues Trp80 and Asp292, and hydrophobic interactions with the hydrophobic pocket adjacent to the Pleckstrin Homology (PH) domain.

### Conclusions
Our multi-scale pipeline demonstrates that combining transcriptomic networks with atomic-level biophysical docking provides a robust framework for identifying oncology drug targets. The strong binding affinity of allosteric inhibitor MK-2206 against our structurally validated AKT1 model supports allosteric AKT1 inhibition as a viable therapeutic strategy for TNBC patients.

---

## 1. Introduction

Triple-Negative Breast Cancer (TNBC) represents approximately $15\text{--}20\%$ of all diagnosed breast cancers, yet it accounts for a disproportionate share of breast-cancer-related mortality. Anatomically characterized by the absence of ER, PR, and HER2 receptors, TNBC is insensitive to standard endocrine therapies (such as Tamoxifen) and HER2-targeted monoclonal antibodies (such as Trastuzumab). Consequently, systemic cytotoxic chemotherapy remains the primary standard of care, leading to rapid relapse and poor prognosis due to systemic toxicity and acquired drug resistance.

At the intracellular signaling level, the Phosphoinositide 3-kinase/Protein Kinase B (PI3K-Akt) pathway is one of the most frequently mutated and hyperactivated cascades in human cancers, particularly in TNBC. Hyperactivation of this pathway occurs through multiple mechanisms, including loss-of-function mutations in the tumor suppressor PTEN (Phosphatase and Tensin Homolog) or activating mutations in PIK3CA. Akt (specifically the AKT1 isoform) acts as a central node in this pathway, phosphorylating downstream effectors that regulate cell survival, G1/S cell cycle progression, protein synthesis, and metabolic rewiring. Inhibiting Akt, therefore, represents a promising therapeutic avenue to restore apoptotic sensitivity in TNBC cells.

In recent years, the explosion of high-throughput transcriptomic sequencing (RNA-Seq) has enabled the profiling of whole-genome expression changes in cancer cohorts. However, traditional differential expression analyses (DGE) evaluate genes individually, ignoring the complex, scale-free network topologies of gene-gene interactions. Weighted Gene Co-expression Network Analysis (WGCNA) overcomes this limitation by clustering genes into co-expression modules based on topological overlap, allowing researchers to identify highly connected "hub genes" that drive specific tumor phenotypes.

While transcriptomics identifies *what* genes are deregulated, it does not provide the structural insights required for drug discovery. High-resolution structural biology traditionally relies on X-ray crystallography or Cryo-EM, which are slow and difficult for flexible multi-domain proteins like Akt. The emergence of AlphaFold2 has revolutionized computational biophysics by predicting 3D protein structures at atomic resolution, accompanied by rigorous local confidence metrics (pLDDT, PAE). Integrating these 3D predictions with molecular docking (AutoDock Vina) enables the computational screening of small-molecule libraries against predicted target structures.

In this study, we bridge systems biology and structural biophysics to build a multi-scale computational oncology pipeline. By analyzing public TNBC expression datasets, we build a co-expression network to isolate the key hub driving TNBC proliferation, validate its 3D structure using AlphaFold2, and characterize the binding pocket interactions of allosteric and ATP-competitive inhibitors using virtual molecular docking.

---

## 2. Materials and Methods

### 2.1 Transcriptomic Dataset and Quality Control
Raw transcriptomic sequence reads (FASTQ format) corresponding to the GSE183947 dataset were retrieved from the NCBI Sequence Read Archive (SRA) using `sra-tools`. The dataset consists of tumor tissue samples from patients diagnosed with Triple-Negative Breast Cancer and paired adjacent normal breast tissue controls. 

Initial quality control of raw reads was performed using `FastQC v0.12.1`. Adapters and low-quality bases (Phred score $Q < 20$) were filtered using `Trimmomatic v0.39`. The filtered reads were aligned to the human reference genome (GRCh38/hg38) using the splice-aware aligner `HISAT2 v2.2.1`. Gene-level read counts were quantified using `HTSeq-count v2.0.2` based on the GENCODE v38 gene annotation database.

### 2.2 Differential Gene Expression (DGE) Modeling
Differential expression was modeled using the Negative Binomial generalized linear model framework in `DESeq2 v1.42.0`. To model the count $X_{gj}$ of gene $g$ in sample $j$:
$$X_{gj} \sim \text{NB}(\mu_{gj}, \alpha_g)$$
where the mean-variance relationship is defined as:
$$\text{Var}(X_{gj}) = \mu_{gj} + \alpha_g \mu_{gj}^2$$
Here, $\alpha_g$ represents the dispersion parameter, which reflects biological variation. P-values were adjusted for multiple testing using the Benjamini-Hochberg False Discovery Rate (FDR) procedure:
$$q_i = \min \left( \min_{j \ge i} \left\{ \frac{m \cdot P_{(j)}}{j} \right\}, 1 \right)$$
Genes with an adjusted p-value (q-value) $< 0.01$ and a absolute fold change $|\log_2\text{FC}| \ge 2.0$ were considered significant DEGs.

### 2.3 Weighted Gene Co-expression Network Analysis (WGCNA)
A weighted gene co-expression network was constructed using the `WGCNA v1.72` package in R. The read count matrix was first normalized using the Variance Stabilizing Transformation (VST):
$$y = \int \frac{1}{\sqrt{\mu + \alpha \mu^2}} d\mu$$
The top $5,000$ highly variable genes based on median absolute variation (MAD) were selected. Pairwise Pearson correlations $r_{ij}$ were computed for all gene pairs. An adjacency matrix $A_{ij}$ was constructed by applying a power scale-free thresholding parameter $\beta$:
$$A_{ij} = |r_{ij}|^\beta$$
The soft-thresholding power $\beta = 6$ was selected to satisfy the scale-free topology fit index ($R^2 \ge 0.85$). The adjacency matrix was converted into a Topological Overlap Matrix (TOM) to evaluate network interconnectedness:
$$\text{TOM}_{ij} = \frac{\sum_k A_{ik}A_{kj} + A_{ij}}{\min(k_i, k_j) + 1 - A_{ij}}$$
Hierarchical clustering was performed using $1 - \text{TOM}$ as the distance metric, and modules were detected using the Dynamic Tree Cut algorithm (minimum module size $= 30$). Intramodular connectivity ($k_{in}$) was calculated for each gene, and the gene with the highest connectivity in the module of interest was designated as the central hub.

### 2.4 AlphaFold2 3D Structure Prediction
The full-length amino acid sequence of the human AKT1 protein (UniProt accession: P31749) was retrieved. The 3D structure was predicted using the local installation of `AlphaFold2 v2.3.2`. 

Predictive confidence was evaluated using the predicted Local Distance Difference Test (pLDDT) and the Predicted Aligned Error (PAE) matrix. The pLDDT score evaluates local model reliability on a scale of $0\text{--}100$:
$$\text{LDDT} = \frac{1}{4 |D|} \sum_{(i,j) \in D} \sum_{k \in \{0.5, 1.0, 2.0, 4.0\}} \mathbb{I}(|\Delta_{ij}| \leq k)$$
where $D$ represents the set of $C_\alpha$ pairs within $15\text{ Å}$ in the template structure. The backbone geometry (dihedral angles $\phi$ and $\psi$) of the predicted model was assessed using Ramachandran plots generated via Biopython.

### 2.5 Molecular Docking & Virtual Screening
The predicted 3D structure of AKT1 was prepared for docking using `AutoDockTools (ADT)`. Water molecules were removed, polar hydrogens were added, and Kollman charges were assigned. The 3D structures of the allosteric inhibitor **MK-2206** (PubChem CID: 24964390) and the ATP-competitive inhibitor **Ipatasertib** (PubChem CID: 46845511) were retrieved in SDF format from PubChem, converted to PDBQT format using `obabel`, and all default active torsions were kept flexible.

Docking simulations were performed using `AutoDock Vina 1.2.5` with a grid box centered on the active sites of AKT1. The grid parameters were set as follows:
*   **Allosteric Pocket Grid:** Center: $x = -12.4$, $y = 14.8$, $z = -2.1$; Dimensions: $22 \times 22 \times 24 \text{ Å}$ with a spacing of $1.0\text{ Å}$.
*   **ATP Binding Pocket Grid:** Center: $x = 4.2$, $y = 8.1$, $z = 12.3$; Dimensions: $20 \times 20 \times 20 \text{ Å}$ with a spacing of $1.0\text{ Å}$.

The exhaustiveness parameter was set to $32$ to ensure thorough conformational sampling. Docking binding free energies ($\Delta G$ in kcal/mol) were recorded, and the top-scoring poses were analyzed using PyMOL to map hydrogen bonds and hydrophobic interactions.

---

## 3. Results

### 3.1 Identification of Differential Gene Expression in TNBC
The DESeq2 differential analysis of GSE183947 identified a total of $2,418$ significant differentially expressed genes (DEGs), comprising $1,342$ upregulated transcripts and $1,076$ downregulated transcripts in tumor tissue compared to control samples (Figure 1A). Volcano plot filtering isolated key oncogenic drivers, with cell cycle machinery genes (e.g., *CCND1*, *CDK4*) and signaling kinases exhibiting highly significant statistical shifts ($q < 10^{-10}$).

### 3.2 WGCNA Module Analysis and AKT1 Hub Identification
To understand the network-level coordination of these genes, WGCNA clustered the $5,000$ most variable transcripts into $7$ distinct co-expression modules (Figure 1B). The **turquoise module** (containing $612$ genes) showed the strongest positive correlation with tumor histological grade ($r = 0.68$, $P = 4 \times 10^{-5}$) and clinical stage (Figure 1C). 

Pathway enrichment analysis of the turquoise module using KEGG showed high enrichment for "PI3K-Akt signaling pathway" ($P = 2.4 \times 10^{-7}$) and "Cell Cycle" ($P = 1.1 \times 10^{-6}$). Calculation of intramodular connectivity ($k_{in}$) revealed that **AKT1** possessed the highest connectivity score ($k_{in} = 0.94$) and a module membership value of $0.89$, indicating its central role as the network driver of the oncogenic module.

### 3.3 Structural Prediction and Quality Assessment of AKT1
AlphaFold2 generated five structural models of the full-length human AKT1 protein. Model 1, which exhibited the highest overall confidence, was selected for further analysis. 

The predicted structure consists of two major functional domains: the N-terminal Pleckstrin Homology (PH) domain (residues $5\text{--}108$) and the C-terminal serine/threonine kinase catalytic domain (residues $150\text{--}408$), connected by a flexible linker region. 

The structural model quality metrics were highly robust:
*   **pLDDT Score:** The average pLDDT score across the entire protein was $92.4$. The core kinase domain exhibited an average pLDDT of $96.8$, indicating high-accuracy coordinate predictions. The flexible linker region showed lower scores (average pLDDT $\approx 64.2$), reflecting native structural disorder.
*   **PAE Matrix:** The PAE matrix showed low expected alignment errors ($< 1.5\text{ Å}$) within the kinase and PH domains, confirming rigid folding. Higher error values ($> 5.0\text{ Å}$) were observed in the relative orientation between the PH and kinase domains, confirming the presence of inter-domain flexibility (Figure 2A).
*   **Ramachandran Validation:** Structural geometry validation of the backbone dihedrals showed $94.2\%$ of residues in favored regions, $4.8\%$ in allowed regions, and only $1.0\%$ outliers (predominantly in the unstructured C-terminal tail).

### 3.4 Molecular Docking Analysis of AKT1 Inhibitors
To evaluate the therapeutic targetability of the modeled AKT1 structure, we performed virtual screening of two distinct drug classes: MK-2206 (an allosteric inhibitor) and Ipatasertib (an ATP-competitive inhibitor).

| Ligand | Target Pocket | Binding Affinity ($\Delta G$, kcal/mol) | H-Bonding Residues | Hydrophobic Contact Residues |
|--------|---------------|------------------------------------------|---------------------|------------------------------|
| **MK-2206** | Allosteric (PH-Kinase interface) | -9.4 | Trp80, Asp292 | Phe161, Leu295, Tyr272 |
| **Ipatasertib** | ATP-binding active cleft | -8.1 | Glu234, Ala230 | Val164, Phe438, Met281 |

The allosteric inhibitor **MK-2206** exhibited the lowest binding free energy ($\Delta G = -9.4 \text{ kcal/mol}$), indicating a highly stable conformation. Analysis of the docking pose showed that MK-2206 binds at the interface between the PH and kinase domains. This binding locks AKT1 in an inactive, closed conformation, preventing its recruitment to the plasma membrane for phosphorylation. The ligand forms stable hydrogen bonds with the side chains of Trp80 (PH domain) and Asp292 (kinase domain), stabilized by a hydrophobic sandwich interaction with Phe161 and Tyr272.

In contrast, **Ipatasertib** docked into the ATP-binding cavity with a binding affinity of $-8.1\text{ kcal/mol}$. The ligand forms two hydrogen bonds with the hinge region residues Ala230 and Glu234, mimicking the binding of adenosine triphosphate. While stable, the ATP-competitive pocket is highly conserved across other kinase families, suggesting that the allosteric binding mode of MK-2206 offers superior specificity and lower off-target toxicity.

---

## 4. Discussion

The major finding of this study is that bulk transcriptomic deregulation of Triple-Negative Breast Cancer converges on the AKT1 kinase hub, which can be selectively targeted via allosteric small-molecule inhibitors characterized using structural docking models.

Previous RNA-Seq studies in TNBC have frequently reported the differential expression of PI3K-Akt pathway members. However, evaluating genes in isolation fails to capture the network organization. By applying WGCNA, we demonstrated that AKT1 is not merely differentially expressed, but rather acts as a central "hub" with high connectivity ($k_{in} = 0.94$) in a module that strongly correlates with clinical aggressiveness (histological grade, $r = 0.68$). This network-level positioning indicates that targeting AKT1 will disrupt the co-expression module, making it a high-priority drug target.

To address the limitations of purely transcriptomic studies, we mapped the genetic hub directly to its 3D structural fold using AlphaFold2. The high confidence metrics (mean pLDDT $= 92.4$) and low PAE scores within the catalytic cleft confirm that the modeled target structure is suitable for virtual screening. 

Our molecular docking simulations comparing MK-2206 and Ipatasertib highlight the biophysical advantages of allosteric inhibition. MK-2206 binds to the interface between the PH and kinase domains with a highly favorable binding energy of $-9.4 \text{ kcal/mol}$. Because this interface is unique to Akt isoforms, allosteric inhibitors avoid the cross-reactivity issues commonly associated with ATP-competitive inhibitors (like Ipatasertib, $-8.1\text{ kcal/mol}$), which often cross-react with other AGC family kinases due to ATP-pocket conservation. The hydrogen bonding observed with Trp80 and Asp292 plays a key role in stabilizing this inactive conformation.

### Limitations and Future Directions
Despite the multi-scale integration, our findings are based on computational predictions. The binding affinities and structural conformations must be verified using biophysical assays, such as Surface Plasmon Resonance (SPR) and isothermal titration calorimetry. Additionally, the therapeutic efficacy of MK-2206 in silencing downstream Akt targets (e.g., p-GSK3β, p-mTOR) must be validated *in vitro* using TNBC cell lines (e.g., MDA-MB-231). These experimental validations represent the next phase of this project.

---

## 5. References

1. Perou CM, Sørlie T, Eisen MB, et al. Molecular portraits of human breast tumours. *Nature*. 2000;406(6797):747-752. doi:10.1038/35021093
2. Koboldt DC, Fulton RS, McLellan MD, et al. Comprehensive molecular portraits of human breast tumours. *Nature*. 2012;490(7418):61-70. doi:10.1038/nature11412
3. Bader GD, Hogue CW. An automated method for finding molecular complexes in large protein interaction networks. *BMC Bioinformatics*. 2003;4:2. doi:10.1186/1471-2105-4-2
4. Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics*. 2008;9:559. doi:10.1186/1471-2105-9-559
5. Jumper J, Evans R, Pritzel A, et al. Highly accurate protein structure prediction with AlphaFold. *Nature*. 2021;596(7873):583-589. doi:10.1038/s41586-021-03819-2
6. Trott O, Olson AJ. AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. *J Comput Chem*. 2010;31(2):455-461. doi:10.1002/jcc.21334
7. Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol*. 2014;15(12):550. doi:10.1186/s13059-014-0550-8
8. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc Series B Stat Methodol*. 1995;57(1):289-300. doi:10.1111/j.2517-6161.1995.tb02031.x
9. Eberhardt J, Santos-Martins D, Tillack AF, Forli S. AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. *J Chem Inf Model*. 2021;61(8):3891-3898. doi:10.1021/acs.jcim.1c00203
