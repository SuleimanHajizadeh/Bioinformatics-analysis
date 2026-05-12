# 🧬 Personal Genomics & Clinical Annotation Pipeline

## 📌 Project Overview
This repository contains a full bioinformatics pipeline for processing, filtering, and medically annotating Direct-to-Consumer (DTC) raw DNA data (MyHeritage, 23andMe, AncestryDNA). The pipeline cross-references raw genotypic data containing over 600,000 SNPs with the NIH ClinVar Database, accurately highlighting drug responses, pathogenic risks, and protective traits.

## 🛠️ Tools & Technologies Used
* **Python**: Core logic and pipeline construction.
* **Pandas**: Used to manipulate, merge, and filter massive genomic datasets efficiently.
* **NIH ClinVar**: Large-scale public medical genetics database integration.
* **Bash/Linux**: Data handling and gzip extraction.

## 🔬 Files and Workflow
1. `analyze_dna.py`: 
   - A fast lookup script that scans the user's raw DNA against a hand-curated dictionary of high-impact traits (e.g., lactose tolerance, COMT "Worrier" gene, caffeine metabolism, PRNP resistance).
2. `clinvar_annotator.py`:
   - An automated tool that downloads the latest version of the NIH ClinVar database (variant_summary.txt.gz).
   - Uses Pandas `merge` to intersect the user's ~609,000 SNPs with ClinVar's millions of records.
   - Cleans the data and extracts variants classified as *Drug Response*, *Protective*, *Pathogenic*, and *Risk Factors*.
3. `MyHeritage_ClinVar_Matches.csv`:
   - The final output showing 21,894 intersecting clinical SNPs.

## 🎓 Academic Goal
This repository serves as a practical demonstration of **Data Science in Structural Biology & Computational Genomics**, aligning with the entry requirements for advanced graduate programs (MPhil in Computational Biology) at institutions such as the University of Cambridge.
