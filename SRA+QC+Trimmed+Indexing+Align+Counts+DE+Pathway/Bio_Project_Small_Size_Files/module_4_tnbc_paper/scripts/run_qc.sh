#!/bin/bash
cd /bioinformatics/module_4_tnbc_paper/reads/

echo "====== [1/2] 20 TNBC Human Samples - Initial FastQC ======"
mkdir -p fastqc_out

# Run FastQC on all 20 files using 20 threads to finish fast
/bioinformatics/miniconda3/envs/rna_seq_env/bin/fastqc -t 20 -o fastqc_out/ *.fastq.gz

echo "====== [2/2] Aggregating Raw MultiQC Report ======"
cd fastqc_out/
/bioinformatics/miniconda3/envs/rna_seq_env/bin/multiqc . -n raw_multiqc_report.html

echo "Initial QC completed. Report is at fastqc_out/raw_multiqc_report.html"
