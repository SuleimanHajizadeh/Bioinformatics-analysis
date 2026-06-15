#!/bin/bash
cd /bioinformatics/module_4_tnbc_paper/reads/

mkdir -p trimmed_reads
mkdir -p trimmed_fastqc

echo "====== TRIMMOMATIC (GNU PARALLEL 20 SAMPLES) ======"

# Trimmomatic environment variables
export TRIMMOMATIC_BIN="/bioinformatics/miniconda3/envs/rna_seq_env/bin/trimmomatic"

# Custom bash function to process individual FASTQ files
trim_file() {
    FILE=$1
    BASENAME=$(basename "$FILE" .fastq.gz)
    echo "[$(date +'%T')] Qayçıya Verildi: $BASENAME"
    
    $TRIMMOMATIC_BIN SE -threads 2 -phred33 \
    "$FILE" "trimmed_reads/${BASENAME}_trimmed.fastq.gz" \
    ILLUMINACLIP:adapters.fa:2:30:10 \
    SLIDINGWINDOW:4:20 MINLEN:36 > "trimmed_reads/${BASENAME}_trim.log" 2>&1
    
    echo "[$(date +'%T')] Tamamlandı: $BASENAME"
}
export -f trim_file

# Execute using GNU Parallel (10 files simultaneously, roughly 20 threads total)
ls *.fastq.gz | /bioinformatics/miniconda3/envs/rna_seq_env/bin/parallel -j 10 trim_file {}

echo "====== POST-TRIMMING FastQC & MultiQC ======"
/bioinformatics/miniconda3/envs/rna_seq_env/bin/fastqc -t 20 -o trimmed_fastqc/ trimmed_reads/*_trimmed.fastq.gz

cd trimmed_fastqc/
/bioinformatics/miniconda3/envs/rna_seq_env/bin/multiqc . -n post_trim_multiqc_report.html

echo "====== Bütün Mərhələ 2 Tamamlandı! ======"
