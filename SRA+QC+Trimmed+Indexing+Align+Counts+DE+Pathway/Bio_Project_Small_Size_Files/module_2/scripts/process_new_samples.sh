#!/bin/bash
set -e

# Directories
DIR="/bioinformatics/module_2"
READS="${DIR}/reads"
TRIMMED="${DIR}/trimmed_reads"
HISAT2_DIR="${DIR}/hcc1395_hisat2"
REF_DIR="${DIR}/references/full_genome"

mkdir -p $READS $TRIMMED $HISAT2_DIR

SAMPLES=("SRR1553606" "SRR1553608")

echo "=== Mühit Yüklənir ==="
source /bioinformatics/miniconda3/etc/profile.d/conda.sh
conda activate rna_seq_env

for SAMPLE in "${SAMPLES[@]}"; do
    echo "======================================"
    echo "        İşlənən Nümunə: $SAMPLE       "
    echo "======================================"
    
    cd $READS
    echo "1. Downloading $SAMPLE (prefetch & fasterq-dump)..."
    prefetch $SAMPLE
    fasterq-dump --split-files -e 30 $SAMPLE
    
    echo "2. Trimming $SAMPLE (Trimmomatic)..."
    cd $TRIMMED
    trimmomatic PE -threads 30 \
        ${READS}/${SAMPLE}_1.fastq ${READS}/${SAMPLE}_2.fastq \
        ${SAMPLE}_1_paired.fastq ${SAMPLE}_1_unpaired.fastq \
        ${SAMPLE}_2_paired.fastq ${SAMPLE}_2_unpaired.fastq \
        ILLUMINACLIP:/bioinformatics/miniconda3/envs/rna_seq_env/share/trimmomatic-0.40-0/adapters/TruSeq3-PE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        
    echo "3. HISAT2 Alignment $SAMPLE..."
    cd $HISAT2_DIR
    hisat2 -p 30 -x ${REF_DIR}/hg38/genome \
        -1 ${TRIMMED}/${SAMPLE}_1_paired.fastq \
        -2 ${TRIMMED}/${SAMPLE}_2_paired.fastq \
        -S ${SAMPLE}.sam
        
    echo "4. Samtools Sort & Index $SAMPLE..."
    samtools sort -@ 30 -o ${SAMPLE}.bam ${SAMPLE}.sam
    samtools index -@ 30 ${SAMPLE}.bam
    
    # Clean up to save space
    rm ${TRIMMED}/${SAMPLE}_*unpaired.fastq ${SAMPLE}.sam
done

echo "=== Bütün 3 Nümunə Üzərində featureCounts İcrası ==="
cd $HISAT2_DIR
featureCounts -T 30 -p --countReadPairs \
    -a ${REF_DIR}/gencode.v44.annotation.gtf \
    -g gene_name \
    -o all_samples_expression.txt \
    SRR1553607.bam SRR1553606.bam SRR1553608.bam

echo "=== Yekun 'final_counts.csv' Faylının Təmizlənməsi ==="
# SRR1553607.bam, SRR1553606.bam, SRR1553608.bam faylları 7,8,9-cu sütunlardır
tail -n +2 all_samples_expression.txt | cut -f1,7,8,9 | sed 's/\t/,/g' > final_counts.csv

echo "======================================"
echo "    BÜTÜN PROSESLƏR UĞURLA BİTDİ!     "
echo "    Nəticə: final_counts.csv hazırdır."
echo "======================================"
