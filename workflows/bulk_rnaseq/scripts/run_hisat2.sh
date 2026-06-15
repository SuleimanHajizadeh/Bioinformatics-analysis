#!/bin/bash
cd /bioinformatics/module_4_tnbc_paper/

mkdir -p alignment

echo "====== MƏRHƏLƏ 3: HISAT2 ALIGNMENT BAŞLAYIR (GNU PARALLEL) ======"

export SUMMARY_LOG="/bioinformatics/module_4_tnbc_paper/alignment/hisat2_summary.log"
> $SUMMARY_LOG
echo "HISAT2 Master Başlanğıc: $(date)" > $SUMMARY_LOG

export HISAT2_BIN="/bioinformatics/miniconda3/envs/rna_seq_env/bin/hisat2"
export SAMTOOLS_BIN="/bioinformatics/miniconda3/envs/rna_seq_env/bin/samtools"
# The 6GB hg38 human index we built in module 2!
export INDEX="/bioinformatics/module_2/references/full_genome/hg38/genome"

align_sample() {
    FILE=$1
    BASENAME=$(basename "$FILE" _trimmed.fastq.gz)
    
    if [ -s "/bioinformatics/module_4_tnbc_paper/alignment/${BASENAME}.bam" ]; then
        echo "[$(date +'%T')] KESİLDİ (Artıq Hazırdır): $BASENAME"
        return 0
    fi
    
    echo "=================================================" >> $SUMMARY_LOG
    echo "[$(date +'%T')] HISAT2 Başladı: $BASENAME" >> $SUMMARY_LOG
    echo "[$(date +'%T')] Xəritələnir: $BASENAME"
    
    # Run HISAT2 (-U for single-end), stream output directly to Samtools View (to BAM) and Samtools Sort to save disk space
    $HISAT2_BIN -p 2 -x $INDEX -U /bioinformatics/module_4_tnbc_paper/reads/trimmed_reads/$FILE 2> "/bioinformatics/module_4_tnbc_paper/alignment/${BASENAME}_align.log" | \
    $SAMTOOLS_BIN view -bS - | \
    $SAMTOOLS_BIN sort -@ 2 -o "/bioinformatics/module_4_tnbc_paper/alignment/${BASENAME}.bam"
    
    # Index the sorted BAM file (.bam.bai)
    $SAMTOOLS_BIN index -@ 2 "/bioinformatics/module_4_tnbc_paper/alignment/${BASENAME}.bam"
    
    echo "[$(date +'%T')] HISAT2 Tamamlandı: $BASENAME (BAM Indexed)" >> $SUMMARY_LOG
    echo "[$(date +'%T')] Tamamlandı: $BASENAME"
}
export -f align_sample

cd /bioinformatics/module_4_tnbc_paper/reads/trimmed_reads/
# GNU Parallel limits to 3 parallel jobs (Golden ratio for 32GB RAM: 3x6GB + Samtools = 22.5GB RAM)
ls *_trimmed.fastq.gz | /bioinformatics/miniconda3/envs/rna_seq_env/bin/parallel -j 3 align_sample {}

echo "====== BÜTÜN MƏRHƏLƏ 3 BİTDİ! BÜTÜN 20 XƏSTƏ HG38-ə QÜSURSUZ DÜZÜLDÜ ======"
