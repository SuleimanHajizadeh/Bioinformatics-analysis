#!/bin/bash
cd /bioinformatics/module_4_tnbc_paper/

echo "====== MƏRHƏLƏ 4: KƏMİYYƏTLƏNDİRMƏ (FEATURECOUNTS) BAŞLAYIR ======"

# 1. Modullar və Dəyişənlər
export FEATURECOUNTS_BIN="/bioinformatics/miniconda3/envs/rna_seq_env/bin/featureCounts"
export GTF_FILE="/bioinformatics/module_4_tnbc_paper/references/gencode.v44.annotation.gtf"

# 2. Gen Sayımı
# DİQQƏT: Pasiyent məlumatları Single-End (Tək-oxunuş) olduğu üçün '-p --countReadPairs' ləğv edildi!
# -T 15 (sürət üçün 15 nüvə) | -t exon | -g gene_name (Rəsmi gen adları üçün)
echo "1. FeatureCounts 20 BAM faylı üzrə saymağa başlayır..."
$FEATURECOUNTS_BIN -T 15 -a $GTF_FILE -o temp_counts.txt -t exon -g gene_name alignment/*.bam

# 3. Matrisin Təmizlənməsi (Stream Processing)
echo "2. Matris avtomatik (Pipe) ilə təmizlənir..."

# İlk proqramlama zibili və lazımsız koordinatlar sütunlarını ləğv edirik
# cut -f 1,7-26: 1-ci sütun (Gene_name) və 7-dən 26-a qədər (20 pasiyent)
cat temp_counts.txt | sed '1d' | cut -f 1,7-26 | tr '\t' ',' > tnbc_gene_counts.csv

# 4. Təmizlik (İstifadəçi əmri)
echo "3. Müvəqqəti hesablamalar təmizlənir..."
rm -f temp_counts.txt temp_counts.txt.summary

echo "====== MƏRHƏLƏ 4 TAMAMLANDI! tnbc_gene_counts.csv HAZIRDIR ======"
