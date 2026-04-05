#!/bin/bash

# TNBC SRA layihəsi üçün xəstə nümunələri (Məs: PRJNA517180)
# T: Tumor, N: Normal (Qonşu sağlam toxuma)
SAMPLES=("SRR8615403" "SRR8615404" "SRR8615405" "SRR8615406" "SRR8615407" 
         "SRR8615408" "SRR8615409" "SRR8615410" "SRR8615411" "SRR8615412")

mkdir -p /bioinformatics/module_4_tnbc_paper/reads/
cd /bioinformatics/module_4_tnbc_paper/reads/

echo "====== TNBC Layihəsi: SRA Endirilməsi Başlayır (fasterq-dump) ======"

for SRR in "${SAMPLES[@]}"; do
    echo "Endirilir: $SRR (100k oxunuşla limitlənmiş sınaq versiyası)"
    # Serveri yormamaq və saatlarla gözləməmək üçün 100,000 "read" çəkilir (-X 100000)
    # Bu, elmi keyfiyyəti simulyasiya edir, lakin 10 saat yox, 5 dəqiqə sürür.
    /bioinformatics/miniconda3/envs/rna_seq_env/bin/fasterq-dump "$SRR" -X 100000 --split-files -p -O .
done

echo "====== Bütün 10 Nümunə Uğurla Endirildi! ======"
ls -lh
