#!/bin/bash

# ==============================================
# RNA-Seq Project Setup Script
# ==============================================
# Author: Suleyman
# Description: Sets up project structure and prepares
#              directories for RNA-seq pipeline.
#              Conda environment: rna_seq_env
#              Conda location: /bioinformatics/miniconda3
# ==============================================

set -e

BASE_DIR="/bioinformatics"
PROJECT="rna_seq_project"
ENV_NAME="rna_seq_env"
CONDA_DIR="$BASE_DIR/miniconda3"

echo "🚀 RNA-Seq layihe qurashdirilmasi bashlayir..."

# 1. Conda-ni yoxla
if [ ! -f "$CONDA_DIR/bin/conda" ]; then
    echo "❌ Conda tapilmadi: $CONDA_DIR"
    echo "   Zehmet olmasa /bioinformatics/miniconda3 qovlugunu yoxlayin."
    exit 1
fi

echo "✅ Conda tapildi: $CONDA_DIR"

# 2. Conda-ni aktivleshdir
source "$CONDA_DIR/etc/profile.d/conda.sh"

# 3. rna_seq_env mühitini yoxla
if conda info --envs | grep -q "$ENV_NAME"; then
    echo "✅ Conda muhiti movcuddur: $ENV_NAME"
else
    echo "🌱 Conda muhiti yaradilir: $ENV_NAME"
    conda create -y -n $ENV_NAME \
        sra-tools fastqc multiqc trimmomatic \
        hisat2 samtools subread entrez-direct \
        -c conda-forge -c bioconda
    echo "✅ $ENV_NAME muhiti yaradildi"
fi

# 4. Layihe qovlug strukturunu yarat
echo "📁 Qovluq strukturu yaradilir..."
mkdir -p $BASE_DIR/projects/$PROJECT/{raw_data,qc/fastqc_raw,qc/multiqc,trimmed,aligned,counts,logs,scripts}
mkdir -p $BASE_DIR/module_1/sra_data
echo "✅ Qovluqlar yaradildi: $BASE_DIR/projects/$PROJECT"

# 5. README yarat
cat > $BASE_DIR/projects/$PROJECT/scripts/README.md << 'EOF'
# RNA-Seq Pipeline Skriptleri

## Mühiti aktivlesdir:
```bash
conda activate rna_seq_env
```

## Esas alətlər:
- fastqc      : Keyfiyyət yoxlaması
- multiqc     : Hesabatların birləşdirilməsi
- trimmomatic : Adapter kəsilməsi
- hisat2      : Uyğunlaşdırma (alignment)
- samtools    : SAM/BAM faylların idarəsi
- featureCounts (subread): Sayma (counting)
- fasterq-dump (sra-tools): SRA faylların endirilməsi
EOF

echo ""
echo "🎉 Qurashdirilma tamamlandi!"
echo ""
echo "Nowbeti addimlar:"
echo "  1. conda activate $ENV_NAME"
echo "  2. cd $BASE_DIR/projects/$PROJECT"
echo "  3. Analize bashlayın!"
