#!/bin/bash
mkdir -p /bioinformatics/module_4_tnbc_paper/reads/
cd /bioinformatics/module_4_tnbc_paper/reads/
rm -f *.fastq

# List of 10 real sample FTP links (PRJNA578488)
LINKS=(
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/058/SRR10314058/SRR10314058.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/059/SRR10314059/SRR10314059.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/067/SRR10314067/SRR10314067.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/090/SRR10314090/SRR10314090.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/093/SRR10314093/SRR10314093.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/097/SRR10314097/SRR10314097.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/007/SRR10314107/SRR10314107.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/042/SRR10314042/SRR10314042.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/049/SRR10314049/SRR10314049.fastq.gz"
"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/051/SRR10314051/SRR10314051.fastq.gz"
)

for LINK in "${LINKS[@]}"; do
    FILE=$(basename "$LINK")
    echo "Downloading Real File: $FILE"
    wget -c -q --show-progress "$LINK"
done

echo "GSE138914 Real Data Download Completed!"
/bioinformatics/miniconda3/envs/rna_seq_env/bin/seqkit stats *.fastq.gz > seqkit_report.txt
cat seqkit_report.txt
