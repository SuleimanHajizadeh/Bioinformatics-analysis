#!/bin/bash
# 20 Nümunə (Homo Sapiens - TNBC vs Normal) Endirmə Skripti

mkdir -p /bioinformatics/module_4_tnbc_paper/reads/
cd /bioinformatics/module_4_tnbc_paper/reads/

# TNBC Tumor
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR131/000/SRR1313090/SRR1313090_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR131/001/SRR1313091/SRR1313091_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR131/002/SRR1313092/SRR1313092_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR131/004/SRR1313094/SRR1313094_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR131/006/SRR1313096/SRR1313096_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR131/009/SRR1313099/SRR1313099_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR131/001/SRR1313101/SRR1313101_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR131/001/SRR1313111/SRR1313111_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR131/002/SRR1313112/SRR1313112_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR131/003/SRR1313113/SRR1313113_1.fastq.gz

# Normal (Uninvolved)
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR131/005/SRR1313115/SRR1313115_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR131/008/SRR1313118/SRR1313118_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR131/000/SRR1313120/SRR1313120_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR131/003/SRR1313123/SRR1313123_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR131/007/SRR1313127/SRR1313127_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR131/008/SRR1313128/SRR1313128_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR131/007/SRR1313067/SRR1313067_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR131/008/SRR1313068/SRR1313068_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR131/008/SRR1313078/SRR1313078_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR131/006/SRR1313066/SRR1313066_1.fastq.gz

echo "20 FASTQ endirməsi bitdi!"
