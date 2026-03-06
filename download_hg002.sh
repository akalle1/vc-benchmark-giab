#!/bin/bash
#SBATCH --job-name=download_hg002
#SBATCH --partition=shared
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=logs/download_%j.log

echo "Starting download at $(date)"
echo "Current directory: $(pwd)"

cd /home/akalle1/scr4_cbradbu3/data

# Download with auto-retry
wget \
  --continue \
  --timeout=30 \
  --tries=0 \
  --retry-connrefused \
  --waitretry=10 \
  --progress=dot:giga \
  https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam

echo "Download complete at $(date)!"
ls -lh HG002.GRCh38.300x.bam
