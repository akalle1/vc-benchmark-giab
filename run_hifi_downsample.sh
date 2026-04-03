#!/bin/bash
#SBATCH --job-name=hifi_downsample
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#SBATCH --account=cbradbu3
#SBATCH --output=logs/hifi_downsample_%j.log
#SBATCH --error=logs/hifi_downsample_%j.err

cd /home/akalle1/scr4_cbradbu3
module load samtools/1.15.1

# Delete any corrupted output first
rm -f results/preprocessing/HG002_hifi_30x_acmg84.bam
rm -f results/preprocessing/HG002_hifi_30x_acmg84.bam.bai

echo "Starting downsampling at $(date)"

samtools view \
    -b \
    -@ 4 \
    -s 42.625 \
    -o results/preprocessing/HG002_hifi_30x_acmg84.bam \
    results/preprocessing/HG002_hifi_acmg84.bam

echo "Indexing at $(date)"
samtools index -@ 4 results/preprocessing/HG002_hifi_30x_acmg84.bam

echo "Checking coverage at $(date)"
samtools depth \
    -b config/acmg84_hg38.bed \
    results/preprocessing/HG002_hifi_30x_acmg84.bam \
    | awk '{sum+=$3; count++} END {print "Mean coverage:", sum/count}'

echo "Done at $(date)"
