#!/bin/bash
#SBATCH --job-name=chr21_pipeline
#SBATCH --time=35:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --partition=shared
#SBATCH --output=logs/chr21_%j.log

cd /home/akalle1/scr4_cbradbu3

# Load modules (use snakemake/7.6.0 NOT 7.3.8!)
module load snakemake/7.6.0
module load samtools/1.15.1
module load bcftools/1.15.1
module load gatk/4.2.6.1-Java-11
module load singularity/3.8.7

echo "Starting chr21 pipeline at $(date)"

# Verify modules loaded
echo "Loaded modules:"
module list

echo "Tool paths:"
which snakemake
which samtools
which bcftools

snakemake \
  --cores 16 \
  --use-singularity \
  --printshellcmds \
  --keep-going

echo "Pipeline complete at $(date)"
