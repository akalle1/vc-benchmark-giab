#!/bin/bash
#SBATCH --job-name=hifi_happy
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
#SBATCH --account=cbradbu3
#SBATCH --output=logs/hifi_happy_%j.log


cd /home/akalle1/scr4_cbradbu3

# Load modules
module load singularity/3.8.7

mkdir -p results/benchmarking/deepvariant_hifi_acmg84

echo "Starting HiFi hap.py benchmarking at $(date)"


singularity exec \
  --bind reference:/ref:ro \
  --bind giab_truth:/truth:ro \
  --bind config:/conf:ro \
  --bind results/deepvariant_hifi:/query:ro \
  --bind results/benchmarking/deepvariant_hifi_acmg84:/output \
  docker://pkrusche/hap.py:latest \
  /opt/hap.py/bin/hap.py \
  /truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  /query/HG002_hifi_30x_acmg84.vcf.gz \
  -o /output/happy_hifi_acmg84 \
  -r /ref/GRCh38.fa \
  -f /truth/HG002_acmg84_confident.bed \
  --engine=vcfeval \
  --threads=8

echo "hap.py complete at $(date)"


