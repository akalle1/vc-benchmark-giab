#!/bin/bash
#SBATCH --job-name=hifi_deepvariant
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
#SBATCH --account=cbradbu3
#SBATCH --output=logs/hifi_deepvariant_%j.log

cd /home/akalle1/scr4_cbradbu3
module load singularity/3.8.7

mkdir -p results/deepvariant_hifi

echo "Starting HiFi DeepVariant at $(date)"

{ /usr/bin/time -v \
singularity exec \
  --bind results/preprocessing:/input:ro \
  --bind reference:/ref:ro \
  --bind results/deepvariant_hifi:/output \
  docker://google/deepvariant:1.6.1 \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=PACBIO \
  --ref=/ref/GRCh38.fa \
  --reads=/input/HG002_hifi_30x_acmg84.bam \
  --output_vcf=/output/HG002_hifi_30x_acmg84.vcf.gz \
  --output_gvcf=/output/HG002_hifi_30x_acmg84.g.vcf.gz \
  --num_shards=8 ; } \
2> logs/hifi_deepvariant_time.log

echo "DeepVariant complete at $(date)"
