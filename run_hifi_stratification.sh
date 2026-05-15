#!/bin/bash
#SBATCH --job-name=difficult_stratification
#SBATCH --time=12:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=12
#SBATCH --partition=shared
#SBATCH --account=cbradbu3
#SBATCH --output=logs/all_stratification_%j.log
cd /home/akalle1/scr4_cbradbu3
module load singularity/3.8.7
CONTAINER="docker://pkrusche/hap.py:latest"
TRUTH_VCF="giab_truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
REF="reference/GRCh38.fa"
for TOOL in gatk4 bcftools deepvariant deepvariant_hifi; do
    if [ "${TOOL}" == "deepvariant_hifi" ]; then
        QUERY_VCF="results/deepvariant_hifi/HG002_hifi_30x_acmg84.vcf.gz"
    else
        QUERY_VCF="results/${TOOL}/HG002_30x_acmg84.vcf.gz"
    fi
    for STRAT in acmg84_difficult; do
        STRAT_BED="stratifications/${STRAT}_confident_corrected.bed"
        OUT_PREFIX="results/benchmark/${TOOL}_HG002_30x_${STRAT}"
        echo "Running hap.py: ${TOOL} vs ${STRAT} at $(date)"
        singularity exec \
          --bind $(dirname ${QUERY_VCF}):/query:ro \
          --bind $(dirname ${TRUTH_VCF}):/truth:ro \
          --bind $(dirname ${STRAT_BED}):/truthbed:ro \
          --bind $(dirname ${REF}):/ref:ro \
          --bind results/benchmark:/output \
          ${CONTAINER} \
          /opt/hap.py/bin/hap.py \
          /truth/$(basename ${TRUTH_VCF}) \
          /query/$(basename ${QUERY_VCF}) \
          -f /truthbed/$(basename ${STRAT_BED}) \
          -r /ref/$(basename ${REF}) \
          -o /output/$(basename ${OUT_PREFIX}) \
          --engine=vcfeval \
          --threads=12
        echo "Finished: ${TOOL} vs ${STRAT} at $(date)"
    done
done
echo "All stratification benchmarks complete at $(date)"
