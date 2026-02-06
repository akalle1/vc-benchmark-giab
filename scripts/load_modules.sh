#!/bin/bash
#SBATCH --job-name=test_modules
#SBATCH --account=cbradbu3
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

#Running on a computer node

module load samtools/1.15.1
module load bcftools/1.15.1
module load gatk/4.2.6.1-Java-11
module load singularity/3.8.7

#test versions
samtools --version
bcftools --version
gatk --version
singularity --version
