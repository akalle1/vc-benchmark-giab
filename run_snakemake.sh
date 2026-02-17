#!/bin/bash
#SBATCH --job-name=smk_test
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=parallel

module load snakemake/7.3.8

snakemake -n -p

