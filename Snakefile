# Variant Calling Benchmarking Pipeline
# Author: Anagha Kalle
# Date: February 2026

import pandas as pd
from pathlib import Path
import os

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Load configuration
configfile: "config/config.yaml"

# Set working directory
workdir: config["workdir"]

# Load sample information
samples_df = pd.read_csv("config/samples.tsv", sep="\t")
#Strip whitespace
samples_df.columns = samples_df.columns.str.strip()
samples_df = samples_df.apply(
    lambda x: x.str.strip() if x.dtype == "object" else x
)

# Extract unique values
SAMPLES = samples_df["sample"].unique().tolist()
COVERAGES = samples_df["coverage"].unique().tolist()
REGIONS = samples_df["region"].unique().tolist()
TOOLS = config["tools"]

wildcard_constraints:
  sample = "[^_]+",
  coverage = "\d+",
  region = "chr21|acmg84",
  tool = "gatk4|bcftools|deepvariant"

#helper function
#get bedfile for the region
def get_region_bed(wildcards):
    region_info = config["regions"][wildcards.region]
    return region_info["bed"]	

#get bed file for the hap.py benchmarking tool
def get_benchmark_region(wildcards):
    """Get GIAB confidence BED for hap.py benchmarking"""
    if wildcards.region == "chr21":
        return config["giab"]["bed"]
    elif wildcards.region == "acmg84":
        return config["regions"]["acmg84"]["bed"]
    else:
        return config["giab"]["bed"]

#contract rule- target workflow
rule all:
  input:
    #Downsampled BAM
    expand(
    "results/preprocessing/{sample}_{coverage}x_{region}.bam",
    sample=SAMPLES, coverage=COVERAGES, region=REGIONS),
    #Variant caller
    expand("results/{tool}/{sample}_{coverage}x_{region}.vcf.gz",
      tool=TOOLS, sample=SAMPLES, coverage=COVERAGES, region=REGIONS),
    #Benchmarking results
    expand("results/benchmark/{tool}_{sample}_{coverage}x_{region}.summary.csv",
      tool=TOOLS, sample=SAMPLES, coverage=COVERAGES, region=REGIONS),     
    #Summary result
    "results/benchmark_summary.txt"



#Include (modular pipeline design)
#splt the files logical files so that that Snakemakefile is readable
include: "workflow/rules/preprocessing.smk"
include: "workflow/rules/variant_calling.smk"
include: "workflow/rules/benchmarking.smk"
