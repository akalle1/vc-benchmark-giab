# Variant Calling Benchmarking Pipeline

Systematic benchmarking of GATK4, bcftools, and DeepVariant for germline variant calling.

**Author:** Anagha Kalle  
**Advisor:** Prof. Christopher Bradburne  
**Institution:** Johns Hopkins University, MS Bioinformatics  
**Date:** February 2026

---

## Overview

This pipeline benchmarks three variant calling tools against GIAB v4.2.1 standards:
- GATK4 HaplotypeCaller (probabilistic)
- bcftools mpileup + call (statistical)
- DeepVariant v1.6.1 (deep learning)

**Coverage:** 30x clinical WGS  
**Reference:** HG002 (Ashkenazi trio)  
**Focus:** ACMG 84 secondary findings genes

---

## Quick Start
```bash
# 1. Test workflow
snakemake -n

# 2. Run on SLURM
snakemake -j 10 --cluster "sbatch -p shared -A cbradbu3"

# 3. View results
cat results/benchmark_summary.txt
```

---

## Project Structure
```
.
├── Snakefile                 # Main workflow
├── config/
│   ├── config.yaml          # Parameters
│   └── samples.tsv          # Sample info
├── workflow/rules/          # Snakemake rules
├── data/                    # Input BAM files
├── reference/               # GRCh38 genome
├── giab_truth/             # GIAB truth sets
└── results/                # Output files
```

---

## Configuration

Edit `config/config.yaml` for your system:
```yaml
workdir: "/home/akalle1/scr4_cbradbu3/akalle1"
input:
  bam: "data/HG002_GRCh38_300x.bam"
```

---

## Tools Used

- samtools 1.15.1
- bcftools 1.15.1
- GATK 4.2.6.1
- DeepVariant 1.6.1
- hap.py (benchmarking)

---

## Output

**Results location:** `results/`

Key files:
- `benchmark_summary.txt` - Overall comparison
- `benchmark/*.summary.csv` - Per-tool metrics
- `gatk4/`, `bcftools/`, `deepvariant/` - Called variants

---

## Citation

If you use this workflow, please cite:
- GIAB: Wagner et al. (2022) Cell Genomics
- GATK: Van der Auwera et al. (2013)
- DeepVariant: Poplin et al. (2018) Nature Biotech

---

## Contact

Questions? Email: akalle1@jhu.edu
