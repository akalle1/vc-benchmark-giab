# Germline Variant Calling Benchmarking on ACMG84 Secondary Finding Genes

**Author:** Anagha Kalle  
**Institution:** Johns Hopkins University, MS Bioinformatics  
---

## Overview

Systematic benchmarking of open-source germline variant calling pipelines against GIAB v4.2.1 on clinically actionable ACMG84 secondary finding genes. This project bridges the gap between theoretical accuracy benchmarking and clinical feasibility by integrating diagnostic accuracy with normalized computational cost into a **Resource Efficiency Matrix** — a per-gene decision guide for clinical laboratories operating under resource constraints.

### Key Findings

- All three short-read tools collapse to identical INDEL recall of **0.733** in Lowmap/Segdup regions — confirming the bottleneck is read technology, not algorithm
- PacBio HiFi DeepVariant rescues INDEL F1 to **97.04%** in complex regions at only **1.67x** the per-base compute cost of short-read DeepVariant
- Exon-guided routing identifies **12 high-complexity genes** requiring HiFi and **72 standard genes** appropriate for short-read pipelines
- bcftools delivers **50.34 F1/core-hour** — 430x more efficient than DeepVariant in accessible regions

---

## Tools Benchmarked

| Tool | Version | Algorithmic Paradigm |
|------|---------|---------------------|
| GATK4 HaplotypeCaller | v4.2.6.1 | Probabilistic Hidden Markov Model |
| bcftools mpileup/call | v1.15.1 | Statistical/Bayesian |
| DeepVariant (Illumina) | v1.6.1 | Deep Learning (CNN) |
| DeepVariant (PacBio HiFi) | v1.6.1 | Deep Learning (CCS) |

---

## Methods Summary

### Sample and Reference
- **Sample:** HG002 Ashkenazi Jewish reference genome (300x WGS BAM, GRCh38)
- **Coverage:** Downsampled to 30x clinical standard (CAP/CLIA aligned)
- **Reference:** GRCh38 (UCSC/NCBI)
- **Truth set:** GIAB v4.2.1 (Wagner et al. 2022, Cell Genomics)
- **Cluster:** JHU Rockfish HPC (MARCC), account cbradbu3

### Benchmarking
- **Tool:** hap.py with vcfeval engine (GA4GH recommended haplotype-aware benchmarking)
- **Regions:** chr21 (validation baseline) → ACMG84 gene panel (81 assessable genes)
- **Stratification:** GA4GH v3.6 BED files (Easy, Difficult, Lowmap/Segdup)
- **Confidence:** All stratification BEDs intersected with GIAB v4.2.1 high-confidence regions

### Exon-Guided Routing
Coding exon coordinates obtained from UCSC Table Browser (GENCODE v49, knownGene table, hg38). Transcript-to-gene-symbol mapping via kgXref table. Exon coordinates intersected with GA4GH Lowmap/Segdup BED using bedtools intersect. Any gene with one or more coding exons in Lowmap/Segdup regions is assigned to the High-Complexity tier (HiFi required).

**Result:** 12 high-complexity genes / 72 standard genes (3 chrX genes excluded from analysis)

### Computational Profiling
Resource consumption measured using `/usr/bin/time -v` on Rockfish HPC. Normalized to **CPU-hours per million reads** (universal exchange rate for cross-depth comparison). HiFi DeepVariant normalized to **CPU-hours per Mbp** due to variable long-read length.

---

## Resource Efficiency Matrix

| Tier | Target Genes | Recommended Pipeline | INDEL F1 | Normalized Cost |
|------|-------------|---------------------|----------|-----------------|
| Low-Cost | 72 standard genes | bcftools | 0.755 | 0.012 CPU-hrs/M reads |
| Standard | Moderate complexity | GATK4 | 0.804 | 0.183 CPU-hrs/M reads |
| High-Complexity | 12 complex genes | HiFi DeepVariant | 0.970 | 0.055 CPU-hrs/Mbp |

**12 High-Complexity genes:** BMPR1A, BRCA1, CALM1, FLNC, MYH11, MYH7, PKP2, PMS2, PTEN, SCN5A, SDHC, SDHD

---

## Project Structure
vc-benchmark-giab/
├── Snakefile                                      # Main workflow
├── config/
│   ├── config.yaml                                # Pipeline parameters
│   ├── acmg84_exons_hg38.bed                     # ACMG84 coding exon coordinates (GENCODE v49)
│   └── knownGene_to_symbol.tsv                   # Transcript-to-gene symbol mapping (UCSC)
├── workflow/
│   ├── preprocessing.smk                          # Read alignment and preprocessing
│   ├── variant_calling.smk                        # GATK4, bcftools, DeepVariant rules
│   └── benchmarking.smk                           # hap.py benchmarking rules
├── stratifications/
│   ├── acmg84_easy.bed                            # Easy region intersections
│   ├── acmg84_difficult.bed                       # Difficult region intersections
│   ├── acmg84_lowmap_segdup.bed                   # Lowmap/Segdup intersections
│   ├── acmg84_easy_confident.bed                  # Easy + GIAB confidence
│   ├── acmg84_difficult_confident.bed             # Difficult + GIAB confidence
│   ├── acmg84_difficult_confident_corrected.bed   # Corrected BED (April 2025)
│   ├── acmg84_lowmap_segdup_confident.bed         # Lowmap/Segdup + GIAB confidence
│   ├── gene_complexity_breakdown.txt              # Per-gene complexity classification
│   ├── gene_easy_bp.txt                           # Base pairs per gene in Easy regions
│   ├── gene_difficult_bp.txt                      # Base pairs per gene in Difficult regions
│   └── gene_lowmap_bp.txt                         # Base pairs per gene in Lowmap/Segdup
├── run_stratifications.sh                         # Stratification benchmarking script
├── run_hifi_stratification.sh                     # HiFi stratification benchmarking script
└── README.md

---

## Quick Start

```bash
# 1. Test workflow (dry run)
snakemake -n

# 2. Run on SLURM
snakemake -j 10 --cluster "sbatch -p shared -A cbradbu3"

# 3. Run stratification benchmarking
bash run_stratifications.sh

# 4. Run HiFi stratification
bash run_hifi_stratification.sh

# 5. View results
cat results/benchmark_summary.txt
```

---

## Data Sources (not included — download separately)

The following large files must be downloaded before running the pipeline:

| File | Source |
|------|--------|
| HG002 300x BAM (GRCh38) | NCBI GIAB FTP |
| GIAB v4.2.1 truth VCF + BED | NIST GIAB FTP |
| GRCh38 reference genome | UCSC/NCBI |
| GRCh38_alldifficultregions.bed | NIST GIAB FTP (GA4GH v3.6) |
| GRCh38_alllowmapandsegdupregions.bed | NIST GIAB FTP (GA4GH v3.6) |
| GRCh38_notinalldifficultregions.bed | NIST GIAB FTP (GA4GH v3.6) |

---

## Dependencies

- Snakemake
- Singularity >= 3.8.7
- GATK4 v4.2.6.1
- bcftools v1.15.1
- DeepVariant v1.6.1
- hap.py + vcfeval
- bedtools
- BWA-MEM
- samtools v1.15.1

---

## Citations

1. Wagner J, et al. Benchmarking challenging small variants with linked and long reads. *Cell Genomics*. 2022;2(5):100128.
2. Krusche P, et al. Best practices for benchmarking germline small-variant calls in human genomes. *Nat Biotechnol*. 2019;37(5):555-560.
3. Poplin R, et al. A universal SNP and small-indel variant caller using deep neural networks. *Nature Biotechnology*. 2018;36:983-987.
4. Van der Auwera GA, et al. From FastQ data to high-confidence variant calls. *Curr Protoc Bioinformatics*. 2013;43:11.10.1-33.

---

## Contact

Questions? Email: akalle1@jh.edu
