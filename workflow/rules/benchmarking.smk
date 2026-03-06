
def get_truth_bed(wildcards):
    if wildcards.region == "chr21":
        return config["giab"]["bed_chr21"]
    else:
        return config["giab"]["bed"]


rule benchmark_variants:
    input:
        query_vcf = "results/{tool}/{sample}_{coverage}x_{region}.vcf.gz",
        query_idx = "results/{tool}/{sample}_{coverage}x_{region}.vcf.gz.tbi",
        truth_vcf = config["giab"]["vcf"],
        truth_idx = config["giab"]["vcf_index"],
        truth_bed = get_truth_bed,
        ref = config["reference"]["genome"],
        ref_fai = config["reference"]["fai"]
    output:
        summary   = "results/benchmark/{tool}_{sample}_{coverage}x_{region}.summary.csv",
        extended  = "results/benchmark/{tool}_{sample}_{coverage}x_{region}.extended.csv",
        vcf       = "results/benchmark/{tool}_{sample}_{coverage}x_{region}.vcf.gz"
    params:
        singularity_module = config["modules"]["singularity"],
        container = "docker://pkrusche/hap.py:latest",
        prefix = lambda wc: f"results/benchmark/{wc.tool}_{wc.sample}_{wc.coverage}x_{wc.region}"
    threads: config["resources"]["benchmark"]["cpus"]
    resources:
        mem_mb = config["resources"]["benchmark"]["mem_mb"],
        time_min = config["resources"]["benchmark"]["time_min"],
        partition = config["resources"]["default"]["partition"],
        account = config["resources"]["default"]["account"]
    log:
        "logs/benchmark/{tool}_{sample}_{coverage}x_{region}.log"
    shell:
        r"""
        module load {params.singularity_module}

        singularity exec \
          --bind $(dirname {input.query_vcf}):/query:ro \
          --bind $(dirname {input.truth_vcf}):/truth:ro \
          --bind $(dirname {input.truth_bed}):/truthbed:ro \
          --bind $(dirname {input.ref}):/ref:ro \
          --bind $(dirname {output.summary}):/output \
          {params.container} \
          /opt/hap.py/bin/hap.py \
          /truth/$(basename {input.truth_vcf}) \
          /query/$(basename {input.query_vcf}) \
          -f /truthbed/$(basename {input.truth_bed}) \
          -r /ref/$(basename {input.ref}) \
          -o /output/$(basename {params.prefix}) \
          --threads {threads} \
          --engine=vcfeval \
          2>> {log}
        """


# AGGREGATE RESULTS

rule aggregate_benchmarks:
    """Combine all benchmark results into summary report"""
    input:
        summaries = expand(
            "results/benchmark/{tool}_{sample}_{coverage}x_{region}.summary.csv",
            tool=TOOLS,
            sample=SAMPLES,
            coverage=COVERAGES,
            region=REGIONS
        )
    output:
        report = "results/benchmark_summary.txt"
    log:
        "logs/benchmark/aggregate.log"
    run:
        with open(output.report, 'w') as f:
            f.write("VARIANT CALLING BENCHMARK SUMMARY\n")  
            f.write(f"Tools benchmarked: {', '.join(TOOLS)}\n")
            f.write(f"Sample: {', '.join(SAMPLES)}\n")
            f.write(f"Regions: {', '.join(REGIONS)}\n")
            f.write(f"Total comparisons: {len(input.summaries)}\n\n")
            
            for summary_file in input.summaries:
                f.write(f"\n{summary_file}:\n")
                try:
                    with open(summary_file) as sf:
                        content = sf.read()
                        f.write(content)
                        if not content.endswith('\n'):
                            f.write('\n')
                except Exception as e:
                    f.write(f"Error reading file: {e}\n")
            
            f.write("END OF SUMMARY\n")
