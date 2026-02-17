rule benchmarking_variants:
    input:
        query_vcf = "results/{tool}/{sample}_{coverage}x_{region}.vcf.gz",
        query_idx = "results/{tool}/{sample}_{coverage}x_{region}.vcf.gz.tbi",
        truth_vcf = config["giab"]["vcf"],
        truth_idx = config["giab"]["vcf_index"],
        truth_bed = config["giab"]["bed"],
        ref = config["reference"]["genome"],
        ref_fai = config["reference"]["fai"]
    output:
        summary   = "results/benchmark/{tool}__{sample}__{coverage}x__{region}.summary.csv",
        extended  = "results/benchmark/{tool}__{sample}__{coverage}x__{region}.extended.csv",
        vcf       = "results/benchmark/{tool}__{sample}__{coverage}x__{region}.vcf.gz"
    params:
        singularity_module = config["modules"]["singularity"],
        container = "docker://pkrusche/hap.py:latest",
        prefix = lambda wc: f"results/benchmark/{wc.tool}__{wc.sample}__{wc.coverage}x__{wc.region}"
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

