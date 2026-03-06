# Variant Calling Rules

rule gatk4_haplotypecaller:
    input:
        bam = "results/preprocessing/{sample}_{coverage}x_{region}.rg.bam",
        bai = "results/preprocessing/{sample}_{coverage}x_{region}.rg.bam.bai",
        ref = config["reference"]["genome"],
        ref_dict = config["reference"]["dict"]
    output:
        vcf = "results/gatk4/{sample}_{coverage}x_{region}.vcf.gz",
        idx = "results/gatk4/{sample}_{coverage}x_{region}.vcf.gz.tbi"
    params:
        java_mem = config["gatk4"]["java_mem"],
        extra = config["gatk4"]["extra_args"]
    threads: config["resources"]["gatk4"]["cpus"]
    envmodules:
        "gatk/4.2.6.1-Java-11"
    log:
        "logs/gatk4/{sample}_{coverage}x_{region}.log"
    shell:
        """
        echo "Starting GATK4 HaplotypeCaller" > {log}

        gatk --java-options "-Xmx{params.java_mem}" HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.vcf} \
            --native-pair-hmm-threads {threads} \
            {params.extra} \
            2>> {log}

        echo "GATK4 complete" >> {log}
        """

rule bcftools_mpileup_call:
    input:
        bam = "results/preprocessing/{sample}_{coverage}x_{region}.rg.bam",
        bai = "results/preprocessing/{sample}_{coverage}x_{region}.rg.bam.bai",
        ref = config["reference"]["genome"]
    output:
        vcf = "results/bcftools/{sample}_{coverage}x_{region}.vcf.gz",
        idx = "results/bcftools/{sample}_{coverage}x_{region}.vcf.gz.tbi"
    params:
        max_depth = config["bcftools"]["max_depth"],
        min_MQ = config["bcftools"]["min_MQ"],
        min_BQ = config["bcftools"]["min_BQ"]
    threads: config["resources"]["bcftools"]["cpus"]
    log:
        "logs/bcftools/{sample}_{coverage}x_{region}.log"
    shell:
        """
        echo "Starting bcftools" > {log}

        module load bcftools/1.15.1

        bcftools mpileup \
            --threads {threads} \
            --max-depth {params.max_depth} \
            --min-MQ {params.min_MQ} \
            --min-BQ {params.min_BQ} \
            -Ou -f {input.ref} {input.bam} 2>> {log} \
        | bcftools call \
            --threads {threads} \
            --ploidy 2 \
            --multiallelic-caller \
            --variants-only \
            -Oz -o {output.vcf} 2>> {log}

        bcftools index --threads {threads} --tbi {output.vcf} 2>> {log}

        echo "bcftools complete" >> {log}
        """

rule deepvariant_call:
    input:
        bam = "results/preprocessing/{sample}_{coverage}x_{region}.rg.bam",
        bai = "results/preprocessing/{sample}_{coverage}x_{region}.rg.bam.bai",
        ref = config["reference"]["genome"],
        ref_fai = config["reference"]["fai"]
    output:
        vcf = "results/deepvariant/{sample}_{coverage}x_{region}.vcf.gz",
        idx = "results/deepvariant/{sample}_{coverage}x_{region}.vcf.gz.tbi",
        gvcf = "results/deepvariant/{sample}_{coverage}x_{region}.g.vcf.gz"
    params:
        model_type = config["deepvariant"]["model_type"],
        container = config["deepvariant"]["container"],
        outdir = lambda wildcards: f"results/deepvariant/tmp_{wildcards.sample}_{wildcards.coverage}x_{wildcards.region}"
    threads: config["resources"]["deepvariant"]["cpus"]
    envmodules:
        "singularity/3.8.7"
    log:
        "logs/deepvariant/{sample}_{coverage}x_{region}.log"
    shell:
        """
        echo "Starting DeepVariant" > {log}

        mkdir -p {params.outdir}

        singularity exec \
            --bind $(dirname {input.bam}):/input:ro \
            --bind $(dirname {input.ref}):/ref:ro \
            --bind {params.outdir}:/output \
            {params.container} \
            /opt/deepvariant/bin/run_deepvariant \
            --model_type={params.model_type} \
            --ref=/ref/$(basename {input.ref}) \
            --reads=/input/$(basename {input.bam}) \
            --output_vcf=/output/$(basename {output.vcf}) \
            --output_gvcf=/output/$(basename {output.gvcf}) \
            --num_shards={threads} \
            2>> {log}

        mv {params.outdir}/$(basename {output.vcf}) {output.vcf}
        mv {params.outdir}/$(basename {output.vcf}).tbi {output.idx}
        mv {params.outdir}/$(basename {output.gvcf}) {output.gvcf}

        rm -rf {params.outdir}

        echo "DeepVariant complete" >> {log}
        """
