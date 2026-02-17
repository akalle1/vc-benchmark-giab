# ==============================================================================
# GATK4 HAPLOTYPECALLER
# ==============================================================================

rule gatk4_haplotypecaller:
    """Call variants using GATK4 HaplotypeCaller"""
    input:
        bam = "results/preprocessing/{sample}_{coverage}x_{region}.bam",
        bai = "results/preprocessing/{sample}_{coverage}x_{region}.bam.bai",
        ref = config["reference"]["genome"],
        ref_dict = config["reference"]["dict"]
    output:
        vcf = "results/gatk4/{sample}_{coverage}x_{region}.vcf.gz",
        idx = "results/gatk4/{sample}_{coverage}x_{region}.vcf.gz.tbi"
    params:
        gatk_module = config["modules"]["gatk"],
        java_mem = config["gatk4"]["java_mem"],
        extra = config["gatk4"]["extra_args"]
    threads: config["resources"]["gatk4"]["cpus"]
    resources:
        mem_mb = config["resources"]["gatk4"]["mem_mb"],
        time_min = config["resources"]["gatk4"]["time_min"],
        partition = config["resources"]["default"]["partition"],
        account = config["resources"]["default"]["account"]
    log:
        "logs/gatk4/{sample}_{coverage}x_{region}.log"
    benchmark:
        "logs/gatk4/{sample}_{coverage}x_{region}.benchmark.txt"
    shell:
        """
        module load {params.gatk_module}
        
        echo "Starting GATK4 HaplotypeCaller" > {log}
        echo "Sample: {wildcards.sample}" >> {log}
        echo "Coverage: {wildcards.coverage}x" >> {log}
        echo "Region: {wildcards.region}" >> {log}
        echo "Threads: {threads}" >> {log}
        echo "Memory: {params.java_mem}" >> {log}
        
        gatk --java-options "-Xmx{params.java_mem}" HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.vcf} \
            --native-pair-hmm-threads {threads} \
            {params.extra} \
            2>> {log}
        
        echo "GATK4 HaplotypeCaller complete" >> {log}
        """

# ==============================================================================
# BCFTOOLS MPILEUP + CALL
# ==============================================================================

rule bcftools_mpileup_call:
    """Call variants using bcftools mpileup + call"""
    input:
        bam = "results/preprocessing/{sample}_{coverage}x_{region}.bam",
        bai = "results/preprocessing/{sample}_{coverage}x_{region}.bam.bai",
        ref = config["reference"]["genome"]
    output:
        vcf = "results/bcftools/{sample}_{coverage}x_{region}.vcf.gz",
        idx = "results/bcftools/{sample}_{coverage}x_{region}.vcf.gz.tbi"
    params:
        bcftools_module = config["modules"]["bcftools"],
        max_depth = config["bcftools"]["max_depth"],
        min_MQ = config["bcftools"]["min_MQ"],
        min_BQ = config["bcftools"]["min_BQ"]
    threads: config["resources"]["bcftools"]["cpus"]
    resources:
        mem_mb = config["resources"]["bcftools"]["mem_mb"],
        time_min = config["resources"]["bcftools"]["time_min"],
        partition = config["resources"]["default"]["partition"],
        account = config["resources"]["default"]["account"]
    log:
        "logs/bcftools/{sample}_{coverage}x_{region}.log"
    benchmark:
        "logs/bcftools/{sample}_{coverage}x_{region}.benchmark.txt"
    shell:
        """
        module load {params.bcftools_module}
        
        echo "Starting bcftools mpileup + call" > {log}
        echo "Sample: {wildcards.sample}" >> {log}
        echo "Coverage: {wildcards.coverage}x" >> {log}
        echo "Region: {wildcards.region}" >> {log}
        echo "Max depth: {params.max_depth}" >> {log}
        echo "Min MQ: {params.min_MQ}" >> {log}
        echo "Min BQ: {params.min_BQ}" >> {log}
        
        # mpileup generates likelihoods, call makes the actual calls
        bcftools mpileup \
            --threads {threads} \
            --max-depth {params.max_depth} \
            --min-MQ {params.min_MQ} \
            --min-BQ {params.min_BQ} \
            -Ou \
            -f {input.ref} \
            {input.bam} \
            2>> {log} \
        | bcftools call \
            --threads {threads} \
            --ploidy 2 \
            --multiallelic-caller \
            --variants-only \
            -Oz \
            -o {output.vcf} \
            2>> {log}
        
        # Index the output
        bcftools index --threads {threads} --tbi {output.vcf} 2>> {log}
        
        echo "bcftools variant calling complete" >> {log}
        """

# ==============================================================================
# DEEPVARIANT (via Singularity)
# ==============================================================================

rule deepvariant_call:
    """Call variants using DeepVariant"""
    input:
        bam = "results/preprocessing/{sample}_{coverage}x_{region}.bam",
        bai = "results/preprocessing/{sample}_{coverage}x_{region}.bam.bai",
        ref = config["reference"]["genome"],
        ref_fai = config["reference"]["fai"]
    output:
        vcf = "results/deepvariant/{sample}_{coverage}x_{region}.vcf.gz",
        idx = "results/deepvariant/{sample}_{coverage}x_{region}.vcf.gz.tbi",
        gvcf = "results/deepvariant/{sample}_{coverage}x_{region}.g.vcf.gz"
    params:
        singularity_module = config["modules"]["singularity"],
        model_type = config["deepvariant"]["model_type"],
        container = config["deepvariant"]["container"],
        # DeepVariant needs output directory for intermediate files
        outdir = lambda wildcards: f"results/deepvariant/tmp_{wildcards.sample}_{wildcards.coverage}x_{wildcards.region}"
    threads: config["resources"]["deepvariant"]["cpus"]
    resources:
        mem_mb = config["resources"]["deepvariant"]["mem_mb"],
        time_min = config["resources"]["deepvariant"]["time_min"],
        partition = config["resources"]["default"]["partition"],
        account = config["resources"]["default"]["account"]
    log:
        "logs/deepvariant/{sample}_{coverage}x_{region}.log"
    benchmark:
        "logs/deepvariant/{sample}_{coverage}x_{region}.benchmark.txt"
    shell:
        """
        module load {params.singularity_module}
        
        echo "Starting DeepVariant" > {log}
        echo "Sample: {wildcards.sample}" >> {log}
        echo "Coverage: {wildcards.coverage}x" >> {log}
        echo "Region: {wildcards.region}" >> {log}
        echo "Model: {params.model_type}" >> {log}
        echo "Threads: {threads}" >> {log}
        
        # Create output directory
        mkdir -p {params.outdir}
        
        # Run DeepVariant via Singularity
        # Note: DeepVariant paths inside container are /input, /output, /ref
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
        
        # Move outputs to final location
        mv {params.outdir}/$(basename {output.vcf}) {output.vcf}
        mv {params.outdir}/$(basename {output.vcf}).tbi {output.idx}
        mv {params.outdir}/$(basename {output.gvcf}) {output.gvcf}
        
        # Clean up intermediate files
        rm -rf {params.outdir}
        
        echo "DeepVariant complete" >> {log}
        """

