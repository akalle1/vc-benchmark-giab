# Preprocessing Rules
# Downsample and region extraction

rule downsample_bam:
    input:
        bam = config["input"]["bam"],
        bai = config["input"]["bam_index"]
    output:
        bam = temp("results/preprocessing/{sample}_{coverage}x_full.bam"),
        bai = temp("results/preprocessing/{sample}_{coverage}x_full.bam.bai")
    params:
        module = config["modules"]["samtools"],
        seed = config["resources"]["downsample"]["seed"],
        frac_digits = config["resources"]["downsample"]["fraction"]
    threads: config["resources"]["downsample"]["cpus"]
    log:
        "logs/preprocessing/downsample_{sample}_{coverage}x.log"
    shell:
        r"""
        module load {params.module}

        echo "Downsampling with samtools -s {params.seed}.{params.frac_digits}" > {log}

        samtools view -@ {threads} -b -s {params.seed}.{params.frac_digits} {input.bam} \
          > {output.bam} 2>> {log}

        samtools index -@ {threads} {output.bam} 2>> {log}
        """

rule extract_region:
    input:
        bam = "results/preprocessing/{sample}_{coverage}x_full.bam",
        bai = "results/preprocessing/{sample}_{coverage}x_full.bam.bai"
    output:
        bam = "results/preprocessing/{sample}_{coverage}x_{region}.bam",
        bai = "results/preprocessing/{sample}_{coverage}x_{region}.bam.bai"
    params:
        module = config["modules"]["samtools"],
        region_bed = get_region_bed
    threads: 4
    log:
        "logs/preprocessing/extract_{sample}_{coverage}x_{region}.log"
    shell:
        r"""
        module load {params.module}
        if [ "{params.region_bed}" = "chr21" ]; then
	    samtools view -b {input.bam} chr21 > {output.bam}
        else
            samtools view -b -L {params.region_bed} {input.bam} > {output.bam}
     	fi
        """

