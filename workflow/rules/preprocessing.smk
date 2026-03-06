# Preprocessing Rules
 
rule downsample_bam:
    input:
        bam = config["input"]["bam"],
        bai = config["input"]["bam_index"]
    output:
        bam = temp("results/preprocessing/{sample}_{coverage}x_full.bam"),
        bai = temp("results/preprocessing/{sample}_{coverage}x_full.bam.bai")
    params:
        sampling = lambda wc: f"{config['resources']['downsample']['seed']}.{str(config['resources']['downsample']['fraction']).split('.')[1]}"
    threads: config["resources"]["downsample"]["cpus"]
    envmodules:
        "samtools/1.15.1"
    log:
        "logs/preprocessing/downsample_{sample}_{coverage}x.log"
    shell:
        r"""
        echo "Downsampling with samtools -s {params.sampling}" > {log}

        samtools view -@ {threads} -b -s {params.sampling} {input.bam} \
          > {output.bam} 2>> {log}

        samtools index -@ {threads} {output.bam} 2>> {log}

        echo "Downsampling complete" >> {log}
        """

rule extract_region:
    input:
        bam = "results/preprocessing/{sample}_{coverage}x_full.bam",
        bai = "results/preprocessing/{sample}_{coverage}x_full.bam.bai"
    output:
        bam = "results/preprocessing/{sample}_{coverage}x_{region}.bam",
        bai = "results/preprocessing/{sample}_{coverage}x_{region}.bam.bai"
    params:
        region_bed = get_region_bed
    threads: 4
    envmodules:
        "samtools/1.15.1"
    log:
        "logs/preprocessing/extract_{sample}_{coverage}x_{region}.log"
    shell:
        r"""
        echo "Extracting region: {wildcards.region}" > {log}
        echo "Using: {params.region_bed}" >> {log}

        if [ "{params.region_bed}" = "chr21" ]; then
            echo "Using chromosome name directly" >> {log}
            samtools view -b -@ {threads} {input.bam} chr21 > {output.bam} 2>> {log}
        else
            echo "Using BED file" >> {log}
            samtools view -b -@ {threads} -L {params.region_bed} {input.bam} > {output.bam} 2>> {log}
        fi

        samtools index -@ {threads} {output.bam} 2>> {log}
        
        echo "Extraction complete" >> {log}
        """

    
rule add_readgroup:
    input:
        bam = "results/preprocessing/{sample}_{coverage}x_{region}.bam",
        bai = "results/preprocessing/{sample}_{coverage}x_{region}.bam.bai"
    output:
        bam = "results/preprocessing/{sample}_{coverage}x_{region}.rg.bam",
        bai = "results/preprocessing/{sample}_{coverage}x_{region}.rg.bam.bai"
    threads: 4
    log:
	"logs/preprocessing/addRG_{sample}_{coverage}x_{region}.log"
    shell:
        r"""
        samtools addreplacerg -w \
          -r '@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}\tPL:ILLUMINA' \
          -@ {threads} \
          -o {output.bam} \
          {input.bam} 2> {log}

        samtools index -@ {threads} {output.bam} 2>> {log}
        """
