RAW_DATA_PATH = get_data_path()


# copy files to local
### TODO: change to copying via tar or remove & expect local path 
'''rule copy_fastq:
    input:
        fastqs=get_fastqs,
    output:
        raw1=f"{RAW_DATA_PATH}{{project}}/{{sample}}_R1.fastq.gz",
        raw2=f"{RAW_DATA_PATH}{{project}}/{{sample}}_R2.fastq.gz",
    params:
        outdir=lambda wildcards, output: Path(output.raw1).parent,
    log:
        "logs/{project}/copy_data/{sample}.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "(mkdir -p {params.outdir} && "
        "cp -v {input.fastqs[0]} {output.raw1} && "
        "cp -v {input.fastqs[1]} {output.raw2}) > {log} 2>&1"'''


# fastp in paired-end mode for Illumina paired-end data
rule fastp:
    input:
        sample=get_local_fastqs,
    output:
        trimmed=temp(
            [
                "results/{project}/trimmed/fastp/{sample}.1.fastq.gz",
                "results/{project}/trimmed/fastp/{sample}.2.fastq.gz",
            ]
        ),
        html=temp("results/{project}/trimmed/fastp/{sample}.html"),
        json="results/{project}/trimmed/fastp/{sample}.fastp.json",
    params:
        adapters=get_adapters,
        extra="--qualified_quality_phred {phred} --length_required {minlen}".format(
            phred=(config["quality-criteria"]["min-PHRED"]),
            minlen=(config["quality-criteria"]["min-length-reads"]),
        ),
    log:
        "logs/{project}/fastp/{sample}.log",
    threads: 2
    wrapper:
        "v2.6.0/bio/fastp"
