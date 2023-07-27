RAW_DATA_PATH = get_data_path()


# copy files to local
rule copy_fastp:
    input:
        get_fastqs,
    output:
        raw1="{data}{{project}}/{{sample}}_R1.fastq.gz".format(data=RAW_DATA_PATH),
        raw2="{data}{{project}}/{{sample}}_R2.fastq.gz".format(data=RAW_DATA_PATH),
    params:
        outdir=lambda wildcards, output: Path(output.raw1).parent,
    log:
        "logs/{project}/copy_data/{sample}.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "(mkdir -p {params.outdir} && "
        "cp -v -t {params.outdir} {input}) > {log} 2>&1"


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
        html="results/{project}/trimmed/fastp/{sample}.html",
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
        "v1.23.5/bio/fastp"
