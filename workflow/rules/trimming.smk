RAW_DATA_PATH = get_data_path()

"""
# copy files to local
rule copy_fastq:
    output:
        raw1=f"{RAW_DATA_PATH}{{project}}/{{sample}}_R1.fastq.gz",
        raw2=f"{RAW_DATA_PATH}{{project}}/{{sample}}_R2.fastq.gz",
    params:
        outdir=lambda wildcards, output: Path(output.raw1).parent,
        # returns a list of folder and the filenames for R1 and R2 reads
        fastqs=get_fastqs,
    threads: 20
    log:
        "logs/{project}/copy_data/{sample}.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "(mkdir -p {params.outdir} && "
        "tar cpfz - -C {params.fastqs[0]}/ {params.fastqs[1]} {params.fastqs[2]} | "
        "(cd {params.outdir} ; tar xpfz -)) > {log} 2>&1"
"""


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
    threads: 30
    wrapper:
        "v2.6.0/bio/fastp"
