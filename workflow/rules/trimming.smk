# fastp in paired-end mode for Illumina paired-end data
rule fastp:
    input:
        sample=get_fastqs,
    output:
        trimmed=
            temp([
                "results/{project}/trimmed/fastp/{sample}.1.fastq.gz",
                "results/{project}/trimmed/fastp/{sample}.2.fastq.gz",
            ]),
        html="results/{project}/trimmed/fastp/{sample}.html",
        json="results/{project}/trimmed/fastp/{sample}.fastp.json",
    params:
        adapters=get_adapters,
        extra="--qualified_quality_phred {phred} --length_required {minlen}".format(
            phred=(config["quality-criteria"]["min-PHRED"]),
            minlen=(config["quality-criteria"]["min-length-reads"])
        ),
    log:
        "logs/{project}/fastp/{sample}.log",
    threads: 2
    wrapper:
        "v1.23.5/bio/fastp"