## read QC
rule fastqc:
    input:
        get_trimmed_fastqs,
    output:
        html="results/{project}/qc/fastqc/{sample}_trimmed.html",
        zip="results/{project}/qc/fastqc/{sample}_trimmed_fastqc.zip",
    log:
        "logs/{project}/fastqc/{sample}.log",
    wrapper:
        "v1.23.5/bio/fastqc"


rule multiqc:
    input:
        expand(
            [
                "results/{{project}}/qc/fastqc/{sample}_trimmed_fastqc.zip",
                "results/{{project}}/trimmed/fastp/{sample}.fastp.json",
            ],
            sample=get_samples(),
        ),
    output:
        report(
            "results/{project}/qc/multiqc.html",
            htmlindex="multiqc.html",
            category="3. Sequencing Details",
        ),
    params:
        extra=(
            "--config config/multiqc_config.yaml --title 'Results for data from {project} project'"
        ),
    log:
        "logs/{project}/multiqc.log",
    wrapper:
        "v1.23.5/bio/multiqc"


## bin QC
rule checkm2_DB_download:
    output:
        dbfile=str(config["checkm2"]),
    params:
        direct=lambda wildcards, output: Path(output.dbfile).parent.parent,
    log:
        "logs/checkm2_DB_download.log",
    conda:
        "../envs/checkm2.yaml"
    shell:
        "checkm2 database --download --path {params.direct} > {log} 2>&1"


rule checkm2_run:
    input:
        bins="results/{project}/das_tool/{sample}/{sample}_DASTool_bins/",
        dbfile=str(config["checkm2"]),
    output:
        stats="results/{project}/qc/checkm2/{sample}/quality_report.tsv",
    params:
        outdir=lambda wildcards, output: Path(output.stats).parent,
    log:
        "logs/{project}/checkm2/{sample}.log",
    threads: 4
    conda:
        "../envs/checkm2.yaml"
    shell:
        "checkm2 predict -x fa --threads {threads} "
        "--input {input.bins}/ --output-directory {params.outdir}/ "
        "> {log} 2>&1"
        # --remove_intermediates-x fa.gz
