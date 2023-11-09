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
        dbfile=get_checkm2_db(),  #"{}/{}".format(config["data-handling"]["resources"], config["checkm2"]),
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
        dbfile=get_checkm2_db(),  #"{}/{}".format(config["data-handling"]["resources"], config["checkm2"]),
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
        "checkm2 predict -x fa --threads {threads} --force "
        "--input {input.bins}/ --output-directory {params.outdir}/ "
        "> {log} 2>&1"
        # --remove_intermediates-x fa.gz


rule bin_summary:
    input:
        tool="results/{project}/das_tool/{sample}/{sample}_DASTool_summary.tsv",
        checkm="results/{project}/qc/checkm2/{sample}/quality_report.tsv",
        gtdb="results/{project}/classification/{sample}/{sample}.bac120.summary.tsv",
    output:
        csv_mags="results/{project}/report/{sample}/mags_summary.csv",
        csv_bins="results/{project}/report/{sample}/bin_summary.csv",
        csv_tax="results/{project}/report/{sample}/bin_taxonomy.csv",
        csv_checkm="results/{project}/report/{sample}/checkm_summary.csv",
        csv_dastool="results/{project}/report/{sample}/DASTool_summary.csv",
        html_mags=report(
            "results/{project}/report/{sample}/mags_summary.html",
            htmlindex="mags_summary.html",
            category="MAGs",
        ),
        html_bins=report(
            "results/{project}/report/{sample}/bin_summary.html",
            htmlindex="bin_summary.html",
            category="Bins",
        ),
        html_tax=report(
            "results/{project}/report/{sample}/bin_taxonomy.html",
            htmlindex="bin_taxonomy.html",
            category="Bins",
        ),
        html_checkm=report(
            "results/{project}/report/{sample}/checkm_summary.html",
            htmlindex="checkm_summary.html",
            category="Bins",
        ),
        html_dastool=report(
            "results/{project}/report/{sample}/DASTool_summary.html",
            htmlindex="DASTool_summary.html",
            category="Bins",
        ),
    params:
        max_cont=config["MAG-criteria"]["max-contamination"],  #snakemake.params.contamination
        min_comp=config["MAG-criteria"]["min-completeness"],
    log:
        "logs/{project}/bin_summary/{sample}.log",
    threads: 4
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/bin_summary.py"
