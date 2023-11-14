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
            category="1. Quality control",
            labels={
              "sample": "all samples",
              "type": "view"
            }
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
        html_bins=report(
            "results/{project}/report/{sample}/bin_summary.html",
            htmlindex="bin_summary.html",
            category="4. Binning results",
            subcategory="4.1 Summary",
            labels={
              "sample": "{sample}",
              "type": "1. view"
            }
        ),
        csv_bins=report(
            "results/{project}/report/{sample}/bin_summary.csv",
            category="4. Binning results",
            subcategory="4.1 Summary",
            labels={
              "sample": "{sample}",
              "type": "2. download"
            }
        ),
        html_checkm=report(
            "results/{project}/report/{sample}/checkm_summary.html",
            htmlindex="checkm_summary.html",
            category="4. Binning results",
            subcategory="4.2 Quality control",
            labels={
              "sample": "{sample}",
              "tool": "CheckM 2",
              "type": "1. view"
            }
        ),
        csv_checkm=report(
            "results/{project}/report/{sample}/checkm_summary.csv",
            category="4. Binning results",
            subcategory="4.2 Quality control",
            labels={
              "sample": "{sample}",
              "tool": "CheckM 2",
              "type": "2. download"
            }
        ),
        html_dastool=report(
            "results/{project}/report/{sample}/DASTool_summary.html",
            htmlindex="DASTool_summary.html",
            category="4. Binning results",
            subcategory="4.2 Quality control",
            labels={
              "sample": "{sample}",
              "tool": "DAS Tool",
              "type": "1. view"
            }
        ),
        csv_dastool=report(
            "results/{project}/report/{sample}/DASTool_summary.csv",
            category="4. Binning results",
            subcategory="4.2 Quality control",
            labels={
              "sample": "{sample}",
              "tool": "DAS Tool",
              "type": "2. download"
            }
        ),
        html_tax=report(
            "results/{project}/report/{sample}/bin_taxonomy.html",
            htmlindex="bin_taxonomy.html",
            category="4. Binning results",
            subcategory="4.3 Taxonomy classification",
            labels={
              "sample": "{sample}",
              "type": "1. view"
            }
        ),
        csv_tax=report(
            "results/{project}/report/{sample}/bin_taxonomy.csv",
            category="4. Binning results",
            subcategory="4.3 Taxonomy classification",
            labels={
              "sample": "{sample}",
              "type": "2. download"
            }
        ),
        html_mags=report(
            "results/{project}/report/{sample}/mags_summary.html",
            htmlindex="mags_summary.html",
            category="5. Taxonomic classification",
            subcategory="5.1 MAGs classification",
            labels={
              "sample": "{sample}",
              "type": "1. view"
            }
        ),
        csv_mags=report(
            "results/{project}/report/{sample}/mags_summary.csv",
            category="5. Taxonomic classification",
            subcategory="5.1 MAGs classification",
            labels={
              "sample": "{sample}",
              "type": "2. download"
            }
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
