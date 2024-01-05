include: "host_filtering.smk"


## read QC
rule fastqc:
    input:
        get_trimmed_fastqs,
    output:
        html=temp("results/{project}/qc/fastqc/{sample}_trimmed.html"),
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
            labels={"sample": "all samples"},
        ),
        "results/{project}/qc/multiqc_data.zip",
    params:
        extra=(
            "--zip-data-dir "
            "--config config/multiqc_config.yaml "
            "--title 'Results for data from {project} project'"
        ),
        use_input_files_only=True,
    log:
        "logs/{project}/multiqc.log",
    wrapper:
        "v3.3.1/bio/multiqc"


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


# change to new folder of gz bins
rule checkm2_run:
    input:
        bins="results/{project}/das_tool/{sample}/{sample}_DASTool_bins/",
        dbfile=get_checkm2_db(),  #"{}/{}".format(config["data-handling"]["resources"], config["checkm2"]),
    output:
        stats=temp("results/{project}/qc/checkm2/{sample}/quality_report.tsv"),
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


rule bin_summary_sample:
    input:
        tool="results/{project}/das_tool/{sample}/{sample}_DASTool_summary.tsv",
        checkm="results/{project}/qc/checkm2/{sample}/quality_report.tsv",
        gtdb="results/{project}/classification/{sample}/{sample}.bac120.summary.tsv",
    output:
        csv_bins="results/{project}/output/report/{sample}/{sample}_bin_summary.csv",
        csv_checkm="results/{project}/output/report/{sample}/{sample}_checkm2_summary.csv",
        csv_dastool="results/{project}/output/report/{sample}/{sample}_DASTool_summary.csv",
        csv_tax="results/{project}/output/report/{sample}/{sample}_bin_taxonomy.csv",
        csv_mags="results/{project}/output/report/{sample}/{sample}_mags_summary.csv",
    params:
        max_cont=config["MAG-criteria"]["max-contamination"],  #snakemake.params.contamination
        min_comp=config["MAG-criteria"]["min-completeness"],
    log:
        "logs/{project}/bin_summary/{sample}.log",
    threads: 4
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/bin_summary_sample.py"


use rule kraken2_report as bin_sample_report with:
    input:
        "results/{project}/output/report/{sample}/{sample}_bin_summary.csv",
    output:
        report(
            directory("results/{project}/output/report/{sample}/bin/"),
            htmlindex="index.html",
            category="4. Binning results",
            subcategory="4.1 Summary",
            labels={"sample": "{sample}"},
        ),
    params:
        pin_until="bin",
        styles="resources/report/tables/",
        name="{sample}_bin_summary",
        header="Bin summary for sample {sample}",
    log:
        "logs/{project}/report/{sample}/bin_rbt_csv.log",


use rule kraken2_report as dastool_report with:
    input:
        "results/{project}/output/report/{sample}/{sample}_DASTool_summary.csv",
    output:
        report(
            directory("results/{project}/output/report/{sample}/dastool/"),
            htmlindex="index.html",
            category="4. Binning results",
            subcategory="4.2 Quality control",
            labels={
                "sample": "{sample}",
                "tool": "DAS Tool",
            },
        ),
    params:
        pin_until="bin",
        styles="resources/report/tables/",
        name="{sample}_DASTool_summary",
        header="DAS Tool summary for sample {sample}",
    log:
        "logs/{project}/report/{sample}/dastool_rbt_csv.log",


use rule kraken2_report as checkm2_report with:
    input:
        "results/{project}/output/report/{sample}/{sample}_checkm2_summary.csv",
    output:
        report(
            directory("results/{project}/output/report/{sample}/checkm2/"),
            htmlindex="index.html",
            category="4. Binning results",
            subcategory="4.2 Quality control",
            labels={
                "sample": "{sample}",
                "tool": "CheckM 2",
            },
        ),
    params:
        pin_until="bin",
        styles="resources/report/tables/",
        name="{sample}_CheckM2_summary",
        header="CheckM2 summary for sample {sample}",
    log:
        "logs/{project}/report/{sample}/checkm2_rbt_csv.log",


use rule kraken2_report as taxonomy_report with:
    input:
        "results/{project}/output/report/{sample}/{sample}_bin_taxonomy.csv",
    output:
        report(
            directory("results/{project}/output/report/{sample}/taxonomy/"),
            htmlindex="index.html",
            category="4. Binning results",
            subcategory="4.3 Taxonomy classification",
            labels={"sample": "{sample}"},
        ),
    params:
        pin_until="bin",
        styles="resources/report/tables/",
        name="{sample}_taxonomy_summary",
        header="Taxonomy summary for sample {sample}",
    log:
        "logs/{project}/report/{sample}/taxonomy_rbt_csv.log",


use rule kraken2_report as mag_report with:
    input:
        "results/{project}/output/report/{sample}/{sample}_mags_summary.csv",
    output:
        report(
            directory("results/{project}/output/report/{sample}/mags/"),
            htmlindex="index.html",
            category="5. Taxonomic classification",
            subcategory="5.1 MAGs classification",
            labels={"sample": "{sample}"},
        ),
    params:
        pin_until="MAG",
        styles="resources/report/tables/",
        name="{sample}_MAG_summary",
        header="MAG summary for sample {sample}",
    log:
        "logs/{project}/report/{sample}/mag_rbt_csv.log",


rule bin_summary_all:
    input:
        csv_mags=expand(
            "results/{{project}}/output/report/{sample}/{sample}_mags_summary.csv",
            sample=get_samples(),
        ),
        csv_bins=expand(
            "results/{{project}}/output/report/{sample}/{sample}_bin_summary.csv",
            sample=get_samples(),
        ),
    output:
        "results/{project}/output/report/all/binnning_summary_all.csv",
    log:
        "logs/{project}/bin_summary/all.log",
    threads: 4
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/bin_summary_all.py"


use rule kraken2_report as bin_all_report with:
    input:
        "results/{project}/output/report/all/binnning_summary_all.csv",
    output:
        report(
            directory("results/{project}/output/report/all/binning/"),
            htmlindex="index.html",
            category="4. Binning results",
            subcategory="4.1 Summary",
            labels={"sample": "all"},
        ),
    params:
        pin_until="sample",
        styles="resources/report/tables/",
        name="bin_summary",
        header="Bin summary for all samples",
    log:
        "logs/{project}/report/all_bin_rbt_csv.log",
