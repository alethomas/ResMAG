from pathlib import Path


rule download_kraken_db:
    output:
        hfile=get_kraken_db_file(),
    params:
        download=config["kraken"]["download-path"],
        db_folder=lambda wildcards, output: Path(output.hfile).parent,
    log:
        "logs/kraken2_DB_download.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "(mkdir -p {params.db_folder} && "
        "wget -c {params.download} -O - | "
        "tar -zxv -C {params.db_folder}) > {log} 2>&1"


if not config["testing"]:

    rule kraken2:
        input:
            hfile=get_kraken_db_file(),
            fastqs=get_trimmed_fastqs,
        output:
            clf1=temp("results/{project}/filtered/{sample}_clf_1.fastq"),
            clf2=temp("results/{project}/filtered/{sample}_clf_2.fastq"),
            unclf1=temp("results/{project}/filtered/{sample}_unclf_1.fastq"),
            unclf2=temp("results/{project}/filtered/{sample}_unclf_2.fastq"),
            report="results/{project}/filtered/kraken/{sample}_report.tsv",
            outfile=temp("results/{project}/filtered/kraken/{sample}_outfile.tsv"),
        params:
            clf=lambda wildcards, output: output.clf1.replace("_clf_1", "_clf#"),
            unclf=lambda wildcards, output: output.unclf1.replace("_unclf_1", "_unclf#"),
            db=lambda wildcards, input: Path(input.hfile).parent,
        threads: 32
        log:
            "logs/{project}/kraken2/run/{sample}.log",
        conda:
            "../envs/kraken2.yaml"
        shell:
            "kraken2 --db {params.db} --threads {threads} --paired --classified-out {params.clf} "
            "--unclassified-out {params.unclf} --output {output.outfile} "
            "--report {output.report} --gzip-compressed {input.fastqs} > {log} 2>&1"


if config["testing"]:

    rule gunzip_kraken_output:
        input:
            clf_gz_1="results/{project}/filtered/{sample}_clf_1.fastq.gz",
            clf_gz_2="results/{project}/filtered/{sample}_clf_2.fastq.gz",
        output:
            clf1=temp("results/{project}/filtered/{sample}_clf_1.fastq"),
            clf2=temp("results/{project}/filtered/{sample}_clf_2.fastq"),
        shell:
            "gunzip -k {input.clf_gz_1} {input.clf_gz_1}"


rule remove_human_reads:
    input:
        clf1="results/{project}/filtered/{sample}_clf_1.fastq",
        clf2="results/{project}/filtered/{sample}_clf_2.fastq",
        report="results/{project}/filtered/kraken/{sample}_report.tsv",
        outfile="results/{project}/filtered/kraken/{sample}_outfile.tsv",
    output:
        out1=temp("results/{project}/filtered/{sample}_non_human_1.fastq"),
        out2=temp("results/{project}/filtered/{sample}_non_human_2.fastq"),
    log:
        "logs/{project}/kraken2/remove_human/{sample}.log",
    params:
        human_tax=9606,
    threads: 4
    conda:
        "../envs/kraken2.yaml"
    shell:
        "extract_kraken_reads.py -s1 {input.clf1} -s2 {input.clf2} -k {input.outfile} "
        "-r {input.report} -t {params.human_tax}  --exclude "
        "-o {output.out1} -o2 {output.out2} --fastq-output > {log} 2>&1"


## Add possibility to also remove other hosts eg pig


rule add_unclf_reads:
    input:
        filt="results/{project}/filtered/{sample}_non_human_{read}.fastq",
        unclf="results/{project}/filtered/{sample}_unclf_{read}.fastq",
    output:
        temp("results/{project}/filtered/fastqs/{sample}_{read}.fastq"),
    log:
        "logs/{project}/kraken2/add_unclf/{sample}_{read}.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "cat {input.filt} {input.unclf} > {output} 2> {log}"


rule gzip_kraken_output:
    input:
        "results/{project}/filtered/fastqs/{sample}_{read}.fastq",
    output:
        "results/{project}/filtered/fastqs/{sample}_{read}.fastq.gz",
    log:
        "logs/{project}/kraken2/gzip_fastq/{sample}_{read}.log",
    threads: 4
    conda:
        "../envs/unix.yaml"
    shell:
        "gzip -k {input} > {log} 2>&1"


rule kraken_summary:
    input:
        reports=expand(
            "results/{{project}}/filtered/kraken/{sample}_report.tsv",
            sample=get_samples(),
        ),
    output:
        csv="results/{project}/report/kraken2_summary.csv",
        html=report(
            "results/{project}/report/kraken2_summary.html",
            category="1. Species diversity",
            subcategory="before host read filtering",
            caption="../report/kraken.rst",
        ),
    log:
        "logs/{project}/kraken2/summary.log",
    params:
        taxid_dict=get_taxID_dict(),
    threads: 2
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/kraken_summary.py"


rule kraken2_postfilt:
    input:
        hfile=get_kraken_db_file(),
        fastqs=get_filtered_gz_fastqs,
    output:
        report="results/{project}/filtered/kraken_postfilt/{sample}_report.tsv",
        outfile=temp("results/{project}/filtered/kraken_postfilt/{sample}_outfile.tsv"),
    params:
        db=lambda wildcards, input: Path(input.hfile).parent,
    threads: 32
    log:
        "logs/{project}/kraken2/run_postfilt/{sample}.log",
    conda:
        "../envs/kraken2.yaml"
    shell:
        "kraken2 --db {params.db} --threads {threads} --paired "
        "--output {output.outfile} --report {output.report} "
        "--gzip-compressed {input.fastqs} > {log} 2>&1"


rule bracken_analysis:
    input:
        hfile=get_kraken_db_file(),
        kreport="results/{project}/filtered/kraken_postfilt/{sample}_report.tsv",
    output:
        breport="results/{project}/filtered/bracken/reports/{sample}.breport",
        bfile="results/{project}/filtered/bracken/files/{sample}.bracken",
    params:
        db=lambda wildcards, input: Path(input.hfile).parent,
    log:
        "logs/{project}/bracken/{sample}.log",
    resources:
        mem_mb=100,
    threads: 1
    params:
        threads=1,
    conda:
        "../envs/bracken.yaml"
    shell:
        "bracken -d {params.db} -i {input.kreport} -l F -o {output.bfile} -w {output.breport} > {log} 2>&1"


rule merge_bracken:
    input:
        expand(
            "results/{{project}}/filtered/bracken/files/{sample}.bracken",
            sample=get_samples(),
        ),
    output:
        "results/{project}/filtered/bracken/merged.bracken.txt",
    log:
        "logs/{project}/bracken/merge_bracken.log",
    resources:
        mem_mb=100,
    threads: 1
    params:
        threads=1,
    conda:
        "../envs/bracken.yaml"
    shell:
        "(python $CONDA_PREFIX/bin/combine_bracken_outputs.py --files {input} --output {output}) > {log} 2>&1"


rule create_bracken_plot:
    input:
        "results/{project}/filtered/bracken/merged.bracken.txt",
    output:
        report(
            "results/{project}/report/bracken_plot.png",
            category="1. Species diversity",
            subcategory="after host read filtering",
        ),
    params:
        threshold=0.001,
    log:
        "logs/{project}/report/bracken_plot.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/brackenplot.py"


rule kraken2krona:
    input:
        "results/{project}/filtered/kraken_postfilt/{sample}_report.tsv",
    output:
        temp("results/{project}/filtered/krona/{sample}.krona"),
    threads: 4
    log:
        "logs/{project}/kraken2/postfilt_krona/{sample}.log",
    conda:
        "../envs/kraken2.yaml"
    shell:
        "(python $CONDA_PREFIX/bin/kreport2krona.py -r {input} "
        "-o {output}) > {log} 2>&1"


rule krona_html:
    input:
        "results/{project}/filtered/krona/{sample}.krona",
    output:
        report(
            "results/{project}/report/{sample}/kraken.krona.html",
            category="1. Species diversity",
            subcategory="after host read filtering",
            caption="../report/kraken.rst",
        ),
    threads: 1
    log:
        "logs/{project}/krona/{sample}.log",
    conda:
        "../envs/kaiju.yaml"
    shell:
        "(ktImportText -o {output} {input}) > {log} 2>&1"


#  results/autobrewer/filtered/krona/ABS_24_07.krona.html results/autobrewer/filtered/krona/ABS_24_07.krona
