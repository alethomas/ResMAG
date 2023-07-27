from pathlib import Path


## copy DB to local
rule copy_kraken_db:
    input:
        get_kraken_db(),
    output:
        directory("{}krakenDB/".format(get_resource_path())),
    threads: 2
    log:
        "logs/copy_kraken_db.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "cp -vr {input} {output} > {log} 2>&1"


if not config["testing"]:

    rule kraken2:
        input:
            fastqs=get_trimmed_fastqs,
            db=get_local_krakenDB(),
        output:
            clf1=temp("results/{project}/filtered/{sample}_clf_1.fastq"),
            clf2=temp("results/{project}/filtered/{sample}_clf_2.fastq"),
            unclf1=temp("results/{project}/filtered/{sample}_unclf_1.fastq"),
            unclf2=temp("results/{project}/filtered/{sample}_unclf_2.fastq"),
            report="results/{project}/filtered/kraken/{sample}_report.tsv",
            outfile="results/{project}/filtered/kraken/{sample}_outfile.tsv",
        params:
            clf=lambda wildcards, output: output.clf1.replace("_clf_1", "_clf#"),
            unclf=lambda wildcards, output: output.unclf1.replace("_unclf_1", "_unclf#"),
        threads: 8
        log:
            "logs/{project}/kraken2/{sample}.log",
        conda:
            "../envs/kraken2.yaml"
        shell:
            "kraken2 --db {input.db} --threads {threads} --paired --classified-out {params.clf} "
            "--unclassified-out {params.unclf} --output {output.outfile} "
            "--report {output.report} --gzip-compressed {input.fastqs} > {log} 2>&1"


if config["testing"]:

    rule gzip_kraken_output:
        input:
            clf_gz_1="results/{project}/filtered/{sample}_clf_1.fastq.gz",
            clf_gz_2="results/{project}/filtered/{sample}_clf_2.fastq.gz",
        output:
            clf1=temp("results/{project}/filtered/{sample}_clf_1.fastq"),
            clf2=temp("results/{project}/filtered/{sample}_clf_2.fastq"),
        shell:
            "gunzip -k {input.clf_gz_1} {input.clf_gz_1}"


rule extract_kraken_reads:
    input:
        clf1="results/{project}/filtered/{sample}_clf_1.fastq",
        clf2="results/{project}/filtered/{sample}_clf_2.fastq",
        report="results/{project}/filtered/kraken/{sample}_report.tsv",
        outfile="results/{project}/filtered/kraken/{sample}_outfile.tsv",
    output:
        out1=temp(
            "results/{project}/filtered/{kraken_ref}/{sample}_{kraken_ref}_1.fastq"
        ),
        out2=temp(
            "results/{project}/filtered/{kraken_ref}/{sample}_{kraken_ref}_2.fastq"
        ),
    log:
        "logs/{project}/kraken2_extract_reads/{sample}_{kraken_ref}.log",
    params:
        taxid=get_taxID,
    threads: 2
    conda:
        "../envs/kraken2.yaml"
    shell:
        "extract_kraken_reads.py -s1 {input.clf1} -s2 {input.clf2} -k {input.outfile} "
        "-r {input.report} -t {params.taxid} --include-children "
        "-o {output.out1} -o2 {output.out2} --fastq-output > {log} 2>&1"


rule add_unclf_reads:
    input:
        clf="results/{project}/filtered/{kraken_ref}/{sample}_{kraken_ref}_{read}.fastq",
        unclf="results/{project}/filtered/{sample}_unclf_{read}.fastq",
    output:
        out="results/{project}/filtered/{kraken_ref}/{sample}_{kraken_ref}_all_{read}.fastq.gz",
    log:
        "logs/{project}/kraken2_add_unclf_reads/{sample}_{kraken_ref}_{read}.log",
    params:
        taxid=get_taxID,
        inter=lambda wildcards, output: output.out.split(".gz")[0],
    conda:
        "../envs/unix.yaml"
    shell:
        "(cat {input.clf} {input.unclf} > {params.inter} "
        "&& gzip {params.inter}) > {log} 2>&1"


rule kraken_summary:
    input:
        report="results/{project}/filtered/kraken/{sample}_report.tsv",
        json="results/{project}/trimmed/fastp/{sample}.fastp.json",
    output:
        ensure("results/{project}/filtered/{sample}_summary.csv", non_empty=True),
    log:
        "logs/{project}/kraken2/{sample}_summary.log",
    params:
        taxid_dict=get_taxID_dict(),
    threads: 2
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/kraken_summary.py"
