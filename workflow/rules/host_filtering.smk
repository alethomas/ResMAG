from pathlib import Path

REFERENCE_DIR = get_reference_dir()
KRAKEN_DB = get_kraken_db()

'''## indexing
rule minimap2_index:
    input:
        target = get_reference_file
    output:
        "{0}{{ref}}.mmi".format(REFERENCE_DIR)
    log:
        "logs/minimap2/{ref}.log"
    params:
        extra=""  # optional additional args
    threads: 3
    wrapper:
        "v1.25.0/bio/minimap2/index"

## minimap2 alignment & filtering
rule filter_host_reads:
    input:
        ref_genomes = get_reference_index,
        query = get_trimmed_fastqs,
    output:
        filt_fastq1 = 'results/{project}/filtered/{sample}_final.1.fastq.gz',
        filt_fastq2 = 'results/{project}/filtered/{sample}_final.2.fastq.gz',
    log:
        "logs/{project}/host_filtering/{sample}.log",
    params:
        filt_path = lambda wildcards, output: Path(output.filt_fastq1).parent,
    threads: 3
    conda:
        "../envs/minimap2.yaml"
    script:
        "../scripts/host_filtering.py"'''

if not config["testing"]:

    rule kraken2:
        input:
            get_trimmed_fastqs,
        output:
            clf1=temp("results/{project}/filtered/{sample}_clf_1.fastq"),
            clf2=temp("results/{project}/filtered/{sample}_clf_2.fastq"),
            report="results/{project}/filtered/kraken/{sample}_report.tsv",
            outfile="results/{project}/filtered/kraken/{sample}_outfile.tsv",
        params:
            clf="results/{project}/filtered/{sample}_clf#.fastq",
            db=KRAKEN_DB,
        threads: 3
        log:
            "logs/{project}/kraken2/{sample}.log",
        conda:
            "../envs/kraken2.yaml"
        shell:
            "kraken2 --db {params.db} --threads {threads} --paired "
            "--classified-out {params.clf} --output {output.outfile} "
            "--report {output.report} --gzip-compressed {input} 2>> {log}"


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
        out1="results/{project}/filtered/{kraken_ref}/{sample}_{kraken_ref}_1.fastq.gz",
        out2="results/{project}/filtered/{kraken_ref}/{sample}_{kraken_ref}_2.fastq.gz",
    log:
        "logs/{project}/extract_kraken_reads/{sample}_{kraken_ref}.log",
    params:
        taxid=get_taxID,
        int1=lambda wildcards, output: output.out1.split(".gz")[0],
        int2=lambda wildcards, output: output.out2.split(".gz")[0],
    threads: 2
    conda:
        "../envs/kraken2.yaml"
    shell:
        "extract_kraken_reads.py -s1 {input.clf1} -s2 {input.clf2} -k {input.outfile} "
        "-r {input.report} -t {params.taxid} --include-children "
        "-o {params.int1} -o2 {params.int2} --fastq-output && "
        "gzip {params.int1} {params.int2} 2>> {log}"


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
        "../envs/kraken2.yaml"
    script:
        "../scripts/kraken_summary.py"
