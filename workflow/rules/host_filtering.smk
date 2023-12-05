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
            fastqs=get_prefiltered_fastqs,
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
        human_tax=get_human_tax_ID(),
    threads: 4
    conda:
        "../envs/kraken2.yaml"
    shell:
        "extract_kraken_reads.py -s1 {input.clf1} -s2 {input.clf2} -k {input.outfile} "
        "-r {input.report} -t {params.human_tax}  --exclude "
        "-o {output.out1} -o2 {output.out2} --fastq-output > {log} 2>&1"


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
        jsons=expand(
            "results/{{project}}/trimmed/fastp/{sample}.fastp.json",
            sample=get_samples(),
        ),
    output:
        csv="results/{project}/report/kraken2_summary.csv",
    log:
        "logs/{project}/kraken2/summary.log",
    params:
        taxid_dict=get_taxID_dict(),
    threads: 2
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/kraken_summary.py"


rule kraken2_report:
    input:
        "results/{project}/report/kraken2_summary.csv",
    output:
        report(
            directory("results/{project}/report/kraken2/"),
            htmlindex="index.html",
            caption="../report/kraken.rst",
            category="2. Species diversity",
            subcategory="2.1 pre-filtering human reads",
            labels={
                "sample": "all samples",
            },
        ),
    params:
        pin_until="sample",
        styles="resources/report/tables/",
        name="kraken2_summary",
        header="Kraken2 summary",
    log:
        "logs/{project}/report/kraken2_rbt_csv.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt csv-report {input} --pin-until {params.pin_until} {output} && "
        "(sed -i '/>github<\/a>/a \\\\t\\t\\t</li>\\n\\t\\t\\t<li class=\"nav-item\">"
        "\\n\\t\\t\\t\\t<a class=\"nav-link\" href=\"#\">{params.header}</a>' "
        "{output}/indexes/index1.html && "
        "sed -i 's/report.xlsx/{params.name}_report.xlsx/g' {output}/indexes/index1.html) && "
        "mv {output}/report.xlsx {output}/{params.name}_report.xlsx && "
        "cp {params.styles}* {output}/css/ > {log} 2>&1"


if config["host_filtering"]["do_host_filtering"]:
    rule map_to_host:
        input:
            fastqs=get_trimmed_fastqs,
        output:
            temp("results/{project}/host_filtering/alignments/{sample}.bam"),
        params:
            ref=config["host_filtering"]["ref_genome"],
        threads: 10
        log:
            "logs/{project}/host_filtering/map_to_host_{sample}.log",
        conda:
            "../envs/minimap2.yaml"
        shell:
            "(minimap2 -a -xsr -t {threads} {params.ref} {input.fastqs} | "
            "samtools view -bh | "
            "samtools sort --threads {threads} -o {output}) > {log} 2>&1"


    rule index_host_alignment:
        input:
            "results/{project}/host_filtering/alignments/{sample}.bam",
        output:
            temp("results/{project}/host_filtering/alignments/{sample}.bam.bai"),
        threads: 3
        log:
            "logs/{project}/host_filtering/index_host_alignment_{sample}.log",
        conda:
            "../envs/minimap2.yaml"
        shell:
            "samtools index {input} > {log} 2>&1"


    rule filter_host:
        input:
            bam="results/{project}/host_filtering/alignments/{sample}.bam",
            bai="results/{project}/host_filtering/alignments/{sample}.bam.bai",
        output:
            non_host=expand("results/{{project}}/host_filtering/non_host/{{sample}}_{reads}.fastq.gz",reads=["R1","R2"]),
        threads: 3
        log:
            "logs/{project}/host_filtering/filter_host_{sample}.log", 
        conda:
            "../envs/minimap2.yaml"
        shell:
            "(samtools fastq -F 3584 -f 77 {input.bam} | gzip -c > {output.non_host[0]} && "
            "samtools fastq -F 3584 -f 141 {input.bam} | gzip -c > {output.non_host[1]}) > {log} 2>&1"
    

    rule host_filtering_summary:
        input:
            csv="results/{project}/report/kraken2_summary.csv",
            jsons=expand(
                "results/{{project}}/trimmed/fastp/{sample}.fastp.json",
                sample=get_samples(),
            ),
        output:
            csv="results/{project}/report/host_filtering_summary.csv",
        params:
            host_name=config["host_filtering"]["host_name"]
        log:
            "logs/{project}/host_filtering/summary.log",
        threads: 2
        conda:
            "../envs/python.yaml"
        script:
            "../scripts/host_filtering_summary.py"

        
    use rule kraken2_report as host_filtering_report with:
        input:
            "results/{project}/report/host_filtering_summary.csv",
        output:
            report(
                directory("results/{project}/report/host_filtering/"),
                htmlindex="index.html",
                category="1. Quality control",
                labels={"sample": "all samples"},
            ),
        params:
            pin_until="sample",
            styles="resources/report/tables/",
            name="host_filtering",
            header="Filtering out {} reads".format(config["host_filtering"]["host_name"]),
        log:
            "logs/{project}/report/host_filtering_rbt_csv.log",


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


rule bracken_genus:
    input:
        hfile=get_kraken_db_file(),
        kreport="results/{project}/filtered/kraken_postfilt/{sample}_report.tsv",
    output:
        breport="results/{project}/report/bracken/reports_genus/{sample}.breport",
        bfile="results/{project}/report/bracken/files_genus/{sample}.bracken",
    params:
        db=lambda wildcards, input: Path(input.hfile).parent,
        level= "G",
    log:
        "logs/{project}/bracken/{sample}_genus.log",
    resources:
        mem_mb=100,
    threads: 1
    conda:
        "../envs/bracken.yaml"
    shell:
        "bracken -d {params.db} -i {input.kreport} -l {params.level} -o {output.bfile} -w {output.breport} > {log} 2>&1"


use rule bracken_genus as bracken_family with:
    input:
        hfile=get_kraken_db_file(),
        kreport="results/{project}/filtered/kraken_postfilt/{sample}_report.tsv",
    output:
        breport="results/{project}/report/bracken/reports_family/{sample}.breport",
        bfile="results/{project}/report/bracken/files_family/{sample}.bracken",
    params:
        db=lambda wildcards, input: Path(input.hfile).parent,
        level= "F",
    log:
        "logs/{project}/bracken/{sample}_family.log",


use rule bracken_genus as bracken_phylum with:
    input:
        hfile=get_kraken_db_file(),
        kreport="results/{project}/filtered/kraken_postfilt/{sample}_report.tsv",
    output:
        breport="results/{project}/report/bracken/reports_phylum/{sample}.breport",
        bfile="results/{project}/report/bracken/files_phylum/{sample}.bracken",
    params:
        db=lambda wildcards, input: Path(input.hfile).parent,
        level= "P",
    log:
        "logs/{project}/bracken/{sample}_phylum.log",


use rule bracken_genus as bracken_class with:
    input:
        hfile=get_kraken_db_file(),
        kreport="results/{project}/filtered/kraken_postfilt/{sample}_report.tsv",
    output:
        breport="results/{project}/report/bracken/reports_class/{sample}.breport",
        bfile="results/{project}/report/bracken/files_class/{sample}.bracken",
    params:
        db=lambda wildcards, input: Path(input.hfile).parent,
        level= "C",
    log:
        "logs/{project}/bracken/{sample}_class.log",


rule merge_bracken:
    input:
        expand(
            "results/{{project}}/report/bracken/files_{{level}}/{sample}.bracken",
            sample=get_samples(),
        ),
    output:
        "results/{project}/report/bracken/merged.bracken_{level}.txt",
    log:
        "logs/{project}/bracken/merge_bracken_{level}.log",
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
        "results/{project}/report/bracken/merged.bracken_{level}.txt",
    output:
        report(
            "results/{project}/report/bracken_{level}_plot.png",
            caption="../report/bracken_plot.rst",
            category="2. Species diversity",
            subcategory="2.2 post-filtering human reads",
            labels={"sample": "all samples", "level":"{level}"},
        ),
    params:
        threshold=0.001,
    log:
        "logs/{project}/report/bracken_{level}_plot.log",
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
            caption="../report/kraken_krona.rst",
            category="2. Species diversity",
            subcategory="2.2 post-filtering human reads",
            labels={"sample": "{sample}"},
        ),
    threads: 1
    log:
        "logs/{project}/krona/{sample}.log",
    conda:
        "../envs/kaiju.yaml"
    shell:
        "(ktImportText -o {output} {input}) > {log} 2>&1"
