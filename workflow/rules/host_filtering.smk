from pathlib import Path


if config["human-filtering"]["use-local"]:

    rule copy_local_human_ref:
        output:
            fasta=get_human_ref(),
        params:
            local=get_human_local_folder(),
            folder=lambda wildcards, output: Path(output.fasta).parent,
            file=lambda wildcards, output: Path(output.fasta).name,
        log:
            "logs/human_ref_local_copy.log",
        group:
            "refGenome_depended"
        conda:
            "../envs/unix.yaml"
        shell:
            "(mkdir -p {params.folder} && "
            "tar cpfz - -C {params.local} {params.file} | "
            "(cd {params.folder} ; tar xpfz -)) > {log} 2>&1"

else:

    rule download_human_ref:
        output:
            fasta=get_human_ref(),
        params:
            download=config["human-filtering"]["download-path"],
            folder=lambda wildcards, output: Path(output.fasta).parent,  #get_resource_path(),
        log:
            "logs/human_ref_download.log",
        group:
            "refGenome_depended"
        conda:
            "../envs/unix.yaml"
        shell:
            "(mkdir -p {params.folder} && "
            "cd {params.folder} && "
            "wget {params.download}) > {log} 2>&1"


rule map_to_human:
    input:
        fastqs=get_prefiltered_fastqs,
        ref=get_human_ref(),
    output:
        bam=temp("results/{project}/human_filtering/alignments/{sample}.bam"),
    threads: 64
    log:
        "logs/{project}/human_filtering/map_to_human_{sample}.log",
    conda:
        "../envs/minimap2.yaml"
    shell:
        "(minimap2 -a -xsr -t {threads} {input.ref} {input.fastqs} | "
        "samtools view -bh | "
        "samtools sort --threads {threads} -o {output.bam}) > {log} 2>&1"


rule index_human_alignment:
    input:
        rules.map_to_human.output.bam,
    output:
        bai=temp("results/{project}/human_filtering/alignments/{sample}.bam.bai"),
    threads: 20
    log:
        "logs/{project}/human_filtering/index_human_alignment_{sample}.log",
    conda:
        "../envs/minimap2.yaml"
    shell:
        "samtools index {input} > {log} 2>&1"


rule filter_human:
    input:
        bam=rules.map_to_human.output.bam,
        bai=rules.index_human_alignment.output.bai,
    output:
        filtered=temp(
            expand(
                "results/{{project}}/filtered/fastqs/{{sample}}_{read}.fastq",
                read=["R1", "R2"],
            )
        ),
    threads: 64
    log:
        "results/{project}/report_prerequisites/qc/filter_human_{sample}.log",
    conda:
        "../envs/minimap2.yaml"
    shell:
        "(samtools fastq --threads {threads} -F 3584 -f 77 "
        "-o {output.filtered[0]} {input.bam} && "
        "samtools fastq --threads {threads} -F 3584 -f 141 "
        "-o {output.filtered[1]} {input.bam}) > {log} 2>&1"


rule gzip_filtered_reads:
    input:
        "results/{project}/filtered/fastqs/{sample}_{read}.fastq",
    output:
        "results/{project}/filtered/fastqs/{sample}_{read}.fastq.gz",
    log:
        "logs/{project}/human_filtering/gzip_{sample}_{read}.log",
    threads: 64
    conda:
        "../envs/unix.yaml"
    shell:
        "pigz -k {input} > {log} 2>&1"


if config["host_filtering"]["do_host_filtering"]:
    # if there is a different host than human,
    # this is run before filtering human reads

    use rule map_to_human as map_to_host with:
        input:
            fastqs=get_trimmed_fastqs,
            ref=config["host_filtering"]["ref_genome"],
        output:
            bam=temp("results/{project}/host_filtering/alignments/{sample}.bam"),
        params:
            ref=config["host_filtering"]["ref_genome"],
        threads: 20
        log:
            "logs/{project}/host_filtering/map_to_host_{sample}.log",

    use rule index_human_alignment as index_host_alignment with:
        input:
            rules.map_to_host.output.bam,
        output:
            bai=temp("results/{project}/host_filtering/alignments/{sample}.bam.bai"),
        threads: 3
        log:
            "logs/{project}/host_filtering/index_host_alignment_{sample}.log",

    rule filter_host:
        input:
            bam=rules.map_to_host.output.bam,
            bai=rules.index_host_alignment.output.bai,
        output:
            filtered=temp(
                expand(
                    "results/{{project}}/host_filtering/non_host/{{sample}}_{read}.fastq.gz",
                    read=["R1", "R2"],
                )
            ),
        threads: 64
        log:
            "results/{project}/report_prerequisites/qc/filter_host_{sample}.log",
        conda:
            "../envs/minimap2.yaml"
        shell:
            "(samtools fastq -F 3584 -f 77 {input.bam} | "
            "pigz -c > {output.filtered[0]} && "
            "samtools fastq -F 3584 -f 141 {input.bam} | "
            "pigz -c > {output.filtered[1]}) > {log} 2>&1"


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
            fastqs=get_filtered_gz_fastqs,
        output:
            report="results/{project}/output/classification/reads/{sample}/{sample}_kraken2_report.tsv",
            outfile=temp(
                "results/{project}/output/classification/reads/{sample}/{sample}_kraken2_outfile.tsv"
            ),
        params:
            db=lambda wildcards, input: Path(input.hfile).parent,
        threads: 64
        log:
            "logs/{project}/kraken2/run/{sample}.log",
        conda:
            "../envs/kraken2.yaml"
        shell:
            "kraken2 --db {params.db} --threads {threads} --paired "
            "--output {output.outfile} --report {output.report} "
            "--gzip-compressed {input.fastqs} > {log} 2>&1"


if config["testing"]:

    rule gunzip_kraken_output:
        input:
            clf_gz_1="results/{project}/filtered/{sample}_clf_1.fastq.gz",
            clf_gz_2="results/{project}/filtered/{sample}_clf_2.fastq.gz",
        output:
            clf1=temp("results/{project}/filtered/{sample}_clf_1.fastq"),
            clf2=temp("results/{project}/filtered/{sample}_clf_2.fastq"),
        threads: 64
        shell:
            "pigz -dk {input.clf_gz_1} {input.clf_gz_1}"


## all reports are done on human (+ optional other different host) filtered
rule diversity_summary:
    input:
        reports=expand(
            "results/{{project}}/output/classification/reads/{sample}/{sample}_kraken2_report.tsv",
            sample=get_samples(),
        ),
        jsons=expand(
            "results/{{project}}/trimmed/fastp/{sample}.fastp.json",
            sample=get_samples(),
        ),
        human_logs=expand(
            "results/{{project}}/report_prerequisites/qc/filter_human_{sample}.log",
            sample=get_samples(),
        ),
        host_logs=get_host_map_statistics,
    output:
        csv="results/{project}/output/report/all/diversity_summary.csv",
    log:
        "logs/{project}/kraken2/summary.log",
    params:
        taxid_dict=get_taxID_dict(),
        other_host=config["host_filtering"]["do_host_filtering"],
        hostname=config["host_filtering"]["host_name"],
    threads: 2
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/diversity_summary.py"


use rule qc_summary_report as diversity_summary_report with:
    input:
        rules.diversity_summary.output.csv,
    output:
        report(
            directory("results/{project}/output/report/all/diversity_summary/"),
            htmlindex="index.html",
            caption="../report/kraken.rst",
            category="2. Species diversity",
            labels={
                "sample": "all samples",
            },
        ),
    params:
        pin_until="sample",
        styles="resources/report/tables/",
        name="diversity_summary",
        header="Diversity summary based on mapping to host genome(s) and Kraken2",
        pattern=config["tablular-config"],
    log:
        "logs/{project}/report/kraken2_rbt_csv.log"


rule create_host_plot:
    input:
        csv=rules.diversity_summary.output.csv,
    output:
        html=report(
            "results/{project}/output/report/all/host_contamination.html",
            caption="../report/host_plot.rst",
            category="2. Species diversity",
            labels={"sample": "all samples"},
        ),
    params:
        other_host=config["host_filtering"]["do_host_filtering"],
        hostname=config["host_filtering"]["host_name"],
    log:
        "logs/{project}/report/host_plot.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot_host.py"


rule bracken_genus:
    input:
        hfile=get_kraken_db_file(),
        kreport=rules.kraken2.output.report,
    output:
        breport=temp(
            "results/{project}/output/report/all/diversity_abundance/reports_genus/{sample}.breport"
        ),
        bfile=temp(
            "results/{project}/output/report/all/diversity_abundance/files_genus/{sample}.bracken"
        ),
    params:
        db=lambda wildcards, input: Path(input.hfile).parent,
        level="G",
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
        kreport=rules.kraken2.output.report,
    output:
        breport=temp(
            "results/{project}/output/report/all/diversity_abundance/reports_family/{sample}.breport"
        ),
        bfile=temp(
            "results/{project}/output/report/all/diversity_abundance/files_family/{sample}.bracken"
        ),
    params:
        db=lambda wildcards, input: Path(input.hfile).parent,
        level="F",
    log:
        "logs/{project}/bracken/{sample}_family.log",


use rule bracken_genus as bracken_phylum with:
    input:
        hfile=get_kraken_db_file(),
        kreport=rules.kraken2.output.report,
    output:
        breport=temp(
            "results/{project}/output/report/all/diversity_abundance/reports_phylum/{sample}.breport"
        ),
        bfile=temp(
            "results/{project}/output/report/all/diversity_abundance/files_phylum/{sample}.bracken"
        ),
    params:
        db=lambda wildcards, input: Path(input.hfile).parent,
        level="P",
    log:
        "logs/{project}/bracken/{sample}_phylum.log",


use rule bracken_genus as bracken_class with:
    input:
        hfile=get_kraken_db_file(),
        kreport=rules.kraken2.output.report,
    output:
        breport=temp(
            "results/{project}/output/report/all/diversity_abundance/reports_class/{sample}.breport"
        ),
        bfile=temp(
            "results/{project}/output/report/all/diversity_abundance/files_class/{sample}.bracken"
        ),
    params:
        db=lambda wildcards, input: Path(input.hfile).parent,
        level="C",
    log:
        "logs/{project}/bracken/{sample}_class.log",


rule merge_bracken:
    input:
        expand(
            "results/{{project}}/output/report/all/diversity_abundance/files_{{level}}/{sample}.bracken",
            sample=get_samples(),
        ),
    output:
        "results/{project}/output/report/all/diversity_abundance/merged.bracken_{level}.txt",
    log:
        "logs/{project}/bracken/merge_bracken_{level}.log",
    resources:
        mem_mb=100,
    threads: 1
    conda:
        "../envs/bracken.yaml"
    shell:
        "(python $CONDA_PREFIX/bin/combine_bracken_outputs.py --files {input} --output {output}) > {log} 2>&1"


rule create_bracken_plot:
    input:
        "results/{project}/output/report/all/diversity_abundance/merged.bracken_{level}.txt",
    output:
        report(
            "results/{project}/output/report/all/abundance_{level}.html",
            caption="../report/bracken_plot.rst",
            category="2. Species diversity",
            labels={"sample": "all samples", "level": "{level}"},
        ),
    params:
        # all level values below 1% will be summed up as other
        threshold=0.01,
    log:
        "logs/{project}/report/bracken_{level}_plot.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/brackenplot.py"
