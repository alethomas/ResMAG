from pathlib import Path

# MIN_CONTIG_LEN=get_contig_length_threshold()


rule megahit:
    input:
        fastqs=get_filtered_gz_fastqs,
    output:
        contigs="results/{project}/megahit/{sample}/final.contigs.fa",
        folder=temp(directory("results/{project}/megahit/{sample}/")),
        #contigs_len=f"results/{{project}}/megahit/{{sample}}/final.contigs_{MIN_CONTIG_LEN}.fa",
    params:
        outdir=lambda wildcards, output: Path(output.contigs).parent,
        threshold=get_contig_length_threshold(),
    threads: 64
    log:
        "logs/{project}/assembly/{sample}_megahit.log",
    conda:
        "../envs/megahit.yaml"
    shell:
        "megahit -1 {input.fastqs[0]} -2 {input.fastqs[1]} "
        "--min-contig-len {params.threshold} "
        "--out-dir {params.outdir} -f > {log} 2>&1"
        #"&& cp {output.contigs} {output.contigs_len}) "


rule gzip_assembly:
    input:
        contigs=get_assembly,
    output:
        "results/{project}/assembly/{sample}_final.contigs.fa.gz",
    threads: 4
    log:
        "logs/{project}/assembly/{sample}_gzip.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "gzip -c {input.contigs[0]} > {output} 2> {log}"


rule assembly_summary:
    input:
        expand(
            "logs/{{project}}/assembly/{sample}_megahit.log",
            sample=get_samples(),
        ),
    output:
        html=report(
            "results/{project}/report/assembly_summary.html",
            category="3. Assembly results",
            labels={"sample": "all samples", "type": "view"},
        ),
        csv=report(
            "results/{project}/report/assembly_summary.csv",
            category="3. Assembly results",
            labels={"sample": "all samples", "type": "download"},
        ),
    log:
        "logs/{project}/assembly/summary.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/assembly_summary.py"


rule csv_report:
    input:
        # a csv formatted file containing the data for the report
        "results/{project}/report/assembly_summary.csv",
    output:
        # path to the resulting report directory
        directory("results/{project}/report/assembly_test/"),
    log:
        "logs/{project}/assembly/rbt-csv-report.log",
    wrapper:
        "v2.6.0/bio/rbt/csvreport"
