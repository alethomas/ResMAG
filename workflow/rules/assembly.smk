from pathlib import Path


include: "host_filtering.smk"


rule megahit:
    input:
        fastqs=get_filtered_gz_fastqs,
    output:
        contigs=temp("results/{project}/megahit/{sample}/final.contigs.fa"),
        outdir=temp(directory("results/{project}/megahit/{sample}/")),
    params:
        threshold=get_contig_length_threshold(),
    threads: 64
    log:
        "logs/{project}/assembly/{sample}_megahit.log",
    conda:
        "../envs/megahit.yaml"
    shell:
        "megahit -1 {input.fastqs[0]} -2 {input.fastqs[1]} "
        "--min-contig-len {params.threshold} -t {threads} "
        "--out-dir {output.outdir} -f > {log} 2>&1"


rule remove_megahit_intermediates:
    input:
        contigs=get_assembly,
    output:
        touch("logs/{project}/assembly/{sample}_intermediate_removal.done"),
    params:
        outdir=lambda wildcards, input: Path(input.contigs).parent,
    log:
        "logs/{project}/assembly/{sample}_intermediate_removal.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "(find {params.outdir}/ -mindepth 1 -type d -exec rm -rf {{}} +) > {log} 2>&1"


rule gzip_assembly:
    input:
        contigs=get_assembly,
    output:
        "results/{project}/output/fastas/{sample}/{sample}_final.contigs.fa.gz",
    threads: 4
    log:
        "logs/{project}/assembly/{sample}_gzip.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "gzip -c {input.contigs} > {output} 2> {log}"


rule assembly_summary:
    input:
        expand(
            "logs/{{project}}/assembly/{sample}_megahit.log",
            sample=get_samples(),
        ),
    output:
        csv="results/{project}/output/report/all/assembly_summary.csv",
    log:
        "logs/{project}/assembly/summary.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/assembly_summary.py"


use rule kraken2_report as assembly_report with:
    input:
        "results/{project}/output/report/all/assembly_summary.csv",
    output:
        report(
            directory("results/{project}/output/report/all/assembly/"),
            htmlindex="index.html",
            category="3. Assembly results",
            labels={
                "sample": "all samples",
            },
        ),
    params:
        pin_until="sample",
        styles="resources/report/tables/",
        name="assembly_summary",
        header="Assembly summary",
    log:
        "logs/{project}/report/assembly_rbt_csv.log",
