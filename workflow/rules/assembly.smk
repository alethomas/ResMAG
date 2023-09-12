from pathlib import Path


rule megahit:
    input:
        fastqs=get_filtered_gz_fastqs,
    output:
        contigs="results/{project}/megahit/{sample}/final.contigs.fa",
    params:
        outdir=lambda wildcards, output: Path(output.contigs).parent,
    threads: 64
    log:
        "logs/{project}/megahit/{sample}.log",
    conda:
        "../envs/megahit.yaml"
    shell:
        "megahit -1 {input.fastqs[0]} -2 {input.fastqs[1]} "
        "--out-dir {params.outdir} -f > {log} 2>&1"


rule gzip_assembly:
    input:
        "results/{project}/megahit/{sample}/final.contigs.fa",
    output:
        "results/{project}/assembly/{sample}_final.contigs.fa.gz",
    threads: 4
    log:
        "logs/{project}/assembly/{sample}_gzip.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "gzip -c {input} > {output} 2> {log}"


rule assembly_summary:
    input:
        "logs/{project}/megahit/{sample}.log",
    output:
        "results/{project}/assembly/{sample}_summary.csv",
    log:
        "logs/{project}/assembly/{sample}_summary.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/assembly_summary.py"
