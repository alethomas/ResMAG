from pathlib import Path


rule megahit:
    input:
        bact_fq1="results/{project}/filtered/kraken/{sample}_bacteria_1.fastq",
        bact_fq2="results/{project}/filtered/kraken/{sample}_bacteria_2.fastq",
    output:
        contigs="results/{project}/assembly_megahit/{sample}/final.contigs.fa",
    params:
        outdir=lambda wildcards, output: Path(output.contigs).parent,
    threads: 8
    log:
        "logs/{project}/assembly_megahit/{sample}.log",
    conda:
        "../envs/megahit.yaml"
    shell:
        "megahit -1 {input.bact_fq1} -2 {input.bact_fq2} --out-dir {params.outdir} -f > {log} 2>&1"


rule assembly_summary:
    input:
        "logs/{project}/assembly_megahit/{sample}.log",
    output:
        "results/{project}/assembly_megahit/{sample}_summary.csv",
    log:
        "logs/{project}/assembly_summary/{sample}.log",
    conda:
        "../envs/kraken2.yaml"
    script:
        "../scripts/assembly_summary.py"
