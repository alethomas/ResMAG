
from pathlib import Path
'''
rule metaspades:
    input:
        bact_fq1="results/{project}/filtered/kraken/{sample}_bacteria_1.fastq",
        bact_fq2="results/{project}/filtered/kraken/{sample}_bacteria_2.fastq",
    output:
        contigs = "results/{project}/assembly/{sample}/contigs.fasta",
        scaffolds = "results/{project}/assembly/{sample}/scaffolds.fasta",
    params:
        outdir = lambda wildcards, output: Path(output.scaffolds).parent,
        mem = "--memory 180", #180GB ma
    threads: 8
    log:
        "logs/{project}/assembly/{sample}.log"
    conda:
        "../envs/metaspades.yaml"
    shell:
        "spades.py --meta --threads {threads} {params.mem} "
        "-1 {input.bact_fq1} -2 {input.bact_fq2} -o {params.outdir} > {log} 2>&1"


rule metaspades:
    input:
        reads=["results/{project}/filtered/kraken/{sample}_bacteria_1.fastq", 
            "results/{project}/filtered/kraken/{sample}_bacteria_2.fastq"],
    output:
        contigs="results/{project}/assembly/{sample}_contigs.fasta",
        scaffolds="results/{project}/assembly/{sample}_scaffolds.fasta",
        dir=directory("results/{project}/assembly/{sample}_intermediate_files"),
    log:
        "logs/{project}/assembly/{sample}.log",
    threads: 16
    resources:
        mem_mb=180000,
    wrapper:
        "v1.27.0/bio/spades/metaspades"
'''

rule megahit:
    input:
        bact_fq1="results/{project}/filtered/kraken/{sample}_bacteria_1.fastq",
        bact_fq2="results/{project}/filtered/kraken/{sample}_bacteria_2.fastq",
    output:
        contigs = "results/{project}/assembly_megahit/{sample}/final.contigs.fa",
    params:
        outdir = lambda wildcards, output: Path(output.contigs).parent,
    threads: 8
    log:
        "logs/{project}/assembly_megahit/{sample}.log"
    conda:
        "../envs/megahit.yaml"
    shell:
        "megahit -1 {input.bact_fq1} -2 {input.bact_fq2} --out-dir {params.outdir} -f > {log} 2>&1"
