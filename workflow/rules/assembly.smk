from pathlib import Path


rule megahit:
    input:
        fastqs=get_bacterial_reads,
    output:
        contigs="results/{project}/assembly/{sample}/final.contigs.fa",
    params:
        outdir=lambda wildcards, output: Path(output.contigs).parent,
    threads: 8
    log:
        "logs/{project}/megahit/{sample}.log",
    conda:
        "../envs/megahit.yaml"
    shell:
        "megahit -1 {input.fastqs[0]} -2 {input.fastqs[1]} --out-dir {params.outdir} -f > {log} 2>&1"


rule assembly_summary:
    input:
        "logs/{project}/megahit/{sample}.log",
    output:
        "results/{project}/assembly/{sample}_summary.csv",
    log:
        "logs/{project}/assembly_summary/{sample}.log",
    conda:
        "../envs/kraken2.yaml"
    script:
        "../scripts/assembly_summary.py"
