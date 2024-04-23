from pathlib import Path


include: "host_filtering.smk"


rule megahit:
    input:
        fastqs=get_filtered_gz_fastqs,
    output:
        contigs="results/{project}/megahit/{sample}_to_{pathogen}/final.contigs.fa",
    params:
        outdir=lambda wildcards, output: Path(output.contigs).parent,
        threshold=get_contig_length_threshold(),
    threads: 64
    log:
        "logs/{project}/assembly/{sample}_to_{pathogen}/megahit.log",
    conda:
        "../envs/megahit.yaml"
    shell:
        "megahit -1 {input.fastqs[0]} -2 {input.fastqs[1]} "
        "--min-contig-len {params.threshold} -t {threads} "
        "--out-dir {params.outdir} -f > {log} 2>&1"


rule remove_megahit_intermediates:
    input:
        contigs=get_assembly,
    output:
        touch(
            "logs/{project}/assembly/{sample}_to_{pathogen}_intermediate_removal.done"
        ),
    params:
        outdir=lambda wildcards, input: Path(input.contigs).parent,
    log:
        "logs/{project}/assembly/{sample}_to_{pathogen}_intermediate_removal.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "(find {params.outdir}/ -mindepth 1 -type d -exec rm -rf {{}} +) > {log} 2>&1"


rule gzip_assembly:
    input:
        contigs=get_assembly,
    output:
        "results/{project}/output/fastas/{sample}/{sample}_to_{pathogen}_final.contigs.fa.gz",
    threads: 4
    log:
        "logs/{project}/assembly/{sample}_to_{pathogen}_gzip.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "gzip -c {input.contigs} > {output} 2> {log}"


rule assembly_summary:
    input:
        expand(
            "logs/{{project}}/assembly/{sample}_to_{pathogen}_megahit.log",
            sample=get_samples(),
            pathogen=get_pathogens(),
        ),
    output:
        csv="results/{project}/output/report/all/assembly_summary.csv",
    log:
        "logs/{project}/assembly/summary.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/assembly_summary.py"
