from pathlib import Path


rule metabinner_unzip_fastqs:
    input:
        fq1=get_fastqs,
        fq2=get_fastqs,
    output:
        fq1_unzipped=temp("results/{project}/data/intermediate/{sample}_1.fastq"),
        fq2_unzipped=temp("results/{project}/data/intermediate/{sample}_2.fastq"),
    log:
        "logs/{project}/metabinner/{sample}/unzip_fastqs.log",
    conda:
        "../envs/metabinner_env.yaml"
    shell:
        "gunzip -c {input.fq1} > {output.fq1_unzipped};"
        "gunzip -c {input.fq2} > {output.fq2_unzipped}"


rule metabinner_final_contig:
    input:
        contig_file="results/{project}/assembly/{sample}/final.contigs.fa",
    output:
        contig_file="results/{project}/assembly/{sample}/final.contigs_{threshold}.fa",
    log:
        "logs/{project}/metabinner/{sample}/final_contig_{threshold}.log",
    conda:
        "../envs/metabinner_env.yaml"
    shell:
        "python $CONDA_PREFIX/bin/scripts/Filter_tooshort.py {input.contig_file} {wildcards.threshold}"


rule metabinner_coverage_profile:
    input:
        contig_file=expand(
            "results/{{project}}/assembly/{{sample}}/final.contigs_{threshold}.fa",
            threshold=get_threshold(),
        ),
        fastq1="results/{project}/data/intermediate/{sample}_1.fastq",
        fastq2="results/{project}/data/intermediate/{sample}_2.fastq",
    output:
        outfile="results/{project}/metabinner/{sample}/coverage_profile/coverage_profile.tsv",
    threads: 8
    params:
        threads=config["binning"]["threads"],
        threshold=config["binning"]["min_contig_length"],
        outdir=lambda wildcards, output: Path(output.outfile).parent,
    log:
        "logs/{project}/metabinner/{sample}/coverage_profile.log",
    conda:
        "../envs/metabinner_env.yaml"
    shell:
        "bash $CONDA_PREFIX/bin/scripts/gen_coverage_file.sh "
        "-t {params.threads} "
        "-a {input.contig_file} "
        "-o {params.outdir} "
        "-l {params.threshold} "
        "{input.fastq1} {input.fastq2} 2>{log}"


rule metabinner_composition_profile:
    input:
        contig_file="results/{project}/assembly/{sample}/final.contigs_{threshold}.fa",
    output:
        outfile="results/{project}/metabinner/{sample}/composition_profile/final.contigs_{threshold}_kmer_{kmer_size}_f{threshold}.csv",
    threads: 8
    params:
        root=get_root(),
        outdir=lambda wildcards, output: Path(output.outfile).parent,
    log:
        "logs/{project}/metabinner/{sample}/composition_profile_{kmer_size}_{threshold}.log",
    conda:
        "../envs/metabinner_env.yaml"
    shell:
        "python {params.root}/workflow/scripts/gen_kmer.py {input.contig_file} {wildcards.threshold} {wildcards.kmer_size} 2>{log}; "
        "mv results/{wildcards.project}/assembly/{wildcards.sample}/final.contigs_{wildcards.threshold}_kmer_{wildcards.kmer_size}_f{wildcards.threshold}.csv {params.outdir} 2>>{log}"


rule metabinner_run:
    input:
        contig_file=expand(
            "results/{{project}}/assembly/{{sample}}/final.contigs_{threshold}.fa",
            threshold=get_threshold(),
        ),
        coverage_profile="results/{project}/metabinner/{sample}/coverage_profile/coverage_profile.tsv",
        kmer_profile=expand(
            "results/{{project}}/metabinner/{{sample}}/composition_profile/final.contigs_{threshold}_kmer_{kmer_size}_f{threshold}.csv",
            kmer_size=get_kmersize(),
            threshold=get_threshold(),
        ),
    output:
        outfile="results/{project}/metabinner/{sample}/metabinner_res/metabinner_result.tsv",
    threads: 8
    params:
        threads=config["binning"]["threads"],
        outdir=lambda wildcards, output: Path(output.outfile).parent.parent,
        root=get_root(),
    log:
        "logs/{project}/metabinner/{sample}/metabinner.log",
    conda:
        "../envs/metabinner_env.yaml"
    shell:
        "run_metabinner.sh "
        "-a {params.root}/{input.contig_file} "
        "-d {params.root}/{input.coverage_profile} "
        "-k {params.root}/{input.kmer_profile} "
        "-o {params.root}/{params.outdir} "
        "-p $CONDA_PREFIX/bin/ "
        "-t {params.threads} > {log} 2>&1"
