from pathlib import Path


'''rule metabinner_unzip_fastqs:
    input:
        fastqs=get_filtered_reads,
    output:
        fq1_unzip=temp("results/{project}/temp/{sample}_1.fastq"),
        fq2_unzip=temp("results/{project}/temp/{sample}_2.fastq"),
    log:
        "logs/{project}/metabinner/{sample}/unzip_fastqs.log",
    threads: 2
    conda:
        "../envs/metabinner_env.yaml"
    shell:
        "(gunzip -c {input.fastqs[0]} > {output.fq1_unzip} && "
        "gunzip -c {input.fastqs[1]} > {output.fq2_unzip}) 2> {log}"'''


rule metabinner_filter_contigs:
    input:
        "results/{project}/assembly/{sample}/final.contigs.fa",
    output:
        "results/{project}/assembly/{sample}/final.contigs_{threshold}.fa",
    log:
        "logs/{project}/metabinner/{sample}/filter_contigs_{threshold}.log",
    threads: 2
    conda:
        "../envs/metabinner_env.yaml"
    shell:
        "python $CONDA_PREFIX/bin/scripts/Filter_tooshort.py "
        "{input} {wildcards.threshold} > {log} 2>&1"


rule metabinner_coverage_profile:
    input:
        contig_file=expand(
            "results/{{project}}/assembly/{{sample}}/final.contigs_{threshold}.fa",
            threshold=get_threshold(),
        ),
        fastqs=get_filtered_reads,
    output:
        outfile="results/{project}/metabinner/{sample}/coverage_profile/coverage_profile.tsv",
    threads: 32
    params:
        threshold=get_threshold(),
        threads=config["binning"]["threads"],
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
        "{input.fastqs[0]} {input.fastqs[1]} > {log} 2>&1"


rule metabinner_composition_profile:
    input:
        contigs="results/{project}/assembly/{sample}/final.contigs_{threshold}.fa",
    output:
        outfile="results/{project}/metabinner/{sample}/composition_profile/final.contigs_{threshold}_kmer_{kmer_size}_f{threshold}.csv",
    threads: 32
    params:
        #root=get_root(),
        script_path="{}/workflow/scripts/gen_kmer.py".format(get_root()),
        outdir=lambda wildcards, output: Path(output.outfile).parent,
        contigs_path=lambda wildcards, input: Path(input.contigs).parent,
        filename=lambda wildcards, output: Path(output.outfile).name,
    log:
        "logs/{project}/metabinner/{sample}/composition_profile_{kmer_size}_{threshold}.log",
    conda:
        "../envs/metabinner_env.yaml"
    shell:
        "(python {params.script_path} {input} {wildcards.threshold} {wildcards.kmer_size} && "
        "mv {params.contigs_path}/{params.filename} {params.outdir}) > {log} 2>&1"


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
    threads: 64
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
