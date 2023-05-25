rule metabinner_unzip_fastqs:
    input:
        fq1=get_fastqs,
        fq2=get_fastqs
    output:
        fq1_unzipped=temp("data/intermediate/{sample}_1.fastq"),
        fq2_unzipped=temp("data/intermediate/{sample}_2.fastq")
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
        contig_file=expand("results/{{project}}/assembly/{{sample}}/final.contigs_{threshold}.fa", threshold=get_threshold()),
        fastq1="data/intermediate/{sample}_1.fastq",
        fastq2="data/intermediate/{sample}_2.fastq",
    output:
        "results/{project}/metabinner/{sample}/coverage_profile/coverage_profile.tsv"
    threads: 8
    params:
        threads=config["binning"]["threads"],
        threshold=config["binning"]["min_contig_length"],
        outdir="results/{project}/metabinner/{sample}/coverage_profile/",
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
        outfile="results/{project}/metabinner/{sample}/composition_profile/final.contigs_{threshold}_kmer_{kmer_size}_f{threshold}.csv"
    threads: 8
    params:
        root_path=config["root_path"],
        outdir="results/{project}/metabinner/{sample}/composition_profile"
    log:
        "logs/{project}/metabinner/{sample}/composition_profile_{kmer_size}_{threshold}.log",        
    conda:
        "../envs/metabinner_env.yaml"
    shell:
        "python {params.root_path}/workflow/scripts/gen_kmer.py {input.contig_file} {wildcards.threshold} {wildcards.kmer_size} 2>{log}; "
        "mv results/{wildcards.project}/assembly/{wildcards.sample}/final.contigs_{wildcards.threshold}_kmer_{wildcards.kmer_size}_f{wildcards.threshold}.csv {params.outdir} 2>>{log}"

rule metabinner_run:
    input:
        contig_file=expand("results/{{project}}/assembly/{{sample}}/final.contigs_{threshold}.fa",
        threshold=get_threshold()),
        coverage_profile="results/{project}/metabinner/{sample}/coverage_profile/coverage_profile.tsv",
        kmer_profile=expand(
            "results/{{project}}/metabinner/{{sample}}/composition_profile/final.contigs_{threshold}_kmer_{kmer_size}_f{threshold}.csv",
            kmer_size=get_kmersize(),
            threshold=get_threshold(),
            )
    output:
        "results/{project}/metabinner/{sample}/metabinner_res/metabinner_result.tsv"
    threads: 8
    params:
        threads=config["binning"]["threads"],
        outdir="results/{project}/metabinner/{sample}",
        root_path=config["root_path"]
    log:
        "logs/{project}/metabinner/{sample}/metabinner.log"
    conda:
        "../envs/metabinner_env.yaml"
    shell:
        "run_metabinner.sh "
        "-a {params.root_path}/{input.contig_file} "
        "-d {params.root_path}/{input.coverage_profile} "
        "-k {params.root_path}/{input.kmer_profile} "
        "-o {params.root_path}/{params.outdir} "
        "-p $CONDA_PREFIX/bin/ "
        "-t {params.threads} 2>{log}"