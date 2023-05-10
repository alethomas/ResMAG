rule vamb_contig_catalogue:
    input:
        contigs="results/{project}/assembly_megahit/{sample}/final.contigs.fa",
    output:
        "results/{project}/vamb/catalogue.fna.gz"
    script:
        "../vamb/src/concatenate.py"

rule metabinner_repo:
    output:
        directory("MetaBinner")
    shell:
        "git clone https://github.com/ziyewang/MetaBinner.git"

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

rule metabinner_coverage_profile:
    input:
        contig_file="results/{project}/assembly/{sample}/final.contigs.fa",
        fastq1=expand("data/intermediate/{sample}_1.fastq", sample=get_samples()),
        fastq2=expand("data/intermediate/{sample}_2.fastq", sample=get_samples()),
        git_repo="MetaBinner"
    output:
        "results/{project}/metabinner/{sample}/coverage_profile/coverage_profile_f1000.tsv"
    params:
        threads = 60,
        outdir="results/{project}/metabinner/{sample}/coverage_profile/"
    conda:
        "../envs/metabinner_env.yaml"
    shell:
        "bash {input.git_repo}/scripts/gen_coverage_file.sh -t {params.threads} -a {input.contig_file} -o {params.outdir} {input.fastq1} {input.fastq2}"

rule metabinner_composition_profile:
    input:
        contig_file="results/{project}/assembly/{sample}/final.contigs.fa",
        git_repo="MetaBinner"
    output:
        "results/{project}/metabinner/{sample}/composition_profile/final.contigs_kmer_4_f1000.csv"
    params:
        contig_length=1000,
        kmer_size=4,
        outdir="results/{project}/metabinner/{sample}/composition_profile"
    conda:
        "../envs/metabinner_env.yaml"
    shell:
        "cp results/{wildcards.project}/assembly/{wildcards.sample}/final.contigs.fa {params.outdir};"
        "python {input.git_repo}/scripts/gen_kmer.py {params.outdir}/final.contigs.fa {params.contig_length} {params.kmer_size}"

rule metabinner_run:
    input:
        contig_file="results/{project}/assembly/{sample}/final.contigs.fa",
        coverage_profile="results/{project}/metabinner/{sample}/coverage_profile/coverage_profile_f1000.tsv",
        kmer_profile="results/{project}/metabinner/{sample}/composition_profile/final.contigs_kmer_4_f1000.csv",
        git_repo="MetaBinner"
    output:
        "results/{project}/metabinner/{sample}/metabinner_res/metabinner_result.tsv"
    params:
        threads=60,
        outdir="results/{project}/metabinner/{sample}/",
        root_path=config["root_path"]
    log:
        "logs/{project}/metabinner/{sample}/metabinner.log"
    conda:
        "../envs/metabinner_env.yaml"
    shell:
        "bash {params.root_path}/{input.git_repo}/run_metabinner.sh "
        "-a {params.root_path}/{input.contig_file} "
        "-d {params.root_path}/{input.coverage_profile} "
        "-k {params.root_path}/{input.kmer_profile} "
        "-o {params.root_path}/{params.outdir} "
        "-p {params.root_path}/{input.git_repo} "
        "-t {params.threads}"