from pathlib import Path

# ROSELLA_FOLDER=get_rosella_folder()


rule clone_rosella:
    output:
        yml=get_rosella_yaml(),
        install=get_rosella_install(),
    log:
        "logs/rosella/rosella_clone.log",
    params:
        git_url=get_rosella_git(),
        folder=lambda wildcards, output: Path(output.install).parent,
    conda:
        "../envs/unix.yaml"
    shell:
        "git clone --recursive {params.git_url} {params.folder} > {log} 2>&1"


rule install_rosella:
    input:
        install=get_rosella_install(),
    output:
        "logs/{project}/rosella/rosella_install.log",
    log:
        "logs/{project}/rosella/rosella_install.log",
    params:
        folder=lambda wildcards, input: Path(input.install).parent,
        root=get_root(),
    conda:
        "../envs/rosella.yml" #get_rosella_yaml()
    shell:
        "cd {params.folder} && bash install.sh && cd {params.root} > {params.root}/{log} 2>&1"


rule rosella_run:
    input:
        fasta="results/{project}/assembly/{sample}/final.contigs.fa",
        coverage="results/{project}/assembly/{sample}/abundance.tsv",
        fq1="results/{project}/filtered/bacteria/{sample}_bacteria_1.fastq.gz",
        fq2="results/{project}/filtered/bacteria/{sample}_bacteria_2.fastq.gz",
        log_install="logs/{project}/rosella/rosella_install.log",
    output:
        outfile="results/{project}/rosella/{sample}/rosella_res/rosella_kmer_table.tsv",
    params:
        outdir=lambda wildcards, output: Path(output.outfile).parent,
        threads=config["binning"]["threads"],
    threads: 8
    log:
        "logs/{project}/rosella/{sample}/rosella_run.log",
    conda:
        "../envs/rosella.yml" #get_rosella_yaml()
    shell:
        "rosella recover -r {input.fasta} -1 {input.fq1} -2 {input.fq2}  -o {params.outdir} -t {params.threads} 2> {log} "
