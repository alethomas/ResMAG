from pathlib import Path


rule clone_rosella:
    output:
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
        "../envs/rosella.yml"
    shell:
        "cd {params.folder} && bash install.sh && cd {params.root} > {params.root}/{log} 2>&1"


rule rosella_run:
    input:
        contigs="results/{project}/assembly/{sample}/final.contigs.fa",
        abd="results/{project}/binning_prep/{sample}/abundance_metabat.tsv",
        log_install="logs/{project}/rosella/rosella_install.log",
    output:
        outfile="results/{project}/rosella/{sample}/rosella_kmer_table.tsv",
    params:
        outdir=lambda wildcards, output: Path(output.outfile).parent,
    threads: 16
    log:
        "logs/{project}/rosella/{sample}/rosella_run.log",
    conda:
        "../envs/rosella.yml"
    shell:
        "rosella recover -r {input.contigs} --coverage-values {input.abd} -o {params.outdir} -t {threads} > {log} 2>&1"
