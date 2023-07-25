from pathlib import Path


rule clone_rosella:
    output:
        install=get_rosella_install(),
    log:
        "logs/rosella_clone.log",
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
        "logs/rosella_install.log",
    log:
        "logs/rosella_install.log",
    params:
        folder=lambda wildcards, input: Path(input.install).parent,
        root=get_root(),
    conda:
        "../envs/rosella.yml"
    shell:
        "(cd {params.folder} && "
        "bash -v install.sh && "
        "cd {params.root}) > {params.root}/{log} 2>&1"


rule checkm_resources:
    output:
        marker="resources/checkm/taxon_marker_sets.tsv",
        tar="resources/checkm/checkm_data.tar.gz",
    params:
        outdir=lambda wildcards, output: Path(output.tar).parent,
        url="https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz",
    log:
        "logs/checkm_download.log",
    threads: 2
    conda:
        "../envs/unix.yaml"
    shell:
        "(curl --create-dirs -o {output.tar} {params.url} && "
        "tar -zvxf {output.tar} -C {params.outdir}) > {log} 2>&1"


rule checkm_init:
    input:
        marker="resources/checkm/taxon_marker_sets.tsv",
    output:
        "logs/checkm_init.log",
    params:
        direct=lambda wildcards, input: Path(input.marker).parent,
    log:
        "logs/checkm_init.log",
    conda:
        "../envs/rosella.yml"
    shell:
        "checkm data setRoot {params.direct} > {log} 2>&1"


rule rosella_run:
    input:
        contigs="results/{project}/assembly/{sample}/final.contigs.fa",
        abd="results/{project}/binning_prep/{sample}/abundance_metabat.tsv",
        log_install="logs/rosella_install.log",
    output:
        outfile="results/{project}/rosella/{sample}/rosella_kmer_table.tsv",
        outdir=directory("results/{project}/rosella/{sample}/"),
    threads: 16
    log:
        "logs/{project}/rosella/{sample}/rosella_run.log",
    conda:
        "../envs/rosella.yml"
    shell:
        "rosella recover -r {input.contigs} --coverage-values {input.abd} "
        "-o {output.outdir} -t {threads} > {log} 2>&1"
