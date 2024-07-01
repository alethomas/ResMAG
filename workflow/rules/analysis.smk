rule load_genomad_DB:
    output:
        file=get_genomad_DB_file(),
    params:
        folder=lambda wildcards, output: Path(output.file).parent.parent,
    log:
        "logs/load_genomad_DB.log",
    conda:
        "../envs/genomad.yaml"
    shell:
        "genomad download-database {params.folder}/ > {log} 2>&1"


rule genomad_run:
    input:
        db=rules.load_genomad_DB.output.file,
        asmbl=rules.gzip_assembly.output,
    output:
        plasmid_tsv="results/{project}/output/plasmids/{sample}/{sample}_final.contigs_summary/{sample}_final.contigs_plasmid_summary.tsv",
        virus_tsv="results/{project}/output/plasmids/{sample}/{sample}_final.contigs_summary/{sample}_final.contigs_virus_summary.tsv",
    params:
        db_folder=lambda wildcards, input: Path(input.db).parent,
        outdir=lambda wildcards, output: Path(output.plasmid_tsv).parent.parent,
    log:
        "logs/{project}/plasmids/{sample}.log",
    threads: 64
    conda:
        "../envs/genomad.yaml"
    shell:
        "genomad end-to-end --cleanup -t {threads} "
        "{input.asmbl} {params.outdir}/ "
        "{params.db_folder}/ > {log} 2>&1"


rule download_CARD_data:
    output:
        json=get_card_db_file(),
    params:
        download=config["card"]["url"],
        folder=lambda wildcards, output: Path(output.json).parent,
        filename=lambda wildcards, output: Path(output.json).name,
    log:
        "logs/CARD_data_download.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "(cd {params.folder} && "
        "wget {params.download} && "
        "tar -xvf data ./{params.filename}) > {log} 2>&1"


rule CARD_load_DB:
    input:
        get_card_db_file(),
    output:
        touch("results/CARD_load_DB.done"),
    log:
        "logs/CARD_load_DB.log",
    conda:
        "../envs/card.yaml"
    shell:
        "rgi clean --local && "
        "rgi load --card_json {input} --local > {log} 2>&1"


rule CARD_assembly_run:
    input:
        fa=rules.gzip_assembly.output,
        db=rules.CARD_load_DB.output,
    output:
        txt="results/{project}/output/ARGs/assembly/{sample}/{sample}.txt",
    params:
        path_wo_ext=lambda wildcards, output: Path(output.txt).with_suffix(""),
    log:
        "logs/{project}/ARGs/assembly/{sample}.log",
    threads: 64
    conda:
        "../envs/card.yaml"
    shell:
        "rgi main -i {input.fa} -o {params.path_wo_ext} "
        "-t contig -a DIAMOND --low_quality --local "
        "-n {threads} --clean > {log} 2>&1"


use rule CARD_assembly_run as CARD_bin_run with:
    input:
        fa="results/{project}/output/fastas/{sample}/bins/{binID}.fa.gz",
        db=rules.CARD_load_DB.output,
    output:
        txt="results/{project}/output/ARGs/bins/{sample}/{binID}.txt",
    params:
        path_wo_ext=lambda wildcards, output: Path(output.txt).with_suffix(""),
    log:
        "logs/{project}/ARGs/bins/{sample}/{binID}.log",


rule wrap_bin_ARGs:
    input:
        get_bin_ARGs,
        #expand("results/{{project}}/output/ARGs/bins/{{sample}}/{binID}.txt", binID=get_binIDs_for_sample),
    output:
        "results/{project}/output/ARGs/bins/{sample}/all_bins.done",
    log:
        "logs/{project}/ARGs/bins/{sample}/all_bins.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "touch {output} > {log} 2>&1"
