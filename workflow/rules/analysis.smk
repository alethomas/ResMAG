# Plasmid analysis
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


# Resistance analysis
rule download_CARD_data:
    output:
        json=get_card_db_file(),
    params:
        download=config["card"]["url"],
        folder=lambda wildcards, output: Path(output.json).parent,
        filename=lambda wildcards, output: Path(output.json).name,
    log:
        "logs/CARD_data_download.log",
    threads: 30
    conda:
        "../envs/unix.yaml"
    shell:
        "(cd {params.folder} && "
        "wget {params.download} && "
        "tar -xvf data ./{params.filename}) > {log} 2>&1"


rule CARD_load_DB_for_reads:
rule CARD_load_DB_for_reads:
    input:
        get_card_db_file(),
    output:
        temp(touch("results/CARD_load_DB_for_reads.done")),
    log:
        "logs/CARD_load_DB_for_reads.log",
    conda:
        "../envs/card.yaml"
    shell:
        "rgi clean --local && "
        "rgi load --card_json {input} --local > {log} 2>&1"


rule CARD_annotation:
    input:
        json=get_card_db_file(),
        load=rules.CARD_load_DB_for_reads.output,
    output:
        temp(touch("results/CARD_annotation.done")),
    params:
        folder=lambda wildcards, input: Path(input.json).parent,
        file=lambda wildcards, input: Path(input.json).name,
        ann=get_card_annotation_file(),
    log:
        "logs/CARD_annotation.log",
    conda:
        "../envs/card.yaml"
    shell:
        "(cd {params.folder}/ && "
        "rgi card_annotation -i {params.file}) && "
        "rgi load -i {input.json} --card_annotation {params.ann} --local > {log} 2>&1"


rule CARD_read_run:
    input:
        fa=get_filtered_gz_fastqs,
        db=rules.CARD_annotation.output,
    output:
        txt="results/{project}/output/ARGs/reads/{sample}/{sample}.gene_mapping_data.txt",
    params:
        folder=lambda wildcards, output: Path(output.txt).parent,
    log:
        "logs/{project}/ARGs/reads/{sample}.log",
    threads: 64
    conda:
        "../envs/card.yaml"
    shell:
        "rgi bwt -1 {input.fa[0]} -2 {input.fa[1]} "
        "-o {params.folder}/{wildcards.sample} --local "
        "--clean -n {threads} > {log} 2>&1"


rule CARD_read_sample_summary:
    input:
        txt=rules.CARD_read_run.output.txt,
    output:
        csv="results/{project}/output/ARGs/reads/{sample}/{sample}_read_ARGs.csv",
    params:
        case="reads",
    log:
        "logs/{project}/ARGs/reads/{sample}.log",
    threads: 2
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/arg_summary_sample.py"


# updates CARD database to use for contigs instead of reads
# read based classification must be finished before
rule CARD_load_DB:
    input:
        db=get_card_db_file(),
        read_args=expand(
            "results/{project}/output/ARGs/reads/{sample}/{sample}.gene_mapping_data.txt",
            sample=get_samples(),
            project=get_project(),
        ),
    output:
        touch("results/CARD_load_DB.done"),
    log:
        "logs/CARD_load_DB.log",
    conda:
        "../envs/card.yaml"
    shell:
        "rgi clean --local && "
        "rgi load --card_json {input.db} --local > {log} 2>&1"
        "rgi load --card_json {input.db} --local > {log} 2>&1"


rule CARD_assembly_run:
    input:
        fa=rules.gzip_assembly.output,
        db=rules.CARD_load_DB.output,
    output:
        txt="results/{project}/output/ARGs/assembly/{sample}/{sample}.txt",
        json="results/{project}/output/ARGs/assembly/{sample}/{sample}.json",
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


use rule CARD_read_sample_summary as CARD_assembly_sample_summary with:
    input:
        txt=rules.CARD_assembly_run.output.txt,
        json=rules.CARD_assembly_run.output.json,
    output:
        csv="results/{project}/output/ARGs/assembly/{sample}/{sample}_assembly_ARGs.csv",
    params:
        case="assembly",
    log:
        "logs/{project}/ARGs/assembly/{sample}.log",


use rule CARD_assembly_run as CARD_mag_run with:
    input:
        fa="results/{project}/output/fastas/{sample}/mags/{binID}.fa.gz",
        fa="results/{project}/output/fastas/{sample}/mags/{binID}.fa.gz",
        db=rules.CARD_load_DB.output,
    output:
        txt="results/{project}/output/ARGs/mags/{sample}/{binID}.txt",
        txt="results/{project}/output/ARGs/mags/{sample}/{binID}.txt",
    params:
        path_wo_ext=lambda wildcards, output: Path(output.txt).with_suffix(""),
    log:
        "logs/{project}/ARGs/mags/{sample}/{binID}.log",
        "logs/{project}/ARGs/mags/{sample}/{binID}.log",


rule wrap_mag_ARGs:
    input:
        get_mag_ARGs,
    output:
        "results/{project}/output/ARGs/mags/{sample}/all_mags.done",
        "results/{project}/output/ARGs/mags/{sample}/all_mags.done",
    log:
        "logs/{project}/ARGs/mags/{sample}/all_mags.log",
        "logs/{project}/ARGs/mags/{sample}/all_mags.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "touch {output} > {log} 2>&1"
