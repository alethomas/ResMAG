rule download_kaiju:
    output:
        db_files=get_kaiju_files(),
    params:
        download=config["kaiju"]["download"],
        db_folder=lambda wildcards, output: Path(output.db_files[0]).parent,
    log:
        "logs/kaiju_DB_download.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "(mkdir -p {params.db_folder} && "
        "wget -c {params.download} -O - | "
        "tar -zxv -C {params.db_folder}) > {log} 2>&1"


rule run_kaiju:
    input:
        db_files=get_kaiju_files(),
        fastqs=get_filtered_gz_fastqs,
    output:
        "results/{project}/kaiju/{sample}/kaiju.out",
    threads: 64
    log:
        "logs/{project}/kaiju/{sample}_run.log",
    conda:
        "../envs/kaiju.yaml"
    shell:
        "kaiju -z {threads} -t {input.db_files[0]} -f {input.db_files[1]} "
        "-i {input.fastqs[0]} -j {input.fastqs[1]} -o {output} > {log} 2>&1"


rule kaiju2krona:
    input:
        db_files=get_kaiju_files(),
        kaiju_report="results/{project}/kaiju/{sample}/kaiju.out",
    output:
        krona=temp("results/{project}/kaiju/{sample}/kaiju.out.krona"),
        html=report(
            "results/{project}/report/{sample}/kaiju.out.html",
            htmlindex="index.html",
            category="5. Taxonomic classification",
            subcategory="5.2 Read classification",
            labels={"sample": "{sample}", "type": "view"},
        ),
        #Kaiju read classification krona plot
    log:
        "logs/{project}/kaiju/{sample}_2krona.log",
    conda:
        "../envs/kaiju.yaml"
    shell:
        "(kaiju2krona -t {input.db_files[0]} -n {input.db_files[2]} "
        "-i {input.kaiju_report} -o {output.krona} && "
        "ktImportText -o {output.html} {output.krona}) > {log} 2>&1"


rule download_GTDB:
    output:
        "logs/gtdbtk_download_DB.log",
    log:
        "logs/gtdbtk_download_DB.log",
    conda:
        "../envs/gtdbtk.yaml"
    shell:
        "download-db.sh > {log} 2>&1"


##remove output here, log can be used as input too


rule gtdbtk_classify_wf:
    input:
        bins="results/{project}/das_tool/{sample}/{sample}_DASTool_bins/",
        log="logs/gtdbtk_download_DB.log",
    output:
        summary="results/{project}/classification/{sample}/{sample}.bac120.summary.tsv",
    params:
        outdir=lambda wildcards, output: Path(output.summary).parent,
    threads: 12
    log:
        "logs/{project}/gtdbtk/{sample}_classify.log",
    conda:
        "../envs/gtdbtk.yaml"
    shell:
        "gtdbtk classify_wf --skip_ani_screen --prefix {wildcards.sample} "
        "-x fa --cpus {threads} --genome_dir {input.bins} --out_dir {params.outdir} "
        "> {log} 2>&1"
