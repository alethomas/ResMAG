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
        report="results/{project}/output/classification/reads/{sample}/kaiju.out",
    threads: 64
    log:
        "logs/{project}/kaiju/{sample}_run.log",
    conda:
        "../envs/kaiju.yaml"
    shell:
        "kaiju -z {threads} -t {input.db_files[0]} -f {input.db_files[1]} "
        "-i {input.fastqs[0]} -j {input.fastqs[1]} -o {output.report} > {log} 2>&1"


rule kaiju_table:
    input:
        kaiju_report=rules.run_kaiju.output.report,
        db_files=get_kaiju_files(),
    output:
        summary="results/{project}/output/classification/reads/{sample}/{sample}_{level}_kaiju_summary.tsv",
    threads: 2
    log:
        "logs/{project}/kaiju/{sample}_{level}_table.log",
    conda:
        "../envs/kaiju.yaml"
    shell:
        "kaiju2table -t {input.db_files[0]} -n {input.db_files[2]} "
        "-r {wildcards.level} -o {output.summary} {input.kaiju_report} "
        "> {log} 2>&1"


rule kaiju2krona:
    input:
        db_files=get_kaiju_files(),
        kaiju_report=rules.run_kaiju.output.report,
    output:
        krona=temp(
            "results/{project}/output/classification/reads/{sample}/kaiju.out.krona"
        ),
        html=report(
            "results/{project}/output/report/{sample}/{sample}_kaiju.out.html",
            htmlindex="index.html",
            category="5. Taxonomic classification",
            subcategory="5.2 Read classification",
            labels={"sample": "{sample}"},
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
    log:
        "logs/gtdbtk_download_DB.log",
    conda:
        "../envs/gtdbtk.yaml"
    shell:
        "download-db.sh > {log} 2>&1"


rule gtdbtk_classify_wf:
    input:
        bins=rules.dastool_run.output.bins,
        dastool=rules.dastool_run.output.outdir,
        log="logs/gtdbtk_download_DB.log",
    output:
        summary="results/{project}/output/classification/bins/{sample}/{sample}.bac120.summary.tsv",
        json="results/{project}/output/classification/bins/{sample}/gtdbtk.json",
        outdir=temp(
            directory("results/{project}/output/classification/bins/{sample}/gtdbtk/")
        ),
    params:
        clf_outdir=lambda wildcards, output: Path(output.summary).parent,
        sum_name=lambda wildcards, output: Path(output.summary).name,
        json=lambda wildcards, output: Path(output.json).name,
    threads: 30
    log:
        "logs/{project}/gtdbtk/{sample}_classify.log",
    conda:
        "../envs/gtdbtk.yaml"
    shell:
        "(gtdbtk classify_wf --skip_ani_screen --prefix {wildcards.sample} "
        "-x fa --cpus {threads} --genome_dir {input.bins}/ --out_dir {output.outdir}/ && "
        "cp {output.outdir}/classify/{params.sum_name} {params.clf_outdir}/ && "
        "cp {output.outdir}/{params.json} {params.clf_outdir}/) > {log} 2>&1"
