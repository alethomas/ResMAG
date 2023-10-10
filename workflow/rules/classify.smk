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