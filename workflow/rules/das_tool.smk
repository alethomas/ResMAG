rule postprocess_vamb:
    input:
        "results/{project}/vamb/{sample}/vamb_res/clusters.tsv",
    output:
        "results/{project}/output/contig2bins/{sample}/vamb_contig2bin.tsv",
    log:
        "logs/{project}/contig2bins/{sample}/postprocess_vamb.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "(awk '{{print $2 \"\tvamb_bin_\" $1}}' {input} > {output} && "
        "sed -i 's/S1C//g' {output}) > {log} 2>&1"


rule postprocess_metabat:
    input:
        "results/{project}/metabat2/{sample}/",
    output:
        "results/{project}/output/contig2bins/{sample}/metabat2_contig2bin.tsv",
    params:
        binner="metabat2",
    log:
        "logs/{project}/contig2bins/{sample}/postprocess_metabat2.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/binner_postprocess.py"


use rule postprocess_metabat as postprocess_rosella with:
    input:
        "results/{project}/rosella/{sample}/",
    output:
        "results/{project}/output/contig2bins/{sample}/rosella_contig2bin.tsv",
    params:
        binner="rosella",
    log:
        "logs/{project}/contig2bins/{sample}/postprocess_rosella.log",


rule postprocess_metacoag:
    input:
        c2bin="results/{project}/metacoag/{sample}/contig_to_bin.tsv",
        folder="results/{project}/metacoag/{sample}/",
    output:
        "results/{project}/output/contig2bins/{sample}/metacoag_contig2bin.tsv",
    log:
        "logs/{project}/contig2bins/{sample}/postprocess_metacoag.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "(awk '{{print $1 \"\t\" $NF}}' {input.c2bin} | "
        "sed 's/len=[0-9]*,/metacoag_/g' > {output}) > {log} 2>&1"


## move metabinner output to contig2bin folder
rule postprocess_metabinner:
    input:
        "results/{project}/metabinner/{sample}/metabinner_res/metabinner_result.tsv",
    output:
        "results/{project}/output/contig2bins/{sample}/metabinner_contig2bin.tsv",
    log:
        "logs/{project}/contig2bins/{sample}/postprocess_metabinner.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "awk '{{print $1 \"\tmetabinner_bin_\" $2}}' {input} > {output} 2> {log}"


## tests if binner created bins and writes file for DAS Tool
rule binner_control:
    input:
        get_all_contig2bin_files,
    output:
        temp("results/{project}/das_tool/{sample}_binner_control.csv"),
    params:
        get_binners(),
    log:
        "logs/{project}/das_tool/{sample}/binner_control.log",
    conda:
        "../envs/unix.yaml"
    script:
        "../scripts/binner_control.py"


rule dastool_run:
    input:
        binc="results/{project}/das_tool/{sample}_binner_control.csv",
        contigs=get_assembly,
        asml_folder=rules.megahit.output.outdir,
    output:
        summary="results/{project}/das_tool/{sample}/{sample}_DASTool_summary.tsv",
        contig2bin="results/{project}/das_tool/{sample}/{sample}_DASTool_contig2bin.tsv",
        bins=temp(
            directory("results/{project}/das_tool/{sample}/{sample}_DASTool_bins/")
        ),
        outdir=temp(directory("results/{project}/das_tool/{sample}/")),
    params:
        path_bin_list=get_paths_binner,
    threads: get_DAS_Tool_threads()
    log:
        "logs/{project}/das_tool/{sample}/das_tool_run.log",
    conda:
        "../envs/das_tool.yaml"
    shell:
        "DAS_Tool -t {threads} "
        "--debug --write_bins "
        "-i {params.path_bin_list[0]} "
        "-l {params.path_bin_list[1]} "
        "-c {input.contigs} "
        "-o {output.outdir}/{wildcards.sample} "
        "> {log} 2>&1 "


if bins_for_sample:

    rule move_dastool_output:
        input:
            contig2bin="results/{project}/das_tool/{sample}/{sample}_DASTool_contig2bin.tsv",
            summary="results/{project}/das_tool/{sample}/{sample}_DASTool_summary.tsv",
            # so folder is not removed before copying
            dastool=rules.dastool_run.output.outdir,
        output:
            contig2bin="results/{project}/output/contig2bins/{sample}/DASTool_contig2bin.tsv",
            summary="results/{project}/output/report/{sample}/{sample}_DASTool_summary.tsv",
        threads: 2
        log:
            "logs/{project}/das_tool/{sample}/move_output.log",
        conda:
            "../envs/unix.yaml"
        shell:
            "(cp {input.contig2bin} {output.contig2bin} && "
            "cp {input.summary} {output.summary}) > {log} 2>&1 "

    rule gzip_bins:
        input:
            bins=rules.dastool_run.output.bins,
            dastool=rules.dastool_run.output.outdir,
        output:
            bins=directory("results/{project}/output/fastas/{sample}/bins/"),
        threads: 2
        log:
            "logs/{project}/bins/{sample}/gz_bins.log",
        conda:
            "../envs/unix.yaml"
        shell:
            "(gzip -k {input.bins}/*.fa && "
            "mkdir -p {output.bins}/ && "
            "mv {input.bins}/*.fa.gz {output.bins}/ ) > {log} 2>&1"

    rule move_MAGs:
        input:
            bin_dir=rules.gzip_bins.output.bins,
            csv="results/{project}/output/report/{sample}/{sample}_mags_summary.csv",
        output:
            outdir=directory("results/{project}/output/fastas/{sample}/mags/"),
        threads: 2
        log:
            "logs/{project}/bins/{sample}/move_MAGs.log",
        conda:
            "../envs/python.yaml"
        script:
            "../scripts/move_MAGs.py"
