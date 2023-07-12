rule binner_control:
    input:
        get_all_contig2bin,
    output:
        "results/{project}/das_tool/{sample}/binner_control.csv",
    log:
        "logs/{project}/das_tool/{sample}/binner_control.log",
    conda:
        "../envs/das_tool.yaml"
    script:
        "../scripts/binner_control.py"


rule postprocess_vamb:
    input:
        "results/{project}/vamb/{sample}/vamb_res/model.pt",
    output:
        "results/{project}/binning_rev/{sample}/vamb_contig2bin.tsv",
    params:
        cluster="results/{project}/vamb/{sample}/vamb_res/clusters.tsv",
        tmp="results/{project}/vamb/{sample}/vamb_res/temp",
        root=get_root(),
    log:
        "logs/{project}/binning_rev/{sample}/postprocess_vamb.log",
    conda:
        "../envs/das_tool.yaml"
    shell:
        "bash {params.root}/workflow/scripts/postprocess_vamb.sh {params.cluster} {output} {params.tmp} 2> {log}"


rule postprocess_metabat:
    input:
        "results/strawpigs/metabat2/{sample}/",
    output:
        "results/{project}/binning_rev/{sample}/metabat_contig2bin.tsv",
    log:
        "logs/{project}/binning_rev/{sample}/postprocess_metabat.log",
    conda:
        "../envs/das_tool.yaml"
    script:
        "../scripts/postprocess_metabat.py"


rule postprocess_rosella:
    input:
        "results/{project}/rosella/{sample}/",
    output:
        "results/{project}/binning_rev/{sample}/rosella_contig2bin.tsv",
    log:
        "logs/{project}/binning_rev/{sample}/postprocess_rosella.log",
    conda:
        "../envs/das_tool.yaml"
    script:
        "../scripts/postprocess_metabat.py"


rule dastool_run:
    input:
        all_bins=get_all_contig2bin,
        binc="results/{project}/das_tool/{sample}/binner_control.csv",
        contigs="results/{project}/assembly/{sample}/final.contigs.fa",
    output:
        outfile="results/{project}/das_tool/{sample}/{sample}_DASTool_summary.tsv",
    params:
        path_bin_list=get_paths_binner,
        outdir=lambda wildcards, output: Path(output.outfile).parent,
    log:
        "logs/{project}/das_tool/{sample}/das_tool_run.log",
    conda:
        "../envs/das_tool.yaml"
    shell:
        "DAS_Tool -t 64 "
        "--debug "
        "-i {params.path_bin_list[0]} "
        "-l {params.path_bin_list[1]} "
        "-c {input.contigs} "
        "-o {params.outdir}/{wildcards.sample} "
        "> {log} 2>&1 "
