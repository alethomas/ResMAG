rule binner_control:
    input:
        bin1="results/{project}/metabinner/{sample}/metabinner_res/metabinner_result.tsv",
        bin2="results/{project}/vamb/{sample}/vamb_res/swaped_clusters.tsv",
        bin3="results/{project}/metabat2/{sample}/metabat_contig2bin.tsv",
        bin4="results/{project}/metacoag/{sample}/contig_to_bin.tsv",
        bin5="results/{project}/rosella/{sample}/rosella_contig2bin.tsv",
    output:
        "results/{project}/das_tool/{sample}/binner_control.csv",
    log:
        "logs/{project}/das_tool/{sample}/binner_control.log",
    conda:
        "../envs/das_tool.yaml"
    script:
        "workflow/scripts/binner_control.py"

rule postprocess_vamb:
    input:
        "results/{project}/vamb/{sample}/vamb_res/model.pt",
    output:
        "results/{project}/vamb/{sample}/vamb_res/swaped_clusters.tsv",
    params:
        cluster="results/{project}/vamb/{sample}/vamb_res/clusters.tsv",
        tmp="results/{project}/vamb/{sample}/vamb_res/temp",
        root_path=get_root(),
    log:
        "logs/{project}/das_tool/{sample}/vamb_postprocessing.log",
    conda:
        "../envs/das_tool.yaml"
    shell:
        "bash {params.root_path}/workflow/scripts/swap_vamb_clusters.sh {params.cluster} {output} {params.tmp} 2> {log}"

rule postprocess_metabat:
    input:
        "results/strawpigs/metabat2/{sample}/"
    output:
        "results/{project}/binning_rev/{sample}/metabat_contig2bin.tsv",
    log:
        "logs/{project}/binning_rev/{sample}/metabat_postprocessing.log",
    conda:
        "../envs/das_tool.yaml"
    script:
        "workflow/scripts/postprocess_metabat.py"

rule postprocess_rosella:
    input:
        "results/{project}/rosella/{sample}/rosella_bins.json",
    output:
        "results/{project}/rosella/{sample}/rosella_contig2bin.tsv",
    log:
        "logs/{project}/das_tool/{sample}/rosella_postprocessing.log",
    conda:
        "../envs/das_tool.yaml"
    script:
        "workflow/scripts/postprocess_rosella.py"

rule dastool_run:
    input:
        bin1="results/{project}/metabinner/{sample}/metabinner_res/metabinner_result.tsv",
        bin2="results/{project}/vamb/{sample}/vamb_res/swaped_clusters.tsv",
        bin3="results/{project}/metabat2/{sample}/metabat_contig2bin.tsv",
        bin4="results/{project}/metacoag/{sample}/contig_to_bin.tsv",
        bin5="results/{project}/rosella/{sample}/rosella_contig2bin.tsv",
        binc="results/{project}/das_tool/{sample}/binner_control.csv",
        contigs="results/{project}/assembly/{sample}/final.contigs.fa",
    output:
        outfile="results/{project}/das_tool/{sample}/{sample}_DASTool_summary.tsv",
    params:
        paths=lambda wildcards, input: get_paths_binner(input.binc)[0],
        binner=lambda wildcards, input: get_paths_binner(input.binc)[1],
        outdir=lambda wildcards, output: Path(output.outfile).parent,
    log:
        "logs/{project}/das_tool/{sample}/das_tool_run.log",
    conda:
        "../envs/das_tool.yaml"
    shell:
        "DAS_Tool -t 64 "
        "--debug "
        "-i {params.paths} "
        "-l {params.binner} "
        "-c {input.contigs} "
        "-o {params.outdir}/{wildcard.sample} "
        " 2> {log} "
