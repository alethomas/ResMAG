rule swap_vamb_res:
    input:
        "results/{project}/vamb/{sample}/vamb_res/model.pt",
    output:
        "results/{project}/vamb/{sample}/vamb_res/swaped_clusters.tsv",
    params:
        cluster="results/{project}/vamb/{sample}/vamb_res/clusters.tsv",
        tmp="results/{project}/vamb/{sample}/vamb_res/temp",
        root_path=config["root_path"],
    log:
        "logs/{project}/das_tool/{sample}/swap_vamb_res.log",
    conda:
        "../envs/das_tool.yaml"
    shell:
        "bash {params.root_path}/workflow/scripts/swap_vamb_clusters.sh {params.cluster} {output} {params.tmp} 2> {log}"


## TODO -l adding additional binner
rule dastool_run:
    input:
        bin1="results/{project}/metabinner/{sample}/metabinner_res/metabinner_result.tsv",
        bin2="results/{project}/vamb/{sample}/vamb_res/swaped_clusters.tsv",
    output:
        outfile="results/{project}/das_tool/{sample}/{sample}_DASTool_summary.tsv",
    params:
        contigs="results/strawpigs/assembly/{sample}/final.contigs.fa",
        threads=config["binning"]["threads"],
        prefix="{sample}",
        outdir=lambda wildcards, output: Path(output.outfile).parent,
    log:
        "logs/{project}/das_tool/{sample}/das_tool_run.log",
    conda:
        "../envs/das_tool.yaml"
    shell:
        "DAS_Tool --debug "
        "-i {input.bin1},{input.bin2} "
        "-l metabinner,vamb "
        "-c {params.contigs} "
        "-o {params.outdir}/{params.prefix} "
        "-t {params.threads} "
        "--write_bins 2> {log} "
