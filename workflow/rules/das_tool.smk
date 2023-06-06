rule swap_vamb_res:
    input: 
        "results/{project}/vamb/{sample}/vamb_res/model.pt"
    output: 
        "results/{project}/vamb/{sample}/vamb_res/swaped_clusters.tsv"
    params: 
        cluster="results/{project}/vamb/{sample}/vamb_res/clusters.tsv"
    log:
        "logs/{project}/das_tool/{sample}/swap_vamb_res.log",
    shell:
        "awk '{{ print $2 " " $1}}' {params.cluster} > {output} 2>{log}"

rule rename_vamb_res:
    input: 
        "results/{project}/vamb/{sample}/vamb_res/swaped_clusters.tsv"
    output: 
        "results/{project}/vamb/{sample}/vamb_res/renamed_swaped_clusters.tsv"
    log:
        "logs/{project}/das_tool/{sample}/rename_vamb_res.log",
    shell:
        "sed 's/S[12]C//g' {input} > {output} 2>{log}"

rule dastool_run:
    input:
        bin1="results/{project}/metabinner/{sample}/metabinner_res/metabinner_result.tsv",
        bin2="results/{project}/vamb/{sample}/vamb_res/swaped_clusters.tsv",
        # bins=get_binners_bins,
        # contigs=get_binners_contigs,
    output:
        "results/{project}/das_tool/{sample}/{sample}_DASTool_summary.tsv"
    params:
        contigs="results/strawpigs/assembly/I15546-L1/final.contigs_500.fa",
        threads=config["binning"]["threads"],
        outdir="results/{project}/das_tool/{sample}/{sample}",
    log:
        "logs/{project}/das_tool/{sample}/das_tool_run.log",       
    conda:
        "../envs/das_tool.yaml"    
    shell:
        "DAS_Tool --debug "
        "-i {input.bin1},{input.bin2} "
        "-l metabinner,vamb "   #, metacoag  
        "-c {params.contigs} "
        "-o {params.outdir} "
        "-t {params.threads} "
        "--write_bins 2> {log} "