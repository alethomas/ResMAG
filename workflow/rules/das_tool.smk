rule dastool_run:
    input:
        bins="results/{project}/metabinner/{sample}/metabinner_res/metabinner_result.tsv",
        contigs="results/{project}/metabinner/{sample}/coverage_profile/work_files/assembly.fa",
        # contigs=expand("results/{{project}}/das_tool/{sample}/{sample}.bam",sample=get_samples()),
    output:
        "results/{project}/das_tool/{sample}/{sample}_DASTool_summary.tsv"
    params:
        threads=config["binning"]["threads"],
        outdir="results/{project}/das_tool/{sample}/{sample}"
    log:
        "logs/{project}/das_tool/{sample}/das_tool_run.log",       
    conda:
        "../envs/das_tool.yaml"    
    shell:
        "DAS_Tool -i {input.bins} "
        "-l metabinner " # vamb, metacoag  
        "-c {input.contigs} "
        "-o {params.outdir} "
        "-t {params.threads} "
        "--write_bins --debug 2>{log}"