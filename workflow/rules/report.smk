rule snakemake_report:
    input:
        # 1. Quality control
        "results/{project}/qc/multiqc.html",
        # 2. Species diversity
        "results/{project}/report/kraken2/",
        "results/{project}/report/bracken_plot.png",
        expand(
            "results/{{project}}/report/{sample}/kraken.krona.html",
            sample=get_samples(),
        ),
        # 3. Assembly results
        "results/{project}/report/assembly/",
        # 4. Binning results
        expand(
            "results/{{project}}/report/{sample}/bin/",
            sample=get_samples(),
        ),
        expand(
            "results/{{project}}/report/{sample}/checkm2/",
            sample=get_samples(),
        ),
        expand(
            "results/{{project}}/report/{sample}/dastool/",
            sample=get_samples(),
        ),
        expand(
            "results/{{project}}/report/{sample}/taxonomy/",
            sample=get_samples(),
        ),
        # 5. Taxonomic classification
        expand(
            "results/{{project}}/report/{sample}/mags/",
            sample=get_samples(),
        ),
        expand(
            "results/{{project}}/report/{sample}/kaiju.out.html",
            sample=get_samples(),
        ),
    output:
        "results/{project}/report/report.zip",  #html",
    params:
        style="resources/report/custom-stylesheet.css",
    #    for_testing=get_if_testing("--snakefile ../workflow/Snakefile"),
    log:
        "logs/{project}/snakemake-report.log",
    conda:
        "../envs/snakemake.yaml"
    shell:
        "snakemake --nolock --report {output} --report-stylesheet {params.style} "

        "> {log} 2>&1"
        #"{params.for_testing} "


"""
expand(
            "results/{{project}}/report/{sample}/bin_summary.html",
            sample=get_samples(),
        ),
        expand(
            "results/{{project}}/report/{sample}/bin_summary.csv",
            sample=get_samples(),
        ),
        expand(
            "results/{{project}}/report/{sample}/checkm2_summary.html",
            sample=get_samples(),
        ),
        expand(
            "results/{{project}}/report/{sample}/DASTool_summary.html",
            sample=get_samples(),
        ),
        expand(
            "results/{{project}}/report/{sample}/checkm2_summary.csv",
            sample=get_samples(),
        ),
        expand(
            "results/{{project}}/report/{sample}/DASTool_summary.csv",
            sample=get_samples(),
        ),
        expand(
            "results/{{project}}/report/{sample}/bin_taxonomy.html",
            sample=get_samples(),
        ),
        expand(
            "results/{{project}}/report/{sample}/bin_taxonomy.csv",
            sample=get_samples(),
        ),expand(
            "results/{{project}}/report/{sample}/mags_summary.html",
            sample=get_samples(),
        ),
        expand(
            "results/{{project}}/report/{sample}/mags_summary.csv",
            sample=get_samples(),
        ),


        "results/{project}/qc/multiqc.html",
        expand("results/{{project}}/report/{sample}/kaiju.out.html",sample=get_samples(),
        ),
        "results/{project}/report/assembly_summary.csv",
        "results/{project}/report/assembly_summary.html",
        expand("results/{{project}}/report/{sample}/mags_summary.csv",sample=get_samples(),
        ),
        expand("results/{{project}}/report/{sample}/bin_summary.csv",sample=get_samples(),
        ),
        expand("results/{{project}}/report/{sample}/bin_taxonomy.csv",sample=get_samples(),
        ),
        expand("results/{{project}}/report/{sample}/checkm2_summary.csv",sample=get_samples(),
        ),
        expand("results/{{project}}/report/{sample}/DASTool_summary.csv",sample=get_samples(),
        ),
        expand("results/{{project}}/report/{sample}/mags_summary.html",sample=get_samples(),
        ),
        expand("results/{{project}}/report/{sample}/bin_summary.html",sample=get_samples(),
        ),
        expand("results/{{project}}/report/{sample}/bin_taxonomy.html",sample=get_samples(),
        ),
        expand("results/{{project}}/report/{sample}/checkm2_summary.html",sample=get_samples(),
        ),
        expand("results/{{project}}/report/{sample}/DASTool_summary.html",sample=get_samples(),
        ),"""
