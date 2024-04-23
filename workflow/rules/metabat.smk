rule coverm_metabat:
    input:
        bact_reads=get_filtered_gz_fastqs,
        contigs=get_assembly,
    output:
        "results/{project}/binning_prep/{sample}/{sample}_to_{pathogen}_abundance.tsv",
    threads: 2
    log:
        "logs/{project}/coverm/{sample}_to_{pathogen}.log",
    conda:
        "../envs/coverm.yaml"
    shell:
        "coverm contig -1 {input.bact_reads[0]} -2 {input.bact_reads[1]} "
        "-r {input.contigs} -o {output} -t {threads} > {log} 2>&1"


rule metabat:
    input:
        contigs=rules.gzip_assembly.output,
        abd=rules.coverm_metabat.output,
    output:
        directory("results/{project}/metabat/{sample}/{sample}_to_{pathogen}/"),
    threads: 64
    params:
        prefix="bin",
    log:
        "logs/{project}/metabat/{sample}_to_{pathogen}.log",
    conda:
        "../envs/metabat2.yaml"
    shell:
        "metabat1 -i {input.contigs} -a {input.abd} "
        "--seed 1 --p1 95 --p2 90 "
        "-B 20 --pB 50 --minProb 80 --minBinned 20 "
        "-m 1500 -s 100000 --fuzzy "
        "--minCorr 95 --minContigByCorr 300 --minSamples 1 "
        "-t {threads} -v -l -o {output}/{params.prefix} > {log} 2>&1"
