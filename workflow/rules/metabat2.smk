rule coverm_metabat:
    input:
        bact_reads=get_bacterial_reads,
        contigs="results/{project}/assembly/{sample}/final.contigs.fa",
    output:
        "results/{project}/binning_prep/{sample}/abundance_metabat.tsv",
    threads: 2
    log:
        "logs/{project}/coverm_metabat/{sample}.log",
    conda:
        "../envs/coverm.yaml"
    shell:
        "coverm contig -m metabat -1 {input.bact_reads[0]} "
        "-2 {input.bact_reads[1]} -r {input.contigs} "
        "-t {threads} -o {output} > {log} 2>&1"


rule metabat2:
    input:
        contigs="results/{project}/assembly/{sample}/final.contigs.fa",
        abd="results/{project}/binning_prep/{sample}/abundance_metabat.tsv",
    output:
        directory("results/{project}/metabat2/{sample}"),
    threads: 8
    params:
        prefix="bin",
    log:
        "logs/{project}/metabat2/{sample}.log",
    conda:
        "../envs/metabat2.yaml"
    shell:
        "metabat2 -i {input.contigs} -a {input.abd} "
        "-t {threads} -v -o {output}/{params.prefix} > {log} 2>&1"