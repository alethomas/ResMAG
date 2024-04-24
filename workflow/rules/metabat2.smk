rule coverm_metabat2:
    input:
        bact_reads=get_filtered_gz_fastqs,
        contigs=get_assembly,
    output:
        abd=temp("results/{project}/binning_prep/{sample}/abundance_metabat.tsv"),
    threads: 2
    log:
        "logs/{project}/coverm_metabat/{sample}.log",
    conda:
        "../envs/coverm.yaml"
    shell:
        "coverm contig -m metabat -1 {input.bact_reads[0]} "
        "-2 {input.bact_reads[1]} -r {input.contigs} "
        "-t {threads} -o {output.abd} > {log} 2>&1"


rule metabat2:
    input:
        contigs=get_assembly,
        abd=rules.coverm_metabat2.output.abd,
    output:
        temp(directory("results/{project}/metabat2/{sample}/")),
    threads: 32
    params:
        prefix="bin",
    log:
        "logs/{project}/metabat2/{sample}.log",
    conda:
        "../envs/metabat2.yaml"
    shell:
        "metabat2 -i {input.contigs} -a {input.abd} "
        "-t {threads} -v -o {output}/{params.prefix} > {log} 2>&1"
