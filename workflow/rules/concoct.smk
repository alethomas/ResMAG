rule concoct_cut_fa:
    input:
        contigs=get_assembly,
    output:
        bed="results/{project}/binning_prep/{sample}/contigs_10K.bed",
        fa="results/{project}/binning_prep/{sample}/contigs_10K.fa",
    threads: 30
    log:
        "logs/{project}/concoct/{sample}/cut_fa.log",
    conda:
        "../envs/concoct.yaml"
    shell:
        "(cut_up_fasta.py {input.contigs} -c 10000 --merge_last "
        "-o 0 -b {output.bed} > {output.fa}) > {log} 2>&1"


rule concoct_coverage:
    input:
        bed=rules.concoct_cut_fa.output.bed,
        bam=rules.map_to_assembly.output.bam,
        bai=rules.reads_mapped_assembly.output.bai,
    output:
        tsv="results/{project}/binning_prep/{sample}/coverage_table.tsv",
    threads: 30
    log:
        "logs/{project}/concoct/{sample}/coverage.log",
    conda:
        "../envs/concoct.yaml"
    shell:
        "(concoct_coverage_table.py {input.bed} "
        "{input.bam} > {output.tsv}) > {log} 2>&1"


rule concoct_run:
    input:
        fa=rules.concoct_cut_fa.output.fa,
        tsv=rules.concoct_coverage.output.tsv,
    output:
        csv="results/{project}/concoct/{sample}/clustering_gt1000.csv",
    params:
        outfolder=lambda wildcards, output: Path(output.csv).parent,
    threads: 64
    log:
        "logs/{project}/concoct/{sample}/run.log",
    conda:
        "../envs/concoct.yaml"
    shell:
        "concoct --composition_file {input.fa} --coverage_file {input.tsv} "
        "--threads {threads} -b {params.outfolder} > {log} 2>&1"


rule concoct_merge:
    input:
        csv=rules.concoct_run.output.csv,
    output:
        merge="results/{project}/concoct/{sample}/clustering_merged.csv",
    threads: 30
    log:
        "logs/{project}/concoct/{sample}/merge.log",
    conda:
        "../envs/concoct.yaml"
    shell:
        "(merge_cutup_clustering.py {input.csv} "
        "> {output.merge}) > {log} 2>&1"


rule concoct_bins_fastas:
    input:
        csv=rules.concoct_merge.output.merge,
        contigs=get_assembly,
    output:
        bins=directory("results/{project}/concoct/{sample}/bins/"),
    threads: 64
    log:
        "logs/{project}/concoct/{sample}/bins_fasta.log",
    conda:
        "../envs/concoct.yaml"
    shell:
        "extract_fasta_bins.py {input.contigs} {input.csv} "
        "--output_path {output.bins} > {log} 2>&1"
