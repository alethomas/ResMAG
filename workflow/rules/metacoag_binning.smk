rule coverm_metacoag:
    input:
        bact_reads=get_bacterial_reads,
        contigs="results/{project}/assembly/{sample}/final.contigs.fa",
    output:
        "results/{project}/binning_prep/{sample}/abundance.tsv",
    threads: 8
    log:
        "logs/{project}/coverm/{sample}.log",
    conda:
        "../envs/coverm.yaml"
    shell:
        "coverm contig -1 {input.bact_reads[0]} -2 {input.bact_reads[1]} "
        "-r {input.contigs} -o {output} -t {threads}"


rule edit_abundance_file:
    input:
        "results/{project}/binning_prep/{sample}/abundance.tsv",
    output:
        "results/{project}/binning_prep/{sample}/abundance_metacoag.tsv",
    threads: 8
    log:
        "logs/{project}/coverm/{sample}.log",
    conda:
        "../envs/metacoag.yaml"
    shell:
        "cp {input} {output} && sed -i '1d' {output} > {log} 2>&1"


rule fastg_assembly_tree:
    input:
        "results/{project}/assembly/{sample}/final.contigs.fa",
    output:
        "results/{project}/binning_prep/{sample}/assembly_tree.fastg",
    threads: 2
    log:
        "logs/{project}/fastg_assembly_tree/{sample}.log",
    conda:
        "../envs/megahit.yaml"
    script:
        "../scripts/fastg_assembly_tree.py"


rule fastg2gfa:
    input:
        "results/{project}/binning_prep/{sample}/assembly_tree.fastg",
    output:
        "results/{project}/binning_prep/{sample}/assembly_tree.gfa",
    params:
        fastg2gfa_program="workflow/scripts/fastg2gfa",
    threads: 2
    log:
        "logs/{project}/fastg2gfa/{sample}.log",
    conda:
        "../envs/metacoag.yaml"
    shell:
        "{params.fastg2gfa_program} {input} > {output} 2>> {log}"


rule metacoag:
    input:
        contigs="results/{project}/assembly/{sample}/final.contigs.fa",
        gfa="results/{project}/binning_prep/{sample}/assembly_tree.gfa",
        abundance="results/{project}/binning_prep/{sample}/abundance_metacoag.tsv",
    output:
        dir=directory("results/{project}/metacoag/{sample}/"),
        out="results/{project}/metacoag/{sample}/contig_to_bin.tsv",
    params:
        fastg2gfa_program="resources/gfaview/misc/fastg2gfa",
    threads: 2
    log:
        "logs/{project}/metacoag/{sample}.log",
    conda:
        "../envs/metacoag.yaml"
    shell:
        "metacoag --assembler megahit --graph {input.gfa} "
        "--contigs {input.contigs} --abundance {input.abundance} "
        "--output {output.dir} > {log} 2>&1"
