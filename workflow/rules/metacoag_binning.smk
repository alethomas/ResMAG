rule coverm:
    input:
        bact_fq1="results/{project}/filtered/bacteria/{sample}_bacteria_1.fastq.gz",
        bact_fq2="results/{project}/filtered/bacteria/{sample}_bacteria_2.fastq.gz",
        contigs="results/{project}/assembly/{sample}/final.contigs.fa",
    output:
        "results/{project}/assembly/{sample}/abundance.tsv",
    threads: 8
    log:
        "logs/{project}/coverm/{sample}.log",
    conda:
        "../envs/metacoag.yaml"
    shell:
        "coverm contig -1 {input.bact_fq1} -2 {input.bact_fq2} "
        "-r {input.contigs} -o {output} -t {threads} && "
        "sed -i '1d' {output} > {log} 2>&1"
