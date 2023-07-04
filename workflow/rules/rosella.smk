rule rosella_run:
    input:
<<<<<<< HEAD
        fasta="results/{project}/assembly/{sample}/final.contigs_300.fa",
        coverage="results/{project}/assembly/{sample}/abundance.tsv",
        fq1="results/{project}/filtered/bacteria/{sample}_bacteria_1.fastq.gz",
        fq2="results/{project}/filtered/bacteria/{sample}_bacteria_2.fastq.gz",
    output:
        outfile="results/{project}/rosella/{sample}/rosella_res/rosella_kmer_table.tsv",
=======
        catalogue="results/{project}/rosella/{sample}/catalogue.fna.gz",
        bam="results/{project}/rosella/{sample}/{sample}.bam",
    output:
        outfile="results/{project}/rosella/{sample}/rosella_res/model.pt",
>>>>>>> 01415cf53a2b888037da8d7ee38e55670005b80c
    params:
        outdir=lambda wildcards, output: Path(output.outfile).parent,
        threads=config["binning"]["threads"],
    log:
        "logs/{project}/rosella/{sample}/rosella_run.log",
    conda:
        "../envs/rosella.yaml"
    shell:
<<<<<<< HEAD
        "bash clone rosella"
        "rosella recover -r {input.fasta} -1 {input.fq1} -2 {input.fq2}  -o {params.outdir} -t {params.threads} 2> {log} "
=======
        "rosella bin -r {scaffolds.fasta} --coverage-values {coverm.cov} -o {params.outdir} -t {params.threads} "
>>>>>>> 01415cf53a2b888037da8d7ee38e55670005b80c
