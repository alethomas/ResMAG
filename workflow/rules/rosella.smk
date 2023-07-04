rule rosella_run:
    input:
        fasta="results/{project}/assembly/{sample}/final.contigs_300.fa",
        coverage="results/{project}/assembly/{sample}/abundance.tsv",
        fq1="results/{project}/filtered/bacteria/{sample}_bacteria_1.fastq.gz",
        fq2="results/{project}/filtered/bacteria/{sample}_bacteria_2.fastq.gz",
    output:
        outfile="results/{project}/rosella/{sample}/rosella_res/rosella_kmer_table.tsv",
    params:
        outdir=lambda wildcards, output: Path(output.outfile).parent,
        threads=config["binning"]["threads"],
    log:
        "logs/{project}/rosella/{sample}/rosella_run.log",
    conda:
        "../envs/rosella.yaml"
    shell:
        "bash clone rosella"
        "rosella recover -r {input.fasta} -1 {input.fq1} -2 {input.fq2}  -o {params.outdir} -t {params.threads} 2> {log} "
