rule rosella_run:
    input:
        catalogue="results/{project}/rosella/{sample}/catalogue.fna.gz",
        bam="results/{project}/rosella/{sample}/{sample}.bam",
    output:
        outfile="results/{project}/rosella/{sample}/rosella_res/model.pt",
    params:
        outdir=lambda wildcards, output: Path(output.outfile).parent,
        threads=config["binning"]["threads"],
    log:
        "logs/{project}/rosella/{sample}/rosella_run.log",
    conda:
        "../envs/rosella.yaml"
    shell:
        "rosella bin -r {scaffolds.fasta} --coverage-values {coverm.cov} -o {params.outdir} -t {params.threads} "
