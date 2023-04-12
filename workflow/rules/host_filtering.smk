from pathlib import Path

REFERENCE_DIR = get_reference_dir()

## indexing
rule minimap2_index:
    input:
        target = get_reference_file
    output:
        "{0}{{ref}}.mmi".format(REFERENCE_DIR)
    log:
        "logs/minimap2/{ref}.log"
    params:
        extra=""  # optional additional args
    threads: 3
    wrapper:
        "v1.25.0/bio/minimap2/index"

## minimap2 alignment & filtering
rule filter_host_reads:
    input:
        ref_genomes = get_reference_index,
        query = get_trimmed_fastqs,
    output:
        filt_fastq1 = 'results/{project}/filtered/{sample}_final.1.fastq.gz',
        filt_fastq2 = 'results/{project}/filtered/{sample}_final.2.fastq.gz',
    log:
        "logs/{project}/host_filtering/{sample}.log",
    params:
        filt_path = lambda wildcards, output: Path(output.filt_fastq1).parent,
    threads: 3
    conda:
        "../envs/minimap2.yaml"
    script:
        "../scripts/host_filtering.py"


