rule fastqc:
    input:
        get_fastqs
    output:
        html="results/{project}/qc/fastqc/{sample}.html",
        zip="results/{project}/qc/fastqc/{sample}_fastqc.zip",
    log:
        "logs/{project}/fastqc/{sample}.log"
    wrapper:
        "v1.23.5/bio/fastqc"

rule multiqc:
    input:
        expand_samples_for_project(
            [
                "results/{{project}}/qc/fastqc/{sample}_fastqc.zip",
                "results/{{project}}/trimmed/fastp/{sample}.fastp.json",
            ]
        ),
    output:
        "results/{project}/qc/multiqc.html",
    params:
        extra=(
            "--config config/multiqc_config.yaml --title 'Results for data from {project} project'"
        ),
    log:
        "logs/{project}/multiqc.log",
    wrapper:
        "v1.23.5/bio/multiqc"