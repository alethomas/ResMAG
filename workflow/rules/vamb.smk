rule vamb_contig_catalogue:
    input:
        contigs=expand("results/{{project}}/assembly/{sample}/final.contigs.fa",sample=get_samples()),
    output:
        catalogue="results/{project}/vamb/catalogue.fna.gz"
    params:
        threshold=config["binning"]["min_contig_length"]
    log:
        "logs/{project}/vamb/contig_catalogue.log",        
    conda:
        "../envs/vamb.yaml"
    shell:
        "python $CONDA_PREFIX/bin/concatenate.py -m {params.threshold} {output.catalogue} {input.contigs}"

rule vamb_catalogue_index:
    input:
        catalogue="results/{project}/vamb/catalogue.fna.gz"
    output:
        index="results/{project}/vamb/catalogue.mmi"
    log:
        "logs/{project}/vamb/contig_catalogue_index.log",       
    conda:
        "../envs/vamb.yaml"
    shell:
        "minimap2 -d {output.index} {input.catalogue} 2>{log}"        

rule vamb_map_reads:
    input:
        f1="results/{project}/trimmed/fastp/{sample}.1.fastq.gz",
        f2="results/{project}/trimmed/fastp/{sample}.2.fastq.gz",        
        index="results/{project}/vamb/catalogue.mmi",
    output:
        bam=temp("results/{project}/vamb/{sample}/{sample}.bam")
    params:
        threads=config["binning"]["threads"]
    log:
        "logs/{project}/vamb/{sample}/map_reads.log",       
    conda:
        "../envs/vamb.yaml"
    shell:
        "minimap2 -t {params.threads} -N 5 -ax sr {input.index} "
        "{input.f1} {input.f2} | samtools view -F 3584 -b --threads 8 > {output.bam} 2>{log}"

rule vamb_run:
    input:
        catalogue="results/{project}/vamb/catalogue.fna.gz",
        bam=expand("results/{{project}}/vamb/{sample}/{sample}.bam",sample=get_samples()),
    output:
        "results/{project}/vamb/{sample}/vamb_res/model.pt"
    params:
        outdir="results/{project}/vamb/{sample}/vamb_res"
    log:
        "logs/{project}/vamb/{sample}/vamb_run.log",       
    conda:
        "../envs/vamb.yaml"    
    shell:
        "rm -r {params.outdir}; "
        "vamb --outdir {params.outdir} --fasta {input.catalogue} --bamfiles {input.bam} 2>{log}"
