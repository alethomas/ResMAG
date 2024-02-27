from pathlib import Path


rule metabinner_coverage_profile:
    input:
        contig_file=get_assembly,
        fastqs=get_filtered_fastqs,
    output:
        outfile=temp(
            "results/{project}/metabinner/{sample}/coverage_profile_{threshold}/coverage_profile.tsv"
        ),
    threads: 32
    params:
        threads=config["binning"]["threads"],
        outdir=lambda wildcards, output: Path(output.outfile).parent,
    log:
        "logs/{project}/metabinner/{sample}/coverage_profile_{threshold}.log",
    conda:
        "../envs/metabinner_env.yaml"
    shell:
        "(bash $CONDA_PREFIX/bin/scripts/gen_coverage_file.sh "
        "-t {params.threads} "
        "-a {input.contig_file} "
        "-o {params.outdir} "
        "-l {wildcards.threshold} "
        "{input.fastqs[0]} {input.fastqs[1]}) > {log} 2>&1"


rule remove_metabinner_cov_overload:
    input:
        file="results/{project}/metabinner/{sample}/coverage_profile_{threshold}/coverage_profile.tsv",
    output:
        touch("logs/{project}/metabinner/{sample}/coverage_profile_{threshold}_removal.done"),
    params:
        outdir=lambda wildcards, input: Path(input.file).parent,
        filename=lambda wildcards, input: Path(input.file).name,
    log:
        "logs/{project}/metabinner/{sample}/coverage_profile_{threshold}_removal.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "(find {params.outdir}/ -mindepth 1 -type d -exec rm -rf {{}} + && "
        "find {params.outdir}/ -type f ! -name '{params.filename}' "
        "-exec rm -f {{}} +) > {log} 2>&1"


rule metabinner_composition_profile:
    input:
        contigs=get_assembly,
    output:
        outfile=temp(
            "results/{project}/metabinner/{sample}/composition_profile/final.contigs_kmer_{kmer_size}_f{threshold}.csv"
        ),
    threads: 32
    params:
        script_path="{}/workflow/scripts/gen_kmer.py".format(get_root()),  ## does it work when we change this to ../scripts?, AD what is changed?
        indir=lambda wildcards, input: Path(input.contigs).parent,
        outdir=lambda wildcards, output: Path(output.outfile).parent,
        filename=lambda wildcards, output: Path(output.outfile).name,
    log:
        "logs/{project}/metabinner/{sample}/composition_profile_kmer_{kmer_size}_f{threshold}.log",
    conda:
        "../envs/metabinner_env.yaml"
    shell:
        "(python {params.script_path} {input.contigs} {wildcards.threshold} {wildcards.kmer_size} && "
        "mv {params.indir}/{params.filename} {params.outdir}) > {log} 2>&1"


use rule remove_metabinner_cov_overload as remove_metabinner_comp_overload with:
    input:
        file="results/{project}/metabinner/{sample}/composition_profile/final.contigs_kmer_{kmer_size}_f{threshold}.csv",
    output:
        touch("logs/{project}/metabinner/{sample}/composition_profile_kmer_{kmer_size}_f{threshold}_removal.done"),
    log:
        "logs/{project}/metabinner/{sample}/composition_profile_kmer_{kmer_size}_f{threshold}_removal.log",


rule metabinner_run:
    input:
        contig_file=get_assembly,
        coverage_profile=expand(
            "results/{{project}}/metabinner/{{sample}}/coverage_profile_{threshold}/coverage_profile.tsv",
            threshold=get_contig_length_filter(),
        ),
        kmer_profile=expand(
            "results/{{project}}/metabinner/{{sample}}/composition_profile/final.contigs_kmer_{kmer_size}_f{threshold}.csv",
            kmer_size=get_kmersize(),
            threshold=get_contig_length_filter(),
        ),
    output:
        outfile=temp(
            "results/{project}/metabinner/{sample}/metabinner_res/metabinner_result.tsv"
        ),
        outfolder=
            directory("results/{project}/metabinner/{sample}/metabinner_res/")
        ,
    threads: 64
    params:
        threads=config["binning"]["threads"],
        outdir=lambda wildcards, output: Path(output.outfile).parent.parent,
        root=get_root(),
    log:
        "logs/{project}/metabinner/{sample}/metabinner.log",
    conda:
        "../envs/metabinner_env.yaml"
    shell:
        "run_metabinner.sh "
        "-a {params.root}/{input.contig_file} "
        "-d {params.root}/{input.coverage_profile} "
        "-k {params.root}/{input.kmer_profile} "
        "-o {params.root}/{params.outdir} "
        "-p $CONDA_PREFIX/bin/ "
        "-t {params.threads} > {log} 2>&1"


use rule remove_metabinner_cov_overload as remove_metabinner_overload with:
    input:
        file="results/{project}/metabinner/{sample}/metabinner_res/metabinner_result.tsv",
    output:
        touch("logs/{project}/metabinner/{sample}/metabinner_removal.done"),
    log:
        "logs/{project}/metabinner/{sample}/metabinner_removal.log",