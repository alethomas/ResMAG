# when assembly.gz & "results/{project}/das_tool/{sample}/{sample}_DASTool_summary.tsv" are there
# rm results/strawpigs/megahit/
# rm results/strawpigs/metacoag/
# rm results/strawpigs/metabat2
# rm results/strawpigs/rosella
# rm results/strawpigs/metabinner
# temp(directory(...))
# mv dastool contig2bin file into contgi2bin folder
# save dastool summary & binner control file; remove rest of das tool folder
# rm binning_prep folder

rule remove_megahit_intermediates:
    input:
        contigs=get_assembly,
    output:
        touch("logs/{project}/assembly/{sample}_intermediate_removal.done"),
    params:
        outdir=lambda wildcards, input: Path(input.contigs).parent,
    log:
        "logs/{project}/assembly/{sample}_intermediate_removal.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "(find {params.outdir}/ -mindepth 1 -type d -exec rm -rf {{}} +) > {log} 2>&1"

