import os
import sys

sys.stderr = open(snakemake.log[0], "w")

contigs = snakemake.input.contigs[0]
fastg = snakemake.output[0]


with open(contigs, "r") as fcont:
    line = fcont.readline()
    k = line[line.find("k") + 1 : line.find("_")]

os.system(f"megahit_core contig2fastg {k} {contigs} > {fastg}")
