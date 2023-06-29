import os
from pathlib import Path
from Bio import SeqIO

sample="ERR970400"
contigs = f"results/metacock_test/assembly/{sample}/final.contigs.fa"

with open(contigs, "r") as fcont:
    for record in SeqIO.parse(fcont, "fasta"):
        name = record.name
        k = name[name.find("k") +1 : name.find("_")]
        break
   
out_folder = f"{Path(contigs).parent}/metacock/"
fastg = f"{out_folder}final.fastg"
os.system(f"mkdir {out_folder}")
os.system(f"megahit_core contig2fastg {k} {contigs} > {fastg}")

gfa = f"{out_folder}final.gfa"
print("starting fastg to gfa conversion:")
convert_path = "resources/gfaview/misc/fastg2gfa"
os.system(f"{convert_path} {fastg} > {gfa}")
print("done")

filt_r1=f"results/metacock_test/filtered/bacteria/{sample}_bacteria_1.fastq.gz"
filt_r2=f"results/metacock_test/filtered/bacteria/{sample}_bacteria_2.fastq.gz"
abundance=f"{out_folder}abundance.tsv"
os.system(f"coverm contig -1 {filt_r1} -2 {filt_r2} -r {contigs} -o {abundance} -t 8 && sed -i '1d' {abundance}")

bin_folder = f"results/metacock_test/metacock/{sample}/"
os.system(f"metacoag --assembler megahit --graph {gfa} --contigs {contigs} --abundance {abundance} --output {bin_folder}")