#import command
import sys
import os
from os import path

sys.stderr = open(snakemake.log[0], "w")

def run_command(cmd):
    cmd = cmd.split(" ")
    command.run(cmd)


sample = snakemake.wildcards.sample #"I15546-L1" #
project = snakemake.wildcards.project #'strawpigs' #
threads = snakemake.threads #"3" #
log = snakemake.log

ref_indeces = snakemake.input.ref_genomes #['resources/reference/pig.mmi', 'resources/reference/human.mmi'] #
#ref = ["pig", 'human'] #snakemake.wildcards.ref

filt_path = snakemake.params.filt_path #f"results/{project}/filtered" #
final_filtered_1 = snakemake.output.filt_fastq1 #f"{filt_path}{sample}_final.1.fastq.gz" #from 
final_filtered_2 = snakemake.output.filt_fastq2 #f"{filt_path}{sample}_final.2.fastq.gz" #from 

if not path.exists(f'{filt_path}/'):
    cmd = f'mkdir {filt_path}/'
    os.system(cmd)

to_delete_ls = []

first_ref = True
for ref_index in ref_indeces:
    ref_name = ref_index.split('/')[-1].split('.')[0]

    if first_ref:
        fastq_in = snakemake.input.query #[f"results/{project}/trimmed/fastp/{sample}.1.fastq.gz", f"results/{project}/trimmed/fastp/{sample}.2.fastq.gz"] #snakemake.input.query
        first_ref = False
    else:
        fastq_in = [filt_fastq_1, filt_fastq_2]

    alignment = f"{filt_path}/{sample}_{ref_name}.sam"
    sorted_alignment = f"{filt_path}/{sample}_{ref_name}.sorted.sam"

    filt_fastq_1 = f"{filt_path}/{sample}_{ref_name}.1.fastq.gz"
    filt_fastq_2 = f"{filt_path}/{sample}_{ref_name}.2.fastq.gz"

    ## files that are discarded afterwards, reads will otherwise be printed to stdout
    single_reads = f"{filt_path}/{sample}_{ref_name}.single.fastq.gz"
    supp_reads = f"{filt_path}/{sample}_{ref_name}.supp.fastq.gz"

    to_delete_ls.extend([filt_fastq_1, filt_fastq_2, single_reads, supp_reads, alignment, sorted_alignment])

    fastqs = ' '.join(fastq_in)
    cmd = f"minimap2 -t {threads} -o {alignment} -a {ref_index} {fastqs} 2>> {log}"
    os.system(cmd)

    cmd = f"samtools sort -n -o {sorted_alignment} {alignment} 2>> {log}"
    os.system(cmd)

    cmd = f'samtools fastq -f 4 -n -s {single_reads} -1 {filt_fastq_1} -2 {filt_fastq_2} -0 {supp_reads} {sorted_alignment} 2>> {log}'
    os.system(cmd)

cmd = f'cp {filt_fastq_1} {final_filtered_1}'
os.system(cmd)

cmd = f'cp {filt_fastq_2} {final_filtered_2}'
os.system(cmd)

to_delete = ' '.join(to_delete_ls)
cmd = f'rm {to_delete}'
os.system(cmd)
    
# minimap gegen verschiedene references nacheinander
# bei mehreren werden unmapped output reads (filtered fastq) als input für minimap genutzt
# letztendlich filtered final file -> als Hauptoutput, alle ziwschen filter files werden gelöscht
# samtools fastq -F 4 -n -s {singletons_file} -1 {read1_file} -2 {read2_file} {output_align} ## host reads -> don't save!!
# out_filt_fastq_1 = "results/{project}/out_filtered/{sample}_{ref}.1.fastq.gz"
# out_filt_fastq_2 = "results/{project}/out_filtered/{sample}_{ref}.2.fastq.gz"