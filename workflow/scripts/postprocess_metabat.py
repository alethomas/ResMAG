import os
import sys

sys.stderr = open(snakemake.log[0], "w")

def fasta_to_tsv(directory, output_tsv):
    with open(output_tsv, 'w') as f_out:
        f_out.write("Sequence\tFile\n")
        
        for file_name in os.listdir(directory):
            if file_name.endswith(".fasta") or file_name.endswith(".fa"):
                file_path = os.path.join(directory, file_name)
                
                with open(file_path, 'r') as f_in:
                    for line in f_in:
                        if line.startswith(">"):
                            sequence_name = line[1:].strip()
                            f_out.write(f"{sequence_name}\t{file_name}\n")

# Get the input and output file paths from command-line arguments
directory = snakemake.input[0]
output_tsv = snakemake.output[0]

# Call the function with the provided file paths

fasta_to_tsv(directory, output_tsv)