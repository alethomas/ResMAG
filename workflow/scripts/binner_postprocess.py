import os
import sys
import pandas as pd
from Bio import SeqIO


sys.stderr = open(snakemake.log[0], "w")


def fasta_to_tsv(directory, output_tsv):
    results_df = pd.DataFrame()
    for file_name in os.listdir(directory):
        if (
            file_name.endswith(".fasta")
            or file_name.endswith(".fa")
            or file_name.endswith(".fna")
        ):
            record_dict = SeqIO.to_dict(
                SeqIO.parse(f"{directory}/{file_name}", "fasta")
            )
            for record in record_dict:
                record_dict[record] = file_name

            if results_df.empty:
                results_df = pd.DataFrame.from_dict(record_dict, orient="index")
            else:
                df = pd.DataFrame.from_dict(record_dict, orient="index")
                results_df = pd.concat([results_df, df])

    results_df.to_csv(output_tsv, sep="\t", header=False)


# Get the input and output file paths from command-line arguments
directory = snakemake.input[0]
output_tsv = snakemake.output[0]

# Call the function with the provided file paths

fasta_to_tsv(directory, output_tsv)
