import json
import sys

sys.stderr = open(snakemake.log[0], "w")

def json_to_tsv(json_file, tsv_file):
    with open(json_file, 'r') as f:
        data = json.load(f)

    with open(tsv_file, 'w') as f:
        for key, value in data.items():
            f.write(f"{value}\t{key}\n")


# Get the input and output file paths from command-line arguments
json_file = snakemake.input[0]
tsv_file = snakemake.output[0]

# Call the function with the provided file paths
json_to_tsv(json_file, tsv_file)