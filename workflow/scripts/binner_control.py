import os
import sys
import pathlib

sys.stderr = open(snakemake.log[0], "w")

def check_files(files):
    paths = []
    binner= []
    for file_path in files:
        if os.path.isfile(file_path) and os.path.getsize(file_path) == 0:
            continue
        else:
            paths.append(f"{file_path}")
            files = ",".join(paths)
            binner.append(file_path.split("/")[2]) 
            binners = ",".join(binner)
    with open(output_csv, 'w') as f_out:
        f_out.write(f"{files}\n{binners}\n")


# Get the input and output file paths from command-line arguments
files = snakemake.input
output_csv = snakemake.output[0]

# Call the function with the provided file paths
check_files(files)