import os
import sys

sys.stderr = open(snakemake.log[0], "w")

def check_files(files, binner_ls):
    path_ls = []
    results_binner = binner_ls.copy()

    for binner in binner_ls:
        for file_path in files:
            if file_path.find(binner) >= 0:
                # check if file exists and is not empty
                if os.path.isfile(file_path) and os.path.getsize(file_path) == 0:
                    results_binner.remove(binner)
                    break
                else:
                    path_ls.append(file_path)
                    break
                
    with open(output_csv, 'w') as f_out:
        f_out.write("{paths}\n{binners}\n".format(paths=(",".join(path_ls)),binners=(",".join(results_binner))))



# Get the input and output file paths from command-line arguments
files = snakemake.input
binner_ls = snakemake.params[0]
output_csv = snakemake.output[0]

# Call the function with the provided file paths
check_files(files, binner_ls)