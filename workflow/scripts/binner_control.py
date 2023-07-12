import os
import sys
import pathlib

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

# Check if the correct number of command-line arguments is provided
if len(sys.argv) != 7:
    print("Usage: python script.py <Binner 1> <Binner 2> <Binner 3> <Binner 4> <Binner 5> <Output>")
    sys.exit(1)

# Get the input and output file paths from command-line arguments
files = [
    sys.argv[1],
    sys.argv[2],
    sys.argv[3],
    sys.argv[4],
    sys.argv[5]
]
output_csv = sys.argv[6]

# Call the function with the provided file paths
check_files(files)