import json
import sys

def json_to_tsv(json_file, tsv_file):
    with open(json_file, 'r') as f:
        data = json.load(f)

    with open(tsv_file, 'w') as f:
        for key, value in data.items():
            f.write(f"{value}\t{key}\n")

# Check if the correct number of command-line arguments is provided
if len(sys.argv) != 3:
    print("Usage: python script.py <input_json_file> <output_tsv_file>")
    sys.exit(1)

# Get the input and output file paths from command-line arguments
json_file = sys.argv[1]
tsv_file = sys.argv[2]

# Call the function with the provided file paths
json_to_tsv(json_file, tsv_file)