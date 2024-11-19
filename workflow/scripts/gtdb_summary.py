import json
import os
from pathlib import Path

sys.stderr = open(snakemake.log[0], "w")

folder = snakemake.params.outdir
json_infile = snakemake.input.json
output_file = snakemake.output.summary

# list of summary outfiles produced by GTDBTk
summary_files = []

with open(json_infile, "r") as f:
    data = json.load(f)

# Find the output for the step with name 'classify'
classify_step = next(
    (step for step in data["steps"] if step["name"] == "classify"), None
)
if classify_step:
    output_files = classify_step["output_files"]
    for key in output_files.keys():
        for file in output_files[key]:
            if file.endswith(f".{key}.summary.tsv"):
                summary_files.append(file)
                filename = Path(file).name
                new_file = f"{folder}/{filename}"
                cmd = f"cp {file} {new_file}"
                os.system(cmd)

with open(output_file, "w") as outfile:
    for i, file in enumerate(summary_files):
        with open(file, "r") as infile:
            lines = infile.readlines()
            if i > 0:  # Skip header for all but the first file
                lines = lines[1:]
            outfile.writelines(lines)
