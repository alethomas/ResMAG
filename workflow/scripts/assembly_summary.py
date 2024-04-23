import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")
log_files = snakemake.input

summary_dict = {}
for log_file in log_files:
    summary_sample_dict = {}

    with open(log_file, "r") as f_log:
        for line in f_log:
            if line.find("total") >= 0:
                info = line.split(" - ")[-1]
                seps = info.strip().split(", ")
                for sep in seps:
                    sep = sep.split()
                    if "contigs" in sep:
                        colname = "# contigs"
                        value = int(sep[0])
                        summary_sample_dict[colname] = f"{value:,}"
                    elif "total" in sep:
                        colname = "total bp"
                        value = int(sep[1])
                        summary_sample_dict[colname] = f"{value:,}"
                    else:
                        colname = "contigs {} bp".format(sep[0])
                        value = int(sep[1])
                        summary_sample_dict[colname] = f"{value:,}"

    sample = log_file[log_file.rfind("/") + 1 : log_file.rfind("_megahit")]
    summary_dict[sample] = summary_sample_dict

summary_df = pd.DataFrame.from_dict(summary_dict, orient="index")
summary_df.index.name = "sample"
summary_df.to_csv(snakemake.output.csv)
