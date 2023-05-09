import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")
f_log = open(snakemake.input[0], "r")

summary_dict = {}
summary_sample_dict = {}

for line in f_log:
    if line.find("total") >= 0:
        info = line.split(" - ")[-1]
        seps = info.strip().split(", ")
        for sep in seps:
            sep = sep.split()
            if "contigs" in sep:
                summary_sample_dict[sep[1]] = sep[0]
            else:
                colname = sep[0] + "_bp"
                summary_sample_dict[colname] = sep[1]
f_log.close()
summary_dict[snakemake.wildcards.sample] = summary_sample_dict

summary_df = pd.DataFrame.from_dict(summary_dict, orient="index")
summary_df.index.name = "sample"
summary_df.to_csv(snakemake.output[0])
