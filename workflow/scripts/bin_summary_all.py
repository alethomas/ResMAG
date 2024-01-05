import pandas as pd
import sys
import os

sys.stderr = open(snakemake.log[0], "w")
csv_mags = snakemake.input.csv_mags
csv_bins = snakemake.input.csv_bins

summary_dict = {}
for bin_file in csv_bins:
    sample = os.path.basename(os.path.dirname(bin_file))

    summary_sample_dict = {}

    mag_file=[file for file in csv_mags if os.path.basename(os.path.dirname(file)) == sample][0]
    mag_df=pd.read_csv(mag_file)
    summary_sample_dict["# MAGs"] = len(mag_df)

    bin_df=pd.read_csv(bin_file)
    summary_sample_dict["# bins"] = len(bin_df)

    summary_sample_dict["least contigs"] = bin_df["contigs"].min()

    summary_sample_dict["highest N50"] = f'{(bin_df["contig_N50"].max()):,}'

    summary_dict[sample]=summary_sample_dict

summary_df = pd.DataFrame.from_dict(summary_dict, orient="index")
summary_df.index.name = "sample"
summary_df.sort_index(inplace=True)

summary_df.to_csv(snakemake.output[0])
