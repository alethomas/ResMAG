import pandas as pd
import json
import sys

sys.stderr = open(snakemake.log[0], "w")


csv = snakemake.input.csv
jsons = snakemake.input.jsons
host=snakemake.params.host_name
outfile = snakemake.output.csv

results_dict = {}

kraken_df= pd.read_csv(csv, index_col=0)

for sample in kraken_df.index.to_list():
    sample_results_dict = {}
    total_reads=int((kraken_df.loc[sample]["total reads (#)"]).replace(",", ""))

    for json_file in jsons:
        if json_file.find(sample) >= 0:
            with open(json_file, "r") as read_file:
                jdata = json.load(read_file)
                # fastp json shows reads instead of read pairs as in the kraken reports
                total_pre_host_filt = int((jdata["summary"]["after_filtering"]["total_reads"])/2)

            host_reads = int(total_pre_host_filt - total_reads)
            host_perc = (host_reads/total_pre_host_filt)*100
            non_host_perc = (total_reads/total_pre_host_filt)*100

            sample_results_dict["total reads before (#)"] = f'{total_pre_host_filt:,}'
            sample_results_dict[f"{host} (%)"] = "%.2f" % host_perc
            sample_results_dict[f"{host} (#)"] = f'{host_reads:,}'
            sample_results_dict[f"non-{host} (%)"] = "%.2f" % non_host_perc
            sample_results_dict[f"non-{host} (#)"] = f'{total_reads:,}'

            jsons.remove(json_file)

    results_dict[sample] = sample_results_dict


results_df = pd.DataFrame.from_dict(results_dict, orient="index")
results_df.index.name = "sample"

results_df.to_csv(outfile)
