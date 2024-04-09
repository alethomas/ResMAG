import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")


reports = snakemake.input.reports
tax_ids = snakemake.params.taxid_dict
outfile = snakemake.output.csv


header_ls = ["prct", "total_reads", "lvl_reads", "lvl", "tax_id", "name"]
results_dict = {}

for report in reports:
    report_df = pd.read_table(report, names=header_ls)
    sample_results_dict = {}

    for ref in tax_ids:
        line = report_df.loc[report_df["tax_id"] == tax_ids[ref]]
        if line.iloc[0]["lvl"] == "S":
            no_reads = line.iloc[0]["lvl_reads"]
        else:
            no_reads = line.iloc[0]["total_reads"]
        prct_reads = line.iloc[0]["prct"]

        sample_results_dict[f"{ref} (%)"] = prct_reads
        sample_results_dict[f"{ref} (#)"] = f'{no_reads:,}'
    
    line = report_df.loc[report_df["tax_id"] == 0]
    if not line.empty:
        prct_reads = line.iloc[0]["prct"]
        no_reads = line.iloc[0]["total_reads"]

    else:
        prct_reads = 0
        no_reads = 0

    sample_results_dict["unclassified (%)"] = prct_reads
    sample_results_dict["unclassified (#)"] = f'{no_reads:,}'

    line = report_df.loc[report_df["tax_id"] == 1]
    total_reads = no_reads + line.iloc[0]["total_reads"]

    sample_results_dict["total reads (#)"] = f'{total_reads:,}'

    sample = report[report.rfind("/") + 1 : report.rfind("_report")]
    results_dict[sample] = sample_results_dict


results_df = pd.DataFrame.from_dict(results_dict, orient="index")
results_df.index.name = "sample"

# move total reads column to the front of the table
first_column = results_df.pop("total reads (#)") 
results_df.insert(0, "total reads (#)", first_column)

results_df.to_csv(outfile)
