import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")


reports = snakemake.input.reports
tax_ids = snakemake.params.taxid_dict
outfile = snakemake.output.csv
html_path = snakemake.output.html


header_ls = ["prct", "total_reads", "lvl_reads", "lvl", "tax_id", "name"]
results_dict = {}

for report in reports:
    report_df = pd.read_table(report, names=header_ls)
    sample_results_dict = {}

    line = report_df.loc[report_df["tax_id"] == 0]
    if not line.empty:
        prct_reads = line.iloc[0]["prct"]
        no_reads = line.iloc[0]["total_reads"]

    else:
        prct_reads = 0
        no_reads = 0

    sample_results_dict["% unclassified"] = prct_reads
    sample_results_dict["# reads unclassified"] = no_reads

    for ref in tax_ids:
        line = report_df.loc[report_df["tax_id"] == tax_ids[ref]]
        if line.iloc[0]["lvl"] == "S":
            no_reads = line.iloc[0]["lvl_reads"]
        else:
            no_reads = line.iloc[0]["total_reads"]
        prct_reads = line.iloc[0]["prct"]

        sample_results_dict[f"% {ref}"] = prct_reads
        sample_results_dict[f"# reads {ref}"] = no_reads
    sample = report[report.rfind("/") + 1 : report.rfind("_report")]
    results_dict[sample] = sample_results_dict


results_df = pd.DataFrame.from_dict(results_dict, orient="index")

results_df.to_csv(outfile)


tabtitle = "Kraken2 Summary"

tablestyle = [
    dict(selector="tr:hover", props=[("background", "#9ec83e")]),
    dict(
        selector="td",
        props=[
            ("background", "white"),
            ("border-top", "2px solid #c0e078"),
            ("padding", "8px"),
            ("font-family", "Arial"),
            ("font-size", "15px"),
        ],
    ),
    dict(
        selector="th",
        props=[
            ("color", "black"),
            ("padding", "8px 13px"),
            ("background", "#c0e078"),
            ("border", "2px solid white"),
            ("font-size", "15px"),
            ("font-family", "Arial"),
        ],
    ),
    dict(
        selector="caption",
        props=[
            ("padding", "8px"),
            ("font-size", "18px"),
            ("background", "#c0e078"),
            ("border", "2px solid white"),
            ("caption-side", "top"),
            ("font-family", "Arial"),
        ],
    ),
]


dfall = (
    results_df.style.set_table_attributes('style="border-collapse:collapse"')
    .set_properties(**{"text-align": "center"})
    .format(precision=2)
    .set_caption(tabtitle)
    .set_table_styles(tablestyle)
)

with open(html_path, "w") as outhtml:
    dfall.to_html(buf=outhtml)
