import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")
log_files = snakemake.input
html_path = snakemake.output.html

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
                        summary_sample_dict[colname] = sep[0]
                    elif "total" in sep:
                        colname = "total bp"
                        summary_sample_dict[colname] = sep[1]
                    else:
                        colname = "contigs {} bp".format(sep[0])
                        summary_sample_dict[colname] = sep[1]

    sample = log_file[log_file.rfind("/") + 1 : log_file.rfind("_megahit")]
    summary_dict[sample] = summary_sample_dict

summary_df = pd.DataFrame.from_dict(summary_dict, orient="index")
summary_df.to_csv(snakemake.output.csv)

tabtitle = "Assembly Summary"

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
    summary_df.style.set_table_attributes('style="border-collapse:collapse"')
    .set_properties(**{"text-align": "center"})
    .format(precision=2)
    .set_caption(tabtitle)
    .set_table_styles(tablestyle)
)

with open(html_path, "w") as outhtml:
    dfall.to_html(buf=outhtml)
