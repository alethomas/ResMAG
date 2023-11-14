import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")

## input files
in_dastool= snakemake.input.tool #"ResMAG/results/autobrewer/das_tool/ABS_24_07/ABS_24_07_DASTool_summary.tsv"
in_checkm = snakemake.input.checkm #"ResMAG/results/autobrewer/qc/checkm2/ABS_24_07/quality_report.tsv"
in_gtdb=snakemake.input.gtdb #"ResMAG/results/autobrewer/classification/ABS_24_07/ABS_24_07.bac120.summary.tsv"

## output files
### csv
csv_path_mags = snakemake.output.csv_mags #"ResMAG/results/autobrewer/report/ABS_24_07/mags_summary.csv" #
csv_path_bins = snakemake.output.csv_bins #"ResMAG/results/autobrewer/report/ABS_24_07/bin_summary.csv" 
csv_path_tax = snakemake.output.csv_tax #"ResMAG/results/autobrewer/report/ABS_24_07/bin_taxonomy.csv"
csv_path_checkm = snakemake.output.csv_checkm #"ResMAG/results/autobrewer/report/ABS_24_07/checkm_summary.csv"
csv_path_dastool = snakemake.output.csv_dastool #"ResMAG/results/autobrewer/report/ABS_24_07/DASTool_summary.csv"
### html
html_path_mags = snakemake.output.html_mags #"ResMAG/results/autobrewer/report/ABS_24_07/mags_summary.html"
html_path_bins = snakemake.output.html_bins #"ResMAG/results/autobrewer/report/ABS_24_07/bin_summary.html" 
html_path_tax = snakemake.output.html_tax #"ResMAG/results/autobrewer/report/ABS_24_07/bin_taxonomy.html" 
html_path_checkm = snakemake.output.html_checkm #"ResMAG/results/autobrewer/report/ABS_24_07/checkm_summary.html"
html_path_dastool = snakemake.output.html_dastool #"ResMAG/results/autobrewer/report/ABS_24_07/DASTool_summary.html"

## params
max_cont= snakemake.params.max_cont
min_comp= snakemake.params.min_comp


def save_html_csv_table(csv_path, html_path, summary_df, title):
    summary_df.columns.name='bin_name'
    summary_df.index.name = None
    summary_df.to_csv(csv_path)

    tabtitle = title

    tablestyle = [
        dict(selector="tr:hover", props=[("background-color", "#9ec83e")]),
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


tool_df = pd.read_table(in_dastool)
tool_df.drop(["bin_set"], axis=1, inplace=True)
tool_df.rename({"bin":"bin_name"},axis=1,inplace=True)
tool_df.set_index("bin_name",inplace=True)
col_order=["bin_score","contigs","size","N50","SCG_set","unique_SCGs","SCG_completeness","redundant_SCGs","SCG_redundancy"]
tool_df=tool_df[col_order]

del col_order[:2]
tool_red_df= tool_df.drop(col_order, axis=1)

tab_title="DAS Tool bin results"
save_html_csv_table(csv_path_dastool, html_path_dastool, tool_red_df, tab_title)


checkm_df = pd.read_table(in_checkm)
rm_list=["Translation_Table_Used", "Additional_Notes"]
checkm_df.drop(rm_list, axis=1, inplace=True)
checkm_df.columns=checkm_df.columns.str.lower()
checkm_df.rename({"name":"bin_name", "gc_content":"GC_content","contig_n50":"contig_N50"},axis=1,inplace=True)
checkm_df.set_index("bin_name",inplace=True)
col_order=["completeness","contamination","genome_size","GC_content","contig_N50",
           "total_coding_sequences","coding_density","average_gene_length","completeness_model_used"]
checkm_df=checkm_df[col_order]

del col_order[:5]
checkm_red_df=checkm_df.drop(col_order, axis=1)

tab_title="CheckM2 bin results"
save_html_csv_table(csv_path_checkm, html_path_checkm, checkm_red_df, tab_title)


gtdb_df=pd.read_table(in_gtdb)
gtdb_df.rename({"user_genome":"bin_name"},axis=1,inplace=True)
gtdb_df.set_index("bin_name",inplace=True)
cols=['classification','fastani_reference','fastani_reference_radius','fastani_ani','fastani_af','classification_method']
gtdb_red_df=gtdb_df[cols]

tab_title="Bins taxonomy classification"
save_html_csv_table(csv_path_tax, html_path_tax, gtdb_red_df, tab_title)


bins_df = pd.concat([tool_red_df,checkm_red_df,gtdb_red_df[['classification','fastani_reference','fastani_ani']]], axis=1)
col_order = ["completeness","contamination","bin_score","contigs","genome_size","contig_N50","GC_content",'classification','fastani_reference']
bins_df=bins_df[col_order]

tab_title="Bins summary"
save_html_csv_table(csv_path_bins, html_path_bins, bins_df, tab_title)


mags_df = bins_df.loc[(bins_df["completeness"] >= min_comp) & (bins_df["contamination"]<= max_cont),:]
tab_title="MAGs summary"
save_html_csv_table(csv_path_mags, html_path_mags, mags_df, tab_title)


"""
if we want to sort by mean of deviation from min completeness & max contamination

mean_dict={}
for binid in bins_df.index:
    cont=bins_df.at[binid,"contamination"]
    cont_per=round(((max_contamination - cont)/max_contamination), 4)
    comp=bins_df.at[binid,"completeness"]
    comp_per=round((comp/100),4)
    mean_dict[binid]=round(((cont_per + comp_per)/2),4)

to_sort=pd.DataFrame.from_dict(mean_dict,orient="index",columns=['mean'])

# index list sorted by mean used to reindex original df
new_ind=to_sort.sort_values(["mean"],ascending = [False]).index.to_list()
sorted_all=bins_df.reindex(new_ind)

"""