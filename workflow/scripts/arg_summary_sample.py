import pandas as pd
import json
import sys

sys.stderr = open(snakemake.log[0], "w")

# one of reads, assembly, mags
in_case = snakemake.params.case

# minimal average coverage of reference sequence
min_coverage = 15

def do_card_main_postprocessing(infile,json_infile, min_coverage):
    to_keep=["Best_Identities","Percentage Length of Reference Sequence","Cut_Off","Best_Hit_ARO","ARO",
             "Pass_Bitscore","Best_Hit_Bitscore","AMR Gene Family","Drug Class","Resistance Mechanism","Antibiotic"]

    df=pd.read_table(infile, index_col="ORF_ID")
    df_red=df[to_keep]

    with open(json_infile, 'r') as f:
        data = json.load(f)

    eval_dict={}
    for orfid in df.index.to_list():
        ident=df.at[orfid,"ID"]
        evalue=data[orfid][ident]["evalue"]
        eval_dict[orfid]=evalue

    eval_df=pd.DataFrame.from_dict(eval_dict,orient='index',columns=["evalue"])

    df_red=pd.concat([eval_df,df_red],axis=1)
    df_red=df_red.sort_values(["Cut_Off","Best_Identities","Percentage Length of Reference Sequence"], ascending=[True,False,False])
    df_red = df_red[df_red['Percentage Length of Reference Sequence'] > min_coverage]
    df_red.index.name = "ORF_ID"

    return df_red


if in_case == "reads":
    infile=snakemake.input.txt
    outfile=snakemake.output.csv

    to_keep=["ARO Term","ARO Accession","All Mapped Reads","Average Percent Coverage",
            "Average Length Coverage (bp)","Reference Length","Average MAPQ (Completely Mapped Reads)",
            "Resistomes & Variants: Observed Pathogen(s)","AMR Gene Family","Drug Class","Resistance Mechanism"]
    
    df=pd.read_table(infile)
    df_red=df[to_keep]
    df_red=df_red.set_index("ARO Accession")
    df_red=df_red.sort_values(["Average Percent Coverage","All Mapped Reads"],ascending=False)
    df_red = df_red[df_red['Average Percent Coverage'] > min_coverage]

    df_red.to_csv(outfile)


elif in_case == "assembly":
    infile=snakemake.input.txt
    json_f = snakemake.input.json
    outfile=snakemake.output.csv
    
    asbl_df = do_card_main_postprocessing(infile,json_f, min_coverage)
    asbl_df.to_csv(outfile)


'''
if in_case == "mags":
    infile=f"/local/work/josefa/ResMAG/results/clinic_june_240903/output/ARGs/mags/{sample}/{sample}.txt"
    json_f= f"/local/work/josefa/ResMAG/results/clinic_june_240903/output/ARGs/mags/{sample}/{sample}.json"
    outfile=f"/local/work/josefa/ResMAG/results/clinic_june_240903/output/ARGs/mags/{sample}/{sample}_mags_ARGs.csv"



    
    for sample in samples:
    sample_path=f"{inpath}{sample}/"
    outfile=f"{sample_path}/ARGs_mags_{sample}.csv"
    list_outfile=f"{sample_path}/ARGs_mags_list_summary_{sample}.csv"

    mags=[f for f in os.listdir(sample_path) if f.endswith(".txt")]
    first_MAG=True
    
    for mag in mags:
        infile=f"{sample_path}{mag}"
        json_f= infile.replace(".txt",".json")

        df=pd.read_table(infile)
        if df.empty:
            continue
        df.set_index("ORF_ID",inplace=True)
        df_red=df[cols]

        with open(json_f, 'r') as f:
            data = json.load(f)
        eval_dict={}
        for orfid in df.index.to_list():
            ident=df.at[orfid,"ID"]
            evalue=data[orfid][ident]["evalue"]
            eval_dict[orfid]=evalue
        eval_df=pd.DataFrame.from_dict(eval_dict,orient='index',columns=["evalue"])
        df_red=pd.concat([eval_df,df_red],axis=1)
        if first_MAG:
            df_all=df_red
            first_MAG=False
        else:
            df_all=pd.concat([df_all,df_red])
    if len(mags)<1:
        df_all=df_red[0:0]
    df_all=df_all.sort_values(["Cut_Off","Best_Hit_Bitscore"], ascending=[True,False])
    df_all.to_csv(outfile)
'''

