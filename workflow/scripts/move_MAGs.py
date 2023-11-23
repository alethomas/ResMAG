import pandas as pd
import sys
import os

sys.stderr = open(snakemake.log[0], "w")


summary = snakemake.input.csv
bin_dir=snakemake.input.bin_dir
outdir= snakemake.output.outdir

mags_df=pd.read_csv(summary, index_col=0)
mags=mags_df.index.to_list()

if len(mags) == 0:
    cmd=f"mkdir -p {outdir}"
    os.system(cmd)
else:
    for mag in mags:
        cmd=f"mkdir -p {outdir} && mv {bin_dir}/{mag}.fa.gz {outdir}/"
        os.system(cmd)