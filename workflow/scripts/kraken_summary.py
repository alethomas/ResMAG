import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")

fastp_json = snakemake.input.json #f'results/{project}/trimmed/fastp/{sample}.fastp.json'
report = snakemake.input.report #f'results/{project}/filtered/kraken/{sample}_report.tsv'
tax_ids = snakemake.params.taxid_dict #{'bacteria': 2, 'fungi': 4751, 'pig': 9823, 'human': 9606, 'plant': 33090}
outfile = snakemake.output[0] #f'results/{project}/filtered/kraken/{sample}_summary.csv'


header_ls = ['prct', 'total_reads', 'lvl_reads', 'lvl', 'tax_id', 'name']
report_df = pd.read_table(report, names=header_ls)
results_dict = {}

line = report_df.loc[report_df['tax_id']== 0]
if not line.empty:
    no_reads = line.iloc[0]['total_reads']
    prct_reads = line.iloc[0]['prct']
    results_dict['unclassified'] = [no_reads, prct_reads]
else:
    results_dict['unclassified'] = [0, 0]

for ref in tax_ids:
    line = report_df.loc[report_df['tax_id']== tax_ids[ref]]
    if line.iloc[0]['lvl'] == 'S':
        no_reads = line.iloc[0]['lvl_reads']
    else:
        no_reads = line.iloc[0]['total_reads']
    prct_reads = line.iloc[0]['prct']
    results_dict[ref] = [no_reads, prct_reads]

header = ['#read_pairs_per_reference', '%reference_reads_of_total']
results_df = pd.DataFrame.from_dict(results_dict, orient='index', columns=header)
results_df.index.name = 'reference'
results_df.to_csv(outfile)
