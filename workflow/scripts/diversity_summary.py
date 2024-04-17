import pandas as pd
import json
import sys

sys.stderr = open(snakemake.log[0], "w")

jsons = snakemake.input.jsons
human_logs = snakemake.input.human_logs
if snakemake.params.other_host:
    hostname=snakemake.params.hostname
    host_logs = snakemake.input.host_logs
reports = snakemake.input.reports
tax_ids = snakemake.params.taxid_dict
outfile = snakemake.output.csv

# header for kraken report file
header_ls = ["prct", "total_reads", "lvl_reads", "lvl", "tax_id", "name"]
results_dict = {}

for report in reports:
    sample = report[report.rfind("/") + 1 : report.rfind("_kraken2_report")]
    sample_results_dict = {}

    for json_file in jsons:
        if json_file.find(sample) >= 0:
            with open(json_file, "r") as read_file:
                jdata = json.load(read_file)
                # fastp json shows reads instead of read pairs as in the kraken reports
                total_pre_host_filt = int((jdata["summary"]["after_filtering"]["total_reads"])/2)
                print(f"after fastp pairs: {total_pre_host_filt}")
            break
    sample_results_dict["#total_reads"] = total_pre_host_filt

    if snakemake.params.other_host:
        for host_log in host_logs:
            if host_log.find(sample) >= 0:
                with open(host_log, "r") as file:
                    for line in file.readlines():
                        if line.find("processed") >= 0:
                            non_host_reads = int(line.split()[2])
                            print(f"non pig pairs: {non_host_reads}")
                            no_reads = total_pre_host_filt - non_host_reads
                            print(f"pig pairs: {no_reads}")
                            prct_reads = (int(no_reads) / int(total_pre_host_filt)) * 100
                            break

                sample_results_dict[f"%{hostname}"] = "%.3f" %prct_reads
                sample_results_dict[f"#reads_{hostname}"] = no_reads
                break


    for human_log in human_logs:
        if human_log.find(sample) >= 0:
            with open(human_log, "r") as file:
                for line in file.readlines():
                    if line.find("processed") >= 0:
                        non_human_reads = int(line.split()[2])
                        print(f"non human pairs: {non_human_reads}")
                        if snakemake.params.other_host:
                            no_reads = non_host_reads - non_human_reads
                        else:
                            no_reads = total_pre_host_filt - non_human_reads
                        print(f"human pairs: {no_reads}")
                        prct_reads = (int(no_reads) / int(total_pre_host_filt)) * 100
                        break

            sample_results_dict[f"%human"] = "%.3f" %prct_reads
            sample_results_dict[f"#reads_human"] = no_reads
            break
        

    report_df = pd.read_table(report, names=header_ls)

    for ref in tax_ids:
        line = report_df.loc[report_df["tax_id"] == tax_ids[ref]]

        no_reads = line.iloc[0]["total_reads"]
        prct_reads = (int(no_reads) / int(total_pre_host_filt)) * 100

        sample_results_dict[f"%{ref}"] = "%.3f" %prct_reads
        sample_results_dict[f"#reads_{ref}"] = no_reads
    
    # add unclassified reads
    line = report_df.loc[report_df["tax_id"] == 0]
    if not line.empty:
        no_reads = line.iloc[0]["total_reads"]
        prct_reads = (int(no_reads) / int(total_pre_host_filt)) * 100

    else:
        no_reads = 0
        prct_reads = 0

    sample_results_dict["%unclassified"] = "%.3f" %prct_reads
    sample_results_dict["#reads_unclassified"] = no_reads


    results_dict[sample] = sample_results_dict


results_df = pd.DataFrame.from_dict(results_dict, orient="index")
results_df.index.name = "sample"

results_df=results_df.sort_index()
results_df.to_csv(outfile)
