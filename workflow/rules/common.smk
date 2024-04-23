import os
import pandas as pd


configfile: "config/config.yaml"


pathogen_file = pd.read_csv(config["pathogens"])


def get_root():
    return os.getcwd()


def get_data_path():
    return config["data-handling"]["data"]


def get_resource_path():
    return config["data-handling"]["resources"]


def get_project():
    return config["project-name"]


def get_samples():
    return list(pep.sample_table["sample_name"].values)


def get_pathogens():
    return list(pathogen_file["tax_abbr"].values)


def expand_samples_for_project(paths, **kwargs):
    def inner(wildcards):
        return expand(
            paths,
            sample=get_samples(),
            **kwargs,
        )

    return inner


def get_fastqs(wildcards):
    return (
        pep.sample_sheet.loc[wildcards.sample]["fq1"],
        pep.sample_sheet.loc[wildcards.sample]["fq2"],
    )


def get_local_fastqs(wildcards):
    path = get_data_path()
    return (
        "{data}{{project}}/{{sample}}_R1.fastq.gz".format(data=path),
        "{data}{{project}}/{{sample}}_R2.fastq.gz".format(data=path),
    )


def get_adapters(wildcards):
    return config["adapter_seqs"]


def get_trimmed_fastqs(wildcards):
    return [
        "results/{project}/trimmed/fastp/{sample}.1.fastq.gz",
        "results/{project}/trimmed/fastp/{sample}.2.fastq.gz",
    ]


def get_prefiltered_fastqs(wildcards):
    if config["host_filtering"]["do_host_filtering"]:
        return [
            "results/{project}/host_filtering/non_host/{sample}_R1.fastq.gz",
            "results/{project}/host_filtering/non_host/{sample}_R2.fastq.gz",
        ]
    else:
        return get_trimmed_fastqs(wildcards)


def get_host_map_statistics(wildcards):
    if config["host_filtering"]["do_host_filtering"]:
        logs = (
            expand_samples_for_project(
                "logs/{project}/pathogen_filtering/filter_{sample}_to_{pathogen}.log",
            ),
        )
        return logs
    else:
        return []


# this function only works for locally available genomes right now, since download link is link not path
def get_pathogen_ref(wildcards):
    if config["pathogen-recruitment"]["use-local"]:
        return "{}ref_genome/{}".format(
            get_resource_path(),
            pathogen_file.loc[
                pathogen_file["tax_abbr"] == wildcards.pathogen, "local_path"
            ]
            .values[0]
            .split("/")[-1],
        )
    else:
        return pathogen_file.loc[
            pathogen_file["tax_abbr"] == wildcards.pathogen, "download_link"
        ].values[0]


## TODO: rename function
def get_pathogen_local_folder(wildcards):
    path = pathogen_file.loc[
        pathogen_file["tax_abbr"] == wildcards.pathogen, "local_path"
    ].values[0]
    folder = Path(path).parent
    filename = Path(path).name
    return [folder, filename]


def get_checkm2_db():
    file = "{}{}".format(get_resource_path(), config["checkm2"])
    return file


def get_filtered_fastqs(wildcards):
    return [
        "results/{project}/filtered/fastqs/{sample}_to_{pathogen}_R1.fastq",
        "results/{project}/filtered/fastqs/{sample}_to_{pathogen}_R2.fastq",
    ]


def get_filtered_gz_fastqs(wildcards):
    return [
        "results/{project}/filtered/fastqs/{sample}_to_{pathogen}_R1.fastq.gz",
        "results/{project}/filtered/fastqs/{sample}_to_{pathogen}_R2.fastq.gz",
    ]


def get_assembly(wildcards):
    return "results/{project}/megahit/{sample}_to_{pathogen}/final.contigs.fa"


## binning parameters
def get_contig_length_threshold():
    return config["binning"]["min_contig_length"]


def get_contig_length_filter():
    filt = int(config["binning"]["min_contig_length"]) - 1
    return filt


# def get_kmersize():
#     return config["binning"]["kmer_length"]


# def get_binners():
#     return config["das_tool"]["binner-list"]


# def get_rosella_install():
#     folder = get_resource_path()
#     script = f"{folder}rosella/install.sh"
#     return script


# def get_rosella_git():
#     return config["rosella"]["gitURL"]


# def get_all_contig2bin_files(wildcards):
#     binners = get_binners()
#     file_list = [
#         f"results/{{project}}/output/contig2bins/{{sample}}/{binner}_contig2bin.tsv"
#         for binner in binners
#     ]
#     return file_list


# ## reads in binner control file and returns list with paths to contig2bin files
# ## and a list with name of the binners that produced results
# def get_paths_binner(wildcards):
#     file = f"results/{wildcards.project}/das_tool/{wildcards.sample}_binner_control.csv"
#     lines = open(file).readlines()
#     paths = str(lines[0].rstrip("\n"))
#     binner = str(lines[1].rstrip("\n"))
#     return paths, binner


# def bins_for_sample(wildcards):
#     if len(get_paths_binner[0]) > 0:
#         return True
#     else:
#         return False
# def get_DAS_Tool_threads():
#     return config["das_tool"]["threads"]
