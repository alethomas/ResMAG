import os


configfile: "config/config.yaml"


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
        pep.sample_table.loc[wildcards.sample]["fq1"],
        pep.sample_table.loc[wildcards.sample]["fq2"],
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


def get_kraken_db():
    return config["kraken"]["kraken-db"]


def get_local_krakenDB():
    path = "{}kraken2_plusPF/".format(get_resource_path())
    return path


def get_kraken_ref():
    return config["kraken"]["ref-analysis"]

def get_kraken_db_path():
    download = config["kraken"]["download-path"]
    db_path = config["kraken"]["kraken-db"]
    return db_path + download.split("/")[-1]

def get_taxID_dict():
    return config["kraken"]["taxIDs-ref"]


def get_taxID(wildcards):
    return config["kraken"]["taxIDs-ref"][wildcards.kraken_ref]


def get_filtered_fastqs(wildcards):
    return [
        "results/{project}/filtered/non_human/{sample}_all_1.fastq",
        "results/{project}/filtered/non_human/{sample}_all_2.fastq",
    ]


def get_filtered_gz_fastqs(wildcards):
    return [
        "results/{project}/filtered/non_human/{sample}_all_1.fastq.gz",
        "results/{project}/filtered/non_human/{sample}_all_2.fastq.gz",
    ]


## binning parameters
def get_threshold():
    return config["binning"]["min_contig_length"]


def get_kmersize():
    return config["binning"]["kmer_length"]


def get_binners():
    return config["das_tool"]["binner-list"]


def get_rosella_install():
    folder = get_resource_path()
    script = f"{folder}rosella/install.sh"
    return script


def get_rosella_git():
    return config["rosella"]["gitURL"]


def get_all_contig2bin_files(wildcards):
    binners = get_binners()
    file_list = [
        f"results/{{project}}/contig2bins/{{sample}}/{binner}_contig2bin.tsv"
        for binner in binners
    ]
    return file_list


## reads in binner control file and returns list with paths to contig2bin files
## and a list with name of the binners that produced results
def get_paths_binner(wildcards):
    file = f"results/{wildcards.project}/das_tool/{wildcards.sample}/binner_control.csv"
    lines = open(file).readlines()
    paths = str(lines[0].rstrip("\n"))
    binner = str(lines[1].rstrip("\n"))
    return paths, binner


def get_DAS_Tool_threads():
    return config["das_tool"]["threads"]
