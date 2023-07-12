import json
import os


configfile: "config/config.yaml"


def get_samples():
    return list(pep.sample_table["sample_name"].values)


def get_fastqs(wildcards):
    return (
        pep.sample_table.loc[wildcards.sample]["fq1"],
        pep.sample_table.loc[wildcards.sample]["fq2"],
    )


def get_project():
    return config["project-name"]


def get_adapters(wildcards):
    return config["adapter_seqs"]


def get_executable_dir(executable):
    import os

    path = os.system(f"which run_metabinner.sh")
    print(path)
    return os.path.dirname(path)


def get_threshold():
    return config["binning"]["min_contig_length"]


def get_kmersize():
    return config["binning"]["kmer_length"]


"""def get_samples_for_project(project):
    # select samples for given project
    df = pep.sample_table
    df = df[df["project"] == project]

    samples_of_run = list(df["sample_name"].values)
    return samples_of_run"""


def expand_samples_for_project(paths, **kwargs):
    def inner(wildcards):
        return expand(
            paths,
            sample=get_samples(),  # _for_project(wildcards.project),
            **kwargs,
        )

    return inner


def get_trimmed_fastqs(wildcards):
    return [
        "results/{project}/trimmed/fastp/{sample}.1.fastq.gz",
        "results/{project}/trimmed/fastp/{sample}.2.fastq.gz",
    ]


def get_reference_species():
    return config["reference-genomes"]["reference-species"]


def get_reference_dir():
    return config["reference-genomes"]["reference-dir"]


def get_reference_file(wildcards):
    ext = config["reference-genomes"]["reference-file-ext"]
    ref_path = "{}{}{}".format(get_reference_dir(), wildcards.ref, ext)
    return ref_path


def get_reference_index(wildcards):
    ind_path = "{}{{ref}}.mmi".format(get_reference_dir())
    return expand(ind_path, ref=get_reference_species())
    """ind_list = []
    for ref in get_reference_species():
        ind_list.append("{}{}.mmi".format(get_reference_dir(),ref))
    return ind_list"""


def get_bacterial_reads(wildcards):
    return [
        "results/{project}/filtered/bacteria/{sample}_bacteria_1.fastq.gz",
        "results/{project}/filtered/bacteria/{sample}_bacteria_2.fastq.gz",
    ]


def get_kraken_db():
    return config["kraken"]["kraken-db"]


def get_kraken_ref():
    return config["kraken"]["ref-analysis"]


def get_taxID_dict():
    return config["kraken"]["taxIDs-ref"]


def get_taxID(wildcards):
    return config["kraken"]["taxIDs-ref"][wildcards.kraken_ref]

def get_rosella_install():
    folder = get_rosella_folder()
    script = f"{folder}/install.sh"
    return script


def get_rosella_git():
    return config["rosella"]["gitURL"]


def get_rosella_folder():
    return config["rosella"]["rosella_dir"]


def get_root():
    return os.getcwd()

def get_paths_binner(file):
    lines = open(file).readlines()
    paths = str(lines[0].rstrip('\n'))
    binner= str(lines[1].rstrip('\n'))
    return paths, binner
