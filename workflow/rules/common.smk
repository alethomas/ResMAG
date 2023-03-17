configfile: "config/config.yaml"

def get_samples():
    return list(pep.sample_table["sample_name"].values)

def get_fastqs(wildcards):
    return pep.sample_table.loc[wildcards.sample][["fq1", "fq2"]]

def get_adapters(wildcards):
    return config["nextera_transposase"]#"adapter_seqs"]

def get_samples_for_project(project):
    # select samples for given project
    df = pep.sample_table
    df = df[df["project"] == project]

    samples_of_run = list(df["sample_name"].values)
    return samples_of_run

def expand_samples_for_project(paths, **kwargs):
    def inner(wildcards):
        return expand(
            paths,
            sample=get_samples_for_project(wildcards.project),
            **kwargs,
        )

    return inner
