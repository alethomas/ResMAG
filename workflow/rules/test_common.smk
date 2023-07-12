rule test_rule:
    input:
        "results/dummyfile.csv",
    output:
        "results/test.txt",
    params:
        paths=get_paths_binner("results/dummyfile.csv")[0],
        binner=get_paths_binner("results/dummyfile.csv")[1],
    log:
        "results/test.log"
    shell:
        "echo {params.paths} {params.binner} > {output} "
