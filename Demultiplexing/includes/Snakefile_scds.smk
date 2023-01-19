#!/usr/bin/env python
import os
import pandas as pd
from glob import glob

##############################
############ SCDS ############
##############################
rule scds:
    input:
        matrix_dir = lambda wildcards: scrnaseq_libs_df["Matrix_Directories"][wildcards.pool],
    output: 
        doublets= "results/{pool}/scds/scds_doublets.txt",
        variables = "results/{pool}/scds/scds_variables.txt"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * scds_dict["scds_memory"],
        disk_mb=lambda wildcards, attempt: attempt * scds_dict["scds_memory"]
    threads: scds_dict["scds_threads"]
    params:
        script = "/opt/WG1-pipeline-QC_GoogleCloud/Demultiplexing/scripts/scds.R",
        out = workflow.default_remote_prefix + "/results/{pool}/scds/",
    log: "results/logs/scds.{pool}.log"
    shell:
        """
        echo {wildcards.pool} > {output.variables}
        echo {params.out} >> {output.variables}
        echo {input.matrix_dir} >> {output.variables}
        Rscript {params.script} {output.variables}

        """


rule scds_results_temp:
    input:
        "results/{pool}/scds/scds_doublets.txt"
    output:
        "results/{pool}/CombinedResults/scds_results.txt"
    resources:
        mem_mb=15000,
        disk_mb=15000
    threads: 1
    shell:
        """
        awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' {input}  > {output}
        """
