#!/usr/bin/env python
import os
import pandas as pd
from glob import glob

##################################
############ SCRUBLET ############
##################################
rule make_scrublet_selection_df:
    input:
        input_dict["samplesheet_filepath"]
    output:
        "results/manual_selections/scrublet/scrublet_percentile_manual_selection.tsv"
    resources:
        mem_mb = 15000,
        disk_mb = 15000
    threads: 1
    log: "results/logs/make_scrublet_selection_df.log"
    shell:
        """
        awk 'BEGIN{{OFS=FS="\t"}}{{print $1 "\t"}}' {input} | sed "1s/.*/Pool\tscrublet_Percentile/" > {output} 2> {log}
        """

if bucket.blob("results/manual_selections/scrublet/scrublet_percentile_manual_selection.tsv").exists():

    blob_scrub = bucket.get_blob("results/manual_selections/scrublet/scrublet_percentile_manual_selection.tsv")
    downloaded_blob_scrub = blob_scrub.download_as_string()
    scrublet_selection = pd.read_csv(BytesIO(downloaded_blob_scrub), sep = "\t")

    if scrublet_selection["scrublet_Percentile"].count() != len(scrublet_selection):
        ready = False
    elif scrublet_selection["scrublet_Percentile"].count() == len(scrublet_selection):
        ready = True
        step = "ready"
    else:
        sys.exit()

    if (scrublet_manual_dict["run_scrublet_manual"] == False) or (scrublet_manual_dict["run_scrublet_manual"] == True and scrublet_selection["scrublet_Percentile"].count() == len(scrublet_selection)):
        log = workflow.default_remote_prefix + "results/{pool}/scrublet_{pctl}/default_run_variables.txt"
        sim_dbl = scrublet_extra_dict["sim_dbl"]
        min_counts = scrublet_extra_dict["min_counts"]
        min_cells = scrublet_extra_dict["min_cells"]
        n_prin_comps = scrublet_extra_dict["n_prin_comps"]
        step = "default"
        scrublet_doublet_threshold = None
    elif (scrublet_manual_dict["run_scrublet_manual"] == True): ### This deals with the possibility that the user still indicated that defaults need to be run but have completed the dataframe 
        log = workflow.default_remote_prefix + "/results/{pool}/scrublet_{pctl}/manual_rerun_variables.txt"
        step = "manual"
        sim_dbl = scrublet_extra_dict["sim_dbl"]
        min_counts = scrublet_extra_dict["min_counts"]
        min_cells = scrublet_extra_dict["min_cells"]
        n_prin_comps = scrublet_extra_dict["n_prin_comps"]
        scrublet_doublet_threshold = lambda wildcards: scrublet_manual_dict["scrublet_manual_threshold_thresholds"][scrublet_manual_dict["scrublet_manual_threshold_pools"].index(wildcards.pool)]
    else:
        sys.exit()

    rule scrublet:
        input:
            matrix = lambda wildcards: scrnaseq_libs_df["Matrix_Directories"][wildcards.pool],
            barcodes = lambda wildcards: scrnaseq_libs_df["Barcode_Files"][wildcards.pool],
            df = ancient("results/manual_selections/scrublet/scrublet_percentile_manual_selection.tsv")
        output:
            results = "results/{pool}/scrublet_{pctl}/scrublet_results.txt",
            log = log,
            figure = report("results/{pool}/scrublet_{pctl}/doublet_score_histogram.png", category = "Scrublet", caption = "../report_captions/scrublet.rst", subcategory = "{pool}")
        resources:
            mem_mb = lambda wildcards, attempt: attempt * scrublet_dict["scrublet_memory"],
            disk_mb = lambda wildcards, attempt: attempt * scrublet_dict["scrublet_memory"],
        threads: scrublet_dict["scrublet_threads"]
        params:
            out = workflow.default_remote_prefix + "/results/{pool}/scrublet_{pctl}/",
            script = "/opt/WG1-pipeline-QC_GoogleCloud/Demultiplexing/scripts/scrublet_pipeline.py",
            sim_dbl = sim_dbl,
            min_counts = min_counts,
            min_cells = min_cells,
            n_prin_comps = n_prin_comps,
            scrublet_doublet_threshold = scrublet_doublet_threshold,
            step = step,
            ready = ready
        shell:
            """
            if [ {params.ready} == "True" ]
            then 
                echo "No need to rerun scrublet since the parameters have alreeady been chosen. Will move on to sorting results and merging with results from all other softwares"
            elif [ {params.ready} == "False" ]
            then 
                if [ {params.step} == "default" ]
                then
                    python {params.script} \
                        --counts_matrix {input.matrix} \
                        --barcodes {input.barcodes} \
                        --sim_doublet_ratio {params.sim_dbl} \
                        --min_counts {params.min_counts} \
                        --min_cells {params.min_cells} \
                        --n_prin_comps {params.n_prin_comps} \
                        --min_gene_variability_pctl {wildcards.pctl} \
                        -o {params.out}
                elif [ {params.step} == "manual" ]
                then
                    python {params.script} \
                        --counts_matrix {input.matrix} \
                        --barcodes {input.barcodes} \
                        --sim_doublet_ratio {params.sim_dbl} \
                        --min_counts {params.min_counts} \
                        --min_cells {params.min_cells} \
                        --n_prin_comps {params.n_prin_comps} \
                        --min_gene_variability_pctl {wildcards.pctl} \
                        -o {params.out} \
                        --scrublet_doublet_threshold {params.scrublet_doublet_threshold}
                fi
                echo "The pool:" {wildcards.pool} >> {output.log}
                echo "This was a" {params.step} "run" >> {output.log}
                echo "The number of doublets simulated per droplet:" {params.sim_dbl} >> {output.log}
                echo "The min number of counts used for filtering cells prior to PCA:" {params.min_counts} >> {output.log}
                echo "The number of cells for a gene to be expressed in for filtering cells prior to PCA:" {params.min_cells} >> {output.log}
                echo "The number of principle components used to embed the trnscriptomes prior to k-nearest-neighbor graph:" {params.n_prin_comps} >> {output.log}
                echo "The manual doublet threshold set:" {params.scrublet_doublet_threshold} >> {output.log}
            fi
            [[ -s {output.results} ]]
            echo $?
            """


    rule scrublet_check_user_input:
        input:
            df = ancient("results/manual_selections/scrublet/scrublet_percentile_manual_selection.tsv"),
            results = "results/{pool}/scrublet_{pctl}/scrublet_results.txt"
        output:
            "results/{pool}/CombinedResults/{pctl}_scrublet_results.txt"
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 15000,
            disk_mb=lambda wildcards, attempt: attempt * 15000
        threads: 1
        params:
            ready = ready
        log: "results/logs/scrublet_check_user_input.{pool}_{pctl}.log"
        shell:
            """
            if [ {params.ready} == "True" ]
            then 
                echo "Looks like you put percentile selections into the scrublet_gene_pctl.txt file." 2> {log}
                echo "The scrublet check is done and the next step of the pipeline will proceed" 2>> {log}
                awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' {input.results} > {output}
            elif [ {params.ready} == "False" ]
            then
                echo "You haven't put the scrublet gene percentile selection for each of the pools into the scrublet_gene_pctl.txt" 2> {log}
                echo "Please check the scrublet outputs and choose the best variable genes percentile - rerun any of the pools where the thresholding failed (see the docs) to choose a manual threshold"
                echo "Once you are happy with the thresholding, input the correct gene percentiles (as numbers between 0 and 100) into the second column of the scrublet_gene_pctl.txt file and restart the snakemake pipeline"
                echo 1
            fi
            """
