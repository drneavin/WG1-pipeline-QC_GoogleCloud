#!/usr/bin/env python
import os
import pandas as pd
from glob import glob

###########################################
############ DOUBLET DETECTION ############
###########################################
rule make_DoubletDetection_selection_df:
    input:
        input_dict["samplesheet_filepath"]
    output:
        output_dict["output_dir"] + "/manual_selections/DoubletDetection/DoubletDetection_manual_selection.tsv"
    resources:
        mem_per_thread_gb = 1,
        disk_per_thread_gb = 1
    threads: 1
    log: output_dict["output_dir"] + "/logs/make_DoubletDetection_selection_df.log"
    shell:
        """
        awk 'BEGIN{{OFS=FS="\t"}}{{print $1 "\t"}}' {input} | sed "1s/.*/Pool\tDoubletDetection_PASS_FAIL/" > {output} 2> {log}
        """

if os.path.exists(output_dict["output_dir"] + "/manual_selections/DoubletDetection/DoubletDetection_manual_selection.tsv"):
    DoubletDetection_selection = pd.read_csv(output_dict["output_dir"] + "/manual_selections/DoubletDetection/DoubletDetection_manual_selection.tsv", sep = "\t")
    if len(DoubletDetection_selection[DoubletDetection_selection['DoubletDetection_PASS_FAIL'].astype(str).str.contains('PASS', na=False)]) != len(DoubletDetection_selection):
        ready = False
    elif len(DoubletDetection_selection[DoubletDetection_selection['DoubletDetection_PASS_FAIL'].astype(str).str.contains('PASS', na=False)]) == len(DoubletDetection_selection):
        ready = True
        step = "ready"
    else:
        sys.exit()
    if DoubletDetection_manual_dict["run_DoubletDetection_manual"] == False or (DoubletDetection_manual_dict["run_DoubletDetection_manual"] == True and len(DoubletDetection_selection[DoubletDetection_selection['DoubletDetection_PASS_FAIL'].astype(str).str.contains('PASS', na=False)]) == len(DoubletDetection_selection)):
        step = "default"
        log = output_dict["output_dir"] + "/{pool}/DoubletDetection/default_run_variables.txt"
        n_iterations = DoubletDetection_extra_dict["n_iterations"]
        phenograph = DoubletDetection_extra_dict["phenograph"]
        standard_scaling = DoubletDetection_extra_dict["standard_scaling"]
        p_thresh = DoubletDetection_extra_dict["p_thresh"]
        voter_thresh = DoubletDetection_extra_dict["voter_thresh"]
    elif DoubletDetection_manual_dict["run_DoubletDetection_manual"] == True:
        step = "manual"
        log = output_dict["output_dir"] + "/{pool}/DoubletDetection/manual_rerun_variables.txt"
        n_iterations = DoubletDetection_manual_dict["n_iterations"]
        phenograph = DoubletDetection_manual_dict["phenograph"]
        standard_scaling = DoubletDetection_manual_dict["standard_scaling"]
        p_thresh = DoubletDetection_manual_dict["p_thresh"]
        voter_thresh = DoubletDetection_manual_dict["voter_thresh"]
    else:
        sys.exit()
        
    rule DoubletDetection:
        input:
            barcodes = lambda wildcards: scrnaseq_libs_df["Barcode_Files"][wildcards.pool],
            matrix = lambda wildcards: scrnaseq_libs_df["Matrix_Files"][wildcards.pool],
            df = ancient(output_dict["output_dir"] + "/manual_selections/DoubletDetection/DoubletDetection_manual_selection.tsv")
        output:
            doublets = output_dict["output_dir"] + "/{pool}/DoubletDetection/DoubletDetection_results.txt",
            figure = report(output_dict["output_dir"] + "/{pool}/DoubletDetection/convergence_test.pdf", category = "DoubletDetection", subcategory = "{pool}", caption = "../report_captions/DoubletDetection.rst"),
            log = log
        resources:
            mem_per_thread_gb = lambda wildcards, attempt: attempt * DoubletDetection_dict["DoubletDetection_memory"],
            disk_per_thread_gb = lambda wildcards, attempt: attempt * DoubletDetection_dict["DoubletDetection_memory"]
        threads: DoubletDetection_dict["DoubletDetection_threads"]
        params:
            script = "/opt/WG1-pipeline-QC_GoogleCloud/Demultiplexing/scripts/DoubletDetection.py",
            out = output_dict["output_dir"] + "/{pool}/DoubletDetection/",
            n_iterations = n_iterations,
            phenograph = phenograph,
            standard_scaling = standard_scaling,
            p_thresh = p_thresh,
            voter_thresh = voter_thresh,
            ready = ready,
            step = step
        shell:
            """
            if [ {params.ready} == "True" ]
            then 
                echo "No need to rerun DoubletDetection since the parameters have alreeady been chosen. Will move on to sorting results and merging with results from all other softwares"
            elif [ {params.ready} == "False" ]
            then 
                python {params.script} \
                    --counts_matrix {input.matrix} \
                    --barcodes {input.barcodes} \
                    --n_iterations {params.n_iterations} \
                    --phenograph {params.phenograph} \
                    --standard_scaling {params.standard_scaling} \
                    --p_thresh {params.p_thresh} \
                    --voter_thresh {params.voter_thresh} \
                    -o {params.out}
				echo "The pool:" {wildcards.pool} >> {output.log}
				echo "This was a" {params.step} "run" >> {output.log}
				echo "The number of iterations used to determine doublets:" {params.n_iterations} >> {output.log}
				echo "The phenograph was was used:" {params.phenograph} >> {output.log}
				echo "The standard scaling was used:" {params.standard_scaling} >> {output.log}
				echo "The p threshold was used:" {params.p_thresh} >> {output.log}
				echo "The voter threshold is:" {params.voter_thresh} >> {output.log}
            fi
            """
            
    rule DoubletDetection_check_user_input:
        input:
            results = output_dict["output_dir"] + "/{pool}/DoubletDetection/DoubletDetection_results.txt",
            df = ancient(output_dict["output_dir"] + "/manual_selections/DoubletDetection/DoubletDetection_manual_selection.tsv")
        output:
            output_dict["output_dir"] + "/{pool}/CombinedResults/DoubletDetection_results.txt"
        resources:
            mem_per_thread_gb=lambda wildcards, attempt: attempt * 1,
            disk_per_thread_gb=lambda wildcards, attempt: attempt * 1
        threads: 1
        params:
            ready = ready
        log: output_dict["output_dir"] + "/logs/DoubletDetection_check_user_input.{pool}.log"
        shell:
            """
            if [ {params.ready} == "True" ]
            then 
                echo "Looks like you put PASS into all the rows of DoubletDetection_manual_selection.tsv file." > {log}
                echo "The DoubletDetection check is done and the next step of the pipeline will proceed" >> {log}
                awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' {input.results} > {output}
            elif [ {params.ready} == "False" ]
            then
                echo "You haven't put PASS/FAIL values into the DoubletDetection_manual_selection.tsv file." > {log}
                echo "Please check the DoubletDetection outputs and decide if the pools passed - rerun any of the pools where the doublet numbers don't reach convergence using the manual selections (see the docs)" >> {log}
                echo "Once you are happy with the results, input PASS into the second column of the DoubletDetection_manual_selection.tsv file and restart the snakemake pipeline" >> {log}
                echo 1
            fi
            """