#!/usr/bin/env python
import os
import pandas as pd
from glob import glob

#################################
######## COMBINE RESULTS ########
#################################
if os.path.exists("results/manual_selections/scrublet/scrublet_percentile_manual_selection.tsv"):
    scrublet_selection = pd.read_csv("results/manual_selections/scrublet/scrublet_percentile_manual_selection.tsv", sep = "\t")

rule join_results:
    input:
        demuxlet = "results/{pool}/CombinedResults/demuxlet_results.txt",
        souporcell = "results/{pool}/CombinedResults/souporcell_results.txt",
        scrublet = lambda wildcards: expand("results/{pool}/CombinedResults/{pctl}_scrublet_results.txt", zip, pool = wildcards.pool, pctl = scrublet_selection.scrublet_Percentile[scrublet_selection.Pool == wildcards.pool]),
        scds = "results/{pool}/CombinedResults/scds_results.txt",
        DoubletDetection = "results/{pool}/CombinedResults/DoubletDetection_results.txt"
    output:
        "results/{pool}/CombinedResults/CombinedDropletAssignments.tsv"
    resources:
        mem_mb=15000,
        disk_mb=15000
    threads: 1
    log: "results/logs/join_results.{pool}.log"
    shell:
        """
        join -a1 -a2 -1 1 -2 1 -t "\t" -e"NA" -o "0,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.5" {input.demuxlet} {input.souporcell} | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1"}}' | \
            join -a1 -a2 -1 1 -2 1 -t "\t" -e"NA" -o "0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,2.2,2.3" - {input.scrublet} | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1"}}' | \
            join -a1 -a2 -1 1 -2 1 -t "\t" -e"NA" -o "0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,2.2,2.3" - {input.scds} | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1"}}' | \
            join -a1 -a2 -1 1 -2 1 -t "\t" -e"NA" -o "0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,2.2" - {input.DoubletDetection} > {output}
        """


#####################################################################
############ SCRIPT FOR CALLING FINAL BARCODE ASSIGNMENT ############
#####################################################################
rule final_assignments:
    input:
        assignments = "results/{pool}/CombinedResults/CombinedDropletAssignments_w_genotypeIDs.tsv"
    output:
        figure = report("results/{pool}/CombinedResults/DropletType_Assignment_BarPlot.png", category = "Number Individuals Summary", caption =  "../report_captions/final_assignments.rst"),
        table = "results/{pool}/CombinedResults/Final_Assignments_demultiplexing_doublets.txt",
        variables = "results/{pool}/CombinedResults/variables.txt"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * CombineResults_dict["FinalAssignments_memory"],
        disk_mb = lambda wildcards, attempt: attempt * CombineResults_dict["FinalAssignments_memory"]
    threads: CombineResults_dict["FinalAssignments_threads"]
    params:
        out = workflow.default_remote_prefix + "/results/",
        script = "/opt/WG1-pipeline-QC_GoogleCloud/Demultiplexing/scripts/FinalBarcodeAssignments.R"
    log: "results/logs/final_assignments.{pool}.log"
    shell:
        """
        echo {params.out} > {output.variables}
        echo {wildcards.pool} >> {output.variables}
        Rscript {params.script} {output.variables}
        [[ -s {output.figure} ]]
        echo $?
        """

####################################################
############ SCRIPT TO PRODUCE QC PLOTS ############
####################################################
rule echo_final_assignments:
    input:
        expand("results/{pool}/CombinedResults/Final_Assignments_demultiplexing_doublets.txt", pool = samples.Pool)
    output:
        "results/QC_figures/final_assignments.txt"
    resources:
        mem_mb=15000,
        disk_mb=15000
    threads: 1
    log: "results/logs/echo_final_assignments.log"
    shell:
        """
        echo {input} | tr ' ' '\n' >> {output}
        """

rule final_assignments_check:
    input:
        assignment_list = "results/QC_figures/final_assignments.txt",
        meta=input_dict["samplesheet_filepath"]
    output:
        assignment_list = "results/QC_figures/final_assignments_comparison.txt",
        meta = "results/QC_figures/meta_comparison.txt"
    resources:
        mem_mb=15000,
        disk_mb=15000
    threads: 1
    log: "results/logs/final_assignments_check.log"
    shell:
        """
        cat {input.assignment_list} > {output.assignment_list}
        awk 'BEGIN{{FS=OFS="\t"}}{{print $1}}' {input.meta} | tail -n+2 > {output.meta}
        if [ "$(wc -l < {output.meta})" -eq "$(wc -l < {output.assignment_list})" ]
        then 
            echo 0
        else 
            echo "The number of pools in the final_assignments.txt file don't match the number of pools in your sample sheet" > {log}
            rm {input.assignment_list}
            rm {output.assignment_list}
            rm {output.meta}
            echo 1
        fi
        """


rule expected_observed_numbers:
    input:
        final_assignment = "results/QC_figures/final_assignments.txt",
        sample_sheet = input_dict["samplesheet_filepath"],
        assignemnts = 'results/test_dataset/CombinedResults/Final_Assignments_demultiplexing_doublets.txt',
    output:
        report("results/QC_figures/expected_observed_individuals_classifications.png", category = "Number Individuals Summary", caption = "../report_captions/expected_observed_numbers.rst")
    resources:
        mem_mb = lambda wildcards, attempt: attempt * CombineResults_dict["expected_observed_numbers_memory"],
        disk_mb = lambda wildcards, attempt: attempt * CombineResults_dict["expected_observed_numbers_memory"]
    threads: CombineResults_dict["expected_observed_numbers_threads"]
    params:
        script = "/opt/WG1-pipeline-QC_GoogleCloud/Demultiplexing/scripts/expected_observed_individuals_doublets.R",
        out = workflow.default_remote_prefix + "/results/QC_figures/",
        basedir = workflow.default_remote_prefix + "/results/"
    log: "results/logs/expected_observed_numbers.log"
    shell:
        """
        Rscript {params.script} {input.sample_sheet} {params.out} {params.basedir}
        """


rule QC_plots:
    input:
        assignment_list = "results/QC_figures/final_assignments_comparison.txt",
        pools = "results/QC_figures/meta_comparison.txt",
        script = 'Singlet_QC_Figures.R',
        assignments = "results/test_dataset/CombinedResults/CombinedDropletAssignments_w_genotypeIDs.tsv",
        finals = 'results/test_dataset/CombinedResults/Final_Assignments_demultiplexing_doublets.txt',
        files = 'file_directories.txt',
        barcode_files = barcodes,
        gene_files = genes,
        matrix_files = matrices
    output:
        variables = temp("results/QC_figures/R_variables.txt"),
        fig1 = report("results/QC_figures/nCount_RNA_violin_MAD_All.png", category = "QC", caption = "../report_captions/QC_plots_nCount.rst"),
        fig2 = report("results/QC_figures/nCount_RNA_violin_MADper_Pool.png", category = "QC", caption = "../report_captions/QC_plots_nCount_RNA_MADall.rst"),
        fig3 = report("results/QC_figures/nCount_RNA_violin_noMADlines.png", category = "QC", caption = "../report_captions/QC_plots_nCount.rst"),
        fig4 = report("results/QC_figures/nFeature_RNA_violin_MAD_All.png", category = "QC", caption = "../report_captions/QC_plots_feature_MADall.rst"),
        fig5 = report("results/QC_figures/nFeature_RNA_violin_MADper_Pool.png", category = "QC", caption = "../report_captions/QC_plots_feature_MADperPool.rst"),
        fig6 = report("results/QC_figures/nFeature_RNA_violin_noMADlines.png", category = "QC", caption = "../report_captions/QC_plots_feature.rst"),
        fig7 = report("results/QC_figures/nFeatures_vs_percentMT_QC_scatter_colorPool.png", category = "QC", caption = "../report_captions/QC_plots_features_mt_pool.rst"),
        fig8 = report("results/QC_figures/nFeatures_vs_percentMT_QC_scatter_w_MADlines.png", category = "QC", caption = "../report_captions/QC_plots_features_mt_MAD.rst"),
        fig9 = report("results/QC_figures/percent.mt_violin_MAD_All.png", category = "QC", caption = "../report_captions/QC_plots_mt_MADall.rst"),
        fig10 = report("results/QC_figures/percent.mt_violin_MADper_Pool.png", category = "QC", caption = "../report_captions/QC_plots_mt_MADpool.rst"),
        fig11 = report("results/QC_figures/percent.mt_violin_noMADlines.png", category = "QC", caption = "../report_captions/QC_plots_mt.rst"),
        fig12 = report("results/QC_figures/percent.rb_violin_MAD_All.png", category = "QC", caption = "../report_captions/QC_plots_rb_MADall.rst"),
        fig13 = report("results/QC_figures/percent.rb_violin_MADper_Pool.png", category = "QC", caption = "../report_captions/QC_plots_rb_MADperPool.rst"),
        fig14 = report("results/QC_figures/percent.rb_violin_noMADlines.png", category = "QC", caption = "../report_captions/QC_plots_rb.rst"),
        fig15 = report("results/QC_figures/UMI_vs_Genes_QC_scatter.png", category = "QC", caption = "../report_captions/QC_plots_UMI_features.rst"),
        fig16 = report("results/QC_figures/UMI_vs_Genes_QC_scatter_w_MADlines.png", category = "QC", caption = "../report_captions/QC_plots_UMI_features_MADall.rst"),
        fig17 = report("results/QC_figures/UMI_vs_percentMT_QC_scatter_colorPool.png", category = "QC", caption = "../report_captions/QC_plots_UMI_mt_pool.rst"),
        fig18 = report("results/QC_figures/UMI_vs_percentMT_QC_scatter_w_MADlines.png", category = "QC", caption = "../report_captions/QC_plots_UMI_mt_MDA_all.rst"),
    resources:
        mem_mb = lambda wildcards, attempt: attempt * CombineResults_dict["FinalQC_memory"],
        disk_mb = lambda wildcards, attempt: attempt * CombineResults_dict["FinalQC_memory"]
    threads: CombineResults_dict["FinalQC_threads"]
    params:
        script = "/opt/WG1-pipeline-QC_GoogleCloud/Demultiplexing/scripts/Singlet_QC_Figures.R",
        main_dir = workflow.default_remote_prefix + "/results/",
        out = workflow.default_remote_prefix + "/results/QC_figures/",
        rb_genes = "/opt/WG1-pipeline-QC_GoogleCloud/Demultiplexing/Ribosomal_genes.txt",
        mt_genes = "/opt/WG1-pipeline-QC_GoogleCloud/Demultiplexing/Mitochondrial_genes.txt",
        remote_prefix = workflow.default_remote_prefix
    log: "results/logs/QC_plots.log"
    shell:
        """
        echo {params.main_dir} > {output.variables}
        echo {input.pools} >> {output.variables}
        echo {input.files} >> {output.variables}
        echo {params.out} >> {output.variables}
        echo {params.rb_genes} >> {output.variables}
        echo {params.mt_genes} >> {output.variables}
        echo {params.remote_prefix} >> {output.variables}
        Rscript {params.script} {output.variables}
        [[ -s {output.fig18} ]]
        echo $?
        """
