#!/usr/bin/env python
import os
import pandas as pd
from glob import glob


####################################
############ SOUPORCELL ############
####################################
rule souporcell_unzip_barcodes:
    input:
        barcodes = lambda wildcards: scrnaseq_libs_df["Barcode_Files"][wildcards.pool]
    threads: 1
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * 5,
        disk_per_thread_gb = lambda wildcards, attempt: attempt * 5
    output:
        output_dict["output_dir"] + "/{pool}/souporcell/barcodes.tsv"
    log: output_dict["output_dir"] + "/logs/souporcell_unzip_barcodes.{pool}.log"
    shell:
        """
        if [[ {input.barcodes} == *".gz"* ]]
        then
            gunzip < {input.barcodes} > {output}
        else 
            cp {input.barcodes} {output}
        fi
        """

rule souporcell:
    input:
        bam = lambda wildcards: scrnaseq_libs_df["Bam_Files"][wildcards.pool],
        barcodes = "results/{pool}/souporcell/barcodes.tsv",
        fasta = fasta,
        snps = input_dict["snp_genotypes_filepath"]
    threads: souporcell_dict["souporcell_threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * souporcell_dict["souporcell_memory"],
        disk_mb = lambda wildcards, attempt: attempt * souporcell_dict["souporcell_memory"]
    output:
        clusters = "results/{pool}/souporcell/clusters.tsv",
        genotypes = "results/{pool}/souporcell/cluster_genotypes.vcf",
        troublet = "results/{pool}/souporcell/troublet.done",
        concensus = "results/{pool}/souporcell/consensus.done",
        clustering = "results/{pool}/souporcell/clustering.done"
    params:
        out = workflow.default_remote_prefix + "/results/{pool}/souporcell/",
        N = lambda wildcards: samples.N[samples.Pool == wildcards.pool].iloc[0],
        min_alt = souporcell_extra_dict["min_alt"],
        min_ref = souporcell_extra_dict["min_ref"],
        max_loci = souporcell_extra_dict["max_loci"] 
    log: "results/logs/souporcell.{pool}.log"
    shell:
        """
        souporcell_pipeline.py \
            -i {input.bam} \
            -b {input.barcodes} \
            -f {input.fasta} \
            -t {threads} \
            -o {params.out} \
            -k {params.N} \
            --common_variants {input.snps} \
            --min_alt {params.min_alt} \
            --min_ref {params.min_ref} \
            --max_loci {params.max_loci} 2> {log}
        [[ -s {output.genotypes} ]]
        echo $?
        """

#####################################################
############ REFORMAT SOUPORCELL RESULTS ############
#####################################################
rule souporcell_results_temp:
    input:
        souporcell = "results/{pool}/souporcell/clusters.tsv"
    output:
        "results/{pool}/CombinedResults/souporcell_results.txt"
    resources:
        mem_mb=15000,
        disk_mb=15000
    threads: 1
    log: "results/logs/souporcell_results_temp.{pool}.log"
    shell:
        """
        awk 'BEGIN{{OFS=FS="\\t"}}{{print $1,$2,$3,$4,$5}}' {input.souporcell} | \
            awk 'BEGIN{{FS=OFS="\t"}} $2=="doublet" {{$3="doublet"}}1' | \
            awk 'BEGIN{{FS=OFS="\t"}} $2=="unassigned" {{$4="unassigned"}}1' | \
            sed "s/status/DropletType/g" | sed "s/assignment/Assignment/g" | \
            sed "s/log_prob_singleton/LogProbSinglet/g" | \
            sed "s/log_prob_doublet/LogProbDoublet/g" | \
            sed "s/barcode/Barcode/g" | \
            sed "1s/\t/\tsouporcell_/g" | \
            awk 'NR<2{{print $0;next}}{{print $0| "sort -k1"}}' > {output} 2> {log}
        """

#####################################
############ SUBSET VCFS ############
#####################################
rule souporcell_pool_vcf:
    input:
        genotypes = input_dict["snp_genotypes_filepath"], 
        cluster_geno = "results/{pool}/souporcell/cluster_genotypes.vcf",
        individuals = lambda wildcards: scrnaseq_libs_df["Individuals_Files"][wildcards.pool]
    output:
        filtered_refs = "results/{pool}/souporcell/Individual_genotypes_subset.vcf.gz",
        filtered_refs_temp = "results/{pool}/souporcell/Individual_genotypes_subset.vcf"
    resources:
        mem_mb=15000,
        disk_mb=15000
    threads: 1
    log: "results/logs/souporcell_pool_vcf.{pool}.log"
    shell:
        """
        bedtools intersect -a {input.genotypes} -b {input.cluster_geno} -f 1.0 -r -wa -header > {output.filtered_refs_temp} 2> {log}
        echo "bedtools complete"
        /opt/bcftools-1.10.2/bcftools view -S {input.individuals} -Oz -o {output.filtered_refs} {output.filtered_refs_temp} 2>> {log}
        """


###############################################################################
############ CORRELATE INDIVIDUAL GENOTYPES WITH CLUSTER GENOTYPES ############
###############################################################################
## To take in souporcell_pool_vcf output and souporcell cluster vcf output
rule souporcell_correlate_genotypes:
    input:
        genotypes = "results/{pool}/souporcell/Individual_genotypes_subset.vcf.gz",
        ref = "results/{pool}/souporcell/cluster_genotypes.vcf",
        assignments = "results/{pool}/CombinedResults/CombinedDropletAssignments.tsv"
    output:
        assignments = "results/{pool}/CombinedResults/CombinedDropletAssignments_w_genotypeIDs.tsv",
        variables = temp("results/{pool}/souporcell/souporcel_genotypes_variables"),
        correlation = report("results/{pool}/souporcell/genotype_correlations/pearson_correlation.png", category = "Souporcell Genotype Correlations", subcategory = "{pool}", caption = "../report_captions/souporcell.rst")
    resources:
        mem_mb = souporcell_dict["souporcell_correlations_memory"],
        disk_mb = souporcell_dict["souporcell_correlations_memory"]
    threads: souporcell_dict["souporcell_correlations_threads"]
    params:
        out = workflow.default_remote_prefix + "/results/",
        script = "/opt/WG1-pipeline-QC_GoogleCloud/Demultiplexing/scripts/Assign_Indiv_by_Geno.R",
        cor_thresh = souporcell_dict["souporcell_genotype_correlation_threshold"]
    log: "results/logs/souporcell_correlate_genotypes.{pool}.log"
    shell:
        """
        echo {params.out} > {output.variables}
        echo {wildcards.pool} >> {output.variables}
        echo {input.assignments} >> {output.variables}
        echo {params.cor_thresh} >> {output.variables}
        Rscript {params.script} {output.variables} 2> {log}
        [[ -s {output.assignments} ]]
        echo $?
        """
