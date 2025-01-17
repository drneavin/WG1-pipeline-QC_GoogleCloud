#!/usr/bin/env python
import os
import pandas as pd
from glob import glob

###################################
############# POPSCLE #############
###################################
###### popscle Preprocessing ######
rule popscle_pileup:
    input:
        vcf = input_dict["snp_genotypes_filepath"],
        barcodes = lambda wildcards: scrnaseq_libs_df["Barcode_Files"][wildcards.pool],
        bam = lambda wildcards: scrnaseq_libs_df["Bam_Files"][wildcards.pool],
        files = 'file_directories.txt',
        individuals = lambda wildcards: scrnaseq_libs_df["Individuals_Files"][wildcards.pool]
    output:
        "results/{pool}/popscle/pileup/pileup.var.gz",
        "results/{pool}/popscle/pileup/pileup.cel.gz",
        "results/{pool}/popscle/pileup/pileup.plp.gz",
        "results/{pool}/popscle/pileup/pileup.umi.gz"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * popscle_dict["pileup_memory"],
        disk_mb=lambda wildcards, attempt: attempt * popscle_dict["pileup_memory"]
    threads: popscle_dict["pileup_threads"]
    params:
        out = workflow.default_remote_prefix + "/results/{pool}/popscle/pileup/pileup",
        tag_group = popscle_dict["tag_group"],
        tag_UMI = popscle_dict["tag_UMI"],
        cap_BQ = popscle_extra_dict["cap_BQ"],
        min_BQ = popscle_extra_dict["min_BQ"],
        min_MQ = popscle_extra_dict["min_MQ"],
        min_TD = popscle_extra_dict["min_TD"],
        excl_flag = popscle_extra_dict["excl_flag"],
        min_total = popscle_extra_dict["min_total"],
        min_snp = popscle_extra_dict["min_snp"]
    log: "results/logs/popscle_pileup.{pool}.log"
    shell:
        """
        popscle dsc-pileup \
            --sam {input.bam} \
            --tag-group {params.tag_group} \
            --tag-UMI {params.tag_UMI} \
            --sm-list {input.individuals} \
            --vcf {input.vcf} \
            --cap-BQ {params.cap_BQ} \
            --min-BQ {params.min_BQ} \
            --min-MQ {params.min_MQ} \
            --min-TD {params.min_TD} \
            --excl-flag {params.excl_flag} \
            --group-list {input.barcodes} \
            --min-total {params.min_total} \
            --min-snp {params.min_snp} \
            --out {params.out}
        """

##### Popscle Demuxlet Demultiplexing #####
rule popscle_demuxlet:
    input:
        pileup_cel = "results/{pool}/popscle/pileup/pileup.cel.gz",
        pileup_plp = "results/{pool}/popscle/pileup/pileup.plp.gz",
        pileup_umi = "results/{pool}/popscle/pileup/pileup.umi.gz",
        pileup_var = "results/{pool}/popscle/pileup/pileup.var.gz",
        snps = input_dict["snp_genotypes_filepath"],
        barcodes = lambda wildcards: scrnaseq_libs_df["Barcode_Files"][wildcards.pool],
        individuals = lambda wildcards: scrnaseq_libs_df["Individuals_Files"][wildcards.pool]
    output:
        "results/{pool}/popscle/demuxlet/demuxletOUT.best"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * popscle_dict["demuxlet_memory"],
        disk_mb = lambda wildcards, attempt: attempt * popscle_dict["demuxlet_memory"]
    threads: popscle_dict["demuxlet_threads"]
    params:
        pileup = workflow.default_remote_prefix + "/results/{pool}/popscle/pileup/pileup",
        out = workflow.default_remote_prefix + "/results/{pool}/popscle/demuxlet/",
        field = popscle_dict["genotype_field"],
        geno_error_offset = popscle_extra_dict["geno_error_offset"],
        geno_error_coeff = popscle_extra_dict["geno_error_coeff"],
        r2_info = popscle_extra_dict["r2_info"],
        min_mac = popscle_extra_dict["min_mac"],
        min_callrate = popscle_extra_dict["min_callrate"],
        doublet_prior = popscle_extra_dict["doublet_prior"],
        cap_BQ = popscle_extra_dict["cap_BQ"],
        min_BQ = popscle_extra_dict["min_BQ"],
        min_MQ = popscle_extra_dict["min_MQ"],
        min_TD = popscle_extra_dict["min_TD"],
        excl_flag = popscle_extra_dict["excl_flag"],
        min_total = popscle_extra_dict["min_total"],
        min_snp = popscle_extra_dict["min_snp"]
    log: "results/logs/popscle_demuxlet.{pool}.log"
    shell:
        """
        popscle demuxlet \
            --plp {params.pileup} \
            --vcf {input.snps} \
            --field {params.field} \
            --geno-error-offset {params.geno_error_offset} \
            --geno-error-coeff {params.geno_error_coeff} \
            --r2-info {params.r2_info} \
            --min-mac {params.min_mac} \
            --min-callrate {params.min_callrate} \
            --group-list {input.barcodes} \
            --sm-list {input.individuals} \
            --doublet-prior {params.doublet_prior} \
            --cap-BQ {params.cap_BQ} \
            --min-BQ {params.min_BQ} \
            --min-MQ {params.min_MQ} \
            --min-TD {params.min_TD} \
            --excl-flag {params.excl_flag} \
            --min-total {params.min_total} \
            --min-snp {params.min_snp} \
            --out {params.out}demuxletOUT 
        [[ -s {output} ]]
        echo $?
        """

###################################################
############ REFORMAT DEMUXLET RESULTS ############
###################################################
rule demuxlet_results_temp:
    input:
        demuxlet = "results/{pool}/popscle/demuxlet/demuxletOUT.best"
    output:
        "results/{pool}/CombinedResults/demuxlet_results.txt"
    resources:
        mem_mb=15000,
        disk_mb=15000
    threads: 1
    log: "results/logs/demuxlet_results_temp.{pool}.log"
    shell:
        """
            awk 'BEGIN{{OFS=FS="\\t"}}{{print $2,$3,$5,$13,$14,$19,$20}}' {input.demuxlet} | \
            sed "s/SNG/singlet/g" | \
            sed "s/DBL/doublet/g" | \
            sed "s/AMB/unassigned/g" | \
            awk 'BEGIN{{FS=OFS="\t"}} $3=="doublet" {{$4="doublet"}}1' | \
            sed "s/NUM.SNPS/nSNP/g" | \
            sed "s/DROPLET.TYPE/DropletType/g" | \
            sed "s/singlet.BEST.GUESS/Assignment/g" | \
            sed "s/singlet.BEST.LLK/SingletLLK/g" | \
            sed "s/doublet.BEST.LLK/DoulbetLLK/g" | \
            sed "s/DIFF.LLK.singlet.doublet/DiffLLK/g" | \
            sed "1s/\t/\tdemuxlet_/g" | \
            sed "s/BARCODE/Barcode/g" | \
            awk 'NR<2{{print $0;next}}{{print $0 | "sort -k1"}}'  > {output}
        """