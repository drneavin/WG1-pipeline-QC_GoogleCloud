#!/usr/bin/env python
shell.executable('bash')


rule indiv_missingness:
    input:
        pgen = pgen,
        pvar = pvar,
        psam = psam,
    output:
        bed = "results/indiv_missingness/indiv_missingness.pgen",
        bim = "results/indiv_missingness/indiv_missingness.pvar",
        fam = "results/indiv_missingness/indiv_missingness.psam",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["indiv_missingness_memory"],
        disk_mb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["indiv_missingness_memory"]
    threads: plink_gender_ancestry_QC_dict["indiv_missingness_threads"]
    params:
       infile = workflow.default_remote_prefix + "/" + re.sub(".psam", "", psam),
       out = workflow.default_remote_prefix + "/results/indiv_missingness/indiv_missingness",
       mind = plink_gender_ancestry_QC_dict["indiv_missingness_mind"],
    shell:
        """
        echo {params.infile}
        plink2 --threads {threads} \
            --pfile {params.infile} \
            --make-pgen 'psam-cols='fid,parents,sex,phenos \
            --mind {params.mind} \
            --out {params.out}
        """

rule check_sex:
    input:
        bed = "results/indiv_missingness/indiv_missingness.pgen",
        bim = "results/indiv_missingness/indiv_missingness.pvar",
        fam = "results/indiv_missingness/indiv_missingness.psam",
    output:
        bed = "results/check_sex/check_sex.bed",
        bim = "results/check_sex/check_sex.bim",
        fam = "results/check_sex/check_sex.fam",
        hh = "results/check_sex/check_sex.hh",
        log = "results/check_sex/check_sex.log",
        nosex = "results/check_sex/check_sex.nosex",
        sexcheck = "results/check_sex/check_sex.sexcheck",
        sexcheck_tab = "results/check_sex/check_sex.sexcheck.tsv"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["check_sex_memory"],
        disk_mb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["check_sex_memory"]
    threads: plink_gender_ancestry_QC_dict["check_sex_threads"]
    params:
        infile = workflow.default_remote_prefix + "/results/indiv_missingness/indiv_missingness",
        out = workflow.default_remote_prefix + "/results/check_sex/check_sex",
    shell:
        """
        plink2 --threads {threads} --pfile {params.infile} --make-bed --max-alleles 2 --out {params.out}
        plink --threads {threads} --bfile {params.out} --check-sex --out {params.out}
        touch {output.nosex}
        sed 's/^ \+//g' {output.sexcheck} | sed 's/ \+/\t/g' > {output.sexcheck_tab}
        """

### Pull just common SNPs between two groups ###
rule common_snps:
    input:
        bed = "results/indiv_missingness/indiv_missingness.pgen",
        bim = "results/indiv_missingness/indiv_missingness.pvar",
        fam = "results/indiv_missingness/indiv_missingness.psam",
    output:
        snps_data = "results/common_snps/snps_data.tsv",
        snps_1000g = "results/common_snps/snps_1000g.tsv",
        bed = "results/common_snps/subset_data.pgen",
        bim = "results/common_snps/subset_data.pvar",
        fam = "results/common_snps/subset_data.psam",
        bed_1000g = "results/common_snps/subset_1000g.pgen",
        bim_1000g = "results/common_snps/subset_1000g.pvar",
        fam_1000g = "results/common_snps/subset_1000g.psam",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["common_snps_memory"],
        disk_mb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["common_snps_memory"]
    threads: plink_gender_ancestry_QC_dict["common_snps_threads"]
    params:
        bim_1000 = "/opt/1000G/all_phase3_filtered.pvar",
        infile = workflow.default_remote_prefix + "/results/indiv_missingness/indiv_missingness",
        infile_1000g = "/opt/1000G/all_phase3_filtered",
        out = workflow.default_remote_prefix + "/results/common_snps/subset_data",
        out_1000g = workflow.default_remote_prefix + "/results/common_snps/subset_1000g",
    shell:
        """
        awk 'NR==FNR{{a[$1,$2,$4,$5];next}} ($1,$2,$4,$5) in a{{print $3}}' {input.bim} {params.bim_1000} > {output.snps_1000g}
        awk 'NR==FNR{{a[$1,$2,$4,$5];next}} ($1,$2,$4,$5) in a{{print $3}}' {params.bim_1000} {input.bim} > {output.snps_data}
        plink2 --threads {threads} --pfile {params.infile} --extract {output.snps_data} --make-pgen 'psam-cols='fid,parents,sex,phenos --out {params.out}
        plink2 --threads {threads} --pfile {params.infile_1000g} --extract {output.snps_1000g} --make-pgen --out {params.out_1000g}
        """

### Prune with --indep,
rule prune_1000g:
    input:
        bed_1000g = "results/common_snps/subset_1000g.pgen",
        bim_1000g = "results/common_snps/subset_1000g.pvar",
        fam_1000g = "results/common_snps/subset_1000g.psam",
        bim = "results/common_snps/subset_data.pvar",
        bed = "results/common_snps/subset_data.pgen",
        fam = "results/common_snps/subset_data.psam",
    output:
        prune_out_1000g = "results/common_snps/subset_pruned_1000g.prune.out",
        prune_out = "results/common_snps/subset_data.prune.out",
        bed_1000g = "results/common_snps/subset_pruned_1000g.pgen",
        bim_1000g = "results/common_snps/subset_pruned_1000g.pvar",
        fam_1000g = "results/common_snps/subset_pruned_1000g.psam",
        bed = "results/common_snps/subset_pruned_data.pgen",
        bim = "results/common_snps/subset_pruned_data.pvar",
        bim_temp = "results/common_snps/subset_pruned_data_temp.pvar",
        bim_old = "results/common_snps/subset_pruned_data_original.pvar",
        fam = "results/common_snps/subset_pruned_data.psam",
        data_1000g_key = "results/common_snps/subset_pruned_data_1000g_key.txt",
        SNPs2keep = "results/common_snps/SNPs2keep.txt"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["prune_1000g_memory"],
        disk_mb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["prune_1000g_memory"]
    threads: plink_gender_ancestry_QC_dict["prune_1000g_threads"]
    params:
        out_1000g = workflow.default_remote_prefix + "/results/common_snps/subset_pruned_1000g",
        infile_1000g = workflow.default_remote_prefix + "/results/common_snps/subset_1000g",
        infile = workflow.default_remote_prefix + "/results/common_snps/subset_data",
        out = workflow.default_remote_prefix + "/results/common_snps/subset_pruned_data"
    shell:
        """
        plink2 --threads {threads} --pfile {params.infile_1000g} \
            --indep-pairwise 50 5 0.5 \
            --out {params.out_1000g}
        plink2 --threads {threads} --pfile {params.infile_1000g} --extract {output.prune_out_1000g} --make-pgen --out {params.out_1000g}
        if [[ $(grep "##" {input.bim} | wc -l) > 0 ]]
        then
            grep "##" {input.bim} > {output.data_1000g_key}
        fi
        awk -F"\\t" 'BEGIN{{OFS=FS = "\\t"}} NR==FNR{{a[$1 FS $2 FS $4 FS $5] = $0; next}} {{ind = $1 FS $2 FS $4 FS $5}} ind in a {{print a[ind], $3}}' {output.bim_1000g} {input.bim} | grep -v "##" >> {output.data_1000g_key}
        grep -v "##" {output.data_1000g_key} | awk 'BEGIN{{FS=OFS="\t"}}{{print $NF}}' > {output.prune_out}
        plink2 --threads {threads} --pfile {params.infile} --extract {output.prune_out} --make-pgen 'psam-cols='fid,parents,sex,phenos --out {params.out}
        cp {output.bim} {output.bim_old}
        grep -v "#" {output.bim_old} | awk 'BEGIN{{FS=OFS="\t"}}{{print($3)}}' > {output.SNPs2keep}
        grep "#CHROM" {output.data_1000g_key} > {output.bim}
        grep -Ff {output.SNPs2keep} {output.data_1000g_key} >> {output.bim}
        awk 'BEGIN{{FS=OFS="\t"}}NF{{NF-=1}};1' < {output.bim} > {output.bim_temp}
        grep "##" {output.bim_1000g} > {output.bim}
        cat {output.bim_temp} >> {output.bim}
        """
        
rule final_pruning: ### put in contingency for duplicated snps - remove from both 1000G and your dataset
    input:
        bed = "results/common_snps/subset_pruned_data.pgen",
        bim = "results/common_snps/subset_pruned_data.pvar",
        fam = "results/common_snps/subset_pruned_data.psam",
    output:
        bed = "results/common_snps/final_subset_pruned_data.pgen",
        bim = "results/common_snps/final_subset_pruned_data.pvar",
        fam = "results/common_snps/final_subset_pruned_data.psam",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["prune_1000g_memory"],
        disk_mb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["prune_1000g_memory"]
    threads: plink_gender_ancestry_QC_dict["prune_1000g_threads"]
    params:
        infile = workflow.default_remote_prefix + "/results/common_snps/subset_pruned_data",
        out = workflow.default_remote_prefix + "/results/common_snps/final_subset_pruned_data"
    shell:
        """
        plink2 --rm-dup 'force-first' -threads {threads} --pfile {params.infile} --make-pgen 'psam-cols='fid,parents,sex,phenos --out {params.out}
        """


### use PCA from plink for PCA and projection
rule pca_1000g:
    input:
        bed_1000g = "results/common_snps/subset_pruned_1000g.pgen",
        bim_1000g = "results/common_snps/subset_pruned_1000g.pvar",
        fam_1000g = "results/common_snps/subset_pruned_1000g.psam",
        bed = "results/common_snps/subset_pruned_data.pgen" 
    output:
        out = "results/pca_projection/subset_pruned_1000g_pcs.acount",
        eig_all = "results/pca_projection/subset_pruned_1000g_pcs.eigenvec.allele",
        eig_vec = "results/pca_projection/subset_pruned_1000g_pcs.eigenvec",
        eig = "results/pca_projection/subset_pruned_1000g_pcs.eigenval",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["pca_1000g_memory"],
        disk_mb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["pca_1000g_memory"]
    threads: plink_gender_ancestry_QC_dict["pca_1000g_threads"]
    params:
        infile = workflow.default_remote_prefix + "/results/common_snps/subset_pruned_1000g",
        out = workflow.default_remote_prefix + "/results/pca_projection/subset_pruned_1000g_pcs"
    shell:
        """
        plink2 --threads {threads} --pfile {params.infile} \
            --freq counts \
            --pca allele-wts \
            --out {params.out}
        """


### use plink pca results to plot with R ###
rule pca_project:
    input:
        bed = "results/common_snps/final_subset_pruned_data.pgen",
        bim = "results/common_snps/final_subset_pruned_data.pvar",
        fam = "results/common_snps/final_subset_pruned_data.psam",
        frq = "results/pca_projection/subset_pruned_1000g_pcs.acount",
        scores = "results/pca_projection/subset_pruned_1000g_pcs.eigenvec.allele",
        bed_1000g = "results/common_snps/subset_pruned_1000g.pgen",
        bim_1000g = "results/common_snps/subset_pruned_1000g.pvar",
        fam_1000g = "results/common_snps/subset_pruned_1000g.psam"
    output:
        projected_scores = "results/pca_projection/final_subset_pruned_data_pcs.sscore",
        projected_1000g_scores = "results/pca_projection/subset_pruned_1000g_pcs_projected.sscore"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["pca_project_memory"],
        disk_mb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["pca_project_memory"]
    threads: plink_gender_ancestry_QC_dict["pca_project_threads"]
    params:
        infile = workflow.default_remote_prefix + "/results/common_snps/final_subset_pruned_data",
        infile_1000g = workflow.default_remote_prefix + "/results/common_snps/subset_pruned_1000g",
        out = workflow.default_remote_prefix + "/results/pca_projection/final_subset_pruned_data_pcs",
        out_1000g = workflow.default_remote_prefix + "/results/pca_projection/subset_pruned_1000g_pcs_projected"
    shell:
        """
        export OMP_NUM_THREADS={threads}
        plink2 --threads {threads} --pfile {params.infile} \
            --read-freq {input.frq} \
            --score {input.scores} 2 5 header-read no-mean-imputation \
                    variance-standardize \
            --score-col-nums 6-15 \
            --out {params.out}
        plink2 --threads {threads} --pfile {params.infile_1000g} \
            --read-freq {input.frq} \
            --score {input.scores} 2 5 header-read no-mean-imputation \
                    variance-standardize \
            --score-col-nums 6-15 \
            --out {params.out_1000g}
       """

rule pca_projection_assign:
    input:
        projected_scores = "results/pca_projection/final_subset_pruned_data_pcs.sscore",
        projected_1000g_scores = "results/pca_projection/subset_pruned_1000g_pcs_projected.sscore",
        fam_1000g = "results/common_snps/subset_1000g.psam",
        psam = psam,
        sexcheck = "results/check_sex/check_sex.sexcheck.tsv",
    output:
        sexcheck = "results/pca_sex_checks/check_sex_update_remove.tsv",
        anc_check = "results/pca_sex_checks/ancestry_update_remove.tsv",
        anc_fig = report("results/pca_sex_checks/Ancestry_PCAs.png", category = "Ancestry", caption = "../report_captions/ancestry_pca.rst")
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 15000,
        disk_mb=lambda wildcards, attempt: attempt * 15000
    threads: 2
    params:
        variables = workflow.default_remote_prefix + "/results/pca_sex_checks/variables.tsv",
        outdir = workflow.default_remote_prefix + "/results/pca_sex_checks/",
        script = "/opt/WG1-pipeline-QC/Imputation/scripts/PCA_Projection_Plotting.R"
    shell:
        """
        echo {params.outdir} > {params.variables}
        echo {input.projected_scores} >> {params.variables}
        echo {input.projected_1000g_scores} >> {params.variables}
        echo {input.fam_1000g} >> {params.variables}
        echo {input.psam} >> {params.variables}
        echo {input.sexcheck} >> {params.variables}
        Rscript {params.script} {params.variables}
        """


rule summary_ancestry_sex:
    input:
        sexcheck = "results/check_sex/check_sex.sexcheck.tsv",
        sexcheck_tsv = "results/pca_sex_checks/check_sex_update_remove.tsv",
        fam = "results/indiv_missingness/indiv_missingness.psam",
        anc_check = "results/pca_sex_checks/ancestry_update_remove.tsv"
    output:
        report("results/metrics/sex_summary.png", category = "Ancestry and Sex Summary", caption = "../report_captions/sex_summary.rst"),
        report("results/metrics/ancestry_summary.png", category = "Ancestry and Sex Summary", caption = "../report_captions/ancestry_summary.rst")
    resources:
        mem_mb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["summary_ancestry_sex_memory"],
        disk_mb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["summary_ancestry_sex_memory"]
    threads: plink_gender_ancestry_QC_dict["summary_ancestry_sex_threads"]
    params:
        outdir = workflow.default_remote_prefix + "/results/metrics/",
        basedir = workflow.default_remote_prefix + "/results",
        script = "/opt/WG1-pipeline-QC/Imputation/scripts/sex_ancestry_summaries.R"
    shell:
        """
        Rscript {params.script} {params.basedir} {params.outdir}
        """



rule separate_indivs:
    input:
        sexcheck = "results/pca_sex_checks/check_sex_update_remove.tsv",
        anc_check = "results/pca_sex_checks/ancestry_update_remove.tsv"
    output:
        update_sex = "results/separate_indivs/sex_update_indivs.tsv",
        remove_indiv = "results/separate_indivs/remove_indivs.tsv",
        remove_indiv_temp = "results/separate_indivs/remove_indivs_temp.tsv",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 15000,
        disk_mb=lambda wildcards, attempt: attempt * 15000
    threads: 1
    shell:
        """
        grep "UPDATE" {input.sexcheck} | awk 'BEGIN{{FS=OFS="\t"}}{{print($1,$2,$4)}}' | sed 's/SNPSEX/SEX/g' > {output.update_sex}
        grep "REMOVE" {input.sexcheck} | awk 'BEGIN{{FS=OFS="\t"}}{{print($1,$2)}}'> {output.remove_indiv_temp}
        grep "REMOVE" {input.anc_check} | awk 'BEGIN{{FS=OFS="\t"}}{{print($1,$2)}}' >> {output.remove_indiv_temp}
        sort -u {output.remove_indiv_temp} > {output.remove_indiv}
        """


rule update_sex_ancestry:
    input:
        bim = "results/indiv_missingness/indiv_missingness.pgen",
        psam_indiv = "results/indiv_missingness/indiv_missingness.psam",
        pvar = "results/indiv_missingness/indiv_missingness.pvar",
        psam = psam,
        update_sex = "results/separate_indivs/sex_update_indivs.tsv",
        remove_indiv = "results/separate_indivs/remove_indivs.tsv",
        psam_updated = "results/pca_sex_checks/updated_psam.psam",
    output:
        bed = "results/update_sex_ancestry/update_sex.pgen",
        bim = "results/update_sex_ancestry/update_sex.pvar",
        psam = "results/update_sex_ancestry/update_sex.psam",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["update_sex_memory"],
        disk_mb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["update_sex_memory"]
    threads: plink_gender_ancestry_QC_dict["update_sex_threads"]
    params:
        anc_updated_psam = workflow.default_remote_prefix + "/results/pca_sex_checks/updated_psam.psam",
        infile = workflow.default_remote_prefix + "/results/indiv_missingness/indiv_missingness",
        psam_temp = workflow.default_remote_prefix + "/results/update_sex_ancestry/temp/indiv_missingness.psam_temp",
        tdir = workflow.default_remote_prefix + "/results/update_sex_ancestry/temp/",
        out = workflow.default_remote_prefix + "/results/update_sex_ancestry/update_sex",
    shell:
        """
        mkdir -p {params.tdir}
        cp {params.infile}* {params.tdir}
        cp {params.anc_updated_psam} {params.tdir}/indiv_missingness.psam 
        plink2 --threads {threads} --pfile {params.tdir}/indiv_missingness --update-sex {input.update_sex} --remove {input.remove_indiv} --make-pgen 'psam-cols='fid,parents,sex,phenos --out {params.out}
        """


rule subset_ancestry:
    input:
        psam = "results/update_sex_ancestry/update_sex.psam",
        pgen = "results/update_sex_ancestry/update_sex.pgen",
        pvar = "results/update_sex_ancestry/update_sex.pvar"
    output:
        keep = "results/subset_ancestry/{ancestry}_individuals.psam",
        pgen = "results/subset_ancestry/{ancestry}_subset.pgen",
        psam = "results/subset_ancestry/{ancestry}_subset.psam",
        pvar = "results/subset_ancestry/{ancestry}_subset.pvar"
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 15000,
        disk_mb = lambda wildcards, attempt: attempt * 15000
    threads: 1
    params:
        infile = workflow.default_remote_prefix + "/results/update_sex_ancestry/update_sex",
        out = workflow.default_remote_prefix + "/results/subset_ancestry/{ancestry}_subset"
    shell:
        """
        grep {wildcards.ancestry} {input.psam} > {output.keep}
        plink2 --threads {threads} --pfile {params.infile} --keep {output.keep} --max-alleles 2 --make-pgen 'psam-cols='fid,parents,sex,phenos --out {params.out}
        """


