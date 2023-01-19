
#!/usr/bin/env python
shell.executable('bash')



# Converts BIM to BED and converts the BED file via CrossMap.
# Finds excluded SNPs and removes them from the original plink file.
# Then replaces the BIM with CrossMap's output.
rule crossmap:
    input:
        pgen = "results/subset_ancestry/{ancestry}_subset.pgen",
        psam = "results/subset_ancestry/{ancestry}_subset.psam",
        pvar = "results/subset_ancestry/{ancestry}_subset.pvar"
    output:
        bed = "results/crossmapped/{ancestry}_crossmapped_plink.bed",
        bim = "results/crossmapped/{ancestry}_crossmapped_plink.bim",
        fam = "results/crossmapped/{ancestry}_crossmapped_plink.fam",
        inbed = "results/crossmapped/{ancestry}_crossmap_input.bed",
        outbed = "results/crossmapped/{ancestry}_crossmap_output.bed",
        excluded_ids = "results/crossmapped/{ancestry}_excluded_ids.txt"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * imputation_dict["crossmap_memory"],
        disk_mb=lambda wildcards, attempt: attempt * imputation_dict["crossmap_memory"]
    threads:
        imputation_dict["crossmap_threads"]
    params:
        in_plink = workflow.default_remote_prefix + "/results/subset_ancestry/{ancestry}_subset",
        out = workflow.default_remote_prefix + "/results/crossmapped/{ancestry}_crossmapped_plink",
        chain_file = "/opt/GRCh37_to_GRCh38.chain"
    shell:
        """
        awk 'BEGIN{{FS=OFS="\t"}}{{print $1,$2,$2+1,$3,$4,$5}}' {input.pvar} > {output.inbed}
        CrossMap.py bed {params.chain_file} {output.inbed} {output.outbed}
        awk '{{print $4}}' {output.outbed}.unmap > {output.excluded_ids}
        plink2 --pfile {params.in_plink} --exclude {output.excluded_ids} --make-bed --output-chr MT --out {params.out}
        awk -F'\t' 'BEGIN {{OFS=FS}} {{print $1,$4,0,$2,$6,$5}}' {output.outbed} > {output.bim}
        """

rule sort_bed:
    input:
        pgen = "results/crossmapped/{ancestry}_crossmapped_plink.bed",
        psam = "results/crossmapped/{ancestry}_crossmapped_plink.bim",
        pvar = "results/crossmapped/{ancestry}_crossmapped_plink.fam"
    output:
        bed = "results/crossmapped_sorted/{ancestry}_crossmapped_sorted.bed",
        bim = "results/crossmapped_sorted/{ancestry}_crossmapped_sorted.bim",
        fam = "results/crossmapped_sorted/{ancestry}_crossmapped_sorted.fam"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * imputation_dict["sort_bed_memory"],
        disk_mb=lambda wildcards, attempt: attempt * imputation_dict["sort_bed_memory"]
    threads:
        imputation_dict["sort_bed_threads"]
    params:
        infile = workflow.default_remote_prefix + "/results/crossmapped/{ancestry}_crossmapped_plink",
        out = workflow.default_remote_prefix + "/results/crossmapped_sorted/{ancestry}_crossmapped_sorted"
    shell:
        """
        plink2 --bfile {params.infile} --make-bed --max-alleles 2 --output-chr MT --out {params.out}
        """


rule harmonize_hg38:
    input:
        bed = "results/crossmapped_sorted/{ancestry}_crossmapped_sorted.bed",
        bim = "results/crossmapped_sorted/{ancestry}_crossmapped_sorted.bim",
        fam = "results/crossmapped_sorted/{ancestry}_crossmapped_sorted.fam",
        vcf = vcf_dir,
        index = vcf_dir + ".tbi"
    output:
        bed = "results/harmonize_hg38/{ancestry}.bed",
        bim = "results/harmonize_hg38/{ancestry}.bim",
        fam = "results/harmonize_hg38/{ancestry}.fam"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 15000,
        java_mem = lambda wildcards, attempt: attempt * imputation_dict["harmonize_hg38_java_memory"],
        disk_mb=lambda wildcards, attempt: attempt * imputation_dict["harmonize_hg38_memory"]
    threads:
        imputation_dict["harmonize_hg38_threads"]
    params:
        infile = workflow.default_remote_prefix + "/results/crossmapped_sorted/{ancestry}_crossmapped_sorted",
        out = workflow.default_remote_prefix + "/results/harmonize_hg38/{ancestry}",
        jar = "/opt/GenotypeHarmonizer-1.4.23/GenotypeHarmonizer.jar"
    shell:
        """
        java -Xmx{resources.java_mem}M -jar {params.jar}\
            --input {params.infile}\
            --inputType PLINK_BED\
            --ref {input.vcf}\
            --refType VCF\
            --update-id\
            --output {params.out}
        """


rule plink_to_vcf:
    input:
        bed = "results/harmonize_hg38/{ancestry}.bed",
        bim = "results/harmonize_hg38/{ancestry}.bim",
        fam = "results/harmonize_hg38/{ancestry}.fam"
    output:
        data_vcf_gz = "results/harmonize_hg38/{ancestry}_harmonised_hg38.vcf.gz",
        index = "results/harmonize_hg38/{ancestry}_harmonised_hg38.vcf.gz.csi"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * imputation_dict["plink_to_vcf_memory"],
        disk_mb=lambda wildcards, attempt: attempt * imputation_dict["plink_to_vcf_memory"]
    threads:
        imputation_dict["plink_to_vcf_threads"]
    params:
        infile = workflow.default_remote_prefix + "/results/harmonize_hg38/{ancestry}",
        out = workflow.default_remote_prefix + "/results/harmonize_hg38/{ancestry}_harmonised_hg38"
    shell:
        """
        plink2 --bfile {params.infile} --recode vcf id-paste=iid --chr 1-22 --out {params.out}

        bgzip {params.out}.vcf
        bcftools index {output.data_vcf_gz}
        """


rule vcf_fixref_hg38:
    input:
        fasta = fasta,
        vcf = vcf_dir,
        index = vcf_dir + ".tbi",
        data_vcf = "results/harmonize_hg38/{ancestry}_harmonised_hg38.vcf.gz"
    output:
        vcf = "results/vcf_fixref_hg38/{ancestry}_fixref_hg38.vcf.gz",
        index = "results/vcf_fixref_hg38/{ancestry}_fixref_hg38.vcf.gz.csi"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * imputation_dict["vcf_fixref_hg38_memory"],
        disk_mb=lambda wildcards, attempt: attempt * imputation_dict["vcf_fixref_hg38_memory"]
    threads:
        imputation_dict["vcf_fixref_hg38_threads"]
    shell:
        """
        bcftools +fixref {input.data_vcf} -- -f {input.fasta} -i {input.vcf} | \
        bcftools norm --check-ref x -f {input.fasta} -Oz -o {output.vcf}

        #Index
        bcftools index {output.vcf}
        """


rule filter_preimpute_vcf:
    input:
        vcf = "results/vcf_fixref_hg38/{ancestry}_fixref_hg38.vcf.gz"
    output:
        tagged_vcf = "results/filter_preimpute_vcf/{ancestry}_tagged.vcf.gz",
        filtered_vcf = "results/filter_preimpute_vcf/{ancestry}_filtered.vcf.gz",
        filtered_index = "results/filter_preimpute_vcf/{ancestry}_filtered.vcf.gz.csi"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * imputation_dict["filter_preimpute_vcf_memory"],
        disk_mb=lambda wildcards, attempt: attempt * imputation_dict["filter_preimpute_vcf_memory"]
    threads:
        imputation_dict["filter_preimpute_vcf_threads"]
    params:
        maf = lambda wildcards: float(maf_df["MAF"][maf_df.Ancestry == wildcards.ancestry].values),
        missing = imputation_dict["snp_missing_pct"],
        hwe = imputation_dict["snp_hwe"]
    shell:
        """
        #Add tags
        bcftools +fill-tags {input.vcf} -Oz -o {output.tagged_vcf}

        #Filter rare and non-HWE variants and those with abnormal alleles and duplicates
        bcftools filter -i 'INFO/HWE > {params.hwe} & F_MISSING < {params.missing} & MAF[0] > {params.maf}' {output.tagged_vcf} |\
        bcftools filter -e 'REF="N" | REF="I" | REF="D"' |\
        bcftools filter -e "ALT='.'" |\
        bcftools norm -d all |\
        bcftools norm -m+any |\
        bcftools view -m2 -M2 -Oz -o {output.filtered_vcf}

        #Index the output file
        bcftools index {output.filtered_vcf}
        """

rule het:
    input:
        vcf = "results/filter_preimpute_vcf/{ancestry}_filtered.vcf.gz",
    output:
        tmp_vcf = temp("results/het/{ancestry}_filtered_temp.vcf"),
        inds = "results/het/{ancestry}_het_failed.inds",
        het = "results/het/{ancestry}_het.het",
        passed = "results/het/{ancestry}_het_passed.inds",
        passed_list = "results/het/{ancestry}_het_passed.txt"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * imputation_dict["het_memory"],
        disk_mb=lambda wildcards, attempt: attempt * imputation_dict["het_memory"]
    threads: imputation_dict["het_threads"]
    params:
        het_base = workflow.default_remote_prefix + "/results/het/{ancestry}_het",
        script = "/opt/WG1-pipeline-QC/Imputation/scripts/filter_het.R",
        hwe = workflow.default_remote_prefix + "/results/hwe/{ancestry}_hwe",
        out = workflow.default_remote_prefix + "/results/het/{ancestry}_het"
    shell:
        """
        gunzip -c {input.vcf} | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' > {output.tmp_vcf}
        vcftools --vcf {output.tmp_vcf} --het --out {params.het_base}
        Rscript {params.script} {output.het} {output.inds} {output.passed} {output.passed_list}
        """

rule het_filter:
    input:
        passed_list = "results/het/{ancestry}_het_passed.txt",
        vcf = "results/filter_preimpute_vcf/{ancestry}_filtered.vcf.gz"
    output:
        vcf = "results/het_filter/{ancestry}_het_filtered.vcf.gz",
        index = "results/het_filter/{ancestry}_het_filtered.vcf.gz.csi"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * imputation_dict["het_filter_memory"],
        disk_mb=lambda wildcards, attempt: attempt * imputation_dict["het_filter_memory"]
    threads: imputation_dict["het_filter_threads"]
    params:
        hwe = workflow.default_remote_prefix + "/results/hwe/{ancestry}_hwe",
        out = workflow.default_remote_prefix + "/results/het_filter/{ancestry}_het_filter",
    shell:
        """
        bcftools view -S {input.passed_list} {input.vcf} -Oz -o {output.vcf}

        #Index the output file
        bcftools index {output.vcf}
        """


rule calculate_missingness:
    input:
        filtered_vcf = "results/het_filter/{ancestry}_het_filtered.vcf.gz",
        filtered_index = "results/het_filter/{ancestry}_het_filtered.vcf.gz.csi"
    output:
        tmp_vcf = temp("results/filter_preimpute_vcf/{ancestry}_het_filtered.vcf"),
        miss = "results/filter_preimpute_vcf/{ancestry}_genotypes.imiss",
        individuals = "results/genotype_donor_annotation/{ancestry}_individuals.tsv"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * imputation_dict["calculate_missingness_memory"],
        disk_mb=lambda wildcards, attempt: attempt * imputation_dict["calculate_missingness_memory"]
    threads:
        imputation_dict["calculate_missingness_threads"]
    params:
        out = workflow.default_remote_prefix + "/results/filter_preimpute_vcf/{ancestry}_genotypes"
    shell:
        """
        gunzip -c {input.filtered_vcf} | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' > {output.tmp_vcf}

        vcftools --gzvcf {output.tmp_vcf} --missing-indv --out {params.out}

        bcftools query -l {input.filtered_vcf} >> {output.individuals}
        """


rule split_by_chr:
    input:
        filtered_vcf = "results/het_filter/{ancestry}_het_filtered.vcf.gz",
        filtered_index = "results/het_filter/{ancestry}_het_filtered.vcf.gz.csi"
    output:
        vcf = "results/split_by_chr/{ancestry}_chr_{chr}.vcf.gz",
        index = "results/split_by_chr/{ancestry}_chr_{chr}.vcf.gz.csi"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * imputation_dict["split_by_chr_memory"],
        disk_mb=lambda wildcards, attempt: attempt * imputation_dict["split_by_chr_memory"]
    threads:
        imputation_dict["split_by_chr_threads"]
    shell:
        """
        bcftools view -r {wildcards.chr} {input.filtered_vcf} -Oz -o {output.vcf}
        bcftools index {output.vcf}
        """


rule eagle_prephasing:
    input:
        vcf = "results/split_by_chr/{ancestry}_chr_{chr}.vcf.gz",
        vcf_index = "results/split_by_chr/{ancestry}_chr_{chr}.vcf.gz.csi",
        map_file = genetic_map,
        phasing_file = phasing_dir + "/chr{chr}.bcf",
        phasing_index = phasing_dir + "/chr{chr}.bcf.csi"
    output:
        vcf = "results/eagle_prephasing/{ancestry}_chr{chr}_phased.vcf.gz"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * imputation_dict["eagle_prephasing_memory"],
        disk_mb=lambda wildcards, attempt: attempt * imputation_dict["eagle_prephasing_memory"]
    threads: imputation_dict["eagle_prephasing_threads"]
    params:
        out = workflow.default_remote_prefix + "/results/eagle_prephasing/{ancestry}_chr{chr}_phased"
    shell:
        """
        eagle --vcfTarget={input.vcf} \
            --vcfRef={input.phasing_file} \
            --geneticMapFile={input.map_file} \
            --chrom={wildcards.chr} \
            --outPrefix={params.out} \
            --numThreads={threads}
        """


rule minimac_imputation:
    input:
        vcf = "results/eagle_prephasing/{ancestry}_chr{chr}_phased.vcf.gz",
        impute_file = impute_dir + "/chr{chr}.m3vcf.gz"
    output:
        "results/minimac_imputed/{ancestry}_chr{chr}.dose.vcf.gz"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * imputation_dict["minimac_imputation_memory"],
        disk_mb=lambda wildcards, attempt: attempt * imputation_dict["minimac_imputation_memory"]
    threads:
        imputation_dict["minimac_imputation_threads"]
    params:
        out = workflow.default_remote_prefix + "/results/minimac_imputed/{ancestry}_chr{chr}",
        minimac4 = "/opt/bin/minimac4",
        chunk_length = imputation_dict["chunk_length"]
    shell:
        """
        {params.minimac4} --refHaps {input.impute_file} \
            --haps {input.vcf} \
            --prefix {params.out} \
            --format GT,DS,GP \
            --noPhoneHome \
            --cpus {threads} \
            --ChunkLengthMb {params.chunk_length}
        """


rule combine_vcfs_ancestry:
    input:
        vcfs = lambda wildcards: expand("results/minimac_imputed/{ancestry}_chr{chr}.dose.vcf.gz", chr = chromosomes, ancestry = ancestry_subsets)
    output:
        combined = "results/vcf_merged_by_ancestries/{ancestry}_imputed_hg38.vcf.gz",
        ind = "results/vcf_merged_by_ancestries/{ancestry}_imputed_hg38.vcf.gz.csi"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * imputation_dict["combine_vcfs_ancestry_memory"],
        disk_mb=lambda wildcards, attempt: attempt * imputation_dict["combine_vcfs_ancestry_memory"]
    threads: imputation_dict["combine_vcfs_ancestry_threads"]
    params:
        files_begin = workflow.default_remote_prefix + "/results/minimac_imputed/{ancestry}_chr*.dose.vcf.gz"
    shell:
        """
        bcftools concat -Oz {params.files_begin} > {output.combined}
        bcftools index {output.combined}
        """


rule combine_vcfs_all:
    input:
        vcfs = lambda wildcards: expand("results/vcf_merged_by_ancestries/{ancestry}_imputed_hg38.vcf.gz", ancestry = ancestry_subsets)
    output:
        combined = "results/vcf_all_merged/imputed_hg38.vcf.gz",
        ind = "results/vcf_all_merged/imputed_hg38.vcf.gz.csi"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * imputation_dict["combine_vcfs_memory"],
        disk_mb=lambda wildcards, attempt: attempt * imputation_dict["combine_vcfs_memory"]
    threads: imputation_dict["combine_vcfs_threads"]
    shell:
        """
        if [[ $(ls -l {input.vcfs} | wc -l) > 1 ]]
        then
            bcftools merge -Oz {input.vcfs} > {output.combined}
        else
            cp {input.vcfs} {output.combined}
        fi
        bcftools index {output.combined}
        """

rule filter4demultiplexing:
    input:
        "results/vcf_all_merged/imputed_hg38.vcf.gz"
    output:
        info_filled = "results/vcf_all_merged/imputed_hg38_info_filled.vcf.gz",
        qc_filtered = "results/vcf_all_merged/imputed_hg38_R2_0.3_MAF0.05.vcf.gz",
        location_filtered = temp("results/vcf_all_merged/imputed_hg38_R2_0.3_MAF0.05_exons.recode.vcf"),
        complete_cases = "results/vcf_all_merged/imputed_hg38_R2_0.3_MAF0.05_exons_complete_cases.recode.vcf"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * imputation_dict["filter4demultiplexing_memory"],
        disk_mb=lambda wildcards, attempt: attempt * imputation_dict["filter4demultiplexing_memory"]
    threads: imputation_dict["filter4demultiplexing_threads"]
    params:
        out = workflow.default_remote_prefix + "/results/vcf_all_merged/imputed_hg38_R2_0.3_MAF0.05_exons",
        complete_out = workflow.default_remote_prefix + "/results/vcf_all_merged/imputed_hg38_R2_0.3_MAF0.05_exons_complete_cases",
        bed = "/opt/hg38exonsUCSC.bed"
    shell:
        """
        ##### Add all the info fields
        bcftools +fill-tags -Oz --output {output.info_filled} {input}

        ##### Filter the Imputed SNP Genotype by Minor Allele Frequency (MAF) and INFO scores #####
        bcftools filter --include 'MAF>=0.05 & R2>=0.3' -Oz --output {output.qc_filtered} {output.info_filled}

        vcftools \
            --gzvcf {output.qc_filtered} \
            --max-alleles 2 \
            --remove-indels \
            --bed {params.bed} \
            --recode \
            --recode-INFO-all \
            --out {params.out}

        vcftools --recode --recode-INFO-all --vcf {output.location_filtered} --max-missing 1 --out {params.complete_out}
        """

rule sort4demultiplexing:
    input:
        complete_cases = "results/vcf_all_merged/imputed_hg38_R2_0.3_MAF0.05_exons_complete_cases.recode.vcf"
    output:
        complete_cases_sorted = "results/vcf_4_demultiplex/imputed_hg38_R2_0.3_MAF0.05_exons_sorted.vcf"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * imputation_dict["sort4demultiplexing_memory"],
        java_mem = lambda wildcards, attempt: attempt * imputation_dict["sort4demultiplexing_java_memory"],
        disk_mb=lambda wildcards, attempt: attempt * imputation_dict["sort4demultiplexing_memory"]
    threads: imputation_dict["sort4demultiplexing_threads"]
    shell:
        """
        java -Xmx{resources.java_mem}M -Xms{resources.java_mem}M -jar /opt/picard/build/libs/picard.jar SortVcf \
            -I {input.complete_cases} \
            -O {output.complete_cases_sorted}
        """


rule count_snps:
    input:
        info_filled = "results/vcf_all_merged/imputed_hg38_info_filled.vcf.gz",
        qc_filtered = "results/vcf_all_merged/imputed_hg38_R2_0.3_MAF0.05.vcf.gz",
        complete_cases_sorted = "results/vcf_4_demultiplex/imputed_hg38_R2_0.3_MAF0.05_exons_sorted.vcf"
    output:
        report("results/metrics/Number_SNPs.png", category = "SNP Numbers", caption = "../report_captions/counts_snps.rst")
    resources:
        mem_mb=lambda wildcards, attempt: attempt * imputation_dict["count_snps_memory"],
        disk_mb=lambda wildcards, attempt: attempt * imputation_dict["count_snps_memory"]
    threads: imputation_dict["count_snps_threads"]
    params:
        basedir = workflow.default_remote_prefix + "/results",
        outdir = workflow.default_remote_prefix + "/results/metrics/",
        script = "/opt/WG1-pipeline-QC/Imputation/scripts/SNP_numbers.R"
    shell:
        """
        Rscript {params.script} {params.basedir} {params.outdir}
        """



rule genotype_donor_annotation:
    input:
        individuals = lambda wildcards: expand("results/genotype_donor_annotation/{ancestry}_individuals.tsv", ancestry = ancestry_subsets),
        updated_psam = "results/update_sex_ancestry/update_sex.psam"
    output:
        out_temp = temp("results/genotype_donor_annotation/genotype_donor_annotation_temp.tsv"),
        combined_individuals = "results/genotype_donor_annotation/combined_individuals.tsv",
        final = "results/genotype_donor_annotation/genotype_donor_annotation.tsv",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * imputation_dict["genotype_donor_annotation_memory"],
        disk_mb=lambda wildcards, attempt: attempt * imputation_dict["genotype_donor_annotation_memory"]
    threads: imputation_dict["genotype_donor_annotation_threads"]
    shell:
        """
        cut -f2,5- {input.updated_psam} | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' > {output.out_temp}
        cat {input.individuals} >> {output.combined_individuals}
        sed -i '1 i\IID' {output.combined_individuals}
        awk -F"\t" 'NR==FNR {{a[$1]; next}} $1 in a' {output.combined_individuals} {output.out_temp} | awk 'BEGIN{{FS=OFS="\t"}}{{sub("1","M",$2);print}}' | awk 'BEGIN{{FS=OFS="\t"}}{{sub("2","F",$2);print}}' > {output.final}
        sed -i 's/^IID\tSEX\tProvided_Ancestry/donor_id\tsex\tethnicity_super_population/g' {output.final}
        sed -i 's/$/\tsceQTL-Gen_hg38_imputation_pipeline\t1000g_30x_GRCh38_ref/' {output.final}
        sed -i '1s/SangerImputationServer/imputation_method/g' {output.final}
        sed -i '1s/1000g_30x_GRCh38_ref/imputation_reference/g' {output.final}
        """