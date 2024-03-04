#!/usr/bin/bash

bgen_fp=/disk/genetics/data/ukb/private/v3/raw/imputed/ukb_imp_chr22_v3.bgen
sample_fp=/disk/genetics2/ukb/orig/UKBv3/sample/ukb11425_imp_chr1_22_v3_s487395.sample
sibs_fp=/disk/genetics/ukb/mahdimir/imputed_genotype_corr_Tammys_analysis_replication/sibs.txt
snps_30_fp=/disk/genetics/ukb/mahdimir/imputed_genotype_corr_Tammys_analysis_replication/snps_2_keep/snps_30.txt
out_dir=/disk/genetics/ukb/mahdimir/imputed_genotype_corr_Tammys_analysis_replication/plink_out_3

/disk/genetics/tools/plink2/plink2 --bgen $bgen_fp ref-first --export bgen-1.2 --sample $sample_fp --keep $sibs_fp --extract $snps_30_fp --out "$out_dir"/bgen_30