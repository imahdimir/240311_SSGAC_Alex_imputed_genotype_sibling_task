#!/usr/bin/bash

pgen_fp_pat=/disk/genetics/ukb/mahdimir/imputed_genotype_corr_Tammys_analysis_replication/plink_out_Tammys_way/bgen_*.pgen
out_dir=/disk/genetics/ukb/mahdimir/imputed_genotype_corr_Tammys_analysis_replication/plink_out_Tammys_way/bed

for pfile in $pgen_fp_pat
do 

    prefix=${pfile%.*}
    out=$(basename $prefix)

    /disk/genetics/tools/plink2/plink2 --pfile ${prefix} --make-bed --out "$out_dir"/${out}

done