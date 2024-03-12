#!/usr/bin/zsh

bgen_fp=/disk/genetics/data/ukb/private/v3/raw/imputed/ukb_imp_chr22_v3.bgen
sample_fp=/disk/genetics2/ukb/orig/UKBv3/sample/ukb11425_imp_chr1_22_v3_s487395.sample
po_fp=/var/genetics/ws/mahdimir/DropBox/2-GPro/1-med-GPro/imputed-genotype-sibling-task-240311/parent_offspring.txt
snps_30_fp=/var/genetics/ws/mahdimir/DropBox/2-GPro/1-med-GPro/imputed-genotype-sibling-task-240311/snps_30.txt
snps_99_fp=/var/genetics/ws/mahdimir/DropBox/2-GPro/1-med-GPro/imputed-genotype-sibling-task-240311/snps_99.txt
out_dir=/var/genetics/ws/mahdimir/DropBox/2-GPro/1-med-GPro/imputed-genotype-sibling-task-240311/plink_out

plink2 --bgen $bgen_fp ref-first --export bgen-1.2 --sample $sample_fp --keep $po_fp --extract $snps_30_fp --out "$out_dir"/snps_30_po
plink2 --bgen $bgen_fp ref-first --export bgen-1.2 --sample $sample_fp --keep $po_fp --extract $snps_99_fp --out "$out_dir"/snps_99_po
