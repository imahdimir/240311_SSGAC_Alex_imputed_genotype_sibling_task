#!/usr/bin/env Rscript

install.packages("data.table")
install.packages("dplyr")
install.packages("glue")
install.packages("readxl")

library("data.table")
library("dplyr")
library("plinkFile")
library("genio")
library("ggplot2")
library("glue")
library('readxl')


plot_corr <- function(pair_type, pn, gtype, gn, info, dn) {

  st <- paste('corr', pair_type, gtype, info, sep = '_')
  fp <- glue('{dta_dir}/{st}.xlsx')
  
  ost <- paste('corr_hist', pair_type, gtype, info, sep = '_')
  fpo <- glue('{fig_dir}/{ost}.png')
  
  df <- read_excel(fp, col_names = F)
  
  n_snps <- nrow(df)
  
  ggplot(data.frame(df), aes(x = ...1)) +
    geom_histogram(binwidth = 0.02, fill = "#56B4E9", colour = "#56B4E9", alpha = 0.5) +
    labs(title = paste0(pn," Corrs (",gn , ") for SNPS with INFO = ", dn, "-", dn + 0.01, " (n = ", n_snps,")"),
         x = "Genotype Correlation",
         y = "Count") +
    geom_vline(xintercept=0.5, linetype="dotted") +
    geom_vline(xintercept=mean(df$...1), linetype="dashed", color = 'red') +
    theme_classic() +
    theme(legend.position="none")
  
  ggsave(fpo)
  
}


out_dir <- '/Users/mahdi/Library/CloudStorage/Dropbox/0-all/1-out-all/imputed_genotype-sibling-task-240311'
dta_dir <- glue('{out_dir}/dta')
fig_dir <- glue('{out_dir}/fig')


pair_type <- c('sibs', 'parent_offspring')
gtype <- c('dosages', 'hard_calls')
info <- c(30, 99)

prd <- expand.grid(pair_type, gtype, info)

prd$pn <- with(prd, ifelse(Var1=='sibs', 'Sibs', 'Parent-Offspring'))
prd$gn <- with(prd, ifelse(Var2=='dosages', 'Dosages', 'Hard Calls'))
prd$dn <- with(prd, ifelse(Var3==30, .3, .99))


for (x in 1:nrow(prd)){
  pair_type = prd[x, 'Var1']
  pair_name = prd[x, 'pn']
  
  gtype = prd[x, 'Var2']
  gn = prd[x, 'gn']
  
  info = prd[x, 'Var3']
  dn = prd[x, 'dn']
  
  plot_corr(pair_type,pair_name, gtype,gn, info, dn)
}





