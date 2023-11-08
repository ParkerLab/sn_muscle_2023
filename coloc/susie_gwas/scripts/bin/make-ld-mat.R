#!/usr/bin/env Rscript

library(susieR)
library(glue)
suppressPackageStartupMessages(library(dplyr))
library(optparse)


process_dosage = function(f, snplist){
    ld = read.csv(f, sep='\t', check.names = F)
    dups = ld[ (duplicated(ld$ID) | duplicated(ld$ID, fromLast = TRUE)),]
    print(glue("N duplicates = {nrow(dups)}"))
    ld = ld[! (duplicated(ld$ID) | duplicated(ld$ID, fromLast = TRUE)),]
    row.names(ld) = ld$ID
    ld$ID = NULL
    idlist = intersect(snplist, row.names(ld))
    ld = ld[idlist,]
    print(ld[1:5, 1:10])
    ld = cor(t(ld))
    return(ld)
}

option_list <- list(
    make_option(c("--trait1"), type = "character", help = "[Required] trait1 dataframe"),
    make_option(c("--trait1_ld"), type = "character", help = "[Required] genotype dosage for trait1 reference panel"),
    make_option(c("--prefix"), type = "character", help = "output prefix")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

## read vars
prefix = opts$prefix

## runsusie GWAS only technically needs beta and se since it is case control, but for calculating the s metric, we also need Z scores. i.e using p values column.
## read trait1 
input1 = read.csv(opts$trait1, sep='\t', header=T, as.is=T, check.names=F)
snps = input1$snp
print(head(input1))

t1.ld = process_dosage(opts$trait1_ld, snps)
print("dosage done")
## LD mat
write.table(t1.ld, sep='\t', quote=F, file=glue("{prefix}.ld.tsv"))
