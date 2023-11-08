#!/usr/bin/env Rscript

library(qvalue)
library(tidyr)
library(optparse)

option_list <- list(
    make_option(c("--df"), type = "character", help = "[Required] tsv file with p value column."),
    make_option(c("--p"), type = "character", help = "[Required] p value column name "),
    make_option(c("--select"), type = "character", help = "[Required] file with selected features"),
    make_option(c("--q"), type = "character", help = "[Required] intended qvalue column name. [Default = qvalue]"),
    make_option(c("--output"), type = "character", help = "output file name. Rows with NAs in the p value column will not be printed.")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

s = read.csv(opts$select, sep='\t', header=F, as.is=T)$V1
## read data
d = read.csv(opts$df, sep='\t', header=T, as.is=T)
d = d[d$gene_name %in% s,]

d = d %>% drop_na(opts$p)
d[[opts$q]] = qvalue(d[[opts$p]])$qvalue
write.table(d, opts$output, sep='\t', quote=FALSE, row.names=FALSE)
