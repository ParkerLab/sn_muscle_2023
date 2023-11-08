#!/usr/bin/env Rscript

library(glue)
library(cit)
library(bedr)
library(optparse)
library(dplyr)


option_list <- list(
    make_option(c("--cit"), type = "character", help = "[Required] comma separated rds files for cit runs"),
    make_option(c("--prefix"), type = "character", help = "[Required] output prefix"),
    make_option(c("--citdir"), type = "character", help = "[Required] dataframe to concat results to")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

prefix = opts$prefix
df = data.frame()
l1 = c()
l2 = c()
l3 = c()
l4 = c()

for (f in unlist(strsplit(opts$cit, ","))){
  print(f)
  result = readRDS(f)
  pref = gsub(".citres.rds", ".cit.tsv", basename(f))
  dfname = glue("{opts$citdir}/{pref}")
  df = rbind(df, read.csv(dfname, sep='\t', header=T))
  l1 = c(l1, lapply(seq(1, length(result)), function(i){result[[i]][[1]]}))
  l2 = c(l2, lapply(seq(1, length(result)), function(i){result[[i]][[2]]}))
  l3 = c(l3, lapply(seq(1, length(result)), function(i){result[[i]][[3]]}))
  l4 = c(l4, lapply(seq(1, length(result)), function(i){result[[i]][[4]]}))
  
}

dl1 = cbind(df, fdr.cit(l1))
dl2 = cbind(df, fdr.cit(l2))
dl3 = cbind(df, fdr.cit(l3))
dl4 = cbind(df, fdr.cit(l4))

write.table(dl1, file=glue("{prefix}.causal_ca_eq_hit.tsv"), sep='\t', quote=F, row.names=F)
write.table(dl2, file=glue("{prefix}.revcausal_ca_eq_hit.tsv"), sep='\t', quote=F, row.names=F)
write.table(dl3, file=glue("{prefix}.causal_eq_ca_hit.tsv"), sep='\t', quote=F, row.names=F)
write.table(dl4, file=glue("{prefix}.revcausal_eq_ca_hit.tsv"), sep='\t', quote=F, row.names=F)
