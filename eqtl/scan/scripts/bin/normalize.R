#!/usr/bin/env Rscript

args = commandArgs(TRUE)
library(tidyr)
library(dplyr)
library(ggplot2)
library(RNOmni)
library(glue)

counts = args[1]
prefix = args[2]
## samplelist = args[3]
## read data

## slist = read.csv(samplelist, sep='\t', header=F, as.is=T)$V1
d = read.csv(counts, sep='\t', header=T, as.is=T, check.names=F)
d = d %>% rename( "X.Chr" = "#Chr")
d$length = d$length/1000
# get gene TSS from gene coordinates
d$start = ifelse(d$strand == "+", d$start, d$end-1)
d$end = d$start + 1
d = d[with(d, order(X.Chr, start, end)), ]
                                        # Get TPM and inverse normalize across samples
infocols = c('X.Chr', 'start', 'end', 'pid', 'gid', 'strand', 'length')
## keep.columns = c(infocols, slist)
## print(keep.columns)
## d = d[, keep.columns]

d1 = d[, ! names(d) %in% infocols]
d1 = d1[,apply(d1,2,function(x) !all(x==0))] 
d1 = d1/d$length

d1= d1 %>% mutate_if(is.numeric, funs(./(sum(.)/1e6)))
head(d1)

                                        # get top 10k median tpm genes
d1$median_tpm = apply(d1, 1, median, na.rm = T)
dt = cbind(d[, infocols], d1)
dt =  dt[order(-dt$median_tpm),]
write.table(dt, glue("{prefix}.median_tpm.txt"), quote=FALSE, row.names=FALSE)
dt = head(dt, 10000)
write.table(dt$pid, glue("{prefix}.top10kgenes.txt"), col.names=FALSE, quote=FALSE, row.names=FALSE)
d1$median_tpm = NULL

                                        # Inverse normalize across samples
d1 = t(apply(d1, 1, function(row){rankNorm(row)} ) )
do = cbind(d[, infocols], d1)
do$length = NULL
write.table(do, glue("{prefix}.bed"), sep='\t', quote=FALSE, row.names=FALSE)
