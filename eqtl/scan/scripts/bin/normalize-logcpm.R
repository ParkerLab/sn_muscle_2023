#!/usr/bin/env Rscript

args = commandArgs(TRUE)
library(tidyr)
library(dplyr)
library(ggplot2)


                                        # read data
d = read.csv(args[1], sep='\t', header=T)
                                        # Get TPM and inverse normalize across samples
infocols = c('X.Chr', 'start', 'end', 'pid', 'gid', 'strand', 'length')

d1 = d[, ! names(d) %in% infocols]
d1 = d1[,apply(d1,2,function(x) !all(x==0))] 

d1= d1 %>% mutate_if(is.numeric, funs(./(sum(.)/1000000)))
head(d1)

                                        # get top 10k median cpm genes
d1$median_cpm = apply(d1, 1, median, na.rm = T)
dt = cbind(d[, infocols], d1)
dt =  dt[order(-dt$median_cpm),]
write.table(dt, paste(args[3], ".t.txt", sep=""), quote=FALSE, row.names=FALSE)
dt = head(dt, 10000)
write.table(dt$pid, args[3], col.names=FALSE, quote=FALSE, row.names=FALSE)
d1$median_cpm = NULL

                                        # log trasnform counts
d1 = log10(d1 + 1)
do = cbind(d[, infocols], d1)
write.table(do, args[2], sep='\t', quote=FALSE, row.names=FALSE)
