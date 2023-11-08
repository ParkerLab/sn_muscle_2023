#!/usr/bin/env Rscript

args = commandArgs(TRUE)
library(tidyr)
library(dplyr)
library(ggplot2)
library(RNOmni)
library(optparse)
library(glue)

option_list <- list(
    make_option(c("--peak_counts"), type = "character", help = "[Required] tsv file with peak fragment counts by sample."),
    make_option(c("--summit_list"), type = "character", help = "[Required] List of selected summits"),
    make_option(c("--sample_list"), type = "character", help = "[Required] headerless, one colum file with the subset of samples to retain"),
    make_option(c("--output_bed"), type = "character", help = "output formatted bed file")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

## read data
d = read.csv(opts$peak_counts, sep='\t', header=T, as.is=T, row.names=1, check.names=F)
samples_in = colnames(d)
samples_select = read.csv(opts$sample_list, sep='\t', header=F, as.is=T)$V1
samples_filter = intersect(samples_in, samples_select)

summits_in = rownames(d)
summits_select = read.csv(opts$summit_list, sep='\t', header=F, as.is=T)$V1
summits_filter = intersect(summits_in, summits_select)

print(glue("Samples being filtered = {length(samples_filter)}; Summits being filtered = {length(summits_select)}"))

d = d[summits_filter, samples_filter]

# there should be no samples with 0 counts but check anyways
d = d[,apply(d, 2, function(x) !all(x==0))]

peaks = data.frame(rownames(d))
peaks$pid = peaks$rownames.d.
peaks$gid = peaks$pid
peaks = peaks %>% separate(rownames.d., c("chrom", "start", "end"), ":")
peaks$start = as.numeric(peaks$start)
peaks$end = as.numeric(peaks$end)
peaks$length = (peaks$end - peaks$start)/1000
peaks$strand = "."
## get peak midpoint from coordinates

peaks$start = peaks$start + round((peaks$end - peaks$start)/2, 0)
peaks$end = peaks$start + 1

## Get TPM and inverse normalize across samples
infocols = c('chrom', 'start', 'end', 'pid', 'gid', 'strand')

d = d/peaks$length

d = d %>% mutate_if(is.numeric, funs(./(sum(.)/1e6)))

                                        # Inverse normalize across samples
d = t(apply(d, 1, function(row){rankNorm(row)} ) )
do = cbind(peaks[, infocols], d)

do = do[with(do, order(chrom, start, end)), ]
colnames(do) = gsub("--.*", "", colnames(do))
options(scipen=10)
write.table(do, opts$output_bed, sep='\t', quote=FALSE, row.names=FALSE)
