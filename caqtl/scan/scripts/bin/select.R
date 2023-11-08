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
    make_option(c("--sample_list"), type = "character", help = "[Required] one column file with sample names to select."),
    make_option(c("--min_count"), type = "numeric", help = "[Required] Min count in a sample for a peak to consider that sample"),
    make_option(c("--min_samples"), type = "numeric", help = "[Required] Min number of sample satisfying the min_count threshold for a peak to be considered"),
    make_option(c("--output_bed"), type = "character", help = "output formatted bed file")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

## read data
d = read.csv(opts$peak_counts, sep='\t', header=T, as.is=T, row.names=1, check.names=F)

## selected sample:
samples = read.csv(opts$sample_list, sep='\t', header=T, as.is=T)$sample
d = d[, samples]


d = d[,apply(d, 2, function(x) !all(x==0))]

# select peaks with minct and minsample 

d1 = data.frame(apply(d, 2, function(x){ifelse(x < opts$min_count, 0, x)}))
row.names(d1) = row.names(d)
d1$nonzero = apply(d1, 1, function(c)sum(c!=0))

selected_peaks = row.names(d1[d1$nonzero >= opts$min_samples,])
d = d[selected_peaks, ]

peaks = data.frame(rownames(d))
write.table(peaks, opts$output_bed, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
