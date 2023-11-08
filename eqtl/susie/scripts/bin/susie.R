#!/usr/bin/env Rscript

library(locuscomparer)
library(glue)
library(susieR)
suppressPackageStartupMessages(library(eQTLUtils))
suppressPackageStartupMessages(library(dplyr))
library(optparse)

option_list <- list(
    make_option(c("--phenotype"), type = "character", help = "[Required] QTLtools format phenotype bed"),
    make_option(c("--gene_name"), type = "character", help = "File with list of gene names to select using the gid field of the phenotype bed."),
    make_option(c("--dosage"), type = "character", help = "[Required] dosage file"),
    make_option(c("--min_abs_cor"), type = "numeric", help = "[Required] Min abs correlation b/w any pair of variants in the set"),
    make_option(c("--window"), type = "numeric", help = "[Required] window over phenotype end to consider variants "),
    make_option(c("--prefix"), type = "character", help = "output prefix")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)


finemap = function(b, opts, gene_name){
    prefix = glue("{opts$prefix}.{gene_name}")
    print(prefix)
    
    Yi = b[b$gid == gene_name,]
    chrom = Yi[['#chr']]
    window = opts$window
    start = Yi[['start']] - window
    end = Yi[['end']] + window
    Yi[,c("#chr", "start", "end", "id", "gid", "strd")] = NULL
    samples = colnames(Yi)

    ## actual dosages from DS
    genotype_matrix_ds = eQTLUtils::extractGenotypeMatrixFromDosage(
        chr = chrom, 
        start = start, 
        end = end, 
        dosage_file = opts$dosage)
    
    genotype_matrix_ds = t(genotype_matrix_ds[,samples])
    
    out_ds = susie(genotype_matrix_ds, t(Yi), L = 10,
                   estimate_residual_variance = TRUE, 
                   estimate_prior_variance = TRUE,
                   verbose = TRUE,
                   min_abs_cor = opts$min_abs_cor)
    print(out_ds$sets)

    save(out_ds, file=glue("{prefix}.susie.Rda"))

    ## save:
    ## pip dataframe
    pipv = susie_get_pip(out_ds)
    pdf = data.frame(list("pip" = pipv, "var" = colnames(genotype_matrix_ds)))
    pdf$var = gsub("chrchr", "chr", pdf$var)
    
    write.table(pdf, file=glue("{prefix}.susie-pip.tsv"), sep='\t', quote=F)
    
    ## pdf = pdf[order(pdf[,c("pip")], decreasing = TRUE),]
    
    ## write.table(pdf, file=glue("{prefix}.susie-pip.tsv"), sep='\t', row.names=F, quote=F)
    
    ## how many 95% csets?
    n95 = length(out_ds$sets$cs)
    
    ## best set
    best = row.names(out_ds$sets$purity[which.max(out_ds$sets$purity$min.abs.corr),])
    
    ## size of the main set
    len_1_95 = length(out_ds$sets$cs[[best]])
    
    ## save the main 95% cset
    write.table(pdf[out_ds$sets$cs[[best]],], file=glue("{prefix}.susie-cset95.tsv"), sep='\t', quote=F)

    print(glue("{n95} 95% credible sets found. Length of the first 95% credible set = {len_1_95}"))
    ## print(glue("{n99} 99% credible sets found. Length of the first 99% credible set = {len_1_99}"))


    png(glue("{prefix}.susie.png"), height=2.5, width=4, units="in", res=150)
    susie_plot(out_ds, y="PIP")
    dev.off()
}


phenotype = read.csv(opts$phenotype, sep='\t', header=T, as.is=T, check.names = F)

genenamelist = read.csv(opts$gene_name, sep='\t', header=F, as.is=T)$V1

lapply(genenamelist, function(i){finemap(phenotype, opts, i)})

## for (i in genenamelist){
##     finemap(phenotype, opts$prefix, i)
## }
