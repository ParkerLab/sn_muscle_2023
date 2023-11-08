#!/usr/bin/env Rscript

library(locuscomparer)
library(glue)
library(susieR)
suppressPackageStartupMessages(library(eQTLUtils))
suppressPackageStartupMessages(library(tidyr))
library(optparse)
library(coloc)

option_list <- list(
    make_option(c("--eqtl"), type = "character", help = "[Required] eqtl SuSiE output to run coloc with"),
    make_option(c("--caqtl"), type = "character", help = "[Required] caqtl SuSiE outputs to run coloc with"),
    make_option(c("--prefix"), type = "character", help = "[Required] prefix")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

loadRData <- function(fileName){
    ##loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}

snpincset = function(object, snp, index){
    cset.snp.index = object$sets$cs[[glue("L{index}")]]
    cset.snps = colnames(object$lbf_variable)[cset.snp.index]
      
    if (snp %in% cset.snps){
        return("yes")
    } else {
        return("no")
    }
}

run.coloc = function(eq, cafile, prefix){
    ca = loadRData(cafile)
    if(length(ca$sets$cs)){
        res = coloc.susie(eq, ca)
        ## check if the hits found by coloc are present in the orignial credible set. New hits could be assigned by SuSiE if not all SNPs overlap between the two signals
        res$summary$hit1.in.eqtl = unlist(apply(res$summary, 1, function(i){snpincset(eq, i['hit1'], i['idx1'])}))
        res$summary$hit2.in.caqtl = unlist(apply(res$summary, 1, function(i){snpincset(ca, i['hit2'], i['idx2'])}))

        write.table(res$summary, file = glue("{prefix}.coloc.tsv"), sep='\t', quote=F)
    }
}

## load eQTL SuSiE object
eq = loadRData(opts$eqtl) ## loads eq

ca.file.list = unlist(strsplit(opts$caqtl, ","))

if(length(eq$sets$cs)){
    lapply(ca.file.list, function(cafile){
        caqtl.locus = unlist(strsplit(basename(cafile), "\\."))[6]
        prefix = glue("{opts$prefix}--{caqtl.locus}")
        run.coloc(eq, cafile, prefix)

    })
}
