#!/usr/bin/env Rscript

library(locuscomparer)
library(glue)
library(susieR)
suppressPackageStartupMessages(library(eQTLUtils))
suppressPackageStartupMessages(library(tidyr))
library(optparse)
library(coloc)

option_list <- list(
    make_option(c("--gwas"), type = "character", help = "[Required] GWAS SuSiE output"),
    make_option(c("--caqtl"), type = "character", help = "space-separated caQTL SuSiE outputs to run coloc with"),
    make_option(c("--ukbb"), type = "character", help = "[Required] ukbb dosage file to match rsids")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

matchval = function(df, colname, value, ret.colname){
  out = value
  sub = df[df[[colname]] == value,]
  if (nrow(sub)==1){
    out = sub[[ret.colname]]
  }
  return(out)
}

snpincset = function(object, snp, index){
    if (snp %in% names(object$sets$cs[[glue("L{index}")]])){
        return("yes")
    } else {
        return("no")
    }
}

run.coloc = function(gwas.res, caqtl, rsdf, prefix){

    load(caqtl) ## loads out_ds
    ## caQTL snp names are not rsids, fix that.
    snp = colnames(out_ds$lbf_variable)
    snp = gsub("chrchr", "chr", snp)
    d = data.frame(list("index" = seq(1:length(snp)), "snp" = snp))
    d$snp1 = d$snp
    d = tidyr::separate(data = d, col = snp1, into = c("chrom", "pos", "ref_f", "alt_f"), sep = "_")
    
    dx = merge(d, rsdf, on=c("chrom", "pos"))
    dx = dx[((dx$ref_f == dx$ref_u) & (dx$alt_f == dx$alt_u)) | ((dx$ref_f == dx$alt_u) & (dx$alt_f == dx$ref_u)),]
    print(head(dx))
    colnames(out_ds$lbf_variable) = unlist(lapply(gsub("chrchr", "chr", colnames(out_ds$lbf_variable)), function(i){matchval(dx, "snp", i, "id")}))
    colnames(out_ds$alpha) = unlist(lapply(gsub("chrchr", "chr", colnames(out_ds$alpha)), function(i){return(matchval(dx, "snp", i, "id"))}))
    names(out_ds$pip) = unlist(lapply(gsub("chrchr", "chr", names(out_ds$pip)), function(i){matchval(dx, "snp", i, "id")}))
    
    if(length(out_ds$sets$cs)){
        out_ds$sets$cs = lapply(out_ds$sets$cs, function(cset){
            names(cset) = unlist(lapply(cset, function(i){
                matchval(dx, "index", i, "id")
            }))
            return(cset)
        })
        print(head(colnames(out_ds$lbf_variable)))
        print(intersect(colnames(gwas.res$lbf_variable),colnames(out_ds$lbf_variable)))

        res = coloc.susie(gwas.res, out_ds)
        # record coloc sensitivity 
        for (row in 1:nrow(res$summary)){
            hit1 = res$summary[row, "hit1"]
            hit2 = res$summary[row, "hit2"]
            sensitivity.prefix = glue("{prefix}.coloc-sensitivity.hit1{hit1}-hit2{hit2}.row{row}")
            png(glue("{sensitivity.prefix}.png"), width=7, height=3, units="in", res=300)
            sens = sensitivity(res, rule="H4 > 0.5", plot.manhattans = FALSE, row=row)
            dev.off()
            write.table(sens, file = glue("{sensitivity.prefix}.tsv"), sep='\t', quote=F)
            ## get the min p12 where
            sens = sens[sens$pass,]
            if (nrow(sens) > 0 ) {
                p12min = head(sens, 1)$p12
            } else {
                p12min = 1
            }
            res$summary[row, "p12min"] = p12min
        }
        
        ## check if the hits found by coloc are present in the orignial credible set. New hits could be assigned by SuSiE if not all SNPs overlap between the two signals
        res$summary$hit1.in.gwas = unlist(apply(res$summary, 1, function(i){snpincset(gwas.res, i['hit1'], i['idx1'])}))
        res$summary$hit2.in.qtl = unlist(apply(res$summary, 1, function(i){snpincset(out_ds, i['hit2'], i['idx2'])}))
        write.table(res$summary, file = glue("{prefix}.coloc.tsv"), sep='\t', quote=F)
        
    } else {
        print(glue("{prefix} - no caQTL sets"))
    }
}

## load GWAS SuSiE object
load(opts$gwas) ## loads S1
gwas.locus = gsub(".selected.Rda", "", basename(opts$gwas))
gwas.locus = gsub(".Rda", "", gwas.locus)

dos = read.csv(opts$ukbb, sep='\t', header = F, as.is = T)
colnames(dos) = c("chrom", "pos", "id", "ref_u", "alt_u", "score", "filter", "info")
dos$id = paste(dos$id, dos$ref_u, dos$alt_u, sep="-")
ca.file.list = unlist(strsplit(opts$caqtl, ","))

lapply(ca.file.list, function(cafile){
    caqtl.locus = gsub(".susie.Rda", "", basename(cafile))
    prefix = glue("{gwas.locus}__{caqtl.locus}")
    run.coloc(S1, cafile, dos, prefix)
})
