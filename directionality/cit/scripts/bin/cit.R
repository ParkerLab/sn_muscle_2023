#!/usr/bin/env Rscript

library(glue)
library(cit)
library(bedr)
library(optparse)
library(dplyr)

option_list <- list(
    make_option(c("--cit"), type = "character", help = "[Required] dataframe listing capeak, egene, calead snp, eqtl lead snp etc to run CIT on"),
    make_option(c("--vcf"), type = "character", help = "[Required] dosage dataframe"),
    make_option(c("--eqpheno"), type = "character", help = "[Required] eqtl phenotypes"),
    make_option(c("--capheno"), type = "character", help = "[Required] caqtl phenotypes"),
    make_option(c("--covs"), type = "character", help = "[Required] covariates dataframe"),
    make_option(c("--nperm"), type = "numeric", default=0, help = "N permutations"),
    make_option(c("--seed"), type = "numeric", default=123456, help = "Seed"),
    ## make_option(c("--cahit"), type = "character", help = "caqtl hit from coloc"),
    make_option(c("--output"), type = "character", help = "[Required] output file")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

select = function(df, feature){
  df1 = df[df$gid == feature,]
  df1[, c('#Chr', 'start', 'end', 'pid', 'gid', 'strand')] = NULL
  df1 = t(df1)
  return(df1)
  
}

get_dosages = function(vcf, chrom, pos, ref=NULL, alt=NULL, id=NULL){
    start =  pos - 1
    region = paste(chrom, ":", start, "-", pos, sep="")
    print(region)
    o = bedr::tabix(region, vcf, check.zero.based=F, verbose=F)
    if (! is.null(ref)){
        o = o[o$REF == ref,]
    }
    if (! is.null(alt)){
        o = o[o$ALT == alt,]
    }
    if (! is.null(id)){
        o = o[o$ID == id,]
    }
    stopifnot(nrow(o)==1)
    
    o[, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")] = NULL
    o = o %>% mutate_all(function(i){as.numeric(unlist(strsplit(i, ":"))[2])})
    return(o)
}

greprs = function(df, rs){colnames(df)[grep(rs, colnames(df))]}

nperm = opts$nperm
seed = opts$seed
vcf = opts$vcf

d = read.csv(opts$cit, sep='\t', header=T, stringsAsFactors=F)
head(d)

covs = read.csv(opts$covs, sep='\t', header=T, check.names=F, row.names=1, as.is=F)
rownames(covs) = as.character(rownames(covs))
print(head(covs))

ca = read.csv(opts$capheno, sep='\t', header=T, check.names=F)
print("ca phenotypes read")
print(ca[1:5,1:8])

eq = read.csv(opts$eqpheno, sep='\t', header=T, check.names=F)
print("eq phenotypes read")
print(eq[1:5,1:8])

make.dos.df = function(vcf, variant.s){
    varlist = unlist(strsplit(variant.s, ","))
    doslist = lapply(varlist, function(var){
        h = unlist(strsplit(var, "_"))
        dosage = t(get_dosages(vcf, h[1], as.numeric(h[2]), ref = h[3], alt = h[4]))
        return(dosage)
    })

    dos = do.call("cbind", doslist)
    colnames(dos) = 1:ncol(dos)
    return(dos)
}

result = lapply(seq(1, nrow(d)), function(i){
    peak = d[i, 'capeak']
    gene = d[i, 'egene']
    chrom = d[i, 'chrom']
    ca1 = select(ca, peak)
    eq1 = select(eq, gene)

    ## cahit(s) -> capeak -> egene model
    dos = make.dos.df(vcf, d[i, "cahit"])
    samples = intersect(intersect(rownames(ca1), rownames(eq1)), rownames(dos))
    causal_ca_eq_hit = cit.cp(L=dos[samples,], G=ca1[samples,], T=eq1[samples,], C=covs[samples,], rseed=seed, n.perm=nperm)
    revcausal_ca_eq_hit = cit.cp(L=dos[samples,], G=eq1[samples,], T=ca1[samples,], C=covs[samples,], rseed=seed, n.perm=nperm)

    ## eqhit -> egene -> capeak model
    dos = make.dos.df(vcf, d[i, "eqhit"])
    samples = intersect(intersect(rownames(ca1), rownames(eq1)), rownames(dos))
    causal_eq_ca_hit = cit.cp(L=dos[samples,], G=eq1[samples,], T=ca1[samples,], C=covs[samples,], rseed=seed, n.perm=nperm)
    revcausal_eq_ca_hit = cit.cp(L=dos[samples,], G=ca1[samples,], T=eq1[samples,], C=covs[samples,], rseed=seed, n.perm=nperm)

    out = list(causal_ca_eq_hit, revcausal_ca_eq_hit, causal_eq_ca_hit, revcausal_eq_ca_hit)
    return(out)

})

print(result)
saveRDS(result, file=opts$output)
## d = cbind(d, t(result))
## write.table(d, file=opts$output, sep='\t', quote=F)
