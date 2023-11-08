#!/usr/bin/env Rscript
library(TwoSampleMR)
library(glue)
library(tidyr)
args = commandArgs(TRUE)

mrsfile = args[1]
qtlvars = args[2]
eqtl_samples = as.numeric(args[3])
caqtl_samples = as.numeric(args[4])
prefix = args[5]


qtl = read.csv(qtlvars, sep='\t', header=F, 
               col.names=c('chr', 'SNP', 'gd', 'position', 'effect_allele', 'other_allele', 'eaf'))
head(qtl)


mrs = read.csv(mrsfile, sep='\t', header=T, as.is=T)
stopifnot(  length(unique(mrs$eqpath))==1, length(unique(mrs$capath))==1)

eqtl = mrs[1, "eqpath"]
caqtl = mrs[1, "capath"]

## eqtl
format.qtl = function(filename, qtl, samples){
    d = read.csv(filename, sep='\t', header=T)
    d = d[, c("chrom", "snp_end", "snp", "gene_name", "p_nominal", "slope", "se")]
    colnames(d) = c("chr", "snp_end", "SNP", "Phenotype", "pval", "beta", "se")
    d$chr = as.numeric(gsub("chr", "", d$chr))
    d = merge(d, qtl, by=c("chr", "SNP")) ## don't match QTLtools output with just chrom:position becuase indel positions are different.
    d$samplesize = samples
    d$id = d$gene = d$Phenotype
    d$var = glue("chr{d$chr}_{d$position}_{d$other_allele}_{d$effect_allele}") ## position is from the original vcf, use this instead of snp_end from 
    return(d)
}

e = format.qtl(eqtl, qtl, eqtl_samples)
ca = format.qtl(caqtl, qtl, caqtl_samples)

runmr = function(exp.df, out.df, exp, out, hits, direction, outprefix){
    print(glue("{exp} on {out}"))
    vars = unlist(strsplit(hits, ","))
    ## don't match QTLtools output with just chrom:position becuase indel positions are different.
    ## for the hit from SuSiE of the form chr_pos_ref_alt, get the resid, then use the rsid to fetch qtl summary stats
    exposure =  format_data(exp.df[exp.df$Phenotype == exp & exp.df$var %in% vars,], type="exposure")
    outcome = format_data(out.df[out.df$Phenotype == out & out.df$var %in% vars,], type="outcome")
    
    dat <- harmonise_data(exposure, outcome)

    dat$rout = NA
    dat$rexp = NA

    o = mr_steiger(p_exp = dat$pval.exposure, p_out = dat$pval.outcome, n_exp = dat$samplesize.exposure, n_out = dat$samplesize.outcome,
                      r_xxo = 1, r_yyo = 1, r_exp = dat$rexp, r_out = dat$rout)
    print(o)
    odf = data.frame(o[c("r2_exp", "r2_out", "r2_exp_adj", "r2_out_adj", "correct_causal_direction", "steiger_test",
                         "correct_causal_direction_adj", "steiger_test_adj", "vz", "vz0",  "vz1", "sensitivity_ratio")])
    odf$exposure = exp
    odf$outcome = out
    odf$mrs_direction = direction
    write.table(odf, file = glue("{outprefix}.sensitivity.tsv"), sep='\t', quote=F)
                       
    png(glue("{outprefix}.sensitivity.png"))
    print(o$sensitivity_plot)
    dev.off()
}



if (nrow(mrs) > 0){
    print("caqtl on eqtl")    
    ## instead of clumping, select the actual signal hits from SuSiE as instruments
    apply(mrs, 1, function(row){runmr(ca, e, row['capeak'], row['egene'], row['cahit'], "ca-to-e", glue("{prefix}.{row['capeak']}-to-{row['egene']}"))})
}


if (nrow(mrs) > 0){
    print("eqtl on caqtl")
    apply(mrs, 1, function(row){runmr(e, ca, row['egene'], row['capeak'], row['eqhit'], "e-to-ca", glue("{prefix}.{row['egene']}-to-{row['capeak']}"))})
}

