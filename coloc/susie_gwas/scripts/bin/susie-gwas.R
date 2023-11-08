#!/usr/bin/env Rscript

library(susieR)
library(coloc)
library(glue)
suppressPackageStartupMessages(library(dplyr))
library(ggplot2)
library(optparse)
library(gtable)
library(grid)
library(gridExtra)
library(cowplot)
ggplot2::theme_set(theme_cowplot())

plot.set = function(setlist, trait, input, title, marker = NULL){
    print(marker)
    input$set = "none"
    for (i in seq_along(setlist)){
        input$set = ifelse(input$snp %in% names(setlist[[i]]), glue("set{i}"), input$set)
    }
    
    if (trait == "gwas"){
        p = ggplot(input, aes(x=snp_end, y=-log10(pval))) +
            geom_point(aes(color=set))
    } else {
        p = ggplot(input, aes(x=snp_end, y=-log10(p_nominal))) +
            geom_point(aes(color=set))
    }

    if (! is.null(marker)){
        p = p +
            geom_point(data = input[input$snp == marker,], shape=1, color="blue")
    }
    p = p +
        labs(x="position", title=title) +
        theme(legend.position="bottom", plot.title=element_text(size=7))
    return(p)
}

plot.result = function(result.vector, name, df, marker = NULL){
    new = list()
    new[[name]] = result.vector
    df = cbind(df, data.frame(new))
    p = ggplot(df, aes_string(x="snp_end", y=name)) +
        geom_point(shape=16, color="grey")
    if (! is.null(marker)){
        p = p +
            geom_point(data = df[df$snp == marker,], shape=1, color="blue")
    }
    
    return(p)
}

process_dosage = function(f, snplist){
    ld = read.csv(f, sep='\t', check.names = F)
    dups = ld[ (duplicated(ld$ID) | duplicated(ld$ID, fromLast = TRUE)),]
    print(glue("N duplicates = {nrow(dups)}"))
    ld = ld[! (duplicated(ld$ID) | duplicated(ld$ID, fromLast = TRUE)),]
    row.names(ld) = ld$ID
    ld$ID = NULL
    idlist = intersect(snplist, row.names(ld))
    ld = ld[snplist,]
    ld = cor(t(ld))
    out = list("ld" = ld, "idlist" = idlist)
    return(out)
}

get_var = function(df, trait, betacol, secol, pcol){
    vcol = glue("{trait}_var")
    zcol = glue("{trait}_z")

    if (is.null(secol)) {
        df[[zcol]] = (df[[betacol]]/abs(df[[betacol]]))*qnorm(df[[pcol]]/2)
        df[[vcol]] = (df[[betacol]]/df[[zcol]])^2
    } else {
        df[[vcol]] = df[[secol]]^2
        df[[zcol]]=df[[betacol]]/df[[secol]]
    }
    return(df)
}

make_coloc_dataset = function(df, beta, var, type, position, snp, sdY=NULL, maf=NULL, N=NULL, LD=NULL) {
    D = list(beta = df[[beta]], varbeta = df[[var]],
             snp = df[[snp]], position = df[[position]], type = type, LD = LD)

    if (type == "quant"){
        if ( ! is.null(sdY)){
            D = c(D, list(sdY = sdY))
        } else {
            stopifnot((! is.null(maf)) & (! is.null(N)) )
            D = c(D, list(MAF = df[[maf]], N = df[[N]]))
        }
    }
    return(D)
}

susie_dropsets = function(res, drop){
    res$sets$cs = res$sets$cs[-drop]
    res$sets$cs_index = res$sets$cs_index[-drop]
    res$sets$coverage = res$sets$coverage[-drop]
    res$sets$purity = res$sets$purity[-drop,,drop=FALSE]
    return(res)
}

option_list <- list(
    make_option(c("--trait1"), type = "character", help = "[Required] trait1 dataframe"),
    make_option(c("--type"), type = "character", help = "[Required] trait1 type-  cc or quant"),
    make_option(c("--beta"), type = "character", help = "[Required] column name for trait1 beta"),
    make_option(c("--se"), type = "character", help = "column name for trait1 se. If not available, provide p value column in --trait1_p"),
    make_option(c("--p"), type = "character", help = "column name for p value - will calculate z and se if se is not available for trait1"),
    make_option(c("--sdY"), type = "numeric", help = "If trait1 is quant, sdY (phenotype). 1 if the phenotypes were inverse normalized"),
    make_option(c("--maf"), type = "character", help = "if trait1 is quant and sdY is not known, provide the column name for maf"),
    make_option(c("--N"), type = "character", help = "if trait1 is quant and sdY is not known, provide the column name for N"),
    make_option(c("--ld"), type = "character", help = "[Required] genotype dosage for trait1 reference panel"),
    make_option(c("--ld_mat"), type = "character", default=NULL, help = "genotype dosage correlation matrix for trait1 reference panel. Can be provided instead of the --trait1_ld option. [DEFAULT = NULL]"),
    make_option(c("--coverage"), type = "numeric", default=0.95, help = "SuSiE coverage [DEFAULT = 0.95]"),
    make_option(c("--r2_prune"), type = "numeric", default=0.8, help = "SuSiE r2.prune [DEFAULT = 0.8]"),
    make_option(c("--maxit"), type = "numeric", default=10000, help = "SuSiE max iterations [DEFAULT = 10000]"),
    make_option(c("--min_abs_corr"), type = "numeric", default=0.1, help = "SuSiE min_abs_corr [DEFAULT = 0.1]"),
    make_option(c("--number_signals"), type = "numeric", default=10, help = "SuSie L [DEFAULT = 10]"),
    make_option(c("--s_threshold"), type = "numeric", default=0.3, help = "max s threshold above which change L [DEFAULT = 0.3]"),
    make_option(c("--number_signals_high_s"), type = "numeric", default=1, help = "SuSiE L if s > max s_threshold [DEFAULT = 1]"),
    make_option(c("--reduce_corr"), type = "numeric", default=NULL, help = "If no credible sets are found, sequentially reduce the min_abs_corr by 0.05 until this value is reached and see if sets can be found. [Default = NULL]"),
    make_option(c("--reduce_coverage"), type = "numeric", default=NULL, help = "If no credible sets are found, sequentially reduce the coverage by 0.05 until this value is reached and see if sets can be found. [Default = NULL]"),
    make_option(c("--marker"), type = "character", default=NULL, help = "highlight in red an inteded lead SNP by providing rsid"),
    make_option(c("--dropsets_if_not_contain"), type = "character", default=NULL, help = "drop the credible set if it doesn't contain this rsID"),
    make_option(c("--dropsets_if_no_proxy"), type = "numeric", default=NULL, help = "drop the credible set if it doesn't contain this rsID or at least 1 proxy with this r2 or greater. Flag --dropsets_if_not_contain if provided takes precedent over this."),
    make_option(c("--prefix"), type = "character", help = "output prefix")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

## read vars
prefix = opts$prefix
t1.beta = opts$beta
t1.se = opts$se
t1.p = opts$p
t1.type = opts$type
t1.sdY = opts$sdY
t1.maf = opts$maf
t1.N = opts$N
t1.ldmat = opts$ld_mat

## runsusie GWAS only technically needs beta and se since it is case control, but for calculating the s metric, we also need Z scores. i.e using p values column.
## read trait1 
input1 = read.csv(opts$trait1, sep='\t', header=T, as.is=T, check.names=F)
## remove all duplicates
input1 = input1[!duplicated(input1), ]
duplicated.input = input1[duplicated(input1[, c('#snp_chrom', 'snp', 'snp_end')]) | duplicated(input1[, c('#snp_chrom', 'snp', 'snp_end')], fromLast=TRUE), ]
if  (nrow(duplicated.input) > 0){
    print(glue("WARNING: {nrow(duplicated.input)} duplicates found in data. Removing for now. Check data."))
    print(duplicated.input)
    input1 = dplyr::setdiff(input1, duplicated.input)
}
input1 = get_var(input1, "trait1", betacol = t1.beta, secol = t1.se, pcol = t1.p)
snps = input1$snp
print(head(input1))
## 
## read or process and save LD matrix
if (! is.null(t1.ldmat)){
    t1.ld = read.csv(t1.ldmat, sep='\t', header=T, row.names=1, as.is=T, check.names=F)
    idlist = intersect(row.names(t1.ld), snps)
    ldf = t1.ld[idlist, idlist]
    t1.ld = as.matrix(ldf)
    rownames(ldf) = gsub("-.*", "", rownames(ldf))
    colnames(ldf) = gsub("-.*", "", colnames(ldf))
} else {
    t1.ld.out = process_dosage(opts$trait1_ld, snps)
    print("dosage done")
    ## LD mat
    t1.ld = t1.ld.out[['ld']]
    idlist = t1.ld.out[['idlist']]
    write.table(t1.ld, sep='\t', quote=F, file=glue("{prefix}.ld.tsv"))

}
## update input df - in case some variants were missing from the ld matrix.
row.names(input1) = input1$snp
input1 = input1[idlist,]

## make coloc dataset
position = "snp_end"
snp = "snp"
D1 = make_coloc_dataset(input1, beta = t1.beta, var = "trait1_var", type = t1.type,
                        position = position, snp = snp, sdY = t1.sdY, maf = t1.maf, N = t1.N,
                        LD = t1.ld)

## Set default metrics
coverage = opts$coverage
r2.prune = opts$r2_prune
maxit = opts$maxit
min_abs_corr = opts$min_abs_corr

## Calculate diagnostic
diagnostic_s_gwas = estimate_s_rss(input1$trait1_z, t1.ld)
print(glue("s = {diagnostic_s_gwas}"))

## plot diagnostic
png(glue("{prefix}.s.png"), height=4, width=4, units="in", res=150)
condz_in = kriging_rss(input1$trait1_z, t1.ld)
print(condz_in$plot)
dev.off()

## if diagnostic s is very bad, set max signals to 1
if(diagnostic_s_gwas > opts$s_threshold){
    print("Setting numbersignals to 1")
    numbersignals = opts$number_signals_high_s
} else {
    numbersignals = opts$number_signals
}

## runsusie 
S1 = coloc::runsusie(D1, coverage = coverage, r2.prune = r2.prune, maxit = maxit,
                     min_abs_corr = min_abs_corr, repeat_until_convergence = F, L = numbersignals, estimate_prior_variance=FALSE)
    
## if fail to converge, reduce max number of signals
while (S1$converged == F  & (numbersignals > 0)) {
    numbersignals = numbersignals - 1
    S1 = coloc::runsusie(D1, coverage = coverage, r2.prune = r2.prune, maxit = maxit,
                         min_abs_corr = min_abs_corr, repeat_until_convergence = F, L = numbersignals)
    if (S1$converged) {
        print(glue("Converged after considering max signals = {numbersignals}"))
    } else {
        print(glue("Didn't converge with max signals = {numbersignals}"))
    }
} # end max signals while loop
    
## if no signals found, reduce minimum abs correlation if that option was selected
if (!is.null(opts$reduce_corr)){
    while ((is.null(S1$sets$cs)) & (min_abs_corr >= opts$reduce_corr + 0.05) & (S1$converged == T)) {
        min_abs_corr = min_abs_corr - 0.05
        S1 = coloc::runsusie(D1, coverage = coverage, r2.prune = r2.prune, maxit = maxit,
                             min_abs_corr = min_abs_corr, repeat_until_convergence = F, L = numbersignals)
        if (! is.null(S1$sets$cs)){
            print(glue("Csets found after considering min_abs_corr = {min_abs_corr}"))
        } else {        
            print(glue("No csets found with min_abs_corr = {min_abs_corr}"))
        }
    } # end min abs correlation loop
}

## if still no signals found, reduce coverage sequentially if that option was selected
if (!is.null(opts$reduce_coverage)){
    print(glue("will try to find signals after reducing coverage to {opts$reduce_coverage}"))
    step = 0.05
    while ((is.null(S1$sets$cs)) & (round(coverage - step, 2) >= opts$reduce_coverage) & (S1$converged == T)) {
        coverage = coverage - step
        S1 = coloc::runsusie(D1, coverage = coverage, r2.prune = r2.prune, maxit = maxit,
                             min_abs_corr = min_abs_corr, repeat_until_convergence = F, L = numbersignals)
        if (! is.null(S1$sets$cs)){
            print(glue("Csets found after considering coverage = {coverage}"))
        } else {
            print(glue("No csets found with coverage = {coverage}"))
        }
    }
}


## write susie results
sfname_gwas = glue("{prefix}.Rda")
save(S1, file=sfname_gwas) 


results = c(glue("prefix\tdiagnostic_s_gwas\tmin_abs_corr\tnumbersignals\tcoverage\tnsets\tlensets\tn_selected_sets\tlen_selected_sets"))

if (is.null(S1$sets$cs)){
    nsets = 0
    lensets = 0
    n_selected_sets = 0
    len_selected_sets = 0
} else {
    nsets = length(S1$sets$cs)
    lensets = paste(lapply(S1$sets$cs, function(i){length(i)}), collapse = ',')
    n_selected_sets = nsets
    len_selected_sets = lensets
}


## color csets over GWAS P values
## intended lead SNP:
if (! is.null(opts$marker)){
    marker.id = input1[gsub("-.*", "", input1$snp) == opts$marker,]$snp
} else {
    marker.id = NULL
}
title = glue("{prefix}\ns = {diagnostic_s_gwas}, min_corr = {min_abs_corr}\n L = {numbersignals}, coverage = {coverage}\nnsets = {nsets}")
p1 = plot.set(S1$sets$cs, "gwas",  input1, title, marker = marker.id)

## plot PIP
p2 = plot.result(S1$pip, "pip", input1, marker = marker.id)
## g2 = ggplotGrob(p2)

png(glue("{prefix}.png"), height=6, width=8, units="in", res=150)
cowplot::plot_grid(p1, p2, align = 'v', ncol = 1, rel_heights = c(0.6, 0.4)) ## display plotpip
dev.off()

## Select the credible set to use downstream. This can be either that set that contains the intended lead SNP itself or a proxy in high LD with it.
if (length(S1$sets$cs)>0) {
    print("There was at least 1 cset identified")
    drop = c()
    check.proxy = FALSE
    stop = FALSE
    
    if (!is.null(opts$dropsets_if_not_contain)) {
        check.match = opts$dropsets_if_not_contain
        ## see if any sets to drop
        print(glue("Will remove any sets that don't contain {check.match}"))
    } else if (!is.null(opts$dropsets_if_no_proxy)) {
        ## see if each cset either contains the exact SNP or a proxy
        check.match = opts$marker
        check.proxy = opts$dropsets_if_no_proxy
        print(glue("Will look for proxy of lead SNP r2 > {check.proxy}"))
    } else {
        stop = TRUE
        print("no selection option provided, leaving.")
    }
    
    if (! stop) {
        ## see if any sets to drop
        for (i in seq_along(S1$sets$cs)){
            cset.snps = gsub("-.*", "", names(S1$sets$cs[[i]]))
            if (! check.match %in% rownames(ldf)){
                print(glue("SNP {check.match} didn't occur in the LD matrix. Exit."))
                break
            }
            if (! check.match %in% cset.snps){
                print(glue("SNP {check.match} not found in set {i}"))
                if (! check.proxy){
                    drop = c(drop, i)
                    print(glue("removing set {i}"))
                } else {
                    ## check for proxy
                    cset.r2 = ldf[check.match, cset.snps]**2
                    max.cset.r2 = max(cset.r2)
                    if (max.cset.r2 < check.proxy){
                        drop = c(drop, i)
                        print(glue("Max r2 with of marker with cset SNP {max.cset.r2} < required {check.proxy}; removing set {i}"))
                    } else {
                        print(glue(("Max r2 with of marker with cset SNP {max.cset.r2} >= required {check.proxy}; retaining set {i}")))
                    }
                }
            } else {
                print(glue("SNP {check.match} found in set {i}"))
            }
        }
    }
    
    ## drop all such sets
    if (length(drop)>0){
        print(glue("Started with number of credible sets = {length(S1$sets$cs)}"))
        S1 = susie_dropsets(S1, drop)
        print(glue("Final number of selected credible sets = {length(S1$sets$cs)}"))
    }
    
    n_selected_sets = length(S1$sets$cs)
    len_selected_sets = paste(lapply(S1$sets$cs, function(i){length(i)}), collapse = ',')

    ## if any selected csets remain, save the .selected.Rda file and plot the selected gwas, pip and lbf
    if (n_selected_sets > 0){
        fname = glue("{prefix}.selected.Rda")
        save(S1, file = fname)
        
        ## plot lbf
        idx = S1$sets$cs_index
        isnps = colnames(S1$lbf_variable)
        bf = S1$lbf_variable[idx, isnps, drop=FALSE]
        p3 = plot.result(t(bf)[,1], "lbf", input1, marker = marker.id)

        png(glue("{prefix}.selected.png"), height=9, width=8, units="in", res=150)
        print(cowplot::plot_grid(plotlist = list(p1, p2, p3), align = 'v', ncol = 1, nrow=3)) ## display plotpip
        dev.off()

    }
}
results = c(results, glue("{prefix}\t{diagnostic_s_gwas}\t{min_abs_corr}\t{numbersignals}\t{coverage}\t{nsets}\t{lensets}\t{n_selected_sets}\t{len_selected_sets}"))
writeLines(results, glue("{prefix}.results.tsv"))

