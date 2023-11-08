#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(glue)
ggplot2::theme_set(theme_cowplot())
library(data.table)
library(viridis)
library(gridExtra)
library(ggrepel)
args = commandArgs(TRUE)

prefix = "caqtl"
savep = function(p, f, h=8, w=4){
    print(p)
    png(glue("{prefix}.{f}.png"), height=h, width=w, units="in", res = 150)
    print(p)
    dev.off()
}

covlist = list(
    "covlist.0" = "age-sex-bmi-batch",
    "covlist.1" = "5PCs",
    "covlist.2" = "10PCSs", 
    "covlist.3" = "20PCs",
    "covlist.4" = "30PCs",
    "covlist.5" = "50PCs", 
    "covlist.6" = "5PCs-n_nuclei",
    "covlist.7" = "10PCs-n_nuclei",
    "covlist.8" = "20PCs-n_nuclei",
    "covlist.9" = "30PCs-n_nuclei",
    "covlist.10" = "50PCs-n_nuclei",
    "covlist.11" = "age-sex-bmi-batch-5PCs-n_nuclei-hqaa_umi_median-fraction_mitochondrial_median",
    "covlist.12" = "age-sex-bmi-batch-10PCs-n_nuclei-hqaa_umi_median-fraction_mitochondrial_median",
    "covlist.13" = "age-sex-bmi-batch-20PCs-n_nuclei-hqaa_umi_median-fraction_mitochondrial_median",
    "covlist.14" = "age-sex-bmi-batch-30PCs-n_nuclei-hqaa_umi_median-fraction_mitochondrial_median",
    "covlist.15" = "age-sex-bmi-batch-50PCs-n_nuclei-hqaa_umi_median-fraction_mitochondrial_median"
    )
    
ncovs = list('covlist.0'=0,
               'covlist.1'=5,
               'covlist.2'=10,
               'covlist.3'=20,
               'covlist.4'=30,
               'covlist.5'=50,
               'covlist.6'=5,
               'covlist.7'=10,
               'covlist.8'=20,
               'covlist.9'=30,
               'covlist.10'=50,
               'covlist.11'=5,
               'covlist.12'=10,
               'covlist.13'=20,
               'covlist.14'=30,
               'covlist.15'=50
              )

tcovs = list('covlist.0'="age-sex-bmi-batch",
               'covlist.1'="PCs",
               'covlist.2'="PCs",
               'covlist.3'="PCs",
               'covlist.4'="PCs",
               'covlist.5'="PCs",
               'covlist.6'="PCs-n_nuclei",
               'covlist.7'="PCs-n_nuclei",
               'covlist.8'="PCs-n_nuclei",
               'covlist.9'="PCs-n_nuclei",
               'covlist.10'="PCs-n_nuclei",
               'covlist.11'="age-sex-bmi-batch-PCs-nnuclei-hqaa-mito",
               'covlist.12'="age-sex-bmi-batch-PCs-nnuclei-hqaa-mito",
               'covlist.13'="age-sex-bmi-batch-PCs-nnuclei-hqaa-mito",
               'covlist.14'="age-sex-bmi-batch-PCs-nnuclei-hqaa-mito",
               'covlist.15'="age-sex-bmi-batch-PCs-nnuclei-hqaa-mito"               
              )

covsets = list('allvars' = glue("covlist.{c(0,11,12,13,14,15)}"),
              'pcs' = glue("covlist.{c(seq(1,5,1))}"),
              'pcsn' = glue("covlist.{seq(6,10,1)}"))

colors = c("#a6cee3", "#1f78b4", "#b2dlf8a", "#33a02c",
                 "#6a3dl9a", "#fb9a99", "#fdlbf6f", "#cab2dl6",
                 "#b15928", "#ff7f00","#e31a1c")

fix = function(f, filename, nseq){
    d = read.csv(f, sep='\t', header=T, as.is=T)
    d = d %>% separate(name, c("cluster", "samples", "covlist"), "__")
    d$samples = gsub("_samples", "", gsub("top_", "", d$samples))
    d$covlist = gsub("_", "", d$covlist)
    d$covs = plyr::mapvalues(x = d$covlist,
                             from = paste("covlist", nseq, sep="."),
                             to = covlist)
    d$ncovs = plyr::mapvalues(x = d$covlist,
                             from = paste("covlist", nseq, sep="."),
                             to = ncovs)
    d$tcovs = plyr::mapvalues(x = d$covlist,
                             from = paste("covlist", nseq, sep="."),
                             to = tcovs)
    d$covs = as.character(d$covs)
    d$ncovs = as.character(d$ncovs)
    d$tcovs = as.character(d$tcovs)

    d$samples = as.numeric(d$samples)
    write.table(d, filename, sep='\t', row.names = F, quote=F)
    return(d)
    }

df = fix(args[1], args[2], seq(0, 15, 1))
df$quant = "inorm_tpm"
df$ncovs = as.numeric(df$ncovs)

df[df$covlist==args[3] & df$samples==as.numeric(args[4]),]

lapply(unique(df$cluster), function(i){
    d1 = df[df$cluster == i,]
    p = ggplot(d1, aes(x=ncovs, y=nsig)) +
    geom_line(aes(color=tcovs)) +
    geom_point(aes(color=tcovs)) +
    facet_grid(quant~samples) +
    labs(title=i)
    savep(p, glue("{i}.npcs"), 5, 12)
    })




