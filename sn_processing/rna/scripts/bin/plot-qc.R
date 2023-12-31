#!/usr/bin/env Rscript

library(ggplot2)
library(ggpointdensity)
library(viridis)

args = commandArgs(TRUE)

d = read.csv(args[1], sep='\t')

makeplot = function(d){
    p = ggplot(d) +
        geom_pointdensity(aes(log_umis, log_mito), size = 0.5, shape = 16) +
        scale_color_viridis() +
        facet_wrap(~library, nrow = 2) +
        labs(y = 'log10(Mitochondrial fraction)', x = 'log10(# UMIs)') +
        theme_bw()
    
    return(p)
}

png(args[2], height=3.5, width=8, units="in", res=150)
print(makeplot(d))
dev.off()

