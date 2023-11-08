#!/bin/bash

################## running coloc for  ${coloc_locus}, ${trait2_name}, ${trait2_locus}: 

## fetch caQTL variants 

tabix -H ${t2h} > caqtl.header ;
for i in ${t2_summary}; do tabix $$i $fetch ; done | grep ${trait2_locus} | cat caqtl.header -  > caqtl.tsv ;

## Coloc 
coloc.R --trait1 ${coloc_locus}.gwas.tsv --trait2 caqtl.tsv --prefix ${prefix} --trait1_ld ${coloc_locus}.ukbb-dosages.tsv --trait2_ld ${coloc_locus}.fusion-dosages.tsv ${t1_params} ${t2_params} ;
		
# rm caqtl.tsv 
