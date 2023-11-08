#!/bin/bash

################## running SuSiE for  ${susie_locus}:

## Susie 
susie-gwas.R --trait1 ${susie_locus}.gwas.tsv --prefix ${susie_locus} --ld_mat ${susie_locus}.ld.tsv  ${t1_params} ;
		
