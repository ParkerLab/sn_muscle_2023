#!/bin/bash

## fetch variants in the region and intersect UKBB and FUSION vcfs
for i in $t1_vcf; do tabix $$i $ukbb_fetch | awk '{if (($$0 !~ /^#/ && $$0 !~ /^chr/)) print "chr"$$0; else print $$0}' ; done | sort | uniq > ukbb.genotypes
zcat $t1_vcf1 | head -10000 | awk '{if (($$0 ~ /^#/)) print $$0}' > ukbb.header
cat ukbb.header ukbb.genotypes | bgzip -c > ukbb.vcf.gz; tabix ukbb.vcf.gz
rm ukbb.genotypes ukbb.header

## fetch UKBB dosages 
zcat ukbb.vcf.gz | head -10000 | awk -F'\t' '{if (($$0 ~/^#CHROM/)) print $$0}' OFS='\t' | sed -e 's:#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT:ID:g' > ukbb-header.txt 
bcftools query -f "%ID-%REF-%ALT[\t%DS]\n" ukbb.vcf.gz | cat ukbb-header.txt - > ${susie_locus}.ukbb-dosages.tsv 

## fetch GWAS variants 
tabix -h $t1_summary $fetch > gwas.tsv ;

## align GWAS alleles with UKBB reference and have consistent rsids
align-gwas-refpanel-alleles.py ukbb.vcf.gz gwas.tsv ${susie_locus}.gwas.tsv

## cleanup
rm -rf ukbb-header.txt gwas.tsv ukbb.vcf.gz*

