#!/usr/bin/env python3
# coding: utf-8

# In[28]:


import pandas
import numpy
import sys
import os
import glob

# argv1 = ukbb vcf
# argv2 = gwas df
v = pandas.read_csv(sys.argv[1], sep='\t', comment="#", header=None,
                    usecols=[0,1,2,3,4]).rename(columns = {0:"#snp_chrom", 
                                                           1:"snp_end",
                                                           2:"snp",
                                                           3:"REF",
                                                           4:"ALT"})
# to avoid duplicates from UKBB, update snp names to include the alleles
v['snp'] = v.apply(lambda x: f"{x['snp']}-{x['REF']}-{x['ALT']}", axis=1)

g = pandas.read_csv(sys.argv[2], sep='\t')
g = pandas.merge(g, v, how="inner", on = ['#snp_chrom', 'snp_end']).drop_duplicates()

g['EA'] = g['EA'].str.upper()
g['NEA'] = g['NEA'].str.upper()

# if there are duplicates for the chrom pos level (multiallelic so UKBB vcf has multiple entries), we should retain where the alt allele == EA..?
dups = g[g.duplicated(['#snp_chrom', 'snp_end'], keep=False)]
if not dups.empty:
    print("Some duplicates found, likely multiallelic SNPs in the UKBB vcf")
    print(dups)
    guniq = g[~g.duplicated(['#snp_chrom', 'snp_end'], keep=False)]
    dups = dups[dups['EA']==dups['ALT']]
    print("After selecting where the EA is the ALT")
    print(dups)
    g = pandas.concat([guniq, dups])

mismatch = g[~(((g['EA'] == g['REF']) & (g['NEA'] == g['ALT'])) | ((g['EA'] == g['ALT']) & (g['NEA'] == g['REF']) ))]
if not mismatch.empty:
    print(f"mismatches = {len(mismatch.index)}")
    print(mismatch)

prev = len(g.index)
# retain only markers where positions and both alleles match
g = g[(((g['EA'] == g['REF']) & (g['NEA'] == g['ALT'])) | ((g['EA'] == g['ALT']) & (g['NEA'] == g['REF']) ))]

print(f"prev = {prev}; new = {len(g.index)}")
# update beta sign if the effect allele was not the alt allele
g['beta'] = g.apply(lambda x: -1*x['beta'] if x['EA']!=x['ALT'] else x['beta'], axis=1)
g['EA'] = g['ALT']
g['NEA'] = g['REF']

g.to_csv(sys.argv[3], sep='\t', index=False, na_rep="NA")


