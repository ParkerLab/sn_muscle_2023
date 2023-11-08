#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas
import numpy
import math
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os
import glob
import time
import pybedtools
dpi = 150
matplotlib.rcParams['figure.dpi']= dpi
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams.update({'font.size': 12})
import argparse

def getOpts():
    parser = argparse.ArgumentParser(description='LD-based colod b/w caQTL and GWAS')
    parser.add_argument('--permute', required=True, nargs="+", help ="""Formatted permutation scan""")
    parser.add_argument('--conditional', required=True, nargs="+",  help ="""Formatted conditional scan""")
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    
    args = getOpts()
    
    permutes = {os.path.basename(f).split('--')[0] : f for f in args.permute}
    conditionals = {os.path.basename(f).split('--')[0] : f for f in args.conditional}
    
    clusters = [c for c in permutes.keys() if c in conditionals.keys()]

    for cluster in clusters: 
    
        # assemble primary and conditional signals
        # primary
        perm = pandas.read_csv(permutes[cluster], sep='\t',
                               usecols=['gene_name', 'distance_var_pheno', 'snp', 'snp_chrom', 'snp_end', 'p_nominal', 'slope', 'qvalue_storey'])
        perm = perm[perm['qvalue_storey']<0.05]
        perm['association_rank'] = 0
        
        # conditional 
        cond = pandas.read_csv(conditionals[cluster], sep='\t',
                               usecols=['gene_name', 'distance_var_pheno', 'association_rank', 'snp', 'snp_chrom', 'snp_end', 'backward_p', 'backward_slope'])
        cond = cond[cond['association_rank']>0]
        
        cond.rename(columns= {'backward_p': "p_nominal",
                              'backward_slope': "slope"}, inplace=True)

        s = pandas.concat([perm.drop('qvalue_storey', axis=1), cond])

        s.to_csv(f"{cluster}.caqtl.tsv", sep='\t', index=False, na_rep="NA")
        
        
