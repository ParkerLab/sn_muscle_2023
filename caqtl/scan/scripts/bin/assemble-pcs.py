#!/usr/bin/env python
# coding: utf-8

import pandas
import numpy
import sys
import os
import glob
import argparse

def getOpts():
    parser = argparse.ArgumentParser(description='make qtltools covariates file')
    parser.add_argument('--pheno', required=True,  help ="""Phenotype pcs""")
    parser.add_argument('--geno', required=True, help ="""Genotype pcs""")
    parser.add_argument('--output', required=True, help ="""output file""")
    parser.add_argument('--sample-info',  help ="""Sample info such as SNG.1ST, age, sex, batch, bmi will be fetched from this df""")
    parser.add_argument('--ataqv', help ="""Sample ataqv""")
    parser.add_argument('--cluster', required=True, help ="""cluster name""")

    args = parser.parse_args()
    return args

def uniq(x):
    s = list(set(x))
    assert len(s) == 1
    return s[0]

def get_pc_list(n):
    return [f'PC{i}' for i in range(1, n+1)]
   
if __name__ == '__main__':
    
    args = getOpts()
    phenofile = args.pheno
    genofile = args.geno

    p = pandas.read_csv(phenofile, sep=' ')
    p['SampleID'] = p['SampleID'].map(lambda x: f"pheno-{x.split('_')[-1]}")

    print(p.head())
    
    g = pandas.read_csv(genofile, sep=' ')
    g['SampleID'] = g['SampleID'].map(lambda x: f"geno-{x.split('_')[-1]}")
    print(g.head())
    
    dflist = [p, g]
    if args.sample_info is not None:
        si = pandas.read_csv(args.sample_info, sep='\t', dtype={'SNG.1ST': str})
        si = si[(si['coarse_cluster_name']==args.cluster) & (si['modality']=="atac")].set_index("SNG.1ST").T.reset_index().rename(columns={'index': 'SampleID'})
        if not si.empty:
            dflist.append(si)
        
    if args.ataqv is not None:
        aq = pandas.read_csv(args.ataqv, sep='\t', usecols=['name', 'tss_enrichment']).drop_duplicates()
        aq[['cluster', 'sample']] = aq['name'].str.split('-', n=1, expand=True)
        aq[['SNG.1ST', 'batch']] = aq['sample'].str.split('--', expand=True)
        aq = aq[aq['cluster'] == args.cluster][['SNG.1ST','tss_enrichment']].drop_duplicates().set_index('SNG.1ST').T.reset_index().rename(columns={'index': 'SampleID'})
        if not aq.empty:
            dflist.append(aq)
    
    
    d = pandas.concat(dflist).dropna(axis='columns')
    d.rename(columns={'SampleID': 'id'}).to_csv(args.output, sep='\t', index=False)
