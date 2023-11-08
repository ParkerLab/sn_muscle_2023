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
    parser.add_argument('--sample-covs',  help ="""tab separated file with cluster-sample info such as age sex bmi n_nuclei etc """)
    parser.add_argument('--sample-col', default="SNG.1ST",  help ="""Column for sample name in sample-covs file. [Default = SNG.1ST]""")    
    parser.add_argument('--pheno-pca', default="coarse_cluster_name", help ="""phenotype pca output""")
    parser.add_argument('--geno-pca', default="SNG.1ST", help ="""genotype pca output""")
    parser.add_argument('--output', required=True, help ="""output pca""")
    args = parser.parse_args()
    return args

def add_pcs(pcfile, df, sample_col, prefix):
    pc = pandas.read_csv(pcfile, delim_whitespace=True, index_col=0).T.reset_index()
    pc.rename(columns = {c: f"{prefix}-{c.split('_')[-1]}" for c in pc.columns if c != "index"}, inplace=True)
    if not df.empty:
        pc.rename(columns = {"index": sample_col}, inplace=True)
        df = pandas.merge(df, pc, how="inner", on=sample_col)
    else:
        df = pc
    return df

if __name__ == '__main__':
    
    args = getOpts()

    if args.sample_covs is not None:
        sample_col = args.sample_col
        d = pandas.read_csv(args.sample_covs, sep='\t', dtype={sample_col: str})

    else:
        sample_col = "index"
        d = pandas.DataFrame()
        
    d = add_pcs(args.pheno_pca, d, sample_col, 'pheno')
    d = add_pcs(args.geno_pca, d, sample_col, 'geno')
    
    dc = d.set_index(sample_col).T.reset_index()
    dc.rename(columns={'index':'id'}, inplace=True)
    dc.to_csv(args.output, sep='\t', index=False)
