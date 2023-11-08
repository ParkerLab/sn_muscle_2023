#!/usr/bin/env python
# coding: utf-8

import pandas
import numpy
import sys
import os
import glob
from scipy import stats
from statsmodels.stats.multitest import multipletests
import argparse

def getOpts():
    parser = argparse.ArgumentParser(description='compile eQTL results across chunks')
    parser.add_argument('--qtl', required=True, nargs='+', help ="""Outputs from qtltools permute scan""")
    parser.add_argument('--output',  help ="""Output prefix.""")
    args = parser.parse_args()
    return args

def fixdf(f):
    p = pandas.read_csv(f, delim_whitespace=True, header=None, names=['gene_name', 'chrom', 'start_pheno', 'end_pheno', 'strand',
                                                                      'n_variants_tested', 'distance_var_pheno', 'snp', 'snp_chrom',
                                                                      'snp_start', 'snp_end', 'p_nominal', 'slope', 'is_top'],
                        usecols = ['snp_chrom', 'snp_start', 'snp_end', 'snp', 'gene_name', 'start_pheno', 'end_pheno', 'strand',
                                   'n_variants_tested', 'distance_var_pheno', 'p_nominal', 'slope', 'is_top'])
    p.dropna(axis=0, inplace=True)
    return p
if __name__ == '__main__':
    
    args = getOpts()

    # concat results across chunks and calculate FDR
    d = pandas.concat([fixdf(f) for f in args.qtl])

    d.to_csv(args.output, sep='\t', index=False, na_rep="NA")
