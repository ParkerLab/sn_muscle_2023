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
    parser.add_argument('--features', required=True, help ="""Outputs from qtltools permute scan""")
    parser.add_argument('--output',  help ="""Output prefix.""")
    args = parser.parse_args()
    return args

def getfeatures(f):
    p = pandas.read_csv(f, delim_whitespace=True, header=None, names=['gene_name', 'chrom', 'start_pheno', 'end_pheno', 'strand',
                                                                      'n_variants_tested', 'distance_var_pheno', 'snp', 'snp_chrom',
                                                                      'snp_start', 'snp_end', 'degrees_of_freedom', 'estimated_degress_of_freedom', 'beta_first_param',
                                                                      'n_effective_tests', 'p_nominal', 'r2', 'slope', 'se', 'p_empirical', 'p_beta'],
                        usecols = ['gene_name'])
    return p


if __name__ == '__main__':
    
    args = getOpts()

    # concat results across chunks and calculate FDR
    tested = pandas.concat([getfeatures(f) for f in args.qtl])['gene_name'].drop_duplicates().tolist()
    features = pandas.read_csv(args.features, sep='\t', usecols=['pid'], squeeze=True).tolist()
    untested = [ft for ft in features if ft not in tested]

    with open(args.output, 'w') as f:
        for ft in untested:
            f.write(f"{ft}\n")
