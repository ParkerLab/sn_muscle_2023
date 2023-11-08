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
    parser.add_argument('--pheno-pcs', nargs='+', default=[30], type=int, help ="""Number of top phenotype PCs to take""")
    parser.add_argument('--geno-pcs', nargs='+', default=[4], type=int, help ="""Number of top genotype PCs to take""")
    parser.add_argument('--other-vars', nargs='+', type=str, help ="""Other covariate names""")

    args = parser.parse_args()
    return args

def uniq(x):
    s = list(set(x))
    assert len(s) == 1
    return s[0]

def get_pc_list(n, prefix):
    return [f'{prefix}-PC{i}' for i in range(1, n+1)]
   
if __name__ == '__main__':

    args = getOpts()
    
    # MAke covariate lists
    n_pheno_pcs = args.pheno_pcs
    n_geno_pcs = args.geno_pcs
    other_vars = args.other_vars
    
    # covs = [get_pc_list(n, 'pheno') + get_pc_list(ng, 'geno') for  n in n_pheno_pcs for ng in n_geno_pcs] + \ ## pheno + geno only
    covs = [other_vars + get_pc_list(ng, 'geno') for ng in n_geno_pcs] + [get_pc_list(n, 'pheno') + get_pc_list(ng, 'geno') for n in n_pheno_pcs for ng in n_geno_pcs] + [get_pc_list(n, 'pheno') + get_pc_list(ng, 'geno') + other_vars for n in n_pheno_pcs for ng in n_geno_pcs]

    print(covs)
    for i, clist in enumerate(covs):
        filename = f"covlist.{i}.txt"
        with open(filename, 'w') as f:
            f.write('\n'.join(clist))

    
