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
    parser.add_argument('--info', required=True,  help ="""tab separated file with sample info such as age sex bmi etc """)
    parser.add_argument('--vars', required=True, nargs='+', help ="""space separated list of variables to include - should be column names in info file""")
    parser.add_argument('--sample-col', required=True, help ="""Column name for samples (should match vcf..). This column name will be renamed to 'sample'.""")
    parser.add_argument('--cluster-col', required=True, help ="""Column name for cluster""")
    parser.add_argument('--cluster-name', required=True, help ="""Name of cluster for which covariates are being made""")
    parser.add_argument('--pca', required=True, help ="""qtltools pca output""")
    parser.add_argument('--pcs', type=int, help ="""Number of PCs to include""")
    parser.add_argument('--aggregate-sum', nargs='+', help = """On these columns, get sum of values per sample per cluster """)
    parser.add_argument('--aggregate-median', nargs='+', help = """On these columns, get median values per sample per cluster """)
    parser.add_argument('--replace', help ="""Replace this string from PC headers""")
    parser.add_argument('--output',  help ="""Output file""")
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    
    args = getOpts()
    
    cols = [args.sample_col, args.cluster_col] + args.vars + args.aggregate_sum + args.aggregate_median
    df = pandas.read_csv(args.info, sep='\t')

    for c in cols:
        assert c in df.columns
        
    df = df[cols].drop_duplicates()
    df = df[df[args.cluster_col] == args.cluster_name]
    df.rename(columns = {args.sample_col: "sample"}, inplace=True)

    functions = {}
    if args.aggregate_sum is not None:
        for c in args.aggregate_sum:
            functions[c] = 'sum'
    if args.aggregate_median is not None:
        for c in args.aggregate_median:
            functions[c] = 'median'

    if not bool(functions):
        df = df.groupby([args.cluster_col, args.sample_col]).agg(functions)

    df.drop_duplicates(inplace = True)
    t = pandas.read_csv(args.pca, delim_whitespace=True, index_col=0).T.reset_index()
    rename = {}
    for c in t.columns:
        if c == "index":
            rename[c] = "sample"
        else:
            rename[c] = c.replace(args.replace, "")
    print(rename)
    t.rename(columns=rename, inplace=True)

    print(t.head())
    d = pandas.merge(df, t[['sample'] + [f"PC{x}" for x in list(range(1, args.pcs))]], how="inner", on="sample")
    d['batch'] = d['batch'].map(lambda x: f"b{x}")

    print(d.head())
    dc = d.set_index('sample').T.reset_index()
    dc.rename(columns={'index':'id'}, inplace=True)
    dc.to_csv(args.output, sep='\t', index=False)
