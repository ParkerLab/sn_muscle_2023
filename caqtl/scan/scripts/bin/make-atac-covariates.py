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
    parser.add_argument('--sample-ataqv', required=True,  help ="""cluster-sample level ataqv metrics""")
    parser.add_argument('--modality', default="atac",  help ="""modality to select. Default = atac""")
    parser.add_argument('--cohort', default=['FUSION'], nargs='+', help ="""cohort to select. Default = FUSION""")
    parser.add_argument('--output',  help ="""Output file""")
    parser.add_argument('--min-nuclei', nargs='+', type=int, help ="""Number of min nuclei in a sample for it to be selected""")
    args = parser.parse_args()
    return args

def uniq(x):
    s = list(set(x))
    assert len(s) == 1
    return s[0]

if __name__ == '__main__':
    
    args = getOpts()

    d = pandas.read_csv(args.info, sep='\t',
                        usecols=['index','SNG.1ST', 'sample', 'cohort', 'modality', 'age', 'sex', 'bmi', 'batch', 'coarse_cluster_name', 'hqaa_umi', 'fraction_mitochondrial'],
                        dtype={'SNG.1ST': str}).drop_duplicates()

    d = d[(d['cohort'].isin(args.cohort)) & (d['modality']== args.modality)]

    cluster_col = "coarse_cluster_name"
    sample_col = "sample"
    d1 = d.groupby([cluster_col, sample_col]).agg({'hqaa_umi': [numpy.mean, numpy.median, sum],
                                                   'fraction_mitochondrial': numpy.median,
                                                   'index': len,
                                                   'age': uniq,
                                                   'sex': uniq,
                                                   'bmi': uniq,
                                                   'batch': uniq}).reset_index()
    d1.columns = ["_".join(c) if "" not in c else "".join(c) for c in d1.columns]
    d1.columns = [x.replace("_uniq", "") for x in d1.columns]
    d1.rename(columns = {'index_len': 'n_nuclei'}, inplace=True)
    
    i = pandas.read_csv(args.sample_ataqv, sep='\t', usecols=['name', 'tss_enrichment']).drop_duplicates()
    
    d1['name'] = d1.apply(lambda x: f"{x['coarse_cluster_name']}-{x['sample']}", axis=1)
    
    d1 = pandas.merge(d1, i, how="inner", on="name")
    d1.drop(['name'], axis=1, inplace=True)
    d1['SNG.1ST'] = d1['sample'].map(lambda x: x.split("--")[0])
    d1.to_csv(args.output, sep='\t', index=False, na_rep="NA")

    
    for i, grp in d1.groupby("coarse_cluster_name"):
        for n_nuclei in args.min_nuclei:
            filename = f"{i}__sampleslist.min_{n_nuclei}_nuc_samples.txt"
            grpdf = grp[grp['n_nuclei'] >= n_nuclei][['sample','SNG.1ST']].drop_duplicates()
            assert len(grpdf.index) == len(grpdf['sample'].drop_duplicates().tolist()) == len(grpdf['SNG.1ST'].drop_duplicates().tolist())
            grpdf[['sample']].to_csv(filename, sep='\t', header=False, index=False)
