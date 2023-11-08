#!/usr/bin/env python
# coding: utf-8

# In[1]:

import pandas
import numpy
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import os
dpi = 150
matplotlib.rcParams['figure.dpi']= dpi
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams.update({'font.size': 12})
import argparse
import sys

# Filter features now:
def fix_index(x):
    end = x.split(':')[-1]
    new = x.replace(f":{end}", f"-{end}")
    return new

        
def getOpts():
    parser = argparse.ArgumentParser(description='Select ATAC samples with min 10 or 5 counts in clusters')
    parser.add_argument('--cluster-info', required=True, help ="""input cluster info file""")
    parser.add_argument('--cohort', default="FUSION", help ="""cohort to select""")
    parser.add_argument('--min-nuclei', type=int, default=10, help ="""min nuclei in sample to select sample for that cluster""")
    parser.add_argument('--exception-clusters', nargs='+', help ="""Clusters with a different min nuclei threhshold that --min-nuclei""")
    parser.add_argument('--exception-min-nuclei', type=int, help ="""min nuclei in sample to select sample for exception clusters""")
    parser.add_argument('--modality', default="atac", help ="""modality to select""")
    parser.add_argument('--prefix', default="atac", help ="""prefix""")
    args = parser.parse_args()
    return args


if __name__ == '__main__':

    args = getOpts()

    prefix = args.prefix
    d = pandas.read_csv(args.cluster_info, sep='\t', usecols=['index', 'sample', 'SNG.1ST', 'cohort', 'modality', 'coarse_cluster_name'], dtype={'SNG.1ST': str})
    d = d[(d['cohort']==args.cohort) & (d['modality']==args.modality)]
    
    dx = d[d['SNG.1ST'].isin(['12004', '22011', '32071'])]
    
    if dx.empty:
        print("No need to worry about the non Finnish sample and the related samples")
    else:
        d = d[~ d['SNG.1ST'].isin(['12004', '22011', '32071'])]

    d1 = d.groupby(['coarse_cluster_name', 'sample', 'SNG.1ST']).size().to_frame(name="n_nuclei").reset_index()

    min_nuclei = args.min_nuclei
    exception_clusters = args.exception_clusters # ['Neuromuscular_junction', 'Adipocyte']
    exception_min_nuclei = args.exception_min_nuclei

    print(d1.tail())

    if exception_clusters is not None and exception_min_nuclei is not None:
        d1 = d1[((d1['coarse_cluster_name'].isin(exception_clusters)) & (d1['n_nuclei'] >= exception_min_nuclei)) | ( (~d1['coarse_cluster_name'].isin(exception_clusters)) & (d1['n_nuclei']>= min_nuclei))]
    else:
        d1 = d1[d1['n_nuclei']>= min_nuclei]
    assert not d1.empty
    
    ds = d1.groupby('coarse_cluster_name').agg({
        'n_nuclei': [min, numpy.mean, numpy.median],
        'sample': len}).reset_index()
    
    ds.columns = ["_".join(c) if "" not in c else "".join(c) for c in ds.columns]
    ds.sort_values("n_nuclei_mean", inplace=True)

    assert not ds.empty
    ds.to_csv(f"{prefix}.cluster_selected_sample_counts.tsv", sep='\t', index=False)


    for cluster, grp in d1.groupby("coarse_cluster_name"):
        grp[['sample']].to_csv(f"{prefix}.{cluster}.sample_list.txt", sep='\t', index=False)
        grp[['SNG.1ST']].to_csv(f"{prefix}.{cluster}.sample_id_list.txt", sep='\t', index=False)


    

