#!/usr/bin/env python
# coding: utf-8

# In[18]:


import pandas
import numpy
import math
import sys
import os
import glob
import argparse

def getOpts():
    parser = argparse.ArgumentParser(description='compile caQTL results across scans')
    parser.add_argument('--info', required=True, help ="""cluster info to calculate median nuclei per sample""")
    parser.add_argument('--sig', required=True, help ="""significant caQTL across scans""")
    parser.add_argument('--tpm', required=True, help ="""directory with *.bed.gz TPM matrices to count total genes tested""")
    parser.add_argument('--select', required=True, help ="""selected scan eg __top_280_samples___covlist.4""")
    parser.add_argument('--output',  help ="""Output prefix.""")
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    
    args = getOpts()

    d = pandas.read_csv(args.info,  sep='\t', dtype={'SNG.1ST': str})
    d = d[(d['cohort']=="FUSION") & (d['modality']=="atac")]
    print(d.shape)


    # median nuclei per clust, sampl
    d1 = d.groupby(['SNG.1ST', 'coarse_cluster_name']).size().to_frame(name="N").reset_index()
    d2 = d1.groupby('coarse_cluster_name').agg({'N':numpy.median}).reset_index()


    # tested genes
    dlist= []
    for f in glob.glob(f"{args.tpm}*bed.gz"):
        n_tested = len(pandas.read_csv(f, sep='\t', usecols=['gid']).index)
        cluster = os.path.basename(f).replace(".bed.gz", "")
        dlist.append([cluster, n_tested])
        
    d3 = pandas.DataFrame(dlist, columns=['coarse_cluster_name', 'n_tested'])

    # 10% FDR sig
    x = pandas.read_csv(args.sig, sep='\t')
    x = x[x['name'].str.contains(args.select)]
    x['coarse_cluster_name'] = x['name'].str.replace(args.select, "")
    x = x[['coarse_cluster_name', 'nsig']]
    
    x1 = pandas.merge(pandas.merge(d2, d3, how="inner", on="coarse_cluster_name"), x,  how="inner", on="coarse_cluster_name")
    x1['percent'] = x1['nsig']/x1['n_tested']*100
    
    x1.to_csv(args.output, sep='\t', index=False, na_rep="NA")




