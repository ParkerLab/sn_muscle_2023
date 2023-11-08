#!/usr/bin/env python
# coding: utf-8

# In[ ]:


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
import subprocess as sp
dpi = 150
matplotlib.rcParams['figure.dpi']= dpi
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams.update({'font.size': 12})
import pandas_extra
from scipy import stats
from statsmodels.stats.multitest import multipletests
import argparse

def getOpts():
    parser = argparse.ArgumentParser(description='Plot PC scan results')
    parser.add_argument('--max-signals', help ="""n pcs that maximize signals file""")
    parser.add_argument('--permute-dir', help ="""formatted permutation scan output dir""")
    parser.add_argument('--matrix-dir', help ="""atac counts matrix dir""")
    parser.add_argument('--sample-dir', help ="""sample list file dir""")
    args = parser.parse_args()
    return args


args = getOpts()


"""Fix the directory where all figures will be saved"""

fdir = "."
if not os.path.exists(fdir):
    os.mkdir(fdir)

pe = pandas_extra.ExtraFunctions(fdir)
        
ms = pandas.read_csv(args.max_signals, sep='\t')
ms[['cluster', 'cov']] = ms['name'].str.split("--", expand=True)
ms.head()

flist = {cluster: f"{args.permute_dir}/{cluster}--{ms[ms['cluster']==cluster].iloc[0]['cov']}.permute.tsv" for cluster in ms['cluster'].drop_duplicates().tolist()}
flist

def getfeatures(m, minct, minsamp):
    m1 = m.copy()
    m1[m1<minct] = 0
    m1['nsamp'] = m1.astype(bool).sum(axis=1)
    features = m1[m1['nsamp']> minsamp].index.tolist()
    return features
    
def getsig(cluster):
    f = flist[cluster]
    d = pandas.read_csv(f, sep='\t')
    d.head()

    mfile = f"{args.matrix_dir}/atac.{cluster.replace('fusion.', '')}.tsv"
    sfile = f"{args.sample_dir}/{cluster}.sample_list.txt"

    m = pandas.read_csv(mfile, sep='\t', index_col=0)
    samples = pandas.read_csv(sfile, sep='\t', squeeze=True).tolist()

    selected_samples = [s for s in samples if s in m.columns]
    missing_samples = [s for s in samples if s not in m.columns]
    assert len(missing_samples) == 0
    m = m.loc[d['gene_name'].tolist(), selected_samples]
    assert not m.empty

    mincounts = [2,5,6,7,8,9,10,12,15]
    minsamples = [10]

    siglist = []
    for minct in mincounts:
        for minsamp in minsamples:
            features = getfeatures(m, minct, minsamp)
            ntested = len(features)
            d1 = d[d['gene_name'].isin(features)]
            d1['qval'] = multipletests(d1['p_beta'], method="fdr_bh")[1]
            nsig = len(d1[d1['qval']<0.05].index)
            siglist.append([cluster, minct, minsamp, ntested, nsig])
    print(siglist)
    return siglist


# In[6]:


out = [getsig(cluster) for cluster in flist.keys()]
out1 = [j for i in out for j in i]
out1


# In[7]:


df = pandas.DataFrame(out1, columns = ['cluster', 'mincount', 'minsample', 'ntested', 'nsig'])
df['cluster'] = df['cluster'].str.replace("fusion.", "")
df.to_csv("multiple_testing_titration.tsv", sep='\t', index=False)
df.head()


clusters = ['Type_1', 'Type_2a', 'Type_2x', 'Mesenchymal_Stem_Cell', 'Endothelial']
clusters2 = [c for c in df['cluster'].drop_duplicates() if c not in clusters]


# In[9]:


clist = clusters
ncols = len(clist)
fig, ax = plt.subplots(2, ncols, figsize=(ncols*3, 6))

for i, cluster in enumerate(clist):
    ax1 = ax[0, i]
    ax2 = ax[1, i]
    pdf = df[df['cluster']==cluster]
    sns.lineplot(data=pdf, x="mincount", y="ntested", ax = ax1)
    sns.scatterplot(data=pdf, x="mincount", y="ntested", ax = ax1)
    ax1.set_title(cluster)
    sns.lineplot(data=pdf, x="mincount", y="nsig", ax = ax2)
    sns.scatterplot(data=pdf, x="mincount", y="nsig", ax = ax2)


plt.tight_layout()
pe.save("fig.top5.nsig_ntested.png")


# In[10]:


clist = clusters2
ncols = len(clist)
fig, ax = plt.subplots(2, ncols, figsize=(ncols*3, 6))

for i, cluster in enumerate(clist):
    ax1 = ax[0, i]
    ax2 = ax[1, i]
    pdf = df[df['cluster']==cluster]
    sns.lineplot(data=pdf, x="mincount", y="ntested", ax = ax1)
    sns.scatterplot(data=pdf, x="mincount", y="ntested", ax = ax1)
    ax1.set_title(cluster)
    sns.lineplot(data=pdf, x="mincount", y="nsig", ax = ax2)
    sns.scatterplot(data=pdf, x="mincount", y="nsig", ax = ax2)


plt.tight_layout()
pe.save("fig.otherclusters.nsig_ntested.png")


# In[ ]:




