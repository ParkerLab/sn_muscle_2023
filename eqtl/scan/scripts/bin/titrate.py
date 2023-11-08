#!/usr/bin/env python
# coding: utf-8

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
import shutil
from upsetplot import plot
import argparse
from rpy2 import robjects as ro
ro.r("options(rlib_downstream_check = FALSE)")
from rpy2.robjects.packages import importr
qvalue = importr('qvalue')


def getOpts():
    parser = argparse.ArgumentParser(description='get optimal list of features to test eQTL')
    parser.add_argument('--max-signals', help ="""best nPCs""")
    parser.add_argument('--permute-dir', help ="""permute dir""")
    parser.add_argument('--rna-counts-dir', help ="""rna counts dir glob""")
    parser.add_argument('--sample-id-dir', help ="""sample id list glob""")
    
    args = parser.parse_args()
    return args



args = getOpts()

"""Fix the directory where all figures will be saved"""
fdir = "."
pe = pandas_extra.ExtraFunctions(fdir)

max_signals = args.max_signals
permute_dir = args.permute_dir
rna_counts_dir = args.rna_counts_dir
sample_id_dir = args.sample_id_dir

ms = pandas.read_csv(max_signals,  sep='\t')
ms[['cluster', 'cov']] = ms['name'].str.split("--", expand=True)
ms.head()

flist = {cluster: f"{permute_dir}/{cluster}--{ms[ms['cluster']==cluster].iloc[0]['cov']}.permute.tsv" for cluster in ms['cluster'].drop_duplicates().tolist()}
flist


# In[11]:


def getfeatures(m, minct, minsamp):
    m1 = m.copy()
    m1[m1<minct] = 0
    m1['nsamp'] = m1.astype(bool).sum(axis=1)
    features = m1[m1['nsamp']> minsamp].index.tolist()
    return features
    
def getsig(cluster):
    f = flist[cluster]
    covtype = os.path.basename(f).replace(".permute.tsv", "")
    d = pandas.read_csv(f, sep='\t')
    mfile = f"{rna_counts_dir}/{cluster}.tsv"
    sfile = f"{sample_id_dir}/{cluster}.sample_id_list.txt"
   
    m = pandas.read_csv(mfile, sep='\t').set_index("pid")
    m.drop(['#Chr', 'start', 'end', 'gid', 'strand', 'length'], axis=1, inplace=True)
    samples = pandas.read_csv(sfile, sep='\t', squeeze=True, dtype={'SNG.1ST': str}).tolist()

    selected_samples = [s for s in samples if s in m.columns]
    missing_samples = [s for s in samples if s not in m.columns]
    
    assert len(missing_samples) == 0
    m = m.loc[d['gene_name'].tolist(), selected_samples]
    assert not m.empty

    mincounts = [1,2,3,4,5,6,7,8,9,10,15]
    minsamples = [10, 20, 30]

    siglist = []
    
    for minct in mincounts:
        for minsamp in minsamples:
            features = getfeatures(m, minct, minsamp)
            filename = f"{covtype}__{minct}__{minsamp}.tested"
            with open(filename, 'w') as f1:
                f1.write("\n".join(features))
            ntested = len(features)
            d1 = d[d['gene_name'].isin(features)]
            try:
                p = ro.vectors.FloatVector(d1['p_beta'])
                d1.loc[:, 'qval'] = qvalue.qvalue(p)[2]
                sig = d1[d1['qval']<0.05]
                phenofilename = f"{covtype}__{minct}__{minsamp}.significant"
                with open(phenofilename, 'w') as f1:
                    f1.write("\n".join(sig['gene_name'].tolist()))
                nsig = len(sig.index)
                siglist.append([cluster, minct, minsamp, ntested, nsig, covtype])
            # try:
            #     d1.loc[:, 'qval'] = multipletests(d1['p_beta'], method="fdr_bh")[1]
            #     sig = d1[d1['qval']<0.05]
            #     phenofilename = f"{covtype}__{minct}__{minsamp}.significant"
            #     with open(phenofilename, 'w') as f1:
            #         f1.write("\n".join(sig['gene_name'].tolist()))
            #     nsig = len(sig.index)
            #     siglist.append([cluster, minct, minsamp, ntested, nsig, covtype])

            except:
                pass
            
    print(siglist)
    return siglist


out = [getsig(cluster) for cluster in flist.keys()]
out1 = [j for i in out for j in i]

df = pandas.DataFrame(out1, columns = ['cluster', 'mincount', 'minsample', 'ntested', 'nsig', 'covtype'])
# df['cluster'] = df['cluster'].str.replace("fusion.", "")
df.to_csv("multiple_testing_res_pbeta.tsv", sep='\t', index=False)
df.head()

def plot_titration(df, cluster_list, fname):
    ncols = len(cluster_list)
    fig, ax = plt.subplots(2, ncols, figsize=(ncols*3, 6))

    for i, cluster in enumerate(cluster_list):
        ax1 = ax[0, i]
        ax2 = ax[1, i]
        pdf = df[df['cluster']==cluster]
        pdf['minsample'] = pdf['minsample'].astype(int)
        print(pdf.head())
        sns.lineplot(data=pdf, x="mincount", y="ntested", hue="minsample", ax = ax1, legend=False)
        sns.scatterplot(data=pdf, x="mincount", y="ntested", hue="minsample", ax = ax1, legend=False)
        ax1.set_title(cluster)
        sns.lineplot(data=pdf, x="mincount", y="nsig", hue="minsample", ax = ax2, legend=False)
        sns.scatterplot(data=pdf, x="mincount", y="nsig", hue="minsample", ax = ax2, legend=False)


    plt.tight_layout()
    pe.save(fname)


clusters = ['fusion.Type_1', 'fusion.Type_2a', 'fusion.Type_2x', 'fusion.Mesenchymal_Stem_Cell', 'fusion.Endothelial', 'fusion.aggregate']
clusters2 = [c for c in df['cluster'].drop_duplicates() if c not in clusters]

plot_titration(df, clusters, "fig.top5.nsig_ntested.png")
plot_titration(df, clusters2, "fig.otherclusters.nsig_ntested.png")

# In[9]:
## get final list of tested features
df = df[df['minsample'] == 10]
df = df.sort_values("nsig", ascending=False).drop_duplicates("cluster")

for i, r in df.iterrows():
    shutil.copy(f"{r['covtype']}__{r['mincount']}__{r['minsample']}.tested", f"{r['covtype']}.tested-features")
    shutil.copy(f"{r['covtype']}__{r['mincount']}__{r['minsample']}.significant", f"{r['covtype']}.significant-features")


# plot upset
sig_features = glob.glob("*significant-features")
def get(f):
    name = os.path.basename(f).replace(".significant-features", "").replace("fusion.", "")
    d = pandas.read_csv(f, sep='\t', header=None, names=['gene_name'])
    d['cluster'] = name
    return d

df = pandas.concat([get(f) for f in sig_features])
df[['cluster', 'cov']] = df['cluster'].str.split("--", expand=True)
df['sig'] = True
dn = pandas.pivot_table(data=df, index="gene_name", columns="cluster", values="sig").fillna(False)

def plotupset(dn, cplot, filename):
    dplot = dn[cplot]
    dplot = dplot.loc[dplot.any(axis=1),:]
    dnx = dplot.groupby(cplot).size()
    plot(dnx, sort_by='cardinality', show_counts='%d')
    pe.saveb(filename)

# top 5
cplot = ['Endothelial', 'Mesenchymal_Stem_Cell', 'Type_1', 'Type_2a', 'Type_2x', 'aggregate']
plotupset(dn, cplot, "fig.upset_top5-with-agg.png")

# top5 no aggregate
cplot = ['Endothelial', 'Mesenchymal_Stem_Cell', 'Type_1', 'Type_2a', 'Type_2x']
plotupset(dn, cplot, "fig.upset_top5.png")

# fibers
cplot = ['Muscle_Fiber_Mixed', 'Neuromuscular_junction', 'Satellite_Cell', 'Type_1', 'Type_2a', 'Type_2x', 'aggregate']
plotupset(dn, cplot, "fig.upset_fibers-with-agg.png")

# others
cplot = ['Endothelial', 'Mesenchymal_Stem_Cell', 'Satellite_Cell', 'Smooth_Muscle', 'aggregate']
plotupset(dn, cplot, "fig.upset_others-with-agg.png")
