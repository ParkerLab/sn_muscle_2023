#!/usr/bin/env python
# coding: utf-8

# In[1]:

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
import argparse

def getOpts():
    parser = argparse.ArgumentParser(description='Plot PC scan results')
    parser.add_argument('--n-pheno-pcs', nargs='+', type=int, help ="""space-separated list of number of phenotype PCs tested""")
    parser.add_argument('--n-geno-pcs', nargs='+', type=int, help ="""space-separated list of number of genotype PCs tested""")
    parser.add_argument('--other-covs', nargs='+', help ="""space-separated list of number of other covariates tested""")
    parser.add_argument('--scan-result', required=True, help ="""Number of significant QTL scan results from QTLtools""")
    parser.add_argument('--prefix', required=True, help ="""output file prefix""")

    args = parser.parse_args()
    return args


if __name__ == '__main__':

    args = getOpts()

    n_pheno_pcs = args.n_pheno_pcs
    n_geno_pcs= args.n_geno_pcs
    other_covs: args.other_covs
    covlist = ['vars+geno-3'] + [f"geno-3+pheno-{n}" for n in n_pheno_pcs] + [f"vars+geno-3+pheno-{n}" for n in n_pheno_pcs]

    covmap = {f"covlist.{k}": cov for k, cov in enumerate(covlist)}

    res = args.scan_result
    o = pandas.read_csv(res, sep='\t')
    o[['cluster', 'cov']] = o['name'].str.split("--", expand=True)
    o.cluster = o.cluster.str.replace("fusion.", "")
    o['cov'] = o['cov'].map(covmap)
    o['nPCs'] = o['cov'].map(lambda x: int(x.split('-')[-1]) if "pheno" in x else 0)
    o['covtype'] = o['cov'].map(lambda x: '-'.join(x.split('-')[0:-1]))
    o.head()
    


    idx = o.groupby(['cluster'])['nsig'].transform(max) == o['nsig']
    o1 = o[idx].sort_values('nsig', ascending=False).drop_duplicates("cluster")
    order = o1['cluster'].tolist()

    o1.to_csv(f"{args.prefix}.max_signals.tsv", sep='\t', index=False)

    # In[8]:
    o.sort_values('covtype', inplace=True)
    g = sns.relplot(data=o, x="nPCs", col="cluster", col_wrap=5, height=2.5, aspect=1,
                    hue="covtype", kind="line", y="nsig", facet_kws={'sharey': False},
                    col_order=order)
    
    for clust, ax in zip(order, g.axes):
        sns.scatterplot(data=o[o['cluster']==clust], x="nPCs", y="nsig", ax=ax, hue="covtype", legend=False)
        x = o1[o1['cluster']==clust].iloc[0]['nPCs']
        y = o1[o1['cluster']==clust].iloc[0]['nsig']
        ax.text(x, y, y)
        #ax.plot([x], [y], linestyle="dashed", color="black")

    plt.savefig(f"{args.prefix}.max_nsig_pcs.png", bbox_inches="tight", dpi=200)
