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
import pandas_extra
import argparse

# In[2]:
parser = argparse.ArgumentParser(description="Select pairs to test coloc ")
parser.add_argument("--coloc", required=True, type=str, help="""[Required]. coloc output directory, will glob all files ending with .coloc.tsv""")
parser.add_argument("--pheno", type=str, help="""[Required]. UKBB phenotype description df""")
args = parser.parse_args()

"""Fix the directory where all figures will be saved"""

fdir = "."
if not os.path.exists(fdir):
    os.mkdir(fdir)
    
pe = pandas_extra.ExtraFunctions(fdir)



def getf(f):
    try:
        d = pandas.read_csv(f, sep='\t')
        n = os.path.basename(f).replace(".coloc.tsv", "")
        d['n'] = n
    except pandas.errors.EmptyDataError:
        d = pandas.DataFrame()
        
    return d    
    
d = pandas.concat([getf(f) for f in glob.glob(f"{args.coloc}/*.coloc.tsv")])
d[['gwas', 'eqtl']] = d['n'].str.replace(".Rda", "").str.split("__fusion.", expand=True)
d['cluster']= d['eqtl'].map(lambda x: x.split("--")[0].replace("fusion.", ""))

print(f"total tested = {len(d.index)}")
print(f"total coloc'd = {len(d[d['PP.H4.abf']>0.5].index)}")
print(f"total GWAS loci coloc'd = {len(d[d['PP.H4.abf']>0.5]['gwas'].drop_duplicates())}")

def getl(x):
    s= x.split("__")
    out = f"{s[1]}__{s[2]}__{s[3]}"
    return out

d['trait'] = d['gwas'].map(lambda x: x.split("__")[0])
d['gwas_locus'] = d['gwas'].map(getl)
p = pandas.read_csv(args.pheno, sep='\t')

def replacena(x):
    try:
        if math.isna(float(x['traitname'])):
            out = x['trait']
        else:
            out = x['traitname']
    except:
        out = x['traitname']
    return out

d = pandas.merge(d, p, how="left")
d['traitname'] = d.apply(lambda x: replacena(x), axis=1) 
d.to_csv("summary.tsv", sep='\t', index=False)

# only plot where the coloc-reported hits were in the actual credible set
d = d[(d['hit1.in.gwas']=="yes") & (d['hit2.in.qtl']=="yes")]

def heatmap_filtered(d, selected, trait, minppa = 0.5):
    print(selected)
    grp = d[d['PP.H4.abf'] >= minppa].sort_values('PP.H4.abf', ascending=False).drop_duplicates(subset=['trait', 'hit1', 'cluster'])
    corder = [c for c in selected if c in grp['cluster'].drop_duplicates().tolist()]
    print(f"Found colocs for {len(corder)} clusters")
        
    if len(corder) > 0:
        dp = grp[grp['cluster'].isin(corder)].pivot_table(index="cluster", columns = "hit1", values = "PP.H4.abf")

        # get order of loci to show
        go = dp.T.fillna(0).reset_index().set_index("hit1")
        go['n'] = go.astype(bool).sum(axis=1)
        go.sort_values(['n']+corder, ascending=False, inplace=True)
        lorder = go.index.tolist()

        dp = dp.loc[corder, lorder]
        total = len(lorder)
        specific = len(go[go['n'] == 1])
        perc = round(float(specific/total*100), 2)

        xlen = 40 if total > 10 else total * 2
        plt.figure(figsize=(xlen, 6))
        g = sns.heatmap(data=dp, cmap="viridis", vmin = minppa, vmax = 1, xticklabels=False,
                        cbar_kws={'label': 'PP shared causal variant\n(Coloc PP H4)',
                                  "shrink": 0.8})

        plt.xlabel(f"{trait} GWAS signals\n(N = {total}); Specific = {specific} ({perc}%)", size=22)
        plt.ylabel("Cluster", size=22)
        g.set_yticklabels(g.get_yticklabels(), size = 20)
        g.set_xticks([])
        g.set_xticklabels([])
        pe.saveb(f"fig.{trait.replace('/', '_').replace(' ', '_')}.filtered-ppa-{minppa}.viridis.png")
        return lorder
    else:
        return []

    
def heatmap_full(grp, selected, trait, lorder=None):
    corder = [c for c in selected if c in grp['cluster'].drop_duplicates().tolist()]
    print(f"Found colocs for {len(corder)} clusters")
    if len(corder) > 0:
        dp = grp[grp['cluster'].isin(corder)].pivot_table(index="cluster", columns = "hit1", values = "PP.H4.abf")
        dp = dp.loc[corder, lorder]
    
        total = len(lorder)
        matplotlib.rcParams.update({'font.size': 12})
        xlen = 40 if total > 10 else total * 2
        plt.figure(figsize=(xlen, 6))
        g = sns.heatmap(data=dp, cmap=sns.color_palette("Reds", as_cmap=True), vmin = 0.0, vmax = 1, xticklabels=False,
                        cbar_kws={'label': 'PP shared causal variant\n(Coloc PP H4)',
                                  "shrink": 0.8})
        plt.ylabel("")
        plt.xlabel("")
        g.set_yticklabels(g.get_yticklabels(), size = 20)
        g.set_xticks([])
        g.set_xticklabels([])

        pe.saveb(f"fig.{trait.replace('/', '_').replace(' ', '_')}.full-ppa-reds.png")
        return g
    else:
        return []
  


# In[17]:



rename = {"Mesenchymal_Stem_Cell": "FAP",
         "_": " "}
def renamec(df, col):
    for key in rename.keys():
        df[col] = df[col].str.replace(key, rename[key])
    return df
corder = ['Type 1', 'Type 2a', 'Type 2x', 'FAP', 'Endothelial']
d = renamec(d, "cluster")

for trait, grp in d.groupby("traitname"):
    lorder = heatmap_filtered(grp, corder, trait)
    if len(lorder) > 0:
        dx = grp.sort_values('PP.H4.abf', ascending=False).drop_duplicates(subset=['trait', 'hit1', 'cluster'])
        heatmap_full(dx, corder, trait, lorder=lorder)


     



