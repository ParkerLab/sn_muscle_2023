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
parser.add_argument("--coloc", required=True, type=str, help="""[Required]. coloc output dir, will glob all files ending with .coloc.tsv""")
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
d[['cluster', 'egene', 'capeak']] = d['n'].str.split("--", expand=True)

print(f"total runs = {len(d.index)}")

d = d[(d['hit1.in.eqtl']=="yes") & (d['hit2.in.caqtl']=="yes")]
print(f"total tests = {len(d.index)}")

d = d[d['PP.H4.abf']>0.5]
if d.empty:
    print("No colocs. Exiting")
    sys.exit(0)
    
print(f"total coloc'd = {len(d.index)}")
print(f"total egenes coloc'd = {len(d['egene'].drop_duplicates())}")
d.head()

d.to_csv("summary.tsv", sep='\t', index=False)

d = d[d['cluster']!="aggregate"]
d1 = d[['cluster', 'egene', 'capeak', 'PP.H4.abf']].sort_values('PP.H4.abf', ascending=False).drop_duplicates(subset=['cluster', 'egene'])
d1.head()



corder = ['Type_1', 'Type_2a', 'Type_2x', 'Mesenchymal_Stem_Cell', 'Endothelial']
corder = [c for c in corder if c in d1['cluster'].drop_duplicates().tolist()]
print(f"Found colocs for {len(corder)} clusters")
dp = d1.pivot_table(index="cluster", columns = "egene", values = "PP.H4.abf")

go = dp.T.fillna(0).reset_index().set_index("egene")
go['n'] = go.astype(bool).sum(axis=1)
#print(go['n'].drop_duplicates())
print(go.head())
go.sort_values(['n']+corder, ascending=False, inplace=True)
lorder = go.index.tolist()
go.head()

dp = dp.loc[corder, lorder]

total = len(lorder)
specific = len(go[go['n'] == 1])
perc = round(float(specific/total*100), 2)

matplotlib.rcParams.update({'font.size': 12})
xlen = 40 if total > 10 else total * 2
plt.figure(figsize=(xlen, 6))
g = sns.heatmap(data=dp, cmap="viridis", vmin = 0.5, vmax = 1,
                cbar_kws={'label': 'PP shared causal variant\n(Coloc PP H4)',
                          "shrink": 0.8})
plt.xlabel(f"N eGenes\n(N = {total}); Specific = {specific} ({perc}%)", size=22)
plt.ylabel("Cluster", size=22)
g.set_yticklabels(g.get_yticklabels(), size = 20)
pe.saveb(f"fig.caqtl-eqtl.coloc.png")

    


# In[ ]:




