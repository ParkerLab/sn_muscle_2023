#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas
import numpy
import sys
import os
import glob
import argparse

parser = argparse.ArgumentParser(description="Select pairs to test CIT ")
parser.add_argument("--coloc", required=True, type=str, help="""[Required]. coloc output summary df""")
parser.add_argument("--chunk", default=100, type=int, help="""[Required]. how many rows per output""")
args = parser.parse_args()

# for e-ca pairs that coloc, prepare to run CIT. Use signal hits from SuSiE
usecols = ['hit1', 'hit2', 'PP.H4.abf', 'cluster', 'egene', 'capeak']
c = pandas.read_csv(args.coloc, sep='\t', usecols = usecols).rename(columns = {'hit1': "eqhit",
                                                                               'hit2': "cahit"})
c = c[(c['cluster']!="aggregate") & (c['PP.H4.abf']>0.5)]
c['eqhit'] = c['eqhit'].str.replace("chrchr", "chr")
c['cahit'] = c['cahit'].str.replace("chrchr", "chr")
c = c.groupby(['cluster', 'egene', 'capeak']).agg({'eqhit': lambda x: ",".join(set(x)), 'cahit': lambda x: ",".join(set(x))}).reset_index()

for cluster, grp in c.groupby("cluster"):
    n = len(grp.index)
    starts = list(numpy.arange(0, n, args.chunk))
    for s in starts:
        grp.iloc[s:s+args.chunk].to_csv(f"fusion.{cluster}--{s}.cit.tsv", sep='\t', index=False)
