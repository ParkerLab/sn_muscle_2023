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
parser.add_argument("--eqtl-glob", required=True, nargs='+', help="""[Required]. egene-chunk mapping glob""")
parser.add_argument("--caqtl-glob", required=True, nargs='+', help="""[Required]. capeak-chunk mapping glob""")
args = parser.parse_args()

# for e-ca pairs that coloc, prepare to run CIT. Use signal hits from SuSiE
usecols = ['hit1', 'hit2', 'PP.H4.abf', 'cluster', 'egene', 'capeak']
c = pandas.read_csv(args.coloc, sep='\t', usecols = usecols).rename(columns = {'hit1': "eqhit",
                                                                               'hit2': "cahit"})
c = c[(c['cluster']!="aggregate") & (c['PP.H4.abf']>0.5)]
c['eqhit'] = c['eqhit'].str.replace("chrchr", "chr")
c['cahit'] = c['cahit'].str.replace("chrchr", "chr")
c['chrom'] = c['cahit'].map(lambda x: x.split("_")[0])
c = c[['egene', 'capeak', 'chrom', 'cluster', 'cahit', 'eqhit']].drop_duplicates()
c = c.groupby(['chrom', 'cluster', 'egene', 'capeak']).agg({'eqhit': lambda x: ",".join(set(x)), 'cahit': lambda x: ",".join(set(x))}).reset_index()

# cis indexed files for egene and capeak
cafiles = pandas.concat([pandas.read_csv(f, header=None, names=['capeak']).assign(capath=(os.path.abspath(f)).replace('.txt', '.bed.gz').replace("cafiles", "cis-indexed")) for f in args.caqtl_glob])
eqfiles = pandas.concat([pandas.read_csv(f, header=None, names=['egene']).assign(eqpath=(os.path.abspath(f)).replace('.txt', '.bed.gz').replace("eqfiles", "cis-indexed")) for f in args.eqtl_glob])

cafiles['cluster'] = cafiles['capath'].map(lambda x: os.path.basename(x).split("--")[0].replace("fusion.", ""))
eqfiles['cluster'] = eqfiles['eqpath'].map(lambda x: os.path.basename(x).split("--")[0].replace("fusion.", ""))


c = pandas.merge(pandas.merge(c, cafiles, how="inner"), eqfiles, how="inner")

for (chrom, cluster, capath, eqpath), grp in c.groupby(['chrom', 'cluster', 'capath', 'eqpath']):
    cachunk = os.path.basename(capath).replace(".bed.gz", "").split("--")[-1]
    eqchunk = os.path.basename(eqpath).replace(".bed.gz", "").split("--")[-1]
    fname = f"{chrom}@fusion.{cluster}@{cachunk}@{eqchunk}.mrs.tsv"
    grp.to_csv(fname, sep='\t', index=False)

# for cluster, grp in c.groupby("cluster"):
#     n = len(grp.index)
#     starts = list(numpy.arange(0, n, args.chunk))
#     for s in starts:
#         grp.iloc[s:s+args.chunk].to_csv(f"fusion.{cluster}--{s}.cit.tsv", sep='\t', index=False)
