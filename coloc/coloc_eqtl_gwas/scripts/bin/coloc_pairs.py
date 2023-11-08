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
import pandas_extra
import argparse


# In[2]:

parser = argparse.ArgumentParser(description="Select pairs to test coloc ")
parser.add_argument("--eqtl", required=True, nargs='+', help="""[Required]. SuSiE output directories for eqtl""")
parser.add_argument("--gene-info", required=True, type=str, help="""[Required]. Gene info for TSS""")
parser.add_argument("--gwas", required=True, nargs='+', help="""[Required]. SuSiE output directories for GWAS""")
args = parser.parse_args()

"""Fix the directory where all figures will be saved"""

fdir = "."
pe = pandas_extra.ExtraFunctions(fdir)

g = pandas.read_csv(args.gene_info, sep='\t')
g['end'] = g.apply(lambda x: x['s'] if x['strand']=='+' else x['e'], axis=1)
g['start'] = g['end'] - 1 - 250000
g['end'] = g['end'] + 250000
g['start'] = numpy.where(g['start']<0, 0, g['start'])

# loci to test coloc for
fl = []
for  dirname in args.eqtl:
    fl.extend(glob.glob(f"{dirname}/*/*.Rda"))
flist = sorted(fl)
assert len(flist) > 0
def get_gene_name(f):
    s = os.path.basename(f).replace(".susie.Rda", "").split(".")
    gene = ".".join(s[5:])
    return gene

f = pandas.DataFrame([[get_gene_name(f), os.path.abspath(f)] for f in flist], columns = ["gene_name", "f"])
f = pandas.merge(f, g, how="inner", on="gene_name")
fbed = pybedtools.BedTool().from_dataframe(f[['chrom', 'start', 'end', 'f']])


fl = []
for  dirname in args.gwas:
    fl.extend(glob.glob(f"{dirname}/*.results.tsv"))
dlist = sorted(fl)

assert len(dlist) > 0

def getf(f):
    d = pandas.read_csv(f, sep='\t', usecols = ['nsets'])
    d = d[d['nsets']>0]
    if not d.empty:
        region = os.path.basename(f).split("__")[4]
        path = os.path.abspath(f).replace(".results.tsv",".Rda")
        d['region_gw'] = region
        d['f'] = path
    return d

d = pandas.concat([getf(f) for f in dlist])
d[['chrom', 'start', 'end']] = d['region_gw'].str.split("-", expand=True)
dbed = pybedtools.BedTool().from_dataframe(d[['chrom', 'start', 'end', 'f']])


o = dbed.intersect(fbed, wa=True, wb=True).to_dataframe(names=['c', 'gs', 'ge', 'gwas', 'c2', 'cs', 'ce', 'eqtl'])
o['region'] = o.apply(lambda x: f"{x['c']}:{x['gs']}-{x['ge']}", axis=1)


for i, grp in o.groupby("region"):
    gwaspath = grp.iloc[0]['gwas']
    gwasname = os.path.basename(gwaspath).split("__")[0]
    filename = f"{gwasname}.{i}.coloc.sh"
    cmd = f"coloc.R --gwas {gwaspath} --eqtl {','.join(grp['eqtl'].tolist())} --ukbb {gwasname}.{i}.rsid \n"
    with open(filename, 'w') as f:
        f.write("#!/bin/bash\n\n")
        f.write(cmd)        




