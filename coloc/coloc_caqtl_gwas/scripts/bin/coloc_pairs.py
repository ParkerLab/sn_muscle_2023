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
parser.add_argument("--caqtl", required=True, type=str, help="""[Required]. Glob for caqtl SuSiE""")
parser.add_argument("--gwas", required=True, nargs='+', help="""[Required]. Glob for GWAS SuSiE""")
args = parser.parse_args()

"""Fix the directory where all figures will be saved"""

fdir = "."
if not os.path.exists(fdir):
    os.mkdir(fdir)
    
pe = pandas_extra.ExtraFunctions(fdir)


# loci to test coloc for
flist = sorted(glob.glob(f"{args.caqtl}/*/*.Rda"))
assert len(flist) > 0
# don't run for aggregate caQTL
flist = [f for f in flist if "aggregate" not in f]  
f = pandas.DataFrame([[os.path.basename(f).split(".")[5], os.path.abspath(f)] for f in flist], columns = ["region_ca", "f"])
f[['chrom', 'start', 'end']] = f['region_ca'].str.split(":", expand=True)
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


o = dbed.intersect(fbed, wa=True, wb=True).to_dataframe(names=['c', 'gs', 'ge', 'gwas', 'c2', 'cs', 'ce', 'caqtl'])
o['region'] = o.apply(lambda x: f"{x['c']}:{x['gs']}-{x['ge']}", axis=1)
# o['gwas'] = o['gwas'].str.replace("/lab/work/", "/gpfs/accounts/scjp_root/scjp99/")
# o['caqtl'] = o['caqtl'].str.replace("/lab/work/", "/gpfs/accounts/scjp_root/scjp99/")

o1 = o.groupby(['gwas']).size().to_frame(name="N").reset_index()
o1.sort_values('N', ascending=False)

    
for i, grp in o.groupby("region"):
    gwaspath = grp.iloc[0]['gwas']
    gwasname = os.path.basename(gwaspath).split("__")[0]
    filename = f"{gwasname}.{i}.coloc.sh"
    cmd = f"coloc.R --gwas {gwaspath} --caqtl {','.join(grp['caqtl'].tolist())} --ukbb {gwasname}.{i}.rsid \n"
    with open(filename, 'w') as f:
        f.write("#!/bin/bash\n\n")
        f.write(cmd)
        




