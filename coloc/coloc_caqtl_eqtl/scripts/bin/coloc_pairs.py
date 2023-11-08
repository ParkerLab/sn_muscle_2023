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
parser.add_argument("--eqtl", required=True, type=str, help="""[Required]. Glob for eqtl SuSiE""")
parser.add_argument("--gene-info", required=True, type=str, help="""[Required]. Gene info for TSS""")
parser.add_argument("--caqtl", required=True, help="""[Required]. Glob for caqtl SuSiE""")
args = parser.parse_args()

"""Fix the directory where all figures will be saved"""

fdir = "."
if not os.path.exists(fdir):
    os.mkdir(fdir)
    
pe = pandas_extra.ExtraFunctions(fdir)

g = pandas.read_csv(args.gene_info, sep='\t')
g['end'] = g.apply(lambda x: x['s'] if x['strand']=='+' else x['e'], axis=1)
g['start'] = g['end'] - 1 - 250000
g['end'] = g['end'] + 250000
g['start'] = numpy.where(g['start']<0, 0, g['start'])

# loci to test coloc for
elist = sorted(glob.glob(f"{args.eqtl}/*/*.Rda"))
assert len(elist) > 0

def get_gene_name(f):
    s = os.path.basename(f).replace(".susie.Rda", "").split(".")
    gene = ".".join(s[5:])
    return gene

fe = pandas.DataFrame([[os.path.basename(f).split("--")[0].replace("fusion.", ""), get_gene_name(f), os.path.abspath(f)] for f in elist], columns = ["cluster", "gene_name", "fe"])
fe = pandas.merge(fe, g, how="inner", on="gene_name")
fe.rename(columns = {'gene_name': "egene"}, inplace=True)
febed = pybedtools.BedTool().from_dataframe(fe[['chrom', 'start', 'end', 'cluster', 'egene', 'fe']])

calist = sorted(glob.glob(f"{args.caqtl}/*/*.Rda"))
assert len(calist) > 0
fc = pandas.DataFrame([[os.path.basename(f).split("--")[0].replace("fusion.", ""), os.path.basename(f).split(".")[5], os.path.abspath(f)] for f in calist], columns = ["cluster", "gene_name_caqtl", "fc"])
fc[['chrom', 'start', 'end']] = fc['gene_name_caqtl'].str.split(":", expand=True)
cabed = pybedtools.BedTool().from_dataframe(fc[['chrom', 'start', 'end', 'cluster', 'gene_name_caqtl', 'fc']])


o = febed.intersect(cabed, wa=True, wb=True).to_dataframe(names=['c', 'es', 'ee', 'cluster', 'egene', 'fe', 'c2', 'cs', 'ce', 'cluster1', 'capeak', 'fc'])
o = o[o['cluster'] == o['cluster1']][['cluster', 'egene', 'fe', 'capeak', 'fc']]
o.head()

for i, grp in o.groupby(['cluster', 'egene']):
    filename = f"{i[0]}--{i[1]}.coloc.sh"
    cmd = f"coloc.R --eqtl {grp.iloc[0]['fe']} --caqtl {','.join(grp['fc'].tolist())} --prefix {i[0]}--{i[1]} \n"
    with open(filename, 'w') as f:
        f.write("#!/bin/bash\n\n")
        f.write(cmd)
        




