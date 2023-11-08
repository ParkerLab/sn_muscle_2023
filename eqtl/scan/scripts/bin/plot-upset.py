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


def getOpts():
    parser = argparse.ArgumentParser(description='make upset plot for eQTL')
    parser.add_argument('--sig-dir', help ="""dir with list of significant eGenes""")
    
    args = parser.parse_args()
    return args



args = getOpts()

"""Fix the directory where all figures will be saved"""
fdir = "."
pe = pandas_extra.ExtraFunctions(fdir)



# plot upset
eqtls = glob.glob(f"{args.sig_dir}/*.significant-features")

def get(f):
    name = os.path.basename(f).split("--")[0].replace("fusion.", "")
    d = pandas.read_csv(f, sep='\t', header=None, names=['gene_name'])
    d['cluster'] = name
    return d

df = pandas.concat([get(f) for f in eqtls])
df['sig'] = True
df['cluster'] = df['cluster'].str.replace('Mesenchymal_Stem_Cell', 'Fibro-adipogenic\nprogenitors')
dn = pandas.pivot_table(data=df, index="gene_name", columns="cluster", values="sig").fillna(False)

# facecolor=["#a6cee3", "#1f78b4", "#8dd3c7", "#6a3d9a","#fb9a99"]
# facecolor = [tuple(int(h.lstrip("#")[i:i+2], 16) for i in (0, 2, 4)) for h in facecolor]
# print(facecolor)
def plotupset(dn, cplot, filename):
    dplot = dn[cplot]
    dplot = dplot.loc[dplot.any(axis=1),:]
    dnx = dplot.groupby(cplot).size()
    plot(dnx, sort_by='cardinality', show_counts='%d', min_subset_size=30)
    pe.saveb(filename)

# top 5
cplot = ['Endothelial', 'Fibro-adipogenic\nprogenitors', 'Type_1', 'Type_2a', 'Type_2x', 'aggregate']
plotupset(dn, cplot, "fig.upset_top5-with-agg.png")

# top5 no aggregate
cplot = ['Endothelial', 'Fibro-adipogenic\nprogenitors', 'Type_1', 'Type_2a', 'Type_2x']
plotupset(dn, cplot, "fig.upset_top5.png")

# fibers
cplot = ['Muscle_Fiber_Mixed', 'Neuromuscular_junction', 'Satellite_Cell', 'Type_1', 'Type_2a', 'Type_2x', 'aggregate']
plotupset(dn, cplot, "fig.upset_fibers-with-agg.png")

# others
cplot = ['Endothelial', 'Fibro-adipogenic\nprogenitors', 'Satellite_Cell', 'Smooth_Muscle', 'aggregate']
plotupset(dn, cplot, "fig.upset_others-with-agg.png")
