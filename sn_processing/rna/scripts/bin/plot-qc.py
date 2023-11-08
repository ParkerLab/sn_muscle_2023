#!/usr/bin/env python

import pandas
import numpy
import glob
import os
import argparse
import re
import sys
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import math
dpi=75
matplotlib.rcParams['figure.dpi']= dpi
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams.update({'font.size': 12})
from utils_customgrid import makeGrid
import pandas_extra


# In[41]:
def getOpts():
    parser = argparse.ArgumentParser(description='Make RNA qc plot for channels')
    parser.add_argument('--input', required=True, nargs='+', help ="""QC input files.""")
    parser.add_argument('--data', required=True, help ="""Output data file name.""")
    parser.add_argument('--fig-dir', default='.', help ="""Directory to output figures.""")
    parser.add_argument('--facet-wrap', type=int, default=4, help ="""Facet wrap""")
    args = parser.parse_args()
    return args


if __name__ == '__main__':

    args = getOpts()

    pe = pandas_extra.ExtraFunctions(args.fig_dir)
    def read(f):
        df = pandas.read_csv(f, sep='\t')
        library = os.path.basename(f).replace("-hg19.txt", "").replace("-hg38.txt", "")
        df['library'] = library
        return df

    dlist = [read(f) for f in args.input]
    d = pandas.concat(dlist, ignore_index=True)
    d = d[d['number_umis']>100]
    d['log_umis'] = d['number_umis'].map(lambda x: math.log(x, 10))
    d['log_mito'] = d['fraction_mitochondrial'].map(lambda x: math.log(x + 0.0001, 10))
    d.sort_values('library', inplace=True)
    d.to_csv(args.data, sep='\t', index=False, na_rep="NA")

    libraries = [os.path.basename(f).replace("-hg19.txt", "").replace("-hg38.txt", "") for f in args.input]
    libraries.sort()
    d['row'] = d['library'].map(lambda x: "A" if x in libraries[0:args.facet_wrap] else "B")

    # N umis jointplot
    fig_width = args.facet_wrap * 3
    fig_height = (len(libraries)/args.facet_wrap) * 3

    g = makeGrid(d, row='row', col='library', x="log_umis", y="log_mito", kind="scatter", figsize=(fig_width, fig_height),
                 figname="fig.umis_joint.png", facet_wrap=args.facet_wrap, **{'main_x_label':"log10(# UMIs)", 'main_y_label':"log10(Mitochondrial fraction)", 'y_label':"",
                                                  'line': [[3.5, 3.5], [-5, 0]] })

    # N umis cumulative fraction
    plt.figure(figsize=(4,4))
    for i, grp in d.groupby('library'):
        g = pe.get_cdf_column(grp, 'number_umis', label=i)
        plt.xscale('log')
    plt.xlabel("# UMIs")
    plt.ylabel("Cumulative fraction of droplets")
    plt.legend()
    pe.legend_out()
    plt.savefig("fig.umis_cumulative_fraction.pdf", bbox_inches="tight")
    
    # Number of nuclei with  UMIs greater than:

    plt.figure(figsize=(4,4))
    for i, grp in d.groupby('library'):
        pe.get_rows_greater_than(grp, 'number_umis', start=1000, step=100, end=10000, label=i)
    plt.legend()
    plt.xlabel("X")
    plt.ylabel("# droplets with UMIx >= X")
    pe.legend_out()
    plt.savefig("fig.n_nuclei_umis.pdf", bbox_inches="tight")
    
    
    
    
