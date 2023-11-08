#!/usr/bin/env python
# coding: utf-8

# In[3]:


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
    parser = argparse.ArgumentParser(description='Plot genotype PCs')
    parser.add_argument('--pc', required=True, help ="""PCs from QTLtools""")
    parser.add_argument('--cluster-info', required=True, help ="""input cluster info file""")
    parser.add_argument('--output', required=True, help ="""output file""")

    args = parser.parse_args()
    return args

def pls(df, p1, ax):
    p2 = p1 + 1
    a = f"scaled-centered-genotype_1_1_svd_PC{p1}"
    b = f"scaled-centered-genotype_1_1_svd_PC{p2}"
    sns.scatterplot(data=df, x=a, y=b, hue="sex", ax=ax, legend=False, **{'alpha': 0.5, 'linewidth':0})
    ax.set_title(f"PC {p1} vs PC {p2}")
    ax.set_xlabel("")
    ax.set_ylabel("")

if __name__ == '__main__':

    args = getOpts()

    d = pandas.read_csv(args.pc, sep=' ', index_col=[0]).T.reset_index()
    d.head()


    c = pandas.read_csv(args.cluster_info, sep='\t', usecols=['SNG.1ST', 'age', 'sex'], dtype={'SNG.1ST': str}).drop_duplicates()
    c.head()
    d = pandas.merge(d, c, how="inner", left_on="index", right_on="SNG.1ST")




    fig, axs = plt.subplots(2, 2)  #plt.figure(figsize=(6,6))

    pls(d, 1, axs[0,0])
    pls(d, 2, axs[0,1])
    pls(d, 3, axs[1,0])
    pls(d, 4, axs[1,1])

    for ax in fig.get_axes():
        ax.label_outer()

    plt.savefig(args.output, bbox_inches="tight", dpi=200)
