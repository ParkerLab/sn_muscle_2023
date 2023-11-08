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
dpi = 150
matplotlib.rcParams['figure.dpi']= dpi
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams.update({'font.size': 12})
from scipy import stats
from statsmodels.stats.multitest import multipletests
from scipy.stats import mannwhitneyu
import pandas_extra
import argparse

def getOpts():
    parser = argparse.ArgumentParser(description='compile eQTL results across scans')
    parser.add_argument('--qtl', required=True, nargs='+', help ="""Outputs from qtltools permute scan""")
    parser.add_argument('--threshold', required=True, type=float, help ="""qvalue significance threshold to count caQTL""")
    parser.add_argument('--qvalue', required=True, type=str, help ="""column name with the q values""")
    parser.add_argument('--output',  help ="""Output prefix.""")
    args = parser.parse_args()
    return args

def fixdf(f, fdrthresh, col):
    name = os.path.basename(f).replace(".permute.tsv", "")
    p = pandas.read_csv(f, sep='\t')
    n = len(p[p[col] < fdrthresh].index)
    return [name, n]

if __name__ == '__main__':
    
    args = getOpts()

    # concat results across chunks and calculate FDR
    d = pandas.DataFrame([fixdf(f, args.threshold, args.qvalue) for f in args.qtl], columns=['name', 'nsig'])
    d.to_csv(args.output, sep='\t', index=False)
