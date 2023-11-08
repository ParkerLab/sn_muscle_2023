#!/usr/bin/env python
# coding: utf-8

# In[1]:

import pandas
import numpy
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import os
dpi = 150
matplotlib.rcParams['figure.dpi']= dpi
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams.update({'font.size': 12})
import argparse
import sys

# Filter features now:
def fix_index(x):
    end = x.split(':')[-1]
    new = x.replace(f":{end}", f"-{end}")
    return new

def get_feature_counts(m, minct, annot_count, annot_sample):
    df = m.copy()
    df[df<minct] = 0
    nsamp = df.astype(bool).sum(axis=1)
    s = pandas.DataFrame(list(range(1, 100)), columns=['n'])
    s['features'] = s['n'].map(lambda x: len(nsamp[nsamp>=x].index))

    sns.lineplot(data=s, x='n', y='features', label=minct)    
    if minct == annot_count:
        xannot = annot_sample
        yannot = s[s['n'] == annot_sample].iloc[0]['features']
        plt.plot([xannot, xannot], [0, yannot], color="black", linestyle="dashed")
        plt.plot([0, xannot], [yannot, yannot], color="black", linestyle="dashed")
        plt.text(xannot + 1, yannot + 1, yannot)

    s['minct'] = minct
    return s[['n','minct','features']]
        
def getOpts():
    parser = argparse.ArgumentParser(description='Make N peak at various min count and min sample thresholds')
    parser.add_argument('--peak-counts', required=True, help ="""input peak count matrix""")
    parser.add_argument('--samples', required=True, help ="""One column file with sample as header containing samples to select""")
    parser.add_argument('--cluster', required=True, help ="""Cluster name""")
    parser.add_argument('--min-counts', required=True, nargs='+', type=int, help ="""Space-separated list of min counts to check""")
    parser.add_argument('--annot-count', required=True, type=int, help ="""Annotate curve for N features passing this min count threshold""")
    parser.add_argument('--annot-sample', required=True, type=int, help ="""Annotate curve for N features passing this min sample threshold""")
    parser.add_argument('--prefix', required=True, help ="""output file prefix""")

    args = parser.parse_args()
    return args


if __name__ == '__main__':

    args = getOpts()


    min_counts = args.min_counts 
    
    m = pandas.read_csv(args.peak_counts, sep='\t', index_col=0)

    samples = pandas.read_csv(args.samples, sep='\t', squeeze=True).tolist()
    
    selected_samples = [s for s in samples if s in m.columns]
    missing_samples = [s for s in samples if s not in m.columns]
    
    assert len(missing_samples) == 0
        # print(f"WARNING: samples {','.join(missing_samples)} were not found in the counts matrix.")
    
    m = m[selected_samples]
    
    if m.empty:
        print("Empty dataframe, check script and params")
        sys.exit(1)
        
    plt.figure(figsize=(4,3))

    dflist = []
    for c in min_counts:
        dflist.append(get_feature_counts(m, c, args.annot_count, args.annot_sample))
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    plt.title(args.cluster)
    plt.savefig(f"{args.prefix}.png", bbox_inches="tight")


    df = pandas.concat(dflist)
    df['cluster'] = args.cluster
    df.to_csv(f"{args.prefix}.tsv", sep='\t', index=False, na_rep="NA")





