#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas
import numpy
import sys
import os
import glob


def get_info(x, col):
    out = [y.replace(f"{col} ", "").replace('"', '').replace(" ", "") for y in x.split(';') if col in y][0]
    return out

biotypes = ['protein_coding']

annot = pandas.read_csv(sys.argv[1], sep='\t', header=None, 
                        names=['chrom', 'i', 'type', 'start', 'end', 'a', 'strand', 'b', 'info'],
                       usecols=['chrom', 'start', 'end', 'type', 'strand', 'info'])
annot = annot[annot['type']=="gene"]

for col in ['gene_id', 'gene_type', 'gene_name']:
    annot[col] = annot['info'].map(lambda x: get_info(x, col))

# select biotype
annot = annot[annot['gene_type'].isin(biotypes)]

# Get TSS based on strand and gene lengths
annot['tss_start'] = annot.apply(lambda x: x['start']-1 if x['strand']=="+" else x['end']-1, axis=1)
annot['tss_end'] = annot.apply(lambda x: x['start'] if x['strand']=="+" else x['end'], axis=1)
annot['length'] = annot['end'] - annot['start']

# Duplicate gene names..
# Update duplicate gene names:
dup_names = annot[annot['gene_name'].duplicated()]['gene_name'].tolist()
annot['gene_name'] = annot.apply(lambda x: f"{x['gene_name']}--{x['gene_id']}" if x['gene_name'] in dup_names else x['gene_name'], axis=1)

annot['pid'] = annot['gene_name']

annot = annot[['chrom', 'tss_start', 'tss_end', 'pid', 'gene_name', 'strand', 'length']]
annot.rename(columns = {'chrom': '#Chr',
                   'tss_start': "start",
                   'tss_end': "end",
                   "gene_name": "gid"}, inplace=True)

annot.to_csv(sys.argv[2], sep='\t', index=False)
