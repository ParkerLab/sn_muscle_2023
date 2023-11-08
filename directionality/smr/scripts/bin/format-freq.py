#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas
import numpy
import sys
import os


d = pandas.read_csv(sys.argv[1], sep='\t', skiprows=[0], header=None,
                                        names=['chrom', 'pos', 'N_ALLELES', 'ext', 'refinfo', 'altinfo'])
d = d[~ d['altinfo'].str.contains("<")]
d[['EA', 'EAF']] = d['altinfo'].str.split(":", expand=True, n=2)
d['NEA'] = d['refinfo'].map(lambda x: x.split(":")[0])
d['chr'] = d['chrom'].str.replace("chr", "")
d['gd'] = 0
d = d[['chrom', 'gd', 'pos', 'EA', 'NEA', 'EAF']]
len1 = len(d.index)

di = pandas.read_csv(sys.argv[2], sep='\t', comment="#",
                     usecols=[0,1,2,3,4], names=['chrom', 'pos', 'id', 'NEA', 'EA']).drop_duplicates(['chrom','pos','NEA','EA'])

d = pandas.merge(d, di, how="inner", on=['chrom', 'pos', 'EA', 'NEA'])
len2 = len(d.index)

print(f"{len1}, new = {len2}")
#assert len1 == len2

d['chrom'] = d['chrom'].str.replace("chr", "")

d = d[['chrom', 'id', 'gd', 'pos', 'EA', 'NEA', 'EAF']]

d.to_csv(sys.argv[3], sep='\t', index=False, header=False)
