
# coding: utf-8

# In[15]:


import json
import pandas
import glob
import os
import argparse
import re
import sys

# In[41]:
def getOpts():
    parser = argparse.ArgumentParser(description='Make config json')
    parser.add_argument('--batch', required=True, help="""Batch name, eg. 499-NM """)
    parser.add_argument('--seqdir',  required=True, help="""path to seqcore results dir with folders  """)
    parser.add_argument('--genome',  default="hg38", help="""genome assembly - hg19 or hg38. [Default = hg38]""")
    parser.add_argument('--output', required=True, help ="""Output config file.""")
    args = parser.parse_args()
    return args

def data_path(x):
    return os.path.join(args.datadir, x)

if __name__ == '__main__':

    args = getOpts()

    batch = args.batch
    ldict = {
        'libraries': {}
    }

    readtypes = {"1": "1",
                 "2": "2"}
    readgroups = ['a','b','c','d']
    
    for lib in list(range(1, 9)):
        library = f"{batch}-{lib}"
        
        ldict['libraries'][library] = {}
        ldict['libraries'][library]['readgroups'] = {}
        ldict['libraries'][library]['genome'] = [args.genome]
            
        for rg in readgroups:
            lib_readgroup = f"{library}-3GEX_{rg}"
            rdict = {}

            for r in readtypes.keys():
                fastq = glob.glob(os.path.join(args.seqdir, f"Sample_{lib_readgroup}", f"*_R{readtypes[r]}_001.fastq.gz"))
                assert len(fastq) == 1
                print(fastq)
                rdict[r] = fastq[0]
                
            ldict['libraries'][library]['readgroups'][lib_readgroup] = rdict

    with open(args.output, 'w', encoding='utf-8') as f:
        json.dump(ldict, f, ensure_ascii=False, indent=4)
    

