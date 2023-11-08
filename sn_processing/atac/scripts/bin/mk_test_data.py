
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
    parser.add_argument('--src', required=True, help="""script directory. Should also contain environment.yml """)
    parser.add_argument('--email', default="arushiv@umich.edu", help="""email """)
    parser.add_argument('--results', required=True, help="""results dir """)
    parser.add_argument('--seqdir',  required=True, help="""path to seqcore results dir parker/  """)
    parser.add_argument('--datadir', required=True, help ="""Path to data directory""")
    parser.add_argument('--wtype', required=True, help ="""type of workflow - provide atac or rna""")
    parser.add_argument('--output', required=True, help ="""Output config file.""")
    args = parser.parse_args()
    return args

def data_path(x):
    return os.path.join(args.datadir, x)

if __name__ == '__main__':

    args = getOpts()
    
    ldict = {'results': args.results,
             'src': args.src,
             'email': args.email,
             "whitelist_10X": data_path("barcode-whitelist/3M-february-2018.txt"),
             "barcode-whitelist": data_path("barcode-whitelist/737K-cratac-v1.txt"),
             "blacklist": {
                 "hg19": [
                     data_path("mappability/wgEncodeDacMapabilityConsensusExcludable.bed.gz"),
                     data_path("mappability/wgEncodeDukeMapabilityRegionsExcludable.bed.gz")
                 ]
             },
             "bwa_index": {
                 "hg19": data_path("bwa/hg19/hg19")
             },
             "tss": {
                 "hg19": data_path("tss/hg19.tss.bed.gz")
             },
             "whitelist": {
                 "hg19": []
             },
             'vcf': "/lab/work/arushiv/muscle-sn/work/Run_3094/data_vcfs/batch1-all.vcf.gz",
             'vcf_index': "/lab/work/arushiv/muscle-sn/work/Run_3094/data_vcfs/batch1-all.vcf.gz.tbi",
             "reorder_bam": "/lab/work/arushiv/muscle-sn/data/reorder-bam/reorder.dict"
    }

    
    if args.wtype == "rna":
        readtypes = {"1": "1",
                     "2": "2"}
    elif args.wtype == "atac":
        readtypes = {"1": "1",
                          "2": "3",
                          "index": "2"}
    else:
        print("Specify either atac or rna")
        sys.exit(1)

    if args.seqdir:
        parkerdir = os.path.join(args.seqdir, "parker/")
        ldict['libraries'] = {}
        dirs = os.listdir(parkerdir)
        libraries = list(set([x.split('_')[1] for x in dirs]))

        for lib in libraries:
            ldict['libraries'][lib] = {}
            ldict['libraries'][lib]['readgroups'] = {}
            ldict['libraries'][lib]['genome'] = ["hg19"]

            lib_readgroups = [d for d in dirs if re.match(f"Sample_{lib}_", d)]
            if args.wtype == "atac":
                cellstring = os.listdir(os.path.join(parkerdir, lib_readgroups[0]))[0].split("_S")[0]
                cellranger = glob.glob(os.path.join(args.seqdir, f"Sample_{cellstring}/outs/possorted_*bam.bam"))
                assert len(cellranger) == 1
                cellranger = cellranger[0]
                ldict['libraries'][lib]['cellranger'] = cellranger
                

            for rg in lib_readgroups:
                rdict = {}
                for r in readtypes.keys():
                    fastq = glob.glob(os.path.join(parkerdir, rg, f"*_R{readtypes[r]}_001.fastq.gz"))[0]
                    print(fastq)
                    rdict[r] = fastq
                ldict['libraries'][lib]['readgroups'][rg] = rdict


    with open(args.output, 'w', encoding='utf-8') as f:
        json.dump(ldict, f, ensure_ascii=False, indent=4)
    

