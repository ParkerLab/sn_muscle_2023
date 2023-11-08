#!/usr/bin/env python3
# coding: utf-8

# In[28]:


import pandas
import numpy
import sys
import os
import glob
import time
import argparse
from string import Template
import yaml
import configargparse

def getOpts():
    parser =  configargparse.ArgParser(config_file_parser_class=configargparse.YAMLConfigFileParser,
                                       description='From a list of GWAS lead SNPs, identify regions to run susie on and make a bash script for each region')
    parser.add('--config', required=True, is_config_file=True, type=yaml.safe_load, help='config file path')
    parser.add('--base', help="""base data dir if paths are not absolute""")
    parser.add('--trait1-leads', required=True, help="""glob of headerless 4-column trait 1 bed files with lead SNPs. All should have the same --trait1-info etc.""")
    parser.add('--trait1-ref', required=True, help="""genotype dosage vcf for trait1""")
    parser.add('--trait1-ref-format', required=True, help="""Is the ref vcf chrom is of the format chr10 (chr) or 10 (int)""")
    parser.add_argument('--trait1-info', required=True, type=yaml.safe_load, help="""for each trait1, dictionary of column names for beta, se, maf etc.""")
    parser.add_argument('--trait1-type', required=True, type=yaml.safe_load, help="""for each trait1, dictionary of trait type.""")
    parser.add('--trait1-dir', required=True, help="""directory with trait1 summary stats indexed .bed.gz files""")
    parser.add('--susie-window', type=int, default=[250000], action="append", help="""Flank window on trait1 lead SNP for SuSiE""")
    parser.add('--dropsets-if-not-contain-lead', action='store_true', default=False, help="""fetch out the signal lead SNP id from the locus name, drop SuSiE sets that don't contain this snp if this parameter is set as True. Default False.""")
    parser.add('--susie-template', required=True, help="""bash template file for susie run script""")
    parser.add('--prep-template', required=True, help="""bash template file for susie prep script """)
    args = parser.parse_args()
    return args

def susie_region(chrom, pos, cflank):
    start= pos - 1 - cflank
    start = 0 if start < 0 else start
    end = pos + cflank
    fetch = f"{chrom}:{start}-{end}"
    return fetch

def get_group_val(g, col):
    val = g.iloc[0][col]
    return val

def make_susie_sh(grp, prep_template, susie_template, chrom_format):
    fetch = grp.name[0]
    susie_locus = grp.name[1]
    
    prep_sh_filename = f"{susie_locus}.susieprep.sh"
    # variables to populate the susieprep.sh template
    sub_dict = {
        "fetch": fetch,
        "ukbb_fetch": fetch.replace("chr","") if chrom_format == "int" else fetch,
        "t1_vcf": " ".join(get_group_val(grp, "t1_vcf")),
        "t1_vcf1": get_group_val(grp, "t1_vcf")[0],
        "t1_summary": get_group_val(grp, 't1_summary'),
        "susie_locus": susie_locus
    }

    with open(prep_template, 'r') as f:
        txt = Template(f.read())
        cmds = txt.substitute(sub_dict)
        
    with open(prep_sh_filename, 'w+') as f:
        f.write(cmds)

    susie_sh_filename = f"{susie_locus}.susie.sh"
    # variables to populate the susie.sh template
    all_params = get_group_val(grp, 't1_params')        
    sub_dict.update({'t1_params': all_params})
        
    with open(susie_template, 'r') as f:
        txt = Template(f.read())
        cmds = txt.substitute(sub_dict)
        
    with open(susie_sh_filename, 'w+') as f:
        f.write(cmds)
        
def checkpath(x):
    x = x.format(base = base)
    if os.path.exists(x):
        return(x)
    else:
        print(f"path {x} does not exist")
        sys.exit(1)

def get_ukbb_vcf(ukbbdir, chrom):
    """ can be multiple vcfs per chrom now"""
    vcf = sorted(glob.glob(f"{ukbbdir}/{chrom}.*.vcf.gz"))
    assert len(vcf) > 0
    return vcf

if __name__ == '__main__':
    
    args = getOpts()

    base = args.base
    trait1_leads = args.trait1_leads.format(base = base)
    trait1_ref = checkpath(args.trait1_ref)
    trait1_info = args.trait1_info
    trait1_type = args.trait1_type

    trait1_dir = checkpath(args.trait1_dir)

    susie_window = args.susie_window

    prep_template = checkpath(args.prep_template)
    susie_template = checkpath(args.susie_template)

    trait1_lead_files = glob.glob(trait1_leads)
 
    assert len(trait1_lead_files) > 0 
    
    trait1_replace = os.path.basename(trait1_leads).replace("*", "")

    for trait1 in trait1_lead_files:
        trait1_name = os.path.basename(trait1).replace(trait1_replace, "")
        ncols = len(pandas.read_csv(trait1, sep='\t', header=None).columns)
        if ncols == 4:
            colnames = ['chrom', 'start', 'snp_end', 'locus']
        elif ncols == 5:
            colnames = ['chrom', 'start', 'snp_end', 'locus', 'window']
        else:
            print(f"signal index bed file {trait1} does not have 4 or 5 columns")
            break
            
        d = pandas.read_csv(trait1, sep='\t', header=None, names=colnames,
                            dtype={'snp_end': int, 'start': int, 'window': int})
        specialChars = ", !#$%^&*();[]"
        def replaceall(string, specialChars, replaceWith):
            for c in specialChars:
                string = string.replace(c, replaceWith)
            return string
        d['trait1_locus'] = d['locus'].map(lambda x: replaceall(x, specialChars, ""))
        d['trait1_name'] = trait1_name
        d['marker'] = d['trait1_locus'].map(lambda x: x.split('__')[1])
        d.drop(['locus'], axis=1, inplace=True)
        
        print(trait1_name)
        t1_summary = glob.glob(f"{trait1_dir}/{trait1_name}*.bed.gz")
        assert len(t1_summary) == 1
        t1_summary = t1_summary[0]

        d['chrom'] = numpy.where(d['chrom']=="chr23", "chrX", d['chrom'])
        d['t1_vcf'] = d['chrom'].map(lambda x: get_ukbb_vcf(trait1_ref, x))
        d['t1_summary'] = t1_summary
        d['t1_params'] = d['trait1_name'].map(lambda x: " ".join([f'--{key} {trait1_info[key]}' for key in trait1_info.keys()]))
        d['t1_params'] = d['t1_params'] + d['trait1_name'].map(lambda x: f' --type {trait1_type[x]}')
        d['t1_params'] = d['t1_params'] + " --marker " + d['marker']

        if args.dropsets_if_not_contain_lead:
            d['t1_params'] = d['t1_params'] + " --dropsets_if_not_contain " + d['marker']
            
        if "window" in d.columns:
            d['fetch'] = d.apply(lambda x: susie_region(x['chrom'], x['snp_end'], x['window']), axis=1)
            d['susie_locus'] = d.apply(lambda x: f"{x['trait1_name']}__{x['trait1_locus']}__{x['fetch'].replace(':', '-')}__{int(x['window']/1000)}kb", axis=1)
        
            print(d.iloc[1]['t1_params'])
            d.groupby(['fetch', 'susie_locus']).apply(make_susie_sh, prep_template, susie_template, args.trait1_ref_format)
            
        else:
            for window in susie_window:
                d['fetch'] = d.apply(lambda x: susie_region(x['chrom'], x['snp_end'], window), axis=1)
                d['susie_locus'] = d.apply(lambda x: f"{x['trait1_name']}__{x['trait1_locus']}__{x['fetch'].replace(':', '-')}__{int(window/1000)}kb", axis=1)
        
                print(d.head())
                d.groupby(['fetch', 'susie_locus']).apply(make_susie_sh, prep_template, susie_template, args.trait1_ref_format)

