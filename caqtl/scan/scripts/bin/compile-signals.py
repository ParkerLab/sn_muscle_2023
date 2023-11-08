#!/usr/bin/env python
# coding: utf-8

import pandas
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
dpi = 150
matplotlib.rcParams['figure.dpi']= dpi
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams.update({'font.size': 12})
import argparse

def getOpts():
    parser = argparse.ArgumentParser(description='compile primary, secondary etc eQTLs across scans')
    parser.add_argument('--conditional', required=True, help ="""Formatted output from qtltools conditional scan""")
    parser.add_argument('--threshold', type=float, default=0.05, help ="""qvalue significance threshold to count caQTL. [Default = 0.05]""")
    parser.add_argument('--qvalue', type=str, default="qvalue_storey", help ="""column name with the q values. [Default = qvalue_storey]""")
    parser.add_argument('--permute', required=True, help ="""Formatted output from the permutation scan""")
    parser.add_argument('--prefix',  help ="""Output prefix.""")
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    
    args = getOpts()

    """Fix the directory where all figures will be saved"""

    # Notes
    # For some features, association rank 0 (primary association) is missing from the conditional results. 
    ## Turns out that the top variant from the permutation scan is identified as best in the forward scan with rank 0, but in the backward scan,
    ## it fails to pass the feature level threshold.. Just throw out such cases ?
    # for some features, the best SNPs don't match the best from permutation - throw out
    # Some duplicate SNP for some features - remove all these


    d = pandas.read_csv(args.conditional, delim_whitespace=True, usecols = ['gene_name', 'n_variants_tested', 'distance_var_pheno', 'snp', 'snp_chrom',
       'snp_end', 'association_rank', 'forward_p', 'is_forward_top', 'is_forward_significant', 'backward_p', 'backward_slope', 'is_backward_top',
       'is_backward_significant'])
    d1 = pandas.read_csv(args.permute, sep='\t')
    
    primary = d[d['association_rank']==0]
    permute = d1[d1[args.qvalue] < args.threshold]

    fconditional = d['gene_name'].drop_duplicates().tolist()
    fprimary = primary['gene_name'].drop_duplicates().tolist()
    fpermute = permute['gene_name'].drop_duplicates().tolist()

    print(f"N features with primary association from conditional scan = {len(fprimary)}")
    print(f"N significant features from permute scan = {len(fpermute)}")

    missing = [x for x in fpermute if x not in fprimary]
    print(f"{len(missing)} missing caQTL features in the conditional data")
    
    # remove features for which there are no primary snps
    remove = [x for x in fconditional if x not in fprimary]
    print(f"Removing {len(remove)} features in the conditional data that don't have primary SNPs")
    d = d[d['gene_name'].isin(fprimary)]

    # Remove all features with duplicate signals
    print(f"Removing {len(d[d.duplicated(['gene_name', 'snp'])]['gene_name'])} features that have duplicate SNPs in the conditional data")
    d = d[~ d.duplicated(['gene_name', 'snp'], keep=False)]

    # How many best primary SNPs are the same from the permutation scan?
    d2 = pandas.merge(primary, permute, how="inner", on=["gene_name"], suffixes=("_cond", "_perm"))
    primary_concordant = d2[d2['snp_cond'] == d2['snp_perm']]['gene_name'].to_list()
    n_primary_concordant = len(primary_concordant)
    print(f"N features with concordant primary signals = {n_primary_concordant}, out of {len(fprimary)} ({n_primary_concordant/len(fprimary)*100}) ")

    # Cases where nonprimary signals were the best from the permutation scan..
    d3 = pandas.merge(d, d1, how="inner", on=["gene_name"], suffixes=("_cond", "_perm"))
    d4 = d3[(d3['association_rank']>0) & (d3['snp_cond'] == d3['snp_perm'])]
    print(f"{len(d4.index)} features have best permute scan SNPs as the nonprimary signals")

    # Only keep features where the primary SNPs are concordant with the permutation scan
    d = d[d['gene_name'].isin(primary_concordant)]
    features = len(d[['gene_name']].drop_duplicates().index)
    permute_features = len(permute.index)
    secondary = len(d.index) - features
    title = f"{args.prefix} {permute_features} primary caQTL and {secondary} secondary caQTL"
    print(title)
    
    # plot
    # plt.figure(figsize=(4,3))
    # d['abs_distance'] = d['distance_var_pheno'].map(abs)
    # sns.ecdfplot(data=d[d['association_rank']<3], x="abs_distance", hue="association_rank")
    # plt.title(title)
    # plt.savefig(f"fig.{args.prefix}.conditional_distance.png", bbox_inches="tight", dpi=200)
    
    # save filtered
    d.to_csv(f"{args.prefix}.conditional.tsv", sep='\t', index=False)




