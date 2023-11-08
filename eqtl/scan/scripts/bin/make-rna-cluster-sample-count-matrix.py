#!/usr/bin/env python
# coding: utf-8

import pandas
import numpy
import sys
import os
import glob
import pybedtools
import argparse

def getOpts():
    parser = argparse.ArgumentParser(description='make-rna-cluster-sample-count-matrix')
    parser.add_argument('--counts', required=True, nargs='+', help ="""input gene by nuclei counts matrix. Column name syntax eg: <modality>.499.NM.1.<barcode> """)
    parser.add_argument('--cluster-info', required=True, type=str, help ="""Tab separated map of nuclei to cluster and sample. Must have columns modality, batch, sample, library, SNG.1ST, cluster_name, barcode""")
    parser.add_argument('--cluster-col', required=True, help ="""Name of the cluster column. Will output files per cluster""")
    parser.add_argument('--gene-info', help ="""If making qtltools input format file, provide gene info tsv with chrom, unique gene_name, start, end, strand""")    
    parser.add_argument('--cols', required=True, nargs='+', help ="""Space separated column names from --cluster-info file to be used""")
    parser.add_argument('--batch-col', required=True, help ="""Name of the batch column. This will be used to confirm that nuclei for one sample come from only one batch""")    
    parser.add_argument('--sample-col', required=True, help ="""Name of the sample column.""")    
    parser.add_argument('--min-samples', type=int, default=1, help ="""No. of samples in which a gene should have a nonzero value to keep that gene""")
    parser.add_argument('--min-fusion-samples', type=int, default=1, help ="""No. of samples in which a gene should have a nonzero value to keep that gene for the eQTL subset file""")
    parser.add_argument('--separate-fusion-for-eqtl', action="store_true", help ="""Separate FUSION samples and organize for sample names to match vcf. Format output per qtltools requirements.""")
    parser.add_argument('--output',  help ="""Output prefix.""")
    args = parser.parse_args()
    return args

def subset_genes(d, min_samples):
    assert "clist" not in d.columns
    d['clist'] = d[[c for c in d.columns if c!="gene_name"]].astype(bool).sum(axis=1)
    d = d[d['clist'] >= min_samples]
    assert not d.empty
    d.drop('clist', axis=1, inplace=True)
    return d


if __name__ == '__main__':
    
    args = getOpts()

    cols = args.cols
    cluster_col = args.cluster_col
    sample_col = args.sample_col
    batch_col = args.batch_col
    
    metadata = pandas.read_csv(args.cluster_info, sep='\t', usecols=cols)
    metadata = metadata[metadata['modality']=="rna"]

    clusters = set(metadata[cluster_col].tolist())

    if args.separate_fusion_for_eqtl:
        chroms = [f"chr{i}" for i in list(range(1, 23)) ]
        geneinfo = pandas.read_csv(args.gene_info, sep='\t', usecols=['chrom', 'start', 'end', 'gene_name', 'strand']).drop_duplicates()
        geneinfo = geneinfo[geneinfo['chrom'].isin(chroms)]
        geneinfo.rename(columns = {'gene_name': "gid",
                                   'chrom': "#Chr"}, inplace=True)
        geneinfo['pid'] = geneinfo['gid']
        geneinfo['length'] = geneinfo.apply(lambda x: abs(x['end'] - x['start']), axis=1)
        geneinfo = geneinfo[['#Chr', 'start', 'end', 'pid', 'gid', 'strand', 'length']]
        
    maindict = {}
    for clust in clusters:
        maindict[clust] = []

    # get sum of counts for sample in each batch file
    for f in args.counts:
        mtx = pandas.read_csv(f, sep='\t')
        mtx.rename(columns = {'Unnamed: 0':'gene_name'}, inplace=True)
        meta = metadata[metadata['index'].isin(mtx.columns)]

        for clust, cgrp in meta.groupby(cluster_col):
            o = mtx[['gene_name']]
            for samp, grp in cgrp.groupby(sample_col):
                assert len(grp[batch_col].drop_duplicates().tolist()) == 1 # one sample should come from only 1 batch
                nuclei = grp['index'].tolist()
                o[samp] = mtx[nuclei].sum(axis=1)

            o.set_index('gene_name', inplace=True)
            maindict[clust].append(o)

    # concat samples across batches for each cluster
    # Subset genes to min nonzero exp in at least x samples 
    for clust in maindict.keys():
        df = pandas.concat(maindict[clust], axis=1)
        df.fillna(0, inplace=True)
        df = df.astype(int)
        df = subset_genes(df, args.min_samples)
        
        samples = df.columns
        df.reset_index(inplace=True)
        print(df.head())
        df.rename(columns={'index':"gene_name"}, inplace=True)
        print(f"cluster {clust} had {len(samples)} samples and {len(df.index)} genes")
        df.to_csv(f"{args.output}.allsamples.{clust}.tsv", sep='\t', na_rep="NA", index=False)


        if args.separate_fusion_for_eqtl:

            fusion = metadata[(metadata['cohort']=="FUSION") & (metadata[cluster_col]==clust)]['sample'].drop_duplicates().tolist()
            fusion_samples = [s for s in fusion if s in df.columns]
            df = df[['gene_name'] + fusion_samples]
            if len(fusion_samples) == len(set([s.split('--')[0] for s in fusion_samples] )):
                rename = {s: s.split('--')[0] for s in fusion_samples}
                df.rename(columns = rename, inplace=True)

            df = subset_genes(df, args.min_fusion_samples)
            df.rename(columns = {'gene_name': "gid"}, inplace=True)
            df = pandas.merge(geneinfo, df, how="inner", on=['gid'])
                    
            df.sort_values(['#Chr','start','end'], inplace=True)

            print(f"cluster {clust} had {len(samples)} FUSION samples and {len(df.index)} genes")
            df.to_csv(f"{args.output}.fusion.{clust}.tsv", sep='\t', na_rep="NA", index=False)


