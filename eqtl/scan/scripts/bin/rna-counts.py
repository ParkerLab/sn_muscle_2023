#!/usr/bin/env python
# coding: utf-8

import pandas
import numpy
import math
import os
import glob
import argparse

def getOpts():
    parser = argparse.ArgumentParser(description='make-rna-cluster-sample-count-matrix')
    parser.add_argument('--counts-dir', required=True, help ="""directory containing all decontx outputs across libraries. Will look for 10X sparse matrix folders named <library>.decont/ """)
    parser.add_argument('--clusters', nargs='+', help ="""clusters to subset""")
    parser.add_argument('--aggregate', action="store_true", help ="""If this flag is provided, will consider all nuclei in aggregate irrespective of cluster""")    
    parser.add_argument('--cluster-info', required=True, type=str, help ="""Tab separated map of nuclei to cluster and sample. Must have columns modality, batch, sample, library, SNG.1ST, cluster_name, barcode""")
    parser.add_argument('--gene-info', help ="""If making qtltools input format file, provide gene info tsv with chrom, unique gene_name, start, end, strand""")
    parser.add_argument('--gene-types', nargs='+', help ="""gene types to fetch using the "gene_type" column in the --gene-info file.""")
    parser.add_argument('--exclude-chr', nargs='+', help ="""Exclude genes on these chromosomes""")    
    parser.add_argument('--min-fusion-samples', type=int, default=1, help ="""No. of samples in which a gene should have a nonzero value to keep that gene for the eQTL subset file""")
    parser.add_argument('--prefix',  help ="""Output prefix.""")
    args = parser.parse_args()
    return args


def get_counts(matrix_dir, genedf, metadata):
    library = os.path.basename(matrix_dir).split(".")[0]
    mat = pandas.read_csv(os.path.join(matrix_dir, "matrix.mtx"), sep=' ', skiprows=[0,1,2], header=None, names=['gid', 'bid', 'count'])

    features_path = os.path.join(matrix_dir, "features.tsv")
    fdf = pandas.read_csv(features_path, sep="\t", header=None, names=['gene_id', 'gene_name', 'gt'], usecols = ['gene_id']).reset_index().rename(columns = {'index': "gid"})
    fdf['gid'] = fdf['gid'] + 1
    fdf = pandas.merge(fdf, genedf[['gene_id', 'gene_name']], how="inner", on="gene_id")
    fi = fdf['gid']
    
    barcodes_path = os.path.join(matrix_dir, "barcodes.tsv")
    bdf = pandas.read_csv(barcodes_path, sep="\t", header=None, names=['bc']).reset_index().rename(columns = {'index': "bid"})
    bdf['bc'] = bdf['bc'] + "." + library
    bdf['bid'] = bdf['bid'] + 1
    bdf = pandas.merge(bdf, metadata[['bc','coarse_cluster_name','SNG.1ST']], how="inner", on="bc")
    bi = bdf['bid']

    mat = mat[(mat['gid'].isin(fi)) & (mat['bid'].isin(bi))]
    if not mat.empty:
        mat = pandas.merge(pandas.merge(mat, fdf, how='inner', on='gid'),
                           bdf, how='inner', on='bid')
        print(mat.head())
        print(mat[mat['count'].isna()])
        mat['count'] = mat['count'].map(lambda x: round(x, 0))
        mat = mat[mat['count']>0]
        mat = mat.drop(['gid','bid','gene_id'], axis=1).groupby(['coarse_cluster_name', 'SNG.1ST', 'gene_name']).sum().reset_index()
    else:
        mat = pandas.DataFrame()
    return mat


if __name__ == '__main__':
    
    args = getOpts()

    # format barcodes
    metadata = pandas.read_csv(args.cluster_info, sep='\t', dtype={'SNG.1ST':str})
    metadata = metadata[(metadata['modality']=="rna") & (metadata['cohort'] == "FUSION")]
    metadata['bc'] = metadata['barcode'] + '.' + metadata['library']
    metadata = metadata[['bc', 'coarse_cluster_name', 'SNG.1ST']]

    if args.aggregate:
        meta1 = metadata.copy()
        meta1['coarse_cluster_name'] = "aggregate"
        metadata = pandas.concat([metadata, meta1]).drop_duplicates()
        
    if args.clusters is not None:
        clusters = args.clusters
        if args.aggregate:
            clusters += "aggregate"
            clusters = list(set(clusters))
        metadata = metadata[metadata['coarse_cluster_name'].isin(clusters)]
    else:
        clusters = metadata['coarse_cluster_name'].drop_duplicates().tolist()

    # select genes
    g = pandas.read_csv(args.gene_info, sep='\t')
    if args.gene_types is not None:
        g = g[(g['gene_type'].isin(args.gene_types))]
    if args.exclude_chr is not None:
        g = g[(~g['chrom'].isin(args.exclude_chr))]
    # update duplicate gene names:
    g1 = g[~g['gene_name'].duplicated(keep=False)]
    g2 = g[g['gene_name'].duplicated(keep=False)].sort_values("gene_name")
    g2['gene_name'] = g2['gene_name'] + "_" + g2['gene_id']
    g = pandas.concat([g1, g2])
    print(g.head())

    # iterate over decontx folders, per batch
    l = pandas.DataFrame({'path': glob.glob(f"{args.counts_dir}/*.decont")})
    l['library'] = l['path'].map(lambda x: os.path.basename(x).split(".")[0])
    l['batch'] = l['library'].map(lambda x: x.split('-')[0])
    l.sort_values(['batch', 'library'], inplace=True)
    l.head()

    
    m = {}
    for b, grp in l.groupby(['batch']):
        print(f"starting {b}")
        m[b] = pandas.concat([get_counts(dirname, g, metadata) for dirname in grp['path'].tolist()])
        print(f"nrows = {len(m[b].index)}")
        if not m[b].empty:
            m[b] = m[b].groupby(['coarse_cluster_name', 'SNG.1ST', 'gene_name']).sum().reset_index()
            print(f"rows after adding up {len(m[b].index)}")

    mtx = pandas.concat([m[b] for b in m.keys()])

    # gene for where at least n samples
    g = g[['chrom', 's', 'e', 'gene_name', 'strand', 'length']].rename(columns = {'chrom': '#Chr',
                                                                        's': 'start',
                                                                        'e': 'end',
                                                                        'gene_name': 'gid'})
    g['pid'] = g['gid']
    g = g[['#Chr', 'start', 'end', 'gid', 'pid', 'strand', 'length']]
    
    for c, grp in mtx.groupby("coarse_cluster_name"):
        nsamp = grp[['gene_name', 'SNG.1ST']].groupby('gene_name').size().to_frame(name="n").reset_index()
        genes = nsamp[nsamp['n']>=args.min_fusion_samples]['gene_name']
        grp = grp[grp['gene_name'].isin(genes)]
        out = grp[['SNG.1ST', 'gene_name', 'count']].pivot_table(index="gene_name", columns="SNG.1ST", values="count").fillna(0).astype(int).reset_index().rename(columns = {'gene_name': 'gid'})
        out = pandas.merge(g, out, how="inner", on="gid")
        print(f"cluster {c} had {len(out.columns)-1} FUSION samples and {len(out.index)} genes")
        out.to_csv(f"{args.prefix}.{c}.tsv", sep='\t', index=False)
        




