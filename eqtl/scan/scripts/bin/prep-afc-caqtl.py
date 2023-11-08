#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import numpy as np
from qtl import genotype, norm
import glob
import re
import os
from plutils import general
import sys
import logging

#CELL_TYPE = 'Endothelial'
#PREFIX = 'test.'
CELL_TYPE, PREFIX = sys.argv[1:]

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

EQTL_RESULTS = glob.glob(f'/lab/work/arushiv/muscle-sn/analyses_hg38/caqtl/consensus_0327_freeze/results/caqtl-signals/fusion.{CELL_TYPE}--covlist.*.caqtl.tsv')
assert(len(EQTL_RESULTS) == 1)
EQTL_RESULTS = EQTL_RESULTS[0]
COVLIST_NUMBER = re.match('fusion.*--covlist.(\d+).caqtl.tsv', os.path.basename(EQTL_RESULTS)).group(1)
COVLIST = f'/lab/work/arushiv/muscle-sn/analyses_hg38/caqtl/consensus_0327_freeze/results/covariate-lists/covlist.{COVLIST_NUMBER}.txt'
COVARIATE_MATRIX = f'/lab/work/arushiv/muscle-sn/analyses_hg38/caqtl/consensus_0327_freeze/results/sample-info/fusion.{CELL_TYPE}-covariates.tsv'
COUNT_MATRIX = f'/lab/work/arushiv/muscle-sn/analyses_hg38/files_for_peter/atac/quantified_consensus_features/atac.{CELL_TYPE}.tsv'
VCF = '/lab/work/arushiv/muscle-sn/analyses_hg38/eqtl/eqtl_11/results/vcf/filtered-vcf.maf0.05-hwe1e6.recode.vcf.gz'
SAMPLE_LIST = f'/lab/work/arushiv/muscle-sn/analyses_hg38/caqtl/consensus_0327_selected/data/sample-lists/fusion.{CELL_TYPE}.sample_list.txt'
# TODO: taking sample list from covariate matrix...which is maybe not what we want

GENOTYPE_TO_DOSAGE = {
    '1|1': 2,
    '0|1': 1,
    '1|0': 1,
    '0|0': 0
}


USE_SAMPLES = pd.read_csv(SAMPLE_LIST).iloc[:,0].to_list()
USE_SAMPLES = set([x.split('--')[0] for x in USE_SAMPLES])

# read in caQTL results
caqtl = pd.read_csv(EQTL_RESULTS, sep='\t')

regions = caqtl.snp_chrom + ':' + caqtl.snp_end.astype(str) + '-' + caqtl.snp_end.astype(str)
regions = list(regions.unique())

genotypes = general.tabix(VCF, regions).drop_duplicates()

VCF_SAMPLES = genotypes.columns.to_list()[9:]
genotype_dosages = genotypes[['ID'] + VCF_SAMPLES].drop_duplicates().set_index('ID')
#genotype_dosages = genotypes[['ID'] + USE_SAMPLES].set_index('ID')
genotype_dosages = genotype_dosages.transform(lambda x: x.str.split(':', expand=True)[0])
assert(genotype_dosages.isin(GENOTYPE_TO_DOSAGE.keys()).all().all())
genotype_dosages = genotype_dosages.transform(lambda x: x.map(GENOTYPE_TO_DOSAGE))
genotype_dosages = genotype_dosages.loc[list(caqtl.snp.unique())]

# make variant df
variant_df = genotypes[['CHROM', 'POS', 'ID']].drop_duplicates().set_index('ID').rename(columns={'CHROM': 'chrom', 'POS': 'pos'}).loc[genotype_dosages.index.to_list()]

covariates = pd.read_csv(COVARIATE_MATRIX, sep='\t')
covariate_list = pd.read_csv(COVLIST, header=None)[0].to_list()
covariates = covariates.set_index('id').loc[covariate_list,:]
COVARIATE_SAMPLES = set(covariates.columns)

counts = pd.read_csv(COUNT_MATRIX, sep='\t')
counts = counts.set_index('feature')
counts.columns = [i.split('--')[0] for i in counts.columns]
#counts = counts[USE_SAMPLES]
COUNTS_SAMPLES = set(counts.columns)

logging.info('{} samples in samples list'.format(len(USE_SAMPLES)))
USE_SAMPLES = list(USE_SAMPLES.intersection(COUNTS_SAMPLES).intersection(COVARIATE_SAMPLES).intersection(VCF_SAMPLES))
logging.info('{} samples in VCF'.format(len(VCF_SAMPLES)))
logging.info('{} samples in covariates file'.format(len(COVARIATE_SAMPLES)))
logging.info('{} samples in counts file'.format(len(COUNTS_SAMPLES)))
logging.info('{} samples in intersection of the above'.format(len(USE_SAMPLES)))
covariates = covariates[USE_SAMPLES]
genotype_dosages = genotype_dosages[USE_SAMPLES]
counts = counts[USE_SAMPLES]

normalized_counts = norm.deseq2_normalized_counts(counts)

# output everything
caqtl[['gene_name', 'snp', 'slope']].rename(columns={'gene_name': 'gene_id', 'snp': 'variant_id'}).to_csv(f'{PREFIX}assoc-df.tsv', sep='\t', index=False)
normalized_counts.to_csv(f'{PREFIX}normalized-counts.tsv', sep='\t')
genotype_dosages.to_csv(f'{PREFIX}genotype-dosages.tsv', sep='\t')
variant_df.to_csv(f'{PREFIX}variant-df.tsv', sep='\t')
covariates.to_csv(f'{PREFIX}covariates.tsv', sep='\t')