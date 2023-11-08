#!/usr/bin/env python

import os
import sys
import pysam
import json
import argparse
import logging
import numpy as np

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

parser = argparse.ArgumentParser(description='Gather stats from a snRNA-seq bam file.', add_help = True)
parser.add_argument('bam', type = str,  help = 'BAM file, sorted by read name.')
parser.add_argument('--cell-tag', dest='cell_tag', type = str, default = 'CB', help = 'Tag denoting the cell/nucleus (default: CB)')
parser.add_argument('--gene-tag', dest='gene_tag', type = str, default = 'GX', help = 'Tag denoting the gene assignment (default: GX)')
parser.add_argument('--umi-tag', dest='umi_tag', type = str, default = 'UB', help = 'Tag denoting the gene assignment (default: UB)')
parser.add_argument('--min-reads', dest='min_reads', type = int, default = 0, help = 'Suppress output for cells with fewer than this many reads (default: 0).')
args = parser.parse_args()

CELL_TAG = args.cell_tag
GENE_TAG = args.gene_tag
UMI_TAG = args.umi_tag
MIN_READS_FOR_OUTPUT = args.min_reads

from cell import Cell

if __name__ == '__main__':
    cells = dict()
    no_cell_tag = 0
    total_reads = 0

    with pysam.AlignmentFile(args.bam, 'rb') as f:
        for read in f.fetch(until_eof=True):
            total_reads += 1
            if total_reads % 1000000 == 0:
                logging.info('Processed {} reads'.format(total_reads))
            if not read.has_tag(CELL_TAG):
                no_cell_tag += 1
                continue
            barcode = read.get_tag(CELL_TAG)
            if barcode not in cells:
                cells[barcode] = Cell(barcode, GENE_TAG, UMI_TAG)
            cells[barcode].record_alignment(read)

    logging.info('Finished reading bam file.')
    JSON = [cell.gather_metrics() for  barcode, cell in cells.items() if cell.total_reads >= MIN_READS_FOR_OUTPUT ]
    print(json.dumps(JSON, indent=4, sort_keys=True))
    logging.info('Encountered {} reads (of {}) with no cell tag'.format(no_cell_tag, total_reads))
