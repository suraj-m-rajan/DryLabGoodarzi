#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 10:57:37 2020

@author: mitch
Platform: Unix

Given bed file, genome fasta, and phyloP bigwig, creates input for 
CNN model
"""

import pandas as pd
import numpy as np
from scipy.stats import entropy
import pyBigWig
import subprocess
import os
import argparse
import logging
import sys
from Bio import SeqIO
from itertools import product

def parse_arguments():
    """Read arguments from the command line"""    
    parser = argparse.ArgumentParser(description='Process genomics data to create '
                                     + 'input for CNN m6A prediction model')
    parser.add_argument('--bed', required=True,
                        help='Path to bed file of potential m6A sites\n'
                        + 'Width of sequences is 1 (0-indexed). BED6 Format\n' 
                        + 'chrom\tstart\tend\tname\tscore\tstrand')
    parser.add_argument('--labels', default='unknown', 
                        choices=['unknown','pos','neg','score'],
                        help='m6A status for sequences\n'
                        + 'Use unknown for prediction.'
                        + 'Use pos or neg if all sequences are '
                        + 'methylated (pos) or not(neg). '
                        + 'Use score if score column contains label '
                        + '1 for methylated and 0 for unmethylated')
    parser.add_argument('--fasta', required=True,
                        help='Path to genome fasta file')
    parser.add_argument('--phylop', required=True,
                        help='Path to phyloP bigwig file')
    parser.add_argument('--phylop_min', default=-20, type=float,
                        help='Genome phyloP min')
    parser.add_argument('--phylop_max', default=9.873, type=float,
                        help='Genome phyloP max')
    parser.add_argument('--outdir', default='./')
    parser.add_argument('--outfile', default='cnn_input.csv.gz')
    args = parser.parse_args()
    
    logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s', 
                        level=logging.INFO, datefmt='%m.%d/%Y %I:%M:%S %p',
                        stream=sys.stdout)
    logging.info('\n' + ' '.join(sys.argv) + '\n')
    return args


def check_args(args):
    """Check that bed file is in expected format"""
    bed = args.bed
    df_bed = pd.read_table(bed, header=None)
    # check formatting
    try:
        assert(df_bed.shape[1] == 6)
    except AssertionError:
        logging.error('Bed file does not contain 6 columns. Aborting')
        sys.exit(1)
        
    df_bed.columns = ['chrom', 'start','end','name','score','strand']
    
    if not 'chr' in df_bed['chrom'].iloc[0]:
        df_bed['chrom'] = df_bed['chrom'].map('chr{}'.format)
    
    if args.labels == 'unknown':
        df_bed['score'] = 'NA'
    elif args.labels == 'pos':
        df_bed['score'] = 1
    elif args.labels == 'neg':
        df_bed['score'] = 0
    else:
        try:
            assert(set(df_bed['score'].values).difference([0,1]) == set())
        except AssertionError:
            logging.error('Score column of bed files contain values other '
                          + 'than 0 (unmethylated) and 1 (methylated). Aborting')
            sys.exit(1)
    
    return df_bed

def expandBed(df, size, outfile):
    df_bed = df.copy(deep=True)
    df_bed['start'] -= size
    df_bed['end'] += size
    df_bed['name'] = df_bed.apply(lambda row:'{}:{}-{}({}):{}'.format(
    row['chrom'], row['start'], row['end'], row['strand'], row['score']),
    axis=1)
    cols = ['chrom','start','end','name','score','strand']
    df_bed.to_csv(outfile, sep='\t', header=False, columns=cols, index=False)


def writeBed(df_bed, outdir):
    """Write bed file for bedtools getFasta"""
    #expand bed    
    cnn_m6a_bed = os.path.join(outdir, 'cnn201.bed')
    expandBed(df_bed, 100, cnn_m6a_bed)

    logging.info('Successfully expanded bed files')

    
def getFasta(outdir, fasta):
    """Call bedtools getfasta"""
    cnn_bed = os.path.join(outdir, 'cnn201.bed')
    
    try:
        assert(os.path.exists(cnn_bed))
    except AssertionError:
        logging.error('Expanded bed files are not found. Aborting.')
        sys.exit(1)
        
    cnn_fa = os.path.join(outdir, 'cnn201.fa')

    cmd_cnn = 'bedtools getfasta -s -name -fo {} -fi {} -bed {}'.format(
        cnn_fa, fasta, cnn_bed)
    try:
        logging.debug(cmd_cnn)
        assert(subprocess.call(cmd_cnn, shell=True) == 0)
    except AssertionError:
        logging.error('Could not obtain fasta files from bed files. Aborting.')
        sys.exit(1)
    
    logging.info('Successfuly extracted fasta sequences')

def retrieveBigWig(row, bw, p_max, p_min):
    scores = bw.intervals(row['chr'], int(row['start']), int(row['end']))
    if scores:
        scores = [(x[2] - p_min)/(p_max - p_min) for x in scores]
        if row['strand'] == '-':
            scores = scores[::-1]
        return pd.Series(scores)
    else: return pd.Series([np.nan] * int(row['end'] - row['start']))

def createCNNBaseInput(outdir, outfile, phylo, p_max, p_min):
    """Create input (dataframe) with 201 bp, kmers and phyloP"""
    cnn_fa = os.path.join(outdir, 'cnn201.fa')
    seqrecords = list(SeqIO.parse(cnn_fa, 'fasta'))
    seqs = [str(x.seq).upper() for x in seqrecords]
    names = [x.name for x in seqrecords]
    seq_len = len(seqs[0])
    df_seq = pd.DataFrame(list(zip(seqs,names)), columns=['seq','loc'])
    bases = ['A','C','G','T']
    twomers = [''.join(p) for p in product(bases, repeat=2)]
    threemers = [''.join(p) for p in product(bases, repeat=3)]
    
    logging.debug('Calculating kmer frequencies')
    kmer_cols = bases + twomers + threemers 
    for kmer in kmer_cols:
        regex_pat = '(?={})'.format(kmer)
        df_seq[kmer] = df_seq['seq'].str.count(regex_pat) / (seq_len // len(kmer))
    logging.info('Succssfully calculated kmer frequencies')
    
    df_seq['gc_content'] = df_seq['G'] + df_seq['C']
    df_seq['amino_content'] = df_seq['A'] + df_seq['C']
    df_seq['purine_content'] = df_seq['A'] + df_seq['G']
    
    df_seq['s_entropy'] = df_seq[bases].apply(entropy, axis=1)
    property_cols = ['gc_content','amino_content','purine_content','s_entropy']
    bp_cols = ['bp_{}'.format(i) for i in range(seq_len)]
    df_seq[bp_cols] = df_seq['seq'].apply(lambda x: pd.Series(list(x)))
    
    bw = pyBigWig.open(phylo)
    try:
        assert(bw.isBigWig())
    except AssertionError:
        logging.error('PhyloP bigwig file is not correctly formatted. Aborting')
        bw.close()
        sys.exit(1)
    
    logging.debug('Retrieving conservation scores')
    chr_regex_pat = '(chr\w+):(\d+)-(\d+)\(([+|-])\):([0|1|NA])'
    loc_cols = ['chr','start','end','strand','group']
    df_seq[loc_cols] = df_seq['loc'].str.extract(chr_regex_pat)
    df_seq[['start','end']] = df_seq[['start','end']].astype(float)
    phylo_cols = ['phylo_{}'.format(i) for i in range(seq_len)]
    df_seq[phylo_cols] = df_seq[loc_cols].apply(lambda row: retrieveBigWig(row, bw, p_max, p_min),
                                                axis=1)
    logging.info('Successfully retrieved conservation scores.')
    df_seq = df_seq[bp_cols + phylo_cols + kmer_cols + property_cols + ['group']]
    df_seq = df_seq.dropna(subset=phylo_cols)
    df_seq.to_csv(os.path.join(outdir, outfile), index=False)
    
def main():
    args = parse_arguments()
    df_bed = check_args(args)
    writeBed(df_bed, args.outdir)
    getFasta(args.outdir, args.fasta)
    createCNNBaseInput(args.outdir, args.outfile, args.phylop, args.phylop_max, args.phylop_min)
    return 0

if __name__ == '__main__':
    main()