# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 10:32:47 2020

@author: mitch
Platform Unix
"""
import pandas as pd
import subprocess
import os

def createBaseInput(infile, outfile):
    homer_cols = ['name','chr','start','end','strand']
    df_m6a = pd.read_table(infile, header=None,
                           names=homer_cols)
    df_m6a['score'] = 1
    infile_back = infile.replace('.homer', '.background.homer')
    df_back = pd.read_table(infile_back, header=None,
                            names=homer_cols)
    df_back['score'] = 0
    
    df = pd.concat([df_m6a, df_back])
    
    df['start'] += 99
    df['end'] -= 100
    
    bed6_cols = ['chr','start','end','name','score','strand']
    df = df[bed6_cols]
    
    df.to_csv(outfile, sep='\t', header=False, index=False)
    
    
lst_df = []
base_dir = '../'

infiles = [os.path.join(base_dir,'processed_data','{}_m6a.201.homer'.format(x)) \
           for x in ['a549','cd8t','hek293']]
outfiles = [x.replace('_m6a.201.homer','_training.bed') for x in infiles]

for infile, outfile in zip(infiles,outfiles):
    createBaseInput(infile, outfile)
    py_script = os.path.join(base_dir,'python','20200827_generate_model_input.py')
    genome_fa = '/mnt/e/Genomes/hg19/hg19.fa'
    phylop = '/mnt/e/Genomes/hg19/hg19.100way.phyloP100way.bw'
    outdir = os.path.join(base_dir,'processed_data')
    base_outfile = os.path.basename(outfile)
    sample = base_outfile.split('_')[0]
    model_input = '{}_training_input.csv.gz'.format(sample)
    cmd = 'python {} --bed {} --labels score --fasta {} --phylop {} --outdir {} --outfile {}'.format(
        py_script, outfile, genome_fa, phylop, outdir, model_input)
    subprocess.call(cmd, shell=True)
    df = pd.read_csv(os.path.join(outdir, model_input))
    df['sample'] = sample
    lst_df.append(df.copy(deep=True))
    
df_final = pd.concat(lst_df)
df_final.to_csv(os.path.join(base_dir,'processed_data','20200918_training.csv.gz'),
                index=False)