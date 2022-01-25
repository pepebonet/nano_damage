#!/usr/bin/envs python3 

import os
import sys
import click
from itertools import islice
import pandas as pd
from collections import Counter

from bgreference import refseq

sys.path.append('../')
import plots as pl
import utils as ut

names_all_1 = ['chrom', 'pos', 'base', 'strand', 'mer', \
    'min_coverage', 'untreated_freq', 'treated_freq']

names_all_2 = ['index', 'chrom', 'pos', 'base', 'strand', 'mer', 'min coverage',
       'untreated_freq', 'treated_freq', 'diff', 'group',
       'motif_1 DRWGGDD P-value', 'motif']

names_all_3 = ['index', 'chrom', 'pos', 'base', 'strand', 'mer', 'min coverage',
       'untreated_freq', 'treated_freq', 'diff', 'group',
       'motif_1 P-value', 'motif_2 P-value', 'motif']


def get_contexts(df, chrom_len): 
    penta_ref, triplet_ref = ut.counts_reference_genome(chrom_len)
    
    df['PENTAMER'] = df.apply(obtain_context, args=(5,2), axis=1)
    df['TRIPLET'] = df.apply(obtain_context, args=(3,1), axis=1)
    
    penta_exp = get_context_counts(df, 'PENTAMER')
    triplet_exp = get_context_counts(df, 'TRIPLET')

    triplet_context = get_context_norm(triplet_ref, triplet_exp)
    penta_context = get_context_norm(penta_ref, penta_exp)  

    return df, triplet_context, penta_context 


#Otain the windows of size n of the yeast genome
def slicing_window(seq, n):

    it = iter(seq)
    result = ''.join(islice(it, n))

    if len(result) == n:
        yield result

    for elem in it:
        result = result[1:] + elem
        yield result


def obtain_context(df, cont, cent):
    if cont == 5:
        return df['mer'][8:13]
    else: 
        return df['mer'][9:12]


def get_context_counts(df, context):
    return Counter(list(df[context]))


#normalize counts by the reference genome
def get_context_norm(ref, exp):
    new_d = {}
    for el in exp:
        if el != '-':
            try:
                norm = exp.get(el) / ref.get(el)
                new_d.update({el : norm})
            except:
                #May produce this exception at NNN or NNNNN contexts
                import pdb;pdb.set_trace()
    
    df = pd.DataFrame(new_d.items(), columns=['CONTEXT', 'REF_NORM'])
    df['TOTAL_NORM'] = df['REF_NORM'] / df['REF_NORM'].sum()

    return df


def generate_plots(df_t, df_p, td, pd, out_dir):
    pl.obtain_plots(df_t, td, 'triplet', 16)
    pl.obtain_plots(df_p, pd, 'pentamer', 256)


@click.command(short_help='Get enrichment')
@click.option('-d', '--data', required=True)
@click.option('-cl', '--chrom_len', required=True)
@click.option('-bf', '--base_filter', default='', multiple=True)
@click.option('-o', '--output', required=True)
def main(data, chrom_len, base_filter, output):
    df = pd.read_csv(data, sep='\t', compression='gzip')
    
    try:
        df.columns = names_all_1
    except:
        try:
            df.columns = names_all_2
        except:
            df.columns = names_all_3
            df.drop('motif_2 P-value', axis=1, inplace=True)

    if base_filter: 
        new_df = pd.DataFrame()
        for el in base_filter:
            new_df = pd.concat([new_df, df[df['base'] == el]])
        df = new_df
    
    #optional to take only those where the treated has higher frequency
    df['value'] = df['treated_freq'] - df['untreated_freq']
    df = df[df['value'] > 0]
    
    chrom_lens = pd.read_csv(
        chrom_len, sep='\t', names=['CHROM', 'LEN']
    )

    df, triplet_context, penta_context = get_contexts(df, chrom_lens)
    
    triplet_dir = os.path.join(output, 'triplet')
    penta_dir = os.path.join(output, 'pentamer')
    
    generate_plots(
        triplet_context, penta_context, triplet_dir, 
        penta_dir, output
    )

    #write output
    ut.write_outputs(
        df, triplet_context, penta_context, output, 'Nanopore'
    )
    

if __name__ == "__main__":
    main()