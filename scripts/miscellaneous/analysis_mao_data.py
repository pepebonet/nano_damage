#!/usr/bin/env python3
import os 
import sys
import click
import pandas as pd

sys.path.append('./scripts')
import utils as ut
import plots as pl


def get_base_percentage(df):
    only_met =  df[df['READS'] != 0]

    q1 = only_met['READS'].quantile([0.90]).tolist()
    new_df = only_met[only_met['READS'] > q1[0]]

    return new_df


def load_data(plus_s, minus_s, chrom_len): 

    B_p = pd.read_csv(plus_s, sep='\t', \
        names=['CHROM', 'P0', 'pos', 'ID', 'READS'])
    B_p['strand'] = '+'
    B_m = pd.read_csv(minus_s, sep='\t', \
        names=['CHROM', 'P0', 'pos', 'ID', 'READS'])
    B_m['strand'] = '-'
    B = pd.concat([B_p, B_m])

    chrom_len = pd.read_csv(
        chrom_len, sep='\t', names=['CHROM', 'LEN']
    )

    return B, chrom_len


@click.command(short_help='Get data from Mao experiment')
@click.option('-ms', '--minus_strand', required=True)
@click.option('-ps', '--plus_strand', required=True)
@click.option('-cl', '--chrom_len', required=True)
@click.option('-o', '--output', required=True)
def main(minus_strand, plus_strand, chrom_len, output):

    df, chrom_len = load_data(plus_strand, minus_strand, chrom_len)

    penta_ref, triplet_ref = ut.counts_reference_genome(chrom_len)

    df['SEQ'] = df.apply(ut.annot, axis = 1)

    df = get_base_percentage(df)

    df['PENTAMER'] = df.apply(ut.obtain_contex, args=(5,2), axis=1)
    df['TRIPLET'] = df.apply(ut.obtain_contex, args=(3,1), axis=1)

    penta_exp = ut.get_context_counts(df, 'PENTAMER')
    triplet_exp = ut.get_context_counts(df, 'TRIPLET')

    triplet_context = ut.get_context_norm(triplet_ref, triplet_exp)
    penta_context = ut.get_context_norm(penta_ref, penta_exp)

    triplet_dir = os.path.join(output, 'triplet')
    penta_dir = os.path.join(output, 'pentamer')

    #Obtain plots
    pl.obtain_plots(triplet_context, triplet_dir, 'triplet', 16)
    pl.obtain_plots(penta_context, penta_dir, 'pentamer', 256)

    ut.write_outputs(
        df, triplet_context, penta_context, output, 'Nanopore'
    )
    
    


if __name__ == '__main__':
    main()