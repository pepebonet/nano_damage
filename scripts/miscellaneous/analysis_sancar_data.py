#!/usr/bin/env python3
import os 
import sys
import click
import pandas as pd

sys.path.append('./scripts')
import utils as ut
import plots as pl

col_names = ['CHROM', 'Start', 'pos', 'SEQ', 'strand']

def load_data(rep1_path, rep2_path, chrom_len):
    
    rep1 = pd.read_csv(rep1_path, sep='\t', header=None, names=col_names)
    rep1['ID'] = rep1['CHROM'] + '_' + rep1['pos'].astype(str) + '_' + rep1['strand']

    rep2 = pd.read_csv(rep2_path, sep='\t', header=None, names=col_names)
    rep2['ID'] = rep2['CHROM'] + '_' + rep2['pos'].astype(str) + '_' + rep2['strand']

    data = pd.merge(rep1, rep2, how='outer', on='ID')
    import pdb;pdb.set_trace()
    chrom_len = pd.read_csv(
        chrom_len, sep='\t', names=['CHROM', 'LEN']
    )
    
    return data, chrom_len


@click.command(short_help='Get data from Sancar experiment')
@click.option('-r1', '--replicate_1', required=True)
@click.option('-r2', '--replicate_2', required=True)
@click.option('-cl', '--chrom_len', required=True)
@click.option('-o', '--output', required=True)
def main(replicate_1, replicate_2, chrom_len, output):

    df, chrom_len = load_data(replicate_2, replicate_1, chrom_len)

    penta_ref, triplet_ref = ut.counts_reference_genome(chrom_len)
    import pdb;pdb.set_trace()
    # df['SEQ'] = df.apply(ut.annot, axis = 1)

    df['PENTAMER'] = df.apply(ut.obtain_contex, args=(5,2), axis=1)
    df['TRIPLET'] = df.apply(ut.obtain_contex, args=(3,1), axis=1)

    import pdb;pdb.set_trace()
    penta_exp = ut.get_context_counts(df, 'PENTAMER')
    triplet_exp = ut.get_context_counts(df, 'TRIPLET')

    triplet_context = ut.get_context_norm(triplet_ref, triplet_exp)
    penta_context = ut.get_context_norm(penta_ref, penta_exp)

    triplet_dir = os.path.join(output, 'triplet')
    penta_dir = os.path.join(output, 'pentamer')
    import pdb;pdb.set_trace()
    #Obtain plots
    pl.obtain_plots(triplet_context, triplet_dir, 'triplet', 16)
    pl.obtain_plots(penta_context, penta_dir, 'pentamer', 256)

    ut.write_outputs(
        df, triplet_context, penta_context, output, 'Nanopore'
    )
    

if __name__ == '__main__':
    main()