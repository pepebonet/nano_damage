#!/usr/bin/env python3
import os 
import sys
import click
import pandas as pd
from tqdm import tqdm
from collections import Counter

from bgreference import hg19

sys.path.append('./scripts')
import utils as ut
import plots as pl

col_names = ['CHROM', 'Start', 'pos', 'SEQ', 'strand']


#Obtain context for every base
def obtain_contex(df, cont, cent):
    try:
        seq = hg19(df['CHROM'], df['pos'] - cent, cont)
        
        if len(seq) < cont:
            return '-'
        else:
            if df['strand'] == '-':
                return ut.comp_seq(seq)
            else:
                return seq

    except:
        seq = '-'
    return seq


def counts_reference_genome(chrom_len):
    from Bio import SeqIO
    all_penta = {}; all_triplet = {}
    for file in tqdm(os.listdir(chrom_len)):
        fastas = SeqIO.parse(open(os.path.join(chrom_len, file)),'fasta')

        for fasta in fastas:
            _, sequence = fasta.id, str(fasta.seq)
            seq = sequence.upper()
            
            penta = Counter(list(ut.slicing_window(seq, 5)))
            triplet = Counter(list(ut.slicing_window(seq, 3)))

            all_triplet = dict(Counter(triplet) + Counter(all_triplet))
            all_penta = dict(Counter(penta) + Counter(all_penta))

    return all_penta, all_triplet


def load_data(rep1_path, rep2_path):
    
    rep1 = pd.read_csv(rep1_path, sep='\t', header=None, names=col_names)
    rep1['ID'] = rep1['CHROM'] + '_' + rep1['pos'].astype(str) + '_' + rep1['strand']

    rep2 = pd.read_csv(rep2_path, sep='\t', header=None, names=col_names)
    rep2['ID'] = rep2['CHROM'] + '_' + rep2['pos'].astype(str) + '_' + rep2['strand']

    return pd.concat([rep1, rep2], axis=0).drop_duplicates()



@click.command(short_help='Get data from Sancar experiment')
@click.option('-r1', '--replicate_1', required=True)
@click.option('-r2', '--replicate_2', required=True)
@click.option('-cl', '--chrom_len', required=True)
@click.option('-o', '--output', required=True)
def main(replicate_1, replicate_2, chrom_len, output):

    df = load_data(replicate_2, replicate_1)

    penta_ref, triplet_ref = counts_reference_genome(chrom_len)

    df['PENTAMER'] = df.apply(obtain_contex, args=(5,2), axis=1)
    df['TRIPLET'] = df.apply(obtain_contex, args=(3,1), axis=1)

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