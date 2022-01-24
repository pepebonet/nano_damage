#!/usr/bin/env python3
import os
import click
import pandas as pd
from scipy import spatial


def load_data(triplet_mao, pentamer_mao, triplet_nanopore, pentamer_nanopore):
    tn = pd.read_csv(triplet_nanopore, sep='\t')
    tn.columns = ['CONTEXT', 'REF_NORM_NANOPORE', 'Norm Nanopore']
    tm = pd.read_csv(triplet_mao, sep='\t')
    tm.columns = ['CONTEXT', 'REF_NORM_MAO', 'Norm Mao']
    
    pn = pd.read_csv(pentamer_nanopore, sep='\t')
    pn.columns = ['CONTEXT', 'REF_NORM_NANOPORE', 'Norm Nanopore']
    pm = pd.read_csv(pentamer_mao, sep='\t')
    pm.columns = ['CONTEXT', 'REF_NORM_MAO', 'Norm Mao']
    
    return tn, tm, pn, pm


#obtain comparison df joining mao and nanopore data
def obtain_df(n, m):
    
    n = n.set_index('CONTEXT')
    m = m.set_index('CONTEXT')
    result = pd.concat([m, n], axis=1, join='inner')

    result['Nanopore - Mao'] = result['Norm Nanopore'] - result['Norm Mao']
    
    cos_sim = 1 - spatial.distance.cosine(
        result['Norm Nanopore'], result['Norm Mao'])
    
    return result.reset_index(), cos_sim



@click.command(short_help='Get comparison Nanopore and Mao')
@click.option('-mt', '--mao_triplet', required=True)
@click.option('-mp', '--mao_pentamer', required=True)
@click.option('-tn', '--triplet_nanopore', required=True)
@click.option('-pn', '--pentamer_nanopore', required=True)
@click.option('-o', '--output', required=True)
def main(mao_triplet, mao_pentamer, triplet_nanopore, pentamer_nanopore, output):

    tn, tm, pn, pm = load_data(
        mao_triplet, mao_pentamer, triplet_nanopore, pentamer_nanopore
    )
    
    triplet, cos_sim_tri = obtain_df(tn, tm)
    pentamer, cos_sim_pent = obtain_df(pn, pm)
    import pdb;pdb.set_trace()
    



if __name__ == '__main__':
    main()