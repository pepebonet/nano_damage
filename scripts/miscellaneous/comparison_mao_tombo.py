#!/usr/bin/env python3
import os
import click
import pandas as pd
from scipy import spatial
import matplotlib.pyplot as plt

import sys
sys.path.append('./scripts/')
import utils as ut


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

    result = pd.merge(n, m, on='CONTEXT', how='outer').fillna(0)

    result['Nanopore - Mao'] = result['Norm Nanopore'] - result['Norm Mao']
    
    cos_sim = 1 - spatial.distance.cosine(
        result['Norm Nanopore'], result['Norm Mao'])
    
    return result.reset_index(), cos_sim


def order_plots(data, context, cent=None):
    if context == 'triplet':
        order_plot = ut.order_muts_triplet(cent)
    else:
        order_plot = ut.order_muts_pentamers(cent)

    df = data.set_index('CONTEXT')
    df = df.loc[order_plot].reset_index()

    return df[df['REF_NORM_MAO'] >= 0]


def do_plots(df, cos_sim, output, label):

    fig, ax = plt.subplots(figsize=(10, 5))
    df = order_plots(df, label)
    custom_lines = []
    
    plt.bar(df['CONTEXT'], df['Norm Mao'] * -1, 
        color='#59abff', alpha=0.3, label='Negative Mao')
    plt.bar(df['CONTEXT'], df['Norm Nanopore'], 
        color='lightblue', alpha=0.7, label='Nanopore')
    plt.bar(df['CONTEXT'], df['Nanopore - Mao'], 
        color='#16416e', alpha=0.9, label='Nanopore - Mao')
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    for el in [('Mao Enrichment (Negative)', '#59abff'), 
        ('Nanopore Enrichment', 'lightblue'), ('Nanopore - Mao', '#16416e')]:
        custom_lines.append(
            plt.plot([],[], marker="o", ms=8, ls="", mec='black', 
            mew=0, color=el[1], label=el[0])[0] 
        )

    ax.legend()
    if label == 'triplet':
        plt.xticks(fontsize=12, rotation=90)
    else:
        ax.get_xaxis().set_visible(False)

    plt.hlines(y=0, xmin=-0.5, xmax=len(df) - 0.5, color='black', 
        alpha=0.8, ls='dashed')

    plt.text(0.5, 1.1, 'Cosine Similariy: {}'.format(round(cos_sim, 2)), 
        horizontalalignment='center', verticalalignment='center', 
        transform=ax.transAxes, fontsize=14)

    ax.legend(
        bbox_to_anchor=(0., 1.0, 1., .102),
        handles=custom_lines, loc='upper right', 
        facecolor='white', ncol=1, fontsize=10, frameon=False
    )

    out_file = os.path.join(output, 'comp_mao_nanopore_{}.pdf'.format(label))
    fig.tight_layout()
    plt.savefig(out_file)
    plt.close()



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

    do_plots(triplet, cos_sim_tri, output, 'triplet')
    do_plots(pentamer, cos_sim_pent, output, 'pentamer')

    

if __name__ == '__main__':
    main()