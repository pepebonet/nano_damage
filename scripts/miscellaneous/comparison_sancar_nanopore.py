#!/usr/bin/env python3
import os
import click
import pandas as pd
from scipy import spatial
import matplotlib.pyplot as plt

import sys
sys.path.append('./scripts/')
import utils as ut


def load_data(triplet_sancar, pentamer_sancar, triplet_nanopore, pentamer_nanopore):
    tn = pd.read_csv(triplet_nanopore, sep='\t')
    tn.columns = ['CONTEXT', 'REF_NORM_NANOPORE', 'Norm Nanopore']
    tm = pd.read_csv(triplet_sancar, sep='\t')
    tm.columns = ['CONTEXT', 'REF_NORM_sancar', 'Norm sancar']
    
    pn = pd.read_csv(pentamer_nanopore, sep='\t')
    pn.columns = ['CONTEXT', 'REF_NORM_NANOPORE', 'Norm Nanopore']
    pm = pd.read_csv(pentamer_sancar, sep='\t')
    pm.columns = ['CONTEXT', 'REF_NORM_sancar', 'Norm sancar']
    
    return tn, tm, pn, pm


#obtain comparison df joining sancar and nanopore data
def obtain_df(n, m):
    m = m[m['CONTEXT' ].str.contains('N') == False]
    m['Norm sancar'] = m['REF_NORM_sancar'] / m['REF_NORM_sancar'].sum()

    result = pd.merge(n, m, on='CONTEXT', how='outer').fillna(0)

    result['Nanopore - sancar'] = result['Norm Nanopore'] - result['Norm sancar']
    
    cos_sim = 1 - spatial.distance.cosine(
        result['Norm Nanopore'], result['Norm sancar'])
    
    return result.reset_index(), cos_sim


def order_plots(data, context, cent=None):
    if context == 'triplet':
        order_plot = ut.order_muts_triplet(cent)
    else:
        order_plot = ut.order_muts_pentamers(cent)

    df = data.set_index('CONTEXT')
    df = df.loc[order_plot].reset_index()

    return df[df['REF_NORM_sancar'] >= 0]


def do_plots(df, cos_sim, output, label_method, label):

    fig, ax = plt.subplots(figsize=(10, 5))
    df = order_plots(df, label)
    custom_lines = []
    
    plt.bar(df['CONTEXT'], df['Norm sancar'] * -1, 
        color='#6baed6', alpha=1, label='Negative sancar')
    plt.bar(df['CONTEXT'], df['Norm Nanopore'], 
        color='#2171b5', alpha=1, label='Nanopore')
    plt.bar(df['CONTEXT'], df['Nanopore - sancar'], 
        color='#08306b', alpha=1, label='Nanopore - sancar')
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    for el in [('Sancar Enrichment (Negative)', '#6baed6'), 
        (f'{label_method} Enrichment', '#2171b5'), 
        (f'{label_method} - Sancar', '#08306b')]:
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

    out_file1 = os.path.join(output, f'comp_sancar_{label_method}_{label}.pdf')
    out_file2 = os.path.join(output, f'comp_sancar_{label_method}_{label}.png')
    fig.tight_layout()
    plt.savefig(out_file1)
    plt.savefig(out_file2)
    plt.close()



@click.command(short_help='Get comparison Nanopore and sancar')
@click.option('-st', '--sancar_triplet', required=True)
@click.option('-sp', '--sancar_pentamer', required=True)
@click.option('-tn', '--triplet_nanopore', required=True)
@click.option('-pn', '--pentamer_nanopore', required=True)
@click.option(
    '-lm', '--label_method', default='Tombo', type=click.Choice(['Tombo', 'SemiSup'])
)
@click.option('-o', '--output', required=True)
def main(sancar_triplet, sancar_pentamer, triplet_nanopore, pentamer_nanopore, 
    label_method, output):

    tn, tm, pn, pm = load_data(
        sancar_triplet, sancar_pentamer, triplet_nanopore, pentamer_nanopore
    )

    triplet, cos_sim_tri = obtain_df(tn, tm)
    pentamer, cos_sim_pent = obtain_df(pn, pm)

    do_plots(triplet, cos_sim_tri, output, label_method, 'triplet')
    do_plots(pentamer, cos_sim_pent, output, label_method, 'pentamer')

    

if __name__ == '__main__':
    main()