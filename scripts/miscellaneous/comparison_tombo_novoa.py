#!/usr/bin/env python3
import os
import click
import pandas as pd
from scipy import spatial
import matplotlib.pyplot as plt

import sys
sys.path.append('./scripts/')
import utils as ut


def load_data(triplet_novoa, pentamer_novoa, triplet_tombo, pentamer_tombo):
    tt = pd.read_csv(triplet_tombo, sep='\t')
    tt.columns = ['CONTEXT', 'REF_NORM_TOMBO', 'Norm Tombo']
    tn = pd.read_csv(triplet_novoa, sep='\t')
    tn.columns = ['CONTEXT', 'REF_NORM_NOVOA', 'Norm Novoa']
    
    pt = pd.read_csv(pentamer_tombo, sep='\t')
    pt.columns = ['CONTEXT', 'REF_NORM_TOMBO', 'Norm Tombo']
    pn = pd.read_csv(pentamer_novoa, sep='\t')
    pn.columns = ['CONTEXT', 'REF_NORM_NOVOA', 'Norm Novoa']
    
    return tt, tn, pt, pn


#obtain comparison df joining mao and nanopore data
def obtain_df(n, m):

    result = pd.merge(n, m, on='CONTEXT', how='outer').fillna(0)

    result['Tombo - Novoa'] = result['Norm Tombo'] - result['Norm Novoa']
    
    cos_sim = 1 - spatial.distance.cosine(
        result['Norm Tombo'], result['Norm Novoa'])
    
    return result.reset_index(), cos_sim


def order_plots(data, context, cent=None):
    if context == 'triplet':
        order_plot = ut.order_muts_triplet(cent)
    else:
        order_plot = ut.order_muts_pentamers(cent)

    df = data.set_index('CONTEXT')
    df = df.loc[order_plot].reset_index()

    return df[df['REF_NORM_TOMBO'] >= 0]


def do_plots(df, cos_sim, output, label):

    fig, ax = plt.subplots(figsize=(10, 5))
    df = order_plots(df, label)
    custom_lines = []
    
    plt.bar(df['CONTEXT'], df['Norm Novoa'] * -1, 
        color='#6baed6', alpha=1, label='Negative Novoa')
    plt.bar(df['CONTEXT'], df['Norm Tombo'], 
        color='#2171b5', alpha=1, label='Tombo')
    plt.bar(df['CONTEXT'], df['Tombo - Novoa'], 
        color='#08306b', alpha=1, label='Tombo - Novoa')
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    for el in [('SemiSup Enrichment (Negative)', '#6baed6'), 
        ('Tombo Enrichment', '#2171b5'), ('Tombo - SemiSup', '#08306b')]:
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

    out_file = os.path.join(output, 'comp_Novoa_Tombo_{}.pdf'.format(label))
    fig.tight_layout()
    plt.savefig(out_file)
    plt.close()



@click.command(short_help='Get comparison Tombo and Novoas method')
@click.option('-tn', '--triplet_novoa', required=True)
@click.option('-pn', '--pentamer_novoa', required=True)
@click.option('-tt', '--triplet_tombo', required=True)
@click.option('-pt', '--pentamer_tombo', required=True)
@click.option('-o', '--output', required=True)
def main(triplet_novoa, pentamer_novoa, triplet_tombo, pentamer_tombo, output):

    tt, tn, pt, pn = load_data(
        triplet_novoa, pentamer_novoa, triplet_tombo, pentamer_tombo
    )

    triplet, cos_sim_tri = obtain_df(tt, tn)
    pentamer, cos_sim_pent = obtain_df(pt, pn)

    do_plots(triplet, cos_sim_tri, output, 'triplet')
    do_plots(pentamer, cos_sim_pent, output, 'pentamer')

    

if __name__ == '__main__':
    main()