#!/usr/bin/env python3
import os
import click
import pandas as pd
import matplotlib.pyplot as plt

import sys
sys.path.append('./scripts/')
import utils as ut



def do_plots(t, n, output):
    fig, ax = plt.subplots(figsize=(7, 4))
    
    nt = n[n['Context'] == 'Triplet']
    np = n[n['Context'] == 'Pentamer']
    tt = t[t['Context'] == 'Triplet']
    tp = t[t['Context'] == 'Pentamer']

    x = list(range(len(nt['DNA structure'])))

    plt.errorbar(
        [a + 0.075 for a in x], np['Cosine-similarity'], 
        ls='none', marker='o', mfc='#6baed6', mec='black', ms=10, mew=1, 
        ecolor='#6baed6', capsize=2.5, elinewidth=0.7, capthick=0.7
    )

    plt.errorbar(
        [a + 0.225 for a in x], nt['Cosine-similarity'], 
        ls='none', marker='o', mfc='#08519c', mec='black', ms=10, mew=1, 
        ecolor='#08519c', capsize=2.5, elinewidth=0.7, capthick=0.7
    )

    plt.errorbar(
        [a - 0.075 for a in x], tt['Cosine-similarity'], 
        ls='none', marker='o', mfc='#006d2c', mec='black', ms=10, mew=1, 
        ecolor='#006d2c', capsize=2.5, elinewidth=0.7, capthick=0.7
    )

    plt.errorbar(
        [a - 0.225 for a in x], tp['Cosine-similarity'], 
        ls='none', marker='o', mfc='#74c476', mec='black', ms=10, mew=1, 
        ecolor='#74c476', capsize=2.5, elinewidth=0.7, capthick=0.7
    )

    custom_lines = []
    for el in [('Pentamer Tombo', '#74c476'), ('Triplet Tombo', '#006d2c'), 
        ('Pentamer SemiSup', '#6baed6'), ('Triplet SemiSup', '#08519c')]:
        custom_lines.append(
            plt.plot([],[], marker="o", ms=8, ls="", mec='black', 
            mew=0, color=el[1], label=el[0])[0] 
        )

    plt.ylim(round(t['Cosine-similarity'].min() - 0.04, 2), 1.03)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_xlabel("Genomic regions", fontsize=12)
    ax.set_ylabel("Cosine similarity", fontsize=12)
    ax.set_xticklabels([''] + nt['DNA structure'].tolist())
    from matplotlib.offsetbox import AnchoredText
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    title_text = 'Similarity on Genomic Regions'

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top", size="7.5%", pad=0)
    cax.get_xaxis().set_visible(False)
    cax.get_yaxis().set_visible(False)

    cax.set_facecolor('#c7cbd4')
    at = AnchoredText(title_text, loc=10, frameon=False,
            prop=dict(backgroundcolor='#c7cbd4',
                    size=8, color='black'))
    cax.add_artist(at)

    ax.legend(
        bbox_to_anchor=(0., 0.3, 1.0, .102),
        handles=custom_lines, loc='upper right', 
        facecolor='white', ncol=1, fontsize=10, frameon=False
    )

    fig.tight_layout()
    out_file = os.path.join(output, 'cosine_similarity_novao_tombo.pdf')
    plt.savefig(out_file)
    plt.close()


@click.command(short_help='Comparison Tombo and Novoas method for chromatin features')
@click.option('-td', '--tombo_data', required=True)
@click.option('-nd', '--novoa_data', required=True)
@click.option('-o', '--output', required=True)
def main(tombo_data, novoa_data, output):

    tombo = pd.read_csv(tombo_data, sep='\t') 
    novoa = pd.read_csv(novoa_data, sep='\t') 

    do_plots(tombo, novoa, output)



    

if __name__ == '__main__':
    main()