#!/usr/bin/env python3
import os
import click
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter


def plot_norm_nucleosomes(mms, cisplatin, outdir):
    # import pdb;pdb.set_trace()
    x = mms['Position'].values.flatten()

    yhat_mms = mms['Relative Increase'].values.flatten()
    yhat_mms = savgol_filter(yhat_mms, 11, 3)

    yhat_cisplatin = cisplatin['Relative Increase'].values.flatten()
    yhat_cisplatin = savgol_filter(yhat_cisplatin, 11, 3)
    
    fig, ax = plt.subplots(figsize=(10, 4))

    plt.sca(ax)
    ax.plot(x, yhat_mms, linewidth=4, color='#a6cee3', alpha=0.8)
    ax.plot(x, yhat_cisplatin, linewidth=4, color='#08519c', alpha=0.8)

    ax.set_ylabel('Relative Increase',fontsize=14)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(-1, 148)
    ax.set_ylim(-0.3, 0.3)
    plt.setp([ax.get_xticklines(), ax.get_yticklines()], color='grey')

    ax.xaxis.set_ticks_position('none')
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(0.2)

    ax.xaxis.set_tick_params(pad=0.5)
    ax.yaxis.set_tick_params(pad=0.5, width=0.5)

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tick_params(axis='both', which='both', bottom=False, left = False)
    ax.set_xlabel('Position', fontsize=14)

    custom_lines = []

    custom_lines.append(
        plt.plot([],[], marker="o", ms=10, ls="", mec='black', 
        mew=0, color='#a6cee3', label='MMS Damage')[0] 
    )

    custom_lines.append(
        plt.plot([],[], marker="o", ms=10, ls="", mec='black', 
        mew=0, color='#08519c', label='Cisplatin Damage')[0] 
    )

    ax.legend(
        bbox_to_anchor=(0., 0.9, 1., .102),
        handles=custom_lines, loc='upper center', 
        facecolor='white', ncol=1, fontsize=12, frameon=False
    )
    
    fig.tight_layout()

    plt.savefig(outdir)
    plt.close()

# ------------------------------------------------------------------------------
# CLICK
# ------------------------------------------------------------------------------

@click.command(short_help='script to plot damage differences '
'between cisplatin and mms')
@click.option(
    '-dm', '--damage_mms', default='', 
    help='data containing the damage from the mms experiment'
)
@click.option(
    '-dc', '--damage_cisplatin', default='', 
    help='data containing the damage from the cisplatin experiment '
)
@click.option(
    '-o', '--output', default='', help='output folder'
)
def main(damage_mms, damage_cisplatin, output):
    #Obtain data
    mms = pd.read_csv(damage_mms, sep='\t')
    cisplatin = pd.read_csv(damage_cisplatin, sep='\t')

    fig_name = os.path.join(output, 'nuc_mms_vs_cisplatin')
    plot_norm_nucleosomes(mms, cisplatin, fig_name)


if __name__ == '__main__':
    main()