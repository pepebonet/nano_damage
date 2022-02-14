#!/usr/bin/envs/ python3
import os
import click 
import pandas as pd
import matplotlib.pyplot as plt


palette_dict = { 
        'm7G':'#e62725',
        'm3A':'#1ebff0'
}

def plot_rep_damage(df, output):
    fig, ax = plt.subplots(figsize=(5, 5), facecolor='white')
    custom_lines = []

    for dam in ['m7G', 'm3A']:
        vals = df[df['Residue'] == dam]

        plt.errorbar([dam], vals['Representation of Damage'].unique() * 100,
            marker='o', mfc=palette_dict[dam], mec='black', ms=8, 
            mew=1,  ecolor='black', capsize=1, capthick=0.5, elinewidth=1, 
            c=palette_dict[dam]
        )

        custom_lines.append(
            plt.plot([],[], marker="o", ms=7, ls="", mec='black', 
            mew=0, color=palette_dict[dam], label=vals['Residue'].unique()[0])[0] 
        )
    

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.legend().set_visible(False)
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(0, 100)
    ax.legend(
        bbox_to_anchor=(0., 0.9, 1., .102),
        handles=custom_lines, loc='upper right', 
        facecolor='white', ncol=1, fontsize=8, frameon=False
    )

    ax.set_xlabel("", fontsize=12)
    ax.set_ylabel("Percentage, %", fontsize=12)

    out_fig = os.path.join(output, 'representation_damage_plot.pdf')
    plt.tight_layout()
    plt.savefig(out_fig)
    plt.close()


def plot_prop_damage(df, output):
    
    fig, ax = plt.subplots(figsize=(5, 5), facecolor='white')
    custom_lines = []

    for dam in ['m7G', 'm3A']:
        vals = df[df['Residue'] == dam]
        std_vals = vals['Damage Proportion'].std()

        plt.errorbar([dam], vals['Mean'].unique(), yerr=std_vals,
            marker='o', mfc=palette_dict[dam], mec='black', ms=8, 
            mew=1,  ecolor='black', capsize=1, capthick=0.5, elinewidth=1, 
            c=palette_dict[dam]
        )

        custom_lines.append(
            plt.plot([],[], marker="o", ms=7, ls="", mec='black', 
            mew=0, color=palette_dict[dam], label=vals['Residue'].unique()[0])[0] 
        )
    

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.legend().set_visible(False)
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(0, 5.6)
    
    ax.legend(
        bbox_to_anchor=(0., 0.9, 1., .102),
        handles=custom_lines, loc='upper right', 
        facecolor='white', ncol=1, fontsize=8, frameon=False
    )

    ax.set_xlabel("", fontsize=12)
    ax.set_ylabel("Percentage, %", fontsize=12)

    out_fig = os.path.join(output, 'proportion_damage_plot.pdf')
    plt.tight_layout()
    plt.savefig(out_fig)
    plt.close()
    

@click.command(short_help='Script to get figures from mass spec')
@click.option('-df', '--dataframe', required=True)
@click.option('-o', '--output', required=True)
def main(dataframe, output):

    df = pd.read_csv(dataframe, sep='\t')
    plot_prop_damage(df, output)
    plot_rep_damage(df, output)


if __name__ == '__main__':
    main()
