#!/usr/bin/envs python3

import click
import pandas as pd

from bgreference import refseq



#Obtain base for the position in the file
def annot(df):
    seq_pos = refseq('saccer3', df['chrom'], df['Start'], 1)
    return seq_pos


# ------------------------------------------------------------------------------
# CLICK
# ------------------------------------------------------------------------------

@click.command(short_help='script to get the strand information')
@click.option(
    '-da', '--damage_path', default='', help='damaged bases'
)
def main(damage_path):
    damage = pd.read_csv(damage_path, sep='\t')
    print(damage.iloc[2])
    print(refseq('saccer3', damage.iloc[2]['chrom'], damage.iloc[2]['pos'] - 10, 21))
    print('you clearly have the negative strand sequence here')
    import pdb;pdb.set_trace()


if __name__ == '__main__':
    main()