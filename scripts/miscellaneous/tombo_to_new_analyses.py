#!/usr/bin/envs python3

import os
import click
import pandas as pd


# ------------------------------------------------------------------------------
# CLICK
# ------------------------------------------------------------------------------

@click.command(short_help='script to get the nucleosome damage')
@click.option(
    '-nd', '--new_damage', default='', 
    help='dataframe as an example of the new damage structure from novoa'
)
@click.option(
    '-td', '--tombo_damage', default='', 
    help='damage dataframe from old tombo analysis'
)
@click.option(
    '-o', '--output', default='', help='output folder'
)
def main(new_damage, tombo_damage, output):
    new = pd.read_csv(new_damage, sep='\t')
    tombo = pd.read_csv(tombo_damage, sep='\t')

    tombo['untreated_freq'] = 0
    tombo['treated_freq'] = 0
    tombo['min_coverage'] = 0
    tombo['mer'] = '-'
    tombo.drop(['diff'], inplace=True, axis=1)

    tombo = tombo.rename(
        columns={'CHROM': 'chrom', 'p-value':'value', 'SEQ':'base'}
    )
    
    tombo = tombo[new.columns]

    out_file = os.path.join(output, 'damage_most_sig_Nanopore.tsv')
    tombo.to_csv(out_file, sep='\t', index=None)


if __name__ == '__main__':
    main()