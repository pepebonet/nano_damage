#!/usr/bin/envs/ python3

import click
import pandas as pd

names_all_1 = ['chrom', 'pos', 'base', 'strand', 'mer', \
    'min_coverage', 'untreated_freq', 'treated_freq']

names_all_2 = ['index', 'chrom', 'pos', 'base', 'strand', 'mer', 'min coverage',
       'untreated_freq', 'treated_freq', 'diff', 'group',
       'motif_1 DRWGGDD P-value', 'motif']


def load_datasets(data_one, data_two):
    """ Load datasets 
    Args:
        data_one: dataset of the first batch
        data_two: dataset of the second batch
    Returns:
        df1, df2: both datasets to merge
    """

    df1 = pd.read_csv(data_one, sep='\t', compression='gzip')
    df2 = pd.read_csv(data_two, sep='\t', compression='gzip')

    try:
        df1.columns = names_all_1
        df2.columns = names_all_1
    except:
        df1.columns = names_all_2
        df2.columns = names_all_2

    df1['ID'] = df1['chrom'] + '_' + df1['pos'].astype(str)
    df2['ID'] = df2['chrom'] + '_' + df2['pos'].astype(str)

    return df1, df2


@click.command(short_help='Script to merge datasets from experiment 2')
@click.option('-do', '--data_one', required=True, help='first dataset to merge')
@click.option('-dt', '--data_two', required=True, help='second dataset to merge')
@click.option('-o', '--output', default='', help='output path')
def main(data_one, data_two, output):

    df1, df2 = load_datasets(data_one, data_two)

    

    import pdb;pdb.set_trace()


if __name__ == '__main__':
    main()