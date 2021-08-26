#!/usr/bin/envs/ python3

import os
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


def get_non_overlapping_pos(df1, df2):
    """ Get positions of damage not overlapping
    Args:
        df1: dataset of the first batch
        df2: dataset of the second batch
    Returns:
        df: dataset containing non-overlapping positions
    """

    non_overlap = list(set(df1['ID'].tolist()) ^ set(df2['ID'].tolist()))

    df1_non_overlap = df1[df1['ID'].isin(non_overlap)]
    df2_non_overlap = df2[df2['ID'].isin(non_overlap)]

    df = pd.concat([df1_non_overlap, df2_non_overlap])

    return df.drop(['ID'], axis=1)



def get_overlapping_pos(df1, df2):
    """ Get positions of damage overlapping
    Args:
        df1: dataset of the first batch
        df2: dataset of the second batch
    Returns:
        df: dataset containing overlapping positions
    """
    
    overlap = pd.merge(df1, df2, on='ID', how='inner')

    df = overlap[overlap.columns[:6]]
    df.columns = names_all_1[:6]

    df['untreated_freq'] = (overlap['untreated_freq_x'] + \
        overlap['untreated_freq_y']) / 2
    
    df['treated_freq'] = (overlap['treated_freq_x'] + \
        overlap['treated_freq_y']) / 2

    return df


@click.command(short_help='Script to merge datasets from experiment 2')
@click.option('-do', '--data_one', required=True, help='first dataset to merge')
@click.option('-dt', '--data_two', required=True, help='second dataset to merge')
@click.option('-o', '--output', default='', help='output path')
def main(data_one, data_two, output):

    df1, df2 = load_datasets(data_one, data_two)

    df_non_overlap = get_non_overlapping_pos(df1, df2)
    
    df_overlap = get_overlapping_pos(df1, df2)

    df_final = pd.concat([df_non_overlap, df_overlap])

    filename = data_one.rsplit('/', 1)[1]
    out_file = os.path.join(output, filename)

    df_final.to_csv(out_file, sep='\t', index=None, compression='gzip')


if __name__ == '__main__':
    main()