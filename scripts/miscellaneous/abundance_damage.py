#!/usr/bin/envs python3
import click
import pandas as pd


def get_triplet_info(path):
    df = pd.read_csv(path, sep='\t')
    
    try:
        df[['0', '1', '2', '3', '4']] = df.CONTEXT.str.split('', expand=True)
    except:
        df[['1', '2', '3']] = pd.DataFrame(df['CONTEXT'].apply(lambda x: list(x)).tolist())

    df_gs = df[df['2'] == 'G']
    prop_gs = round(df_gs['TOTAL_NORM'].sum(), 2)
    df_as = df[df['2'] == 'A']
    prop_as = round(df_as['TOTAL_NORM'].sum(), 2)
    print('Triplet info:')
    print(f'Proportion of Gs: {prop_gs}')
    print(f'Proportion of As: {prop_as}')


def get_pentamer_info(path):
    df = pd.read_csv(path, sep='\t')
    try:
        df[['-1', '0', '1', '2', '3', '4', '5']] = df.CONTEXT.str.split('', expand=True)
    except:
        df[['0', '1', '2', '3', '4']] = pd.DataFrame(df['CONTEXT'].apply(lambda x: list(x)).tolist())
    
    df_gs = df[df['2'] == 'G']
    prop_gs = round(df_gs['TOTAL_NORM'].sum(), 2)
    df_as = df[df['2'] == 'A']
    prop_as = round(df_as['TOTAL_NORM'].sum(), 2)

    print('Pentamer info:')
    print(f'Proportion of Gs: {prop_gs}')
    print(f'Proportion of As: {prop_as}')


@click.command(short_help='Get stats from cosines')
@click.option('-td', '--triplet_data', required=True)
@click.option('-pd', '--pentamer_data', required=True)
def main(triplet_data, pentamer_data):

    get_triplet_info(triplet_data)
    get_pentamer_info(pentamer_data)
    


if __name__ == '__main__':
    main()

