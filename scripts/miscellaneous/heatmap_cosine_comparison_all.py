#!/usr/bin/env python3

import click
import numpy as np
import pandas as pd
from scipy import spatial
import seaborn as sns

import sys
sys.path.append('./scripts/')
import utils as ut

names = ['Novoa Cisplatin', 'Novoa MMS', 'Tombo MMS', 'Tombo Cisplatin',
        'Mao MMS', 'Sancar Cisplatin']


def get_cosines(nc, nm, tm, tc, mm, sc):

    df = pd.DataFrame(columns=names, index=names)

    for el1 in [(nc, 'Novoa Cisplatin'), (nm, 'Novoa MMS'), 
        (tm, 'Tombo MMS'), (tc, 'Tombo Cisplatin'), (mm, 'Mao MMS'), 
            (sc, 'Sancar Cisplatin')]:
        for el2 in  [(nc, 'Novoa Cisplatin'), (nm, 'Novoa MMS'), 
            (tm, 'Tombo MMS'), (tc, 'Tombo Cisplatin'), (mm, 'Mao MMS'), 
                (sc, 'Sancar Cisplatin')]:
            cos_sim = round(1 - spatial.distance.cosine(
                el1[0]['TOTAL_NORM'].tolist(), el2[0]['TOTAL_NORM'].tolist()), 2)
            print(f'Cosine {el1[1]} - {el2[1]}: {cos_sim}')

            df[el1[1]][el2[1]] = cos_sim
    return df


def fix_data(df, good_df):

    aa = list(set(df['CONTEXT']) ^ set(good_df['CONTEXT']))
    bb = pd.DataFrame([aa, [0]*len(aa), [0]*len(aa)]).T
    bb.set_axis(list(df.columns), axis=1, inplace=True)
    return pd.concat([df, bb], axis=0).reset_index(drop=True).sort_values('CONTEXT')


def fix_sancar(df):
    return df[ df['CONTEXT' ].str.contains('N') == False]


def load_data(novoa_cis, novoa_mms, tombo_mms, tombo_cis, mao_mms, sancar_cis):
    tm = pd.read_csv(tombo_mms, sep='\t').sort_values('CONTEXT')
    tc = pd.read_csv(tombo_cis, sep='\t').sort_values('CONTEXT')

    nc = fix_data(pd.read_csv(novoa_cis, sep='\t'), tm).sort_values('CONTEXT')
    nm = fix_data(pd.read_csv(novoa_mms, sep='\t'), tm).sort_values('CONTEXT')
    
    mm = fix_data(pd.read_csv(mao_mms, sep='\t'), tm)
    sc = fix_sancar(pd.read_csv(sancar_cis, sep='\t')).sort_values('CONTEXT')

    return nc, nm, tm, tc, mm, sc


def plot_heatmap(df, output):

    corr = np.corrcoef(df.values)
    mask = np.zeros_like(corr)
    mask[np.triu_indices_from(mask)] = True
    with sns.axes_style("white"):
        f, ax = plt.subplots(figsize=(7, 5))
        ax = sns.heatmap(corr, mask=mask, vmax=.3, square=True)


@click.command(short_help='Get comparison all and heatmap')
@click.option('-nm', '--novoa_mms', required=True)
@click.option('-nc', '--novoa_cis', required=True)
@click.option('-mm', '--mao_mms', required=True)
@click.option('-sc', '--sancar_cis', required=True)
@click.option('-tm', '--tombo_mms', required=True)
@click.option('-tc', '--tombo_cis', required=True)
@click.option('-o', '--output', required=True)
def main(novoa_cis, novoa_mms, mao_mms, sancar_cis, tombo_mms, 
    tombo_cis, output):

    nc, nm, tm, tc, mm, sc = load_data(
        novoa_cis, novoa_mms, tombo_mms, tombo_cis, mao_mms, sancar_cis
    )
     
    df = get_cosines(nc, nm, tm, tc, mm, sc)


    plot_heatmap(df, output)

    

if __name__ == '__main__':
    main()