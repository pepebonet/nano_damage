#!/usr/bin/env python3
import os
import click
import numpy as np
import pandas as pd
from scipy import spatial
import matplotlib.pyplot as plt
import seaborn as sns

import sys
sys.path.append('./scripts/')
import utils as ut

names = ['SemiSup Cisplatin', 'SemiSup MMS', 'Tombo MMS', 'Tombo Cisplatin',
        'Mao MMS', 'Sancar Cisplatin']

names_extended = ['SemiSup Cisplatin', 'SemiSup Cisplatin Mao', 'SemiSup MMS', 
    'SemiSup MMS Mao', 'Tombo MMS', 'Tombo MMS Mao', 'Tombo Cisplatin',
    'Tombo Cisplatin Mao', 'Mao MMS', 'Sancar Cisplatin', 'Sancar Cisplatin Mao']


def get_cosines(nc, nm, tm, tc, mm, sc):

    df = pd.DataFrame(columns=names, index=names, dtype='float64')

    for el1 in [(nc, 'SemiSup Cisplatin'), (nm, 'SemiSup MMS'), 
        (tm, 'Tombo MMS'), (tc, 'Tombo Cisplatin'), (mm, 'Mao MMS'), 
            (sc, 'Sancar Cisplatin')]:
        for el2 in  [(nc, 'SemiSup Cisplatin'), (nm, 'SemiSup MMS'), 
            (tm, 'Tombo MMS'), (tc, 'Tombo Cisplatin'), (mm, 'Mao MMS'), 
                (sc, 'Sancar Cisplatin')]:

            if el1[1] == 'Mao MMS' or el2[1] == 'Mao MMS':
                cos_sim = round(1 - spatial.distance.cosine(
                    el1[0][1]['TOTAL_NORM'].tolist(), el2[0][1]['TOTAL_NORM'].tolist()), 2)
                
            else:
                
                cos_sim = round(1 - spatial.distance.cosine(
                    el1[0][0]['TOTAL_NORM'].tolist(), el2[0][0]['TOTAL_NORM'].tolist()), 2)

            df[el1[1]][el2[1]] = cos_sim

    return df


def fix_data(df, good_df):

    aa = list(set(df['CONTEXT']) ^ set(good_df['CONTEXT']))
    bb = pd.DataFrame([aa, [0]*len(aa), [0]*len(aa)]).T
    bb.set_axis(list(df.columns), axis=1, inplace=True)
    return pd.concat([df, bb], axis=0).reset_index(drop=True).sort_values('CONTEXT')


def fix_sancar(df):
    return df[ df['CONTEXT' ].str.contains('N') == False]


def fix_for_mao(df, mao):
    return pd.merge(mao['CONTEXT'], df, how='inner', on='CONTEXT')


def get_signature_df(mm, tm, tmm, tc, tcm, nc, ncm, nm, nmm, sc, scm):
    sig = pd.DataFrame(columns=['SemiSup Cisplatin', 'SemiSup MMS', 'Tombo MMS', 
        'Tombo Cisplatin', 'Sancar Cisplatin'])
    sig_mao = pd.DataFrame(columns=names)
    
    sig['SemiSup Cisplatin'] = nc['TOTAL_NORM']
    sig['SemiSup MMS'] = nm['TOTAL_NORM']
    sig['Tombo MMS'] = tm['TOTAL_NORM']
    sig['Tombo Cisplatin'] = tc['TOTAL_NORM']
    sig['Sancar Cisplatin'] = sc['TOTAL_NORM']

    sig_mao['SemiSup Cisplatin'] = ncm['TOTAL_NORM']
    sig_mao['SemiSup MMS'] = nmm['TOTAL_NORM']
    sig_mao['Tombo MMS'] = tmm['TOTAL_NORM']
    sig_mao['Tombo Cisplatin'] = tcm['TOTAL_NORM']
    sig_mao['Sancar Cisplatin'] = scm['TOTAL_NORM']
    sig_mao['Mao MMS'] = mm['TOTAL_NORM']

    return sig, sig_mao


def load_data(novoa_cis, novoa_mms, tombo_mms, tombo_cis, mao_mms, sancar_cis):
    mm = pd.read_csv(mao_mms, sep='\t').sort_values('CONTEXT')\
        .reset_index(drop=True)

    tm = pd.read_csv(tombo_mms, sep='\t').sort_values('CONTEXT')\
        .reset_index(drop=True)
    tmm = fix_for_mao(tm, mm).sort_values('CONTEXT').reset_index(drop=True)
    
    tc = pd.read_csv(tombo_cis, sep='\t').sort_values('CONTEXT')\
        .reset_index(drop=True)
    tcm = fix_for_mao(tc, mm).sort_values('CONTEXT').reset_index(drop=True)

    nc = fix_data(pd.read_csv(novoa_cis, sep='\t'), tm)\
        .sort_values('CONTEXT').reset_index(drop=True)
    ncm = fix_for_mao(nc, mm).sort_values('CONTEXT').reset_index(drop=True)

    nm = fix_data(pd.read_csv(novoa_mms, sep='\t'), tm)\
        .sort_values('CONTEXT').reset_index(drop=True)
    nmm = fix_for_mao(nm, mm).sort_values('CONTEXT').reset_index(drop=True)
    
    
    sc = fix_sancar(pd.read_csv(sancar_cis, sep='\t'))\
        .sort_values('CONTEXT').reset_index(drop=True)
    scm = fix_for_mao(sc, mm).sort_values('CONTEXT').reset_index(drop=True)

    sig, sig_mao = get_signature_df(
        mm, tm, tmm, tc, tcm, nc, ncm, nm, nmm, sc, scm
    )

    return (nc, ncm), (nm, nmm), (tm, tmm), (tc, tcm), (mm, mm), \
        (sc, scm), sig, sig_mao


def plot_heatmap(df, output):

    mask = np.zeros_like(df.values)
    mask[np.triu_indices_from(mask)] = True

    with sns.axes_style("white"):
        fig, ax = plt.subplots(figsize=(7, 5))
        ax = sns.heatmap(df, mask=mask, square=True, cmap='Blues', 
            cbar_kws={'label': 'Cosine Similarity'}, annot=True)
    
    outdir = os.path.join(output, "heatmap_all.pdf")
    fig.tight_layout()

    plt.savefig(outdir)
    plt.close()


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

    nc, nm, tm, tc, mm, sc, sig, sig_mao = load_data(
        novoa_cis, novoa_mms, tombo_mms, tombo_cis, mao_mms, sancar_cis
    )
    sig.to_csv(os.path.join(output, 'signatures.tsv'), sep='\t', index=None)    
    sig_mao.to_csv(os.path.join(output, 'signatures_mao.tsv'), sep='\t', index=None)

    df = get_cosines(nc, nm, tm, tc, mm, sc)
    df.to_csv(os.path.join(output, 'data_cosines.tsv'), sep='\t')

    plot_heatmap(df, output)

    

if __name__ == '__main__':
    main()