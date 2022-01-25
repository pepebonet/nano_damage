#!/usr/bin/envs python3

import os 
import click 
import pybedtools
import numpy as np
import pandas as pd
from bgreference import refseq
import matplotlib.pyplot as plt
from collections import Counter


chr_dict = {'chrI': 1, 'chrII': 2, 'chrIII': 3, 'chrIV': 4, 'chrV': 5, 
    'chrVI': 6, 'chrVII': 7, 'chrVIII': 8, 'chrIX': 9, 'chrX': 10, 'chrXI': 11, 
    'chrXII': 12, 'chrXIII': 13, 'chrXIV': 14, 'chrXV': 15, 'chrXVI': 16, 
    'chrM': 17, 'Mito' : 17}

chr_dict2 = {'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5, 
    'VI': 6, 'VII': 7, 'VIII': 8, 'IX': 9, 'X': 10, 'XI': 11, 
    'XII': 12, 'XIII': 13, 'XIV': 14, 'XV': 15, 'XVI': 16, 
    'M': 17, 'Mito' : 17}

base_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

# import plots as pl

# ------------------------------------------------------------------------------
# FUNCTIONS
# ------------------------------------------------------------------------------


def chr2num(df):
    df['CHROM'] = df.CHROM.apply(lambda x : chr_dict[x])
    return df


def chr2num2(df):
    df['CHROM'] = df.CHROM.apply(lambda x : chr_dict2[x])
    return df


def num2chr(df):
    df['CHROM'] = df.CHROM.apply(lambda x : find_chr(x))
    return df


def find_chr(x):
    for k, v in chr_dict.items():
        if v == x: 
            return k


def comp_seq(seq):
    s = ''
    for el in seq:
        s = base_dict[el] + s
    return s


def get_damage(damage):
    df = pd.read_csv(damage, sep='\t')

    df = df[df['base'] == 'G']

    df['Start'] = df['pos'] - 1
    df['End'] = df['pos']
    df = df.rename(columns={'chrom': 'CHROM'})

    cols = df.columns.tolist()
    cols = cols[:1] + cols[-2:] + cols[1:-2]
    df = df[cols]

    df = chr2num(df)

    pos_strand = df[df['strand'] == '+']
    neg_strand = df[df['strand'] == '-']

    return pos_strand, neg_strand


def get_transcription(df, flag):

    df = df.rename(columns={"Chromosome/scaffold name": "CHROM"})
    pos_strand = df[df['Strand'] == 1]
    neg_strand = df[df['Strand'] == -1]

    if flag == 'tss':
        pos_strand['Start'] = pos_strand['Transcription start site (TSS)'] - 650
        pos_strand['End'] = pos_strand['Transcription start site (TSS)'] + 650
        neg_strand['Start'] = neg_strand['Transcription start site (TSS)'] - 650
        neg_strand['End'] = neg_strand['Transcription start site (TSS)'] + 650
    else:
        pos_strand['Start'] = pos_strand['Transcript end (bp)'] - 650
        pos_strand['End'] = pos_strand['Transcript end (bp)'] + 650
        neg_strand['Start'] = neg_strand['Transcript start (bp)'] - 650
        neg_strand['End'] = neg_strand['Transcript start (bp)'] + 650
    
    df = pd.concat([pos_strand, neg_strand])

    cols = df.columns.tolist()
    cols.remove('CHROM')
    
    cols = ['CHROM'] + cols[-2:] + cols[:-2]
    df = df[cols]

    df = df[df['Start'] >= 0]
    
    df = chr2num2(df)

    pos_strand = df[df['Strand'] == 1]
    neg_strand = df[df['Strand'] == -1]

    return pos_strand, neg_strand


def intersect(damage, transcription):

    a = pybedtools.BedTool.from_dataframe(damage)
    b = pybedtools.BedTool.from_dataframe(transcription)
    result = a.intersect(b, wao = True)
    cols = damage.columns.tolist() + ['chrom', 'start', 'end'] + \
        transcription.columns.tolist()[3:] + ['Overlap'] 

    df = pd.read_csv(result.fn, sep='\t', names=cols)

    df = df[df['Overlap'] == 1]

    return df
        
#TODO check with negative strand
def get_expected(df):
    df = num2chr(df)
    genome = df.apply(
        lambda x: refseq('saccer3', x[0], x[1], x[2] - x[1]), axis=1
    )
    
    all_seqs = []
    for seq in genome:
        if df['Strand'].unique()[0] == -1:
            all_seqs.append(comp_seq(seq))
        else:
            all_seqs.append(seq)
    
    aa = list(genome)
    my = np.empty([1300, len(genome)], dtype='O')
    for i in range(len(genome)):
        try:
            my[:,i] = list(aa[i])
        except:
            print(len(aa[i]))
    
    total_Gs = []
    for el in my:
        # total_Gs.append(Counter(el)['G'])
        if Counter(el)['G'] < 40:
            total_Gs.append(500)
        else:
            total_Gs.append(Counter(el)['G'])
            # import pdb;pdb.set_trace()
    
    Gs = pd.DataFrame([range(-649, 651), total_Gs]).T
    Gs.columns = ['Relative Position', 'Normalizing Counts']

    df = chr2num(df)
      
    return Gs
        


def get_df_positions(df, expected, flagT, flagS):

    if flagT == 'tss':
        df['Relative Position'] = df['pos'] - \
            df['Transcription start site (TSS)'].astype(int)
    else:
        if flagS == 'pos':
            df['Relative Position'] = df['pos'] - \
                df['Transcript end (bp)'].astype(int)
        else:
            df['Relative Position'] = df['pos'] - \
                df['Transcript start (bp)'].astype(int)
    
    
    rel_pos = df.groupby('Relative Position').apply(
        lambda x: x.shape[0]).reset_index()

    norm =  pd.merge(rel_pos, expected, on='Relative Position', how='outer').fillna(0)

    norm['Normalized Counts'] = norm[0] / norm['Normalizing Counts']
    norm = norm.sort_values(by='Relative Position')

    return df, norm


def do_plots(df, output, flagT, flagS):

    fig, ax = plt.subplots(figsize=(10, 5))

    custom_lines = []
    
    plt.plot(df['Relative Position'], df['Normalized Counts'], 
        color='#16416e', alpha=1, label='Negative Mao')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    for el in [('Relative Increase', '#16416e')]:
        custom_lines.append(
            plt.plot([],[], marker="o", ms=8, ls="", mec='black', 
            mew=0, color=el[1], label=el[0])[0] 
        )

    ax.legend(
        bbox_to_anchor=(0., 1.0, 1., .102),
        handles=custom_lines, loc='upper right', 
        facecolor='white', ncol=1, fontsize=10, frameon=False
    )

    out_file = os.path.join(output, '{}_{}.pdf'.format(flagT, flagS))

    fig.tight_layout()
    plt.savefig(out_file)
    plt.close()



@click.command(short_help='Analyse trancription sites')
@click.option(
    '-dp', '--damaged-positions', required=True, 
    help='damaged positions output by DeepMP or tombo'
)
@click.option(
    '-ts', '--transcription-sites', required=True, 
    help='transcription sites from ensemble. Deafault: tss_tts_yeast.txt'
)
@click.option(
    '-o', '--output', default='', help='output folder'
)
def main(damaged_positions, transcription_sites, output):

    pos_damage, neg_damage = get_damage(damaged_positions)
    
    transcription = pd.read_csv(transcription_sites, sep='\t')
    transcription = transcription[~transcription['Gene name'].isnull()]
    
    pos_tss, neg_tss = get_transcription(transcription, 'tss')
    pos_tts, neg_tts = get_transcription(transcription, 'tts')

    exp_pos_tss = get_expected(pos_tss)
    exp_pos_tts = get_expected(pos_tts)
    exp_neg_tss = get_expected(neg_tss)
    exp_neg_tts = get_expected(neg_tts)

    in_pos_tss = intersect(pos_damage, pos_tss)
    in_pos_tts = intersect(pos_damage, pos_tts)

    in_neg_tss = intersect(neg_damage, neg_tss)
    in_neg_tts = intersect(neg_damage, neg_tts)

    df_pos_tss, rel_pos_tss = get_df_positions(in_pos_tss, exp_pos_tss, 'tss', 'pos')
    df_pos_tts, rel_pos_tts = get_df_positions(in_pos_tts, exp_pos_tts, 'tts', 'pos')
    df_neg_tss, rel_neg_tss = get_df_positions(in_neg_tss, exp_neg_tss, 'tss', 'neg')
    df_neg_tts, rel_neg_tts = get_df_positions(in_neg_tts, exp_neg_tts, 'tts', 'neg')

    do_plots(rel_pos_tss, output, 'tss', 'pos')
    do_plots(rel_pos_tts, output, 'tts', 'pos')
    do_plots(rel_neg_tss, output, 'tss', 'neg')
    do_plots(rel_neg_tts, output, 'tts', 'neg')
    
    



if __name__ == '__main__':
    main()