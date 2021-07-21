#!/usr/bin/envs python3

import os
import click
import pybedtools
import pandas as pd 
from operator import add
from bgreference import refseq
from collections import Counter
from scipy.signal import savgol_filter

import utils as ut
import plots as pl


# ------------------------------------------------------------------------------
# FUNCTIONS
# ------------------------------------------------------------------------------

#Add start and end positions of the nucleosomes to the df
def add_nuc_se(df):
    df['Nuc_center'] = df['End']
    df['Start'] = df['End'] - 200
    df['End'] = df['End'] + 200
    df = df[df['Start'] > 0].reset_index(drop=True)
    return df


#Add start and end positions of the damage to the df
def add_dam_se(df):
    df['Start'] = df['pos'] - 1
    df = df.rename(columns={"pos": "End"})
    cols = df.columns.tolist()
    cols = cols[:1] + cols[-1:] + cols[1:-1]
    df = df[cols]
    return df


# Change chr format for pybedtools
def chr2num(df):
    df['CHROM'] = df.CHROM.apply(lambda x : ut.chr_dict[x])
    return df


#Obtain base for the position in the file
def annot(df):
    seq_pos = refseq('saccer3', df['CHROM'], df['Start'], 400)
    return seq_pos


#Obtain number of nucleotides within each nucleosome position for normalization
def obtain_nuc_norm(df):

    seq = df.apply(annot, axis = 1)
    #Init sequence to counts
    seq2count = list(seq[0])

    triplet_nuc = Counter(ut.chunks(seq[0], 3))

    for i in range(1, len(seq)):
        #Ignore any nucleosome at the end of the chromosome 
        #Number of bases lower than 147 
        if len(list(seq[i])) != 400:
            continue
        # add sequence of every nucleosome for each on of the 147 bases
        seq2count = list(map(add, seq2count, seq[i]))
        #Nucleosome triple context
        triplet_nuc += Counter(ut.chunks(seq[i], 3))

    norm = []
    #Count nucleotides for each base in a nucleosome
    for el in seq2count:
        norm.append(Counter(list(el)))

    return norm, triplet_nuc


#Intersect damage and nucleosomes through bedtools
def intersect(damage, nucleosomes):
    a = pybedtools.BedTool.from_dataframe(damage)
    b = pybedtools.BedTool.from_dataframe(nucleosomes)
    result = a.intersect(b, wao = True)

    if damage.shape[1] == 15:
        df = pd.read_csv(
            result.fn, sep='\t', names=['CHROM', 'Start', 'End', 'SEQ', 'strand', 
            'mer', 'min coverage', 'untreated_freq', 'treated_freq', 'group', 
            'motif_1 DRWGGDD P-value', 'motif', 'value', 'PENTAMER', 'TRIPLET', 
            'Chr_nuc', 'Start_nuc', 'End_nuc', 'Val_1', 'Val_2', 'Center_nuc', 
            'Overlapped'])
    else: 
        df = pd.read_csv(
            result.fn, sep='\t', names=['CHROM', 'Start', 'End', 'SEQ', 'strand', 
            'mer', 'min coverage', 'untreated_freq', 'treated_freq', 'value', 
            'PENTAMER', 'TRIPLET', 'Chr_nuc', 'Start_nuc', 'End_nuc',
            'Val_1', 'Val_2', 'Center_nuc', 'Overlapped'])

    return df


# #Per base enrichment
# def obtain_per_base_enrichment(df_nuc, norm, base):
#     counts_position = Counter(df_nuc['End_nuc'] - df_nuc['End'])
    
#     df = pd.DataFrame(counts_position.items(), columns=['Position', 'Damage'])
#     df.sort_values(by=['Position'], inplace=True)

#     return df.reset_index(drop=True)

#Per base enrichment
def obtain_per_base_enrichment(df_nuc, norm, base):
    d = {}
    #Add the nucleosome postion of the damage
    counts_position = Counter(df_nuc['End_nuc'] - df_nuc['End'])
    for k, v in counts_position.items():
        n = v / norm[k][base]
        d.update({k : n})
    
    df = pd.DataFrame(d.items(), columns=['Position', 'NORM'])
    df['NORM_2'] = df['NORM'] / df['NORM'].sum()
    df.sort_values(by=['Position'], inplace=True)

    return df.reset_index(drop=True)


def load_data(damage_path, nucleosome_info, base_study):
    damage = pd.read_csv(damage_path, sep='\t')
    if 'diff' in damage:
        damage = damage.drop(columns=['diff'])

    nucleosomes = pd.read_csv(
        nucleosome_info, sep='\t', #header=None,
    )
    nucleosomes.columns = ["CHROM", "Start", "End", "Val_1", "Val_2"]
    
    #Filter data for only the target nucleotide
    if 'index' in damage.columns: 
        damage = damage.drop('index', axis=1)
    damage = damage.rename(columns={"chrom": "CHROM", "base": "SEQ"})

    return damage, nucleosomes


# ------------------------------------------------------------------------------
# CLICK
# ------------------------------------------------------------------------------

@click.command(short_help='script to get the nucleosome damage')
@click.option(
    '-da', '--damage', default='', help='damaged bases'
)
@click.option(
    '-ni', '--nucleosome_information', default='', 
    help='nucleosome informationd bases'
)
@click.option(
    '-bs', '--base_study', default='G', help='damaged bases'
)
@click.option(
    '-o', '--output', default='', help='output folder'
)
def main(damage, nucleosome_information, base_study, output):
    #Obtain data
    damage, nucleosomes = load_data(damage, nucleosome_information, base_study)

    #obtain start and end positions
    damage = add_dam_se(damage)
    nucleosomes = add_nuc_se(nucleosomes)
    
    #obtain norm nucleotide counts
    norms, triplet_norm = obtain_nuc_norm(nucleosomes)

    #chr in int format for pybedtools
    damage_bed = chr2num(damage)
    nucleosomes_bed = chr2num(nucleosomes)

    #interset bedtools
    df = intersect(damage_bed, nucleosomes_bed)

    #Select only damage in the nucleosome regions
    df_nuc = df[df['Overlapped'] != 0]

    #obtain per-base enrichment within the nucleosome
    enrichment_df = obtain_per_base_enrichment(df_nuc, norms, base_study)
    
    enrichment_df['smooth_damage']= savgol_filter(
        enrichment_df['NORM_2'].values, 31, 3
    )

    pl.plot_damage_nuc_linker(enrichment_df, output)
    


if __name__ == '__main__':
    main()