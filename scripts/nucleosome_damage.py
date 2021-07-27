#!/usr/bin/envs python3

import os
import click
import pybedtools
import numpy as np
import pandas as pd 
from operator import add
from bgreference import refseq
from collections import Counter

import utils as ut
import plots as pl


# ------------------------------------------------------------------------------
# FUNCTIONS
# ------------------------------------------------------------------------------

#Add start and end positions of the nucleosomes to the df
def add_nuc_se(df):
    df['Nuc_center'] = df['End']
    df['Start'] = df['End'] - 73
    df['End'] = df['End'] + 74
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
    seq_pos = refseq('saccer3', df['CHROM'], df['Start'], 147)
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
        if len(list(seq[i])) != 147:
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


#Per base enrichment
def obtain_per_base_enrichment(df_nuc, norm, base):
    d = {}
    #Add the nucleosome postion of the damage
    #TODO <JB> Fix a possible problem with the strand 
    # neg = df_nuc[df_nuc['strand'] == '-']
    # counts_neg = Counter(neg['End_nuc'] - neg['Start'])
    # pos = df_nuc[df_nuc['strand'] == '+']
    # counts_pos = Counter(pos['End'] - pos['Start_nuc'])
    # counts_position = counts_pos + counts_neg
    counts_position = Counter(df_nuc['End'] - df_nuc['Start_nuc'])
    for k, v in counts_position.items():
        n = v / norm[k - 1][base[int(np.floor(len(base) / 2))]]
        d.update({k - 1 : n})
    
    df = pd.DataFrame(d.items(), columns=['POSITION', 'NORM'])
    df['NORM_2'] = df['NORM'] / df['NORM'].sum()
    df.sort_values(by=['POSITION'], inplace=True)

    return df.reset_index(drop=True)


#Triplet enrichment within the nucleosomes
def obtain_nucleosome_enrichment(df, triplet_norm, out_dir):
    triplet_exp = ut.get_context_counts(df, 'TRIPLET')
    
    results = ut.get_context_norm(triplet_norm, triplet_exp)

    triplet_dir = os.path.join(out_dir, 'triplet_norm_nuc')
    pl.obtain_plots(results, triplet_dir, 'triplet', 16)


#Employ EWA to smooth a function. Related to AdamOptimzer of AI. 
def exp_weighted_averages(df):

    for el in ['NORM_2', 'Random Model']:
        to_s = df[el].tolist()
        B = 0.8; V = [to_s[0] - to_s[0] * (1 - B)]
        for i in range(len(to_s)):
            V.append(B * V[i] + (1 - B) * to_s[i])

        name = '{}_smooth'.format(el)
        df[name] = V[1:]

    return df


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

    if len(base_study) == 1:
        damage = damage[damage['SEQ'] == base_study]
    else:
        damage['Motif'] = damage['mer'].apply(obtain_context, cent=base_study)
        damage = damage[damage['Motif'] == base_study]
        damage = damage.drop(['Motif'], axis=1)

    print('Initializing the analysis on nucleotide: {}'.format(base_study))

    return damage, nucleosomes


def obtain_context(df, cent):

    mid = int(np.floor(len(df) / 2))

    if len(cent) % 2 == 0:
        plus = int(len(cent) / 2)
        less = int(len(cent) / 2 - 1)
        return df[mid - less: mid + plus + 1]
    
    else:
        plus = int(np.floor(len(cent) / 2))
        less = int(np.floor(len(cent) / 2))
        return df[mid - less: mid + plus + 1]
        


def get_random_damage(norms, base):

    tot_bases = 0
    for el in norms: 
        tot_bases += el[base[int(np.floor(len(base) / 2))]]

    random_damage = [i[base[int(np.floor(len(base) / 2))]] / \
        tot_bases for i in norms]

    return random_damage
    


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

    #obtain random model of damage
    random_model = get_random_damage(norms, base_study)

    #chr in int format for pybedtools
    damage_bed = chr2num(damage)
    nucleosomes_bed = chr2num(nucleosomes)

    #interset bedtools
    df = intersect(damage_bed, nucleosomes_bed)

    #Select only damage in the nucleosome regions
    df_nuc = df[df['Overlapped'] != 0]
    print(df_nuc.shape)
    #Nucleosome enrichment as enrichment.py
    obtain_nucleosome_enrichment(df_nuc, triplet_norm, output)

    #obtain per-base enrichment within the nucleosome
    enrichment_df = obtain_per_base_enrichment(df_nuc, norms, base_study)

    #add random model to the enrichment. Be careful with dimensions here 
    enrichment_df['Random Model'] = random_model
       
    #Exponentially weighted averages for smoothing
    enrichment_df_smooth = exp_weighted_averages(enrichment_df)

    per_base_dir = os.path.join(
        output,'study_per_base_smooth_{}'.format(base_study)
    )
    
    pl.plot_per_base_enrichment(enrichment_df_smooth, per_base_dir, label='smooth')

    per_base_dir = os.path.join(
        output,'study_per_base_{}'.format(base_study)
    )
    pl.plot_per_base_enrichment(enrichment_df, per_base_dir, label='norm')


if __name__ == '__main__':
    main()