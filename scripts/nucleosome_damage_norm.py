#!/usr/bin/envs python3

import os
import click
import pybedtools
import numpy as np
import pandas as pd 
from operator import add
from bgreference import refseq
from collections import Counter

from tqdm import tqdm
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


def num2chr(df):
    df['CHROM'] = df.CHROM.apply(lambda x : ut.find_chr(x))
    return df


#Obtain base for the position in the file
def annot(df):
    seq_pos = refseq('saccer3', df['CHROM'].unique()[0], 
        df['Start_nuc'].unique()[0] - 1, 149)
    return seq_pos


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
            'Chr_nuc', 'Start_nuc', 'End_nuc', 'Val_1', 'Val_2', 'ID', 'Center_nuc', 
            'Overlapped'])
    else: 
        df = pd.read_csv(
            result.fn, sep='\t', names=['CHROM', 'Start', 'End', 'SEQ', 'strand', 
            'mer', 'min coverage', 'untreated_freq', 'treated_freq', 'value', 
            'PENTAMER', 'TRIPLET', 'Chr_nuc', 'Start_nuc', 'End_nuc',
            'Val_1', 'Val_2', 'ID', 'Center_nuc', 'Overlapped'])

    return df


def load_data(damage_path, nucleosome_info, enrichment_path):
    """ Load datasets for further analysis 
    Args:
        damage_path: path to damage dataset
        nucleosome_info: path to nucleosome information
    Returns:
        damage, nucleosomes: both datasets
    """

    damage = pd.read_csv(damage_path, sep='\t')
    if 'diff' in damage:
        damage = damage.drop(columns=['diff'])

    nucleosomes = pd.read_csv(
        nucleosome_info, sep='\t', #header=None,
    )
    nucleosomes.columns = ["CHROM", "Start", "End", "Val_1", "Val_2"]
    nucleosomes['ID'] = range(len(nucleosomes))
    
    #Clean-up damage data
    if 'index' in damage.columns: 
        damage = damage.drop('index', axis=1)
    damage = damage.rename(columns={"chrom": "CHROM", "base": "SEQ"})

    enrichment = pd.read_csv(enrichment_path, sep='\t')

    return damage, nucleosomes, enrichment


def obtain_observed_damage(df_nuc):
    """ Get the observed damage per position 
    Args:
        df_nuc: dataframe of damage in the nucleosomes
    Returns:
        observed: observed damage per position in the nucleosomes
    """
    
    counts_position = Counter(df_nuc['End'] - df_nuc['Start_nuc'])

    observed = pd.DataFrame(
        counts_position.items(), columns=['Position', 'Observed_damage']
    ).sort_values(by='Position', ascending=True)
    
    return observed


def get_expected_damage(df, enrichment):
    """ Get expected damage based on trinucleotide context 
    Args:
        df: dataframe with damage in nucleosomes
        enrichment: dataframe containing the triple enrichment probabilities
    Returns:
        expected: expected damage per position in the nucleosome
    """

    enrichment = enrichment.set_index('CONTEXT').to_dict()['TOTAL_NORM']

    df = num2chr(df)
    prob_all_nuc = []

    for i, j in tqdm(df.groupby(['ID', 'strand'])):
        N = j.shape[0]
        seq = annot(j)

        if j['strand'].unique()[0] == '-':
            seq = ut.comp_seq(seq)
        
        prob_nuc = []

        for k in ut.slicing_window(seq, 3):
            try: 
                prob_nuc.append(enrichment[k])
            except:
                prob_nuc.append(0)
        prob_nuc_norm = [x / sum(prob_nuc) * N for x in prob_nuc]
        prob_all_nuc.append(prob_nuc_norm)
    
    expected = [sum(x) for x in zip(*prob_all_nuc)]
    
    return expected


def arrange_df_and_save(final_dam, expected_damage, output):
    """ Get final df and save to tsv 
    Args:
        final_dam: dataframe containing the observed damage
        expected_damage: list of expected damage per nucleosome site
    Returns:
        final_dam: final damage df containing relative increase of damage
    """
    exp = pd.DataFrame(expected_damage).reset_index()
    exp['index'] = exp['index'] + 1
    exp = exp.rename(columns={'index': 'Position', 0: 'Expected_damage'})
    final_dam = final_dam.merge(exp, on='Position')

    final_dam['Relative Increase'] = (final_dam['Observed_damage'].values - \
        final_dam['Expected_damage'].values) / final_dam['Expected_damage'].values
    
    out_file = os.path.join(output, 'exp_obs_relinc.tsv')
    final_dam.to_csv(out_file, sep='\t', index=None)

    return final_dam

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
    '-ed', '--enrichment_data', default='', help='data containing \
        enrichment probabilities normalized to one'
)
@click.option(
    '-o', '--output', default='', help='output folder'
)
def main(damage, nucleosome_information, enrichment_data, output):
    #Obtain data
    damage, nucleosomes, enrichment = load_data(
        damage, nucleosome_information, enrichment_data
    )
    
    #obtain start and end positions
    damage = add_dam_se(damage)
    nucleosomes = add_nuc_se(nucleosomes)
    
    #chr in int format for pybedtools
    damage_bed = chr2num(damage)
    nucleosomes_bed = chr2num(nucleosomes)

    #interset bedtools
    df = intersect(damage_bed, nucleosomes_bed)

    #Select only damage in the nucleosome regions
    df_nuc = df[df['Overlapped'] != 0]
    print(df_nuc.shape)

    out_nuc_dam = os.path.join(output, 'damage_in_nucleosomes.tsv')
    df_nuc.to_csv(out_nuc_dam, sep='\t', index=None)

    #Obtain observed damage in the nucleosomes
    final_dam = obtain_observed_damage(df_nuc)
    
    #compute expected damage (Needs revision)
    expected_damage = get_expected_damage(df_nuc, enrichment)
    
    #save damage dataframe
    final_dam = arrange_df_and_save(final_dam, expected_damage, output)
    
    #plot relative damage increase in the nucleosome
    per_base_dir = os.path.join(output,'study_norm_all_bases')
    pl.plot_norm_nucleosomes(final_dam, per_base_dir)


if __name__ == '__main__':
    main()