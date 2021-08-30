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
            'Chr_nuc', 'Start_nuc', 'End_nuc', 'Val_1', 'Val_2', 'ID', 'Center_nuc', 
            'Overlapped'])
    else: 
        df = pd.read_csv(
            result.fn, sep='\t', names=['CHROM', 'Start', 'End', 'SEQ', 'strand', 
            'mer', 'min coverage', 'untreated_freq', 'treated_freq', 'value', 
            'PENTAMER', 'TRIPLET', 'Chr_nuc', 'Start_nuc', 'End_nuc',
            'Val_1', 'Val_2', 'ID', 'Center_nuc', 'Overlapped'])

    return df


#Per base enrichment
def obtain_observed_damage(df_nuc):

    import pdb;pdb.set_trace()
    counts_position = Counter(df_nuc['End'] - df_nuc['Start_nuc'])
    import pdb; pdb.set_trace()
    return df_nuc


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
        #TODO <JB> review this bit of the normalization it is not working properly
        for k in ut.slicing_window(seq, 3):
            try: 
                prob_nuc.append(enrichment[k] * N)
            except:
                prob_nuc.append(0)
        
        prob_all_nuc.append(prob_nuc)

    pos_expected_damage = sum([ sum(x) for x in zip(*prob_all_nuc) ])
    
    return pos_expected_damage



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
    # import pdb;pdb.set_trace()
    #obtain start and end positions
    damage = add_dam_se(damage)
    nucleosomes = add_nuc_se(nucleosomes)
    # import pdb;pdb.set_trace()
    #chr in int format for pybedtools
    damage_bed = chr2num(damage)
    nucleosomes_bed = chr2num(nucleosomes)

    #interset bedtools
    df = intersect(damage_bed, nucleosomes_bed)

    #Select only damage in the nucleosome regions
    df_nuc = df[df['Overlapped'] != 0]
    print(df_nuc.shape)

    #compute expected damage (Needs revision)
    expected_damage = get_expected_damage(df_nuc, enrichment)
    import pdb;pdb.set_trace()
    #Obtain observed damage in the nucleosomes
    observed_damage = obtain_observed_damage(df_nuc)

    # per_base_dir = os.path.join(output,'study_norm_all_bases')
    # pl.plot_norm_nucleosomes(enrichment_df, per_base_dir, label='norm')


if __name__ == '__main__':
    main()