#!/usr/bin/env python3
import click
import pandas as pd
from tqdm import tqdm
from bgreference import refseq
from collections import Counter

import utils as ut
import nucleosome_damage_norm as ndn


def load_data(damage_nucleosomes, enrichment_path):
    """ Load datasets for further analysis 
    Args:
        damage_nucleosomes: path to damage within the nucleosomes
        enrichment_path: data of the sequence context probability 
    Returns:
        dam_nuc, enrichment: datasets
    """

    dam_nuc = pd.read_csv(damage_nucleosomes, sep='\t')
    enrichment = pd.read_csv(enrichment_path, sep='\t')

    return dam_nuc, enrichment


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


#Obtain base for the position in the file
def annot(df):
    seq_pos = refseq('saccer3', df['CHROM'], df['Start_nuc'] - 1, 149)
    return seq_pos


def get_info_damage(df):
    """ Get information to compute the expected damage 
    Args:
        df: dataframe with damage in nucleosomes
    Returns:
        sequences: sequence of every nucleosome
        strand: strand of the sequence
        N_damage: amount of the damage in that sequence
    """

    df = ndn.num2chr(df)
    
    df_id = pd.DataFrame(
        df.groupby(['ID', 'strand']).apply(lambda x: x['CHROM'].unique()[0]),
            columns=['CHROM']).reset_index()
    
    df_id['Start_nuc'] = df.groupby(
        ['ID', 'strand']).apply(lambda x: x['Start_nuc'].unique()[0]
    ).tolist()

    strand = df_id['strand'].tolist()
    N_damage = df.groupby(['ID', 'strand']).apply(lambda x: x.shape[0]).tolist()

    sequences = df_id.apply(annot, axis=1)

    return sequences, strand, N_damage


def get_expected_damage(df, enrichment):
    """ Get expected damage based on trinucleotide context 
    Args:
        df: dataframe with damage in nucleosomes
        enrichment: dataframe containing the triple enrichment probabilities
    Returns:
        expected: expected damage per position in the nucleosome
    """

    enrichment = enrichment.set_index('CONTEXT').to_dict()['TOTAL_NORM']
    prob_all_nuc = []

    sequences, strand, N_damage = get_info_damage(df)

    import pdb;pdb.set_trace()
    for i, seq in tqdm(enumerate(sequences.tolist())):

        if  strand[i] == '-':
            seq = ut.comp_seq(seq)
        
        prob_nuc = []

        for k in ut.slicing_window(seq, 3):
            try: 
                prob_nuc.append(enrichment[k])
            except:
                prob_nuc.append(0)
        prob_nuc_norm = [x / sum(prob_nuc) for x in prob_nuc]
        prob_all_nuc.append(prob_nuc_norm)
    import pdb;pdb.set_trace()
    

# ------------------------------------------------------------------------------
# CLICK
# ------------------------------------------------------------------------------

@click.command(short_help='script to radomized damage in nucleosomes')
@click.option(
    '-dn', '--damage_nucleosomes', default='', 
    help='nucleosome informationd bases'
)
@click.option(
    '-ed', '--enrichment_data', default='', help='data containing \
        enrichment probabilities normalized to one'
)
@click.option(
    '-o', '--output', default='', help='output folder'
)
def main(damage_nucleosomes, enrichment_data, output):

    dam_nuc, enrichment = load_data(damage_nucleosomes, enrichment_data)

    #Obtain observed damage in the nucleosomes
    final_dam = obtain_observed_damage(dam_nuc)

    #compute expected damage (Needs revision)
    expected_damage = get_expected_damage(dam_nuc, enrichment)

    import pdb;pdb.set_trace()

if __name__ == '__main__':
    main()