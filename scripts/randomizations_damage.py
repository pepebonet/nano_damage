#!/usr/bin/env python3
import click
import numpy as np
import pandas as pd
from tqdm import tqdm
from bgreference import refseq
from collections import Counter

import utils as ut
import plots as pl
import nucleosome_damage_norm as ndn
import miscellaneous.new_spectral as ns
import miscellaneous.plots_periodicity as pp


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



def annot(df):
    """ Obtain nucleosome sequence
    Args:
        df: dataframe with the position information
    Returns:
        seq_pos: nucleosome sequence for the position
    """
    seq_pos = refseq('saccer3', df['CHROM'], df['Start_nuc'] - 1, 149)
    return seq_pos


def get_info_damage(df):
    """ Get information to compute the expected damage 
    Args:
        df: dataframe with damage in nucleosomes
    Returns:
        sequences: sequence of every nucleosome
        strand: strand of the sequence
        N_damage: vector with the amount of damage
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
        mean_expected: mean expected damage per position in the nucleosome
        randoms_expected: randomizations of damage per position in the nucleosome
    """

    enrichment = enrichment.set_index('CONTEXT').to_dict()['TOTAL_NORM']
    prob_all_nuc = []; prob_all_nuc_mean = []

    sequences, strand, N_damage = get_info_damage(df)

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
        prob_nuc_norm_mean = [x / sum(prob_nuc) * N_damage[i] for x in prob_nuc]

        prob_all_nuc.append(prob_nuc_norm)
        prob_all_nuc_mean.append(prob_nuc_norm_mean)
    
    mean_expected = [sum(x) for x in zip(*prob_all_nuc_mean)]
    randoms_expected = do_randomizations(prob_all_nuc, N_damage)

    return mean_expected, randoms_expected


def do_randomizations(probs, N_damage):
    """ Get randomizations of expected damage
    Args:
        probs: normalized probabilities for every  
        N_damage: amount of the damage for every sequence
    Returns:
        expecteds: 1000 randomizations of expected damage
    """

    expecteds = []; 
    for i in tqdm(range(1000)):
        randoms = np.zeros([len(probs), 147])
        for j in range(len(probs)):
            random_draws = np.random.choice(
                range(147), N_damage[j], p=probs[j], replace=True
            )
            for el in random_draws: 
                randoms[j, el] += 1

        expecteds.append(sum(randoms))    

    return expecteds


def get_rel_increase(mean_expected, randoms_expected):
    """ Get relative increase of the randoms against the mean expected
    Args:
        mean_expected: normalized probabilities for every  
        randoms_expected: 1000 randomizations of damage
    Returns:
        rel_increases: relative increases for every randomization
    """

    rel_increases = []
    for el in randoms_expected:
        rel_increases.append((el - mean_expected) / mean_expected)
    
    return rel_increases


def compute_snr(rel_increase, peak_obs):
    """ Get relative increase of the randoms against the mean expected
    Args:
        rel_increase: relative increases for every randomization
        preak_obs: peak from the observed data 
    Returns:
        peaks: peak of periodicity
        snrs: signal to noise ratio calculated
    """

    peaks, snrs = [], []

    for signal in tqdm(rel_increase):
        x, y, snr, peak = ns.compute_spectrum(
            signal, norm=True, low_p=5, high_p=15, low_t=0, high_t=len(signal)-2,
            center=peak_obs
        )
        peaks.append(peak); snrs.append(snr)
    
    return peaks, snrs


def get_snr_observed(obs, exp):
    """ Get snr of the observed damage and the spectrum information
    Args:
        obs: dataframe with the observed damage
        exp: expected damage
    Returns:
        obs: dataframe containing the relative increase
        peak: peak of periodicity
        snr: signal to noise ratio calculated
        x: spectrum output x-axis
        y: spectrum output y-axis
    """

    exp = pd.DataFrame(exp).reset_index()
    exp['index'] = exp['index'] + 1
    exp = exp.rename(columns={'index': 'Position', 0: 'Expected_damage'})
    obs = obs.merge(exp, on='Position')

    rel_inc = (obs['Observed_damage'].values - \
        obs['Expected_damage'].values) / obs['Expected_damage'].values

    obs['rel_inc'] = rel_inc

    x, y, snr, peak = ns.compute_spectrum(
        rel_inc, norm=True, low_p=5, high_p=15, low_t=0, high_t=len(rel_inc)-2
    )
    
    return obs, peak, snr, x, y


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
    mean_expected, randoms_expected = get_expected_damage(dam_nuc, enrichment)

    rel_increases = get_rel_increase(mean_expected, randoms_expected)

    rel_inc_obs, peak_obs, snr_obs, x, y = get_snr_observed(final_dam, mean_expected)

    peaks, snrs = compute_snr(rel_increases, peak_obs)

    #TODO <JB> Maybe remove
    pl.plot_peaks(peaks, peak_obs, output)
    pl.plot_snrs(snrs, snr_obs, output)

    p_val = (len(np.argwhere(np.asarray(snrs) > snr_obs)) + 1) / (len(snrs) + 1)
    print(p_val)

    pp.plot(rel_inc_obs, x, y, snr_obs, peak_obs, p_val, output)

    #TODO <JB> 
    #   1.- Write information for functions properly
    #   2.- Try to expand width of plots to 147
    #   3.- Clean up and organize spectrum and plot scripts 
    #   4.- Extract damage at minor-in/out 
    

if __name__ == '__main__':
    main()