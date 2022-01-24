#!/usr/bin/env python3
import sys
import pandas as pd
from collections import Counter

from bgreference import refseq

import utils as ut


def polII_context1(df, sub, add, strand='+'):
    try:
        seq = refseq('saccer3', df['Chromosome'], 
            df['Start'] - sub, df['diff'] + add)
        if strand == '-':
            return ut.comp_seq(seq)
        else:
            return seq
    except:
        seq = '-'
    return seq


def polII_context2(genome, window):
    seq = Counter({});
    for el in genome:
        seq = seq + Counter(list(ut.slicing_window(el, window)))
    return seq


def get_context_polII(df):
    df['diff'] = df['End'] - df['Start']
    df = ut.num2chr(df)

    df['triplet'] = df.apply(polII_context1, args=(0, 2), axis=1)
    df['pentamer'] = df.apply(polII_context1, args=(1, 4), axis=1)

    penta = polII_context2(df['pentamer'].values, 5)
    triplet = polII_context2(df['triplet'].values, 3)

    return penta, triplet


def get_polII_damage(polii, damage, gen_triplet_prob, 
    output, gen_tri, gen_pent, label):
    polii['Start'] = polii['Start'].astype(int) + 1
    polii['End'] = polii['End'].astype(int) + 1

    trans = ut.chr2num(polii[polii['occupancy'] > 0])
    if label == 'both':
        non_trans = ut.chr2num(polii[polii['occupancy'] < 0])
    else: 
        non_trans = ut.chr2num(polii[polii['occupancy'] == 0])

    names = ['Chr', 'Start', 'End', 'base', 'strand', 'mer', 'min_coverage',
        'untreated_freq', 'treated_freq', 'value', 'PENTAMER', 'TRIPLET', 
        'Chromosome_t', 'Start_t', 'End_t', 'occupancy', 'Overlapped']

    damage_trans = ut.intersect(trans, damage, names)
    damage_non_trans = ut.intersect(non_trans, damage, names)
    
    trans_pent, trans_tri = get_context_polII(trans)
    non_trans_pent, non_trans_tri = get_context_polII(non_trans)

    ut.get_expected(damage_trans, gen_triplet_prob, trans_tri, 
        output, label='obs_exp_polii_trans')
    _, _ = ut.pre_enrichment_step(
        damage_trans, trans_tri, trans_pent, output, 
        label='polii_trans'
    )

    ut.get_expected(damage_non_trans, gen_triplet_prob, non_trans_tri, 
        output, label='obs_exp_polii_non_trans')
    _, _ = ut.pre_enrichment_step(
        damage_non_trans, non_trans_tri, non_trans_pent, output, 
        label='polii_non_trans'
    )


def select_polii_analysis(polii_occupancy, polii_type, damage, gen_triplet_prob,
    output, tri_counts, pent_counts):
    polii = pd.read_csv(polii_occupancy, header=None, sep='\t', 
        names=['Chromosome', 'Start', 'End', 'occupancy']).dropna()
    if polii_type == 'plus':
        damage = damage[damage['strand'] == '+']
    elif polii_type == 'minus':
        damage = damage[damage['strand'] == '-']
    get_polII_damage(
        polii, damage, gen_triplet_prob, output, tri_counts, pent_counts, 
        polii_type
    )


if __name__ == '__main__':
    pass