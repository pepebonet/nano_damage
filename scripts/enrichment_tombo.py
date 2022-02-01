#!/usr/bin/env python3
# coding=utf8
import os
import sys
import math
import tombo
import argparse
import numpy as np
import pandas as pd
import statsmodels.api
from scipy import stats
from itertools import islice
from bgreference import refseq
from collections import Counter

from itertools import groupby
from operator import itemgetter

from tqdm import tqdm

from tombo import tombo_helper, tombo_stats, resquiggle
from tombo.tombo_stats import ModelStats, LevelStats

sys.path.append('../')
import plots as pl
import utils as ut


# ------------------------------------------------------------------------------
# ARGPARSER
# ------------------------------------------------------------------------------

def parse_args():

    parser = argparse.ArgumentParser(
          description='*** Enrichment analysis of the context of modified ' \
            'bases detected through tombo detection methods. ***'
    )
    parser.add_argument(
        'input', help='Absolute or relative path to stats data. ' \
           'Statiscical data is a binary file obtained after running tombo ' \
           'detect modification methods.'
    )
    parser.add_argument(
        '-chl', '--chrom_len', type=str, default='',
        help='Path to saccharomyces cerevisiae chromosome lenghts. '
    )
    parser.add_argument(
        '-o', '--output', type=str, default='',
        help='Path to output directory. '
    )
    parser.add_argument(
        '-ls', '--level_stats', action='store_true', default=False,
        help='Level type stats file '
    )
    parser.add_argument(
        '-not', '--norm_control_triplet', type=str, default='',
        help='Absolute or relative path to control levels. Usage: Normalize \
        results'
    )
    parser.add_argument(
        '-nop', '--norm_control_pentamer', type=str, default='',
        help='Absolute or relative path to control levels. Usage: Normalize \
        results'
    )
    parser.add_argument(
        '-ms', '--most_significant', type=int, default=100000,
        help='Number of most siginificant bases to be obtained \
        from the stats file'
    )
    parser.add_argument(
        '-bn', '--bases_number', type=int, default=1,
        help='Number of bases to analyze for significance. Default = 1'
    )
    parser.add_argument(
        '-bf', '--base_filter', nargs="+", default=[],
        help='Whether to filter for a single base'
    )
    parser.add_argument(
        '-mt', '--multiple_testing', action='store_true', default=False,
        help='Perform multiple testing correction. Default = False'
    )
    parser.add_argument(
        '-f', '--filters', action='store_true', default=False,
        help='Filter chr XII and M, and adjust adjacent sites.'
    )
    parser.add_argument(
        '-pv', '--pvalue', type=float, default=0.05,
        help='p-value threshold based on control comparison. default = 0.05.'
    )

    args = parser.parse_args()

    return args


# ------------------------------------------------------------------------------
# FUNCTIONS
# ------------------------------------------------------------------------------

#Tranform tombo interval data into dataframe
def intervalData2df(ido):
    c = []; p = []; snd = []; pv = []
    for el in ido: 
        c.append('chr' + el.chrm); p.append(el.end); snd.append(el.strand)
        pv.append(10 ** - float(el.reg_text.split(': ')[-1]))

    d = {'CHROM': c, 'pos': p, 'strand': snd, 'p-value':pv}
    df = pd.DataFrame(d)
    df['CHROM'] = df['CHROM'].replace('chrMito', 'chrM')

    return df


def filter_data(data):
    new_data = []
    for i in data:
        if i.chrm != 'XII' and i.chrm != 'Mito':
            new_data.append(i)
    return new_data


def ranges(nums):
    nums = sorted(set(nums))
    gaps = [[s, e] for s, e in zip(nums, nums[1:]) if s+1 < e]
    edges = iter(nums[:1] + sum(gaps, []) + nums[-1:])
    return list(zip(edges, edges))


def get_pos_info(ordered, el):
    pos1 = ordered[ordered['pos_str']\
                    .str.contains(str(el[0]))].index[0]
    pos2 = ordered[ordered['pos_str']\
        .str.contains(str(el[1]))].index[0]
    return ordered.loc[pos1:pos2]


#Parse tombo stats file
def obtain_dfs(args):
    #Init tombo class depending on the type of detection method employed:
    #Level sample compare
    if args.level_stats:
        ts = LevelStats(args.input)

    #Model sample compare
    else: 
        ts = ModelStats(args.input)
    #Extract most significant regions
    most_sig = ts.get_most_signif_regions(
        args.bases_number, args.most_significant
    )
    #Filter chromsome XII and Mito
    # if args.filters:
    most_sig_filt = filter_data(most_sig)

    return intervalData2df(most_sig_filt)


#Order control equally as the input df 
def order_df(df, control):
    order = df['CONTEXT'].tolist()
    control = control.set_index('CONTEXT')
    control = control.loc[order]
    return control.reset_index()


#Perform multiple testing correction
def mtc(df, alpha=0.01, method='fdr_bh'):
    _, corr, _, _ = statsmodels.stats.multitest.multipletests(
        df['p-value'].values, alpha=alpha, method=method
    )
    df['corrected_p'] = corr
    df = df[df['corrected_p'] <= alpha]
    print(df.shape)
    return df


#Obtain plots
def generate_plots(df_t, df_p, td, pd, out_dir):
    pl.obtain_plots(df_t, td, 'triplet', 16)
    pl.obtain_plots(df_p, pd, 'pentamer', 256)
    # pl.significant_bases_distance(sig_bases, out_dir)


def get_output_dir(args, damage):
    if args.output:
        out_dir = os.path.join(
            args.output, 'pval_{}_damage_{}'.format(args.pvalue, damage)
        )
    else:
        out_dir = os.path.join(
            os.path.dirname(args.input), 'pval_{}_damage_{}'.format(
                args.pvalue, damage
            )
        )
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    return out_dir


def df_for_pipeline(df):

    df = df.rename(columns={'CHROM': 'chrom', 'SEQ':'base', 'p-value': 'value'})

    if 'corrected_p' in df.columns:
        df.drop(['corrected_p'], inplace=True, axis=1)

    df['mer'] = df['PENTAMER']
    df['min_coverage'] = 0
    df['untreated_freq'] = 0
    df['treated_freq'] = 0

    df = df[['chrom', 'pos', 'base', 'strand', 'mer', 'min_coverage', \
        'untreated_freq', 'treated_freq', 'value', 'PENTAMER', 'TRIPLET']]

    return df


def main(args):
    chrom_len = pd.read_csv(
        args.chrom_len, sep='\t', names=['CHROM', 'LEN']
    )

    penta_ref, triplet_ref = ut.counts_reference_genome(chrom_len)

    df = obtain_dfs(args)
    df['SEQ'] = df.apply(ut.annot, axis = 1)

    if args.base_filter:
        new_df = pd.DataFrame()
        for el in args.base_filter:
            new_df = pd.concat([new_df, df[df['SEQ'] == el]])
        df = new_df

    df = df[df['p-value'] < args.pvalue]
    print(df.shape)

    if args.multiple_testing:
        df = mtc(df)

    df['PENTAMER'] = df.apply(ut.obtain_contex, args=(5,2), axis=1)
    df['TRIPLET'] = df.apply(ut.obtain_contex, args=(3,1), axis=1)

    penta_exp = ut.get_context_counts(df, 'PENTAMER')
    triplet_exp = ut.get_context_counts(df, 'TRIPLET')

    triplet_context = ut.get_context_norm(triplet_ref, triplet_exp)
    penta_context = ut.get_context_norm(penta_ref, penta_exp)

    out_dir = get_output_dir(args, df.shape[0])

    triplet_dir = os.path.join(out_dir, 'triplet.pdf')
    penta_dir = os.path.join(out_dir, 'pentamer.pdf')

    #Obtain plots
    generate_plots(
        triplet_context, penta_context, triplet_dir, 
        penta_dir, out_dir
        )

    df = df_for_pipeline(df)

    #write output
    ut.write_outputs(
        df, triplet_context, penta_context, out_dir, 'Nanopore'
    )


if __name__ == '__main__':
    args = parse_args()
    main(args)