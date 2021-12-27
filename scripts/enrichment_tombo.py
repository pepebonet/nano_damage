#!/usr/bin/env python3
# coding=utf8
import os
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


def decision_adjacent(pos_info, triplet, pentamer):
    trip_merge = pd.merge(
        pos_info, triplet, left_on='TRIPLET', right_on='CONTEXT'
        )
    pent_merge = pd.merge(
        pos_info, pentamer, left_on='PENTAMER', right_on='CONTEXT'
        )
    pent_merge['new_p-value'] = pent_merge['p-value'] / pent_merge['TOTAL_NORM'] / \
        trip_merge['TOTAL_NORM']

    return pos_info.loc[pent_merge['new_p-value'].idxmin() : \
        pent_merge['new_p-value'].idxmin()]


def filter_adjacent(df, triplet, pentamer):
    df_out = pd.DataFrame()
    for i, j in df.groupby('CHROM'):
        ordered = j.sort_values(by='pos')
        positions = ordered['pos'].tolist()

        rang = ranges(positions)
        
        ordered = ordered.reset_index(drop=True)
        ordered['pos_str'] = ordered['pos'].apply(str)

        for el in tqdm(rang):
            if el[0] == el[1]:
                p = ordered[ordered['pos_str']\
                    .str.contains(str(el[0]))].index[0]
                selected_pos = ordered.loc[p:p]
            else:
                pos_info = get_pos_info(ordered, el)
                selected_pos = decision_adjacent(pos_info, triplet, pentamer)

            df_out = pd.concat([df_out, selected_pos])
    print(df_out.shape)
    return df_out


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


#normalize results by the control vs control levels (Background noise)
def normalize_by_control(df, control):
    #First delete general norm
    df = df.drop(['NORM_2'], axis=1)

    #Order equally to substract context enrichment results
    control = order_df(df, control)

    #substract and make 0 in case the control is higher
    df['control_norm'] = df['REF_NORM'] - control['REF_NORM']
    df.loc[(df['control_norm'] < 0), 'control_norm'] = 0
    df['TOTAL_NORM'] = df['control_norm'] / df['control_norm'].sum()

    return df.rename(
        columns={"REF_NORM": "No_control_norm", "control_norm": "REF_NORM"}
    )


#Perform multiple testing correction
def mtc(df, alpha=0.01, method='fdr_bh'):
    _, corr, _, _ = statsmodels.stats.multitest.multipletests(
        df['p-value'].values, alpha=alpha, method=method
    )
    df['corrected_p'] = corr
    df = df[df['corrected_p'] <= alpha]
    print(df.shape)
    return df


#Extract control file and normalize treatment by it
def extract_control_norm(args, df, out_dir, path_extension):
    control = pd.read_csv(args.norm_control_triplet, sep='\t')
    
    if control.empty:
        raise Exception('Control file is empty')
    
    outfile = os.path.join(out_dir, path_extension)

    return normalize_by_control(df, control), outfile


def analysis_adjacent(df):
    sig_bases_dist = []
    final_df = pd.DataFrame()
    for i, j in df.groupby('CHROM'):
        chrom_df = j.sort_values(by='pos')
        chrom_df['diff'] = chrom_df.pos.diff()

        num = 0
        for el in chrom_df['diff'].tolist():
            if el == 1.0:
                num += 1
            elif math.isnan(el):
                pass
            else: 
                sig_bases_dist.append(num)
                num = 0
            
        final_df = pd.concat([final_df, chrom_df])
    
    return final_df, dict(Counter(sig_bases_dist))


#Obtain plots
def generate_plots(df_t, df_p, td, pd, sig_bases, out_dir):
    pl.obtain_plots(df_t, td, 'triplet', 16)
    pl.obtain_plots(df_p, pd, 'pentamer', 256)
    pl.significant_bases_distance(sig_bases, out_dir)


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


def main(args):
    chrom_len = pd.read_csv(
        args.chrom_len, sep='\t', names=['CHROM', 'LEN']
    )

    penta_ref, triplet_ref = ut.counts_reference_genome(chrom_len)

    df = obtain_dfs(args)

    df['SEQ'] = df.apply(ut.annot, axis = 1)
    df = df[df['p-value'] < args.pvalue]
    print(df.shape)
    import pdb;pdb.set_trace()
    if args.multiple_testing:
        df = mtc(df)

    df['PENTAMER'] = df.apply(ut.obtain_contex, args=(5,2), axis=1)
    df['TRIPLET'] = df.apply(ut.obtain_contex, args=(3,1), axis=1)

    penta_exp = ut.get_context_counts(df, 'PENTAMER')
    triplet_exp = ut.get_context_counts(df, 'TRIPLET')

    triplet_context = ut.get_context_norm(triplet_ref, triplet_exp)
    penta_context = ut.get_context_norm(penta_ref, penta_exp)
    
    if args.filters:
        df = filter_adjacent(df, triplet_context, penta_context)
        penta_exp = ut.get_context_counts(df, 'PENTAMER')
        triplet_exp = ut.get_context_counts(df, 'TRIPLET')

        triplet_context = ut.get_context_norm(triplet_ref, triplet_exp)
        penta_context = ut.get_context_norm(penta_ref, penta_exp)

    df, sig_bases_dist = analysis_adjacent(df)

    out_dir = get_output_dir(args, df.shape[0])

    #Use control to normalize triplet and pentamer context
    if args.norm_control_triplet:
        triplet_context, triplet_dir = extract_control_norm(
            args, triplet_context, out_dir, 'triplet_norm_control'
        )
    elif args.norm_control_pentamer:
        penta_context, penta_dir = extract_control_norm(
            args, penta_context, out_dir, 'penta_norm_control'
        )
    else:
        triplet_dir = os.path.join(out_dir, 'triplet')
        penta_dir = os.path.join(out_dir, 'pentamer')
    import pdb;pdb.set_trace()
    #Obtain plots
    generate_plots(
        triplet_context, penta_context, triplet_dir, 
        penta_dir, sig_bases_dist, out_dir
        )

    #write output
    ut.write_outputs(
        df, triplet_context, penta_context, out_dir, 'Nanopore'
    )


if __name__ == '__main__':
    args = parse_args()
    main(args)