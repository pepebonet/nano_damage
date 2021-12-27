#!/usr/bin/env python3

import os
import pybedtools
import pandas as pd
from itertools import islice
from bgreference import refseq
from collections import Counter
from scipy.stats import power_divergence
from sklearn.metrics.pairwise import cosine_similarity

import plots as pl


base = ['A', 'C', 'G', 'T']

base_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

chr_dict = {'chrI': 1, 'chrII': 2, 'chrIII': 3, 'chrIV': 4, 'chrV': 5, 
    'chrVI': 6, 'chrVII': 7, 'chrVIII': 8, 'chrIX': 9, 'chrX': 10, 'chrXI': 11, 
    'chrXII': 12, 'chrXIII': 13, 'chrXIV': 14, 'chrXV': 15, 'chrXVI': 16, 
    'chrM': 17, 'chrmt' : 17}


def chr2num(df):
    df['Chromosome'] = df.Chromosome.apply(lambda x : chr_dict[x])
    return df


def num2chr(df):
    df['Chromosome'] = df.Chromosome.apply(lambda x : find_chr(x))
    return df


def find_chr(x):
    for k, v in chr_dict.items():
        if v == x: 
            return k


#Obtain reverse complementary sequence
def comp_seq(seq):
    s = ''
    for el in seq:
        s = base_dict[el] + s
    return s


def counts_reference_genome(chrom_len):
    genome = chrom_len.apply(lambda x: refseq('saccer3', x[0], 1, x[1]), axis=1)
    
    seq = ''
    for el in genome:
        seq = seq + el

    penta = Counter(list(slicing_window(seq, 5)))
    triplet = Counter(list(slicing_window(seq, 3)))

    return penta, triplet


#Otain the windows of size n of the yeast genome
def slicing_window(seq, n):

    it = iter(seq)
    result = ''.join(islice(it, n))

    if len(result) == n:
        yield result

    for elem in it:
        result = result[1:] + elem
        yield result



#Obtain base for the position in the file
def annot(df):
    seq = refseq('saccer3', df['CHROM'], df['pos'], 1)
    if df['strand'] == '-':
        return base_dict[seq]
    else:
        return seq


#Obtain context for every base
def obtain_contex(df, cont, cent):
    try:
        seq = refseq('saccer3', df['CHROM'], df['pos'] - cent, cont)

        if len(seq) < cont:
            return '-'
        else:
            if df['strand'] == '-':
                return comp_seq(seq)
            else:
                return seq

    except:
        seq = '-'
    return seq


#Intersect damage and replication through bedtools
def intersect(damage, replication, names):
    a = pybedtools.BedTool.from_dataframe(damage)
    b = pybedtools.BedTool.from_dataframe(replication)
    result = b.intersect(a, wao=True)

    df = pd.read_csv(
        result.fn, sep='\t', names=names, low_memory=False)

    df = df[df['Overlapped'] > 0]
    
    return df


def pre_enrichment_step(damage, reg_tri, reg_pent, output, label):
    penta_exp = get_context_counts(damage, 'PENTAMER')
    triplet_exp = get_context_counts(damage, 'TRIPLET')
    get_enrichment(
        triplet_exp, penta_exp, reg_tri, reg_pent, output, label=label
    )
    return triplet_exp, penta_exp


def get_enrichment(triplet_exp, penta_exp, genome_tri, genome_pent,  
    output, label=''):
    triplet_context = get_context_norm(genome_tri, triplet_exp)
    penta_context = get_context_norm(genome_pent, penta_exp)

    triplet_dir = os.path.join(output, 'triplet_{}'.format(label))
    penta_dir = os.path.join(output, 'pentamer_{}'.format(label))

    generate_plots(
        triplet_context, penta_context, triplet_dir, 
        penta_dir
    )


def generate_plots(df_t, df_p, td, pd):
    pl.obtain_plots(df_t, td, 'triplet', 16)
    pl.obtain_plots(df_p, pd, 'pentamer', 256)


def calc_cosine_sim(genome, subset, output, label):
    
    df = pd.merge(genome, subset, on='CONTEXT', how='outer').fillna(0)
    score = cosine_similarity(df['TOTAL_NORM_x'].values.reshape(1, -1), 
        df['TOTAL_NORM_y'].values.reshape(1, -1))[0][0]
    
    df_score = pd.DataFrame([[score, label]], 
        columns=['Cosine-similarity', 'DNA structure'])
    
    return df_score


#normalize counts by the reference genome
def get_context_norm(ref, exp):
    new_d = {}
    for el in exp:
        if el != '-':
            try:
                norm = exp.get(el) / ref.get(el)
                new_d.update({el : norm})
            except:
                norm = 0
                new_d.update({el : norm})
                print('why?')
    
    df = pd.DataFrame(new_d.items(), columns=['CONTEXT', 'REF_NORM'])
    df['TOTAL_NORM'] = df['REF_NORM'] / df['REF_NORM'].sum()

    return df


def get_context_counts(df, context):
    return Counter(list(df[context]))


def get_expected(df, triplet_prob, reg_counts_triplet, 
    output, label=''):
    expected = triplet_prob.apply(
        lambda x: round(x[1] * reg_counts_triplet[x[0]], 2), axis=1
    )
    dict_pos = {'Observed' : df.shape[0], 'Expected' : 
         round(sum(expected))}

    stat, pval = power_divergence(
        [df.shape[0], round(sum(expected))], lambda_='log-likelihood'
    )

    pl.obs_exp_simple_plot(dict_pos, output, pval,label=label)


#get counts pentamers and triplets from a given genome region
def counts_segment(df, label=''):

    if label == 'genome':
        genome = df.apply(lambda x: refseq('saccer3', x[0], 1, x[1]), axis=1)
    else:
        if label == 'neg':
            genome = df.apply(
                lambda x: refseq('saccer3', x[0], x[1], x[2] - x[1]), axis=1
            )
            genome = genome.apply(lambda x: comp_seq(x)) 
        else:
            genome = df.apply(
                lambda x: refseq('saccer3', x[0], x[1], x[2] - x[1]), axis=1
            )

    penta = Counter({}); triplet = Counter({}); seq = ''
    for el in genome:
        seq = seq + el
        penta = penta + Counter(list(slicing_window(el, 5)))
        triplet = triplet + Counter(list(slicing_window(el, 3)))

    bases_dict = Counter(seq)
    total_bases = sum(bases_dict.values())

    for k, v in bases_dict.items():
        bases_dict[k] = round(v / total_bases, 4)
    
    return penta, triplet, bases_dict


def context_counts(df):
    df['Triplet'] = df.apply(counts_segment_open_close, args=(3,1), axis=1)
    df['Pentamer'] = df.apply(counts_segment_open_close, args=(5,2), axis=1)

    return Counter(df['Triplet']), Counter(df['Pentamer'])


def counts_segment_open_close(df, cont, cent):
    try:
        seq = refseq('saccer3', df['Chromosome'], df['End'] - cent, cont)
        if len(seq) < cont:
            return '-'
        else:
            return seq
    except:
        return '-'


#Yield successive n-sized chunks from l.
def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]


def order_plots(data, context, cent=None):
    if context == 'triplet':
        order_plot = order_muts_triplet(cent)
    else:
        order_plot = order_muts_pentamers(cent)

    df = data.set_index('CONTEXT')
    df = df.loc[order_plot]
    y = df['TOTAL_NORM'].tolist()

    return order_plot, y


#Obtain list of ordered triplets
def order_muts_triplet(cent=None):
    total = []
    if cent:
        for f in base:
            for last in base:
                a = '{}{}{}'.format(f, cent, last)
                total.append(a)
    else:
        for ref in base:
            for f in base:
                for last in base:
                    a = '{}{}{}'.format(f, ref, last)
                    total.append(a)

    return total


#Obtain list of ordered pentamers
def order_muts_pentamers(cent=None):
    total = []
    if cent:
        for f in base:
                for s in base:
                    for four in base:
                        for last in base:
                            a = '{}{}{}{}{}'.format(f, s, cent, four, last)
                            total.append(a)
    else:
        for ref in base:
                for f in base:
                    for s in base:
                        for four in base:
                            for last in base:
                                a = '{}{}{}{}{}'.format(f, s, ref, four, last)
                                total.append(a)

    return total


# ------------------------------------------------------------------------------
# OUTPUT FUNCTIONS
# ------------------------------------------------------------------------------


#Write outputs from enrichment analysis for Mao and Nanopore
def write_outputs(df, triplet, pentamer, out_dir, method, base=''):
    df.to_csv(os.path.join(
        out_dir, 'damage_most_sig_{}.tsv'.format(method)), 
        sep='\t', index=False
    )
    triplet.to_csv(os.path.join(
        out_dir, 'triplet_normalized_{}.tsv'.format(method)), 
        sep='\t', index=False
    )
    pentamer.to_csv(os.path.join(
        out_dir, 'pentamer_normalized_{}.tsv'.format(method)), 
        sep='\t', index=False
    )
