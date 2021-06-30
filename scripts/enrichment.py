#!/usr/bin/envs python3 

import os
import click
from itertools import islice
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt

from bgreference import refseq

names_all_1 = ['chrom', 'pos', 'base', 'strand', 'mer', \
    'min_coverage', 'untreated_freq', 'treated_freq']

names_all_2 = ['index', 'chrom', 'pos', 'base', 'strand', 'mer', 'min coverage',
       'untreated_freq', 'treated_freq', 'diff', 'group',
       'motif_1 DRWGGDD P-value', 'motif']

base_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

base = ['A', 'C', 'G', 'T']


def get_contexts(df, chrom_len): 
    penta_ref, triplet_ref = counts_reference_genome(chrom_len)
    
    df['PENTAMER'] = df.apply(obtain_context, args=(5,2), axis=1)
    df['TRIPLET'] = df.apply(obtain_context, args=(3,1), axis=1)
    
    penta_exp = get_context_counts(df, 'PENTAMER')
    triplet_exp = get_context_counts(df, 'TRIPLET')
    
    triplet_context = get_context_norm(triplet_ref, triplet_exp)
    penta_context = get_context_norm(penta_ref, penta_exp)  

    return df, triplet_context, penta_context 


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


def obtain_context(df, cont, cent):
    if cont == 5:
        return df['mer'][5:10]
    else: 
        return df['mer'][6:9]


def get_context_counts(df, context):
    return Counter(list(df[context]))


#normalize counts by the reference genome
def get_context_norm(ref, exp):
    new_d = {}
    for el in exp:
        if el != '-':
            norm = exp.get(el) / ref.get(el)
            new_d.update({el : norm})
    
    df = pd.DataFrame(new_d.items(), columns=['CONTEXT', 'REF_NORM'])
    df['TOTAL_NORM'] = df['REF_NORM'] / df['REF_NORM'].sum()

    return df


def generate_plots(df_t, df_p, td, pd, out_dir):
    obtain_plots(df_t, td, 'triplet', 16)
    obtain_plots(df_p, pd, 'pentamer', 256)


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



def order_plots(data, context, cent=None):
    if context == 'triplet':
        order_plot = order_muts_triplet(cent)
    else:
        order_plot = order_muts_pentamers(cent)
 
    df = data.set_index('CONTEXT')
    df = df.loc[order_plot]
    y = df['TOTAL_NORM'].tolist()

    return order_plot, y


def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]


#Obtain context plots for triplets and pentamers
def obtain_plots(data, outdir, context, chunk):
    order_plot, y = order_plots(data, context)

    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(20, 5),
                            gridspec_kw={'height_ratios': [1, 30]})

    #1st axis. Upper color plot
    colors = []
    colors_base = ['#1ebff0', '#050708', '#e62725', '#cbcacb']
    bot = 0.5
    for ix, c in enumerate(chunks(y, chunk)):
        colors.extend([colors_base[ix] for s in c])
        axs[0].barh(1, chunk, left=bot, color=colors_base[ix])
        bot += chunk


    axs[0].set_xlim(0, 4 * chunk)
    axs[0].spines['top'].set_visible(False)
    axs[0].spines['bottom'].set_visible(False)
    axs[0].spines['left'].set_visible(False)
    axs[0].spines['right'].set_visible(False)

    axs[0].get_yaxis().set_visible(False)
    axs[0].get_xaxis().set_visible(False)

    #2nd axis. Bar plot
    plt.sca(axs[1])

    axs[1].bar(order_plot, y, color=colors)
    axs[1].set_ylabel('Relative Probability',fontsize=24)
    
    if context == 'triplet':
        axs[1].set_xticks(order_plot)
        axs[1].set_xlabel('Triplet Context', fontsize=24)
    else: 
        axs[1].tick_params(labelbottom=False)
        axs[1].set_xlabel('Pentamer Context', fontsize=24)

    axs[1].spines['top'].set_visible(False)
    axs[1].spines['right'].set_visible(False)

    plt.setp([axs[1].get_xticklines(), axs[1].get_yticklines()], color='grey')

    axs[1].xaxis.set_ticks_position('none')
    for axis in ['top', 'bottom', 'left', 'right']:
        axs[1].spines[axis].set_linewidth(0.2)

    axs[1].xaxis.set_tick_params(pad=0.5)
    axs[1].yaxis.set_tick_params(pad=0.5, width=0.5)

    plt.xticks(fontsize=20, rotation=90)
    plt.yticks(fontsize=16)
    plt.tick_params(axis='both', which='both', bottom=False, left = False)

    fig.tight_layout()

    plt.savefig(outdir)
    plt.close()


@click.command(short_help='Get enrichment')
@click.option('-d', '--data', required=True)
@click.option('-cl', '--chrom_len', required=True)
@click.option('-o', '--output', required=True)
def main(data, chrom_len, output):
    df = pd.read_csv(data, sep='\t', compression='gzip')

    try:
        df.columns = names_all_1
    except:
        df.columns = names_all_2

    #optional to take only those where the treated has higher frequency
    df['value'] = df['treated_freq'] - df['untreated_freq']
    df = df[df['value'] > 0]

    chrom_lens = pd.read_csv(
        chrom_len, sep='\t', names=['CHROM', 'LEN']
    )

    df, triplet_context, penta_context = get_contexts(df, chrom_lens)

    triplet_dir = os.path.join(output, 'triplet')
    penta_dir = os.path.join(output, 'pentamer')

    generate_plots(
        triplet_context, penta_context, triplet_dir, 
        penta_dir, output
    )
    



if __name__ == "__main__":
    main()