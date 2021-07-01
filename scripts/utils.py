#!/usr/bin/env python3

import os


base_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

base = ['A', 'C', 'G', 'T']


def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]


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
