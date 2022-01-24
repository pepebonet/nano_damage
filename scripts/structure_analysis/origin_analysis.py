#!/usr/bin/env python3
import sys
import pandas as pd

sys.path.append('../')
import utils as ut


def get_origins_damaged(damage, replication):
    #chr in int format for pybedtools
    damage_bed = ut.chr2num(damage)
    replication_bed = ut.chr2num(replication.copy())
    names_replication = ['Chromosome_Or', 
        'Start_Or', 'End_Or', 'ACS', 'Strand_Or', 'Class', 'Timing', 'IGR', 
        'WT G1 NDR (bp )', 'rpd3 G1 NDR (bp )', 'Chr', 'Start', 'End', 'base',
        'strand', 'mer', 'min_coverage', 'untreated_freq', 'treated_freq',
        'value', 'PENTAMER', 'TRIPLET', 'Overlapped']
    
    df = ut.intersect(damage_bed, replication_bed, names_replication)
    return df, replication_bed


def get_expected_origins(intersect, rep, ndr, gen_triplet_prob, output):
    for i in intersect.Timing.unique():
        df = intersect[intersect['Timing'] == i]
        if i == '.':
            label = 'obs_exp_origins_NDR'; enrich = 'origins_NDR'
            pent, tri, or_prop_bases = ut.counts_segment(ndr)
        else:
            label = 'obs_exp_origins_{}'.format(i)
            enrich = 'origins_{}'.format(i)
            pent, tri, or_prop_bases = ut.counts_segment(rep[rep['Timing'] == i])
        
        ut.get_expected(df, gen_triplet_prob, tri, output, label)
        _, _ = ut.pre_enrichment_step(df, tri, pent, output, enrich)


def do_origin_analysis(damage, replication, only_or, ndr, 
    gen_triplet_prob, output):
    or_intersect, replications = get_origins_damaged(damage, replication)

    or_pent, or_tri, or_prop_bases = ut.counts_segment(replication)
    ut.get_expected(or_intersect, gen_triplet_prob, or_tri, 
        output, label='obs_exp_origins')
    triplet_exp, penta_exp = ut.pre_enrichment_step(
        or_intersect, or_tri, or_pent, output, label='origins'
    )

    triplet_origins = ut.get_context_norm(or_tri, triplet_exp)
    cosine_rep_or = ut.calc_cosine_sim(
        gen_triplet_prob, triplet_origins, output, 'Replication Origins'
    )
    
    get_expected_origins(or_intersect, only_or, ndr, gen_triplet_prob, output)

    return cosine_rep_or