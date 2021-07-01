#!/usr/bin/env python3
import pandas as pd
import bottleneck as bn

import utils_structure as ut

def add_dam_se(df):
    df['Start'] = df['pos'] - 1
    df = df.rename(columns={"pos": "End"})
    cols = df.columns.tolist()
    cols = cols[:1] + cols[-1:] + cols[1:-1]
    df = df[cols]
    return df


def get_chrom_info(path):
    chrom_info = pd.read_csv(
        path, sep='\t', names=['Chromosome', 'LEN']
    )
    return chrom_info, bn.nansum(chrom_info['LEN'].values)


def do_genomewise_analysis(chrom_len, damage, output):
    chrom_info, genome = get_chrom_info(chrom_len)
    genome_pent, genome_tri, gen_prop_bases = ut.counts_segment(
        chrom_info, 'genome'
    )

    damage = add_dam_se(damage.copy())

    triplet_exp, penta_exp = ut.pre_enrichment_step(
        damage, genome_tri, genome_pent, output, label='genome'
    )
    gen_triplet_prob = ut.get_context_norm(genome_tri, triplet_exp)

    return gen_triplet_prob, chrom_info, damage, genome_tri, genome_pent