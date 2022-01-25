#!/usr/bin/env python3
import sys
import pandas as pd

sys.path.append('../')
import utils as ut


def get_df_non_telomeres(df):
    df['Start'] = 30002; df['End'] = df['LEN'] - 30001
    df = df.drop(columns=['LEN'])
    return ut.chr2num(df)


def get_non_telemores_damaged(chrom_info, damage):
    
    non_telomeres = get_df_non_telomeres(chrom_info)
    names_tel = ['Chromosome_Or', 
            'Start_Or', 'End_Or', 'Chr', 'Start', 'End', 'strand', 
            'stat', 'damaged_base', 'PENTAMER', 'TRIPLET', 'Overlapped']

    return ut.intersect(damage, non_telomeres, names_tel), non_telomeres


def get_df_telomeres(df):
    start_tel = df.copy(); end_tel = df.copy()

    start_tel['Start'] = 1;  start_tel['End'] = 30001
    start_tel = start_tel.drop(columns='LEN')
    end_tel['Start'] = end_tel['LEN'] - 30000
    end_tel = end_tel.rename(columns={'LEN': 'End'}) 
    end_tel = end_tel[['Chromosome', 'Start', 'End']]
    
    return ut.chr2num(pd.concat([start_tel, end_tel]))


def get_telemores_damaged(chrom_info, damage):
    
    telomeres = get_df_telomeres(chrom_info)
    names_tel = ['Chromosome_Or', 
            'Start_Or', 'End_Or', 'Chr', 'Start', 'End', 'strand', 
            'stat', 'damaged_base', 'PENTAMER', 'TRIPLET', 'Overlapped']

    return ut.intersect(damage, telomeres, names_tel), telomeres


def do_telomere_analysis(chrom_info, damage, gen_triplet_prob, 
    gen_penta_prob, output):
    tel_intersect, telomeres = get_telemores_damaged(chrom_info, damage)
    
    tel_chr = ut.num2chr(telomeres.copy())
    tel_pent, tel_tri, tel_prop_bases = ut.counts_segment(tel_chr)
    ut.get_expected(
        tel_intersect, gen_triplet_prob, tel_tri, output, label='obs_exp_telomeres'
    )
    triplet_exp, penta_exp = ut.pre_enrichment_step(
        tel_intersect, tel_tri, tel_pent, output, label='telomere'
    )

    triplet_telomeres = ut.get_context_norm(tel_tri, triplet_exp)
    cosine_tel_tri = ut.calc_cosine_sim(
        gen_triplet_prob, triplet_telomeres, 'Triplet', 'Subtelomeric'
    )

    penta_telomeres = ut.get_context_norm(tel_pent, penta_exp)
    cosine_tel_pent = ut.calc_cosine_sim(
        gen_penta_prob, penta_telomeres, 'Pentamer', 'Subtelomeric'
    )

    damage_random_subset = damage.sample(n=tel_intersect.shape[0])
    _, _ = ut.pre_enrichment_step(
        damage_random_subset, tel_tri, tel_pent, output, label='telomere'
    )

    return telomeres, tel_intersect, cosine_tel_tri, cosine_tel_pent


def do_non_telomere_analysis(chrom_info, damage, gen_triplet_prob, 
    gen_penta_prob, output):
    non_tel_intersect, non_telomeres = get_non_telemores_damaged(
        chrom_info, damage
    )
    non_tel_chr = ut.num2chr(non_telomeres.copy())
    non_tel_pent, non_tel_tri, non_tel_prop_bases = ut.counts_segment(non_tel_chr)
    ut.get_expected(
        non_tel_intersect, gen_triplet_prob, non_tel_tri, 
        output, label='obs_exp_non_telomeres'
    )   
    triplet_exp, penta_exp = ut.pre_enrichment_step(
        non_tel_intersect, non_tel_tri, non_tel_pent, output, label='non_telomere'
    )

    triplet_non_telomeres = ut.get_context_norm(non_tel_tri, triplet_exp)
    cosine_non_tel_tri = ut.calc_cosine_sim(
        gen_triplet_prob, triplet_non_telomeres, 'Triplet', 'Non subtelomeric'
    )

    penta_non_telomeres = ut.get_context_norm(non_tel_pent, penta_exp)
    cosine_non_tel_pent = ut.calc_cosine_sim(
        gen_penta_prob, penta_non_telomeres, 'Pentamer', 'Non subtelomeric'
    )

    return cosine_non_tel_tri, cosine_non_tel_pent