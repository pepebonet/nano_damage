#!/usr/bin/env python3
import pandas as pd

import utils as ut


def get_open_close_damage(open_pos, close_pos, damage):

    open_damage = pd.merge(
        damage, open_pos, how='inner', on=['Chromosome', 'End']
    ).drop(columns=['Start_y']).rename(columns={'Start_x':'Start'})

    close_damage = pd.merge(
        damage, close_pos, how='inner', on=['Chromosome', 'End']
    ).drop(columns=['Start_y']).rename(columns={'Start_x':'Start'})

    return open_damage, close_damage


def do_accessibility_analysis(access_data, damage, gen_triplet_prob, output):
    open_pos = access_data[access_data['counts'] > 0].reset_index(drop=True)
    close_pos = access_data[access_data['counts'] <= 0].reset_index(drop=True)
    
    close_tri, close_pent = ut.context_counts(close_pos)
    open_tri, open_pent = ut.context_counts(open_pos)
    
    open_pos = ut.chr2num(open_pos); close_pos = ut.chr2num(close_pos)
    open_damage, close_damage = get_open_close_damage(
        open_pos, close_pos, damage
    )
    ut.get_expected(
        open_damage, gen_triplet_prob, open_tri, output, label='obs_exp_open'
    )
    ut.get_expected(
        close_damage, gen_triplet_prob, close_tri, output, label='obs_exp_close'
    )

    triplet_exp_open, penta_exp = ut.pre_enrichment_step(
        open_damage, open_tri, open_pent, output, label='open'
    )

    triplet_open = ut.get_context_norm(open_tri, triplet_exp_open)
    cosine_open = ut.calc_cosine_sim(
        gen_triplet_prob, triplet_open, output, 'Open chromatin'
    )

    triplet_exp_close, penta_exp = ut.pre_enrichment_step(
        close_damage, close_tri, close_pent, output, label='close'
    )

    triplet_close = ut.get_context_norm(close_tri, triplet_exp_close)
    cosine_close = ut.calc_cosine_sim(
        gen_triplet_prob, triplet_close, output, 'Close chromatin'
    )

    return open_damage, close_damage, cosine_open, cosine_close


def do_telomere_accessibility_analysis(open_damage, close_damage, 
    telomeres, access_data, gen_triplet_prob, output):
    names_access = ['Chromosome', 'Start', 'End', 'strand', 'p-value', 
        'SEQ', 'PENTAMER', 'TRIPLET', 'counts', 'Chromosome_tel', 'Start_tel', 
        'End_tel', 'Overlapped']

    telomeres['Start'] = telomeres['Start'] - 1
    access_data = ut.chr2num(access_data)
    open_tel_damage = ut.intersect(telomeres, open_damage, names_access)
    close_tel_damage = ut.intersect(telomeres, close_damage, names_access)

    names_access_tel = ['Chromosome', 'Start', 'End', 'counts', 
        'Chromosome_2', 'Start_2', 'End_2',  'Overlapped']

    access_tel = ut.num2chr(
        ut.intersect(telomeres, access_data, names_access_tel)
    )
    access_tel['counts'] = access_tel['counts'].astype('int64')
    open_tel = access_tel[access_tel['counts'] > 0]
    close_tel = access_tel[access_tel['counts'] <= 0]

    close_tel_tri, close_tel_pent = ut.context_counts(close_tel)
    open_tel_tri, open_tel_pent = ut.context_counts(open_tel)

    ut.get_expected(
        open_tel_damage, gen_triplet_prob, open_tel_tri, 
        output, label='obs_exp_open_tel'
    )
    ut.get_expected(
        close_tel_damage, gen_triplet_prob, close_tel_tri, 
        output, label='obs_exp_close_tel'
    )

    _, _ = ut.pre_enrichment_step(
        open_damage, open_tel_tri, open_tel_pent, output, label='open_telomere'
    )
    _, _ = ut.pre_enrichment_step(
        close_damage, close_tel_tri, close_tel_pent, output, 
        label='close_telomere'
    )

    return open_tel_damage, close_tel_damage