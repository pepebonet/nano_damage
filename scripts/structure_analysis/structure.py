#!/usr/bin/env python3
import os
import sys
import click
import numpy as np
import pandas as pd

sys.path.append('/home/jbonet/Desktop/nano_damage/scripts/')
import plots as pl

import polii_analysis as pa
import origin_analysis as oa
import telomere_analysis as ta
import genomewide_analysis as ga
import accessibility_analysis as aa


def get_data(damage, replication_origins, ndr, accessibility):
    replication_data = pd.read_csv(
        replication_origins, sep ='\t').drop(columns=['Origin']
    )
    damage_data = pd.read_csv(
        damage, sep='\t').rename(columns={"chrom": "Chromosome"})

    nucl_dep = pd.read_csv(ndr, sep='\t') 
    
    replication_tot = pd.concat([replication_data, nucl_dep], sort=False)
    access_data = get_accessibility_data(accessibility)

    return damage_data, replication_tot, access_data, replication_data, nucl_dep


def get_accessibility_data(accessibility):
    return pd.read_csv(accessibility, 
        names=['Chromosome', 'Start', 'End', 'counts'], sep ='\t', 
        compression='gzip')
    

def get_mappability(mappability):
    df = pd.DataFrame(columns=['Chromosome', 'pos', 'mappability'])

    for file in os.listdir(mappability):
        df_tmp = pd.DataFrame(columns=['Chromosome', 'pos', 'mappability'])
        map_chr = np.fromfile(
            os.path.join(mappability, file), dtype='uint8'
        )
        pos = np.arange(1, len(map_chr) + 1)
        chrom = np.asarray([file.split('.', 1)[0]] * len(map_chr))
        df_tmp['Chromosome'] = chrom; df_tmp['pos'] = pos
        df_tmp['mappability'] = map_chr
        df = pd.concat([df, df_tmp])

    return df


def get_mappable_damage(mappability, damage):
    map_data = get_mappability(mappability)
    damage_map = pd.merge(
        damage, map_data, how='inner', on=['Chromosome', 'pos']
    )
    
    return damage_map[damage_map['mappability'] > 0].drop(
        columns=['mappability'])



@click.command(short_help='Parser of fast5 files')
@click.option(
    '-sd', '--significant-damage', required=True,
    help='Table containing all significant damaged bases'
)
@click.option(
    '-ro', '--replication-origins', required=True, 
    help='Yeast chromatin accessibility file'
)
@click.option(
    '-nd', '--nucleosome-depleated', required=True, 
    help='nucleosome depleated regions '
)
@click.option(
    '-chl', '--chrom_len', default='',
    help='Path to saccharomyces cerevisiae chromosome lenghts. '
)
@click.option(
    '-a', '--accessibility', 
    help='Yeast chromatin accessibility file. DNase-seq_saccer3.bed.gz'
)
@click.option(
    '-m', '--mappability', default='',
    help='Yeast genome mappability folder. globalmap_k28tok103'
)
@click.option(
    '-ts', '--transcribed-sites', default='',
    help='Yeast genome transcribed sites per chromosome.'
)
@click.option(
    '-po', '--polii_occupancy', default='',
    help='RNA polii occupancy (transcribed analysis).'
)
@click.option(
    '-pt', '--polii_type', type=click.Choice(['minus', 'plus', 'both']),
    help='Select the polii anaylsis (minus or plus strand or both).'
)
@click.option(
    '-o', '--output', default='', help='output folder'
)
def main(significant_damage, replication_origins, nucleosome_depleated, 
    chrom_len, accessibility, mappability, transcribed_sites, 
    polii_occupancy, polii_type, output):
    damage, replication, access_data, rep, ndr = get_data(
        significant_damage, replication_origins, 
        nucleosome_depleated, accessibility
    )

    #mappability
    if mappability:
        damage = get_mappable_damage(mappability, damage) 
    

    #genome
    gen_triplet_prob, chrom_info, damage, tri_counts, pent_counts, \
        gen_penta_prob = ga.do_genomewise_analysis(
            chrom_len, damage, output)
    #Origins and timing 
    cosine_tri_or, cosine_pent_or = oa.do_origin_analysis(
        damage, replication, rep, ndr, gen_triplet_prob, gen_penta_prob,
        output
    )
    
    #telomeres
    telomeres, tel_intersect, cosine_tel_tri, cosine_tel_pent \
        = ta.do_telomere_analysis(
        chrom_info, damage, gen_triplet_prob, gen_penta_prob, output
    )
    #polII
    if polii_occupancy:
        pa.select_polii_analysis(
            polii_occupancy, polii_type, damage.copy(), gen_triplet_prob, output, 
            tri_counts, pent_counts
        )
    #non telomeres
    cosine_non_tel_tri, cosine_non_tel_pent = ta.do_non_telomere_analysis(
        chrom_info, damage, gen_triplet_prob, gen_penta_prob, output
    )

    # accessibility
    open_damage, close_damage, cosine_open_tri, cosine_close_tri, cosine_open_pent,  \
        cosine_close_pent = aa.do_accessibility_analysis(
        access_data, damage, gen_triplet_prob, gen_penta_prob, output
    )

    #Damage in open and close telomeres
    open_tel_damage, close_tel_damage = aa.do_telomere_accessibility_analysis(
        open_damage, close_damage, telomeres, access_data, gen_triplet_prob, output
    )

    cosine = pd.concat(
        [cosine_tri_or, cosine_pent_or, cosine_tel_tri, cosine_tel_pent, 
        cosine_non_tel_tri, cosine_non_tel_pent, cosine_open_tri, cosine_close_tri,
        cosine_open_pent, cosine_close_pent]
    )
    import pdb;pdb.set_trace()
    pl.plot_cosine(cosine, output)
    print('All done!')


if __name__ == '__main__':
    main()
