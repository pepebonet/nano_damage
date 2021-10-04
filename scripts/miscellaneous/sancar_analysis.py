#!/usr/bin/envs/ python3

import click
import numpy as np
from Bio import SeqIO
from collections import Counter



#Parse fastq files to extract last base enrichment and per read length
def parse_fastq(df):
    a = ''; read_length_dist = []

    for record in SeqIO.parse(df, "fastq"):
        import pdb;pdb.set_trace()
        read_length_dist.append(np.log(len(record.seq)))
        a += record.seq[-1]
    return Counter(a), read_length_dist

# ------------------------------------------------------------------------------
# CLICK
# ------------------------------------------------------------------------------

@click.command(short_help='script to get the nucleosome damage')
@click.option(
    '-ds', '--damage_sancar', default='', help='damaged fastq from sancar'
)
@click.option(
    '-o', '--output', default='', help='output folder'
)
def main(damage_sancar, output):
    #Obtain data
    parse_fastq(damage_sancar)
    import pdb;pdb.set_trace()

    #TODO <JB>:
    #   1.- Remove PCR artifacts 
    #       1.1. Use samtools and bedtools to delete reads where 
    #               pairs are mapped to the same location
    #   2.- Determine damage position and nucleotide composition
    #       2.1. 2 nt upstream of the first read start using samtools and bedtools



if __name__ == '__main__':
    main()

