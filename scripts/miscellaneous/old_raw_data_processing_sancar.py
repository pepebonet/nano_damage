#!/usr/bin/envs/ python3

import click
import pysam
import numpy as np
from collections import defaultdict



def read_pair_generator(bam, region_string=None):
    """ Generate read pairs in a BAM file or within a region string.
        Reads are added to read_dict until a pair is found.
    Args:
        bam: bam file to process
        region_string: Region to process
    Returns:
        generator object containig pairs of reads
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]



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
    bam = pysam.AlignmentFile(damage_sancar, 'rb')

    reads_to_delete = []; counter = 0
    for read1, read2 in read_pair_generator(bam):
        # import pdb;pdb.set_trace()
        if (read1.pos == read2.pos) and \
            (read1.query_alignment_end == read2.query_alignment_end): 

            reads_to_delete.append(read1)
            reads_to_delete.append(read2)
            # import pdb;pdb.set_trace()
        
        elif read1.is_duplicate or read2.is_duplicate:
            import pdb;pdb.set_trace()

        elif read1.is_unmapped or read2.is_unmapped:
            import pdb;pdb.set_trace()

        else:
            counter += 1

            if  counter % 500 == 0: 
                print(counter)
    import pdb;pdb.set_trace()
    #TODO <JB>:
    #   1.- Remove PCR artifacts 
    #       1.1. Use samtools and bedtools to delete reads where 
    #               pairs are mapped to the same location
    #   2.- Determine damage position and nucleotide composition
    #       2.1. 2 nt upstream of the first read start using samtools and bedtools



if __name__ == '__main__':
    main()

