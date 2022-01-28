import os
import sys
import subprocess
from os.path import dirname, exists
from bgreference import hg19
import pysam
import pandas as pd

MAIN_DIR = '/workspace/projects/repair_states/data/UV_damage/sancar/reads/'
BOWTIE = 'bowtie2'
SAMTOOLS = 'samtools'
JAVA = 'java'
REF = '/workspace/datasets/genomes/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome'
C = 12
GB = 200

REMOVE_DUPLICATES_PAIRS = True
DISCARD_MISLIGATED_ADAPTER = True

#23
ANALYSIS_DIR = MAIN_DIR.rstrip('/')
FILENAMES = ['SRR5461431', 'SRR5461432']
#FILENAMES = ['SRR5461431']
FS = ['{}_1.fastq.gz'.format(filename) for filename in FILENAMES] + ['{}_2.fastq.gz'.format(filename) for filename in FILENAMES]
FILE_PAIR_FORMATS = ['{}_{}.fastq.gz'.format(filename, '{}') for filename in FILENAMES]
discard_adapter_5p = 'GACTGGTTCCAATTGAAAGTGCTCTTCCGATCT'
# TODO try reverse complement and on R2 as well?

if not REMOVE_DUPLICATES_PAIRS:
    FILE_PAIR_FORMATS = []


def _get_proc_result(proc_out):
    return proc_out.read().decode("utf-8").strip()


def get_bedtools_path():
    proc = subprocess.Popen("which closestBed", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    proc.wait()
    bedtools_path = dirname(_get_proc_result(proc.stdout))
    return bedtools_path


def check_install():
    bow = subprocess.Popen('which {}'.format(BOWTIE), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    if not bow.stdout.read():
        print('{} NOT INSTALLED'.format(BOWTIE))
        sys.exit(1)

    sam = subprocess.Popen('which {}'.format(SAMTOOLS), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    if not sam.stdout.read():
        print('{} NOT INSTALLED'.format(SAMTOOLS))
        sys.exit(1)

    jav = subprocess.Popen('which {}'.format(JAVA), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    if not jav.stdout.read():
        print('{} NOT INSTALLED'.format(JAVA))
        sys.exit(1)


def align_paired(file_name_1, file_name_2, stats_file):
    sam_path = '{}/mapping_hg19/{}_both.sam'.format(ANALYSIS_DIR, file_name_1.split('.fastq')[0])
    if not exists(sam_path.replace('.sam','.sort.bam')):
        command = '{} -x {} -1 {}/{} -2 {}/{} -S {} -p {} --phred33 --maxins 1000 --seed 123'.format(BOWTIE, REF, ANALYSIS_DIR, file_name_1, ANALYSIS_DIR, file_name_2, sam_path, C)
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        stats_file.write('---- Alignment\n')
        stats_file.write('command: ' + command + '\n\n')
        stats_file.write(get_proc_result(proc) + '\n\n')
        err = proc.stderr.read()
        if err != b'':
            print(err)
        proc.wait()

        explore_samtools(sam_path, stats_file)
    return sam_path


def get_proc_result(proc):
    return proc.stdout.read().decode("utf-8").strip()+proc.stderr.read().decode("utf-8").strip()


def explore_samtools(sam_path, stats_file):
    command_qual_15 = '{} view -c -q 15 {} -@ {}'.format(SAMTOOLS, sam_path, C)
    proc_15 = subprocess.Popen(command_qual_15, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    command_qual_20 = '{} view -c -q 20 {} -@ {}'.format(SAMTOOLS, sam_path, C)
    proc_20 = subprocess.Popen(command_qual_20, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    command_qual_30 = '{} view -c -q 30 {} -@ {}'.format(SAMTOOLS, sam_path, C)
    proc_30 = subprocess.Popen(command_qual_30, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

    proc_15.wait()
    proc_20.wait()
    proc_30.wait()

    stats_file.write('---- MAPQ\n')
    stats_file.write('command: {}\n'.format(command_qual_15))
    stats_file.write('command: {}\n'.format(command_qual_20))
    stats_file.write('command: {}\n'.format(command_qual_30))
    stats_file.write('{} reads with MAPQ >=15\n'.format(get_proc_result(proc_15))) # todo percentage?
    stats_file.write('{} reads with MAPQ >=20; 0.99 probability of a correct match\n'.format(get_proc_result(proc_20))) # todo percentage?
    stats_file.write('{} reads with MAPQ >=30; 0.999 prob.\n\n'.format(get_proc_result(proc_30)))

    command_forward = '{} view -c -F 20 {} -@ {}'.format(SAMTOOLS, sam_path, C)
    proc_fwd = subprocess.Popen(command_forward, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    command_reverse = '{} view -c -F 4 -f 16 {} -@ {}'.format(SAMTOOLS, sam_path, C)
    proc_rev = subprocess.Popen(command_reverse, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

    proc_fwd.wait()
    proc_rev.wait()

    stats_file.write('---- strands (mapped reads only)\n')
    stats_file.write('command forward: {}\n'.format(command_forward))
    stats_file.write('command reverse: {}\n'.format(command_reverse))
    stats_file.write('forward: {}\n'.format(get_proc_result(proc_fwd)))
    stats_file.write('reverse: {}\n\n'.format(get_proc_result(proc_rev)))


def _is_non_zero_file(fpath):
    """ https://stackoverflow.com/a/15924160 """
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0


def generate_bam(file_name, stats_file):
    """
    command: samtools view -F 4 -b ./2018-11-20_BCCW53ANXX/mapping_hg19/18_3757_2_30459_ATGTCA.sam > ./2018-11-20_BCCW53ANXX/mapping_hg19/18_3757_2_30459_ATGTCA.bam -@ 45
    sort command: samtools sort ./2018-11-20_BCCW53ANXX/mapping_hg19/18_3757_2_30459_ATGTCA.bam -@ 45 -o ./2018-11-20_BCCW53ANXX/mapping_hg19/18_3757_2_30459_ATGTCA.sort.bam
    index command: samtools index ./2018-11-20_BCCW53ANXX/mapping_hg19/18_3757_2_30459_ATGTCA.sort.bam -@ 45
    :param file_name:
    :return:
    """
    sam_path = file_name
    bam_path = sam_path.replace('.sam','.bam')
    sort_path = bam_path.replace('.bam','.sort.bam')

    if not exists(sort_path):
        stats_file.write('---- bam file (mapped reads only)\n')
        command_bam = '{} view -f 2 -@ {} -b {} > {} '.format(SAMTOOLS, C, sam_path, bam_path)
        bam_proc = subprocess.Popen(command_bam, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        stats_file.write('command: {}\n'.format(command_bam))
        bam_proc.wait()

        command_sort = '{} sort {} -@ {} -m {}G -o {}'.format(SAMTOOLS, bam_path, C, int(GB/C),sort_path)
        proc_sort = subprocess.Popen(command_sort, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        stats_file.write('command sort: {}\n'.format(command_sort))
        err = proc_sort.stderr.read()
        if err != b'':
            print(err)
        proc_sort.wait()

    if not exists(sort_path+'.bai'):
        command_index = '{} index {} -@ {}'.format(SAMTOOLS, sort_path, C)
        proc_ind = subprocess.Popen(command_index, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        stats_file.write('command index: {}\n'.format(command_index))
        err = proc_ind.stderr.read()
        if err != b'':
            print(err)
        proc_ind.wait()

    for fpath in [sam_path, bam_path]:
        if _is_non_zero_file(sort_path) and exists(fpath):
            os.remove(fpath)
    return sort_path


def discard_misligated_adapters_at_start_paired(stats_file, filename1, filename2, adapter_seq):
    filename1_discarded_misligated = filename1.replace('.fastq.gz','.discarded_misligated.fastq.gz')
    filename2_discarded_misligated = filename2.replace('.fastq.gz','.discarded_misligated.fastq.gz')
    f1_t_path = '{}/{}'.format(ANALYSIS_DIR, filename1_discarded_misligated)
    f2_t_path = '{}/{}'.format(ANALYSIS_DIR, filename2_discarded_misligated)
    if not os.path.exists(f1_t_path) and not os.path.exists(f2_t_path):
        command_cutadapt = 'cutadapt --times 5 --overlap 9 -j 0 -g {} -G {} -o {} -p {} --pair-filter=any --discard ' \
                           '{} {}'.format(adapter_seq, adapter_seq,
                                          '{}/{}'.format(ANALYSIS_DIR, filename1_discarded_misligated),
                                          '{}/{}'.format(ANALYSIS_DIR, filename2_discarded_misligated),
                                          '{}/{}'.format(ANALYSIS_DIR, filename1),
                                          '{}/{}'.format(ANALYSIS_DIR, filename2))
        proc = subprocess.Popen(command_cutadapt, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        stats_file.write('---- Discarding adapters\n')
        stats_file.write('command: ' + command_cutadapt + '\n\n')
        err = proc.stderr.read()
        stats_file.write(get_proc_result(proc) + '\n\n')
        if err != b'':
            print(err)
        proc.wait()
    return filename1_discarded_misligated, filename2_discarded_misligated


# HERE A SECTION OF FIGURING OUT THE POSITON OF THE DAMAGE


def reverse_complement(bases):
    """
    Reverse complements given string of bases
    :param bases:
    :return:
    """
    complementary_base = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    reverse_complemented = [complementary_base[base] for base in list(bases)][::-1]

    return ''.join(reverse_complemented)


def get_preceding(chrom, pos, length=4):
    try:
        context = hg19(chrom, pos, length)  # bgreference is 1 based; bed is 0 based
        if len(context) != length:
            context += 'N' * (length - len(context))
    except ValueError:  # something does not map to bgreference
        if pos == 0:  # 5' neighbour sticks out of the assembly
            context = 'N' + hg19(chrom, pos + 1, length - 1)
        elif pos < 0:
            context = 'N' * (-(pos - 1)) + hg19(chrom, 1, (length + (pos - 1)))
        else:
            context = ''
            print('sth wrong with getting the bases')
            print(chrom, pos)
            import sys
            sys.exit()
    return context


def extract_pos_from_sam(samfile):
    allowed_chrs = ['chr' + str(x) for x in range(1, 23)] + ['chrX', 'chrY', 'chrM']

    bed_lines = []
    num = 0
    for c, chrom in enumerate(samfile.references):
        if chrom not in allowed_chrs:
            # print(chrom)
            continue
        for read in samfile.fetch(reference=chrom):
            if not read.is_secondary:  # only primary alignments
                if read.is_read2:
                    continue
                if read.is_read1:
                    s, e = read.reference_start, read.reference_end  # 0 based start; reference_end points to one past the last aligned residue
                    if not read.is_reverse:  # strand + of the read1
                        if s != 0:  # exclude reads mapping to beginning of assembly & damage is outside of the reference
                            damage_pos = s -1 # 1 downstream of 3' end;

                            ref = get_preceding(chrom, damage_pos, 2)
                            bed_lines.append((chrom, damage_pos - 1, damage_pos+1, ref, '+'))

                    else:  # strand -
                        if e != samfile.lengths[c]:  # exclude reads mapping to end of assembly & damage is outside of the reference
                            damage_pos = e + 1  # "reference_end points to one past the last aligned residue"

                            ref = reverse_complement(get_preceding(chrom, damage_pos, 2))
                            bed_lines.append((chrom, damage_pos - 1, damage_pos+1, ref, '-'))

    samfile.close()
    return bed_lines


def _run_bedtool(command, save_path):
    print(command + '\n')
    proc = subprocess.Popen(command, stdout=open(save_path, 'w'), stderr=subprocess.PIPE, shell=True)
    proc.wait()


def save_damage_pos_bed(bam_sort_path, bedtools_path, genome_path):
    outfile = bam_sort_path.replace('mapping_hg19','analysis_files_hg19').replace('_both.sort.bam','.BED')

    if not exists(outfile):
        os.makedirs(dirname(outfile), exist_ok=True)

        samfile = pysam.AlignmentFile(bam_sort_path, 'rb')
        read_bed_lines = extract_pos_from_sam(samfile)
        print('bed lines before and after deduplication:')
        print(len(read_bed_lines))
        read_bed_lines = list(set(read_bed_lines))
        print(len(read_bed_lines))
        bed_df = pd.DataFrame(read_bed_lines, columns=['chr', 'start', 'end', 'dipyr', 'strand'])
        bed_df = bed_df[bed_df.dipyr.isin(['TT','CT'])].copy() # TODO this is for CPD! for 64-PP is TC.
        bed_df = bed_df.sort_values(by=['chr', 'start'], ascending=[True, True])
        print('and after getting only the wanted dipyrimidine types')
        print(len(bed_df))
        #print(bed_df.chr.value_counts())
        #print(bed_df.strand.value_counts())
        #print(bed_df.dipyr.value_counts())

        bed_df.to_csv(outfile + '_unsorted', sep='\t', index=False, header=False)

        from time import sleep
        sleep(2)

        command = '{}/sortBed -i {} -faidx {}'.format(bedtools_path, outfile + '_unsorted', genome_path)
        print(command)
        _run_bedtool(command, outfile)

    '''if _is_non_zero_file(outfile) and exists(outfile + '_unsorted'):
        os.remove(outfile + '_unsorted')'''

    return outfile




if __name__ == '__main__':
    check_install()
    genome_path = '/workspace/projects/repair_states/data/genomechunks/human.hg19.genome'
    bedtools_path = get_bedtools_path()


    for filename in FILENAMES:
        ANALYSIS_DIR = MAIN_DIR + filename
        print(ANALYSIS_DIR)
        if REMOVE_DUPLICATES_PAIRS:
            os.makedirs(ANALYSIS_DIR + '/dedup', exist_ok=True)
            for pair in FILE_PAIR_FORMATS: #becasue we want to removed duplicates by pairs
                if pair.startswith(filename):
                    print(pair)
                    stats_dedup = open(ANALYSIS_DIR+'/dedup/stats.txt','a')
                    stats_dedup.write('Pair format: {}\n'.format(pair))
                    file1 = '{}/{}/{}'.format(ANALYSIS_DIR,'{}',pair.format('1'))
                    file2 = '{}/{}/{}'.format(ANALYSIS_DIR,'{}',pair.format('2'))
                    if not exists(file1.format('dedup/')) and not exists(file2.format('dedup/')):
                        clump_path = os.path.abspath('/workspace/projects/damage_maps/damage_maps/software/bbmap_38.73/clumpify.sh')
                        command_clumpify = '{} threads={} in1={} in2={} out1={} out2={} dedupe'.format(clump_path,C, file1.format(''), file2.format(''), file1.format('dedup/'),file2.format('dedup/'))
                        print(command_clumpify)
                        clump_proc = subprocess.Popen(command_clumpify, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                        clump_proc.wait()
                        stats_dedup.write(get_proc_result(clump_proc)+'\n\n')
                        stats_dedup.close()
            ANALYSIS_DIR = ANALYSIS_DIR + '/dedup'
        os.chdir(MAIN_DIR)

        os.makedirs('{}/stats/'.format(ANALYSIS_DIR), exist_ok=True)
        os.makedirs('{}/mapping_hg19/unmapped/'.format(ANALYSIS_DIR), exist_ok=True)


        for pair in FILE_PAIR_FORMATS:  # becasue we want to removed duplicates by pairs
            if pair.startswith(filename):
                # ANALYSIS_DIR, 'FASTQ', file_name,
                s = open('{}/stats/{}-stats.txt'.format(ANALYSIS_DIR, pair.replace('{}.fastq.gz', 'paired.discarded_misligated')), 'a')
                s.write('Pair format: {}\n'.format(pair))
                file1 = pair.format('1')
                file2 = pair.format('2')
                if DISCARD_MISLIGATED_ADAPTER:
                    print('testing discarding adapters at start...')
                    file1, file2 = discard_misligated_adapters_at_start_paired(s, file1, file2, discard_adapter_5p)
                    s = open('{}/stats/{}-stats.txt'.format(ANALYSIS_DIR, pair.replace('{}.fastq.gz', 'paired.discarded_misligated')),
                             'a')

                file_sam = align_paired(file1, file2, s)
        bam_sort_path = generate_bam(file_sam,s)

        #sys.exit()
        print(filename, 'generating bed')
        #read_bed_file = save_damage_pos_bed(sam_folder, filename, bedtools_path, genome_path)
        read_bed_file = save_damage_pos_bed(bam_sort_path, bedtools_path, genome_path)
        #print(read_bed_file)
    sys.exit()
