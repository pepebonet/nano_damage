# Nano analysis
Recently, Nanopore long-read sequencing has emerged as a novel technology to enable direct detection of DNA modifications. Here, the idea is to explore the capabilities of Nanopore sequencing data to directly map a variety of DNA damages to the yeast genome in an unsupervised manner. The objective of the scripts is to shed light on the challenging question of whether it is possible to detect damage without previous knowledge and with a low damage proportion. For that we generate the pipeline shown in the following figure: 

<img src="figures/Figure.png" alt="alt text" width=1000 height="whatever">

# Contents
- [Installation](#Installation)
- [Usage](#Usage)
- [Example data](#Example-data)         

# Installation
## Clone repository
First download the repository:

        git clone https://github.com/pepebonet/nano_damage.git

## Install dependencies
We highly recommend to use a virtual environment to run the scripts: 

`Create environment and install nano_damage:`

        conda create --name nanodamage python=3.8
        conda activate nanodamage
        pip install -r requirements.txt

# Usage

## Handle Raw Fast5 files

The objective here is to go from multi-fast5 to resquiggled fast5 files. To do so three things need to be done. 

1.- Get the subset of reads of your treatment:

        fast5_subset --input path_to_input/ --save_path path_to_output/ --read_id_list ../sequencing_summary_treatment.txt --batch_size 2000 -t 56

2.- Get single reads: 

        multi_to_single_fast5 --input_path multi_reads/ --save_path single_reads/ --threads 56

3.- Resquiggle with Tombo: 
        tombo resquiggle single_read/ path_to_reference_genome/ --processes 56 --num-most-common-errors 5

These can also be done by adapting the Snakemake pipeline in the `simulation` folder. Changes to both the Snakefile and config.yaml will be needed. 
        
        snakemake --cores 4


## Enrichment

For the folder of the available data `Desktop/samples_nanopore_bbg_novoa` the enrichment script can be run for the different treatments as follows: 

        python scripts/enrichment.py -o outputs/retest_comnds_Dec2021/mms/ -cl  sacCer3.chrom.sizes -d ../samples_nanopore_bbg_novoa/first_batch_repeated/merged_treated_mms/denovo.tsv.gz

        python scripts/enrichment.py -o outputs/retest_comnds_Dec2021/aag/ -cl  sacCer3.chrom.sizes -d ../samples_nanopore_bbg_novoa/first_batch_repeated/merged_treated_aag/denovo.tsv.gz

        python scripts/enrichment.py -o outputs/retest_commands_Dec2021/ems/ -cl  sacCer3.chrom.sizes -d ../samples_nanopore_bbg_novoa/FAL40053.EMS.mf01.P01/denovo.tsv.gz

        python scripts/enrichment.py -o outputs/retest_comnds_Dec2021/cisplatin/ -cl  sacCer3.chrom.sizes -d ../samples_nanopore_bbg_novoa/FAL40053.cisplatin.mf01.P01/denovo.tsv.gz


## Nucleosomes Norm

After the damage and the enrichment is obtained we want to focus only in the nucleosome. Therefore, the commands below are run: 

        python scripts/nucleosome_damage_norm.py -da outputs/retest_commands_Dec2021/mms/damage_most_sig_Nanopore.tsv -ni nucleosome_info/saccer3/brogaard_saccer3.bed.nooverlapp.bed.gz -ed outputs/retest_commands_Dec2021/mms/triplet_normalized_Nanopore.tsv -o outputs/retest_commands_Dec2021/mms/

        python scripts/nucleosome_damage_norm.py -da outputs/retest_commands_Dec2021/aag/damage_most_sig_Nanopore.tsv -ni nucleosome_info/saccer3/brogaard_saccer3.bed.nooverlapp.bed.gz -ed outputs/retest_commands_Dec2021/aag/triplet_normalized_Nanopore.tsv -o outputs/retest_commands_Dec2021/aag/

        python scripts/nucleosome_damage_norm.py -da outputs/retest_commands_Dec2021/ems/damage_most_sig_Nanopore.tsv -ni nucleosome_info/saccer3/brogaard_saccer3.bed.nooverlapp.bed.gz -ed outputs/retest_commands_Dec2021/ems/triplet_normalized_Nanopore.tsv -o outputs/retest_commands_Dec2021/ems/

        python scripts/nucleosome_damage_norm.py -da outputs/retest_commands_Dec2021/cisplatin/damage_most_sig_Nanopore.tsv -ni nucleosome_info/saccer3/brogaard_saccer3.bed.nooverlapp.bed.gz -ed outputs/retest_commands_Dec2021/cisplatin/triplet_normalized_Nanopore.tsv -o outputs/retest_commands_Dec2021/cisplatin/

## Randomizations

Once the damage in the nucleosomes is obtained the following commands need to be run to obtain the randomizations of the damage (remember the -nr parameter for the number of randomizations to reduce the running time): 

        python scripts/randomizations_damage.py -dn outputs/retest_commands_Dec2021/mms/damage_in_nucleosomes.tsv -ed outputs/retest_commands_Dec2021/mms/triplet_normalized_Nanopore.tsv -p -o outputs/retest_commands_Dec2021/mms/

        python scripts/randomizations_damage.py -dn outputs/retest_commands_Dec2021/aag/damage_in_nucleosomes.tsv -ed outputs/retest_commands_Dec2021/aag/triplet_normalized_Nanopore.tsv -p -o outputs/retest_commands_Dec2021/aag/

        python scripts/randomizations_damage.py -dn outputs/retest_commands_Dec2021/ems/damage_in_nucleosomes.tsv -ed outputs/retest_commands_Dec2021/ems/triplet_normalized_Nanopore.tsv -p -o outputs/retest_commands_Dec2021/ems/ 

        python scripts/randomizations_damage.py -dn outputs/retest_commands_Dec2021/cisplatin/damage_in_nucleosomes.tsv -ed outputs/retest_commands_Dec2021/cisplatin/triplet_normalized_Nanopore.tsv -p -o outputs/retest_commands_Dec2021/cisplatin/

## Zoom out 

On top of looking at the nucleosomes, we can also look at differences between nucleosome and linkers and check for the periodicity. To do so, the following commands are needed: 

        python scripts/zoom_out_nuc.py -da outputs/retest_commands_Dec2021/mms/damage_most_sig_Nanopore.tsv -ni nucleosome_info/saccer3/brogaard_saccer3.bed.nooverlapp.bed.gz -o outputs/retest_commands_Dec2021/mms/ -ed outputs/retest_commands_Dec2021/mms/triplet_normalized_Nanopore.tsv

        python scripts/zoom_out_nuc.py -da outputs/retest_commands_Dec2021/aag/damage_most_sig_Nanopore.tsv -ni nucleosome_info/saccer3/brogaard_saccer3.bed.nooverlapp.bed.gz -o outputs/retest_commands_Dec2021/aag/ -ed outputs/retest_commands_Dec2021/aag/triplet_normalized_Nanopore.tsv

        python scripts/zoom_out_nuc.py -da outputs/retest_commands_Dec2021/ems/damage_most_sig_Nanopore.tsv -ni nucleosome_info/saccer3/brogaard_saccer3.bed.nooverlapp.bed.gz -o outputs/retest_commands_Dec2021/ems/ -ed outputs/retest_commands_Dec2021/ems/triplet_normalized_Nanopore.tsv

        python scripts/zoom_out_nuc.py -da outputs/retest_commands_Dec2021/cisplatin/damage_most_sig_Nanopore.tsv -ni nucleosome_info/saccer3/brogaard_saccer3.bed.nooverlapp.bed.gz -o outputs/retest_commands_Dec2021/cisplatin/ -ed outputs/retest_commands_Dec2021/cisplatin/triplet_normalized_Nanopore.tsv


## Sancar Data

Commands to run the sancar data properly: 

        python scripts/miscellaneous/sancar_analysis.py -ds /workspace/projects/nanopore/sancar_data/SRR3623538.1.mapped.bam.sort_normal -o . 

# Example data