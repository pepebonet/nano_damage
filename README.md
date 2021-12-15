# Nano analysis
Scripts to detect patterns of DNA methylation in yeast. 


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
        pip install -e .

# Usage

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

Once the damage in the nucleosomes is obtained the following commands need to be run to obtain the randomizations of the damage: 

        python scripts/randomizations_damage.py -dn outputs/retest_commands_Dec2021/cisplatin/damage_in_nucleosomes.tsv -ed outputs/retest_commands_Dec2021/cisplatin/triplet_normalized_Nanopore.tsv -p -o outputs/retest_commands_Dec2021/cisplatin/

## Zoom out 

# Example data