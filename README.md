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




# Example data