#!/usr/bin/envs python3

import os
import click
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.spatial.distance import cosine

def dirichlet_generator(signature, n_draws, n_channels, alpha=1, method='sp'):
    
    labels = []; means  = []; pvals  = []

    sig = signature[''].values
    s_pool = np.random.dirichlet(alpha=[alpha for _ in range(n_channels)], size=n_draws)
    cosines = list(map(lambda x: 1-cosine(x, sig), s_pool))
    means.append(np.mean(cosines))
    pvals.append(sum((np.array(cosines) >= 0.8)) / n_draws)
    labels.append(c)
    logpvals = list(map(lambda x: -np.log10(x) if x > 0 else -np.log10(1/n_draws), pvals))
    return labels, logpvals, means


@click.command(short_help='Get comparison all and heatmap')
@click.option('-nm', '--novoa_mms', required=True)
@click.option('-nc', '--novoa_cis', required=True)
@click.option('-mm', '--mao_mms', required=True)
@click.option('-sc', '--sancar_cis', required=True)
@click.option('-tm', '--tombo_mms', required=True)
@click.option('-tc', '--tombo_cis', required=True)
@click.option('-dc', '--data_cosines', required=True)
@click.option('-o', '--output', required=True)
def main(novoa_cis, novoa_mms, mao_mms, sancar_cis, tombo_mms, 
    tombo_cis, data_cosines, output):

    df = pd.read_csv(data_cosines, sep='\t', index_col=0)


if __name__ == '__main__':
    main()