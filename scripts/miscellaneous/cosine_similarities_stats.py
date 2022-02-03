#!/usr/bin/envs python3

import os
import click
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm
import matplotlib.pyplot as plt
from scipy.spatial.distance import cosine


names = ['SemiSup Cisplatin', 'SemiSup MMS', 'Tombo MMS', 'Tombo Cisplatin',
        'Mao MMS', 'Sancar Cisplatin']

def dirichlet_generator_old(signature, cosines_df, n_draws, n_channels, alpha=1):
    
    signatures = signature.columns
    labels = []; means  = []; pvals  = []; cos = []
    

    for c in tqdm(signatures):
        sig = signature[c].values
        s_pool = np.random.dirichlet(
            alpha=[alpha for _ in range(n_channels)], size=n_draws
        )
        cosines = list(map(lambda x: 1-cosine(x, sig), s_pool))
        means.append(np.mean(cosines))
        # for el in tqdm(cosines_df[c].values):
        import pdb;pdb.set_trace()
        pvals.append(sum((np.array(cosines) >= 0.8)) / n_draws)
        labels.append(c); cos.append(0.8)
        

    logpvals = list(map(
        lambda x: -np.log10(x) if x > 0 else -np.log10(1/n_draws), pvals)
    )
    print(labels, means, logpvals, pvals)
    # import pdb;pdb.set_trace()
    return labels, logpvals, means



def dirichlet_generator(signature, cosines_df, n_draws, n_channels, alpha=1):
    
    signatures = signature.columns
    df = pd.DataFrame(columns=names, index=names, dtype='float64')

    for c in tqdm(signatures):
        sig = signature[c].values
        s_pool = np.random.dirichlet(
            alpha=[alpha for _ in range(n_channels)], size=n_draws
        )
        cosines = list(map(lambda x: 1-cosine(x, sig), s_pool))
        for el in cosines_df.index:
            if el != 'Mao MMS':
                pval = sum((np.array(cosines) >= cosines_df[c][el])) / n_draws
                logpval = [-np.log10(pval) if pval > 0 else -np.log10(1/n_draws)][0]
                df[c][el] = logpval
    import pdb;pdb.set_trace()

    # import pdb;pdb.set_trace()
    return df


def dirichlet_vs_dirichlet(n_draws, n_channels, alpha=1):
    
    cosines = []
    sig_pool_1 = np.random.dirichlet(alpha=[alpha for _ in range(n_channels)], size=n_draws)
    sig_pool_2 = np.random.dirichlet(alpha=[alpha for _ in range(n_channels)], size=n_draws)
    cosines = list(map(lambda x, y: 1-cosine(x, y), sig_pool_1, sig_pool_2))
    import pdb;pdb.set_trace()
    return np.array(cosines)


def test_alphas(channels, label, output):

    means = []
    pvals = []
    N = 10000
    x = np.linspace(0.1, 3.5, num=100)
    for alpha in tqdm(x):
        cosines = dirichlet_vs_dirichlet(N, channels, alpha=alpha)
        pvals.append(np.sum(cosines >= 0.85) / N)
        means.append(np.mean(cosines))

    logpvals = list(map(lambda x: -np.log10(x) if x > 0 else -np.log10(1/N), pvals))

    fig, ax1 = plt.subplots()

    ax2 = ax1.twinx()
    ax1.plot(x, logpvals, 'r-')
    ax2.plot(x, means, 'b-')

    ax1.set_xlabel('alpha')
    ax1.set_ylabel('-log10 pvals', color='r')
    ax2.set_ylabel('mean cosine similarity', color='b')

    outdir = os.path.join(output, f'alphas_{label}.png')
    fig.tight_layout()

    plt.savefig(outdir)
    plt.close()


def do_plots(labels, logpvals, means, alpha, label, output):

    fig, ax1 = plt.subplots(figsize=(20, 5))

    ax2 = ax1.twinx()

    x = range(len(labels))

    ax1.plot(x, logpvals, 'ro-')
    ax2.plot(x, means, 'bo-')

    ax1.set_xlabel('signature')
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, rotation=90)
    ax1.hlines(2, 0, len(logpvals), linestyles='dashed')

    for t in x:
        ax1.vlines(t, 0, 4, linestyles='dashed', alpha=0.1)

    ax1.set_ylabel('-log10 pvals', color='r')
    ax1.set_ylim(-0.1, 4.1)
    ax2.set_ylabel('mean cosine similarity', color='b')
    ax2.set_ylim(-0.1, 1.1)
    plt.title(f'Prob(cosine$\geq$0.8)\n Symmetric Dirichlet: a={alpha}')
    
    outdir = os.path.join(output, f'distributions_{label}.png')
    fig.tight_layout()

    plt.savefig(outdir)
    plt.close()


@click.command(short_help='Get stats from cosines')
@click.option('-sa', '--signatures_all', required=True)
@click.option('-sm', '--signatures_mao', required=True)
@click.option('-dc', '--data_cosines', required=True)
@click.option('-o', '--output', required=True)
def main(signatures_all, signatures_mao, data_cosines, output):

    sig = pd.read_csv(signatures_all, sep='\t')
    sig_mao = pd.read_csv(signatures_mao, sep='\t')
    df = pd.read_csv(data_cosines, sep='\t', index_col=0)

    alpha = 2
    n = 10000

    labels, logpvals, means = dirichlet_generator(
        sig, df, n, 64, alpha=alpha
    )
    do_plots(labels, logpvals, means, alpha, 'all', output)

    labels_mao, logpvals_mao, means_mao = dirichlet_generator_mao(
        sig_mao, df, n, 32, alpha=alpha
    )
    do_plots(labels_mao, logpvals_mao, means_mao, alpha, 'mao', output)


    test_alphas(64, 'all', output)
    test_alphas(32, 'mao', output)
    import pdb;pdb.set_trace()


if __name__ == '__main__':
    main()