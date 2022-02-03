"""
Given a treatment, it gives a collection of signatures
that may be related to the treatment after adjusting by the tumor type.
"""

import os
import glob
import click
import json
from tqdm import tqdm
from contextlib import suppress

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from scipy.stats import zscore

import warnings
warnings.filterwarnings('ignore')


np.random.seed(42)


def parse_treat(raw_regression_table, randomization, norm=True):

    samples = list(zip(*randomization))[0]
    treat_df = raw_regression_table.loc[samples, :]

    # reponse variable
    y = treat_df['resp'].values.astype(np.float64)
    treat_df.drop(columns=['resp'], inplace=True)
    treat_df.drop(columns=['tumor_type'], inplace=True)

    if norm:
        # log + z-score transformation
        treat_df.loc[:, :] = np.log(treat_df.values + 1)
        treat_df.loc[:, :] = zscore(treat_df.values, axis=0)

    return treat_df, y


def run_logistic(raw_regression_table, randomization):

    treat_df, y = parse_treat(raw_regression_table, randomization)
    X = treat_df.values.astype(np.float64)
    logit = LogisticRegression(random_state=42, fit_intercept=True,
                               solver='lbfgs', penalty='l2', C=0.1).fit(X, y)
    dict_res = dict((zip(list(treat_df.columns), logit.coef_[0])))
    return dict_res


class Urn:

    """
    Urn class provides with the container and methods to randomly sample
    balanced sets of treated vs untreated samples on a tumor-type by tumor-type basis.
    """

    def __init__(self, toregression):

        dg = toregression
        # tumor_type labels for each sample
        def label_tumor_type(x):
            cancer_label = x.name.split('_')[-2]
            if cancer_label.endswith('AdenoCA'):
                return '_'.join(x.name.split('_')[-3:-1])
            else:
                return 'Breast_AdenoCA'
        dg['tumor_type'] = dg.apply(label_tumor_type, axis=1)
        # response label
        dg['resp'] = toregression['resp'].values
        # bag of tumor-type and response values
        self.bag = dg[['tumor_type', 'resp']]

    @staticmethod
    def balanced_sample(df):
        """sample a balanced set of indices with replacement"""

        df_one  = df[df["resp"] == 1]
        df_zero = df[df["resp"] == 0]

        # split samples according to the number of treated patients.
        parts = 3
        if len(df_one) < 50:
            parts = 2
        m1, m2 = len(df_one), len(df_zero)
        n = min(m1, m2) // parts
        if n > 0:
            ones    = np.random.choice(list(df_one.index), size=n, replace=False)
            zeros   = np.random.choice(list(df_zero.index), size=n, replace=False)
            samples = list(ones) + list(zeros)
        else:
            samples = []

        return [(s, df.loc[s, 'resp']) for s in samples]

    def stratified_random_sample(self, tumor_type='Pan'):
        """
        Generate a treated vs untreated balanced set of indices
        on a tumor-type by tumor-type basis
        """

        pool = []
        if tumor_type == 'Pan':
            for ttype in self.bag['tumor_type'].unique():
                pool += self.balanced_sample(self.bag[self.bag['tumor_type'] == ttype])
        else:
            return self.balanced_sample(self.bag[self.bag['tumor_type'] == tumor_type])
        return pool

    def randomize(self, n=1000, tumor_type='Pan'):

        np.random.seed(42)
        return [self.stratified_random_sample(tumor_type=tumor_type) for _ in range(n)]


def merge_dicts(dict_iterable):

    d = {}
    for k in dict_iterable[0].keys():
        d[k] = list(dict_[k] for dict_ in dict_iterable)
    return d


def effect_size(rand, to_regression):

    zeros, ones = {}, {}
    effect_size_list = []

    for sample in rand:

        treat_df, y = parse_treat(to_regression, sample, norm=False)

        df = treat_df.iloc[(y == 1), :]
        for sign in df.columns:
            ones[sign] = ones.get(sign, []) + df.loc[:, sign].values.tolist()

        dg = treat_df.iloc[(y == 0), :]
        for sign in dg.columns:
            zeros[sign] = zeros.get(sign, []) + dg.loc[:, sign].values.tolist()

        effect = {sign: fold_change(zeros[sign], ones[sign]) for sign in zeros}
        effect_size_list.append(effect)

    effect_dict = merge_dicts(effect_size_list)
    effect_dict = {k: np.mean(v) for k, v in effect_dict.items()}

    return effect_dict


def fold_change(zeros, ones):

    eps = 1
    m_zeros = np.mean(zeros)
    m_ones  = np.mean(ones)  # rule out outlier values

    if m_zeros > eps:
        return m_ones / m_zeros
    else:
        return 0


def empirical_pvalue(values):

    assert (len(values) > 0)
    return len([v for v in values if v <= 0]) / len(values)


def return_metadata(path_metadata):

    # this comes after manually fixing all the data
    metadata = pd.read_csv(path_metadata, sep='\t')
    metadata['primaryTumorLocation_fixed'] = metadata['primaryTumorLocation'].apply(lambda x: str(x).replace(' ', '-').replace('/', '-'))
    forbidden_primary = []  # ['Unknown', 'nan', 'Double_primary', np.nan]
    metadata = metadata[~metadata['primaryTumorLocation_fixed'].isin(forbidden_primary)]
    dic_primary_full = dict(zip(metadata.sampleId, metadata.primaryTumorLocation_fixed))

    return dic_primary_full


def create_matrix_treatments(etiology, n_treated, replicate, deconstruction_folder, shuffle=None):

    if shuffle is None:
        fol_name = f'{etiology}.{n_treated}.{replicate}.catalogue.tsv'
        exposure = pd.read_csv(os.path.join(deconstruction_folder, f'{fol_name}/{fol_name}.exposures.tsv'), sep='\t')
    else:
        fol_name = f'{etiology}.{n_treated}.{replicate}.catalogue.tsv'
        exposure = pd.read_csv(os.path.join(deconstruction_folder, f'{fol_name}/shuffle.{shuffle}.{fol_name}.exposures.tsv'), sep='\t')

    toregression = exposure.transpose()
    # add response
    toregression['resp'] = toregression.apply(lambda r: int(r.name.endswith('_t')), axis=1)
    return toregression


@click.group()
def cli(): pass


def ensemble_regression(etiology, n_treated, replicate, output_folder, deconstruction_folder, shuffle=None):

    # return the matrix, the signatures list and the tumor type list
    raw_regression_table = create_matrix_treatments(etiology, n_treated, replicate, deconstruction_folder, shuffle=shuffle)

    # shuffle treatment labels

    # randomize
    urn = Urn(raw_regression_table)
    rand = urn.randomize(n=1000)

    # run logistic regression for each of the randomizations
    pool = []
    for randomization in rand:
        with suppress(Exception):
            dict_res = run_logistic(raw_regression_table, randomization)
            pool.append(dict_res)

    # get results
    merge_res = merge_dicts(pool)
    if shuffle is None:
        with open(os.path.join(output_folder, f'regression_lor.{etiology}.{n_treated}.{replicate}.json'), 'wt') as f:
            json.dump(merge_res, f)
    else:
        with open(os.path.join(output_folder, f'shuffle.{shuffle}.regression_lor.{etiology}.{n_treated}.{replicate}.json'), 'wt') as f:
            json.dump(merge_res, f)
    # get the effect size
    effect = effect_size(rand, raw_regression_table)  # fold change
    pvals  = {k: empirical_pvalue(v) for k, v in merge_res.items()}

    res = []
    for s, e in effect.items():
        d = {'signature': [s],
             'effect_size': [e],
             'pvals': [pvals[s]]}
        res.append(pd.DataFrame(d))

    res = pd.concat(res, axis=0)
    res.reset_index(drop=True, inplace=True)
    if shuffle is None:
        res.to_csv(os.path.join(output_folder, f'regression_summary.{etiology}.{n_treated}.{replicate}.tsv'),
                   sep='\t', index=False)
    else:
        res.to_csv(os.path.join(output_folder, f'shuffle.{shuffle}.regression_summary.{etiology}.{n_treated}.{replicate}.tsv'),
                   sep='\t', index=False)


def treatment_shuffle(etiology, deconstruction_folder):

    np.random.seed(42)
    for fn in tqdm(glob.glob(os.path.join(deconstruction_folder, f'{etiology}.*.catalogue.tsv/*.exposures.tsv'))):
        for falsenotreat in [5, 10, 20]:
            folder = os.path.dirname(os.path.normpath(fn))
            df = pd.read_csv(fn, sep='\t')
            samples = list(df.columns)
            new_samples = []
            for sample in samples:
                if sample.endswith('_t'):
                    a = np.random.choice([0, 1], p=[falsenotreat/100, 1-(falsenotreat/100)])
                    if a == 0:
                        sample = sample[:-2]
                new_samples.append(sample)
            df.columns = new_samples
            output = os.path.join(folder, f'shuffle.{falsenotreat}' + '.' + os.path.basename(fn))
            df.to_csv(output, sep='\t', index=False)


@cli.command()
@click.option('--folder', type=click.Path())
def merge(folder):

    for mut in tqdm(['snv', 'dbs', 'indels']):
        pool = []
        for fn in glob.glob(os.path.join(folder, f'results_regression_*_{mut}.tsv')):
            df = pd.read_csv(fn, sep='\t')
            pool.append(df)
        df = pd.concat(pool)
        df.to_csv(os.path.join(folder, f'regression_results_treatments_{mut}.tsv'), sep='\t', index=False)


def run_ensemble_regression(deconstruction_folder, regression_folder, shuffle=False, etiology=None):

    etiologies = {'capecitabine': 'SBS17b', 'platinum': 'SBS31', 'flat': 'SBS3',
                  'benzo': 'SBS4', 'poleta': 'SBS9', 'apobec': 'SBS13'}

    if etiology is not None:
        etiologies = [etiology]

    if shuffle:
        shuffle_iterable = [5, 10, 20]
    else:
        shuffle_iterable = [None]

    for etiology in tqdm(etiologies):
        for n_treated in tqdm([10, 25, 50, 100, 150]):
            for replicate in tqdm(range(0, 25)):
                for shuffle in shuffle_iterable:
                    try:
                        ensemble_regression(etiology, n_treated, replicate, regression_folder,
                                            deconstruction_folder, shuffle=shuffle)
                    except:
                        pass


if __name__ == '__main__':

    deconstruction_folder = '/workspace/users/fmuinos/treatment_effects/simulations/deconstruction/breastlungcolon/'
    regression_folder = '/home/fmuinos/projects/oriol/treatment_effects/simulations/regression_results/breastlungcolon/'

    run_ensemble_regression(deconstruction_folder, regression_folder, etiology='poleta', shuffle=True)

    # treatment_shuffle('poleta', deconstruction_folder)
