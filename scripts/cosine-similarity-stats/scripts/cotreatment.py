"""
Explains each signature as a mixture of co-administered matrix_treatments.
For each signature and the collection of co-administered matrix_treatments that are related to it,
we deem "efficiency" the capacity for a treatment to generate exposure of a specific signature.
"""

import os
import tqdm
import click
import numpy as np
import pandas as pd
import scipy.optimize
import dill as pickle
import gzip

from non_negative_least_squares import nnls_input
from synthetic import signature_recovery


class EfficiencyDrugSpecific:

    def __init__(self, exposures_path, matrix_treatments_path, signature_index):

        self.exposures = pd.read_csv(exposures_path, sep='\t')
        self.exposures = self.exposures.transpose()
        self.matrix_treatments = pd.read_csv(matrix_treatments_path, sep='\t', index_col=0)
        self.signature_index = signature_index

    def parse_exposure(self):

        return self.exposures.loc[:, [self.signature_index]]

    def treatfit(self, etiology, placebo):

        # exposures dataframe
        exposure = self.parse_exposure()
        treatments = self.matrix_treatments[[etiology, placebo]]

        # NNLS: transform the input
        E = exposure.values
        T = treatments.values

        A, b = nnls_input(T, E)

        # NNLS: custom bounds -- some values must be zero
        lb = np.zeros(A.shape[1])  # lower bound
        ub = lb + np.inf

        # NNLS: solver
        res = scipy.optimize.lsq_linear(A, b, bounds=(lb, ub))

        # collapse near-zero coefficients
        return dict(zip([etiology, placebo], list(res.x)))

    def treatfit_many_at_once(self, etiology, competitors):

        # exposures dataframe
        exposure = self.parse_exposure()
        treatments = self.matrix_treatments[[etiology] + competitors]

        # NNLS: transform the input
        E = exposure.values
        T = treatments.values

        A, b = nnls_input(T, E)

        # NNLS: custom bounds -- some values must be zero
        lb = np.zeros(A.shape[1])  # lower bound
        ub = lb + np.inf

        # NNLS: solver
        res = scipy.optimize.lsq_linear(A, b, bounds=(lb, ub))

        # collapse near-zero coefficients
        return dict(zip([etiology] + competitors, list(res.x)))

def origin(treatment, signature, exposures_path, discriminate_path, matrix_treatment_path):

    # show progress during execution
    print(treatment, signature)

    # determine efficiencies of competing matrix_treatments
    e = EfficiencyDrugSpecific(exposures_path, discriminate_path, matrix_treatment_path)

    # corner case 1
    if signature not in e.exposures.columns:
        return None, None

    # corner case 2
    if treatment not in e.matrix_treatments.columns:
        return None, None

    dftreatpool = e.select_pools(treatment, signature)

    # corner case 3
    if len(dftreatpool) == 0:
        return None, None

    x = e.treatfit(dftreatpool, [signature])
    x = dict(zip(list(x.index), list(x.iloc[:, 0].values)))

    # most likely origin
    origin = sorted(x.keys(), key=lambda k: x[k], reverse=True)[0]

    # efficiency rate between origin and reference treatment
    if (origin != treatment) and (x[treatment] == 0):
        score = np.inf
    else:
        score = x[origin] / x[treatment]
    return origin, score


def dump_efficiency_dict(exposures_path, discriminate_path, matrix_treatment_path, output):

    eff = EfficiencyDrugSpecific(exposures_path, discriminate_path, matrix_treatment_path)
    with gzip.open(output, 'wb') as fn:
        pickle.dump(eff, fn)


def create_random_matrix_treatments(etiology, n_treated, replicate, overlap, exposure_folder, output_folder, symmetric=True):

    exposure_path = os.path.join(exposure_folder,
                                 f'{etiology}.{n_treated}.{replicate}.catalogue.tsv',
                                 f'{etiology}.{n_treated}.{replicate}.catalogue.tsv.exposures.tsv')
    df = pd.read_csv(exposure_path, sep='\t')
    dg = df.transpose()
    dg[etiology] = dg.apply(lambda x: 1 if x.name.endswith('_t') else 0, axis=1)
    index_treated = dg.loc[dg.index.str.endswith('_t')].index.to_list()
    index_nontreated = dg.loc[~dg.index.str.endswith('_t')].index.to_list()
    for i in range(10):
        placebo_treated = np.random.choice(index_treated, size=len(index_treated) * overlap // 100, replace=False)
        if symmetric:
            placebo_nontreated = np.random.choice(index_nontreated, size=len(index_nontreated) * overlap // 100, replace=False)
        else:
            placebo_nontreated = np.random.choice(index_nontreated, size=n_treated - len(placebo_treated), replace=False)
        placebo = list(placebo_treated) + list(placebo_nontreated)
        dg[f'placebo_{i+1}'] = dg.apply(lambda x: int(x.name in placebo), axis=1)
    dg = dg[[etiology] + [f'placebo_{i+1}' for i in range(10)]]
    output_path = os.path.join(output_folder, f'{etiology}.{n_treated}.{replicate}.{overlap}.treatment_matrix.tsv')
    dg.to_csv(output_path, sep='\t', index=True)


@click.group()
def cli(): pass


def attribution(etiology, n_treated, replicate, deconstruction_folder, matrix_treatment_folder, overlap, signatures_path):

    exposures_path = os.path.join(deconstruction_folder,
                                  f'{etiology}.{n_treated}.{replicate}.catalogue.tsv',
                                  f'{etiology}.{n_treated}.{replicate}.catalogue.tsv.exposures.tsv')
    matrix_treatment_path = os.path.join(matrix_treatment_folder,
                                         f'{etiology}.{n_treated}.{replicate}.{overlap}.treatment_matrix.tsv')

    sig_index, cosine = signature_recovery(etiology, n_treated, replicate, signatures_path, deconstruction_folder)

    # determine efficiencies of competing matrix_treatments
    e = EfficiencyDrugSpecific(exposures_path, matrix_treatment_path, sig_index)
    res_dict = {etiology: [], 'placebo': []}
    for i in range(10):
        fit = e.treatfit(etiology, f'placebo_{i+1}')
        res_dict[etiology].append(fit[etiology])
        res_dict['placebo'].append(fit[f'placebo_{i+1}'])
    return res_dict


def attribution_many_at_once(etiology, n_treated, replicate, deconstruction_folder, matrix_treatment_folder, overlap, signatures_path):

    exposures_path = os.path.join(deconstruction_folder,
                                  f'{etiology}.{n_treated}.{replicate}.catalogue.tsv',
                                  f'{etiology}.{n_treated}.{replicate}.catalogue.tsv.exposures.tsv')
    matrix_treatment_path = os.path.join(matrix_treatment_folder,
                                         f'{etiology}.{n_treated}.{replicate}.{overlap}.treatment_matrix.tsv')

    sig_index, cosine = signature_recovery(etiology, n_treated, replicate, signatures_path, deconstruction_folder)

    # determine efficiencies of competing matrix_treatments
    e = EfficiencyDrugSpecific(exposures_path, matrix_treatment_path, sig_index)
    competitors = pd.read_csv(matrix_treatment_path, sep='\t', index_col=None)
    competitors = list(competitors.columns[2:])
    res_dict = {etiology: [], 'placebo': []}
    for i, c in enumerate(competitors):
        fit = e.treatfit_many_at_once(etiology, competitors[: i+1])  # must define a new treatfit that takes list arg
        res_dict[etiology].append(fit[etiology])
        del fit[etiology]
        max = np.max(list(fit.values()))
        res_dict['placebo'].append(max)
    return res_dict


def run_create_random_matrix_treatments():

    etiologies = {'capecitabine': 'SBS17b', 'platinum': 'SBS31', 'flat': 'SBS3',
                  'benzo': 'SBS4', 'poleta': 'SBS9', 'apobec': 'SBS13'}
    deconstruction_folder = '/workspace/users/fmuinos/treatment_effects/simulations/deconstruction/breastlungcolon'
    treatment_matrix_folder = '/home/fmuinos/projects/oriol/treatment_effects/simulations/treatment_matrix/asymmetric/breastlungcolon'

    for etiology in tqdm.tqdm(etiologies):
        for n_treated in tqdm.tqdm([10, 25, 50, 100, 150]):
            for replicate in range(1, 11):
                for overlap in [25, 50, 75]:
                    try:
                        create_random_matrix_treatments(etiology, n_treated, replicate, overlap, deconstruction_folder,
                                                        treatment_matrix_folder, symmetric=False)
                    except FileNotFoundError:
                        print(etiology, n_treated, replicate)


if __name__ == '__main__':

    deconstruction_folder = '/workspace/users/fmuinos/treatment_effects/simulations/deconstruction/'
    matrix_treatment_folder = '/home/fmuinos/projects/oriol/treatment_effects/simulations/treatment_matrix/asymmetric'
    signatures_path = '/home/fmuinos/projects/oriol/treatment_effects/simulations/signatures.pickle.gz'
    output = '/home/fmuinos/projects/oriol/treatment_effects/simulations/cotreatments/manyatonce/asymmetric'

    # run_create_random_matrix_treatments()

    d = attribution_many_at_once('capecitabine', 100, 1, deconstruction_folder, matrix_treatment_folder, 75, signatures_path)
    # d = attribution('capecitabine', 100, 1, deconstruction_folder, matrix_treatment_folder, 75, signatures_path)
    print(d)
