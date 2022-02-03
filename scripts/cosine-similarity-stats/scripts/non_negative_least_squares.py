"""Util functions for conducting non-negative least squares regression"""


import scipy
import numpy as np
import pandas as pd


def create_block(row, n_repeats):
    """
    Returns:
        matrix: n_repeats x [length(row) * n_repeats]
    """

    block = np.zeros((n_repeats, len(row) * n_repeats))
    for i in range(block.shape[0]):
        block[i, i * len(row): (i + 1) * len(row)] = row
    return block


def nnls_input(T, E):
    """
    Args:
        T: samples x matrix_treatments matrix with binary values
        E: samples x signatures matrix with float values
    Returns:
        A: design matrix for non-negative least squares setting
        b: target vector (length = samples * signatures) for non-negative least squares setting
    """

    # create design matrix
    A = []
    n_signatures = E.shape[1]
    for i in range(T.shape[0]):
        row = T[i, :]
        block = create_block(row, n_signatures)
        A.append(block)
    A = np.concatenate(tuple(A), axis=0)
    A[~np.isfinite(A)] = 0

    # create target vector
    b = np.reshape(E, (E.shape[0] * E.shape[1], 1), order='F')
    b = b.flatten()
    b[~np.isfinite(b)] = 0

    return A, b


def nnls_input_bootstrap(T, E, repeat_size=1000, cv_fold=3):
    """
    Args:
        T: samples x matrix_treatments matrix with binary values
        E: samples x signatures matrix with float values
    Returns:
        [A]: list of design matrices for non-negative least squares setting
        [b]: list of target vectors (length = samples * signatures) for non-negative least squares setting
    """

    design_matrix_list = []
    target_vector_list = []

    n_signatures = E.shape[1]
    n_samples = T.shape[0]

    # bootstrap parameters: cross-validation size; number of repeats
    cv_size = n_samples // cv_fold

    for _ in range(repeat_size):

        A = []
        randomization = np.random.choice(n_samples, size=cv_size, replace=False)

        # create the design matrix
        for i in randomization:
            row = T[i, :]
            block = create_block(row, n_signatures)
            A.append(block)
        A = np.concatenate(tuple(A), axis=0)
        A[~np.isfinite(A)] = 0

        # create target vector
        E_subset = E[randomization, :]
        b = np.reshape(E_subset, (E_subset.shape[0] * E_subset.shape[1], 1), order='F')
        b = b.flatten()
        b[~np.isfinite(b)] = 0

        # output lists
        design_matrix_list.append(A)
        target_vector_list.append(b)

    return design_matrix_list, target_vector_list


def matrix_test():

    I = pd.DataFrame({'sig1': [True, False],
                      'sig2': [False, True],
                      'sig3': [True, True]},
                     index=['treat1', 'treat2'])
    T = np.array([[1, 0], [1, 1], [0, 1], [1, 1]])
    E = np.array([[5, 5, 5], [3, 2, 5], [3, 5, 2], [2, 2, 0]])
    A, b = nnls_input(T, E)

    lb = np.zeros(A.shape[1])  # lower bound

    mask = np.reshape(I.values, (len(I.index) * len(I.columns), 1), order='F')
    mask = list(mask.flatten())

    ub = lb + np.inf
    ub[np.logical_not(mask)] = 0.1  # upper bound

    # NNLS: solver
    res = scipy.optimize.lsq_linear(A, b, bounds=(lb, ub))

    # format solution as a dataframe
    x_matrix = np.reshape(res.x, (T.shape[1], E.shape[1]), order='F')

    print(x_matrix)