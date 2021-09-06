import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import cmath
import click
import os
import json
import gzip
from tqdm import tqdm
from functools import partial, reduce
from multiprocessing import Pool
from miscellaneous.non_linear import create_quadratic_model, NonLinearFitting


""" Utils """


def sliding_window(seq, window, stepsize=1):
    """
    Args:
        seq: np.array
        window: int: window size
        stepsize: int: step size
    Returns:
        np.array of 1-step windows produced from seq
    """
    sliding_list = []
    rolled = seq
    n_steps = (len(seq) - window) // stepsize
    for _ in range(n_steps + 1):
        sliding_list.append(rolled[: window])
        rolled = np.roll(rolled, -stepsize)
    return np.array(sliding_list)


def mean_smooth(seq, window):
    sliding_array = sliding_window(seq, window)
    return list(map(np.mean, sliding_array))


"""Discrete-time Fourier Transform"""


class FourierAnalysis(object):

    def __init__(self, signal):
        self.signal = (signal - np.mean(signal)).flatten()
        self.L = len(signal)

    def dtft(self, smin, smax):
        def power(period):
            """reference: https://www.youtube.com/watch?v=Qs-Zai0F2Pw"""
            return (1 / (smax - smin)) * abs(sum(self.signal[s] * cmath.exp(-1j * s * 2 * np.pi / period) for s in range(smin, smax)))**2

        def phase(period):
            return cmath.phase(sum(self.signal[s] * cmath.exp(-1j * s * 2 * np.pi / period) for s in range(smin, smax)))
        return power, phase

    @staticmethod
    def normalize(func, pmin, pmax):
        def normalized_func(period):
            return (pmax - pmin) * func(period) / sum(func(p) for p in range(pmin, pmax))
        return normalized_func


def normalize(func, pmin, pmax):
    def normalized_func(period):
        return (pmax - pmin) * func(period) / sum(func(p) for p in range(pmin, pmax))
    return normalized_func


def compute_fourier(signal, **kwargs):
    """
    Args:
        signal: numpy array
    Returns:
        function: normalized DTFT (discrete-time Fourier transform) power spectrum
    """

    # quadratic least-squares de-trending
    initial_values_dict = {'a_0': np.mean(signal), 'a_1': 0., 'a_2': 0.}
    params, obj_func = create_quadratic_model(initial_values_dict)
    x = np.arange(len(signal))
    non_linear_fitter = NonLinearFitting(obj_func, params, x, signal)
    result, predicted = non_linear_fitter.least_squares()
    signal_detrended = signal - predicted

    # compute DTFT
    fourier = FourierAnalysis(signal_detrended[:])
    power, phase = fourier.dtft(kwargs['low_t'], kwargs['high_t'])
    return power


def compute_spectrum(signal, norm=True, center=None, **kwargs):
    """
    Args:
        signal: numpy array with discrete-time signal values
        norm: whether to provide normalized spectrogram
        kwargs: parameters to run compute_fourier
    Returns:
        tuple:
            x-values
            y-values: DTFT of x-values.
            snr: float: signal-to-noise ratio
            peak: float: period of maximum power
    """
    signal = mean_smooth(signal, 3)  # 3-bp average smoothing reduces effects of 3-bp periodic signal
    power = compute_fourier(signal, **kwargs)  # compute fourier encloses
    if norm:
        power = normalize(power, kwargs['low_p'], kwargs['high_p'])
    x = np.linspace(kwargs['low_p'], kwargs['high_p'], 100)
    y = list(map(power, x))
    snr, peak = compute_signal_to_noise(x, y, center=center)
    return x, y, snr, peak


def closest_index(x, value):
    """
    Args:
        x: array: x-values
        value: float:
    Returns:
        index i of x such at x[i] is closest to value
    """
    return np.argmin(np.abs(x-value))


def compute_signal_to_noise(x, y, center=None):
    """
    Args:
        x: array: x-values
        y: array: DTFT of x-values
        center=None: center period for signal-to-noise analysis
                     default value gives maximum of power spectrum
    Returns:
        signal-to-noise considering foreground interval and its complementary
        period where the spectrogram peaks
    Source:
        See Lehmann, Machne and Herzel NAR 2014:
        ~/home/fmuinos/bibliography/supplement_methods_structural_code_cyanobacterial_genomes.pdf
        We compute the SNR by centering at the period peaking in the spectrogram.
    """
    x_i, y_i = [], []
    background = []
    # x_i, y_i will keep the values in the peak window
    # background will keep the y-values outside the peak window
    
    assert(len(x) == len(y))
    
    if center is None:
        peak = x[np.argmax(y)]
    else:
        peak = center
    
    for i, v in enumerate(x):
        if peak - 0.5 <= v <= peak + 0.5:
            x_i.append(x[i])
            y_i.append(y[i])
        else:
            background.append(y[i])
    
    x_i, y_i = np.array(x_i), np.array(y_i)
    y_argmax = np.argmax(y_i)
    
    if (y_argmax == 0) or (y_argmax == len(y_i) - 1):
        return y_i[closest_index(x_i, peak)] / np.median(background), x_i[closest_index(x_i, peak)]
    else:
        return max(y_i) / np.median(background), x_i[np.argmax(y_i)]


def build_species_tables(species, folder, center=None):
    motifs = ['aa', 'ac', 'ag', 'at', 'at2',
              'ca', 'cc', 'cg', 'ct',
              'ga', 'gc', 'gg', 'gt',
              'ta', 'tc', 'tg', 'tt',
              'ss', 'ww']
    df = pd.DataFrame(columns=motifs)
    dg = pd.DataFrame(columns=motifs)
    for motif in motifs:
        try:
            file_name = os.path.join(folder, 'pair_counts.{0}.{1}.json'.format(motif, species))
            with open(file_name, 'rt') as f:
                mydict = json.load(f)
            pair_count = np.array(mydict['pair_count'])
            motif_count = mydict['motif_count']
            len_chunk = mydict['chunk_len']
            signal = np.array([(v / len_chunk) / (motif_count / len_chunk) ** 2 for v in pair_count])
            x, y, snr, peak = compute_spectrum(signal, center=center, low_t=30, high_t=100, low_p=5, high_p=20)
            df.loc[species, motif] = snr
            dg.loc[species, motif] = peak
        except:
            df.loc[species, motif] = np.nan
            dg.loc[species, motif] = np.nan
    return df, dg


def motif_tables(chunk_id, folder, center=None, norm=True):
    """
    Args:
        chunk_id: chunk id that is part of the file name
        center: either None or float value
    """
    motifs = ['ww', 'at2', 'ss']
    df = pd.DataFrame(columns=motifs, index=[chunk_id])
    for motif in motifs:
        file_name = os.path.join(folder, '.'.join([chunk_id, motif, 'json.gz']))
        try:
            with gzip.open(file_name, 'rb') as f:
                mydict = json.load(f)
            pair_count = np.array(mydict['pair_count'])
            motif_count = mydict['motif_count']
            len_chunk = mydict['chunk_len']
            signal = np.array([(v / len_chunk) / (motif_count / len_chunk) ** 2 for v in pair_count])
            x, y, snr, peak = compute_spectrum(signal, norm=norm, center=center, low_t=30, high_t=100, low_p=5, high_p=20)
            df.loc[chunk_id, motif] = snr
        except:
            df.loc[chunk_id, motif] = np.nan
    return df


@click.command(context_settings={'help_option_names': ['-h', '--help'], 'max_content_width': 200})
@click.option('--folder', type=click.STRING, help='path to folder with chunk counts', required=True)
@click.option('--cores', type=click.INT, help='number of workers in multiprocessing', required=False, default=1)
@click.option('--label', type=click.STRING, help='label for the output tables', required=True)
@click.option('--center', type=click.FLOAT, help='center for SNR', required=False, default=None)
@click.option('--norm', is_flag=True, help='whether the Fourier transform is normalized', required=False)
def build_tables(folder, cores, label, center, norm):
    """
    \b
    bgqmap run -m 20G -c 50 \\
    'python spectral.py --folder /workspace/projects/periodicity/data/species_genomes/Eukarya/hg19/json_10k/chr1/ \\
    --cores 50 \\
    --label homo_chr1_10k_amplitude \\
    --center 10.1'
    """
    chunk_id_list = []
    for file in os.listdir(folder):
        chunk_id = '.'.join(file.split('.')[:2])
        if chunk_id not in chunk_id_list:
            chunk_id_list.append(chunk_id)
    print(chunk_id_list)
    func = partial(motif_tables, folder=folder, center=center, norm=norm)
    df_list, dg_list = [], []
    with Pool(cores) as pool:
        for a in tqdm(pool.imap(func, chunk_id_list[:100])):
            df_list.append(a)
    df = pd.concat(df_list, axis=0)
    df.to_csv('/workspace/users/fmuinos/periodicity/{0}_snr_{1}.tsv'.format(label, str(center)), sep='\t')


def dtft_spectrum_plot(signal, ax, norm=True, low_t=30, high_t=100, low_p=5, high_p=20, title='Title'):
    """
    Args:
        signal: numpy array with discrete-time signal values
        ax: matplotlib axes
        low_t: int: lower bound for t in DTFT
        high_t: int: upper bound for t in DTFT
        low_p: float: lower bound for p in DTFT
        high_p: float: upper bound for p in DTFT
        title: str: plot title
    Returns:
        plot in axes
    """
    x, y, snr, peak = compute_spectrum(signal, norm=norm, low_t=low_t, high_t=high_t, low_p=low_p, high_p=high_p)
    plt.rc('grid', linestyle='dashed', color='black', alpha=0.1)
    ax.grid(True)
    ax.plot(x, y)
    ax.set_title(f'{title} SNR={snr:.2f}: Peak={peak:.2f}')
    ax.set_xticks(np.arange(low_p, high_p))
    ax.set_ylim(0, 10)
    ax.set_ylabel('Power (a.u.)')
    ax.set_xlabel('Period (bp)')


if __name__ == '__main__':
    build_tables()
    pass
