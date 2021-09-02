import cmath

import numpy as np
from rpy2 import robjects
from rpy2.robjects import pandas2ri

from nucperiod.non_linear import create_quadratic_model, NonLinearFitting

pandas2ri.activate()


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
            return (1 / (smax - smin)) * abs(
                sum(self.signal[s] * cmath.exp(-1j * (s-smin+1) * 2 * np.pi / period) for s in range(smin, smax))) ** 2

        def phase(period):
            return cmath.phase(sum(self.signal[s] * cmath.exp(-1j * (s-smin+1) * 2 * np.pi / period) for s in range(smin, smax)))

        return power, phase


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
    return power, phase


def spectrum(signal, x, norm=True, **kwargs):
    """
    Args:
        signal: numpy array with discrete-time signal values
        x: array-like: grid of periods
        norm: whether to provide normalized spectrogram
        kwargs: parameters to run compute_fourier
    Returns:
        tuple:
            y-values: DTFT of x-values.
    """
    signal = mean_smooth(signal, 3)  # 3-bp average smoothing reduces effects of 3-bp periodic signal
    power, _ = compute_fourier(signal, **kwargs)
    norm_factor = 1
    if norm:
        norm_factor = (kwargs['high_p'] - kwargs['low_p']) * (
                    1 / sum(power(p) for p in range(kwargs['low_p'], kwargs['high_p'])))
    y = list(map(lambda a: norm_factor * power(a), x))
    return y


def closest_index(x, value):
    return np.abs(x-value).argmin()


def signal_to_noise(x, y, center=None):
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
     See Lehmann, Machne and Herzel NAR 2014.

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


def compute(signal, norm=True, center=None, **kwargs):
    """
    Args:
        signal: numpy array with discrete-time signal values
        norm: whether to provide normalized spectrogram
        kwargs: keyword arguments compatible with the function compute_fourier
    Returns:
        tuple:
            x-values
            y-values: DTFT of x-values.
            snr: float: signal-to-noise ratio
            peak: max peak relative to the center
    """
    signal = mean_smooth(signal, 3)  # 3-bp average smoothing reduces effects of 3-bp periodic signal
    power, _ = compute_fourier(signal, **kwargs)
    norm_factor = 1
    if norm:
        norm_factor = (kwargs['high_p'] - kwargs['low_p']) * (
                    1 / sum(power(p) for p in range(kwargs['low_p'], kwargs['high_p'])))
    x = np.linspace(kwargs['low_p'], kwargs['high_p'], 100)
    y = list(map(lambda a: norm_factor * power(a), x))
    snr, peak = signal_to_noise(x, y, center=center)

    return x, y, snr, peak


def synthetic_signal_out(N, P, bounds=None):

    if bounds is None:
        bounds = [-61.8, 61.8]

    x = np.linspace(bounds[0], bounds[1], N)
    return np.cos(2 * np.pi * x / P)


def assess_phase(xvals, signal, N, P, bounds, zoomin):
    """

    :param xvals:
    :param signal:
    :param N:
    :param P:
    :param bounds:
    :return:
    """

    # we will do the crosscorrelation with the smoothing data
    spline = robjects.r["smooth.spline"]
    val = spline(xvals, signal, df=60)

    # create synthetic signals
    synth = synthetic_signal_out(N, P, bounds=bounds)

    # Phase analysis

    _, phase = compute_fourier(synth, low_t=0, high_t=len(synth))
    phase_synt = phase(P)

    _, phase = compute_fourier(list(val[1]), low_t=0, high_t=len(signal))
    phase_real = phase(P)

    # phase differences
    dif_val = np.abs(phase_synt - phase_real)

    # we should make the value less than 2pi
    while dif_val >= np.pi * 2:
        dif_val -= np.pi * 2

    # if this is smaller than pi/2, this means it is out of phase compared to melanoma-like
    if np.abs(np.pi - dif_val) < (np.pi / 2):
        phase = -1
    # this means in phase
    else:
        phase = 1
    if zoomin:
        phase = phase*-1

    return phase
