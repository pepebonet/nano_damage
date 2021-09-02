import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import ticker
from rpy2.robjects import pandas2ri
from rpy2 import robjects
pandas2ri.activate()

from nucperiod import plot as nucperiod_plot, positions, spectral

# for y-axis formating
FORMATTER = ticker.ScalarFormatter(useMathText=True)
FORMATTER.set_scientific(True)
FORMATTER.set_powerlimits((-3, 3))


def get_dict(f):
    df2 = pd.read_csv(f, sep='\t', names=['id', 'pos', 'sig'])
    d_pos2 = df2['pos'].value_counts().to_dict()
    win = 117

    half_win = int((win - 1) / 2)
    xvals2 = []
    yvals2 = []

    for ix, i in enumerate(range(-half_win, half_win + 1)):
        yvals2.append(d_pos2.get(i, 0))
        xvals2.append(i)

    return xvals2, yvals2


spline = robjects.r["smooth.spline"]

def plot(time0, time1):
    nucperiod_plot.config_params_full()
    xvals, yvals = get_dict(time0)
    xvals2, yvals2 = get_dict(time1)
    div = (np.array(yvals) - np.array(yvals2)) / np.array(yvals)
    val = spline(xvals2, div)

    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(4.88 / 2.13, 7.5 / 1.83))
    plt.subplots_adjust(hspace=0.001)

    for xd in positions.DYAD_X_SMALL:
        ax[0].axvline(xd, color='black', linestyle='--', lw=0.2, alpha=0.5)
        ax[1].axvline(xd, color='black', linestyle='--', lw=0.2, alpha=0.5)

    ax[0].set_title('Repair', fontsize = 7.5)
    ax[0].plot(xvals, yvals, c='#1f78b4ff', linewidth=0.8,)
    ax[0].plot(xvals2, yvals2, c='#1f78b4ff', linewidth=0.8, alpha=0.5)

    ax[0].yaxis.set_major_formatter(FORMATTER)

    ax[0].spines['right'].set_visible(False)
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['bottom'].set_visible(False)
    ax[0].get_xaxis().set_visible(False)
    ax[0].set_xlim(-58, 59)
    ax[0].set_ylabel('Damage reads')

    ax[1].plot(xvals, div, linewidth=0.3, color='lightblue', alpha=0.7)

    nucperiod_plot.wave_painting_zoomin(list(val[0]), list(val[1]), ax[1], '#515e19ff', '#f3bb00', '#AFC6E9', 1.2)
    nucperiod_plot.rectangles_drawing(ax[1], 58, '#f3bb00', '#515e19ff')

    x, y, snr, peak = spectral.compute(list(val[1]), center=10.3,
                                       low_t=0, high_t=115, low_p=5, high_p=20, norm=True)

    ax[2].plot(x, y, color='#31a354', lw=1.2)
    ax[2].set_xlim(4, 20)

    ax[1].set_title('n = 1, SNR = {}\nMP = {}, p-value = n1n, phase = -1'.format(round(snr, 2), round(peak, 2),
                                                                                 ), fontsize=7.5)
    ax[2].set_title('SNR = {}, MP = {}, phase = -1'.format(round(snr, 2), round(peak, 2)), fontsize=7.5)
    ax[1].set_ylabel('Repair rate')
    ax[1].set_xlabel('Distance from dyad (bp)')
    ax[2].set_ylabel('Power')
    ax[2].set_xlabel('Period (bp)')
    ax[2].spines['right'].set_visible(False)
    ax[2].spines['top'].set_visible(False)
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['top'].set_visible(False)
    ax[1].spines['left'].set_visible(False)

    plt.tight_layout()