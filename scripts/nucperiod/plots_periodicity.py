import os
from rpy2 import robjects
from matplotlib import ticker
import matplotlib.pyplot as plt
from rpy2.robjects import pandas2ri

import nucperiod.positions as po
from nucperiod.subplots_periodicity import wave_painting_zoomin, \
    rectangles_drawing, config_params_full, wave_painting_zoomout


pandas2ri.activate()

# for y-axis formating
FORMATTER = ticker.ScalarFormatter(useMathText=True)
FORMATTER.set_scientific(True)
FORMATTER.set_powerlimits((-3, 3))


def get_dict(f):
    f['Position'] = range(-73,74)
    d_pos2 = f.set_index('Position').to_dict()['rel_inc']
    win = 147

    half_win = int((win - 1) / 2)
    xvals2 = []
    yvals2 = []

    for ix, i in enumerate(range(-half_win, half_win + 1)):
        yvals2.append(d_pos2.get(i, 0))
        xvals2.append(i)

    return xvals2, yvals2


spline = robjects.r["smooth.spline"]


def plot(time, x, y, snr, peak, pval, output):
    config_params_full()
    xvals, yvals = get_dict(time)

    val = spline(xvals, yvals)

    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7.5 / 2.13, 7.5 / 1.83))
    plt.subplots_adjust(hspace=0.001)

    for xd in po.DYAD_X:
        ax[0].axvline(xd, color='black', linestyle='--', lw=0.2, alpha=0.5)

    ax[0].plot(xvals, yvals, linewidth=0.3, color='lightblue', alpha=0.7)

    wave_painting_zoomin(
        list(val[0]), list(val[1]), ax[0], '#515e19ff', '#f3bb00', '#AFC6E9', 1.2
    )
    rectangles_drawing(ax[0], 58, '#f3bb00', '#515e19ff')
    ax[0].set_xlim(-75, 75)

    ax[1].plot(x, y, color='#31a354', lw=1.2)
    ax[1].set_xlim(4, 50)

    ax[0].set_title('n = 1, SNR = {}\nMP = {}, p-value = {}, phase = -1'.format(
        round(snr, 2), round(peak, 2), round(pval, 2)), fontsize=7.5)
    ax[1].set_title('SNR = {}, MP = {}, phase = -1'.format(
        round(snr, 2), round(peak, 2)), fontsize=7.5
    )
    ax[0].set_ylabel('Relative increase damage')
    ax[0].set_xlabel('Distance from dyad (bp)')
    ax[1].set_ylabel('Power')
    ax[1].set_xlabel('Period (bp)')

    ax[1].spines['right'].set_visible(False)
    ax[1].spines['top'].set_visible(False)
    ax[0].spines['right'].set_visible(False)
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['left'].set_visible(False)

    plt.tight_layout()
    out_file = os.path.join(output, 'nucleosome_periodicity.pdf')
    plt.savefig(out_file)
    plt.close()



def plot_zoomout(time, snr, peak, x, y, signal, output):
    config_params_full()

    xvals = time['Position'].tolist()
    yvals = time['rel_inc'].tolist()

    val = spline(xvals, yvals, spar=0.55)
    
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7.5 / 2.13, 7.5 / 1.83))
    plt.subplots_adjust(hspace=0.001)

    ax[0].plot(xvals, yvals, linewidth=0.3, color='lightblue', alpha=0.7)

    # wave_painting_zoomout(list(val[0]), list(val[1]), ax[0], [353, 500, 647])
    ax[0].plot(list(val[0]), list(val[1]), linewidth=1, color='#08519c', alpha=0.9)

    ax[0].set_xlim(-520, 520)

    ax[1].plot(x, y, color='#31a354', lw=1.2)
    ax[1].set_xlim(100, 250)

    ax[0].set_title('n = 1, SNR = {}\nMP = {}'.format(
        round(snr, 2), round(peak, 2)), fontsize=7.5)
    ax[1].set_title('SNR = {}, MP = {}'.format(
        round(snr, 2), round(peak, 2)), fontsize=7.5
    )
    ax[0].set_ylabel('Relative increase damage')
    ax[0].set_xlabel('Distance from dyad (bp)')
    ax[1].set_ylabel('Power')
    ax[1].set_xlabel('Period (bp)')

    ax[1].spines['right'].set_visible(False)
    ax[1].spines['top'].set_visible(False)
    ax[0].spines['right'].set_visible(False)
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['left'].set_visible(False)

    plt.tight_layout()
    out_file = os.path.join(output, 'zoomout_periodicity.pdf')
    plt.savefig(out_file)
    plt.close()