
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from nucperiod.positions import MINORS_IN, MINORS_OUT, DYAD_X_SMALL


def config_params(font_size=7):
    mpl.rcParams.update(mpl.rcParamsDefault)
    plt.rcParams['font.sans-serif'] = ['arial']
    plt.rcParams['font.size'] = font_size
    plt.rcParams['font.family'] = ['sans-serif']
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.cal'] = 'arial'
    plt.rcParams['mathtext.rm'] = 'arial'


def config_params_full(font_size=7):
    config_params(font_size)
    plt.rcParams['axes.linewidth'] = 0.4
    plt.rcParams['xtick.major.width'] = 0.4
    plt.rcParams['xtick.minor.width'] = 0.4
    plt.rcParams['ytick.major.width'] = 0.4
    plt.rcParams['ytick.minor.width'] = 0.4


def wave_painting_zoomin(xvals, yvals, axs, color_in, color_out, color_other, line_width):

    colors = []
    number_to_split = 10

    all_xs = []
    all_ys = []
    for ix in range(len(xvals) - 1):
        ypointA = yvals[ix]
        ypointB = yvals[ix + 1]
        xpointA = xvals[ix]
        xpointB = xvals[ix + 1]
        increase = (ypointB - ypointA) / number_to_split
        xvals_small = np.arange(xpointA, xpointB, 1 / number_to_split)
        yvals_small = np.arange(ypointA, ypointB, increase)
        all_xs.extend(xvals_small)
        all_ys.extend(yvals_small)

        for x in xvals_small:

            matched_color = 0
            for (s, e) in MINORS_OUT:
                if s <= x < e:
                    colors.append(color_out)
                    matched_color = 1
                    break
            if matched_color == 0:
                for (s, e) in MINORS_IN:
                    if s <= x < e:
                        colors.append(color_in)
                        matched_color = 1
                        break
            if matched_color == 0:
                colors.append(color_other)

    to_plot = []
    for ix, d in enumerate(zip(all_xs, all_ys, colors)):
        to_plot.append(d)

    for i in range(len(to_plot) - 1):
        x, v, c = to_plot[i]
        x2, v2, c2 = to_plot[i + 1]
        axs.plot((x, x2), (v, v2), c=c, lw=line_width, solid_capstyle='round', solid_joinstyle='round')


def wave_painting_zoomout(xvals, yvals, axs, good_local_maxima):

    nucleosome_zone = []
    linker_zone = []

    for ix in range(len(good_local_maxima) - 1):

        maxima = good_local_maxima[ix]

        # this covers half of the nucleosome to the right
        start_linker = maxima + 73

        # this gets the putative length of the linker
        end_linker = (good_local_maxima[ix + 1] - 73)

        for p in range(int(start_linker), int(end_linker) + 1):
            linker_zone.append(p)

        for p in range(int(maxima) - 73, int(maxima) + 74):
            nucleosome_zone.append(p)

    for p in range(int(good_local_maxima[ix + 1]) - 73, int(good_local_maxima[ix + 1]) + 74):  # add last nucleosome
        nucleosome_zone.append(p)

    colors = []

    for ix, i in enumerate(xvals):

        if i in linker_zone:
            colors.append('#7570b3')
        elif i in nucleosome_zone:
            colors.append('#d95f02')
        else:
            colors.append('#AFC6E9')

    # array to plot
    to_plot = []
    for ix, d in enumerate(zip(xvals, yvals, colors)):
        to_plot.append(d)

    for i in range(len(to_plot) - 1):
        x, v, c = to_plot[i]
        x2, v2, c2 = to_plot[i + 1]
        axs.plot((x, x2), (v, v2), c=c, lw=1.2, solid_capstyle='round')


def rectangles_drawing(axs, half_win, color1, color2):

    ax2 = axs.twinx()
    ax2.set_xticks(axs.get_xticks())

    for xv in DYAD_X_SMALL:

        ax2.add_patch(plt.Rectangle((xv - 2.575, -0.105), 5.15, 0.06, color=color1,
                                    clip_on=False, linewidth=0))

        if xv != -10.3:
            ax2.add_patch(plt.Rectangle((xv - 2.575 + 5.15, -0.105), 5.15, 0.06, color=color2,
                                        clip_on=False, linewidth=0))

    ax2.add_patch(plt.Rectangle((-61.8 - 2.575 + 5.15, -0.105), 5.15, 0.06, color=color2,
                                clip_on=False, linewidth=0, ))

    axs.spines["bottom"].set_position(("axes", -0.1093))

    ax2.set_xlim(-half_win - 1, half_win + 1)
    axs.set_xlim(-half_win - 1, half_win + 1)

    ax2.get_yaxis().set_ticks([])

    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)
