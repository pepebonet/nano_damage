#!/usr/bin/envs python3
import matplotlib.pyplot as plt

import utils as ut


#Obtain context plots for triplets and pentamers
def obtain_plots(data, outdir, context, chunk):
    order_plot, y = ut.order_plots(data, context)

    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(20, 5),
                            gridspec_kw={'height_ratios': [1, 30]})

    #1st axis. Upper color plot
    colors = []
    colors_base = ['#1ebff0', '#050708', '#e62725', '#cbcacb']
    bot = 0.5
    for ix, c in enumerate(ut.chunks(y, chunk)):
        colors.extend([colors_base[ix] for s in c])
        axs[0].barh(1, chunk, left=bot, color=colors_base[ix])
        bot += chunk


    axs[0].set_xlim(0, 4 * chunk)
    axs[0].spines['top'].set_visible(False)
    axs[0].spines['bottom'].set_visible(False)
    axs[0].spines['left'].set_visible(False)
    axs[0].spines['right'].set_visible(False)

    axs[0].get_yaxis().set_visible(False)
    axs[0].get_xaxis().set_visible(False)

    #2nd axis. Bar plot
    plt.sca(axs[1])

    axs[1].bar(order_plot, y, color=colors)
    axs[1].set_ylabel('Relative Probability',fontsize=24)
    
    if context == 'triplet':
        axs[1].set_xticks(order_plot)
        axs[1].set_xlabel('Triplet Context', fontsize=24)
    else: 
        axs[1].tick_params(labelbottom=False)
        axs[1].set_xlabel('Pentamer Context', fontsize=24)

    axs[1].spines['top'].set_visible(False)
    axs[1].spines['right'].set_visible(False)

    plt.setp([axs[1].get_xticklines(), axs[1].get_yticklines()], color='grey')

    axs[1].xaxis.set_ticks_position('none')
    for axis in ['top', 'bottom', 'left', 'right']:
        axs[1].spines[axis].set_linewidth(0.2)

    axs[1].xaxis.set_tick_params(pad=0.5)
    axs[1].yaxis.set_tick_params(pad=0.5, width=0.5)

    plt.xticks(fontsize=20, rotation=90)
    plt.yticks(fontsize=16)
    plt.tick_params(axis='both', which='both', bottom=False, left = False)

    fig.tight_layout()

    plt.savefig(outdir)
    plt.close()