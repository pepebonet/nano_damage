#!/usr/bin/env python3
import os
import seaborn as sns
import matplotlib.pyplot as plt

import utils as ut
from scipy.signal import savgol_filter


# ------------------------------------------------------------------------------
# CHROMATIN ANALYSIS
# ------------------------------------------------------------------------------

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



def obs_exp_simple_plot(data, output, pval, label=''):
    fig, ax = plt.subplots(figsize=(5, 5))
    out_file = os.path.join(output, '{}.png'.format(label))

    sns.barplot(x=list(data.keys()), y=list(data.values()), palette=['#08519c', '#a6cee3'])
    
    ax.set_title('p-value: {}'.format(pval))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.tight_layout()
    plt.savefig(out_file)
    plt.close()


def obs_exp_multi_plots(data, output, label=''):
    fig, ax = plt.subplots(nrows=1, ncols=4, figsize=(15, 5)) 
 
    out_file = os.path.join(output, '{}.png'.format(label))  
    i = 0
    for k, v in data.items():
        sns.barplot(x=list(v.keys()), y=list(v.values()), ax=ax[i])
        ax[i].title.set_text(k)
        ax[i].spines['top'].set_visible(False)
        ax[i].spines['right'].set_visible(False)
        ax[i].set_xlabel('')
        ax[i].tick_params(labelrotation=20)
        i += 1 
    
    plt.tight_layout()

    plt.savefig(out_file)
    plt.close()


def plot_cosine(df, output):
    fig, ax = plt.subplots(figsize=(7, 4))
    
    df_triplet = df[df['Context'] == 'Triplet']
    df_penta = df[df['Context'] == 'Pentamer']
    x = list(range(len(df_triplet['DNA structure'])))
    import pdb;pdb.set_trace()
    plt.errorbar(
        [a + 0.1 for a in x], df_triplet['Cosine-similarity'], 
        ls='none', marker='o', mfc='#08519c', mec='black', ms=10, mew=1, 
        ecolor='#08519c', capsize=2.5, elinewidth=0.7, capthick=0.7
    )

    plt.errorbar(
        [a - 0.1 for a in x], df_penta['Cosine-similarity'], 
        ls='none', marker='o', mfc='#de971d', mec='black', ms=10, mew=1, 
        ecolor='#de971d', capsize=2.5, elinewidth=0.7, capthick=0.7
    )

    custom_lines = []
    for el in [('Triplet', '#08519c'), ('Pentamer', '#de971d')]:
        custom_lines.append(
            plt.plot([],[], marker="o", ms=8, ls="", mec='black', 
            mew=0, color=el[1], label=el[0])[0] 
        )

    plt.ylim(round(df['Cosine-similarity'].min() - 0.04, 2), 1.03)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_xlabel("Genomic regions", fontsize=12)
    ax.set_ylabel("Cosine similarity", fontsize=12)
    ax.set_xticklabels([''] + df_penta['DNA structure'].tolist())
    from matplotlib.offsetbox import AnchoredText
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    title_text = 'Similarity on Genomic Regions'

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top", size="7.5%", pad=0)
    cax.get_xaxis().set_visible(False)
    cax.get_yaxis().set_visible(False)

    cax.set_facecolor('#c7cbd4')
    at = AnchoredText(title_text, loc=10, frameon=False,
            prop=dict(backgroundcolor='#c7cbd4',
                    size=8, color='black'))
    cax.add_artist(at)

    ax.legend(
        bbox_to_anchor=(0., 0.2, 1.0, .102),
        handles=custom_lines, loc='upper right', 
        facecolor='white', ncol=1, fontsize=10, frameon=False
    )

    fig.tight_layout()
    out_file = os.path.join(output, 'cosine_similarity.pdf')
    plt.savefig(out_file)
    plt.close()


def significant_bases_distance(distance, output, title=''):
    fig, ax = plt.subplots(figsize=(5, 5))

    sns.barplot(x=list(distance.keys()), y=list(distance.values()))
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # plt.title('Cisplatin', fontsize=16)
    fig.tight_layout()
    out_file = os.path.join(output, 'distribution.png')
    plt.savefig(out_file)
    plt.close()


# ------------------------------------------------------------------------------
# NUCLEOSOME ANALYSIS
# ------------------------------------------------------------------------------

def plot_per_base_enrichment(df, outdir, label):
    x = df['POSITION'].tolist()

    if label == 'smooth':
        yhat = df['NORM_2_smooth'].tolist()
        yrandom = df['Random Model_smooth'].tolist()
    else:
        yhat = df['NORM_2'].tolist()
        yrandom = df['Random Model'].tolist()
    
    # import pdb;pdb.set_trace()
    # Smoothing using: window size 51, polynomial order 3
    yhat = savgol_filter(yhat, 9, 3)
    yrandom = savgol_filter(yrandom, 9, 3)

    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(20, 5),
                            gridspec_kw={'height_ratios': [30, 1]})

    #1st axis. Lower color plot

    aa = [el * 146 / 23 for el in list(range(23))]
    for ix, j in enumerate(aa):
        if ix in [10, 11, 12]:
                color = 'white'
        else:
            if ix % 2 == 0:
                color = '#006400'
            else:
                color = '#FFB90F'
        
        axs[1].barh(1, aa[1], left=j, color=color)
        

    axs[1].set_xlim(-1, 148)
    axs[1].spines['top'].set_visible(False)
    axs[1].spines['bottom'].set_visible(False)
    axs[1].spines['left'].set_visible(False)
    axs[1].spines['right'].set_visible(False)

    axs[1].get_yaxis().set_visible(False)
    axs[1].get_xaxis().set_visible(False)

    #2nd axis. Plot
    plt.sca(axs[0])
    axs[0].plot(x, yhat, linewidth=4)
    axs[0].plot(x, yrandom, linewidth=2, color='grey', ls='--')
    # axs.set_xticks(order_plot)
    axs[0].set_ylabel('Relative Probability',fontsize=24)
    axs[0].spines['top'].set_visible(False)
    axs[0].spines['right'].set_visible(False)
    axs[0].set_xlim(-1, 148)
    plt.setp([axs[0].get_xticklines(), axs[0].get_yticklines()], color='grey')

    axs[0].xaxis.set_ticks_position('none')
    for axis in ['top', 'bottom', 'left', 'right']:
        axs[0].spines[axis].set_linewidth(0.2)

    axs[0].xaxis.set_tick_params(pad=0.5)
    axs[0].yaxis.set_tick_params(pad=0.5, width=0.5)

    plt.xticks(fontsize=20, rotation=90)
    plt.yticks(fontsize=16)
    plt.tick_params(axis='both', which='both', bottom=False, left = False)
    axs[1].set_xlabel('Position', fontsize=24)

    fig.tight_layout()

    plt.savefig(outdir)
    plt.close()


def plot_norm_nucleosomes(df, outdir):
    # import pdb;pdb.set_trace()
    x = df['Position'].values.flatten()

    yhat = df['Relative Increase'].values.flatten()
    yhat = savgol_filter(yhat, 11, 3)
    
    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(20, 5),
                            gridspec_kw={'height_ratios': [30, 1]})

    #1st axis. Lower color plot

    aa = [el * 146 / 23 for el in list(range(23))]
    for ix, j in enumerate(aa):
        if ix in [10, 11, 12]:
                color = 'white'
        else:
            if ix % 2 == 0:
                color = '#006400'
            else:
                color = '#FFB90F'
        
        axs[1].barh(1, aa[1], left=j, color=color)
        

    axs[1].set_xlim(-1, 148)
    axs[1].spines['top'].set_visible(False)
    axs[1].spines['bottom'].set_visible(False)
    axs[1].spines['left'].set_visible(False)
    axs[1].spines['right'].set_visible(False)

    axs[1].get_yaxis().set_visible(False)
    axs[1].get_xaxis().set_visible(False)

    #2nd axis. Plot
    plt.sca(axs[0])
    axs[0].plot(x, yhat, linewidth=4)

    # axs.set_xticks(order_plot)
    axs[0].set_ylabel('Relative Increase',fontsize=24)
    axs[0].spines['top'].set_visible(False)
    axs[0].spines['right'].set_visible(False)
    axs[0].set_xlim(-1, 148)
    axs[0].set_ylim(-0.3, 0.3)
    plt.setp([axs[0].get_xticklines(), axs[0].get_yticklines()], color='grey')

    axs[0].xaxis.set_ticks_position('none')
    for axis in ['top', 'bottom', 'left', 'right']:
        axs[0].spines[axis].set_linewidth(0.2)

    axs[0].xaxis.set_tick_params(pad=0.5)
    axs[0].yaxis.set_tick_params(pad=0.5, width=0.5)

    plt.xticks(fontsize=20, rotation=90)
    plt.yticks(fontsize=16)
    plt.tick_params(axis='both', which='both', bottom=False, left = False)
    axs[1].set_xlabel('Position', fontsize=24)

    fig.tight_layout()

    plt.savefig(outdir)
    plt.close()


def plot_zoom_out_nucleosomes(df, outdir):
    # import pdb;pdb.set_trace()
    x = df['Position'].values.flatten()

    yhat = df['Relative Increase'].values.flatten()
    yhat = savgol_filter(yhat, 21, 3)
    
    fig, ax = plt.subplots(figsize=(6, 5))

    ax.plot(x, yhat, linewidth=2, alpha=0.9)

    ax.set_ylabel('Relative Increase',fontsize=14)
    ax.set_xlabel('Position',fontsize=14)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    ax.vlines(0, ymin=-0.2, ymax=0.2, linestyle='dashed', color='grey', linewidth=1)
    ax.vlines(73, ymin=-0.2, ymax=0.2, linestyle='dashed', color='grey', linewidth=1)
    ax.vlines(-73, ymin=-0.2, ymax=0.2, linestyle='dashed', color='grey', linewidth=1)

    plt.setp([ax.get_xticklines(), ax.get_yticklines()], color='grey')

    ax.xaxis.set_ticks_position('none')
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(0.2)

    ax.xaxis.set_tick_params(pad=0.5)
    ax.yaxis.set_tick_params(pad=0.5, width=0.5)

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tick_params(axis='both', which='both', bottom=False, left = False)

    fig.tight_layout()

    plt.savefig(outdir)
    plt.close()



def plot_damage_nuc_linker(df, output, nuc_signal):
    fig, ax = plt.subplots(figsize=(5, 5))

    plt.plot(df['Position'].values, df['smooth_damage'].values, linewidth=2)
    plt.plot(df['Position'].values, df['smooth_random'].values,
        linewidth=1, color='black', ls='--')

    plt.vlines(200, ymin=0.002, ymax=0.003, linestyle='dashed', color='grey', linewidth=0.5)
    plt.vlines(273, ymin=0.002, ymax=0.003, linestyle='dashed', color='grey', linewidth=0.5)
    plt.vlines(127, ymin=0.002, ymax=0.003, linestyle='dashed', color='grey', linewidth=0.5)
    
    ax.set_ylabel('Relative Probability')
    
    ax2 = ax.twinx()

    color = 'red'
    ax2.set_ylabel('Nucleosome Probability', color=color)
    ax2.plot(nuc_signal['index'].values, nuc_signal[0].values, color='red', alpha=0.7)
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.set_ylim(0, 2)

    fig.tight_layout()
    out_file = os.path.join(output, 'damage_nuc_linker.pdf')
    plt.savefig(out_file)
    plt.close()



def plot_snrs(snrs, snr_obs, output):
    fig, ax = plt.subplots(figsize=(5, 5))

    sns.displot(snrs, kind="kde")

    plt.vlines(
        snr_obs, ymin=0, ymax=0.2, linestyle='dashed', color='red', linewidth=0.5
    )

    fig.tight_layout()
    out_file = os.path.join(output, 'snrs_distribution.pdf')
    plt.savefig(out_file)
    plt.close()


def plot_peaks(peaks, peak_obs, output):
    fig, ax = plt.subplots(figsize=(5, 5))

    sns.displot(peaks, kind="kde")

    plt.vlines(
        peak_obs, ymin=0, ymax=0.2, linestyle='dashed', color='red', linewidth=0.5
    )

    fig.tight_layout()
    out_file = os.path.join(output, 'peaks_distribution.pdf')
    plt.savefig(out_file)
    plt.close()