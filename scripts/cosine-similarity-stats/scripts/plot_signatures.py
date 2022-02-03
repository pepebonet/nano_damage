import itertools
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def mut_key_generator():
    """

    Returns:
        Generates all possible lex sortable mutation keys
        1st component: substitution;
        2nd component: flanks

    """
    subs = ['CA', 'CG', 'CT', 'TA', 'TC', 'TG']
    for s in sorted(subs):
        for c in sorted(itertools.product(set('ACGT'), repeat=2)):
            yield tuple([s, ''.join(c)])


def minor_labels():
    major_labels = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
    flanks = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT',
              'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    minor_labels = []
    for subs in major_labels:
        for flank in flanks:
            minor_labels.append(flank[0] + subs[0] + flank[1])
    return minor_labels


def plot_signature(profile, title=None, ax=None, ymax=0.3):
    """
    Args:
        profile: signature-like object in lexicographic order
        title: string
        ymax: float

    Returns:
        produces the signature bar plot
    """

    total = sum(profile.values())
    if abs(total - 1) > 0.01:
        profile = defaultdict(int, {k: v / total for k, v in profile.items()})
    sns.set(font_scale=1.5)
    sns.set_style('white')
    vector = np.array([profile[k] for k in sorted(mut_key_generator())])

    # bar plot
    barlist = ax.bar(range(96), vector)
    color_list = ['#72bcd4', 'k', 'r', '#7e7e7e', 'g', '#e6add8']
    for category in range(6):
        for i in range(16):
            barlist[category * 16 + i].set_color(color_list[category])
    ax.set_xlim([-0.5, 96])
    ax.set_ylim([0, ymax])

    # ax.set_ylabel('subs rel freq')
    major_labels = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
    major_ticks = np.arange(8, 8 + 16 * 5 + 1, 16)
    minor_ticks = np.arange(0.2, 96.2, 1)
    ax.tick_params(length=0, which='major', pad=30, labelsize=14)
    ax.tick_params(length=0, which='minor', pad=5, labelsize=8)
    ax.set_xticks(major_ticks, minor=False)
    ax.set_xticklabels(major_labels, minor=False)
    ax.set_xticks(minor_ticks, minor=True)
    ax.set_xticklabels(minor_labels(), minor=True, rotation=90)
    ax.set_title(title, fontsize=24)


if __name__ == '__main__':

    fig, ax = plt.subplots()
    uniform = {k: 1 / 96 for k in mut_key_generator()}
    plot_signature(uniform, ax=ax)
    plt.show()
