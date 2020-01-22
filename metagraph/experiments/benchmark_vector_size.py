#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib as mpl
from scipy.special import xlogy
from scipy.ndimage.filters import gaussian_filter1d
import seaborn as sns
mpl.use('agg')
from matplotlib import pyplot as plt

plot_ratio = (1 + np.sqrt(5)) / 2
markers = ['^', 'v', '<', '>', 'o', 's', 'P', '*']


if len(sys.argv) != 2:
    print("Usage:\n{} <results.txt>".format(sys.argv[0]))
    print(
"""\n\n#First, run this script to generate vectors and measure their size

for type in small smart stat sd rrr63 rrr127; do
    for d in $(seq 0 0.02 1); do
        ./run_experiments vectors $type -d $d -l 100000000 >> results.txt
    done
done

or, in parallel:

for type in {"stat","sd","rrr63","small","smart"}; do
    seq 0 0.05 1 | xargs -n 1 -P 50 -I % ./run_experiments vectors $type -d % -l 2000000 >> results.txt;
done

""")
    exit(1)

results_filename = sys.argv[1]

with open(results_filename, "r") as rfile:
    lines = np.vstack([ np.array(a.rstrip("\n").split("\t")) for a in rfile ])

method = lines[:, 0]
vector_size = lines[:, 1].astype(int)
density = lines[:, 3].astype(float)
bits_per_entry = lines[:, 4].astype(float)
RAM_per_entry = lines[:, 5].astype(float) / vector_size * 8
bits_per_entry_expected = lines[:, 6].astype(float)

plt.rcParams.update({ 'font.size': 22 })

figsize = 8

def plot_feature(feature=bits_per_entry, name="Bits per entry"):
    for vs in np.unique(vector_size):
        plt.figure(figsize=(figsize * 1.2, figsize))

        ms = np.linspace(10, 4, len(np.unique(method)))

        for i, m in enumerate(np.unique(method)):
            dens = density[(method == m) & (vector_size == vs)]
            sizes = feature[(method == m) & (vector_size == vs)]
            idx = dens.argsort()
            plt.plot(dens[idx], sizes[idx],
                     label=m, marker='o', ms=ms[i])

        for i, m in enumerate(np.unique(method)):
            dens = density[(method == m) & (vector_size == vs)]
            sizes = bits_per_entry_expected[(method == m) & (vector_size == vs)]
            idx = dens.argsort()
            plt.plot(dens[idx], sizes[idx],
                     label=m + ' predicted', marker='*', ms=ms[i])

        plt.title('{}: {:.1e}'.format(name, vs))
        plt.xlabel('Density')
        plt.ylabel('Bits per entry')

        plt.grid(True)
        plt.legend(fontsize=18)
        plt.tight_layout()
        plt.savefig('vectors_{}_{}.pdf'.format('_'.join(name.split(' ')).lower(), vs), fmt='pdf')
        # plt.show()

plot_feature(bits_per_entry, 'Serialized Size')
plot_feature(RAM_per_entry, 'RAM')


access_time = lines[:, 7].astype(float)
access_word_time = lines[:, 8].astype(float)
rank_time = lines[:, 9].astype(float)
select_time = lines[:, 10].astype(float)

seq_access_time = lines[:, 11].astype(float)
seq_access_word_time = lines[:, 12].astype(float)
seq_rank_time = lines[:, 13].astype(float)
seq_select_time = lines[:, 14].astype(float)

for vs in np.unique(vector_size):
    fig, ax = plt.subplots(4, 2, figsize=(figsize * 1.2 * 2, figsize * 3))

    ms = np.linspace(10, 4, len(np.unique(method)))

    for i, m in enumerate(np.unique(method)):
        dens = density[(method == m) & (vector_size == vs)]
        idx = dens.argsort()
        sns.lineplot(dens[idx], access_time[(method == m) & (vector_size == vs)][idx],
            ax=ax[0, 0], label=m, marker='o', ms=ms[i])
        sns.lineplot(dens[idx], access_word_time[(method == m) & (vector_size == vs)][idx],
            ax=ax[1, 0], label=m, marker='o', ms=ms[i])
        sns.lineplot(dens[idx], rank_time[(method == m) & (vector_size == vs)][idx],
            ax=ax[2, 0], label=m, marker='o', ms=ms[i])
        sns.lineplot(dens[idx], select_time[(method == m) & (vector_size == vs)][idx],
            ax=ax[3, 0], label=m, marker='o', ms=ms[i])

        sns.lineplot(dens[idx], seq_access_time[(method == m) & (vector_size == vs)][idx],
            ax=ax[0, 1], label=m, marker='o', ms=ms[i])
        sns.lineplot(dens[idx], seq_access_word_time[(method == m) & (vector_size == vs)][idx],
            ax=ax[1, 1], label=m, marker='o', ms=ms[i])
        sns.lineplot(dens[idx], seq_rank_time[(method == m) & (vector_size == vs)][idx],
            ax=ax[2, 1], label=m, marker='o', ms=ms[i])
        sns.lineplot(dens[idx], seq_select_time[(method == m) & (vector_size == vs)][idx],
            ax=ax[3, 1], label=m, marker='o', ms=ms[i])

    ax[0, 0].set_title('Random Access: {:.1e}'.format(vs))
    ax[1, 0].set_title('Random Access 64 Bits Word: {:.1e}'.format(vs))
    ax[2, 0].set_title('Random Rank: {:.1e}'.format(vs))
    ax[3, 0].set_title('Random Select: {:.1e}'.format(vs))

    ax[0, 1].set_title('Sequential Access: {:.1e}'.format(vs))
    ax[1, 1].set_title('Sequential Access 64 Bits Word: {:.1e}'.format(vs))
    ax[2, 1].set_title('Sequential Rank: {:.1e}'.format(vs))
    ax[3, 1].set_title('Sequential Select: {:.1e}'.format(vs))

    for axis_ in ax:
        for axis in axis_:
            axis.set_xlabel('Density')
            axis.set_ylabel('Time, s')
            axis.grid(True)
            axis.legend(loc='upper left', fontsize=18)
            axis.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    plt.tight_layout()
    plt.savefig('vectors_time_{}.pdf'.format(vs), fmt='pdf')
    # plt.show()
