#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib as mpl
from scipy.special import xlogy
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
                     label=m,
                     marker='o',
                     ms=ms[i])

        for i, m in enumerate(np.unique(method)):
            dens = density[(method == m) & (vector_size == vs)]
            sizes = bits_per_entry_expected[(method == m) & (vector_size == vs)]
            idx = dens.argsort()
            plt.plot(dens[idx], sizes[idx],
                     label=m + ' predicted',
                     marker='*',
                     ms=ms[i])

        plt.title('{}: {:.1e}'.format(name, vs))
        plt.xlabel('Density')
        plt.ylabel('Bits per entry')

        plt.grid(True)
        plt.legend(fontsize=18)
        plt.tight_layout()
        plt.savefig('vectors_{}_{}.pdf'.format('_'.join(name.split(' ')).lower(), vs), fmt='pdf')
        plt.show()

plot_feature(bits_per_entry, 'Serialized Size')
plot_feature(RAM_per_entry, 'RAM')
