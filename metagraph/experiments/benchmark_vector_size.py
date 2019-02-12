#!/usr/bin/env python3

import sys
import numpy as np
from matplotlib import pyplot as plt

plot_ratio = (1 + np.sqrt(5)) / 2
markers = ['^', 'v', '<', '>', 'o', 's', 'P', '*']


if len(sys.argv) != 2:
    print("Usage:\n{} plot <results.txt>".format(sys.argv[0]))
    print(
"""\n\n#First, run this script to generate vectors and measure their size

for type in small stat sd rrr63 rrr127; do
    for d in $(seq 0 0.02 1); do
        ./run_experiments vectors $type -d $d -l 100000000 >> results.txt
    done
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

plt.rcParams.update({ 'font.size': 22 })

figsize = 8

for vs in np.unique(vector_size):
    plt.figure(figsize=(figsize * 1.2, figsize))

    ms = np.linspace(20, 4, len(np.unique(method)))

    for i, m in enumerate(np.unique(method)):
        dens = density[(method == m) & (vector_size == vs)]
        sizes = bits_per_entry[(method == m) & (vector_size == vs)]
        idx = dens.argsort()
        plt.plot(dens[idx], sizes[idx],
                 label=m,
                 marker=markers[i],
                 ms=ms[i])

    plt.title('Vector size: {:.1e}'.format(vs))
    plt.xlabel('Density')
    plt.ylabel('Bits per entry')

    plt.grid(True)
    plt.legend(fontsize=20)
    plt.tight_layout()
    plt.savefig('vectors_size_{}.pdf'.format(vs), fmt='pdf')
    plt.show()
