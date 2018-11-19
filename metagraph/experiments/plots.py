#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
from scipy.special import comb

mpl.rcParams['figure.dpi'] = 150
rcParams.update({'errorbar.capsize' : 2})
rcParams['text.antialiased'] = True


def remove_padding(plt):
    plt.tight_layout(pad=0, w_pad=0, h_pad=0)


def remove_legend_bars(plt, ax, fontsize=18):
    handles, labels = ax.get_legend_handles_labels()
    handles = [h[0] for h in handles]
    plt.legend(handles, labels, fontsize=16)


plot_ratio = (1 + 5 ** 0.5) / 2

matrix = np.vstack([ np.array(a.rstrip().split(" "))
                            for a in open(sys.argv[1]).readlines() ])
matrix_label = matrix[:, 0]
matrix = matrix[:, 1:].astype(float)

#bytes_per_element = matrix[:,5] / matrix[:,2]
bytes_per_element = matrix[:, 5] / (matrix[:, 0] * matrix[:, 1])

labels = np.unique(matrix_label)
label_dict = {
    "rbfish": "Rainbowfish",
    "brwt": "BRWT",
    "binrel": "BinRelWT",
    "column": "Column compressed",
}

masks = [ (label_dict[label], matrix_label == label) for label in labels ]
#x = comb(matrix[:,1].astype(int), (matrix[:,1] * matrix[:,4]).astype(int)) / matrix[:,0]
x = matrix[:,4]

f, ax = plt.subplots(figsize=(6 * plot_ratio, 6))
for label, mask in masks:
    argsort = np.argsort(x[mask])
    curx = x[mask][argsort]
    y = bytes_per_element[mask][argsort]
    ymean = np.array([np.mean(y[a == curx]) for a in np.unique(curx)])
    ystd = np.array([np.std(y[a == curx]) for a in np.unique(curx)])
    plt.errorbar(np.unique(curx),
                 ymean,
                 yerr=ystd,
                 #linestyle="",
                 marker="o",
                 label=label)


plt.legend()
plt.xlabel("Mean row density")
plt.ylabel("Bytes per matrix element")
ax.grid(linestyle=":", which="both")
#remove_padding(plt)
#plt.xscale("log")
#plt.yscale("log")

plt.show()
