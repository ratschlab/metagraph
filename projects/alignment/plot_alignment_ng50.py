#!/usr/bin/env python3

"""
Jointly plot the scores, NG50, and NG90 for each input alignment JSON file.
JSON files should be in the alignment format of metagraph, vg, etc.

Usage:
./plot_alignment_ng50.py FILE1.json [FILE2.json [...]]
"""

import sys
import json
import numpy as np
from matplotlib import pyplot as plt

def align_stats(file):
    with open(file, "r") as f:
        scores = []
        ng50s = []
        ng90s = []
        unmapped = 0
        for line in f:
            lengths = []
            line_json = json.loads(line)

            if "path" not in line_json:
                unmapped += 1
                continue

            path = line_json["path"]
            if "mapping" not in path:
                unmapped += 1
                continue

            assert("score" in line_json)

            scores.append(line_json["score"])

            length = 0
            for pos in path["mapping"]:
                for edit in pos["edit"]:
                    if "from_length" not in edit \
                      or "to_length" not in edit \
                      or edit["from_length"] != edit["to_length"] \
                      or "sequence" in edit:
                        if length > 0:
                            lengths.append(length)

                        length = 0
                        continue

                    length += edit["from_length"]

            lengths.append(length)
            lengths = np.sort(lengths)[::-1]
            cumsum = np.cumsum(lengths)

            ng50_where = np.where(cumsum > len(line_json["sequence"]) * 0.5)[0]
            ng90_where = np.where(cumsum > len(line_json["sequence"]) * 0.9)[0]

            ng50s.append(lengths[ng50_where[0]] if len(ng50_where) else 0)
            ng90s.append(lengths[ng90_where[0]] if len(ng90_where) else 0)
            assert(ng50s[-1] >= ng90s[-1])

        assert(len(ng50s) == len(ng90s))
        return np.array(scores), np.array(ng50s), np.array(ng90s), unmapped

assert(len(sys.argv) > 1)

stats = [align_stats(a) for a in sys.argv[1:]]
ylabels = ["Score", "NG50", "NG90"]
number_of_bins = 50
fig, axs = plt.subplots(len(ylabels))

# adapted from https://matplotlib.org/3.1.0/gallery/statistics/multiple_histograms_side_by_side.html#sphx-glr-gallery-statistics-multiple-histograms-side-by-side-py
for i, [ax, label] in enumerate(zip(axs, ylabels)):
    means = np.array([np.mean(a[i]) for a in stats])
    medians = np.array([np.median(a[i]) for a in stats])

    hist_range = (np.min([np.min(a[i]) for a in stats]), np.max([np.max(a[i]) for a in stats]))
    binned_data_sets = [np.histogram(a[i], range=hist_range, bins=number_of_bins)[0] for a in stats]
    binned_maximums = np.max(binned_data_sets, axis=1)
    x_locations = np.arange(0, sum(binned_maximums), np.max(binned_maximums))

    # The bin_edges are the same for all of the histograms
    bin_edges = np.linspace(hist_range[0], hist_range[1], number_of_bins + 1)
    centers = 0.5 * (bin_edges + np.roll(bin_edges, 1))[1:]
    heights = np.diff(bin_edges)

    # Cycle through and plot each histogram
    for x_loc, binned_data in zip(x_locations, binned_data_sets):
        lefts = x_loc - 0.5 * binned_data
        ax.barh(centers, binned_data, height=heights, left=lefts)

    ax.set_xticks(x_locations)
    ax.set_xticklabels(sys.argv[1:])

    ax.set_ylabel(label)

    #ax.plot(x_locations, means, marker='o')
    #ax.plot(x_locations, medians, marker='o')

for ax in fig.get_axes():
    ax.label_outer()

plt.show()
