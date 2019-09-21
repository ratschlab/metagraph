#!/usr/bin/env python3

import gzip
import glob
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

ks = [19, 27, 31, 41]

for k in ks:
    for f in sorted(glob.glob('k' + str(k) + "/*.diff*.fasta.gz")):
        lines = [a.decode().strip() for a in gzip.open(f, "rb")]
        seqs = lines[slice(1, len(lines), 2)]
        lens = sorted([len(a) for a in seqs])
        assert(lens[0] > 0)
        plt.semilogx(
            lens,
            np.arange(1, len(lens) + 1) / len(lens),
            label=f + " " + str(scipy.stats.entropy(lens)))

plt.xlabel("Unitig length")
plt.ylabel("CDF")
plt.legend()
plt.show()

