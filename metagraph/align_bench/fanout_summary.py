#!/usr/bin/env python3
"""Summarize align_engine_bench fanout histograms.

usage: fanout_summary.py FANOUT.tsv [...]
For each file and each side (suffix = k-mers ending with the seed, prefix =
k-mers starting with it, i.e. forward resolution), over seed windows with at
least one match on that side: number of windows, p50, p90, p99, p99.9, max,
and the fraction of windows with fan-out <= 16 and > 64.
"""
import os, sys
print('\t'.join(['set', 'side', 'windows', 'matched', 'p50', 'p90', 'p99', 'p99.9', 'max', 'le16', 'gt64']))
for path in sys.argv[1:]:
    hist = {'suffix': {}, 'prefix': {}}
    for line in open(path):
        f = line.split()
        if f[0] in hist: hist[f[0]][int(f[1])] = int(f[2])
    for side, h in hist.items():
        total = sum(h.values())
        m = sorted((v, c) for v, c in h.items() if v > 0)
        n = sum(c for _, c in m)
        def q(p):
            acc = 0
            for v, c in m:
                acc += c
                if acc >= p * n: return v
        le16 = sum(c for v, c in m if v <= 16) / n
        gt64 = sum(c for v, c in m if v > 64) / n
        print('\t'.join(map(str, [os.path.basename(path).split('.')[0], side, total, n,
                                  q(.5), q(.9), q(.99), q(.999), m[-1][0], f'{le16:.4f}', f'{gt64:.5f}'])))
