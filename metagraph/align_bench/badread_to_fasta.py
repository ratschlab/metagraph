#!/usr/bin/env python3
"""Convert badread FASTQ to FASTA with ground truth in the read name.

usage: badread_to_fasta.py CONTIG_MAP.tsv MIN_LEN MAX_LEN N_READS N_NEGATIVES SEED < reads.fq > reads.fa

CONTIG_MAP.tsv: contig<TAB>genome<TAB>contig length. Names: r<i>|<genome>|<contig>|<start>|<strand>,
start 0-based on the forward contig (badread reports '-' strand spans on the
reverse-complemented contig, so they are mirrored here). Reads outside [MIN_LEN, MAX_LEN] are dropped,
then the first N_READS are kept. Negatives: random sequence with lengths of
the kept reads, n<i>|NEG.
"""
import random, sys
cmap, clen = {}, {}
for l in open(sys.argv[1]):
    c, g, n_ = l.split()
    cmap[c], clen[c] = g, int(n_)
lo, hi, n, nneg, seed = map(int, sys.argv[2:7])
rng = random.Random(seed)
lines = sys.stdin.read().split('\n')
kept = []
for j in range(0, len(lines) - 3, 4):
    h, s = lines[j], lines[j + 1]
    f = h[1:].split()
    contig, strand, span = f[1].split(',')
    start, end = map(int, span.split('-'))
    if strand[0] == '-':
        start = clen[contig] - end
    if not (lo <= len(s) <= hi): continue
    kept.append((cmap[contig], contig, start, strand[0], s))
    if len(kept) == n: break
for i, (g, c, p, st, s) in enumerate(kept):
    print(f'>r{i}|{g}|{c}|{p}|{st}\n{s}')
for i in range(nneg):
    L = len(kept[i % len(kept)][4])
    print(f'>n{i}|NEG\n' + ''.join(rng.choice('ACGT') for _ in range(L)))
if len(kept) < n:
    sys.stderr.write(f'only {len(kept)} reads kept\n')
