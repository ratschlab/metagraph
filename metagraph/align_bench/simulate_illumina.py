#!/usr/bin/env python3
"""Simulate Illumina-like reads with ground truth.

usage: simulate_illumina.py GENOME_DIR N_READS LENGTH SUB_RATE N_NEGATIVES SEED > reads.fa

Reads are sampled uniformly over bases of GENOME_DIR/*.fa (one genome per file),
substitution errors only, random strand. Names: r<i>|<genome>|<contig>|<start>|<strand>,
start is 0-based on the forward contig. Negatives are uniform random sequence: n<i>|NEG.
"""
import glob, os, random, sys

def read_fasta(path):
    name, seq = None, []
    for line in open(path):
        if line.startswith('>'):
            if name: yield name, ''.join(seq)
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip().upper())
    if name: yield name, ''.join(seq)

COMP = str.maketrans('ACGTN', 'TGCAN')
def rc(s): return s.translate(COMP)[::-1]

gdir, n, length, sub, nneg, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
rng = random.Random(seed)
contigs = []
for f in sorted(glob.glob(os.path.join(gdir, '*.fa'))):
    g = os.path.basename(f)[:-3]
    for c, s in read_fasta(f):
        if len(s) >= length: contigs.append((g, c, s))
weights = [len(s) - length + 1 for _, _, s in contigs]
for i in range(n):
    while True:
        g, c, s = rng.choices(contigs, weights)[0]
        p = rng.randrange(len(s) - length + 1)
        r = s[p:p + length]
        if 'N' not in r: break
    strand = rng.choice('+-')
    if strand == '-': r = rc(r)
    r = ''.join(rng.choice([b for b in 'ACGT' if b != x]) if rng.random() < sub else x for x in r)
    print(f'>r{i}|{g}|{c}|{p}|{strand}\n{r}')
for i in range(nneg):
    print(f'>n{i}|NEG\n' + ''.join(rng.choice('ACGT') for _ in range(length)))
