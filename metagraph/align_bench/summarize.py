#!/usr/bin/env python3
"""Run metrics.py over every finished arm under runs/ and write runs/summary.tsv.

usage: summarize.py [RUNS_DIR]   (from the work directory holding data/; default runs)
Columns: mode (labeled_<graph> / unlabeled_canonical), seeds, set, engine, then
the metrics.py columns.
"""
import glob, os, subprocess, sys
here = os.path.dirname(os.path.abspath(__file__))
root = sys.argv[1] if len(sys.argv) > 1 else 'runs'
rows, header = [], None
for tj in sorted(glob.glob(f'{root}/*/*/*.time.json')):
    prefix = tj[:-len('.time.json')]
    mode, seeds, base = prefix.split('/')[-3:]
    engine, rset = base.split('.', 1)
    out = subprocess.run([sys.executable, os.path.join(here, 'metrics.py'), 'data/genomes',
                          f'data/reads/{rset}.fa', prefix], capture_output=True, text=True, check=True)
    h, r = out.stdout.strip().split('\n')
    header = ['mode', 'seeds', 'set', 'engine'] + h.split('\t')[1:]
    rows.append([mode, seeds, rset, engine] + r.split('\t')[1:])
with open(f'{root}/summary.tsv', 'w') as f:
    f.write('\t'.join(header) + '\n')
    for r in rows: f.write('\t'.join(r) + '\n')
print(len(rows), 'arms')
