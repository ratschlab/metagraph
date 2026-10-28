#!/usr/bin/env python3
"""Score align_engine_bench output against the simulated ground truth.

usage: metrics.py GENOME_DIR READS.fa OUT_PREFIX [OUT_PREFIX ...]  -> one TSV row per prefix

Truth comes from the read names (r<i>|<genome>|<contig>|<start>|<strand>, or n<i>|NEG).
Per prefix (reads with truth = positives P, negatives N):
  aligned      aligned positives / P
  neg_hit      aligned negatives / N
  top1_prec    positives whose first label is the true genome / aligned positives
  top1_recall  same count / P
  set_recall   positives whose label set contains the true genome / P
  set_size     mean label-set size over aligned positives
  pos_ok       aligned positives whose path spelling starts with a 31-mer found
               (either strand) in the true contig within [start-500, start+1.2*len+500],
               / aligned positives
  identity     mean (=)/(= + X + I + D) of the cigar over aligned positives
  qcov         mean aligned query fraction over aligned positives
  p50/p99/max latency (ms) over all reads that finished (OK or UNMAPPED)
  timeout      reads that hit the watchdog; error: reads that threw; crash: reads that crashed the process;
               skipped: reads not run after the arm hit 25 timeouts
  wall_s, rss_gb from OUT_PREFIX.time.json
"""
import glob, json, os, re, sys
import statistics as st

COMP = str.maketrans('ACGTN', 'TGCAN')
def rc(s): return s.translate(COMP)[::-1]

contigs = {}
for f in glob.glob(os.path.join(sys.argv[1], '*.fa')):
    name, seq = None, []
    for line in open(f):
        if line.startswith('>'):
            if name: contigs[name] = ''.join(seq)
            name, seq = line[1:].split()[0], []
        else: seq.append(line.strip().upper())
    contigs[name] = ''.join(seq)

lengths, name = {}, None
for line in open(sys.argv[2]):
    if line.startswith('>'): name = line[1:].split()[0]; lengths[name] = 0
    else: lengths[name] += len(line.strip())

def genome_of(label):
    return os.path.basename(label).removesuffix('.gz').removesuffix('.fa')

def pct(xs, q):
    xs = sorted(xs)
    return xs[min(len(xs) - 1, int(q * len(xs)))] if xs else float('nan')

cols = ['run', 'P', 'N', 'aligned', 'neg_hit', 'top1_prec', 'top1_recall', 'set_recall', 'set_size',
        'pos_ok', 'identity', 'qcov', 'lat_p50_ms', 'lat_p99_ms', 'lat_max_ms', 'timeout', 'error', 'crash', 'skipped',
        'wall_s', 'rss_gb', 'segments']
print('\t'.join(cols))
for prefix in sys.argv[3:]:
    P = N = al = neg = top1 = setr = pos = to = err = crash = skipped = 0
    sizes, ident, qcov, lat = [], [], [], []
    for line in open(prefix + '.tsv'):
        f = line.rstrip('\n').split('\t')
        name, status = f[0], f[1]
        truth = name.split('|')
        is_neg = truth[1] == 'NEG'
        if is_neg: N += 1
        else: P += 1
        if status == 'TIMEOUT': to += 1; continue
        if status == 'CRASH': crash += 1; continue
        if status == 'SKIPPED': skipped += 1; continue
        if status == 'ERROR': err += 1; continue
        lat.append(int(f[2]) / 1000)
        if status != 'OK': continue
        if is_neg: neg += 1; continue
        al += 1
        g, contig, start = truth[1], truth[2], int(truth[3])
        cigar, qb, qe, spelling, labels = f[4], int(f[6]), int(f[7]), f[8], f[9]
        labels = [genome_of(l) for l in labels.split(';') if l]
        if labels and labels[0] == g: top1 += 1
        if g in labels: setr += 1
        sizes.append(len(labels))
        ops = dict.fromkeys('=XID', 0)
        for n_, op in re.findall(r'(\d+)([=XIDSM])', cigar):
            if op in ops: ops[op] += int(n_)
        tot = sum(ops.values())
        ident.append(ops['='] / tot if tot else 0)
        L = lengths[name]
        qcov.append((qe - qb) / L)
        seed = spelling[:31]
        s = contigs[contig]
        lo = max(0, start - 500)
        win = s[lo:start + int(1.2 * L) + 500]
        if seed in win or rc(seed) in win: pos += 1
    row = dict(run=os.path.basename(prefix), P=P, N=N, aligned=al / P, neg_hit=neg / N if N else 0,
               top1_prec=top1 / al if al else 0, top1_recall=top1 / P, set_recall=setr / P,
               set_size=st.mean(sizes) if sizes else 0, pos_ok=pos / al if al else 0,
               identity=st.mean(ident) if ident else 0, qcov=st.mean(qcov) if qcov else 0,
               lat_p50_ms=pct(lat, .5), lat_p99_ms=pct(lat, .99), lat_max_ms=max(lat) if lat else 0,
               timeout=to, error=err, crash=crash, skipped=skipped)
    t = json.load(open(prefix + '.time.json'))
    row.update(wall_s=t['wall_s'], rss_gb=t['max_rss_kb'] / 1e6, segments=t['segments'])
    print('\t'.join(f'{row[c]:.3f}' if isinstance(row[c], float) else str(row[c]) for c in cols))
