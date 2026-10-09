#!/usr/bin/env python3
"""bench_pattern.py -- a repeatable, paced benchmark of a MetaGraph /pattern server.

Two suites:
  degenerate  how pattern search scales with increasingly degenerate periodic patterns: every p-th base is N
              (XXXXN, XXXN, XXN, XN), at growing lengths up to 50 bp (lengths > k as paths). Per cell, patterns
              cut from real RefSeq records (the mini index's FASTA, so at least one hit) and random ones.
  supported   graph paths vs record-supported paths for long peptides: conserved bacterial protein segments
              (found by their motifs in the FASTA, translated in six frames). Phase 1 slides a 20-aa window over
              each 50-aa segment and records the anchors per strand (admission against max_anchors) and the
              graph paths; phase 2 takes up to three admitted starts per protein at 20/30/40 aa and asks for the
              paths with their labels (mode partial, labels all): label-supported and record-verified paths among
              those returned. Plus the insulin A chain as a eukaryotic reference.
Stdlib only, Python >= 3.10. One request at a time.

USAGE
  # the plan only (no network)
  bench_pattern.py --base URL --suite degenerate --dry-run
  # a run; send it to the server's OWN port: the public ingress cuts requests after about 120 s, and the
  # extended (600 s budget) re-runs of time-stopped patterns need longer. On mex:
  bench_pattern.py --base http://localhost:42010 --suite degenerate,supported --label after-fix \\
      --fasta ~/pattern-bench/fasta
  # summary tables of a run (also written at the end of every run)
  bench_pattern.py --report OUT_DIR
  # before/after per cell (exit code 1 when an exact count differs)
  bench_pattern.py --compare OUT_BEFORE OUT_AFTER
  # pause between requests (e.g. for a quiet e2e run of another session): touch OUT_DIR/PAUSE; rm to resume
  # resume an interrupted run: the same command and --out; measured patterns are reused

COLD AND WARM
  Every pattern is sent twice in mode count: "first" (cold or partly cold page cache) and "repeat". A first
  call that stops on time is re-sent once per (cell, kind) with --extended-budget-ms instead of being
  repeated. On refseq33m (338 GiB graph, a 128 GiB container) first calls were 10-60x slower than repeats:
  compare runs only cell by cell, and only after a similar history (e.g. both right after a restart).

OUTPUT (OUT_DIR, default bench_out/pattern_<label>_<timestamp>)
  meta.json          base URL, arguments, release, index_fp, start and end time
  capabilities.json  GET /traverse/capabilities at the start
  results.jsonl      one row per request
  summary.md         the tables of --report
"""
import argparse
import glob
import json
import os
import random
import statistics
import sys
import time
import urllib.error
import urllib.request

FAMILIES = {5: [15, 20, 25, 30, 40, 50], 4: [16, 20, 24, 28, 40, 48],
            3: [18, 21, 24, 27, 30, 36, 45], 2: [24, 26, 28, 30, 32, 40, 50]}
FAMILY_NAME = {5: 'XXXXN', 4: 'XXXN', 3: 'XXN', 2: 'XN'}
# conserved bacterial protein motifs; the first record holding one (six frames, standard code) gives the
# 50-aa segment starting at the motif
MOTIFS = {'EF-Tu': 'GHVDHGKTTLTAAIT', 'GroEL': 'GDGTTTATVLA', 'DnaK': 'IDLGTTNS', 'RecA': 'GPESSGKTT',
          'RpoB': 'GDKMAGRHGNKG', 'ATPsynB': 'GGAGVGKTV', 'RplB': 'GRNNQGRITTRH'}
INSULIN_A = 'GIVEQCCTSICSLYQLENYCN'
SEED = 20261008
DEFAULT_FASTA = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'build', 'mini_refseq', 'fasta')

_B = 'TCAG'
_AA = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODE = {a + b + c: _AA[16 * i + 4 * j + m] for i, a in enumerate(_B) for j, b in enumerate(_B) for m, c in enumerate(_B)}


def rc(s):
    return s[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def translate(s):
    return ''.join(CODE.get(s[i:i + 3], 'X') for i in range(0, len(s) - 2, 3))


def read_fasta(directory):
    records = []
    for f in sorted(glob.glob(os.path.join(os.path.expanduser(directory), '*.fa'))):
        cur = None
        for line in open(f):
            if line.startswith('>'):
                cur = [os.path.basename(f), line[1:].split()[0], []]
                records.append(cur)
            elif cur is not None:
                cur[2].append(line.strip().upper())
    return [(f, name, ''.join(parts)) for f, name, parts in records]


class Client:
    def __init__(self, base, out_dir, timeout_s, dry_run):
        self.base, self.timeout_s, self.dry_run = base.rstrip('/'), timeout_s, dry_run
        self.pause = os.path.join(out_dir, 'PAUSE')
        self.sent = 0

    def get(self, path):
        with urllib.request.urlopen(self.base + path, timeout=60) as r:
            return json.load(r)

    def pattern(self, body):
        while os.path.exists(self.pause):
            time.sleep(10)
        if self.dry_run:
            print('  would send', json.dumps(body))
            return 0, {}, 0.0
        self.sent += 1
        t = time.time()
        req = urllib.request.Request(self.base + '/pattern', data=json.dumps(body).encode(),
                                     headers={'Content-Type': 'application/json'})
        try:
            with urllib.request.urlopen(req, timeout=self.timeout_s) as r:
                d, code = json.load(r), r.status
        except urllib.error.HTTPError as e:
            raw, code = e.read() or b'{}', e.code
            try:
                d = json.loads(raw)
            except ValueError:  # e.g. a proxy's HTML gateway timeout
                d = {'code': f'non_json_{code}', 'body': raw[:200].decode('utf-8', 'replace')}
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            d, code = {'code': f'client_{type(e).__name__}'}, 0
        return code, d, time.time() - t


def count_of(counts, key):
    x = (counts or {}).get(key) or {}
    return {'value': x.get('value'), 'relation': x.get('relation'), 'lower': x.get('lower'), 'upper': x.get('upper'),
            'by_strand': {s: (v or {}).get('value') for s, v in (x.get('by_strand') or {}).items()}}


def summarize(code, d, wall):
    row = {'http': code, 'wall_s': round(wall, 3)}
    if code != 200 or 'patterns' not in d:
        row['error'] = d.get('code') or d.get('error')
        return row
    e = d['patterns'][0]
    if e.get('error'):
        row['error'] = e['error'].get('code')
        row['bits'] = e.get('information_bits')
        return row
    c = e['counts']
    unit = 'paths' if c.get('paths') else 'contexts'
    main = count_of(c, unit)
    row.update({'unit': unit, 'count': main['value'], 'relation': main['relation'], 'lower': main['lower'],
                'upper': main['upper'], 'by_strand': main['by_strand'], 'bits': e.get('information_bits'),
                'anchors': count_of(c, 'anchors') if unit == 'paths' else None,
                'labels': count_of(c, 'labels'), 'occurrences': count_of(c, 'occurrences'),
                'work': e.get('work'), 'timing': e.get('timing'), 'server_ms': (e.get('timing') or {}).get('elapsed_ms'),
                'steps': (e.get('work') or {}).get('steps'), 'stop': e.get('stop'),
                'withheld': (e.get('withheld') or {}).get('reason'), 'cut': (e.get('cut') or {}).get('reason'),
                'returned': e.get('returned'), 'determinism': e.get('determinism'), 'notes': e.get('notes')})
    res = e.get('results') or []
    if res and any('labels' in r for r in res):
        row['label_supported'] = sum(1 for r in res if r.get('labels'))
        row['record_verified'] = sum(1 for r in res if any(l.get('support') == 'record_verified'
                                                           for l in (r.get('labels') or [])))
        row['distinct_labels_on_results'] = len({l.get('column') for r in res for l in (r.get('labels') or [])})
    return row


def stop_reason(r):
    return (r.get('stop') or {}).get('reason')


class Log:
    def __init__(self, path):
        self.path = path
        self.done = {}
        if os.path.exists(path):
            for line in open(path):
                r = json.loads(line)
                self.done.setdefault(r.get('key'), []).append(r)
        self.f = open(path, 'a')

    def write(self, rows):
        for r in rows:
            self.f.write(json.dumps(r) + '\n')
        self.f.flush()


def emit(r):
    extra = ''
    if r.get('label_supported') is not None:
        extra = f" label_sup {r['label_supported']} rec_ver {r['record_verified']}"
    anchors = (r.get('anchors') or {}).get('by_strand') if r.get('anchors') else None
    print(f"{r['key']:42s} {r['mode']:12s} {r['run']:8s} http {r['http']} {r['wall_s']:7.2f}s {r.get('unit', '')} "
          f"{r.get('count')} {r.get('relation')} anchors {anchors} steps {r.get('steps')} stop {stop_reason(r)} "
          f"withheld {r.get('withheld')} returned {r.get('returned')}{extra} err {r.get('error')}", flush=True)


# ------------------------------------------------------------------------------------------------ degenerate

def degenerate(seq, p):
    return ''.join('N' if i % p == p - 1 else c for i, c in enumerate(seq))


def suite_degenerate(cl, log, a, k, retrieve_max):
    rnd = random.Random(SEED)
    # the windows are cut from each FASTA file's records concatenated, so that the same seed gives the same
    # patterns and runs compare pattern by pattern
    files = {}
    for f, name, s in read_fasta(a.fasta):
        files.setdefault(f, []).append(s)
    records = [(f, '', ''.join(parts)) for f, parts in files.items()]
    if not records:
        sys.exit(f'no FASTA records in {a.fasta} (--fasta)')
    for p in [int(x) for x in a.families.split(',')]:
        limited_run = 0
        for L in FAMILIES[p]:
            if L > a.max_length:
                break
            cell = []
            for kind in ('real', 'random'):
                for j in range(a.per_cell):
                    if kind == 'real':
                        while True:
                            f, _, s = records[rnd.randrange(len(records))]
                            i = rnd.randrange(0, len(s) - L)
                            w = s[i:i + L]
                            if set(w) <= set('ACGT'):
                                break
                        source = f'{f}:{i}'
                    else:
                        w = ''.join(rnd.choice('ACGT') for _ in range(L))
                        source = 'random'
                    pat = degenerate(w, p)
                    key = f'degenerate.p{p}.L{L}.{kind}.{j}'
                    if key in log.done:
                        cell.append(next(r for r in log.done[key] if r['run'] == 'first'))
                        continue
                    base = {'patterns': [{'iupac': pat}]}
                    if L > k:
                        base['long_search'] = 'paths'
                    meta = {'suite': 'degenerate', 'key': key, 'period': p, 'family': FAMILY_NAME[p], 'length': L,
                            'n_count': pat.count('N'), 'kind': kind, 'source': source, 'pattern': pat}
                    rows = []

                    def run(mode, name, extra=None):
                        code, d, wall = cl.pattern({**base, 'mode': mode, **(extra or {})})
                        r = {**summarize(code, d, wall), **meta, 'mode': mode, 'run': name}
                        rows.append(r)
                        return r

                    first = run('count', 'first')
                    if not a.dry_run:
                        if stop_reason(first) != 'time':
                            run('count', 'repeat')
                        elif j == 0 and a.extended_budget_ms > 0:
                            run('count', 'extended', {'time_budget_ms': a.extended_budget_ms})
                        if first.get('relation') == 'exact' and (first.get('count') or 0) <= retrieve_max:
                            run('all_or_count', 'list')
                    log.write(rows)
                    for r in rows:
                        emit(r)
                    cell.append(first)
            limited = [c for c in cell if not c.get('error')]
            limited = bool(limited) and all(stop_reason(c) or c.get('relation') in ('at_least', 'unknown') for c in limited)
            limited_run = limited_run + 1 if limited else 0
            if limited_run >= 2:
                print(f'{FAMILY_NAME[p]}: two consecutive lengths hit a limit; the family stops at {L} bp', flush=True)
                break


# ------------------------------------------------------------------------------------------------- supported

def find_segments(records, length=50):
    found = {}
    for f, name, s in records:
        for strand, seq in (('+', s), ('-', rc(s))):
            for frame in range(3):
                prot = translate(seq[frame:])
                for m, motif in MOTIFS.items():
                    if m in found:
                        continue
                    i = prot.find(motif)
                    if i >= 0 and len(prot) >= i + length and '*' not in prot[i:i + length]:
                        found[m] = {'file': f, 'record': name, 'strand': strand, 'frame': frame, 'aa_start': i,
                                    'segment': prot[i:i + length]}
    return found


def suite_supported(cl, log, a, max_anchors, out_dir):
    segments = find_segments(read_fasta(a.fasta))
    json.dump(segments, open(os.path.join(out_dir, 'segments.json'), 'w'), indent=1)
    print('protein segments:', ', '.join(f"{m} ({v['file']} {v['record']})" for m, v in segments.items()), flush=True)

    def measure(tag, prot, off, L, pep, mode, extra=None):
        key = f'supported.{tag}.{prot}.{off}.{L}.{mode}'
        if key in log.done:
            return log.done[key][0]
        body = {'patterns': [{'protein': pep}], 'mode': mode, 'long_search': 'paths',
                'time_budget_ms': a.extended_budget_ms or 60000, **(extra or {})}
        code, d, wall = cl.pattern(body)
        r = {**summarize(code, d, wall), 'suite': 'supported', 'key': key, 'tag': tag, 'protein': prot,
             'offset': off, 'length_aa': L, 'pattern': pep, 'mode': mode, 'run': tag}
        log.write([r])
        emit(r)
        return r

    admitted = {}
    for prot, info in segments.items():
        seg = info['segment']
        for off in range(0, 31):
            r = measure('scan', prot, off, 20, seg[off:off + 20], 'count')
            anch = r.get('anchors') or {}
            if anch.get('relation') == 'exact' and anch.get('value') is not None and anch['value'] <= max_anchors:
                admitted.setdefault(prot, []).append((anch['value'], off))
    for prot, offs in admitted.items():
        seg = segments[prot]['segment']
        for _, off in sorted(offs)[:3]:
            for L in (20, 30, 40):
                if off + L <= len(seg):
                    measure('paths', prot, off, L, seg[off:off + L], 'count')
                    measure('support', prot, off, L, seg[off:off + L], 'partial', {'output': {'labels': 'all'}})
    measure('paths', 'insulinA', 0, len(INSULIN_A), INSULIN_A, 'count')
    measure('support', 'insulinA', 0, len(INSULIN_A), INSULIN_A, 'partial', {'output': {'labels': 'all'}})


# ---------------------------------------------------------------------------------------------------- report

def load(run_dir):
    return [json.loads(line) for line in open(os.path.join(run_dir, 'results.jsonl'))]


def med(xs):
    xs = [x for x in xs if x is not None]
    return statistics.median(xs) if xs else None


def fmt(x, spec='{:.2f}'):
    return '-' if x is None else spec.format(x)


def degenerate_cells(rows):
    cells = {}
    for r in rows:
        if r.get('suite') == 'degenerate':
            cells.setdefault((r['period'], r['length']), []).append(r)
    return dict(sorted(cells.items(), key=lambda kv: (-kv[0][0], kv[0][1])))


def report(run_dir):
    rows = load(run_dir)
    meta = json.load(open(os.path.join(run_dir, 'meta.json'))) if os.path.exists(os.path.join(run_dir, 'meta.json')) else {}
    out = [f"# Pattern benchmark {meta.get('label', '')}", '',
           f"release {meta.get('release')}, index_fp {str(meta.get('index_fp'))[:16]}, base {meta.get('base')}, "
           f"{meta.get('started')} .. {meta.get('finished')}", '']
    cells = degenerate_cells(rows)
    if cells:
        out += ['## Degenerate periodic patterns (medians per cell; relation e exact, a at_least, u unknown)', '',
                '| family | L | N | unit | count real | count random | first s (med / max) | repeat s | steps | '
                'time stops | extended runs | list s |', '|' + '---|' * 12]
        for (p, L), rs in cells.items():
            first = [r for r in rs if r['run'] == 'first']
            rel = lambda kind: ''.join(sorted({(r.get('relation') or 'x')[0] for r in first if r['kind'] == kind}))
            ext = ' '.join(f"{r['wall_s']:.0f}s:{r.get('count')}{(r.get('relation') or 'x')[0]}" for r in rs if r['run'] == 'extended')
            lst = [r['wall_s'] for r in rs if r['run'] == 'list']
            out.append(f"| {FAMILY_NAME[p]} | {L} | {first[0]['n_count']} | {first[0].get('unit') or first[0].get('error')} | "
                       f"{fmt(med([r.get('count') for r in first if r['kind'] == 'real']), '{:.0f}')} {rel('real')} | "
                       f"{fmt(med([r.get('count') for r in first if r['kind'] == 'random']), '{:.0f}')} {rel('random')} | "
                       f"{fmt(med([r['wall_s'] for r in first]))} / {fmt(max(r['wall_s'] for r in first))} | "
                       f"{fmt(med([r['wall_s'] for r in rs if r['run'] == 'repeat']))} | "
                       f"{fmt(med([r.get('steps') for r in first]), '{:.0f}')} | "
                       f"{sum(1 for r in first if stop_reason(r))} | {ext or '-'} | {fmt(med(lst))} |")
        out.append('')
    sup = [r for r in rows if r.get('suite') == 'supported']
    if sup:
        out += ['## Long peptides: anchors per start window (20 aa)', '',
                '| protein | windows admitted | anchors + (min..max) | anchors - (min..max) |', '|---|---|---|---|']
        for prot in dict.fromkeys(r['protein'] for r in sup if r['tag'] == 'scan'):
            sc = [r for r in sup if r['tag'] == 'scan' and r['protein'] == prot]
            ap = [((r.get('anchors') or {}).get('by_strand') or {}).get('+') or 0 for r in sc]
            am = [((r.get('anchors') or {}).get('by_strand') or {}).get('-') or 0 for r in sc]
            adm = sum(1 for r in sc if r.get('withheld') is None and (r.get('anchors') or {}).get('relation') == 'exact'
                      and r.get('count') is not None)
            out.append(f'| {prot} | {adm} of {len(sc)} | {min(ap)}..{max(ap)} | {min(am)}..{max(am)} |')
        out += ['', '## Long peptides: graph paths vs supported paths (support: among the paths returned, at most max_paths)', '',
                '| protein | start | aa | anchors | graph paths | returned | label-supported | record-verified | '
                'count s | labelled s | rows read |', '|' + '---|' * 11]
        for r in sup:
            if r['tag'] != 'paths':
                continue
            s = next((x for x in sup if x['tag'] == 'support' and x['protein'] == r['protein'] and x['offset'] == r['offset']
                      and x['length_aa'] == r['length_aa']), {})
            out.append(f"| {r['protein']} | {r['offset']} | {r['length_aa']} | {(r.get('anchors') or {}).get('value')} | "
                       f"{r.get('count')} {(r.get('relation') or '')[:1]} | {s.get('returned')} | {s.get('label_supported')} | "
                       f"{s.get('record_verified')} | {fmt(r['wall_s'])} | {fmt(s.get('wall_s'))} | "
                       f"{(s.get('work') or {}).get('annotation_rows')} |")
        out.append('')
    errors = [r for r in rows if r.get('http') != 200 or r.get('error')]
    if errors:
        out += ['## Errors', ''] + [f"- {r['key']} {r['mode']} {r['run']}: http {r['http']} {r.get('error')} "
                                    f"after {r['wall_s']:.1f} s" for r in errors] + ['']
    text = '\n'.join(out)
    open(os.path.join(run_dir, 'summary.md'), 'w').write(text)
    return text


def compare(dir_a, dir_b):
    ca, cb = degenerate_cells(load(dir_a)), degenerate_cells(load(dir_b))
    print(f'| family | L | first s A | first s B | repeat s A | repeat s B | steps A | steps B |\n' + '|' + '---|' * 8)
    changed = 0
    for key in ca:
        if key not in cb:
            continue
        fa = {r['key']: r for r in ca[key] if r['run'] == 'first'}
        fb = {r['key']: r for r in cb[key] if r['run'] == 'first'}
        for k2 in fa.keys() & fb.keys():
            if fa[k2].get('relation') == fb[k2].get('relation') == 'exact' and fa[k2].get('count') != fb[k2].get('count'):
                changed += 1
                print(f"EXACT COUNT CHANGED {k2}: {fa[k2].get('count')} -> {fb[k2].get('count')}", file=sys.stderr)
        g = lambda rows, run, f: med([f(r) for r in rows if r['run'] == run])
        print(f"| {FAMILY_NAME[key[0]]} | {key[1]} | {fmt(g(ca[key], 'first', lambda r: r['wall_s']))} | "
              f"{fmt(g(cb[key], 'first', lambda r: r['wall_s']))} | {fmt(g(ca[key], 'repeat', lambda r: r['wall_s']))} | "
              f"{fmt(g(cb[key], 'repeat', lambda r: r['wall_s']))} | {fmt(g(ca[key], 'first', lambda r: r.get('steps')), '{:.0f}')} | "
              f"{fmt(g(cb[key], 'first', lambda r: r.get('steps')), '{:.0f}')} |")
    return 1 if changed else 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument('--base', help='server base URL (prefer the server port, see USAGE)')
    ap.add_argument('--suite', default='degenerate,supported', help='comma list: degenerate, supported')
    ap.add_argument('--out', help='run directory (default bench_out/pattern_<label>_<timestamp>)')
    ap.add_argument('--label', default='run')
    ap.add_argument('--fasta', default=DEFAULT_FASTA, help='directory of *.fa records present in the index '
                    '(default: build/mini_refseq/fasta, whose records are RefSeq records)')
    ap.add_argument('--families', default='5,4,3,2', help='periods to run (5 XXXXN, 4 XXXN, 3 XXN, 2 XN)')
    ap.add_argument('--max-length', type=int, default=999, help='longest pattern of the degenerate suite')
    ap.add_argument('--per-cell', type=int, default=3, help='patterns per kind (real, random) and cell')
    ap.add_argument('--extended-budget-ms', type=int, default=600000,
                    help='budget of the re-run of a time-stopped pattern, and of the supported suite (0: none)')
    ap.add_argument('--timeout-s', type=float, default=900.0, help='client timeout per request')
    ap.add_argument('--dry-run', action='store_true', help='print the requests of the first pattern of each cell')
    ap.add_argument('--report', metavar='RUN_DIR', help='print and write RUN_DIR/summary.md')
    ap.add_argument('--compare', nargs=2, metavar=('A', 'B'), help='degenerate cells of two runs side by side')
    a = ap.parse_args(argv)
    if a.report:
        print(report(a.report))
        return 0
    if a.compare:
        return compare(*a.compare)
    if not a.base:
        ap.error('--base is required for a run')
    out_dir = a.out or os.path.join('bench_out', f"pattern_{a.label}_{time.strftime('%Y%m%d-%H%M%S')}")
    os.makedirs(out_dir, exist_ok=True)
    cl = Client(a.base, out_dir, a.timeout_s, a.dry_run)
    caps = cl.get('/traverse/capabilities')
    json.dump(caps, open(os.path.join(out_dir, 'capabilities.json'), 'w'))
    pcaps = (caps.get('pattern') or {})
    k = pcaps.get('k') or 31
    limits = pcaps.get('caps') or {}
    meta = {'label': a.label, 'base': a.base, 'args': vars(a), 'release': caps.get('release'),
            'index_fp': caps.get('index_fp'), 'k': k, 'started': time.strftime('%Y-%m-%d %H:%M:%S')}
    print(f"release {meta['release']} index_fp {str(meta['index_fp'])[:16]} k {k} out {out_dir}", flush=True)
    log = Log(os.path.join(out_dir, 'results.jsonl'))
    suites = a.suite.split(',')
    if 'degenerate' in suites:
        suite_degenerate(cl, log, a, k, limits.get('max_contexts', 10000))
    if 'supported' in suites:
        suite_supported(cl, log, a, limits.get('max_anchors', 1000), out_dir)
    meta['finished'] = time.strftime('%Y-%m-%d %H:%M:%S')
    meta['requests_sent'] = cl.sent
    json.dump(meta, open(os.path.join(out_dir, 'meta.json'), 'w'), indent=1)
    if not a.dry_run:
        print(report(out_dir))
    return 0


if __name__ == '__main__':
    sys.exit(main())
