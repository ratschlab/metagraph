#!/usr/bin/env python3
"""bench_library.py -- local benchmark of the metagraph.traverse Python library on recorded responses.

No network. It times what an agent does with one retrieval -- parse, summary, walks(top=5), claims,
label_walks / routes for one label and for every label, label() by name x1000, compare (pairs of the
same seed), to_fasta, to_json, to_gfa, memory_bytes -- on recorded /traverse responses, so that two
library versions can be compared on the same machine and the same inputs. Stdlib only, Python >= 3.10.

USAGE
  # the repository fixtures (tests/data/traverse/**.graphlet.json[.gz]) with the repository's library
  bench_library.py --label repo
  # recorded responses (bench_traverse raw/ dirs, the refseq-exp traverse/raw dir, any *.json[.gz])
  bench_library.py RUN_DIR/raw ../other/raw --max-items 40 --label before
  # another library version (the directory that holds the `metagraph` package)
  bench_library.py RUN_DIR/raw --lib /path/to/old/api/python --label old
  # before/after: per-op timing ratios + result fingerprints (exit code 1 when a result changed)
  bench_library.py --compare OUT_OLD/library_results.json OUT_NEW/library_results.json
  # two versions timed interleaved per (item, op) -- the reliable A/B on a busy machine
  bench_library.py RUN_DIR/raw --fixtures --interleave /path/to/old/api/python /path/to/new/api/python --rounds 3

INPUTS
  Files or directories (directories: *.json, *.json.gz, *.mgt, non-recursive unless --recursive). A JSON
  file is either a recorded request {meta, request, response} (bench_traverse raw/, the experiment's raw/)
  or a /traverse response {results: [...]} (the repository fixtures *.graphlet.json); every result with a
  `graphlet` is one item (Graphlet.from_response). An .mgt file is parsed with parse() (body only unless
  it carries a J line: ops that need the envelope then report MissingEnvelope). Identical graphlet bodies
  are timed once (--no-dedupe keeps them); compare() pairs are formed before deduplication from items of
  the same seed (validated seed id + sequence in the S line). With no inputs the repository fixtures are
  used (always this checkout's tests/data/traverse, whatever --lib is); --fixtures adds them to explicit
  inputs. Fixtures change with the tree: for a before/after pair over days, pass a fixed copy of the
  inputs (e.g. a bench_traverse raw/ directory) to both runs.

METHOD
  Timing: time.perf_counter, garbage collector off during the timed call. COLD = a fresh parse before
  every repetition (only the op is timed), so the library's per-graphlet caches start empty; min and
  median over the repetitions (7 / 5 / 3 for bodies < 20 KB / < 80 KB / larger; the all-label loops and
  label() x1000: 3 / 1). --warm also times the same op again on the same object (caches built).
  Memory: tracemalloc peak (and held, for parse: the parsed model) in a separate pass that is never timed.
  Fingerprints: every op's result is reduced to a stable digest (dataclasses by field, sets sorted, no
  object addresses) in another untimed pass; --compare flags (item, op) pairs whose digests differ.
  memory_bytes is informational (it changes whenever the model's layout does) and never flagged.
  Work: with a library that has stage L budgets (f667d775 and later), one more untimed call per (item, op)
  under local_budget() without limits records its usage: lwu (local work units at the cold price,
  deterministic) and model_kb (the modelled memory account's peak). It is the call's consumed work in the
  library's own units, not CPU time.
  --compare compares TIME ONLY ON THE SAME COMPLETED WORK: an (item, op) whose result fingerprints differ is
  flagged and left out of the ratios and sums (its times measure different work), and pairs with equal
  results whose lwu differ are listed (the work model or its charges changed).

INTERLEAVED A/B (--interleave LIB_A LIB_B)
  Two separate runs compare the machine's state as much as the code: on a loaded machine whole runs of one
  library differed by +-25 % on every op alike (A1 -> A2 1.20, A2 -> A3 0.77, 2026-10-04). --interleave starts
  one worker process per library (the same inputs, items and ops) and times each (item, op) in --rounds rounds
  of A B B A, each request the cold repetitions above; per (item, op) the minimum of each side, B/A per op over
  the pairs with equal result fingerprints (others listed and left out; exit code 1 when any changed), and a
  NOISE CONTROL: A's first request of each round against its last (A2/A1, the same code milliseconds apart),
  as a geometric mean and as the 90th percentile of its spread -- a B/A inside that spread is no difference.
  Output: ab_results.json (one row per (item, op): a_ms, b_ms, a1_ms, a2_ms, fp_a, fp_b, err) and summary.md.

OUTPUT (--out DIR, default bench_library_out/<label>_<timestamp>)
  library_results.json   meta (library path, module file, source digest, git head + dirty flag, python),
                         items, one row per (item, op): t_min_ms, t_med_ms, reps, peak_kb, held_kb, err, fp,
                         usage {lwu, model_kb} (null without stage L); compare rows
  summary.md             per-op table (median / p90 / max over items), size bins, the slowest calls, compare,
                         errors

The library is imported read-only: bytecode writing is switched off, so nothing is written into the tree.
"""

import argparse
import dataclasses
import datetime
import gc
import glob
import gzip
import hashlib
import json
import math
import os
import platform
import re
import statistics
import subprocess
import sys
import time
import tracemalloc

sys.dont_write_bytecode = True

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_LIB = os.path.normpath(os.path.join(HERE, '..', '..', 'api', 'python'))
FIXTURE_GLOBS = ['tests/data/traverse/documents/*.graphlet.json',
                 'tests/data/traverse/documents/*/*.graphlet.json',
                 'tests/data/traverse/review3/*.graphlet.json',
                 'tests/data/traverse/real/*/*.graphlet.json.gz']
ALL_OPS = ['parse', 'summary', 'walks_top5', 'claims', 'label_walks_one', 'routes_one', 'label_walks_all',
           'routes_all', 'label_by_name_x1000', 'to_fasta', 'to_json', 'to_gfa', 'memory_bytes']
LOOP_OPS = {'label_walks_all', 'routes_all', 'label_by_name_x1000'}
INFORMATIONAL = {'memory_bytes'}
COMPARE_MODES = ['claims', 'walks', 'labels', 'prefix_subset']

T = None          # metagraph.traverse, imported by load_library()
parser_mod = None


# ============================================================================ library

def load_library(path):
    global T, parser_mod
    path = os.path.abspath(path)
    if not os.path.isdir(os.path.join(path, 'metagraph', 'traverse')):
        raise SystemExit('--lib %s: no metagraph/traverse package there' % path)
    for m in [m for m in sys.modules if m == 'metagraph' or m.startswith('metagraph.')]:
        del sys.modules[m]
    sys.path.insert(0, path)
    import metagraph.traverse as traverse  # noqa: E402
    from metagraph.traverse import parser  # noqa: E402
    if not os.path.abspath(traverse.__file__).startswith(path + os.sep):
        raise SystemExit('imported %s, not the library under %s' % (traverse.__file__, path))
    T, parser_mod = traverse, parser
    return lib_meta(path, traverse.__file__)


def lib_meta(path, module_file):
    pkg = os.path.dirname(module_file)
    h = hashlib.sha256()
    for f in sorted(glob.glob(os.path.join(pkg, '*.py'))):
        h.update(os.path.basename(f).encode())
        with open(f, 'rb') as fh:
            h.update(fh.read())
    meta = {'path': path, 'module': module_file, 'source_sha': h.hexdigest()[:16], 'git_head': None,
            'git_dirty': None}
    try:
        head = subprocess.run(['git', '-C', pkg, 'rev-parse', '--short', 'HEAD'], capture_output=True,
                              text=True, timeout=20)
        if head.returncode == 0:
            meta['git_head'] = head.stdout.strip()
            st = subprocess.run(['git', '-C', pkg, 'status', '--porcelain', '--', '.'], capture_output=True,
                                text=True, timeout=20)
            meta['git_dirty'] = bool(st.stdout.strip()) if st.returncode == 0 else None
    except (OSError, subprocess.SubprocessError):
        pass
    return meta


# ============================================================================ inputs

def expand_inputs(paths, recursive):
    files = []
    for p in paths:
        if os.path.isdir(p):
            pats = ['*.json', '*.json.gz', '*.mgt']
            for pat in pats:
                files += glob.glob(os.path.join(p, '**', pat) if recursive else os.path.join(p, pat),
                                   recursive=recursive)
        elif os.path.exists(p):
            files.append(p)
        else:
            files += glob.glob(p)
    out, seen = [], set()
    for f in sorted(files):
        a = os.path.abspath(f)
        if a not in seen and not f.endswith(('results.json', 'library_results.json')):
            seen.add(a)
            out.append(a)
    return out


def fixture_files(root=REPO_LIB):
    """The repository's fixtures (always this checkout's, whatever --lib is, so two library versions are
    timed on the same inputs)."""
    files = []
    for g in FIXTURE_GLOBS:
        files += glob.glob(os.path.join(root, g))
    return sorted(files)


def read_json(path):
    if path.endswith('.gz'):
        with gzip.open(path, 'rt', encoding='utf-8') as f:
            return json.load(f)
    with open(path, encoding='utf-8') as f:
        return json.load(f)


def seed_key(body):
    for line in body.split('\n', 3)[:3]:
        if line.startswith('S '):
            f = line.split(' ')
            return f[-1] if len(f) >= 6 else line   # the sequence: equal across releases
    return None


def load_items(files, base_dir=None):
    items = []
    for path in files:
        rel = os.path.relpath(path, base_dir) if base_dir else os.path.basename(path)
        try:
            if path.endswith('.mgt'):
                with open(path, encoding='utf-8', newline='') as f:
                    text = f.read()
                items.append({'file': rel, 'idx': 0, 'kind': 'mgt', 'text': text, 'result': None, 'response': None,
                              'body': text})
                continue
            d = read_json(path)
        except (OSError, ValueError) as e:
            print('skip %s: %s' % (rel, e), file=sys.stderr)
            continue
        if not isinstance(d, dict):
            continue
        resp = d['response'] if isinstance(d.get('response'), dict) else d
        for i, r in enumerate(resp.get('results') or []):
            if isinstance(r, dict) and isinstance(r.get('graphlet'), str):
                items.append({'file': rel, 'idx': i, 'kind': 'response', 'result': r, 'response': resp,
                              'body': r['graphlet']})
    for it in items:
        it['body_bytes'] = len(it['body'].encode('utf-8'))
        it['lines'] = it['body'].count('\n')
        it['sha'] = hashlib.sha256(it['body'].encode('utf-8')).hexdigest()[:12]
        it['key'] = '%s#%d@%s' % (it['file'], it['idx'], it['sha'][:8])
        it['seed'] = seed_key(it['body'])
    return items


def parse_item(it):
    if it['kind'] == 'mgt':
        return parser_mod.parse(it['text'])
    return parser_mod.from_response(it['result'], it['response'])


# ============================================================================ ops

def sides(g):
    return list(g.arms.keys())


def unique_names(g):
    cnt = {}
    for l in g.labels:
        cnt[l.name] = cnt.get(l.name, 0) + 1
    return [n for n, c in cnt.items() if c == 1]


def op_table():
    def label_by_name(g):
        names = unique_names(g)
        if not names:
            raise LookupError('no uniquely named label')
        return [g.label({'name': names[j % len(names)]}).id for j in range(1000)]

    def first_label(g):
        if not g.labels:
            raise LookupError('no labels')
        return {'id': 0}

    return {
        'summary': lambda g: g.summary(),
        'walks_top5': lambda g: [g.walks(s, top=5) for s in sides(g)],
        'claims': lambda g: g.claims(),
        'label_walks_one': lambda g: g.label_walks(first_label(g)),
        'routes_one': lambda g: [g.routes(first_label(g), s, spell=True) for s in sides(g)],
        'label_walks_all': lambda g: [g.label_walks({'id': i}) for i in range(len(g.labels))],
        'routes_all': lambda g: [g.routes({'id': i}, s, spell=True) for i in range(len(g.labels)) for s in sides(g)],
        'label_by_name_x1000': label_by_name,
        'to_fasta': lambda g: g.to_fasta(),
        'to_json': lambda g: json.dumps(g.to_json(), separators=(',', ':'), sort_keys=True),
        'to_gfa': lambda g: g.to_gfa(),
        'memory_bytes': lambda g: g.memory_bytes(),
    }


def call(fn, g):
    try:
        return fn(g), None
    except Exception as e:  # an op that does not apply (annotate mode, body-only graphlet, old library)
        return None, '%s: %s' % (type(e).__name__, str(e)[:160])


# ============================================================================ fingerprints

def canonical(x, depth=0, seen=None):
    if seen is None:
        seen = set()
    if depth > 12:
        return '<deep>'
    if x is None or isinstance(x, (bool, int, str)):
        return x
    if isinstance(x, float):
        return round(x, 9) if math.isfinite(x) else str(x)
    if isinstance(x, (bytes, bytearray)):
        return hashlib.sha256(bytes(x)).hexdigest()
    if isinstance(x, dict):
        return {str(k): canonical(v, depth + 1, seen) for k, v in sorted(x.items(), key=lambda kv: str(kv[0]))}
    if isinstance(x, (list, tuple)):
        return [canonical(v, depth + 1, seen) for v in x]
    if isinstance(x, (set, frozenset)):
        return sorted((canonical(v, depth + 1, seen) for v in x), key=lambda v: json.dumps(v, sort_keys=True,
                                                                                          default=str))
    oid = id(x)
    if oid in seen:
        return '<cycle %s>' % type(x).__name__
    seen = seen | {oid}
    if hasattr(x, 'as_dict') and callable(x.as_dict):
        try:
            return {'_t': type(x).__name__, 'v': canonical(x.as_dict(), depth + 1, seen)}
        except Exception:
            pass
    if dataclasses.is_dataclass(x) and not isinstance(x, type):
        return {'_t': type(x).__name__,
                'v': {f.name: canonical(getattr(x, f.name, None), depth + 1, seen)
                      for f in dataclasses.fields(x) if f.name not in ('cache',) and f.compare is not False}}
    if hasattr(x, '__dict__'):
        return {'_t': type(x).__name__, 'v': canonical({k: v for k, v in vars(x).items() if not k.startswith('_')},
                                                       depth + 1, seen)}
    return '<%s>' % type(x).__name__


def fingerprint(x):
    return hashlib.sha256(json.dumps(canonical(x), sort_keys=True, default=str).encode()).hexdigest()[:16]


# ============================================================================ measurement

def reps_for(size, op, cap):
    if op in LOOP_OPS:
        r = 3 if size < 20000 else 1
    else:
        r = 7 if size < 20000 else 5 if size < 80000 else 3
    return max(1, min(r, cap)) if cap else r


def time_parse(it, reps):
    ts = []
    for _ in range(reps):
        gc.collect()
        gc.disable()
        t = time.perf_counter()
        try:
            parse_item(it)
        finally:
            ts.append(time.perf_counter() - t)
            gc.enable()
    return ts


def time_op(it, fn, reps, warm):
    cold, hot, err = [], [], None
    for _ in range(reps):
        g = parse_item(it)
        gc.collect()
        gc.disable()
        t = time.perf_counter()
        _, err = call(fn, g)
        cold.append(time.perf_counter() - t)
        if warm:
            t = time.perf_counter()
            call(fn, g)
            hot.append(time.perf_counter() - t)
        gc.enable()
        if err:
            break
    return cold, hot, err


def peak_of(fn):
    gc.collect()
    tracemalloc.start()
    try:
        tracemalloc.reset_peak()
        base = tracemalloc.get_traced_memory()[0]
        res = fn()
        cur, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    return res, (peak - base) / 1024.0, (cur - base) / 1024.0


def usage_of(fn):
    """Stage L's account of one call of |fn| under a budget without limits: {lwu (local work units, cold price:
    deterministic, the same in any process), model_kb (the modelled memory account's peak)}; None for a library
    without stage L (before f667d775) or a call that raised. Never timed."""
    lb = getattr(T, 'local_budget', None)
    if lb is None:
        return None
    try:
        with lb() as b:
            fn()
    except Exception:  # noqa -- an op that does not apply raises the same without a budget
        return None
    u = b.usage()
    return {'lwu': u.get('work_units'), 'model_kb': (u.get('memory_bytes') or 0) / 1024.0}


def measure(items, ops, args):
    table = op_table()
    rows = []
    t_start = time.time()
    for k, it in enumerate(items):
        try:
            g0 = parse_item(it)
        except Exception as e:
            rows.append({'key': it['key'], 'op': 'parse', 'err': '%s: %s' % (type(e).__name__, str(e)[:160])})
            continue
        it['labels'] = len(g0.labels)
        it['mode'] = g0.mode
        it['segments'] = sum(len(a.segments) for a in g0.arms.values())
        it['bases'] = sum(a.counts.bases for a in g0.arms.values())
        it['leaves'] = sum(a.counts.leaves for a in g0.arms.values())
        size = it['body_bytes']
        if 'parse' in ops:
            r = reps_for(size, 'parse', args.reps)
            ts = time_parse(it, r)
            g, peak, held = peak_of(lambda: parse_item(it))
            try:   # the canonical body text the parsed model writes back
                fp = fingerprint(parser_mod.dump(g, envelope=False))
            except Exception:  # an older library without dump(envelope=...)
                fp = None
            rows.append({'key': it['key'], 'op': 'parse', 't_min_ms': min(ts) * 1e3,
                         't_med_ms': statistics.median(ts) * 1e3, 'reps': r, 'peak_kb': peak, 'held_kb': held,
                         'err': None, 'fp': fp, 'usage': usage_of(lambda: parse_item(it))})
        for name in ops:
            if name == 'parse':
                continue
            fn = table[name]
            r = reps_for(size, name, args.reps)
            cold, hot, err = time_op(it, fn, r, args.warm)
            g = parse_item(it)
            (res, err2), peak, held = peak_of(lambda: call(fn, g))
            fp = None
            if err2 is None:
                try:
                    fp = fingerprint(res)
                except Exception as e:  # noqa
                    fp = 'unhashable: %s' % type(e).__name__
            g = parse_item(it)
            row = {'key': it['key'], 'op': name, 't_min_ms': min(cold) * 1e3, 't_med_ms': statistics.median(cold) * 1e3,
                   'reps': len(cold), 'peak_kb': peak, 'held_kb': held, 'err': err or err2, 'fp': fp,
                   'usage': usage_of(lambda: fn(g)) if not (err or err2) else None}
            if hot:
                row['t_warm_ms'] = min(hot) * 1e3
            rows.append(row)
        if args.progress:
            print('%3d/%d %-60s %7d B %4d labels  %.1f s' % (k + 1, len(items), it['key'][:60], size,
                                                           it.get('labels', 0), time.time() - t_start),
                  file=sys.stderr, flush=True)
    return rows


def compare_pairs(all_items, max_pairs):
    by_seed = {}
    for it in all_items:
        if it['seed']:
            by_seed.setdefault(it['seed'], []).append(it)
    pairs = []
    for its in by_seed.values():
        its = sorted(its, key=lambda x: (x['file'], x['idx']))
        for a, b in zip(its, its[1:], strict=False):
            pairs.append((a, b))
    pairs.sort(key=lambda p: -(p[0]['body_bytes'] + p[1]['body_bytes']))
    return pairs[:max_pairs]


def measure_compare(pairs, args):
    rows = []
    for a, b in pairs:
        for mode in COMPARE_MODES:
            ts, res, err = [], None, None
            for _ in range(3):
                try:
                    ga, gb = parse_item(a), parse_item(b)
                except Exception as e:
                    err = '%s: %s' % (type(e).__name__, str(e)[:160])
                    break
                gc.collect()
                gc.disable()
                t = time.perf_counter()
                res, err = call(lambda g: g.compare(gb, mode=mode), ga)
                ts.append(time.perf_counter() - t)
                gc.enable()
                if err:
                    break
            peak = None
            if not err:
                ga, gb = parse_item(a), parse_item(b)
                _, peak, _ = peak_of(lambda: call(lambda g: g.compare(gb, mode=mode), ga))
            rows.append({'a': a['key'], 'b': b['key'], 'mode': mode, 'bytes': a['body_bytes'] + b['body_bytes'],
                         't_min_ms': min(ts) * 1e3 if ts else None, 'peak_kb': peak, 'err': err,
                         'comparable': None if res is None else str(getattr(res, 'comparable', None)),
                         'equal': None if res is None else getattr(res, 'equal', None),
                         'fp': fingerprint(res) if res is not None else None})
    return rows


# ============================================================================ reporting

def fmt(x, nd=2):
    if x is None:
        return '-'
    if isinstance(x, float):
        if abs(x) >= 100:
            return '%.0f' % x
        if 0 < abs(x) < 0.1:
            return '%.2g' % x
        return ('%.' + str(nd) + 'f') % x
    return str(x)


def md_table(header, rows):
    out = ['| ' + ' | '.join(header) + ' |', '|' + '|'.join('---' for _ in header) + '|']
    out += ['| ' + ' | '.join(str(c).replace('|', '/') for c in r) + ' |' for r in rows]
    return '\n'.join(out)


def pct(xs, q):
    xs = sorted(xs)
    if not xs:
        return None
    i = min(len(xs) - 1, max(0, int(math.ceil(q * len(xs))) - 1))
    return xs[i]


def summary_md(res):
    meta, items, rows = res['meta'], res['items'], res['rows']
    by_key = {it['key']: it for it in items}
    L = ['# metagraph.traverse library benchmark: %s\n' % meta['label']]
    lib = meta['library']
    L.append('- library: `%s` (source sha %s, git %s%s)' % (lib['path'], lib['source_sha'], lib['git_head'],
                                                             ' dirty' if lib.get('git_dirty') else ''))
    L.append('- python %s on %s; started %s, took %.0f s' % (meta['python'], meta['platform'], meta['started'],
                                                            meta.get('seconds') or 0))
    sizes = [it['body_bytes'] for it in items]
    L.append('- %d graphlets (%s deduplicated) from %d files, bodies %s-%s bytes (median %s), labels up to %s' % (
        len(items), meta.get('deduplicated', 0), meta.get('files', 0), min(sizes) if sizes else '-',
        max(sizes) if sizes else '-', int(statistics.median(sizes)) if sizes else '-',
        max((it.get('labels') or 0) for it in items) if items else '-'))
    L.append('- method: %s\n' % meta['method'])
    ops = []
    for r in rows:
        if r['op'] not in ops:
            ops.append(r['op'])
    trows = []
    for op in ops:
        rs = [r for r in rows if r['op'] == op and r.get('t_min_ms') is not None and not r.get('err')]
        errs = sum(1 for r in rows if r['op'] == op and r.get('err'))
        if not rs:
            trows.append([op, 0, '-', '-', '-', '-', '-', '-', '-', errs])
            continue
        ts = [r['t_min_ms'] for r in rs]
        worst = max(rs, key=lambda r: r['t_min_ms'])
        lwu = [r['usage']['lwu'] for r in rs if (r.get('usage') or {}).get('lwu') is not None]
        trows.append([op, len(rs), fmt(statistics.median(ts)), fmt(pct(ts, 0.9)), fmt(worst['t_min_ms']),
                      '%s (%d B)' % (worst['key'][:48], by_key.get(worst['key'], {}).get('body_bytes', 0)),
                      fmt(statistics.median([r['peak_kb'] for r in rs if r.get('peak_kb') is not None] or [0]), 1),
                      fmt(max([r['peak_kb'] for r in rs if r.get('peak_kb') is not None] or [0]), 0),
                      fmt(statistics.median(lwu)) if lwu else '-', errs])
    L.append('## per op (cold: fresh parse before every repetition; t_min over repetitions)\n')
    L.append('median lwu = stage L local work units of one call (deterministic; "-" without stage L).\n')
    L.append(md_table(['op', 'items', 'median ms', 'p90 ms', 'max ms', 'slowest item', 'median peak KB',
                       'max peak KB', 'median lwu', 'errors'], trows))
    L.append('')
    bins = [(0, 5000, '< 5 KB'), (5000, 20000, '5-20 KB'), (20000, 80000, '20-80 KB'), (80000, 10 ** 12, '>= 80 KB')]
    brows = []
    for op in ops:
        cells = [op]
        for lo, hi, _ in bins:
            ts = [r['t_min_ms'] for r in rows if r['op'] == op and r.get('t_min_ms') is not None and not r.get('err')
                  and lo <= by_key.get(r['key'], {}).get('body_bytes', -1) < hi]
            cells.append('%s (%d)' % (fmt(statistics.median(ts)), len(ts)) if ts else '-')
        brows.append(cells)
    L.append('## median ms by body size (items)\n')
    L.append(md_table(['op'] + [b[2] for b in bins], brows))
    L.append('')
    slow = sorted([r for r in rows if r.get('t_min_ms') is not None and not r.get('err')],
                  key=lambda r: -r['t_min_ms'])[:12]
    L.append('## slowest calls\n')
    L.append(md_table(['op', 'item', 'bytes', 'labels', 't_min ms', 't_med ms', 'peak KB'],
                      [[r['op'], r['key'], by_key.get(r['key'], {}).get('body_bytes'),
                        by_key.get(r['key'], {}).get('labels'), fmt(r['t_min_ms']), fmt(r['t_med_ms']),
                        fmt(r.get('peak_kb'), 0)] for r in slow]))
    L.append('')
    if res.get('compare'):
        L.append('## compare() on pairs of the same seed\n')
        L.append(md_table(['a', 'b', 'mode', 'bytes', 't_min ms', 'peak KB', 'comparable', 'equal', 'error'],
                          [[c['a'][:40], c['b'][:40], c['mode'], c['bytes'], fmt(c['t_min_ms']),
                            fmt(c['peak_kb'], 0), c['comparable'], c['equal'], (c['err'] or '')[:60]]
                           for c in res['compare']]))
        L.append('')
    errs = [r for r in rows if r.get('err')]
    if errs:
        L.append('## errors (an op that does not apply to an item, or a library that lacks it)\n')
        seen = {}
        for r in errs:
            seen.setdefault((r['op'], r['err'][:80]), []).append(r['key'])
        L.append(md_table(['op', 'error', 'items', 'example'],
                          [[op, e, len(ks), ks[0][:50]] for (op, e), ks in sorted(seen.items())]))
        L.append('')
    return '\n'.join(L)


def compare_runs(pa, pb):
    with open(pa) as f:
        A = json.load(f)
    with open(pb) as f:
        B = json.load(f)
    ra = {(r['key'], r['op']): r for r in A['rows']}
    rb = {(r['key'], r['op']): r for r in B['rows']}
    common = [k for k in ra if k in rb]
    L = ['# library benchmark comparison: %s -> %s\n' % (A['meta']['label'], B['meta']['label'])]
    for tag, R in (('A', A), ('B', B)):
        lib = R['meta']['library']
        L.append('- %s: %s, `%s` (source sha %s, git %s%s), %d items' % (
            tag, R['meta']['label'], lib['path'], lib['source_sha'], lib['git_head'],
            ' dirty' if lib.get('git_dirty') else '', len(R['items'])))
    L.append('')
    ops = []
    for k in common:
        if k[1] not in ops:
            ops.append(k[1])
    rows, mism, lwu_changed = [], [], []
    for op in ops:
        rat, peak = [], []
        worst = None
        nfp = 0
        tot_a = tot_b = 0.0
        wa = wb = 0
        nw = 0
        for k in common:
            if k[1] != op:
                continue
            a, b = ra[k], rb[k]
            if a.get('err') or b.get('err') or not a.get('t_min_ms') or not b.get('t_min_ms'):
                if bool(a.get('err')) != bool(b.get('err')):
                    mism.append((k, 'error %s -> %s' % (a.get('err'), b.get('err'))))
                continue
            if op not in INFORMATIONAL and a.get('fp') and b.get('fp') and a['fp'] != b['fp']:
                # not the same completed work: its time is not compared
                nfp += 1
                mism.append((k, 'result fingerprint %s -> %s' % (a['fp'], b['fp'])))
                continue
            x = b['t_min_ms'] / a['t_min_ms']
            rat.append(x)
            tot_a += a['t_min_ms']
            tot_b += b['t_min_ms']
            if worst is None or x > worst[0]:
                worst = (x, k[0], a['t_min_ms'], b['t_min_ms'])
            if a.get('peak_kb') and b.get('peak_kb'):
                peak.append(b['peak_kb'] / a['peak_kb'])
            ua, ub = a.get('usage') or {}, b.get('usage') or {}
            if ua.get('lwu') is not None and ub.get('lwu') is not None:
                nw += 1
                wa += ua['lwu']
                wb += ub['lwu']
                if ua['lwu'] != ub['lwu']:
                    lwu_changed.append((k, ua['lwu'], ub['lwu']))
        if not rat:
            continue
        gm = math.exp(sum(math.log(x) for x in rat) / len(rat))
        rows.append([op, len(rat), fmt(gm, 3), fmt(statistics.median(rat), 3), fmt(tot_a), fmt(tot_b),
                     fmt(tot_b / tot_a if tot_a else None, 3),
                     '%s (%s -> %s ms) %s' % (fmt(worst[0], 2), fmt(worst[2]), fmt(worst[3]), worst[1][:40]),
                     fmt(statistics.median(peak), 2) if peak else '-', nfp if op not in INFORMATIONAL else 'n/a',
                     '%s -> %s' % (fmt(wa), fmt(wb)) if nw else '-'])
    L.append('B/A of t_min per (item, op), over the pairs that did the same completed work (equal result '
             'fingerprints); < 1 = B faster. "sum" = total over those items; "results changed" = pairs left out '
             'because their results differ. lwu = stage L local work units of the call (cold price, deterministic; '
             '"-" for a library without stage L), summed over the compared pairs.\n')
    L.append(md_table(['op', 'items', 'geomean B/A', 'median B/A', 'sum A ms', 'sum B ms', 'sum B/A',
                       'worst B/A (item)', 'median peak B/A', 'results changed', 'lwu A -> B'], rows))
    ca = {(c['a'], c['b'], c['mode']): c for c in A.get('compare') or []}
    cb = {(c['a'], c['b'], c['mode']): c for c in B.get('compare') or []}
    crow = []
    for k in [k for k in ca if k in cb]:
        a, b = ca[k], cb[k]
        same = a.get('fp') == b.get('fp')
        if not same:
            mism.append(((k[0] + ' vs ' + k[1]), 'compare(%s) %s -> %s' % (k[2], a.get('fp'), b.get('fp'))))
        crow.append([k[0][:36], k[1][:36], k[2], fmt(a.get('t_min_ms')), fmt(b.get('t_min_ms')),
                     fmt(b['t_min_ms'] / a['t_min_ms'] if a.get('t_min_ms') and b.get('t_min_ms') else None, 2),
                     'same' if same else 'CHANGED'])
    if crow:
        L.append('\n## compare()\n')
        L.append(md_table(['a', 'b', 'mode', 'A ms', 'B ms', 'B/A', 'result'], crow))
    L.append('')
    if mism:
        L.append('**%d results changed** (they must be equal unless the library\'s behaviour changed on purpose):\n'
                 % len(mism))
        for k, why in mism[:60]:
            L.append('- %s: %s' % (k if isinstance(k, str) else '%s %s' % (k[1], k[0]), why))
    else:
        L.append('No result changed (memory_bytes is informational and not compared).')
    if lwu_changed:
        L.append('\n%d (item, op) pairs with equal results charged different local work units (the work model '
                 'or the charges changed; informational):' % len(lwu_changed))
        for k, x, y in lwu_changed[:20]:
            L.append('- %s %s: %s -> %s lwu' % (k[1], k[0], x, y))
    only = len(set(ra) ^ set(rb))
    if only:
        L.append('\n(%d (item, op) rows are in only one of the two runs)' % only)
    return '\n'.join(L), mism


def select_items(args):
    """-> (files, all items, the items timed, duplicates left out) as the arguments choose them."""
    files = expand_inputs(args.inputs, args.recursive)
    if not args.inputs or args.fixtures:
        files += fixture_files()
    if not files:
        raise SystemExit('no input files')
    all_items = load_items(files)
    if not all_items:
        raise SystemExit('no graphlets in %d input files' % len(files))
    items, seen, dup = [], set(), 0
    for it in all_items:
        if not args.no_dedupe and it['sha'] in seen:
            dup += 1
            continue
        seen.add(it['sha'])
        items.append(it)
    if args.max_items and len(items) > args.max_items:
        srt = sorted(items, key=lambda x: x['body_bytes'])
        step = (len(srt) - 1) / float(args.max_items - 1) if args.max_items > 1 else 0
        pick = sorted({int(round(i * step)) for i in range(args.max_items)})
        items = [srt[i] for i in pick]
    return files, all_items, items, dup


# ============================================================================ interleaved A/B

def worker(args):
    """--worker LIB: one library version in its own process, timing on request. Writes a header line
    {library, items: [[key, bytes]...]}, then answers each line {i, op, reps} with {t: [ms...], fp, err}."""
    libm = load_library(args.worker)
    _, _, items, _ = select_items(args)
    table = op_table()
    print(json.dumps({'library': libm, 'items': [[it['key'], it['body_bytes']] for it in items]}), flush=True)
    fps = {}
    for line in sys.stdin:
        q = json.loads(line)
        it, op = items[q['i']], q['op']
        err = None
        if op == 'parse':
            ts = time_parse(it, q['reps'])
        else:
            ts, _, err = time_op(it, table[op], q['reps'], False)
        k = (q['i'], op)
        if k not in fps and not err and op not in INFORMATIONAL:
            try:
                if op == 'parse':
                    fps[k] = fingerprint(parser_mod.dump(parse_item(it), envelope=False))
                else:
                    res, e2 = call(table[op], parse_item(it))
                    fps[k] = None if e2 else fingerprint(res)
            except Exception as e:  # noqa
                fps[k] = 'unhashable: %s' % type(e).__name__
        print(json.dumps({'t': [x * 1e3 for x in ts], 'fp': fps.get(k), 'err': err}), flush=True)
    return 0


def interleave(args, ops):
    """--interleave LIB_A LIB_B: both versions in two worker processes, timed per (item, op) in rounds of
    A B B A (each request a few cold repetitions), so that the machine's drift between whole runs (+-25 % on a
    loaded machine) cancels; per (item, op) the minimum of each side. A noise control per (item, op): A's first
    request of each round against its last (A1 / A2, the same code a few milliseconds apart)."""
    passthrough = list(args.inputs)
    for flag in ('fixtures', 'recursive', 'no_dedupe'):
        if getattr(args, flag):
            passthrough.append('--' + flag.replace('_', '-'))
    if args.max_items:
        passthrough += ['--max-items', str(args.max_items)]
    procs, heads = [], []
    for lib in args.interleave:
        p = subprocess.Popen([sys.executable, os.path.abspath(__file__), '--worker', os.path.abspath(lib)]
                             + passthrough, stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, bufsize=1)
        procs.append(p)
        heads.append(json.loads(p.stdout.readline()))
    if heads[0]['items'] != heads[1]['items']:
        raise SystemExit('the two workers chose different items')

    def ask(p, i, op, reps):
        p.stdin.write(json.dumps({'i': i, 'op': op, 'reps': reps}) + '\n')
        p.stdin.flush()
        return json.loads(p.stdout.readline())
    started = datetime.datetime.now().astimezone().isoformat(timespec='seconds')
    t0 = time.time()
    rows = []
    for i, (key, size) in enumerate(heads[0]['items']):
        for op in ops:
            reps = reps_for(size, op, args.reps)
            a, b, a1, a2, fa, fb, err = [], [], [], [], None, None, None
            for _ in range(max(1, args.rounds)):
                for j, tag in enumerate('ABBA'):
                    x = ask(procs[0] if tag == 'A' else procs[1], i, op, reps)
                    if x['err']:
                        err = '%s: %s' % (tag, x['err'])
                        break
                    if tag == 'A':
                        a += x['t']
                        (a1 if j == 0 else a2).extend(x['t'])
                        fa = x['fp']
                    else:
                        b += x['t']
                        fb = x['fp']
                if err:
                    break
            rows.append({'key': key, 'op': op, 'bytes': size, 'err': err,
                         'a_ms': min(a) if a else None, 'b_ms': min(b) if b else None,
                         'a1_ms': min(a1) if a1 else None, 'a2_ms': min(a2) if a2 else None,
                         'fp_a': fa, 'fp_b': fb})
        if args.progress:
            print('%3d/%d %s  %.0f s' % (i + 1, len(heads[0]['items']), key[:60], time.time() - t0),
                  file=sys.stderr, flush=True)
    for p in procs:
        p.stdin.close()
        p.wait()
        p.stdout.close()
    meta = {'label': args.label, 'library_a': heads[0]['library'], 'library_b': heads[1]['library'],
            'python': sys.version.split()[0], 'platform': platform.platform(), 'started': started,
            'seconds': time.time() - t0, 'ops': ops, 'rounds': args.rounds, 'argv': sys.argv[1:],
            'method': 'two worker processes, one per library; per (item, op) rounds of A B B A, each request '
                      'cold repetitions by body size (fresh parse, gc off during the timed call); the minimum '
                      'per side; noise control A1/A2 = the first A request of each round against the last'}
    res = {'format': 'metagraph-bench-library-ab', 'version': 1, 'meta': meta, 'rows': rows}
    text, changed = interleave_md(res)
    stamp = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    out = args.out or os.path.join('bench_library_out', '%s_%s' % (re.sub(r'[^A-Za-z0-9._-]+', '_', args.label), stamp))
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, 'ab_results.json'), 'w') as f:
        json.dump(res, f, indent=1, default=str)
        f.write('\n')
    with open(os.path.join(out, 'summary.md'), 'w') as f:
        f.write(text + '\n')
    print(text)
    print('wrote %s/{ab_results.json,summary.md}: %d rows, %.0f s' % (out, len(rows), time.time() - t0))
    return 1 if changed else 0


def interleave_md(res):
    """The interleaved comparison as markdown -> (text, the rows whose results differ)."""
    meta, rows = res['meta'], res['rows']
    L = ['# library benchmark, interleaved A/B: %s\n' % meta['label']]
    for tag in ('a', 'b'):
        lib = meta['library_' + tag]
        L.append('- %s: `%s` (source sha %s, git %s%s)' % (tag.upper(), lib['path'], lib['source_sha'],
                                                            lib['git_head'], ' dirty' if lib.get('git_dirty') else ''))
    L.append('- %s; %d rounds of A B B A per (item, op); python %s; started %s, took %.0f s\n' % (
        meta['method'], meta['rounds'], meta['python'], meta['started'], meta.get('seconds') or 0))
    changed = [r for r in rows if not r['err'] and r['op'] not in INFORMATIONAL and r['fp_a'] != r['fp_b']]
    trows = []
    for op in meta['ops']:
        rs = [r for r in rows if r['op'] == op and not r['err'] and r['a_ms'] and r['b_ms'] and r not in changed]
        if not rs:
            continue
        x = [r['b_ms'] / r['a_ms'] for r in rs]
        noise = [r['a2_ms'] / r['a1_ms'] for r in rs if r['a1_ms'] and r['a2_ms']]
        # the spread of A2/A1 as a factor: 1.15 = the same code differed by up to 15 % (90th percentile)
        spread = math.exp(pct([abs(math.log(v)) for v in noise], 0.9)) if noise else None
        w = max(rs, key=lambda r: r['b_ms'] / r['a_ms'])
        trows.append([op, len(rs), fmt(geomean(x), 3), fmt(statistics.median(x), 3),
                      fmt(sum(r['a_ms'] for r in rs)), fmt(sum(r['b_ms'] for r in rs)),
                      fmt(sum(r['b_ms'] for r in rs) / sum(r['a_ms'] for r in rs), 3),
                      fmt(geomean(noise), 3) if noise else '-', fmt(spread, 2) if noise else '-',
                      '%s (%s -> %s ms) %s' % (fmt(w['b_ms'] / w['a_ms'], 2), fmt(w['a_ms']), fmt(w['b_ms']),
                                               w['key'][:40])])
    L.append('B/A of the minimum per (item, op) over the pairs with equal results; < 1 = B faster. Noise: A2/A1, '
             'the same code a few ms apart (its geometric mean, and the 90th percentile of its spread as a factor): '
             'a B/A within that spread is not a difference.\n')
    L.append(md_table(['op', 'items', 'geomean B/A', 'median B/A', 'sum A ms', 'sum B ms', 'sum B/A',
                       'noise A2/A1', 'noise p90 x', 'worst B/A (item)'], trows))
    errs = [r for r in rows if r['err']]
    if errs:
        L.append('\n%d (item, op) rows not timed (an op that does not apply): %s' % (
            len(errs), ', '.join(sorted({r['op'] for r in errs}))))
    if changed:
        L.append('\n**%d results changed** (left out of the ratios):\n' % len(changed))
        for r in changed[:60]:
            L.append('- %s %s: %s -> %s' % (r['op'], r['key'], r['fp_a'], r['fp_b']))
    else:
        L.append('\nNo result changed (memory_bytes is informational and not compared).')
    return '\n'.join(L), changed


def geomean(xs):
    return math.exp(sum(math.log(x) for x in xs) / len(xs)) if xs else None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog='See python3 bench_library.py --help-long for the inputs, the method, '
                                        'the interleaved A/B and the outputs.')
    ap.add_argument('inputs', nargs='*', help='files or directories of recorded responses (default: the '
                                              'repository fixtures)')
    ap.add_argument('--lib', default=REPO_LIB, help='directory holding the metagraph package (default: %s)' % REPO_LIB)
    ap.add_argument('--fixtures', action='store_true', help='add this checkout\'s fixtures (tests/data/traverse) '
                                                           'to the inputs')
    ap.add_argument('--recursive', action='store_true', help='descend into input directories')
    ap.add_argument('--max-items', type=int, default=0, help='at most N graphlets, spread evenly over the size '
                                                             'range (smallest and largest included)')
    ap.add_argument('--no-dedupe', action='store_true', help='time identical graphlet bodies more than once')
    ap.add_argument('--ops', default=','.join(ALL_OPS), help='comma list of ops (%s)' % ','.join(ALL_OPS))
    ap.add_argument('--no-compare', action='store_true', help='skip compare() on pairs of the same seed')
    ap.add_argument('--max-pairs', type=int, default=12)
    ap.add_argument('--reps', type=int, default=0, help='cap the repetitions per op (default: by body size)')
    ap.add_argument('--warm', action='store_true', help='also time the op again on the same object')
    ap.add_argument('--out', help='output directory (default bench_library_out/<label>_<timestamp>)')
    ap.add_argument('--label', default='lib')
    ap.add_argument('--progress', action='store_true', help='one line per item on stderr')
    ap.add_argument('--compare', nargs=2, metavar=('A', 'B'), help='compare two library_results.json')
    ap.add_argument('--compare-out', help='also write the comparison markdown here')
    ap.add_argument('--interleave', nargs=2, metavar=('LIB_A', 'LIB_B'),
                    help='time two library versions interleaved per (item, op) (see INTERLEAVED A/B)')
    ap.add_argument('--rounds', type=int, default=2, help='with --interleave: A B B A rounds per (item, op)')
    ap.add_argument('--worker', help=argparse.SUPPRESS)
    ap.add_argument('--help-long', action='store_true', help='print the full documentation')
    args = ap.parse_args(argv)

    if args.help_long:
        print(__doc__)
        return 0

    if args.compare:
        text, mism = compare_runs(*args.compare)
        print(text)
        if args.compare_out:
            with open(args.compare_out, 'w') as f:
                f.write(text + '\n')
        return 1 if mism else 0

    ops = [o.strip() for o in args.ops.split(',') if o.strip()]
    bad = [o for o in ops if o not in ALL_OPS]
    if bad:
        ap.error('unknown ops: %s' % ', '.join(bad))
    if args.worker:
        return worker(args)
    if args.interleave:
        return interleave(args, ops)
    libm = load_library(args.lib)
    files, all_items, items, dup = select_items(args)
    started = datetime.datetime.now().astimezone().isoformat(timespec='seconds')
    t0 = time.time()
    rows = measure(items, ops, args)
    cmp_rows = []
    if not args.no_compare:
        keys = {it['key'] for it in items}
        pool = [it for it in all_items if it['key'] in keys or not args.max_items]
        cmp_rows = measure_compare(compare_pairs(pool, args.max_pairs), args)
    meta = {'label': args.label, 'library': libm, 'python': sys.version.split()[0], 'platform': platform.platform(),
            'started': started, 'seconds': time.time() - t0, 'files': len(files), 'deduplicated': dup,
            'ops': ops, 'argv': sys.argv[1:],
            'method': 'perf_counter, gc off during the timed call; cold = a fresh parse before each repetition '
                      '(op only timed); min/median over 7/5/3 repetitions by body size (<20 KB/<80 KB/larger), '
                      'all-label loops and label() x1000 3/1; tracemalloc peak in a separate untimed pass'}
    out_items = [{k: it.get(k) for k in ('key', 'file', 'idx', 'kind', 'body_bytes', 'lines', 'sha', 'labels', 'mode',
                                          'segments', 'bases', 'leaves')} for it in items]
    res = {'format': 'metagraph-bench-library-results', 'version': 1, 'meta': meta, 'items': out_items,
           'rows': rows, 'compare': cmp_rows}
    stamp = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    out = args.out or os.path.join('bench_library_out', '%s_%s' % (re.sub(r'[^A-Za-z0-9._-]+', '_', args.label), stamp))
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, 'library_results.json'), 'w') as f:
        json.dump(res, f, indent=1, default=str)
        f.write('\n')
    with open(os.path.join(out, 'summary.md'), 'w') as f:
        f.write(summary_md(res))
    print('wrote %s/{library_results.json,summary.md}: %d graphlets x %d ops, %d compare rows, %.0f s' % (
        out, len(items), len(ops), len(cmp_rows), time.time() - t0))
    return 0


if __name__ == '__main__':
    sys.exit(main())
