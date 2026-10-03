"""SCALE: the largest cached retrievals through the library, the store and the MCP layer
(DESIGN-traverse-graphlet.md §5, §5.1, §6, §7, §14). SLOW: every test skips unless
METAGRAPH_REAL_SLOW=1 (the whole module takes a few minutes and a few GB of RAM).

    cd api/python/tests/real && METAGRAPH_REAL_SLOW=1 python3 -m unittest -v test_real_scale

What is measured, on the largest retrievals of the cache (SRA annotate beam 20 at 10 kb,
the design §7 reference case, ~3.6 MB of body; exhaustive constrain and label-free
tries; the 508-carrier hub seed), and asserted against generous bounds that catch a
blow-up, not a slow machine:

  * parse time, the model's retained heap and the parse peak (tracemalloc), the
    library's own memory_bytes() estimate against the traced heap;
  * dump (canonical: equal to the body), to_json (equal to the server's detail-full
    JSON after the §2.5 normalisation on the reference case), walks()/claims() over
    every arm, label_summary, summary (<= 2 KB), to_fasta / to_gfa;
  * compare (claims, walks, labels, prefix_subset) of the two largest comparable
    retrievals, with verified index identity (mini_refseq) and without (SRA);
  * the store's spool round trip (put -> restart -> parse from the spool -> save ->
    load, bodies deduplicated) and a deliberately small store (max_ram_mb) holding the
    beam / exhaustive cells: spooled delivery, LRU eviction, oversize never resident;
  * a cold MCP session (a fresh store and GraphletTools over a spool written by an
    earlier session), offline through a replaying client, and once live (SRA);
  * SUPER-LINEAR behaviour: the same operation on retrievals of increasing size, the
    exponent of time ~ size fitted on a log-log scale, against the work the operation
    has to do (W = MGT lines + the decoded label-set cardinalities + runs + the walks'
    chain lengths, + output bytes / 8 for the exports; the exponent against the body
    size is reported too: the body compresses label sets, so W grows faster than it).

Numbers go to <scratch>/scale_report.json (merged per section, so a partial run keeps
the rest). Product bugs found here are expected failures, each with a comment naming
the bug (see PRODUCT BUGS below); when one is fixed its test becomes an unexpected
success.

PRODUCT BUGS (measured on the cached retrievals, 2026-10-02; all three FIXED on
2026-10-03, their tests are plain tests again -- see the 'Fixed' note after each):

  B1 compare(mode='walks'|'prefix_subset') costs the whole trie, whatever the depth it
     compares: _walk_keys / _route_pairs call support_profile() (every P run of every
     segment of every walk's chain, as frozensets and sorted Label lists) and
     walk_bases() (the whole chain) for every walk before cutting at the common depth.
     sra_rand50_06__beam20 vs sra_rand50_06__beam1 at depth_used = 1 bp: walks 8.2 s,
     claims 0.05 s; sra_hub__beam20 vs __beam1 (depth 1): walks 107 s, prefix_subset
     195 s; time ~ walks^2.1 on the SRA beam ladder. The MCP graphlet_compare blocks as
     long, and §14's local-library budget (an interrupted comparison returns
     comparable: unknown) does not exist yet.
     Fixed: ops._cuts keys every walk by its cut [0, m), computed once per segment
     holding base m - 1 and reading the chain only up to m (memoized per comparison),
     and the walks sharing a cut are keyed once: the depth-1 pair now takes walks
     0.003 s, prefix_subset 0.005 s; on the ladder (mean of cold runs) the exponent
     is 1.08. The §14 budget is still missing (a design question).
  B2 GraphletStore(max_ram_mb=...) is not a bound on the heap: the store accounts each
     resident model by Graphlet.memory_bytes(), which counts 35-67 % of what
     tracemalloc traces for the model (the per-run frozenset `stars`, the ints/floats
     of R and T records, T extras tuples, EndLabel caches are not counted). Part of it
     is waste: parser.py _run_record builds a new frozenset per run (`frozenset()` is
     no singleton since Python 3.10: 216 B each), 3.3 MB of the 13.2 MB model of
     sra_hub__exhaustive (15,922 runs). A store with max_ram_mb=8 holds ~2-3x that.
     Fixed: shared star sets in the parser, memory_bytes() a deduplicating deep
     sizeof (worst ratio to the traced heap 1.00), and the store re-measures a
     resident model's caches on access and after every tool call (an 8 MB store
     now holds 7.8 MB of heap).
  B3 graphlet_walks recomputes the whole ranked walk list, spelling every walk, on every
     page (it calls Graphlet.walks(arm) without top): a 2 KB page costs as much as
     walks() over the whole arm (~0.5 s on the reference case), and the SRA rows (8
     labels with ~120-byte path names) fill a page with one row, so listing the 10,000
     walks of one arm would take ~10,000 x 0.5 s: quadratic in the walks.
     Fixed: graphlet_walks ranks the arm once on cheap keys (ops.rank_walks, cached
     per handle and arguments) and builds only its page's rows (ops.walks_at): a
     page costs 0.2 ms against walks() 0.26 s on the reference arm.
"""

import copy
import datetime
import gc
import json
import math
import os
import platform
import shutil
import sys
import tempfile
import time
import tracemalloc
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import realdata as R  # noqa: E402  (puts api/python on sys.path)

from metagraph.traverse import (  # noqa: E402
    GraphletStore, TraverseClient, derive, from_response, parse,
)
from metagraph.traverse.mcp_tools import (  # noqa: E402
    DEFAULT_MAX_BYTES, SEQUENCE_MAX_BYTES, GraphletTools,
)

SLOW = os.environ.get('METAGRAPH_REAL_SLOW') == '1'
# the 31.75 MB body (SRA 16S x beam 20, 403 MB of full JSON) joins the size ladder
# unless METAGRAPH_REAL_HEAVY=0
HEAVY = os.environ.get('METAGRAPH_REAL_HEAVY', '1') != '0'
REPORT_PATH = os.path.join(R.SCRATCH, 'scale_report.json')
TMP_ROOT = os.path.join(R.SCRATCH, 'scale_tmp')

slow = unittest.skipUnless(SLOW, 'slow scale suite: set METAGRAPH_REAL_SLOW=1')
# METAGRAPH_REAL_NO_XFAIL=1 runs the product-bug tests as plain tests (their real
# assertion messages, and proof that they fail on the assertion, not on an error)
xfail = (lambda f: f) if os.environ.get('METAGRAPH_REAL_NO_XFAIL') == '1' \
    else unittest.expectedFailure

MB = float(1 << 20)

# the largest retrievals of the cache, one per class (a cell missing from the cache
# is skipped, not failed)
REFERENCE = ('sra', 'sra_rand50_00__beam20_10k')   # design §7: beam 20, 10 kb, ~3.6 MB
TARGETS = [
    REFERENCE + ('SRA annotate beam 20, 10 kb (design §7 reference case)',),
    ('sra', 'sra_rand50_01__beam20_10k', 'SRA annotate beam 20, 10 kb, 10,000-path cap'),
    ('uhgg', 'uhgg_16s__beam20', 'UHGG annotate beam 20, 5 kb, 16S window'),
    ('sra', 'sra_hub__beam20', '508-carrier hub, annotate beam 20, 3 kb'),
    ('sra', 'sra_hub__beam1', '508-carrier hub, annotate beam 1, 10 kb'),
    ('sra', 'sra_hub__exhaustive', '508-carrier hub, exhaustive constrain trie, 300 bp'),
    ('sra', 'sra_hub__limit2', '508-carrier hub, constrain, branch limit 2, merge'),
    ('sra', 'sra_hub__cut_events', '508-carrier hub, exhaustive constrain, events cut'),
    ('sra', 'sra_rand50_00__exhaustive', 'exhaustive constrain trie, 35 carriers'),
    ('sra', 'sra_16s_PZ326290__exhaustive', 'exhaustive constrain trie, 16S, 6 carriers'),
    ('sra', 'sra_16s_OK510341__annotate_exh', 'largest label-free (annotate) exhaustive trie'),
    ('sra', 'sra_hub__annotate_exh', '508-carrier hub, label-free exhaustive trie'),
    ('mini_refseq', 'mini_ndm1_rc__beam20', 'largest mini_refseq retrieval (verified index)'),
]
HEAVIEST = ('sra', 'sra_16s_PZ326290__beam20')

# size ladders (same index, same strategy class), increasing body
LADDER_BEAM = ['sra_rand50_05__beam20', 'sra_rand50_04__beam20', 'sra_rand50_02__beam20',
               'sra_rand50_08__beam20', 'sra_rand50_09__beam20', 'sra_rand50_07__beam20',
               'sra_rand50_00__beam20', 'sra_rand50_06__beam20', 'sra_16s_OK510341__beam20',
               'sra_rand50_00__beam20_10k']
LADDER_EXH = ['sra_rand50_06__exhaustive', 'sra_rand50_04__exhaustive',
              'sra_rand50_01__exhaustive', 'sra_16s_OK510341__exhaustive',
              'sra_16s_PZ326290__exhaustive', 'sra_rand50_00__exhaustive',
              'sra_hub__exhaustive']
# compare ladder: beam 20 vs beam 1 of one seed (comparison depth <= 190 bp)
LADDER_CMP = ['sra_rand50_04', 'sra_rand50_02', 'sra_rand50_08', 'sra_rand50_09',
              'sra_rand50_07', 'sra_rand50_00', 'sra_rand50_06']
DEPTH1_PAIR = ('sra', 'sra_rand50_06__beam20', 'sra_rand50_06__beam1')

TO_JSON_MAX_FULL = 160 << 20     # to_json is skipped where the server's full JSON is larger
SUPERLINEAR = 1.35               # an exponent above this (vs the work W) is flagged
MIN_FIT_T = 1e-3                 # points faster than this are noise for the fit
OPS = ('parse', 'dump', 'memory', 'to_json', 'walks', 'claims', 'label_summary',
       'summary', 'to_fasta', 'to_gfa')

# ------------------------------------------------------------------ the report

REPORT = {}


def record(section, key, value):
    REPORT.setdefault(section, {})[key] = value


def _git_head():
    try:
        with open(os.path.join(R.REPO, '..', '.git', 'HEAD')) as f:
            ref = f.read().strip()
        if ref.startswith('ref: '):
            with open(os.path.join(R.REPO, '..', '.git', ref[5:])) as f:
                return f.read().strip()
        return ref
    except OSError:
        return None


def tearDownModule():
    if not REPORT:
        return
    try:
        with open(REPORT_PATH) as f:
            out = json.load(f)
    except (OSError, ValueError):
        out = {}
    for section, rows in REPORT.items():
        out.setdefault(section, {}).update(rows)
    out['meta'] = {
        'updated': datetime.datetime.now().isoformat(timespec='seconds'),
        'python': sys.version.split()[0], 'platform': platform.platform(),
        'cpu_count': os.cpu_count(), 'git_head': _git_head(),
        'cache_updated': R.load_manifest().get('updated'),
        'superlinear_threshold': SUPERLINEAR,
        'work_measure': 'W = MGT lines + decoded label-set cardinalities (G entry/end, '
                        'P runs) + R runs + the walks\' chain lengths (+ output bytes / 8 '
                        'for to_fasta / to_gfa)',
    }
    os.makedirs(os.path.dirname(REPORT_PATH), exist_ok=True)
    tmp = REPORT_PATH + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(out, f, indent=1, sort_keys=True)
    os.replace(tmp, REPORT_PATH)


# ------------------------------------------------------------------ helpers

def _cell(index, cell):
    R.require_cache(index, cell)
    return R.load_cell(index, cell)


def _key(index, cell):
    return '%s/%s' % (index, cell)


def _timed(fn, reps=1):
    """(best wall time of |reps| runs, the last result)."""
    best = out = None
    for _ in range(reps):
        out = None
        gc.collect()
        t0 = time.perf_counter()
        out = fn()
        dt = time.perf_counter() - t0
        best = dt if best is None else min(best, dt)
    return best, out


class _Cold:
    """Resets a graphlet's derived caches to what the parser left (the per-arm caches
    hold derivations; the parser itself fills a few): each operation is timed cold."""

    def __init__(self, g):
        self.g = g
        self.base = (dict(g.cache), {s: dict(a.cache) for s, a in g.arms.items()})

    def reset(self):
        self.g.cache.clear()
        self.g.cache.update(self.base[0])
        for s, a in self.g.arms.items():
            a.cache.clear()
            a.cache.update(self.base[1][s])
        return self.g


def _cold_mean(fn, colds, at_least=0.05, max_reps=200):
    """(mean wall time of fn(), the last result, the runs): fn is run cold -- every _Cold
    of |colds| reset before each run, untimed -- and repeated until the runs total
    |at_least| seconds (or |max_reps| runs), so that an operation of a millisecond is
    timed above the timer's noise without warming its caches."""
    total, reps, out = 0.0, 0, None
    gc.collect()
    while reps == 0 or (total < at_least and reps < max_reps):
        for c in colds:
            c.reset()
        t0 = time.perf_counter()
        out = fn()
        total += time.perf_counter() - t0
        reps += 1
    return total / reps, out, reps


def _traced_parse(body):
    """-> (graphlet, retained bytes, peak bytes) of parse(body) under tracemalloc."""
    gc.collect()
    tracemalloc.start()
    try:
        base = tracemalloc.get_traced_memory()[0]
        tracemalloc.reset_peak()
        g = parse(body)
        gc.collect()
        cur, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    return g, cur - base, peak - base


def _content(g, body):
    """Size measures of one model: E (lines + decoded set cardinalities + runs) and the
    walks' chain lengths."""
    sets = 0
    runs = 0
    for a in g.arms.values():
        runs += len(a.runs)
        for s in a.segments:
            sets += len(s.entry) + len(s.end) + sum(len(p.labels) for p in s.presence)
    chain = sum(len(derive.chain(a, p.leaf)) for a in g.arms.values()
                for p in derive.paths(a))
    return {'lines': body.count('\n'), 'set_cardinality': sets, 'runs': runs,
            'E': body.count('\n') + sets + runs, 'chain': chain,
            'segments': sum(len(a.segments) for a in g.arms.values()),
            'walks': sum(len(derive.paths(a)) for a in g.arms.values()),
            'labels': len(g.labels), 'mode': g.mode}


_MEASURED = {}


def measure(index, cell):
    """Every operation of OPS on results[0] of |cell|, each on a cold model; memoised.
    Keeps numbers only (no model) between calls. On the reference case to_json() is
    also checked against the server's detail-full JSON."""
    key = _key(index, cell)
    if key in _MEASURED:
        return _MEASURED[key]
    to_json_check = (index, cell) == REFERENCE
    c = _cell(index, cell)
    resp = c.graphlet_response
    r = resp['results'][0]
    body = r['graphlet']
    full_bytes = (c.meta.get('full') or {}).get('json_bytes') or 0
    reps = 3 if len(body) < 300000 else 1
    m = {'cell': key, 'body_bytes': len(body), 'full_json_bytes': full_bytes,
         'graphlet_bytes_server': r.get('graphlet_bytes')}
    m['parse_s'], g = _timed(lambda: parse(body), reps)
    m['from_response_s'], g = _timed(lambda: from_response(r, resp), reps)
    m.update(_content(g, body))
    m['memory_bytes'] = g.memory_bytes()
    cold = _Cold(g)
    for a in g.arms.values():
        m.setdefault('complete_to_bp', {})[a.side] = a.complete_to_bp

    def op(name, fn):
        best = None
        out = None
        for _ in range(reps):
            cold.reset()
            out = None
            dt, out = _timed(fn)
            best = dt if best is None else min(best, dt)
        m[name + '_s'] = best
        return out

    text = op('dump', lambda: g.dump(envelope=False))
    m['dump_equal'] = text == body
    del text
    ws = op('walks', lambda: {s: g.walks(s) for s in g.arms})
    m['walks_n'] = {s: len(v) for s, v in ws.items()}
    m['walks_bp'] = sum(w.length_bp for v in ws.values() for w in v)
    m['walks_segments'] = sum(len(w.segments) for v in ws.values() for w in v)
    m['walks_leaves_ok'] = all(sorted(w.path_id for w in ws[s]) ==
                               list(range(len(derive.paths(g.arms[s])))) for s in ws)
    del ws
    cl = op('claims', lambda: {s: g.claims(s, strict=False) for s in g.arms})
    m['claims_n'] = {s: len(v) for s, v in cl.items()}
    del cl
    ls = op('label_summary', lambda: g.label_summary())
    m['label_summary_rows'] = len(ls)
    m['label_summary_ok'] = len(ls) == len({l.ref for l in g.labels}) and all(
        set(g.arms) <= set(row) for row in ls.values())
    del ls
    sm = op('summary', lambda: g.summary())
    m['summary_bytes'] = len(json.dumps(sm, separators=(',', ':'), ensure_ascii=False)
                             .encode('utf-8'))
    del sm
    fa = op('to_fasta', lambda: g.to_fasta())
    m['fasta_bytes'] = len(fa)
    m['fasta_records'] = fa.count('>')
    del fa
    gfa = op('to_gfa', lambda: g.to_gfa())
    m['gfa_bytes'] = len(gfa)
    kinds = {}
    for line in gfa.split('\n'):
        if line:
            kinds[line[0]] = kinds.get(line[0], 0) + 1
    m['gfa_records'] = kinds
    m['gfa_links_expected'] = sum(max(1, len(s.parents)) for a in g.arms.values()
                                  for s in a.segments)
    del gfa
    if full_bytes <= TO_JSON_MAX_FULL:
        j = op('to_json', lambda: g.to_json())
        m['to_json_paths'] = {s: len(a['paths']) for s, a in j['arms'].items()}
        if to_json_check:
            t0 = time.perf_counter()
            full = json.loads(c.raw('full'))['results'][0]
            mine = R.normalize_full_result(j)
            theirs = R.normalize_full_result(full)
            m['to_json_equal'] = mine == theirs
            if not m['to_json_equal']:
                m['to_json_diff'] = R.json_diff(mine, theirs)
            m['to_json_check_s'] = time.perf_counter() - t0
            del full, mine, theirs
        del j
    else:
        m['to_json_s'] = None
        m['to_json_skipped'] = 'full JSON %.0f MB > %.0f MB' % (full_bytes / MB,
                                                                 TO_JSON_MAX_FULL / MB)
    del cold, g
    gc.collect()
    g, retained, peak = _traced_parse(body)
    m['model_retained_bytes'] = retained
    m['parse_peak_bytes'] = peak
    del g
    gc.collect()
    m['W'] = m['E'] + m['chain']
    _MEASURED[key] = m
    return m


def fit_exponent(xs, ys, min_y=MIN_FIT_T):
    """Least-squares slope of log y ~ log x over the points with x > 0, y >= min_y ->
    (slope, r2, n) or None (fewer than 4 points, or x spanning less than a decade)."""
    pts = [(math.log(x), math.log(y)) for x, y in zip(xs, ys)
           if x and y is not None and x > 0 and y >= min_y]
    if len(pts) < 4:
        return None
    lx = [p[0] for p in pts]
    if max(lx) - min(lx) < math.log(10):
        return None
    mx = sum(lx) / len(pts)
    my = sum(p[1] for p in pts) / len(pts)
    sxx = sum((p[0] - mx) ** 2 for p in pts)
    sxy = sum((p[0] - mx) * (p[1] - my) for p in pts)
    syy = sum((p[1] - my) ** 2 for p in pts)
    slope = sxy / sxx
    r2 = (sxy * sxy) / (sxx * syy) if syy else 1.0
    return round(slope, 3), round(r2, 3), len(pts)


def _round(x, n=4):
    return None if x is None else round(x, n)


def _brief(m):
    """The report row of one measured cell."""
    keys = ('body_bytes', 'full_json_bytes', 'lines', 'segments', 'walks', 'labels', 'mode',
            'runs', 'set_cardinality', 'chain', 'W', 'parse_s', 'from_response_s', 'dump_s',
            'to_json_s', 'walks_s', 'claims_s', 'label_summary_s', 'summary_s',
            'to_fasta_s', 'to_gfa_s', 'model_retained_bytes', 'parse_peak_bytes',
            'memory_bytes', 'fasta_bytes', 'gfa_bytes', 'summary_bytes', 'claims_n',
            'complete_to_bp', 'to_json_equal', 'to_json_check_s', 'to_json_skipped')
    out = {}
    for k in keys:
        if k in m:
            v = m[k]
            out[k] = _round(v) if isinstance(v, float) else v
    if m.get('model_retained_bytes'):
        out['memory_bytes_ratio'] = round(m['memory_bytes'] / m['model_retained_bytes'], 3)
    return out


# ------------------------------------------------------------------ a replaying client

class ReplayClient(TraverseClient):
    """A TraverseClient that answers /traverse from the retrieval cache (the request
    as built by the tool layer, timing ignored) and capabilities from the recorded ones:
    the MCP layer runs offline, deterministically, on the real responses."""

    def __init__(self, index, cells):
        super().__init__('127.0.0.1', 1, timeout=1)
        self.index = index
        self.calls = 0
        self._by_request = {}
        for cell in cells:
            c = R.load_cell(index, cell)
            self._by_request[self._norm(c.request)] = c

    @staticmethod
    def _norm(req):
        r = copy.deepcopy(req)
        r['strategy'].setdefault('output', {})['timing'] = False
        return json.dumps(r, sort_keys=True, separators=(',', ':'))

    def capabilities(self):
        rec = R.recorded_capabilities(self.index) or {}
        return copy.deepcopy(rec.get('capabilities') or {})

    def traverse_raw(self, request):
        self.calls += 1
        c = self._by_request.get(self._norm(request))
        if c is None:
            raise AssertionError('ReplayClient: no cached cell for %s'
                                 % json.dumps(request)[:300])
        return json.loads(c.raw('graphlet'))


def _strategy_of(c):
    """The cell's strategy as an agent passes it (output.detail / timing are the tool
    layer's)."""
    st = copy.deepcopy(c.request['strategy'])
    out = {k: v for k, v in (st.pop('output', None) or {}).items()
           if k not in ('detail', 'timing')}
    if out:
        st['output'] = out
    return st


def _as_delivered(body, delivery):
    """|body| with the O record's delivery dimension set to |delivery|: what the model of
    a spooled entry dumps (store.set_delivery makes the outcome say what happened; the
    spool keeps the body as received)."""
    code = {'inline': 'i', 'spooled': 's', 'paged': 'p'}[delivery]
    lines = body.split('\n')
    for i, line in enumerate(lines):
        if line.startswith('O '):
            f = line.split(' ')
            f[-1] = code
            lines[i] = ' '.join(f)
            break
    return '\n'.join(lines)


def _spooled_body(store, handle):
    with open(store._body_path(store.get(handle).digest), 'r', encoding='utf-8',
              newline='') as f:
        return f.read()


def _size(obj):
    return len(json.dumps(obj, separators=(',', ':'), ensure_ascii=False).encode('utf-8'))


def _tmpdir(tag):
    os.makedirs(TMP_ROOT, exist_ok=True)
    return tempfile.mkdtemp(prefix=tag + '-', dir=TMP_ROOT)


# ====================================================================== the largest

@slow
class TestLargestRetrievals(unittest.TestCase):
    """Every operation on the largest retrievals, against bounds relative to their size
    (generous: they flag a 10x regression, not a slow machine)."""

    def _targets(self):
        out = []
        for index, cell, what in TARGETS:
            meta = R.cell_meta(index, cell)
            if meta is None or meta.get('status') != 'ok':
                continue
            out.append((index, cell, what))
        if not out:
            self.skipTest('none of the scale cells is cached (run realdata.py fill)')
        return out

    def _m(self, index, cell, what):
        m = measure(index, cell)
        row = _brief(m)
        row['what'] = what
        record('targets', _key(index, cell), row)
        return m

    def test_parse_and_canonical_dump(self):
        for index, cell, what in self._targets():
            m = self._m(index, cell, what)
            body_mb = m['body_bytes'] / MB
            with self.subTest(cell=cell):
                self.assertTrue(m['dump_equal'], 'dump() differs from the body: %s' % cell)
                self.assertEqual(m['graphlet_bytes_server'], m['body_bytes'], cell)
                self.assertLess(m['parse_s'], 1.0 + 3.0 * body_mb, m)
                self.assertLess(m['from_response_s'], 1.5 + 2 * m['parse_s'], m)
                self.assertLess(m['dump_s'], 1.0 + 3.0 * m['parse_s'], m)

    def test_model_memory(self):
        """The retained heap of the parsed model and the parse peak (tracemalloc)."""
        for index, cell, what in self._targets():
            m = self._m(index, cell, what)
            body = m['body_bytes']
            with self.subTest(cell=cell):
                # measured 1.4-17x the body (17x: the constrain hub, 16k runs)
                self.assertLess(m['model_retained_bytes'], 16 * MB + 40 * body, m)
                self.assertLess(m['parse_peak_bytes'],
                                16 * MB + 2 * m['model_retained_bytes'] + 4 * body, m)
                # memory_bytes() is the store's RAM accounting: it must at least be
                # the right order of magnitude (B2 is the strict version)
                ratio = m['memory_bytes'] / m['model_retained_bytes']
                self.assertGreater(ratio, 0.25, m)
                self.assertLess(ratio, 2.0, m)
        ref = _key(*REFERENCE)
        if ref in _MEASURED:
            # design §7: "model 30 MB" for the reference case (memory_bytes); the traced
            # heap is what a process pays
            self.assertLess(_MEASURED[ref]['model_retained_bytes'], 128 * MB)

    # B2 (fixed): GraphletStore's RAM budget is enforced with memory_bytes(), which
    # counted 35-67 % of the traced heap of the model (see the module docstring)
    def test_memory_bytes_tracks_the_traced_heap(self):
        worst = None
        for index, cell, what in self._targets():
            m = self._m(index, cell, what)
            ratio = m['memory_bytes'] / m['model_retained_bytes']
            worst = ratio if worst is None else min(worst, ratio)
        record('findings', 'B2_memory_bytes_worst_ratio', round(worst, 3))
        self.assertGreaterEqual(worst, 0.75, 'memory_bytes() counts %.0f %% of the traced '
                                'heap on the worst target' % (100 * worst))

    def test_to_json(self):
        for index, cell, what in self._targets():
            m = self._m(index, cell, what)
            if m.get('to_json_s') is None:
                continue
            with self.subTest(cell=cell):
                self.assertEqual(m['to_json_paths'], m['walks_n'], cell)
                self.assertLess(m['to_json_s'],
                                2.0 + 0.25 * m['full_json_bytes'] / MB + 3 * m['parse_s'], m)
        ref = _MEASURED.get(_key(*REFERENCE))
        if ref is None:
            self.skipTest('the reference cell is not cached')
        # the reference case's full JSON (72 MB) is above the harness self-check's 40 MB
        # cut: to_json() must equal it here
        self.assertTrue(ref.get('to_json_equal'), ref.get('to_json_diff'))

    def test_walks_and_claims_over_everything(self):
        for index, cell, what in self._targets():
            m = self._m(index, cell, what)
            with self.subTest(cell=cell):
                self.assertTrue(m['walks_leaves_ok'], cell)
                self.assertEqual(sum(m['walks_n'].values()), m['walks'], cell)
                self.assertEqual(m['walks_segments'], m['chain'], cell)
                rec = R.cell_meta(index, cell)['results'][0]['arms']
                self.assertEqual({s: a['leaves'] for s, a in rec.items()}, m['walks_n'], cell)
                self.assertLess(m['walks_s'], 2.0 + 5 * m['parse_s'] + 2e-6 * m['W'], m)
                self.assertLess(m['claims_s'], 1.0 + 2 * m['parse_s'], m)

    def test_label_summary_and_summary(self):
        for index, cell, what in self._targets():
            m = self._m(index, cell, what)
            with self.subTest(cell=cell):
                self.assertTrue(m['label_summary_ok'], cell)
                self.assertLess(m['label_summary_s'], 1.0 + m['parse_s'], m)
                self.assertLessEqual(m['summary_bytes'], 2048, m)
                self.assertLess(m['summary_s'], 1.0 + m['parse_s'], m)

    def test_fasta_and_gfa(self):
        for index, cell, what in self._targets():
            m = self._m(index, cell, what)
            with self.subTest(cell=cell):
                self.assertEqual(m['fasta_records'], m['walks'], cell)
                recs = m['gfa_records']
                self.assertEqual(recs.get('S', 0), m['segments'] + 1, cell)    # + the seed
                self.assertEqual(recs.get('P', 0), m['walks'], cell)
                self.assertEqual(recs.get('L', 0), m['gfa_links_expected'], cell)
                self.assertLess(m['to_fasta_s'], 1.0 + 0.1 * m['fasta_bytes'] / MB, m)
                self.assertLess(m['to_gfa_s'], 1.0 + 0.2 * m['gfa_bytes'] / MB, m)


# ====================================================================== super-linear

def _ladder_points(index, cells):
    pts = []
    for cell in cells:
        meta = R.cell_meta(index, cell)
        if meta is None or meta.get('status') != 'ok':
            continue
        pts.append(measure(index, cell))
    return pts


def _work(m, op):
    if op in ('to_fasta', 'to_gfa'):
        return m['W'] + m[{'to_fasta': 'fasta_bytes', 'to_gfa': 'gfa_bytes'}[op]] / 8
    return m['W']


def _y(m, op):
    if op == 'memory':
        return m['model_retained_bytes']
    if op == 'parse':
        return m['parse_s']
    return m.get(op + '_s')


@slow
class TestSuperLinear(unittest.TestCase):
    """The same operation on retrievals of increasing size: the fitted exponent against
    the work W must stay <= SUPERLINEAR (1 = linear); the exponent against the body is
    reported (the body compresses label sets: W/body grows from 0.2 to 2)."""

    def _ladder(self, name, index, cells):
        pts = _ladder_points(index, cells)
        if len(pts) < 4:
            self.skipTest('fewer than 4 ladder cells cached for %s' % name)
        out = {'cells': [p['cell'] for p in pts],
               'points': [{k: _brief(p).get(k) for k in
                           ('body_bytes', 'W', 'E', 'chain', 'parse_s', 'dump_s', 'to_json_s',
                            'walks_s', 'claims_s', 'label_summary_s', 'summary_s',
                            'to_fasta_s', 'to_gfa_s', 'model_retained_bytes', 'fasta_bytes',
                            'gfa_bytes')} for p in pts],
               'exponents': {}}
        flagged = []
        for op in OPS:
            ys = [_y(p, op) for p in pts]
            vs_w = fit_exponent([_work(p, op) for p in pts], ys,
                                min_y=1 if op == 'memory' else MIN_FIT_T)
            vs_body = fit_exponent([p['body_bytes'] for p in pts], ys,
                                   min_y=1 if op == 'memory' else MIN_FIT_T)
            out['exponents'][op] = {'vs_W': vs_w, 'vs_body': vs_body}
            big = max((y for y in ys if y is not None), default=0)
            meaningful = op == 'memory' or big >= 0.05
            if vs_w is not None and meaningful and vs_w[0] > SUPERLINEAR:
                flagged.append('%s: time ~ W^%.2f (r2 %.2f, %d points)' % ((op,) + vs_w))
        out['flagged'] = flagged
        record('scaling', name, out)
        self.assertEqual([], flagged)

    def test_annotate_beam_ladder(self):
        cells = list(LADDER_BEAM)
        if HEAVY:
            cells.append(HEAVIEST[1])
        self._ladder('sra_annotate_beam20', 'sra', cells)

    def test_constrain_exhaustive_ladder(self):
        self._ladder('sra_constrain_exhaustive', 'sra', LADDER_EXH)


# ====================================================================== compare

def _pairs_by_seed(index):
    """Per seed, the pair (a, b) of its two largest retrievals with DIFFERENT bodies
    (deterministic, single-seed, body <= 8 MB: the 31.75 MB body is left out); sorted by
    the pair's combined body size, largest first."""
    by_seed = {}
    for row in R.cells(index, heavy=None, max_body_bytes=R.HEAVY_BODY_BYTES):
        if row['seed_kind'] == 'batch' or row.get('n_results') != 1:
            continue
        if not row.get('deterministic') or (row.get('n_errors') or 0):
            continue
        by_seed.setdefault(row['seed'], []).append(row)
    out = []
    for rows in by_seed.values():
        rows.sort(key=lambda r: (-r['graphlet']['body_bytes'], r['cell']))
        a = rows[0]
        text_a = R.load_cell(index, a).text(0)
        b = next((r for r in rows[1:] if R.load_cell(index, r).text(0) != text_a), None)
        if b is not None:
            out.append((a['graphlet']['body_bytes'] + b['graphlet']['body_bytes'], a, b))
    out.sort(key=lambda x: (-x[0], x[1]['cell']))
    return [(a, b) for _, a, b in out]


@slow
class TestCompareAtScale(unittest.TestCase):

    def _compare(self, ga, gb, mode):
        return _timed(lambda: ga.compare(gb, mode=mode))

    def _run_pair(self, tag, index, a, b, modes, allowed):
        ca, cb = R.load_cell(index, a), R.load_cell(index, b)
        ga, gb = ca.graphlet(0), cb.graphlet(0)
        row = {'a': _key(index, a), 'b': _key(index, b),
               'a_body': len(ca.text(0)), 'b_body': len(cb.text(0)), 'modes': {}}
        for mode in modes:
            dt, cmp = self._compare(ga, gb, mode)
            row['modes'][mode] = {'s': round(dt, 4), 'comparable': cmp.comparable,
                                  'equal': cmp.equal, 'depth_used': cmp.depth_used,
                                  'only_in_a': len(cmp.only_in_a),
                                  'only_in_b': len(cmp.only_in_b),
                                  'differ': len(cmp.differ)}
            with self.subTest(pair=tag, mode=mode):
                self.assertIn(cmp.comparable, allowed, (cmp.reason, cmp.notes))
                # equality is decided only for a comparable True pair
                self.assertEqual(cmp.comparable is True, cmp.equal is not None, cmp.notes)
                self.assertEqual(cmp.depth_used, min(
                    min(ga.arms[s].complete_to_bp, gb.arms[s].complete_to_bp)
                    for s in ga.arms if s in gb.arms))
        record('compare', tag, row)
        return row

    def test_largest_verified_pair(self):
        """mini_refseq carries index_fp: comparable True, equality decided."""
        R.require_cache('mini_refseq')
        pairs = _pairs_by_seed('mini_refseq')
        if not pairs:
            self.skipTest('no comparable mini_refseq pair cached')
        a, b = pairs[0]
        row = self._run_pair('largest_verified', 'mini_refseq', a['cell'], b['cell'],
                             ('claims', 'labels', 'walks', 'prefix_subset'),
                             (True, 'qualified', 'unknown'))
        for mode, v in row['modes'].items():
            self.assertLess(v['s'], 30.0, (mode, v))

    def test_largest_unverifiable_pair(self):
        """SRA has no index_fp: comparable 'unverifiable', equal None, the lists still
        computed. claims and labels are cheap; walks / prefix_subset are timed (B1)."""
        R.require_cache('sra')
        pairs = _pairs_by_seed('sra')
        if not pairs:
            self.skipTest('no comparable SRA pair cached')
        a, b = pairs[0]
        row = self._run_pair('largest_unverifiable', 'sra', a['cell'], b['cell'],
                             ('claims', 'labels', 'walks', 'prefix_subset'),
                             ('unverifiable', 'unknown'))
        self.assertLess(row['modes']['claims']['s'], 10.0, row)
        self.assertLess(row['modes']['labels']['s'], 10.0, row)
        # a sanity cap only: B1 makes these minutes on other pairs
        self.assertLess(row['modes']['walks']['s'], 300.0, row)
        self.assertLess(row['modes']['prefix_subset']['s'], 300.0, row)

    def test_identical_bodies_compare_without_differences(self):
        """sra_rand50_01 x beam20 and x beam20_10k: the 10,000-path cap trips before 3
        kb, the bodies are byte-identical (different requests)."""
        a, b = 'sra_rand50_01__beam20', 'sra_rand50_01__beam20_10k'
        ca, cb = _cell('sra', a), _cell('sra', b)
        if ca.text(0) != cb.text(0):
            self.skipTest('the two bodies are no longer identical')
        row = self._run_pair('identical_bodies', 'sra', a, b, ('claims', 'labels'),
                             ('unverifiable',))
        for mode, v in row['modes'].items():
            self.assertEqual((0, 0, 0), (v['only_in_a'], v['only_in_b'], v['differ']),
                             (mode, v))
            self.assertLess(v['s'], 10.0, (mode, v))

    # B1 (fixed): compare(mode='walks'|'prefix_subset') built every walk's whole
    # support profile and spelling before cutting at the common depth (see the module
    # docstring)
    def test_compare_walks_at_depth_1_costs_like_claims(self):
        index, a, b = DEPTH1_PAIR
        ga, gb = _cell(index, a).graphlet(0), _cell(index, b).graphlet(0)
        t_claims, cmp = self._compare(ga, gb, 'claims')
        if cmp.depth_used != 1:
            self.skipTest('the pair no longer compares at depth 1 (%s)' % cmp.depth_used)
        t_walks, _ = self._compare(ga, gb, 'walks')
        t_prefix, _ = self._compare(ga, gb, 'prefix_subset')
        bound = max(1.0, 20 * t_claims)
        record('findings', 'B1_depth1_pair', {
            'a': _key(index, a), 'b': _key(index, b), 'depth_used': 1,
            'walks_a': sum(len(derive.paths(x)) for x in ga.arms.values()),
            'claims_s': round(t_claims, 4), 'walks_s': round(t_walks, 3),
            'prefix_subset_s': round(t_prefix, 3), 'bound_s': round(bound, 3)})
        self.assertLess(t_walks, bound)
        self.assertLess(t_prefix, bound)

    # B1 again (fixed), as an exponent: compare(mode='walks') time was ~ walks^2 on the
    # beam ladder (comparison depth <= 190 bp throughout); the work is linear in the walks
    def test_compare_walks_scales_linearly_in_walks(self):
        xs, ys, rows = [], [], []
        for seed in LADDER_CMP:
            a, b = R.cell_id(seed, 'beam20'), R.cell_id(seed, 'beam1')
            if not all((R.cell_meta('sra', c) or {}).get('status') == 'ok' for c in (a, b)):
                continue
            ga, gb = R.load_cell('sra', a).graphlet(0), R.load_cell('sra', b).graphlet(0)
            colds = (_Cold(ga), _Cold(gb))
            n = sum(len(derive.paths(x)) for g in (ga, gb) for x in g.arms.values())
            # since the fix a comparison takes about a millisecond: one cold run is below
            # MIN_FIT_T, so the mean of repeated cold runs is fitted
            dt, cmp, reps = _cold_mean(lambda: ga.compare(gb, mode='walks'), colds)
            xs.append(n)
            ys.append(dt)
            rows.append({'seed': seed, 'walks': n, 'depth_used': cmp.depth_used,
                         'walks_s': round(dt, 5), 'runs': reps})
            del ga, gb, colds
        fit = fit_exponent(xs, ys, min_y=1e-5)
        if fit is None:
            self.skipTest('too few compare ladder points')
        record('findings', 'B1_compare_walks_ladder', {'points': rows, 'exponent_vs_walks': fit})
        self.assertLessEqual(fit[0], SUPERLINEAR, fit)


# ====================================================================== the store

STORE_CELLS = [  # (index, cell): the beam and exhaustive cells the store is driven with
    ('sra', 'sra_hub__exhaustive'), ('sra', 'sra_rand50_00__exhaustive'),
    ('sra', 'sra_16s_PZ326290__exhaustive'), ('sra', 'sra_16s_OK510341__exhaustive'),
    ('sra', 'sra_16s_OK510341__annotate_exh'), ('sra', 'sra_hub__annotate_exh'),
    ('sra', 'sra_hub__beam1'), ('sra', 'sra_rand50_00__beam20'),
    ('sra', 'sra_hub__beam20'), REFERENCE,
]


def _store_cells():
    return [(i, c) for i, c in STORE_CELLS
            if (R.cell_meta(i, c) or {}).get('status') == 'ok']


@slow
class TestStoreAtScale(unittest.TestCase):

    def setUp(self):
        self.dir = _tmpdir('store')

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def test_spool_round_trip(self):
        """put -> a restarted store parses it from the spool -> save -> load: the body
        and the envelope survive unchanged, at every size."""
        cells = _store_cells()
        if not cells:
            self.skipTest('no store cell cached')
        spool = os.path.join(self.dir, 'spool')
        st = GraphletStore(spool)
        handles = {}
        rows = {}
        for index, cell in cells:
            c = R.load_cell(index, cell)
            resp, req = c.graphlet_response, c.request
            t, h = _timed(lambda: st.put(resp, req, source=index))
            handles[(index, cell)] = h
            rows[_key(index, cell)] = {'body_bytes': len(c.text(0)), 'put_s': round(t, 4)}
        del st
        gc.collect()
        st2 = GraphletStore(spool)                     # a restart: nothing in RAM
        for (index, cell), h in handles.items():
            c = R.load_cell(index, cell)
            body = c.text(0)
            row = rows[_key(index, cell)]
            m = measure(index, cell)
            with self.subTest(cell=cell):
                self.assertIn(h, st2)
                self.assertFalse(st2.in_ram(h))
                t, g = _timed(lambda: st2.graphlet(h))
                row['get_cold_s'] = round(t, 4)
                self.assertEqual(body, g.dump(envelope=False))
                self.assertEqual(c.graphlet_result(0).get('seed'),
                                 (g.seed_summary or {}).get('seed'))
                with open(st2._body_path(st2.get(h).digest), 'rb') as f:
                    self.assertEqual(body.encode('utf-8'), f.read())
                path = os.path.join(self.dir, cell + '.mgt')
                t, n = _timed(lambda: st2.save(h, path))
                row['save_s'] = round(t, 4)
                row['mgt_file_bytes'] = n
                t, h2 = _timed(lambda: st2.load(path))
                row['load_s'] = round(t, 4)
                g2 = st2.graphlet(h2)
                self.assertEqual(body, g2.dump(envelope=False))
                self.assertEqual(g.envelope, g2.envelope)
                self.assertEqual(st2.get(h).digest, st2.get(h2).digest)   # one body file
                self.assertEqual(st2.request_of(g2)['seeds'][0]['sequence'],
                                 c.seed_sequence(0).upper())
                bound = 2.0 + 3 * m['parse_s']
                for k in ('put_s', 'get_cold_s', 'load_s'):
                    self.assertLess(row[k], bound + m['dump_s'], (k, row))
                self.assertLess(row['save_s'], 2.0 + 3 * m['dump_s'], row)
                st2.free(h2)
        bodies = os.listdir(os.path.join(spool, 'bodies'))
        self.assertEqual(len({st2.get(h).digest for h in handles.values()}), len(bodies))
        record('store', 'spool_round_trip', rows)

    def test_identical_bodies_share_one_spool_file(self):
        a, b = ('sra', 'sra_rand50_01__beam20'), ('sra', 'sra_rand50_01__beam20_10k')
        ca, cb = _cell(*a), _cell(*b)
        if ca.text(0) != cb.text(0):
            self.skipTest('the two bodies are no longer identical')
        st = GraphletStore(os.path.join(self.dir, 'spool'), max_ram_mb=1)
        ha = st.put(ca.graphlet_response, ca.request)
        hb = st.put(cb.graphlet_response, cb.request)
        self.assertNotEqual(ha, hb)                       # entry identity != body identity
        self.assertEqual(st.get(ha).digest, st.get(hb).digest)
        self.assertNotEqual(st.get(ha).request, st.get(hb).request)
        bodies = os.path.join(self.dir, 'spool', 'bodies')
        self.assertEqual(1, len(os.listdir(bodies)))
        st.free(ha)
        self.assertEqual(1, len(os.listdir(bodies)))      # still held by hb
        self.assertEqual(ca.text(0), st.graphlet(hb).dump(envelope=False))
        st.free(hb)
        self.assertEqual(0, len(os.listdir(bodies)))

    def test_small_store_evicts_lru_and_never_holds_an_oversize_model(self):
        cells = _store_cells()
        if len(cells) < 4:
            self.skipTest('too few store cells cached')
        budget_mb = 16
        st = GraphletStore(os.path.join(self.dir, 'spool'), max_ram_mb=budget_mb,
                           max_handles=64)
        budget = st.max_ram_bytes
        sizes, order, handles, rows = {}, [], {}, {}
        for index, cell in cells:
            c = R.load_cell(index, cell)
            t, h = _timed(lambda: st.put(c.graphlet_response, c.request, source=index))
            handles[h] = (index, cell)
            # what the store charged (memory_bytes() of its own model); a model never
            # admitted is charged nothing, its size is the same estimate
            sizes[h] = st._ram[h][1] if st.in_ram(h) else measure(index, cell)['memory_bytes']
            order.append(h)
            # the LRU model: the longest recency suffix of the admissible handles that
            # fits the budget (an oversize model is never admitted)
            expect, total = [], 0
            for x in reversed(order):
                if sizes[x] > budget:
                    continue
                if total + sizes[x] > budget:
                    break
                expect.append(x)
                total += sizes[x]
            with self.subTest(cell=cell):
                self.assertEqual(sorted(expect), sorted(x for x in order if st.in_ram(x)))
                self.assertLessEqual(st._ram_bytes, budget)
                self.assertEqual(sizes[h] <= budget, st.in_ram(h))
            rows[_key(index, cell)] = {'put_s': round(t, 4), 'memory_bytes': sizes[h],
                                       'resident_after_put': st.in_ram(h)}
        # touching the oldest resident model makes it the most recent: the next put
        # evicts the one after it instead
        resident = [x for x in order if st.in_ram(x)]
        if len(resident) >= 2:
            oldest, second = resident[0], resident[1]
            st.graphlet(oldest)
            order.remove(oldest)
            order.append(oldest)
            index, cell = handles[second]
            c = R.load_cell(index, cell)
            h = st.put(c.graphlet_response, c.request)
            handles[h] = (index, cell)
            sizes[h] = sizes[second]
            order.append(h)
            if sizes[h] + sizes[oldest] <= budget:
                self.assertTrue(st.in_ram(oldest))
        # every evicted or oversize handle answers from the spool, identical
        for h in order:
            if st.in_ram(h):
                continue
            index, cell = handles[h]
            t, g = _timed(lambda: st.graphlet(h))
            self.assertEqual(R.load_cell(index, cell).text(0), g.dump(envelope=False))
            rows.setdefault(_key(index, cell), {})['reparse_s'] = round(t, 4)
            if sizes[h] > budget:
                self.assertFalse(st.in_ram(h), cell)       # parsed on demand, not kept
        self.assertLessEqual(st._ram_bytes, budget)
        record('store', 'small_store_%dmb' % budget_mb, rows)

    # B2 (fixed): the store's resident models took 2-3x max_ram_mb of heap, because the
    # budget was charged with the shallow memory_bytes() (see the module docstring)
    def test_small_store_heap_stays_within_max_ram_mb(self):
        # the constrain tries (most runs per byte: the worst-accounted models), smallest
        # first, so that the store ends up holding the largest ones that fit
        cells = sorted((x for x in _store_cells() if measure(*x)['mode'] == 'constrain'
                        and measure(*x)['memory_bytes'] < 8 * MB),
                       key=lambda x: measure(*x)['memory_bytes'])
        if len(cells) < 2:
            self.skipTest('too few constrain store cells cached')
        budget_mb = 8
        payload = [(R.load_cell(i, c).graphlet_response, R.load_cell(i, c).request)
                   for i, c in cells]
        gc.collect()
        tracemalloc.start()
        try:
            base = tracemalloc.get_traced_memory()[0]
            st = GraphletStore(os.path.join(self.dir, 'spool'), max_ram_mb=budget_mb)
            for resp, req in payload:
                st.put(resp, req)
            gc.collect()
            held = tracemalloc.get_traced_memory()[0] - base
        finally:
            tracemalloc.stop()
        record('findings', 'B2_small_store_heap', {
            'cells': [_key(*x) for x in cells],
            'max_ram_mb': budget_mb, 'accounted_bytes': st._ram_bytes,
            'traced_heap_bytes': held, 'resident': sum(st.in_ram(h) for h in
                                                       [e['handle'] for e in st.list()])})
        self.assertLessEqual(st._ram_bytes, st.max_ram_bytes)
        self.assertLessEqual(held, 1.25 * st.max_ram_bytes + 2 * MB,
                             'a %d MB store holds %.1f MB of heap' % (budget_mb, held / MB))


# ====================================================================== MCP sessions

@slow
class TestMcpAtScale(unittest.TestCase):
    """A cold MCP session over big retrievals: a fresh store and GraphletTools over a
    spool an earlier session wrote, and spooled delivery through traverse_fetch."""

    def setUp(self):
        self.dir = _tmpdir('mcp')

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def _check(self, out, tool, ceiling=DEFAULT_MAX_BYTES, allow_error=()):
        self.assertLessEqual(_size(out), ceiling, (tool, _size(out)))
        if 'error' in out:
            self.assertIn(out['error'], allow_error, (tool, out))
        return out

    def test_cold_session_over_a_written_spool(self):
        index, cell = REFERENCE
        c = _cell(index, cell)
        other = R.cell_id(c.seed_name, 'beam20')
        spool = os.path.join(self.dir, 'spool')
        st = GraphletStore(spool)
        h = st.put(c.graphlet_response, c.request, source=index)
        hb = None
        if (R.cell_meta(index, other) or {}).get('status') == 'ok':
            cb = R.load_cell(index, other)
            hb = st.put(cb.graphlet_response, cb.request, source=index)
        del st
        gc.collect()
        m = measure(index, cell)
        times = {}

        def call(name, fn, ceiling=DEFAULT_MAX_BYTES, allow_error=()):
            t, out = _timed(fn)
            times[name] = round(t, 4)
            return self._check(out, name, ceiling, allow_error)

        t0 = time.perf_counter()
        store = GraphletStore(spool)
        tools = GraphletTools(store, {index: ReplayClient(index, [cell])})
        times['session_start'] = round(time.perf_counter() - t0, 4)
        s = call('graphlet_summary (cold)', lambda: tools.graphlet_summary(h))
        self.assertEqual(s['summary']['outcome']['walks'], c.graphlet_result(0)['outcome']['walks'])
        call('graphlet_summary (warm)', lambda: tools.graphlet_summary(h))
        call('graphlet_summary limitations', lambda: tools.graphlet_summary(
            h, detail='limitations'))
        g = store.graphlet(h)
        for side in g.arms:
            w = call('graphlet_walks %s p1' % side, lambda: tools.graphlet_walks(h, side))
            self.assertEqual(len(derive.paths(g.arms[side])), w['total'])
            self.assertGreaterEqual(len(w['rows']), 1)
            if w.get('next_cursor'):
                call('graphlet_walks %s p2' % side, lambda: tools.graphlet_walks(
                    h, side, cursor=w['next_cursor']))
            top = w['rows'][0]['walk']
            call('graphlet_walk %s' % side, lambda: tools.graphlet_walk(h, side, top))
            call('graphlet_support %s' % side, lambda: tools.graphlet_support(h, side, top))
            call('graphlet_claims %s' % side, lambda: tools.graphlet_claims(h, side))
            call('graphlet_splits %s' % side, lambda: tools.graphlet_splits(h, side))
            sq = call('graphlet_sequence %s' % side, lambda: tools.graphlet_sequence(
                h, side, walk=top, with_seed=True), SEQUENCE_MAX_BYTES)
            self.assertEqual(g.spell(side, top, with_seed=True)[:len(sq['sequence'])],
                             sq['sequence'])
        call('graphlet_labels', lambda: tools.graphlet_labels(h))
        call('graphlet_labels at_bp', lambda: tools.graphlet_labels(h, at_bp=10))
        if hb is not None:
            cmp = call('graphlet_compare claims', lambda: tools.graphlet_compare(h, hb))
            self.assertEqual('unverifiable' if index != 'mini_refseq' else True,
                             cmp['comparable'])
        files = {}
        for fmt in ('fasta', 'gfa', 'json', 'mgt'):
            out = call('graphlet_export %s' % fmt, lambda: tools.graphlet_export(h, fmt))
            files[fmt] = out['bytes']
            self.assertTrue(os.path.exists(out['path']))
            self.assertEqual(os.path.getsize(out['path']), out['bytes'])
        self.assertEqual(files['fasta'], m['fasta_bytes'])
        self.assertEqual(files['gfa'], m['gfa_bytes'])
        ld = call('graphlet_load', lambda: tools.graphlet_load('%s.mgt' % h))
        self.assertEqual(c.text(0), store.graphlet(ld['handle']).dump(envelope=False))
        record('mcp', 'cold_session', {'cell': _key(index, cell), 'times_s': times,
                                       'export_bytes': files, 'parse_s': m['parse_s']})
        # bounds: the first call parses from the spool, the others are local
        self.assertLess(times['graphlet_summary (cold)'], 2.0 + 3 * m['parse_s'], times)
        self.assertLess(times['graphlet_summary (warm)'], 0.5, times)
        for name, t in times.items():
            if name.startswith(('graphlet_walks', 'graphlet_compare')):
                self.assertLess(t, 5.0 + 3 * m['walks_s'], (name, times))
            elif name.startswith('graphlet_export'):
                self.assertLess(t, 5.0 + 5 * max(m['to_gfa_s'], m['to_json_s'] or 0),
                                (name, times))
            elif name in ('graphlet_load',):
                self.assertLess(t, 2.0 + 3 * m['parse_s'], (name, times))
            elif not name.startswith(('graphlet_summary (cold)', 'session_start')):
                self.assertLess(t, 2.0, (name, times))

    def test_spooled_delivery_and_eviction_through_traverse_fetch(self):
        """The beam / exhaustive cells fetched through the tool layer (replayed from the
        cache) into a deliberately small store: a body above max_graphlet_mb is spooled
        (delivery 'spooled' in the tool's answer, the summary and the stored model), a
        model above max_ram_mb is never resident, and every handle still answers."""
        cells = _store_cells()
        if not cells:
            self.skipTest('no store cell cached')
        by_index = {}
        for i, c in cells:
            by_index.setdefault(i, []).append(c)
        store = GraphletStore(os.path.join(self.dir, 'spool'), max_ram_mb=8, max_handles=4)
        tools = GraphletTools(store, {i: ReplayClient(i, cs) for i, cs in by_index.items()})
        threshold_mb = 0.25
        rows = {}
        handles = []
        for index, cell in cells:
            c = R.load_cell(index, cell)
            seed = {'sequence': c.seed_sequence(0)}
            t, out = _timed(lambda: tools.traverse_fetch(
                index=index, seed=seed, strategy=_strategy_of(c),
                max_graphlet_mb=threshold_mb))
            self._check(out, 'traverse_fetch')
            body = len(c.text(0))
            spooled = body > threshold_mb * MB
            h = out['handle']
            handles.append((h, index, cell))
            with self.subTest(cell=cell):
                self.assertEqual('spooled' if spooled else 'inline', out['delivery'])
                self.assertEqual(out['delivery'], out['summary']['outcome']['delivery'])
                self.assertEqual(out['delivery'], out['evidence']['outcome']['delivery'])
                if spooled:
                    self.assertFalse(store.in_ram(h))
                    self.assertEqual('spooled', store.get(h).delivery)
                self.assertLessEqual(store._ram_bytes, store.max_ram_bytes)
                self.assertLessEqual(len(store._ram), 4)
            rows[_key(index, cell)] = {'body_bytes': body, 'delivery': out['delivery'],
                                       'fetch_s': round(t, 4)}
        # every handle answers (re-parsed from the spool where evicted / oversize), with
        # the delivery that actually happened
        for h, index, cell in handles:
            c = R.load_cell(index, cell)
            t, s = _timed(lambda: tools.graphlet_summary(h))
            self._check(s, 'graphlet_summary')
            t2, w = _timed(lambda: tools.graphlet_walks(h, 'right'))
            self._check(w, 'graphlet_walks', allow_error=('bad_arm', 'not_in_view'))
            row = rows[_key(index, cell)]
            row['summary_after_s'] = round(t, 4)
            row['walks_after_s'] = round(t2, 4)
            row['resident_at_end'] = store.in_ram(h)
            self.assertEqual(row['delivery'], s['evidence']['outcome']['delivery'], cell)
            # the spool holds the body as received; the model says how it was delivered
            self.assertEqual(c.text(0), _spooled_body(store, h), cell)
            self.assertEqual(_as_delivered(c.text(0), row['delivery']),
                             store.graphlet(h).dump(envelope=False), cell)
            m = measure(index, cell)
            self.assertLess(t, 2.0 + 3 * m['parse_s'], (cell, row))
        self.assertLessEqual(store._ram_bytes, store.max_ram_bytes)
        record('mcp', 'small_store_fetch', {'max_ram_mb': 8, 'max_handles': 4,
                                            'max_graphlet_mb': threshold_mb, 'cells': rows})

    # B3 (fixed): graphlet_walks rebuilt and spelled the whole ranked walk list on every
    # page (see the module docstring)
    def test_walks_page_costs_a_page_not_the_arm(self):
        index, cell = REFERENCE
        c = _cell(index, cell)
        store = GraphletStore(os.path.join(self.dir, 'spool'))
        h = store.put(c.graphlet_response, c.request)
        tools = GraphletTools(store)
        g = store.graphlet(h)
        side = max(g.arms, key=lambda s: len(derive.paths(g.arms[s])))
        g.walks(side)                                         # warm the derivations
        t_all, ws = _timed(lambda: g.walks(side))
        t_top, _ = _timed(lambda: g.walks(side, top=10))
        tools.graphlet_walks(h, side)
        t_p1, p1 = _timed(lambda: tools.graphlet_walks(h, side))
        t_p2, p2 = _timed(lambda: tools.graphlet_walks(h, side, cursor=p1['next_cursor']))
        per_page = max(1, len(p1['rows']))
        projected = p1['total'] / per_page * (t_p1 + t_p2) / 2
        record('findings', 'B3_walks_paging', {
            'cell': _key(index, cell), 'arm': side, 'walks': len(ws),
            'rows_per_page': per_page, 'walks_all_s': round(t_all, 4),
            'walks_top10_s': round(t_top, 4), 'page1_s': round(t_p1, 4),
            'page2_s': round(t_p2, 4), 'projected_full_listing_s': round(projected, 1)})
        self.assertLessEqual(_size(p1), DEFAULT_MAX_BYTES)
        self.assertLess(max(t_p1, t_p2), 0.25 * t_all,
                        'a page costs %.2f s, walks() over the whole arm %.2f s; listing '
                        'all %d walks would take ~%.0f s' % (t_p1, t_all, len(ws), projected))

    @R.skip_unless_server('sra')
    def test_live_cold_session_fetch_and_continue(self):
        """The reference case fetched live through the tool layer (a fresh store, a real
        client): spooled above max_graphlet_mb, the body byte-identical to the cache
        (the strategy is deterministic), then one continuation of the top walk."""
        index, cell = REFERENCE
        c = _cell(index, cell)
        if not R.live_matches_cache(index):
            self.skipTest('the live %s server is not the index the cache came from' % index)
        hp = R.server(index)
        store = GraphletStore(os.path.join(self.dir, 'spool'), max_ram_mb=16)
        t0 = time.perf_counter()
        tools = GraphletTools(store, {index: TraverseClient(hp[0], hp[1], timeout=300)})
        times = {'session_start': time.perf_counter() - t0}
        t, caps = _timed(lambda: tools.traverse_capabilities(index))
        times['traverse_capabilities'] = t
        self._check(caps, 'traverse_capabilities')
        t, out = _timed(lambda: tools.traverse_fetch(
            index=index, seed={'sequence': c.seed_sequence(0)}, strategy=_strategy_of(c),
            max_graphlet_mb=1))
        times['traverse_fetch'] = t
        self._check(out, 'traverse_fetch')
        self.assertEqual('spooled', out['delivery'])
        h = out['handle']
        self.assertFalse(store.in_ram(h))
        self.assertEqual(c.text(0), _spooled_body(store, h))     # deterministic strategy
        self.assertEqual(_as_delivered(c.text(0), 'spooled'),
                         store.graphlet(h).dump(envelope=False))
        w = self._check(tools.graphlet_walks(h, 'right', n=1), 'graphlet_walks')
        top = w['rows'][0]['walk']
        t, cont = _timed(lambda: tools.traverse_continue(
            h, 'right', top, overrides={'bounds': {'max_extension_bp': 200}}))
        times['traverse_continue'] = t
        self._check(cont, 'traverse_continue', allow_error=('seed_failed',))
        if 'error' not in cont:
            self.assertEqual({'handle': h, 'arm': 'right', 'walk': top},
                             {k: cont['parent'][k] for k in ('handle', 'arm', 'walk')})
            self.assertEqual(cont['parent'], store.get(cont['handle']).parent)
        record('mcp', 'live_session', {'cell': _key(index, cell),
                                       'times_s': {k: round(v, 3) for k, v in times.items()},
                                       'continuation': cont.get('error') or 'ok'})
        self.assertLess(times['traverse_fetch'], 120.0, times)


if __name__ == '__main__':
    unittest.main()
