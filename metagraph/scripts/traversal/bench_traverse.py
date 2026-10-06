#!/usr/bin/env python3
"""bench_traverse.py -- a repeatable, paced benchmark of a MetaGraph /traverse + /resolve server.

It re-measures, with as few requests as possible, the effects found on the refseq33m staging
server (decode dominates; first touch vs warm; row width; the row-diff path re-decoded per fetch
call; time-budget overruns in the derivation / depth 0 / lookahead; /resolve discover costs), so
that a run before an efficiency fix and a run after it can be compared request by request.
Stdlib only, Python >= 3.10.

USAGE
  # plan only (no network)
  bench_traverse.py --base URL --panel bench_panel_refseq33m.json --dry-run
  # a full run (54 requests on the refseq33m panel; hard cap 70)
  bench_traverse.py --base https://metagraph.ethz.ch:8080/api/refseq33m-experimental \\
      --panel bench_panel_refseq33m.json --label before-pathreuse --window-offset 0
  # after the fix: the SAME windows as the baseline (same K), hours later so they are cold again
  bench_traverse.py ... --label after-pathreuse --window-offset 2500
  # optional extra: a true first touch on fresh windows (timings confounded by locus width, see WINDOW OFFSET)
  bench_traverse.py ... --label after-pathreuse-fresh --window-offset auto --avoid OUT_BEFORE/results.json
  # before/after table (exit code 1 when a deterministic result changed)
  bench_traverse.py --compare OUT_BEFORE/results.json OUT_AFTER/results.json
  # rebuild results.json + summary.md from a run's raw/ (after changing the analysis)
  bench_traverse.py --reanalyze OUT_DIR
  # record coordinates (feature level 6): the trace walks ask for them (cap 16, or --coordinate-cap)
  bench_traverse.py ... --label level6-coords --coordinates
  # then: elapsed with vs without them on the same windows (the walks' content is the same, so the
  # pairs compare as identical completed work), and the block's size per seed
  bench_traverse.py --compare OUT_OPTOUT/results.json OUT_COORDS/results.json
  # offline (no server): the block size a run's trace walks WOULD carry, from its recorded responses
  bench_traverse.py --coordinates-estimate OUT_DIR [OUT_DIR ...]

THE REQUEST SET (fixed; ids are stable across runs, PLAN_VERSION changes when it does)
  caps      GET /traverse/capabilities (every run starts with it: limits, feature_level).
  budgets   time-budget adherence on FRESH windows (panel roles), each sent twice (first touch, then the
            identical warm repeat), batch_kmers pinned:
              budgets.derived_1s.{first,warm}       derived column labels, r1000, 1 s budget, batch 64
              budgets.named_hdr_1s.{first,warm}     the record header named, r1000, 1 s, batch 64
              budgets.lookahead_250ms.{first,warm}  derived columns, r3000, 250 ms, batch_kmers 1000
  walks     per walk seed (one per width class: sparse, narrow x2, medium, wide, widest), at the seed's radius,
            5 s budget, each sent twice in a row (first touch, then warm repeat):
              walks.<seed>.named.{first,warm}       labels = [record header], support kmer
              walks.<seed>.trace.{first,warm}       the same header, support trace, on_reconverge keep, on the
                                                    seed's TRACE window (offset + trace_delta, so its first touch
                                                    is not warmed by the named walk)
  callsize  per call-size seed (wide): /resolve, explicit headers (up to 3), support trace, over 30, 1, 3, 30
            consecutive k-mers of the seed window: callsize.<seed>.{k30_first,k1,k3,k30_warm}. A call costs about
            one row-diff path decode plus a per-row increment: k1 ~ per-call cost, (k30-k1)/29 ~ per extra row.
  resolve   per resolve seed: resolve.<seed>.{explicit_first,explicit_warm,discover_column,discover_header}
            (explicit = up to 3 column labels, support kmer; discover max_labels 5): the discover overhead.
  annotate  annotate.<seed>.{first,warm}: labels mode annotate (columns), max_labels_per_node 64,
            max_live_paths 64, r1000, 5 s.
  lookahead lookahead.<seed>.{b64_first,b64_warm,b1,b1000}: derived columns, r300, 5 s, batch_kmers 64/1/1000
            (rows fetched per base consumed = the lookahead over-read; results must be identical).
  batch     batch.<n>seeds: the panel's batch seeds (<= max_seeds) in one request, named header, r1000, 5 s; each
            seed's graphlet must equal walks.<seed>.named (body without the H line, which holds the seed index).
  smoke     (not in "all") caps + one explicit /resolve + one named walk r100 on the smoke seed's own window
            (offset + suites.smoke.delta, away from every measured window): 3 requests.
  --coordinates  every support-trace walk also sets output.coordinates true (and
            max_coordinate_occurrences, default 16: --coordinate-cap N | unlimited); the request ids stay the
            same, so --compare pairs a run with an opt-out run of the same windows (see COMPARE).
  On the refseq33m panel "all" is 54 requests. The order puts the fresh-window suites (budgets, walks) first;
  callsize / resolve / annotate / lookahead / batch run on windows the walks already read (warm by design;
  run alone, their *_first request is the first touch).

RULES FOR SHARED SERVERS (non-local --base)
  strictly serial (an exclusive lock per server on this machine), >= --pause-s (min 2 s) between requests,
  >= 6 s after a request over 3 s, >= 10 s after one over 20 s, the pause state kept across runs; a hard
  cap (--max-requests, failed requests count); stops after 3 consecutive 5xx / connection errors, pauses
  60 s after a request over 60 s and stops after a second one; gzip responses. A local server
  (127.0.0.1 / localhost) may use --pause-s 0.

WINDOW OFFSET
  Every seed is a window of a longer source sequence in the panel (source + offset). --window-offset K moves
  every window to offset + K (wrapped inside its source), so FIRST-TOUCH requests read rows the previous
  run did not. Results of different K are different walks: --compare then compares timings only and marks
  equality "request differs". results.json records the source intervals each run read (window +- reach +
  margin); --window-offset auto --avoid PREV/results.json [...] picks the smallest K (step 100 bp) whose
  first-touch windows (+- 100 bp) miss all of them. Rows also go cold by themselves (on staging within
  ~20 min once other loci were read), so the same K hours later is usually cold too, but not guaranteed.
  COMPARE BEFORE/AFTER AT THE SAME K: moving the windows moves them onto loci of another width (the ndm1
  window is narrow at K=0 but 65-column wide at K=2500), which confounds every timing; a fresh K is only
  for an extra first-touch run. The refseq33m staging baseline (2026-10-04, release refseq97-8a98759a)
  used K=2500; it is kept in ~/.cache/metagraph-traverse-bench/staging-baseline. If a first-touch time
  comes out close to its warm repeat, the window was not cold.

OUTPUT (--out DIR, default bench_out/<label>_<timestamp>)
  results.json   run meta, capabilities, one record per request (wall time, HTTP status, bytes, server timing:
                 elapsed_ms, annotation_fetch_ms, rows_fetched, tuple_rows_fetched, coords_mapped, cache_hits,
                 serialize_ms and any other number the server adds, a nested object's as "name.field" (e.g.
                 path_cache.hits); per arm status / complete_to_bp / reach / leaves_by_reason / counters; outcome,
                 limitations; release / algorithm_version / feature_level; an attempt's usage (work_units,
                 peak_admitted_bytes) when the request was ledger-managed; derived: ms per fetched row, budget
                 ratio, bp reached, bp per second, rows per base; result hashes), suite tables, invariants, read
                 intervals. Per seed also:
                   reach     CERTIFIED reach: bp (complete_to_bp, both arms), walked_bp (counts.max_bp), steps (the
                             walker's steps: nodes walked), segments, leaves, complete_walks (leaves whose end reason
                             is the walk's own, not a resource stop: max_steps, max_live_paths, max_paths,
                             max_output_bp, time_budget, beam_pruned, resource_limit) and stopped_walks, stop (failed,
                             resource:<memory|work>, the arms' cap reasons, else the walks outcome)
                   work      CONSUMED work: work_units (the arms' counter, reported under a memory or work budget),
                             usage_work_units (the attempt's usage), and the logical counters every server reports --
                             steps, successor_enumerations, annotation rows_requested and keys_mapped; rows_decoded
                             (rows_fetched + tuple_rows_fetched) is physical
                   physical  the row-diff path cache: reuse (timing.path_cache.hits, paths cut at a cached row),
                             stored_rows_read, peak_cache_bytes (timing.path_cache.peak_bytes); None on a server that
                             reports none, 0 when a server that reports timing.seed_phase_ms omits the object (it is
                             omitted when all its counters are 0)
                 coordinates   (a request with output.coordinates) the block's kind, compact JSON bytes and its
                               share of the seed result, entries (seed labels, runs), occurrences, lists cut,
                               lower-bound runs, the largest list, complete; or the null reason
  summary.md     tables by suite, the invariant checks, certified reach and consumed work per seed, record
                 coordinates (when asked for), all requests
  raw/NNN_<id>.json.gz   {meta, request, response} per request
  progress.jsonl one line per request as it completes

COMPARE (--compare A B; exit code 1 when a deterministic result changed)
  Each request (each seed of a batch) present in both runs is one pair, classified:
    completed         both ran to their end. LATENCY is compared, and only on IDENTICAL COMPLETED WORK: the same
                      result (walks_sha / graphlet_sha; /resolve: the result digest) and the same logical work
                      counters (steps, successor_enumerations, rows_requested, keys_mapped; /resolve: num_kmers).
                      A pair whose result or logical work differs is shown in parentheses and left out of the
                      suite's latency; equal results with different logical work are flagged (the server's
                      counting changed).
    deadline-limited  a time budget, deadline or memory / work budget stopped either side. Their elapsed is about
                      the budget on both sides and hides a gain, so the compare shows CERTIFIED REACH (A/B bp,
                      bp B/A, steps, complete walks / leaves, the stop) and CONSUMED WORK (work units where the
                      server reports them, else the annotation rows requested; B/A only in one unit).
    not compared      an HTTP error on a side, a failed seed outside a deadline, a request that differs.
  A pair whose requests differ only in output.coordinates / max_coordinate_occurrences (a --coordinates run
  against an opt-out one) is the same request: its walks are the same (the MGT body differs at most by the
  coordinates limitation's K record, which walks_sha leaves out), so its latency is the cost of the
  coordinates; a "record coordinates" table then lists the block bytes and share per pair.
  Both tables also show the physical path-cache reuse and peak cache bytes (A and B) where the servers report
  them. Per suite: the geometric mean of B/A latency over the identical completed pairs (< 1 = B faster) and of
  B/A certified bp over the deadline-limited pairs (> 1 = B reached further), how many reached further / the
  same / less, and the geometric mean of B/A consumed work. A results.json written before script version 2
  (no counters, leaves_by_reason or physical fields: the staging baseline, mini_before) is analysed again in
  memory from the raw/ directory next to it; without raw/ its stored fields are used (steps and complete walks
  then show "-"). Nothing is written into either run.

RESULT IDENTITY
  graphlet_sha (whole body), graphlet_sha_noH (without the H line; equal between a batch and its singles),
  walks_sha (only the walk-content lines S X L V G P E T C R; equal when the walks are, even if counters or
  limitations changed). /resolve: the response without timing / capabilities / release. A seed result is
  deterministic when no time budget stopped it (no bounds.time_budget_ms limitation, no time_budget cap
  trigger, no time resource stop); only deterministic results are expected to repeat byte-identically.

COORDINATES ESTIMATE (--coordinates-estimate RUN_DIR ...; offline, no server)
  For every support-trace seed result recorded in the runs' raw/ responses (a run without --coordinates), the
  coordinates block it would carry, built as the server writes it (keys sorted, compact JSON): one entry per
  seed label and one per run of each arm (the R records), each with ONE occurrence whose positions have
  --coordinate-digits digits (default 7: positions within a bacterial record; a column label's global positions
  have more). Occurrences per entry are not in an opt-out response, so the estimate is a LOWER bound where a
  record holds the walked sequence several times (rRNA operons, insertion sequences). Prints a table and writes
  coordinates_estimate.json beside each run's results.json.

COMPATIBILITY
  Works with servers without feature_level / attempts (8a98759a) and with newer ones (0a880475+): every field
  is optional; unknown timing fields are recorded as they come. The path cache's statistics
  (timing.path_cache, and timing.seed_phase_ms / seed_fetch_ms / label_resolve_ms beside it) come with the review
  fixes of pass 5; a renamed reuse or peak field is still found by name (path_cache|dependency ... hits|reuse,
  ... peak ... bytes). Results of script version 1 are read by --compare (see COMPARE).
"""

import argparse
import contextlib
import datetime
import glob
import gzip
import hashlib
import json
import math
import os
import platform
import re
import statistics
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request

try:
    import fcntl
except ImportError:  # pragma: no cover (non-POSIX)
    fcntl = None

# 3: record coordinates (--coordinates, the per-seed coordinates facts, --coordinates-estimate); a
# results.json of version 2 is read as it is (no coordinates facts)
SCRIPT_VERSION = 3
PLAN_VERSION = 1
# the coordinates block of a trace walk with --coordinates (None: opt-out, the requests as they were)
COORDINATES_CAP_DEFAULT = 16
ALL_SUITES = ['caps', 'budgets', 'walks', 'callsize', 'resolve', 'annotate', 'lookahead', 'batch']
EXTRA_SUITES = ['smoke']
WALK_BUDGET_MS = 5000
DEFAULT_BATCH_KMERS = 64
TIMING_FIELDS = ('elapsed_ms', 'annotation_fetch_ms', 'rows_fetched', 'tuple_rows_fetched',
                 'coords_mapped', 'cache_hits', 'serialize_ms')
WALK_LINES = set('SXLVGPETCR')
LOCAL_HOSTS = {'127.0.0.1', 'localhost', '::1'}
# leaf end reasons that stop an exploration short of its requested domain (the walker's resource
# stops, _codec.RESOURCE_CODES): a leaf with any other reason is a COMPLETE walk
RESOURCE_LEAF_REASONS = frozenset(('max_steps', 'max_live_paths', 'max_paths', 'max_output_bp', 'time_budget',
                                   'beam_pruned', 'resource_limit'))
# the logical (deterministic) work counters of a seed: equal for equal results on any server build of
# one algorithm, whatever its caches did; a completed pair is "identical completed work" only when they are
COMPLETED_WORK_KEYS = ('steps', 'successor_enumerations', 'rows_requested', 'keys_mapped')
# the row-diff path cache's physical statistics: timing.path_cache {hits, stored_rows_read, rows_kept,
# bytes_kept, peak_bytes} per seed (omitted when all are 0). Older servers have none; a renamed field is
# still found by these patterns over the flattened timing names
PATH_CACHE_REUSE = re.compile(r'^path_cache\.hits$|(path_cache|dependenc)\w*[._](hits|reuse)')
PATH_CACHE_PEAK = re.compile(r'^path_cache\.peak_bytes$|(path_cache|dependenc)\w*[._]peak\w*bytes')


# ============================================================================ small helpers

def now_iso():
    return datetime.datetime.now().astimezone().isoformat(timespec='milliseconds')


def sha(obj):
    if isinstance(obj, str):
        b = obj.encode('utf-8')
    elif isinstance(obj, (bytes, bytearray)):
        b = bytes(obj)
    else:
        b = json.dumps(obj, sort_keys=True, separators=(',', ':')).encode('utf-8')
    return hashlib.sha256(b).hexdigest()[:16]


def r1(x, nd=1):
    if x is None:
        return None
    try:
        return round(float(x), nd)
    except (TypeError, ValueError):
        return None


def div(a, b):
    try:
        if a is None or b is None or b == 0:
            return None
        return a / b
    except TypeError:
        return None


def fmt(x, nd=1, dash='-'):
    if x is None:
        return dash
    if isinstance(x, bool):
        return 'yes' if x else 'no'
    if isinstance(x, int):
        return str(x)
    if isinstance(x, float):
        if math.isnan(x):
            return dash
        if abs(x) >= 1000:
            return '%.0f' % x
        if 0 < abs(x) < 0.1:
            return '%.2g' % x      # small values (ms per row on a small index): two significant digits
        return ('%.' + str(nd) + 'f') % x
    return str(x)


def md_table(header, rows):
    out = ['| ' + ' | '.join(header) + ' |', '|' + '|'.join('---' for _ in header) + '|']
    for r in rows:
        out.append('| ' + ' | '.join(str(c).replace('|', '/') for c in r) + ' |')
    return '\n'.join(out)


def write_json_atomic(path, obj):
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(obj, f, indent=1, default=str)
        f.write('\n')
    os.replace(tmp, path)


def is_local(base):
    host = urllib.parse.urlparse(base).hostname or ''
    return host in LOCAL_HOSTS


def geomean(xs):
    xs = [x for x in xs if x and x > 0]
    if not xs:
        return None
    return math.exp(sum(math.log(x) for x in xs) / len(xs))


# ============================================================================ panel

class Panel:
    def __init__(self, path):
        with open(path, 'rb') as f:
            raw = f.read()
        d = json.loads(raw)
        if d.get('format') != 'metagraph-bench-traverse-panel':
            raise SystemExit('%s: not a bench panel (format %r)' % (path, d.get('format')))
        self.path = os.path.abspath(path)
        self.sha = hashlib.sha256(raw).hexdigest()[:16]
        self.meta = d.get('meta', {})
        self.name = self.meta.get('name', os.path.basename(path))
        self.sources = d['sources']
        self.seeds = d['seeds']
        self.suites = d.get('suites', {})
        for sid, s in self.seeds.items():
            src = self.sources[s['source']]['sequence']
            L = s.get('length', 60)
            if s['offset'] < 0 or s['offset'] + L > len(src):
                raise SystemExit('panel %s: seed %s lies outside its source' % (self.name, sid))
            if s.get('sequence') and src[s['offset']:s['offset'] + L] != s['sequence']:
                raise SystemExit('panel %s: seed %s does not match its source at offset %d'
                                 % (self.name, sid, s['offset']))

    def window(self, sid, K, delta=0):
        """-> dict(source, start, end, sequence, wrapped, delta) of the seed window moved by K (+delta)."""
        s = self.seeds[sid]
        src = self.sources[s['source']]['sequence']
        L = s.get('length', 60)
        n = len(src) - L + 1
        start = s['offset'] + delta + K
        wrapped = not 0 <= start < n
        start %= n
        return {'seed': sid, 'source': s['source'], 'start': start, 'end': start + L,
                'sequence': src[start:start + L], 'wrapped': wrapped, 'delta': delta, 'K': K}

    def trace_delta(self, sid):
        """The trace walk's own window: far enough from the named walk that its first touch is cold."""
        s = self.seeds[sid]
        if 'trace_delta' in s:
            return s['trace_delta']
        src_len = len(self.sources[s['source']]['sequence'])
        L = s.get('length', 60)
        d = 2 * s.get('radius', 1000) + 200
        hi = src_len - L
        ok = [x for x in (d, -d) if 0 <= s['offset'] + x <= hi]
        if ok:
            # the direction that keeps the trace window farther from the source's ends
            return max(ok, key=lambda x: min(s['offset'] + x, hi - s['offset'] - x))
        # short source: the farthest position from the seed
        return hi - s['offset'] if hi - s['offset'] >= s['offset'] else -s['offset']


# ============================================================================ the request plan

def strategy(mode, radius, budget_ms, batch_kmers=DEFAULT_BATCH_KMERS, mlp=None, mlpn=None,
             detail='graphlet'):
    st = {'direction': 'both',
          'bounds': {'max_extension_bp': radius, 'time_budget_ms': budget_ms},
          'output': {'detail': detail, 'timing': True},
          'annotation': {'batch_kmers': batch_kmers}}
    if mode == 'derived_column':
        st['labels'] = {'seed_label_kind': 'column'}
    elif mode == 'named':
        pass
    elif mode == 'trace':
        st['support'] = 'trace'
        st['branching'] = {'on_reconverge': 'keep'}
    elif mode == 'annotate':
        st['labels'] = {'mode': 'annotate', 'seed_label_kind': 'column'}
        if mlpn is not None:
            st['labels']['max_labels_per_node'] = mlpn
    else:
        raise ValueError(mode)
    if mlp is not None:
        st['bounds']['max_live_paths'] = mlp
    return st


class Planner:
    def __init__(self, panel, K, coordinates=None):
        self.p = panel
        self.K = K
        self.reqs = []
        # --coordinates: the cap (an int or "unlimited") every support-trace walk asks for; None: none asked
        self.coordinates = coordinates

    def add(self, rid, suite, route, body=None, method='POST', windows=(), params=None):
        self.reqs.append({'id': rid, 'suite': suite, 'route': route, 'method': method, 'body': body,
                          'windows': [{k: v for k, v in w.items() if k != 'sequence'} for w in windows],
                          'params': params or {}})

    def win(self, sid, delta=0):
        return self.p.window(sid, self.K, delta)

    def seed_obj(self, w, labels=None):
        o = {'seed_id': w['seed'], 'sequence': w['sequence']}
        if labels is not None:
            o['labels'] = list(labels)
        return o

    def traverse(self, rid, suite, windows, mode, radius, budget_ms, labels=None, batch_kmers=DEFAULT_BATCH_KMERS,
                 mlp=None, mlpn=None, extra=None):
        seeds = []
        for i, w in enumerate(windows):
            seeds.append(self.seed_obj(w, labels[i] if labels is not None else None))
        body = {'seeds': seeds, 'strategy': strategy(mode, radius, budget_ms, batch_kmers, mlp, mlpn)}
        params = {'mode': mode, 'radius': radius, 'budget_ms': budget_ms, 'batch_kmers': batch_kmers,
                  'seeds': [w['seed'] for w in windows]}
        if mode == 'trace' and self.coordinates is not None:
            # record coordinates are reported under support trace only: the other walks stay opt-out
            body['strategy']['output'].update(coordinates=True,
                                              max_coordinate_occurrences=self.coordinates)
            params['coordinates'] = self.coordinates
        if labels is not None:
            params['labels'] = labels
        if mlp is not None:
            params['max_live_paths'] = mlp
        if mlpn is not None:
            params['max_labels_per_node'] = mlpn
        if extra:
            params.update(extra)
        self.add(rid, suite, '/traverse', body, windows=windows, params=params)

    def resolve(self, rid, suite, w, body_extra, params, sub=None):
        seq = w['sequence'] if sub is None else w['sequence'][sub[0]:sub[1]]
        body = dict({'sequence': seq}, **body_extra)
        ww = dict(w)
        if sub is not None:
            ww = dict(w, start=w['start'] + sub[0], end=w['start'] + sub[1])
        p = dict(params, seed=w['seed'], bp=len(seq), kmers=max(0, len(seq) - 30))
        self.add(rid, suite, '/resolve', body, windows=[ww], params=p)

    # ---------------------------------------------------------------- suites
    def caps(self):
        self.add('caps', 'caps', '/traverse/capabilities', None, method='GET', params={})

    def smoke(self):
        cfg = self.p.suites['smoke']
        sid = cfg['seed']
        s = self.p.seeds[sid]
        w = self.win(sid, cfg.get('delta', 0))   # a window of its own: the smoke must not warm a measured one
        self.resolve('smoke.%s.resolve_explicit' % sid, 'smoke', w,
                     {'labels': [s['record_header']], 'support': 'kmer'}, {'resolve': 'explicit'})
        self.traverse('smoke.%s.named_r100' % sid, 'smoke', [w], 'named', 100, WALK_BUDGET_MS,
                      labels=[[s['record_header']]])

    def budgets(self):
        roles = self.p.suites['budgets']
        plan = [('derived_1s', 'derived_column', 1000, 1000, DEFAULT_BATCH_KMERS),
                ('named_hdr_1s', 'named', 1000, 1000, DEFAULT_BATCH_KMERS),
                ('lookahead_250ms', 'derived_column', 3000, 250, 1000)]
        for role, mode, radius, budget, bk in plan:
            if role not in roles:
                continue
            sid = roles[role]
            w = self.win(sid)
            labels = [[self.p.seeds[sid]['record_header']]] if mode == 'named' else None
            for phase in ('first', 'warm'):
                self.traverse('budgets.%s.%s' % (role, phase), 'budgets', [w], mode, radius, budget,
                              labels=labels, batch_kmers=bk, extra={'phase': phase, 'role': role})

    def walks(self):
        for sid in self.p.suites['walks']['seeds']:
            s = self.p.seeds[sid]
            hdr = [[s['record_header']]]
            w = self.win(sid)
            for phase in ('first', 'warm'):
                self.traverse('walks.%s.named.%s' % (sid, phase), 'walks', [w], 'named', s['radius'],
                              WALK_BUDGET_MS, labels=hdr, extra={'phase': phase, 'class': s.get('class')})
            wt = self.win(sid, self.p.trace_delta(sid))
            for phase in ('first', 'warm'):
                self.traverse('walks.%s.trace.%s' % (sid, phase), 'walks', [wt], 'trace', s['radius'],
                              WALK_BUDGET_MS, labels=hdr, extra={'phase': phase, 'class': s.get('class')})

    def callsize(self):
        for sid in self.p.suites['callsize']['seeds']:
            s = self.p.seeds[sid]
            w = self.win(sid)
            ex = {'labels': s['headers'][:3], 'support': 'trace'}
            for tag, sub in (('k30_first', None), ('k1', (0, 31)), ('k3', (10, 43)), ('k30_warm', None)):
                self.resolve('callsize.%s.%s' % (sid, tag), 'callsize', w, ex,
                             {'resolve': 'explicit_trace', 'tag': tag}, sub)

    def resolve_suite(self):
        for sid in self.p.suites['resolve']['seeds']:
            s = self.p.seeds[sid]
            w = self.win(sid)
            ex = {'labels': s['columns'][:3], 'support': 'kmer'}
            self.resolve('resolve.%s.explicit_first' % sid, 'resolve', w, ex, {'resolve': 'explicit'})
            self.resolve('resolve.%s.explicit_warm' % sid, 'resolve', w, ex, {'resolve': 'explicit'})
            for kind in ('column', 'header'):
                self.resolve('resolve.%s.discover_%s' % (sid, kind), 'resolve', w,
                             {'discover': {'kind': kind, 'max_labels': 5}}, {'resolve': 'discover_' + kind})

    def annotate(self):
        cfg = self.p.suites['annotate']
        sid = cfg['seed']
        w = self.win(sid)
        for phase in ('first', 'warm'):
            self.traverse('annotate.%s.%s' % (sid, phase), 'annotate', [w], 'annotate', cfg.get('radius', 1000),
                          WALK_BUDGET_MS, mlp=64, mlpn=64, extra={'phase': phase})

    def lookahead(self):
        cfg = self.p.suites['lookahead']
        sid = cfg['seed']
        w = self.win(sid)
        radius = cfg.get('radius', 300)
        for tag, bk in (('b64_first', 64), ('b64_warm', 64), ('b1', 1), ('b1000', 1000)):
            self.traverse('lookahead.%s.%s' % (sid, tag), 'lookahead', [w], 'derived_column', radius,
                          WALK_BUDGET_MS, batch_kmers=bk, extra={'tag': tag})

    def batch(self):
        cfg = self.p.suites['batch']
        sids = cfg['seeds']
        ws = [self.win(sid) for sid in sids]
        labels = [[self.p.seeds[sid]['record_header']] for sid in sids]
        self.traverse('batch.%dseeds' % len(sids), 'batch', ws, 'named', cfg.get('radius', 1000), WALK_BUDGET_MS,
                      labels=labels)

    def build(self, suites):
        fn = {'caps': self.caps, 'smoke': self.smoke, 'budgets': self.budgets, 'walks': self.walks,
              'callsize': self.callsize, 'resolve': self.resolve_suite, 'annotate': self.annotate,
              'lookahead': self.lookahead, 'batch': self.batch}
        for s in suites:
            if s == 'smoke' or s in self.p.suites or s == 'caps':
                fn[s]()
            else:
                print('note: the panel has no "%s" suite; skipped' % s, file=sys.stderr)
        return self.reqs


def parse_suites(text):
    out = []
    for s in text.split(','):
        s = s.strip()
        if not s:
            continue
        if s == 'all':
            out += ALL_SUITES
        elif s in ALL_SUITES or s in EXTRA_SUITES:
            out.append(s)
        else:
            raise SystemExit('unknown suite %r (all, %s)' % (s, ', '.join(ALL_SUITES + EXTRA_SUITES)))
    if 'caps' not in out:
        out.insert(0, 'caps')
    seen, uniq = set(), []
    for s in out:
        if s not in seen:
            seen.add(s)
            uniq.append(s)
    # smoke already starts with caps
    return uniq


def first_touch_windows(reqs):
    """Windows whose first touch the run measures: budgets/walks/smoke requests (first phase)."""
    out = []
    for r in reqs:
        if r['suite'] in ('budgets', 'walks', 'smoke') and r['params'].get('phase', 'first') == 'first':
            out += r['windows']
    return out


def avoided_intervals(panel, avoid_files):
    """Read intervals of earlier runs on the SAME panel (another panel is another index: ignored)."""
    avoid, used = {}, []
    for f in avoid_files:
        d = load_results(f)
        name = (d.get('meta', {}).get('panel') or {}).get('name')
        if name not in (None, panel.name):
            print('warning: %s was run on panel %s, not %s: ignored' % (f, name, panel.name), file=sys.stderr)
            continue
        for src, ivs in (d.get('read_intervals') or {}).items():
            avoid.setdefault(src, []).extend(ivs)
        used.append(f)
    return avoid, used


def choose_offset(panel, suites, avoid_files, step=100, max_k=40000, margin=100):
    avoid, used = avoided_intervals(panel, avoid_files)
    if not avoid:
        return 0, 0, 'no read intervals to avoid'
    best = None
    for K in range(0, max_k + 1, step):
        reqs = Planner(panel, K).build(suites)
        n = 0
        for w in first_touch_windows(reqs):
            lo, hi = w['start'] - margin, w['end'] + margin
            for a, b in avoid.get(w['source'], []):
                if lo < b and a < hi:
                    n += 1
                    break
        if best is None or n < best[1]:
            best = (K, n)
        if n == 0:
            break
    return best[0], best[1], 'chose K=%d: %d first-touch windows overlap the %d avoided runs' % (
        best[0], best[1], len(used))


# ============================================================================ the runner

class Runner:
    def __init__(self, args, out_dir):
        self.args = args
        self.base = args.base.rstrip('/')
        self.local = is_local(self.base)
        self.out = out_dir
        self.raw = os.path.join(out_dir, 'raw')
        os.makedirs(self.raw, exist_ok=True)
        self.pause = args.pause_s if self.local else max(2.0, args.pause_s)
        host = urllib.parse.urlparse(self.base)
        key = re.sub(r'[^A-Za-z0-9._-]+', '_', '%s_%s' % (host.hostname, host.port or host.scheme))
        tmp = tempfile.gettempdir()
        self.state_path = os.path.join(tmp, 'metagraph-bench-traverse-%s.state.json' % key)
        self.lock_path = os.path.join(tmp, 'metagraph-bench-traverse-%s.lock' % key)
        self.lock_f = None
        self.sent = 0
        self.consec_fail = 0
        self.slow = 0
        self.caps = None

    # ---------------------------------------------------------------- lock + pace
    def lock(self):
        if fcntl is None:
            return
        self.lock_f = open(self.lock_path, 'w')
        try:
            fcntl.flock(self.lock_f, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            raise SystemExit('another bench_traverse run against this server is active (lock %s)'
                             % self.lock_path) from None

    def unlock(self):
        if self.lock_f is not None:
            with contextlib.suppress(OSError):
                fcntl.flock(self.lock_f, fcntl.LOCK_UN)
            self.lock_f.close()
            self.lock_f = None

    def _state(self):
        try:
            with open(self.state_path) as f:
                return json.load(f)
        except (OSError, ValueError):
            return {'last_end': 0.0, 'last_wall': 0.0}

    def _save_state(self, last_end, last_wall):
        with contextlib.suppress(OSError):
            write_json_atomic(self.state_path, {'last_end': last_end, 'last_wall': last_wall, 'base': self.base})

    def pace(self):
        st = self._state()
        need = self.pause
        if not self.local:
            if st['last_wall'] > 3.0:
                need = max(need, 2 * self.pause, 6.0)
            if st['last_wall'] > 20.0:
                need = max(need, 10.0)
        wait = st['last_end'] + need - time.time()
        if wait > 0:
            time.sleep(wait)

    # ---------------------------------------------------------------- one request
    def send(self, req, n):
        body = None if req['body'] is None else json.dumps(req['body']).encode('utf-8')
        hdrs = {'Accept-Encoding': 'gzip', 'User-Agent': 'metagraph-bench-traverse/%d' % SCRIPT_VERSION}
        if body is not None:
            hdrs['Content-Type'] = 'application/json'
        http_req = urllib.request.Request(self.base + req['route'], data=body, method=req['method'], headers=hdrs)
        self.pace()
        started = now_iso()
        t0 = time.time()
        code, wire, enc, err, t_hdr = None, b'', None, None, None
        try:
            with urllib.request.urlopen(http_req, timeout=self.args.timeout_s) as r:
                t_hdr = time.time() - t0
                code, enc = r.status, r.headers.get('Content-Encoding')
                wire = r.read()
        except urllib.error.HTTPError as e:
            t_hdr = time.time() - t0
            code, enc = e.code, e.headers.get('Content-Encoding')
            with contextlib.suppress(Exception):
                wire = e.read()
        except Exception as e:  # timeouts, connection errors
            err = '%s: %s' % (type(e).__name__, e)
        wall = time.time() - t0
        self._save_state(time.time(), wall)
        try:
            text = gzip.decompress(wire) if enc == 'gzip' and wire else wire
        except OSError as e:
            text, err = wire, (err or '') + ' gzip: %s' % e
        try:
            resp = json.loads(text) if text else None
        except ValueError:
            resp = {'_unparsed': text[:4000].decode('utf-8', 'replace')}
        meta = {'n': n, 'id': req['id'], 'suite': req['suite'], 'route': req['route'], 'method': req['method'],
                'params': req['params'], 'windows': req['windows'], 'started_at': started,
                'wall_ms': wall * 1000.0, 'headers_ms': None if t_hdr is None else t_hdr * 1000.0,
                'http': code, 'error': err, 'content_encoding': enc, 'bytes_wire': len(wire),
                'bytes_json': len(text), 'request_bytes': len(body) if body else 0,
                'request_sha': sha(req['body']) if req['body'] is not None else None,
                'request_sha_nocoords': (sha(without_coordinates(req['body']))
                                         if req['body'] is not None else None),
                'base': self.base}
        fname = '%03d_%s.json.gz' % (n, re.sub(r'[^A-Za-z0-9._-]+', '_', req['id']))
        with gzip.open(os.path.join(self.raw, fname), 'wt') as f:
            json.dump({'meta': meta, 'request': req['body'], 'response': resp}, f)
        meta['raw_file'] = 'raw/' + fname
        return meta, resp

    # ---------------------------------------------------------------- the run
    def run(self, reqs, on_result):
        stop = None
        for req in reqs:
            if self.sent >= self.args.max_requests:
                stop = 'hard cap of %d requests reached (%d planned requests not sent)' % (
                    self.args.max_requests, len(reqs) - self.sent)
                break
            if req['suite'] == 'batch' and self.caps and isinstance(self.caps.get('max_seeds'), int):
                ms = self.caps['max_seeds']
                if ms and len(req['body']['seeds']) > ms:
                    req['body']['seeds'] = req['body']['seeds'][:ms]
                    req['windows'] = req['windows'][:ms]
                    req['params']['seeds'] = req['params']['seeds'][:ms]
                    req['params']['note'] = 'batch cut to the server max_seeds %d' % ms
            self.sent += 1
            meta, resp = self.send(req, self.sent)
            if req['route'] == '/traverse/capabilities' and isinstance(resp, dict) and meta['http'] == 200:
                self.caps = resp
            rec = on_result(meta, req['body'], resp)
            print(brief_line(rec), flush=True)
            failed = meta['http'] is None or meta['http'] >= 500
            self.consec_fail = self.consec_fail + 1 if failed else 0
            if self.consec_fail >= 3:
                stop = 'three consecutive server errors / connection failures: stopped'
                break
            if meta['wall_ms'] > 60000 and not self.local:
                self.slow += 1
                if self.slow >= 2:
                    stop = 'a second request over 60 s: stopped'
                    break
                print('     a request took over 60 s: pausing 60 s', flush=True)
                time.sleep(60)
        return stop


# ============================================================================ record coordinates

def without_coordinates(body):
    """|body| without strategy.output.coordinates / max_coordinate_occurrences (a copy): the identity of a
    request apart from the coordinates it asks for, so that a --coordinates run pairs with an opt-out one."""
    if not isinstance(body, dict) or not isinstance((body.get('strategy') or {}).get('output'), dict):
        return body
    out = json.loads(json.dumps(body))
    o = out['strategy']['output']
    o.pop('coordinates', None)
    o.pop('max_coordinate_occurrences', None)
    return out


def coordinates_facts(res):
    """What a seed result's coordinates field says (None when the request asked for none): the block's kind, its
    compact JSON bytes and share of the result, entries, occurrences, lists cut, lower-bound runs, the largest
    list, complete; or the null reason."""
    if 'coordinates' not in res:
        return None
    c = res['coordinates']
    if not isinstance(c, dict):
        return {'present': False, 'reason': res.get('coordinates_reason')}
    block = len(json.dumps(c, separators=(',', ':'), sort_keys=True).encode('utf-8'))
    whole = len(json.dumps(res, separators=(',', ':'), sort_keys=True).encode('utf-8'))
    seed = c.get('seed') or []
    runs = [e for es in (c.get('arms') or {}).values() for e in es]
    lists = seed + runs
    totals = [e.get('occurrences_total', len(e.get('occurrences') or [])) for e in lists]
    return {'present': True, 'kind': c.get('kind'), 'max_occurrences': c.get('max_occurrences'),
            'complete': c.get('complete'), 'block_bytes': block, 'share': div(block, whole),
            'seed_entries': len(seed), 'runs': len(runs),
            'occurrences': sum(len(e.get('occurrences') or []) for e in lists),
            'occurrences_total': sum(totals), 'max_list': max(totals, default=0),
            'lists_cut': sum(1 for e in lists if 'occurrences_total' in e),
            'runs_lower_bound': c.get('runs_lower_bound', 0),
            'runs_chains_ended': sum(1 for e in runs if e.get('chains_ended'))}


def estimate_block(res, digits=7):
    """The coordinates block a support-trace seed result without one would carry, built as the server writes it
    (keys sorted, compact JSON) with ONE occurrence per entry whose positions have |digits| digits: (block, its
    compact bytes, the result's compact bytes). Counts are exact (an entry per seed label, one per R record of
    each arm); the occurrences per entry are not in an opt-out response (a lower bound where a record repeats the
    walked sequence)."""
    body = res.get('graphlet') or ''
    k = None
    nsl = (res.get('seed') or {}).get('num_seed_labels') or 0
    seed_bp = (res.get('seed') or {}).get('length_bp') or 0
    arms = {}
    kinds = set()
    for line in body.split('\n'):
        f = line.split(' ')
        if f[0] == 'H' and len(f) > 3:
            k = int(f[3])
        elif f[0] == 'L' and len(f) > 1:
            kinds.add(f[1])
        elif f[0] == 'A' and len(f) > 1:
            cur = arms.setdefault({'l': 'left', 'r': 'right'}[f[1]], [])
        elif f[0] == 'R' and len(f) > 5 and arms:
            # R <segment> <label> <from_bp> <to_bp> <end> ... (the run's own interval is enough here)
            try:
                a, b = int(f[3]), int(f[4])
            except ValueError:
                a, b = 0, 0
            cur.append((len(cur), int(f[2]) if f[2].isdigit() else 0, a, b))
    pos = 10 ** (digits - 1)
    block = {'complete': True, 'k': k or 31, 'max_occurrences': COORDINATES_CAP_DEFAULT,
             'kind': 'mixed' if len(kinds) > 1 else 'record' if kinds == {'h'} else 'column',
             'seed': [{'label': i, 'occurrences': [[pos, pos + seed_bp]]} for i in range(nsl)],
             'arms': {side: [{'from_bp': a, 'label': lab, 'occurrences': [[pos, pos + b - a]], 'run': i,
                              'to_bp': b} for i, lab, a, b in runs] for side, runs in arms.items()}}
    nb = len(json.dumps(block, separators=(',', ':'), sort_keys=True).encode('utf-8'))
    whole = len(json.dumps(res, separators=(',', ':'), sort_keys=True).encode('utf-8'))
    return block, nb, whole


def coordinates_estimate(run_dirs, digits=7, out=print):
    """--coordinates-estimate: the table over every support-trace seed result in the runs' raw responses."""
    rows = []
    for d in run_dirs:
        found = []
        for p in sorted(glob.glob(os.path.join(d, 'raw', '[0-9]*.json.gz'))):
            with gzip.open(p, 'rt') as f:
                e = json.load(f)
            req, resp = e.get('request') or {}, e.get('response') or {}
            st = (req.get('strategy') or {}) if isinstance(req, dict) else {}
            if st.get('support') != 'trace' or not isinstance(resp, dict):
                continue
            for i, r in enumerate(resp.get('results') or []):
                if not isinstance(r.get('graphlet'), str):
                    continue
                block, nb, whole = estimate_block(r, digits)
                n_runs = sum(len(v) for v in block['arms'].values())
                row = {'run': os.path.basename(os.path.normpath(d)), 'id': e['meta']['id'], 'result': i,
                       'kind': block['kind'], 'seed_labels': len(block['seed']), 'runs': n_runs,
                       'block_bytes_at_least': nb, 'result_bytes': whole,
                       'share_at_least': div(nb, whole + nb),
                       'graphlet_bytes': len(r['graphlet'].encode('utf-8')),
                       'elapsed_ms': _num((r.get('timing') or {}).get('elapsed_ms'))}
                found.append(row)
                rows.append(row)
        if found and os.path.isdir(d):
            write_json_atomic(os.path.join(d, 'coordinates_estimate.json'),
                              {'format': 'metagraph-bench-coordinates-estimate', 'version': SCRIPT_VERSION,
                               'digits': digits, 'rows': found})
    out('coordinates block estimate (one occurrence per entry, %d-digit positions; a lower bound):' % digits)
    out(md_table(['run', 'id', 'kind', 'seed labels', 'runs', 'block B >=', 'result B', 'share >=',
                  'graphlet B', 'elapsed ms'],
                 [[r['run'], '%s#%d' % (r['id'], r['result']), r['kind'], r['seed_labels'], r['runs'],
                   r['block_bytes_at_least'], r['result_bytes'], fmt(r['share_at_least'], 3),
                   r['graphlet_bytes'], fmt(r['elapsed_ms'])] for r in rows]))
    if rows:
        xs = [r['block_bytes_at_least'] for r in rows]
        out('\n%d trace seed results: block >= %s B median (%s .. %s), share >= %s median'
            % (len(rows), fmt(statistics.median(xs), 0), min(xs), max(xs),
               fmt(statistics.median([r['share_at_least'] for r in rows]), 3)))
    return rows


# ============================================================================ analysis of one request

def _num(x):
    return x if isinstance(x, (int, float)) and not isinstance(x, bool) else None


def lim_strings(res):
    out = []
    for l in res.get('limitations') or []:
        out.append('%s:%s' % (l.get('kind'), l.get('knob')))
    for side, a in (res.get('arms') or {}).items():
        for l in a.get('limitations') or []:
            out.append('%s.%s:%s' % (side[0], l.get('kind'), l.get('knob')))
    return out


def time_bound(res):
    for l in res.get('limitations') or []:
        if l.get('knob') == 'bounds.time_budget_ms':
            return True
    for a in (res.get('arms') or {}).values():
        if (a.get('cap_trigger') or {}).get('reason') == 'time_budget':
            return True
        for l in a.get('limitations') or []:
            if l.get('knob') == 'bounds.time_budget_ms':
                return True
    rs = res.get('resource_stop') or {}
    if rs.get('resource') == 'time':
        return True
    err = res.get('error') or ''
    return 'time budget' in err or 'time_budget' in err


def graphlet_hashes(text):
    """-> (whole body, body without the H line, walk content). The walk content leaves out the H/A/B/Q/K/O/Z
    bookkeeping lines and the validated seed id in S (it changes with the server's release / index identity)."""
    lines = text.split('\n')
    noh = '\n'.join(l for l in lines if not l.startswith('H '))
    walk = []
    for l in lines:
        if l[:1] in WALK_LINES and l[1:2] == ' ':
            if l.startswith('S '):
                f = l.split(' ')
                l = ' '.join(['S', '*'] + f[2:])
            walk.append(l)
    return sha(text), sha(noh), sha('\n'.join(walk))


def analyze_seed_result(res, budget_ms, window, n_seeds, top_elapsed, usage_work=None):
    t = res.get('timing') or {}
    timing = {k: _num(t.get(k)) for k in TIMING_FIELDS}
    for k, v in t.items():          # newer servers: any extra numeric timing field
        if k not in timing and _num(v) is not None:
            timing[k] = v
        elif isinstance(v, dict):   # and the numbers of a nested object (timing.path_cache), as "name.field"
            for kk, vv in v.items():
                if _num(vv) is not None:
                    timing['%s.%s' % (k, kk)] = vv
    arms = {}
    for side, a in (res.get('arms') or {}).items():
        c = a.get('counts') or {}
        lbr = c.get('leaves_by_reason')
        arms[side] = {'status': a.get('status'), 'complete_to_bp': _num(a.get('complete_to_bp')),
                      'max_bp': _num(c.get('max_bp')), 'bases': _num(c.get('bases')),
                      'leaves': _num(c.get('leaves')), 'segments': _num(c.get('segments')),
                      'cap_reason': (a.get('cap_trigger') or {}).get('reason'),
                      'leaves_by_reason': ({k: v for k, v in lbr.items() if _num(v) is not None}
                                           if isinstance(lbr, dict) else None),
                      'counters': {k: v for k, v in (a.get('counters') or {}).items() if _num(v) is not None}}
    seed = res.get('seed') or {}
    rs = res.get('resource_stop')
    tb = time_bound(res)
    out = {'seed_id': seed.get('seed_id'), 'window': window, 'outcome': res.get('outcome'),
           'error': res.get('error'), 'timing': timing, 'annotation': res.get('annotation'),
           'arms': arms, 'limitations': lim_strings(res),
           'resource_stop': ({k: rs.get(k) for k in ('resource', 'phase', 'requested', 'effective', 'used')}
                             if isinstance(rs, dict) else None),
           'labels': {k: seed.get(k) for k in ('num_labels', 'num_seed_labels', 'labels_supporting_total',
                                               'labels_dropped', 'labels_from_seed')},
           'label_mode': res.get('label_mode'), 'budget_ms': budget_ms, 'time_bound': tb}
    if 'graphlet' in res and isinstance(res['graphlet'], str):
        g, gn, gw = graphlet_hashes(res['graphlet'])
        out.update(graphlet_bytes=len(res['graphlet'].encode('utf-8')), graphlet_sha=g, graphlet_sha_noH=gn,
                   walks_sha=gw)
    else:
        body = {k: v for k, v in res.items() if k != 'timing'}
        out['result_sha'] = sha(body)
    out['deterministic'] = not tb and not res.get('error')
    coords = coordinates_facts(res)
    if coords is not None:
        out['coordinates'] = coords
    # derived
    el = timing.get('elapsed_ms')
    if el is None and top_elapsed is not None:
        el = top_elapsed / max(1, n_seeds)
    rows = None
    if timing.get('rows_fetched') is not None or timing.get('tuple_rows_fetched') is not None:
        rows = (timing.get('rows_fetched') or 0) + (timing.get('tuple_rows_fetched') or 0)
    bpl = (arms.get('left') or {}).get('complete_to_bp')
    bpr = (arms.get('right') or {}).get('complete_to_bp')
    bp_total = None if bpl is None and bpr is None else (bpl or 0) + (bpr or 0)
    bases = sum((a.get('bases') or 0) for a in arms.values()) if arms else None
    out['derived'] = {
        'elapsed_ms': el, 'rows': rows,
        'ms_per_row': div(timing.get('annotation_fetch_ms'), rows),
        'elapsed_ms_per_row': div(el, rows),
        'fetch_share': div(timing.get('annotation_fetch_ms'), el),
        'budget_ratio': div(el, budget_ms),
        'bp_left': bpl, 'bp_right': bpr, 'bp_total': bp_total,
        'reach_left': (arms.get('left') or {}).get('max_bp'), 'reach_right': (arms.get('right') or {}).get('max_bp'),
        'bp_per_s': div(bp_total, (el or 0) / 1000.0) if el else None,
        'rows_per_base': div(rows, bases) if bases else None,
    }
    reach_and_work(out, usage_work)
    return out


def _arm_total(arms, get):
    vals = [get(a) for a in arms.values()]
    vals = [v for v in vals if v is not None]
    return sum(vals) if vals else None


def stop_reason(s):
    """What ended a seed's exploration, in one word: failed, resource:<r> (a memory / work budget), the arms'
    cap reasons (time_budget, max_steps, ...), else the walks outcome (complete)."""
    if s.get('error'):
        return 'failed'
    rs = s.get('resource_stop') or {}
    if rs.get('resource'):
        return 'resource:%s' % rs['resource']
    caps = sorted({a.get('cap_reason') for a in (s.get('arms') or {}).values() if a.get('cap_reason')})
    if caps:
        return ','.join(caps)
    return (s.get('outcome') or {}).get('walks') or 'complete'


def reach_and_work(s, usage_work=None):
    """Adds to an analysed seed result (also one read back from an older results.json, whose arms may lack the
    counters and leaves_by_reason: those numbers are then None):
      reach     the CERTIFIED reach -- bp (complete_to_bp, both arms), walked_bp (the longest walks, counts.max_bp),
                steps (walker steps: the nodes walked), segments, leaves, complete_walks (leaves that ended for a
                reason of the walk's own, not a resource stop) and stopped_walks, stop (stop_reason)
      work      the CONSUMED work -- work_units (the arms' counter, reported under a memory / work budget),
                usage_work_units (the attempt's usage, ledger-managed requests), and the logical counters every
                server reports: steps, successor_enumerations, annotation.rows_requested / keys_mapped; rows_decoded
                (rows_fetched + tuple_rows_fetched) is physical (lookahead and caches change it)
      physical  the path cache's reuse (timing.path_cache.hits: paths cut at a cached row), stored_rows_read,
                peak_cache_bytes; None where the server reports none (absent on a server that reports
                timing.seed_phase_ms, the same release: 0, since the object is omitted when all counters are 0)"""
    arms = s.get('arms') or {}
    d = s.get('derived') or {}
    known = bool(arms) and all(isinstance(a.get('leaves_by_reason'), dict) for a in arms.values())
    complete = stopped = None
    if known:
        complete = sum(n for a in arms.values() for r, n in a['leaves_by_reason'].items()
                       if r not in RESOURCE_LEAF_REASONS)
        stopped = sum(n for a in arms.values() for r, n in a['leaves_by_reason'].items()
                      if r in RESOURCE_LEAF_REASONS)

    def counter(name):
        return _arm_total(arms, lambda a: _num((a.get('counters') or {}).get(name)))
    s['reach'] = {'bp': d.get('bp_total'), 'bp_left': d.get('bp_left'), 'bp_right': d.get('bp_right'),
                  'walked_bp': _arm_total(arms, lambda a: a.get('max_bp')), 'steps': counter('steps'),
                  'segments': _arm_total(arms, lambda a: a.get('segments')),
                  'leaves': _arm_total(arms, lambda a: a.get('leaves')), 'complete_walks': complete,
                  'stopped_walks': stopped, 'stop': stop_reason(s)}
    ann = s.get('annotation') or {}
    s['work'] = {'work_units': counter('work_units'), 'usage_work_units': usage_work,
                 'steps': counter('steps'), 'successor_enumerations': counter('successor_enumerations'),
                 'rows_requested': _num(ann.get('rows_requested')), 'keys_mapped': _num(ann.get('keys_mapped')),
                 'rows_decoded': d.get('rows')}
    timing = s.get('timing') or {}
    phys = {'reuse': None, 'stored_rows_read': None, 'peak_cache_bytes': None, 'fields': []}
    for k, v in timing.items():
        if v is None:
            continue
        if PATH_CACHE_REUSE.search(k) and phys['reuse'] is None:
            phys['reuse'] = v
            phys['fields'].append(k)
        elif PATH_CACHE_PEAK.search(k) and phys['peak_cache_bytes'] is None:
            phys['peak_cache_bytes'] = v
            phys['fields'].append(k)
        elif k == 'path_cache.stored_rows_read':
            phys['stored_rows_read'] = v
            phys['fields'].append(k)
    if phys['reuse'] is None and 'seed_phase_ms' in timing:
        phys.update(reuse=0, stored_rows_read=0, peak_cache_bytes=0)
    s['physical'] = phys
    return s


def analyze(meta, request, resp):
    rec = {k: meta.get(k) for k in ('n', 'id', 'suite', 'route', 'method', 'params', 'windows', 'started_at',
                                    'wall_ms', 'headers_ms', 'http', 'error', 'content_encoding', 'bytes_wire',
                                    'bytes_json', 'request_bytes', 'request_sha', 'request_sha_nocoords',
                                    'raw_file')}
    if not isinstance(resp, dict):
        return rec
    caps = resp if meta['route'] == '/traverse/capabilities' else (resp.get('capabilities') or {})
    rec['release'] = resp.get('release', caps.get('release') if isinstance(caps, dict) else None)
    rec['algorithm_version'] = resp.get('algorithm_version')
    if isinstance(caps, dict):
        rec['feature_level'] = caps.get('feature_level')
        rec['index_fp'] = caps.get('index_fp')
    if meta['http'] != 200:
        rec['error_body'] = (resp.get('error') or resp.get('_unparsed') or json.dumps(resp))[:800]
        return rec
    if meta['route'] == '/traverse/capabilities':
        rec['capabilities'] = {k: resp.get(k) for k in ('release', 'feature_level', 'max_seeds', 'max_time_ms',
                                                         'max_seed_labels', 'index_fp', 'index_ns', 'regime', 'k',
                                                         'num_labels', 'schema_version')}
        rec['capabilities']['attempts'] = 'attempts' in resp
        return rec
    timing = resp.get('timing') or {}
    rec['server'] = {'elapsed_ms': _num(timing.get('elapsed_ms')), 'serialize_ms': _num(timing.get('serialize_ms'))}
    rec['overhead_ms'] = (meta['wall_ms'] - rec['server']['elapsed_ms']) if rec['server']['elapsed_ms'] else None
    if meta['route'] == '/traverse':
        st = resp.get('strategy') or {}
        budget = _num((st.get('bounds') or {}).get('time_budget_ms'))
        rec['clamped'] = st.get('clamped') or None
        res = resp.get('results') or []
        wins = meta.get('windows') or []
        u = resp.get('usage')
        per_seed = [None] * len(res)
        if isinstance(u, dict):
            # an attempt's usage (ledger-managed requests): the work its seeds consumed
            rec['usage'] = {'work_units': _num(u.get('work_units')), 'reason': u.get('reason'),
                            'elapsed_ms': _num(u.get('elapsed_ms')),
                            'peak_admitted_bytes': _num((u.get('memory') or {}).get('peak_admitted_bytes'))}
            for e in u.get('per_seed') or []:
                i = e.get('index')
                if isinstance(i, int) and 0 <= i < len(res):
                    per_seed[i] = _num(e.get('work_units'))
            if len(res) == 1 and per_seed[0] is None:
                per_seed[0] = rec['usage']['work_units']
        rec['results'] = [analyze_seed_result(r, budget, wins[i] if i < len(wins) else None, len(res),
                                              rec['server']['elapsed_ms'], per_seed[i]) for i, r in enumerate(res)]
        if len(res) > 1 and budget:
            rec['batch_budget_ratio'] = div(rec['server']['elapsed_ms'], budget * len(res))
    elif meta['route'] == '/resolve':
        labels = resp.get('labels') or []
        body = {k: v for k, v in resp.items() if k not in ('timing', 'capabilities', 'release')}
        nk = _num(resp.get('num_kmers'))
        rec['resolve'] = {'elapsed_ms': rec['server']['elapsed_ms'], 'num_kmers': nk, 'n_labels': len(labels),
                          'labels_truncated': resp.get('labels_truncated'), 'graph_runs': resp.get('graph_runs'),
                          'support': resp.get('support'), 'ms_per_kmer': div(rec['server']['elapsed_ms'], nk),
                          'result_sha': sha(body),
                          'labels_full_length': sum(1 for l in labels if nk and l.get('kmers_supported') == nk)}
    return rec


def brief_line(rec):
    head = '#%03d %-42s http=%s wall=%7.1f ms' % (rec.get('n') or 0, rec['id'][:42], rec.get('http'),
                                                 rec.get('wall_ms') or 0)
    if rec.get('error'):
        return head + ' ERROR ' + rec['error'][:200]
    if rec.get('http') != 200:
        return head + ' ' + (rec.get('error_body') or '')[:200]
    if 'capabilities' in rec:
        c = rec['capabilities']
        return head + ' release=%s feature_level=%s max_seeds=%s max_time_ms=%s attempts=%s' % (
            c.get('release'), c.get('feature_level'), c.get('max_seeds'), c.get('max_time_ms'), c.get('attempts'))
    if 'resolve' in rec:
        z = rec['resolve']
        return head + ' server=%s ms kmers=%s labels=%s ms/kmer=%s' % (
            fmt(z['elapsed_ms']), z['num_kmers'], z['n_labels'], fmt(z['ms_per_kmer'], 2))
    parts = []
    for s in rec.get('results') or []:
        d = s['derived']
        arms = ' '.join('%s:%s/%s' % (k[0], (v['status'] or '?')[:4], v['complete_to_bp'])
                        for k, v in s['arms'].items())
        parts.append('el=%s fetch=%s rows=%s ms/row=%s ratio=%s %s %s%s' % (
            fmt(d['elapsed_ms']), fmt(s['timing'].get('annotation_fetch_ms')), d['rows'], fmt(d['ms_per_row'], 2),
            fmt(d['budget_ratio'], 2), (s.get('outcome') or {}).get('walks'), arms,
            ' ERR' if s.get('error') else ''))
    return head + ' ' + ' || '.join(parts)


# ============================================================================ run-level tables

def by_id(records):
    return {r['id']: r for r in records}


def seed_res(rec, i=0):
    if not rec or rec.get('http') != 200:
        return None
    res = rec.get('results') or []
    return res[i] if i < len(res) else None


def same_result(a, b, key='walks'):
    """-> (verdict, detail) for two seed results of the same request."""
    if a is None or b is None:
        return 'n/a', 'missing'
    if not (a.get('deterministic') and b.get('deterministic')):
        return 'time-bound', 'not expected to repeat'
    if 'graphlet_sha' in a and 'graphlet_sha' in b:
        if key == 'noH':
            return ('equal' if a['graphlet_sha_noH'] == b['graphlet_sha_noH'] else 'DIFFERENT'), ''
        if a['graphlet_sha'] == b['graphlet_sha']:
            return 'equal', ''
        if a['walks_sha'] == b['walks_sha']:
            return 'walks equal', 'bookkeeping lines differ'
        return 'DIFFERENT', 'walk content differs'
    return ('equal' if a.get('result_sha') == b.get('result_sha') else 'DIFFERENT'), ''


def invariants(records):
    R = by_id(records)
    out = []

    def check(name, a_id, b_id, ia=0, ib=0, key='walks'):
        a, b = seed_res(R.get(a_id), ia), seed_res(R.get(b_id), ib)
        if a is None or b is None:
            return
        v, d = same_result(a, b, key)
        out.append({'check': name, 'a': a_id, 'b': b_id, 'verdict': v, 'detail': d})

    for rid in R:
        m = re.match(r'^(walks\.[^.]+\.(?:named|trace)|annotate\.[^.]+|budgets\.[^.]+)\.first$', rid)
        if m:
            check('first == warm repeat', rid, m.group(1) + '.warm')
    for rid in R:
        m = re.match(r'^lookahead\.([^.]+)\.b64_first$', rid)
        if m:
            base = 'lookahead.%s.' % m.group(1)
            for other in ('b64_warm', 'b1', 'b1000'):
                check('identical across batch_kmers', rid, base + other)
    for rid, rec in R.items():
        if rid.startswith('batch.') and rec.get('http') == 200:
            for i in range(len(rec.get('results') or [])):
                sid = (rec['params'].get('seeds') or [None] * (i + 1))[i]
                single = 'walks.%s.named.warm' % sid
                if single not in R:
                    single = 'walks.%s.named.first' % sid
                srec = R.get(single)
                if srec and srec['params'].get('radius') == rec['params'].get('radius') \
                        and srec.get('request_sha') and _same_seed_window(srec, rec, sid):
                    check('batch seed == single', rid, single, ia=i, ib=0, key='noH')
    for rid, rec in R.items():
        if rec.get('route') != '/resolve':
            continue
        m = re.match(r'^(callsize\.[^.]+)\.k30_first$|^(resolve\.[^.]+)\.explicit_first$', rid)
        if m:
            other = (m.group(1) + '.k30_warm') if m.group(1) else (m.group(2) + '.explicit_warm')
            a, b = R.get(rid), R.get(other)
            if a and b and a.get('http') == 200 and b.get('http') == 200:
                v = 'equal' if a['resolve']['result_sha'] == b['resolve']['result_sha'] else 'DIFFERENT'
                out.append({'check': 'resolve first == warm', 'a': rid, 'b': other, 'verdict': v, 'detail': ''})
    return out


def _same_seed_window(single, batch, sid):
    ws = [w for w in single.get('windows') or [] if w.get('seed') == sid]
    wb = [w for w in batch.get('windows') or [] if w.get('seed') == sid]
    return bool(ws and wb and ws[0]['start'] == wb[0]['start'] and ws[0]['source'] == wb[0]['source'])


def read_intervals(records):
    """Source intervals each request read: window +- reach + margin (k + batch_kmers lookahead)."""
    iv = {}
    for rec in records:
        if rec.get('http') != 200 or rec.get('route') not in ('/traverse', '/resolve'):
            continue
        res = rec.get('results') or [None] * len(rec.get('windows') or [])
        bk = (rec.get('params') or {}).get('batch_kmers') or 0
        for i, w in enumerate(rec.get('windows') or []):
            s = res[i] if i < len(res) else None
            left = right = 0
            if s:
                left = s['derived'].get('reach_left') or s['derived'].get('bp_left') or 0
                right = s['derived'].get('reach_right') or s['derived'].get('bp_right') or 0
            margin = 31 + bk
            iv.setdefault(w['source'], []).append([w['start'] - left - margin, w['end'] + right + margin])
    merged = {}
    for src, xs in iv.items():
        xs.sort()
        m = []
        for a, b in xs:
            if m and a <= m[-1][1]:
                m[-1][1] = max(m[-1][1], b)
            else:
                m.append([a, b])
        merged[src] = m
    return merged


def suite_tables(records):
    """Derived per-suite numbers (also rendered in summary.md)."""
    R = by_id(records)
    T = {}
    # callsize: per call ~ c + r x (rows - 1)
    cs = {}
    for rid, rec in R.items():
        m = re.match(r'^callsize\.([^.]+)\.(k30_first|k1|k3|k30_warm)$', rid)
        if m and rec.get('http') == 200:
            cs.setdefault(m.group(1), {})[m.group(2)] = rec['resolve']['elapsed_ms']
    for d in cs.values():
        if 'k1' in d and 'k30_warm' in d:
            per_row = (d['k30_warm'] - d['k1']) / 29.0
            d['marginal_ms_per_row'] = per_row
            d['per_call_ms'] = d['k1']
            d['k30_over_k1'] = div(d['k30_warm'], d['k1'])
            d['k3_per_row'] = div(d.get('k3'), 3)
            d['k30_per_row'] = div(d['k30_warm'], 30)
    T['callsize'] = cs
    rs = {}
    for rid, rec in R.items():
        m = re.match(r'^resolve\.([^.]+)\.(explicit_first|explicit_warm|discover_column|discover_header)$', rid)
        if m and rec.get('http') == 200:
            rs.setdefault(m.group(1), {})[m.group(2)] = rec['resolve']['elapsed_ms']
    for d in rs.values():
        x = d.get('explicit_warm')
        d['discover_column_over_explicit'] = div(d.get('discover_column'), x)
        d['discover_header_over_explicit'] = div(d.get('discover_header'), x)
    T['resolve'] = rs
    fw = {}
    for rid, rec in R.items():
        m = re.match(r'^(walks\.[^.]+\.(?:named|trace)|annotate\.[^.]+|budgets\.[^.]+)\.first$', rid)
        if m:
            a, b = seed_res(rec), seed_res(R.get(m.group(1) + '.warm'))
            if a and b:
                fw[m.group(1)] = {'first_ms': a['derived']['elapsed_ms'], 'warm_ms': b['derived']['elapsed_ms'],
                                  'first_over_warm': div(a['derived']['elapsed_ms'], b['derived']['elapsed_ms']),
                                  'first_ms_per_row': a['derived']['ms_per_row'],
                                  'warm_ms_per_row': b['derived']['ms_per_row'],
                                  'first_bp': a['derived']['bp_total'], 'warm_bp': b['derived']['bp_total']}
    T['first_vs_warm'] = fw
    la = {}
    for rid, rec in R.items():
        m = re.match(r'^lookahead\.([^.]+)\.(b64_first|b64_warm|b1|b1000)$', rid)
        s = seed_res(rec)
        if m and s:
            la.setdefault(m.group(1), {})[m.group(2)] = {
                'elapsed_ms': s['derived']['elapsed_ms'], 'rows': s['derived']['rows'],
                'rows_per_base': s['derived']['rows_per_base'], 'bp_total': s['derived']['bp_total']}
    T['lookahead'] = la
    ratios = [s['derived']['budget_ratio'] for rec in records for s in (rec.get('results') or [])
              if s['derived'].get('budget_ratio') is not None and s.get('time_bound')]
    T['budget_ratio_time_bound'] = {
        'n': len(ratios), 'median': statistics.median(ratios) if ratios else None,
        'max': max(ratios) if ratios else None}
    return T


# ============================================================================ summary.md

def seed_cells(s):
    if s is None:
        return ['-'] * 5
    d = s['derived']
    st = '/'.join((a['status'] or '?')[:4] for a in s['arms'].values()) or (s.get('outcome') or {}).get('walks') or '-'
    if s.get('error'):
        st = 'failed'
    return [fmt(d['elapsed_ms']), fmt(d['ms_per_row'], 2), '%s+%s' % (fmt(d['bp_left']), fmt(d['bp_right'])),
            fmt(d['rows']), st]


def summary_md(results):
    meta, recs = results['meta'], results['requests']
    R = by_id(recs)
    L = []
    srv = results.get('server') or {}
    L.append('# Traversal benchmark: %s\n' % meta.get('label'))
    L.append('- base: `%s`' % meta.get('base'))
    L.append('- server: release `%s`, algorithm `%s`, feature_level `%s`, index_fp `%s`' % (
        srv.get('release'), srv.get('algorithm_version'), srv.get('feature_level'), srv.get('index_fp')))
    c = results.get('capabilities') or {}
    if c:
        L.append('- limits: max_seeds %s, max_time_ms %s, attempts %s' % (
            c.get('max_seeds'), c.get('max_time_ms'), c.get('attempts')))
    pnl = meta.get('panel') or {}
    note = (' (' + meta['offset_note'] + ')') if meta.get('offset_note') else ''
    L.append('- panel: %s (sha %s), window offset K=%s%s' % (pnl.get('name'), pnl.get('sha'), meta.get('window_offset'),
                                                           note))
    wrapped = sorted({w['seed'] for r in recs for w in (r.get('windows') or []) if w.get('wrapped')})
    if wrapped:
        L.append('- windows wrapped inside their source: %s' % ', '.join(wrapped))
    L.append('- suites: %s; plan version %s; requests sent %s of %s planned (cap %s)' % (
        ', '.join(meta.get('suites') or []), meta.get('plan_version'), meta.get('requests_sent'),
        meta.get('requests_planned'), meta.get('max_requests')))
    L.append('- started %s, finished %s%s' % (meta.get('started'), meta.get('finished'),
                                             ('; STOPPED: ' + meta['stop']) if meta.get('stop') else ''))
    tot_wall = sum(r.get('wall_ms') or 0 for r in recs)
    tot_srv = sum((r.get('server') or {}).get('elapsed_ms') or 0 for r in recs)
    br = results['tables'].get('budget_ratio_time_bound') or {}
    L.append('- request time: %s s wall in total, %s s server elapsed; time-bound seeds: %s, budget ratio '
             'median %s, max %s' % (fmt(tot_wall / 1000.0, 1), fmt(tot_srv / 1000.0, 1), br.get('n'),
                                    fmt(br.get('median'), 3), fmt(br.get('max'), 2)))
    inv = results.get('invariants') or []
    bad = [i for i in inv if i['verdict'] == 'DIFFERENT']
    L.append('- invariants: %d checked, %d DIFFERENT, %d time-bound (not comparable)\n' % (
        len(inv), len(bad), sum(1 for i in inv if i['verdict'] == 'time-bound')))

    # budgets
    rows = []
    for rid, rec in R.items():
        if rec['suite'] != 'budgets':
            continue
        s = seed_res(rec)
        w = (rec.get('windows') or [{}])[0]
        rows.append([rid, '%s@%s' % (w.get('seed'), w.get('start')), fmt(rec['params'].get('budget_ms')),
                     rec['params'].get('batch_kmers')] + (
            [fmt(s['derived']['elapsed_ms']), fmt(s['derived']['budget_ratio'], 2),
             (s.get('outcome') or {}).get('walks') or '-', '%s+%s' % (fmt(s['derived']['bp_left']),
                                                                    fmt(s['derived']['bp_right'])),
             fmt(s['derived']['rows']), fmt(s['timing'].get('annotation_fetch_ms')),
             (s.get('error') or '')[:70]] if s else ['-', '-', 'HTTP %s' % rec.get('http'), '-', '-', '-',
                                                    (rec.get('error_body') or rec.get('error') or '')[:70]]))
    if rows:
        L.append('## budgets: time-budget adherence on fresh windows\n')
        L.append(md_table(['id', 'seed@start', 'budget ms', 'batch', 'elapsed ms', 'ratio', 'walks', 'bp L+R',
                           'rows', 'fetch ms', 'error'], rows))
        L.append('')

    # first vs warm (walks, annotate)
    rows = []
    for rid, rec in R.items():
        m = re.match(r'^(walks|annotate)\.(.+)\.first$', rid)
        if not m:
            continue
        a, b = seed_res(rec), seed_res(R.get(rid[:-len('first')] + 'warm'))
        v, _ = same_result(a, b)
        p = rec['params']
        rows.append([m.group(2), p.get('class') or '-', p.get('mode'), p.get('radius')] + seed_cells(a)
                    + seed_cells(b) + [fmt(div(a and a['derived']['elapsed_ms'], b and b['derived']['elapsed_ms']), 1),
                                       v])
    if rows:
        L.append('## walks / annotate: first touch, then the identical warm repeat (5 s budget)\n')
        L.append('ms/row = annotation_fetch_ms / (rows_fetched + tuple_rows_fetched); bp = complete_to_bp per arm.\n')
        L.append(md_table(['request', 'class', 'mode', 'radius', 'first ms', 'first ms/row', 'first bp', 'first rows',
                           'first arms', 'warm ms', 'warm ms/row', 'warm bp', 'warm rows', 'warm arms',
                           'first/warm', 'result'], rows))
        L.append('')

    # callsize
    cs = results['tables'].get('callsize') or {}
    if cs:
        rows = [[sid, fmt(d.get('k30_first')), fmt(d.get('k1')), fmt(d.get('k3')), fmt(d.get('k30_warm')),
                 fmt(d.get('per_call_ms')), fmt(d.get('marginal_ms_per_row'), 2), fmt(d.get('k30_per_row'), 2),
                 fmt(d.get('k30_over_k1'), 2)] for sid, d in cs.items()]
        L.append('## callsize: /resolve explicit headers, support trace, 1 / 3 / 30 consecutive k-mers\n')
        L.append('A call costs about one row-diff path decode (~ the k1 time) plus a per-row increment '
                 '((k30 - k1) / 29). Path reuse across calls shows as k1 and per-call cost falling.\n')
        L.append(md_table(['seed', 'k30 first ms', 'k1 ms', 'k3 ms', 'k30 warm ms', 'per call ms',
                           'marginal ms/row', 'k30 ms/row', 'k30/k1'], rows))
        L.append('')

    rs = results['tables'].get('resolve') or {}
    if rs:
        rows = [[sid, fmt(d.get('explicit_first')), fmt(d.get('explicit_warm')), fmt(d.get('discover_column')),
                 fmt(d.get('discover_header')), fmt(d.get('discover_column_over_explicit'), 2),
                 fmt(d.get('discover_header_over_explicit'), 2)] for sid, d in rs.items()]
        L.append('## resolve: explicit labels vs discover\n')
        L.append(md_table(['seed', 'explicit first ms', 'explicit warm ms', 'discover column ms',
                           'discover header ms', 'dcol / explicit', 'dhdr / explicit'], rows))
        L.append('')

    la = results['tables'].get('lookahead') or {}
    if la:
        rows = []
        for sid, d in la.items():
            for tag in ('b64_first', 'b64_warm', 'b1', 'b1000'):
                x = d.get(tag)
                if not x:
                    continue
                v, _ = same_result(seed_res(R.get('lookahead.%s.b64_first' % sid)),
                                   seed_res(R.get('lookahead.%s.%s' % (sid, tag))))
                rows.append([sid, tag, fmt(x['elapsed_ms']), fmt(x['rows']), fmt(x['bp_total']),
                             fmt(x['rows_per_base'], 2), v if tag != 'b64_first' else '(reference)'])
        L.append('## lookahead: derived columns, r300, batch_kmers 64 / 1 / 1000\n')
        L.append('rows/base = rows fetched per distinct base output (the lookahead over-read).\n')
        L.append(md_table(['seed', 'request', 'elapsed ms', 'rows', 'bp L+R', 'rows/base', 'result vs b64_first'],
                          rows))
        L.append('')

    rows = []
    for rid, rec in R.items():
        if rec['suite'] != 'batch':
            continue
        if rec.get('http') != 200:
            rows.append([rid, '-', 'HTTP %s' % rec.get('http')] + ['-'] * 5)
            continue
        for i, s in enumerate(rec.get('results') or []):
            sid = rec['params']['seeds'][i]
            single = R.get('walks.%s.named.warm' % sid) or R.get('walks.%s.named.first' % sid)
            v = '-'
            if single and single['params'].get('radius') == rec['params'].get('radius'):
                v, _ = same_result(s, seed_res(single), key='noH')
            rows.append([rid, sid] + seed_cells(s) + [v])
        rows.append([rid, '(request)', fmt(rec['server']['elapsed_ms']), '-', '-', '-',
                     'batch ratio %s' % fmt(rec.get('batch_budget_ratio'), 2), '-'])
    if rows:
        L.append('## batch: several seeds in one request\n')
        L.append(md_table(['request', 'seed', 'elapsed ms', 'ms/row', 'bp L+R', 'rows', 'arms', 'equal to single'],
                          rows))
        L.append('')

    if inv:
        L.append('## invariants\n')
        L.append(md_table(['check', 'a', 'b', 'verdict', 'detail'],
                          [[i['check'], i['a'], i['b'], i['verdict'], i['detail']] for i in inv]))
        L.append('')

    rows = []
    for rec in recs:
        for i, s in enumerate(rec.get('results') or []):
            if 'reach' not in s:
                reach_and_work(s)
            r, w, ph = s['reach'], s['work'], s.get('physical') or {}
            rows.append([rec['id'] + ('#%d' % i if len(rec['results']) > 1 else ''), r['stop'],
                         '%s+%s' % (fmt(r['bp_left']), fmt(r['bp_right'])), fmt(r['walked_bp']), fmt(r['steps']),
                         '%s/%s' % (fmt(r['complete_walks']), fmt(r['leaves'])), fmt(w['work_units']),
                         fmt(w['usage_work_units']), fmt(w['rows_requested']), fmt(w['rows_decoded']),
                         fmt(ph.get('reuse')), fmt(ph.get('peak_cache_bytes'))])
    if rows:
        L.append('## certified reach and consumed work (per seed)\n')
        L.append('stop = what ended the exploration; bp = complete_to_bp per arm (certified); walked = the longest '
                 'walks (counts.max_bp); walks = complete walks (a leaf reason of the walk\'s own, not a resource '
                 'stop) / leaves; wu = the arms\' work_units (under a budget), usage wu = the attempt\'s usage; rows '
                 'requested is logical, rows decoded physical; reuse = timing.path_cache.hits, peak cache = '
                 'timing.path_cache.peak_bytes.\n')
        L.append(md_table(['id', 'stop', 'bp L+R', 'walked bp', 'steps', 'walks', 'wu', 'usage wu',
                           'rows requested', 'rows decoded', 'reuse', 'peak cache B'], rows))
        L.append('')

    rows = []
    for rec in recs:
        res = rec.get('results') or []
        if res:
            for i, s in enumerate(res):
                d = s['derived']
                rows.append([rec['n'], rec['id'] + ('#%d' % i if len(res) > 1 else ''), rec['http'],
                             fmt(rec['wall_ms']), fmt(d['elapsed_ms']), fmt(s['timing'].get('annotation_fetch_ms')),
                             fmt(d['rows']), fmt(d['ms_per_row'], 2), '%s+%s' % (fmt(d['bp_left']), fmt(d['bp_right'])),
                             '%s+%s' % (fmt(d['reach_left']), fmt(d['reach_right'])),
                             fmt(d['budget_ratio'], 2), (s.get('outcome') or {}).get('walks') or '-',
                             '%s/%s' % (rec.get('bytes_wire'), rec.get('bytes_json'))])
        else:
            z = rec.get('resolve') or {}
            rows.append([rec['n'], rec['id'], rec['http'], fmt(rec['wall_ms']),
                         fmt((rec.get('server') or {}).get('elapsed_ms')), '-', '-',
                         fmt(z.get('ms_per_kmer'), 2) + (' /kmer' if z else ''), '-', '-', '-',
                         (rec.get('error_body') or rec.get('error') or '')[:40] or '-',
                         '%s/%s' % (rec.get('bytes_wire'), rec.get('bytes_json'))])
    crow = []
    for rid, rec in R.items():
        for i, s in enumerate(rec.get('results') or []):
            c = s.get('coordinates')
            if not c:
                continue
            name = rec['id'] + ('#%d' % i if len(rec.get('results') or []) > 1 else '')
            if not c.get('present'):
                crow.append([name, '-', 'null: %s' % c.get('reason'), '-', '-', '-', '-', '-', '-', '-', '-'])
                continue
            crow.append([name, c['kind'], c['block_bytes'], fmt(c['share'], 3), c['seed_entries'], c['runs'],
                         '%s/%s' % (c['occurrences'], c['occurrences_total']), c['max_list'], c['lists_cut'],
                         c['runs_lower_bound'], c['complete']])
    if crow:
        L.append('## record coordinates\n')
        L.append('block = the compact JSON bytes of the coordinates block; share = of the seed result; '
                 'occurrences listed / true; max = the largest true count of a list.\n')
        L.append(md_table(['id', 'kind', 'block B', 'share', 'seed entries', 'runs', 'occurrences', 'max',
                           'lists cut', 'lower-bound runs', 'complete'], crow))
        L.append('')
    L.append('## all requests\n')
    L.append('bp = complete_to_bp per arm (certified); reach = the longest walk per arm (counts.max_bp).\n')
    L.append(md_table(['n', 'id', 'http', 'wall ms', 'server ms', 'fetch ms', 'rows', 'ms/row', 'bp L+R',
                       'reach L+R', 'budget ratio', 'walks', 'bytes wire/json'], rows))
    L.append('')
    return '\n'.join(L)


# ============================================================================ results assembly

def assemble(meta, records):
    caps_rec = next((r for r in records if r.get('capabilities')), None)
    srv = {}
    for r in records:
        for k in ('release', 'algorithm_version', 'feature_level', 'index_fp'):
            if r.get(k) is not None and srv.get(k) is None:
                srv[k] = r[k]
    results = {'format': 'metagraph-bench-traverse-results', 'version': SCRIPT_VERSION, 'meta': meta,
               'capabilities': caps_rec['capabilities'] if caps_rec else None, 'server': srv,
               'requests': records}
    results['invariants'] = invariants(records)
    results['tables'] = suite_tables(records)
    results['read_intervals'] = read_intervals(records)
    return results


def load_results(path):
    if os.path.isdir(path):
        path = os.path.join(path, 'results.json')
    with open(path) as f:
        return json.load(f)


def write_outputs(out_dir, results):
    write_json_atomic(os.path.join(out_dir, 'results.json'), results)
    with open(os.path.join(out_dir, 'summary.md'), 'w') as f:
        f.write(summary_md(results))


# ============================================================================ compare

def keyed(results):
    out = {}
    for rec in results['requests']:
        res = rec.get('results') or []
        if len(res) > 1:
            for i, s in enumerate(res):
                sid = (rec.get('params') or {}).get('seeds', [None] * len(res))[i]
                out['%s#%s' % (rec['id'], sid)] = (rec, s)
        else:
            out[rec['id']] = (rec, res[0] if res else None)
    return out


def metrics(rec, s):
    """The numbers --compare reads from one request (and seed): latency, certified reach, consumed work and the
    physical cache work. A result read from an older results.json gets its reach / work first."""
    if s is not None:
        if 'reach' not in s:
            reach_and_work(s)
        d, r, w, ph = s['derived'], s['reach'], s['work'], s.get('physical') or {}
        return {'elapsed': d['elapsed_ms'], 'ms_row': d['ms_per_row'], 'bp': d['bp_total'],
                'ratio': d['budget_ratio'], 'bp_s': d['bp_per_s'], 'rows': d['rows'], 'time_bound': s.get('time_bound'),
                'stop': r['stop'], 'steps': r['steps'], 'complete_walks': r['complete_walks'], 'leaves': r['leaves'],
                'work': work_amount(w), 'logical': {k: w.get(k) for k in COMPLETED_WORK_KEYS},
                'reuse': ph.get('reuse'), 'peak_cache': ph.get('peak_cache_bytes')}
    base = {'ms_row': None, 'bp': None, 'ratio': None, 'bp_s': None, 'rows': None, 'time_bound': False,
            'stop': None, 'steps': None, 'complete_walks': None, 'leaves': None, 'work': (None, None),
            'logical': {}, 'reuse': None, 'peak_cache': None}
    if rec.get('resolve'):
        z = rec['resolve']
        return dict(base, elapsed=z['elapsed_ms'], ms_row=z['ms_per_kmer'], rows=z['num_kmers'],
                    logical={'num_kmers': z['num_kmers']})
    return dict(base, elapsed=(rec.get('server') or {}).get('elapsed_ms'))


def work_amount(w):
    """-> (amount, unit) of a seed's consumed work: work units when the server reported them (a budget's meter or
    an attempt's usage), else the annotation rows the walk requested (logical, reported by every server)."""
    if w.get('work_units') is not None:
        return w['work_units'], 'wu'
    if w.get('usage_work_units') is not None:
        return w['usage_work_units'], 'wu'
    if w.get('rows_requested') is not None:
        return w['rows_requested'], 'rows'
    return None, None


def _work_cell(x):
    amount, unit = x
    return '-' if amount is None else '%s %s' % (fmt(amount), unit)


def _logical_diff(la, lb):
    """'same', '?' (a side does not report them) or the counters that differ, 'name a->b'."""
    keys = [k for k in la if la.get(k) is not None and lb.get(k) is not None]
    if not keys:
        return '?'
    diff = ['%s %s->%s' % (k, fmt(la[k]), fmt(lb[k])) for k in keys if la[k] != lb[k]]
    return '; '.join(diff) if diff else 'same'


def same_request(ra, rb):
    """The same request, or the same apart from the record coordinates it asks for (a --coordinates run
    against an opt-out one: the same walks, so the pair measures what the coordinates cost). A record of an
    older run has no request_sha_nocoords: its request_sha stands for it (it asked for no coordinates)."""
    if ra.get('request_sha') == rb.get('request_sha'):
        return True
    return (ra.get('request_sha_nocoords') or ra.get('request_sha')) == \
        (rb.get('request_sha_nocoords') or rb.get('request_sha'))


def _pair_class(ra, sa, rb, sb, xa, xb):
    """completed (both ran to their end: latency is comparable on identical completed work), deadline (a time
    budget, deadline or resource budget stopped either: compare certified reach and consumed work), or other
    (an HTTP error on a side, a failed seed, a request that differs)."""
    if ra.get('http') != 200 or rb.get('http') != 200 or not same_request(ra, rb):
        return 'other'
    if sa is None and sb is None:
        return 'completed'          # /resolve: no budget stops it
    if sa is None or sb is None:
        return 'other'
    if xa['time_bound'] or xb['time_bound'] or str(xa['stop']).startswith('resource:') \
            or str(xb['stop']).startswith('resource:'):
        return 'deadline'
    if sa.get('error') or sb.get('error'):
        return 'other'
    return 'completed'


def compare(path_a, path_b):
    A, B = load_for_compare(path_a), load_for_compare(path_b)
    ka, kb = keyed(A), keyed(B)
    L = []
    ma, mb = A['meta'], B['meta']
    L.append('# Benchmark comparison: %s -> %s\n' % (ma.get('label'), mb.get('label')))
    for tag, m, r in (('A', ma, A), ('B', mb, B)):
        s = r.get('server') or {}
        L.append('- %s: %s, release `%s`, feature_level %s, panel %s (%s), K=%s, started %s, %s requests%s' % (
            tag, m.get('label'), s.get('release'), s.get('feature_level'), (m.get('panel') or {}).get('name'),
            (m.get('panel') or {}).get('sha'), m.get('window_offset'), m.get('started'), m.get('requests_sent'),
            ' (analysed again from its raw/ responses: an older results.json)' if r.get('_rederived') else ''))
    if ma.get('plan_version') != mb.get('plan_version'):
        L.append('- WARNING: plan versions differ (%s vs %s): ids may not mean the same request' % (
            ma.get('plan_version'), mb.get('plan_version')))
    L.append('')
    L.append('Completed requests (both ran to their end) compare LATENCY, and only on identical completed work: '
             'the same result and the same logical work counters (%s; /resolve: num_kmers). Deadline-limited '
             'requests (a time budget, deadline or resource budget stopped either side) compare CERTIFIED REACH '
             '(complete_to_bp of both arms, steps, complete walks, the stop) and CONSUMED WORK (work units where '
             'the server reports them, else the annotation rows requested): their elapsed is about the budget '
             'on both sides and hides any gain. reuse = paths cut at a cached row (timing.path_cache.hits), '
             'peak cache = timing.path_cache.peak_bytes; "-" where a server reports none.\n'
             % ', '.join(COMPLETED_WORK_KEYS))
    done_rows, dl_rows, other_rows, flags, work_flags = [], [], [], [], []
    per_suite = {}
    for key in [k for k in ka if k in kb]:
        (ra, sa), (rb, sb) = ka[key], kb[key]
        if ra.get('route') == '/traverse/capabilities':
            continue
        xa, xb = metrics(ra, sa), metrics(rb, sb)
        if ra.get('http') != 200 or rb.get('http') != 200:
            verdict = 'http %s/%s' % (ra.get('http'), rb.get('http'))
        elif not same_request(ra, rb):
            verdict = 'request differs'
        elif sa is not None or sb is not None:
            verdict, _ = same_result(sa, sb)
        else:
            verdict = 'equal' if (ra.get('resolve') or {}).get('result_sha') == (rb.get('resolve') or {}).get(
                'result_sha') else 'DIFFERENT'
        if verdict == 'DIFFERENT':
            flags.append(key)
        cls = _pair_class(ra, sa, rb, sb, xa, xb)
        ps = per_suite.setdefault(ra['suite'], {'lat': [], 'excluded': 0, 'reach': [], 'work': [], 'further': 0,
                                                'same': 0, 'less': 0, 'deadline': 0})
        ratio = div(xb['elapsed'], xa['elapsed'])
        if cls == 'completed':
            logical = _logical_diff(xa['logical'], xb['logical'])
            identical = verdict in ('equal', 'walks equal') and logical in ('same', '?')
            if identical and ratio:
                ps['lat'].append(ratio)
            else:
                ps['excluded'] += 1
                if verdict in ('equal', 'walks equal') and logical not in ('same', '?'):
                    work_flags.append('%s (%s)' % (key, logical))
            done_rows.append([key, fmt(xa['elapsed']), fmt(xb['elapsed']), fmt(ratio, 2) if identical else
                              '(%s)' % fmt(ratio, 2), fmt(xa['ms_row'], 2), fmt(xb['ms_row'], 2), logical,
                              fmt(xa['reuse']), fmt(xb['reuse']), fmt(xa['peak_cache']), fmt(xb['peak_cache']),
                              verdict])
        elif cls == 'deadline':
            ps['deadline'] += 1
            bpa, bpb = xa['bp'], xb['bp']
            if bpa is not None and bpb is not None:
                ps['further' if bpb > bpa else 'less' if bpb < bpa else 'same'] += 1
                if bpa > 0 and bpb > 0:
                    ps['reach'].append(bpb / bpa)
            (wa, ua), (wb, ub) = xa['work'], xb['work']
            wr = div(wb, wa) if ua == ub and wa else None
            if wr:
                ps['work'].append(wr)
            dl_rows.append([key, fmt(xa['elapsed']), fmt(xb['elapsed']), xa['stop'] or '-', xb['stop'] or '-',
                            fmt(bpa), fmt(bpb), fmt(div(bpb, bpa), 2), fmt(xa['steps']), fmt(xb['steps']),
                            '%s/%s' % (fmt(xa['complete_walks']), fmt(xa['leaves'])),
                            '%s/%s' % (fmt(xb['complete_walks']), fmt(xb['leaves'])),
                            _work_cell(xa['work']), _work_cell(xb['work']), fmt(wr, 2),
                            fmt(xa['reuse']), fmt(xb['reuse']), fmt(xa['peak_cache']), fmt(xb['peak_cache']),
                            verdict])
        else:
            other_rows.append([key, fmt(xa['elapsed']), fmt(xb['elapsed']), xa['stop'] or '-', xb['stop'] or '-',
                               verdict])
    only_a = [k for k in ka if k not in kb]
    only_b = [k for k in kb if k not in ka]
    coord_rows = []
    for key in [k for k in ka if k in kb]:
        (ra, sa), (rb, sb) = ka[key], kb[key]
        ca, cb = (sa or {}).get('coordinates'), (sb or {}).get('coordinates')
        if not ca and not cb:
            continue
        xa, xb = metrics(ra, sa), metrics(rb, sb)

        def cell(c, f, nd=None):
            if not c:
                return '-'
            if not c.get('present'):
                return 'null' if f == 'block_bytes' else '-'
            return fmt(c.get(f), nd) if nd is not None else c.get(f)
        coord_rows.append([key, cell(ca, 'block_bytes'), cell(cb, 'block_bytes'), cell(ca, 'share', 3),
                           cell(cb, 'share', 3), cell(ca, 'occurrences'), cell(cb, 'occurrences'),
                           cell(ca, 'lists_cut'), cell(cb, 'lists_cut'), fmt(xa['elapsed']),
                           fmt(xb['elapsed']), fmt(div(xb['elapsed'], xa['elapsed']), 2)])
    if coord_rows:
        L.append('## record coordinates\n')
        L.append('block = the coordinates block\'s compact JSON bytes ("-": not asked for, "null": asked for, '
                 'none reported); share = of the seed result; elapsed B/A is the cost of the coordinates where '
                 'only B asked for them (the same walks: see the completed table).\n')
        L.append(md_table(['id', 'A block B', 'B block B', 'A share', 'B share', 'A occurrences',
                           'B occurrences', 'A cut', 'B cut', 'A ms', 'B ms', 'B/A'], coord_rows))
        L.append('')
    L.append('Per suite. Latency: geometric mean of B/A elapsed over the completed pairs with identical work '
             '(< 1 = B faster; "excluded": completed pairs whose result or logical work differ). Reach: geometric '
             'mean of B/A certified bp over the deadline-limited pairs (> 1 = B reached further; pairs at 0 bp '
             'only counted), B further / same / less, and the geometric mean of B/A consumed work.\n')
    L.append(md_table(['suite', 'completed', 'latency B/A', 'excluded', 'deadline-limited', 'reach B/A',
                       'further/same/less', 'work B/A'],
                      [[s, len(v['lat']), fmt(geomean(v['lat']), 3), v['excluded'], v['deadline'],
                        fmt(geomean(v['reach']), 3), '%d/%d/%d' % (v['further'], v['same'], v['less']),
                        fmt(geomean(v['work']), 3)] for s, v in per_suite.items()]))
    L.append('')
    if done_rows:
        L.append('## completed: latency on identical completed work\n')
        L.append('B/A in parentheses: not comparable (the result or the logical work differ).\n')
        L.append(md_table(['id', 'A ms', 'B ms', 'B/A', 'A ms/row', 'B ms/row', 'logical work', 'A reuse',
                           'B reuse', 'A peak cache B', 'B peak cache B', 'result'], done_rows))
        L.append('')
    if dl_rows:
        L.append('## deadline-limited: certified reach and consumed work\n')
        L.append('bp = complete_to_bp of both arms; walks = complete walks / leaves; work = work units (wu) or '
                 'annotation rows requested (rows); work B/A only in one unit.\n')
        L.append(md_table(['id', 'A ms', 'B ms', 'A stop', 'B stop', 'A bp', 'B bp', 'bp B/A', 'A steps',
                           'B steps', 'A walks', 'B walks', 'A work', 'B work', 'work B/A', 'A reuse', 'B reuse',
                           'A peak cache B', 'B peak cache B', 'result'], dl_rows))
        L.append('')
    if other_rows:
        L.append('## not compared\n')
        L.append(md_table(['id', 'A ms', 'B ms', 'A stop', 'B stop', 'result'], other_rows))
        L.append('')
    if only_a or only_b:
        L.append('- only in A: %s' % (', '.join(only_a) or '-'))
        L.append('- only in B: %s' % (', '.join(only_b) or '-'))
    if work_flags:
        L.append('\n**Completed with equal results but different logical work: %s** -- the server\'s counting '
                 'changed; their latency is not compared.' % ', '.join(work_flags))
    if flags:
        L.append('\n**DIFFERENT deterministic results (same request): %s** -- results must be equal unless the '
                 'server\'s behaviour changed on purpose.' % ', '.join(flags))
    else:
        L.append('\nNo deterministic result changed (rows marked "request differs" or "time-bound" were not compared).')
    return '\n'.join(L), flags


def load_for_compare(path):
    """A results.json for --compare. One written before SCRIPT_VERSION 2 lacks the arms' counters and
    leaves_by_reason (steps, complete walks) and the physical fields: when its raw/ responses are next to it,
    its requests are analysed again from them, in memory (nothing is written); otherwise reach and work are
    filled from what it recorded (the missing numbers stay None)."""
    R = load_results(path)
    run_dir = path if os.path.isdir(path) else os.path.dirname(os.path.abspath(path))
    old = any('reach' not in s for rec in R.get('requests') or [] for s in rec.get('results') or [])
    raws = sorted(glob.glob(os.path.join(run_dir, 'raw', '[0-9]*.json.gz')))
    if old and raws:
        records = []
        for p in raws:
            with gzip.open(p, 'rt') as f:
                e = json.load(f)
            records.append(analyze(dict(e['meta'], raw_file='raw/' + os.path.basename(p)), e['request'],
                                   e['response']))
        R = assemble(R.get('meta') or {}, records)
        R['_rederived'] = True
    return R


# ============================================================================ main

def reanalyze(out_dir):
    results_path = os.path.join(out_dir, 'results.json')
    meta = {}
    if os.path.exists(results_path):
        meta = load_results(results_path).get('meta', {})
    records = []
    for p in sorted(glob.glob(os.path.join(out_dir, 'raw', '[0-9]*.json.gz'))):
        with gzip.open(p, 'rt') as f:
            e = json.load(f)
        m = dict(e['meta'], raw_file='raw/' + os.path.basename(p))
        records.append(analyze(m, e['request'], e['response']))
    meta['reanalyzed'] = now_iso()
    meta['requests_sent'] = len(records)
    results = assemble(meta, records)
    write_outputs(out_dir, results)
    print('rewrote %s/results.json and summary.md from %d raw responses' % (out_dir, len(records)))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog='See the module docstring (python3 bench_traverse.py --help-long) for the '
                                        'request set, the pacing rules and the outputs.')
    ap.add_argument('--base', help='server base URL, e.g. https://metagraph.ethz.ch:8080/api/refseq33m-experimental '
                                   'or http://127.0.0.1:5871')
    ap.add_argument('--panel', default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                    'bench_panel_refseq33m.json'))
    ap.add_argument('--suite', default='all', help='comma list: all | %s' % ' | '.join(ALL_SUITES + EXTRA_SUITES))
    ap.add_argument('--only', help='regex on request ids (applied after --suite)')
    ap.add_argument('--max-requests', type=int, default=70, help='hard cap per run (failed requests count)')
    ap.add_argument('--pause-s', type=float, default=3.0,
                    help='pause between requests (min 2 s on a non-local server; 6 s after one over 3 s)')
    ap.add_argument('--timeout-s', type=float, default=180.0)
    ap.add_argument('--out', help='run directory (default bench_out/<label>_<timestamp>)')
    ap.add_argument('--label', default='run')
    ap.add_argument('--window-offset', default='0', help='K bp (int) or "auto" (with --avoid)')
    ap.add_argument('--avoid', nargs='*', default=[], help='results.json of earlier runs whose read intervals '
                                                           '--window-offset auto avoids (and warns about)')
    ap.add_argument('--dry-run', action='store_true', help='print the request plan, send nothing')
    ap.add_argument('--verbose', action='store_true', help='with --dry-run: print every request body')
    ap.add_argument('--compare', nargs=2, metavar=('A', 'B'), help='before/after table of two results.json')
    ap.add_argument('--compare-out', help='also write the comparison markdown here')
    ap.add_argument('--reanalyze', metavar='RUN_DIR', help='rebuild results.json + summary.md from RUN_DIR/raw')
    ap.add_argument('--help-long', action='store_true', help='print the full documentation')
    ap.add_argument('--coordinates', action='store_true',
                    help='support-trace walks ask for record coordinates (feature level 6)')
    ap.add_argument('--coordinate-cap', default=str(COORDINATES_CAP_DEFAULT),
                    help='with --coordinates: max_coordinate_occurrences, an int >= 1 or "unlimited" '
                         '(default %d)' % COORDINATES_CAP_DEFAULT)
    ap.add_argument('--coordinates-estimate', nargs='+', metavar='RUN_DIR',
                    help='offline: the coordinates block size the runs\' recorded trace walks would carry')
    ap.add_argument('--coordinate-digits', type=int, default=7,
                    help='with --coordinates-estimate: digits per position (default 7)')
    args = ap.parse_args(argv)

    if args.help_long:
        print(__doc__)
        return 0
    if args.compare:
        text, flags = compare(*args.compare)
        print(text)
        if args.compare_out:
            with open(args.compare_out, 'w') as f:
                f.write(text + '\n')
        return 1 if flags else 0
    if args.reanalyze:
        reanalyze(args.reanalyze)
        return 0
    if args.coordinates_estimate:
        coordinates_estimate(args.coordinates_estimate, args.coordinate_digits)
        return 0
    cap = None
    if args.coordinates:
        cap = args.coordinate_cap
        if cap != 'unlimited':
            try:
                cap = int(cap)
            except ValueError:
                ap.error('--coordinate-cap is an integer >= 1 or "unlimited"')
            if cap < 1:
                ap.error('--coordinate-cap is an integer >= 1 or "unlimited"')
    if not args.base and not args.dry_run:
        ap.error('--base is required (except with --dry-run, --compare, --reanalyze)')

    panel = Panel(args.panel)
    suites = parse_suites(args.suite)
    offset_note = None
    if args.window_offset == 'auto':
        K, overlaps, offset_note = choose_offset(panel, suites, args.avoid)
    else:
        K = int(args.window_offset)
    reqs = Planner(panel, K, cap).build(suites)
    if args.only:
        rx = re.compile(args.only)
        reqs = [r for r in reqs if rx.search(r['id'])]
    if args.avoid and args.window_offset != 'auto':
        avoid, _ = avoided_intervals(panel, args.avoid)
        hits = sorted({w['seed'] for w in first_touch_windows(reqs)
                       for a, b in avoid.get(w['source'], []) if w['start'] - 100 < b and a < w['end'] + 100})
        if hits:
            print('warning: first-touch windows overlap rows read by an avoided run: %s' % ', '.join(hits),
                  file=sys.stderr)

    if args.dry_run:
        print('panel %s (sha %s), K=%d%s, suites %s: %d requests (cap %d)' % (
            panel.name, panel.sha, K, (' [' + offset_note + ']') if offset_note else '', ','.join(suites),
            len(reqs), args.max_requests))
        for i, r in enumerate(reqs, 1):
            p = r['params']
            ws = ' '.join('%s@%d%s' % (w['seed'], w['start'], '(wrapped)' if w['wrapped'] else '')
                          for w in r['windows'])
            desc = ' '.join('%s=%s' % (k, v) for k, v in p.items() if k not in ('seeds', 'seed', 'labels'))
            lab = p.get('labels') or (r['body'] or {}).get('labels')
            print('%3d %s%-44s %-22s %s %s%s' % (i, '*' if i > args.max_requests else ' ', r['id'], r['route'], ws,
                                                  desc, (' labels=%s' % lab) if lab else ''))
            if args.verbose and r['body'] is not None:
                print('      ' + json.dumps(r['body']))
        if len(reqs) > args.max_requests:
            print('(* beyond the cap: not sent)')
        local = bool(args.base) and is_local(args.base)
        pause = args.pause_s if local else max(2.0, args.pause_s)
        print('pauses alone: >= %.0f s' % (pause * max(0, min(len(reqs), args.max_requests) - 1)))
        return 0

    stamp = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    out_dir = args.out or os.path.join('bench_out', '%s_%s' % (re.sub(r'[^A-Za-z0-9._-]+', '_', args.label), stamp))
    if os.path.exists(os.path.join(out_dir, 'results.json')):
        raise SystemExit('%s already holds a results.json: choose another --out' % out_dir)
    os.makedirs(out_dir, exist_ok=True)
    meta = {'label': args.label, 'base': args.base, 'script_version': SCRIPT_VERSION, 'plan_version': PLAN_VERSION,
            'panel': {'name': panel.name, 'path': panel.path, 'sha': panel.sha}, 'window_offset': K,
            'offset_note': offset_note, 'suites': suites, 'only': args.only, 'max_requests': args.max_requests,
            'pause_s': args.pause_s, 'requests_planned': len(reqs), 'started': now_iso(),
            'coordinates': cap,
            'python': sys.version.split()[0], 'platform': platform.platform(), 'argv': sys.argv[1:]}
    runner = Runner(args, out_dir)
    records = []
    progress = open(os.path.join(out_dir, 'progress.jsonl'), 'a')

    def on_result(m, body, resp):
        rec = analyze(m, body, resp)
        records.append(rec)
        progress.write(json.dumps(rec, default=str) + '\n')
        progress.flush()
        return rec

    print('bench_traverse: %s -> %s (%d requests planned, cap %d, K=%d)' % (
        panel.name, args.base, len(reqs), args.max_requests, K), flush=True)
    runner.lock()
    stop = None
    try:
        stop = runner.run(reqs, on_result)
    except KeyboardInterrupt:
        stop = 'interrupted'
    finally:
        runner.unlock()
        progress.close()
        meta.update(finished=now_iso(), requests_sent=runner.sent, stop=stop)
        results = assemble(meta, records)
        write_outputs(out_dir, results)
    if stop:
        print('STOP: %s' % stop)
    inv = results['invariants']
    bad = [i for i in inv if i['verdict'] == 'DIFFERENT']
    print('wrote %s/{results.json,summary.md,raw/}: %d requests, invariants %d checked, %d DIFFERENT' % (
        out_dir, runner.sent, len(inv), len(bad)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
