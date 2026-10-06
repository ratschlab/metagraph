"""Conformance of every cached real retrieval (design §2.3, §2.5, §9 T37; spec §7.5).

For every cell of the retrieval cache (realdata.py: SRA, UHGG, mini_refseq; each cell
fetched once with detail "graphlet" and once with detail "full"), and for every seed
result in it:

  roundtrip   parse(text).dump() == text byte-exact; is_canonical(text); the body-only
              dump; graphlet_bytes / graphlet_lines of the summary; the Z line
  counts      the A record's six body counts recounted from the TEXT (independent of
              the parser) == A == the JSON summary's counts; max_bp, label_ends and
              leaves_by_reason of the summary == what the body and the full JSON say
  rules       check_rules(reference=full result) reports nothing (star mismatches,
              non-canonical explicit fields, §4/§9 invariants, the v5.2 outcome rule)
  to_json     Graphlet.from_response(...).to_json() == the cached detail-full result
              after the §2.5/§9 normalisation (events sorted within equal at_bp,
              needed_budgets sorted, timing removed), written independently of the
              library's own normalize_result
  invariants  from the full JSON alone: #label_end events <= #runs per arm; every run
              entered by switch has a prev_run with that label; paths/end_labels
              consistent with runs
  summary     the per-seed JSON summary agrees with the body (O, K, A, S, H) and with
              the detail-full result (seed, outcome, limitations, per-arm fields)
  identity    H index_ns / index_fp / index_meta_fp == the server's capabilities (the
              manifest's record of them); mini_refseq verifiable, SRA/UHGG not
  seed        S == the request's seed (upper case), num_kmers == L - k + 1
  save_load   save() -> load(): the .mgt file round-trips (J line, body, to_json)

A derivation failure ({seed, error}, batch_refuse) is checked to carry no body and the
same error in both detail levels.

Tests are generated per (cell, check) from the cache manifest at import, one class per
index, so a failure names its cell. Without a cache they all skip. Non-deterministic
cells (time_tight: the graphlet and the full request ran separately and stopped at
different points) skip the checks that compare the two fetches. Heavy cells (full JSON
> 64 MB) run the body-only checks; their full-JSON comparisons skip unless
METAGRAPH_REAL_MAX_FULL_MB is raised above their size (default 40 MB: loading the 403 MB
SRA full JSON costs several GB of RAM).

    cd api/python/tests/real && python3 -m unittest -v test_real_conformance
"""

import collections
import os
import re
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import realdata as R  # noqa: E402

from metagraph.traverse import (  # noqa: E402
    GraphletFormatError, dump, from_response, is_canonical, load, parse,
)
from metagraph.traverse._codec import REASON  # noqa: E402
from metagraph.traverse.export import (  # noqa: E402
    limitation_json, normalize_result, value_json,
)

MAX_FULL_BYTES = int(float(os.environ.get('METAGRAPH_REAL_MAX_FULL_MB', '40')) * (1 << 20))

# (cell, check) -> reason: a PRODUCT bug this suite found; the test is kept and marked
# unittest.expectedFailure with the reason below (see the module's product-bug notes)
KNOWN_BUGS = {}

CHECKS = ('counts', 'identity', 'invariants', 'roundtrip', 'rules', 'save_load', 'seed',
          'summary', 'to_json')


# ---------------------------------------------------------------- independent text reads

def text_lines(text):
    assert text.endswith('\n')
    return text.split('\n')[:-1]


def text_arm_counts(text):
    """Per arm, the six A counts recounted from the records alone (no parser code):
    segments = G lines, runs = R lines, leaves = T lines, merges = G with > 1 parent,
    splits = distinct parents of single-parent G records, bases = sum of G lengths;
    plus max_bp = the largest end of a leaf segment, and the A record's own fields."""
    arms = collections.OrderedDict()
    cur = None
    last_g = None
    for line in text_lines(text):
        tag = line[:1]
        if tag == 'A':
            f = line.split(' ')
            cur = {'A': [int(x) for x in f[14:20]], 'segments': 0, 'runs': 0, 'leaves': 0,
                   'merges': 0, 'split_parents': set(), 'bases': 0, 'max_bp': 0,
                   'leaf_reasons': collections.Counter(), 'run_ends': collections.Counter(),
                   'G': []}
            arms[{'l': 'left', 'r': 'right'}[f[1]]] = cur
            last_g = None
        elif cur is None:
            continue
        elif tag == 'G':
            f = line.split(' ')
            parents = [] if f[1] == '*' else [int(x) for x in f[1].split(',')]
            from_bp, length = int(f[2]), int(f[3])
            cur['segments'] += 1
            cur['bases'] += length
            if len(parents) > 1:
                cur['merges'] += 1
            elif len(parents) == 1:
                cur['split_parents'].add(parents[0])
            last_g = (from_bp, length)
            cur['G'].append(last_g)
        elif tag == 'T':
            f = line.split(' ')
            cur['leaves'] += 1
            cur['max_bp'] = max(cur['max_bp'], last_g[0] + last_g[1])
            cur['leaf_reasons'][f[1]] += 1
        elif tag == 'R':
            f = line.split(' ')
            cur['runs'] += 1
            cur['run_ends'][f[5]] += 1
        elif tag == 'Z':
            cur = None
    for a in arms.values():
        a['splits'] = len(a.pop('split_parents'))
    return arms


def z_value(text):
    last = text_lines(text)[-1]
    m = re.match(r'Z (\d+)\Z', last)
    return int(m.group(1)) if m else None


def arm_head_json(g, a):
    """The per-arm summary fields as the body states them (A, K), built from the model
    without to_json() (which would materialise every path of a heavy cell)."""
    j = {'status': a.status, 'complete_to_bp': a.complete_to_bp,
         'completeness_scope': a.scope,
         'frontier_remaining': {'live_paths': a.frontier.live_paths,
                                'live_labels': a.frontier.live_labels,
                                'exact': a.frontier.exact},
         'labels_per_node': {'cap': g.cap, 'max_seen': a.labels_per_node.max_seen,
                             'nodes_truncated': a.labels_per_node.nodes_truncated},
         'counters': dict(a.counters),
         'evidence': {'complete': a.evidence_complete_to_bp is None,
                      'complete_to_bp': a.evidence_complete_to_bp},
         'limitations': [limitation_json(l) for l in g.limitations if l.arm == a.side]}
    c = a.cap_trigger
    if c is not None:
        j['cap_trigger'] = {'reason': REASON[c.reason], 'at_bp': c.at_bp,
                            'segment': c.segment, 'live_paths': c.live_paths,
                            'live_labels': c.live_labels, 'exact': c.exact}
    return j


def full_label_end_events(arm_json):
    return [e for s in arm_json.get('segments') or [] for e in s.get('events') or []
            if e.get('type') == 'label_end']


# ---------------------------------------------------------------- a one-cell cache

class _Loaded:
    """The cell being checked: the CachedCell (its JSON memoised) and its parsed
    graphlets. Checks of one cell run consecutively (method names sort by cell), so one
    entry is enough and memory stays at one cell."""

    key = None
    cell = None
    graphlets = None

    @classmethod
    def get(cls, row):
        key = (row['index'], row['cell'])
        if cls.key != key:
            cls.key, cls.cell, cls.graphlets = None, None, None
            cell = R.load_cell(row['index'], row)
            gs = []
            for i in range(cell.n_results()):
                gs.append(cell.graphlet(i))
            cls.key, cls.cell, cls.graphlets = key, cell, gs
        return cls.cell, cls.graphlets


def _full_ok(row):
    fb = (row.get('full') or {}).get('json_bytes') or 0
    return fb <= MAX_FULL_BYTES


# ---------------------------------------------------------------- the checks

class _ConformanceBase(unittest.TestCase):
    maxDiff = 4000

    def cell(self, row):
        R.require_cache(row['index'], row['cell'])
        return _Loaded.get(row)

    def need_full(self, row, what):
        if not row.get('deterministic'):
            self.skipTest('%s: not deterministic (the graphlet and the full fetch stopped '
                          'at different points): %s not comparable' % (row['cell'], what))
        if not _full_ok(row):
            self.skipTest('%s: full JSON %d MB > METAGRAPH_REAL_MAX_FULL_MB' % (
                row['cell'], ((row.get('full') or {}).get('json_bytes') or 0) >> 20))

    def results(self, row):
        """[(i, graphlet result, graphlet or None)] of the cell."""
        cell, gs = self.cell(row)
        return [(i, cell.graphlet_result(i), gs[i]) for i in range(cell.n_results())]

    # ------------------------------------------------------------ roundtrip

    def check_roundtrip(self, row):
        cell, gs = self.cell(row)
        n = 0
        for i, res, g in self.results(row):
            text = res.get('graphlet')
            if text is None:
                self.assertIn('error', res, (row['cell'], i))
                with self.assertRaises(ValueError):
                    from_response(res, cell.graphlet_response)
                continue
            n += 1
            with self.subTest(result=i):
                self.assertEqual(res['graphlet_bytes'], len(text.encode('utf-8')))
                lines = text_lines(text)
                self.assertEqual(res['graphlet_lines'], len(lines))
                self.assertEqual(z_value(text), len(lines), 'Z line')
                self.assertTrue(lines[0].startswith('H mgt 1 '), lines[0][:40])
                g2 = parse(text)
                out = g2.dump()
                if out != text:
                    a, b = text_lines(text), text_lines(out)
                    first = next((j for j, (x, y) in enumerate(zip(a, b)) if x != y),
                                 min(len(a), len(b)))
                    self.fail('dump(parse(text)) != text at line %d:\n  %r\n  %r' % (
                        first + 1, a[first][:200] if first < len(a) else None,
                        b[first][:200] if first < len(b) else None))
                self.assertEqual(dump(g2), text)          # the module-level writer too
                self.assertTrue(is_canonical(text))
                # the graphlet built from the response: body-only dump equals the body,
                # with the envelope a J line follows H
                self.assertEqual(g.dump(envelope=False), text)
                withj = text_lines(g.dump(envelope=True))
                self.assertTrue(withj[1].startswith('J {'), withj[1][:20])
                self.assertEqual(withj[0], lines[0])
                self.assertEqual(withj[2:-1], lines[1:-1])
                self.assertEqual(withj[-1], 'Z %d' % len(withj))
        self.assertEqual(n, sum(1 for g in gs if g is not None))

    # ------------------------------------------------------------ counts

    def check_counts(self, row):
        cell, _ = self.cell(row)
        full = cell.full_response if (row.get('deterministic') and _full_ok(row)) else None
        for i, res, g in self.results(row):
            if g is None:
                continue
            text = res['graphlet']
            mine = text_arm_counts(text)
            self.assertEqual(list(mine), [s for s in ('left', 'right') if s in g.arms])
            for side, c in mine.items():
                with self.subTest(result=i, arm=side):
                    six = [c['segments'], c['runs'], c['leaves'], c['splits'], c['merges'],
                           c['bases']]
                    self.assertEqual(c['A'], six, 'A counts vs the records')
                    sc = res['arms'][side]['counts']
                    self.assertEqual(
                        [sc['segments'], sc['runs'], sc['leaves'], sc['splits'],
                         sc['merges'], sc['bases']], six, 'summary counts vs the records')
                    self.assertEqual(sc['max_bp'], c['max_bp'], 'max_bp')
                    # leaves_by_reason: path_reason name per T, '*' = semantic
                    lbr = collections.Counter()
                    for code, nn in c['leaf_reasons'].items():
                        lbr['semantic' if code == '*' else REASON[code]] += nn
                    self.assertEqual(sc['leaves_by_reason'], dict(lbr), 'leaves_by_reason')
                    # label_ends: every ended run by enum name -- including the silent Lw
                    # switch-source ends, which growth[].label_ends counts too (design
                    # §0.4) -- but not the merge-closed m runs
                    a = g.arms[side]
                    le = collections.Counter(r.reason for r in a.runs if not r.merged)
                    self.assertEqual(sc['label_ends'], dict(le), 'label_ends')
                    self.assertEqual(sum(le.values()), sum(n for b in a.growth
                                                           for n in b.label_ends.values())
                                     if a.growth else sum(le.values()))
                    n_silent = sum(1 for r in a.runs if r.silent)
                    if full is not None:
                        fa = full['results'][i]['arms'][side]
                        self.assertEqual(len(fa['segments']), c['segments'])
                        self.assertEqual(len(fa['runs']), c['runs'])
                        self.assertEqual(len(fa['paths']), c['leaves'])
                        self.assertEqual(len(fa['splits']), c['splits'])
                        self.assertEqual(sum(1 for s in fa['segments']
                                             if len(s['parents']) > 1), c['merges'])
                        self.assertEqual(sum(s['length_bp'] for s in fa['segments']),
                                         c['bases'])
                        self.assertEqual(max((p['length_bp'] for p in fa['paths']),
                                             default=0), c['max_bp'])
                        fl = collections.Counter(p.get('path_reason', 'semantic')
                                                 for p in fa['paths'])
                        self.assertEqual(dict(fl), sc['leaves_by_reason'])
                        fe = collections.Counter(e['reason'] for e in
                                                 full_label_end_events(fa))
                        # label_end events carry the qualifier text when there is one;
                        # the counts carry the enum name, and a silent Lw end has no
                        # event: totals differ by the silent ends exactly
                        self.assertEqual(sum(fe.values()),
                                         sum(sc['label_ends'].values()) - n_silent)

    # ------------------------------------------------------------ rules

    def check_rules(self, row):
        cell, _ = self.cell(row)
        use_ref = row.get('deterministic') and _full_ok(row)
        for i, res, g in self.results(row):
            if g is None:
                continue
            with self.subTest(result=i):
                ref = cell.full_result(i) if use_ref else None
                found = g.check_rules(ref)
                self.assertEqual(found, [], 'check_rules: %s' % found[:5])

    # ------------------------------------------------------------ to_json

    def check_to_json(self, row):
        self.need_full(row, 'to_json')
        cell, _ = self.cell(row)
        for i, res, g in self.results(row):
            full = cell.full_result(i)
            if g is None:
                self.assertIn('error', full, (row['cell'], i))
                self.assertEqual(full['error'], res['error'])
                continue
            with self.subTest(result=i):
                j = g.to_json()
                mine = R.normalize_full_result(j)
                theirs = R.normalize_full_result(full)
                if mine != theirs:
                    self.fail('to_json() != detail full (normalised): %s'
                              % R.json_diff(mine, theirs, limit=12))
                # the library's own normaliser (the one T37 uses) reaches the same
                # verdict (its tie-break among equal at_bp differs from the harness's,
                # so each normaliser is compared with itself)
                self.assertTrue(normalize_result(j) == normalize_result(full),
                                'equal under the harness normalisation, not under '
                                'export.normalize_result')

    # ------------------------------------------------------------ invariants

    def check_invariants(self, row):
        self.need_full(row, 'full-JSON invariants')
        cell, _ = self.cell(row)
        for i, res, g in self.results(row):
            if g is None:
                continue
            full = cell.full_result(i)
            for side, fa in full['arms'].items():
                with self.subTest(result=i, arm=side):
                    runs = fa['runs']
                    self.assertLessEqual(len(full_label_end_events(fa)), len(runs),
                                         '#label_end events <= #runs')
                    for j, r in enumerate(runs):
                        self.assertLessEqual(r['from_bp'], r['to_bp'], ('run', j))
                        seg = fa['segments'][r['segment']]
                        self.assertLessEqual(r['to_bp'], seg['from_bp'] + seg['length_bp'],
                                             ('run past its anchor', j))
                        self.assertGreaterEqual(r['to_bp'], seg['from_bp'],
                                                ('run before its anchor', j))
                        if r['entered_by'] == 'switch':
                            self.assertIn('from', r, j)
                    # end_labels of a path are runs anchored at its leaf ending there
                    for p in fa['paths']:
                        leaf = p['segments'][-1]
                        end = fa['segments'][leaf]['from_bp'] + fa['segments'][leaf]['length_bp']
                        self.assertEqual(p['length_bp'], end)
                        self.assertEqual(p['n_labels'], len(p['end_labels']))
                        for e in p['end_labels']:
                            r = runs[e['run']]
                            self.assertEqual(r['label'], e['label'])
                            self.assertEqual(r['segment'], leaf)
                            self.assertEqual(r['to_bp'], end)
                        if full['label_mode'] == 'constrain':
                            self.assertEqual(sorted(e['label'] for e in p['end_labels']),
                                             sorted(fa['segments'][leaf]['labels_at_end']))
                        for a, b in zip(p['segments'], p['segments'][1:]):
                            self.assertEqual(fa['segments'][b]['parents'][0], a)

    # ------------------------------------------------------------ summary

    def check_summary(self, row):
        cell, _ = self.cell(row)
        resp = cell.graphlet_response
        self.assertEqual(resp['capabilities'].get('graphlet_format'), 1)
        self.assertIn('graphlet', resp['capabilities'].get('detail_levels') or [])
        self.assertEqual(resp['strategy']['output']['detail'], 'graphlet')
        self.assertIs(resp['strategy']['output']['timing'], False)
        cmp_full = row.get('deterministic') and _full_ok(row)
        for i, res, g in self.results(row):
            if g is None:
                continue
            with self.subTest(result=i):
                self.assertEqual(res['outcome'], g.outcome.as_dict(), 'O vs outcome')
                self.assertEqual(res['label_mode'], g.mode)
                self.assertEqual(res['seed']['validated_seed_id'], g.seed.validated_seed_id)
                self.assertEqual(res['seed']['length_bp'], g.seed.length_bp)
                self.assertEqual(res['seed']['num_kmers'], g.seed.num_kmers)
                self.assertEqual(res['seed']['num_seed_labels'], g.seed.num_seed_labels)
                self.assertEqual(res['seed']['num_labels'], len(g.labels))
                self.assertEqual(res['limitations'], [limitation_json(l) for l in
                                                      g.limitations if l.arm is None],
                                 'seed-level K records')
                if g.resource_stop is not None:
                    self.assertIn('resource_stop', res)
                    q = g.resource_stop
                    self.assertEqual(
                        res['resource_stop'],
                        {'scope': q.scope, 'resource': q.resource, 'phase': q.phase,
                         'requested': value_json(q.requested),
                         'effective': value_json(q.effective), 'used': value_json(q.used),
                         'remaining': value_json(q.remaining), 'actions': list(q.actions),
                         'message': q.message})
                else:
                    self.assertNotIn('resource_stop', res)
                self.assertEqual(sorted(res['arms']), sorted(g.arms))
                for side, a in g.arms.items():
                    sa = res['arms'][side]
                    ja = arm_head_json(g, a)
                    for k in ('status', 'complete_to_bp', 'completeness_scope',
                              'frontier_remaining', 'labels_per_node', 'counters',
                              'evidence', 'limitations'):
                        self.assertEqual(sa[k], ja[k], (side, k))
                    self.assertEqual(sa.get('cap_trigger'), ja.get('cap_trigger'),
                                     (side, 'cap_trigger'))
                    self.assertEqual(sa['branch_events_total'], a.branch_events_total)
                    self.assertLessEqual(len(a.branch_events), a.branch_events_total)
                    # the conservative rule (v5.2): a limitation of a class never sits
                    # beside 'complete' in its dimension
                kinds = {l['kind'] for l in res['limitations']} | {
                    l['kind'] for s in res['arms'].values() for l in s['limitations']}
                if kinds & {'walk_domain', 'seed_labels'}:
                    self.assertNotEqual(res['outcome']['walks'], 'complete', kinds)
                if 'branch_events' in kinds:
                    self.assertEqual(res['outcome']['branch_diagnostics'], 'cut', kinds)
                if kinds & {'label_lists', 'inexact_counts', 'seed_labels', 'switch_sources',
                            'greedy_losses'}:
                    self.assertNotEqual(res['outcome']['label_evidence'], 'complete', kinds)
                if cmp_full:
                    full = cell.full_result(i)
                    for k in ('label_mode', 'outcome', 'limitations', 'annotation'):
                        self.assertEqual(res.get(k), full.get(k), k)
                    for k, v in res['seed'].items():
                        if k in full['seed']:
                            self.assertEqual(v, full['seed'][k], ('seed', k))
                    self.assertEqual(len(full['label_dict']), res['seed']['num_labels'])
                    self.assertEqual(len(full['seed']['labels']),
                                     res['seed']['num_seed_labels'])
                    for side, sa in res['arms'].items():
                        fa = full['arms'][side]
                        for k in ('status', 'complete_to_bp', 'completeness_scope',
                                  'frontier_remaining', 'labels_per_node', 'counters',
                                  'evidence', 'limitations'):
                            self.assertEqual(sa[k], fa[k], (side, k))
                        self.assertEqual(sa.get('cap_trigger'), fa.get('cap_trigger'))
                        self.assertEqual(sa['branch_events_total'],
                                         len(fa['branch_events'])
                                         + fa['branch_events_truncated'])

    # ------------------------------------------------------------ identity

    def check_identity(self, row):
        cell, _ = self.cell(row)
        rec = (R.load_manifest().get('servers') or {}).get(row['index'])
        if not rec:
            self.skipTest('the manifest records no server identity for %s' % row['index'])
        caps = cell.graphlet_response['capabilities']
        for k in ('index_ns', 'index_fp', 'index_meta_fp', 'k', 'regime', 'alphabet'):
            self.assertEqual(caps.get(k), rec.get(k), ('capabilities vs manifest', k))
        for i, res, g in self.results(row):
            if g is None:
                continue
            with self.subTest(result=i):
                self.assertEqual(g.index_ns, rec['index_ns'])
                self.assertEqual(g.index_fp, rec['index_fp'])
                self.assertEqual(g.index_meta_fp, rec['index_meta_fp'])
                self.assertEqual(g.k, rec['k'])
                self.assertEqual(g.regime, rec['regime'])
                self.assertEqual(g.alphabet, rec['alphabet'])
                if row['index'] == 'mini_refseq':
                    # started with --index-manifest: identity verifiable (a sha256 hex)
                    self.assertRegex(g.index_fp or '', r'\A[0-9a-f]{64}\Z')
                else:
                    self.assertIsNone(g.index_fp, 'no manifest: index_fp must be "*"')
                h = text_lines(res['graphlet'])[0].split(' ')
                self.assertEqual(h[13], rec['index_ns'] or '*')
                self.assertEqual(h[14], rec['index_fp'] or '*')
                self.assertEqual(h[15], rec['index_meta_fp'] or '*')

    # ------------------------------------------------------------ seed

    def check_seed(self, row):
        cell, _ = self.cell(row)
        req = cell.request
        for i, res, g in self.results(row):
            if g is None:
                continue
            with self.subTest(result=i):
                seq = req['seeds'][i]['sequence'].upper()
                self.assertEqual(g.seed.sequence, seq)
                self.assertEqual(g.seed.length_bp, len(seq))
                self.assertEqual(g.seed.num_kmers, len(seq) - g.k + 1)
                self.assertEqual(g.seed_index, i)
                self.assertEqual(len({l.ref for l in g.labels}), len(g.labels),
                                 'label refs are unique')
                # (a column-kind cell derives the seed's COLUMNS: the catalog lists the
                # header carriers /resolve reported for it)
                if g.mode == 'constrain' and 'column_kind' not in (row.get('tags') or ()):
                    rec = cell.seed_record
                    if rec is not None and rec['kind'] != 'batch' and \
                            not g.seed_summary['seed'].get('labels_dropped'):
                        # derived seed labels = the catalog's full-length carriers
                        self.assertEqual(sorted(l.name for l in g.labels[
                            :g.seed.num_seed_labels]), sorted(rec['labels']))

    # ------------------------------------------------------------ save / load

    def check_save_load(self, row):
        cell, _ = self.cell(row)
        for i, res, g in self.results(row):
            if g is None:
                continue
            with self.subTest(result=i), tempfile.TemporaryDirectory() as d:
                path = os.path.join(d, 'x.mgt')
                n = g.save(path)
                with open(path, encoding='utf-8', newline='') as f:
                    saved = f.read()
                self.assertEqual(n, len(saved.encode('utf-8')))
                self.assertTrue(is_canonical(saved))
                g2 = load(path)
                self.assertEqual(g2.dump(), saved)
                self.assertEqual(g2.dump(envelope=False), res['graphlet'])
                self.assertTrue(g2.has_envelope)
                if row.get('deterministic') and _full_ok(row):
                    self.assertEqual(R.normalize_full_result(g2.to_json()),
                                     R.normalize_full_result(g.to_json()))
                # a body cut at any record boundary never parses into a graphlet
                lines = text_lines(saved)
                for cut in sorted({2, len(lines) // 2, len(lines) - 1}):
                    if 0 < cut < len(lines):
                        with self.assertRaises(GraphletFormatError):
                            parse('\n'.join(lines[:cut]) + '\n')


# ---------------------------------------------------------------- test generation

def _make(check, row):
    def test(self):
        getattr(self, 'check_' + check)(row)
    test.__doc__ = '%s: %s/%s' % (check, row['index'], row['cell'])
    reason = KNOWN_BUGS.get((row['cell'], check))
    if reason:
        test = unittest.expectedFailure(test)
        test.__doc__ += ' [expected failure: %s]' % reason
    return test


def _rows(index):
    try:
        return R.cells(index, heavy=None)
    except Exception:                              # noqa: BLE001 -- no cache: no tests
        return []


def _build_classes():
    out = {}
    for index in R.INDEXES:
        attrs = {}
        for row in _rows(index):
            for check in CHECKS:
                attrs['test_%s__%s' % (row['cell'], check)] = _make(check, row)
        cls = type('TestConformance_%s' % index, (_ConformanceBase,), attrs)
        cls.__doc__ = 'Every cached %s cell: %d cells x %d checks' % (
            index, len(_rows(index)), len(CHECKS))
        cls.__module__ = __name__
        out[cls.__name__] = cls
    return out


globals().update(_build_classes())


class TestConformanceCoverage(unittest.TestCase):
    """The generated suite covers the cache: every ok cell of every index, and the
    situations the design's T37 list names are present in it."""

    def test_the_cache_is_present_and_every_cell_is_covered(self):
        R.require_cache()
        rows = R.cells(heavy=None)
        names = {n for n in dir(sys.modules[__name__]) if n.startswith('TestConformance_')}
        generated = set()
        for n in names:
            cls = getattr(sys.modules[__name__], n)
            generated |= {m for m in dir(cls) if m.startswith('test_')}
        for row in rows:
            for check in CHECKS:
                self.assertIn('test_%s__%s' % (row['cell'], check), generated)

    def test_the_t37_situations_occur_in_the_cache(self):
        R.require_cache()
        seen = collections.Counter()
        for row in R.cells(max_full_bytes=8 << 20):
            meta = row.get('results') or []
            for r in meta:
                if 'error' in r:
                    seen['derivation_error'] += 1
                    continue
                if r.get('label_mode') == 'annotate':
                    seen['annotate'] += 1
                for side, a in (r.get('arms') or {}).items():
                    if a.get('scope') == 'united_history':
                        seen['merge'] += 1
                    if a.get('cap'):
                        seen['cap_' + a['cap']] += 1
                    if 'branch_events' in (a.get('limitations') or []):
                        seen['cut_events'] += 1
                    if 'label_lists' in (a.get('limitations') or []):
                        seen['cut_lists'] += 1
            if 'sequences_false' in (row.get('tags') or ()):
                seen['sequences_false'] += 1
            if 'switch' in (row.get('tags') or ()):
                seen['switch'] += 1
            if row['strategy'] in ('left_only', 'right_only'):
                seen['one_arm'] += 1
        for want in ('derivation_error', 'annotate', 'merge', 'cut_events', 'cut_lists',
                     'sequences_false', 'switch', 'one_arm'):
            self.assertGreater(seen[want], 0, (want, dict(seen)))


if __name__ == '__main__':
    unittest.main()
