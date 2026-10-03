"""metagraph.traverse library operations, end to end on real retrievals.

Every check here runs the library (walks, claims, label_walks, routes, support_profile,
support_changes, splits, to_fasta, to_gfa, subgraph, compare) on graphlets the real
servers returned (the retrieval cache of realdata.py) and compares the answer with an
expectation computed SEPARATELY: from the server's own detail-full JSON of the same
cell, from the independent spellers of the harness, or from the /resolve oracle (the
search path, which shares no code with the walker). Ranking keys, displayed support,
route choices at merges and the k-1 overlaps are re-derived here from the full JSON,
never by calling the library's derive/ops helpers.

Scope and budgets (deterministic, environment-tunable):
  * cells: every deterministic cached cell with a body <= METAGRAPH_REAL_OPS_MAX_BODY
    (1 MB) and a full JSON <= METAGRAPH_REAL_OPS_MAX_FULL (32 MB); the heavy beam cells
    are left out (they belong to the scale tests);
  * per-leaf checks visit at most METAGRAPH_REAL_OPS_LEAF_CAP leaves per arm (evenly
    spaced, so the sample is the same on every run);
  * oracle checks query /resolve for a few walks per arm of the core seeds (cached on
    disk under the scratch oracle_cache, so a rerun is offline-fast).

Skips cleanly without the cache (cache tests) or without the servers (oracle and live
tests). Run:

    cd api/python/tests/real && python3 -m unittest -v test_real_ops
"""

import collections
import json
import math
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import realdata as R  # noqa: E402

from metagraph.traverse import Graphlet, GraphletView, UnknownLabel  # noqa: E402
from metagraph.traverse._codec import REASON, RESOURCE_CODES  # noqa: E402


def _env_int(name, default):
    try:
        return int(os.environ.get(name, '') or default)
    except ValueError:
        return default


MAX_BODY = _env_int('METAGRAPH_REAL_OPS_MAX_BODY', 1 << 20)
MAX_FULL = _env_int('METAGRAPH_REAL_OPS_MAX_FULL', 32 << 20)
LEAF_CAP = _env_int('METAGRAPH_REAL_OPS_LEAF_CAP', 250)
ORACLE_WALKS = _env_int('METAGRAPH_REAL_OPS_ORACLE_WALKS', 4)
LABEL_CAP = _env_int('METAGRAPH_REAL_OPS_LABEL_CAP', 16)
# annotate leaves carry hundreds to thousands of recorded labels each: fewer per arm
ANNOTATE_LEAF_CAP = _env_int('METAGRAPH_REAL_OPS_ANNOTATE_LEAF_CAP', 40)


def leaf_cap(g):
    return LEAF_CAP if g.mode == 'constrain' else ANNOTATE_LEAF_CAP

RESOURCE_NAMES = frozenset(REASON[c] for c in RESOURCE_CODES)
END_REASON_NAMES = frozenset(REASON.values()) | {
    'minority', 'below_min_labels', 'split_limit', 'hairpin', 'superseded',
    'switch_sources'}

# what each test visited, printed at the end of the module (the report's cell counts)
COVERAGE = collections.defaultdict(set)


def tearDownModule():
    if COVERAGE:
        allc = set().union(*COVERAGE.values())
        sys.stderr.write('\n[test_real_ops] %d distinct cells covered: %s\n' % (
            len(allc), ', '.join('%s=%d' % (k, len(v)) for k, v in sorted(COVERAGE.items()))))


# ------------------------------------------------------------------ cell selection

def sample(seq, n):
    """At most |n| items of |seq|, evenly spaced (deterministic)."""
    seq = list(seq)
    if n is None or len(seq) <= n:
        return seq
    step = len(seq) / float(n)
    return [seq[int(i * step)] for i in range(n)]


def op_rows(index=None, **kw):
    kw.setdefault('deterministic', True)
    kw.setdefault('max_body_bytes', MAX_BODY)
    kw.setdefault('max_full_bytes', MAX_FULL)
    return R.cells(index, **kw)


def results(rows, tag, *, full=True, mode=None, sequences=None):
    """(row, cell, i, graphlet, full result or None) per result with a graphlet."""
    for row in rows:
        c = R.load_cell(row['index'], row)
        for i in range(c.n_results()):
            g = c.graphlet(i)
            if g is None:
                continue
            if mode is not None and g.mode != mode:
                continue
            has_bases = all(s.walk is not None for a in g.arms.values() for s in a.segments)
            if sequences is not None and has_bases != sequences:
                continue
            COVERAGE[tag].add(c.index + '/' + c.cell)
            yield row, c, i, g, (c.full_result(i) if full else None)


def require_rows(rows):
    if not rows:
        raise unittest.SkipTest('no cached cells for this selection (run realdata.py fill)')
    return rows


def live_rows(rows):
    """The rows whose index server answers (oracle tests): skips when none does."""
    up = {i for i in R.INDEXES if R.server(i) is not None and R.server_up(i)}
    got = [r for r in rows if r['index'] in up]
    if not got:
        raise unittest.SkipTest('no index server reachable for the oracle')
    return got


# ------------------------------------------------------------------ independent derivations

def full_children(fa):
    out = collections.defaultdict(list)
    for s in fa['segments']:
        for p in s['parents']:
            out[p].append(s['id'])
    return out


def full_leaf_set(fa):
    return {p['segments'][-1] for p in fa['paths']}


def walking_bases(fa, side, sid):
    """A segment's bases in walking order (the full JSON stores natural orientation)."""
    seq = fa['segments'][sid].get('sequence')
    if seq is None:
        return None
    return seq if side == 'right' else seq[::-1]


def full_flank(fa, side, chain):
    """The natural flank of a segment chain root -> x, from the full JSON."""
    parts = [walking_bases(fa, side, s) for s in chain]
    if any(p is None for p in parts):
        return None
    w = ''.join(parts)
    return w if side == 'right' else w[::-1]


def expected_labels_full(g, fa, path):
    """Labels with displayed support on the whole walk, from the full JSON. Constrain:
    leaf labels whose entry was never routed in through a non-first parent (T route_bp
    0) and whose run starts at the seed boundary; annotate: the labels recorded at the
    seed boundary and on every node of the walk."""
    if g.mode == 'constrain':
        runs = fa['runs']
        return sorted(e['label'] for e in path['end_labels']
                      if e['route_bp'] == 0 and runs[e['run']]['from_bp'] == 0)
    return sorted(annotate_alive_end(fa)[path['segments'][-1]])


def annotate_alive_end(fa):
    """Annotate, per segment: the labels recorded at the seed boundary and on every node
    of its first-parent chain (one pass, parents first; memoised on the full result)."""
    got = fa.get('_alive_end')
    if got is None:
        segs = fa['segments']
        got = [None] * len(segs)
        for s in segs:
            alive = set(s['labels']) if not s['parents'] else set(got[s['parents'][0]])
            for ls in s.get('label_sets') or []:
                if not alive:
                    break
                alive &= set(ls['labels'])
            got[s['id']] = frozenset(alive)
        fa['_alive_end'] = got
    return got


def expected_alive(g, fa, path):
    if g.mode == 'constrain':
        return sorted(e['label'] for e in path['end_labels'])
    return sorted(fa['segments'][path['segments'][-1]]['labels_at_end'])


def expected_route_consistent(g, fa, path):
    if g.mode == 'constrain':
        return {e['label'] for e in path['end_labels'] if e['route_bp'] == 0}
    return set(expected_labels_full(g, fa, path))


def expected_min_loss(g, path):
    if g.mode != 'constrain':
        return 0.0
    return min((e['loss'] for e in path['end_labels']), default=math.inf)


def ranking(g, fa, by):
    """The path ids in the order walks(by=...) must return them (§5: support ranks by
    (-len(labels_full), -n_alive, -length_bp, path_id))."""
    rows = fa.get('_rank_rows')
    if rows is None:
        rows = []
        for p in fa['paths']:
            rows.append((p['id'], len(expected_labels_full(g, fa, p)),
                         len(expected_alive(g, fa, p)), p['length_bp'],
                         expected_min_loss(g, p)))
        fa['_rank_rows'] = rows
    keys = {
        'support': lambda r: (-r[1], -r[2], -r[3], r[0]),
        'length': lambda r: (-r[3], r[0]),
        'loss': lambda r: (r[4], -r[2], -r[3], r[0]),
        'id': lambda r: r[0],
    }
    return [r[0] for r in sorted(rows, key=keys[by])]


def full_displayed_pieces(fa, path):
    """Constrain: maximal (from, to, frozenset) of the labels alive at each base of the
    walk's first-parent chain, from the full JSON alone: a segment's entry labels, minus
    every ended run anchored there with to_bp inside the segment (absent from to_bp on),
    plus every switch-in (present from at_bp on); ends before switch-ins at one position
    (spec §7.1: label sets change only through events)."""
    segs = fa['segments']
    runs = fa['runs']
    by_seg = collections.defaultdict(list)
    for r in runs:
        if 'end_reason' in r:
            by_seg[r['segment']].append(r)
    out = []
    for sid in path['segments']:
        s = segs[sid]
        lo, hi = s['from_bp'], s['from_bp'] + s['length_bp']
        if hi == lo:
            continue
        ops_ = [(r['to_bp'], 0, r['label']) for r in by_seg[sid] if r['to_bp'] < hi]
        ops_ += [(e['at_bp'], 1, e['to']) for e in s['events'] if e['type'] == 'switch']
        ops_.sort()
        cur = set(s['labels'])
        pos, k = lo, 0
        while pos < hi:
            while k < len(ops_) and ops_[k][0] <= pos:
                (cur.discard if ops_[k][1] == 0 else cur.add)(ops_[k][2])
                k += 1
            nxt = min(max(ops_[k][0] if k < len(ops_) else hi, pos + 1), hi)
            out.append((pos, nxt, frozenset(cur)))
            pos = nxt
    return merge_pieces(out)


def full_presence_pieces(fa, path):
    """Annotate: the recorded label sets along the chain, maximal, with their totals."""
    out = []
    for sid in path['segments']:
        for ls in fa['segments'][sid].get('label_sets') or []:
            out.append((ls['from_bp'], ls['to_bp'], frozenset(ls['labels']),
                        ls['labels_total']))
    merged = []
    for f, t, s, tot in out:
        if merged and merged[-1][1] == f and merged[-1][2] == s and merged[-1][3] == tot:
            merged[-1] = (merged[-1][0], t, s, tot)
        else:
            merged.append((f, t, s, tot))
    return merged


def merge_pieces(pieces):
    out = []
    for f, t, s in pieces:
        if out and out[-1][1] == f and out[-1][2] == s:
            out[-1] = (out[-1][0], t, s)
        else:
            out.append((f, t, s))
    return out


def lineage_label_at(g_runs, run, depth):
    """§5.1: the lineage's incoming label at |depth| (before a switch committed there)."""
    r = run
    while r.from_bp >= depth and r.from_label is not None and r.prev_run is not None:
        r = g_runs[r.prev_run]
    return r.label


# ------------------------------------------------------------------ oracle geometry

def node_kmer(side, seed_len, flank_len, k, i):
    """The k-mer start, in the natural molecule (left: flank + seed; right: seed +
    flank), of the node at outward base |i| of a walk."""
    if side == 'right':
        return seed_len + i - k + 1
    return flank_len - 1 - i


def node_span(side, seed_len, flank_len, k, f, t):
    """The k-mer start interval [a, b) of the nodes at outward bases [f, t)."""
    if side == 'right':
        return seed_len + f - k + 1, seed_len + t - k + 1
    return flank_len - t, flank_len - f


def covered(runs, a, b):
    """[a, b) inside one oracle run."""
    if b <= a:
        return True
    return any(x <= a and b <= y for x, y in runs)


def touches(runs, a, b):
    return any(x < b and a < y for x, y in runs)


def has_bases(g):
    return all(s.walk is not None for a in g.arms.values() for s in a.segments)


# ================================================================== walks

def a_parents(g, side, w):
    """The merges (segments with several parents) on a walk's first-parent chain."""
    a = g.arms[side]
    return [x for x in w.segments if len(a.segments[x].parents) > 1]


class TestWalksRankingsAndFilters(unittest.TestCase):
    """walks() under every ranking and filter, against the server's detail-full paths."""

    @classmethod
    def setUpClass(cls):
        R.require_cache()
        cls.rows = require_rows(op_rows())

    def test_every_ranking_returns_every_walk_in_the_documented_order(self):
        n = 0
        for row, c, i, g, full in results(self.rows, 'walks.rankings'):
            for side in g.arms:
                fa = full['arms'][side]
                npaths = len(fa['paths'])
                # every walk is materialised (chain, spelling, claims) by walks(): on a
                # large arm the ranking is checked on its first LEAF_CAP entries
                cap = None if npaths <= LEAF_CAP else LEAF_CAP
                with self.subTest(cell=c.cell, result=i, arm=side):
                    for by in ('support', 'length', 'loss', 'id'):
                        got = [w.path_id for w in g.walks(side, by=by, top=cap)]
                        self.assertEqual(ranking(g, fa, by)[:cap], got, by)
                        if cap is None:
                            self.assertEqual(npaths, len(set(got)))
                        top = [w.path_id for w in g.walks(side, by=by, top=3)]
                        self.assertEqual(got[:3], top, by)
                    self.assertEqual(npaths, g.arms[side].counts.leaves)
                    lengths = sorted(p['length_bp'] for p in fa['paths'])
                    if lengths:
                        m = lengths[len(lengths) // 2]
                        got = {w.path_id for w in g.walks(side, by='id', min_bp=m)}
                        self.assertEqual({p['id'] for p in fa['paths']
                                          if p['length_bp'] >= m}, got)
                    n += 1
        self.assertGreater(n, 0)

    def test_walk_fields_match_the_full_json_and_the_independent_speller(self):
        for row, c, i, g, full in results(self.rows, 'walks.fields'):
            bases = has_bases(g)
            for side in g.arms:
                fa = full['arms'][side]
                ctb = fa['complete_to_bp']
                with self.subTest(cell=c.cell, result=i, arm=side):
                    walks = g.walks(side, by='id', top=LEAF_CAP)
                    self.assertEqual(min(LEAF_CAP, len(fa['paths'])), len(walks))
                    for w, p in zip(walks, fa['paths']):
                        self.assertEqual(p['id'], w.path_id)
                        self.assertEqual(p['segments'], list(w.segments))
                        self.assertEqual(p['segments'][-1], w.leaf)
                        self.assertEqual(p['length_bp'], w.length_bp)
                        self.assertEqual(p.get('path_reason'), w.path_reason)
                        self.assertEqual(p['end_reasons'], w.end_reasons)
                        self.assertEqual(len(expected_alive(g, fa, p)), w.n_alive)
                        resource = p.get('path_reason') in RESOURCE_NAMES
                        self.assertEqual((not resource) and p['length_bp'] <= ctb, w.complete)
                        self.assertEqual(max(0, p['length_bp'] - ctb), w.beyond_certified_bp)
                        if bases:
                            mine = R.spell_with_seed(g, side, w.leaf)
                            seed = g.seed.sequence
                            flank = mine[len(seed):] if side == 'right' else \
                                mine[:len(mine) - len(seed)]
                            self.assertEqual(flank, w.sequence, p['id'])
                            self.assertEqual(full_flank(fa, side, p['segments']), w.sequence)
                        else:
                            self.assertIsNone(w.sequence)

    def test_labels_full_and_claims_are_consistent(self):
        for row, c, i, g, full in results(self.rows, 'walks.labels_full'):
            for side in g.arms:
                fa = full['arms'][side]
                with self.subTest(cell=c.cell, result=i, arm=side):
                    for w, p in zip(g.walks(side, by='id', top=LEAF_CAP), fa['paths']):
                        full_ids = [l.id for l in w.labels_full]
                        self.assertEqual(expected_labels_full(g, fa, p), full_ids, p['id'])
                        alive = set(expected_alive(g, fa, p))
                        self.assertLessEqual(set(full_ids), alive)
                        by_label = {}
                        for cl in w.claims:
                            self.assertEqual(w.leaf, cl.segment)
                            self.assertEqual(w.length_bp, cl.to_bp)
                            self.assertEqual(side, cl.arm)
                            self.assertNotIn(cl.label.id, by_label, 'two claims of one label')
                            by_label[cl.label.id] = cl
                        if g.mode == 'constrain':
                            # one claim per run ending at the leaf: exactly the alive labels
                            self.assertEqual(alive, set(by_label))
                            route = {e['label']: e['route_bp'] for e in p['end_labels']}
                            for lid, cl in by_label.items():
                                run = fa['runs'][cl.run]
                                self.assertEqual(lid, run['label'])
                                if cl.kind == 'boundary':
                                    self.assertEqual(0, w.length_bp)
                                    continue
                                self.assertEqual('alive' if p.get('path_reason') else 'end',
                                                 cl.kind)
                                self.assertEqual(max(route[lid], run['from_bp']),
                                                 cl.evidence_from, lid)
                        else:
                            # annotate: the maximal route ends at the leaf's end (the §5.1
                            # union rule, GPT review finding 3): those displayed from 0 are
                            # labels_full; the others entered the displayed chain through a
                            # non-first merge parent and are displayed from route_bp on
                            self.assertLessEqual(set(by_label), alive)
                            shown = {lid for lid, cl in by_label.items() if not cl.route_bp}
                            self.assertEqual(set(full_ids), shown)
                            for cl in by_label.values():
                                self.assertEqual(0, cl.from_bp)
                                if w.length_bp == 0:
                                    self.assertIsNone(cl.evidence_from)
                                else:
                                    self.assertEqual(cl.route_bp, cl.evidence_from)
                                    self.assertLess(cl.evidence_from, w.length_bp + 1)
                                    if cl.route_bp:
                                        self.assertGreater(len(a_parents(g, side, w)), 0)
                        for lid in full_ids:
                            cl = by_label[lid]
                            self.assertTrue(cl.evidence_from == 0 or
                                            (cl.kind == 'boundary' and w.length_bp == 0), lid)

    def test_label_filters_route_consistent_and_union(self):
        for row, c, i, g, full in results(self.rows, 'walks.filters'):
            for side in g.arms:
                fa = full['arms'][side]
                alive_by_path = {p['id']: set(expected_alive(g, fa, p)) for p in fa['paths']}
                rc_by_path = {p['id']: expected_route_consistent(g, fa, p) for p in fa['paths']}
                seen = sorted(set().union(*alive_by_path.values())) if alive_by_path else []
                labels = sample(seen, LABEL_CAP if len(fa['paths']) <= LEAF_CAP else 4)
                never = [l.id for l in g.labels if l.id not in set(seen)][:1]
                with self.subTest(cell=c.cell, result=i, arm=side):
                    for lid in labels + never:
                        lab = {'id': lid}
                        got_rc = {w.path_id for w in g.walks(side, by='id', labels=[lab])}
                        got_any = {w.path_id for w in g.walks(side, by='id', labels=[lab],
                                                              route_consistent=False)}
                        self.assertEqual({p for p, s in rc_by_path.items() if lid in s},
                                         got_rc, lid)
                        self.assertEqual({p for p, s in alive_by_path.items() if lid in s},
                                         got_any, lid)
                        self.assertLessEqual(got_rc, got_any)
                        # a selector by ref and by name selects the same walks
                        ref = g.labels[lid].ref
                        self.assertEqual(got_rc, {w.path_id for w in g.walks(
                            side, by='id', labels=[{'ref': ref}])})
                    if len(labels) >= 2:
                        a, b = labels[0], labels[-1]
                        both = {w.path_id for w in g.walks(side, by='id', labels=[
                            {'id': a}, {'id': b}], route_consistent=False)}
                        self.assertEqual({p for p, s in alive_by_path.items()
                                          if a in s or b in s}, both)


class TestWalksAgainstTheOracle(unittest.TestCase):
    """labels_full is displayed-path support: /resolve must report every such label on
    every node of the walk. Constrain: the labels are seed carriers, so every k-mer of
    the molecule (seed + flank, natural orientation); annotate: the seed-boundary node
    and every flank node (an annotate label need not carry the whole seed)."""

    def test_labels_full_cover_the_whole_molecule(self):
        rows = live_rows(require_rows(op_rows(tier='core')))
        checked = 0
        for row, c, i, g, full in results(rows, 'oracle.labels_full', full=False,
                                          sequences=True):
            if g.support != 'kmer' and g.mode == 'constrain':
                support = 'trace'
            else:
                support = 'kmer'
            for side in g.arms:
                walks = [w for w in g.walks(side, top=ORACLE_WALKS) if w.labels_full]
                for w in walks:
                    mol = R.spell_with_seed(g, side, w.leaf)
                    if len(mol) < g.k:
                        continue
                    n = R.oracle_num_kmers(c.index, mol)
                    runs = R.oracle_label_runs(c.index, mol, w.labels_full, support=support)
                    if g.mode == 'constrain':
                        lo, hi = 0, n
                    else:
                        seed_len = len(g.seed.sequence)
                        lo, hi = node_span(side, seed_len, w.length_bp, g.k, -1, w.length_bp)
                    with self.subTest(cell=c.cell, result=i, arm=side, walk=w.path_id):
                        for lab in w.labels_full:
                            if g.mode == 'constrain':
                                self.assertEqual([(0, n)], runs[lab.name], lab.name)
                            else:
                                self.assertTrue(covered(runs[lab.name], lo, hi),
                                                (lab.name, lo, hi, runs[lab.name]))
                    checked += 1
        self.assertGreater(checked, 0)


# ================================================================== label walks and routes

class TestLabelWalksAndRoutes(unittest.TestCase):
    """label_walks() and routes() for every label of every constrain retrieval: one entry
    per run, each run's OWN route (chosen at every merge through the partition holding the
    lineage's incoming label), its end as the label_end event states it."""

    @classmethod
    def setUpClass(cls):
        R.require_cache()
        cls.rows = require_rows(op_rows(tag='constrain'))

    def test_one_label_walk_per_run_with_the_runs_fields(self):
        n = 0
        for row, c, i, g, full in results(self.rows, 'labels.label_walks', mode='constrain'):
            bases = has_bases(g)
            for side in g.arms:
                fa = full['arms'][side]
                runs = fa['runs']
                children = full_children(fa)
                leaves = full_leaf_set(fa)
                events = {(s['id'], e['label'], e['at_bp']): e['reason'] for s in fa['segments']
                          for e in s['events'] if e['type'] == 'label_end'}
                with self.subTest(cell=c.cell, result=i, arm=side):
                    for lab in g.labels:
                        lws = g.label_walks(lab, side)
                        mine = [k for k, r in enumerate(runs) if r['label'] == lab.id]
                        self.assertEqual(mine, [lw.run for lw in lws], lab.name)
                        for lw in lws:
                            r = runs[lw.run]
                            self.assertEqual((r['from_bp'], r['to_bp']), (lw.from_bp, lw.to_bp))
                            self.assertEqual(lab.ref, lw.label.ref)
                            self.assertEqual(side, lw.arm)
                            if 'end_reason' not in r:
                                want = 'merged'
                            else:
                                want = events.get((r['segment'], lab.id, r['to_bp']), 'switched')
                                if want == 'switched':
                                    self.assertEqual('label_lost', r['end_reason'])
                            self.assertEqual(want, lw.end, lw.run)
                            # the leaves below the anchor
                            below, stack, seen = set(), [r['segment']], set()
                            while stack:
                                s = stack.pop()
                                if s in seen:
                                    continue
                                seen.add(s)
                                if s in leaves:
                                    below.add(s)
                                stack.extend(children[s])
                            self.assertEqual(sorted(below), lw.leaves_below)
                            if want == 'merged':
                                into = [x for x in children[r['segment']]
                                        if len(fa['segments'][x]['parents']) > 1
                                        and fa['segments'][x]['from_bp'] == r['to_bp']]
                                self.assertIn(lw.merged_into, into)
                            else:
                                self.assertIsNone(lw.merged_into)
                            self.assertGreaterEqual(lw.evidence_from, lw.from_bp)
                            self.assertLessEqual(lw.evidence_from, max(lw.from_bp, lw.to_bp))
                            if not bases:
                                self.assertIsNone(lw.sequence)
                        n += 1
        self.assertGreater(n, 0)

    def test_routes_follow_parent_edges_and_honour_merge_partitions(self):
        merges_seen = 0
        for row, c, i, g, full in results(self.rows, 'labels.routes', mode='constrain'):
            bases = has_bases(g)
            for side in g.arms:
                fa = full['arms'][side]
                segs = fa['segments']
                a = g.arms[side]
                with self.subTest(cell=c.cell, result=i, arm=side):
                    for lab in g.labels:
                        rts = g.routes(lab, side)
                        spelled = g.routes(lab, side, spell=True)
                        lws = g.label_walks(lab, side)
                        own = [r for r in a.runs if r.label == lab.id]
                        self.assertEqual(len(own), len(rts))
                        for run, route, (route2, rbases), lw in zip(own, rts, spelled, lws):
                            self.assertEqual(route, route2)
                            self.assertEqual([], segs[route[0]]['parents'])
                            self.assertEqual(run.segment, route[-1])
                            for x, y in zip(route, route[1:]):
                                parents = segs[y]['parents']
                                self.assertIn(x, parents, (run.id, x, y))
                                if len(parents) > 1:
                                    merges_seen += 1
                                    j = parents.index(x)
                                    part = segs[y]['labels_via_parent']
                                    incoming = lineage_label_at(a.runs, run, segs[y]['from_bp'])
                                    self.assertIn(incoming, part[j], (run.id, y))
                                    # the first partition holding the incoming label
                                    first = next(k for k, p in enumerate(part) if incoming in p)
                                    self.assertEqual(first, j, (run.id, y))
                            # spelled: the route's own bases up to to_bp, natural orientation
                            if bases:
                                flank = full_flank(fa, side, route)
                                walk = flank if side == 'right' else flank[::-1]
                                cut = walk[:run.to_bp]
                                self.assertEqual(cut if side == 'right' else cut[::-1], rbases)
                                part = cut[run.from_bp:run.to_bp]
                                self.assertEqual(part if side == 'right' else part[::-1],
                                                 lw.sequence)
                            else:
                                self.assertIsNone(rbases)
        if merges_seen == 0:
            self.skipTest('no route through a merge in the cached constrain cells')

    def test_a_run_with_route_bp_0_supports_its_displayed_path(self):
        """Spec §7.1 counts runs with entered_by seed, from_bp 0 and route_bp 0 as sharing
        the displayed flank: no such run's displayed chain may pass a merge its lineage
        entered through a non-first parent. Fixed SRV-ROUTE-BP-CLONE (walker.cpp,
        commit_entries): a run cloned at a split AFTER such a merge did not inherit the
        merge's route stamp (the library derives evidence from the partitions and was
        right; the R / runs[].route_bp field was not). Old repro:
        sra_16s_PZ326290__limit1, left arm, runs[7] (label 0, [0, 53), route_bp 0,
        segment 5): chain 0-1-3-5, merge 3 labels_via_parent [[3], [0,1,2,4,5]] -> label
        0 came via parent 2 (see test_the_split_clone_keeps_the_merge_route_stamp). Needs
        the cache filled from the fixed server.
        """
        bad = []
        for row, c, i, g, full in results(self.rows, 'labels.route_bp', mode='constrain'):
            for side in g.arms:
                fa = full['arms'][side]
                a = g.arms[side]
                for k, r in enumerate(fa['runs']):
                    if r['entered_by'] != 'seed' or r['from_bp'] != 0 or r['route_bp'] != 0:
                        continue
                    chain_merge = self._displayed_deviation(fa, a, a.runs[k])
                    if chain_merge is not None:
                        bad.append((c.cell, i, side, k, chain_merge))
        self.assertEqual([], bad[:20], '%d runs with route_bp 0 whose own route leaves the '
                         'displayed chain' % len(bad))

    @staticmethod
    def _displayed_deviation(fa, a, run):
        """The from_bp of a merge on the run's displayed chain (root -> anchor via first
        parents) at or before to_bp whose first partition does not hold the lineage's
        incoming label; None when the displayed chain is the lineage's own route."""
        segs = fa['segments']
        s = run.segment
        while True:
            seg = segs[s]
            if not seg['parents']:
                return None
            if len(seg['parents']) > 1 and seg['from_bp'] <= run.to_bp:
                incoming = lineage_label_at(a.runs, run, seg['from_bp'])
                if incoming not in seg['labels_via_parent'][0]:
                    return seg['from_bp']
            s = seg['parents'][0]


class TestRoutesAgainstTheOracle(unittest.TestCase):
    """Each run's own route is carried by its label: /resolve reports the label on every
    node of [from_bp, to_bp) along the route (k-mer support), whatever the displayed
    chain is. Annotate stretches: the label is recorded on every node of [0, j) of its
    walk and (exact recording) absent at node j when the stretch ends inside the walk."""

    def test_every_sampled_run_is_carried_along_its_own_route(self):
        rows = live_rows(require_rows(op_rows(tier='core', tag='constrain')))
        checked = deviating = 0
        for row, c, i, g, full in results(rows, 'oracle.routes', full=False,
                                          mode='constrain', sequences=True):
            if g.support != 'kmer':
                continue
            for side in g.arms:
                a = g.arms[side]
                per_label = collections.defaultdict(list)
                for lab in g.labels:
                    for run, (route, bases) in zip([r for r in a.runs if r.label == lab.id],
                                                   g.routes(lab, side, spell=True)):
                        if run.to_bp > run.from_bp:
                            per_label[lab.id].append((run, route, bases))
                picked = []
                for lid in sample(sorted(per_label), 12):
                    rs = per_label[lid]
                    dev = [x for x in rs if x[1] != _first_chain(a, x[0].segment)]
                    picked.extend(dev[:2] or rs[:1])
                for run, route, bases in picked:
                    seed = g.seed.sequence
                    mol = seed + bases if side == 'right' else bases + seed
                    lo, hi = node_span(side, len(seed), len(bases), g.k, run.from_bp, run.to_bp)
                    if run.from_bp == 0 and run.from_label is None:
                        # entered by the seed: the label carries the seed too
                        if side == 'right':
                            lo = 0
                        else:
                            hi = R.oracle_num_kmers(c.index, mol)
                    runs = R.oracle_label_runs(c.index, mol, [g.labels[run.label]])
                    with self.subTest(cell=c.cell, result=i, arm=side, run=run.id):
                        self.assertTrue(covered(runs[g.labels[run.label].name], lo, hi),
                                        (lo, hi, runs[g.labels[run.label].name]))
                    checked += 1
                    deviating += route != _first_chain(a, run.segment)
        self.assertGreater(checked, 0)
        sys.stderr.write('\n[routes oracle] %d runs, %d on a route off the displayed chain\n'
                         % (checked, deviating))

    def test_the_split_clone_keeps_the_merge_route_stamp(self):
        """The minimal reproduction of the route_bp clone bug (fixed SRV-ROUTE-BP-CLONE,
        see TestLabelWalksAndRoutes): the run is a clone made at a split after the merge
        (segment 3, from 50) its lineage entered through a non-first parent. Its R and
        runs[].route_bp now carry the merge's stamp (50; the old server wrote 0), equal to
        the library's partition-derived evidence; /resolve: its own route is covered,
        the displayed chain is not."""
        R.require_cache('sra', 'sra_16s_PZ326290__limit1')
        c = R.load_cell('sra', 'sra_16s_PZ326290__limit1')
        g, full = c.graphlet(0), c.full_result(0)
        a = g.arms['left']
        run = a.runs[7]
        self.assertEqual((0, 0, 53), (run.label, run.from_bp, run.to_bp))
        claim = [x for x in g.claims('left', labels=[{'id': 0}]) if x.run == 7][0]
        self.assertEqual(50, claim.evidence_from)
        self.assertEqual(50, run.route_bp)
        self.assertEqual(50, full['arms']['left']['runs'][7]['route_bp'])
        # the run the clone was made from: the same lineage on the other child of the
        # split, stamped by the merge itself
        twins = [r for r in a.runs if r.id != 7 and (r.label, r.from_bp, r.prev_run)
                 == (run.label, run.from_bp, run.prev_run)]
        self.assertTrue(twins)
        self.assertIn(50, [r.route_bp for r in twins])
        lab = g.labels[0]
        displayed = R.spell_with_seed(g, 'left', run.segment)
        n = R.oracle_num_kmers('sra', displayed)
        self.assertNotEqual([(0, n)], R.oracle_label_runs('sra', displayed, [lab])[lab.name])
        route, bases = g.routes(lab, 'left', spell=True)[
            [r.id for r in a.runs if r.label == 0].index(7)]
        own = bases + g.seed.sequence
        m = R.oracle_num_kmers('sra', own)
        self.assertEqual([(0, m)], R.oracle_label_runs('sra', own, [lab])[lab.name])

    def test_annotate_stretches_are_recorded_and_maximal(self):
        rows = live_rows(require_rows(op_rows(tier='core', tag='annotate')))
        checked = 0
        for row, c, i, g, full in results(rows, 'oracle.stretches', full=False,
                                          mode='annotate', sequences=True):
            for side in g.arms:
                a = g.arms[side]
                exact = a.labels_per_node.nodes_truncated == 0
                labels = sample([l for l in g.labels if l.id in set(a.segments[0].entry)], 6)
                for lab in labels:
                    lws = g.label_walks(lab, side)
                    rts = g.routes(lab, side)        # one route per label walk, same order
                    self.assertEqual(len(lws), len(rts))
                    for lw, route in sample(list(zip(lws, rts)), 3):
                        if lw.leaves_below and lw.evidence_from == 0:
                            # a displayed walk spells the label walk's whole route
                            mol = R.spell_with_seed(g, side, lw.leaves_below[0])
                        else:
                            # a route through a non-first merge parent (the union rule
                            # of direct_bp): no walk displays it, so it is spelled here
                            # from its own G records, one base beyond its end when the
                            # end lies inside the anchor
                            anchor = a.segments[route[-1]]
                            hi = lw.to_bp + (1 if lw.to_bp < anchor.end_bp else 0)
                            flank = ''.join(a.segments[x].walk for x in route)[:hi]
                            mol = (g.seed.sequence + flank if side == 'right'
                                   else flank[::-1] + g.seed.sequence)
                        flank_len = len(mol) - len(g.seed.sequence)
                        runs = R.oracle_label_runs(c.index, mol, [lab])[lab.name]
                        with self.subTest(cell=c.cell, arm=side, label=lab.ref, to=lw.to_bp):
                            lo, hi = node_span(side, len(g.seed.sequence), flank_len, g.k,
                                               0, lw.to_bp)
                            self.assertTrue(covered(runs, lo, hi), (lo, hi, runs))
                            if exact and lw.to_bp < flank_len:
                                j = node_kmer(side, len(g.seed.sequence), flank_len, g.k,
                                              lw.to_bp)
                                self.assertFalse(covered(runs, j, j + 1), (j, runs))
                        checked += 1
        self.assertGreater(checked, 0)


def _first_chain(a, seg):
    out = [seg]
    while a.segments[out[-1]].parents:
        out.append(a.segments[out[-1]].parents[0])
    return out[::-1]


# ================================================================== support profiles

class TestSupportProfilesAndChanges(unittest.TestCase):
    """support_profile and support_changes along every leaf (capped per arm): the
    displayed profile re-derived from the full JSON, the route profile containing it,
    and every change explained -- a label leaves by an end (with its reason), a switch or
    a split, and joins by a switch or a merge."""

    @classmethod
    def setUpClass(cls):
        R.require_cache()
        cls.rows = require_rows(op_rows())

    def test_displayed_profile_matches_the_full_json(self):
        n = 0
        for row, c, i, g, full in results(self.rows, 'support.profile'):
            for side in g.arms:
                fa = full['arms'][side]
                segs = fa['segments']
                with self.subTest(cell=c.cell, result=i, arm=side):
                    for p in sample(fa['paths'], leaf_cap(g)):
                        prof = g.support_profile(side, p['id'])
                        L = p['length_bp']
                        if L == 0:
                            self.assertEqual([], prof)
                            continue
                        self.assertEqual(0, prof[0].from_bp)
                        self.assertEqual(L, prof[-1].to_bp)
                        for x, y in zip(prof, prof[1:]):
                            self.assertEqual(x.to_bp, y.from_bp)
                        if g.mode == 'constrain':
                            mine = [(r.from_bp, r.to_bp, frozenset(l.id for l in r.labels))
                                    for r in prof]
                            for x, y in zip(mine, mine[1:]):
                                self.assertNotEqual(x[2], y[2], 'not maximal')
                            self.assertEqual(full_displayed_pieces(fa, p), mine, p['id'])
                            self.assertTrue(all(r.exact for r in prof))
                            # at each segment's last base: the segment's labels_at_end
                            for sid in p['segments']:
                                s = segs[sid]
                                end = s['from_bp'] + s['length_bp']
                                if s['length_bp'] == 0:
                                    continue
                                at = [r for r in mine if r[0] <= end - 1 < r[1]][0]
                                self.assertEqual(set(s['labels_at_end']), set(at[2]), sid)
                            route = g.support_profile(side, p['id'], kind='route')
                            for r in route:
                                for d in prof:
                                    if d.from_bp <= r.from_bp < d.to_bp:
                                        self.assertLessEqual({l.id for l in d.labels},
                                                             {l.id for l in r.labels})
                        else:
                            mine = [(r.from_bp, r.to_bp, frozenset(l.id for l in r.labels),
                                     r.total) for r in prof]
                            self.assertEqual(full_presence_pieces(fa, p), mine, p['id'])
                            for r in prof:
                                self.assertEqual(r.total == len(r.labels), r.exact)
                        n += 1
        self.assertGreater(n, 0)

    def test_support_changes_explain_every_label_leaving_and_joining(self):
        why_seen = collections.Counter()
        for row, c, i, g, full in results(self.rows, 'support.changes'):
            for side in g.arms:
                fa = full['arms'][side]
                ctx = _ChangeContext(g, fa, side)
                with self.subTest(cell=c.cell, result=i, arm=side):
                    for p in sample(fa['paths'], leaf_cap(g)):
                        chain = p['segments']
                        prof = g.support_profile(side, p['id'])
                        changes = g.support_changes(side, p['id'])
                        # where: exactly the boundaries where the set changes (and the
                        # first base, against the root's entry)
                        prev = set(ctx.segs[chain[0]]['labels'])
                        want = []
                        for r in prof:
                            cur = {l.id for l in r.labels}
                            if cur != prev:
                                want.append((r.from_bp, sorted(cur - prev), sorted(prev - cur)))
                            prev = cur
                        self.assertEqual(want, [(x.at_bp, [l.id for l in x.added],
                                                 [l.id for l in x.removed]) for x in changes])
                        walk = ctx.walk(chain)
                        for ch in changes:
                            removed = {l.id for l in ch.removed}
                            added = {l.id for l in ch.added}
                            got = {(r['change'], r['label']['ref']) for r in ch.reasons}
                            self.assertEqual({('removed', g.labels[x].ref) for x in removed}
                                             | {('added', g.labels[x].ref) for x in added}, got)
                            for r in ch.reasons:
                                lid = ctx.by_ref[r['label']['ref']]
                                why_seen[(g.mode, r['change'], r['why'])] += 1
                                problem = ctx.unexplained(walk, ch.at_bp, lid, r)
                                self.assertIsNone(problem, (p['id'], ch.at_bp, r))
        sys.stderr.write('\n[support_changes] reasons: %s\n' % dict(why_seen))
        self.assertTrue(why_seen)


    def test_annotate_split_reasons_name_the_branch_that_took_the_label(self):
        """A label recorded on the parent's last node and on NO child's first node is on
        no successor: 'absent', never {'why': 'split', 'took': []} -- a split that no
        branch took (fixed product bug SPLIT-REASON-EMPTY, ops._change_reasons; it was
        13 % of the annotate split reasons). Repro of the old answer:
        sra_16s_PZ326290__annotate_exh, left arm, support_changes('left', 16): at 16 label
        c:142 removed as split with took [] although none of segment 1's children 3, 4, 5
        records it (nodes_truncated 0)."""
        bad = []
        rows = op_rows(tag='annotate', tier='core', strategy=('annotate_exh',))
        for row, c, i, g, full in results(rows, 'support.annotate_split', full=False):
            for side, a in g.arms.items():
                if a.labels_per_node.nodes_truncated:
                    continue
                for w in sample(g.walks(side, by='id', top=LEAF_CAP), 10):
                    for ch in g.support_changes(side, w.path_id):
                        for r in ch.reasons:
                            if r['why'] == 'split' and not r['took']:
                                bad.append((c.cell, side, w.path_id, ch.at_bp,
                                            r['label']['ref']))
        self.assertEqual([], bad[:10], '%d split reasons with nothing taken' % len(bad))


class _ChangeContext:
    """Per arm, from the full JSON: what explains a label leaving or joining a walk at a
    position (written independently of ops._change_reasons). unexplained() returns None
    when the stated reason holds, else what is wrong."""

    def __init__(self, g, fa, side):
        self.g = g
        self.side = side
        self.segs = fa['segments']
        self.by_ref = {l.ref: l.id for l in g.labels}
        self.children = full_children(fa)
        self.events = {(s['id'], e['label'], e['at_bp']): e['reason'] for s in self.segs
                       for e in s['events'] if e['type'] == 'label_end'}
        self.switches = collections.defaultdict(list)       # at_bp -> [(segment, event)]
        for s in self.segs:
            for e in s['events']:
                if e['type'] == 'switch':
                    self.switches[e['at_bp']].append((s['id'], e))
        self.runs_end = collections.defaultdict(list)       # (label, to_bp) -> [run]
        for r in fa['runs']:
            self.runs_end[(r['label'], r['to_bp'])].append(r)

    def walk(self, chain):
        single, merge = {}, {}
        for sid in chain:
            s = self.segs[sid]
            (merge if len(s['parents']) > 1 else single if s['parents'] else {}) \
                .setdefault(s['from_bp'], sid)
        return set(chain), single, merge

    def first_base(self, sid):
        seq = self.segs[sid].get('sequence')
        if not seq:
            return None
        return seq[0] if self.side == 'right' else seq[-1]

    def unexplained(self, walk, at, lid, r):
        on_chain, single, merge = walk
        segs = self.segs
        why = r['why']
        g = self.g
        if g.mode != 'constrain':
            if r['change'] == 'added':
                return None if why == 'present' else 'annotate addition %r' % why
            if why == 'absent':
                return None
            if why != 'split':
                return 'annotate removal %r' % why
            s = single.get(at)
            if s is None:
                return 'a split reason off a split'
            if lid not in segs[segs[s]['parents'][0]]['labels_at_end']:
                return 'split of a label not at the parent end'
            return None
        if r['change'] == 'removed':
            if why == 'ended':
                return 'a label left without a reason'
            if why == 'switched':
                ok = any(x['segment'] in on_chain and x.get('end_reason') == 'label_lost'
                         and (x['segment'], lid, at) not in self.events
                         for x in self.runs_end[(lid, at)])
                if not ok:
                    return 'no silent end of the label here'
                tos = sorted({e['to'] for s, e in self.switches[at] if s in on_chain
                              and e['from'] == lid})
                if [g.labels[x].ref for x in tos] != [x['ref'] for x in r['to']]:
                    return 'switch targets %r' % tos
                return None
            if why == 'split':
                s = single.get(at)
                if s is None:
                    return 'a split reason off a split'
                parent = segs[s]['parents'][0]
                if lid not in segs[parent]['labels_at_end'] or lid in segs[s]['labels']:
                    return 'not a label that went another way'
                sibs = [x for x in self.children[parent] if x != s]
                if all(self.first_base(x) is not None or segs[x]['length_bp'] == 0
                       for x in sibs):
                    took = sorted({self.first_base(x) for x in sibs
                                   if lid in segs[x]['labels']} - {None})
                    if took != r['took']:
                        return 'took %r' % took
                    if not took:
                        return 'a split that no other branch took'
                return None
            if why not in END_REASON_NAMES:
                return 'unknown reason %r' % why
            if not any(self.events.get((x, lid, at)) == why for x in on_chain):
                return 'no label_end %s here' % why
            return None
        if why == 'entered':
            return 'a label joined without a reason'
        if why == 'switch':
            ok = any(s in on_chain and e['to'] == lid
                     and g.labels[e['from']].ref == r['from']['ref']
                     for s, e in self.switches[at])
            return None if ok else 'no switch into the label here'
        if why != 'merge':
            return 'unknown joining reason %r' % why
        m = merge.get(at)
        if m is None:
            return 'a merge reason off a merge'
        parents = segs[m]['parents']
        if r['via_parent'] not in parents:
            return 'via a non-parent'
        j = parents.index(r['via_parent'])
        if j == 0 or lid not in segs[m]['labels_via_parent'][j]:
            return 'not routed in through a later parent'
        return None


def refused(a, at, ref, g):
    """Whether a stored branch event at |at| refused a successor to the label."""
    lid = g.label({'ref': ref}).id
    return any(be.at_bp == at and any(lid in r.labels for r in be.refused)
               for be in a.branch_events)


class TestSupportAgainstTheOracle(unittest.TestCase):
    """The displayed profile on real k-mers: constrain -- every label displayed at base i
    is on node i, and a label that left by label_lost is not on the node it was lost at;
    annotate with exact recording -- the recorded set at node i IS the set of dictionary
    labels /resolve reports on node i."""

    def test_displayed_support_is_present_and_label_lost_is_absent(self):
        rows = live_rows(require_rows(op_rows(tier='core', tag='constrain')))
        checked = lost = 0
        for row, c, i, g, full in results(rows, 'oracle.support', full=False,
                                          mode='constrain', sequences=True):
            for side in g.arms:
                a = g.arms[side]
                for w in sample(g.walks(side, by='id'), ORACLE_WALKS):
                    prof = g.support_profile(side, w.path_id)
                    names = sorted({l.name for r in prof for l in r.labels}
                                   | {l.name for ch in g.support_changes(side, w.path_id)
                                      for l in ch.removed})
                    if not names or w.length_bp == 0:
                        continue
                    mol = R.spell_with_seed(g, side, w.leaf)
                    seed_len = len(g.seed.sequence)
                    runs = R.oracle_label_runs(c.index, mol, names)
                    with self.subTest(cell=c.cell, result=i, arm=side, walk=w.path_id):
                        for r in prof:
                            lo, hi = node_span(side, seed_len, w.length_bp, g.k,
                                               r.from_bp, r.to_bp)
                            for lab in r.labels:
                                self.assertTrue(covered(runs[lab.name], lo, hi),
                                                (lab.name, r.from_bp, r.to_bp))
                        for ch in g.support_changes(side, w.path_id):
                            for why in ch.reasons:
                                if why['change'] != 'removed' or ch.at_bp >= w.length_bp:
                                    continue
                                # label_lost: the label is not on the next node. A split
                                # of an unlimited (exhaustive) trie: the label went the
                                # other way only, so it is not on this branch's node --
                                # unless a refusal (quorum, split limit) kept it off
                                if why['why'] == 'split':
                                    if not g.envelope['strategy'].get('exhaustive') or \
                                            refused(a, ch.at_bp, why['label']['ref'], g):
                                        continue
                                elif why['why'] != 'label_lost':
                                    continue
                                lab = g.label({'ref': why['label']['ref']})
                                j = node_kmer(side, seed_len, w.length_bp, g.k, ch.at_bp)
                                self.assertFalse(covered(runs[lab.name], j, j + 1),
                                                 (lab.name, ch.at_bp, why['why']))
                                lost += 1
                    checked += 1
        self.assertGreater(checked, 0)
        sys.stderr.write('\n[support oracle] %d walks, %d label_lost ends and exhaustive '
                         'split departures checked absent\n' % (checked, lost))

    def test_annotate_recorded_sets_are_the_node_label_sets(self):
        rows = live_rows(require_rows(op_rows(tier='core', tag='annotate')))
        checked = 0
        for row, c, i, g, full in results(rows, 'oracle.recorded', full=False,
                                          mode='annotate', sequences=True):
            if len(g.labels) > 3000:
                continue
            for side in g.arms:
                for w in sample(g.walks(side, by='id'), ORACLE_WALKS):
                    prof = g.support_profile(side, w.path_id)
                    if not prof or not all(r.exact for r in prof):
                        continue
                    mol = R.spell_with_seed(g, side, w.leaf)
                    runs = R.oracle_label_runs(c.index, mol, g.labels)
                    seed_len = len(g.seed.sequence)
                    with self.subTest(cell=c.cell, result=i, arm=side, walk=w.path_id):
                        for r in prof:
                            lo, hi = node_span(side, seed_len, w.length_bp, g.k,
                                               r.from_bp, r.to_bp)
                            want = {l.name for l in r.labels}
                            cov = {n for n, rr in runs.items() if covered(rr, lo, hi)}
                            touch = {n for n, rr in runs.items() if touches(rr, lo, hi)}
                            self.assertEqual(want, cov, (r.from_bp, r.to_bp))
                            self.assertEqual(want, touch, (r.from_bp, r.to_bp))
                    checked += 1
        self.assertGreater(checked, 0)


# ================================================================== splits

class TestSplits(unittest.TestCase):
    """derive.splits (what graphlet_splits pages) vs the stored G.split and the server's
    splits[]: every split parent and only split parents carry G.split, kind ambiguous
    iff G.split is 1, the A count holds."""

    @classmethod
    def setUpClass(cls):
        R.require_cache()
        cls.rows = require_rows(op_rows())

    def test_splits_match_G_split_and_the_full_json(self):
        from metagraph.traverse import derive
        n = 0
        for row, c, i, g, full in results(self.rows, 'splits'):
            for side, a in g.arms.items():
                fa = full['arms'][side]
                with self.subTest(cell=c.cell, result=i, arm=side):
                    got = derive.splits(a, g.mode)
                    self.assertEqual(a.counts.splits, len(got))
                    self.assertEqual(len(fa['splits']), len(got))
                    single = collections.Counter(s.parents[0] for s in a.segments
                                                 if len(s.parents) == 1)
                    parents = {sp.segment for sp in got}
                    for s in a.segments:
                        is_split = single.get(s.id, 0) >= 2
                        self.assertEqual(is_split, s.split is not None, s.id)
                        self.assertEqual(is_split, s.id in parents, s.id)
                    for sp, fs in zip(got, fa['splits']):
                        self.assertEqual(fs['segment'], sp.segment)
                        self.assertEqual(fs['at_bp'], sp.at_bp)
                        self.assertEqual(fs['children'], list(sp.children))
                        self.assertEqual(fs['labels_before'], sp.labels_before)
                        self.assertEqual('ambiguous' if a.segments[sp.segment].split
                                         else 'divergence', fs['kind'])
                        self.assertEqual(fs['kind'] == 'ambiguous', sp.ambiguous)
                        self.assertGreaterEqual(len(sp.children), 2)
                        seg = a.segments[sp.segment]
                        self.assertEqual(seg.end_bp, sp.at_bp)
                        for ch in sp.children:
                            self.assertEqual((sp.segment,), a.segments[ch].parents)
                            self.assertEqual(sp.at_bp, a.segments[ch].from_bp)
                        if g.mode == 'constrain':
                            # labels follow branches: a child's entry is a subset of the
                            # parent's end set
                            for ch in sp.children:
                                self.assertLessEqual(set(a.segments[ch].entry), set(seg.end))
                        fb = [(b['char'], b['labels_distinct']) for b in fs['branches']]
                        self.assertEqual(fb, [(b['char'], b['labels_distinct']) for b in
                                              derive.split_branches(a, sp, g.cap)])
                    n += 1
        self.assertGreater(n, 0)


# ================================================================== FASTA and GFA

def read_fasta_text(text):
    out = []
    head, parts = None, []
    for line in text.split('\n'):
        if line.startswith('>'):
            if head is not None:
                out.append((head, ''.join(parts)))
            head, parts = line[1:], []
        elif line:
            parts.append(line)
    if head is not None:
        out.append((head, ''.join(parts)))
    return out


def header_fields(head):
    name, *rest = head.split(' ')
    f = dict(x.split('=', 1) for x in rest)
    side, pid = name.rsplit('_', 1)
    return side, int(pid), f


class TestFasta(unittest.TestCase):
    """to_fasta re-read: one record per walk, the sequence equal to spell() and to the
    independent speller, headers carrying the arm, walk id, length, end reason(s), the
    number of labels alive at the leaf and whether the seed is included."""

    @classmethod
    def setUpClass(cls):
        R.require_cache()
        cls.rows = require_rows(op_rows())

    def test_records_reread_equal_the_walks(self):
        n = 0
        for row, c, i, g, full in results(self.rows, 'export.fasta', sequences=True):
            with self.subTest(cell=c.cell, result=i):
                recs = read_fasta_text(g.to_fasta())
                self.assertEqual(sum(len(fa['paths']) for fa in full['arms'].values()),
                                 len(recs))
                seen = set()
                for head, seq in recs:
                    side, pid, f = header_fields(head)
                    fa = full['arms'][side]
                    p = fa['paths'][pid]
                    seen.add((side, pid))
                    self.assertEqual(R.spell_with_seed(g, side, p['segments'][-1]), seq)
                    self.assertEqual(g.spell(side, pid, with_seed=True), seq)
                    self.assertEqual(int(f['length_bp']), p['length_bp'])
                    self.assertEqual(int(f['labels']), len(expected_alive(g, fa, p)))
                    self.assertEqual('yes', f['seed'])
                    reason = p.get('path_reason') or ','.join(sorted(p['end_reasons'])) or '.'
                    self.assertEqual(reason, f['reason'])
                self.assertEqual(len(recs), len(seen), 'duplicate record names')
                # per arm, a subset of walks, wrapped, without the seed, walking order
                for side, fa in full['arms'].items():
                    ids = sample([p['id'] for p in fa['paths']], 5)
                    wrapped = read_fasta_text(g.to_fasta(side, ids, width=60))
                    self.assertEqual([g.spell(side, x, with_seed=True) for x in ids],
                                     [s for _, s in wrapped])
                    text = g.to_fasta(side, ids, width=60)
                    self.assertTrue(all(len(line) <= 60 for line in text.split('\n')
                                        if line and not line.startswith('>')))
                    bare = read_fasta_text(g.to_fasta(side, ids, with_seed=False))
                    self.assertEqual([full_flank(fa, side, fa['paths'][x]['segments'])
                                      for x in ids], [s for _, s in bare])
                    self.assertTrue(all(header_fields(h)[2]['seed'] == 'no' for h, _ in bare))
                    walk = read_fasta_text(g.to_fasta(side, ids, orientation='walk'))
                    self.assertEqual([g.spell(side, x, 'walk', True) for x in ids],
                                     [s for _, s in walk])
                n += 1
        self.assertGreater(n, 0)

    def test_records_are_in_the_graph(self):
        rows = live_rows(require_rows(op_rows(tier='core')))
        checked = 0
        for row, c, i, g, full in results(rows, 'oracle.fasta', full=False, sequences=True):
            for side in g.arms:
                ids = sample([w.path_id for w in g.walks(side, by='id')], 2)
                for head, seq in read_fasta_text(g.to_fasta(side, ids)):
                    if len(seq) < g.k:
                        continue
                    with self.subTest(cell=c.cell, record=head):
                        self.assertTrue(R.oracle_present(c.index, seq), head)
                    checked += 1
        self.assertGreater(checked, 0)


def parse_gfa(text):
    S, L, P = {}, [], []
    for line in text.splitlines():
        f = line.split('\t')
        if f[0] == 'S':
            S[f[1]] = f[2]
        elif f[0] == 'L':
            L.append(f)
        elif f[0] == 'P':
            P.append(f)
    return S, L, P


class TestGfa(unittest.TestCase):
    """to_gfa re-parsed: S per segment (+ the seed), L per parent edge (merge parents
    included) with a k-1 overlap that holds on the real bases, P per walk spelling the
    walk's molecule on both arms, LB = the leaf's labels, ER = its end reason."""

    @classmethod
    def setUpClass(cls):
        R.require_cache()
        cls.rows = require_rows(op_rows())

    def test_gfa_reparses_to_the_walks(self):
        n = merges = 0
        for row, c, i, g, full in results(self.rows, 'export.gfa', sequences=True):
            k1 = g.k - 1
            with self.subTest(cell=c.cell, result=i):
                S, L, P = parse_gfa(g.to_gfa())
                self.assertEqual(1 + sum(len(a.segments) for a in g.arms.values()), len(S))
                self.assertEqual(sum(max(1, len(s.parents)) for a in g.arms.values()
                                     for s in a.segments), len(L))
                self.assertEqual(sum(len(fa['paths']) for fa in full['arms'].values()), len(P))
                for l in L:
                    self.assertEqual('%dM' % k1, l[5])
                    self.assertEqual(S[l[1]][-k1:], S[l[3]][:k1], l[:4])
                for side, a in g.arms.items():
                    for s in a.segments:
                        if len(s.parents) > 1:
                            merges += 1
                for p in P:
                    names = [x[:-1] for x in p[2].split(',')]
                    self.assertTrue(all(x.endswith('+') for x in p[2].split(',')))
                    spelled = S[names[0]] + ''.join(S[x][k1:] for x in names[1:])
                    side, pid = p[1].rsplit('_', 1)
                    fa = full['arms'][side]
                    path = fa['paths'][int(pid)]
                    self.assertEqual(R.spell_with_seed(g, side, path['segments'][-1]), spelled)
                    tags = dict(t.split(':Z:', 1) for t in p[4:])
                    refs = [g.labels[x].ref for x in expected_alive(g, fa, path)]
                    self.assertEqual(','.join(refs) or '.', tags['LB'])
                    self.assertEqual(path.get('path_reason') or 'semantic', tags['ER'])
                # without the seed: no seed record, P spells k-1 bases of context + flank
                S2, L2, P2 = parse_gfa(g.to_gfa(with_seed=False))
                self.assertNotIn('seed', S2)
                self.assertEqual(len(P), len(P2))
                for p in P2:
                    names = [x[:-1] for x in p[2].split(',')]
                    spelled = S2[names[0]] + ''.join(S2[x][k1:] for x in names[1:])
                    side, pid = p[1].rsplit('_', 1)
                    fa = full['arms'][side]
                    flank = full_flank(fa, side, fa['paths'][int(pid)]['segments'])
                    seed = g.seed.sequence
                    want = seed[-k1:] + flank if side == 'right' else flank + seed[:k1]
                    self.assertEqual(want, spelled)
                n += 1
        self.assertGreater(n, 0)
        if merges == 0:
            sys.stderr.write('\n[gfa] no merge in the selection: merge L lines untested\n')

    def test_gfa_refuses_a_retrieval_without_bases(self):
        rows = op_rows(tag='sequences_false')
        if not rows:
            self.skipTest('no sequences:false cell cached')
        for row, c, i, g, full in results(rows[:3], 'export.gfa_nobases', full=False):
            with self.assertRaises(ValueError):
                g.to_gfa()


# ================================================================== subgraph views

class TestSubgraphViews(unittest.TestCase):
    """subgraph(): ids preserved, the selected labels' walks all present (and only
    those), claims restricted to them, completeness qualified, save/load keeps the
    backing body byte-exact."""

    @classmethod
    def setUpClass(cls):
        R.require_cache()
        cls.rows = require_rows(op_rows())

    def _picks(self, g, fa_by_side):
        seen = set()
        for fa in fa_by_side.values():
            for p in fa['paths']:
                seen.update(expected_alive(g, fa, p))
        return sample(sorted(seen), 4)

    def test_views_keep_ids_and_every_walk_of_the_selected_labels(self):
        n = 0
        for row, c, i, g, full in results(self.rows, 'subgraph.views'):
            fas = full['arms']
            picks = self._picks(g, fas)
            with self.subTest(cell=c.cell, result=i):
                for lid in picks:
                    view = g.subgraph([{'id': lid}])
                    self.assertIsInstance(view, GraphletView)
                    self.assertIs(g, view.backing)
                    for side, a in g.arms.items():
                        fa = fas[side]
                        kept = view.segments[side]
                        self.assertTrue(set(kept) <= set(range(len(a.segments))))
                        for sid in kept:     # closed under the first-parent chain
                            if a.segments[sid].parents:
                                self.assertIn(a.segments[sid].parents[0], set(kept))
                        want = sorted(p['id'] for p in fa['paths']
                                      if lid in expected_alive(g, fa, p))
                        self.assertEqual(want, view.path_ids[side], (side, lid))
                        got = view.walks(side, by='id')
                        self.assertEqual(want, [w.path_id for w in got])
                        backing = {w.path_id: w for w in g.walks(side, by='id')}
                        for w in got:
                            b = backing[w.path_id]
                            self.assertEqual((b.leaf, b.segments, b.sequence, b.length_bp),
                                             (w.leaf, w.segments, w.sequence, w.length_bp))
                            self.assertIn(w.leaf, set(kept))
                    self.assertEqual([x.as_dict() for x in g.claims(labels=[{'id': lid}],
                                                                   strict=False)],
                                     [x.as_dict() for x in view.claims(strict=False)])
                    comp = view.completeness()
                    for side, a in g.arms.items():
                        self.assertEqual(a.complete_to_bp, comp[side]['complete_to_bp'])
                        self.assertEqual('for the selected labels', comp[side]['qualified'])
                        self.assertEqual(a.scope, comp[side]['scope'])
                    summ = view.summary()
                    self.assertEqual([g.labels[lid].as_dict()], summ['labels'])
                    self.assertEqual([{'ref': g.labels[lid].ref}], summ['view']['selectors'])
                    if has_bases(g):
                        recs = read_fasta_text(view.to_fasta())
                        self.assertEqual(sorted((s, x) for s in g.arms
                                                for x in view.path_ids[s]),
                                         sorted((header_fields(h)[0], header_fields(h)[1])
                                                for h, _ in recs))
                if len(picks) >= 2:
                    a_, b_ = picks[0], picks[-1]
                    allv = g.subgraph([{'id': a_}, {'id': b_}], mode='all')
                    anyv = g.subgraph([{'id': a_}, {'id': b_}], mode='any')
                    for side, fa in fas.items():
                        both = sorted(p['id'] for p in fa['paths']
                                      if {a_, b_} <= set(expected_alive(g, fa, p)))
                        either = sorted(p['id'] for p in fa['paths']
                                        if {a_, b_} & set(expected_alive(g, fa, p)))
                        self.assertEqual(both, allv.path_ids[side])
                        self.assertEqual(either, anyv.path_ids[side])
                n += 1
        self.assertGreater(n, 0)

    def test_save_load_keeps_the_backing_body_and_the_view(self):
        n = 0
        rows = sample(self.rows, 40)
        with tempfile.TemporaryDirectory() as tmp:
            for row, c, i, g, full in results(rows, 'subgraph.save', full=False):
                body = c.text(i)
                labs = [l for l in g.labels if any(l.id in s.end for a in g.arms.values()
                                                   for s in a.segments)]
                if not labs:
                    continue
                side = sorted(g.arms)[0]
                view = g.subgraph([labs[0]], arm=side)
                path = os.path.join(tmp, '%s_%d.mgt' % (c.cell, i))
                view.save(path)
                with open(path, encoding='utf-8') as f:
                    text = f.read()
                lines = text.split('\n')
                blines = body.split('\n')
                with self.subTest(cell=c.cell, result=i):
                    self.assertTrue(lines[1].startswith('J '))
                    # the body unchanged around the J line; Z counts the document's lines,
                    # J included
                    self.assertEqual(blines[:-2], [lines[0]] + lines[2:-2])
                    self.assertEqual(['Z %d' % (int(blines[-2][2:]) + 1), ''], lines[-2:])
                    j = json.loads(lines[1][2:])
                    self.assertEqual(view.view_spec(), j['view'])
                    back = Graphlet.load(path)
                    self.assertIsInstance(back, GraphletView)
                    self.assertEqual(view.segments, back.segments)
                    self.assertEqual(view.path_ids, back.path_ids)
                    self.assertEqual(view.view_spec(), back.view_spec())
                    self.assertEqual(body, back.backing.dump(envelope=False))
                n += 1
        self.assertGreater(n, 0)


# ================================================================== compare

def _cell_graphlet(index, cell, i=0):
    meta = R.cell_meta(index, cell)
    if meta is None or meta.get('status') != 'ok':
        return None
    return R.load_cell(index, cell).graphlet(i)


def _rank(v):
    return {True: 0, 'qualified': 1, 'unverifiable': 2, 'unknown': 3, False: 4}[v]


def expected_verdict(a, b, sides=None):
    """The comparability compare() must report, from first principles: identity
    (manifest digest on both -> True, else unverifiable), qualified when a side's label
    evidence is not complete or the completeness scopes differ, unknown without a common
    certified depth; the weakest wins."""
    v = True if a.index_fp and b.index_fp else 'unverifiable'
    sides = sides or [s for s in ('left', 'right') if s in a.arms and s in b.arms]
    if a.outcome.label_evidence != 'complete' or b.outcome.label_evidence != 'complete' or \
            any(a.arms[s].labels_per_node.nodes_truncated or b.arms[s].labels_per_node
                .nodes_truncated for s in sides) or \
            any(a.arms[s].scope != b.arms[s].scope for s in sides):
        v = max(v, 'qualified', key=_rank)
    if min(min(a.arms[s].complete_to_bp, b.arms[s].complete_to_bp) for s in sides) == 0:
        v = max(v, 'unknown', key=_rank)
    return v


class TestCompare(unittest.TestCase):
    """compare(): constrain vs annotate on the same seed and family, tuned vs exhaustive
    prefix_subset, the index identity of each deployment (SRA and UHGG have no manifest:
    unverifiable; mini_refseq has one: comparable), lower bounds and differing scopes
    qualified, a self-comparison clean."""

    @classmethod
    def setUpClass(cls):
        R.require_cache()

    def _pairs(self, a_strategy, b_strategy, index=None):
        out = []
        for idx in ([index] if index else R.INDEXES):
            for seed in R.catalog_seeds(idx):
                if seed['kind'] == 'batch':
                    continue
                a = _cell_graphlet(idx, R.cell_id(seed['name'], a_strategy))
                b = _cell_graphlet(idx, R.cell_id(seed['name'], b_strategy))
                if a is not None and b is not None:
                    COVERAGE['compare'].update({idx + '/' + R.cell_id(seed['name'], a_strategy),
                                                idx + '/' + R.cell_id(seed['name'], b_strategy)})
                    out.append((idx, seed['name'], a, b))
        if not out:
            self.skipTest('no cached %s / %s pairs' % (a_strategy, b_strategy))
        return out

    def test_constrain_vs_annotate_claims(self):
        """The constrained exhaustive trie and the label-free annotate trie of one seed
        make the same claims for the seed's labels within the common certified depth;
        equal only when the identity is verifiable (mini_refseq) and both sides' label
        evidence is complete, else no difference and no equality claimed."""
        equal_on_mini = 0
        for idx, seed, a, b in self._pairs('exhaustive', 'annotate_exh'):
            sel = [{'ref': l.ref} for l in a.labels[:a.seed.num_seed_labels]]
            with self.subTest(index=idx, seed=seed):
                cmp = a.compare(b, labels=sel, mode='claims')
                want = expected_verdict(a, b)
                self.assertEqual(want, cmp.comparable, cmp.notes)
                self.assertEqual(min(min(a.arms[s].complete_to_bp, b.arms[s].complete_to_bp)
                                     for s in a.arms), cmp.depth_used)
                self.assertEqual(([], [], []), (cmp.only_in_a, cmp.only_in_b, cmp.differ))
                self.assertEqual(True if want is True else None, cmp.equal)
                if want is not True:
                    self.assertTrue(any('equality is not claimed' in x for x in cmp.notes))
                if idx == 'mini_refseq' and cmp.equal:
                    equal_on_mini += 1
                self.assertEqual(('kmer', 'kmer'), tuple(cmp.support))
                # label-consistent route support agrees between the modes (reach_bp does
                # not: annotate reach is the last recorded presence, constrain reach the
                # lineage's end)
                lab = a.compare(b, labels=sel, mode='labels')
                for d in lab.differ:
                    self.assertEqual(d['a']['direct_bp'], d['b']['direct_bp'], d)
                self.assertEqual(([], []), (lab.only_in_a, lab.only_in_b))
        self.assertGreater(equal_on_mini, 0, 'no equal pair on mini_refseq')

    def test_identity_unverifiable_without_a_manifest_and_comparable_with_one(self):
        for idx, seed, a, b in self._pairs('exhaustive', 'annotate_exh'):
            with self.subTest(index=idx, seed=seed):
                ident = R.recorded_capabilities(idx) or {}
                if idx == 'mini_refseq':
                    self.assertIsNotNone(a.index_fp)
                    self.assertEqual(ident.get('index_fp'), a.index_fp)
                else:
                    self.assertIsNone(a.index_fp)
                    self.assertIsNone(b.index_fp)
                    cmp = a.compare(b, mode='claims')
                    self.assertIn(cmp.comparable, ('unverifiable', 'unknown'))
                    self.assertIsNone(cmp.equal)
                self.assertEqual(a.index_meta_fp, b.index_meta_fp)
                self.assertEqual(a.index_ns, b.index_ns)

    def test_different_indexes_and_different_seeds_are_not_comparable(self):
        pairs = {idx: (a, b) for idx, seed, a, b in self._pairs('strict', 'exhaustive')}
        if len(pairs) < 2:
            self.skipTest('needs cells of two indexes')
        idxs = sorted(pairs)
        a = pairs[idxs[0]][0]
        b = pairs[idxs[1]][0]
        cmp = a.compare(b)
        self.assertIs(False, cmp.comparable)
        self.assertIn('different index', cmp.reason)
        self.assertIsNone(cmp.equal)
        for idx in idxs:
            seeds = [x for x in self._pairs('strict', 'exhaustive', idx)]
            if len(seeds) >= 2:
                cmp = seeds[0][2].compare(seeds[1][2])
                self.assertIs(False, cmp.comparable)
                self.assertEqual('different oriented seed sequences', cmp.reason)

    def test_lower_bound_and_differing_scopes_are_qualified(self):
        checked = collections.Counter()
        for other in ('cut_lists', 'switch1', 'strict', 'limit1', 'limit2', 'annotate_merge'):
            for idx, seed, a, b in self._pairs(other, 'exhaustive'):
                with self.subTest(strategy=other, index=idx, seed=seed):
                    # the seed labels both sides recorded (a ref only one side knows makes
                    # compare() raise: see test_a_ref_known_to_one_side_only_is_compared)
                    known = {l.ref for l in a.labels}
                    sel = [{'ref': l.ref} for l in b.labels[:b.seed.num_seed_labels]
                           if l.ref in known] or None
                    for x, y in ((a, b), (b, a)):
                        cmp = x.compare(y, labels=sel, mode='claims')
                        want = expected_verdict(x, y)
                        self.assertEqual(want, cmp.comparable, cmp.notes)
                        if want is not True:
                            self.assertIsNone(cmp.equal)
                        lower = [t for t, g in (('a', x), ('b', y))
                                 if g.outcome.label_evidence != 'complete']
                        for t in lower:
                            self.assertTrue(any(n.startswith(t + ': ') for n in cmp.notes),
                                            cmp.notes)
                        if any(x.arms[s].scope != y.arms[s].scope for s in x.arms):
                            self.assertTrue(any('completeness scopes differ' in n
                                                for n in cmp.notes), cmp.notes)
                            checked['scope'] += 1
                        if lower:
                            checked['lower'] += 1
                        if idx == 'mini_refseq' and want == 'qualified':
                            checked['qualified_on_mini'] += 1
        self.assertGreater(checked['scope'], 0)
        self.assertGreater(checked['lower'], 0)
        self.assertGreater(checked['qualified_on_mini'], 0)

    def test_prefix_subset_tuned_vs_exhaustive(self):
        """Every walk a tuned retrieval reports under a label is a prefix of an exhaustive
        walk on which the label supports that prefix; the walks the tuned one omits come
        with a reason, and none is unexplained when both keep per-path histories."""
        held = 0
        for tuned in ('limit2_keep', 'strict', 'limit1', 'limit2', 'switch1'):
            for idx, seed, a, b in self._pairs(tuned, 'exhaustive'):
                with self.subTest(tuned=tuned, index=idx, seed=seed):
                    cmp = a.compare(b, mode='prefix_subset')
                    self.assertEqual([], cmp.only_in_a)
                    want = expected_verdict(a, b)
                    self.assertEqual(want, cmp.comparable, cmp.notes)
                    self.assertEqual(True if want is True else None, cmp.equal)
                    same_scope = all(a.arms[s].scope == b.arms[s].scope for s in a.arms)
                    bwalks = {s: [w.sequence for w in b.walks(s, by='id')] for s in b.arms}
                    for row in cmp.only_in_b:
                        self.assertIn('reason', row)
                        self.assertEqual(b.label({'ref': row['ref']}).name, row['name'])
                        if same_scope:
                            self.assertNotEqual('unexplained', row['reason'], row)
                        side = row['arm']
                        if side == 'right':
                            self.assertTrue(any(w.startswith(row['walk']) for w in bwalks[side]))
                        else:
                            self.assertTrue(any(w.endswith(row['walk']) for w in bwalks[side]))
                    held += 1
        self.assertGreater(held, 0)

    def test_self_comparison_finds_no_difference(self):
        rows = sample(op_rows(), 60)
        n = 0
        for row, c, i, g, full in results(rows, 'compare.self', full=False):
            for mode in ('claims', 'walks', 'labels', 'prefix_subset'):
                with self.subTest(cell=c.cell, result=i, mode=mode):
                    cmp = g.compare(g, mode=mode)
                    self.assertEqual(([], []), (cmp.only_in_a, cmp.only_in_b))
                    self.assertEqual([], cmp.differ)
                    want = expected_verdict(g, g)
                    if not has_bases(g) and mode != 'labels':
                        # claims and walks are keyed by their displayed bases: without
                        # them they cannot be matched (unknown, the reason naming
                        # output.sequences; GPT review, finding 5) -- stated, never guessed
                        want = max(want, 'unknown', key=_rank)
                        self.assertIn('output.sequences', cmp.reason)
                    self.assertEqual(want, cmp.comparable, cmp.notes)
                    self.assertEqual(True if want is True else None, cmp.equal)
            n += 1
        self.assertGreater(n, 0)


    def test_a_ref_known_to_one_side_only_is_compared(self):
        """A label selector resolves against |a|, then |b| (fixed product bug
        COMPARE-SELECTORS-A-ONLY): a {'ref': ...} that |b| recorded and |a| did not (a
        cut list, a label absent from one retrieval) is compared -- b's claims of it are
        only_in_b -- instead of raising UnknownLabel; comparisons are keyed by LabelRef
        (§5), and a lower-bound side missing a label is exactly what a comparison must
        show. Was: a = mini_ndm1__cut_lists, b = mini_ndm1__exhaustive, a.compare(b,
        labels=[{'ref': r} for r in b's 19 seed-label refs]) -> UnknownLabel 'no label
        with ref h:1:1' (14 of the 19 are not in a's dictionary). Both directions must
        now give the same rows with the sides swapped."""
        n = 0
        for idx, seed, a, b in self._pairs('cut_lists', 'exhaustive'):
            refs = [{'ref': l.ref} for l in b.labels[:b.seed.num_seed_labels]]
            known_a = {x.ref for x in a.labels}
            if all(l['ref'] in known_a for l in refs):
                continue
            with self.subTest(index=idx, seed=seed):
                ab = a.compare(b, labels=refs)
                ba = b.compare(a, labels=refs)
                self.assertEqual((ab.comparable, ab.depth_used), (ba.comparable, ba.depth_used))
                key = lambda r: json.dumps(r, sort_keys=True)  # noqa: E731
                self.assertEqual(sorted(map(key, ab.only_in_a)), sorted(map(key, ba.only_in_b)))
                self.assertEqual(sorted(map(key, ab.only_in_b)), sorted(map(key, ba.only_in_a)))
                wanted = {l['ref'] for l in refs}
                for row in ab.only_in_a + ab.only_in_b:
                    self.assertIn(row['label']['ref'], wanted)
                # b's claims of a label a never recorded are b's alone
                only_b = {r['label']['ref'] for r in ab.only_in_b}
                self.assertTrue(only_b & (wanted - known_a) or not any(
                    c.label.ref in wanted - known_a
                    for c in b.claims(at_most_bp=ab.depth_used, strict=False)))
                # a ref neither side knows is still an error
                with self.assertRaises(UnknownLabel):
                    a.compare(b, labels=[{'ref': 'h:99999:99999'}])
            n += 1
        if n == 0:
            self.skipTest('no cut_lists / exhaustive pair with a one-sided label')

    def test_a_retrieval_without_bases_compares_to_the_same_one_with_bases(self):
        """Claims and walks are keyed by their displayed bases; a retrieval made with
        output.sequences false has none (fixed product bug COMPARE-NO-BASES; GPT review,
        finding 5). Against the SAME traversal made with bases, claims, walks and
        prefix_subset are 'unknown' (equal None, no rows, the reason naming
        output.sequences), never spurious one-sided rows with comparable True and equal
        False; labels mode needs no bases. Was: mini_ndm1__limit1 vs
        mini_ndm1__no_sequences, mode claims -> (True, False, 53 only_in_a, 48
        only_in_b); prefix_subset raised ValueError."""
        n = 0
        for idx, seed, a, b in self._pairs('limit1', 'no_sequences'):
            self.assertEqual(a.compare(b, mode='labels').only_in_a, [])
            for x, y in ((a, b), (b, a)):
                with self.subTest(index=idx, seed=seed, a=has_bases(x)):
                    for mode in ('claims', 'walks', 'prefix_subset'):
                        cmp = x.compare(y, mode=mode)
                        self.assertEqual(('unknown', None), (cmp.comparable, cmp.equal))
                        self.assertEqual(([], [], []),
                                         (cmp.only_in_a, cmp.only_in_b, cmp.differ))
                        self.assertIn('output.sequences', cmp.reason)
            n += 1
        self.assertGreater(n, 0)


class TestCompareAtDepthZero(unittest.TestCase):
    """A retrieval certified to 0 bp (bounds.max_steps 1 stops the right arm before its
    first base) compares 'unknown' with anything, on every deployment: equality over an
    empty window is not equality."""

    def test_unknown_at_depth_zero(self):
        n = 0
        for idx in R.INDEXES:
            if not R.server_up(idx):
                continue
            seeds = [s for s in R.catalog_seeds(idx, tier='core') if s['kind'] != 'batch']
            if not seeds:
                continue
            seed = seeds[0]
            other = _cell_graphlet(idx, R.cell_id(seed['name'], 'exhaustive'))
            if other is None:
                continue
            req, resp = R.live_traverse(idx, seed['name'], strategy_dict={
                'exhaustive': True, 'bounds': {'max_steps': 1}})
            from metagraph.traverse import from_response
            g = from_response(resp['results'][0], resp)
            if g.arms['right'].complete_to_bp != 0:
                sys.stderr.write('\n[depth 0] max_steps 1 certified %d bp on %s: skipped\n'
                                 % (g.arms['right'].complete_to_bp, idx))
                continue
            for mode in ('claims', 'walks', 'labels', 'prefix_subset'):
                for x, y in ((g, other), (other, g), (g, g)):
                    with self.subTest(index=idx, mode=mode):
                        cmp = x.compare(y, arm='right', mode=mode)
                        self.assertEqual(('unknown', None, 0),
                                         (cmp.comparable, cmp.equal, cmp.depth_used))
                        self.assertIn('no common certified depth (complete_to_bp 0 on a or b)',
                                      cmp.notes)
            n += 1
        if n == 0:
            self.skipTest('no server up')


if __name__ == '__main__':
    unittest.main()
