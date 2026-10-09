"""The third external review of stage 2 (GPT, HEAD 278a53dd), Python side: one class per
finding or design decision.

  4   continuing could exceed the original loss budget: next_request() subtracted C's
      loss_used, the SMALLEST terminal loss of the continued labels;
  5   a switched continuation built an invalid request: its labels.extra still held the
      label the continuation now seeds with (backend 400);
  6   merges beyond the comparison depth caused false differences before it: compare()
      cut the deeper DAG's displayed walks and claims instead of restricting the DAG;
  N2  an argument-binding error ignored the max_bytes the caller supplied;
  N3  the live test oracle accepted a server whose manifest digest differed;
  D1  annotate witness routes, D4 one-sided manifests (unverifiable, opt-in), D5 the
      walk filter's key, D6 receipts of operations that make something, D7 the local
      operations' (missing) budgets.

And the adversarial review of those fixes (round 3):

  A   compare() still differed when the comparison depth equalled a merge position;
  B   a malformed continuation override, an arm that is no string and resolve labels
      that are no list raised instead of answering a bounded error;
  C   traverse_continue made the per-label loss list essential, so a walk with many
      labels was refused at the default ceiling;
  D   a label alive at the leaf but not among the continuation's labels was dropped from
      the continuation unstated.

data/traverse/review3/ holds three responses (detail graphlet, timing off) of the review's
reproducers, produced with build_debug/metagraph on the CLI fixture indexes of
scripts/traversal/graphlet_fixtures.py with an index manifest (so their identity is
verifiable): bubbles_annotate_r30 / _r100 (the annotate fixture's request, right arm,
max_labels_per_node 64, radius 30 / 100) and budgetchain_r80 (the switch_chain request on
the chain index extended by the records E, F, G and H, seed labels A and E, extra
B C D F G H). Round 3 adds, from the same binary and indexes: merge_constrain_r36 /
_r100 / _b0_r100 (the 'merge' CLI fixture walked to 36 and 100, and to 100 with no branch
allowance), bubbles_annotate_r36 / _r62 (the request of bubbles_annotate_r100 walked to the
merges), budgetchain_ae_r31 (the budgetchain request walked to 31: B switched in at 30) and
wide_forbid_r40 (30 records with long names on one walk, forbid, radius 40). The classes
marked CLI also run build_debug/metagraph and skip without it.
"""

import collections
import copy
import gzip
import hashlib
import json
import os
import random
import shutil
import subprocess
import tempfile
import unittest

import traverse_testlib as T

from metagraph.traverse import (  # noqa: E402
    GraphletStore, IncompatibleContinuations, NextRequest, derive, from_response, ops,
    parse,
)
from metagraph.traverse.client import TraverseClient  # noqa: E402
from metagraph.traverse.mcp_tools import MIN_MAX_BYTES, GraphletTools, _size  # noqa: E402
from metagraph.traverse.model import Leaf, PresenceRun  # noqa: E402

DATA3 = os.path.join(T.DATA, 'review3')
REAL = os.path.join(T.DATA, 'real')
MODES = ('claims', 'walks', 'labels', 'prefix_subset')


def response3(name):
    with open(os.path.join(DATA3, name + '.graphlet.json'), encoding='utf-8') as f:
        return json.load(f)


def graphlet3(name):
    o = response3(name)
    return from_response(o['results'][0], o)


def real_graphlets():
    """The offline real fixtures (data/traverse/real): (cell, graphlet) per retrieval."""
    out = []
    for index in sorted(os.listdir(REAL)):
        d = os.path.join(REAL, index)
        if not os.path.isdir(d):
            continue
        for name in sorted(os.listdir(d)):
            if name.endswith('.graphlet.json.gz'):
                with gzip.open(os.path.join(d, name), 'rt', encoding='utf-8') as f:
                    o = json.load(f)
                for i, r in enumerate(o.get('results') or []):
                    if 'graphlet' in r:
                        out.append(('%s/%s#%d' % (index, name[:-len('.graphlet.json.gz')], i),
                                    from_response(r, o)))
    return out


def fixture_graphlets():
    """Every whole-document fixture, the comparison bodies and the review's responses."""
    out = [(name, T.graphlet(name)) for name in sorted(T.DOCUMENTS)]
    for name in sorted(os.listdir(T.COMPARE)):
        if name.endswith('.mgt'):
            out.append((name, parse(T.compare_text(name[:-4]))))
    for name in ('bubbles_annotate_r30', 'bubbles_annotate_r100', 'budgetchain_r80'):
        out.append((name, graphlet3(name)))
    return out


# ------------------------------------------------------------------ the CLI

def _cli():
    path = os.path.join(T.REPO, 'build_debug', 'metagraph')
    return path if os.path.exists(path) else None


class _Indexes:
    """The CLI fixture indexes, built once per test run in a temporary directory, each
    with an index manifest (so that its retrievals' identity is verifiable)."""
    root = None
    built = {}

    @classmethod
    def G(cls):
        if not hasattr(cls, '_G'):
            cls._G = T.fixture_script()
            # the review's budgetchain: the chain index with E (S T1 T2 and the first 20 bp
            # of T3, so it ends inside the walk) and F, G, H chaining 30 bp blocks on
            G = cls._G
            G.INDEXES['budgetchain'] = copy.deepcopy(G.INDEXES['chain'])
            spec = G.INDEXES['budgetchain']
            spec['blocks'].update(T5=30, T6=30, T7=30)
            spec['files']['chain.fa'].extend([
                ('E', ['S', 'T1', 'T2', ('T3', 0, 20)]), ('F', [('T3', 10, 30), 'T5']),
                ('G', [('T5', 10, 30), 'T6']), ('H', [('T6', 10, 30), 'T7'])])
        return cls._G

    @classmethod
    def get(cls, name):
        if name not in cls.built:
            if cls.root is None:
                cls.root = tempfile.mkdtemp(prefix='review3-')
            G = cls.G()
            graph, anno = G.build_index(name, G.resolve_cli(_cli()), cls.root)
            d = os.path.dirname(graph)
            files = []
            loaded = {os.path.basename(graph), os.path.basename(anno)}
            for p in sorted(os.listdir(d)):
                full = os.path.join(d, p)
                if p == 'index.manifest.json' or not os.path.isfile(full):
                    continue
                # one graph with one annotation: the column annotation the column_coord one
                # was transformed from is another annotation, which a manifest must not list
                # (the server refuses it since the review of pass 5)
                if p.endswith('dbg') and p not in loaded:
                    continue
                with open(full, 'rb') as f:
                    digest = hashlib.sha256(f.read()).hexdigest()
                files.append({'path': p, 'size': os.path.getsize(full), 'sha256': digest})
            text = ''.join('%s\t%d\t%s\n' % (x['path'], x['size'], x['sha256']) for x in files)
            manifest = os.path.join(d, 'index.manifest.json')
            with open(manifest, 'w') as f:
                json.dump({'files': files,
                           'index_fp': hashlib.sha256(text.encode()).hexdigest()}, f)
            cls.built[name] = (graph, anno, manifest)
        return cls.built[name]

    @classmethod
    def run(cls, name, request):
        """-> (returncode, parsed stdout or its text)."""
        graph, anno, manifest = cls.get(name)
        with tempfile.NamedTemporaryFile('w', suffix='.json', delete=False) as f:
            json.dump(dict(request), f)
        try:
            p = subprocess.run([cls.G().resolve_cli(_cli()), 'traverse', '--index-name', name,
                                '--index-manifest', manifest, '-i', graph, '-a', anno, f.name],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        finally:
            os.unlink(f.name)
        out = p.stdout.decode('utf-8')
        try:
            return p.returncode, json.loads(out)
        except ValueError:
            return p.returncode, out

    @classmethod
    def graphlet(cls, name, request):
        rc, out = cls.run(name, request)
        if rc:
            raise AssertionError('traverse failed (%d): %s' % (rc, str(out)[-500:]))
        return from_response(out['results'][0], out)


def tearDownModule():
    if _Indexes.root is not None:
        shutil.rmtree(_Indexes.root, ignore_errors=True)
        # another module reusing _Indexes (test_traverse_stage2_recheck) builds them again
        _Indexes.root = None
        _Indexes.built = {}


def cli_requests():
    """(name, index, request with detail graphlet) of every CLI fixture."""
    G = _Indexes.G()
    for name, fx in sorted(G.CLI_FIXTURES.items()):
        r = G.materialize(fx['request'], G.blocks_of(G.INDEXES[fx['index']]))
        r['strategy'].setdefault('output', {})['detail'] = 'graphlet'
        yield name, fx['index'], r


# ------------------------------------------------------------------ finding 4

class TestFinding4ContinuationLossBudget(unittest.TestCase):
    """A continuation must never let a route exceed its original loss budget. The C
    record's loss_used is the SMALLEST terminal loss of the labels the walk continues
    with: the review's budgetchain walk 1 ends with E at loss 0 and C at loss 2 under a
    budget of 3, so subtracting loss_used left the whole budget 3 to C as well, and the
    continuation reached 180 bp accepting a cumulative cost of 4 where the uninterrupted
    walk stops at 150 (loss_budget). The request carries one loss budget and no per-label
    starting loss, so it is now reduced by the LARGEST terminal loss -- exact for C,
    conservative for E -- and next_request(), deepen(), the continuation and
    traverse_continue say so."""

    def setUp(self):
        self.g = graphlet3('budgetchain_r80')

    def test_the_continuation_states_each_labels_loss(self):
        c = self.g.continuation('right', 1)
        self.assertEqual(['E', 'C'], [l.name for l in c.labels])
        self.assertEqual((0.0, 2.0), c.losses)
        self.assertEqual(0.0, c.loss_used)                  # C's: the smallest
        self.assertIn('conservative', c.note)
        self.assertIn("'E' 3", c.note)                      # E's own remaining budget
        one = self.g.continuation('right', 0)                # D alone: nothing to state
        self.assertEqual(((1.0,), None), (one.losses, one.note))

    def test_the_budget_is_reduced_by_the_largest_terminal_loss(self):
        req = self.g.next_request('right', [1], bp=100)
        self.assertIsInstance(req, NextRequest)
        self.assertEqual(1.0, req['strategy']['labels']['loss_budget'])   # 3 - 2, was 3
        self.assertEqual({'original': 3.0, 'effective': 1.0, 'largest_terminal_loss': 2.0},
                         {k: req.loss_budget[k] for k in ('original', 'effective',
                                                          'largest_terminal_loss')})
        self.assertEqual({'E': (0.0, 3.0), 'C': (2.0, 1.0)},
                         {x['name']: (x['loss'], x['remaining'])
                          for x in req.loss_budget['labels']})
        note, = [n for n in req.notes if 'loss_budget' in n and 'conservative' in n]
        self.assertIn("'E' (loss 0, could spend 3)", note)
        # nothing of it reaches the server: json writes the request's fields only
        self.assertEqual({'seeds', 'strategy'}, set(json.loads(json.dumps(req))))

    def test_keeping_the_whole_budget_is_stated(self):
        req = self.g.next_request('right', [1], bp=100, reduce_budget=False)
        self.assertEqual(3.0, req['strategy']['labels']['loss_budget'])
        self.assertTrue(any('may exceed its original budget by up to 2' in n
                            for n in req.notes), req.notes)

    def test_several_walks_take_the_largest_loss_of_all(self):
        # walk 0 continues under D, walk 1 under E and C: their pools differ, so they
        # are not one request; per walk, each takes its own labels' largest loss
        with self.assertRaises(IncompatibleContinuations):
            self.g.next_request('right', [0, 1])
        r0, r1 = self.g.next_requests('right', [0, 1], bp=100)
        self.assertEqual((2.0, 1.0), (r0['strategy']['labels']['loss_budget'],
                                      r1['strategy']['labels']['loss_budget']))

    def test_the_callers_overrides_are_final_and_stated(self):
        # a budget the caller raises is theirs, and said to allow exceeding the original
        req = self.g.next_request('right', [1], labels={'loss_budget': 3})
        self.assertEqual(3, req['strategy']['labels']['loss_budget'])
        self.assertTrue(any("caller's" in n and 'may exceed' in n for n in req.notes))
        self.assertEqual([], T.request_violations(req))
        # labels.extra is rebuilt against the cost model the request will carry
        req = self.g.next_request('right', [1], labels={
            'change_cost': {'model': 'table', 'default': 'forbid',
                            'entries': [['C', 'F', 1], ['E', 'G', 0.5]]}})
        self.assertEqual(['F', 'G'], req['strategy']['labels']['extra'])
        self.assertEqual([], T.request_violations(req))
        # a list the caller gives is sent as given
        req = self.g.next_request('right', [1], labels={'extra': ['F']})
        self.assertEqual(['F'], req['strategy']['labels']['extra'])

    def test_a_next_request_is_a_plain_request_to_its_users(self):
        req = self.g.next_request('right', [1])
        self.assertEqual(dict(req), req)
        again = copy.deepcopy(req)
        self.assertEqual((req.notes, req.loss_budget, req.left_out),
                         (again.notes, again.loss_budget, again.left_out))
        self.assertEqual(json.loads(json.dumps(req)), dict(req))

    def test_switch_chain_is_unchanged(self):
        # one label per walk: the largest loss is loss_used (3 - 2)
        req = T.graphlet('switch_chain').next_request('right', [1])
        self.assertEqual(1.0, req['strategy']['labels']['loss_budget'])
        self.assertFalse([n for n in req.notes if 'conservative' in n])

    def test_deepen_and_the_tool_carry_the_statement(self):
        class Session:
            def __init__(self, response):
                self.response, self.sent = response, []

            def request(self, method, url, data=None, headers=None, timeout=None):
                self.sent.append(json.loads(data) if data else None)
                resp = type('R', (), {})()
                resp.status_code, resp.headers = 200, {}
                resp.content = json.dumps(self.response).encode()
                return resp
        s = Session(response3('budgetchain_r80'))
        out = TraverseClient('h', 1, session=s).deepen(self.g, 'right', [1], bp=100)
        self.assertEqual(1.0, s.sent[0]['strategy']['labels']['loss_budget'])
        self.assertTrue(any('conservative' in n for n in out.notes), out.notes)
        with tempfile.TemporaryDirectory() as d:
            store = GraphletStore(d)
            o = response3('budgetchain_r80')
            h = store.put_parsed(self.g, o['results'][0]['graphlet'], {'seeds': []})
            dry = GraphletTools(store).traverse_continue(h, 'right', 1, execute=False)
            self.assertEqual(1.0, dry['request']['strategy']['labels']['loss_budget'])
            self.assertTrue(any('conservative' in n for n in dry['notes']), dry)
            self.assertEqual(2.0, dry['loss_budget']['largest_terminal_loss'])


@unittest.skipUnless(_cli(), 'needs build_debug/metagraph')
class TestFinding4NeverDeeperThanUninterrupted(unittest.TestCase):
    """CLI: the review's fixture both ways. Uninterrupted (radius 180) walk 1's lineage
    goes E (0) -> C (1, switched at 80) -> F (2) -> G (3) and stops at 150 (loss_budget).
    Continued from 80 with 100 bp more: was 180 (H entered at a cumulative cost of 4);
    now 120 -- 30 bp shallower, conservative for E. Walk 0 (D alone, exact) reaches 90 in
    both."""

    def test_both_walks(self):
        req = copy.deepcopy(response3('budgetchain_r80')['strategy'])
        req.pop('clamped', None)
        base = {'seeds': [{'seed_id': 'switch_chain',
                           'sequence': self.seed(), 'labels': ['A', 'E']}],
                'strategy': req}
        g = _Indexes.graphlet('budgetchain', base)
        full = copy.deepcopy(base)
        full['strategy']['bounds']['max_extension_bp'] = 180
        whole = _Indexes.graphlet('budgetchain', full)
        reach = {w.path_id: w.length_bp for w in whole.walks('right', by='id')}
        self.assertEqual({0: 90, 1: 150}, reach)
        a = g.arms['right']
        got = {}
        for p in derive.paths(a):
            nxt = g.next_request('right', [p.id], bp=100)
            self.assertEqual([], T.request_violations(nxt))
            child = _Indexes.graphlet('budgetchain', nxt)
            got[p.id] = p.length_bp + max(q.length_bp for q in derive.paths(child.arms['right']))
        # walk 0 (D, exact) as far as uninterrupted; walk 1 (E, C) 30 bp shorter (was 180)
        self.assertEqual({0: 90, 1: 120}, got)
        for pid in got:
            self.assertLessEqual(got[pid], reach[pid])

    def seed(self):
        return graphlet3('budgetchain_r80').seed.sequence


# ------------------------------------------------------------------ finding 5

def _constrain_continuations(graphlets):
    for name, g in graphlets:
        if g.mode != 'constrain' or not g.has_envelope:
            continue
        for side, a in g.arms.items():
            for p in derive.paths(a):
                if a.segments[p.leaf].leaf.continuation is not None:
                    yield name, g, side, p.id


class TestFinding5ContinuationRequestsAreValid(unittest.TestCase):
    """A switched continuation promoted its labels into the seed but kept the original
    labels.extra: switch_chain walk 1 gave seed ['C'] with extra ['B', 'C', 'D'], a 400
    ('Extra label C duplicates a seed label'), from next_request(), deepen() and
    traverse_continue alike. labels.extra is now rebuilt around the new seed labels: the
    retrieval's permitted labels minus them (A, an original seed label, is a switch target
    now), each kept that a seed label reaches in one switch within the budget."""

    def test_the_reviewers_case(self):
        req = T.graphlet('switch_chain').next_request('right', [1], bp=50)
        self.assertEqual(['C'], req['seeds'][0]['labels'])
        self.assertEqual(['A', 'B', 'D'], req['strategy']['labels']['extra'])
        self.assertEqual([], T.request_violations(req))

    def test_every_continuation_of_every_fixture_is_valid(self):
        n = 0
        for name, g, side, pid in _constrain_continuations(fixture_graphlets()
                                                           + real_graphlets()):
            with self.subTest(doc=name, arm=side, walk=pid):
                try:
                    req = g.next_request(side, [pid])
                except ValueError as e:              # UnverifiableLabelName, no sequence
                    self.assertNotIsInstance(e, IncompatibleContinuations)
                    continue
                self.assertEqual([], T.request_violations(req))
                seeds = set(req['seeds'][0]['labels'])
                # every label of the retrieval is a seed label, an extra label or listed
                # as left out (the notes name the first few, under a cost that switches)
                left = {x['ref'] for x in req.left_out}
                for l in g.labels:
                    if l.name in seeds or l.name in req['strategy']['labels']['extra']:
                        self.assertNotIn(l.ref, left)
                    else:
                        self.assertIn(l.ref, left, (l, req.left_out))
                model = req['strategy']['labels']['change_cost']['model']
                if left and model != 'forbid':
                    self.assertTrue(any('leaves out' in x for x in req.notes), req.notes)
                n += 1
        self.assertGreater(n, 20)

    def test_walks_share_a_request_only_with_the_same_pool(self):
        sw = T.graphlet('switch_chain')
        with self.assertRaises(IncompatibleContinuations) as cm:
            sw.next_request('right', [0, 1])
        self.assertIn('next_requests()', str(cm.exception))
        self.assertIsInstance(cm.exception, ValueError)
        reqs = sw.next_requests('right', [0, 1])
        self.assertEqual([['D'], ['C']], [r['seeds'][0]['labels'] for r in reqs])
        for r in reqs:
            self.assertEqual([], T.request_violations(r))
        # forbid: no switch target at all, so the walks share one request as before
        two = T.graphlet('fork').next_request('left', [0, 1])
        self.assertEqual(2, len(two['seeds']))
        self.assertEqual([], two['strategy']['labels']['extra'])
        self.assertEqual([], T.request_violations(two))

    def test_an_extra_name_that_cannot_be_verified_is_left_out_and_stated(self):
        sw = T.graphlet('switch_chain')
        sw.labels[1].name = 'B�'                       # B: a replaced name
        req = sw.next_request('right', [1])
        self.assertEqual(['A', 'D'], req['strategy']['labels']['extra'])
        self.assertTrue(any('h:0:1' in n and 'U+FFFD' in n for n in req.notes), req.notes)

    def test_an_unreachable_target_is_left_out_and_stated(self):
        sw = T.graphlet('switch_chain')
        sw.envelope['strategy']['labels']['change_cost'] = {
            'model': 'table', 'default': 'forbid', 'entries': [['C', 'D', 0.5]]}
        req = sw.next_request('right', [1])
        self.assertEqual(['D'], req['strategy']['labels']['extra'])
        self.assertEqual([], T.request_violations(req))
        note, = [n for n in req.notes if 'leaves out' in n]
        self.assertIn("'A'", note)
        self.assertIn("'B'", note)


@unittest.skipUnless(_cli(), 'needs build_debug/metagraph')
class TestFinding5TheServerAcceptsThem(unittest.TestCase):
    """CLI: every continuation of every constrain CLI fixture, built per walk (and per
    group where the pools agree), is accepted by the server's own validation."""

    def test_every_continuation_is_accepted(self):
        n = 0
        for name, index, request in cli_requests():
            rc, out = _Indexes.run(index, request)
            self.assertEqual(0, rc, (name, out))
            g = from_response(out['results'][0], out)
            if g.mode != 'constrain':
                continue
            for side, a in g.arms.items():
                pids = [p.id for p in derive.paths(a)
                        if a.segments[p.leaf].leaf.continuation is not None]
                for pid in pids:
                    with self.subTest(fixture=name, arm=side, walk=pid):
                        req = g.next_request(side, [pid], bp=30)
                        rc, res = _Indexes.run(index, req)
                        self.assertEqual(0, rc, res)
                        self.assertNotIn('error', res['results'][0])
                        n += 1
                if len(pids) > 1:
                    try:
                        req = g.next_request(side, pids, bp=30)
                    except IncompatibleContinuations:
                        continue
                    rc, res = _Indexes.run(index, req)
                    self.assertEqual(0, rc, (name, side, res))
        self.assertGreater(n, 5)


# ------------------------------------------------------------------ finding 6

def truncated(g, depth):
    """|g| as a retrieval walked to |depth| would have recorded it, built from the model
    with the library's PRE-EXISTING accessors -- each run re-anchored on its own route as
    routes() reconstructs it, displayed evidence then by derive.evidence() on the copy --
    and not with the restricted derivations of compare() it is checked against: segments
    starting before the depth, cut there (a segment reaching it is a leaf, ended by the
    radius); runs that exist at the depth, cut there; P runs, events and branch events
    before it."""
    t = copy.deepcopy(g)
    t.cache = {}
    for side, a in t.arms.items():
        segs = a.segments
        route_of = {}
        if g.mode == 'constrain':
            per = collections.defaultdict(list)
            for r in g.arms[side].runs:
                per[r.label].append(r.id)
            for l, ids in per.items():
                for rid, route in zip(ids, ops.routes(g, {'id': l}, side)):
                    route_of[rid] = route
        ends = {s.id: s.end_bp for s in segs}
        keep = [s for s in segs if s.from_bp < depth]
        new = {s.id: i for i, s in enumerate(keep)}
        runs = []
        for r in a.runs:
            zero = r.from_bp == r.to_bp
            if (r.from_bp > depth) if zero else (r.from_bp >= depth):
                continue
            if r.to_bp > depth or segs[r.segment].from_bp >= depth:
                hold = [x for x in route_of[r.id]
                        if segs[x].from_bp <= depth - 1 < ends[x]]
                assert len(hold) == 1, (r, hold)
                r.segment = hold[0]
                r.to_bp = min(r.to_bp, depth)
                r.merged = r.silent = False
                r.end, r.reason, r.qualifier, r.needed_budget = 'X', 'max_extension_bp', None, None
            runs.append(r)
        rid = {r.id: i for i, r in enumerate(runs)}
        for r in runs:
            r.id, r.segment = rid[r.id], new[r.segment]
            r.prev_run = None if r.prev_run is None else rid[r.prev_run]
        for s in keep:
            end = min(ends[s.id], depth)
            s.id = new[s.id]
            s.parents = tuple(new[p] for p in s.parents)
            s.length_bp = end - s.from_bp
            if s.walk is not None:
                s.walk = s.walk[:s.length_bp]
            s.presence = [PresenceRun(p.from_bp, min(p.to_bp, depth), p.total, p.labels)
                          for p in s.presence if p.from_bp < depth]
            s.events = [e for e in s.events if e.at_bp < depth
                        and (e.segment is None or e.segment in new)]
            for e in s.events:
                if e.segment is not None:
                    e.segment = new[e.segment]
            s.children = []
        for s in keep:
            for p in s.parents:
                keep[p].children.append(s.id)
        by_seg = collections.defaultdict(list)
        for r in runs:
            by_seg[r.segment].append(r.id)
        for s, orig in zip(keep, [x for x in sorted(ends) if x in new]):
            if s.children:
                s.leaf = None
            else:
                s.split = None
                if s.leaf is None or ends[orig] > depth:
                    s.leaf = Leaf('X', {}, None)
            if g.mode == 'constrain':
                s.end = derive.rule_end(g.mode, s, s.entry, s.presence, by_seg[s.id], runs)
            else:
                s.end = list(s.presence[-1].labels) if s.presence else list(s.entry)
        a.segments, a.runs = keep, runs
        a.branch_events = [be for be in a.branch_events
                           if be.at_bp < depth and be.segment in new]
        for be in a.branch_events:
            be.segment = new[be.segment]
        a.complete_to_bp = min(a.complete_to_bp, depth)
        if a.evidence_complete_to_bp is not None:
            a.evidence_complete_to_bp = min(a.evidence_complete_to_bp, depth)
        a.cache = {}
    return t


def _differences(c, mode):
    if mode == 'prefix_subset':
        return c.only_in_a or c.only_in_b
    return c.only_in_a or c.only_in_b or c.differ


class TestFinding6CompareOverTheRestrictedDag(unittest.TestCase):
    """compare() derived its keys from the deeper DAG's displayed walks and its claims cut
    at the depth: on the review's bubbles (annotate, radius 30 vs 100, compared at 30) the
    merges at 36 and 62 made claims, walks and prefix_subset report differences (b.fa,
    c.fa and both.fa on the branch that enters the merge at 36 through its non-first
    parent: no displayed walk of the 100 bp retrieval shows it). Both DAGs are now
    restricted to the depth before anything is keyed."""

    def test_the_reviewers_case(self):
        small, big = graphlet3('bubbles_annotate_r30'), graphlet3('bubbles_annotate_r100')
        for mode in MODES:
            for a, b in ((small, big), (big, small)):
                c = a.compare(b, mode=mode)
                with self.subTest(mode=mode, a=a.arms['right'].complete_to_bp):
                    self.assertEqual((True, 30, True), (c.comparable, c.depth_used, c.equal),
                                     c.as_dict())
                    if mode != 'labels':
                        self.assertTrue(any('restricted to [0, 30)' in n for n in c.notes))

    def test_claims_at_a_cut_keep_the_normative_cut(self):
        # what the deeper DAG says about [0, 30) (§5.1: evidence at the anchored endpoint,
        # then clipped) is unchanged; compare() keys the restricted DAG's claims instead
        big = graphlet3('bubbles_annotate_r100')
        cut = [c for c in big.claims('right', at_most_bp=30) if c.label.name == 'c.fa']
        self.assertEqual(['route_only'], [c.kind for c in cut])
        mine = [c for c in ops._restricted_claims(big, 'right', 30) if c.label.name == 'c.fa']
        self.assertEqual([('stretch', 0)], [(c.kind, c.evidence_from) for c in mine])

    def test_the_restricted_leaves(self):
        big = graphlet3('bubbles_annotate_r100')
        a = big.arms['right']
        leaves = ops.restricted_leaves(a, 30)
        shown = {s for p in derive.paths(a) for s in derive.chain(a, p.leaf)}
        # one of them is on no displayed walk of the deeper DAG
        self.assertTrue([x for x in leaves if x not in shown])
        self.assertEqual(sorted(p.leaf for p in derive.paths(a)),
                         ops.restricted_leaves(a, 10 ** 6))

    def test_an_unreconstructible_route_is_qualified_not_a_difference(self):
        g = T.graphlet('merge')
        g.index_fp = 'f' * 64                  # verifiable: only 'qualified' can weaken it
        t = truncated(g, 50)
        self.assertEqual((True, True), (t.compare(g, mode='claims').comparable,
                                        t.compare(g, mode='claims').equal))
        a = g.arms['right']
        m = next(s for s in a.segments if len(s.parents) > 1 and s.from_bp == 62)
        for part in m.partition:
            for i in range(len(part) - 1, -1, -1):
                if g.labels[part[i]].name == 'b.fa':
                    del part[i]
        g.cache.clear()
        a.cache.clear()
        c = t.compare(g, mode='claims')
        self.assertEqual(('qualified', None), (c.comparable, c.equal))
        self.assertTrue(any('could not be reconstructed exactly' in n for n in c.notes),
                        c.notes)

    def depths(self, g, rng):
        out = set()
        for a in g.arms.values():
            top = a.complete_to_bp
            if top < 2:
                continue
            out.update(rng.sample(range(1, top), min(3, top - 1)))
            for s in a.segments:
                if len(s.parents) > 1:
                    out.update(x for x in (s.from_bp - 1, s.from_bp, s.from_bp + 1)
                               if 0 < x < top)
        return sorted(out)

    def check_truncations(self, graphlets, rng):
        n = 0
        for name, g in graphlets:
            if any(s.walk is None for a in g.arms.values() for s in a.segments):
                continue                          # no bases: 'unknown' by design
            for depth in self.depths(g, rng):
                t = truncated(g, depth)
                for mode in MODES:
                    for x, y in ((t, g), (g, t)):
                        c = x.compare(y, mode=mode)
                        with self.subTest(doc=name, depth=depth, mode=mode):
                            self.assertNotEqual(False, c.comparable)
                            self.assertFalse(_differences(c, mode), c.as_dict())
                            if c.comparable is True:
                                self.assertTrue(c.equal)
                        n += 1
        return n

    def test_every_fixture_equals_itself_truncated(self):
        self.assertGreater(self.check_truncations(fixture_graphlets(), random.Random(3)), 500)

    def test_every_offline_real_retrieval_equals_itself_truncated(self):
        self.assertGreater(self.check_truncations(real_graphlets(), random.Random(5)), 300)


@unittest.skipUnless(_cli(), 'needs build_debug/metagraph')
class TestFinding6RealRetrievalsAtRandomRadii(unittest.TestCase):
    """CLI: every CLI fixture's request at a few random smaller radii against the
    fixture's own radius: a retrieval and a deeper one of the same seed never differ at
    the shallower depth (before the fix: 34 of 512 comparisons did, at seed 7)."""

    def test_no_false_difference(self):
        rng = random.Random(11)
        n = 0
        for name, index, request in cli_requests():
            radius = request['strategy'].get('bounds', {}).get('max_extension_bp', 5000)
            deep = _Indexes.graphlet(index, request)
            for depth in sorted({1, max(1, radius // 2), rng.randint(1, max(1, radius - 1))}):
                r = copy.deepcopy(request)
                r['strategy']['bounds']['max_extension_bp'] = depth
                shallow = _Indexes.graphlet(index, r)
                for mode in MODES:
                    for x, y in ((shallow, deep), (deep, shallow)):
                        c = x.compare(y, mode=mode)
                        with self.subTest(fixture=name, depth=depth, mode=mode):
                            if mode == 'prefix_subset' and x is deep:
                                self.assertFalse(c.only_in_a, c.as_dict())
                            else:
                                self.assertFalse(_differences(c, mode), c.as_dict())
                            self.assertNotEqual(False, c.equal)
                        n += 1
        self.assertGreater(n, 100)


# ------------------------------------------------------------------ N2

class TestN2BindingErrorsHoldTheCeiling(unittest.TestCase):
    """graphlet_list(max_bytes=64, foo=1) answered 77 bytes and graphlet_summary(
    max_bytes=64) 74: an error raised while binding the arguments was fitted to the
    default ceiling, not to the one the caller supplied."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.tools = GraphletTools(GraphletStore(self.tmp.name))

    def tearDown(self):
        self.tmp.cleanup()

    def test_the_reviewers_calls(self):
        for name, args in (('graphlet_list', {'max_bytes': 64, 'foo': 1}),
                           ('graphlet_summary', {'max_bytes': 64}),
                           ('graphlet_walks', {'max_bytes': 70, 'handle': 'g_x'}),
                           ('traverse_continue', {'max_bytes': 64, 'bogus': True})):
            out = getattr(self.tools, name)(**args)
            with self.subTest(tool=name):
                self.assertEqual('bad_argument', out['error'])
                self.assertLessEqual(_size(out), args['max_bytes'], out)

    def test_a_positional_ceiling_holds_too(self):
        out = self.tools.graphlet_list(None, MIN_MAX_BYTES, 'extra')
        self.assertEqual('bad_argument', out['error'])
        self.assertLessEqual(_size(out), MIN_MAX_BYTES)

    def test_an_invalid_ceiling_falls_back_to_the_default(self):
        out = self.tools.graphlet_list(max_bytes=8, foo=1)
        self.assertEqual('bad_argument', out['error'])
        self.assertLessEqual(_size(out), self.tools.max_bytes)
        self.assertIn('foo', out['message'])


# ------------------------------------------------------------------ N3

class TestN3LiveOracleComparesTheManifest(unittest.TestCase):
    """realdata.live_matches_cache() compared index_meta_fp only: a live mini_refseq whose
    index_fp differs from the recorded one was accepted, and the /resolve oracle ran
    against it. Every fingerprint both sides state is compared now; a digest on one side
    only cannot be verified; with no manifest on either side equal index_meta_fp is the
    (negative) check there is, stated as unverifiable. Offline: the capability reads are
    replaced, no socket is opened."""

    @classmethod
    def setUpClass(cls):
        import test_traverse_real_offline as off
        cls.h = off.H
        cls.R = cls.h.R

    def patched(self, index, **changes):
        R = self.R
        rec = R.load_manifest()['servers'][index]
        live = dict(copy.deepcopy(rec), **changes)
        saved = R.server_up, R.live_capabilities
        R.server_up = lambda i: True
        R.live_capabilities = lambda i: live
        self.addCleanup(lambda: (setattr(R, 'server_up', saved[0]),
                                 setattr(R, 'live_capabilities', saved[1])))
        return rec

    def test_a_different_digest_is_refused(self):
        rec = self.patched('mini_refseq', index_fp='f' * 64)
        self.assertTrue(rec['index_fp'] and rec['index_fp'] != 'f' * 64)
        self.assertFalse(self.R.live_matches_cache('mini_refseq'))
        self.assertEqual('different', self.R.live_identity('mini_refseq')[0])
        with self.assertRaises(unittest.SkipTest) as cm:
            self.h.O._require_oracle('mini_refseq')
        self.assertIn('index_fp differs', str(cm.exception))

    def test_the_same_digest_is_accepted(self):
        self.patched('mini_refseq')
        self.assertEqual('same', self.R.live_identity('mini_refseq')[0])
        self.assertTrue(self.R.live_matches_cache('mini_refseq'))

    def test_a_digest_on_one_side_only_is_unverifiable(self):
        self.patched('mini_refseq', index_fp=None)
        status, why = self.R.live_identity('mini_refseq')
        self.assertEqual('unverifiable', status)
        self.assertIn('negative check', why)
        self.assertFalse(self.R.live_matches_cache('mini_refseq'))

    def test_without_manifests_equal_metadata_is_the_stated_check(self):
        rec = self.patched('sra')
        self.assertFalse(rec.get('index_fp'))
        status, why = self.R.live_identity('sra')
        self.assertEqual('unverifiable', status)
        self.assertIn('neither', why)
        self.assertTrue(self.R.live_matches_cache('sra'))
        self.patched('sra', index_meta_fp='0' * 16)
        self.assertFalse(self.R.live_matches_cache('sra'))


# ------------------------------------------------------------------ D1, D7

class TestToolDescriptionsStateTheirBudgets(unittest.TestCase):
    """The descriptions of the compare and export tools, which an agent reads: without a
    budget they have no work or allocation budget, and a budget is there to be passed."""

    def test_compare_and_export(self):
        from metagraph.traverse import mcp_tools
        for tool in (mcp_tools.GraphletTools.graphlet_compare,
                     mcp_tools.GraphletTools.graphlet_export):
            flat = ' '.join(tool.__doc__.split())
            self.assertRegex(flat, r'(?i)no work (or|and) allocation budget', tool.__name__)
            self.assertRegex(flat, r'(?i)without (a budget|local limits|local_limits)',
                             tool.__name__)


# ------------------------------------------------------------------ D4

class _Caps:
    """A backend double: capabilities with a settable index_fp, and the switch_chain
    fixture's graphlet response carrying the same digest."""

    def __init__(self, fp):
        self.fp = fp
        self.caps = {'index_ns': None, 'index_fp': fp, 'index_meta_fp': None, 'release': None}
        self.requests = []

    def capabilities(self):
        return dict(self.caps)

    def build_request(self, seeds, strategy, detail='graphlet'):
        return {'seeds': seeds, 'strategy': strategy}

    def traverse_raw(self, req):
        self.requests.append(req)
        resp = copy.deepcopy(T.doc_json('switch_chain', 'graphlet'))
        body = resp['results'][0]['graphlet']
        head, rest = body.split('\n', 1)
        f = head.split(' ')
        f[-2] = self.fp or '*'                       # H: ... <index_fp> <index_meta_fp>
        resp['results'][0]['graphlet'] = ' '.join(f) + '\n' + rest
        resp['results'][0].pop('graphlet_bytes', None)
        return resp


class TestD4OneSidedManifestIsUnverifiable(unittest.TestCase):
    """A manifest digest on one side only proves neither identity nor difference: it was
    refused as index_mismatch (proven different). It is now index_unverifiable, still
    refused automatically, and run only when the caller passes allow_unverified_index=
    true -- the result then states the identity unverified and accepted."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.store = GraphletStore(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def fetch(self, client):
        tools = GraphletTools(self.store, {'mini': client})
        out = tools.traverse_fetch(index='mini',
                                   seed={'sequence': T.graphlet('switch_chain').seed.sequence})
        self.assertIn('handle', out, out)
        return tools, out['handle']

    def head_fp(self):
        body = T.doc_json('switch_chain', 'graphlet')['results'][0]['graphlet']
        return body.split('\n', 1)[0].split(' ')

    def test_unverifiable_is_refused_without_the_opt_in(self):
        client = _Caps(fp='a' * 64)
        tools, parent = self.fetch(client)
        client.fp = client.caps['index_fp'] = None
        for out in (tools.traverse_continue(parent, 'right', 1),
                    tools.traverse_fetch(replay=parent)):
            self.assertEqual('index_unverifiable', out['error'], out)
            self.assertIn('not proven different', out['message'])
            self.assertIn('allow_unverified_index', out['hint'])
        self.assertEqual(1, len(client.requests))           # nothing sent

    def test_the_opt_in_runs_it_and_says_so(self):
        client = _Caps(fp='a' * 64)
        tools, parent = self.fetch(client)
        client.fp = client.caps['index_fp'] = None
        out = tools.traverse_continue(parent, 'right', 1, allow_unverified_index=True)
        self.assertIn('handle', out, out)
        self.assertEqual((False, 'allow_unverified_index'),
                         (out['identity']['verified'], out['identity']['accepted']))
        rep = tools.traverse_fetch(replay=parent, allow_unverified_index=True)
        self.assertFalse(rep['identity']['verified'])

    def test_a_proven_difference_is_refused_even_with_the_opt_in(self):
        client = _Caps(fp='a' * 64)
        tools, parent = self.fetch(client)
        client.fp = client.caps['index_fp'] = 'b' * 64
        out = tools.traverse_continue(parent, 'right', 1, allow_unverified_index=True)
        self.assertEqual('index_mismatch', out['error'])

    def test_the_opt_in_is_a_boolean(self):
        client = _Caps(fp=None)
        tools, parent = self.fetch(client)
        out = tools.traverse_continue(parent, 'right', 1, allow_unverified_index='yes')
        self.assertEqual('bad_argument', out['error'])


# ------------------------------------------------------------------ D5

class TestD5WalksCountMergeEnteredWalks(unittest.TestCase):
    """graphlet_walks counted the walks its route_consistent filter removed as
    route_only (a claim's kind at a cut) while graphlet_claims counted the same situation
    as merge_entered. Both say merge_entered now; an annotate walk whose label is recorded
    at its end but on no route from the seed boundary is not_from_seed."""

    def tools(self, name):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        store = GraphletStore(tmp.name)
        resp = T.doc_json(name, 'graphlet')
        h = store.put_parsed(T.graphlet(name), resp['results'][0]['graphlet'], None)
        return GraphletTools(store), h

    def test_constrain(self):
        tools, h = self.tools('merge')
        w = tools.graphlet_walks(h, 'right', label={'ref': 'c:1'})
        self.assertEqual({'merge_entered': 1}, w['filtered'])
        self.assertIn('non-first merge parent', w['hint'])
        c = tools.graphlet_claims(h, 'right', max_bytes=8192)
        self.assertEqual({'merge_entered': 2}, c['filtered'])
        self.assertIn('non-first merge parent', c['hint'])
        # the claims it removed are alive (displayed support from the merge on), not
        # route_only, the kind of a claim with no displayed support at a cut
        removed = [x for x in T.graphlet('merge').claims('right') if x.route_bp]
        self.assertEqual({'alive'}, {x.kind for x in removed})
        self.assertTrue(all(x.evidence_from < x.to_bp for x in removed))

    def test_annotate(self):
        g = T.graphlet('annotate')
        seen = collections.Counter()
        for side, a in g.arms.items():
            for l in g.labels:
                every = set(ops.rank_walks(g, side, labels=[{'id': l.id}],
                                           route_consistent=False))
                kept = set(ops.rank_walks(g, side, labels=[{'id': l.id}]))
                for pid in every - kept:
                    why = ops.walk_filter_reason(g, side, pid, [{'id': l.id}])
                    seen[why] += 1
                    alive_end, _ = ops._annotate_route_ends(a)
                    self.assertEqual(why == 'merge_entered',
                                     l.id in alive_end[derive.paths(a)[pid].leaf])
        self.assertTrue(seen, 'the annotate fixture filters some walk')


# ------------------------------------------------------------------ D6

class TestD6ReceiptsAreNeverReplaced(unittest.TestCase):
    """An operation that made something answered result_too_large when its answer did not
    fit, and the caller never learned the handle or the path. Its optional fields go
    first now (named in fields_cut); a receipt that could not fit even then is refused
    before anything is stored or written."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.store = GraphletStore(self.tmp.name)
        resp = T.doc_json('merge', 'graphlet')
        self.h = self.store.put_parsed(T.graphlet('merge'), resp['results'][0]['graphlet'],
                                       {'seeds': [{'sequence': 'ACGT'}]})

    def tools(self, max_bytes, export_dir=None):
        return GraphletTools(self.store, max_bytes=max_bytes, export_dir=export_dir)

    def test_load_keeps_its_handle(self):
        big = self.tools(4096)
        saved = big.graphlet_save(self.h, 'm.mgt')
        self.assertIn('path', saved)
        small = self.tools(MIN_MAX_BYTES)
        n = len(self.store.list())
        out = small.graphlet_load(os.path.basename(saved['path']))
        self.assertTrue(out['handle'].startswith('g_'), out)
        self.assertEqual(['summary'], out['fields_cut'])
        self.assertLessEqual(_size(out), MIN_MAX_BYTES)
        self.assertEqual(n + 1, len(self.store.list()))

    def test_export_and_save_refuse_before_writing(self):
        with tempfile.TemporaryDirectory() as d:
            root = os.path.realpath(d)
            limit = _size({'path': os.path.join(root, 'a.fa'),
                           'fields_cut': ['records', 'bytes']})
            tools = self.tools(limit, export_dir=d)
            name = 'x' * 40 + '.fa'
            for out in (tools.graphlet_export(self.h, path=name),
                        tools.graphlet_save(self.h, name)):
                self.assertEqual('receipt_too_large', out['error'], out)
            self.assertEqual([], os.listdir(d))
            ok = tools.graphlet_export(self.h, path='a.fa')
            self.assertEqual(os.path.join(root, 'a.fa'), ok['path'], ok)
            self.assertLessEqual(_size(ok), limit)
            self.assertEqual(['a.fa'], os.listdir(d))

    def test_the_receipt_drops_optional_fields_in_order(self):
        from metagraph.traverse.mcp_tools import _Receipt
        r = _Receipt({'path': 'p', 'bytes': 10 ** 15, 'records': 10 ** 15},
                     ('records', 'bytes'))
        full = _size(dict(r))
        self.assertEqual(dict(r), r.fit(full))
        one = r.fit(full - 1)
        self.assertEqual((['records'], 10 ** 15), (one['fields_cut'], one['bytes']))
        both = r.fit(_size(r.minimal()))
        self.assertEqual({'path': 'p', 'fields_cut': ['records', 'bytes']}, both)

    def test_subtrie_keeps_its_handle(self):
        tools = self.tools(MIN_MAX_BYTES * 2)
        out = tools.graphlet_subtrie(self.h, ['a.fa', 'b.fa'])
        self.assertTrue(out['handle'].startswith('g_'), out)
        self.assertIn('summary', out['fields_cut'])
        self.assertLessEqual(_size(out), MIN_MAX_BYTES * 2)

    def test_fetch_keeps_its_handle_or_stores_nothing(self):
        client = _Caps(fp=None)
        tools = GraphletTools(self.store, {'mini': client})
        n = len(self.store.list())
        out = tools.traverse_fetch(index='mini', seed={'sequence': 'ACGT'}, max_bytes=120)
        self.assertTrue(out['handle'].startswith('g_'), out)
        self.assertIn('summary', out['fields_cut'])
        self.assertLessEqual(_size(out), 120)
        self.assertEqual(n + 1, len(self.store.list()))
        out = tools.traverse_fetch(index='mini', seed={'sequence': 'ACGT'},
                                   max_bytes=MIN_MAX_BYTES)
        self.assertEqual('receipt_too_large', out['error'], out)
        self.assertEqual(n + 1, len(self.store.list()))    # nothing stored


# ------------------------------------------------------------------ round 3 (the fixes' review)

class TestRound3AMergeAtTheComparisonDepth(unittest.TestCase):
    """compare() still reported false differences when the comparison depth equalled a
    merge position: a retrieval walked to D records the merge at D as a zero-length
    segment, its runs closed there as merged and its runs ending there anchored on it; the
    restricted derivation kept that merge (outside [0, D)). The 'merge' CLI fixture
    (bubbles, merges at 36 and 62) walked to 36 against the same request to 100 compared
    unequal in claims; on mini_refseq 90 of 2,016 comparisons, all at D == a merge. Data:
    merge_constrain_r36 / _r100, and bubbles_annotate_r36 / _r62 (the annotate fixture's
    request as bubbles_annotate_r100)."""

    def test_a_radius_at_a_merge_compares_equal(self):
        for deep, shallows in (('merge_constrain_r100', ('merge_constrain_r36',)),
                               ('bubbles_annotate_r100', ('bubbles_annotate_r36',
                                                          'bubbles_annotate_r62'))):
            deep = graphlet3(deep)
            for name in shallows:
                shallow = graphlet3(name)
                r = shallow.arms['right'].complete_to_bp
                mode = shallow.mode
                merges = [s for s in shallow.arms['right'].segments
                          if len(s.parents) > 1 and s.from_bp == r]
                self.assertTrue(merges, 'the radius-%d retrieval ends in a merge' % r)
                for m in MODES:
                    for a, b in ((shallow, deep), (deep, shallow)):
                        c = a.compare(b, mode=m)
                        with self.subTest(mode=mode, radius=r, compare=m):
                            self.assertEqual((True, r, True),
                                             (c.comparable, c.depth_used, c.equal), c.as_dict())

    def test_a_run_closed_at_the_depth_is_open_there(self):
        deep = graphlet3('merge_constrain_r100')
        a = deep.arms['right']
        closed = [r for r in a.runs if r.merged and r.to_bp == 36]
        self.assertTrue(closed)
        got = {c.run: c for c in ops._restricted_claims(deep, 'right', 36)}
        for r in closed:
            self.assertEqual('open', got[r.id].end_class)
            self.assertNotEqual('merged', got[r.id].kind)

    def test_true_differences_at_a_merge_depth_are_kept(self):
        # no branch allowance against an unlimited one: a real difference before 36
        shallow, b0 = graphlet3('merge_constrain_r36'), graphlet3('merge_constrain_b0_r100')
        for m in MODES:
            c = shallow.compare(b0, mode=m)
            with self.subTest(compare=m):
                self.assertEqual((True, 36, False), (c.comparable, c.depth_used, c.equal))


class TestRound3BMalformedOverridesAreResults(unittest.TestCase):
    """The rebuilt next_request() read the merged override sections as objects: an
    override such as {'labels': 5} raised AttributeError, which the tool wrapper does not
    catch, where HEAD answered a bounded backend_error. An unhashable arm (graphlet_walks,
    graphlet_summary, ... traverse_continue) and traverse_resolve(labels=5) raised
    TypeError. Each is now a bounded bad_argument or bad_arm."""

    BAD = ({'labels': 5}, {'labels': {'change_cost': 5}}, {'branching': 5},
           {'labels': {'extra': 'B'}}, {'labels': {'loss_budget': 'x'}},
           {'labels': {'change_cost': {'model': 'constant', 'value': None}}},
           {'labels': {'change_cost': {'model': 'table', 'entries': 'x'}}},
           {'labels': {'change_cost': {'model': 'table', 'entries': [['A', 'B', 1]],
                                       'default': 'x'}}})

    def setUp(self):
        self.g = graphlet3('budgetchain_r80')
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.store = GraphletStore(self.tmp.name)
        resp = response3('budgetchain_r80')
        self.h = self.store.put_parsed(self.g, resp['results'][0]['graphlet'],
                                       {'seeds': [{'sequence': 'ACGT'}]})

    def test_next_request_refuses_with_a_value_error(self):
        for ov in self.BAD:
            with self.subTest(overrides=ov):
                with self.assertRaises(ValueError):
                    self.g.next_request('right', [1], **copy.deepcopy(ov))

    def test_the_tools_answer_within_the_ceiling(self):
        tools = GraphletTools(self.store, {'x': _Caps(fp=None)})
        for ov in self.BAD:
            for execute in (False, True):
                out = tools.traverse_continue(self.h, 'right', 1, overrides=copy.deepcopy(ov),
                                              execute=execute, max_bytes=MIN_MAX_BYTES)
                with self.subTest(overrides=ov, execute=execute):
                    self.assertEqual('bad_argument', out['error'], out)
                    self.assertLessEqual(_size(out), MIN_MAX_BYTES)
        for arm in ({'a': 1}, ['right'], 3):
            for name, kw in (('graphlet_walks', {}), ('graphlet_summary', {}),
                             ('graphlet_claims', {}), ('graphlet_splits', {}),
                             ('graphlet_labels', {}), ('graphlet_support', {'walk': 0}),
                             ('traverse_continue', {'walk': 0})):
                out = getattr(tools, name)(self.h, arm=arm, max_bytes=MIN_MAX_BYTES, **kw)
                with self.subTest(tool=name, arm=arm):
                    self.assertEqual('bad_arm', out['error'], out)
                    self.assertLessEqual(_size(out), MIN_MAX_BYTES)
        out = tools.traverse_resolve(index='x', sequence='ACGT', labels=5,
                                     max_bytes=MIN_MAX_BYTES)
        self.assertEqual('bad_argument', out['error'], out)
        self.assertLessEqual(_size(out), MIN_MAX_BYTES)


class _SameIndex(_Caps):
    """A backend double answering a continuation with |name|'s response itself, under its
    own index identity (so that a stored entry of that response verifies)."""

    def __init__(self, name):
        self.resp = response3(name)
        super().__init__(self.resp['capabilities'].get('index_fp'))
        self.caps = {k: self.resp['capabilities'].get(k) for k in
                     ('index_ns', 'index_fp', 'index_meta_fp', 'release')}

    def traverse_raw(self, req):
        self.requests.append(req)
        return copy.deepcopy(self.resp)


class TestRound3CContinueReceiptsFitWideWalks(unittest.TestCase):
    """traverse_continue made the per-label loss budget an essential part of its receipt,
    for every constrain request, even without a loss budget (each label at loss 0,
    remaining 0): a continuation of a walk carrying 30 labels was refused with
    receipt_too_large at the default 2,048-byte ceiling (3,393 bytes, 3,121 of them the
    list), where HEAD returned the handle; and the refusal told it to use a shorter file
    name. The list is now stated only where a loss budget applies, it is the first
    optional field cut (fields_cut), and each tool names its own lever."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.store = GraphletStore(self.tmp.name)

    def entry(self, name):
        resp = response3(name)
        g = from_response(resp['results'][0], resp)
        return self.store.put_parsed(g, resp['results'][0]['graphlet'],
                                     {'seeds': [{'sequence': 'ACGT'}]}, source='x'), g

    def test_a_walk_of_thirty_labels_keeps_its_handle(self):
        h, g = self.entry('wide_forbid_r40')
        self.assertEqual(30, len(g.continuation('right', 0).labels))
        self.assertIsNone(g.next_request('right', [0]).loss_budget)
        tools = GraphletTools(self.store, {'x': _SameIndex('wide_forbid_r40')})
        out = tools.traverse_continue(h, 'right', 0)
        self.assertTrue(out.get('handle', '').startswith('g_'), out)
        self.assertNotIn('loss_budget', out)
        self.assertNotIn('loss_budget_labels', out)
        self.assertLessEqual(_size(out), 2048)

    def test_the_per_label_list_is_cut_first(self):
        h, g = self.entry('budgetchain_r80')
        req = g.next_request('right', [1])
        self.assertEqual(3.0, req.loss_budget['original'])
        tools = GraphletTools(self.store, {'x': _SameIndex('budgetchain_r80')})
        full = tools.traverse_continue(h, 'right', 1, max_bytes=1 << 20)
        self.assertIn('loss_budget_labels', full)
        self.assertNotIn('labels', full['loss_budget'])
        cut_first = 0
        for limit in range(200, _size(full) + 1, 7):
            out = tools.traverse_continue(h, 'right', 1, max_bytes=limit)
            if 'error' in out:
                self.assertEqual('receipt_too_large', out['error'], out)
                continue
            with self.subTest(limit=limit):
                self.assertTrue(out['handle'].startswith('g_'), out)
                self.assertLessEqual(_size(out), limit)
                # the essential statement stays whole; the per-label list goes first
                self.assertEqual(full['loss_budget'], out['loss_budget'])
                self.assertEqual(full['notes'], out['notes'])
                if 'fields_cut' in out:
                    self.assertEqual('loss_budget_labels', out['fields_cut'][0])
                    cut_first += 1
        self.assertGreater(cut_first, 0)
        # execute=False carries both, the request under the sequence ceiling
        dry = tools.traverse_continue(h, 'right', 1, execute=False)
        self.assertEqual(full['loss_budget_labels'], dry['loss_budget_labels'])

    def test_each_tool_names_its_own_lever(self):
        h, g = self.entry('budgetchain_r80')
        tools = GraphletTools(self.store, {'x': _SameIndex('budgetchain_r80')})
        out = tools.traverse_continue(h, 'right', 1, max_bytes=MIN_MAX_BYTES)
        self.assertEqual('receipt_too_large', out['error'], out)
        long = tools.traverse_continue(h, 'right', 1, max_bytes=300)
        self.assertEqual('receipt_too_large', long['error'], long)
        self.assertNotIn('file name', long['message'])
        self.assertIn('raise max_bytes', long['message'])


class TestRound3DAliveLabelsThatAreNotSeededAreStated(unittest.TestCase):
    """A label alive at the leaf that switched in within the last bases does not cover
    the continuation's whole tail, so it is not among C's labels and the request does not
    seed it: under switch_on: loss the continuation can enter it only where a seed label
    is lost, and may lack the lineage one uninterrupted walk keeps. next_request() said
    nothing (notes [], left_out without it). Data: budgetchain_ae_r31 (seed labels A and
    E, extra B C D F G H, radius 31: B switched in at 30 and is alive at the leaf at loss
    1, C's labels are [E])."""

    def test_the_label_is_stated(self):
        g = graphlet3('budgetchain_ae_r31')
        a = g.arms['right']
        leaf = next(p for p in derive.paths(a) if a.segments[p.leaf].leaf.continuation)
        alive = {g.labels[e.label].name for e in derive.end_labels(a, leaf.leaf)}
        cont = {l.name for l in g.continuation('right', leaf.id).labels}
        self.assertEqual({'B'}, alive - cont)
        req = g.next_request('right', [leaf.id])
        self.assertEqual(['E'], req['seeds'][0]['labels'])
        (entry,) = [x for x in req.left_out if x['why'] == 'alive_not_seeded']
        self.assertEqual(('B', leaf.id, 1.0), (entry['name'], entry['walk'], entry['loss']))
        self.assertTrue(any("'B'" in n and 'not seeded' in n for n in req.notes), req.notes)
        # the request itself is unchanged by the statement: B stays a switch target
        self.assertIn('B', req['strategy']['labels']['extra'])

    def test_a_walk_whose_labels_all_cover_the_tail_states_nothing_of_it(self):
        g = graphlet3('budgetchain_r80')
        for walk in (0, 1):
            req = g.next_request('right', [walk])
            self.assertFalse([x for x in req.left_out if x['why'] == 'alive_not_seeded'])
            self.assertFalse([n for n in req.notes if 'not seeded' in n])


if __name__ == '__main__':
    unittest.main()
