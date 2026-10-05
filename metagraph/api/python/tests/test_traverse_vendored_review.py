"""The search service's code review of 2026-10-05, its part on the vendored library (31
findings VOP1-*, VOP2-*, VPC-*, VMT-*, VMD-*, reviewed at e7d99e0c): a regression test per
finding that still held, its offline reproduction restated here (the service's repository
is not read at test time), and a guard for the ones fixed before.

  VOP2-03  a comparison's cut (annotate) starts from the seed boundary's set
  VOP2-01  a view's walks are ranked, then cut to top
  VMD-06   a replay and a continuation keep the graph of a multi-graph server
  VMD-01   a handle in steady use survives a restart (the entry's last use persisted)
  VPC-04   the RANGES writer refuses every negative id
  VPC-05   a J line with one empty side reads back as written
  VMT-04   graphlet_sequence(segment=...) takes no with_seed
  VOP1-01  an infinite or NaN branch allowance is a named ValueError (dcc0cebd)
  VOP1-02  {'id': True} is no label id
  VOP1-03  walks() builds the claims of its own leaves only (f667d775)
  VOP1-04  _switch_cost, never called, is gone
  VOP2-02  compare(arm=X) of a pair where b lacks X answers instead of a KeyError
  VOP2-04  a view's backing digest is computed only when none is given
  VOP2-05  a view's segments in one pass (no chain per hit)
  VOP2-07  a view needs a label; a saved view without one is a format error
  VOP2-08  compare()'s early answers state both support kinds, after the mode check
  VOP2-09  one restoration of a saved view (parser.load and the store)
  VPC-01   a non-integer transport count is a GraphletFormatError
  VPC-02   from_response and standalone_text share the transport checks
  VPC-03   standalone_text makes no extra copy of the body
  VMD-02   processes sharing a spool: no body deleted under another's entry; free, `in`
           and list() see the other process's entries (review of these fixes, finding 5)
  VMD-03   an entry file of another version does not stop the store opening
  VMD-04   a missing body expires its entry, also on a copy without a parse (review of
           these fixes, finding 4)
  VMD-05   to_gfa's context read back as far as needed (f667d775)
  VMD-07   a body that is not UTF-8, or does not decode, is a TraverseError
  VMD-08   one statement of the label changes inside a segment
  VMT-01   traverse_continue's own arguments are no overrides (f667d775)
  VMT-02   a page is filled in linear time; max_bytes has a ceiling
"""

import copy
import gzip
import json
import math
import os
import random
import sys
import tempfile
import time
import tracemalloc
import gc
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402
import traverse_golden as TG  # noqa: E402
from test_traverse_mcp import CONTINUATION, SEED, FakeClient  # noqa: E402
from test_traverse_scale import comb_annotate  # noqa: E402

import metagraph.traverse as MT  # noqa: E402
from metagraph.traverse import (  # noqa: E402
    GraphletStore, GraphletView, TraverseClient, derive, ops, parse, parser,
)
from metagraph.traverse import store as store_mod  # noqa: E402
from metagraph.traverse._codec import (  # noqa: E402
    CodecError, GraphletFormatError, decode_ranges, encode_ranges,
)
from metagraph.traverse.client import ServerInitializing, TraverseError  # noqa: E402
from metagraph.traverse.mcp_tools import (  # noqa: E402
    MAX_MAX_BYTES, GraphletTools, ToolLimits, _size,
)
from metagraph.traverse.model import BadSelector  # noqa: E402
from metagraph.traverse.store import UnknownHandle  # noqa: E402

DOC_NAMES = ('ambiguous_split', 'annotate', 'caps', 'clipped_merge', 'fork', 'limits',
             'linear', 'merge', 'no_bins', 'quorum', 'reentry', 'same_name', 'switch_chain')


def _graphlets():
    """(name, a fresh model) of every CLI document and committed retrieval."""
    out = [(n, (lambda n=n: parse(T.doc_text(n)))) for n in DOC_NAMES]
    for key, r, resp in TG.items_of(TG.golden_files()):
        out.append((key, (lambda r=r, resp=resp: MT.from_response(r, resp))))
    return out


GRAPHLETS = _graphlets()


# ======================================================================= VOP2-03

ORIG_ROOT = 'G * 0 20 0-1 4 * * 0 * TTTCCTATTTAGCCTCTGTC\nP 0 20 4 !\n'
# b.fa (label 1) absent on the seed boundary node, recorded on every node from step 1 on
EDITED_ROOT = 'G * 0 20 0 4 * * 0 * TTTCCTATTTAGCCTCTGTC\nP 0 20 4 0-1\n'


def _cut_reference(g, a, anchor, m):
    """_cut_info()'s labels by definition: the displayed support of every base [0, m) of
    the anchor's chain, from the seed boundary's set (annotate: the root's entry, then
    each P run; constrain: the alive pieces)."""
    if anchor is None:
        return frozenset(a.segments[0].entry)
    inter = None if g.mode == 'constrain' else frozenset(a.segments[0].entry)
    for sid in derive.chain(a, anchor):
        s = a.segments[sid]
        pieces = ops._alive_pieces(a, s) if g.mode == 'constrain' else \
            [(p.from_bp, p.to_bp, frozenset(p.labels)) for p in s.presence]
        for f, t, labs in pieces:
            if f >= m:
                break
            inter = labs if inter is None else inter & labs
    return inter or frozenset()


class TestCutFromTheBoundary(unittest.TestCase):
    """VOP2-03."""

    def setUp(self):
        text = T.doc_text('annotate')
        self.assertIn(ORIG_ROOT, text)
        self.orig = text
        self.edited = text.replace(ORIG_ROOT, EDITED_ROOT)

    def test_the_reviewers_case(self):
        g = parse(self.edited)
        names = {l.id: l.name for l in g.labels}
        for depth in (1, 10, 70):
            seen = set()
            for leaf, m, info in ops._cuts(g, 'right', depth, {}):
                if id(info) in seen:
                    continue
                seen.add(id(info))
                self.assertNotIn('b.fa', {names[l] for l in info[1]}, depth)
                self.assertLessEqual(set(info[2]), set(g.arms['right'].segments[0].entry))
        # the label set at a deeper cut is within the boundary's (monotone in depth)
        a = g.arms['right']
        self.assertLessEqual(ops._cut_info(g, a, 0, 1)[1], ops._cut_info(g, a, None, 0)[1])

    def test_compare_walks_agrees_with_claims_and_labels(self):
        a, b = parse(self.orig), parse(self.edited)
        modes = {m: a.compare(b, mode=m) for m in ('claims', 'labels')}
        self.assertTrue(modes['claims'].only_in_a)
        self.assertTrue(modes['labels'].differ)
        # the walks' keys at a cut inside the boundary label's stretch: b.fa (c:1)
        # supports the original's walk there, not the edited one's (it was on both)
        keys = [ops._walk_keys(g, ['right'], 10, None, {}) for g in (a, b)]
        self.assertEqual([{'c:0', 'c:1'}, {'c:0'}], [set().union(*k.values()) for k in keys])
        pairs = [{x[2] for x in ops._route_pairs(g, ['right'], 10, {})} for g in (a, b)]
        self.assertEqual([{'c:0', 'c:1'}, {'c:0'}], pairs)
        sup = [ops._supported_prefixes(g, ['right'], 70, {}) for g in (a, b)]
        self.assertIn(('right', 'c:1'), sup[0])
        self.assertNotIn(('right', 'c:1'), sup[1])
        # and compare(mode=walks) of retrievals certified to 10 bp says so
        for g in (a, b):
            g.arms['right'].complete_to_bp = 10
        walks = a.compare(b, mode='walks', arm='right')
        self.assertEqual(10, walks.depth_used)
        self.assertEqual(1, len(walks.differ), walks)
        self.assertEqual(['c:0', 'c:1'], [x['ref'] for x in walks.differ[0]['a']])
        self.assertEqual(['c:0'], [x['ref'] for x in walks.differ[0]['b']])

    def test_the_cut_is_the_displayed_support_on_every_fixture(self):
        for name, make in GRAPHLETS:
            g = make()
            for side, a in g.arms.items():
                if not a.segments:
                    continue
                for depth in sorted({1, a.complete_to_bp // 2, a.complete_to_bp}):
                    if depth <= 0:
                        continue
                    for leaf in ops.restricted_leaves(a, depth):
                        m = min(a.segments[leaf].end_bp, depth)
                        k = None
                        if m:
                            k = leaf
                            while a.segments[k].from_bp >= m:
                                k = a.segments[k].parents[0]
                        self.assertEqual(_cut_reference(g, a, k, m),
                                         ops._cut_info(g, a, k, m)[1], (name, side, depth))


# ======================================================================= VOP2-01, 05

def _view_reference(g, labels, mode='any', arm=None):
    """A view's segments and walks as GraphletView built them before (every hit's whole
    first-parent chain)."""
    ids = {l.id for l in labels}
    segments, path_ids = {}, {}
    for side, a in g.arms.items():
        if arm is not None and side != arm:
            continue
        rbs = derive.runs_by_segment(a)
        keep = set()
        for s in a.segments:
            present = set(s.entry) | set(s.end) | {a.runs[r].label for r in rbs[s.id]}
            for p in s.presence:
                present.update(p.labels)
            hit = ids & present
            if (mode == 'any' and hit) or (mode == 'all' and hit == ids):
                keep.update(derive.chain(a, s.id))
        segments[side] = sorted(keep)
        out = []
        for p in derive.paths(a):
            if p.leaf not in keep:
                continue
            have = {e.label for e in derive.end_labels(a, p.leaf)} if g.mode == 'constrain' \
                else set(a.segments[p.leaf].end)
            hit = ids & have
            if (bool(hit) if mode == 'any' else hit == ids):
                out.append(p.id)
        path_ids[side] = out
    return segments, path_ids


class TestViews(unittest.TestCase):
    def test_top_is_cut_after_the_view(self):
        # VOP2-01: the backing arm's top walk is outside the view
        g = parse(T.doc_text('ambiguous_split'))
        v = g.subgraph([{'ref': 'c:3'}])
        self.assertEqual([1], v.path_ids['right'])
        self.assertEqual([1], [w.path_id for w in v.walks('right', top=1)])

    def test_every_view_ranks_then_cuts(self):
        n = 0
        for name, make in GRAPHLETS[:13]:
            g = make()
            for l in g.labels[:6]:
                v = g.subgraph([{'ref': l.ref}])
                for side in v.arms():
                    keep = set(v.path_ids[side])
                    # without top: the backing walks, filtered, object for object
                    self.assertEqual([w for w in g.walks(side) if w.path_id in keep],
                                     v.walks(side), (name, l.ref, side))
                    for by in ('support', 'length', 'id'):
                        ranked = [p for p in ops.rank_walks(g, side, by=by) if p in keep]
                        for k in range(1, len(keep) + 1):
                            self.assertEqual(ranked[:k], [w.path_id for w in v.walks(
                                side, top=k, by=by)], (name, l.ref, side, by, k))
                            n += 1
        self.assertGreater(n, 100)

    def test_a_budgeted_view_list_states_its_stop(self):
        g = parse(comb_annotate(200))
        v = g.subgraph([{'id': 0}])
        want = [w.path_id for w in v.walks('right', by='id')]
        full = MT.LocalBudget()
        v.walks('right', by='id', budget=full)
        resumed = 0
        for w in range(0, full.used_work, max(1, full.used_work // 40)):
            try:
                v.walks('right', by='id', budget=MT.LocalBudget(work_units=w))
                continue
            except MT.LocalBudgetExceeded as e:
                part = e.partial
            if part is None or not part.rows:
                continue                          # stopped while ranking: no rows yet
            rest = v.walks('right', by='id', resume=part.resume)
            self.assertEqual(want, [x.path_id for x in part.rows + rest])
            resumed += 1
            with self.assertRaises(ValueError):
                v.walks('right', resume=part.resume)
        self.assertGreater(resumed, 3)
        for w in range(0, full.used_work, max(1, full.used_work // 20)):
            try:
                v.walks('right', budget=MT.LocalBudget(work_units=w))
            except MT.LocalBudgetExceeded as e:
                self.assertIsNone(e.partial)       # a ranked list is never partial

    def test_the_segments_are_the_chains_of_the_hits(self):
        # VOP2-05: one pass instead of a chain per hit, the same view
        for name, make in GRAPHLETS:
            g = make()
            for l in g.labels[:5]:
                for mode in ('any', 'all'):
                    sel = [l] if mode == 'any' else g.labels[:2]
                    v = g.subgraph([{'ref': x.ref} for x in sel], mode=mode)
                    self.assertEqual(_view_reference(g, sel, mode), (v.segments, v.path_ids),
                                     (name, l.ref, mode))

    def test_a_view_of_a_wide_comb_is_linear(self):
        g = parse(comb_annotate(20000))
        derive.paths(g.arms['right'])
        t0 = time.perf_counter()
        v = g.subgraph([{'id': 1}], of='x')
        took = time.perf_counter() - t0
        self.assertEqual(20001, len(v.segments['right']))     # the root and the spine
        self.assertLess(took, 2.0)

    def test_a_given_of_skips_the_body_digest(self):
        # VOP2-04
        calls = []
        real = parser.dump

        def counting(*a, **k):
            calls.append(1)
            return real(*a, **k)
        g = parse(T.doc_text('fork'))
        parser.dump = counting
        try:
            v = GraphletView(g, [{'id': 0}], of='x')
            self.assertEqual('x', v.of)
            self.assertEqual('x', g.subgraph([{'id': 0}], of='x').of)
            self.assertEqual([], calls)
            self.assertTrue(g.subgraph([{'id': 0}]).of.startswith('sha256:'))
            self.assertEqual(1, len(calls))
        finally:
            parser.dump = real

    def test_subtrie_dumps_the_body_once(self):
        calls = []
        real_p, real_s = parser.dump, store_mod.dump

        def counting(real):
            def f(*a, **k):
                calls.append(1)
                return real(*a, **k)
            return f
        with tempfile.TemporaryDirectory() as d:
            store = GraphletStore(d)
            tools = GraphletTools(store, secret=b'k' * 32)
            g = T.graphlet('fork')
            h = store.put_parsed(g, T.doc_json('fork', 'graphlet')['results'][0]['graphlet'],
                                 {'seeds': []})
            parser.dump, store_mod.dump = counting(real_p), counting(real_s)
            try:
                out = tools.graphlet_subtrie(h, [{'ref': g.labels[0].ref}])
            finally:
                parser.dump, store_mod.dump = real_p, real_s
            self.assertNotIn('error', out, out)
            self.assertEqual(h, out['view']['of'])
            self.assertEqual(1, len(calls))          # the stored body of the new entry

    def test_a_view_needs_a_label(self):
        # VOP2-07
        g = parse(T.doc_text('fork'))
        for sel in (None, []):
            for mode in ('any', 'all'):
                with self.assertRaises(BadSelector):
                    g.subgraph(sel, mode=mode)
                with self.assertRaises(BadSelector):
                    GraphletView(g, sel, mode=mode, budget=MT.LocalBudget())


def _saved_view(d, spec_edit=None):
    """A saved view of fork (label 0), its J line's view edited by |spec_edit|."""
    g = parse(T.doc_text('fork'))
    p = os.path.join(d, 'v.mgt')
    g.subgraph([{'id': 0}]).save(p)
    if spec_edit is not None:
        lines = T.read(p).split('\n')
        j = json.loads(lines[1][2:])
        j['view'] = spec_edit(j['view'])
        lines[1] = 'J ' + json.dumps(j, sort_keys=True, separators=(',', ':'))
        with open(p, 'w', encoding='utf-8', newline='') as f:
            f.write('\n'.join(lines))
    return p


class TestSavedViews(unittest.TestCase):
    """VOP2-07, VOP2-09."""

    def test_one_restoration(self):
        with tempfile.TemporaryDirectory() as d:
            p = _saved_view(d)
            loaded = parser.load(p)
            store = GraphletStore(os.path.join(d, 'sp'))
            h = store.load(p)
            views = [store.view(h), store.view(h)]
            for v in views:
                self.assertIsInstance(v, GraphletView)
                self.assertEqual(loaded.view_spec(), v.view_spec())
                self.assertEqual(loaded.labels, v.labels)
                self.assertEqual((loaded.segments, loaded.path_ids), (v.segments, v.path_ids))
                self.assertEqual(loaded.dump(), v.dump())

    def test_a_saved_view_without_labels_is_refused_on_both_paths(self):
        for edit in (lambda v: dict(v, selectors=[]), lambda v: dict(v, selectors=None),
                     lambda v: {}):
            with tempfile.TemporaryDirectory() as d:
                p = _saved_view(d, edit)
                with self.assertRaises(GraphletFormatError):
                    parser.load(p)
                store = GraphletStore(os.path.join(d, 'sp'))
                h = store.load(p)
                with self.assertRaises(GraphletFormatError):
                    store.view(h)
                tools = GraphletTools(store, secret=b'k' * 32)
                self.assertEqual('format_error', tools.graphlet_summary(h)['error'])


# ======================================================================= VMD-06

class TestReplayOnItsGraph(unittest.TestCase):
    GRAPH = 'sra-logan-chunks/007'

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.store = GraphletStore(self.tmp.name)
        self.client = FakeClient()
        self.client.graph = self.GRAPH
        self.tools = GraphletTools(self.store, {'mini': self.client})

    def tearDown(self):
        self.tmp.cleanup()

    def test_request_of_keeps_the_graph(self):
        resp = copy.deepcopy(T.doc_json('linear', 'graphlet'))
        g = MT.from_response(resp['results'][0], resp)
        self.assertNotIn('graph', GraphletStore.request_of(g))
        self.assertEqual(self.GRAPH, GraphletStore.request_of(g, graph=self.GRAPH)['graph'])
        resp['graph'], resp['graph_path'] = self.GRAPH, '/data/chunk007.dbg'
        g = MT.from_response(resp['results'][0], resp)
        req = GraphletStore.request_of(g)
        self.assertEqual((self.GRAPH, '/data/chunk007.dbg'), (req['graph'], req['graph_path']))

    def test_a_loaded_entry_replays_on_the_clients_graph(self):
        path = os.path.join(self.tmp.name, 'x.mgt')
        T.graphlet('fork').save(path)
        h = self.tools.graphlet_load(path)['handle']
        self.assertNotIn('graph', self.store.get(h).request)
        out = self.tools.traverse_fetch(replay=h, index='mini')
        self.assertNotIn('error', out, out)
        self.assertEqual(self.GRAPH, self.client.requests[-1]['graph'])
        self.assertEqual(self.GRAPH, self.store.get(out['handle']).request['graph'])
        # the stored request itself is not changed in place
        self.assertNotIn('graph', self.store.get(h).request)

    def test_a_continuation_runs_and_is_stored_on_the_graph(self):
        h = self.tools.traverse_fetch(seed={'sequence': SEED['switch_chain']})['handle']
        self.assertEqual(self.GRAPH, self.client.requests[-1]['graph'])
        done = self.tools.traverse_continue(h, 'right', 1)
        self.assertNotIn('error', done, done)
        self.assertEqual(CONTINUATION, self.client.requests[-1]['seeds'][0]['sequence'])
        self.assertEqual(self.GRAPH, self.client.requests[-1]['graph'])
        self.assertEqual(self.GRAPH, self.store.get(done['handle']).request['graph'])
        # a graph the overrides name wins
        self.tools.traverse_continue(h, 'right', 1, overrides={'graph': 'other'})
        self.assertEqual('other', self.client.requests[-1]['graph'])


# ======================================================================= VMD-01

class Clock:
    t = 1_000_000.0

    def __call__(self):
        return self.t


class TestSteadyUseSurvivesARestart(unittest.TestCase):
    TTL = 7 * 86400

    def test_a_handle_used_every_twentieth_of_the_ttl(self):
        with tempfile.TemporaryDirectory() as d:
            clock = Clock()
            store = GraphletStore(d, clock=clock, ttl_disk_s=self.TTL)
            h = store.load(_doc_file(d, 'linear'))
            for _ in range(30):                        # 1.5 TTL of use, each gap TTL / 20
                clock.t += self.TTL / 20
                store.get(h)
            restarted = GraphletStore(d, clock=clock, ttl_disk_s=self.TTL)
            self.assertEqual([], restarted.sweep()['expired'])
            restarted.get(h)
            # the gaps just under a tenth of the TTL, used across a restart
            for _ in range(15):
                clock.t += self.TTL / 10 - 1
                restarted.get(h)
            self.assertEqual([], GraphletStore(d, clock=clock,
                                               ttl_disk_s=self.TTL).sweep()['expired'])

    def test_freed_and_expired_handles_are_forgotten(self):
        with tempfile.TemporaryDirectory() as d:
            clock = Clock()
            store = GraphletStore(d, clock=clock, ttl_disk_s=100)
            a, b = store.load(_doc_file(d, 'linear')), store.load(_doc_file(d, 'fork'))
            store.free(a)
            self.assertNotIn(a, store._persisted)
            clock.t += 1000
            self.assertEqual([b], store.sweep()['expired'])
            self.assertNotIn(b, store._persisted)


def _doc_file(d, name):
    p = os.path.join(d, name + '.mgt')
    if not os.path.exists(p):
        T.graphlet(name).save(p)
    return p


# ======================================================================= VPC-04, 05

class TestCodecAndJLine(unittest.TestCase):
    def test_every_negative_id_is_refused(self):
        for ids in ([-5, 3], [-5, -4, 7], [-1], [-3, -2], [0, 1 << 64]):
            with self.assertRaises(CodecError, msg=ids):
                encode_ranges(ids)
        rng = random.Random(7)
        for _ in range(500):
            ids = sorted(rng.sample(range(200), rng.randint(0, 40)))
            self.assertEqual(ids, list(decode_ranges(encode_ranges(ids))))

    def cases(self):
        resp = T.doc_json('linear', 'graphlet')
        body = resp['results'][0]['graphlet']
        counts = {k: resp['results'][0][k] for k in ('graphlet_bytes', 'graphlet_lines')}
        return [('real', resp['results'][0], resp),
                ('empty summary, empty envelope', {'graphlet': body}, {'results': [None]}),
                ('empty summary', {'graphlet': body}, resp),
                ('transport counts only', dict(counts, graphlet=body), resp),
                ('empty envelope', resp['results'][0], {'results': resp['results']})]

    def test_a_j_line_reads_back_as_written(self):
        # VPC-05
        for name, result, response in self.cases():
            with self.subTest(name):
                g = MT.from_response(result, response)
                text = parser.dump(g, envelope=True)
                self.assertTrue(parser.is_canonical(text))
                self.assertEqual(text, parser.dump(parse(text)))
                with tempfile.TemporaryDirectory() as d:
                    p1, p2 = os.path.join(d, 'a.mgt'), os.path.join(d, 'b.mgt')
                    parser.save(g, p1)
                    back = parser.load(p1)
                    self.assertEqual(g.has_envelope, back.has_envelope)
                    parser.save(back, p2)
                    self.assertEqual(T.read(p1), T.read(p2))
                st = parser.standalone_text(result['graphlet'], result, response)
                self.assertEqual(text, st)
                self.assertTrue(parser.is_canonical(st))

    def test_the_store_keeps_an_envelope_with_one_empty_side(self):
        for name, result, response in self.cases()[2:]:
            with self.subTest(name), tempfile.TemporaryDirectory() as d:
                g = MT.from_response(result, response)
                store = GraphletStore(d, max_ram_mb=0)
                h = store.put_parsed(g, result['graphlet'], {'seeds': []})
                self.assertTrue(store.graphlet(h).has_envelope)
                self.assertEqual(parser.dump(g, envelope=True), store.standalone_text(h))


# ======================================================================= VPC-01..03

class TestTransportChecks(unittest.TestCase):
    def setUp(self):
        self.resp = T.doc_json('linear', 'graphlet')
        self.result = self.resp['results'][0]
        self.body = self.result['graphlet']

    def both(self, result):
        """(from_response's error, standalone_text's error): (type, message)."""
        out = []
        for fn in (lambda: MT.from_response(result, self.resp),
                   lambda: parser.standalone_text(result['graphlet'], result, self.resp)):
            try:
                fn()
                out.append(None)
            except Exception as e:
                out.append((type(e), str(e)))
        return out

    def test_a_count_that_is_no_integer(self):
        # VPC-01
        for key in ('graphlet_bytes', 'graphlet_lines'):
            for bad in (None, '12', 12.5, True, [1]):
                got = self.both(dict(self.result, **{key: bad}))
                for e in got:
                    self.assertIsNotNone(e, (key, bad))
                    self.assertIs(GraphletFormatError, e[0], (key, bad, e))

    def test_the_same_checks_in_the_same_order(self):
        # VPC-02: a body with a broken record and a wrong line count is the count's error
        # on both paths (from_response checked it after the parse)
        lines = self.body.split('\n')
        broken = '\n'.join(lines[:3] + ['Q nonsense'] + lines[3:])
        n = self.result['graphlet_lines']
        for result in (dict(self.result, graphlet_bytes=self.result['graphlet_bytes'] + 1),
                       dict(self.result, graphlet_lines=n + 1),
                       dict(self.result, graphlet=broken, graphlet_lines=n + 7,
                            graphlet_bytes=len(broken.encode()))):
            a, b = self.both(result)
            self.assertEqual(a, b)
            self.assertIs(GraphletFormatError, a[0])

    def test_no_copy_of_the_body_for_its_j_test(self):
        # VPC-03: the splice holds the body, its middle and the text (3.0x the body before)
        resp = T.doc_json('annotate', 'graphlet')
        r = {k: v for k, v in resp['results'][0].items()
             if k not in ('graphlet_bytes', 'graphlet_lines')}
        r['graphlet'] = body = comb_annotate(20000)
        gc.collect()
        tracemalloc.start()
        base = tracemalloc.get_traced_memory()[0]
        try:
            out = parser.standalone_text(body, r, resp)
            peak = tracemalloc.get_traced_memory()[1] - base
        finally:
            tracemalloc.stop()
        self.assertTrue(out.startswith(body[:body.index('\n') + 1] + 'J '))
        self.assertLess(peak, 2.3 * len(body))


# ======================================================================= VMT-*, VOP1-*

class ToolsCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.store = GraphletStore(self.tmp.name)
        self.client = FakeClient()
        self.tools = GraphletTools(self.store, {'mini': self.client})

    def tearDown(self):
        self.tmp.cleanup()

    def fetch(self, name):
        out = self.tools.traverse_fetch(seed={'sequence': SEED[name]})
        self.assertNotIn('error', out, out)
        return out['handle']


class TestToolArguments(ToolsCase):
    def test_a_segment_takes_no_seed(self):
        # VMT-04
        h = self.fetch('fork')
        out = self.tools.graphlet_sequence(h, 'right', segment=0, with_seed=True)
        self.assertEqual('bad_argument', out['error'])
        out = self.tools.graphlet_sequence(h, 'right', segment=0)
        self.assertFalse(out['with_seed'])
        full = self.tools.graphlet_sequence(h, 'right', walk=0, with_seed=True)
        self.assertTrue(full['with_seed'])

    def test_continue_arguments_are_no_overrides(self):
        # VMT-01 (fixed at f667d775)
        h = self.fetch('switch_chain')
        for ov in ({'reset_branches': True}, {'arm': 'left'}, {'leaves': [0]}):
            out = self.tools.traverse_continue(h, 'right', 1, overrides=ov, execute=False)
            self.assertEqual('bad_argument', out['error'], ov)

    def test_an_infinite_branch_allowance(self):
        # VOP1-01 (fixed at dcc0cebd)
        h = self.fetch('switch_chain')
        g = self.store.graphlet(h)
        for v in (math.inf, -math.inf, math.nan):
            with self.assertRaises(ValueError) as e:
                g.next_request('right', [1], branching={'max_label_branches': v})
            self.assertIn('branching.max_label_branches', str(e.exception))
            out = self.tools.traverse_continue(h, 'right', 1, execute=False,
                                               overrides={'branching':
                                                          {'max_label_branches': v}})
            self.assertEqual('bad_argument', out['error'])

    def test_a_bool_is_no_label_id(self):
        # VOP1-02
        g = parse(T.doc_text('fork'))
        for v in (True, False):
            with self.assertRaises(BadSelector):
                ops.label(g, {'id': v})
            with self.assertRaises(BadSelector):
                ops.labels_matching(g, [{'id': v}])
        self.assertEqual(g.labels[1], ops.label(g, {'id': 1}))
        h = self.fetch('fork')
        out = self.tools.graphlet_claims(h, 'right', labels={'id': True})
        self.assertEqual('bad_argument', out['error'])
        self.assertNotIn('error', self.tools.graphlet_claims(h, 'right', labels={'id': 1}))


class TestPaging(ToolsCase):
    """VMT-02."""

    def test_a_large_page_is_linear(self):
        resp = T.doc_json('annotate', 'graphlet')
        r = {k: v for k, v in resp['results'][0].items()
             if k not in ('graphlet_bytes', 'graphlet_lines')}
        r['graphlet'] = comb_annotate(4000)
        g = MT.from_response(r, resp)
        h = self.store.put_parsed(g, r['graphlet'], {'seeds': []})
        t0 = time.perf_counter()
        out = self.tools.graphlet_claims(h, 'right', max_bytes=1 << 20)
        took = time.perf_counter() - t0
        self.assertGreater(len(out['rows']), 3000)
        self.assertLessEqual(_size(out), 1 << 20)
        self.assertLess(took, 5.0)          # minutes when each row re-measured the page
        # the page is what one measured as a whole: one row more would not fit
        nxt = self.tools.graphlet_claims(h, 'right', max_bytes=1 << 20,
                                         cursor=out.get('next_cursor'))
        if out.get('next_cursor'):
            more = dict(out, rows=out['rows'] + nxt['rows'][:1])
            self.assertGreater(_size(more), 1 << 20)

    def test_every_page_is_the_measured_one(self):
        h = self.fetch('annotate')
        for mb in (64, 200, 512, 1000, 4096):
            cur = None
            while True:
                out = self.tools.graphlet_claims(h, 'right', max_bytes=mb, cursor=cur)
                self.assertLessEqual(_size(out), mb)
                cur = out.get('next_cursor')
                if not cur or 'error' in out:
                    break
                if out['rows'] and not out.get('row_truncated'):
                    nxt = self.tools.graphlet_claims(h, 'right', max_bytes=mb, cursor=cur)
                    if nxt.get('rows'):
                        # the next row did not fit this page
                        self.assertGreater(_size(dict(out, rows=out['rows']
                                                      + nxt['rows'][:1])), mb)

    def test_max_bytes_has_a_ceiling(self):
        h = self.fetch('fork')
        out = self.tools.graphlet_claims(h, 'right', max_bytes=MAX_MAX_BYTES + 1)
        self.assertEqual('bad_argument', out['error'])
        self.assertNotIn('error', self.tools.graphlet_claims(h, 'right',
                                                             max_bytes=MAX_MAX_BYTES))

    def test_resolve_is_trimmed_in_linear_time(self):
        class Wide(FakeClient):
            def resolve(self, sequence, **kw):
                return {'k': 31, 'num_kmers': 9,
                        'labels': [{'label': 'sample_%06d_%s' % (i, 'x' * (i % 40)),
                                    'kind': 'column', 'kmers_supported': i}
                                   for i in range(6000)],
                        'candidates': [{}],
                        'selection': {'seeds': [{'seed_id': 's%d' % i, 'sequence': 'ACGT' * 3,
                                                 'labels': ['a'] * (i % 5)}
                                                for i in range(3000)]}}
        tools = GraphletTools(self.store, {'mini': Wide()})
        # the rows kept are those the old trimming (a row off, then the whole answer
        # measured again) kept, where that is affordable
        for limit, mb in ((300, None), (300, 20_000), (10, None), (200, 2_000), (300, 200)):
            out = tools.traverse_resolve(sequence='ACGT', limit=limit, max_bytes=mb)
            self.assertEqual(_resolve_reference(Wide().resolve('ACGT'), limit,
                                                mb or tools.max_bytes), out, (limit, mb))
        # and where it is not (46 s at limit 5,000): fast, within the ceiling, and no row
        # more would fit -- the seeds' last ones off first, then the labels'
        for limit, mb in ((5000, None), (5000, 100_000), (3000, 300_000), (3000, 150_000)):
            t0 = time.perf_counter()
            out = tools.traverse_resolve(sequence='ACGT', limit=limit, max_bytes=mb)
            self.assertLess(time.perf_counter() - t0, 2.0)
            lim = mb or tools.max_bytes
            self.assertLessEqual(_size(out), lim)
            self.assertTrue(out['labels_cut_to_fit'])
            full = Wide().resolve('ACGT')
            if len(out['labels']) < limit:
                self.assertEqual([], out['seeds'])
                more = dict(out, labels=out['labels'] + [
                    {'name': full['labels'][len(out['labels'])]['label'], 'kind': 'column',
                     'kmers_supported': len(out['labels'])}])
            else:
                s = full['selection']['seeds'][len(out['seeds'])]
                more = dict(out, seeds=out['seeds'] + [
                    {'seed_id': s['seed_id'], 'sequence_bp': 12, 'labels': len(s['labels'])}])
            self.assertGreater(_size(more), lim)


def _resolve_reference(out, limit, lim):
    """traverse_resolve's answer as it was trimmed before: a row off, then the whole
    answer measured again (quadratic: for small limits only)."""
    labels_ = out.get('labels', [])
    res = {'index': 'mini', 'k': out.get('k'), 'num_kmers': out.get('num_kmers'),
           'labels_total': len(labels_),
           'labels': [{'name': l['label'], 'kind': l.get('kind'),
                       'kmers_supported': l.get('kmers_supported')} for l in labels_[:limit]],
           'candidates': len(out.get('candidates', []))}
    if 'selection' in out:
        res['seeds'] = [{'seed_id': s['seed_id'], 'sequence_bp': len(s['sequence']),
                         'labels': len(s['labels'])}
                        for s in out['selection']['seeds'][:limit]]
    while _size(res) > lim and (res['labels'] or res.get('seeds')):
        res['labels_cut_to_fit'] = True
        if res.get('seeds'):
            res['seeds'].pop()
        else:
            res['labels'].pop()
    return res


# ======================================================================= VOP1-03, 04

class TestWalksAndDeadCode(unittest.TestCase):
    def test_walk_claims_are_their_leaves_claims(self):
        # VOP1-03 (fixed at f667d775): the claims at each returned leaf, built for it
        for name, make in GRAPHLETS:
            g = make()
            for side in g.arms:
                every = g.claims(side, strict=False)
                for w in g.walks(side, top=5):
                    end = g.arms[side].segments[w.leaf].end_bp
                    self.assertEqual([c for c in every if c.segment == w.leaf
                                      and c.to_bp == end], w.claims, (name, side, w.path_id))

    def test_one_switch_pricing(self):
        # VOP1-04
        self.assertFalse(hasattr(ops, '_switch_cost'))


# ======================================================================= VOP2-02, 08

def _one_arm(name, side):
    g = T.graphlet(name)
    for s in list(g.arms):
        if s != side:
            del g.arms[s]
    return g


class TestCompareEarlyAnswers(unittest.TestCase):
    def test_an_arm_b_lacks(self):
        # VOP2-02
        a, b = T.graphlet('fork'), _one_arm('fork', 'right')
        for mode in ('claims', 'walks', 'labels', 'prefix_subset'):
            c = a.compare(b, arm='left', mode=mode)
            self.assertFalse(c.comparable)
            self.assertIn('left arm was not retrieved in b', c.reason)
            self.assertEqual((a.support, b.support), c.support)
            self.assertTrue(ops.compare_cost(a, b, arm='left', mode=mode)['exact'])
            self.assertTrue(a.compare(b, arm='right', mode=mode).comparable)
        with tempfile.TemporaryDirectory() as d:
            store = GraphletStore(d)
            tools = GraphletTools(store, secret=b'k' * 32, export_dir=d)
            body = T.doc_json('fork', 'graphlet')['results'][0]['graphlet']
            ha = store.put_parsed(a, body, {'seeds': []})
            hb = store.put_parsed(b, body, {'seeds': []})
            out = tools.graphlet_compare(ha, hb, arm='left')
            self.assertNotIn('error', out, out)
            self.assertFalse(out['comparable'])
            out = tools.graphlet_export(ha, format='json', what='compare', other=hb,
                                        arm='left', path='c.json')
            self.assertNotIn('error', out, out)

    def test_no_arm_in_common_states_both_supports(self):
        # VOP2-08
        a, b = _one_arm('fork', 'left'), _one_arm('fork', 'right')
        b.support = 'trace'
        c = a.compare(b)
        self.assertEqual((False, 'no arm in common'), (c.comparable, c.reason))
        self.assertEqual((a.support, 'trace'), c.support)
        self.assertEqual([a.support, 'trace'], list(c.as_dict()['support']))

    def test_the_mode_is_checked_first(self):
        fork, other = T.graphlet('fork'), T.graphlet('fork')
        other.index_meta_fp = 'f' * 16                               # another index
        self.assertIs(False, fork.compare(other).comparable)
        for a, b in ((fork, other), (_one_arm('fork', 'left'), _one_arm('fork', 'right')),
                     (fork, T.graphlet('fork'))):
            with self.assertRaises(ValueError):
                a.compare(b, mode='bogus')
            with self.assertRaises(ValueError):
                a.compare(b, mode='bogus', budget=MT.LocalBudget())


# ======================================================================= VMD-02..04

class TestSharedSpool(unittest.TestCase):
    """VMD-02: two stores (processes) on one spool."""

    def test_no_body_is_deleted_under_another_entry(self):
        with tempfile.TemporaryDirectory() as d:
            p = _doc_file(d, 'fork')
            a = GraphletStore(os.path.join(d, 'spool'))
            h1 = a.load(p)
            b = GraphletStore(os.path.join(d, 'spool'), max_ram_mb=0)
            h2 = b.load(p)                     # the same body, deduplicated
            self.assertEqual(a.get(h1).digest, b.get(h2).digest)
            a.free(h1)                         # a does not know h2
            self.assertEqual(MT.dump(T.graphlet('fork')), MT.dump(b.graphlet(h2)))

    def test_entries_of_the_other_store(self):
        with tempfile.TemporaryDirectory() as d:
            p = _doc_file(d, 'fork')
            clock = Clock()
            a = GraphletStore(os.path.join(d, 'spool'), clock=clock, ttl_disk_s=1000)
            b = GraphletStore(os.path.join(d, 'spool'), clock=clock, ttl_disk_s=1000)
            h = a.load(p)
            self.assertEqual(MT.dump(a.graphlet(h), envelope=True),
                             MT.dump(b.graphlet(h), envelope=True))   # adopted
            # a keeps using it, b never does: b's sweep does not expire it
            for _ in range(5):
                clock.t += 400
                a.get(h)
            self.assertEqual([], b.sweep()['expired'])
            self.assertIn(h, b)
            a.free(h)
            with self.assertRaises(UnknownHandle):
                b.get(h)
            self.assertNotIn(h, b)
            self.assertEqual([], b.list())

    def test_free_in_and_list_before_any_get(self):
        # the review of these fixes, finding 5: free() and `in` read this process's index
        # only -- free of a live handle another process stored answered unknown_handle and
        # left its entry and body on disk, though a summary of it answered
        with tempfile.TemporaryDirectory() as d:
            p = _doc_file(d, 'fork')
            spool = os.path.join(d, 'spool')
            a = GraphletStore(spool)
            b = GraphletStore(spool)
            h = a.load(p)
            digest = a.get(h).digest
            self.assertIn(h, b)
            self.assertEqual([h], [e['handle'] for e in b.list()])
            b.free(h)
            self.assertFalse(os.path.exists(a._entry_path(h)))
            self.assertFalse(os.path.exists(a._body_path(digest)))
            self.assertNotIn(h, a)
            self.assertEqual([], a.list())
            with self.assertRaises(UnknownHandle):
                a.get(h)
            # through the tools of the second store, before it used the handle
            h = a.load(p)
            tools = GraphletTools(GraphletStore(spool), secret=b'k' * 32)
            self.assertEqual({'freed': h}, tools.graphlet_free(h))
            self.assertEqual('unknown_handle', tools.graphlet_summary(h)['error'])
            self.assertEqual('unknown_handle',
                             GraphletTools(a, secret=b'k' * 32).graphlet_summary(h)['error'])
            # a handle no process stored is still unknown everywhere
            with self.assertRaises(UnknownHandle):
                b.free('g_' + '0' * 12)
            self.assertNotIn('g_' + '0' * 12, b)


class TestStoreFiles(unittest.TestCase):
    def test_entry_files_of_other_versions(self):
        # VMD-03
        with tempfile.TemporaryDirectory() as d:
            store = GraphletStore(d)
            h = store.load(_doc_file(d, 'fork'))
            entries = os.path.join(d, 'entries')
            with open(os.path.join(entries, h + '.json')) as f:
                j = json.load(f)
            newer = dict(j, handle='g_' + 'a' * 12, future_field={'x': 1})
            with open(os.path.join(entries, newer['handle'] + '.json'), 'w') as f:
                json.dump(newer, f)
            for name, content in (('g_' + 'b' * 12, [1, 2]),
                                  ('g_' + 'c' * 12, {'handle': 'g_' + 'c' * 12}),
                                  ('g_' + 'd' * 12, None)):
                with open(os.path.join(entries, name + '.json'), 'w') as f:
                    json.dump(content, f)
            reopened = GraphletStore(d)
            self.assertEqual(sorted([h, newer['handle']]),
                             sorted(e['handle'] for e in reopened.list()))
            self.assertEqual(MT.dump(reopened.graphlet(h)),
                             MT.dump(reopened.graphlet(newer['handle'])))
            # a tombstone that is no object
            with open(os.path.join(d, 'tombstones', 'g_' + 'e' * 12 + '.json'), 'w') as f:
                json.dump(['x'], f)
            with self.assertRaises(UnknownHandle) as e:
                reopened.get('g_' + 'e' * 12)
            self.assertFalse(e.exception.replayable)

    def test_a_missing_body_expires_its_entry(self):
        # VMD-04
        with tempfile.TemporaryDirectory() as d:
            store = GraphletStore(os.path.join(d, 'spool'), max_ram_mb=0)
            tools = GraphletTools(store, secret=b'k' * 32)
            with_req = store.load(_doc_file(d, 'fork'))
            self.assertTrue(store.get(with_req).request)
            os.unlink(store._body_path(store.get(with_req).digest))
            with self.assertRaises(UnknownHandle) as e:
                store.graphlet(with_req)
            self.assertTrue(e.exception.replayable)
            self.assertNotIn(with_req, store)
            out = tools.graphlet_summary(with_req)
            self.assertEqual('unknown_handle', out['error'])
            self.assertTrue(out['replayable'])
            self.assertIn('replay', out['hint'])
            # an entry without a request: not offered for replay
            g = T.graphlet('linear')
            bare = store.put_graphlet(g, None)
            os.unlink(store._body_path(store.get(bare).digest))
            with self.assertRaises(UnknownHandle) as e:
                store.graphlet(bare)
            self.assertFalse(e.exception.replayable)

    def test_a_missing_body_expires_its_entry_on_a_copy(self):
        # the review of these fixes, finding 4: the copy without a parse (graphlet_save
        # and graphlet_export(format=mgt) under local limits; body_text, standalone_text
        # and save_body) answered io_error and left the entry listed
        with tempfile.TemporaryDirectory() as d:
            calls = [
                ('save', lambda t, s, h: t.graphlet_save(h, 'x.mgt')),
                ('export_mgt', lambda t, s, h: t.graphlet_export(h, format='mgt',
                                                                 path='y.mgt')),
                ('body_text', lambda t, s, h: s.body_text(h)),
                ('standalone_text', lambda t, s, h: s.standalone_text(
                    h, budget=MT.LocalBudget())),
                ('save_body', lambda t, s, h: s.save_body(h, os.path.join(d, 'z.mgt'))),
            ]
            for name, call in calls:
                store = GraphletStore(os.path.join(d, name), max_ram_mb=0)
                tools = GraphletTools(store, secret=b'k' * 32, export_dir=d,
                                      local_limits=ToolLimits())
                h = store.load(_doc_file(d, 'fork'))
                os.unlink(store._body_path(store.get(h).digest))
                try:
                    out = call(tools, store, h)
                except UnknownHandle as e:
                    self.assertTrue(e.replayable, name)
                else:
                    self.assertEqual('unknown_handle', out['error'], (name, out))
                    self.assertTrue(out['replayable'], name)
                self.assertNotIn(h, store, name)
                self.assertEqual([], store.list(), name)
                self.assertTrue(os.path.exists(store._tomb_path(h)), name)


# ======================================================================= VMD-07, 08

class _Answer:
    def __init__(self, status, content, headers=None, raw_encoded=False):
        self.status_code = status
        self.content = content
        self.headers = headers or {}
        self.raw_encoded = raw_encoded


class _Session:
    def __init__(self, answer):
        self.answer = answer

    def request(self, method, url, data=None, headers=None, timeout=None):
        return self.answer


class TestUndecodableBodies(unittest.TestCase):
    """VMD-07."""

    def post(self, answer):
        return TraverseClient('h', 1, session=_Session(answer)).traverse_raw({'seeds': []})

    def test_bodies(self):
        latin1 = 'Serveur indisponible: d\xe9marrage'.encode('latin-1')
        with self.assertRaises(ServerInitializing) as e:
            self.post(_Answer(503, latin1, {'Retry-After': '30'}))
        self.assertEqual(503, e.exception.status)
        with self.assertRaises(TraverseError) as e:
            self.post(_Answer(502, latin1))
        self.assertEqual(502, e.exception.status)
        self.assertNotIsInstance(e.exception, ServerInitializing)
        self.assertIn('Serveur', e.exception.message)
        with self.assertRaises(TraverseError) as e:
            self.post(_Answer(200, b'\xff\xfe{}'))
        self.assertEqual('the response is not JSON', e.exception.message)
        cut = gzip.compress(json.dumps({'results': []}).encode())[:-8]
        for status in (200, 500):
            with self.assertRaises(TraverseError) as e:
                self.post(_Answer(status, cut, {'Content-Encoding': 'gzip'}, raw_encoded=True))
            self.assertEqual(status, e.exception.status)
            self.assertIn('could not be decoded', e.exception.message)
        self.assertEqual({'results': []}, self.post(_Answer(
            200, gzip.compress(b'{"results": []}'), {'Content-Encoding': 'gzip'},
            raw_encoded=True)))


class TestOneStatementOfLabelChanges(unittest.TestCase):
    def test_the_end_rule_and_the_displayed_support_read_one_list(self):
        # VMD-08
        n = 0
        for name, make in GRAPHLETS:
            g = make()
            if g.mode != 'constrain':
                continue
            for side, a in g.arms.items():
                rbs = derive.runs_by_segment(a)
                for s in a.segments:
                    changes = derive.label_changes(s, rbs[s.id], a.runs)
                    self.assertEqual(changes, ops._segment_ops(a, s))
                    cur = set(s.entry)
                    for _, kind, lab in changes:
                        if kind == 0:
                            cur.discard(lab)
                        else:
                            cur.add(lab)
                    self.assertEqual(sorted(cur), sorted(derive.rule_end(
                        g.mode, s, s.entry, s.presence, rbs[s.id], a.runs)))
                    self.assertEqual(sorted(cur), sorted(s.end), (name, side, s.id))
                    n += 1
        self.assertGreater(n, 100)


if __name__ == '__main__':
    unittest.main()
