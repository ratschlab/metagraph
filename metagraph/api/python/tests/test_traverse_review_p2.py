"""The P2 items of the review of 2026-10-06 in the library's tools, store and exports:

  * O1: a replay or continuation whose parse stops on its local budget is checked against
    its entry's index like a parsed one -- the identity its H record states, read without a
    parse -- before its entry and parent link are stored; and the capabilities read before
    a continuation are those of the graph it is sent to (overrides={'graph': ...});
  * O2: store.standalone_text() builds its text in one join of slices of the body (about
    3 bodies held, the account's 3), and charges a call's base (CALL_BASE);
  * O3: to_json(budget=) charges the seed block (dropped labels and their runs), the
    seed-level and arm-level limitations and label_dict before building them;
  * X1 (library part): the derivation limitation states k-mers read on a walked (D3)
    result and elapsed ms on a failed seed; the library reads the first as k-mers and
    passes the second on untouched.

And the review of those fixes:

  * next_request() with a change_cost table: the copies of its entries, the switch
    search's table, edges and heap are in the account and the work (work model 2), so the
    account bounds the traced peak (it stayed at the 0-entry figure, 9.8x short at 900
    entries);
  * the texts: the derivation row states D3 at feature level 6, and no text says an
    attempt_id is "used up" -- the server refuses it while it holds it and runs it after.
"""

import copy
import gc
import os
import sys
import tempfile
import tracemalloc
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402
import test_traverse_mcp as M  # noqa: E402

from metagraph.traverse import (  # noqa: E402
    GraphletStore, LocalBudget, LocalBudgetExceeded, LocalLimits, from_response,
)
from metagraph.traverse import budget as B  # noqa: E402
from metagraph.traverse.mcp_tools import GraphletTools, ToolLimits  # noqa: E402

FP_A, FP_B = 'a' * 64, 'b' * 64


def _traced(fn):
    gc.collect()
    tracemalloc.start()
    base = tracemalloc.get_traced_memory()[0]
    try:
        out = fn()
        return tracemalloc.get_traced_memory()[1] - base, out
    finally:
        tracemalloc.stop()


# ======================================================================= O1

class _Stamped(M.FakeClient):
    """The fake backend answering with index_fp |stamp| in H (what answered), while its
    capabilities state |caps_fp| per graph (what was asked before)."""

    def __init__(self, caps_fp=None, graph=None):
        super().__init__()
        self.stamp = FP_A
        self.caps_fp = caps_fp or {None: FP_A}
        self.graph = graph
        self.caps_asked = []

    def traverse_raw(self, request):
        self.fp = self.stamp
        return super().traverse_raw(request)

    def capabilities(self, graph=None, graph_path=None):
        self.caps_asked.append(graph)
        return {'release': '', 'index_fp': self.caps_fp.get(graph or self.graph)}


class TestUnparsedIdentity(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.store = GraphletStore(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def tools(self, client):
        return GraphletTools(self.store, {'mini': client}, local_limits=ToolLimits())

    def parent(self, t):
        out = t.traverse_fetch(seed={'sequence': M.SEED['switch_chain']}, max_bytes=8192)
        return out['handle']

    def stopping_parse(self, t):
        # the continuation's parse stops on this (the review's work_units 500)
        t.local_limits = ToolLimits(parse=LocalLimits(500, 1024))

    def entries(self):
        return sorted(m['handle'] for m in self.store.list())

    def test_a_continuation_answered_from_another_index(self):
        c = _Stamped()
        t = self.tools(c)
        h = self.parent(t)
        before = self.entries()
        c.stamp = FP_B                 # the capabilities still say A: only H tells
        self.stopping_parse(t)
        out = t.traverse_continue(h, 'right', 1, max_bytes=8192)
        self.assertEqual('index_mismatch', out.get('error'), out)
        self.assertEqual(before, self.entries())      # no child, no parent link

    def test_a_replay_answered_from_another_index(self):
        c = _Stamped()
        t = self.tools(c)
        h = self.parent(t)
        before = self.entries()
        c.stamp = FP_B
        self.stopping_parse(t)
        out = t.traverse_fetch(replay=h, max_bytes=8192)
        self.assertEqual('index_mismatch', out.get('error'), out)
        self.assertEqual(before, self.entries())

    def test_an_unparsed_answer_without_a_manifest_digest(self):
        # the parent states fp A; the answer's H none: not proven the same index
        c = _Stamped()
        t = self.tools(c)
        h = self.parent(t)
        c.stamp = None
        self.stopping_parse(t)
        out = t.traverse_continue(h, 'right', 1, max_bytes=8192)
        self.assertEqual('index_unverifiable', out.get('error'), out)
        out = t.traverse_continue(h, 'right', 1, max_bytes=8192, allow_unverified_index=True)
        self.assertNotIn('error', out)
        self.assertIs(False, out['parsed'])
        self.assertEqual((False, 'allow_unverified_index'),
                         (out['identity']['verified'], out['identity']['accepted']))
        self.assertEqual(h, self.store.get(out['handle']).parent['handle'])

    def test_the_same_index_unparsed_is_stored_as_before(self):
        c = _Stamped()
        t = self.tools(c)
        h = self.parent(t)
        self.stopping_parse(t)
        out = t.traverse_continue(h, 'right', 1, max_bytes=8192)
        self.assertNotIn('error', out)
        self.assertEqual({'verified': True, 'index_fp': FP_A}, out['identity'])
        self.assertEqual(h, self.store.get(out['handle']).parent['handle'])

    def test_the_capabilities_of_the_graph_the_continuation_is_sent_to(self):
        c = _Stamped({'g1': FP_A, 'g2': FP_B}, graph='g1')
        t = self.tools(c)
        h = self.parent(t)
        n = len(c.requests)
        c.caps_asked.clear()
        for limits in (None, ToolLimits(parse=LocalLimits(500, 1024))):
            t.local_limits = limits
            out = t.traverse_continue(h, 'right', 1, overrides={'graph': 'g2'},
                                      max_bytes=8192)
            self.assertEqual('index_mismatch', out.get('error'), out)
        self.assertEqual(['g2', 'g2'], c.caps_asked)
        self.assertEqual(n, len(c.requests))          # refused before it was sent
        # on the client's own graph the capabilities are asked as before (no argument)
        c.caps_asked.clear()
        t.local_limits = ToolLimits()
        self.assertNotIn('error', t.traverse_continue(h, 'right', 1, max_bytes=8192))
        self.assertEqual([None], c.caps_asked)


# ======================================================================= O2

def _items():
    import glob
    import traverse_golden as TG
    files = sorted(glob.glob(os.path.join(T.DOCS, '**', '*.graphlet.json'), recursive=True))
    return TG.items_of(files + TG.golden_files())


class TestStandaloneText(unittest.TestCase):
    """O2: the account of standalone_text() bounds its traced peak, on every fixture, with
    and without a delivery written into the O record (the review: 26 of 90 calls peaked
    above their account, mini_ndm1__trace by 30%)."""

    def test_the_account_bounds_the_peak(self):
        n = 0
        with tempfile.TemporaryDirectory() as d:
            st = GraphletStore(os.path.join(d, 'sp'))
            for key, r, resp in _items():
                body = r.get('graphlet')
                if not body:
                    continue
                for delivery in (None, 'spooled'):
                    h = st.put_parsed(from_response(r, resp), body, {'seeds': []},
                                      resident=False, delivery=delivery)
                    b = LocalBudget()
                    want = st.standalone_text(h, budget=b)
                    peak, got = _traced(lambda: st.standalone_text(h))
                    self.assertEqual(want, got)
                    self.assertGreaterEqual(b.usage()['memory_bytes'], peak,
                                            (key, delivery, st.get(h).bytes))
                    # a call's base is part of it (a direct library call opened no scope)
                    self.assertGreaterEqual(b.usage()['memory_bytes'], B.CALL_BASE)
                    n += 1
        self.assertGreater(n, 80)

    def test_the_text_is_the_saved_graphlet(self):
        with tempfile.TemporaryDirectory() as d:
            st = GraphletStore(os.path.join(d, 'sp'))
            for key, r, resp in _items()[:20]:
                body = r.get('graphlet')
                if not body:
                    continue
                for delivery in (None, 'spooled', 'paged'):
                    g = from_response(r, resp)
                    h = st.put_parsed(g, body, {'seeds': []}, resident=False,
                                      delivery=delivery)
                    path = os.path.join(d, 'x.mgt')
                    g2 = from_response(r, resp)
                    if delivery is not None:
                        from metagraph.traverse.store import set_delivery
                        set_delivery(g2, delivery)
                    g2.save(path)
                    with open(path, encoding='utf-8', newline='') as f:
                        self.assertEqual(f.read(), st.standalone_text(h), (key, delivery))

    def test_a_limit_between_the_old_peak_and_the_account(self):
        # the review's library case: a limit equal to the call's account completes, and
        # the traced peak stays under it
        key, r, resp = next(x for x in _items() if 'mini_ndm1__trace' in x[0])
        with tempfile.TemporaryDirectory() as d:
            st = GraphletStore(os.path.join(d, 'sp'))
            h = st.put_parsed(from_response(r, resp), r['graphlet'], {'seeds': []},
                              resident=False)
            b = LocalBudget()
            st.standalone_text(h, budget=b)
            acc = b.usage()['memory_bytes']
            limit = LocalBudget(memory_mb=acc / float(1 << 20))
            peak, txt = _traced(lambda: st.standalone_text(h, budget=limit))
            self.assertLessEqual(peak, acc)
            with self.assertRaises(LocalBudgetExceeded):
                st.standalone_text(h, budget=LocalBudget(memory_mb=(acc - 1) / float(1 << 20)))


# ======================================================================= O3

def _dropped_body(dropped, runs, nseed=1000):
    """The review's realistic case: a 1,000-bp seed, |dropped| X records with the server's
    reason (seed_unsupported) and |runs| presence runs each."""
    from test_traverse_scale import comb_annotate
    t = comb_annotate(1)
    seq = ('ACGTTGCA' * (nseed // 8 + 1))[:nseed]
    t = t.replace('S aaaaaaaaaaaaaaaa 10 6 0 ACGTACGTAC',
                  'S aaaaaaaaaaaaaaaa %d %d 0 %s' % (nseed, nseed - 4, seq))
    t = t.rstrip('\n').split('\n')[:-1]
    i_s = next(i for i, line in enumerate(t) if line.startswith('S '))
    xs = []
    for d in range(dropped):
        pairs = ','.join('%d-%d' % (300 + 16 * j, 300 + 16 * j + 7) for j in range(runs))
        xs.append('X seed_unsupported %s requested_label_%06d' % (pairs, d))
    t[i_s + 1:i_s + 1] = xs
    t.append('Z %d' % (len(t) + 1))
    return '\n'.join(t) + '\n'


def _labels_body(labels):
    from test_traverse_scale import comb_annotate
    t = comb_annotate(1).rstrip('\n').split('\n')[:-1]
    i_o = next(i for i, line in enumerate(t) if line.startswith('O '))
    t[i_o:i_o] = ['L c %d * 0 extra_label_%07d' % (2 + i, i) for i in range(labels)]
    t.append('Z %d' % (len(t) + 1))
    return '\n'.join(t) + '\n'


class TestToJsonEnvelope(unittest.TestCase):
    def mk(self, text):
        return from_response({'graphlet': text}, {})

    def test_many_dropped_labels_stop_under_one_mebibyte(self):
        # 10,000 X records of 40 runs: a traced peak of ~34 MB that completed under
        # memory_mb=1 with an account of 30 KB
        text = _dropped_body(10_000, 40)
        with self.assertRaises(LocalBudgetExceeded) as cm:
            self.mk(text).to_json(budget=LocalBudget(memory_mb=1))
        self.assertEqual('memory', cm.exception.stop.resource)

    def test_the_account_bounds_the_peak_and_the_work_follows(self):
        for name, text in (('dropped', _dropped_body(2000, 40)),
                           ('dropped_1run', _dropped_body(2000, 1)),
                           ('labels', _labels_body(20_000)),
                           ('plain', _dropped_body(0, 0))):
            b = LocalBudget()
            want = self.mk(text).to_json(budget=b)
            g = self.mk(text)
            peak, got = _traced(g.to_json)
            self.assertEqual(want, got, name)
            self.assertGreaterEqual(b.peak_bytes, peak, name)
        # the work grows with the records it builds (one per dropped label and run)
        small, large = LocalBudget(), LocalBudget()
        self.mk(_dropped_body(100, 40)).to_json(budget=small)
        self.mk(_dropped_body(2000, 40)).to_json(budget=large)
        self.assertGreater(large.used_work - small.used_work, 1900 * 40)

    def test_unbudgeted_unchanged(self):
        text = _dropped_body(50, 3)
        self.assertEqual(self.mk(text).to_json(), self.mk(text).to_json(budget=LocalBudget()))


# ======================================================================= X1

class TestDerivationUnits(unittest.TestCase):
    """The derivation limitation's observed: k-mers read on a result walked from a partial
    derivation (D3), elapsed ms on a seed whose derivation failed (no graphlet)."""

    def test_a_walked_result_states_kmers(self):
        from test_traverse_coordinates import fresh
        from metagraph.traverse import derive
        g = fresh('d3_kmer')
        lim = derive.partial_derivation(g)
        self.assertEqual([], g.check_rules())
        self.assertLessEqual(lim.observed, g.seed.num_kmers)
        self.assertEqual({'partial': True, 'kmers_read': lim.observed,
                          'of': g.seed.num_kmers}, g.summary()['seed']['derivation'])

    def test_a_failed_seed_passes_its_ms_on(self):
        failed = {'seed': {'seed_id': 's'}, 'error': 'the time budget ran out while deriving '
                  'the permitted set', 'outcome': {'walks': 'failed'},
                  'limitations': [{'kind': 'derivation', 'knob': 'bounds.time_budget_ms',
                                   'limit': 20, 'observed': 12032,
                                   'effect': 'no result: elapsed ms'}]}
        from metagraph.traverse import TraverseClient
        resp = TraverseClient.response({'results': [failed]})
        self.assertEqual([None], resp.graphlets)
        self.assertEqual(failed['limitations'], resp.errors[0]['result']['limitations'])

        class Failing(M.FakeClient):
            def traverse_raw(self, request):
                return {'release': '', 'results': [copy.deepcopy(failed)]}
        with tempfile.TemporaryDirectory() as d:
            t = GraphletTools(GraphletStore(d), {'mini': Failing()})
            out = t.traverse_fetch(seed={'sequence': 'ACGTACGT'})
        self.assertEqual('seed_failed', out['error'])
        self.assertEqual(failed['limitations'], out['limitations'])


# ================================================= next_request() with a change_cost table

def _wide_forbid():
    import traverse_golden as TG
    for key, r, resp in TG.items_of(TG.golden_files()):
        if key == 'wide_forbid_r40.graphlet.json#0':
            return r, resp
    raise AssertionError('the golden wide_forbid_r40 retrieval is missing')


class TestChangeCostTable(unittest.TestCase):
    """The review of the P2 fixes (U15): with labels.change_cost a table, the account of
    next_request() did not grow with its entries -- neither the copies of the entries
    (the merged strategy, the request's) nor the switch search's table were charged, and 2
    lwu per entry were: 39,531 B against a traced peak of 387,409 B at 900 entries."""

    def setUp(self):
        self.r, self.resp = _wide_forbid()
        self.names = sorted({l.name for l in from_response(self.r, self.resp).labels})

    def table(self, extra=0, cost=9):
        entries = [['X%d' % i, 'Y%d' % i, 1] for i in range(extra)]
        entries += [[a, b, cost] for a in self.names for b in self.names]
        return {'model': 'table', 'default': 0.25, 'entries': entries}

    def call(self, many, **kw):
        g = from_response(self.r, self.resp)
        if many:
            return g.next_requests('right', [0], **kw)
        return g.next_request('right', [0], **kw)

    def test_the_account_bounds_the_peak(self):
        # (a call without a table peaks near 50 KB, where the interpreter's free lists move
        # the traced figure by a few KB either way: TestMemoryModel looks above 50 KB too)
        for many in (False, True):
            for table in (self.table(), self.table(cost=0.1), self.table(extra=5000),
                          self.table(extra=20000)):
                kw = {'labels': {'change_cost': table}}
                with self.subTest(many=many, entries=len(table['entries'])):
                    lb = LocalBudget()
                    want = self.call(many, budget=lb, **kw)
                    peak, got = _traced(lambda: self.call(many, **kw))
                    self.assertEqual(want, got)
                    self.assertGreater(peak, 300_000)
                    self.assertLessEqual(peak, lb.peak_bytes, (peak, lb.peak_bytes))

    def test_the_charge_grows_with_the_entries(self):
        def used(table):
            lb = LocalBudget()
            self.call(False, budget=lb, labels={'change_cost': table})
            return lb.used_work, lb.peak_bytes
        w0, m0 = used({'model': 'table', 'default': 0.25, 'entries': []})
        w1, m1 = used(self.table())
        k = len(self.names) ** 2
        # per entry: two copies (6 lwu each), two checks and the searches' check, count
        # and build (work model 1: 3.1 lwu per entry, the memory not at all)
        self.assertGreaterEqual((w1 - w0) / k, 10)
        self.assertGreaterEqual((m1 - m0) / k, 250)
        # a memory limit between the two stops the call with the entries, not without
        mb = (m0 + m1) / 2 / (1 << 20)
        self.call(False, budget=LocalBudget(memory_mb=mb),
                  labels={'change_cost': {'model': 'table', 'default': 0.25, 'entries': []}})
        with self.assertRaises(LocalBudgetExceeded) as cm:
            self.call(False, budget=LocalBudget(memory_mb=mb),
                      labels={'change_cost': self.table()})
        self.assertEqual(('memory', 'next_request'),
                         (cm.exception.stop.resource, cm.exception.stop.op))

    def test_the_switch_search_alone(self):
        # the reviewer's direct probe: 20 sources with an entry to each of n targets
        from metagraph.traverse import ops
        for n in (50, 500, 2000):
            src = ['S%05d' % i for i in range(n)]
            tgt = ['T%05d' % i for i in range(n)]
            cost = {'model': 'table', 'default': 1,
                    'entries': [[a, b, 9] for a in src[:20] for b in tgt]}
            lb = LocalBudget()
            with lb.scope('probe'):
                want = ops._switch_reach(cost, src, tgt, 3.0, lb=lb)
            peak, got = _traced(lambda: ops._switch_reach(cost, src, tgt, 3.0))
            self.assertEqual(want, got)
            self.assertLessEqual(peak, lb.peak_bytes - B.CALL_BASE, n)


# ================================================================ texts

class TestTexts(unittest.TestCase):
    """The library's texts follow the contract: D3 is feature level 6 (SPEC §7.0; builds
    7aaee760..67bef367 state 5 and are never deployed), and a registered attempt_id is
    refused while the server holds it, not "used up" (a 400 from a server below level 5
    registers it; it runs again once that server no longer retains it)."""

    def read(self, *parts):
        with open(os.path.join(*parts)) as f:
            return ' '.join(f.read().split())

    def test_the_texts(self):
        import re
        lib = os.path.join(T.API, 'metagraph', 'traverse')
        rst = self.read(T.REPO, 'docs', 'source', 'graphlets.rst')
        texts = {'graphlets.rst': rst}
        for name in ('client.py', 'attempts.py', 'mcp_tools.py', 'store.py'):
            texts[name] = self.read(lib, name)
        used_up = re.compile(r'\bus(?:e|es|ed|ing) (?:that |the |its |an |it )?'
                             r'(?:attempt_id |id )?up\b|\bid is used up|\bruns once\b')
        for name, text in texts.items():
            self.assertIsNone(used_up.search(text), (name, used_up.search(text)))
        self.assertNotIn('(feature level 5+)', rst)
        self.assertIn('after ``observed`` of its k-mers (feature level 6; builds '
                      '``7aaee760`` to ``67bef367`` carry it while stating level 5 and are '
                      'never deployed)', rst)
        self.assertIn('A level-5 server fails such a seed instead', rst)
        self.assertIn('registers its attempt_id before it answers the 400 (the id is '
                      'refused while that server retains it, and runs again after that)',
                      texts['client.py'])


if __name__ == '__main__':
    unittest.main()
