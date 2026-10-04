"""The external recheck of the stage-2 fixes (GPT, 278a53dd..d8d3ef86), Python side, with
the design answers adopted from it:

  P3a  a spooled fetch stored its entry and answered result_too_large without the handle:
       the receipt was checked with delivery 'inline', one byte shorter than 'spooled';
  P3b  next_request() over several walks took left_out from the first walk only, naming a
       label another walk seeded as unreachable;
  D1   the continuation's branch allowance restarted at the seed (a continuation from 40
       reached 100 where one uninterrupted walk stopped at 46): next_request() reduces
       branching.max_label_branches by the largest terminal branch count, keeps
       "unlimited", states it (branch_budget, notes) and offers reset_branches;
  D2   labels.extra kept only the labels one switch from a seed label reaches: a chain of
       switches within the loss budget is enough (the walk enforces the cumulative loss),
       as the server now accepts.

The CLI classes run build_debug/metagraph on the review3 indexes and skip without it.
"""

import copy
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402
import test_traverse_review3 as R  # noqa: E402

from metagraph.traverse import GraphletStore  # noqa: E402
from metagraph.traverse import ops  # noqa: E402
from metagraph.traverse.mcp_tools import GraphletTools, _Receipt, _size  # noqa: E402
from test_traverse_mcp import FakeClient, SEED  # noqa: E402


def tearDownModule():
    R.tearDownModule()


# ------------------------------------------------------------------ P3a

class TestSpooledReceiptIsCheckedAsStored(unittest.TestCase):
    """The receipt of a spooled fetch says delivery 'spooled': checked as it will be
    stated, a ceiling one byte too small for it stores nothing, and one that fits stores
    the entry and returns its handle."""

    def fetch(self, max_bytes):
        with tempfile.TemporaryDirectory() as root:
            store = GraphletStore(root)
            tools = GraphletTools(store, {'mini': FakeClient()})
            out = tools.traverse_fetch(index='mini', seed={'sequence': SEED['switch_chain']},
                                       max_graphlet_mb=0, max_bytes=max_bytes)
            return dict(out), [(e['handle'], e['delivery']) for e in store.list()]

    def test_the_boundary(self):
        inline = _size(_Receipt({'index': 'mini', 'delivery': 'inline',
                                 'handle': 'g_' + '0' * 12},
                                ('summary', 'evidence', 'graphlet_bytes')).minimal())
        out, stored = self.fetch(inline)
        self.assertIn('error', out)
        self.assertNotIn('handle', out)
        self.assertEqual([], stored)           # no entry the caller is not told about
        out, stored = self.fetch(inline + 1)
        self.assertNotIn('error', out, out)
        self.assertEqual('spooled', out['delivery'])
        self.assertEqual([(out['handle'], 'spooled')], stored)


# ------------------------------------------------------------------ P3b

class TestLeftOutIsPerWalk(unittest.TestCase):
    """fork, left arm: walk 0 continues under acc1 and acc2, walk 1 under acc3, and forbid
    lets no walk switch. Each walk's left_out names what IT leaves out, with the walk."""

    def test_the_order_of_the_walks_does_not_matter(self):
        g = T.graphlet('fork')
        a = g.next_request('left', [0, 1])
        b = g.next_request('left', [1, 0])
        key = lambda r: sorted((x['ref'], x['why'], x['walk']) for x in r.left_out)
        self.assertEqual(key(a), key(b))
        self.assertTrue(a.left_out)

    def test_a_label_a_walk_seeds_is_not_left_out_for_it(self):
        g = T.graphlet('fork')
        req = g.next_request('left', [0, 1])
        seeded = {c: set(s['labels']) for c, s in zip((0, 1), req['seeds'])}
        for x in req.left_out:
            self.assertIn('walk', x)
            self.assertNotIn(x['name'], seeded[x['walk']], x)
        # every label of the retrieval is, for each walk, a seed label, an extra label or
        # left out for that walk
        for w in (0, 1):
            left = {x['ref'] for x in req.left_out if x['walk'] == w}
            for l in g.labels:
                if l.name not in seeded[w] and l.name not in req['strategy']['labels']['extra']:
                    self.assertIn(l.ref, left, (w, l))

    def test_one_walk_is_unchanged_but_for_the_walk(self):
        g = T.graphlet('fork')
        req = g.next_request('left', [1])
        self.assertEqual({1}, {x['walk'] for x in req.left_out})


# ------------------------------------------------------------------ D1

class TestBranchAllowanceIsReduced(unittest.TestCase):
    """switch_chain, walk 1: its lineage used one branch of max_label_branches 1."""

    def test_reduced_by_the_largest_terminal_count(self):
        req = T.graphlet('switch_chain').next_request('right', [1])
        self.assertEqual(0, req['strategy']['branching']['max_label_branches'])
        self.assertEqual({'original': 1, 'effective': 0, 'largest_terminal_branches': 1,
                          'reset': False}, req.branch_budget)
        # one lineage: exact, so no conservative note (stated in branch_budget)
        self.assertFalse([n for n in req.notes if 'max_label_branches' in n], req.notes)

    def test_reset_keeps_the_original_and_says_so(self):
        req = T.graphlet('switch_chain').next_request('right', [1], reset_branches=True)
        self.assertEqual(1, req['strategy']['branching']['max_label_branches'])
        self.assertEqual(True, req.branch_budget['reset'])
        self.assertEqual(1, req.branch_budget['effective'])
        note, = [n for n in req.notes if 'reset_branches=True' in n]
        self.assertIn('may branch further', note)

    def test_unlimited_stays_unlimited(self):
        g = T.graphlet('switch_chain')
        g.envelope['strategy']['branching']['max_label_branches'] = 'unlimited'
        req = g.next_request('right', [1])
        self.assertEqual('unlimited', req['strategy']['branching']['max_label_branches'])
        self.assertEqual('unlimited', req.branch_budget['effective'])

    def test_no_allowance_states_nothing(self):
        g = T.graphlet('switch_chain')
        g.envelope['strategy']['branching']['max_label_branches'] = 0
        req = g.next_request('right', [1])
        self.assertEqual(0, req['strategy']['branching']['max_label_branches'])
        self.assertIsNone(req.branch_budget)

    def test_the_callers_value_is_final_and_stated(self):
        g = T.graphlet('switch_chain')
        req = g.next_request('right', [1], branching={'max_label_branches': 5})
        self.assertEqual(5, req['strategy']['branching']['max_label_branches'])
        self.assertEqual(5, req.branch_budget['effective'])
        self.assertTrue(any("is the caller's" in n for n in req.notes), req.notes)
        with self.assertRaises(ValueError):
            g.next_request('right', [1], branching={'max_label_branches': 'many'})
        with self.assertRaises(ValueError):
            g.next_request('right', [1], reset_branches='yes')

    def test_the_tool_states_it_and_passes_reset(self):
        with tempfile.TemporaryDirectory() as root:
            tools = GraphletTools(GraphletStore(root), {'mini': FakeClient()})
            h = tools.traverse_fetch(index='mini', seed={'sequence': SEED['switch_chain']})
            h = h['handle']
            out = tools.traverse_continue(h, 'right', 1, execute=False)
            self.assertEqual({'original': 1, 'effective': 0, 'largest_terminal_branches': 1,
                              'reset': False}, out['branch_budget'])
            self.assertEqual(0, out['request']['strategy']['branching']['max_label_branches'])
            out = tools.traverse_continue(h, 'right', 1, execute=False, reset_branches=True)
            self.assertTrue(out['branch_budget']['reset'])
            self.assertEqual(1, out['request']['strategy']['branching']['max_label_branches'])
            self.assertEqual('bad_argument', tools.traverse_continue(
                h, 'right', 1, execute=False, reset_branches=1)['error'])
            # executed: the receipt carries it, and it is never cut
            done = tools.traverse_continue(h, 'right', 1)
            self.assertNotIn('error', done, done)
            self.assertIn('branch_budget', done)
            self.assertNotIn('branch_budget', done.get('fields_cut', []))


@unittest.skipUnless(R._cli(), 'needs build_debug/metagraph')
class TestBranchAllowanceAgainstOneWalk(unittest.TestCase):
    """The recheck's case: bubbles, the 'merge' request with max_label_branches 1. One
    uninterrupted walk to 100 ends both.fa at 46 (branch); the continuation of walk 0 from
    40 reached 100 when the allowance restarted, and now ends it where the walk did."""

    def test_the_continuation_stops_where_one_walk_does(self):
        G = R._Indexes.G()
        req = G.materialize(G.CLI_FIXTURES['merge']['request'], G.blocks_of(G.INDEXES['bubbles']))
        req['strategy']['direction'] = 'right'
        req['strategy']['output']['detail'] = 'graphlet'
        req['strategy']['branching']['max_label_branches'] = 1
        req['strategy']['bounds']['max_extension_bp'] = 40
        g = R._Indexes.graphlet('bubbles', req)
        req['strategy']['bounds']['max_extension_bp'] = 100
        whole = R._Indexes.graphlet('bubbles', req)
        q = g.next_request('right', [0], bp=60)
        self.assertEqual(0, q['strategy']['branching']['max_label_branches'])
        child = R._Indexes.graphlet('bubbles', q)

        def ends(x, offset):
            return {rr.to_bp + offset for rr in x.arms['right'].runs
                    if x.labels[rr.label].name == 'both.fa' and rr.reason == 'branch'}
        self.assertEqual({46}, ends(whole, 0))
        self.assertEqual({46}, ends(child, 40))
        reached = max(rr.to_bp + 40 for rr in child.arms['right'].runs
                      if child.labels[rr.label].name == 'both.fa')
        self.assertLessEqual(reached, 46)


# ------------------------------------------------------------------ D2

class TestExtraLabelsByChains(unittest.TestCase):
    """switch_chain walk 1 continues under C, with loss_budget 1 left: C -> D 0.5 and
    D -> A 0.5 reach A through D within it, which one switch did not."""

    def chain(self, entries):
        g = T.graphlet('switch_chain')
        g.envelope['strategy']['labels']['change_cost'] = {
            'model': 'table', 'default': 'forbid', 'entries': entries}
        return g.next_request('right', [1])

    def test_a_chain_within_the_budget_is_kept(self):
        req = self.chain([['C', 'D', 0.5], ['D', 'A', 0.5]])
        self.assertEqual(1.0, req['strategy']['labels']['loss_budget'])
        self.assertEqual(['A', 'D'], req['strategy']['labels']['extra'])
        self.assertEqual([], T.request_violations(req))
        left = [(x['name'], x['why']) for x in req.left_out]
        self.assertEqual([('B', 'unreachable')], left)
        note, = [n for n in req.notes if 'leaves out' in n]
        self.assertIn('no chain of switches', note)

    def test_a_chain_beyond_the_budget_is_left_out(self):
        req = self.chain([['C', 'D', 0.5], ['D', 'A', 0.75]])
        self.assertEqual(['D'], req['strategy']['labels']['extra'])
        self.assertEqual([], T.request_violations(req))
        self.assertEqual({'A', 'B'}, {x['name'] for x in req.left_out})

    def test_the_default_reaches_through_another_label(self):
        # C -> A is explicitly too dear, the default reaches A through B or D
        g = T.graphlet('switch_chain')
        g.envelope['strategy']['labels']['change_cost'] = {
            'model': 'table', 'default': 0.5, 'entries': [['C', 'A', 9]]}
        req = g.next_request('right', [1])
        self.assertEqual(['A', 'B', 'D'], req['strategy']['labels']['extra'])
        self.assertEqual([], T.request_violations(req))

    def test_switch_reach(self):
        table = {'model': 'table', 'default': 'forbid',
                 'entries': [['A', 'B', 1], ['B', 'C', 1], ['C', 'D', 5], ['A', 'B', 0.5]]}
        # the last entry for a pair wins, as on the server
        self.assertEqual({'B': 0.5, 'C': 1.5}, ops._switch_reach(table, ['A'], ['B', 'C', 'D'], 2))
        self.assertEqual({}, ops._switch_reach({'model': 'forbid'}, ['A'], ['B'], 9))
        self.assertEqual({'B': 1.0}, ops._switch_reach({'model': 'constant', 'value': 1},
                                                       ['A'], ['B'], 1))
        self.assertIsNone(ops._switch_reach({'model': 'hll'}, ['A'], ['B'], 1))
        # a sink ends a chain but does not continue one
        self.assertEqual({'B': 0.5}, ops._switch_reach(table, ['A'], ['B', 'C'], 2,
                                                       sinks={'B'}))
        with self.assertRaises(ValueError):
            ops._switch_reach({'model': 'table', 'entries': [['A', 'B', 'x']]}, ['A'], ['B'], 1)


@unittest.skipUnless(R._cli(), 'needs build_debug/metagraph')
class TestTheServerAcceptsAChain(unittest.TestCase):
    """The recheck's case on the chain index: extra B and C, A -> B 1, B -> C 1, loss
    budget 2 — refused before (C is two switches away), accepted now; within 1.5 C is
    refused, by name."""

    def request(self, budget):
        req = T.doc_json('switch_chain', 'request')
        req.pop('fixture_index', None)
        lab = req['strategy']['labels']
        lab['extra'] = ['B', 'C']
        lab['loss_budget'] = budget
        lab['change_cost'] = {'model': 'table', 'default': 'forbid',
                              'entries': [['A', 'B', 1], ['B', 'C', 1]]}
        req['strategy'].setdefault('output', {})['detail'] = 'graphlet'
        return req

    def test_accepted_within_the_budget(self):
        rc, out = R._Indexes.run('chain', self.request(2))
        self.assertEqual(0, rc, out)
        self.assertNotIn('error', out['results'][0])
        self.assertEqual([], T.request_violations(self.request(2)))

    def test_refused_by_name_beyond_it(self):
        rc, out = R._Indexes.run('chain', self.request(1.5))
        self.assertNotEqual(0, rc)
        text = out['error'] if isinstance(out, dict) else str(out)
        self.assertIn("'C'", text)
        self.assertIn('no chain of switches', text)
        self.assertTrue(T.request_violations(self.request(1.5)))


if __name__ == '__main__':
    unittest.main()
