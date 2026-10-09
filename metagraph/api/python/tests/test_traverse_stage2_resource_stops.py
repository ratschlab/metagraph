"""The request budgets during the walk (DESIGN-traverse-graphlet.md §14.1) as a client sees
them: the request budgets (bounds.max_memory_mb, bounds.max_work_units), the resource stop
they produce (JSON `resource_stop`, the MGT `Q` record, `resource_limit` ends) and the
memory_bound_soft limitation, read by the library exactly as it is (MGT v1 is frozen: the
records, fields and tokens are all MGT v1's).

The documents are CLI fixtures under documents/budgets/ (graphlet_fixtures.py
--from-cli; --check keeps them current). A work stop is the same walk in every detail, so
its graphlet must rebuild its full response; a memory stop is charged in the requested
detail, so the same request stops at another depth in detail full.
"""

import copy
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402

from metagraph.traverse import (  # noqa: E402
    GraphletStore, dump, from_response, is_canonical, parse,
)
from metagraph.traverse import ops  # noqa: E402
from metagraph.traverse.export import normalize_result  # noqa: E402
from metagraph.traverse.mcp_tools import GraphletTools  # noqa: E402

BUDGETS = os.path.join(T.DOCS, 'budgets')
# name -> the details written for it
FIXTURES = {
    'work_stop': ('full', 'graphlet'),
    'memory_soft': ('full', 'graphlet'),
    'memory_stop': ('graphlet',),
    'memory_stop_full': ('full',),
    # the budget-aware annotation reads (test_traverse_stage3_decode.py): a stop of a budget-aware annotation read
    'decode_stop': ('graphlet',),
}


def text(name):
    return T.read(os.path.join(BUDGETS, name + '.mgt'))


def response(name, detail):
    with open(os.path.join(BUDGETS, '%s.%s.json' % (name, detail)), encoding='utf-8') as f:
        return json.load(f)


def graphlet(name):
    resp = response(name, 'graphlet')
    return from_response(resp['results'][0], resp)


class TestBudgetDocuments(unittest.TestCase):

    def test_all_budget_fixtures_are_present(self):
        names = {f.split('.')[0] for f in os.listdir(BUDGETS)}
        self.assertEqual(set(FIXTURES), names)
        for name, details in FIXTURES.items():
            for detail in details:
                self.assertTrue(os.path.exists(os.path.join(BUDGETS, '%s.%s.json'
                                                            % (name, detail))), name)

    def test_bodies_round_trip_and_are_canonical(self):
        for name, details in FIXTURES.items():
            if 'graphlet' not in details:
                continue
            with self.subTest(name):
                body = text(name)
                self.assertEqual(body, dump(parse(body)))
                self.assertTrue(is_canonical(body))
                self.assertEqual(body, response(name, 'graphlet')['results'][0]['graphlet'])

    def test_detail_independent_stops_rebuild_the_full_result(self):
        # work is charged on the walk alone, and a memory budget that admits every head
        # changes nothing: the graphlet reproduces the detail full response
        for name in ('work_stop', 'memory_soft'):
            with self.subTest(name):
                full = response(name, 'full')['results'][0]
                mine = normalize_result(graphlet(name).to_json())
                self.assertEqual(normalize_result(full), mine,
                                 '\n'.join(T.json_diff(mine, normalize_result(full))[:20]))
                self.assertEqual([], graphlet(name).check_rules(full))

    def test_the_q_record_is_decoded_with_typed_values(self):
        q = parse(text('work_stop')).resource_stop
        self.assertEqual(('locus', 'work', 'traversal'), (q.scope, q.resource, q.phase))
        self.assertEqual((1500, 1500), (q.requested, q.effective))
        self.assertIsInstance(q.used, int)
        self.assertGreater(q.used, q.effective)
        # the check interval bounds the overshoot (W in GET /traverse/capabilities)
        self.assertLessEqual(q.used - q.effective, 65536)
        self.assertEqual(0, q.remaining)
        self.assertEqual(['raise_work_budget', 'continue_from_leaves'], q.actions)
        self.assertIn('bounds.max_work_units', q.message)
        m = parse(text('memory_stop')).resource_stop
        self.assertEqual(('locus', 'memory', 'traversal'), (m.scope, m.resource, m.phase))
        self.assertEqual((1, 1), (m.requested, m.effective))
        self.assertLessEqual(m.used, m.effective)
        self.assertIn('raise_memory_budget', m.actions)
        self.assertIn('use_graphlet', m.actions)
        self.assertIsNone(parse(text('memory_soft')).resource_stop)

    def test_the_json_states_the_same_stop(self):
        for name, detail in (('work_stop', 'full'), ('work_stop', 'graphlet'),
                             ('memory_stop', 'graphlet'), ('memory_stop_full', 'full')):
            with self.subTest(name + ' ' + detail):
                r = response(name, detail)['results'][0]
                stop = r['resource_stop']
                self.assertEqual({'scope', 'resource', 'phase', 'requested', 'effective', 'used',
                                  'remaining', 'actions', 'message'}, set(stop))
                if 'graphlet' in r:
                    self.assertEqual(ops.resource_stop_dict(parse(r['graphlet']).resource_stop),
                                     stop)

    def test_a_stop_is_a_partial_walk_with_its_limitations(self):
        for name, knob in (('work_stop', 'bounds.max_work_units'),
                           ('memory_stop', 'bounds.max_memory_mb')):
            with self.subTest(name):
                g = parse(text(name))
                self.assertEqual('partial', g.outcome.walks)
                self.assertEqual('inline', g.outcome.delivery)
                walk = [l for l in g.limitations if l.kind == 'walk_domain']
                self.assertEqual(1, len(walk))
                self.assertEqual(('right', knob), (walk[0].arm, walk[0].knob))
                self.assertGreater(walk[0].observed, walk[0].limit)
                arm = g.arms['right']
                self.assertEqual('truncated', arm.status)
                self.assertEqual(walk[0].complete_to_bp, arm.complete_to_bp)
                # the library keeps the MGT code of the reason ('Y' = resource_limit)
                self.assertEqual('Y', arm.cap_trigger.reason)
                # the heads not expanded end with resource_limit (T Y)
                self.assertIn('\nT Y', text(name))
                # the conservative outcome rule finds nothing to object to
                self.assertEqual([], [f for f in g.check_rules() if f['field'].startswith('outcome')])

    def test_memory_bound_soft_is_stated_under_a_memory_budget(self):
        for name in ('memory_soft', 'memory_stop'):
            g = parse(text(name))
            soft = [l for l in g.limitations if l.kind == 'memory_bound_soft']
            self.assertEqual(1, len(soft), name)
            self.assertEqual((None, 'bounds.max_memory_mb', 1), (soft[0].arm, soft[0].knob,
                                                                 soft[0].limit))
            self.assertIsInstance(soft[0].observed, int)
        # it belongs to no outcome class: a budget that admits everything is complete
        self.assertEqual('complete', parse(text('memory_soft')).outcome.walks)
        self.assertEqual([], [l for l in parse(text('work_stop')).limitations
                              if l.kind == 'memory_bound_soft'])

    def test_the_output_detail_is_charged(self):
        # the same request and budget: detail full spells every leaf's chain (quadratic on
        # a comb) and is refused much earlier than the graphlet, which is the stated action
        g = parse(text('memory_stop'))
        full = response('memory_stop_full', 'full')['results'][0]
        self.assertLess(full['arms']['right']['complete_to_bp'], g.arms['right'].complete_to_bp)
        self.assertEqual('memory', full['resource_stop']['resource'])
        self.assertIn('use_graphlet', full['resource_stop']['actions'])
        # the delivered prefix is exact: every path's chain is its leaf's first-parent chain
        arm = full['arms']['right']
        for p in arm['paths']:
            chain = p['segments']
            for parent, child in zip(chain, chain[1:]):
                self.assertEqual(parent, arm['segments'][child]['parents'][0])
            self.assertEqual([], arm['segments'][chain[0]]['parents'])

    def test_work_units_round_trip_through_a(self):
        g = parse(text('work_stop'))
        counters = g.arms['right'].counters
        full = response('work_stop', 'full')['results'][0]['arms']['right']['counters']
        self.assertEqual(full['work_units'], counters['work_units'])
        # the stop's `used` is the seed's: the seed phase (the derivation's reads, 8 per seed
        # k-mer and 1 per annotation entry, as the walk charges a fetched row) and the arm's
        # own units. Every seed k-mer of this fixture carries exactly its four labels.
        seed = response('work_stop', 'full')['results'][0]['seed']
        self.assertEqual(seed['num_kmers'] * (8 + seed['labels_supporting_total']),
                         g.resource_stop.used - counters['work_units'])
        # a request without a budget has no such counter (its responses are unchanged)
        self.assertNotIn('work_units', T.body('merge').arms['right'].counters)

    def test_the_budgets_are_echoed_only_when_given(self):
        self.assertEqual(1500, response('work_stop', 'full')['strategy']['bounds']['max_work_units'])
        self.assertEqual(1, response('memory_soft', 'full')['strategy']['bounds']['max_memory_mb'])
        self.assertNotIn('max_work_units', response('memory_soft', 'full')['strategy']['bounds'])
        merge = T.doc_json('merge', 'full')['strategy']['bounds']
        self.assertNotIn('max_memory_mb', merge)
        self.assertNotIn('max_work_units', merge)


class TestAgentView(unittest.TestCase):
    """What an agent reads: the evidence block and the summary carry the stop and its
    actions; nothing reads as complete."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.store = GraphletStore(self.tmp.name)
        self.tools = GraphletTools(self.store, {})

    def tearDown(self):
        self.tmp.cleanup()

    def test_evidence_block_lists_the_resource_stop(self):
        for name in ('work_stop', 'memory_stop'):
            ev = ops.evidence_block(graphlet(name))
            self.assertIn('resource_stop', ev['limitations']['seed'], name)
            self.assertIn('walk_domain', ev['limitations']['right'], name)
            self.assertEqual('partial', ev['outcome']['walks'])

    def test_summary_caveats_carry_the_stop_and_its_actions(self):
        h = self.store.put(response('work_stop', 'graphlet'), None)
        caveats = self.tools.graphlet_summary(h, max_bytes=8192)['summary']['caveats']
        stop = [c for c in caveats if c['kind'] == 'resource_stop']
        self.assertEqual(1, len(stop))
        self.assertEqual(('work', 'traversal'), (stop[0]['resource'], stop[0]['phase']))
        self.assertEqual(['raise_work_budget', 'continue_from_leaves'], stop[0]['actions'])
        walk = [c for c in caveats if c['kind'] == 'walk_domain'][0]
        self.assertEqual('bounds.max_work_units', walk['knob'])
        full = self.tools.graphlet_summary(h, detail='limitations', max_bytes=8192)
        self.assertEqual('work', full['resource_stop']['resource'])
        m = self.store.put(response('memory_stop', 'graphlet'), None)
        kinds = [c['kind'] for c in self.tools.graphlet_summary(m, max_bytes=8192)['summary']['caveats']]
        self.assertIn('memory_bound_soft', kinds)
        self.assertIn('resource_stop', kinds)


def _cli():
    path = os.path.join(T.REPO, 'build_debug', 'metagraph')
    return path if os.path.exists(path) else None


@unittest.skipUnless(_cli(), 'needs build_debug/metagraph')
class TestLiveCli(unittest.TestCase):
    """The work budget against the CLI: every budget gives a certified prefix, and more
    budget never certifies less."""

    @classmethod
    def setUpClass(cls):
        cls.G = T.fixture_script()
        cls.root = tempfile.mkdtemp(prefix='stage2-')
        cls.graph, cls.anno = cls.G.build_index('bubbles', cls.G.resolve_cli(_cli()), cls.root)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.root, ignore_errors=True)

    def run_request(self, budget, detail):
        fx = self.G.BUDGET_FIXTURES['work_stop']
        r = self.G.materialize(fx['request'], self.G.blocks_of(self.G.INDEXES['bubbles']))
        r = copy.deepcopy(r)
        r['strategy']['bounds']['max_work_units'] = budget
        r['strategy']['output']['detail'] = detail
        return self.G.run_traverse(self.G.resolve_cli(_cli()), self.graph, self.anno, r)

    def test_more_budget_certifies_at_least_as_far(self):
        last = -1
        for budget in (100, 400, 800, 1600, 3200, 6400):
            out = self.run_request(budget, 'graphlet')
            g = from_response(out['results'][0], out)
            arm = g.arms['right']
            self.assertGreaterEqual(arm.complete_to_bp, last, budget)
            last = arm.complete_to_bp
            self.assertEqual(g.dump(), out['results'][0]['graphlet'])
            full = self.run_request(budget, 'full')['results'][0]
            self.assertEqual(normalize_result(full), normalize_result(g.to_json()), budget)
            if g.resource_stop is not None:
                self.assertEqual('partial', g.outcome.walks)
                self.assertGreater(g.resource_stop.used, budget)
        self.assertEqual(100, last, 'the largest budget completes the walk')


if __name__ == '__main__':
    unittest.main()
