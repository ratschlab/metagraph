"""Stage 3 of DESIGN-traverse-graphlet.md §14.1 as a client sees it: on an annotation whose
reads are budget-aware (a row-diff annotation over BRWT or ColumnMajor), a read that does
not fit the memory the request has left stops the walk in phase `annotation_decode`. The
library reads it exactly as it is: MGT v1 is frozen, and the Q record's phase and actions
are free `[a-z_]` tokens, the K effects free text, so no record, field or token is new.

The document is the CLI fixture documents/budgets/decode_stop (graphlet_fixtures.py
--from-cli, index `dense`: the row before a k-mer repeated 120,000 times carries one label,
but its row-diff dependency rows hold 120,000 coordinates).
"""

import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402

from metagraph.traverse import GraphletStore, dump, from_response, is_canonical, parse  # noqa: E402
from metagraph.traverse import ops  # noqa: E402
from metagraph.traverse.mcp_tools import GraphletTools  # noqa: E402

BUDGETS = os.path.join(T.DOCS, 'budgets')


def text():
    return T.read(os.path.join(BUDGETS, 'decode_stop.mgt'))


def response():
    with open(os.path.join(BUDGETS, 'decode_stop.graphlet.json'), encoding='utf-8') as f:
        return json.load(f)


class TestDecodeStop(unittest.TestCase):

    def test_the_body_round_trips_and_is_canonical(self):
        body = text()
        self.assertEqual(body, dump(parse(body)))
        self.assertTrue(is_canonical(body))
        self.assertEqual(body, response()['results'][0]['graphlet'])

    def test_the_q_record_carries_the_phase_and_the_levers(self):
        q = parse(text()).resource_stop
        self.assertEqual(('locus', 'memory', 'annotation_decode'), (q.scope, q.resource, q.phase))
        self.assertEqual((1, 1), (q.requested, q.effective))
        # a read that does not fit is refused whole: the walk never held more than the budget
        self.assertLessEqual(q.used, q.effective)
        # on a row-diff annotation a read decodes a whole row with its dependency rows: what
        # reads fewer rows is the lever, a smaller radius or fewer named labels is not
        self.assertIn('more_selective_seed', q.actions)
        self.assertIn('raise_memory_budget', q.actions)
        self.assertEqual('continue_from_leaves', q.actions[-1])
        for lever in ('lower_max_seed_labels', 'name_fewer_labels', 'label_constrained_query'):
            self.assertNotIn(lever, q.actions)    # the last only in annotate mode
        self.assertIn('row-diff dependency rows', q.message)
        self.assertIn('annotation.batch_kmers', q.message)

    def test_the_json_states_the_same_stop(self):
        r = response()['results'][0]
        stop = r['resource_stop']
        self.assertEqual('annotation_decode', stop['phase'])
        self.assertEqual(ops.resource_stop_dict(parse(r['graphlet']).resource_stop), stop)

    def test_a_decode_stop_is_a_partial_walk_with_its_limitations(self):
        g = parse(text())
        self.assertEqual('partial', g.outcome.walks)
        walk = [l for l in g.limitations if l.kind == 'walk_domain']
        self.assertEqual(1, len(walk))
        self.assertEqual(('right', 'bounds.max_memory_mb'), (walk[0].arm, walk[0].knob))
        # observed: what admitting the refused row needed beside what the walk held, which is
        # more than the budget
        self.assertGreater(walk[0].observed, walk[0].limit)
        self.assertIn('row-diff dependency rows', walk[0].effect)
        arm = g.arms['right']
        self.assertEqual('truncated', arm.status)
        self.assertEqual(walk[0].complete_to_bp, arm.complete_to_bp)
        self.assertEqual('Y', arm.cap_trigger.reason)
        self.assertIn('\nT Y', text())
        self.assertEqual([], [f for f in g.check_rules() if f['field'].startswith('outcome')])

    def test_memory_bound_soft_names_only_what_is_uncharged(self):
        g = parse(text())
        soft = [l for l in g.limitations if l.kind == 'memory_bound_soft']
        self.assertEqual(1, len(soft))
        self.assertEqual((None, 'bounds.max_memory_mb', 1), (soft[0].arm, soft[0].knob, soft[0].limit))
        self.assertIn('charged before they are held', soft[0].effect)
        self.assertIn("fetched rows until its heads are processed", soft[0].effect)

    def test_an_agent_sees_the_phase(self):
        with tempfile.TemporaryDirectory() as tmp:
            store = GraphletStore(tmp)
            tools = GraphletTools(store, {})
            h = store.put(response(), None)
            caveats = tools.graphlet_summary(h, max_bytes=8192)['summary']['caveats']
            stop = [c for c in caveats if c['kind'] == 'resource_stop']
            self.assertEqual(1, len(stop))
            self.assertEqual(('memory', 'annotation_decode'), (stop[0]['resource'], stop[0]['phase']))
            self.assertIn('more_selective_seed', stop[0]['actions'])
            ev = ops.evidence_block(from_response(response()['results'][0], response()))
            self.assertIn('resource_stop', ev['limitations']['seed'])


if __name__ == '__main__':
    unittest.main()
