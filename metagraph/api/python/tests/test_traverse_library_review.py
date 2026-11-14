"""The review of the Python graphlet library made at 6d97d6e8 (2026-11-14), re-checked at
9726cb420 (train 6, lane W5): each finding that held, with the case that failed before.

  1. the MCP tools read a multi-graph server's /traverse answer (the envelope of one pair,
     SPEC-multigraph-envelope.md section 1) as TraverseClient does: results[0] and the
     envelope the graphlet keeps are the pair's, in traverse_fetch and traverse_continue.
"""

import copy
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402
from test_traverse_mcp import BY_SEED, SEED, FakeClient  # noqa: E402

from metagraph.traverse import GraphletStore  # noqa: E402
from metagraph.traverse.mcp_tools import GraphletTools  # noqa: E402

PAIR = {'graph': 'E', 'graph_path': '/g/E', 'annotation_path': '/a/E'}


def envelope(resp, pair=PAIR, refusal=None):
    """|resp| as a release-10 multi-graph server answers /traverse for one pair: the
    answer an entry tagged with the pair, its usage moved to the envelope; a refused
    entry instead with |refusal|."""
    resp = copy.deepcopy(resp)
    usage = resp.pop('usage', None)
    if refusal is None:
        entry = dict(resp, index_fp=None, outcome='answered', **pair)
    else:
        entry = dict(pair, index_fp=None, outcome='refused', refusal=refusal)
    out = {'envelope_version': 1, 'route': '/traverse', 'graphs': [pair['graph']],
           'answered': int(refusal is None), 'refused': int(refusal is not None),
           'not_reached': [], 'unfinished': [], 'request_id': None,
           'server_instance': 'f' * 16,
           'deadline': {'not_after_ms': None, 'finish_by_ms': None,
                        'effective_deadline_at': None, 'cutoff_at': None},
           'answers': [entry], 'timing': {'elapsed_ms': 1.0}, 'work': {'memory_bytes': 0}}
    if usage is not None:
        out['usage'] = usage
    return out


class MultiGraphFake(FakeClient):
    """The MCP tests' fake backend as a multi-graph server: every answer in the envelope of
    one pair (|refusal|: the pair refused); the capabilities of any pair named."""

    def __init__(self, refusal=None):
        super().__init__()
        self.refusal = refusal
        self.caps_asked = []

    def traverse_raw(self, request):
        return envelope(super().traverse_raw(request), refusal=self.refusal)

    def capabilities(self, graph=None, graph_path=None, annotation_path=None):
        self.caps_asked.append((graph, graph_path, annotation_path))
        return super().capabilities()


class ToolsPair(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.single = GraphletStore(os.path.join(self.tmp.name, 'single'))
        self.multi = GraphletStore(os.path.join(self.tmp.name, 'multi'))

    def tearDown(self):
        self.tmp.cleanup()

    def tools(self, client, store=None):
        return GraphletTools(store or self.multi, {'mini': client})


# ======================================================================= finding 1

class TestToolsReadTheMultiGraphEnvelope(ToolsPair):
    def test_a_fetch_stores_the_pairs_graphlet(self):
        one = self.tools(FakeClient(), self.single).traverse_fetch(
            seed={'sequence': SEED['merge']})
        out = self.tools(MultiGraphFake()).traverse_fetch(seed={'sequence': SEED['merge']})
        # before: {"error": "empty_response"} -- results[] is inside answers[0]
        self.assertNotIn('error', out, out)
        self.assertEqual(one['evidence'], out['evidence'])
        g = self.multi.graphlet(out['handle'])
        self.assertEqual(self.single.graphlet(one['handle']).dump(envelope=False),
                         g.dump(envelope=False))
        # the graphlet keeps the pair's answer as its envelope, not the envelope's members
        self.assertEqual(PAIR, {k: g.envelope.get(k) for k in PAIR})
        self.assertNotIn('answers', g.envelope)
        self.assertEqual(T.doc_json('merge', 'graphlet')['strategy'], g.envelope['strategy'])

    def test_a_continuation_reads_it_too(self):
        tools = self.tools(MultiGraphFake())
        h = tools.traverse_fetch(seed={'sequence': SEED['switch_chain']})['handle']
        done = tools.traverse_continue(h, 'right', 1)
        self.assertNotIn('error', done, done)
        self.assertEqual(h, done['parent']['handle'])
        self.assertEqual(PAIR, {k: self.multi.graphlet(done['handle']).envelope.get(k)
                                for k in PAIR})

    def test_a_refused_pair_is_the_backends_error(self):
        refusal = {'http_status': 400, 'error': 'Bad request: the pair is not served'}
        tools = self.tools(MultiGraphFake(refusal))
        out = tools.traverse_fetch(seed={'sequence': SEED['merge']})
        self.assertEqual(('backend_error', 400, refusal['error']),
                         (out.get('error'), out.get('status'), out.get('message')), out)
        self.assertEqual([], self.multi.list())

    def test_a_failed_seed_of_a_pair_is_seed_failed(self):
        out = self.tools(MultiGraphFake()).traverse_fetch(seed={'sequence': 'ACGTACGTAC'})
        self.assertEqual('seed_failed', out.get('error'), out)
        self.assertEqual('no carrier', out['message'])


if __name__ == '__main__':
    unittest.main()
