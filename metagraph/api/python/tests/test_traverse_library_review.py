"""The review of the Python graphlet library made at 6d97d6e8 (2026-11-14), re-checked at
9726cb420 (train 6, lane W5): each finding that held, with the case that failed before.

  1. the MCP tools read a multi-graph server's /traverse answer (the envelope of one pair,
     SPEC-multigraph-envelope.md section 1) as TraverseClient does: results[0] and the
     envelope the graphlet keeps are the pair's, in traverse_fetch and traverse_continue;
  2. a continuation (next_request(), deepen(), traverse_continue), a replay (request_of(),
     traverse_fetch(replay=)) and the client keep the pair a retrieval was answered on --
     graph, graph_path, annotation_path -- and annotation_path is a selector (in a `pairs`
     object, the only place /traverse takes it), never a strategy field;
  3. a view of one arm answers claims() for that arm only, as walks() and to_fasta() do;
  4. frames.records() builds its walk rows from the per-leaf records, never whole walks
     (their chains and bases are quadratic on a comb): the same rows.
"""

import copy
import gc
import glob
import json
import os
import sys
import tempfile
import tracemalloc
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402
import traverse_golden as TG  # noqa: E402
from test_traverse_scale import comb_annotate  # noqa: E402
from test_traverse_mcp import SEED, FakeClient  # noqa: E402

import metagraph.traverse as MT  # noqa: E402
from metagraph.traverse import GraphletStore, LocalBudget, TraverseClient, ops  # noqa: E402
from metagraph.traverse import frames, parse  # noqa: E402
from metagraph.traverse.mcp_tools import GraphletTools  # noqa: E402

PAIR = {'graph': 'E', 'graph_path': '/g/E', 'annotation_path': '/a/E'}


def _items():
    """Every committed retrieval with a graphlet: the documents and the golden inputs."""
    files = sorted(glob.glob(os.path.join(T.DOCS, '**', '*.graphlet.json'), recursive=True))
    return TG.items_of(files + TG.golden_files())


def _peak(fn):
    gc.collect()
    tracemalloc.start()
    base = tracemalloc.get_traced_memory()[0]
    try:
        out = fn()
        return tracemalloc.get_traced_memory()[1] - base, out
    finally:
        tracemalloc.stop()


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


# ======================================================================= finding 2

class _Session:
    """A requests-like session: a multi-graph server answering /traverse for the pair
    the request names (the merge fixture's answer) and the probe of one pair."""

    def __init__(self):
        self.calls = []

    def request(self, method, url, data=None, headers=None, timeout=None):
        body = json.loads(data) if data else None
        self.calls.append((method, url, body))
        if url.endswith('/traverse'):
            pairs = body.get('pairs')
            if pairs:
                pair = dict(pairs[0])
            else:
                pair = {'graph': body.get('graph', 'E'),
                        'graph_path': body.get('graph_path', '/g/E'),
                        'annotation_path': '/a/E'}
            pair.setdefault('graph_path', '/g/' + pair['graph'])
            pair.setdefault('annotation_path', '/a/' + pair['graph'])
            return _Response(200, envelope(T.doc_json('merge', 'graphlet'), pair))
        if '/traverse/capabilities' in url:
            return _Response(200, {'release': '', 'index_fp': None})
        return _Response(404, {'error': 'no route'})

    def bodies(self):
        return [b for m, u, b in self.calls if m == 'POST']


class _Response:
    def __init__(self, status, obj):
        self.status_code = status
        self.content = json.dumps(obj).encode()
        self.headers = {'Content-Type': 'application/json'}


def _selection(req):
    return {k: req[k] for k in ('graph', 'graph_path', 'annotation_path', 'pairs', 'graphs')
            if k in req}


PAIRS = {'pairs': [PAIR]}


class TestTheSourceSelection(unittest.TestCase):
    def client(self, **kw):
        s = _Session()
        return TraverseClient('h', 1, session=s, **kw), s

    def retrieval(self, client, **kw):
        return client.traverse([T.doc_json('merge', 'request')['seeds'][0]],
                               T.doc_json('merge', 'request')['strategy'], **kw).graphlets[0]

    def test_a_continuation_keeps_the_pair_it_was_answered_on(self):
        c, s = self.client()
        g = self.retrieval(c, graph='E')
        # before: no selection at all (a multi-graph server refuses it, or runs it
        # on another pair)
        self.assertEqual(PAIRS, _selection(g.next_request('right', [0])))
        self.assertEqual(PAIRS, _selection(g.next_requests('right', [0])[0]))

    def test_deepen_continues_where_the_graphlet_came_from(self):
        c, s = self.client(graph='B')
        g = self.retrieval(c, graph='E')
        self.assertEqual({'graph': 'E'}, _selection(s.bodies()[-1]))
        c.deepen(g, 'right', [0])
        # before: the client's own graph, B (another index)
        self.assertEqual(PAIRS, _selection(s.bodies()[-1]))
        # a selection the caller gives wins
        c.deepen(g, 'right', [0], graph='F')
        self.assertEqual({'graph': 'F'}, _selection(s.bodies()[-1]))

    def test_a_single_graph_retrieval_continues_on_the_clients_graph(self):
        g = T.graphlet('switch_chain')
        self.assertEqual({}, _selection(g.next_request('right', [1])))
        c, s = self.client(graph='B')
        c.deepen(g, 'right', [1])
        self.assertEqual({'graph': 'B'}, _selection(s.bodies()[-1]))
        c, s = self.client(graph='B', graph_path='/g/B', annotation_path='/a/B')
        c.deepen(g, 'right', [1])
        self.assertEqual({'pairs': [{'graph': 'B', 'graph_path': '/g/B',
                                     'annotation_path': '/a/B'}]},
                         _selection(s.bodies()[-1]))

    def test_annotation_path_is_a_selector_not_a_strategy_field(self):
        g = T.graphlet('switch_chain')
        req = g.next_request('right', [1], graph='E', annotation_path='/a/E')
        # before: strategy.annotation_path, which the server refuses
        self.assertNotIn('annotation_path', req['strategy'])
        self.assertEqual({'pairs': [{'graph': 'E', 'annotation_path': '/a/E'}]},
                         _selection(req))
        req = g.next_request('right', [1], graph='E', graph_path='/g/E')
        self.assertEqual({'graph': 'E', 'graph_path': '/g/E'}, _selection(req))
        # the caller's selection replaces the retrieval's whole
        c, _ = self.client()
        m = self.retrieval(c, graph='E')
        self.assertEqual({'graph': 'F'}, _selection(m.next_request('right', [0], graph='F')))
        self.assertEqual({}, _selection(m.next_request('right', [0], graph=None)))
        with self.assertRaises(ValueError):
            g.next_request('right', [1], annotation_path='/a/E')     # names no graph

    def test_the_client_names_a_pair(self):
        c, s = self.client()
        seq = T.doc_json('merge', 'request')['seeds'][0]
        req = c.build_request([seq], graph='E', graph_path='/g/E', annotation_path='/a/E')
        self.assertEqual(PAIRS, _selection(req))
        self.assertEqual({'graph': 'E', 'graph_path': '/g/E'},
                         _selection(c.build_request([seq], graph='E', graph_path='/g/E')))
        own, own_session = self.client(graph='E', graph_path='/g/E', annotation_path='/a/E')
        self.assertEqual(PAIRS, _selection(own.build_request([seq])))
        self.assertEqual({'graph': 'F'}, _selection(own.build_request([seq], graph='F')))
        self.retrieval(c, graph='E', graph_path='/g/E', annotation_path='/a/E')
        self.assertEqual(PAIRS, _selection(s.bodies()[-1]))
        # the probe names the pair too
        own.capabilities()
        self.assertEqual(
            'http://h:1/traverse/capabilities?graph=E&graph_path=%2Fg%2FE&annotation_path=%2Fa%2FE',
            own_session.calls[-1][1])

    def test_request_of_keeps_the_annotation_path(self):
        c, _ = self.client()
        g = self.retrieval(c, graph='E')
        # before: graph and graph_path only (the pair is ambiguous where a graph has
        # several annotations)
        self.assertEqual(PAIRS, _selection(GraphletStore.request_of(g)))
        self.assertEqual({'graph': 'F'}, _selection(GraphletStore.request_of(g, graph='F')))


class TestTheToolsKeepThePair(ToolsPair):
    def test_continue_and_replay_on_the_pair(self):
        client = MultiGraphFake()
        client.graph = 'E'
        tools = self.tools(client)
        h = tools.traverse_fetch(seed={'sequence': SEED['switch_chain']})['handle']
        self.assertEqual({'graph': 'E'}, _selection(client.requests[-1]))
        done = tools.traverse_continue(h, 'right', 1)
        self.assertNotIn('error', done, done)
        # the pair the parent was answered on, and its capabilities checked
        self.assertEqual(PAIRS, _selection(client.requests[-1]))
        self.assertEqual(('E', '/g/E', '/a/E'), client.caps_asked[-1])
        self.assertEqual(PAIRS, _selection(self.multi.get(done['handle']).request))
        again = tools.traverse_fetch(replay=done['handle'])
        self.assertNotIn('error', again, again)
        self.assertEqual(PAIRS, _selection(client.requests[-1]))


# ======================================================================= finding 3

class TestAViewKeepsItsArm(unittest.TestCase):
    def test_claims_of_another_arm(self):
        g = T.graphlet('fork')
        v = g.subgraph(['acc2'], arm='right')
        self.assertEqual((['right'], [], ''), (v.arms(), v.walks('left'), v.to_fasta('left')))
        # before: the left arm's claim of acc2
        self.assertEqual([], v.claims('left'))
        self.assertEqual([], v.claims('l', budget=LocalBudget()))
        right = [c.as_dict() for c in g.claims('right', ['acc2'])]
        self.assertTrue(right)
        for arm in (None, 'right', 'r'):
            self.assertEqual(right, [c.as_dict() for c in v.claims(arm)])
        both = g.subgraph(['acc2'])
        self.assertEqual([c.as_dict() for c in g.claims('left', ['acc2'])],
                         [c.as_dict() for c in both.claims('left')])
        with self.assertRaises((KeyError, ValueError)):
            v.claims('up')


# ======================================================================= finding 4

def _walk_row(side, w):
    """The row frames.records() made from a whole Walk before (the reference)."""
    return {'arm': side, 'walk': w.path_id, 'leaf': w.leaf, 'length_bp': w.length_bp,
            'path_reason': w.path_reason, 'n_alive': w.n_alive,
            'n_full': len(w.labels_full), 'complete': w.complete,
            'beyond_certified_bp': w.beyond_certified_bp}


class TestTheTableRows(unittest.TestCase):
    def test_the_walk_rows_are_the_walks_fields(self):
        n = 0
        for _, result, response in _items():
            g = MT.from_response(result, response)
            want = [_walk_row(side, w) for side in ('left', 'right') if side in g.arms
                    for w in g.walks(side, by='id')]
            got = frames.records(MT.from_response(result, response))['walks']
            self.assertEqual(want, got)
            n += len(got)
        self.assertGreater(n, 100)

    def test_no_walk_is_built(self):
        g = parse(comb_annotate(60))
        want = [_walk_row('right', w) for w in g.walks('right', by='id')]
        g = parse(comb_annotate(60))
        # before: every Walk was built (its chain, its spelling, its claims)
        with mock.patch.object(ops, '_walk_objects', side_effect=AssertionError('a Walk')), \
                mock.patch.object(MT.derive, 'walk_batch', side_effect=AssertionError('a chain')):
            rec = frames.records(g)
        self.assertEqual(want, rec['walks'])

    def test_the_peak_is_what_the_call_keeps(self):
        # on a comb of 2,000 splits (Python 3.11): before, a traced peak of 25.5 MB for
        # 7.9 MB kept (the rows and the derivations cached), the walks' chains and bases
        # quadratic in the splits; now 8.0 MB for 7.6 MB kept
        g = parse(comb_annotate(2000))
        gc.collect()
        tracemalloc.start()
        try:
            base = tracemalloc.get_traced_memory()[0]
            rec = frames.records(g)
            kept, peak = (x - base for x in tracemalloc.get_traced_memory())
        finally:
            tracemalloc.stop()
        self.assertEqual(2001, len(rec['walks']))
        self.assertLess(peak, 1.5 * kept)


if __name__ == '__main__':
    unittest.main()
