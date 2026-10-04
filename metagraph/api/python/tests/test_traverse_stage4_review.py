"""The reviews of the stage-3 fixes and of the stage-4 backend, Python side:

  P3 (stage-3 fixes)  next_request()'s note on a label left out searched chains only through
       the continued label and the labels left out: a chain through a KEPT extra label, which
       one uninterrupted walk could take, was missed, and the note understated what the
       continuation may miss;
  F4   the attempt-bound 503 (`{error, usage}`, no Retry-After: the attempt ran and its id is
       used up) was raised as ServerInitializing (a loading index, nothing run), whose
       documented recovery -- retry -- meets a 409; an error's usage was reachable only
       through .body;
  MCP  traverse_capabilities() answered result_too_large at its default ceiling (2048 bytes)
       once the server described its attempts (2.4 KB);
  F6   a cancel the server could not tombstone for its whole retention period is a 429
       answer (tombstone: false), returned like the 404s.
"""

import copy
import dataclasses
import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402

from metagraph.traverse import (AttemptAtBound, GraphletStore, ServerInitializing,  # noqa: E402
                                TraverseClient, TraverseError, ops)
from metagraph.traverse.mcp_tools import (CAPABILITIES_MAX_BYTES, DEFAULT_MAX_BYTES,  # noqa: E402
                                          GraphletTools, _size)
from test_traverse_mcp import FakeClient  # noqa: E402


# ------------------------------------------------------------------ P3 (stage-3 fixes)

class TestTheNoteSearchesEveryChain(unittest.TestCase):
    """switch_chain, walk 1, continued under C (loss 2) and D (loss 0) with loss_budget 3:
    the request's budget is 1, so A (D -> B 0.8, B -> A 0.8) is left out, but one
    uninterrupted walk could still have entered it from D through B, a kept extra label."""

    def next_request(self, entries):
        g = T.graphlet('switch_chain')
        names = {l.name: l for l in g.labels}
        original = ops._continuation_seeds

        def two_labels(g_, a_, leaves):
            conts, _ = original(g_, a_, leaves)
            c = dataclasses.replace(conts[0], labels=[names['C'], names['D']],
                                    losses=(2.0, 0.0))
            return [c], [c.as_seed()]

        g.envelope['strategy']['labels']['loss_budget'] = 3
        g.envelope['strategy']['labels']['change_cost'] = {
            'model': 'table', 'default': 'forbid', 'entries': entries}
        ops._continuation_seeds = two_labels
        try:
            return g.next_request('right', [1])
        finally:
            ops._continuation_seeds = original

    def test_a_chain_through_a_kept_extra_label(self):
        req = self.next_request([['D', 'B', 0.8], ['B', 'A', 0.8]])
        self.assertEqual(1.0, req['strategy']['labels']['loss_budget'])
        self.assertEqual(['B'], req['strategy']['labels']['extra'])
        self.assertEqual([('A', 'unreachable', 1)],
                         [(x['name'], x['why'], x['walk']) for x in req.left_out])
        note, = [n for n in req.notes if 'leaves out' in n]
        self.assertIn("'A' (reachable for 'D' within its own remaining budget)", note)
        self.assertIn('one uninterrupted walk could still have entered it', note)
        # the library's chain search agrees: D reaches A through B within its own 3
        self.assertEqual({'A': 1.6, 'B': 0.8},
                         ops._switch_reach(req['strategy']['labels']['change_cost'], ['D'],
                                           ['A', 'B', 'C'], 3.0))

    def test_no_chain_within_any_remaining_budget(self):
        req = self.next_request([['D', 'B', 0.8], ['B', 'A', 2.5]])
        note, = [n for n in req.notes if 'leaves out' in n]
        self.assertNotIn('reachable for', note)
        self.assertNotIn('uninterrupted walk could still', note)


# ------------------------------------------------------------------ F4

class _Response:
    def __init__(self, status, body, headers=None):
        self.status_code = status
        self.content = json.dumps(body).encode()
        self.headers = dict({'Content-Type': 'application/json'}, **(headers or {}))


class _Session:
    def __init__(self, status, body, headers=None):
        self.answer = (status, body, headers)

    def request(self, method, url, data=None, headers=None, timeout=None):
        return _Response(*self.answer)


USAGE = {'attempt_id': 'svc-1', 'reason': 'deadline', 'work_units': 7,
         'bound_ms': 51, 'seeds': {'requested': 1, 'started': 1, 'finished': 1, 'abandoned': 0}}


class TestTheBound503IsNotALoadingServer(unittest.TestCase):

    def client(self, status, body, headers=None):
        return TraverseClient('h', 1, session=_Session(status, body, headers))

    def test_the_bound(self):
        body = {'error': 'the attempt reached the duration bound', 'usage': USAGE}
        with self.assertRaises(AttemptAtBound) as cm:
            self.client(503, body).traverse(['ACGT'], attempt_id='svc-1')
        e = cm.exception
        self.assertNotIsInstance(e, ServerInitializing)
        self.assertIsInstance(e, TraverseError)
        self.assertEqual((503, USAGE), (e.status, e.usage))

    def test_a_loading_server(self):
        body = {'error': 'Server is currently initializing, please come back later.'}
        with self.assertRaises(ServerInitializing) as cm:
            self.client(503, body, {'Retry-After': '60'}).traverse(['ACGT'])
        self.assertEqual((60, None), (cm.exception.retry_after, cm.exception.usage))
        # a 503 without usage stays a loading server, also without the header
        with self.assertRaises(ServerInitializing):
            self.client(503, body).traverse(['ACGT'])

    def test_an_error_carries_its_usage(self):
        body = {'error': "unknown field 'bogus'", 'usage': dict(USAGE, reason='error')}
        with self.assertRaises(TraverseError) as cm:
            self.client(400, body).traverse(['ACGT'], attempt_id='svc-2')
        self.assertEqual('error', cm.exception.usage['reason'])
        with self.assertRaises(TraverseError) as cm:
            self.client(400, {'error': 'request.attempt_id: expected a string'}).traverse(['A'])
        self.assertIsNone(cm.exception.usage)

    def test_the_tool_does_not_say_retry(self):
        class Bound(FakeClient):
            def capabilities(self):
                raise AttemptAtBound(503, 'the attempt reached its bound', {'usage': USAGE})
        with tempfile.TemporaryDirectory() as root:
            tools = GraphletTools(GraphletStore(root), {'mini': Bound()})
            out = tools.traverse_capabilities(index='mini')
        self.assertEqual('backend_error', out['error'])
        self.assertNotIn('still loading', json.dumps(out))
        self.assertIn('new attempt_id', out.get('hint', ''))


# ------------------------------------------------------------------ MCP capabilities

def _capabilities(size):
    """a server's capabilities of about |size| bytes, shaped like GET /traverse/capabilities"""
    caps = {'release': 'r1', 'schema_version': 1, 'budgets': ['max_memory_mb', 'max_work_units'],
            'memory_bound': 'soft', 'content_encodings': ['gzip', 'deflate'],
            'attempts': {'fields': ['attempt_id', 'budget_id', 'locus_id'],
                         'cancel': 'POST /traverse/cancel', 'retention_s': 3600}}
    caps['work_bound'] = 'w' * max(0, size - len(json.dumps(caps)) - 20)
    return caps


class TestCapabilitiesFitTheDefaultCeiling(unittest.TestCase):

    def tools(self, caps, root):
        class Big(FakeClient):
            def capabilities(self):
                return copy.deepcopy(caps)
        return GraphletTools(GraphletStore(root), {'mini': Big()})

    def test_the_servers_description_is_returned_whole(self):
        for size in (2400, 12000):
            caps = _capabilities(size)
            with tempfile.TemporaryDirectory() as root:
                out = self.tools(caps, root).traverse_capabilities(index='mini')
            self.assertNotIn('error', out, (size, str(out)[:200]))
            self.assertEqual(caps, out['capabilities'])
            self.assertGreater(_size(out), DEFAULT_MAX_BYTES)
            self.assertLessEqual(_size(out), CAPABILITIES_MAX_BYTES)

    def test_an_explicit_ceiling_still_holds(self):
        with tempfile.TemporaryDirectory() as root:
            out = self.tools(_capabilities(2400), root).traverse_capabilities(
                index='mini', max_bytes=DEFAULT_MAX_BYTES)
        self.assertEqual('result_too_large', out['error'])
        self.assertLessEqual(_size(out), DEFAULT_MAX_BYTES)


# ------------------------------------------------------------------ F6

class TestACancelThatPromisesNothingIsAnAnswer(unittest.TestCase):

    def test_429_is_returned(self):
        body = {'error': "unknown attempt_id 'z', and it was NOT tombstoned", 'attempt_id': 'z',
                'cancelled': False, 'state': 'unknown', 'tombstone': False}
        out = TraverseClient('h', 1, session=_Session(429, body)).cancel('z')
        self.assertEqual((False, False), (out['cancelled'], out['tombstone']))


if __name__ == '__main__':
    unittest.main()
