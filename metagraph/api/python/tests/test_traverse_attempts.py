"""Stage 4 of DESIGN-traverse-graphlet.md §14, backend half: the library's side of a
ledger-managed attempt.

  - TraverseClient sends `attempt_id`, `budget_id` and `locus_id` (the frozen wire contract's
    request fields, §14.1) and nothing when they are not given, so a request without them is
    the one it always was;
  - TraverseResponse.usage reads the response-level usage block (None without attempt_id);
  - cancel() and attempt() return the server's answer, a 404 included (nothing runs under the
    id: finished, unknown or tombstoned), and raise on the other errors;
  - a graphlet stopped from outside the walk carries only new values of free tokens (Q scope
    `attempt`, resources `cancelled` / `attempt_deadline`, the K knob `attempt_id`): it parses
    and dumps byte for byte (MGT v1 unchanged).

The server side (cancel mid-walk, a gone client, the bound) is in the integration tests.
"""

import json
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from metagraph.traverse import TraverseClient, TraverseError, dump, parse  # noqa: E402

DOCS = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'traverse', 'documents')


class _Response:
    def __init__(self, status, body):
        self.status_code = status
        self.content = json.dumps(body).encode()
        self.headers = {'Content-Type': 'application/json'}


class _Session:
    """answers by (method, path suffix); records every request"""

    def __init__(self, answers):
        self.answers = answers
        self.sent = []

    def request(self, method, url, data=None, headers=None, timeout=None):
        self.sent.append((method, url, json.loads(data) if data else None))
        for (m, suffix), (status, body) in self.answers.items():
            if m == method and url.endswith(suffix):
                return _Response(status, body)
        return _Response(404, {'error': 'no route'})


class TestAttemptFields(unittest.TestCase):

    def test_requests_carry_the_attempt_fields_only_when_given(self):
        client = TraverseClient('h', 1)
        plain = client.build_request(['ACGT'])
        for key in ('attempt_id', 'budget_id', 'locus_id'):
            self.assertNotIn(key, plain)
        req = client.build_request(['ACGT'], attempt_id='a-1', budget_id='b-1', locus_id='l-1')
        self.assertEqual(('a-1', 'b-1', 'l-1'),
                         (req['attempt_id'], req['budget_id'], req['locus_id']))

        session = _Session({('POST', '/traverse'): (200, {'results': [],
                                                          'usage': {'attempt_id': 'a-2'}})})
        client = TraverseClient('h', 1, session=session)
        resp = client.traverse(['ACGT'], attempt_id='a-2', locus_id='l-2')
        self.assertEqual('a-2', session.sent[0][2]['attempt_id'])
        self.assertEqual('l-2', session.sent[0][2]['locus_id'])
        self.assertNotIn('budget_id', session.sent[0][2])
        self.assertEqual({'attempt_id': 'a-2'}, resp.usage)

        session = _Session({('POST', '/traverse'): (200, {'results': []})})
        self.assertIsNone(TraverseClient('h', 1, session=session).traverse(['ACGT']).usage)

    def test_cancel_and_attempt_return_the_answer(self):
        stopping = {'attempt_id': 'a-3', 'cancelled': True, 'state': 'stopping'}
        unknown = {'attempt_id': 'zz', 'cancelled': False, 'state': 'unknown', 'tombstone': True,
                   'error': "unknown attempt_id 'zz'"}
        session = _Session({('POST', '/traverse/cancel'): (200, stopping),
                            ('GET', '/traverse/attempt/a-3'): (200, {'state': 'finished',
                                                                    'reason': 'cancelled'}),
                            ('GET', '/traverse/attempt/zz'): (404, {'state': 'unknown'})})
        client = TraverseClient('h', 1, session=session)
        self.assertEqual(stopping, client.cancel('a-3', wait_ms=2000))
        self.assertEqual(('POST', {'attempt_id': 'a-3', 'wait_ms': 2000}),
                         (session.sent[-1][0], session.sent[-1][2]))
        client.cancel('a-3')
        self.assertEqual({'attempt_id': 'a-3'}, session.sent[-1][2])
        self.assertEqual('cancelled', client.attempt('a-3')['reason'])
        self.assertEqual('unknown', client.attempt('zz')['state'])
        # a 404 to a cancel is an answer: nothing runs under the id
        session.answers[('POST', '/traverse/cancel')] = (404, unknown)
        self.assertEqual(unknown, client.cancel('zz'))
        # other errors still raise
        session.answers[('POST', '/traverse/cancel')] = (400, {'error': 'request.wait_ms'})
        with self.assertRaises(TraverseError) as e:
            client.cancel('zz', wait_ms=10**6)
        self.assertEqual(400, e.exception.status)
        # a 404 of /traverse itself is no answer
        session.answers[('POST', '/traverse')] = (404, {'error': 'Could not find path'})
        with self.assertRaises(TraverseError):
            client.traverse(['ACGT'])


class TestAttemptDocuments(unittest.TestCase):

    def test_a_stop_from_outside_the_walk_round_trips(self):
        """A work-stopped fixture restated as a cancelled attempt: the Q and K records hold
        new values of free tokens only, which the reader accepts and the writer reproduces"""
        with open(os.path.join(DOCS, 'budgets', 'work_stop.mgt')) as f:
            lines = f.read().split('\n')
        q = next(i for i, l in enumerate(lines) if l.startswith('Q '))
        k = next(i for i, l in enumerate(lines)
                 if l.startswith('K r walk_domain bounds.max_work_units '))
        lines[q] = ('Q attempt cancelled traversal i:70000 i:70000 i:1004 i:68996 '
                    'continue_from_leaves the attempt was cancelled (POST /traverse/cancel) and '
                    'the walk stopped at the next checkpoint, 1004 ms after the request was '
                    'received (the bound: 70000 ms), at 41 bp on the right arm: every walk up to '
                    "each arm's complete_to_bp is present and the heads not expanded end with "
                    'resource_limit; no budget of the request ran out; a continuation from a leaf '
                    'is a new traversal')
        fields = lines[k].split(' ', 8)
        fields[3:6] = ['attempt_id', 'i:70000', 'i:1004']
        fields[8] = ('every walk of at most complete_to_bp bp is present, longer walks may be '
                     'missing: the attempt was cancelled (POST /traverse/cancel) and the '
                     'exploration stopped at the next checkpoint (resource_limit; see '
                     'resource_stop; limit: the attempt\'s bound, observed: its elapsed time, '
                     'ms); no budget of the request ran out: a new attempt can continue from '
                     'the leaves')
        lines[k] = ' '.join(fields)
        text = '\n'.join(lines)
        g = parse(text)
        self.assertEqual(text, dump(g))
        stop = g.resource_stop
        self.assertEqual(('attempt', 'cancelled', 'traversal', ['continue_from_leaves']),
                         (stop.scope, stop.resource, stop.phase, stop.actions))
        self.assertIn('attempt_id', [lim.knob for lim in g.limitations])


if __name__ == '__main__':
    unittest.main()
