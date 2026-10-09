"""Feature level 5 in TraverseClient (SPEC §5 and §10.3), on answers recorded from a
level-5 server_query (data/traverse/level5/answers.json: dcc0cebd on build/mini_refseq, the
level3_* answers from ce949da5; test_traverse_level5_live.py asks a live server the same
questions):

  - cancel(attempt_id, wait_ms, not_after_ms) sends not_after_ms, and expect_server_instance
    (an explicit one, or 'auto': the server's own) goes out with attempt_id only -- both only
    to a server that states feature_level >= 5, read once from GET /capabilities and never
    for a request without them (an older server refuses them with a 400 after it
    registered the request's attempt_id, refused then while that server retains it);
  - the answers are typed: AttemptAnswer (status, suppression, the 429's reason),
    AttemptConflict (the duplicate's and the tombstoned id's 409, with the suppression),
    InstanceMismatch, AttemptExpired's fields;
  - release_verdict() applies the capabilities' release_rule exactly and names the first
    condition an answer does not meet;
  - what GET /capabilities stated is kept per server process (an answer naming another
    server_instance, or a 400 refusing a level-5 field as unknown -- an older binary now at
    the address -- makes the client read it again, and a cancel refused that way goes out
    again without not_after_ms), traverse() states what it sent on a transport error too,
    and a `finished` release names the assumptions it rests on.
"""

import copy
import json
import os
import socket
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402

from metagraph.traverse import (  # noqa: E402
    AttemptAnswer, AttemptConflict, AttemptExpired, AttemptSent, InstanceMismatch,
    ReleaseVerdict, Suppression, TraverseClient, TraverseError, TraverseResponse,
    UnsupportedFeature, release_verdict,
)

with open(os.path.join(T.DATA, 'level5', 'answers.json'), encoding='utf-8') as _f:
    REC = json.load(_f)
A = REC['answers']
NAF = REC['not_after_ms']
INST = REC['server_instance']
SEED = 'ACGACGAAAGCGATGGTATCCCTGCCCGACATATGCAAGAGGCTCATCGGCTAGCGTGTTGCTGAGTTGA'


class _Response:
    def __init__(self, status, body):
        self.status_code = status
        self.content = json.dumps(body).encode()
        self.headers = {'Content-Type': 'application/json'}


class _Replay:
    """A requests-like session: answers by (method, path) with a recorded answer; records
    every request (the body as sent)."""

    def __init__(self, routes):
        self.routes = dict(routes)
        self.sent = []

    def request(self, method, url, data=None, headers=None, timeout=None):
        path = url.split('://', 1)[1].split('/', 1)[1]
        body = json.loads(data) if data else None
        self.sent.append((method, '/' + path, body))
        got = self.routes.get((method, '/' + path))
        if callable(got):
            got = got(body)
        if got is None:
            return _Response(404, {'error': 'no route %s %s' % (method, path)})
        status, answer = got
        return _Response(status, copy.deepcopy(answer))


def _client(routes, **kw):
    session = _Replay(routes)
    return TraverseClient('h', 1, session=session, timeout=5, **kw), session


def _paths(session):
    return [p for _, p, _ in session.sent]


class TestTheGate(unittest.TestCase):
    def test_cancel_sends_not_after_ms_to_a_level5_server_only(self):
        c, s = _client({('GET', '/capabilities'): A['capabilities'],
                        ('POST', '/traverse/cancel'): A['cancel_unknown_not_after']})
        got = c.cancel('rec-1', not_after_ms=NAF)
        self.assertEqual(['/capabilities', '/traverse/cancel'], _paths(s))
        self.assertEqual({'attempt_id': 'rec-1', 'not_after_ms': NAF}, s.sent[-1][2])
        self.assertEqual({'attempt_id': 'rec-1', 'not_after_ms': NAF}, got.sent)
        # the level is read once
        c.cancel('rec-1', 100, NAF)
        self.assertEqual({'attempt_id': 'rec-1', 'wait_ms': 100, 'not_after_ms': NAF},
                         s.sent[-1][2])
        self.assertEqual(1, _paths(s).count('/capabilities'))

        old, s = _client({('GET', '/capabilities'): A['level3_capabilities'],
                          ('POST', '/traverse/cancel'): A['level3_cancel_unknown']})
        got = old.cancel('rec-10', not_after_ms=NAF)
        # the cancel still goes out (it stops the attempt), without the field
        self.assertEqual({'attempt_id': 'rec-10'}, s.sent[-1][2])
        self.assertEqual((404, True, None), (got.status, got.tombstone, got.suppression))
        self.assertNotIn('not_after_ms', got.sent)

    def test_a_request_without_level5_fields_reads_no_level(self):
        c, s = _client({('POST', '/traverse/cancel'): A['cancel_unknown_plain'],
                        ('GET', '/traverse/attempt/rec-1'): A['get_tombstoned'],
                        ('POST', '/traverse'): A['traverse_completed']})
        c.cancel('rec-2')
        c.attempt('rec-1')
        req = c.build_request([SEED], attempt_id='rec-5', not_after_ms=NAF)
        self.assertNotIn('expect_server_instance', req)
        c.traverse([SEED], attempt_id='rec-5', not_after_ms=NAF)
        self.assertNotIn('/capabilities', _paths(s))

    def test_expect_server_instance(self):
        c, s = _client({('GET', '/capabilities'): A['capabilities'],
                        ('POST', '/traverse'): A['traverse_completed']})
        with self.assertRaises(ValueError):
            c.build_request([SEED], expect_server_instance='auto')
        self.assertEqual([], s.sent)                     # refused before anything is read
        req = c.build_request([SEED], attempt_id='a', expect_server_instance='auto')
        self.assertEqual(INST, req['expect_server_instance'])
        req = c.build_request([SEED], attempt_id='a', expect_server_instance='00ff')
        self.assertEqual('00ff', req['expect_server_instance'])
        resp = c.traverse([SEED], attempt_id='rec-5', not_after_ms=NAF,
                          expect_server_instance='auto')
        self.assertEqual(AttemptSent('rec-5', NAF, INST), resp.sent)
        self.assertEqual(INST, s.sent[-1][2]['expect_server_instance'])
        self.assertEqual(1, _paths(s).count('/capabilities'))
        self.assertIsNone(c.traverse([SEED]).sent)       # not an attempt

        # a level given to the client is not read again; the instance still is
        c, s = _client({('GET', '/capabilities'): A['capabilities']}, feature_level=5)
        self.assertEqual(INST, c.build_request(
            [SEED], attempt_id='a', expect_server_instance='auto')['expect_server_instance'])

    def test_below_level5_auto_is_left_out_and_an_explicit_instance_raises(self):
        c, s = _client({('GET', '/capabilities'): A['level3_capabilities'],
                        ('POST', '/traverse'): A['traverse_completed']})
        req = c.build_request([SEED], attempt_id='a', expect_server_instance='auto')
        self.assertNotIn('expect_server_instance', req)
        resp = c.traverse([SEED], attempt_id='rec-5', expect_server_instance='auto')
        self.assertEqual(AttemptSent('rec-5'), resp.sent)
        with self.assertRaises(UnsupportedFeature) as cm:
            c.traverse([SEED], attempt_id='b', expect_server_instance=INST)
        self.assertEqual((3, 5), (cm.exception.feature_level, cm.exception.needed))
        self.assertIsInstance(cm.exception, ValueError)
        self.assertNotIn('b', [b.get('attempt_id') for _, p, b in s.sent if b])

    def test_a_server_without_the_capabilities_route_is_level_0(self):
        c, s = _client({('POST', '/traverse/cancel'): A['level3_cancel_unknown']})
        self.assertEqual(0, c.server_feature_level())
        c.cancel('rec-10', not_after_ms=NAF)
        self.assertEqual({'attempt_id': 'rec-10'}, s.sent[-1][2])


def _unknown_field(name, usage_of=None):
    """The 400 of a server that does not know |name| (ce949da5's for not_after_ms in a
    cancel, data/traverse/level5); |usage_of|: a request with that attempt_id, whose id the
    server registered before the 400 (it carries the usage)."""
    body = {'error': "request: unknown field '%s'" % name}
    if usage_of is not None:
        body['usage'] = {'attempt_id': usage_of, 'reason': 'error',
                         'server_instance': A['level3_cancel_unknown'][1]['server_instance']}
    return 400, body


class TestOneServerProcess(unittest.TestCase):
    """What the client read is the process's: a rollback (an older binary at the same
    address, which cannot answer instance_mismatch) and a restart are noticed."""

    def rolled_back(self, **kw):
        c, s = _client({('GET', '/capabilities'): A['capabilities']}, **kw)
        self.assertEqual(5, c.server_feature_level())
        self.assertEqual(INST, c.server_instance())
        s.routes[('GET', '/capabilities')] = A['level3_capabilities']
        s.routes[('POST', '/traverse/cancel')] = lambda body: (
            A['level3_cancel_with_not_after'] if 'not_after_ms' in body
            else A['level3_cancel_unknown'])
        return c, s

    def test_a_cancel_refused_for_not_after_ms_goes_out_again_without_it(self):
        self.assertEqual("request: unknown field 'not_after_ms'",
                         A['level3_cancel_with_not_after'][1]['error'])
        c, s = self.rolled_back()
        got = c.cancel('rec-10', 100, NAF)
        # the 400 stopped nothing: the cancel went out again, and stopped (tombstoned) it
        self.assertEqual([{'attempt_id': 'rec-10', 'wait_ms': 100, 'not_after_ms': NAF},
                          {'attempt_id': 'rec-10', 'wait_ms': 100}],
                         [b for _, p, b in s.sent if p == '/traverse/cancel'])
        self.assertEqual((404, True), (got.status, got.tombstone))
        self.assertEqual({'attempt_id': 'rec-10', 'wait_ms': 100}, got.sent)
        # the level is read again, and the next cancel goes out without it at once
        self.assertEqual(3, c.server_feature_level())
        c.cancel('rec-10', not_after_ms=NAF)
        self.assertEqual({'attempt_id': 'rec-10'}, s.sent[-1][2])
        self.assertEqual(2, _paths(s).count('/capabilities'))

    def test_a_given_level_gives_way_to_that_refusal(self):
        c, s = self.rolled_back(feature_level=5)
        got = c.cancel('rec-10', not_after_ms=NAF)
        self.assertEqual((404, {'attempt_id': 'rec-10'}), (got.status, got.sent))
        self.assertIsNone(c.feature_level)
        self.assertEqual(3, c.server_feature_level())
        # a 400 for another reason is raised as it was, nothing sent again or forgotten
        c, s = _client({('GET', '/capabilities'): A['capabilities'],
                        ('POST', '/traverse/cancel'): (400, {'error': 'bad attempt_id'})},
                       feature_level=5)
        with self.assertRaises(TraverseError):
            c.cancel('rec-10', not_after_ms=NAF)
        self.assertEqual((1, 5), (_paths(s).count('/traverse/cancel'), c.feature_level))

    def test_an_auto_request_refused_for_the_field_reads_the_level_again(self):
        c, s = self.rolled_back()
        s.routes[('POST', '/traverse')] = lambda body: (
            _unknown_field('expect_server_instance', body['attempt_id'])
            if 'expect_server_instance' in body else A['traverse_completed'])
        # the first request to reach the older binary is lost (its 400 came after the id was
        # registered: refused while that binary retains it) ...
        with self.assertRaises(TraverseError) as cm:
            c.traverse([SEED], attempt_id='rec-11', not_after_ms=NAF,
                       expect_server_instance='auto')
        self.assertEqual(AttemptSent('rec-11', NAF, INST), cm.exception.sent)
        self.assertEqual('finished', release_verdict(cm.exception, cm.exception.sent).code)
        # ... the next one reads the level again and goes out without the field
        resp = c.traverse([SEED], attempt_id='rec-5', not_after_ms=NAF,
                          expect_server_instance='auto')
        self.assertEqual(AttemptSent('rec-5', NAF), resp.sent)
        self.assertNotIn('expect_server_instance', s.sent[-1][2])

    def test_another_server_instance_forgets_what_was_read(self):
        other = dict(A['get_finished'][1], server_instance='00000000000000aa')
        c, s = _client({('GET', '/capabilities'): A['capabilities'],
                        ('GET', '/traverse/attempt/rec-5'): A['get_finished'],
                        ('GET', '/traverse/attempt/rec-6'): (200, other)})
        c.server_instance()
        c.attempt('rec-5')                               # this process: kept
        self.assertIsNotNone(c._server_doc)
        c.attempt('rec-6')                               # a restarted one
        self.assertIsNone(c._server_doc)
        # a usage block and a 409's attempt object name the process too
        for answer in (A['traverse_completed'], A['traverse_tombstoned']):
            c.server_instance()
            body = copy.deepcopy(answer[1])
            for block in (body, body.get('usage'), body.get('attempt')):
                if isinstance(block, dict) and 'server_instance' in block:
                    block['server_instance'] = '00000000000000aa'
            s.routes[('POST', '/traverse')] = (answer[0], body)
            try:
                c.traverse([SEED], attempt_id='rec-1')
            except AttemptConflict:
                pass
            self.assertIsNone(c._server_doc, answer[0])
        n = _paths(s).count('/capabilities')
        c.build_request([SEED], attempt_id='x', expect_server_instance='auto')
        self.assertEqual(n + 1, _paths(s).count('/capabilities'))

    def test_the_level_read_wins_over_a_given_one(self):
        c, s = _client({('GET', '/capabilities'): A['level3_capabilities']}, feature_level=5)
        req = c.build_request([SEED], attempt_id='a', expect_server_instance='auto')
        self.assertNotIn('expect_server_instance', req)
        self.assertEqual(3, c.server_feature_level())


class _Silent:
    """A socket that accepts connections and never answers."""

    def __enter__(self):
        self.sock = socket.socket()
        self.sock.bind(('127.0.0.1', 0))
        self.sock.listen(4)
        return self.sock.getsockname()[1]

    def __exit__(self, *exc):
        self.sock.close()


class TestTransportErrors(unittest.TestCase):
    def test_an_unanswered_attempt_states_what_was_sent(self):
        def hang(body):
            raise TimeoutError('timed out')
        c, s = _client({('GET', '/capabilities'): A['capabilities'],
                        ('POST', '/traverse'): hang})
        with self.assertRaises(TimeoutError) as cm:
            c.traverse([SEED], attempt_id='rec-5', not_after_ms=NAF,
                       expect_server_instance='auto')
        self.assertEqual(AttemptSent('rec-5', NAF, INST), cm.exception.sent)
        # what the ledger then does: cancel it and judge the answer against what was sent
        self.assertEqual('tombstone', release_verdict(
            _answer('cancel_unknown_not_after', 'cancel'),
            AttemptSent('rec-1', NAF, INST)).code)
        # the stdlib transport: a server that never answers
        with _Silent() as port:
            c = TraverseClient('127.0.0.1', port, timeout=0.3, feature_level=5)
            with self.assertRaises(OSError) as cm:
                c.traverse([SEED], attempt_id='rec-12', not_after_ms=NAF,
                           expect_server_instance=INST)
            self.assertEqual(AttemptSent('rec-12', NAF, INST), cm.exception.sent)
        # nothing was sent: no .sent
        c, s = _client({('GET', '/capabilities'): A['level3_capabilities']})
        with self.assertRaises(UnsupportedFeature) as cm:
            c.traverse([SEED], attempt_id='b', expect_server_instance=INST)
        self.assertFalse(hasattr(cm.exception, 'sent'))


class TestTypedAnswers(unittest.TestCase):
    def setUp(self):
        routes = {('GET', '/capabilities'): A['capabilities']}
        self.c, self.s = _client(routes)
        self.routes = self.s.routes

    def cancel(self, name, attempt_id, **kw):
        self.routes[('POST', '/traverse/cancel')] = A[name]
        return self.c.cancel(attempt_id, **kw)

    def test_tombstone_answers(self):
        got = self.cancel('cancel_unknown_not_after', 'rec-1', not_after_ms=NAF)
        self.assertIsInstance(got, AttemptAnswer)
        self.assertIsInstance(got, dict)                 # the answer, untyped
        self.assertEqual(A['cancel_unknown_not_after'][1], dict(got))
        self.assertEqual((404, 'cancel', 'unknown', False, True),
                         (got.status, got.route, got.state, got.cancelled, got.tombstone))
        sup = got.suppression
        self.assertEqual(Suppression(A['cancel_unknown_not_after'][1]['suppressed_until_ms'],
                                     NAF, True, None), sup)
        # covers_admission: the hold reaches not_after_ms + the skew allowance
        skew = A['capabilities'][1]['attempts']['clock_skew_allowance_ms']
        self.assertGreaterEqual(sup.suppressed_until_ms, NAF + skew)
        far = self.cancel('cancel_unknown_beyond_max', 'rec-far').suppression
        self.assertEqual((False, 'beyond_tombstone_max'),
                         (far.covers_admission, far.covers_admission_reason))
        plain = self.cancel('cancel_unknown_plain', 'rec-2').suppression
        self.assertEqual((False, 'no_not_after_ms', None),
                         (plain.covers_admission, plain.covers_admission_reason,
                          plain.not_after_ms))
        self.routes[('GET', '/traverse/attempt/rec-1')] = A['get_tombstoned']
        got = self.c.attempt('rec-1')
        self.assertEqual((404, 'attempt', True, NAF, True),
                         (got.status, got.route, got.tombstone, got.suppression.not_after_ms,
                          got.suppression.covers_admission))
        self.assertIsNone(got.attempt)
        self.routes[('GET', '/traverse/attempt/rec-never')] = A['get_unknown']
        got = self.c.attempt('rec-never')
        self.assertEqual((404, None, None), (got.status, got.tombstone, got.suppression))

    def test_429_reasons(self):
        full = self.cancel('cancel_tombstones_full', 'rec-7', not_after_ms=NAF)
        self.assertEqual((429, False, 'tombstones_full', True, None),
                         (full.status, full.tombstone, full.reason, full.retryable,
                          full.suppression))
        none = self.cancel('cancel_no_suppression', 'rec-8', not_after_ms=NAF)
        self.assertEqual((429, False, 'no_suppression', False),
                         (none.status, none.tombstone, none.reason, none.retryable))

    def test_finished_answers(self):
        got = self.cancel('cancel_finished', 'rec-5', not_after_ms=NAF)
        self.assertEqual((404, 'finished', True, 'completed'),
                         (got.status, got.state, got.finished, got.attempt['reason']))
        self.routes[('GET', '/traverse/attempt/rec-5')] = A['get_finished']
        got = self.c.attempt('rec-5')
        self.assertEqual((200, True, 'completed'), (got.status, got.finished, got.reason))
        self.assertEqual(NAF, got.attempt['not_after_ms'])

    def traverse_error(self, name, **kw):
        self.routes[('POST', '/traverse')] = A[name]
        with self.assertRaises(TraverseError) as cm:
            self.c.traverse([SEED], **kw)
        return cm.exception

    def test_the_409s(self):
        e = self.traverse_error('traverse_tombstoned', attempt_id='rec-1', not_after_ms=NAF,
                                expect_server_instance='auto')
        self.assertIsInstance(e, AttemptConflict)
        self.assertEqual((409, True, 'unknown', INST),
                         (e.status, e.tombstoned, e.state, e.server_instance))
        self.assertEqual((NAF, True),
                         (e.suppression.not_after_ms, e.suppression.covers_admission))
        self.assertEqual(AttemptSent('rec-1', NAF, INST), e.sent)
        self.assertIsNone(e.usage)
        e = self.traverse_error('traverse_tombstoned_without_not_after', attempt_id='rec-2')
        self.assertEqual((True, False, 'no_not_after_ms', None),
                         (e.tombstoned, e.suppression.covers_admission,
                          e.suppression.covers_admission_reason, e.suppression.not_after_ms))
        e = self.traverse_error('traverse_duplicate_finished', attempt_id='rec-5')
        self.assertIsInstance(e, AttemptConflict)
        self.assertEqual((False, 'finished', None), (e.tombstoned, e.state, e.suppression))

        self.c.server_instance()
        self.assertIsNotNone(self.c._server_doc)
        e = self.traverse_error('traverse_instance_mismatch', attempt_id='rec-3',
                                expect_server_instance='0123456789abcdef')
        self.assertIsInstance(e, InstanceMismatch)
        self.assertEqual(('0123456789abcdef', INST),
                         (e.expect_server_instance, e.server_instance))
        # the instance is read again before the next 'auto'
        self.assertIsNone(self.c._server_doc)
        n = _paths(self.s).count('/capabilities')
        self.c.build_request([SEED], attempt_id='x', expect_server_instance='auto')
        self.assertEqual(n + 1, _paths(self.s).count('/capabilities'))

        e = self.traverse_error('traverse_expired', attempt_id='rec-4', not_after_ms=1)
        self.assertIsInstance(e, AttemptExpired)
        self.assertEqual((1, INST), (e.not_after_ms, e.server_instance))
        self.assertGreater(e.server_time_ms, 1)
        for x in (AttemptConflict, InstanceMismatch, AttemptExpired):
            self.assertTrue(issubclass(x, TraverseError))


def _answer(name, route):
    status, body = A[name]
    return AttemptAnswer(body, status, route)


class TestReleaseRule(unittest.TestCase):
    SENT = AttemptSent('rec-1', NAF, INST)

    def check(self, answer, sent, code, release=False, early=False, **kw):
        v = release_verdict(answer, sent, **kw)
        self.assertIsInstance(v, ReleaseVerdict)
        self.assertEqual((release, early, code), (v.release, v.early, v.code), v.why)
        self.assertEqual(release, bool(v))
        self.assertTrue(v.why)
        return v

    def test_early_release_needs_every_condition(self):
        tomb = _answer('cancel_unknown_not_after', 'cancel')
        self.check(tomb, self.SENT, 'tombstone', True, True)
        # the request dict as sent works too
        self.check(tomb, {'attempt_id': 'rec-1', 'not_after_ms': NAF,
                          'expect_server_instance': INST, 'seeds': []},
                   'tombstone', True, True)
        self.check(_answer('get_tombstoned', 'attempt'), self.SENT, 'tombstone', True, True)
        self.check(tomb, AttemptSent('rec-1', NAF + 1, INST), 'other_not_after_ms')
        self.check(tomb, AttemptSent('rec-1', NAF, 'ffffffffffffffff'),
                   'other_server_instance')
        self.check(tomb, AttemptSent('rec-1', None, INST), 'sent_without_not_after_ms')
        self.check(tomb, AttemptSent('rec-1', NAF, None),
                   'sent_without_expect_server_instance')
        self.check(tomb, AttemptSent('rec-2', NAF, INST), 'answer_for_other_attempt')
        far = A['cancel_unknown_beyond_max'][1]['not_after_ms']
        v = self.check(_answer('cancel_unknown_beyond_max', 'cancel'),
                       AttemptSent('rec-far', far, INST), 'covers_admission_false')
        self.assertIn('beyond_tombstone_max', v.why)
        self.check(_answer('cancel_unknown_plain', 'cancel'), AttemptSent('rec-2', NAF, INST),
                   'covers_admission_false')
        # a level-3 tombstone states no suppression
        old = _answer('level3_cancel_unknown', 'cancel')
        self.check(old, AttemptSent('rec-10', NAF, old['server_instance']),
                   'no_suppression_stated')
        self.check(_answer('get_unknown', 'attempt'), AttemptSent('rec-never', NAF, INST),
                   'no_tombstone')

    def test_nothing_else_releases_early(self):
        for name in ('cancel_tombstones_full', 'cancel_no_suppression'):
            ans = _answer(name, 'cancel')
            v = self.check(ans, AttemptSent(ans.attempt_id, NAF, ans.server_instance),
                           'not_tombstoned')
            self.assertIn(ans.reason, v.why)
        stopping = AttemptAnswer({'attempt_id': 'rec-1', 'server_instance': INST,
                                  'cancelled': True, 'state': 'stopping'}, 200, 'cancel')
        self.check(stopping, self.SENT, 'stopping')
        running = AttemptAnswer({'attempt_id': 'rec-1', 'state': 'running'}, 200, 'attempt')
        self.check(running, self.SENT, 'running')
        # a 409 for a tombstoned id frees nothing; neither does an instance_mismatch (the
        # duplicate's finished state and an expired 409 release: TestReleaseGrounds409)
        for name, kind, sent in (
                ('traverse_tombstoned', AttemptConflict, self.SENT),
                ('traverse_instance_mismatch', InstanceMismatch,
                 AttemptSent('rec-3', NAF, '0123456789abcdef'))):
            status, body = A[name]
            self.check(kind(status, body['error'], body), sent, 'not_a_release_answer')
        # a 409 about another id (the duplicate of rec-5) says nothing about rec-1
        status, body = A['traverse_duplicate_finished']
        self.check(AttemptConflict(status, body['error'], body), self.SENT,
                   'answer_for_other_attempt')
        self.check(TraverseError(400, 'bad', {'error': 'bad'}), self.SENT,
                   'not_a_release_answer')
        self.check(None, self.SENT, 'no_answer')

    def test_finished_states_release(self):
        sent = AttemptSent('rec-5', NAF, INST)
        self.check(_answer('cancel_finished', 'cancel'), sent, 'finished', True)
        self.check(_answer('get_finished', 'attempt'), sent, 'finished', True)
        resp = TraverseClient.response(A['traverse_completed'][1])
        self.check(resp, sent, 'finished', True)
        self.check(resp, AttemptSent('rec-6', NAF, INST), 'answer_for_other_attempt')
        self.check(TraverseResponse({}, []), sent, 'answer_for_other_attempt')
        # an error after the attempt was registered is its response: it ran
        err = TraverseError(400, 'request.strategy', {'error': 'x',
                                                      'usage': {'attempt_id': 'rec-5'}})
        self.check(err, sent, 'finished', True)

    def test_the_clock(self):
        skew = A['capabilities'][1]['attempts']['clock_skew_allowance_ms']
        limit = NAF + skew + 40000
        kw = dict(clock_skew_allowance_ms=skew, bound_ms=40000)
        self.check(None, self.SENT, 'no_answer', now_ms=limit, **kw)
        self.check(None, self.SENT, 'clock', True, now_ms=limit + 1, **kw)
        v = self.check(None, self.SENT, 'clock', True, now_ms=limit + 1, **kw)
        # the clock takes the attempt as stopped at its bound apart from one uninterruptible
        # step, whose length has no stated bound: every clock release says so
        self.assertEqual(('uninterruptible_overrun',), v.assumptions)
        # an answer saying stopping (or running) past the limit is that step: held by
        # default (LRG-R4), released with the assumption stated only when asked to
        stopping = AttemptAnswer({'attempt_id': 'rec-1', 'state': 'stopping'}, 200, 'cancel')
        v = self.check(stopping, self.SENT, 'stopping', now_ms=limit + 1, **kw)
        self.assertIn('hold_clock_while_running', v.why)
        running = AttemptAnswer({'attempt_id': 'rec-1', 'state': 'running'}, 200, 'attempt')
        self.check(running, self.SENT, 'running', now_ms=limit + 10 ** 6, **kw)
        v = self.check(stopping, self.SENT, 'clock', True, now_ms=limit + 1,
                       hold_clock_while_running=False, **kw)
        self.assertEqual(('uninterruptible_overrun',), v.assumptions)
        self.assertIn('the last answer said stopping', v.why)
        # without not_after_ms the clock releases nothing
        self.check(None, AttemptSent('rec-1'), 'no_answer', now_ms=limit + 10 ** 9, **kw)

    def test_a_finished_release_states_its_assumptions(self):
        attempts = A['capabilities'][1]['attempts']
        sent = AttemptSent('rec-5', NAF, INST)
        resp = TraverseClient.response(A['traverse_completed'][1])
        for answer in (_answer('cancel_finished', 'cancel'), _answer('get_finished', 'attempt'),
                       resp):
            # held through not_after_ms + the skew on that process, a copy elsewhere refused
            v = self.check(answer, sent, 'finished', True, attempts=attempts)
            self.assertEqual((), v.assumptions)
            # without the attempts block the server's hold is not judged, and the verdict
            # says so rather than an empty list (() means "held")
            self.assertEqual(('server_hold_unchecked',),
                             self.check(answer, sent, 'finished', True).assumptions)
            # the whole document works as well as its attempts block
            self.assertEqual(v, release_verdict(answer, sent, attempts=A['capabilities'][1]))
        v = self.check(resp, AttemptSent('rec-5'), 'finished', True, attempts=attempts)
        self.assertEqual(('sent_without_not_after_ms', 'sent_without_expect_server_instance'),
                         v.assumptions)
        self.assertIn('no copy of the request arrives after the attempt left retention',
                      v.why)
        none = dict(attempts, retention_s=0)
        self.assertEqual(('no_retention',), self.check(resp, sent, 'finished', True,
                                                       attempts=none).assumptions)
        # not_after_ms + the skew later than tombstone_max_s after the finish
        finish = release_verdict.__globals__['_finish_ms'](resp.usage)
        far = finish + 1000 * attempts['tombstone_max_s'] - attempts['clock_skew_allowance_ms']

        def sent_with(naf):
            # the response of a copy sent with |naf|: its usage echoes it
            body = copy.deepcopy(A['traverse_completed'][1])
            body['usage']['not_after_ms'] = naf
            return TraverseClient.response(body)
        self.assertEqual((), self.check(sent_with(far), AttemptSent('rec-5', far, INST),
                                        'finished', True, attempts=attempts).assumptions)
        self.assertEqual(('beyond_tombstone_max',), self.check(
            sent_with(far + 1), AttemptSent('rec-5', far + 1, INST), 'finished', True,
            attempts=attempts).assumptions)
        # a finish the answer does not state is not taken to be covered
        bare = TraverseResponse({'usage': {'attempt_id': 'rec-5', 'server_instance': INST,
                                           'not_after_ms': NAF}}, [])
        v = self.check(bare, sent, 'finished', True, attempts=attempts)
        self.assertEqual(('beyond_tombstone_max',), v.assumptions)
        self.assertIn('the finish is not stated', v.why)
        # once the clock has passed too, its release (which assumes only that the
        # attempt's run past its bound has ended) is given
        skew, cap = attempts['clock_skew_allowance_ms'], attempts['hard_cap_ms']
        v = self.check(resp, AttemptSent('rec-5', NAF), 'clock', True,
                       now_ms=NAF + skew + cap + 1, attempts=attempts)
        self.assertEqual(('uninterruptible_overrun',), v.assumptions)
        # the attempts block gives the skew allowance and, as the bound, hard_cap_ms
        self.check(None, sent, 'no_answer', now_ms=NAF + skew + cap, attempts=attempts)
        self.check(None, sent, 'clock', True, now_ms=NAF + skew + cap + 1, attempts=attempts)

    def test_bad_arguments(self):
        with self.assertRaises(ValueError):
            release_verdict(None, {'seeds': []})
        with self.assertRaises(TypeError):
            release_verdict({'state': 'finished'}, self.SENT)


if __name__ == '__main__':
    unittest.main()
