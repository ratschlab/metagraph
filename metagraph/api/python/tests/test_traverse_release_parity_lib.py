"""The search service's release-parity cases, ported (its report RELEASE_PARITY_2026-10-05,
tests/test_traversal_release_parity.py at bbca7869, read-only; nothing of that repository is
used here): the same bytes, read by this library's own client and judged by
release_verdict(), with the library's verdict under the release-gap rules the report names
(LRG-G1..G7, LRG-R1..R4) and, per case, whether it agrees with the service's decision or
remains a documented divergence -- a service contract (the service's own stricter switches
and history) or a possible service bug the service's report names.

The bodies are those of the service's level-5 fake (answering as MetaGraph dcc0cebd, and as
f667d775 for level 4) on its fixed clock -- the server's wall clock half a second after the
grant G -- and the one-condition flips of its attempt-wire tests, written out here from the
fake's rules (_fake_*). Times: N = the sent not_after_ms (G + 10 s), bound_ms 40 s, the
stated skew 2000 ms.

Also here: the readers' gaps the service pinned (LRG-G1..G6) and the light import
(LRG-G7), checked in a subprocess through sys.modules."""

import copy
import datetime
import json
import os
import subprocess
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402

from metagraph.traverse import (  # noqa: E402
    AttemptAnswer, AttemptConflict, AttemptExpired, AttemptSent, InstanceMismatch,
    Suppression, TraverseClient, TraverseError, classify_409, release_verdict,
)

G = 1_791_172_800_000            # the grant instant (the worker's wall clock, epoch ms)
N = G + 10_000                   # the sent not_after_ms
SERVER_NOW = G + 500             # the fake server's wall clock
BOUND_MS = 40_000
SKEW = 2000
INST, OTHER = '0123456789abcdef', 'fedcba9876543210'
A = '66666666-7777-4888-9999-aaaaaaaaaaaa:1'
SENT5 = AttemptSent(A, N, INST)
NONE = AttemptSent(A)
# when a decision is taken: right after the exchange; between the library's clock limit
# (N + skew + bound) and the service sweep's (N + 5 s + bound); past both
NOW = G + 1_000
BETWEEN = N + 2_000 + BOUND_MS + 1_000
LATE = N + 5_000 + BOUND_MS + 1_000
RETENTION_MS = 3_600_000         # the fake's retention_s 3600
HELD_TO = SERVER_NOW + RETENTION_MS


def _iso(ms):
    """An epoch-ms instant as the server writes it."""
    t = datetime.datetime.fromtimestamp(ms / 1000, tz=datetime.timezone.utc)
    return t.isoformat(timespec='milliseconds').replace('+00:00', 'Z')


with open(os.path.join(T.DATA, 'level5', 'answers.json'), encoding='utf-8') as _f:
    _L5 = json.load(_f)['answers']['capabilities'][1]['attempts']
with open(os.path.join(T.DATA, 'level5', 'level4_f667d775.json'), encoding='utf-8') as _f:
    _L4 = json.load(_f)['attempts']
# the attempts blocks the fake states (the canonical documents with its instance, retention,
# skew and tombstone cap): level 5 states tombstone_max_s, level 4 does not
CAPS = {'level5': dict(_L5, server_instance=INST, retention_s=3600, retention_count=10_000,
                       clock_skew_allowance_ms=SKEW, tombstone_max_s=86_400),
        'level4': dict(_L4, server_instance=INST, retention_s=3600, retention_count=10_000,
                       clock_skew_allowance_ms=SKEW),
        'none': None}
LEVEL = {'level5': 5, 'level4': 4, 'none': 0}


# ------------------------------------------------------------------ the fake's bodies

def _unknown(instance=INST):
    return {'error': "unknown attempt_id '%s'" % A, 'attempt_id': A, 'state': 'unknown',
            'server_instance': instance}


def _usage():
    """The fake's usage of the finished attempt (one seed, a 30 s budget, not_after_ms N)."""
    return {'attempt_id': A, 'budget_id': None, 'locus_id': None, 'server_instance': INST,
            'reason': 'completed', 'received_at': _iso(SERVER_NOW), 'elapsed_ms': 1,
            'bound_ms': BOUND_MS, 'bound': {'allowance_ms': 10_000, 'time_budget_ms': 30_000.0,
                                            'enforced': True, 'capped_by': None},
            'seeds': {'requested': 1, 'started': 1, 'finished': 1, 'abandoned': 0},
            'work_units': 1000,
            'memory': {'peak_admitted_bytes': 1 << 20, 'soft_excess_bytes': 0,
                       'held_bound_bytes': None},
            'per_seed': [{'index': 0, 'work_units': 1000, 'elapsed_ms': 1, 'refused_bytes': 0}],
            'not_after_ms': N}


def _state(state, usage=None):
    """The fake's state_json of the attempt sent as SENT5."""
    return {'attempt_id': A, 'server_instance': INST, 'not_after_ms': N, 'state': state,
            'reason': 'completed' if state == 'finished' else None,
            'received_at': _iso(SERVER_NOW), 'usage': usage}


def _tomb_state(**suppression):
    return dict({'attempt_id': A, 'server_instance': INST, 'state': 'unknown',
                 'tombstone': True, 'cancel_requested_at': _iso(SERVER_NOW)}, **suppression)


COVERED_KEYS = {'suppressed_until_ms': HELD_TO, 'not_after_ms': N, 'covers_admission': True}

CANCEL = {
    'tomb_covered': (404, dict(_unknown(), cancelled=False, tombstone=True, **COVERED_KEYS)),
    'tomb_no_not_after': (404, dict(_unknown(), cancelled=False, tombstone=True,
                                    suppressed_until_ms=HELD_TO, covers_admission=False,
                                    covers_admission_reason='no_not_after_ms')),
    # retention_s 1, tombstone_max_s 1: the hold cannot reach not_after_ms + skew
    'tomb_beyond_max': (404, dict(_unknown(), cancelled=False, tombstone=True,
                                  suppressed_until_ms=SERVER_NOW + 1000, not_after_ms=N,
                                  covers_admission=False,
                                  covers_admission_reason='beyond_tombstone_max')),
    'tomb_other_instance': (404, dict(_unknown(OTHER), cancelled=False, tombstone=True,
                                      **COVERED_KEYS)),
    'tomb_level4': (404, dict(_unknown(), cancelled=False, tombstone=True)),
    'finished_404': (404, {'attempt_id': A, 'server_instance': INST,
                           'error': "attempt '%s' has finished: nothing to cancel" % A,
                           'state': 'finished', 'cancelled': False,
                           'attempt': _state('finished', _usage())}),
    'stopping_200': (200, {'attempt_id': A, 'server_instance': INST, 'cancelled': True,
                           'state': 'stopping', 'attempt': _state('stopping')}),
    'finished_200': (200, {'attempt_id': A, 'server_instance': INST, 'cancelled': True,
                           'state': 'finished',
                           'attempt': {'attempt_id': A, 'state': 'finished',
                                       'reason': 'cancelled',
                                       'usage': {'attempt_id': A, 'work_units': 5}}}),
    'full_429': (429, dict(_unknown(), cancelled=False, tombstone=False,
                           reason='tombstones_full')),
    'no_suppression_429': (429, dict(_unknown(), cancelled=False, tombstone=False,
                                     reason='no_suppression')),
    'full_429_level4': (429, dict(_unknown(), cancelled=False, tombstone=False)),
    'unknown_404': (404, _unknown()),
}
COVERED = CANCEL['tomb_covered']


def _with(reply, **over):
    """|reply| with keys of its body replaced (a value None drops the key)."""
    status, body = reply
    out = dict(body, **over)
    return status, {k: v for k, v in out.items() if v is not None}


FLIPS = {
    'covers_missing': _with(COVERED, covers_admission=None),
    'covers_false': _with(COVERED, covers_admission=False,
                          covers_admission_reason='beyond_tombstone_max'),
    'covers_string': _with(COVERED, covers_admission='true'),
    'instance_missing': _with(COVERED, server_instance=None),
    'not_after_other': _with(COVERED, not_after_ms=N + 1),
    'not_after_bool': _with(COVERED, not_after_ms=True),
    'not_after_missing': _with(COVERED, not_after_ms=None),
    'suppression_short': _with(COVERED, suppressed_until_ms=N + 1_999),
    'suppressed_until_missing': _with(COVERED, suppressed_until_ms=None),
    'state_finished': _with(COVERED, state='finished'),
    'tombstone_false': _with(COVERED, tombstone=False),
    'attempt_id_other': _with(COVERED, attempt_id='other:1'),
}
_TOMBED = "attempt_id '%s' was cancelled before this request arrived (tombstoned)" % A
_USED = "attempt_id '%s' is running or was used on this server: an attempt runs once" % A
_finished_409 = (409, {'error': _USED, 'attempt': _state('finished', _usage())})
DISPATCH = {  # (reply, sent, caps)
    'tomb_covered': ((409, {'error': _TOMBED, 'attempt': _tomb_state(**COVERED_KEYS)}),
                     SENT5, 'level5'),
    'tomb_unsent': ((409, {'error': _TOMBED, 'attempt': _tomb_state(
        suppressed_until_ms=HELD_TO, covers_admission=False,
        covers_admission_reason='no_not_after_ms')}), NONE, 'level5'),
    'tomb_level4': ((409, {'error': _TOMBED, 'attempt': _tomb_state()}), NONE, 'level4'),
    'running': ((409, {'error': _USED, 'attempt': _state('running')}), SENT5, 'level5'),
    'finished': (_finished_409, SENT5, 'level5'),
    'finished_other_id': ((409, dict(_finished_409[1], attempt=dict(
        _finished_409[1]['attempt'], attempt_id='other:1'))), SENT5, 'level5'),
    # the server's clock 60 s ahead: past not_after_ms on arrival
    'expired': ((409, {'error': 'not_after_ms %d has passed on this server\'s clock (%d)'
                       % (N, SERVER_NOW + 60_000), 'state': 'expired', 'not_after_ms': N,
                       'server_time_ms': SERVER_NOW + 60_000, 'attempt_id': A,
                       'server_instance': INST}), SENT5, 'level5'),
    # restarted as OTHER
    'instance_mismatch': ((409, {'error': "expect_server_instance '%s' is not this server's "
                                 "instance" % INST, 'state': 'instance_mismatch',
                                 'expect_server_instance': INST, 'server_instance': OTHER,
                                 'attempt_id': A, 'not_after_ms': N}), SENT5, 'level5'),
    # a request that named no instance
    'instance_mismatch_unsent': ((409, {'error': 'expect_server_instance mismatch',
                                        'state': 'instance_mismatch',
                                        'expect_server_instance': INST,
                                        'server_instance': OTHER}), NONE, 'level4'),
    'tomb_other_id': ((409, {'error': 'tombstoned', 'attempt': {
        'attempt_id': 'other:1', 'state': 'unknown', 'tombstone': True,
        'server_instance': INST, 'covers_admission': True}}), SENT5, 'level5'),
}
RESPONSE = (200, {'results': [], 'capabilities': {'feature_level': 5}, 'usage': _usage()})


# ------------------------------------------------------------------ the cases

class Case:
    def __init__(self, id, site, reply, sent=SENT5, caps='level5', now_ms=NOW):
        self.id, self.site, self.reply, self.sent = id, site, reply, sent
        self.caps, self.now_ms = caps, now_ms


def _cases():
    cases = []
    for name, reply in CANCEL.items():
        cases += [Case('worker.' + name, 'worker', reply), Case('sweep.' + name, 'sweep', reply)]
    # the service's switch TRAVERSAL_TOMBSTONE_EARLY_RELEASE on: the same bytes for the
    # library, which has no such switch
    cases += [Case('worker.tomb_covered.switch_on', 'worker', COVERED),
              Case('sweep.tomb_covered.switch_on', 'sweep', COVERED)]
    cases += [Case('worker.flip.' + name, 'worker', reply) for name, reply in FLIPS.items()]
    cases += [
        Case('worker.sent_none', 'worker', COVERED, sent=NONE),
        Case('worker.sent_not_after_only', 'worker', COVERED, sent=AttemptSent(A, N)),
        Case('worker.sent_expect_only', 'worker', COVERED, sent=AttemptSent(A, None, INST)),
        Case('sweep.sent_none', 'sweep', COVERED, sent=NONE),
        # the service's instance history (two processes seen): nothing the library keeps
        Case('worker.two_instances', 'worker', COVERED),
        Case('sweep.two_instances', 'sweep', COVERED),
        Case('sweep.caps_level4', 'sweep', COVERED, caps='level4'),
        Case('sweep.caps_unknown', 'sweep', COVERED, caps='none'),
        Case('sweep.finished_404_other_id', 'sweep',
             _with(CANCEL['finished_404'], attempt_id='other:1')),
    ]
    cases += [
        Case('sweep.tomb_no_not_after.between', 'sweep', CANCEL['tomb_no_not_after'],
             now_ms=BETWEEN),
        Case('sweep.tomb_no_not_after.late', 'sweep', CANCEL['tomb_no_not_after'],
             now_ms=LATE),
        Case('sweep.full_429.late', 'sweep', CANCEL['full_429'], now_ms=LATE),
        Case('sweep.stopping_200.late', 'sweep', CANCEL['stopping_200'], now_ms=LATE),
        Case('sweep.tomb_covered.late', 'sweep', COVERED, now_ms=LATE),
        Case('sweep.tomb_covered.late.switch_on', 'sweep', COVERED, now_ms=LATE),
    ]
    for name, (reply, sent, caps) in DISPATCH.items():
        cases.append(Case('dispatch.' + name, 'dispatch', reply, sent=sent, caps=caps))
    reply, sent, caps = DISPATCH['instance_mismatch']
    cases.append(Case('dispatch.instance_mismatch.switch_on', 'dispatch', reply, sent=sent,
                      caps=caps))
    cases += [Case('response.usage', 'response', RESPONSE),
              Case('response.usage_other_attempt', 'response',
                   (200, dict(RESPONSE[1], usage=dict(RESPONSE[1]['usage'],
                                                      attempt_id='other:1'))))]
    return cases


CASES = _cases()

# The library's verdict: (release, early, code). Every case not listed holds (False, False)
# with the code given in HELD.
RELEASED = {
    'worker.tomb_covered': (True, True, 'tombstone'),
    'sweep.tomb_covered': (True, True, 'tombstone'),
    'worker.tomb_covered.switch_on': (True, True, 'tombstone'),
    'sweep.tomb_covered.switch_on': (True, True, 'tombstone'),
    'worker.finished_404': (True, False, 'finished'),
    'sweep.finished_404': (True, False, 'finished'),
    'worker.finished_200': (True, False, 'finished'),
    'sweep.finished_200': (True, False, 'finished'),
    'worker.flip.state_finished': (True, False, 'finished'),
    'worker.two_instances': (True, True, 'tombstone'),
    'sweep.two_instances': (True, True, 'tombstone'),
    'sweep.caps_level4': (True, True, 'tombstone'),
    'sweep.caps_unknown': (True, True, 'tombstone'),
    'sweep.tomb_no_not_after.between': (True, False, 'clock'),
    'sweep.tomb_no_not_after.late': (True, False, 'clock'),
    'sweep.full_429.late': (True, False, 'clock'),
    'sweep.tomb_covered.late': (True, True, 'tombstone'),
    'sweep.tomb_covered.late.switch_on': (True, True, 'tombstone'),
    'dispatch.finished': (True, False, 'finished'),           # LRG-R1
    'dispatch.expired': (True, False, 'expired'),              # LRG-R2
    'response.usage': (True, False, 'finished'),
}
HELD = {
    'tomb_no_not_after': 'covers_admission_false', 'tomb_beyond_max': 'covers_admission_false',
    'tomb_other_instance': 'other_server_instance', 'tomb_level4': 'no_suppression_stated',
    'stopping_200': 'stopping', 'full_429': 'not_tombstoned',
    'no_suppression_429': 'not_tombstoned', 'full_429_level4': 'not_tombstoned',
    'unknown_404': 'no_tombstone',
    'flip.covers_missing': 'covers_admission_unstated',            # LRG-G5
    'flip.covers_false': 'covers_admission_false',
    'flip.covers_string': 'covers_admission_unstated',             # LRG-G5
    'flip.instance_missing': 'other_server_instance',
    'flip.not_after_other': 'other_not_after_ms',
    'flip.not_after_bool': 'other_not_after_ms',                   # LRG-G5: a bool is no instant
    'flip.not_after_missing': 'other_not_after_ms',
    'flip.suppression_short': 'suppression_short',                 # LRG-R3
    'flip.suppressed_until_missing': 'suppression_short',          # LRG-G5 + R3
    'flip.tombstone_false': 'no_tombstone',
    'flip.attempt_id_other': 'answer_for_other_attempt',
    'sent_none': 'sent_without_not_after_ms',
    'sent_not_after_only': 'sent_without_expect_server_instance',
    'sent_expect_only': 'sent_without_not_after_ms',
    'finished_404_other_id': 'answer_for_other_attempt',
    'stopping_200.late': 'stopping',                               # LRG-R4 / X2
    'dispatch.tomb_covered': 'not_a_release_answer',
    'dispatch.tomb_unsent': 'not_a_release_answer',
    'dispatch.tomb_level4': 'not_a_release_answer',
    'dispatch.running': 'running',
    'dispatch.finished_other_id': 'answer_for_other_attempt',
    'dispatch.instance_mismatch': 'not_a_release_answer',
    'dispatch.instance_mismatch_unsent': 'not_a_release_answer',   # LRG-G2
    'dispatch.tomb_other_id': 'answer_for_other_attempt',
    'dispatch.instance_mismatch.switch_on': 'not_a_release_answer',
    'response.usage_other_attempt': 'answer_for_other_attempt',
}

# Where the service decided otherwise at bbca7869 (its DIVERGENCES table): (the service's
# (release, early), classification, why). Every other case agrees.
SERVICE_DIVERGES = {
    'worker.tomb_covered': ((False, False), 'service contract',
                            'D11: the service\'s early release is off by default'),
    'sweep.tomb_covered': ((False, False), 'service contract', 'D11'),
    'sweep.tomb_covered.late': ((True, False), 'service contract',
                                'D11: the sweep releases by its clock fallback, not early'),
    'worker.finished_404': ((False, False), 'service contract',
                            'D7: the follow-up cancel is judged by the tombstone rule only'),
    'worker.finished_200': ((False, False), 'service contract', 'D7'),
    'worker.flip.state_finished': ((False, False), 'service contract', 'D7'),
    'worker.flip.attempt_id_other': ((True, True), 'possible service bug',
                                     'the service\'s early_release never compares the '
                                     'answer\'s attempt_id'),
    'worker.two_instances': ((False, False), 'service contract',
                             'critique A2: two processes seen behind the endpoint'),
    'sweep.two_instances': ((False, False), 'service contract', 'critique A2'),
    'sweep.caps_level4': ((False, False), 'service contract',
                          'host_not_level5: the host re-read at level 4'),
    'sweep.caps_unknown': ((False, False), 'service contract',
                           'host_not_level5: the capabilities unreadable'),
    'sweep.finished_404_other_id': ((True, False), 'possible service bug',
                                    'classify_cancel\'s finished verdict never compares the '
                                    'attempt_id'),
    'sweep.tomb_no_not_after.between': ((False, False), 'service contract',
                                        'the sweep\'s clock fallback adds at least 5 s '
                                        '(TRAVERSAL_NOT_AFTER_MAX_SKEW_S), not the stated 2 s'),
    'dispatch.finished_other_id': ((True, False), 'possible service bug',
                                   '_on_409 does not compare the finished attempt\'s id'),
    'dispatch.instance_mismatch.switch_on': ((True, False), 'service contract',
                                             'D3: TRAVERSAL_INSTANCE_MISMATCH_RELEASE on'),
    'response.usage_other_attempt': ((True, False), 'service contract',
                                     'critique C2: an echo mismatch is a canary; the '
                                     'received response releases'),
}
# Divergences of the service's table at bbca7869 on which the library decides as the service
NOW_AGREE = {
    'worker.flip.suppression_short': 'LRG-R3: suppressed_until_ms >= not_after_ms + skew',
    'sweep.stopping_200.late': 'LRG-R4 / X2: no clock release against running or stopping',
    'dispatch.finished': 'LRG-R1: a 409\'s finished attempt object is a finished state',
    'dispatch.expired': 'LRG-R2: an expired 409 releases (code expired)',
}


class _Wire:
    """A requests-like session: every request answered with one recorded body."""

    def __init__(self, reply):
        self.status, self.body = reply
        self.sent = []

    def request(self, method, url, data=None, headers=None, timeout=None):
        self.sent.append((method, url, json.loads(data) if data else None))

        class R:
            pass
        r = R()
        r.status_code, r.content, r.headers = self.status, json.dumps(self.body).encode(), {}
        return r


def _request(sent):
    req = {'attempt_id': sent.attempt_id, 'seeds': [{'sequence': 'ACGTACGTAC'}],
           'strategy': {'bounds': {'time_budget_ms': 30_000}}}
    if sent.not_after_ms is not None:
        req['not_after_ms'] = sent.not_after_ms
    if sent.expect_server_instance is not None:
        req['expect_server_instance'] = sent.expect_server_instance
    return req


def library_verdict(case):
    """The library's ReleaseVerdict on the case's bytes, read by its own client."""
    client = TraverseClient('library.invalid', 1, session=_Wire(copy.deepcopy(case.reply)),
                            feature_level=LEVEL[case.caps])
    if case.site == 'worker':        # the worker's follow-up cancel names not_after_ms
        answer = client.cancel(A, wait_ms=0, not_after_ms=case.sent.not_after_ms)
    elif case.site == 'sweep':       # the sweep's cancel never does
        answer = client.cancel(A, wait_ms=0)
    elif case.site == 'dispatch':
        try:
            client.traverse_raw(_request(case.sent))
        except TraverseError as e:
            answer = e
        else:
            raise AssertionError('%s: a 409 did not raise' % case.id)
    else:
        answer = TraverseClient.response(copy.deepcopy(case.reply[1]))
    return release_verdict(answer, case.sent, now_ms=case.now_ms, bound_ms=BOUND_MS,
                           attempts=CAPS[case.caps])


def _expected(case):
    if case.id in RELEASED:
        return RELEASED[case.id]
    key = case.id.split('.', 1)[1] if case.site in ('worker', 'sweep') else case.id
    return (False, False, HELD[key])


class TestReleaseParity(unittest.TestCase):
    def test_the_case_table(self):
        ids = [c.id for c in CASES]
        self.assertEqual(66, len(ids))                  # the service's 66 cases
        self.assertEqual(len(ids), len(set(ids)))
        self.assertLessEqual(set(SERVICE_DIVERGES) | set(NOW_AGREE) | set(RELEASED),
                             set(ids))
        self.assertFalse(set(SERVICE_DIVERGES) & set(NOW_AGREE))

    def test_the_librarys_verdict_per_case(self):
        for case in CASES:
            with self.subTest(case.id):
                v = library_verdict(case)
                self.assertEqual(_expected(case), (v.release, v.early, v.code), v.why)

    def test_agreement_with_the_service(self):
        """A case agrees with the service unless it is a documented divergence; a listed
        one still diverges (else it belongs to NOW_AGREE)."""
        for case in CASES:
            with self.subTest(case.id):
                v = library_verdict(case)
                mine = (v.release, v.early)
                if case.id in SERVICE_DIVERGES:
                    service, cls, _ = SERVICE_DIVERGES[case.id]
                    self.assertNotEqual(service, mine)
                    self.assertIn(cls, ('service contract', 'possible service bug'))
                elif case.id in NOW_AGREE:
                    # the service's decision at bbca7869 is the library's
                    service = {'worker.flip.suppression_short': (False, False),
                               'sweep.stopping_200.late': (False, False),
                               'dispatch.finished': (True, False),
                               'dispatch.expired': (True, False)}[case.id]
                    self.assertEqual(service, mine)
        self.assertEqual(16, len(SERVICE_DIVERGES))     # the table's 20 less NOW_AGREE's 4

    def test_the_assumptions_the_releases_state(self):
        by_id = {c.id: c for c in CASES}
        # the fake's finished attempt: pinned, its own not_after_ms, held through N + skew
        self.assertEqual((), library_verdict(by_id['dispatch.finished']).assumptions)
        self.assertEqual((), library_verdict(by_id['response.usage']).assumptions)
        self.assertEqual((), library_verdict(by_id['dispatch.expired']).assumptions)
        # the hand-written 200 states neither the copy's not_after_ms nor its finish
        self.assertEqual(('other_copy_not_held', 'beyond_tombstone_max'),
                         library_verdict(by_id['worker.finished_200']).assumptions)
        for name in ('sweep.tomb_no_not_after.between', 'sweep.full_429.late'):
            self.assertEqual(('uninterruptible_overrun',),
                             library_verdict(by_id[name]).assumptions)

    def test_the_stopping_reply_releases_by_the_clock_only_when_asked(self):
        case = next(c for c in CASES if c.id == 'sweep.stopping_200.late')
        client = TraverseClient('x', 1, session=_Wire(copy.deepcopy(case.reply)),
                                feature_level=5)
        answer = client.cancel(A, wait_ms=0)
        v = release_verdict(answer, SENT5, now_ms=LATE, bound_ms=BOUND_MS,
                            attempts=CAPS['level5'], hold_clock_while_running=False)
        self.assertEqual((True, False, 'clock', ('uninterruptible_overrun',)),
                         (v.release, v.early, v.code, v.assumptions))


# ------------------------------------------------------------------ the readers (LRG-G1..G6)

class TestReaders(unittest.TestCase):
    def test_g1_a_non_object_attempt_reads_as_absent(self):
        for body in ({'attempt': 'not an object'}, {'attempt': 7}, {'attempt': None}):
            e = AttemptConflict(409, 'x', body)
            self.assertEqual(({}, None, False, None, None, None),
                             (e.attempt, e.state, e.tombstoned, e.server_instance,
                              e.suppression, e.attempt_id))
        for body in (None, [], 'text', 42, {'attempt': 'not an object'}):
            e = classify_409(body)
            self.assertIs(type(e), TraverseError)
            self.assertEqual(409, e.status)

    def test_g2_the_409_classification(self):
        conflict = {'error': 'x', 'attempt': {'attempt_id': A, 'state': 'running'}}
        self.assertIs(type(classify_409(conflict)), AttemptConflict)
        # the attempt object's state first, whatever the top level says
        self.assertIs(type(classify_409(dict(conflict, state='expired'))), AttemptConflict)
        expired = DISPATCH['expired'][0][1]
        self.assertIs(type(classify_409(expired)), AttemptExpired)
        mismatch = DISPATCH['instance_mismatch'][0][1]
        # an instance_mismatch only to a request that named an instance
        self.assertIs(type(classify_409(mismatch, SENT5)), InstanceMismatch)
        self.assertIs(type(classify_409(mismatch, _request(SENT5))), InstanceMismatch)
        for sent in (None, NONE, AttemptSent(A, N), {'seeds': []}):
            self.assertIs(type(classify_409(mismatch, sent)), TraverseError)
        # an attempt object without a state is still the id's
        self.assertIs(type(classify_409({'attempt': {'attempt_id': A}})), AttemptConflict)
        self.assertIs(type(classify_409({'error': 'busy'})), TraverseError)
        self.assertEqual('busy', classify_409({'error': 'busy'}).message)
        # the client uses it: an instance_mismatch to a request without one is not typed
        c = TraverseClient('h', 1, session=_Wire((409, mismatch)), feature_level=5)
        with self.assertRaises(TraverseError) as cm:
            c.traverse_raw(_request(NONE))
        self.assertIs(type(cm.exception), TraverseError)
        with self.assertRaises(InstanceMismatch):
            c.traverse_raw(_request(SENT5))

    def test_g3_the_conflicts_fields_fall_back_to_the_top_level(self):
        body = {'attempt': {'attempt_id': A, 'state': 'finished'}, 'server_instance': INST}
        self.assertEqual(INST, AttemptConflict(409, 'x', body).server_instance)
        body = {'attempt': {'attempt_id': A, 'state': 'unknown'}, 'tombstone': True}
        self.assertTrue(AttemptConflict(409, 'x', body).tombstoned)
        # a null in the object is absent (the duplicate 409 writes "tombstone": null)
        body = {'attempt': {'attempt_id': A, 'state': 'unknown', 'tombstone': None,
                            'server_instance': None}, 'tombstone': True,
                'server_instance': INST}
        e = AttemptConflict(409, 'x', body)
        self.assertEqual((True, INST), (e.tombstoned, e.server_instance))
        dup = {'attempt': dict(_state('finished', _usage()), tombstone=None)}
        self.assertFalse(AttemptConflict(409, 'x', dup).tombstoned)
        # the object wins where it states one
        body = {'attempt': {'attempt_id': A, 'state': 'unknown', 'tombstone': False},
                'tombstone': True}
        self.assertFalse(AttemptConflict(409, 'x', body).tombstoned)

    def test_g4_the_refused_ids_attempt_id_and_usage(self):
        e = classify_409(DISPATCH['finished'][0][1])
        self.assertEqual((A, 1000), (e.attempt_id, e.usage['work_units']))
        e = classify_409(DISPATCH['tomb_other_id'][0][1])
        self.assertEqual(('other:1', None), (e.attempt_id, e.usage))
        # a refused copy whose id is running: no usage yet
        self.assertIsNone(classify_409(DISPATCH['running'][0][1]).usage)

    def test_g5_suppression_typed_and_tolerant(self):
        body = dict(COVERED[1])
        self.assertEqual(Suppression(HELD_TO, N, True, None), Suppression.of(body))
        # an absent suppressed_until_ms is None, not a missing suppression
        sup = Suppression.of({k: v for k, v in body.items() if k != 'suppressed_until_ms'})
        self.assertEqual((None, N, True), (sup.suppressed_until_ms, sup.not_after_ms,
                                           sup.covers_admission))
        # covers_admission True only when stated true; absent or not a bool: None
        self.assertIsNone(Suppression.of(dict(body, covers_admission='true')).covers_admission)
        unstated = {k: v for k, v in body.items() if k != 'covers_admission'}
        self.assertIsNone(Suppression.of(unstated).covers_admission)
        self.assertIs(False, Suppression.of(dict(body, covers_admission=False))
                      .covers_admission)
        # instants: integers of the wire's range only
        for bad in (True, False, 1.0e12, float(N), -1, 1 << 53, str(N)):
            got = Suppression.of(dict(body, not_after_ms=bad, suppressed_until_ms=bad))
            self.assertEqual((None, None), (got.not_after_ms, got.suppressed_until_ms), bad)
        # no suppression key at all: none (a level-4 tombstone; a finished state's own
        # not_after_ms is not one)
        self.assertIsNone(Suppression.of(CANCEL['tomb_level4'][1]))
        self.assertIsNone(Suppression.of(_state('finished')))
        self.assertIsNone(Suppression.of(None))

    def test_g6_an_answer_read_from_its_attempt_object(self):
        # a cancel's answer stating its state only inside `attempt` (the service's batch
        # 6b): a block-only stopping keeps the clock off, a block-only finished is finished
        stopping = AttemptAnswer({'attempt': _state('stopping')}, 200, 'cancel')
        self.assertEqual(('stopping', INST, A), (stopping.state, stopping.server_instance,
                                                 stopping.attempt_id))
        v = release_verdict(stopping, SENT5, now_ms=LATE, bound_ms=BOUND_MS,
                            attempts=CAPS['level5'])
        self.assertEqual((False, 'stopping'), (v.release, v.code))
        done = AttemptAnswer({'attempt': _state('finished', _usage())}, 404, 'cancel')
        self.assertTrue(done.finished)
        self.assertEqual(1000, done.usage['work_units'])
        v = release_verdict(done, SENT5, attempts=CAPS['level5'])
        self.assertEqual((True, 'finished', ()), (v.release, v.code, v.assumptions))
        # GET's usage at the top level
        got = AttemptAnswer(dict(_state('finished', _usage())), 200, 'attempt')
        self.assertEqual(1000, got.usage['work_units'])
        self.assertIsNone(AttemptAnswer(CANCEL['tomb_level4'][1], 404, 'cancel').usage)
        # an answer that names no attempt_id holds
        anon = AttemptAnswer({k: v for k, v in COVERED[1].items() if k != 'attempt_id'},
                             404, 'cancel')
        v = release_verdict(anon, SENT5, attempts=CAPS['level5'])
        self.assertEqual((False, 'answer_for_other_attempt'), (v.release, v.code))
        self.assertIn('names no attempt_id', v.why)


# ------------------------------------------------------------------ the release points

class TestReleaseGrounds(unittest.TestCase):
    def verdict(self, answer, sent=SENT5, **kw):
        kw.setdefault('attempts', CAPS['level5'])
        return release_verdict(answer, sent, **kw)

    def test_r1_a_409s_finished_attempt_object(self):
        e = classify_409(DISPATCH['finished'][0][1])
        v = self.verdict(e)
        self.assertEqual((True, False, 'finished', ()),
                         (v.release, v.early, v.code, v.assumptions))
        # the same conditions as a cancel's finished answer: an unpinned or unbounded send
        # states what it assumes
        v = self.verdict(e, NONE)
        self.assertEqual(('sent_without_not_after_ms', 'sent_without_expect_server_instance'),
                         v.assumptions)
        # a finish stated by another process than the pinned one, or of another copy
        other = copy.deepcopy(DISPATCH['finished'][0][1])
        other['attempt']['server_instance'] = OTHER
        self.assertEqual(('other_server_instance',),
                         self.verdict(classify_409(other)).assumptions)
        other['attempt'].pop('server_instance')
        self.assertEqual(('server_instance_not_stated',),
                         self.verdict(classify_409(other)).assumptions)
        copy2 = copy.deepcopy(DISPATCH['finished'][0][1])
        copy2['attempt']['not_after_ms'] = N - 5000
        copy2['attempt']['usage']['not_after_ms'] = N - 5000
        self.assertEqual(('other_copy_not_held',),
                         self.verdict(classify_409(copy2)).assumptions)

    def test_o16_o20_the_servers_hold_judged_only_with_its_block(self):
        e = classify_409(DISPATCH['finished'][0][1])
        self.assertEqual(('server_hold_unchecked',),
                         release_verdict(e, SENT5).assumptions)
        self.assertEqual(('server_hold_unchecked',), self.verdict(
            e, attempts={'clock_skew_allowance_ms': SKEW}).assumptions)
        # a real level-4 block (f667d775): the id is held only while retained
        with open(os.path.join(T.DATA, 'level5', 'level4_f667d775.json'),
                  encoding='utf-8') as f:
            rec = json.load(f)
        self.assertNotIn('tombstone_max_s', rec['attempts'])
        status, body = rec['get_finished']
        got = AttemptAnswer(body, status, 'attempt')
        sent = AttemptSent('rec4-1', rec['not_after_ms'])
        v = release_verdict(got, sent, attempts=rec['attempts'])
        self.assertEqual(('no_finished_hold', 'sent_without_expect_server_instance'),
                         v.assumptions)
        self.assertIn('held only while the attempt is retained', v.why)
        # retention_s 0 at any level
        self.assertEqual(('no_retention',), self.verdict(
            e, attempts=dict(CAPS['level5'], retention_s=0)).assumptions)

    def test_r2_an_expired_409(self):
        e = classify_409(DISPATCH['expired'][0][1])
        v = self.verdict(e)
        self.assertEqual((True, False, 'expired', ()),
                         (v.release, v.early, v.code, v.assumptions))
        self.assertIn('no copy of the attempt was registered', v.why)
        v = self.verdict(e, AttemptSent(A, N))
        self.assertEqual(('sent_without_expect_server_instance',), v.assumptions)
        # not the answer to this send
        self.assertEqual('other_not_after_ms', self.verdict(e, AttemptSent(A, N + 1, INST)).code)
        self.assertEqual('other_server_instance', self.verdict(e, AttemptSent(A, N, OTHER)).code)
        self.assertEqual('not_a_release_answer', self.verdict(e, NONE).code)
        self.assertEqual('answer_for_other_attempt',
                         self.verdict(e, AttemptSent('b', N, INST)).code)

    def test_r2_the_server_clock_step_back(self):
        # Later copies are refused only while the server's clock reads later than
        # not_after_ms (a strict check, no allowance), and the server's release_rule assumes
        # its clock does not step back below it. Every release assumes a step back of at most
        # clock_skew_allowance_ms, so the expired release states the step of its own wherever
        # the 409's margin is within that, or cannot be judged
        body = DISPATCH['expired'][0][1]

        def at(margin):
            if margin is None:
                return classify_409({k: v for k, v in body.items() if k != 'server_time_ms'},
                                    SENT5)
            return classify_409(dict(body, server_time_ms=N + margin), SENT5)

        for margin, stated in ((1, True), (3, True), (SKEW, True), (SKEW + 1, False),
                               (60_000, False), (None, True)):
            with self.subTest(margin=margin):
                v = self.verdict(at(margin))
                self.assertEqual((True, False, 'expired'), (v.release, v.early, v.code))
                self.assertEqual(('server_clock_step_back',) if stated else (),
                                 v.assumptions, v.why)
                if stated:
                    self.assertIn("this server's clock does not step back below "
                                  'not_after_ms', v.why)
        self.assertIn('server_time_ms - not_after_ms (3 ms) is the step it survives, within '
                      'the clock_skew_allowance_ms (2000 ms)', self.verdict(at(3)).why)
        self.assertIn('states no server_time_ms', self.verdict(at(None)).why)
        # the call's allowance overrides the block's, as everywhere in release_verdict()
        self.assertEqual(('server_clock_step_back',), self.verdict(
            at(SKEW + 1), clock_skew_allowance_ms=SKEW + 1).assumptions)
        self.assertEqual((), self.verdict(at(SKEW), clock_skew_allowance_ms=SKEW - 1).assumptions)
        # no allowance at all (no block, none given): the margin cannot be judged
        v = release_verdict(at(60_000), SENT5)
        self.assertEqual(('server_clock_step_back',), v.assumptions)
        self.assertIn('no clock_skew_allowance_ms was given', v.why)
        # unpinned and just late: both assumptions, in that order
        e = classify_409(dict(body, server_time_ms=N + 3), AttemptSent(A, N))
        self.assertEqual(('server_clock_step_back', 'sent_without_expect_server_instance'),
                         self.verdict(e, AttemptSent(A, N)).assumptions)
        # a 1 ms margin against a 2,000 ms allowance
        n = 1_800_000_000_000
        sent = {'attempt_id': 'a1', 'not_after_ms': n,
                'expect_server_instance': 'abcdef0123456789'}
        e = classify_409({'error': 'x', 'state': 'expired', 'not_after_ms': n,
                          'server_time_ms': n + 1, 'attempt_id': 'a1',
                          'server_instance': 'abcdef0123456789'}, sent)
        v = release_verdict(e, sent, attempts={'clock_skew_allowance_ms': 2000,
                                               'retention_s': 3600, 'tombstone_max_s': 86400,
                                               'hard_cap_ms': 899000})
        self.assertEqual((True, False, 'expired', ('server_clock_step_back',)),
                         (v.release, v.early, v.code, v.assumptions))
        # once the clock has passed too, the clock release is returned (it rests only on
        # the bound's overrun and the shared step-back assumption)
        v = self.verdict(at(3), now_ms=N + SKEW + BOUND_MS + 1, bound_ms=BOUND_MS)
        self.assertEqual(('clock', ('uninterruptible_overrun',)), (v.code, v.assumptions))

    def test_the_clock_states_the_whole_overrun(self):
        # Past its bound an attempt runs on until its next delivery check -- the rest of its
        # piece, the walk up to a poll that reads the clock, the stopped seed's finalisation
        # and the building up to that check -- not "one uninterruptible step"; in the
        # server's words (its bound and release_rule)
        v = self.verdict(None, now_ms=N + SKEW + BOUND_MS + 1, bound_ms=BOUND_MS)
        self.assertEqual(('clock', ('uninterruptible_overrun',)), (v.code, v.assumptions))
        for phrase in ("the attempt's run past its bound has ended",
                       'runs on until its next delivery check (then the 503) or its '
                       "handler's return", 'the rest of the piece it was in when the bound '
                       'passed', 'the walk up to its next poll that reads the clock',
                       "the stopped seed's finalisation", 'up to the first delivery check',
                       'a run of no stated length'):
            self.assertIn(phrase, v.why)
        self.assertNotIn('one uninterruptible step', v.why)
        running = AttemptAnswer({'attempt_id': A, 'server_instance': INST,
                                 'state': 'running'}, 200, 'attempt')
        v = self.verdict(running, now_ms=N + SKEW + BOUND_MS + 1, bound_ms=BOUND_MS)
        self.assertEqual((False, 'running'), (v.release, v.code))
        self.assertIn('its run past its bound, up to its next delivery check, is still going',
                      v.why)
        with open(os.path.join(T.REPO, 'docs', 'source', 'graphlets.rst')) as f:
            rst = ' '.join(f.read().split())
        self.assertNotIn('one uninterruptible step', rst)
        self.assertIn('past it the attempt runs on until its next delivery check (then the '
                      "503) or its handler's return", rst)
        self.assertIn('``server_clock_step_back``', rst)

    def test_r3_the_suppression_reaches_not_after_ms_plus_the_skew(self):
        exact = AttemptAnswer(dict(COVERED[1], suppressed_until_ms=N + SKEW), 404, 'cancel')
        self.assertEqual('tombstone', self.verdict(exact).code)
        self.assertEqual('suppression_short', self.verdict(
            exact, attempts=dict(CAPS['level5'], clock_skew_allowance_ms=SKEW + 1)).code)
        # without a stated skew: against not_after_ms itself (the server judged
        # covers_admission with its own)
        self.assertEqual('tombstone', release_verdict(exact, SENT5).code)
        short = AttemptAnswer(dict(COVERED[1], suppressed_until_ms=N - 1), 404, 'cancel')
        self.assertEqual('suppression_short', release_verdict(short, SENT5).code)

    def test_r4_the_clock_against_running_and_stopping(self):
        limit = N + SKEW + BOUND_MS
        for state in ('running', 'stopping'):
            ans = AttemptAnswer({'attempt_id': A, 'server_instance': INST, 'state': state},
                                200, 'attempt')
            v = self.verdict(ans, now_ms=limit + 1, bound_ms=BOUND_MS)
            self.assertEqual((False, state), (v.release, v.code))
            v = self.verdict(ans, now_ms=limit + 1, bound_ms=BOUND_MS,
                             hold_clock_while_running=False)
            self.assertEqual((True, 'clock', ('uninterruptible_overrun',)),
                             (v.release, v.code, v.assumptions))
            # a 409 whose object says so is the same reply
            e = classify_409({'error': 'x', 'attempt': dict(_state(state))})
            self.assertEqual((False, state), (self.verdict(e, now_ms=limit + 1,
                                                           bound_ms=BOUND_MS).release, state))


class TestLightImport(unittest.TestCase):
    """LRG-G7: the client and the release rule load without the parser, the model or the
    operations; every public name of the package still resolves."""

    def run_py(self, code):
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
        out = subprocess.run([sys.executable, '-c', code], cwd=T.API, env=env,
                             capture_output=True, text=True, check=True)
        return json.loads(out.stdout)

    def test_the_client_is_light(self):
        heavy = ('parser', 'model', 'ops', 'derive', 'coords', 'store', 'export', 'budget',
                 'mcp_tools', 'frames', '_codec')
        for stmt in ('import metagraph.traverse.client',
                     'from metagraph.traverse import release_verdict, TraverseClient, '
                     'AttemptConflict, classify_409',
                     'import metagraph.traverse.attempts'):
            mods = self.run_py('import json, sys\n%s\nprint(json.dumps(sorted(m for m in '
                               'sys.modules if m.startswith("metagraph"))))' % stmt)
            self.assertFalse([m for m in mods if m.rsplit('.', 1)[-1] in heavy],
                             (stmt, mods))

    def test_every_public_name_resolves(self):
        got = self.run_py(
            'import json\nimport metagraph.traverse as MT\nfrom metagraph.traverse import *\n'
            'names = list(MT.__all__)\n'
            'print(json.dumps([names, [n for n in names if n not in globals()], '
            '[n for n in names if getattr(MT, n) is None], MT.ops.__name__, '
            'MT.client.TraverseClient is MT.TraverseClient, '
            'MT.client.release_verdict is MT.attempts.release_verdict]))')
        names, missing, nones, ops_name, same1, same2 = got
        self.assertEqual([], missing)
        self.assertEqual([], nones)
        self.assertEqual('metagraph.traverse.ops', ops_name)
        self.assertTrue(same1 and same2)
        self.assertIn('classify_409', names)


if __name__ == '__main__':
    unittest.main()
