"""Feature level 5 against a live server_query on build/mini_refseq: the client's level-5
fields, typed answers and release_verdict() on what the server really answers (the unit
tests, test_traverse_level5_client.py, replay recorded answers).

The server binary: $METAGRAPH_LEVEL5_SERVER, else the benchmark cache's level-5 build
(~/.cache/metagraph-traverse-bench/bin_dcc0cebd/metagraph); the index:
$METAGRAPH_MINI_REFSEQ, else build/mini_refseq (scripts/traversal/build_mini_refseq.sh). Both
absent, or a server below feature level 5: skipped. A level-3 server
($METAGRAPH_LEVEL3_SERVER, else the cache's bin_ce949da5) checks the gate when present.
Three servers are started (the default retention, --traverse-attempt-retention 1 for the 429
tombstones_full and --traverse-attempt-retention-s 0 for no_suppression), each on a free
port, and stopped at the end. With both binaries, one port is then served in turn by the
level-5 server, the level-3 one (a rollback) and the level-5 one again (a restart): the
client notices each replacement from the answers.
"""

import json
import os
import socket
import subprocess
import sys
import time
import unittest
import urllib.request
import uuid

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402

from metagraph.traverse import (  # noqa: E402
    AttemptConflict, AttemptExpired, AttemptSent, InstanceMismatch, TraverseClient,
    UnsupportedFeature, release_verdict,
)

_CACHE = os.path.expanduser('~/.cache/metagraph-traverse-bench')
SERVER = os.environ.get('METAGRAPH_LEVEL5_SERVER',
                        os.path.join(_CACHE, 'bin_dcc0cebd', 'metagraph'))
OLD_SERVER = os.environ.get('METAGRAPH_LEVEL3_SERVER',
                            os.path.join(_CACHE, 'bin_ce949da5', 'metagraph'))
INDEX = os.environ.get('METAGRAPH_MINI_REFSEQ', os.path.join(T.REPO, 'build', 'mini_refseq'))
GRAPH = os.path.join(INDEX, 'graph_k31.dbg')
ANNO = os.path.join(INDEX, 'annotation.relaxed.relabeled.row_diff_brwt_coord.annodbg')
MANIFEST = os.path.join(INDEX, 'annotation.relaxed.relabeled.manifest.json')
# the start of NZ_LPPQ01000025.1 (mini_refseq's 1296536.fa)
SEED = 'ACGACGAAAGCGATGGTATCCCTGCCCGACATATGCAAGAGGCTCATCGGCTAGCGTGTTGCTGAGTTGA'
STRATEGY = {'direction': 'right', 'bounds': {'max_extension_bp': 40}}


def _have():
    return all(os.path.exists(p) for p in (SERVER, GRAPH, ANNO, MANIFEST))


def _free_port():
    with socket.socket() as s:
        s.bind(('127.0.0.1', 0))
        return s.getsockname()[1]


class _Server:
    def __init__(self, binary, *options, port=None):
        self.port = port or _free_port()
        self.log = open(os.devnull, 'w')
        for attempt in range(20):
            self.proc = subprocess.Popen(
                [binary, 'server_query', '-i', GRAPH, '-a', ANNO, '--index-manifest',
                 MANIFEST, '--address', '127.0.0.1', '--port', str(self.port), '-p', '3']
                + list(options), stdout=self.log, stderr=subprocess.STDOUT)
            try:
                self.doc = self._wait()
                return
            except RuntimeError:
                # a port just given up by the previous server may not be free yet
                if port is None or self.proc.poll() is None or attempt == 19:
                    self.stop()
                    raise
                time.sleep(0.5)

    def _wait(self, timeout=120):
        t0 = time.time()
        while time.time() - t0 < timeout:
            if self.proc.poll() is not None:
                raise RuntimeError('server_query exited (%s)' % self.proc.returncode)
            try:
                with urllib.request.urlopen('http://127.0.0.1:%d/capabilities' % self.port,
                                            timeout=2) as r:
                    doc = json.loads(r.read())
                if doc.get('ready', True):
                    return doc
            except OSError:
                pass
            time.sleep(0.2)
        raise RuntimeError('server_query did not answer in %d s' % timeout)

    def client(self, **kw):
        return TraverseClient('127.0.0.1', self.port, timeout=60, **kw)

    def stop(self):
        self.proc.terminate()
        try:
            self.proc.wait(10)
        except subprocess.TimeoutExpired:
            self.proc.kill()
            self.proc.wait()
        if not self.log.closed:
            self.log.close()


def _id(tag):
    return 'live-%s-%s' % (tag, uuid.uuid4().hex[:12])


def _now_ms():
    return int(time.time() * 1000)


@unittest.skipUnless(_have(), 'needs a level-5 server_query and build/mini_refseq')
class TestLevel5Live(unittest.TestCase):
    servers = []

    @classmethod
    def setUpClass(cls):
        cls.main = _Server(SERVER)
        cls.servers = [cls.main]
        if int(cls.main.doc.get('feature_level') or 0) < 5:
            cls.tearDownClass()
            raise unittest.SkipTest('the server states feature_level %s'
                                    % cls.main.doc.get('feature_level'))
        cls.full = _Server(SERVER, '--traverse-attempt-retention', '1')
        cls.servers.append(cls.full)
        cls.nosupp = _Server(SERVER, '--traverse-attempt-retention-s', '0')
        cls.servers.append(cls.nosupp)
        cls.c = cls.main.client()
        cls.skew = cls.main.doc['attempts']['clock_skew_allowance_ms']

    @classmethod
    def tearDownClass(cls):
        for s in cls.servers:
            s.stop()
        cls.servers = []

    def test_the_client_reads_the_level_and_the_instance(self):
        self.assertGreaterEqual(self.c.server_feature_level(), 5)
        inst = self.c.server_instance()
        self.assertRegex(inst, '^[0-9a-f]{16}$')
        self.assertIn('expect_server_instance', self.main.doc['attempts']['fields'])
        self.assertIn('not_after_ms', self.main.doc['attempts']['cancel_fields'])
        self.assertIn('release_rule', self.main.doc['attempts'])

    def test_a_covering_tombstone_releases_early_and_refuses_the_copy(self):
        aid = _id('tomb')
        naf = _now_ms() + 60000
        sent = AttemptSent(aid, naf, self.c.server_instance())
        got = self.c.cancel(aid, not_after_ms=naf)
        self.assertEqual((404, True, 'unknown'), (got.status, got.tombstone, got.state))
        sup = got.suppression
        self.assertEqual((naf, True), (sup.not_after_ms, sup.covers_admission))
        self.assertGreaterEqual(sup.suppressed_until_ms, naf + self.skew)
        v = release_verdict(got, sent)
        self.assertEqual((True, True, 'tombstone'), (v.release, v.early, v.code), v.why)
        other = AttemptSent(aid, naf + 1, sent.expect_server_instance)
        self.assertEqual('other_not_after_ms', release_verdict(got, other).code)
        self.assertEqual('sent_without_expect_server_instance',
                         release_verdict(got, AttemptSent(aid, naf)).code)
        # GET states the same suppression and releases the same
        state = self.c.attempt(aid)
        self.assertEqual((404, True), (state.status, state.suppression.covers_admission))
        self.assertTrue(release_verdict(state, sent))
        # the delayed copy is refused, its 409 carrying the tombstone judged against its own
        # not_after_ms; the 409 itself releases nothing
        with self.assertRaises(AttemptConflict) as cm:
            self.c.traverse([SEED], STRATEGY, attempt_id=aid, not_after_ms=naf,
                            expect_server_instance='auto')
        e = cm.exception
        self.assertTrue(e.tombstoned)
        self.assertEqual((naf, True),
                         (e.suppression.not_after_ms, e.suppression.covers_admission))
        self.assertEqual(sent, e.sent)
        self.assertEqual('not_a_release_answer', release_verdict(e, e.sent).code)

    def test_tombstones_that_do_not_cover(self):
        aid = _id('plain')
        got = self.c.cancel(aid)
        self.assertEqual((False, 'no_not_after_ms'),
                         (got.suppression.covers_admission,
                          got.suppression.covers_admission_reason))
        v = release_verdict(got, AttemptSent(aid, _now_ms() + 60000, self.c.server_instance()))
        self.assertEqual('covers_admission_false', v.code)
        aid = _id('far')
        far = _now_ms() + self.main.doc['attempts']['tombstone_max_s'] * 1000 + 10 ** 6
        got = self.c.cancel(aid, not_after_ms=far)
        self.assertEqual('beyond_tombstone_max', got.suppression.covers_admission_reason)
        self.assertFalse(release_verdict(got, AttemptSent(aid, far, self.c.server_instance())))

    def test_a_completed_attempt_releases_on_its_response_and_state(self):
        aid = _id('done')
        naf = _now_ms() + 60000
        resp = self.c.traverse([SEED], STRATEGY, attempt_id=aid, not_after_ms=naf,
                               expect_server_instance='auto')
        self.assertEqual(AttemptSent(aid, naf, self.c.server_instance()), resp.sent)
        self.assertEqual(('finished', True), (release_verdict(resp, resp.sent).code,
                                              release_verdict(resp, resp.sent).release))
        self.assertTrue(release_verdict(self.c.attempt(aid), resp.sent))
        # the server holds the id against a replay: the release assumes nothing
        self.assertEqual((), release_verdict(resp, resp.sent,
                                             attempts=self.main.doc['attempts']).assumptions)
        got = self.c.cancel(aid, not_after_ms=naf)
        self.assertEqual((404, 'finished'), (got.status, got.state))
        self.assertTrue(release_verdict(got, resp.sent))
        # a replay of the finished request is refused while the id is held
        with self.assertRaises(AttemptConflict) as cm:
            self.c.traverse([SEED], STRATEGY, attempt_id=aid, not_after_ms=naf,
                            expect_server_instance='auto')
        self.assertEqual(('finished', False), (cm.exception.state, cm.exception.tombstoned))

    def test_instance_mismatch_and_expired(self):
        with self.assertRaises(InstanceMismatch) as cm:
            self.c.traverse([SEED], STRATEGY, attempt_id=_id('other'),
                            expect_server_instance='0123456789abcdef')
        self.assertEqual(('0123456789abcdef', self.c.server_instance()),
                         (cm.exception.expect_server_instance, cm.exception.server_instance))
        with self.assertRaises(AttemptExpired) as cm:
            self.c.traverse([SEED], STRATEGY, attempt_id=_id('late'), not_after_ms=1)
        self.assertEqual((1, self.c.server_instance()),
                         (cm.exception.not_after_ms, cm.exception.server_instance))
        self.assertGreater(cm.exception.server_time_ms, 1)

    def test_the_429s(self):
        naf = _now_ms() + 60000
        c = self.full.client()
        first = c.cancel(_id('full'), not_after_ms=naf)
        aid = _id('full')
        got = c.cancel(aid, not_after_ms=naf)
        if first.status == 404:
            self.assertTrue(first.tombstone)
        self.assertEqual((429, False, 'tombstones_full', True),
                         (got.status, got.tombstone, got.reason, got.retryable))
        self.assertEqual('not_tombstoned', release_verdict(
            got, AttemptSent(aid, naf, c.server_instance())).code)
        c = self.nosupp.client()
        aid = _id('none')
        got = c.cancel(aid, not_after_ms=naf)
        self.assertEqual((429, 'no_suppression', False),
                         (got.status, got.reason, got.retryable))
        self.assertFalse(release_verdict(got, AttemptSent(aid, naf, c.server_instance())))

    def test_a_finished_release_on_a_server_that_keeps_nothing_is_assumed(self):
        c = self.nosupp.client()
        aid = _id('keepnone')
        naf = _now_ms() + 60000
        resp = c.traverse([SEED], STRATEGY, attempt_id=aid, not_after_ms=naf,
                          expect_server_instance='auto')
        v = release_verdict(resp, resp.sent, attempts=self.nosupp.doc['attempts'])
        self.assertEqual((True, 'finished', ('no_retention',)),
                         (v.release, v.code, v.assumptions), v.why)
        # what it assumes: the same request arriving again runs again
        again = c.traverse([SEED], STRATEGY, attempt_id=aid, not_after_ms=naf,
                           expect_server_instance=resp.sent.expect_server_instance)
        self.assertEqual(aid, again.usage['attempt_id'])


@unittest.skipUnless(_have() and os.path.exists(OLD_SERVER),
                     'needs a level-3 server_query and build/mini_refseq')
class TestBelowLevel5Live(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.old = _Server(OLD_SERVER)
        if int(cls.old.doc.get('feature_level') or 0) >= 5:
            cls.old.stop()
            raise unittest.SkipTest('not a server below feature level 5')

    @classmethod
    def tearDownClass(cls):
        cls.old.stop()

    def test_the_level5_fields_are_not_sent(self):
        c = self.old.client()
        aid = _id('old')
        got = c.cancel(aid, not_after_ms=_now_ms() + 60000)
        # the cancel went out without not_after_ms (the server would refuse it: 400)
        self.assertEqual((404, True, None), (got.status, got.tombstone, got.suppression))
        self.assertNotIn('not_after_ms', got.sent)
        self.assertEqual('no_suppression_stated', release_verdict(
            got, AttemptSent(aid, 1, got.server_instance)).code)
        aid = _id('old')
        resp = c.traverse([SEED], STRATEGY, attempt_id=aid, expect_server_instance='auto')
        self.assertEqual(AttemptSent(aid), resp.sent)
        self.assertIsNotNone(resp.graphlets[0])
        with self.assertRaises(UnsupportedFeature):
            c.traverse([SEED], STRATEGY, attempt_id=_id('old'),
                       expect_server_instance=got.server_instance)


@unittest.skipUnless(_have() and os.path.exists(OLD_SERVER),
                     'needs a level-5 and a level-3 server_query and build/mini_refseq')
class TestReplacedServerLive(unittest.TestCase):
    """One address, three processes: the client keeps what GET /capabilities stated only
    for the process that stated it."""

    def test_a_rollback_and_a_restart(self):
        first = _Server(SERVER)
        port = first.port
        try:
            if int(first.doc.get('feature_level') or 0) < 5:
                self.skipTest('the server states feature_level %s'
                              % first.doc.get('feature_level'))
            c = first.client()
            self.assertGreaterEqual(c.server_feature_level(), 5)
            inst = c.server_instance()
        finally:
            first.stop()
        old = _Server(OLD_SERVER, port=port)
        try:
            if int(old.doc.get('feature_level') or 0) >= 5:
                self.skipTest('not a server below feature level 5')
            # the rollback: the level kept says 5, the process refuses not_after_ms; the
            # cancel goes out again without it and is carried out (the id tombstoned)
            aid = _id('rolled')
            naf = _now_ms() + 60000
            got = c.cancel(aid, not_after_ms=naf)
            self.assertEqual((404, True), (got.status, got.tombstone))
            self.assertNotIn('not_after_ms', got.sent)
            self.assertEqual(int(old.doc['feature_level']), c.server_feature_level())
            # the request it cancelled is refused, not run
            with self.assertRaises(AttemptConflict):
                c.traverse([SEED], STRATEGY, attempt_id=aid, not_after_ms=naf)
            # and 'auto' is left out for this process
            aid = _id('rolled')
            resp = c.traverse([SEED], STRATEGY, attempt_id=aid, not_after_ms=naf,
                              expect_server_instance='auto')
            self.assertEqual(AttemptSent(aid, naf), resp.sent)
        finally:
            old.stop()
        again = _Server(SERVER, port=port)
        try:
            # the restart at the same level: an answer names the new instance, and the
            # next 'auto' request names it (no instance_mismatch round trip)
            c.server_instance()
            self.assertEqual(int(old.doc['feature_level']), c.server_feature_level())
            c.attempt(_id('probe'))
            aid = _id('restarted')
            resp = c.traverse([SEED], STRATEGY, attempt_id=aid, not_after_ms=_now_ms() + 60000,
                              expect_server_instance='auto')
            new_inst = again.doc['attempts']['server_instance']
            self.assertNotEqual(inst, new_inst)
            self.assertEqual(new_inst, resp.sent.expect_server_instance)
            self.assertEqual(aid, resp.usage['attempt_id'])
        finally:
            again.stop()


if __name__ == '__main__':
    unittest.main()
