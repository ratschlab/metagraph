"""A /search whose client is gone (lane T7-P2; SPEC-multigraph-envelope.md section 3): once the
client has closed, reset or half-closed its connection, or the server shuts down, the request
starts no further record, batch or pair (a record in progress finishes), writes nothing and logs
one line:

  [Server] Request <id>: <client gone|shutdown>, /search stopped: <r> of <n> records started,
  <p> of <P> pairs started; no response written

so that its handler, its threads and its budget accounts are free for the next request long
before the request would have ended. With the result cache on, a client gone stops nothing: the
request computes to its end and offers its answer to the cache, so that the client's retry is a
hit while the cache retains it (Fable review M1); an identical request that attaches to it while
it computes logs

  [Server] Request <id>: cache WAIT for an identical request still computing

before it waits (GPT review of train 7, P2), and a shutdown stops the producer, each identical
request waiting for it logging

  [Server] Request <id>: shutdown, /search stopped while waiting for request <producer>; no
  response written

The fixture is slow by construction: a telomere-like pair (50 records each holding (CCCTAA)n, so
that each repeat k-mer has about 330 coordinates in every record) searched with many records of
the repeat, query_coords. Each test measures the request without a disconnect first, closes at a
fraction of that time and states the numbers it saw (stderr).
"""

import json
import os
import random
import re
import shlex
import signal
import socket
import subprocess
import sys
import threading
import time
import unittest
from tempfile import TemporaryDirectory

import requests

from base import PROTEIN_MODE, TestingBase, METAGRAPH, graph_file_extension
from test_multigraph import Server

ROW_DIFF_BRWT_COORD = '.row_diff_brwt_coord.annodbg'
STOP = re.compile(r'Request (\d+): (client gone|shutdown), /search stopped: (\d+) of (\d+) '
                  r'records started, (\d+) of (\d+) pairs started; no response written\n')
STARTED = re.compile(r'\[Server\] /search request (\d+) from ')
# a request that waited for a stopped identical request (the result cache on a shutdown)
WAITED = re.compile(r'Request (\d+): shutdown, /search stopped while waiting for request (\d+); no '
                    r'response written\n')
# the client closes at this fraction of the undisturbed request's time, which must be at least
# MIN_FULL_S for the stop to be told from the request's end
CLOSE_AT = 0.15
MIN_FULL_S = 1.5


def note(*args):
    print('[test_search_stop]', *args, file=sys.stderr, flush=True)


@unittest.skipIf(PROTEIN_MODE, "the fixture is DNA")
class TestSearchStop(TestingBase):
    k = 31
    single_records = 100
    multi_records = 30
    names = [f't{i:02d}' for i in range(6)]
    started = 0
    # the single-graph request's answer without a disconnect and its wall time (-p 1
    # --threads-each 2), measured once
    reference = None

    @classmethod
    def setUpClass(cls):
        cls.tempdir = TemporaryDirectory()
        d = cls.tempdir.name
        rng = random.Random(23)
        flank = lambda: ''.join(rng.choice('AGT') for _ in range(40))
        cls.fasta = f'{d}/pt.fa'
        with open(cls.fasta, 'w') as f:
            for i in range(50):
                f.write(f'>t{i}\n{flank()}{"CCCTAA" * 333}{flank()}\n')
        cls.graph = f'{d}/pt' + graph_file_extension['succinct']
        anno = f'{d}/pt_anno'
        cls._build_graph(cls.fasta, f'{d}/pt', cls.k, 'succinct', 'basic')
        cls._annotate_graph(cls.fasta, cls.graph, anno, 'row_diff_brwt_coord',
                            anno_type='filename')
        res = subprocess.run(shlex.split(METAGRAPH) + [
            'annotate', '--anno-filename', '--index-header-coords', '-i', cls.graph, '-o', anno,
            cls.fasta], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        assert res.returncode == 0, res.stderr.decode()
        cls.anno = anno + ROW_DIFF_BRWT_COORD
        cls.csv = f'{d}/graphs.csv'
        with open(cls.csv, 'w') as f:
            for name in cls.names:
                f.write(f'{name},{cls.graph},{cls.anno}\n')
        cls.servers = []

    @classmethod
    def tearDownClass(cls):
        for server in cls.servers:
            server.stop()
        cls.tempdir.cleanup()

    def tearDown(self):
        for server in self.servers:
            server.stop()
        self.servers.clear()

    def server(self, extra=(), multi=False, p=1, e=2):
        TestSearchStop.started += 1
        log = os.path.join(self.tempdir.name, f'server_{TestSearchStop.started}.log')
        args = [self.csv] if multi else ['-i', self.graph, '-a', self.anno]
        # staging's --cache-size 0 unless the test sets the cache (a repeated request would
        # otherwise be a cache hit)
        if '--cache-size' not in extra:
            extra = ['--cache-size', '0'] + list(extra)
        server = Server(args + ['--threads-each', str(e)] + list(extra), log, threads=p)
        self.servers.append(server)
        return server

    # -------------------------------------------------------------------- requests

    def body(self, records, multi=False, name='q'):
        fasta = ''.join(f'>{name}{i}\n{"CCCTAA" * 50}\n' for i in range(records))
        body = {'FASTA': fasta, 'top_labels': 1000, 'discovery_fraction': 0.0,
                'query_coords': True}
        if multi:
            body['graphs'] = list(self.names)
        return body

    def timed(self, server, body):
        """the answer of |body| and its wall time"""
        t = time.monotonic()
        r = server.post('search', body)
        elapsed = time.monotonic() - t
        self.assertEqual(200, r.status_code, r.text[:2000])
        return r, elapsed

    def send(self, server, body):
        """|body| POSTed on a raw socket, which the caller closes (the client leaving)"""
        data = json.dumps(body).encode()
        sock = socket.create_connection(('127.0.0.1', server.port))
        sock.sendall((f'POST /search HTTP/1.1\r\nHost: 127.0.0.1:{server.port}\r\n'
                      f'Content-Type: application/json\r\nContent-Length: {len(data)}\r\n\r\n')
                     .encode() + data)
        return sock

    def request_ids(self, server, before):
        return [int(m.group(1)) for m in STARTED.finditer(server.text()) if m.start() >= before]

    def wait_for(self, server, before, pattern, timeout):
        """the first match of |pattern| after |before| bytes of the log, and when it was seen"""
        deadline = time.monotonic() + timeout
        while True:
            m = next((m for m in pattern.finditer(server.text()) if m.start() >= before), None)
            if m:
                return m, time.monotonic()
            if time.monotonic() > deadline:
                self.fail(f'no {pattern.pattern!r} within {timeout} s: {server.text()[before:][-3000:]}')
            time.sleep(0.02)

    def stopped_quietly(self, server, before, request_id):
        """no line of |request_id| after |before| but its start and its stop: no answer, no budget
        line"""
        text = server.text()[before:]
        for line in text.splitlines():
            if re.search(rf'Request {request_id}\b', line) or re.search(rf'request {request_id} ', line):
                self.assertTrue(STOP.search(line + '\n') or STARTED.search(line), line)

    def close_after(self, full_s):
        """when the client closes, after its request started: CLOSE_AT of the undisturbed
        request's |full_s|, which must be long enough to observe a stop"""
        self.assertGreater(full_s, MIN_FULL_S, 'the fixture runs too fast here to observe a stop')
        return CLOSE_AT * full_s

    def disconnect(self, server, body, after_s, then=None):
        """|body| sent and its client gone after |after_s| (|then| called right after the
        close): (request id, the stop line's match, the seconds from the close to the line, the
        log's length before)"""
        before = len(server.text())
        sock = self.send(server, body)
        m, _ = self.wait_for(server, before, STARTED, 30)
        request_id = int(m.group(1))
        time.sleep(after_s)
        sock.close()
        closed = time.monotonic()
        follower = threading.Thread(target=then, args=(closed,)) if then else None
        if follower:
            follower.start()
        m, seen = self.wait_for(server, before, STOP, 120)
        if follower:
            follower.join(timeout=300)
        self.assertEqual(request_id, int(m.group(1)))
        return request_id, m, seen - closed, before

    def single_reference(self):
        """the single-graph request's answer text and wall time without a disconnect"""
        if TestSearchStop.reference is None:
            server = self.server(p=1, e=2)
            r, full_s = self.timed(server, self.body(self.single_records))
            TestSearchStop.reference = (r.text, full_s)
            self.stop_server(server)
        return TestSearchStop.reference

    # -------------------------------------------------------------------- tests

    def test_single_graph_stop_frees_the_handler(self):
        """-p 1: one handler. The stopped request leaves it free: a second /search is answered
        long before the first would have ended"""
        reference, full_s = self.single_reference()
        server = self.server(p=1, e=2)
        body = self.body(self.single_records)
        small = self.body(1, name='s')
        small_ref, _ = self.timed(server, small)
        # a second request sent right after the close: the one handler answers it once the first
        # stopped
        second = {}

        def next_request(closed):
            second['r'] = server.post('search', small)
            second['s'] = time.monotonic() - closed

        close_s = self.close_after(full_s)
        request_id, m, stop_s, before = self.disconnect(server, body, close_s,
                                                        then=next_request)
        self.assertEqual('client gone', m.group(2))
        started, records, pairs_started, pairs = (int(m.group(i)) for i in (3, 4, 5, 6))
        self.assertEqual(self.single_records, records)
        self.assertLess(started, records)
        self.assertEqual((1, 1), (pairs_started, pairs))
        self.assertEqual(200, second['r'].status_code)
        self.assertEqual(small_ref.text, second['r'].text)
        second_s = second['s']
        note(f'single graph -p 1 --threads-each 2: full {full_s:.2f} s; closed at {close_s:.2f} '
             f's, stopped {stop_s:.2f} s after the close ({started} of {records} records '
             f'started); the next request answered {second_s:.2f} s after the close')
        # the stop comes within a record or two of the close, not at the request's end, and the
        # handler is free then: the next request is answered while the first would still run
        self.assertLess(close_s + stop_s, 0.6 * full_s)
        self.assertLess(close_s + second_s, 0.6 * full_s)
        self.stopped_quietly(server, before, request_id)
        # and the undisturbed request is answered as before, byte for byte
        r, _ = self.timed(server, body)
        self.assertEqual(reference, r.text)

    def test_multi_graph_stop_states_the_pairs(self):
        """staging's flags (-p 2 --threads-each 2): the pairs not started stay so"""
        server = self.server(multi=True, p=2, e=2)
        body = self.body(self.multi_records, multi=True)
        reference, full_s = self.timed(server, body)
        close_s = self.close_after(full_s)
        request_id, m, stop_s, before = self.disconnect(server, body, close_s)
        self.assertEqual('client gone', m.group(2))
        started, records, pairs_started, pairs = (int(m.group(i)) for i in (3, 4, 5, 6))
        self.assertEqual(self.multi_records, records)
        self.assertEqual(len(self.names), pairs)
        self.assertLess(pairs_started, pairs)
        self.assertGreaterEqual(pairs_started, 1)
        self.assertLessEqual(started, records)
        note(f'multi graph -p 2 --threads-each 2, {pairs} pairs: full {full_s:.2f} s; closed at '
             f'{close_s:.2f} s, stopped {stop_s:.2f} s after the close ({started} of {records} '
             f'records, {pairs_started} of {pairs} pairs started)')
        self.assertLess(close_s + stop_s, 0.6 * full_s)
        self.stopped_quietly(server, before, request_id)
        r, _ = self.timed(server, body)
        self.assertEqual(reference.text, r.text)

    def test_shutdown_stops_a_search(self):
        """SIGTERM during a /search: it stops at its next record, nothing is written, and the
        handler returns, so the process exits without the 'did not return' warning"""
        _, full_s = self.single_reference()
        server = self.server(p=1, e=2)
        body = self.body(self.single_records)
        before = len(server.text())
        sock = self.send(server, body)
        self.wait_for(server, before, STARTED, 30)
        time.sleep(self.close_after(full_s))
        signalled = time.monotonic()
        server.process.send_signal(signal.SIGTERM)
        self.assertEqual(0, server.process.wait(timeout=60))
        exit_s = time.monotonic() - signalled
        sock.settimeout(10)
        received = b''
        try:
            while True:
                chunk = sock.recv(65536)
                if not chunk:
                    break
                received += chunk
        except OSError:
            pass
        sock.close()
        text = server.text()[before:]
        m = STOP.search(text)
        self.assertTrue(m, text[-3000:])
        self.assertEqual('shutdown', m.group(2))
        self.assertLess(int(m.group(3)), int(m.group(4)))
        self.assertNotIn('did not return within', text)
        self.assertIn('[Server] Stopped', text)
        self.assertEqual(b'', received)
        note(f'shutdown: full {full_s:.2f} s; exited {exit_s:.2f} s after SIGTERM '
             f'({m.group(3)} of {m.group(4)} records started)')

    def attached(self, server, request_id, producer, timeout=30):
        """waits until |request_id| is attached to |producer|'s pending entry, which its handler
        logs before it waits ("cache WAIT"; GPT review of train 7, P2: inferred before from its
        parsed body, after which a descheduled handler could attach once the answer was stored
        and pass as a plain hit), and returns where; fails if |producer| stores first"""
        deadline = time.monotonic() + timeout
        while True:
            text = server.text()
            m = re.search(rf'\[Server\] Request {request_id}: cache WAIT for an identical '
                          r'request still computing', text)
            if m:
                return m.start()
            self.assertIsNone(re.search(rf'Request {producer}: cache STORE', text),
                              f'request {producer} stored before {request_id} attached')
            self.assertLess(time.monotonic(), deadline, text[-3000:])
            time.sleep(0.02)

    def test_cache_waiter_gets_the_whole_answer(self):
        """--cache-size above 0: the producer computes on when its client leaves, and a duplicate
        that waits for its entry gets the answer byte for byte"""
        reference, full_s = self.single_reference()
        # (the fixture runs long enough to observe a stop)
        self.close_after(full_s)
        body = self.body(self.single_records)
        server = self.server(['--cache-size', '1', '-v'], p=2, e=2)
        before = len(server.text())
        sock = self.send(server, body)
        self.wait_for(server, before, STARTED, 30)
        producer = self.request_ids(server, before)[0]
        answer = {}

        def duplicate():
            answer['r'] = server.post('search', body)

        waiter = threading.Thread(target=duplicate)
        waiter.start()
        deadline = time.monotonic() + 30
        while len(self.request_ids(server, before)) < 2:
            self.assertLess(time.monotonic(), deadline)
            time.sleep(0.02)
        dup = self.request_ids(server, before)[1]
        # the duplicate attached to the pending entry, then the producer's client leaves
        attached_at = self.attached(server, dup, producer)
        sock.close()
        at_close = len(server.text())
        waiter.join(timeout=300)
        r = answer['r']
        self.assertEqual(200, r.status_code, r.text[:2000])
        self.assertEqual(reference, r.text)
        text = server.text()[before:]
        self.assertIsNone(STOP.search(text), text[-3000:])
        # stored after the close: the producer ran on past its client's leaving
        store = re.search(rf'Request {producer}: cache STORE', server.text())
        self.assertTrue(store, text[-3000:])
        self.assertGreaterEqual(store.start(), at_close)
        self.assertRegex(text, rf'Request {dup}: cache HIT')
        # and the duplicate had attached before the close, so before STORE: a waiter, not a hit
        # on the stored answer
        self.assertLess(attached_at, at_close)

    def test_shutdown_with_a_waiter(self):
        """--cache-size above 0, SIGTERM while a duplicate waits for the producer: both are
        written nothing; the producer's line has the counts, the waiter's names the producer"""
        _, full_s = self.single_reference()
        # (the fixture runs long enough to observe a stop)
        self.close_after(full_s)
        body = self.body(self.single_records)
        server = self.server(['--cache-size', '1', '-v'], p=2, e=2)
        before = len(server.text())
        socks = [self.send(server, body)]
        self.wait_for(server, before, STARTED, 30)
        socks.append(self.send(server, body))
        deadline = time.monotonic() + 30
        while len(self.request_ids(server, before)) < 2:
            self.assertLess(time.monotonic(), deadline)
            time.sleep(0.02)
        producer, dup = self.request_ids(server, before)[:2]
        attached_at = self.attached(server, dup, producer)
        at_signal = len(server.text())
        self.assertLess(attached_at, at_signal)
        server.process.send_signal(signal.SIGTERM)
        self.assertEqual(0, server.process.wait(timeout=60))
        for sock in socks:
            sock.settimeout(10)
            received = b''
            try:
                while True:
                    chunk = sock.recv(65536)
                    if not chunk:
                        break
                    received += chunk
            except OSError:
                pass
            sock.close()
            self.assertEqual(b'', received)
        text = server.text()[before:]
        stops = list(STOP.finditer(text))
        self.assertEqual([(producer, 'shutdown')], [(int(m.group(1)), m.group(2)) for m in stops],
                         text[-3000:])
        self.assertLess(int(stops[0].group(3)), int(stops[0].group(4)))
        waited = list(WAITED.finditer(text))
        self.assertEqual([(dup, producer)], [(int(m.group(1)), int(m.group(2))) for m in waited],
                         text[-3000:])
        self.assertNotIn('cache STORE', text)
        self.assertNotIn('cache HIT', text)
        note(f'shutdown with a waiter: producer {producer} stopped ({stops[0].group(3)} of '
             f'{stops[0].group(4)} records), waiter {dup} logged its own line')

    def test_cache_retry_gets_the_stored_answer(self):
        """--cache-size above 0 (Fable review M1): a client that closes and sends its request
        again gets the stored answer. The first request computes to its end although its client
        left and stores it; the retry is a hit (as a waiter, or after the store), the answer byte
        for byte the undisturbed one's"""
        reference, full_s = self.single_reference()
        body = self.body(self.single_records)
        server = self.server(['--cache-size', '1'], p=2, e=2)
        close_s = self.close_after(full_s)
        before = len(server.text())
        sock = self.send(server, body)
        self.wait_for(server, before, STARTED, 30)
        first = self.request_ids(server, before)[0]
        time.sleep(close_s)
        sock.close()
        at_close = len(server.text())
        closed = time.monotonic()
        time.sleep(0.05)
        r = server.post('search', body)
        retry_s = time.monotonic() - closed
        self.assertEqual(200, r.status_code, r.text[:2000])
        self.assertEqual(reference, r.text)
        text = server.text()[before:]
        self.assertIsNone(STOP.search(text), text[-3000:])
        store = re.search(rf'Request {first}: cache STORE', server.text())
        self.assertTrue(store, text[-3000:])
        self.assertGreaterEqual(store.start(), at_close)
        retry = self.request_ids(server, before)[1]
        self.assertRegex(text, rf'Request {retry}: cache HIT')
        self.assertNotRegex(text, rf'Request {retry}: cache STORE')
        note(f'cache retry: full {full_s:.2f} s; closed at {close_s:.2f} s, the retry (sent 0.05 '
             f's after the close) answered {retry_s:.2f} s after the close')
        # the retry waited for the rest of the first request, not for a computation of its own
        self.assertLess(close_s + retry_s, 1.5 * full_s)

    def test_half_close_is_a_client_gone(self):
        """A client that half-closes (shutdown SHUT_WR) after its request reads as gone, as on
        /traverse (Fable review M2, accepted): nothing is written to it and the connection is
        closed; the stop line says client gone"""
        _, full_s = self.single_reference()
        self.close_after(full_s)
        server = self.server(p=1, e=2)
        body = self.body(self.single_records)
        before = len(server.text())
        sock = self.send(server, body)
        sock.shutdown(socket.SHUT_WR)
        sent = time.monotonic()
        sock.settimeout(60)
        received = b''
        try:
            while True:
                chunk = sock.recv(65536)
                if not chunk:
                    break
                received += chunk
        except ConnectionResetError:
            pass  # closed (a reset): nothing more to read either
        closed_s = time.monotonic() - sent
        sock.close()
        self.assertEqual(b'', received)
        m, _ = self.wait_for(server, before, STOP, 30)
        self.assertEqual('client gone', m.group(2))
        self.assertLess(int(m.group(3)), int(m.group(4)))
        self.stopped_quietly(server, before, int(m.group(1)))
        note(f'half-close: full {full_s:.2f} s; the connection closed {closed_s:.2f} s after the '
             f'half-close, nothing read ({m.group(3)} of {m.group(4)} records started)')
        self.assertLess(closed_s, 0.6 * full_s)

    def test_in_ram_wait_is_abandoned(self):
        """`in_ram` on a server on mmap whose --mem-cap-gb holds one pair: a request waiting for
        the memory another request's load holds abandons the wait once its client is gone (it
        loads nothing), while the other request runs on to its answer"""
        size = os.path.getsize(self.graph) + os.path.getsize(self.anno)
        server = self.server(['--mmap', '--mem-cap-gb', f'{1.5 * size / 1e9:.9f}', '-v'],
                             multi=True, p=2, e=2)
        holder = dict(self.body(self.single_records), graphs=['t00'], in_ram=True)
        waiting = dict(self.body(self.multi_records), graphs=['t01'], in_ram=True)
        before = len(server.text())
        answer = {}

        def hold():
            answer['r'] = server.post('search', holder)
            answer['at'] = time.monotonic()

        holding = threading.Thread(target=hold)
        holding.start()
        loading = re.compile(r'Request (\d+): Loading graph ')
        m, _ = self.wait_for(server, before, loading, 60)
        holder_id = int(m.group(1))
        after_load = len(server.text())
        sock = self.send(server, waiting)
        m, _ = self.wait_for(server, after_load, STARTED, 30)
        waiting_id = int(m.group(1))
        # its pair's task began (and waits for the memory)
        self.wait_for(server, after_load,
                      re.compile(rf'Request {waiting_id}: Started querying graph '), 30)
        time.sleep(0.3)
        too_fast = ('the holder answered before {}: the fixture runs too fast here to tell an '
                    'abandoned wait from one that ended with the holder\'s release')
        self.assertTrue(holding.is_alive(), too_fast.format('the close'))
        sock.close()
        closed = time.monotonic()
        m, seen = self.wait_for(server, after_load, STOP, 60)
        stop_s = seen - closed
        holding.join(timeout=300)
        # the holder outlasted the stop line, so the memory was never freed for the waiting one
        self.assertLess(seen, answer['at'], too_fast.format('the stop line'))
        self.assertEqual(waiting_id, int(m.group(1)))
        self.assertEqual(('client gone', 0, self.multi_records, 0, 1),
                         (m.group(2), *(int(m.group(i)) for i in (3, 4, 5, 6))))
        r = answer['r']
        self.assertEqual(200, r.status_code, r.text[:2000])
        # the wait was abandoned: no load of the waiting request's
        text = server.text()[before:]
        self.assertNotRegex(text, rf'Request {waiting_id}: Loading graph ')
        self.assertRegex(text, rf'Request {holder_id}: Loaded and searched graph ')
        note(f'in_ram wait: stopped {stop_s:.2f} s after the close, the holder answered '
             f'{answer["at"] - seen:.2f} s later')
        self.assertLess(stop_s, 1.0)

    def test_stop_wins_over_a_refusal(self):
        """A request refused on the answer side (no partial_refusal: once its answer account
        passed, its units' searches run on, building nothing, and it is refused 400
        search_answer_too_large at its end) whose client leaves after the pass and before the end
        is stopped silently: the stop line, no refusal or budget line"""
        # an unbounded run's budget line: the most one unit's values held, and the answer's
        # bytes a record
        server = self.server(['--search-max-answer-memory-mb', '100000'], p=1, e=2)
        before = len(server.text())
        self.timed(server, self.body(self.single_records))
        m = re.compile(r'answer budget: values at most (\d+) B a unit, answer (\d+) B').search(
                server.text(), before)
        self.assertTrue(m, server.text()[before:][-3000:])
        values, answer = int(m.group(1)), int(m.group(2))
        self.stop_server(server)
        # enough records for an answer of 2.5 times a unit's values, and a bound of 1.2 times
        # them: every unit's values fit, the answer passes at about half the records
        records = -(-int(2.5 * values) * self.single_records // answer)
        bound_mb = -(-int(1.2 * values) >> 20)
        # the records answered when the answer account passes the bound (its bytes a record as
        # measured), and the most that can be in progress beside them (--threads-each)
        passed_at = -(-(bound_mb << 20) * self.single_records // answer)
        in_progress = 2
        body = self.body(records)
        server = self.server(['--search-max-answer-memory-mb', str(bound_mb)], p=1, e=2)
        t = time.monotonic()
        r = server.post('search', body)
        refused_s = time.monotonic() - t
        self.assertEqual(400, r.status_code, r.text[:2000])
        self.assertIn('search_answer_too_large', r.text)
        # closed at three quarters of the refused request's time: after the pass, before the end
        close_s = 0.75 * refused_s
        self.assertGreater(refused_s, 4 * MIN_FULL_S, 'too fast here to close between the two')
        request_id, m, stop_s, before = self.disconnect(server, body, close_s)
        self.assertEqual('client gone', m.group(2))
        self.assertEqual(records, int(m.group(4)))
        self.assertLess(int(m.group(3)), records)
        # the bound had passed at the stop: more records had started than the pass needs
        # answered and the threads can hold in progress
        self.assertGreater(int(m.group(3)), passed_at + in_progress,
                           f'the close came before the answer passed the bound (at about '
                           f'{passed_at} of {records} records)')
        self.stopped_quietly(server, before, request_id)
        text = server.text()[before:]
        self.assertNotIn('search_answer_too_large', text)
        self.assertNotIn('answer budget', text)
        self.assertNotIn('unit refused', text)
        note(f'stop over a refusal: values at most {values} B a unit, {records} records, bound '
             f'{bound_mb} MiB (passed at about {passed_at} records); refused after '
             f'{refused_s:.2f} s; closed at {close_s:.2f} s, stopped {stop_s:.2f} s after the '
             f'close ({m.group(3)} of {m.group(4)} records)')
        self.assertLess(stop_s, 0.15 * refused_s)

    def stop_server(self, server):
        server.stop()
        self.servers.remove(server)


if __name__ == '__main__':
    unittest.main()
