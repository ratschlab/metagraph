"""A /search whose client is gone (lane T7-P2; SPEC-multigraph-envelope.md section 3): once the
client has closed its connection, or the server shuts down, the request starts no further record,
batch or pair (a record in progress finishes), writes nothing and logs one line:

  [Server] Request <id>: <client gone|shutdown>, /search stopped: <r> of <n> records started,
  <p> of <P> pairs started; no response written

so that its handler, its threads and its budget accounts are free for the next request long
before the request would have ended. With the result cache, a producer another identical request
waits for computes on for it, and a stopped one leaves the cache: an identical request after it
computes the answer itself.

The fixture is slow by construction: a telomere-like pair (50 records each holding (CCCTAA)n, so
that each repeat k-mer has about 330 coordinates in every record) searched with many records of
the repeat, query_coords. Each test measures the request without a disconnect first and states
the numbers it saw (stderr).
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

    def disconnect(self, server, body, after_s=0.5, then=None):
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

        request_id, m, stop_s, before = self.disconnect(server, body, then=next_request)
        self.assertEqual('client gone', m.group(2))
        started, records, pairs_started, pairs = (int(m.group(i)) for i in (3, 4, 5, 6))
        self.assertEqual(self.single_records, records)
        self.assertLess(started, records)
        self.assertEqual((1, 1), (pairs_started, pairs))
        self.assertEqual(200, second['r'].status_code)
        self.assertEqual(small_ref.text, second['r'].text)
        second_s = second['s']
        note(f'single graph -p 1 --threads-each 2: full {full_s:.2f} s; stopped {stop_s:.2f} s '
             f'after the close ({started} of {records} records started); the next request '
             f'answered {second_s:.2f} s after the close')
        # the stop comes within a record or two of the close, not at the request's end, and the
        # handler is free then: the next request is answered while the first would still run
        self.assertLess(0.5 + stop_s, 0.6 * full_s)
        self.assertLess(0.5 + second_s, 0.6 * full_s)
        self.stopped_quietly(server, before, request_id)
        # and the undisturbed request is answered as before, byte for byte
        r, _ = self.timed(server, body)
        self.assertEqual(reference, r.text)

    def test_multi_graph_stop_states_the_pairs(self):
        """staging's flags (-p 2 --threads-each 2): the pairs not started stay so"""
        server = self.server(multi=True, p=2, e=2)
        body = self.body(self.multi_records, multi=True)
        reference, full_s = self.timed(server, body)
        request_id, m, stop_s, before = self.disconnect(server, body)
        self.assertEqual('client gone', m.group(2))
        started, records, pairs_started, pairs = (int(m.group(i)) for i in (3, 4, 5, 6))
        self.assertEqual(self.multi_records, records)
        self.assertEqual(len(self.names), pairs)
        self.assertLess(pairs_started, pairs)
        self.assertGreaterEqual(pairs_started, 1)
        self.assertLessEqual(started, records)
        note(f'multi graph -p 2 --threads-each 2, {pairs} pairs: full {full_s:.2f} s; stopped '
             f'{stop_s:.2f} s after the close ({started} of {records} records, {pairs_started} '
             f'of {pairs} pairs started)')
        self.assertLess(0.5 + stop_s, 0.6 * full_s)
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
        time.sleep(0.5)
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

    def test_cache_waiter_gets_the_whole_answer(self):
        """--cache-size above 0: a duplicate that waits for the producer's entry keeps it going
        when the producer's client leaves; the duplicate gets the answer byte for byte"""
        reference, _ = self.single_reference()
        body = self.body(self.single_records)
        server = self.server(['--cache-size', '1'], p=2, e=2)
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
        # the duplicate's handler parses its body and attaches to the entry within this
        time.sleep(1.0)
        sock.close()
        waiter.join(timeout=300)
        r = answer['r']
        self.assertEqual(200, r.status_code, r.text[:2000])
        self.assertEqual(reference, r.text)
        text = server.text()[before:]
        self.assertIsNone(STOP.search(text), text[-3000:])
        self.assertRegex(text, rf'Request {producer}: cache STORE')
        dup = self.request_ids(server, before)[1]
        self.assertRegex(text, rf'Request {dup}: cache HIT')

    def test_cache_late_duplicate_computes(self):
        """--cache-size above 0: a stopped producer (nobody waited) leaves the cache; an identical
        request after it computes the answer itself"""
        reference, full_s = self.single_reference()
        body = self.body(self.single_records)
        server = self.server(['--cache-size', '1'], p=2, e=2)
        request_id, m, stop_s, before = self.disconnect(server, body)
        self.assertEqual('client gone', m.group(2))
        self.assertLess(int(m.group(3)), int(m.group(4)))
        after = len(server.text())
        r, late_s = self.timed(server, body)
        self.assertEqual(reference, r.text)
        text = server.text()[after:]
        late = self.request_ids(server, after)[0]
        self.assertRegex(text, rf'Request {late}: cache STORE')
        self.assertNotRegex(text, rf'Request {late}: cache HIT')
        note(f'cache: stopped {stop_s:.2f} s after the close; the late duplicate computed in '
             f'{late_s:.2f} s (undisturbed {full_s:.2f} s)')

    def stop_server(self, server):
        server.stop()
        self.servers.remove(server)


if __name__ == '__main__':
    unittest.main()
