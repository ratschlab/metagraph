"""TraverseClient against a local HTTP server: explicit Accept-Encoding, gzip and
deflate bodies, non-2xx -> TraverseError, 503 -> ServerInitializing, per-seed errors,
deepening through next_request; and an injected requests-like session."""

import copy
import gzip
import http.server
import json
import os
import sys
import threading
import unittest
import zlib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402

from metagraph.traverse import (  # noqa: E402
    AttemptExpired, ServerInitializing, TraverseClient, TraverseError,
)


class _Handler(http.server.BaseHTTPRequestHandler):
    def log_message(self, *args):
        pass

    def _send(self, status, obj, extra=None):
        data = json.dumps(obj, separators=(',', ':')).encode()
        accept = self.headers.get('Accept-Encoding', '')
        headers = {'Content-Type': 'application/json'}
        enc = self.server.encoding if self.server.encoding in accept else None
        if enc == 'gzip':
            data = gzip.compress(data)
        elif enc == 'deflate':
            data = zlib.compress(data)
        if enc:
            headers['Content-Encoding'] = enc
        headers.update(extra or {})
        self.send_response(status)
        for k, v in headers.items():
            self.send_header(k, v)
        self.send_header('Content-Length', str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        self.server.seen.append(('GET', self.path, dict(self.headers), None))
        if self.server.initializing:
            return self._send(503, {'error': 'Server is currently initializing'},
                              {'Retry-After': '60'})
        if self.path == '/api/traverse/capabilities':
            return self._send(200, {'k': 5, 'release': '', 'graphlet_format': 1})
        if self.path.startswith('/api/traverse/capabilities?'):
            # a multi-graph server's probe of one pair
            from urllib.parse import parse_qs, urlsplit
            q = parse_qs(urlsplit(self.path).query)
            return self._send(200, {'graph': q['graph'][0],
                                    'graph_path': q.get('graph_path', ['/g/' + q['graph'][0]])[0]})
        if self.path == '/api/capabilities':
            return self._send(200, {'mode': 'multi', 'graphs': ['A', 'B'], 'feature_level': 3})
        self._send(404, {'error': 'no route'})

    def do_POST(self):
        body = json.loads(self.rfile.read(int(self.headers['Content-Length'])))
        self.server.seen.append(('POST', self.path, dict(self.headers), body))
        if self.path == '/api/resolve':
            return self._send(200, {'k': 5, 'labels': [{'label': 'acc1', 'kind': 'header',
                                                         'kmers_supported': 6}],
                                    'candidates': []})
        if self.path != '/api/traverse':
            return self._send(404, {'error': 'no route'})
        seq = body['seeds'][0]['sequence']
        if seq == 'BAD':
            return self._send(400, {'error': "seed 'x': invalid character 'B'"})
        if 'not_after_ms' in body and body['not_after_ms'] < 1000:
            return self._send(409, {'error': 'not_after_ms 1 has passed', 'state': 'expired',
                                    'not_after_ms': body['not_after_ms'],
                                    'server_time_ms': 1791137002417,
                                    'server_instance': 'f' * 16})
        if body.get('attempt_id') == 'dup':
            return self._send(409, {'error': 'attempt_id dup is running', 'attempt': {}})
        resp = copy.deepcopy(T.doc_json(self.server.fixture, 'graphlet'))
        result = resp['results'][0]
        resp['results'] = []
        for s in body['seeds']:
            if s['sequence'] == 'NOCARRIER':
                resp['results'].append({'seed': {'seed_id': '', 'length_bp': 9,
                                                 'labels_from_seed': True},
                                        'error': 'no label carries every k-mer'})
            else:
                resp['results'].append(result)
        self._send(200, resp)


class TestClient(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.server = http.server.ThreadingHTTPServer(('127.0.0.1', 0), _Handler)
        cls.server.seen = []
        cls.server.encoding = 'gzip'
        cls.server.initializing = False
        cls.server.fixture = 'switch_chain'
        cls.thread = threading.Thread(target=cls.server.serve_forever, daemon=True)
        cls.thread.start()
        cls.client = TraverseClient('127.0.0.1', cls.server.server_address[1], api_path='api',
                                    timeout=10)

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()

    def setUp(self):
        self.server.seen.clear()
        self.server.encoding = 'gzip'
        self.server.initializing = False

    def test_capabilities_asks_for_compression(self):
        self.assertEqual(1, self.client.capabilities()['graphlet_format'])
        self.assertEqual('gzip, deflate', self.server.seen[0][2]['Accept-Encoding'])

    def test_pass5_capabilities_not_after_and_expired(self):
        """the per-graph probe (?graph=, graph_path percent-encoded), the server-wide document,
        not_after_ms passed through, and the expired 409 told apart from a duplicate's"""
        got = self.client.capabilities(graph='A', graph_path='/data/x y/g.dbg')
        self.assertEqual(('A', '/data/x y/g.dbg'), (got['graph'], got['graph_path']))
        self.assertEqual('/api/traverse/capabilities?graph=A&graph_path=%2Fdata%2Fx%20y%2Fg.dbg',
                         self.server.seen[-1][1])
        multi = TraverseClient('127.0.0.1', self.server.server_address[1], api_path='api',
                               timeout=10, graph='B')
        self.assertEqual('B', multi.capabilities()['graph'])
        self.assertEqual(['A', 'B'], self.client.server_capabilities()['graphs'])
        req = self.client.build_request(['ACGT'], attempt_id='a-1', not_after_ms=1791137060000)
        self.assertEqual(1791137060000, req['not_after_ms'])
        self.assertNotIn('not_after_ms', self.client.build_request(['ACGT']))
        with self.assertRaises(AttemptExpired) as cm:
            self.client.traverse(['AAAAACCCCC'], {'direction': 'right'}, not_after_ms=1)
        self.assertEqual((409, 1), (cm.exception.status, cm.exception.body['not_after_ms']))
        self.assertIsNone(cm.exception.usage)
        self.assertEqual(1, self.server.seen[-1][3]['not_after_ms'])
        with self.assertRaises(TraverseError) as cm:
            self.client.traverse(['AAAAACCCCC'], {'direction': 'right'}, attempt_id='dup')
        self.assertNotIsInstance(cm.exception, AttemptExpired)
        self.assertEqual(409, cm.exception.status)

    def test_traverse_builds_a_graphlet_request(self):
        resp = self.client.traverse(['AAAAACCCCC'], {'direction': 'right'})
        _, path, headers, body = self.server.seen[0]
        self.assertEqual('graphlet', body['strategy']['output']['detail'])
        self.assertTrue(body['strategy']['output']['timing'])
        self.assertEqual([{'sequence': 'AAAAACCCCC'}], body['seeds'])
        self.assertEqual('application/json', headers['Content-Type'])
        (g,) = resp.graphlets
        self.assertTrue(g.has_envelope)
        self.assertEqual('traverse-0.2', resp.envelope['algorithm_version'])
        self.assertEqual(T.doc_text('switch_chain'), g.dump())

    def test_deflate(self):
        self.server.encoding = 'deflate'
        resp = self.client.traverse([{'sequence': 'AAAAACCCCC', 'labels': ['chrA']}])
        self.assertEqual(1, len(resp.graphlets))

    def test_errors(self):
        with self.assertRaises(TraverseError) as e:
            self.client.traverse(['BAD'])
        self.assertEqual(400, e.exception.status)
        self.assertIn('invalid character', e.exception.message)
        self.server.initializing = True
        with self.assertRaises(ServerInitializing) as e:
            self.client.capabilities()
        self.assertEqual((503, 60), (e.exception.status, e.exception.retry_after))

    def test_a_failed_seed_is_a_result_not_an_exception(self):
        resp = self.client.traverse(['AAAAACCCCC', 'NOCARRIER'])
        self.assertIsNotNone(resp.graphlets[0])
        self.assertIsNone(resp.graphlets[1])
        self.assertEqual([1], [e['index'] for e in resp.errors])

    def test_deepen_posts_the_next_request(self):
        g = self.client.traverse(['AAAAACCCCC']).graphlets[0]
        self.server.seen.clear()
        self.client.deepen(g, 'right', [1], bp=40)
        body = self.server.seen[0][3]
        # walk 1 of switch_chain ends at the radius under C, 20 bp after its switch
        self.assertEqual([{'sequence': 'TTACTCGTAGCCGGGCGTGA', 'labels': ['C']}], body['seeds'])
        self.assertEqual(40, body['strategy']['bounds']['max_extension_bp'])

    def test_resolve(self):
        out = self.client.resolve('ACGTACGT', labels=['acc1'])
        self.assertEqual('acc1', out['labels'][0]['label'])
        self.assertEqual({'sequence': 'ACGTACGT', 'labels': ['acc1']}, self.server.seen[0][3])


class _FakeResponse:
    def __init__(self, status, obj):
        self.status_code = status
        self.content = json.dumps(obj).encode()
        self.headers = {'Content-Type': 'application/json'}


class _FakeSession:
    def __init__(self):
        self.calls = []

    def request(self, method, url, data=None, headers=None, timeout=None):
        self.calls.append((method, url, headers, timeout))
        if url.endswith('/traverse'):
            return _FakeResponse(200, T.doc_json('fork', 'graphlet'))
        return _FakeResponse(400, {'error': 'nope'})


class TestInjectedSession(unittest.TestCase):
    def test_session(self):
        s = _FakeSession()
        c = TraverseClient('h', 1, session=s, timeout=5, release='r1', graph='uhgg')
        resp = c.traverse(['ACGT'])
        self.assertEqual(2, len(resp.graphlets[0].arms))
        method, url, headers, timeout = s.calls[0]
        self.assertEqual(('POST', 'http://h:1/traverse', 'gzip, deflate', 5),
                         (method, url, headers['Accept-Encoding'], timeout))
        req = c.build_request(['ACGT'])
        self.assertEqual(('r1', 'uhgg'), (req['release'], req['graph']))
        with self.assertRaises(TraverseError):
            c.capabilities()


if __name__ == '__main__':
    unittest.main()
