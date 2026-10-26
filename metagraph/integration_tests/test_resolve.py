"""Resolve/select through the normal executable, HTTP server and Python clients."""
import json
import os
from pathlib import Path
import random
import socket
import subprocess
import time
from unittest.mock import patch
import requests

from parameterized import parameterized_class
from metagraph.client import GraphClientJson, GraphClient
from base import TestingBase, METAGRAPH, PROTEIN_MODE

MODES = ['basic'] if PROTEIN_MODE else ['basic', 'canonical', 'primary']


def unused_port():
    with socket.socket() as sock:
        sock.bind(('127.0.0.1', 0))
        return sock.getsockname()[1]


def start_server(graph, annotation, log_path):
    for attempt in range(3):
        port = unused_port()
        with log_path.with_suffix(f'.{attempt}.log').open('w+') as log:
            server = subprocess.Popen([METAGRAPH, 'server_query', '-i', graph,
                                       '-a', annotation, '--mmap', '--address', '127.0.0.1',
                                       '--port', str(port), '-p', '1'],
                                      stdout=log, stderr=subprocess.STDOUT)
            client = GraphClientJson('127.0.0.1', port)
            deadline = time.monotonic() + 20
            while time.monotonic() < deadline and server.poll() is None:
                try:
                    client.stats(timeout=1)
                    if server.poll() is None:
                        return server, port
                except (RuntimeError, requests.RequestException):
                    pass
                time.sleep(0.05)
            if server.poll() is None:
                server.terminate()
                server.wait(timeout=10)
                raise AssertionError('Local test server did not become ready')
            log.seek(0)
            error = log.read()
            # Another test can bind after unused_port releases its socket.
            if 'bind: Address already in use' not in error or attempt == 2:
                raise AssertionError(error)


def runs(mask):
    result = []
    for i, match in enumerate(mask):
        if match:
            if result and result[-1][1] == i:
                result[-1][1] = i + 1
            else:
                result.append([i, i + 1])
    return result


@parameterized_class(('mode', 'representation'),
                     [(m, a) for m in MODES for a in ['column', 'brwt']])
class TestResolve(TestingBase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.k = 7
        rng = random.Random(670)
        cls.query = ''.join(rng.choice('ACGT') for _ in range(150))
        cls.records = {'A': cls.query[:90], 'B': cls.query[45:], 'C': 'A' * 30}
        files = []
        for label, sequence in cls.records.items():
            path = Path(cls.tempdir.name, label + '.fa')
            path.write_text(f'>{label}\n{sequence}\n')
            files.append(str(path))
        cls.graph = cls.tempdir.name + '/graph.dbg'
        cls.annotation = cls.tempdir.name + '/annotation.' + cls.representation + '.annodbg'
        cls._build_graph(files, cls.graph, cls.k, 'succinct', mode=cls.mode)
        cls._annotate_graph(files, cls.graph, cls.tempdir.name + '/annotation', cls.representation)
        cls.request = {'sequence': cls.query, 'discover': {'max_labels': 100},
                       'select': {'policy': 'max_support', 'max_seeds': 2,
                                  'merge_overlapping': False}}
        cls.request_path = Path(cls.tempdir.name, 'request.json')
        cls.request_path.write_text(json.dumps(cls.request))
        cls.server, cls.port = start_server(cls.graph, cls.annotation,
                                           Path(cls.tempdir.name, 'server.log'))
        cls.client = GraphClientJson('127.0.0.1', cls.port)

    @classmethod
    def tearDownClass(cls):
        cls.server.terminate()
        cls.server.wait(timeout=10)
        super().tearDownClass()

    def cli(self, mmap=True, paths=None):
        command = [METAGRAPH, 'traverse', '--resolve', '-i', self.graph, '-a', self.annotation]
        if mmap:
            command.append('--mmap')
        return subprocess.run(command + [str(p) for p in (paths or [self.request_path])],
                              capture_output=True, text=True, timeout=30)

    def test_cli_http_clients_and_literal_membership(self):
        expected = {}
        for name, sequence in self.records.items():
            kmers = {sequence[i:i+self.k] for i in range(len(sequence)-self.k+1)}
            if self.mode != 'basic':
                kmers |= {k.translate(str.maketrans('ACGT', 'TGCA'))[::-1] for k in list(kmers)}
            mask = [self.query[i:i+self.k] in kmers for i in range(len(self.query)-self.k+1)]
            if any(mask):
                expected[name] = {'count': sum(mask), 'runs': runs(mask)}
        outputs = []
        for mmap in [False, True]:
            proc = self.cli(mmap)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            outputs.append(json.loads(proc.stdout))
        outputs.append(self.client.resolve(**self.request))
        outputs.append(GraphClient('127.0.0.1', self.port).resolve(**self.request))
        for result in outputs:
            self.assertEqual(result, outputs[0])
            self.assertEqual(result['support'], 'kmer')
            self.assertEqual(result['regime'], self.mode)
            got = {x['label']: {'count': x['kmers_supported'], 'runs': x['runs']}
                   for x in result['labels']}
            self.assertEqual(got, expected)
            for seed in result['selection']['seeds']:
                a, b = seed['kmer_interval']
                self.assertEqual(seed['sequence'], self.query[a:b+self.k-1])
                for label in seed['labels']:
                    self.assertTrue(any(x <= a and b <= y for x, y in expected[label]['runs']))

    def test_limits_and_unsupported_requests(self):
        limited = self.client.resolve(self.query, discover={'max_labels': 1})
        self.assertEqual(len(limited['labels']), 1)
        self.assertEqual(limited['labels_truncated']['total'], 2)
        for options in [{'support': 'trace'}, {'discover': {'kind': 'header'}},
                        {'bounds': {'time_budget_ms': 10}}, {'min_block_kmers': -1}]:
            request = dict(self.request, **options)
            with self.assertRaises(RuntimeError):
                self.client.resolve(**request)
        bad = Path(self.tempdir.name, 'bad.json')
        bad.write_text('{"sequence":"A","sequence":"C"}')
        proc = self.cli(paths=[bad, self.request_path])
        self.assertEqual(proc.returncode, 1)
        results = [json.loads(line) for line in proc.stdout.splitlines()]
        self.assertIn('error', results[0])
        self.assertIn('selection', results[1])

    def test_explicit_selection_validation(self):
        seed = {'kmer_interval': [0, 3], 'labels': ['A']}
        for policy in [None, 'max_support', 'longest_first']:
            for seeds in [[], [seed]]:
                select = {'seeds': seeds}
                if policy is not None:
                    select['policy'] = policy
                request = dict(self.request, select=select)
                with self.assertRaisesRegex(RuntimeError, 'seeds requires explicit'):
                    self.client.resolve(**request)
                path = Path(self.tempdir.name, 'invalid-select.json')
                path.write_text(json.dumps(request))
                result = self.cli(paths=[path])
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('seeds requires explicit', json.loads(result.stdout)['error'])
        select = {'policy': 'explicit', 'seeds': [seed], 'min_block_bp': self.k + 3}
        with self.assertRaisesRegex(RuntimeError, 'shorter than min_block_bp'):
            self.client.resolve(self.query, labels=['A'], select=select)
        select['min_block_bp'] -= 1
        result = GraphClient('127.0.0.1', self.port).resolve(
            self.query, labels=['A'], select=select, support='kmer',
            min_block_kmers=2, run_format='rle', bounds={'max_query_bp': len(self.query)})
        self.assertEqual(result['selection']['seeds'][0]['sequence'], self.query[:self.k + 2])
        self.assertIsInstance(result['labels'][0]['runs'], str)

    def test_missing_resolve_exits_without_reading_stdin(self):
        command = [METAGRAPH, 'traverse', '-i', self.graph, '-a', self.annotation]
        with subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, text=True) as process:
            try:
                self.assertNotEqual(process.wait(timeout=5), 0)
                self.assertIn('traverse requires --resolve', process.stderr.read())
            finally:
                if process.poll() is None:
                    process.terminate()
                    process.wait(timeout=10)

    def test_server_start_retries_a_bind_collision(self):
        with socket.socket() as occupied:
            occupied.bind(('127.0.0.1', 0))
            occupied.listen()
            ports = [occupied.getsockname()[1], unused_port()]
            with patch(__name__ + '.unused_port', side_effect=ports) as allocate:
                server, port = start_server(self.graph, self.annotation,
                                            Path(self.tempdir.name, 'collision.log'))
                try:
                    self.assertEqual(allocate.call_count, 2)
                    self.assertEqual(port, ports[1])
                    self.assertIn('graph', GraphClientJson('127.0.0.1', port).stats())
                finally:
                    server.terminate()
                    server.wait(timeout=10)
