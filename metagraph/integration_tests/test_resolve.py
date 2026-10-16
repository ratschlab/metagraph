"""Resolve/select through the normal executable, HTTP server and Python clients."""
import json
import os
from pathlib import Path
import random
import socket
import subprocess
import time
import requests

from parameterized import parameterized_class
from metagraph.client import GraphClientJson, GraphClient
from base import TestingBase, METAGRAPH, PROTEIN_MODE

MODES = ['basic'] if PROTEIN_MODE else ['basic', 'canonical', 'primary']


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
        with socket.socket() as s:
            s.bind(('127.0.0.1', 0))
            cls.port = s.getsockname()[1]
        cls.server_log = open(Path(cls.tempdir.name, 'server.log'), 'w')
        cls.server = subprocess.Popen([METAGRAPH, 'server_query', '-i', cls.graph,
                                      '-a', cls.annotation, '--mmap', '--address', '127.0.0.1',
                                      '--port', str(cls.port), '-p', '1'],
                                     stdout=cls.server_log, stderr=subprocess.STDOUT)
        cls.client = GraphClientJson('127.0.0.1', cls.port)
        deadline = time.monotonic() + 20
        while time.monotonic() < deadline:
            if cls.server.poll() is not None:
                raise AssertionError(Path(cls.tempdir.name, 'server.log').read_text())
            try:
                cls.client.stats(timeout=1)
                return
            except (RuntimeError, requests.RequestException):
                pass
            time.sleep(0.05)
        cls.server.terminate()
        cls.server.wait(timeout=10)
        raise AssertionError('Local test server did not become ready')

    @classmethod
    def tearDownClass(cls):
        cls.server.terminate()
        cls.server.wait(timeout=10)
        cls.server_log.close()
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
