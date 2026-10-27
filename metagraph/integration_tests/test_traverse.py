from contextlib import contextmanager
from copy import deepcopy
from pathlib import Path
import json
import re
import subprocess
import tempfile
import time
from urllib.parse import urlsplit

import requests

from base import METAGRAPH, TestingBase
from metagraph.client import GraphApiError, GraphClient, GraphClientJson


class TestTraversal(TestingBase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.root = Path(cls.tempdir.name)
        cls.fasta = cls.root / 'source.fa'
        cls.fasta.write_text('>A\nGGCTTACGTAAGCATG\n>B\nCCTGGACGTACCTAGG\n')
        cls.graph = cls.root / 'graph.dbg'
        cls.annotation = cls.root / 'annotation.column_coord.annodbg'
        cls.binary = cls.root / 'annotation.column.annodbg'
        cls.command('build', '-k', '5', '-o', cls.graph, cls.fasta)
        cls.command('annotate', '-i', cls.graph, '--anno-header', '--coordinates',
                    '-o', cls.root / 'annotation', cls.fasta)
        cls.command('transform_anno', '--coordinates', '--anno-type', 'column_coord',
                    '-o', cls.root / 'annotation', cls.binary)
        cls.headers = cls.root / 'headers.column_coord.annodbg'
        cls.command('annotate', '-i', cls.graph, '--anno-filename', '--coordinates',
                    '-o', cls.root / 'headers', cls.fasta)
        cls.command('annotate', '-i', cls.graph, '--anno-filename', '--index-header-coords',
                    '-o', cls.root / 'headers', cls.fasta)
        cls.command('transform_anno', '--coordinates', '--anno-type', 'column_coord',
                    '-o', cls.root / 'headers', cls.root / 'headers.column.annodbg')

    @staticmethod
    def command(*args):
        result = subprocess.run([METAGRAPH, *map(str, args)], capture_output=True, text=True, timeout=120)
        if result.returncode:
            raise AssertionError(f'{args}: {result.stdout}\n{result.stderr}')
        return result.stdout

    @contextmanager
    def server(self, *args):
        with tempfile.NamedTemporaryFile(mode='w+') as log, requests.Session() as client:
            client.trust_env = False
            process = subprocess.Popen([METAGRAPH, 'server_query', *map(str, args),
                                        '--address', '127.0.0.1', '--port', '0', '-p', '4', '-v'],
                                       stdout=log, stderr=log)
            try:
                deadline = time.monotonic() + 30
                while time.monotonic() < deadline:
                    text = Path(log.name).read_text()
                    self.assertIsNone(process.poll(), text)
                    match = re.search(r'\[Server\] Listening on port (\d+)', text)
                    if match:
                        url = f'http://127.0.0.1:{match.group(1)}'
                        if client.get(url + '/stats', timeout=5).status_code == 200:
                            self.server_log = Path(log.name)
                            yield client, url
                            return
                    time.sleep(0.05)
                self.fail('Traversal server did not become ready: ' + text)
            finally:
                if process.poll() is None:
                    process.terminate()
                try:
                    process.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait(timeout=10)
                    raise

    @staticmethod
    def request(direction='both', support='trace'):
        return {'seeds': [{'sequence': 'ACGTA', 'labels': ['A', 'B']}],
                'strategy': {'direction': direction, 'support': support,
                             'branching': {'on_reconverge': 'keep'},
                             'bounds': {'max_extension_bp': 5},
                             'output': {'timing': False}}}

    def post(self, client, url, route, body, status=200, **kwargs):
        response = client.post(url + route, json=body, timeout=30, **kwargs)
        self.assertEqual(response.status_code, status, response.text)
        return response.json()

    @staticmethod
    def paths(result, arm):
        labels = [label['name'] for label in result['label_dict']]
        return {path['continuation']['sequence']:
                {labels[label] for label in path['continuation']['labels']}
                for path in result['arms'][arm]['paths']}

    def test_cli_http_and_python_directions(self):
        expected = {'left': {'CCTGG': {'B'}, 'GGCTT': {'A'}},
                    'right': {'AGCAT': {'A'}, 'CCTAG': {'B'}}}
        for mmap in [False, True]:
            with self.server('-i', self.graph, '-a', self.annotation, *(['--mmap'] if mmap else [])) as (http, url):
                address = urlsplit(url)
                clients = [GraphClientJson(address.hostname, address.port),
                           GraphClient(address.hostname, address.port)]
                for client in clients:
                    caps = client.traversal_capabilities(timeout=10)
                    self.assertTrue(caps['supports_trace'])
                    profile = client.resolve('ACGTA', discover={}, select={'policy': 'longest_first'}, timeout=10)
                    self.assertEqual(profile['selection']['seeds'][0]['labels'], ['A', 'B'])
                    profile = client.resolve('ACGTA', labels='A', timeout=10)
                    self.assertEqual([label['label'] for label in profile['labels']], ['A'])
                for direction in ['left', 'right', 'both']:
                    request = self.request(direction)
                    path = self.root / f'{self._testMethodName}-{mmap}-{direction}.json'
                    path.write_text(json.dumps(request))
                    cli = json.loads(self.command('traverse', '-i', self.graph, '-a', self.annotation,
                                                  *(['--mmap'] if mmap else []), path))
                    output = self.post(http, url, '/traverse', request)
                    python = clients[0].traverse(request, timeout=10)
                    wanted = {'left', 'right'} if direction == 'both' else {direction}
                    for answer in [cli, output, python]:
                        result = answer['results'][0]
                        self.assertEqual(set(result['arms']), wanted)
                        for arm in wanted:
                            self.assertEqual(self.paths(result, arm), expected[arm])
                # The other arm still extends when one direction is a graph boundary.
                request = self.request()
                request['seeds'] = [{'sequence': 'GGCTT', 'labels': ['A']}]
                result = self.post(http, url, '/traverse', request)['results'][0]
                self.assertEqual(max(p['length_bp'] for p in result['arms']['left']['paths']), 0)
                self.assertGreater(max(p['length_bp'] for p in result['arms']['right']['paths']), 0)

    def test_headers_graphlet_and_strict_errors(self):
        with self.server('-i', self.graph, '-a', self.headers, '--mmap') as (http, url):
            resolved = self.post(http, url, '/resolve', {'sequence': 'ACGTA', 'support': 'trace',
                                                       'discover': {'kind': 'header'}})
            names = [label['label'] for label in resolved['labels']]
            self.assertEqual(set(names), {'A', 'B'})
            request = self.request()
            request['strategy']['labels'] = {'seed_label_kind': 'header'}
            request['strategy']['output']['coordinates'] = True
            request['strategy']['output']['detail'] = 'graphlet'
            response = http.post(url + '/traverse', json=request, timeout=30,
                                 headers={'Accept-Encoding': 'gzip;q=0, deflate;q=1'})
            self.assertEqual(response.status_code, 200, response.text)
            self.assertEqual(response.headers.get('Content-Encoding'), 'deflate')
            result = response.json()['results'][0]
            self.assertIn('graphlet', result)
            self.assertIn('coordinates', result)
            self.post(http, url, '/resolve', {'sequence': 'ACGTA', 'discover': {}, 'typo': 1}, 400)
            self.post(http, url, '/traverse', dict(request, in_ram=True), 400)
            self.post(http, url, '/resolve', {'sequence': 'ACGTA', 'discover': {},
                                             'bounds': {'time_budget_ms': 1}}, 400)
            malformed = http.post(url + '/resolve', data='{"sequence":"ACGTA","sequence":"C"}', timeout=10)
            self.assertEqual(malformed.status_code, 400)
        with self.server('-i', self.graph, '-a', self.binary) as (http, url):
            self.post(http, url, '/traverse', self.request(), 400)

    def test_attempt_lifecycle(self):
        with self.server('-i', self.graph, '-a', self.annotation, '--mmap') as (http, url):
            caps = http.get(url + '/capabilities', timeout=10).json()
            instance = caps['attempts']['server_instance']
            request = dict(self.request(), attempt_id='completed', expect_server_instance=instance)
            result = self.post(http, url, '/traverse', request)
            self.assertEqual(result['usage']['attempt_id'], 'completed')
            state = http.get(url + '/traverse/attempt/completed', timeout=10)
            self.assertEqual(state.status_code, 200, state.text)
            self.assertEqual(state.json()['state'], 'finished')
            self.post(http, url, '/traverse', request, 409)
            cancelled = self.post(http, url, '/traverse/cancel', {'attempt_id': 'before-start'}, 404)
            self.assertTrue(cancelled['tombstone'])
            self.post(http, url, '/traverse', dict(request, attempt_id='before-start'), 409)
            self.post(http, url, '/traverse', dict(request, attempt_id='expired', not_after_ms=1), 409)
            self.post(http, url, '/traverse', dict(request, attempt_id='other', expect_server_instance='other'), 409)
            bad = deepcopy(request)
            bad['attempt_id'] = 'bad'
            bad['strategy']['unknown'] = 1
            self.assertEqual(self.post(http, url, '/traverse', bad, 400)['usage']['attempt_id'], 'bad')
            address = urlsplit(url)
            client = GraphClient(address.hostname, address.port)
            colon_request = dict(self.request(), attempt_id='run:1')
            client.traverse(colon_request, timeout=10)
            self.assertEqual(client.traversal_attempt('run:1', timeout=10)['state'], 'finished')
            with self.assertRaises(GraphApiError) as error:
                client.cancel_traversal('python-before-start', timeout=10)
            self.assertEqual(error.exception.status_code, 404)
            self.assertTrue(error.exception.response['tombstone'])

    def test_multigraph_selection(self):
        csv = self.root / 'graphs.csv'
        csv.write_text(f'one,{self.graph},{self.annotation}\n'
                       f'one,{self.graph},{self.annotation}\n'
                       f'ambiguous,{self.graph},{self.annotation}\n'
                       f'ambiguous,{self.graph},{self.binary}\n')
        with self.server(csv, '--mmap') as (http, url):
            request = self.request()
            self.post(http, url, '/traverse', request, 400)
            one = self.post(http, url, '/traverse', dict(request, graph='one'))
            same = self.post(http, url, '/traverse', dict(request, graphs=['one']))
            self.assertEqual(one['results'], same['results'])
            self.post(http, url, '/traverse', dict(request, graph='one', graph_path='missing'), 400)
            self.post(http, url, '/traverse', dict(request, graph='ambiguous'), 400)
            self.post(http, url, '/traverse', dict(request, graph='one', graphs=['one']), 400)
            caps = http.get(url + '/traverse/capabilities', params={'graph': 'one'}, timeout=10)
            self.assertEqual(caps.status_code, 200, caps.text)
            self.assertEqual(caps.json()['graph'], 'one')

    def test_header_index_is_prepared_only_for_traversal(self):
        with self.server('-i', self.graph, '-a', self.headers, '--mmap') as (http, url):
            http.get(url + '/traverse/capabilities', timeout=10).raise_for_status()
            self.post(http, url, '/search', {'FASTA': '>query\nACGTA\n'})
            self.assertNotIn('Built header index', self.server_log.read_text())
            request = {'sequence': 'ACGTA', 'discover': {'kind': 'header'}}
            for _ in range(2):
                profile = self.post(http, url, '/resolve', request)
                self.assertGreaterEqual(profile['timing']['prepare_ms'], 0)
            self.assertEqual(self.server_log.read_text().count('Built header index'), 1)

    def test_record_mapping_rejects_duplicate_headers(self):
        fasta = self.root / 'duplicate.fa'
        fasta.write_text('>A\nGGCTTACGTAAGCATG\n>A\nCCTGGACGTACCTAGG\n')
        result = subprocess.run([METAGRAPH, 'annotate', '-i', str(self.graph),
                                 '--anno-filename', '--index-header-coords',
                                 '-o', str(self.root / 'duplicate'), str(fasta)],
                                capture_output=True, text=True, timeout=30)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('repeats record', result.stderr)

    def test_record_mapping_preserves_whole_headers(self):
        fasta = self.root / 'split-header.fa'
        fasta.write_text('>ACC|\nGGCTTACGTAAGCATG\n')
        result = subprocess.run([METAGRAPH, 'annotate', '-i', str(self.graph),
                                 '--anno-filename', '--index-header-coords',
                                 '--header-delimiter', '|',
                                 '-o', str(self.root / 'split-header'), str(fasta)],
                                capture_output=True, text=True, timeout=30)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('requires whole record headers', result.stderr)

    def test_index_names_follow_the_token_grammar(self):
        for name in ['', 'with space', 'with:colon']:
            result = subprocess.run([METAGRAPH, 'traverse', '--index-inventory',
                                     '-i', str(self.graph), '-a', str(self.annotation),
                                     '--index-name', name],
                                    capture_output=True, text=True, timeout=30)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('--index-name must match', result.stderr)
