from contextlib import contextmanager
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import threading
import re
import subprocess
import tempfile
import time

import requests

from base import METAGRAPH, TEST_DATA_DIR, TestingBase


class TestServerCache(TestingBase):
    @contextmanager
    def running_server(self, *args):
        with tempfile.NamedTemporaryFile(mode='w+') as log, requests.Session() as client:
            client.trust_env = False
            process = subprocess.Popen([METAGRAPH, 'server_query', *map(str, args),
                                        '--address', '127.0.0.1', '--port', '0', '-p', '4', '-v'],
                                       stdout=log, stderr=log)
            try:
                deadline = time.monotonic() + 30
                while time.monotonic() < deadline:
                    self.assertIsNone(process.poll(), 'Server exited during startup')
                    match = re.search(r'\[Server\] Listening on port (\d+)', Path(log.name).read_text())
                    if match:
                        url = f'http://127.0.0.1:{match.group(1)}'
                        if client.get(url + '/stats', timeout=5).status_code == 200:
                            break
                    time.sleep(0.05)
                else:
                    self.fail('Server did not become ready')
                yield client, url + '/search', Path(log.name)
            finally:
                if process.poll() is None:
                    process.terminate()
                try:
                    process.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait(timeout=10)
                    raise

    def tiny_index(self):
        root = Path(self.tempdir.name) / self._testMethodName
        root.mkdir()
        fasta = root / 'input.fa'
        fasta.write_text('>sample\nACGTACGTACGT\n')
        graph = root / 'graph.dbg'
        self._build_graph(str(fasta), str(graph), 6, 'succinct')
        self._annotate_graph(str(fasta), str(graph), str(root / 'annotation'), 'column',
                             extra_params='--anno-label sample')
        return graph, root / 'annotation.column.annodbg'

    def test_hit_compression_oversized_response_and_disabled_cache(self):
        graph, annotation = self.tiny_index()
        payload = {'FASTA': '>q\nACGTACGTACGT\n', 'top_labels': 1}
        with self.running_server('-i', graph, '-a', annotation, '--cache-size', '0.000004096') as (client, url, log):
            first = client.post(url, json=payload, headers={'Accept-Encoding': 'identity'}, timeout=30)
            self.assertEqual(first.status_code, 200, first.text)
            hit = client.post(url, json=payload, headers={'Accept-Encoding': 'deflate'}, timeout=30)
            self.assertEqual(hit.status_code, 200, hit.text)
            self.assertEqual(hit.headers.get('Content-Encoding'), 'deflate')
            self.assertEqual(first.json(), hit.json())
            oversized = dict(payload, FASTA='>' + 'x' * 6000 + '\nACGTACGTACGT\n')
            large = client.post(url, json=oversized, timeout=30)
            self.assertEqual(large.status_code, 200, large.text)
            again = client.post(url, json=payload, timeout=30)
            self.assertEqual(again.json(), first.json())
            text = log.read_text()
            self.assertEqual(text.count('cache STORE'), 1, text[-3000:])
            self.assertEqual(text.count('cache HIT'), 2, text[-3000:])
            self.assertEqual(text.count('cache BYPASS'), 1, text[-3000:])
        with self.running_server('-i', graph, '-a', annotation, '--cache-size', '0') as (client, url, log):
            for _ in range(2):
                response = client.post(url, json=payload, timeout=30)
                self.assertEqual(response.status_code, 200)
                self.assertEqual(response.json(), first.json())
            self.assertNotRegex(log.read_text(), r'cache (HIT|STORE|BYPASS)')

    def test_identical_concurrent_requests_share_one_result(self):
        graph, annotation = self.tiny_index()
        payload = {'FASTA': ''.join(f'>q{i}\nACGTACGTACGT\n' for i in range(200)), 'top_labels': 1}
        with self.running_server('-i', graph, '-a', annotation) as (_, url, log):
            barrier = threading.Barrier(4)
            def post():
                with requests.Session() as client:
                    client.trust_env = False
                    barrier.wait(timeout=10)
                    response = client.post(url, json=payload, timeout=30)
                    self.assertEqual(response.status_code, 200)
                    return response.json()
            with ThreadPoolExecutor(max_workers=4) as pool:
                results = list(pool.map(lambda _: post(), range(4)))
            self.assertTrue(all(result == results[0] for result in results))
            self.assertEqual(log.read_text().count('cache STORE'), 1)
            self.assertEqual(log.read_text().count('cache HIT'), 3)

    def test_graph_error_waits_for_other_queries(self):
        graph, annotation = self.tiny_index()
        bad = graph.with_name('bad.dbg')
        bad.symlink_to(graph)
        indexes = graph.with_name('indexes.csv')
        indexes.write_text(f'bad,{bad},{annotation}\n' + f'slow,{graph},{annotation}\n' * 20)
        with self.running_server(indexes, '--mmap', '--mem-cap-gb', '1') as (client, url, log):
            bad.unlink()  # Already mmap-loaded; only this request's in_ram reload fails.
            response = client.post(url, json={
                'FASTA': ''.join(f'>q{i}\nACGTACGTACGT\n' for i in range(100)),
                'graphs': ['bad', 'slow'], 'in_ram': True}, timeout=30)
            self.assertEqual(response.status_code, 400)
            text = log.read_text()
            completed = re.findall(r'Loaded and searched graph .* in [0-9.e+-]+ sec', text)
            self.assertEqual(len(completed), 20, text[-3000:])
            self.assertLess(text.rfind('Loaded and searched graph'), text.rfind('Error on request'))

    def test_graph_paths_with_delimiters_have_distinct_cache_entries(self):
        root = Path(self.tempdir.name)
        fasta = TEST_DATA_DIR + '/transcripts_100.fa'
        graph = root / 'graph.dbg'
        self._build_graph(fasta, str(graph), 6, 'succinct')
        for label in ['left', 'right']:
            self._annotate_graph(fasta, str(graph), str(root / label), 'column',
                                 extra_params=f'--anno-label {label}')

        left_graph = root / 'left.dbg'
        right_annotation = root / 'right-alias.column.annodbg'
        middle = root / 'middle.dbg'
        right_graph = Path(str(left_graph) + '|' + str(middle))
        left_annotation = Path(str(middle) + '|' + str(right_annotation))
        # Both pairs have the same old encoding: graph + '|' + annotation.
        self.assertEqual(f'{left_graph}|{left_annotation}',
                         f'{right_graph}|{right_annotation}')
        for target, source in [(left_graph, graph), (right_graph, graph),
                               (left_annotation, root / 'left.column.annodbg'),
                               (right_annotation, root / 'right.column.annodbg')]:
            target.parent.mkdir(parents=True, exist_ok=True)
            target.symlink_to(source)
        index_list = root / 'indexes.csv'
        index_list.write_text(f'left,{left_graph},{left_annotation}\n'
                              f'right,{right_graph},{right_annotation}\n')

        with tempfile.NamedTemporaryFile(mode='w+') as log:
            process = subprocess.Popen([METAGRAPH, 'server_query', str(index_list),
                                        '--address', '127.0.0.1', '--port', '0', '-p', '2'],
                                       stdout=log, stderr=log)
            try:
                with requests.Session() as client:
                    client.trust_env = False
                    deadline = time.monotonic() + 30
                    while time.monotonic() < deadline:
                        self.assertIsNone(process.poll(), 'Server exited during startup')
                        match = re.search(r'\[Server\] Listening on port (\d+)',
                                          Path(log.name).read_text())
                        if match:
                            break
                        time.sleep(0.05)
                    else:
                        self.fail('Server did not become ready')
                    url = f'http://127.0.0.1:{match.group(1)}/search'
                    for label in ['left', 'right']:
                        response = client.post(url, json={
                            'FASTA': '>query\nCCTCTGTGGAATCCAATCTGTCTTCCATCCTGCGTGGCCGAGGG\n',
                            'discovery_fraction': 0.01, 'graphs': [label]}, timeout=30)
                        self.assertEqual(response.status_code, 200, response.text)
                        samples = {hit['sample'] for hit in response.json()[0]['results']}
                        self.assertEqual(samples & {'left', 'right'}, {label})
            finally:
                if process.poll() is None:
                    process.terminate()
                try:
                    process.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait(timeout=10)
                    raise
