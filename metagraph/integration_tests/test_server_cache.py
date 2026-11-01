from contextlib import contextmanager
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib
import json
import threading
import re
import subprocess
import tempfile
import time

import requests

from base import METAGRAPH, TEST_DATA_DIR, TestingBase


class TestServerCache(TestingBase):
    def test_cached_search_preserves_gzip_and_other_json_routes(self):
        graph, annotation = self.tiny_index()
        payload = {'FASTA': '>q\nACGTACGTACGT\n', 'top_labels': 1}
        with self.running_server('-i', graph, '-a', annotation) as (client, url, log):
            first = client.post(url, json=payload, headers={'Accept-Encoding': 'identity'}, timeout=30)
            compressed = client.post(url, json=payload, headers={'Accept-Encoding': 'gzip'}, timeout=30)
            self.assertEqual(compressed.status_code, 200, compressed.text)
            self.assertEqual(compressed.headers.get('Content-Encoding'), 'gzip')
            self.assertEqual(first.json(), compressed.json())
            self.assertIsInstance(client.get(url.replace('/search', '/stats'), timeout=5).json(), dict)
            self.assertIsInstance(client.get(url.replace('/search', '/column_labels'), timeout=5).json(), list)
            aligned = client.post(url.replace('/search', '/align'), json=payload, timeout=30)
            self.assertEqual(aligned.status_code, 200, aligned.text)
            self.assertIsInstance(aligned.json(), list)
            self.assertEqual(log.read_text().count('cache STORE'), 1)
            self.assertEqual(log.read_text().count('cache HIT'), 1)

    def test_cache_preserves_configured_graph_selection_limit(self):
        graph, annotation = self.tiny_index()
        indexes = graph.with_name('indexes.csv')
        indexes.write_text(f'left,{graph},{annotation}\nright,{graph},{annotation}\n')
        payload = {'FASTA': '>q\nACGTACGTACGT\n', 'top_labels': 1}
        with self.running_server(indexes, '--max-graphs-without-selection', '1') as (client, url, log):
            selected = dict(payload, graphs=['left'])
            for _ in range(2):
                response = client.post(url, json=selected, timeout=30)
                self.assertEqual(response.status_code, 200, response.text)
            response = client.post(url, json=payload, timeout=30)
            self.assertEqual(response.status_code, 400, response.text)
            self.assertEqual(log.read_text().count('cache STORE'), 1)
            self.assertEqual(log.read_text().count('cache HIT'), 1)

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
        # twenty names of one pair: a pair listed twice under one name is searched once
        slow = [f'slow{i:02}' for i in range(20)]
        indexes.write_text(f'bad,{bad},{annotation}\n'
                           + ''.join(f'{name},{graph},{annotation}\n' for name in slow))
        with self.running_server(indexes, '--mmap', '--mem-cap-gb', '1') as (client, url, log):
            bad.unlink()  # Already mmap-loaded; only this request's in_ram reload fails.
            response = client.post(url, json={
                'FASTA': ''.join(f'>q{i}\nACGTACGTACGT\n' for i in range(100)),
                'graphs': ['bad'] + slow, 'in_ram': True}, timeout=30)
            self.assertEqual(response.status_code, 400)
            text = log.read_text()
            completed = re.findall(r'Loaded and searched graph .* in [0-9.e+-]+ sec', text)
            self.assertEqual(len(completed), 20, text[-3000:])
            self.assertLess(text.rfind('Loaded and searched graph'), text.rfind('Error on request'))

    def test_a_pair_under_two_names_is_answered_with_each_name(self):
        """SPEC-multigraph-envelope.md section 3 (E1): a multi-graph /search tags every row and
        every cut with the name its pair was selected under, so the result cache keys the names
        beside their pairs: one pair listed under two names is answered with each name, never
        with the other name's cached rows."""
        graph, annotation = self.tiny_index()
        indexes = graph.with_name('indexes.csv')
        indexes.write_text(f'right,{graph},{annotation}\nleft,{graph},{annotation}\n'
                           f'left,{graph},{annotation}\n')
        payload = {'FASTA': '>q\nACGTACGTACGT\n', 'top_labels': 1}
        tags = {'graph_path': str(graph), 'annotation_path': str(annotation)}

        def cut(name):
            return dict(tags, graph=name, index_fp=None, returned=1, complete=True,
                        nth_score=None, tied=0, tied_returned=0)

        with self.running_server(indexes, '--cache-size', '0.01') as (client, url, log):
            for names in (['left'], ['right'], ['left'], ['right', 'left']):
                with self.subTest(graphs=names):
                    response = client.post(url, json=dict(payload, graphs=names), timeout=30)
                    self.assertEqual(response.status_code, 200, response.text)
                    [answer] = response.json()
                    # in the index's pair order: the names in byte order, a pair listed twice
                    # under one name once
                    self.assertEqual([dict(tags, sample='sample', kmer_count=7, graph=name)
                                      for name in sorted(names)], answer['results'])
                    self.assertEqual([cut(name) for name in sorted(names)], answer['cuts'])
            text = log.read_text()
            self.assertEqual(text.count('cache STORE'), 3, text[-3000:])
            self.assertEqual(text.count('cache HIT'), 1, text[-3000:])

    def test_single_graph_search_states_one_cut(self):
        """SPEC-multigraph-envelope.md section 3 (E1): a single-graph server leaves its rows
        untagged and states one cut per query sequence, its pair tags null and its index_fp the
        server's: null without --index-manifest, the manifest's digest with it."""
        graph, annotation = self.tiny_index()
        payload = {'FASTA': '>q\nACGTACGTACGT\n>none\nTTTTTTTTTTTT\n', 'top_labels': 1}

        def cut(index_fp, returned):
            return {'graph': None, 'graph_path': None, 'annotation_path': None,
                    'index_fp': index_fp, 'returned': returned, 'complete': True,
                    'nth_score': None, 'tied': 0, 'tied_returned': 0}

        # the bundle's manifest, as scripts/traversal/build_mini_refseq.sh writes it
        files = []
        for path in sorted(graph.parent.iterdir()):
            if path.is_file() and path.name.startswith(('graph', 'annotation')):
                data = path.read_bytes()
                files.append({'path': path.name, 'size': len(data),
                              'sha256': hashlib.sha256(data).hexdigest()})
        index_fp = hashlib.sha256(''.join('%s\t%d\t%s\n' % (f['path'], f['size'], f['sha256'])
                                          for f in files).encode()).hexdigest()
        manifest = graph.with_name('index.manifest.json')
        manifest.write_text(json.dumps({'format': 'metagraph-index-manifest', 'version': 1,
                                        'index_fp': index_fp, 'files': files}))

        for flags, fp in (([], None), (['--index-manifest', manifest], index_fp)):
            with self.subTest(flags=flags), \
                    self.running_server('-i', graph, '-a', annotation, *flags) as (client, url, _):
                response = client.post(url, json=payload, timeout=30)
                self.assertEqual(response.status_code, 200, response.text)
                self.assertEqual([{'seq_description': 'q', 'cuts': [cut(fp, 1)],
                                   'results': [{'sample': 'sample', 'kmer_count': 7}]},
                                  {'seq_description': 'none', 'cuts': [cut(fp, 0)],
                                   'results': []}],
                                 response.json())

    def test_search_refuses_pairs_on_every_server(self):
        """/search selects by `graphs` only (SPEC-multigraph-envelope.md section 6): `pairs` is
        refused 400 pairs_unsupported on a single-graph server too, where it used to be
        ignored (GPT B4), cold and with the same request otherwise cached: the check comes
        before the cache. `graphs` keeps its refusal of before."""
        graph, annotation = self.tiny_index()
        payload = {'FASTA': '>q\nACGTACGTACGT\n', 'top_labels': 1}
        refusal = {'error': 'request: /search selects by graphs only (it has no re-ask): '
                            'pairs is not accepted', 'code': 'pairs_unsupported'}
        indexes = graph.with_name('indexes.csv')
        indexes.write_text(f'one,{graph},{annotation}\n')
        for server in (('-i', graph, '-a', annotation), (indexes,)):
            with self.subTest(server=server), \
                    self.running_server(*server) as (client, url, log):
                for pairs in ([{'graph': 'one'}], [], 'x'):
                    # cold
                    response = client.post(url, json=dict(payload, pairs=pairs), timeout=30)
                    self.assertEqual(400, response.status_code, response.text)
                    self.assertEqual(refusal, response.json())
                # warm: the request without pairs is cached, with them it is refused
                for _ in range(2):
                    response = client.post(url, json=payload, timeout=30)
                    self.assertEqual(200, response.status_code, response.text)
                response = client.post(url, json=dict(payload, pairs=[{'graph': 'one'}]),
                                       timeout=30)
                self.assertEqual(400, response.status_code, response.text)
                self.assertEqual(refusal, response.json())
                self.assertEqual(log.read_text().count('cache STORE'), 1)
                self.assertEqual(log.read_text().count('cache HIT'), 1)
        # a single-graph server's refusal of `graphs` comes first, as before
        with self.running_server('-i', graph, '-a', annotation) as (client, url, _):
            response = client.post(url, json=dict(payload, graphs=['one'], pairs=[]),
                                   timeout=30)
            self.assertEqual(400, response.status_code, response.text)
            self.assertEqual({'error': 'Bad request: no support for filtering graphs on this '
                                       'server'}, response.json())

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
