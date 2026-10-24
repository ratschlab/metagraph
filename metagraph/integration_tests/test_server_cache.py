from pathlib import Path
import re
import subprocess
import tempfile
import time

import requests

from base import METAGRAPH, TEST_DATA_DIR, TestingBase


class TestServerCache(TestingBase):
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
