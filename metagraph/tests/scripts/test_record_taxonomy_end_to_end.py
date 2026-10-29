"""Real CLI/preparation/HTTP checks; set METAGRAPH_BINARY to a DNA build to run."""
import contextlib
import gzip
import importlib.util
import json
import os
from pathlib import Path
import shutil
import socket
import subprocess
import sys
import tempfile
import time
import unittest
import urllib.error
import urllib.request


MG = os.environ.get('METAGRAPH_BINARY')
ROOT = Path(__file__).resolve().parents[2]
PREPARE = ROOT / 'scripts' / 'prepare_record_taxids.py'
MANIFEST = ROOT / 'scripts' / 'traversal' / 'index_manifest.py'
_spec = importlib.util.spec_from_file_location(
    'pattern_contract', ROOT / 'api' / 'python' / 'tests' / 'test_pattern_fixtures.py')
CONTRACT = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(CONTRACT)


def command(*args, ok=True):
    result = subprocess.run([str(a) for a in args], capture_output=True, text=True, timeout=60)
    if ok and result.returncode:
        raise AssertionError(f'{args}: {result.returncode}\n{result.stdout}\n{result.stderr}')
    return result


def http(url, body=None):
    payload = None if body is None else json.dumps(body).encode()
    request = urllib.request.Request(url, data=payload, headers={'Content-Type': 'application/json'})
    try:
        with urllib.request.urlopen(request, timeout=10) as response:
            return response.status, json.load(response)
    except urllib.error.HTTPError as error:
        with error:
            return error.code, json.load(error)


@contextlib.contextmanager
def server(graph, annotation, *extra):
    with socket.socket() as listener:
        listener.bind(('127.0.0.1', 0))
        port = listener.getsockname()[1]
    with tempfile.TemporaryFile() as log:
        process = subprocess.Popen([MG, 'server_query', '-i', str(graph), '-a', str(annotation),
                                    '--address', '127.0.0.1', '--port', str(port), '-p', '2',
                                    '--pattern-min-information-bits', '0', *map(str, extra)],
                                   stdout=log, stderr=subprocess.STDOUT)
        url = f'http://127.0.0.1:{port}'
        try:
            deadline = time.monotonic() + 20
            while process.poll() is None and time.monotonic() < deadline:
                try:
                    if http(url + '/stats')[0] == 200:
                        yield url
                        return
                except (OSError, urllib.error.URLError):
                    pass
                time.sleep(0.1)
            log.seek(0)
            raise AssertionError(f'Server did not become ready: {log.read().decode()}')
        finally:
            if process.poll() is None:
                process.terminate()
                process.wait(timeout=10)


@unittest.skipUnless(MG, 'set METAGRAPH_BINARY to a DNA build')
class RecordTaxonomyEndToEnd(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.folder = Path(cls.temporary.name)
        bins = {
            'bin-a.fa': '>NC_000001.1\nAAAAACCCCC\n>NC_000002.1\nGGGGGTTTTT\n>NC_000003.1\nCCCCCGGGGG\n',
            'bin-b.fa': '>NC_000004.1\nAAAAACCCCCGGGGG\n',
        }
        fastas = []
        for name, text in bins.items():
            path = cls.folder / name
            path.write_text(text)
            fastas.append(path)
        cls.graph = cls.folder / 'graph.dbg'
        base = cls.folder / 'annotation'
        cls.annotation = base.with_suffix('.row_diff_brwt_coord.annodbg')
        command(MG, 'build', '-p', '2', '--mode', 'basic', '-k', '5', '-o', cls.graph, *fastas)
        command(MG, 'transform', '--mask-dummy', '-p', '2', cls.graph)
        command(MG, 'annotate', '-p', '1', '-i', cls.graph, '-o', base,
                '--anno-type', 'column', '--anno-filename', '--coordinates', *fastas)
        rd = cls.folder / 'rd'
        rd.mkdir()
        for stage in range(3):
            command(MG, 'transform_anno', '-p', '2', '--anno-type', 'row_diff',
                    '--row-diff-stage', stage, '--coordinates', '-i', cls.graph,
                    '-o', rd / 'out', base.with_suffix('.column.annodbg'))
        command(MG, 'transform_anno', '-p', '2', '--anno-type', 'row_diff_brwt_coord', '--greedy',
                '-i', cls.graph, '-o', base, rd / 'annotation.column.annodbg')
        columns = command(MG, 'stats', '--print-col-names', '-a', cls.annotation).stdout.splitlines()
        start = next(i for i, line in enumerate(columns) if line.startswith('Number of columns: '))
        count = int(columns[start].split(': ')[1])
        ordered = columns[start + 1:start + 1 + count]
        if sorted(ordered) != sorted(map(str, fastas)):
            raise AssertionError('Fixture annotation columns do not match its FASTA files')
        command(MG, 'annotate', '-p', '2', '-i', cls.graph, '-o', base,
                '--anno-filename', '--index-header-coords', *ordered)
        cls.inventory = cls.folder / 'headers.tsv'
        cls.inventory.write_text(command(MG, 'stats', '--print-sequence-headers',
                                         base.with_suffix('.seqs')).stdout)
        catalogue = cls.folder / 'release.catalog.gz'
        with gzip.open(catalogue, 'wt') as stream:
            for accession, taxid in [('NC_000001.1', 100), ('NC_000002.1', 200),
                                     ('NC_000003.1', 100), ('NC_000004.1', 300)]:
                stream.write(f'{taxid}\ttest organism\t{accession}\tbacteria\tVALIDATED\t15\n')
        prepared = cls.folder / 'prepared'
        command(sys.executable, PREPARE, '--catalogue', catalogue,
                '--chunk', f'test={cls.inventory}', '--output', prepared)
        cls.coverage = json.loads((prepared / 'coverage.json').read_text())
        cls.sidecar = base.with_suffix('.taxids.tsv')
        shutil.copyfile(prepared / 'test' / 'annotation.taxids.tsv', cls.sidecar)

    @unittest.skipUnless(os.environ.get('TAXONOMY_PILOT_BINARY'), 'set TAXONOMY_PILOT_BINARY')
    def test_membership_pilot_matches_coordinate_reference(self):
        for mode in ('sample', 'contiguous', 'full'):
            output = self.folder / ('pilot-' + mode)
            command(os.environ['TAXONOMY_PILOT_BINARY'], self.graph, self.annotation,
                    self.folder / 'annotation.seqs', self.sidecar, output, mode, 120, 16 * 2**20)
            report = json.loads((output / 'result.json').read_text())
            self.assertTrue(report['ok'])
            self.assertGreater(report['reference_rows_checked'], 0)
            self.assertGreater(report['membership_brwt_bytes'], 0)
            self.assertGreater(report['native_nonanchor_rows'], 0)
            self.assertEqual(report['rows_planned'], report['native_rows_validated'])
            self.assertGreater(report['max_source_path_rows'], 1)
            self.assertTrue(report['query_measurements'])
            for timing in report['query_measurements']:
                self.assertEqual(timing['ondemand_hash'], timing['precomputed_hash'])
            native = output / 'taxid-coordinates.row_diff_brwt_coord.annodbg'
            command(MG, 'stats', '--check-names', '-a', native)
            request = {'patterns': [{'dna': 'AAAAACCCCCGGGGG'}], 'mode': 'partial',
                       'strands': 'forward', 'predicate': {'any': ['100']},
                       'predicate_strands': 'context', 'long_search': 'supported_paths',
                       'output': {'labels': 'all', 'occurrences': True}}
            with server(output / 'sample.dbg', native) as url:
                status, answer = http(url + '/pattern', request)
                self.assertEqual(200, status, answer)
                self.assertEqual(0, answer['patterns'][0]['counts']['selected']['value'])
                request['predicate'] = {'any': ['300']}
                status, answer = http(url + '/pattern', request)
                self.assertEqual(200, status, answer)
                self.assertEqual(1, answer['patterns'][0]['counts']['selected']['value'])
                self.assertIn('NC_000004.1', json.dumps(answer))
            refused = command(os.environ['TAXONOMY_PILOT_BINARY'], self.graph, self.annotation,
                              self.folder / 'annotation.seqs', self.sidecar, output, mode,
                              120, 16 * 2**20, ok=False)
            self.assertNotEqual(0, refused.returncode)
            self.assertEqual(report, json.loads((output / 'result.json').read_text()))

    def test_header_preparation_and_http_semantics(self):
        self.assertTrue(self.coverage['complete'])
        self.assertEqual(4, self.coverage['chunks']['test']['mapped_records'])
        self.assertEqual(0, self.coverage['missing_accessions'])
        request = {'patterns': [{'dna': 'AAAAA'}], 'mode': 'partial', 'strands': 'forward',
                   'predicate_namespace': 'taxid', 'predicate': {'any': ['100']},
                   'predicate_strands': 'context',
                   'output': {'labels': 'all', 'occurrences': True}}
        with server(self.graph, self.annotation, '--record-taxonomy') as url:
            status, capabilities = http(url + '/pattern/capabilities')
            self.assertEqual(200, status)
            self.assertEqual(['column', 'taxid'], capabilities['predicate']['namespaces'])
            CONTRACT.Checker(CONTRACT.Raising(), 'capabilities').block(capabilities, 'capabilities')
            status, answer = http(url + '/pattern', request)
            self.assertEqual(200, status, answer)
            CONTRACT.Checker(CONTRACT.Raising(), 'taxid').answer(answer, request, capabilities)
            self.assertEqual(1, answer['patterns'][0]['counts']['selected']['value'])
            self.assertIn('NC_000001.1', json.dumps(answer['patterns'][0]['results']))
            request['predicate'] = {'any': ['200']}
            status, answer = http(url + '/pattern', request)
            self.assertEqual(200, status, answer)
            self.assertEqual(0, answer['patterns'][0]['counts']['selected']['value'])
            request.update(patterns=[{'dna': 'AAAAACCCCCGGGGG'}], long_search='supported_paths',
                           predicate={'any': ['100']})
            status, answer = http(url + '/pattern', request)
            self.assertEqual(200, status, answer)
            self.assertEqual(0, answer['patterns'][0]['counts']['selected']['value'])
            request['supported_paths_level'] = 'label_intersection'
            status, answer = http(url + '/pattern', request)
            self.assertEqual(200, status, answer)
            self.assertEqual(1, answer['patterns'][0]['counts']['selected']['value'])
            CONTRACT.Checker(CONTRACT.Raising(), 'long-taxid').answer(answer, request, capabilities)
        with server(self.graph, self.annotation) as url:
            _, capabilities = http(url + '/pattern/capabilities')
            self.assertEqual(['column'], capabilities['predicate']['namespaces'])
            status, answer = http(url + '/pattern', request)
            self.assertEqual(400, status, answer)
            self.assertIn('taxonomy_unavailable', json.dumps(answer))

    def test_manifest_inventory_identity_and_load_profile(self):
        fingerprints = []
        for enabled in (False, True):
            flags = ['--record-taxonomy'] if enabled else []
            native = json.loads(command(MG, 'traverse', '--index-inventory', '-i', self.graph,
                                        '-a', self.annotation, *flags).stdout)
            python = json.loads(command(sys.executable, MANIFEST, '--inventory', '-i', self.graph,
                                        '-a', self.annotation, *flags).stdout)
            self.assertEqual(native['files'], python['files'])
            roles = [entry['role'] for entry in native['files']]
            self.assertEqual(enabled, 'record_taxonomy' in roles)
            manifest = self.folder / f'manifest-{enabled}.json'
            command(sys.executable, MANIFEST, '-i', self.graph, '-a', self.annotation,
                    '-o', manifest, '--metagraph', MG, *flags)
            command(sys.executable, MANIFEST, '--verify', manifest, '--metagraph', MG, *flags)
            fingerprints.append(json.loads(manifest.read_text())['index_fp'])
            with server(self.graph, self.annotation, '--index-manifest', manifest, *flags):
                pass
            wrong_profile = [] if enabled else ['--record-taxonomy']
            refused = command(sys.executable, MANIFEST, '--verify', manifest, *wrong_profile, ok=False)
            self.assertNotEqual(0, refused.returncode)
        self.assertNotEqual(*fingerprints)

    def test_missing_and_misaligned_mapping_are_refused_only_when_enabled(self):
        with tempfile.TemporaryDirectory() as folder:
            base = Path(folder) / 'annotation'
            annotation = base.with_suffix('.row_diff_brwt_coord.annodbg')
            annotation.symlink_to(self.annotation)
            base.with_suffix('.seqs').symlink_to(self.folder / 'annotation.seqs')
            for data in (None, 'METAGRAPH_RECORD_TAXIDS_V1\n0\t0\tWRONG.1\t100\n'):
                if data is not None:
                    base.with_suffix('.taxids.tsv').write_text(data)
                command(MG, 'stats', '--check-names', '-a', annotation)
                refused = command(MG, 'stats', '--check-names', '-a', annotation,
                                  '--record-taxonomy', ok=False)
                self.assertNotEqual(0, refused.returncode)


if __name__ == '__main__':
    unittest.main()
