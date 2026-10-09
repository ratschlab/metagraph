"""A multi-graph server as production runs them (`server_query GRAPHS.csv --mmap --mem-cap-gb`):
three chunks of the mini index's records (scripts/traversal/build_mini_refseq.sh's FASTA), one
masked with a coordinate annotation and its record mapping like refseq33m, two unmasked with
column annotations, each listed under a name of its own and two of them also under one name
together (a name spanning two graphs). What such a server serves, as /search does:

  - POST /pattern takes /search's `graphs`; each selected pair is answered as a single-graph
    server answers it, tagged with its pair, the answers concatenated;
  - POST /traverse and /resolve take `graphs: [name]` beside `graph`; with it a seed the graph
    does not hold is a per-seed result (outcome.walks not_in_graph, its k-mers present);
  - `in_ram` on /pattern, /traverse and /resolve loads the pair into RAM for the request (a
    server on mmap, the pair within --mem-cap-gb), the budgets after the load, timing.load_ms;
  - GET /capabilities: the per-pair summary (graph_summary) and columns_disjoint;
  - GET /pattern/capabilities?graph= and /traverse/capabilities?graph=: the pair's block.

The answers are checked against oracles of their own: the CLI (`metagraph pattern`) on the same
files, the k-mers of the chunk's records (the contexts of a pattern, a seed's k-mers present),
and the server's answers without the new fields.
"""

import json
import os
import shlex
import socket
import subprocess
import sys
import tempfile
import time
import unittest

import requests

from base import PROTEIN_MODE, TestingBase, METAGRAPH

MINI_DIR = os.environ.get('METAGRAPH_MINI_REFSEQ', os.path.join(os.getcwd(), 'mini_refseq'))
FASTA = os.path.join(MINI_DIR, 'fasta')
K = 31
# the chunks: (name, the mini's FASTA files, masked, annotation)
CHUNKS = [
    ('ecoli', ['562.fa'], True, 'row_diff_brwt_coord'),
    ('kleb', ['573.fa'], False, 'column'),
    ('pseudo', ['287.fa', '470.fa'], False, 'column'),
]
# a name over two graphs, as a production list names an index's chunks
SPLIT = ('enterics', ['ecoli', 'kleb'])


def free_port():
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(('127.0.0.1', 0))
        return s.getsockname()[1]


def records_of(paths):
    out = []
    for path in paths:
        seq = []
        with open(path) as f:
            for line in f:
                if line.startswith('>'):
                    if seq:
                        out.append(''.join(seq))
                    seq = []
                else:
                    seq.append(line.strip().upper())
        if seq:
            out.append(''.join(seq))
    return out


def revcomp(s):
    return s[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def pair_names(capabilities):
    return [p['graph'] for p in capabilities['graph_summary']['pairs']]


def untimed(value):
    if isinstance(value, dict):
        return {k: untimed(v) for k, v in value.items() if k != 'timing'}
    if isinstance(value, list):
        return [untimed(v) for v in value]
    return value


class Server:
    """A server_query process on a free port, ready once GET /stats answers 200."""

    def __init__(self, args, log_path, timeout=300):
        os.environ['NO_PROXY'] = '127.0.0.1'
        for _ in range(20):
            self.port = free_port()
            self.log_path = log_path
            self.log = open(log_path, 'ab')
            self.process = subprocess.Popen(
                shlex.split(METAGRAPH) + ['server_query'] + args
                + ['--port', str(self.port), '--address', '127.0.0.1', '-p', '4'],
                stdout=self.log, stderr=subprocess.STDOUT)
            if self._wait_ready(timeout):
                return
            self.stop()
        raise AssertionError(f"Couldn't start the server: {args} (log: {log_path})")

    def _wait_ready(self, timeout):
        deadline = time.time() + timeout
        while time.time() < deadline:
            if self.process.poll() is not None:
                return False
            try:
                if requests.get(self.url('stats'), timeout=10).status_code == 200:
                    return True
            except requests.exceptions.ConnectionError:
                pass
            time.sleep(0.2)
        return False

    def url(self, route):
        return f'http://127.0.0.1:{self.port}/{route}'

    def post(self, route, payload):
        return requests.post(self.url(route), data=json.dumps(payload), timeout=300)

    def get(self, route):
        return requests.get(self.url(route), timeout=60)

    def text(self):
        with open(self.log_path) as f:
            return f.read()

    def stop(self):
        if self.process.poll() is None:
            self.process.kill()
            self.process.wait()
        self.log.close()


@unittest.skipIf(PROTEIN_MODE, "the mini index is DNA")
@unittest.skipUnless(os.path.isdir(FASTA),
                     "the mini index is not built (scripts/traversal/build_mini_refseq.sh)")
class TestMultiGraphServer(TestingBase):

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        d = cls.tempdir.name
        cls.pairs = {}
        cls.records = {}
        cls.columns = {}
        for name, files, masked, anno in CHUNKS:
            out = os.path.join(d, name)
            os.makedirs(out)
            fastas = [os.path.join(FASTA, f) for f in files]
            graph = os.path.join(out, 'graph.dbg')
            cls._run_command(f'{METAGRAPH} build -p 4 --mode basic --graph succinct --state stat '
                             f'-k {K} {"--mask-dummy --in-ram " if masked else ""}'
                             f'-o {out}/graph {" ".join(fastas)}', 'build')
            cls._annotate_graph(fastas, graph, os.path.join(out, 'anno'), anno,
                                anno_type='filename')
            if anno == 'row_diff_brwt_coord':
                # the record mapping (.seqs), as build_mini_refseq.sh writes it
                cls._run_command(f'{METAGRAPH} annotate -p 4 -i {graph} --anno-filename '
                                 f'--index-header-coords -o {out}/anno {" ".join(fastas)}',
                                 'annotate --index-header-coords')
            cls.pairs[name] = (graph, os.path.join(out, f'anno.{anno}.annodbg'))
            assert all(os.path.exists(p) for p in cls.pairs[name]), cls.pairs[name]
            # the longest record first (the tests cut their seeds and patterns from it)
            cls.records[name] = sorted(records_of(fastas), key=len, reverse=True)
            cls.columns[name] = fastas
        cls.csv = os.path.join(d, 'graphs.csv')
        with open(cls.csv, 'w') as f:
            for name, _, _, _ in CHUNKS:
                f.write(f'{name},{cls.pairs[name][0]},{cls.pairs[name][1]}\n')
            for chunk in SPLIT[1]:
                f.write(f'{SPLIT[0]},{cls.pairs[chunk][0]},{cls.pairs[chunk][1]}\n')
        cls.server = Server([cls.csv, '--mmap', '--mem-cap-gb', '1'],
                            os.path.join(d, 'server.log'))

    @classmethod
    def tearDownClass(cls):
        cls.server.stop()
        super().tearDownClass()

    # -------------------------------------------------------------------- helpers

    def kmers(self, chunk):
        out = set()
        for r in self.records[chunk]:
            for i in range(len(r) - K + 1):
                out.add(r[i:i + K])
        return out

    def contexts(self, chunk, pattern):
        """The oracle of a count on a masked BASIC graph: every (strand, k-mer, offset) whose
        k-mer holds the pattern there, on + as written and on - as its reverse complement."""
        out = []
        for kmer in self.kmers(chunk):
            for strand, p in (('+', pattern), ('-', revcomp(pattern))):
                for o in range(K - len(p) + 1):
                    if kmer[o:o + len(p)] == p:
                        out.append((strand, kmer, o))
        return out

    def cli_pattern(self, chunk, request):
        """`metagraph pattern` on the chunk's files: the answer a single-graph server gives"""
        path = os.path.join(self.tempdir.name, f'request_{chunk}.json')
        with open(path, 'w') as f:
            json.dump(request, f)
        graph, anno = self.pairs[chunk]
        res = subprocess.run(shlex.split(METAGRAPH) + ['pattern', '--json', '-i', graph, '-a',
                                                       anno, path],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(0, res.returncode, res.stderr.decode()[-2000:])
        return json.loads(res.stdout.decode().strip().split('\n')[-1])

    def post(self, route, payload, status=200):
        ret = self.server.post(route, payload)
        self.assertEqual(status, ret.status_code, ret.text[:2000])
        return ret.json()

    # -------------------------------------------------------------------- capabilities

    def test_capabilities_summarise_every_pair(self):
        caps = self.server.get('capabilities').json()
        self.assertEqual('multi', caps['mode'])
        self.assertEqual(['ecoli', 'enterics', 'kleb', 'pseudo'], caps['graphs'])
        self.assertIn('pattern', caps['features'])
        self.assertEqual('POST /pattern', caps['routes']['pattern'])
        self.assertEqual('GET /pattern/capabilities?graph={name}[&graph_path={path}]',
                         caps['routes']['pattern_capabilities'])
        self.assertEqual({'routes': ['pattern', 'resolve', 'search', 'traverse'], 'loads': True,
                          'mem_cap_gb': 1, 'budgets_start': 'after_load'}, caps['in_ram'])
        # the server-wide block: the contract and the caps, no graph
        self.assertTrue(caps['pattern']['available'])
        self.assertIsNone(caps['pattern']['k'])
        summary = caps['graph_summary']
        # the chunks' columns are their FASTA files: no file in two chunks
        self.assertEqual((True, 0), (summary['columns_disjoint'], summary['shared_columns']))
        listed = [(p['graph'], p['graph_path']) for p in summary['pairs']]
        self.assertEqual([('ecoli', self.pairs['ecoli'][0]), ('enterics', self.pairs['ecoli'][0]),
                          ('enterics', self.pairs['kleb'][0]), ('kleb', self.pairs['kleb'][0]),
                          ('pseudo', self.pairs['pseudo'][0])], listed)
        for p in summary['pairs']:
            chunk = {v[0]: c for c, v in self.pairs.items()}[p['graph_path']]
            masked = chunk == 'ecoli'
            with self.subTest(pair=(p['graph'], chunk)):
                self.assertEqual(self.pairs[chunk][1], p['annotation_path'])
                self.assertEqual((K, 'basic', True, None), (p['k'], p['graph_mode'],
                                                            p['available'],
                                                            p['unavailable_reason']))
                self.assertEqual(('file', 'exact') if masked else ('absent', 'upper_bound'),
                                 (p['mask'], p['counting']))
                self.assertEqual((None, None), (p['index_ns'], p['index_fp']))
                self.assertEqual({'regime': 'basic', 'num_labels': len(self.columns[chunk]),
                                  'has_coordinates': masked, 'has_coord_to_header': masked,
                                  'supports_trace': masked}, p['traversal'])
                # what the pair's own documents say
                block = self.server.get(f'pattern/capabilities?graph={chunk}').json()
                for key in ('k', 'graph_mode', 'available', 'mask', 'counting'):
                    self.assertEqual(block[key], p[key], key)
                probe = self.server.get(f'traverse/capabilities?graph={chunk}').json()
                for key in ('regime', 'num_labels', 'has_coordinates', 'has_coord_to_header',
                            'supports_trace'):
                    self.assertEqual(probe[key], p['traversal'][key], key)

    def test_columns_shared_by_two_pairs(self):
        """A column of two pairs (a FASTA annotated on two graphs): columns_disjoint false, and
        the shared name counted once."""
        d = self.tempdir.name
        out = os.path.join(d, 'overlap')
        os.makedirs(out, exist_ok=True)
        ecoli = os.path.join(FASTA, '562.fa')
        self._run_command(f'{METAGRAPH} annotate -p 4 -i {self.pairs["kleb"][0]} --anno-filename '
                          f'-o {out}/anno {ecoli}', 'annotate')
        csv = os.path.join(d, 'overlap.csv')
        with open(csv, 'w') as f:
            f.write(f'ecoli,{self.pairs["ecoli"][0]},{self.pairs["ecoli"][1]}\n'
                    f'other,{self.pairs["kleb"][0]},{out}/anno.column.annodbg\n'
                    f'kleb,{self.pairs["kleb"][0]},{self.pairs["kleb"][1]}\n')
        server = Server([csv], os.path.join(d, 'server_overlap.log'))
        try:
            summary = server.get('capabilities').json()['graph_summary']
            self.assertEqual((False, 1), (summary['columns_disjoint'], summary['shared_columns']))
            self.assertIn('columns of more than one', server.text())
        finally:
            server.stop()

    def test_a_pair_the_pattern_search_does_not_serve(self):
        """A graph the engine does not recognise (a hash graph) beside a served one: its summary
        entry and its block say why, a request naming it is refused with its reason (the whole
        request, as /search fails on one graph), and the server-wide block is available while
        one pair is served -- with the first pair's reason when none is."""
        d = self.tempdir.name
        out = os.path.join(d, 'hash')
        os.makedirs(out, exist_ok=True)
        record = os.path.join(FASTA, '1296536.fa')
        self._run_command(f'{METAGRAPH} build -p 4 --graph hash -k {K} -o {out}/graph {record}',
                          'build hash')
        self._run_command(f'{METAGRAPH} annotate -p 4 -i {out}/graph.orhashdbg --anno-filename '
                          f'-o {out}/anno {record}', 'annotate hash')
        hashed = f'hashed,{out}/graph.orhashdbg,{out}/anno.column.annodbg\n'
        for name, lines, available in (
                ('mixed', hashed + f'ecoli,{self.pairs["ecoli"][0]},{self.pairs["ecoli"][1]}\n',
                 True),
                ('only', hashed, False)):
            csv = os.path.join(d, f'hash_{name}.csv')
            with open(csv, 'w') as f:
                f.write(lines)
            server = Server([csv], os.path.join(d, f'server_hash_{name}.log'))
            try:
                caps = server.get('capabilities').json()
                pair = {p['graph']: p for p in caps['graph_summary']['pairs']}['hashed']
                self.assertEqual((False, 'representation_unsupported', None, None),
                                 (pair['available'], pair['unavailable_reason'], pair['mask'],
                                  pair['counting']))
                self.assertEqual(available, caps['pattern']['available'])
                self.assertEqual(None if available else 'representation_unsupported',
                                 caps['pattern']['unavailable_reason'])
                block = server.get('pattern/capabilities?graph=hashed').json()
                self.assertEqual((False, 'representation_unsupported'),
                                 (block['available'], block['unavailable_reason']))
                p = {'patterns': [{'dna': 'ACGTTGCAACGTTGCAAG'}], 'mode': 'count'}
                ret = server.post('pattern', dict(p, graphs=list(pair_names(caps))))
                self.assertEqual(400, ret.status_code, ret.text)
                self.assertEqual('representation_unsupported', ret.json()['code'])
                if available:
                    ret = server.post('pattern', dict(p, graphs=['ecoli']))
                    self.assertEqual(200, ret.status_code, ret.text)
            finally:
                server.stop()

    def test_pattern_capabilities_name_a_pair(self):
        enterics = self.server.get('pattern/capabilities?graph=enterics')
        self.assertEqual(400, enterics.status_code)
        self.assertIn('spans several graphs', enterics.json()['error'])
        kleb = self.server.get('pattern/capabilities?graph=kleb').json()
        quoted = requests.utils.quote(self.pairs['kleb'][0], safe='')
        picked = self.server.get(f'pattern/capabilities?graph=enterics&graph_path={quoted}').json()
        self.assertEqual(('kleb', self.pairs['kleb'][0]), (kleb.pop('graph'), kleb.pop('graph_path')))
        self.assertEqual(('enterics', self.pairs['kleb'][0]),
                         (picked.pop('graph'), picked.pop('graph_path')))
        self.assertEqual(kleb, picked)
        probe = self.server.get('traverse/capabilities?graph=kleb').json()
        self.assertEqual(dict(kleb, details='GET /pattern/capabilities'), probe['pattern'])
        for query, words in (('', 'needs ?graph=<name>'), ('?graph=kleb&x=1', "unknown parameter 'x'"),
                             ('?graph=nope', "unknown graph 'nope'")):
            ret = self.server.get('pattern/capabilities' + query)
            self.assertEqual(400, ret.status_code, query)
            self.assertIn(words, ret.json()['error'], query)

    def test_unmasked_pairs_are_sampled_at_load(self):
        """The dummy fraction of every graph without its mask is sampled at start-up, in its
        loading thread, as a single-graph server samples its graph: the first /pattern on an
        untouched unmasked pair, with a budget just above the finalisation reserve, answers 200
        and samples nothing (the sample's trace line, one per unmasked graph, is in the log
        before the first request and no other follows it)."""
        d = self.tempdir.name
        server = Server([self.csv, '--mmap', '--mem-cap-gb', '1', '-v'],
                        os.path.join(d, 'server_sampled.log'))
        try:
            sampled = 'Dummy fraction sampled for the pattern search'
            unmasked = [c for c, _, masked, _ in CHUNKS if not masked]
            self.assertEqual(len(unmasked), server.text().count(sampled), server.text()[-3000:])
            caps = server.get('pattern/capabilities?graph=kleb').json()
            reserve = caps['finalize_reserve_ms']
            kleb = self.records['kleb'][0]
            ret = server.post('pattern', {'patterns': [{'dna': kleb[3000:3020]}],
                                          'mode': 'count', 'graphs': ['kleb'],
                                          'time_budget_ms': reserve + 10})
            self.assertEqual(200, ret.status_code, ret.text[:2000])
            answer = ret.json()['answers'][0]
            self.assertEqual('upper_bound', answer['index']['counting'])
            self.assertEqual('sampled', answer['index']['dummy_fraction']['source'])
            self.assertEqual(caps['dummy_fraction'], answer['index']['dummy_fraction'])
            self.assertEqual(len(unmasked), server.text().count(sampled))
        finally:
            server.stop()

    # -------------------------------------------------------------------- POST /pattern

    def test_pattern_answers_each_selected_pair(self):
        ecoli = self.records['ecoli'][0]
        pattern = ecoli[123456:123476]
        request = {'patterns': [{'dna': pattern, 'id': 'p'}, {'dna': 'ACGTTGCAACGTTGCAAG'}],
                   'mode': 'count'}
        out = self.post('pattern', dict(request, graphs=['pseudo', 'ecoli', 'ecoli']))
        self.assertEqual(['ecoli', 'pseudo'], out['graphs'])
        self.assertEqual(1, out['pattern_contract_version'])
        self.assertEqual(['ecoli', 'pseudo'], [a['graph'] for a in out['answers']])
        for answer in out['answers']:
            chunk = answer['graph']
            with self.subTest(chunk=chunk):
                self.assertEqual(self.pairs[chunk], (answer['graph_path'],
                                                     answer['annotation_path']))
                self.assertIsNone(answer['index_fp'])
                tagless = {k: v for k, v in answer.items()
                           if k not in ('graph', 'graph_path', 'annotation_path', 'index_fp')}
                # the answer of a single-graph server on the pair's files
                self.assertEqual(untimed(self.cli_pattern(chunk, request)), untimed(tagless))
        # the masked chunk counts exactly: the contexts of its records' k-mers
        counted = out['answers'][0]['patterns'][0]['counts']['contexts']
        expected = self.contexts('ecoli', pattern)
        self.assertGreater(len(expected), 0)
        self.assertEqual(('exact', len(expected)), (counted['relation'], counted['value']))
        # a name over two graphs: one answer per pair, in the list's order
        both = self.post('pattern', dict(request, graphs=['enterics']))
        self.assertEqual([('enterics', self.pairs['ecoli'][0]),
                          ('enterics', self.pairs['kleb'][0])],
                         [(a['graph'], a['graph_path']) for a in both['answers']])
        self.assertEqual(untimed(out['answers'][0]['patterns']),
                         untimed(both['answers'][0]['patterns']))
        # without graphs: every name of this small server
        every = self.post('pattern', request)
        self.assertEqual(['ecoli', 'enterics', 'kleb', 'pseudo'], every['graphs'])
        self.assertEqual(5, len(every['answers']))

    def test_pattern_refusals(self):
        p = {'patterns': [{'dna': 'ACGTTGCAACGTTGCAAG'}], 'mode': 'count'}
        for payload, status, code, words in (
                (dict(p, graphs=['nope']), 400, 'invalid_request', "unknown graph 'nope'"),
                (dict(p, graphs=[]), 400, 'invalid_request', 'non-empty array'),
                (dict(p, graphs='ecoli'), 400, 'invalid_request', 'non-empty array'),
                (dict(p, graphs=['ecoli'], in_ram='yes'), 400, 'invalid_request', 'in_ram'),
                (dict(p, graphs=['ecoli'], budget_split=1), 400, 'later_increment',
                 'budget_split'),
                (dict(p, graphs=['ecoli'], mode='x'), 400, 'invalid_request', 'mode')):
            with self.subTest(payload=payload):
                body = self.post('pattern', payload, status)
                self.assertEqual({'error', 'code'}, set(body))
                self.assertEqual(code, body['code'])
                self.assertIn(words, body['error'])

    # -------------------------------------------------------------------- in_ram

    def test_in_ram_loads_the_pair_for_the_request(self):
        ecoli = self.records['ecoli'][0]
        p = {'patterns': [{'dna': ecoli[5000:5020]}], 'mode': 'count', 'graphs': ['ecoli']}
        resident = self.post('pattern', p)
        loaded = self.post('pattern', dict(p, in_ram=True))
        self.assertEqual(untimed(resident), untimed(loaded))
        self.assertNotIn('load_ms', resident['answers'][0]['timing'])
        self.assertGreater(loaded['answers'][0]['timing']['load_ms'], 0)
        # in_ram false: stated, nothing loaded
        self.assertEqual(0, self.post('pattern', dict(p, in_ram=False))
                         ['answers'][0]['timing']['load_ms'])

        seed = ecoli[20000:20100]
        t = {'graphs': ['ecoli'], 'seeds': [{'sequence': seed}],
             'strategy': {'direction': 'right', 'bounds': {'max_extension_bp': 50},
                          'output': {'detail': 'full'}}}
        resident = self.post('traverse', t)
        loaded = self.post('traverse', dict(t, in_ram=True, attempt_id='in-ram-1'))
        self.assertEqual(untimed(resident['results']), untimed(loaded['results']))
        self.assertGreater(loaded['timing']['load_ms'], 0)
        self.assertNotIn('load_ms', resident['timing'])
        # the attempt's bound starts after the load
        self.assertGreaterEqual(loaded['usage']['bound']['load_ms'], 1)
        r = {'graphs': ['ecoli'], 'sequence': seed, 'discover': {'max_labels': 10}}
        resident = self.post('resolve', r)
        loaded = self.post('resolve', dict(r, in_ram=True))
        self.assertEqual(untimed(resident), untimed(loaded))
        self.assertGreater(loaded['timing']['load_ms'], 0)
        self.assertIn('loaded into RAM for it', self.server.text())
        # /search as it always was
        query = {'FASTA': f'>q\n{seed}', 'graphs': ['ecoli'], 'discovery_fraction': 0.5}
        self.assertEqual(self.post('search', query), self.post('search', dict(query, in_ram=True)))

    # -------------------------------------------------------------------- the fan-out of /traverse

    def test_traverse_fans_out_with_not_in_graph(self):
        ecoli, kleb = self.records['ecoli'][0], self.records['kleb'][0]
        seeds = {'e': ecoli[30000:30080], 'k': kleb[40000:40080],
                 'mixed': ecoli[60000:60040] + kleb[70000:70040]}
        request = {'seeds': [{'seed_id': i, 'sequence': s} for i, s in seeds.items()],
                   'strategy': {'direction': 'both', 'bounds': {'max_extension_bp': 30},
                                'output': {'detail': 'summary', 'timing': False}}}
        seen = set()
        for chunk in ('ecoli', 'kleb', 'pseudo'):
            kmers = self.kmers(chunk)
            out = self.post('traverse', dict(request, graphs=[chunk]))
            for result in out['results']:
                seed = seeds[result['seed']['seed_id']]
                present = sum(seed[i:i + K] in kmers for i in range(len(seed) - K + 1))
                seen.add('walked' if present == len(seed) - K + 1
                         else 'partly' if present else 'absent')
                with self.subTest(chunk=chunk, seed=result['seed']['seed_id']):
                    if present == len(seed) - K + 1:
                        self.assertNotEqual('not_in_graph', result['outcome']['walks'])
                        self.assertNotIn('not_in_graph', result)
                    else:
                        self.assertEqual('not_in_graph', result['outcome']['walks'])
                        self.assertEqual({'kmers': len(seed) - K + 1, 'kmers_present': present},
                                         result['not_in_graph'])
            # with graph (not graphs): the request fails on its first absent seed, as before
            absent = any(sum(s[i:i + K] in kmers for i in range(len(s) - K + 1))
                         < len(s) - K + 1 for s in seeds.values())
            ret = self.server.post('traverse', dict(request, graph=chunk))
            self.assertEqual(400 if absent else 200, ret.status_code, chunk)
            if absent:
                self.assertIn('Seed is not fully present in the graph', ret.json()['error'])
        # every case was met: walked, partly present, absent
        self.assertEqual({'walked', 'partly', 'absent'}, seen)
        # the seed of a chunk walks there as in a request that names the graph with graph
        one = dict(request, seeds=[{'seed_id': 'e', 'sequence': seeds['e']}])
        self.assertEqual(self.post('traverse', dict(one, graphs=['ecoli']))['results'],
                         self.post('traverse', dict(one, graph='ecoli'))['results'])
        # a name over two graphs: graphs selects as graph does (graph_path picks the graph)
        ret = self.server.post('traverse', dict(one, graphs=['enterics']))
        self.assertEqual(400, ret.status_code)
        self.assertIn('spans several graphs', ret.json()['error'])
        picked = self.post('traverse', dict(one, graphs=['enterics'],
                                            graph_path=self.pairs['ecoli'][0]))
        self.assertEqual(self.post('traverse', dict(one, graph='ecoli'))['results'],
                         picked['results'])
        for payload, words in ((dict(one, graphs=['ecoli'], graph='ecoli'), 'not both'),
                               (dict(one, graphs=['ecoli', 'kleb']), 'expected [name]')):
            ret = self.server.post('traverse', payload)
            self.assertEqual(400, ret.status_code)
            self.assertIn(words, ret.json()['error'])

    def test_resolve_takes_graphs(self):
        seed = self.records['kleb'][0][10000:10200]
        r = {'sequence': seed, 'discover': {'max_labels': 10}}
        self.assertEqual(untimed(self.post('resolve', dict(r, graph='kleb'))),
                         untimed(self.post('resolve', dict(r, graphs=['kleb']))))


if __name__ == '__main__':
    unittest.main()
