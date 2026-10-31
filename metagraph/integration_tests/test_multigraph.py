"""A multi-graph server as production runs them (`server_query GRAPHS.csv --mmap --mem-cap-gb`):
three chunks of the mini index's records (scripts/traversal/build_mini_refseq.sh's FASTA), one
masked with a coordinate annotation and its record mapping like refseq33m, two unmasked with
column annotations, each listed under a name of its own and two of them also under one name
together (a name spanning two graphs). What such a server serves, as /search does:

  - POST /pattern takes /search's `graphs`; each selected pair is answered as a single-graph
    server answers it, tagged with its pair, the answers concatenated in an envelope; a pair
    the request is refused on is a refused entry of the envelope (the pair's own body under
    `refusal`), the other pairs answered, while what no pair decides stays the request's 400;
  - POST /resolve answers the envelope too (SPEC-multigraph-envelope.md section 1): a named
    graph one entry, `graphs` and a request without a selection every selected pair; both
    envelopes state request_id, server_instance and the deadline block, and a request is cut at
    its work cutoff (not_reached, unfinished), its request_id answered by GET /requests;
  - POST /traverse takes `graphs: [name]` beside `graph`, two spellings of one request; a seed
    the graph does not hold is a per-seed result (outcome.walks not_in_graph, its k-mers
    present), whichever spelling named the graph;
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
import random
import shlex
import shutil
import socket
import subprocess
import sys
import tempfile
import time
import unittest

import requests

from base import PROTEIN_MODE, TestingBase, METAGRAPH
from test_pattern import gate_block

MINI_DIR = os.environ.get('METAGRAPH_MINI_REFSEQ', os.path.join(os.getcwd(), 'mini_refseq'))
FASTA = os.path.join(MINI_DIR, 'fasta')
# $METAGRAPH_REQUIRE_GUARDS=1: a missing mini index fails the class instead of skipping it
REQUIRE_GUARDS = os.environ.get('METAGRAPH_REQUIRE_GUARDS', '') == '1'
K = 31
# the chunks: (name, the mini's FASTA files, masked, annotation)
CHUNKS = [
    ('ecoli', ['562.fa'], True, 'row_diff_brwt_coord'),
    ('kleb', ['573.fa'], False, 'column'),
    ('pseudo', ['287.fa', '470.fa'], False, 'column'),
]
# a name over two graphs, as a production list names an index's chunks
SPLIT = ('enterics', ['ecoli', 'kleb'])
# what an envelope entry carries beside the pair's own answer (SPEC §24.1)
TAGS = ('graph', 'graph_path', 'annotation_path', 'index_fp', 'outcome')
# the members of every POST envelope (SPEC-multigraph-envelope.md section 1); /pattern adds
# pattern_contract_version
ENVELOPE = {'envelope_version', 'route', 'graphs', 'answered', 'refused', 'not_reached',
            'unfinished', 'request_id', 'server_instance', 'deadline', 'answers', 'timing',
            'work'}
# the deadline block of a request that gives neither not_after_ms nor finish_by_ms
NO_DEADLINE = {'not_after_ms': None, 'finish_by_ms': None, 'effective_deadline_at': None,
               'cutoff_at': None}


def free_port():
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(('127.0.0.1', 0))
        return s.getsockname()[1]


def named_records_of(paths):
    """The (header, sequence) of every record of the FASTA files, in their order"""
    out = []
    for path in paths:
        header, seq = None, []
        with open(path) as f:
            for line in f:
                if line.startswith('>'):
                    if seq:
                        out.append((header, ''.join(seq)))
                    header, seq = line[1:].strip(), []
                else:
                    seq.append(line.strip().upper())
        if seq:
            out.append((header, ''.join(seq)))
    return out


def records_of(paths):
    return [seq for _, seq in named_records_of(paths)]


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
    """A server_query process on a free port, ready once GET /stats answers 200; |threads| is
    its -p (the graph pool the pairs of a request run on, and its HTTP threads)."""

    def __init__(self, args, log_path, timeout=300, threads=4):
        os.environ['NO_PROXY'] = '127.0.0.1'
        for _ in range(20):
            self.port = free_port()
            self.log_path = log_path
            self.log = open(log_path, 'ab')
            self.process = subprocess.Popen(
                shlex.split(METAGRAPH) + ['server_query'] + args
                + ['--port', str(self.port), '--address', '127.0.0.1', '-p', str(threads)],
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


def guard(condition, reason):
    """unittest.skipUnless(condition, reason); with $METAGRAPH_REQUIRE_GUARDS=1 a class whose
    input is missing fails instead (test_pattern.py's guard)."""
    if condition or not REQUIRE_GUARDS:
        return unittest.skipUnless(condition, reason)

    def decorate(cls):
        def fail(klass):
            raise AssertionError(f'$METAGRAPH_REQUIRE_GUARDS=1: {cls.__name__} cannot run: '
                                 f'{reason}')
        cls.setUpClass = classmethod(fail)
        return cls
    return decorate


@unittest.skipIf(PROTEIN_MODE, "the mini index is DNA")
@guard(os.path.isdir(FASTA), "the mini index is not built (scripts/traversal/build_mini_refseq.sh)")
class TestMultiGraphServer(TestingBase):

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        d = cls.tempdir.name
        cls.pairs = {}
        cls.records = {}
        cls.named_records = {}
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
            cls.named_records[name] = sorted(named_records_of(fastas), key=lambda r: len(r[1]),
                                             reverse=True)
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

    def cli_refusal(self, chunk, request):
        """`metagraph pattern` on the chunk's files when it refuses the request: the body a
        single-graph server answers (the CLI writes it and exits 1)"""
        path = os.path.join(self.tempdir.name, f'refusal_{chunk}.json')
        with open(path, 'w') as f:
            json.dump(request, f)
        graph, anno = self.pairs[chunk]
        res = subprocess.run(shlex.split(METAGRAPH) + ['pattern', '--json', '-i', graph, '-a',
                                                       anno, path],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(1, res.returncode, res.stderr.decode()[-2000:])
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
        # the threshold of the selection rule (the default of --max-graphs-without-selection)
        self.assertEqual(10, caps['max_graphs_without_selection'])
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
        entry and its block say why, a request naming it answers 200 with that pair as a
        refused entry of the envelope carrying its reason (alone: the envelope with its one
        entry refused, nothing answered), and the server-wide block is available while one
        pair is served -- with the first pair's reason when none is."""
        d = self.tempdir.name
        out_dir = os.path.join(d, 'hash')
        os.makedirs(out_dir, exist_ok=True)
        record = os.path.join(FASTA, '1296536.fa')
        self._run_command(f'{METAGRAPH} build -p 4 --graph hash -k {K} -o {out_dir}/graph {record}',
                          'build hash')
        self._run_command(f'{METAGRAPH} annotate -p 4 -i {out_dir}/graph.orhashdbg --anno-filename '
                          f'-o {out_dir}/anno {record}', 'annotate hash')
        hashed = f'hashed,{out_dir}/graph.orhashdbg,{out_dir}/anno.column.annodbg\n'
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
                self.assertEqual(200, ret.status_code, ret.text)
                out = ret.json()
                self.assertEqual((1 if available else 0, 1), (out['answered'], out['refused']))
                entry = {a['graph']: a for a in out['answers']}['hashed']
                self.assertEqual('refused', entry['outcome'])
                self.assertEqual({'http_status', 'error', 'code'}, set(entry['refusal']))
                self.assertEqual((400, 'representation_unsupported'),
                                 (entry['refusal']['http_status'], entry['refusal']['code']))
                self.assertEqual((f'{out_dir}/graph.orhashdbg', f'{out_dir}/anno.column.annodbg', None),
                                 (entry['graph_path'], entry['annotation_path'], entry['index_fp']))
                if available:
                    self.assertEqual('answered', {a['graph']: a for a in out['answers']}
                                     ['ecoli']['outcome'])
                    ret = server.post('pattern', dict(p, graphs=['ecoli']))
                    self.assertEqual(200, ret.status_code, ret.text)
                    self.assertEqual((1, 0), (ret.json()['answered'], ret.json()['refused']))
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
        # the pair's gate block (SPEC-pattern-search.md §23.2)
        self.assertEqual(gate_block(kleb), probe['pattern'])
        self.assertNotIn('predicate', probe['pattern'])
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
        self.assertEqual(ENVELOPE | {'pattern_contract_version'}, set(out))
        # the members of SPEC-multigraph-envelope.md section 1 (E4): no request_id, no
        # deadline given, every pair reached and finished
        self.assertEqual((1, '/pattern', [], [], None, NO_DEADLINE),
                         (out['envelope_version'], out['route'], out['not_reached'],
                          out['unfinished'], out['request_id'], out['deadline']))
        self.assertRegex(out['server_instance'], '^[0-9a-f]{16}$')
        self.assertEqual((2, 0), (out['answered'], out['refused']))
        self.assertEqual(['ecoli', 'pseudo'], [a['graph'] for a in out['answers']])
        for answer in out['answers']:
            chunk = answer['graph']
            with self.subTest(chunk=chunk):
                self.assertEqual(self.pairs[chunk], (answer['graph_path'],
                                                     answer['annotation_path']))
                self.assertIsNone(answer['index_fp'])
                self.assertEqual('answered', answer['outcome'])
                tagless = {k: v for k, v in answer.items() if k not in TAGS}
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
                (dict(p, graphs=['ecoli'], mode='x'), 400, 'invalid_request', 'mode'),
                # the request's own fault first, although the pair would be refused too (the
                # column annotation under output.labels "all": annotation_unbudgeted)
                (dict(p, graphs=['kleb'], mode='all_or_count', output={'labels': 'all'},
                      budget_split=1), 400, 'later_increment', 'budget_split'),
                (dict(p, graphs=['kleb'], mode='all_or_count', output={'labels': 'all'},
                      time_budget_ms='soon'), 400, 'invalid_request', 'time_budget_ms')):
            with self.subTest(payload=payload):
                body = self.post('pattern', payload, status)
                self.assertEqual({'error', 'code'}, set(body))
                self.assertEqual(code, body['code'])
                self.assertIn(words, body['error'])

    def test_graphs_required_above_the_threshold(self):
        """A request without `graphs` queries every name of a list of at most
        --max-graphs-without-selection names (10 by default) and is refused on a longer list,
        on /pattern as on /search, the refusal naming the threshold; GET /capabilities states
        it, so a client can check it before an added name refuses its requests. A value below 1
        does not start."""
        d = self.tempdir.name
        kleb = self.records['kleb'][0]
        p = {'patterns': [{'dna': kleb[3000:3020]}], 'mode': 'count'}
        s = {'FASTA': '>q\n' + kleb[3000:3100]}
        # the default server lists four names, within the default threshold: every name
        self.assertEqual(['ecoli', 'enterics', 'kleb', 'pseudo'], self.post('pattern', p)['graphs'])
        ret = self.server.post('search', s)
        self.assertEqual(200, ret.status_code, ret.text[:2000])
        server = Server([self.csv, '--mmap', '--mem-cap-gb', '1',
                         '--max-graphs-without-selection', '1'],
                        os.path.join(d, 'server_threshold.log'))
        try:
            self.assertEqual(1, server.get('capabilities').json()['max_graphs_without_selection'])
            ret = server.post('pattern', p)
            self.assertEqual(400, ret.status_code, ret.text[:2000])
            self.assertEqual({'code': 'invalid_request',
                              'error': 'request.graphs: required on this server, which hosts 4 '
                                       'graph names (more than 1; GET /capabilities lists them)'},
                             ret.json())
            ret = server.post('search', s)
            self.assertEqual(400, ret.status_code, ret.text[:2000])
            self.assertEqual({'error': 'Bad request: requests without names (no "graphs" field) '
                                       'are only supported for small indexes (<=1 names)'},
                             ret.json())
            # /resolve reads the selection as /pattern does (SPEC-multigraph-envelope.md
            # section 6)
            ret = server.post('resolve', {'sequence': kleb[3000:3100]})
            self.assertEqual(400, ret.status_code, ret.text[:2000])
            self.assertEqual({'error': 'request.graphs: required on this server, which hosts 4 '
                                       'graph names (more than 1; GET /capabilities lists them)'},
                             ret.json())
            ret = server.post('resolve', {'sequence': kleb[3000:3100], 'graph': 'kleb',
                                          'discover': {'max_labels': 2}})
            self.assertEqual(200, ret.status_code, ret.text[:2000])
            self.assertEqual(['kleb'], ret.json()['graphs'])
            # with `graphs` the threshold plays no part
            ret = server.post('pattern', dict(p, graphs=['kleb']))
            self.assertEqual(200, ret.status_code, ret.text[:2000])
            self.assertEqual(['kleb'], ret.json()['graphs'])
            ret = server.post('search', dict(s, graphs=['kleb']))
            self.assertEqual(200, ret.status_code, ret.text[:2000])
        finally:
            server.stop()
        for value in ('0', '-1', 'abc', '1.5', ''):
            res = subprocess.run(shlex.split(METAGRAPH) + [
                'server_query', self.csv, '--port', str(free_port()), '--address', '127.0.0.1',
                '--max-graphs-without-selection', value],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=60)
            self.assertNotEqual(0, res.returncode, value)
            self.assertIn("--max-graphs-without-selection must be an integer in [1, 2^53 - 1], "
                          f"got '{value}'", res.stderr.decode(), value)

    def test_pattern_and_search_select_the_same_graphs(self):
        """SPEC-pattern-search.md §24.1: a well-formed `graphs` (a list of known names, a name
        given twice, a name over two graphs, or no field) selects the same graphs on /pattern and
        on /search; /pattern reads the field more strictly. A /search answer does not name its
        graphs: they are the chunks whose labels its results name (a FASTA file, or a record's
        accession where the chunk has a record mapping), from a query cut from each chunk's
        records (the chunks' labels are disjoint)."""
        chunks = [name for name, _, _, _ in CHUNKS]
        cut = {c: self.records[c][0][4000:4100] for c in chunks}
        chunk_of_label = {os.path.basename(f): c for c in chunks for f in self.columns[c]}
        chunk_of_label.update({h.split()[0]: c for c in chunks for h, _ in self.named_records[c]})
        chunk_of_graph = {self.pairs[c][0]: c for c in chunks}
        s = {'FASTA': ''.join(f'>{c}\n{cut[c]}\n' for c in chunks), 'discovery_fraction': 0.5}
        p = {'patterns': [{'dna': cut[c][:20]} for c in chunks], 'mode': 'count'}

        def searched(selection):
            out = self.post('search', dict(s, **selection))
            return {chunk_of_label[os.path.basename(r['sample'])]
                    for query in out for r in query['results']}

        def patterned(selection):
            out = self.post('pattern', dict(p, **selection))
            return {chunk_of_graph[a['graph_path']] for a in out['answers']}

        for selection, expected in (({'graphs': ['ecoli']}, {'ecoli'}),
                                    ({'graphs': ['kleb']}, {'kleb'}),
                                    ({'graphs': ['pseudo', 'ecoli', 'ecoli']}, {'ecoli', 'pseudo'}),
                                    ({'graphs': ['enterics']}, {'ecoli', 'kleb'}),
                                    ({'graphs': ['kleb', 'enterics']}, {'ecoli', 'kleb'}),
                                    ({}, set(chunks))):
            with self.subTest(selection=selection):
                self.assertEqual(expected, patterned(selection))
                self.assertEqual(expected, searched(selection))

        # where /pattern is stricter: an empty list, a value that is not a list and an item that
        # is not a string are its 400; /search answers [] with null (it queries nothing), reads a
        # value that is not a list as no field, and an item as its text
        for graphs in ([], 'ecoli', None, [1], ['ecoli', 1]):
            with self.subTest(graphs=graphs):
                body = self.post('pattern', dict(p, graphs=graphs), 400)
                self.assertEqual({'code': 'invalid_request',
                                  'error': 'request.graphs: expected a non-empty array of graph '
                                           'names (GET /capabilities lists them)'}, body)
        self.assertIsNone(self.post('search', dict(s, graphs=[])))
        self.assertEqual(set(chunks), searched({'graphs': 'ecoli'}))
        self.assertEqual(set(chunks), searched({'graphs': None}))
        self.assertEqual({'error': 'Request with an uninitialized graph 1'},
                         self.post('search', dict(s, graphs=[1]), 400))
        # an unknown name: 400 on both, each in its own words
        self.assertEqual({'error': 'Request with an uninitialized graph nope'},
                         self.post('search', dict(s, graphs=['nope']), 400))
        self.assertEqual({'code': 'invalid_request',
                          'error': "request.graphs: unknown graph 'nope' (GET /capabilities "
                                   "lists the graphs)"},
                         self.post('pattern', dict(p, graphs=['nope']), 400))

    def test_a_refused_pair_is_an_entry_of_the_envelope(self):
        """A pair the request cannot be answered on -- here output.labels "all" on a column
        annotation, annotation_unbudgeted -- is a refused entry of the envelope: its tags,
        outcome "refused", and under `refusal` the body a single-graph server answers for it
        (the CLI's, byte for byte) with its HTTP status; the other pairs answer, the request is
        200. A name over two graphs lists each pair's own outcome. Every pair refused is still
        the envelope, answered 0."""
        ecoli = self.records['ecoli'][0]
        request = {'patterns': [{'dna': ecoli[123456:123476], 'id': 'p'}],
                   'mode': 'all_or_count', 'output': {'labels': 'all'}}
        out = self.post('pattern', dict(request, graphs=['kleb', 'ecoli']))
        self.assertEqual((['ecoli', 'kleb'], 1, 1),
                         (out['graphs'], out['answered'], out['refused']))
        answered, refused = out['answers']
        self.assertEqual(('ecoli', 'answered'), (answered['graph'], answered['outcome']))
        self.assertNotIn('refusal', answered)
        tagless = {k: v for k, v in answered.items() if k not in TAGS}
        self.assertEqual(untimed(self.cli_pattern('ecoli', request)), untimed(tagless))
        self.assertEqual({'graph', 'graph_path', 'annotation_path', 'index_fp', 'outcome',
                          'refusal'}, set(refused))
        self.assertEqual(('kleb', self.pairs['kleb'][0], self.pairs['kleb'][1], None, 'refused'),
                         (refused['graph'], refused['graph_path'], refused['annotation_path'],
                          refused['index_fp'], refused['outcome']))
        body = self.cli_refusal('kleb', request)
        self.assertEqual({'error', 'code'}, set(body))
        self.assertEqual('annotation_unbudgeted', body['code'])
        self.assertEqual(dict(body, http_status=400), refused['refusal'])
        self.assertIn("refused (400 annotation_unbudgeted)", self.server.text())
        # a name over two graphs: each pair's outcome, in the list's order
        both = self.post('pattern', dict(request, graphs=['enterics']))
        self.assertEqual([('enterics', self.pairs['ecoli'][0], 'answered'),
                          ('enterics', self.pairs['kleb'][0], 'refused')],
                         [(a['graph'], a['graph_path'], a['outcome']) for a in both['answers']])
        self.assertEqual((1, 1), (both['answered'], both['refused']))
        self.assertEqual(untimed(answered['patterns']), untimed(both['answers'][0]['patterns']))
        self.assertEqual(refused['refusal'], both['answers'][1]['refusal'])
        # every pair refused
        none = self.post('pattern', dict(request, graphs=['pseudo', 'kleb']))
        self.assertEqual((['kleb', 'pseudo'], 0, 2),
                         (none['graphs'], none['answered'], none['refused']))
        self.assertEqual([('kleb', 'refused', 400, 'annotation_unbudgeted'),
                          ('pseudo', 'refused', 400, 'annotation_unbudgeted')],
                         [(a['graph'], a['outcome'], a['refusal']['http_status'],
                           a['refusal']['code']) for a in none['answers']])

    # -------------------------------------------------------------------- deadlines, the registry

    def _cutoff_csv(self):
        """A list of three names in byte order: a_kleb, a chunk of 0.7 Mbp on which the slow
        pattern below completes in about a tenth of a second (and the query of the slow resolve
        below, absent from it, is resolved fast); b_big, the whole mini index (its coordinate
        annotation), on which that pattern works for over a second, and the query, 2 Mbp of
        one of its records, takes over a second to resolve (its k-mers' rows decoded, an
        uninterruptible step of /resolve); c_tiny, a graph of 2 kbp. On a server with one pool
        thread the pairs run one after the other in that order, so b_big starts a tenth of a
        second after the request."""
        mini = (os.path.join(MINI_DIR, 'graph_k31.dbg'),
                os.path.join(MINI_DIR, 'annotation.relaxed.relabeled.row_diff_brwt_coord.annodbg'))
        if not all(os.path.exists(p) for p in mini):
            if REQUIRE_GUARDS:
                self.fail(f'$METAGRAPH_REQUIRE_GUARDS=1: the mini index is missing: {mini}')
            self.skipTest(f'the mini index is not built: {mini}')
        d = self.tempdir.name
        csv = os.path.join(d, 'cutoff.csv')
        if os.path.exists(csv):
            return csv
        tiny = os.path.join(d, 'tiny')
        os.makedirs(tiny, exist_ok=True)
        fasta = os.path.join(tiny, 'tiny.fa')
        rng = random.Random(7)
        with open(fasta, 'w') as f:
            f.write('>tiny\n' + ''.join(rng.choice('ACGT') for _ in range(2000)) + '\n')
        self._run_command(f'{METAGRAPH} build -p 4 --mode basic --graph succinct --state stat '
                          f'-k {K} -o {tiny}/graph {fasta}', 'build tiny')
        self._annotate_graph([fasta], f'{tiny}/graph.dbg', f'{tiny}/anno', 'column',
                             anno_type='filename')
        with open(csv, 'w') as f:
            f.write(f'a_kleb,{self.pairs["kleb"][0]},{self.pairs["kleb"][1]}\n'
                    f'b_big,{mini[0]},{mini[1]}\n'
                    f'c_tiny,{tiny}/graph.dbg,{tiny}/anno.column.annodbg\n')
        return csv

    # the request whose b_big pair works for over a second (see _cutoff_csv)
    SLOW_PATTERN = {'patterns': [{'dna': 'ACGT'}] * 16, 'mode': 'partial', 'scope': 'suffix',
                    'strands': 'both', 'max_contexts': 1000000, 'output': {'labels': 'none'}}

    def _slow_resolve(self, bp=2000000):
        seq = self.records['pseudo'][0][:bp]
        return {'sequence': seq, 'discover': {'max_labels': 2},
                'bounds': {'max_query_bp': len(seq)}}

    def _cutoff_server(self, name, *args):
        """A server on _cutoff_csv with one pool thread, a finalisation reserve of 500 ms on
        /pattern (250 ms on /resolve) and the delivery estimate made negligible"""
        return Server([self._cutoff_csv(), '--pattern-finalize-ms', '500',
                       '--pattern-delivery-build-mbps', '1e9',
                       '--pattern-delivery-compress-mbps', '1e9',
                       '--pattern-max-contexts', '1000000'] + list(args),
                      os.path.join(self.tempdir.name, f'server_cutoff_{name}.log'),
                      threads=1)

    def _wait_ended(self, server, request_id, timeout=120):
        """GET /requests/{request_id} until it states the request ended"""
        deadline = time.time() + timeout
        while True:
            ret = server.get(f'requests/{request_id}')
            self.assertEqual(200, ret.status_code, ret.text[:2000])
            if ret.json()['state'] == 'ended' or time.time() > deadline:
                return ret.json()
            time.sleep(0.05)

    def post_to(self, server, route, payload, status=200):
        """POST to |server| (self.post posts to the class's server)"""
        ret = server.post(route, payload)
        self.assertEqual(status, ret.status_code, ret.text[:2000])
        return ret.json()

    @staticmethod
    def _ref(entry):
        return {k: entry[k] for k in ('graph', 'graph_path', 'annotation_path')}

    def test_the_work_cutoff_keeps_the_completed_pair(self):
        """SPEC-multigraph-envelope.md section 4 (E4): the request's work cutoff is its
        effective deadline (arrival + time_budget_ms, 1500 ms, without finish_by_ms) less the
        finalisation reserve (500 ms). On one pool thread a_kleb completes; b_big starts after
        it with its own budget, which ends after the cutoff: it is still working there, so it
        is unfinished, with its usage so far, and told to stop; c_tiny never started:
        not_reached. The completed answer is written beside them, by the effective deadline,
        as a single-graph server answers it (each pair keeps its own budget)."""
        server = self._cutoff_server('cutoff')
        try:
            started = time.time()
            ret = server.post('pattern', dict(self.SLOW_PATTERN, time_budget_ms=1500))
            took = time.time() - started
            self.assertEqual(200, ret.status_code, ret.text[:2000])
            out = ret.json()
            self.assertEqual((['a_kleb', 'b_big', 'c_tiny'], 1, 2),
                             (out['graphs'], out['answered'], out['refused']),
                             [(a['graph'], a['outcome'], a.get('refusal', {}).get('code'))
                              for a in out['answers']])
            kleb, big, tiny = out['answers']
            self.assertEqual(('a_kleb', 'answered'), (kleb['graph'], kleb['outcome']))
            self.assertEqual(16, len(kleb['patterns']))
            self.assertEqual((1500, []), (kleb['limits']['time_budget_ms'],
                                          kleb['limits']['clamped']))
            self.assertEqual(('b_big', 'refused'), (big['graph'], big['outcome']))
            self.assertEqual({'http_status': 503, 'code': 'unfinished', 'reason': 'cutoff'},
                             {k: big['refusal'][k] for k in ('http_status', 'code', 'reason')})
            self.assertIn("had not finished at the request's work cutoff",
                          big['refusal']['error'])
            self.assertIn('the request has no request_id', big['refusal']['error'])
            usage = big['refusal']['usage']
            self.assertEqual({'state', 'started_at', 'elapsed_ms', 'cpu_ms', 'memory_bytes'},
                             set(usage))
            self.assertEqual('stopping', usage['state'])
            self.assertGreater(usage['elapsed_ms'], 100)
            self.assertEqual(('c_tiny', 'refused'), (tiny['graph'], tiny['outcome']))
            self.assertEqual({'http_status': 503, 'code': 'not_reached', 'reason': 'cutoff',
                              'error': tiny['refusal']['error']}, tiny['refusal'])
            self.assertIn('waited for a pool thread', tiny['refusal']['error'])
            # the top-level lists name the pairs
            self.assertEqual([self._ref(big)], out['unfinished'])
            self.assertEqual([self._ref(tiny)], out['not_reached'])
            # no deadline field given: the block states none (the cutoff applies all the same)
            self.assertEqual(NO_DEADLINE, out['deadline'])
            self.assertLess(took, 30)
            # work.memory_bytes: the answered pair's account
            self.assertEqual(kleb['work']['memory_bytes'], out['work']['memory_bytes'])
        finally:
            server.stop()

    def test_finish_by_ms_cuts_a_pair_and_the_registry_states_its_end(self):
        """finish_by_ms (E20) on /resolve without a budget: the effective deadline is
        finish_by_ms less the clock skew allowance (--traverse-clock-skew-ms 500), the cutoff
        250 ms (the reserve) before it; the deadline block states all four, and a pair's
        deadline is the minimum of its own (the server's maximum) and finish_by_ms - the skew.
        On one pool thread b_big resolves the query (over a second of work) under that
        deadline: it is told to stop at the cutoff, so it is cut there (unfinished) or, when
        its own stop by time wins the race to the cutoff, answered with what it found, its
        clamped budget stated; c_tiny, queued behind it, is cut at the cutoff. (A pair cut
        for certain: test_the_work_cutoff_keeps_the_completed_pair, under a budget.) The
        request's request_id (E4): GET /requests/{request_id} states it ended
        once the cut worker ended, each pair with its state, outcome and final usage, the
        response written (status, bytes); a repeat of the id starts no work (409, the entry
        with an error); an id never sent, or past the retention count (--request-retention
        1), is 404."""
        server = self._cutoff_server('finish', '--traverse-clock-skew-ms', '500',
                                     '--request-retention', '1')
        try:
            # effective = finish_by - 500 (skew), the cutoff 250 ms (reserve) before it
            finish_by = int(time.time() * 1000) + 1300
            body = dict(self._slow_resolve(), graphs=['b_big', 'c_tiny'], finish_by_ms=finish_by,
                        request_id='it-finish-1')
            started = time.time()
            ret = server.post('resolve', body)
            took = time.time() - started
            self.assertEqual(200, ret.status_code, ret.text[:2000])
            out = ret.json()
            # by the effective deadline (the writer refuses a later envelope, 503), not after
            # b_big's whole work
            self.assertLess(took, 3)
            outcomes = [(a['graph'], a['outcome'], a.get('refusal', {}).get('code'))
                        for a in out['answers']]
            self.assertEqual(('/resolve', 'it-finish-1', 2), (out['route'], out['request_id'],
                                                              out['answered'] + out['refused']))
            big, tiny = out['answers']
            # c_tiny is cut: never started, or started a moment before the cutoff
            self.assertEqual(('c_tiny', 'refused', 'cutoff'),
                             (tiny['graph'], tiny['outcome'], tiny['refusal']['reason']),
                             outcomes)
            self.assertIn(tiny['refusal']['code'], ('not_reached', 'unfinished'))
            if big['outcome'] == 'refused':
                self.assertEqual(('unfinished', 'cutoff'), (big['refusal']['code'],
                                                            big['refusal']['reason']), outcomes)
                self.assertIn('GET /requests/it-finish-1', big['refusal']['error'])
            else:
                # its deadline: finish_by_ms - the skew, stated as a clamp of its budget (none
                # requested: the server's maximum)
                clamped = {c['field']: c for c in big['limits']['clamped']}
                self.assertIsNone(clamped['bounds.time_budget_ms']['requested'])
                self.assertLess(clamped['bounds.time_budget_ms']['effective'], 800)
            for field in ('not_reached', 'unfinished'):
                self.assertEqual([self._ref(a) for a in out['answers']
                                  if a.get('refusal', {}).get('code') == field], out[field])
            iso = lambda ms: time.strftime('%Y-%m-%dT%H:%M:%S', time.gmtime(ms // 1000)) \
                + f'.{ms % 1000:03d}Z'
            self.assertEqual({'not_after_ms': None, 'finish_by_ms': finish_by,
                              'effective_deadline_at': iso(finish_by - 500),
                              'cutoff_at': iso(finish_by - 750)}, out['deadline'])
            # the registry: the request ends when the cut worker ends
            entry = self._wait_ended(server, 'it-finish-1')
            self.assertEqual(('it-finish-1', '/resolve', 'ended', 'answered'),
                             (entry['request_id'], entry['route'], entry['state'],
                              entry['reason']))
            self.assertEqual(out['server_instance'], entry['server_instance'])
            # the bytes written: the body as sent (compressed, the client accepting gzip)
            self.assertEqual({'status': 200, 'bytes': int(ret.headers['Content-Length']),
                              'written': True}, entry['response'])
            self.assertIsNotNone(entry['retained_until'])
            p_big, p_tiny = entry['pairs']
            self.assertEqual(('b_big', 'ended', big.get('refusal', {}).get('code', 'answered')),
                             (p_big['graph'], p_big['state'], p_big['outcome']))
            self.assertEqual({'started_at', 'finished_at', 'elapsed_ms', 'cpu_ms',
                              'memory_bytes'}, set(p_big['usage']))
            self.assertIsNotNone(p_big['usage']['finished_at'])
            if big['outcome'] == 'refused':
                # the cut worker's final usage: at least what the envelope stated
                self.assertGreaterEqual(p_big['usage']['elapsed_ms'],
                                        big['refusal']['usage']['elapsed_ms'])
            self.assertEqual(('c_tiny', tiny['refusal']['code']),
                             (p_tiny['graph'], p_tiny['outcome']))
            if tiny['refusal']['code'] == 'not_reached':
                self.assertEqual(('never_started', None), (p_tiny['state'], p_tiny['usage']))
            else:
                self.assertEqual('ended', p_tiny['state'])
            # a repeat starts no work: 409, the entry and an error
            ret = server.post('resolve', body)
            self.assertEqual(409, ret.status_code, ret.text[:2000])
            repeat = ret.json()
            self.assertIn("request_id 'it-finish-1' is running, or retained", repeat['error'])
            self.assertEqual(entry, {k: v for k, v in repeat.items() if k != 'error'})
            self.assertEqual(1, server.text().count('graph name(s), request_id it-finish-1'))
            # an id never sent
            ret = server.get('requests/it-never-sent')
            self.assertEqual(404, ret.status_code, ret.text[:2000])
            self.assertEqual({'error', 'request_id', 'server_instance'}, set(ret.json()))
            self.assertEqual(('it-never-sent', out['server_instance']),
                             (ret.json()['request_id'], ret.json()['server_instance']))
            self.assertIn('not known on this server', ret.json()['error'])
            # retention by count (1): the next request with an id drops this one. One pair
            # (inline) under a finish_by_ms 10 s ahead: its deadline is finish_by_ms - the skew,
            # below the server's maximum (--traverse-max-time-ms, 30 s), which a budget-less
            # request gets under finish_by_ms, stated as a clamp of its budget
            out = self.post_to(server, 'resolve', dict(
                    body, graphs=['c_tiny'], request_id='it-finish-2',
                    finish_by_ms=int(time.time() * 1000) + 10000))
            self.assertEqual((1, 0), (out['answered'], out['refused']))
            clamped = {c['field']: c for c in out['answers'][0]['limits']['clamped']}
            self.assertIsNone(clamped['bounds.time_budget_ms']['requested'])
            self.assertLess(clamped['bounds.time_budget_ms']['effective'], 9500)
            self.assertEqual('ended', self._wait_ended(server, 'it-finish-2')['state'])
            ret = server.get('requests/it-finish-1')
            self.assertEqual(404, ret.status_code, ret.text[:2000])
            self.assertIn('at most 1', ret.json()['error'])
            # a malformed request_id is the request's 400 (nothing registered)
            ret = server.post('resolve', dict(body, request_id='a b'))
            self.assertEqual(400, ret.status_code, ret.text[:2000])
            self.assertIn('request.request_id', ret.json()['error'])
        finally:
            server.stop()

    def test_a_short_budget_bounds_the_cutoff_before_a_later_finish_by_ms(self):
        """The cutoff comes from the EFFECTIVE deadline, never from finish_by_ms alone: a short
        time_budget_ms (1500 ms) under a finish_by_ms a minute later cuts b_big at arrival +
        1000 ms, the completed a_kleb written beside it within the budget"""
        server = self._cutoff_server('short')
        try:
            finish_by = int(time.time() * 1000) + 60000
            started = time.time()
            out = self.post_to(server, 'pattern', dict(self.SLOW_PATTERN, time_budget_ms=1500,
                                                       finish_by_ms=finish_by))
            self.assertLess(time.time() - started, 5)
            self.assertEqual([('a_kleb', 'answered'), ('b_big', 'unfinished'),
                              ('c_tiny', 'not_reached')],
                             [(a['graph'], a.get('refusal', {}).get('code', a['outcome']))
                              for a in out['answers']])
            self.assertEqual(finish_by, out['deadline']['finish_by_ms'])
            effective = out['deadline']['effective_deadline_at']
            cutoff = out['deadline']['cutoff_at']
            self.assertLess(cutoff, effective)
            # the effective deadline is the budget's, a minute before finish_by_ms
            self.assertLess(effective, time.strftime('%Y-%m-%dT%H:%M:%S',
                                                     time.gmtime(finish_by / 1000 - 30)))
        finally:
            server.stop()

    def test_a_client_gone_is_recorded_ended(self):
        """A request with a request_id whose client leaves before the answer (it half-closes
        its connection right after the request): the work stops as before, nothing is
        written, and GET /requests states the end: reason client_gone, response not written,
        the started pair with its usage. A repeat of the id starts no work (409)."""
        server = self._cutoff_server('gone')
        try:
            from test_pattern import half_closed_request
            body = dict(self.SLOW_PATTERN, time_budget_ms=60000, graphs=['b_big'],
                        request_id='it-gone-1')
            self.assertEqual(b'', half_closed_request(server.port, 'pattern', body))
            entry = self._wait_ended(server, 'it-gone-1')
            self.assertEqual(('ended', 'client_gone'), (entry['state'], entry['reason']))
            self.assertEqual({'status': None, 'bytes': None, 'written': False},
                             entry['response'])
            (pair,) = entry['pairs']
            self.assertEqual(('b_big', 'ended'), (pair['graph'], pair['state']))
            self.assertIsNotNone(pair['usage']['finished_at'])
            ret = server.post('pattern', body)
            self.assertEqual(409, ret.status_code, ret.text[:2000])
            self.assertEqual('client_gone', ret.json()['reason'])
        finally:
            server.stop()

    def test_not_after_ms_admits_resolve_and_pattern(self):
        """not_after_ms (the latest START, E20) on /resolve and /pattern: passed when the
        handler starts, 409 {error, state: expired, not_after_ms, server_time_ms} and nothing
        runs (nothing registered under its request_id); not passed, the request runs (one
        that runs past it: test_an_admitted_request_runs_past_not_after_ms). A malformed
        request_id, not_after_ms or finish_by_ms is the request's 400."""
        past = int(time.time() * 1000) - 1000
        kleb = self.records['kleb'][0]
        for route, body in (
                ('resolve', {'graph': 'kleb', 'sequence': kleb[3000:3100],
                             'discover': {'max_labels': 2}}),
                ('pattern', {'patterns': [{'dna': kleb[3000:3020]}], 'mode': 'count',
                             'graphs': ['kleb', 'pseudo']})):
            with self.subTest(route=route):
                ret = self.server.post(route, dict(body, not_after_ms=past,
                                                   request_id=f'it-expired-{route}'))
                self.assertEqual(409, ret.status_code, ret.text[:2000])
                out = ret.json()
                self.assertEqual(('expired', past, f'it-expired-{route}'),
                                 (out['state'], out['not_after_ms'], out['request_id']))
                self.assertGreater(out['server_time_ms'], past)
                self.assertIn('has passed on this server\'s clock', out['error'])
                self.assertEqual(404, self.server.get(f'requests/it-expired-{route}').status_code)
        # not passed: the request runs, the deadline block stating it
        later = int(time.time() * 1000) + 60000
        for route, body in (
                ('resolve', {'graphs': ['kleb', 'pseudo'], 'sequence': kleb[3000:3100],
                             'discover': {'max_labels': 2}}),
                ('pattern', {'patterns': [{'dna': kleb[3000:3020]}], 'mode': 'count',
                             'graphs': ['kleb', 'pseudo']})):
            out = self.post(route, dict(body, not_after_ms=later))
            self.assertEqual((2, 0, later), (out['answered'], out['refused'],
                                             out['deadline']['not_after_ms']))
        # malformed: the request's 400
        for field, value in (('not_after_ms', -1), ('not_after_ms', 1.5),
                             ('finish_by_ms', 'soon'), ('request_id', '')):
            with self.subTest(field=field, value=value):
                ret = self.server.post('resolve', {'graph': 'kleb', 'sequence': kleb[:100],
                                                   'discover': {'max_labels': 2}, field: value})
                self.assertEqual(400, ret.status_code, ret.text[:2000])
                self.assertIn(f'request.{field}', ret.json()['error'])
                ret = self.server.post('pattern', {'patterns': [{'dna': kleb[:20]}],
                                                   'graphs': ['kleb'], field: value})
                self.assertEqual(400, ret.status_code, ret.text[:2000])
                self.assertIn(f'request.{field}', ret.json()['error'])

    def test_an_admitted_request_runs_past_not_after_ms(self):
        """not_after_ms is admission only (E20): a request admitted before it runs past it and
        finishes -- b_big's resolve of 3 Mbp, over two seconds, admitted 800 ms before
        not_after_ms, alone (inline) and beside c_tiny; c_tiny, queued behind b_big on the one
        pool thread, was not started by not_after_ms: not_reached, reason not_after_ms (decided
        when not_after_ms passed, not at the cutoff)."""
        server = self._cutoff_server('admission')
        try:
            request = self._slow_resolve(3000000)
            for graphs in (['b_big'], ['b_big', 'c_tiny']):
                with self.subTest(graphs=graphs):
                    not_after = int(time.time() * 1000) + 800
                    out = self.post_to(server, 'resolve', dict(request, graphs=graphs,
                                                               not_after_ms=not_after))
                    self.assertGreater(time.time() * 1000, not_after + 500)
                    big = out['answers'][0]
                    self.assertEqual(('b_big', 'answered'), (big['graph'], big['outcome']))
                    self.assertGreater(big['timing']['elapsed_ms'], 1000)
                    self.assertEqual((not_after, None),
                                     (out['deadline']['not_after_ms'],
                                      out['deadline']['finish_by_ms']))
                    self.assertIsNotNone(out['deadline']['effective_deadline_at'])
                    self.assertEqual([], out['unfinished'])
                    if len(graphs) == 1:
                        self.assertEqual((1, 0, []), (out['answered'], out['refused'],
                                                      out['not_reached']))
                        continue
                    tiny = out['answers'][1]
                    self.assertEqual((1, 1, [self._ref(tiny)]),
                                     (out['answered'], out['refused'], out['not_reached']))
                    self.assertEqual({'http_status': 503, 'code': 'not_reached',
                                      'reason': 'not_after_ms', 'error': tiny['refusal']['error']},
                                     tiny['refusal'])
                    self.assertIn(f'not_after_ms {not_after} had passed', tiny['refusal']['error'])
        finally:
            server.stop()

    def test_traverse_takes_finish_by_ms_and_answers_get_requests(self):
        """/traverse reads finish_by_ms (E20): the attempt's bound is at most finish_by_ms less
        the clock skew allowance (usage.bound.capped_by "finish_by_ms"); its attempt_id is
        the request's handle: GET /requests/{attempt_id} answers from the attempt registry
        (route /traverse, attempt: GET /traverse/attempt/{attempt_id}'s object, the pair)."""
        seed = self.records['ecoli'][0][20000:20100]
        t = {'graphs': ['ecoli'], 'seeds': [{'sequence': seed}], 'attempt_id': 'it-trv-1',
             'finish_by_ms': int(time.time() * 1000) + 30000,
             'strategy': {'direction': 'right', 'bounds': {'max_extension_bp': 50},
                          'output': {'detail': 'summary'}}}
        out = self.post('traverse', t)
        self.assertEqual('finish_by_ms', out['usage']['bound']['capped_by'])
        self.assertLessEqual(out['usage']['bound_ms'], 30000)
        deadline = time.time() + 30
        while True:
            entry = self.server.get('requests/it-trv-1').json()
            if entry['state'] == 'ended' or time.time() > deadline:
                break
            time.sleep(0.05)
        self.assertEqual(('it-trv-1', '/traverse', 'ended'),
                         (entry['request_id'], entry['route'], entry['state']))
        self.assertEqual('it-trv-1', entry['attempt']['attempt_id'])
        self.assertEqual(self.server.get('traverse/attempt/it-trv-1').json()['usage'],
                         entry['attempt']['usage'])
        self.assertEqual([('ecoli', self.pairs['ecoli'][0], self.pairs['ecoli'][1], 'answered')],
                         [(p['graph'], p['graph_path'], p['annotation_path'], p['outcome'])
                          for p in entry['pairs']])
        self.assertEqual({'status': 200, 'written': True},
                         {k: entry['response'][k] for k in ('status', 'written')})

    # -------------------------------------------------------------------- per-task figures

    def test_pattern_envelope_states_the_cpu_of_its_pairs(self):
        """The envelope's timing.cpu_ms (SPEC-pattern-search.md §24.1, §8.4): its own thread's
        part and every pair's work on its pool thread, so at least the answered pairs' own
        cpu_ms together; page_faults on Linux only, likewise summed; work.memory_bytes the
        answered pairs' together."""
        ecoli = self.records['ecoli'][0]
        out = self.post('pattern', {'patterns': [{'dna': ecoli[5000:5020]}], 'mode': 'count',
                                    'graphs': ['enterics', 'pseudo']})
        answered = [a for a in out['answers'] if a['outcome'] == 'answered']
        self.assertEqual(3, len(answered))
        linux = sys.platform.startswith('linux')
        for t in [out['timing']] + [a['timing'] for a in answered]:
            self.assertGreater(t['cpu_ms'], 0)
            self.assertEqual(linux, 'page_faults' in t)
        self.assertGreaterEqual(out['timing']['cpu_ms'],
                                sum(a['timing']['cpu_ms'] for a in answered))
        # work.memory_bytes: the answered pairs' accounts together (each its own)
        self.assertEqual(['memory_bytes'], list(out['work']))
        self.assertEqual(out['work']['memory_bytes'],
                         sum(a['work']['memory_bytes'] for a in answered))
        if linux:
            for k in ('minor', 'major'):
                self.assertGreaterEqual(out['timing']['page_faults'][k],
                                        sum(a['timing']['page_faults'][k] for a in answered))

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
        # the envelope's one entry (SPEC-multigraph-envelope.md section 1)
        self.assertNotIn('load_ms', resident['answers'][0]['timing'])
        self.assertGreater(loaded['answers'][0]['timing']['load_ms'], 0)
        self.assertIn('loaded into RAM for it', self.server.text())
        # /search as it always was
        query = {'FASTA': f'>q\n{seed}', 'graphs': ['ecoli'], 'discovery_fraction': 0.5}
        self.assertEqual(self.post('search', query), self.post('search', dict(query, in_ram=True)))

    def test_in_ram_shares_the_resident_header_index(self):
        """The copy an `in_ram` load makes is built with the resident pair's record mapping (the
        .seqs: the sequence headers and their reverse index, built once at start-up), so a load
        builds no header index: the log's one build line (ecoli, the pair with a mapping, at
        start-up) stays the only one, and requests with header labels on /resolve, /traverse and
        /pattern answer with in_ram as without, timing apart."""
        built = 'Sequence header index built'
        log = self.server.text()
        self.assertEqual(1, log.count(built), log[-3000:])
        self.assertLess(log.index(built), log.index('All graphs were loaded'))
        header, record = self.named_records['ecoli'][0]
        seed = record[20000:20100]

        r = {'graphs': ['ecoli'], 'sequence': seed, 'labels': [header]}
        resident = self.post('resolve', r)
        loaded = self.post('resolve', dict(r, in_ram=True))
        self.assertEqual([(header, 'header')],
                         [(l['label'], l['kind']) for l in resident['answers'][0]['labels']])
        self.assertEqual(untimed(resident), untimed(loaded))
        self.assertGreater(loaded['answers'][0]['timing']['load_ms'], 0)
        d = {'graphs': ['ecoli'], 'sequence': seed,
             'discover': {'max_labels': 10, 'kind': 'header'}}
        resident = self.post('resolve', d)
        self.assertIn(header, [l['label'] for l in resident['answers'][0]['labels']])
        self.assertEqual(untimed(resident), untimed(self.post('resolve', dict(d, in_ram=True))))

        t = {'graphs': ['ecoli'], 'seeds': [{'sequence': seed, 'labels': [header]}],
             'strategy': {'direction': 'right', 'bounds': {'max_extension_bp': 50},
                          'output': {'detail': 'full'}}}
        resident = self.post('traverse', t)
        loaded = self.post('traverse', dict(t, in_ram=True))
        self.assertEqual([(header, 'header')],
                         [(l['name'], l['kind']) for l in resident['results'][0]['label_dict']])
        self.assertEqual(untimed(resident['results']), untimed(loaded['results']))
        self.assertGreater(loaded['timing']['load_ms'], 0)

        # /pattern's copy shares it too: the record labels of its paths
        p = {'patterns': [{'dna': record[30000:30040]}], 'long_search': 'paths',
             'output': {'labels': 'all'}, 'graphs': ['ecoli']}
        resident = self.post('pattern', p)
        loaded = self.post('pattern', dict(p, in_ram=True))
        self.assertEqual(untimed(resident), untimed(loaded))
        self.assertGreater(loaded['answers'][0]['timing']['load_ms'], 0)
        results = resident['answers'][0]['patterns'][0]['results']
        self.assertTrue(results, resident)
        self.assertEqual('record_verified', results[0]['support'], results[0])

        # every load above shared the mapping: nothing was built after start-up
        log = self.server.text()
        self.assertEqual(1, log.count(built), log[-3000:])
        self.assertGreaterEqual(log.count('loaded into RAM for it'), 4)

    def test_in_ram_on_a_pair_without_its_record_mapping(self):
        """A coordinate-annotated pair served without its .seqs (the ecoli chunk's files by
        symlinks, none to the .seqs; a served configuration): the copy an `in_ram` load makes
        has the mapping the resident has -- none --, so the operator's notice that the mapping
        is missing is the start-up's, logged once, not once per request; the answers with
        in_ram are the resident's, timing apart."""
        d = os.path.join(self.tempdir.name, 'bare')
        os.makedirs(d)
        src = os.path.dirname(self.pairs['ecoli'][0])
        for f in os.listdir(src):
            if not f.endswith('.seqs'):
                os.symlink(os.path.join(src, f), os.path.join(d, f))
        graph, anno = (os.path.join(d, os.path.basename(p)) for p in self.pairs['ecoli'])
        csv = os.path.join(d, 'graphs.csv')
        with open(csv, 'w') as f:
            f.write(f'bare,{graph},{anno}\n')
        server = Server([csv, '--mmap', '--mem-cap-gb', '1'], os.path.join(d, 'server.log'))
        try:
            notice = 'No CoordToHeader mapping found'
            self.assertEqual(1, server.text().count(notice), server.text()[-3000:])
            pair = server.get('capabilities').json()['graph_summary']['pairs'][0]
            self.assertEqual(('bare', True, False),
                             (pair['graph'], pair['traversal']['has_coordinates'],
                              pair['traversal']['has_coord_to_header']))
            ecoli = self.records['ecoli'][0]
            seed = ecoli[20000:20100]
            bodies = [
                ('resolve', {'graphs': ['bare'], 'sequence': seed,
                             'discover': {'max_labels': 10}}),
                ('pattern', {'patterns': [{'dna': ecoli[5000:5020]}], 'mode': 'count',
                             'graphs': ['bare']}),
                ('traverse', {'graphs': ['bare'], 'seeds': [{'sequence': seed}],
                              'strategy': {'direction': 'right',
                                           'bounds': {'max_extension_bp': 50},
                                           'output': {'detail': 'full'}}}),
            ]
            for route, body in bodies:
                for _ in range(2):
                    ret = server.post(route, body)
                    self.assertEqual(200, ret.status_code, ret.text[:2000])
                    resident = ret.json()
                    ret = server.post(route, dict(body, in_ram=True))
                    self.assertEqual(200, ret.status_code, ret.text[:2000])
                    loaded = ret.json()
                    if route == 'traverse':
                        resident, loaded = resident['results'], loaded['results']
                    self.assertEqual(untimed(resident), untimed(loaded), route)
            log = server.text()
            self.assertEqual(6, log.count('loaded into RAM for it'), log[-3000:])
            self.assertEqual(1, log.count(notice), log[-3000:])
        finally:
            server.stop()

    @unittest.skipIf(hasattr(os, 'geteuid') and os.geteuid() == 0,
                     "root reads a file without permissions")
    def test_a_failed_in_ram_load_is_answered_and_the_server_lives(self):
        """A pair whose annotation cannot be read any more after start-up (chmod 000; the resident
        index, mapped at start-up, still answers): an `in_ram` request on it is answered, not the
        process ended -- /pattern's pair refused with the loader's message (400, no code, the
        request 200, and so /resolve's), /traverse 400 with it --, its memory released, and the server
        lives: the pair without in_ram, the other pair with it (above the cap: served resident)
        and, once readable again, the pair with in_ram answer, the last proving the release (the
        cap holds one copy of it: a reservation a failed load kept would make it wait). A malformed
        body with in_ram: true naming the pair is refused before any load (SPEC §24.2): the
        server's line of a failed load, which every attempt to load the unreadable file logs,
        not logged again."""
        d = os.path.join(self.tempdir.name, 'unreadable')
        pairs = {}
        for name in ('kleb', 'pseudo'):
            dst = os.path.join(d, name)
            shutil.copytree(os.path.dirname(self.pairs[name][0]), dst)
            pairs[name] = tuple(os.path.join(dst, os.path.basename(p)) for p in self.pairs[name])
        csv = os.path.join(d, 'graphs.csv')
        with open(csv, 'w') as f:
            for name, (graph, anno) in pairs.items():
                f.write(f'{name},{graph},{anno}\n')
        kleb_bytes = sum(os.path.getsize(p) for p in pairs['kleb'])
        pseudo_bytes = sum(os.path.getsize(p) for p in pairs['pseudo'])
        # one copy of kleb fits the cap, pseudo does not
        cap_gb = kleb_bytes * 1.2 / 1e9
        self.assertGreater(pseudo_bytes, cap_gb * 1e9)
        server = Server([csv, '--mmap', '--mem-cap-gb', f'{cap_gb:.9f}'],
                        os.path.join(d, 'server.log'))
        anno = pairs['kleb'][1]
        failed = f"Cannot load annotation from '{anno}'"
        kleb = self.records['kleb'][0]
        seed = kleb[20000:20100]
        p = {'patterns': [{'dna': kleb[5000:5020]}], 'mode': 'count', 'graphs': ['kleb']}
        t = {'graphs': ['kleb'], 'seeds': [{'sequence': seed}],
             'strategy': {'direction': 'right', 'bounds': {'max_extension_bp': 50},
                          'output': {'detail': 'full'}}}
        r = {'graphs': ['kleb'], 'sequence': seed, 'discover': {'max_labels': 10}}
        try:
            os.chmod(anno, 0)
            # the resident index answers
            self.assertEqual(1, server.post('pattern', p).json()['answered'])
            # /pattern: the pair refused with the loader's message, the request answered
            ret = server.post('pattern', dict(p, in_ram=True))
            self.assertEqual(200, ret.status_code, ret.text[:2000])
            body = ret.json()
            self.assertEqual((0, 1), (body['answered'], body['refused']), body)
            entry = body['answers'][0]
            self.assertEqual(('kleb', 'refused'), (entry['graph'], entry['outcome']), entry)
            self.assertEqual({'http_status', 'error'}, set(entry['refusal']), entry)
            self.assertEqual(400, entry['refusal']['http_status'])
            self.assertIn(failed, entry['refusal']['error'])
            self.assertIn('read permission', entry['refusal']['error'])
            # /traverse: 400 with the message; with usage on a managed attempt
            ret = server.post('traverse', dict(t, in_ram=True))
            self.assertEqual(400, ret.status_code, ret.text[:2000])
            self.assertEqual({'error'}, set(ret.json()))
            self.assertIn(failed, ret.json()['error'])
            ret = server.post('traverse', dict(t, in_ram=True, attempt_id='in-ram-unreadable'))
            self.assertEqual(400, ret.status_code, ret.text[:2000])
            self.assertEqual({'error', 'usage'}, set(ret.json()))
            self.assertIn(failed, ret.json()['error'])
            self.assertEqual('error', ret.json()['usage']['reason'])
            # /resolve: the pair refused with the message (400, no code), the envelope 200
            ret = server.post('resolve', dict(r, in_ram=True))
            self.assertEqual(200, ret.status_code, ret.text[:2000])
            body = ret.json()
            self.assertEqual((0, 1, '/resolve'), (body['answered'], body['refused'],
                                                  body['route']), body)
            entry = body['answers'][0]
            self.assertEqual(('kleb', 'refused'), (entry['graph'], entry['outcome']), entry)
            self.assertEqual({'http_status', 'error'}, set(entry['refusal']), entry)
            self.assertEqual(400, entry['refusal']['http_status'])
            self.assertIn(failed, entry['refusal']['error'])
            log = server.text()
            # one line per load attempt: the server's, with the message (the annotation
            # does not open, so the load fails before any loader or the graph's load starts)
            self.assertEqual(4, log.count(failed), log[-3000:])
            self.assertEqual(0, log.count('loaded into RAM for it'), log[-3000:])

            # a malformed body with in_ram naming the pair: refused before any load
            for route, payload, words in (
                    ('pattern', dict(p, in_ram=True, mode='x'), 'mode'),
                    ('traverse', dict(t, in_ram=True, seeds=[]), 'request.seeds'),
                    ('resolve', dict(r, in_ram=True, sequence=''), 'request.sequence')):
                with self.subTest(route=route):
                    ret = server.post(route, payload)
                    self.assertEqual(400, ret.status_code, ret.text[:2000])
                    self.assertIn(words, ret.json()['error'])
                    self.assertNotIn(failed, ret.text)
            log = server.text()
            self.assertEqual(4, log.count(failed), log[-3000:])
            self.assertEqual(0, log.count('loaded into RAM for it'), log[-3000:])

            # the server lives: the pair without in_ram, the other pair with it (above the cap)
            self.assertEqual(1, server.post('pattern', p).json()['answered'])
            pseudo = self.records['pseudo'][0]
            q = {'patterns': [{'dna': pseudo[5000:5020]}], 'mode': 'count', 'graphs': ['pseudo'],
                 'in_ram': True}
            ret = server.post('pattern', q)
            self.assertEqual(200, ret.status_code, ret.text[:2000])
            self.assertEqual(1, ret.json()['answered'], ret.text[:2000])
            self.assertEqual(0, ret.json()['answers'][0]['timing']['load_ms'])
            self.assertIn('too large to fit into RAM', server.text())
        finally:
            os.chmod(anno, 0o644)
        try:
            # readable again: the load succeeds, the memory of the failed loads was released
            # (the cap holds one copy; a kept reservation would make this request wait)
            ret = requests.post(server.url('pattern'), data=json.dumps(dict(p, in_ram=True)),
                                timeout=120)
            self.assertEqual(200, ret.status_code, ret.text[:2000])
            self.assertEqual(1, ret.json()['answered'], ret.text[:2000])
            self.assertGreater(ret.json()['answers'][0]['timing']['load_ms'], 0)
            self.assertEqual(1, server.text().count('loaded into RAM for it'))
            self.assertIsNone(server.process.poll())
        finally:
            server.stop()

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
            via_graphs = self.server.post('traverse', dict(request, graphs=[chunk]))
            self.assertEqual(200, via_graphs.status_code, via_graphs.text[:2000])
            out = via_graphs.json()
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
            # graph, the older spelling: the same bytes, the spelling not echoed (the request
            # was a 400 on its first absent seed before round D)
            via_graph = self.server.post('traverse', dict(request, graph=chunk))
            self.assertEqual(200, via_graph.status_code, via_graph.text[:2000])
            self.assertEqual(via_graphs.content, via_graph.content, chunk)
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

    # -------------------------------------------------------------------- POST /resolve

    def cli_resolve(self, chunk, request):
        """`metagraph traverse --resolve` on the chunk's files: the answer a single-graph server
        gives"""
        path = os.path.join(self.tempdir.name, f'resolve_{chunk}.json')
        with open(path, 'w') as f:
            json.dump(request, f)
        graph, anno = self.pairs[chunk]
        res = subprocess.run(shlex.split(METAGRAPH) + ['traverse', '--resolve', '--json', '-i',
                                                       graph, '-a', anno, path],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(0, res.returncode, res.stderr.decode()[-2000:])
        return json.loads(res.stdout.decode().strip().split('\n')[-1])

    def test_resolve_answers_the_envelope(self):
        """POST /resolve on a multi-graph server (SPEC-multigraph-envelope.md section 1, E2):
        the envelope, also for a named graph (one entry, `graph` and `graphs: [name]` one
        request); `graphs` naming several, their pairs in the index's order; a name over two
        graphs each pair (graph_path still picks one); a request without a selection every
        name of this small server (no 'graph required' refusal any more). Each entry is the
        pair's single-graph answer (the CLI's on the pair's files), tagged; work.memory_bytes
        the answered entries' together; request_id stated and its request registered."""
        seed = self.records['kleb'][0][10000:10200]
        r = {'sequence': seed, 'discover': {'max_labels': 10}}
        named = self.post('resolve', dict(r, graph='kleb'))
        self.assertEqual(untimed(named), untimed(self.post('resolve', dict(r, graphs=['kleb']))))
        self.assertEqual(ENVELOPE, set(named))
        self.assertEqual((1, '/resolve', ['kleb'], 1, 0, [], [], None, NO_DEADLINE),
                         (named['envelope_version'], named['route'], named['graphs'],
                          named['answered'], named['refused'], named['not_reached'],
                          named['unfinished'], named['request_id'], named['deadline']))
        self.assertRegex(named['server_instance'], '^[0-9a-f]{16}$')
        (entry,) = named['answers']
        self.assertEqual(('kleb', *self.pairs['kleb'], None, 'answered'),
                         tuple(entry[k] for k in TAGS))
        self.assertTrue(entry['labels'], entry)
        self.assertEqual(entry['work']['memory_bytes'], named['work']['memory_bytes'])
        tagless = {k: v for k, v in entry.items() if k not in TAGS}
        self.assertEqual(untimed(self.cli_resolve('kleb', r)), untimed(tagless))
        # several names: their pairs in byte order, each its named request's entry
        several = self.post('resolve', dict(r, graphs=['kleb', 'ecoli'],
                                            request_id='it-resolve-1'))
        self.assertEqual((['ecoli', 'kleb'], 2, 0, 'it-resolve-1'),
                         (several['graphs'], several['answered'], several['refused'],
                          several['request_id']))
        self.assertEqual(untimed(self.post('resolve', dict(r, graph='ecoli'))['answers']),
                         untimed(several['answers'][:1]))
        self.assertEqual(untimed(entry), untimed(several['answers'][1]))
        self.assertEqual(sum(a['work']['memory_bytes'] for a in several['answers']),
                         several['work']['memory_bytes'])
        self.assertGreater(several['timing']['cpu_ms'], 0)
        state = self.server.get('requests/it-resolve-1').json()
        self.assertEqual(('/resolve', [('ecoli', 'answered'), ('kleb', 'answered')]),
                         (state['route'], [(p['graph'], p['outcome']) for p in state['pairs']]))
        # a name over two graphs: each pair; graph_path picks one
        both = self.post('resolve', dict(r, graph='enterics'))
        self.assertEqual([('enterics', self.pairs['ecoli'][0]),
                          ('enterics', self.pairs['kleb'][0])],
                         [(a['graph'], a['graph_path']) for a in both['answers']])
        self.assertEqual(untimed(tagless),
                         untimed({k: v for k, v in both['answers'][1].items() if k not in TAGS}))
        picked = self.post('resolve', dict(r, graph='enterics', graph_path=self.pairs['kleb'][0]))
        self.assertEqual([self.pairs['kleb'][0]], [a['graph_path'] for a in picked['answers']])
        # no selection: every name of this small server, each name's pairs
        every = self.post('resolve', r)
        self.assertEqual(['ecoli', 'enterics', 'kleb', 'pseudo'], every['graphs'])
        self.assertEqual(['ecoli', 'enterics', 'enterics', 'kleb', 'pseudo'],
                         [a['graph'] for a in every['answers']])
        self.assertEqual((5, 0), (every['answered'], every['refused']))
        # what no pair decides is the request's 400
        for payload, words in ((dict(r, graphs=['nope']), "unknown graph 'nope'"),
                               (dict(r, graph='nope'), "unknown graph 'nope'"),
                               (dict(r, sequence=''), 'request.sequence'),
                               (dict(r, bounds={'time_budget_ms': 1}), 'time_budget_ms')):
            with self.subTest(payload=payload):
                ret = self.server.post('resolve', payload)
                self.assertEqual(400, ret.status_code, ret.text[:2000])
                self.assertIn(words, ret.json()['error'])


if __name__ == '__main__':
    unittest.main()
