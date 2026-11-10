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
  - GET /capabilities: the per-pair summary (graph_summary) and columns_disjoint, each pair's
    capability profile (E17: index_meta_fp at start-up, alphabet, graph_representation,
    annotation_form, decode_charged, traversal.counts and support_fp);
  - GET /pattern/capabilities?graph= and /traverse/capabilities?graph=: the pair's block, in
    the envelope of one entry (SPEC-multigraph-envelope.md section 1), and every pair's without
    ?graph=; a /traverse answers in the envelope of its one pair.

The answers are checked against oracles of their own: the CLI (`metagraph pattern`) on the same
files, the k-mers of the chunk's records (the contexts of a pattern, a seed's k-mers present),
and the server's answers without the new fields.
"""

import json
import os
import random
import re
import shlex
import shutil
import signal
import socket
import subprocess
import sys
import tempfile
import threading
import time
import unittest

import requests

from base import PROTEIN_MODE, TestingBase, METAGRAPH
from test_pattern import gate_block
from test_traverse_forms import CHARGED

MINI_DIR = os.environ.get('METAGRAPH_MINI_REFSEQ', os.path.join(os.getcwd(), 'mini_refseq'))
FASTA = os.path.join(MINI_DIR, 'fasta')
# $METAGRAPH_REQUIRE_GUARDS=1: a missing mini index fails the class instead of skipping it
REQUIRE_GUARDS = os.environ.get('METAGRAPH_REQUIRE_GUARDS', '') == '1'
# the release-9 binary (8fa32e553, before release 10), whose refusal texts a multi-graph
# client matches whole: compared with where given (not in CI)
BASE_BINARY = os.environ.get('METAGRAPH_BASE_BINARY', '')
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


def strip_column_labels(seqs):
    """|seqs| (a .seqs written since H4) rewritten as an earlier binary writes it: its columns
    without the trailer of column labels (the marker "mtg-cols", then the labels), which an H4
    load reads as a mapping that states none (seqs_order "count_only")"""
    with open(seqs, 'rb') as f:
        data = f.read()
    at = data.rfind(b'mtg-cols')
    assert at > 0, seqs
    with open(seqs, 'wb') as f:
        f.write(data[:at])


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


SPEC_ENVELOPE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'docs',
                             'SPEC-multigraph-envelope.md')


def support_rule():
    """The version line and the fields of support_fp, as SPEC-multigraph-envelope.md section 8.1
    states them (its fenced block: the version line, then "<name>=<value>" per field)"""
    with open(SPEC_ENVELOPE) as f:
        lines = f.read().split('\n')
    start = next(i for i, line in enumerate(lines) if line.startswith('mgt-support '))
    end = next(i for i in range(start, len(lines)) if lines[i].startswith('```'))
    fields = [line[:-len('=<value>')] for line in lines[start + 1:end]]
    assert all(line.endswith('=<value>') for line in lines[start + 1:end]), lines[start:end]
    return lines[start], fields


def support_fp(pair):
    """support_fp recomputed from a graph_summary entry (or pair_summary) alone, by the rule of
    SPEC-multigraph-envelope.md section 8.1: FNV-1a-64 over the version line and a line
    "<name>=<compact JSON>" per field, as 16 lowercase hex digits"""
    version, fields = support_rule()
    text = version + '\n'
    for name in fields:
        value = pair
        for part in name.split('.'):
            value = value.get(part) if isinstance(value, dict) else None
        text += f'{name}={json.dumps(value)}\n'
    h = 14695981039346656037
    for byte in text.encode():
        h = ((h ^ byte) * 1099511628211) & 0xFFFFFFFFFFFFFFFF
    return '%016x' % h


def single(out, route=None):
    """The one entry of a multi-graph server's envelope (SPEC-multigraph-envelope.md section
    1): an answered pair, the single-graph answer with its tags; on /traverse the attempt's
    usage is the envelope's, set beside it"""
    assert out['envelope_version'] == 1 and len(out['answers']) == 1, out
    assert route is None or out['route'] == route, out['route']
    entry = dict(out['answers'][0])
    assert entry['outcome'] == 'answered', entry
    if 'usage' in out:
        entry['usage'] = out['usage']
    return entry


def untimed(value):
    if isinstance(value, dict):
        return {k: untimed(v) for k, v in value.items() if k != 'timing'}
    if isinstance(value, list):
        return [untimed(v) for v in value]
    return value


class Server:
    """A server_query process on a free port, ready once GET /stats answers 200; |threads| is
    its -p (the graph pool the pairs of a request run on, and its HTTP threads)."""

    def __init__(self, args, log_path, timeout=300, threads=4, binary=None):
        os.environ['NO_PROXY'] = '127.0.0.1'
        for _ in range(20):
            self.port = free_port()
            self.log_path = log_path
            self.log = open(log_path, 'ab')
            self.process = subprocess.Popen(
                shlex.split(binary or METAGRAPH) + ['server_query'] + args
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
        # no pair has a manifest: no composite fingerprint (SPEC-multigraph-envelope.md
        # section 8), the multigraph block stated
        self.assertIsNone(caps['composite_fp'])
        self.assertEqual(1, caps['multigraph']['envelope_version'])
        self.assertEqual('GET /requests/{request_id}', caps['routes']['requests'])
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
                                  'supports_trace': masked, 'counts': None}, p['traversal'])
                # what the pair's own documents say
                block = single(self.server.get(f'pattern/capabilities?graph={chunk}').json())
                for key in ('k', 'graph_mode', 'available', 'mask', 'counting'):
                    self.assertEqual(block[key], p[key], key)
                probe = single(self.server.get(f'traverse/capabilities?graph={chunk}').json())
                for key in ('regime', 'num_labels', 'has_coordinates', 'has_coord_to_header',
                            'supports_trace'):
                    self.assertEqual(probe[key], p['traversal'][key], key)

    def test_capabilities_state_each_pairs_profile(self):
        """Each pair's capability profile (SPEC-multigraph-envelope.md section 8.1): its
        index_meta_fp computed at start-up (before any request read the pair: the probes below
        come after) and its probe's, the alphabet its probe states, the succinct graph, the
        annotation form, charged reads on both forms (row_diff_brwt_coord; plain column), no
        counts, and support_fp recomputed from the entry: kleb and pseudo, two chunks of one
        build (unmasked column annotations by file name), share it while their index_meta_fp
        differ; ecoli (masked, coordinates, its .seqs) has another; enterics' entries are
        ecoli's and kleb's"""
        caps = self.server.get('capabilities').json()
        self.assertNotIn('pair_summary', caps)
        log = self.server.text()
        self.assertRegex(log, r'index_meta_fp of the 3 pairs \(\d+ columns\) computed in '
                              r'\d+\.\d+ s')
        self.assertIn('2 capability profile(s) (support_fp)', log)
        by_chunk = {}
        for p in caps['graph_summary']['pairs']:
            chunk = {v[0]: c for c, v in self.pairs.items()}[p['graph_path']]
            with self.subTest(pair=(p['graph'], chunk)):
                probe = single(self.server.get(f'traverse/capabilities?graph={chunk}').json())
                self.assertEqual(probe['index_meta_fp'], p['index_meta_fp'])
                self.assertRegex(p['index_meta_fp'], r'^[0-9a-f]{16}$')
                self.assertEqual(probe['alphabet'], p['alphabet'])
                self.assertEqual('succinct', p['graph_representation'])
                self.assertEqual(dict((c[0], c[3]) for c in CHUNKS)[chunk], p['annotation_form'])
                self.assertIs(True, p['decode_charged'])
                self.assertIsNone(p['traversal']['counts'])
                self.assertNotIn('counts', probe)
                self.assertRegex(p['support_fp'], r'^[0-9a-f]{16}$')
                self.assertEqual(support_fp(p), p['support_fp'])
                if p['graph'] == chunk:
                    by_chunk[chunk] = p
        for p in caps['graph_summary']['pairs']:
            if p['graph'] == 'enterics':
                chunk = {v[0]: c for c, v in self.pairs.items()}[p['graph_path']]
                self.assertEqual(by_chunk[chunk], dict(p, graph=chunk))
        self.assertEqual(by_chunk['kleb']['support_fp'], by_chunk['pseudo']['support_fp'])
        self.assertNotEqual(by_chunk['kleb']['index_meta_fp'], by_chunk['pseudo']['index_meta_fp'])
        self.assertNotEqual(by_chunk['kleb']['traversal']['num_labels'],
                            by_chunk['pseudo']['traversal']['num_labels'])
        self.assertNotEqual(by_chunk['ecoli']['support_fp'], by_chunk['kleb']['support_fp'])

    def test_capabilities_probes_without_a_graph(self):
        """GET /traverse/capabilities and GET /pattern/capabilities without ?graph= answer the
        envelope over every pair, in graph_summary's order, each entry the block the probe of
        that pair alone answers (enterics' pairs are ecoli's and kleb's under another name)"""
        summary = self.server.get('capabilities').json()['graph_summary']['pairs']
        for route in ('pattern/capabilities', 'traverse/capabilities'):
            out = self.server.get(route).json()
            self.assertEqual(('/' + route, ['ecoli', 'enterics', 'kleb', 'pseudo'],
                              len(summary), 0),
                             (out['route'], out['graphs'], out['answered'], out['refused']))
            self.assertEqual([(p['graph'], p['graph_path'], p['annotation_path'], p['index_fp'])
                              for p in summary],
                             [(e['graph'], e['graph_path'], e['annotation_path'], e['index_fp'])
                              for e in out['answers']])
            for entry in out['answers']:
                if entry['graph'] == 'enterics':
                    continue
                alone = single(self.server.get(f'{route}?graph={entry["graph"]}').json())
                for key in ('k', 'graph_mode', 'available', 'mask', 'counting', 'regime',
                            'num_labels', 'has_coordinates', 'index_meta_fp'):
                    self.assertEqual(alone.get(key), entry.get(key), (route, key))

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
                block = single(server.get('pattern/capabilities?graph=hashed').json())
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
        # a name of two pairs: a probe names one (SPEC-multigraph-envelope.md section 6)
        # (in the words of the release-9 server, which a client matches whole)
        self.assertEqual(('pair_ambiguous', f"Bad request: index 'enterics' spans several graphs "
                                            f"({self.pairs['ecoli'][0]}, {self.pairs['kleb'][0]}); "
                                            f"pass 'graph_path' to pick one"),
                         (enterics.json()['code'], enterics.json()['error']))
        self.assertEqual([self.pairs['ecoli'][0], self.pairs['kleb'][0]],
                         [c['graph_path'] for c in enterics.json()['candidates']])
        kleb = single(self.server.get('pattern/capabilities?graph=kleb').json(),
                      '/pattern/capabilities')
        quoted = requests.utils.quote(self.pairs['kleb'][0], safe='')
        picked = single(self.server.get(f'pattern/capabilities?graph=enterics&graph_path={quoted}')
                        .json())
        self.assertEqual(('kleb', self.pairs['kleb'][0]), (kleb.pop('graph'), kleb.pop('graph_path')))
        self.assertEqual(('enterics', self.pairs['kleb'][0]),
                         (picked.pop('graph'), picked.pop('graph_path')))
        self.assertEqual(kleb, picked)
        probe = single(self.server.get('traverse/capabilities?graph=kleb').json())
        # the pair's gate block (SPEC-pattern-search.md §23.2)
        self.assertEqual(gate_block(kleb), probe['pattern'])
        self.assertNotIn('predicate', probe['pattern'])
        for query, words in (('?graph_path=x', 'needs ?graph=<name>'),
                             ('?graph=kleb&x=1', "unknown parameter 'x'"),
                             ('?graph=nope', "unknown graph 'nope'")):
            ret = self.server.get('pattern/capabilities' + query)
            self.assertEqual(400, ret.status_code, query)
            self.assertIn(words, ret.json()['error'], query)

    def _refusal_panel(self):
        """Requests the release-9 server (8fa32e553) refused on this list, each with the
        refusal text it gave, byte for byte (a multi-graph client matches these texts whole;
        SPEC-multigraph-envelope.md section 6 adds a code beside them, never a new wording):
        (method, route, payload or query, code, text)"""
        ecoli, kleb = self.pairs['ecoli'][0], self.pairs['kleb'][0]
        seed = self.records['ecoli'][0][1000:1100]
        t = {'seeds': [{'sequence': seed}]}
        r = {'sequence': seed}
        spans = (f"Bad request: index 'enterics' spans several graphs ({ecoli}, {kleb}); pass "
                 f"'graph_path' to pick one")
        not_part = "Bad request: graph_path '/nowhere' is not part of index 'enterics'"
        unknown = "Bad request: unknown graph 'nope'"
        both = "Bad request: give the graph as 'graph' or as 'graphs', not both"
        one_name = "Bad request: 'graphs' names the one graph a traversal reads: expected [name]"
        required = "Bad request: 'graph' (index name) is required in multi-graph mode"
        panel = []
        for route in ('traverse/capabilities', 'pattern/capabilities'):
            needs = (f'Bad request: in multi-graph mode GET /{route} needs ?graph=<name> (and '
                     f'graph_path=<path> when the name spans several graphs); GET /capabilities '
                     f'lists the graphs')
            panel += [
                ('GET', route, '?graph=enterics', 'pair_ambiguous', spans),
                ('GET', route, '?graph=nope', 'graph_unknown', unknown),
                ('GET', route, '?graph=kleb&x=1', None,
                 f"Bad request: unknown parameter 'x' of GET /{route} (graph, graph_path)"),
                ('GET', route, '?graph=kleb&graph=ecoli', None,
                 "Bad request: the parameter 'graph' is given twice"),
                ('GET', route, '?graph_path=x', None, needs),
                ('GET', route, '?graph=enterics&graph_path=/nowhere', 'pair_not_found', not_part),
            ]
        for route, body in (('traverse', t), ('resolve', r)):
            panel += [
                ('POST', route, dict(body, graph='nope'), 'graph_unknown', unknown),
                ('POST', route, dict(body, graphs=['nope']), 'graph_unknown', unknown),
                ('POST', route, dict(body, graph='kleb', graphs=['kleb']), 'selection_conflict',
                 both),
                ('POST', route, dict(body, graphs=[]), None, one_name),
                ('POST', route, dict(body, graphs='kleb'), None, one_name),
                ('POST', route, dict(body, graph_path=kleb), None, required),
                ('POST', route, dict(body, graph=5), None, required),
                ('POST', route, dict(body, graph='enterics', graph_path='/nowhere'),
                 'pair_not_found', not_part),
            ]
        panel += [
            ('POST', 'traverse', dict(t, graph='enterics'), 'traverse_single_pair', spans),
            ('POST', 'traverse', dict(t, graphs=['ecoli', 'kleb']), 'traverse_single_pair',
             one_name),
            ('POST', 'traverse', t, 'traverse_single_pair', required),
            ('POST', 'search', {'FASTA': '>q\n' + seed, 'graphs': ['nope']}, 'graph_unknown',
             'Request with an uninitialized graph nope'),
            ('POST', 'pattern', {'patterns': [{'dna': seed[:20]}], 'graphs': ['nope']},
             'graph_unknown',
             "request.graphs: unknown graph 'nope' (GET /capabilities lists the graphs)"),
            ('POST', 'pattern', {'patterns': [{'dna': seed[:20]}], 'graphs': []},
             'invalid_request',
             'request.graphs: expected a non-empty array of graph names (GET /capabilities '
             'lists them)'),
        ]
        # Two faults (contract F3 of the release-10 review): the release-9 order of the checks
        # decides the words, the selection's shape, then in_ram, then the names, then the
        # fields, the request's handles (unknown fields to it) last
        in_ram = 'request.in_ram: expected a boolean'
        unknown_pattern = "request.graphs: unknown graph 'nope' (GET /capabilities lists the graphs)"
        p = {'patterns': [{'dna': seed[:20]}]}
        panel += [
            ('POST', 'resolve', dict(r, graph='nope', in_ram='x'), None, in_ram),
            ('POST', 'traverse', dict(t, graph='nope', in_ram='x'), None, in_ram),
            ('POST', 'traverse', dict(t, in_ram='x'), 'traverse_single_pair', required),
            ('POST', 'traverse', dict(t, in_ram='x', graphs=5), None, one_name),
            ('POST', 'resolve', dict(r, in_ram='x', graphs=5), None, one_name),
            ('POST', 'resolve', dict(r, graph='nope', request_id='bad id'), 'graph_unknown',
             unknown),
            ('POST', 'resolve', dict(r, graph='nope', not_after_ms='x'), 'graph_unknown',
             unknown),
            ('POST', 'resolve', dict(r, graph='kleb', sequence=5, finish_by_ms='x'), None,
             'request.sequence: expected a string'),
            ('POST', 'traverse', dict(t, graph='nope', finish_by_ms='soon'), 'graph_unknown',
             unknown),
            ('POST', 'traverse', dict(t, graph='kleb', seeds=5, finish_by_ms='soon'), None,
             'request.seeds: expected a non-empty array'),
            ('POST', 'pattern', dict(p, graphs=['nope'], not_after_ms='x'), 'graph_unknown',
             unknown_pattern),
            ('POST', 'pattern', dict(p, graphs=['nope'], request_id='bad id'), 'graph_unknown',
             unknown_pattern),
            ('POST', 'pattern', dict(p, graphs=['kleb'], patterns=[5], request_id='bad id'),
             'invalid_request', 'request.patterns[0]: expected an object'),
            # several names: the base refused every such list as not one name; with a name the
            # server lacks it still is (a list of known names runs on /resolve since release 10)
            ('POST', 'resolve', dict(r, graphs=['kleb', 'nope']), 'graph_unknown', one_name),
        ]
        return panel

    @staticmethod
    def _ask(server, method, route, request):
        return server.get(route + request) if method == 'GET' else server.post(route, request)

    def test_refusal_texts_are_those_before_release_10(self):
        """Every request the release-9 server refused is refused in its words, byte for byte:
        the probes' and /traverse's and /resolve's selection texts, /search's and /pattern's
        unknown name, the code of SPEC-multigraph-envelope.md section 6 added (a malformed
        selector keeps no code, /pattern's every 400 its invalid_request). Pinned here as
        literals; the same requests against the release-9 binary:
        test_refusal_texts_equal_the_base_binarys."""
        for method, route, request, code, text in self._refusal_panel():
            with self.subTest(route=route, request=request):
                ret = self._ask(self.server, method, route, request)
                self.assertEqual(400, ret.status_code, ret.text[:2000])
                body = ret.json()
                self.assertEqual(text, body['error'])
                self.assertEqual(code, body.get('code'))

    def _names_cap_server(self, name, binary=None):
        return Server([self.csv, '--mmap', '--mem-cap-gb', '1', '--max-graphs-without-selection',
                       '1'], os.path.join(self.tempdir.name, f'server_cap_{name}.log'),
                      binary=binary)

    def _names_cap_panel(self):
        """The requests without a selection beyond the names cap (one name, four listed), each
        with its refusal body, byte for byte the release-9 one where that server answered the
        request: (method, route, payload or query, body, as the base)"""
        kleb = self.records['kleb'][0]
        pattern_text = ('request.graphs: required on this server, which hosts 4 graph names '
                        '(more than 1; GET /capabilities lists them)')
        panel = [
            ('POST', 'search', {'FASTA': '>q\n' + kleb[3000:3100]},
             {'error': 'Bad request: requests without names (no "graphs" field) are only '
                       'supported for small indexes (<=1 names)'}, True),
            ('POST', 'pattern', {'patterns': [{'dna': kleb[3000:3020]}], 'mode': 'count'},
             {'code': 'invalid_request', 'error': pattern_text}, True),
            # new in release 10 (the release-9 /resolve needed a name): /pattern's text
            ('POST', 'resolve', {'sequence': kleb[3000:3100]}, {'error': pattern_text}, False),
            ('POST', 'traverse', {'seeds': [{'sequence': kleb[3000:3100]}]},
             {'error': "Bad request: 'graph' (index name) is required in multi-graph mode"},
             True),
        ]
        for route in ('traverse/capabilities', 'pattern/capabilities'):
            panel.append(('GET', route, '',
                          {'error': f'Bad request: in multi-graph mode GET /{route} needs '
                                    f'?graph=<name> (and graph_path=<path> when the name spans '
                                    f'several graphs); GET /capabilities lists the graphs'},
                          True))
        return panel

    def test_the_names_cap_texts_are_those_before_release_10(self):
        """Beyond --max-graphs-without-selection a request without a selection is refused in
        the words of the release-9 server, byte for byte: /search's, /pattern's (which
        /resolve uses since it reads the cap, release 10), the probes' without ?graph=, and
        /traverse's (it always needed a name). SPEC-multigraph-envelope.md section 6."""
        server = self._names_cap_server('new')
        try:
            for method, route, request, body, _ in self._names_cap_panel():
                with self.subTest(route=route):
                    ret = self._ask(server, method, route, request)
                    self.assertEqual((400, body), (ret.status_code, ret.json()))
        finally:
            server.stop()

    @unittest.skipUnless(BASE_BINARY and os.path.isfile(BASE_BINARY),
                         '$METAGRAPH_BASE_BINARY (the release-9 binary, 8fa32e553) is not set')
    def test_refusal_texts_equal_the_base_binarys(self):
        """The panels above against the release-9 binary on the same lists: the same status
        and error text for every request it answered (the code is release 10's addition)."""
        base = Server([self.csv, '--mmap', '--mem-cap-gb', '1'],
                      os.path.join(self.tempdir.name, 'server_base.log'), binary=BASE_BINARY)
        try:
            for method, route, request, code, text in self._refusal_panel():
                with self.subTest(route=route, request=request):
                    a = self._ask(base, method, route, request)
                    b = self._ask(self.server, method, route, request)
                    self.assertEqual((400, text), (a.status_code, a.json()['error']))
                    self.assertEqual((a.status_code, a.json()['error']),
                                     (b.status_code, b.json()['error']))
        finally:
            base.stop()
        base, new = self._names_cap_server('base', BASE_BINARY), self._names_cap_server('new2')
        try:
            for method, route, request, body, as_base in self._names_cap_panel():
                if not as_base:
                    continue
                with self.subTest(route=route):
                    a = self._ask(base, method, route, request)
                    b = self._ask(new, method, route, request)
                    self.assertEqual((400, body), (a.status_code, a.json()))
                    self.assertEqual((a.status_code, a.content), (b.status_code, b.content))
        finally:
            base.stop()
            new.stop()

    # The differences to a release-10 base binary that E23 (the merge of master 33ad588d2)
    # intends, each pinned by test_align_expired_and_404_against_the_base_binary
    INTENDED_DIFFERENCES = {
        # master's #722: /align refuses an invalid max_alternative_alignments (release 10 fell
        # back to 1, or failed on a value that is no integer)
        'max_alternative_alignments':
            {'error': "'max_alternative_alignments' must be a positive integer."},
        # SPEC-multigraph-envelope.md section 4: not_after_ms is checked at admission
        'expired_409': ("when the request's handler started", "at the request's admission"),
        # SPEC-multigraph-envelope.md section 5: the three cases of GET /requests' 404
        # (F3b: the never-admitted case names the instance_mismatch 409 too)
        'requests_404': ("never received, or no longer retained",
                         "not admitted yet, never admitted (an HTTP-level 400, the expired 409 "
                         "or the instance_mismatch 409; or never received), or ended and no "
                         "longer retained"),
        # F3, section 5: the 404 states the server's clock read in its lookup
        'requests_404_member': 'server_time_ms',
    }

    @unittest.skipUnless(BASE_BINARY and os.path.isfile(BASE_BINARY),
                         '$METAGRAPH_BASE_BINARY (a release-10 binary, 5c7b45b74) is not set')
    def test_align_expired_and_404_against_the_base_binary(self):
        """/align (single-graph and multi-graph), the expired 409 of /resolve, /pattern and
        /traverse, and GET /requests' 404 against a release-10 base binary on the same files:
        byte for byte the same answers but for INTENDED_DIFFERENCES, each checked as stated
        there (a release-9 base, which has no request registry, skips this test)."""
        d = self.tempdir.name
        graph, anno = self.pairs['kleb']
        seed = self.records['kleb'][0][1000:1100]
        fasta = '>q\n' + seed
        base_multi = Server([self.csv, '--mmap', '--mem-cap-gb', '1'],
                            os.path.join(d, 'server_base_r10.log'), binary=BASE_BINARY)
        servers = [base_multi]
        try:
            if 'request_id' not in base_multi.get('capabilities').json().get('multigraph', {}):
                self.skipTest('the base binary has no request registry (before release 10)')
            base_single = Server(['-i', graph, '-a', anno],
                                 os.path.join(d, 'server_base_r10_single.log'), binary=BASE_BINARY)
            servers.append(base_single)
            new_single = Server(['-i', graph, '-a', anno],
                                os.path.join(d, 'server_new_single.log'))
            servers.append(new_single)

            def raw(server, route, data):
                return requests.post(server.url(route), data=data, timeout=300)

            # /align on a single-graph server: (what, body, the base's status)
            same = [
                ('bad json', '{"FASTA": ">q\\nACGT",', 400),
                ('min_exact_match above 1', json.dumps({'FASTA': fasta, 'min_exact_match': 1.5}),
                 400),
                ('min_exact_match below 0', json.dumps({'FASTA': fasta, 'min_exact_match': -0.5}),
                 400),
                # no FASTA: no query, an empty answer
                ('no FASTA', json.dumps({'min_exact_match': 0.5}), 200),
            ]
            for what, data, status in same:
                with self.subTest(align=what):
                    a, b = raw(base_single, 'align', data), raw(new_single, 'align', data)
                    self.assertEqual(status, a.status_code, a.text[:500])
                    self.assertEqual((a.status_code, a.content), (b.status_code, b.content))
            refusal = self.INTENDED_DIFFERENCES['max_alternative_alignments']
            for value in (-1, 0, None, True, 1.5, '2', 2147483648):
                with self.subTest(max_alternative_alignments=value):
                    data = json.dumps({'FASTA': fasta, 'max_alternative_alignments': value})
                    a, b = raw(base_single, 'align', data), raw(new_single, 'align', data)
                    self.assertEqual((400, refusal), (b.status_code, b.json()))
                    self.assertNotEqual((400, refusal), (a.status_code, a.json()))
            data = json.dumps({'FASTA': fasta, 'max_alternative_alignments': 2})
            self.assertEqual((200, 200), (raw(base_single, 'align', data).status_code,
                                          raw(new_single, 'align', data).status_code))

            # /align on a multi-graph server: the refusal in release 10's words
            data = json.dumps({'FASTA': fasta})
            a, b = raw(base_multi, 'align', data), raw(self.server, 'align', data)
            self.assertEqual(400, a.status_code, a.text[:500])
            self.assertEqual((a.status_code, a.content), (b.status_code, b.content))

            # the expired 409: the same body but its wording and the server's clock
            old, new = self.INTENDED_DIFFERENCES['expired_409']
            clock = re.compile(r"server's clock \(\d+, [^)]*\)")
            for route, body in (
                    ('resolve', {'sequence': seed, 'graphs': ['kleb'],
                                 'discover': {'max_labels': 2}, 'not_after_ms': 1,
                                 'request_id': 'it-expired-base-resolve'}),
                    ('resolve', {'sequence': seed, 'graph': 'kleb',
                                 'discover': {'max_labels': 2}, 'not_after_ms': 1}),
                    ('pattern', {'patterns': [{'dna': seed[:20]}], 'mode': 'count',
                                 'graphs': ['kleb', 'pseudo'], 'not_after_ms': 1,
                                 'request_id': 'it-expired-base-pattern'}),
                    ('traverse', {'seeds': [{'sequence': seed}], 'graph': 'kleb',
                                  'not_after_ms': 1, 'attempt_id': 'it-expired-base-traverse'}),
                    ('traverse', {'seeds': [{'sequence': seed}], 'graph': 'kleb',
                                  'not_after_ms': 1})):
                with self.subTest(expired=route, body=sorted(body)):
                    a, b = base_multi.post(route, body), self.server.post(route, body)
                    self.assertEqual((409, 409), (a.status_code, b.status_code), b.text[:500])
                    a, b = a.json(), b.json()
                    self.assertEqual('expired', b['state'])
                    self.assertIn(old, a['error'])
                    self.assertNotIn(old, b['error'])
                    self.assertIn(new, b['error'])
                    self.assertEqual(clock.sub('', a['error'].replace(old, new)),
                                     clock.sub('', b['error']))
                    self.assertEqual(set(a), set(b))
                    varying = ('error', 'server_time_ms', 'server_instance')
                    self.assertEqual({k: v for k, v in a.items() if k not in varying},
                                     {k: v for k, v in b.items() if k not in varying})

            # GET /requests' 404: the same members and server_time_ms, the text naming the
            # three cases
            old, new = self.INTENDED_DIFFERENCES['requests_404']
            a = base_multi.get('requests/it-never-sent')
            b = self.server.get('requests/it-never-sent')
            self.assertEqual((404, 404), (a.status_code, b.status_code))
            a, b = a.json(), b.json()
            added = self.INTENDED_DIFFERENCES['requests_404_member']
            self.assertNotIn(added, a)
            self.assertEqual(set(a) | {added}, set(b))
            self.assertIs(type(b[added]), int)
            self.assertEqual(a['request_id'], b['request_id'])
            self.assertIn(old, a['error'])
            self.assertEqual(a['error'].replace(old, new), b['error'])
        finally:
            for server in servers:
                server.stop()

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
            caps = single(server.get('pattern/capabilities?graph=kleb').json())
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
                # the selection's code (SPEC-multigraph-envelope.md section 6)
                (dict(p, graphs=['nope']), 400, 'graph_unknown', "unknown graph 'nope'"),
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
        self.assertEqual({'error': 'Request with an uninitialized graph 1', 'code': 'graph_unknown'},
                         self.post('search', dict(s, graphs=[1]), 400))
        # an unknown name: 400 on both, each in its own words (/search's with the selector
        # refusal's code, SPEC-multigraph-envelope.md section 6)
        self.assertEqual({'error': 'Request with an uninitialized graph nope', 'code': 'graph_unknown'},
                         self.post('search', dict(s, graphs=['nope']), 400))
        self.assertEqual({'code': 'graph_unknown',
                          'error': "request.graphs: unknown graph 'nope' (GET /capabilities "
                                   "lists the graphs)"},
                         self.post('pattern', dict(p, graphs=['nope']), 400))

    def test_search_tags_rows_and_states_a_cut_per_pair(self):
        """SPEC-multigraph-envelope.md section 3 (E1): on a multi-graph server every /search
        row carries its pair (graph, graph_path, annotation_path); the rows of a query sequence
        come in the index's pair order (the names in byte order, each name's pairs in the
        list's order); each query sequence states one cut per pair searched, in that order,
        empty pairs included, index_fp null for a pair without a manifest. The chunks' labels
        are disjoint, so a query cut from one chunk's record leaves the other pairs empty."""
        chunks = [name for name, _, _, _ in CHUNKS]
        cut = {c: self.records[c][0][4000:4100] for c in chunks}
        s = {'FASTA': ''.join(f'>{c}\n{cut[c]}\n' for c in chunks), 'discovery_fraction': 0.5}
        pair_tags = ('graph', 'graph_path', 'annotation_path')

        def tags(name, chunk):
            return dict(zip(pair_tags, (name,) + self.pairs[chunk]))

        # selected in any order: ecoli < enterics < kleb < pseudo, enterics' pairs as listed
        out = self.post('search', dict(s, graphs=['pseudo', 'enterics', 'ecoli']))
        order = [('ecoli', 'ecoli'), ('enterics', 'ecoli'), ('enterics', 'kleb'),
                 ('pseudo', 'pseudo')]
        self.assertEqual(chunks, [query['seq_description'] for query in out])
        empty = 0
        for query in out:
            chunk = query['seq_description']
            with self.subTest(query=chunk):
                rows = query['results']
                # each pair's rows together, the pairs in the index's order
                seen, returned = [], {}
                for row in rows:
                    pair = next(p for p in order if tags(*p) == {t: row[t] for t in pair_tags})
                    if not seen or seen[-1] != pair:
                        seen.append(pair)
                    returned[pair] = returned.get(pair, 0) + 1
                self.assertEqual(sorted(set(seen), key=order.index), seen)
                # the pairs over the query's chunk answer it (the chunks share sequence, so
                # others may as well)
                self.assertTrue({p for p in order if p[1] == chunk} <= set(seen), seen)
                # one cut per pair, the empty ones included
                self.assertEqual([dict(tags(*pair), index_fp=None, returned=returned.get(pair, 0),
                                       complete=True, nth_score=None, tied=0, tied_returned=0)
                                  for pair in order],
                                 query['cuts'])
                empty += sum(1 for cut in query['cuts'] if cut['returned'] == 0)
                # a pair under two names answers alike under each
                untagged = lambda name: [{k: v for k, v in r.items() if k not in pair_tags}
                                         for r in rows if r['graph'] == name
                                         and r['graph_path'] == self.pairs['ecoli'][0]]
                self.assertEqual(untagged('ecoli'), untagged('enterics'))
        self.assertTrue(empty)

        # the cut: the two labels of pseudo (its FASTA files) on a query of 100 bp of each,
        # the windows chosen so that both match as many k-mers (a basic graph: the k-mers as
        # written); the top 1 is the label first in byte order
        files = [os.path.join(FASTA, f) for f in CHUNKS[2][1]]
        longest = [sorted(records_of([f]), key=len, reverse=True)[0] for f in files]
        texts = ['\n'.join(records_of([f])) for f in files]

        def score(query, text):
            return sum(query[i:i + K] in text for i in range(len(query) - K + 1))

        for offset in range(4000, 40000, 500):
            query = ''.join(seq[offset:offset + 100] for seq in longest)
            scores = [score(query, text) for text in texts]
            if scores[0] == scores[1] > 0:
                break
        else:
            self.fail('no window of the two files matches as many k-mers of each')
        q = {'FASTA': f'>q\n{query}\n', 'discovery_fraction': 0.1, 'graphs': ['pseudo']}
        [complete] = self.post('search', dict(q, top_labels=10))
        self.assertEqual([{'sample': f, 'kmer_count': scores[0], **tags('pseudo', 'pseudo')}
                          for f in sorted(files)],
                         complete['results'])
        [top] = self.post('search', dict(q, top_labels=1))
        self.assertEqual(complete['results'][:1], top['results'])
        self.assertEqual([dict(tags('pseudo', 'pseudo'), index_fp=None, returned=1,
                               complete=False, nth_score=scores[0], tied=2, tied_returned=1)],
                         top['cuts'])

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

    def test_a_pairs_400_in_the_envelope_keeps_the_request_id(self):
        """Only an HTTP-level 400 frees a request_id (section 5). A pair's 400 found by the work
        is a refused entry of a 200 envelope: the request was answered, GET /requests/{id}
        states it ended with that pair refused, and the same id is a 409 (a corrected request
        takes a new id). A mixed fan-out is never forgotten for one pair's 400 (GPT review 3,
        finding 4)."""
        ecoli = self.records['ecoli'][0]
        request = {'patterns': [{'dna': ecoli[123456:123476], 'id': 'p'}],
                   'mode': 'all_or_count', 'output': {'labels': 'all'}}
        for request_id, graphs in (('it-pair400-1', ['kleb']),
                                   ('it-pair400-2', ['kleb', 'ecoli'])):
            out = self.post('pattern', dict(request, graphs=graphs, request_id=request_id))
            self.assertEqual(request_id, out['request_id'])
            refused = [a for a in out['answers'] if a['graph'] == 'kleb']
            self.assertEqual([(400, 'annotation_unbudgeted')],
                             [(a['refusal']['http_status'], a['refusal']['code'])
                              for a in refused])
            entry = self._wait_ended(self.server, request_id)
            self.assertEqual(('/pattern', 'ended', 'answered'),
                             (entry['route'], entry['state'], entry['reason']))
            self.assertEqual({('kleb', 'refused')} | ({('ecoli', 'answered')}
                                                      if 'ecoli' in graphs else set()),
                             {(p['graph'], p['outcome']) for p in entry['pairs']})
            repeat = self.server.post('pattern', dict(request, graphs=['ecoli'],
                                                      request_id=request_id))
            self.assertEqual(409, repeat.status_code, repeat.text[:2000])
            self.assertEqual(request_id, repeat.json()['request_id'])

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

    def _cutoff_ms(self, request):
        """The work cutoff (ms after the request's arrival) of a cutoff test with |request|:
        three times the time a one-pair /resolve of it on kleb takes on this machine now (the
        class's server, a_kleb's files), at least 750 ms (an idle machine's). a_kleb completes
        by then also under load (the tests failed at a machine load of 30 and more with a_kleb
        unfinished at a fixed 750 ms), and b_big, whose work on the request is over four times
        a_kleb's, starting after it on the one pool thread, is still working: the margins scale
        with the load on both sides"""
        started = time.time()
        self.post('resolve', dict(request, graphs=['kleb']))
        return max(750, int(3000 * (time.time() - started)))

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
        effective deadline (arrival + the server's maximum, --traverse-max-time-ms, for a
        /resolve without a budget or finish_by_ms) less the finalisation reserve (250 ms): 1000
        ms on an idle machine, scaled with the load (_cutoff_ms). On
        one pool thread a_kleb completes; b_big starts after it and runs without a deadline (a
        /resolve without one has no budget for the cutoff to cap): it is still working at the
        cutoff, over a second from its end, so it is unfinished, with its usage so far, and
        told to stop; c_tiny never started: not_reached. The completed answer is written beside
        them, by the effective deadline, as a single-graph server answers it. The request's
        request_id (section 5): GET /requests/{request_id} states it ended once the cut worker
        ended, each pair with its state, outcome and final usage."""
        request = self._slow_resolve(3000000)
        cutoff_ms = self._cutoff_ms(request)
        server = self._cutoff_server('cutoff', '--traverse-max-time-ms', str(cutoff_ms + 250))
        try:
            started = time.time()
            ret = server.post('resolve', dict(request, request_id='it-cutoff-1'))
            took = time.time() - started
            self.assertEqual(200, ret.status_code, ret.text[:2000])
            out = ret.json()
            self.assertEqual((['a_kleb', 'b_big', 'c_tiny'], 1, 2),
                             (out['graphs'], out['answered'], out['refused']),
                             [(a['graph'], a['outcome'], a.get('refusal', {}).get('code'))
                              for a in out['answers']])
            kleb, big, tiny = out['answers']
            self.assertEqual(('a_kleb', 'answered'), (kleb['graph'], kleb['outcome']))
            # the single-graph answer on the pair's files (no deadline: no limits, no stop),
            # with the unknown_labels of a multi-graph server's pair
            self.assertEqual(untimed(dict(self.cli_resolve('kleb', request), unknown_labels=[])),
                             untimed({k: v for k, v in kleb.items() if k not in TAGS}))
            self.assertNotIn('limits', kleb)
            self.assertEqual(('b_big', 'refused'), (big['graph'], big['outcome']))
            self.assertEqual({'http_status': 503, 'code': 'unfinished', 'reason': 'cutoff'},
                             {k: big['refusal'][k] for k in ('http_status', 'code', 'reason')})
            self.assertIn("had not finished at the request's work cutoff",
                          big['refusal']['error'])
            self.assertIn('GET /requests/it-cutoff-1', big['refusal']['error'])
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
            # the registry: the request ends when the cut worker ends
            entry = self._wait_ended(server, 'it-cutoff-1')
            self.assertEqual(('ended', 'answered'), (entry['state'], entry['reason']))
            p_kleb, p_big, p_tiny = entry['pairs']
            self.assertEqual([('a_kleb', 'ended', 'answered'), ('b_big', 'ended', 'unfinished'),
                              ('c_tiny', 'never_started', 'not_reached')],
                             [(p['graph'], p['state'], p['outcome']) for p in entry['pairs']])
            self.assertEqual({'started_at', 'finished_at', 'elapsed_ms', 'cpu_ms',
                              'memory_bytes'}, set(p_big['usage']))
            self.assertIsNotNone(p_big['usage']['finished_at'])
            # the cut worker's final usage: at least what the envelope stated
            self.assertGreaterEqual(p_big['usage']['elapsed_ms'], usage['elapsed_ms'])
            self.assertIsNone(p_tiny['usage'])
        finally:
            server.stop()

    def test_a_pooled_pair_stopped_by_time_answers_by_the_cutoff(self):
        """SPEC-multigraph-envelope.md section 4 (the owner's decision of 2026-11-01): a pooled
        pair's own budget is capped at the request's work cutoff (arrival + time_budget_ms 1500
        less the finalisation reserve 500 ms), so that a pair stopped by time answers by the
        cutoff instead of racing it. On one pool thread a_kleb completes, its single-graph
        answer exactly (its own budget stated, nothing clamped); b_big, whose work runs for over
        a second, starts after it, stops by time before the cutoff and answers what it found,
        the cap stated as a clamp of its budget (limits.time_budget_ms, the time_budget_ms entry
        of limits.clamped: the time left before the cutoff at its start)."""
        server = self._cutoff_server('cap')
        try:
            request = dict(self.SLOW_PATTERN, time_budget_ms=1500)
            started = time.time()
            out = self.post_to(server, 'pattern', dict(request, graphs=['a_kleb', 'b_big']))
            took = time.time() - started
            self.assertEqual([('a_kleb', 'answered'), ('b_big', 'answered')],
                             [(a['graph'], a['outcome']) for a in out['answers']])
            self.assertEqual(([], []), (out['not_reached'], out['unfinished']))
            kleb, big = out['answers']
            self.assertEqual((1500, []), (kleb['limits']['time_budget_ms'],
                                          kleb['limits']['clamped']))
            # the pair alone (one pair: inline, the single-graph request) answers the same
            alone = self.post_to(server, 'pattern', dict(request, graphs=['a_kleb']))
            self.assertEqual(untimed(alone['answers'][0]), untimed(kleb))
            # b_big: stopped by time, its budget capped at the cutoff and stated so
            budget = big['limits']['time_budget_ms']
            self.assertEqual([{'field': 'time_budget_ms', 'requested': 1500,
                               'effective': budget}], big['limits']['clamped'])
            self.assertLess(budget, 1000)
            self.assertIn('time', {e['stop']['reason'] for e in big['patterns'] if e['stop']})
            self.assertEqual(16, len(big['patterns']))
            # by the effective deadline
            self.assertLess(took, 5)
        finally:
            server.stop()

    def test_sigterm_during_a_multi_pair_resolve_writes_no_complete_answer(self):
        """SIGTERM while a two-pair /resolve runs (class Shutdown in server.cpp, master's #681:
        the signal thread stops the HTTP workers only, the acceptor and the connections are
        closed after their join): the server exits with status 0 within its bounds (the drain,
        then kExitMs for a handler that does not return), and the client receives no
        complete-looking answer, only a closed connection or a body short of its
        Content-Length. b_big's resolve of the whole record runs for seconds; the signal comes
        once the request is logged as started."""
        server = Server([self._cutoff_csv(), '--traverse-max-time-ms', '600000'],
                        os.path.join(self.tempdir.name, 'server_sigterm_resolve.log'), threads=2)
        got = {}

        def post():
            try:
                got['response'] = server.post('resolve', dict(self._slow_resolve(100000000),
                                                              graphs=['a_kleb', 'b_big'],
                                                              request_id='it-sigterm-1'))
            except requests.exceptions.RequestException as e:
                got['error'] = e
        client = threading.Thread(target=post)
        try:
            client.start()
            deadline = time.time() + 60
            while 'request_id it-sigterm-1' not in server.text():
                self.assertLess(time.time(), deadline, server.text()[-3000:])
                self.assertTrue(client.is_alive(), got)
                time.sleep(0.05)
            time.sleep(0.3)
            self.assertTrue(client.is_alive(), 'the request ended before the signal')
            started = time.time()
            server.process.send_signal(signal.SIGTERM)
            self.assertEqual(0, server.process.wait(timeout=60))
            took = time.time() - started
            client.join(timeout=60)
            self.assertFalse(client.is_alive())
            # the drain (3 s, traversals only) and kExitMs (3 s) bound the exit; slack for load
            self.assertLess(took, 20)
            self.assertNotIn('response', got,
                             got.get('response') is not None and
                             (got['response'].status_code, got['response'].text[:300]))
            self.assertIsInstance(got['error'], (requests.exceptions.ConnectionError,
                                                 requests.exceptions.ChunkedEncodingError))
            text = server.text()
            self.assertIn('[Server] Shutting down...', text)
            self.assertIn('Signal 15: shutting down', text)
        finally:
            server.stop()
            client.join(timeout=60)

    def test_finish_by_ms_cuts_a_pair_and_the_registry_states_its_end(self):
        """finish_by_ms (E20) on /resolve without a budget: the effective deadline is
        finish_by_ms less the clock skew allowance (--traverse-clock-skew-ms 500), the cutoff
        250 ms (the reserve) before it (finish_by_ms 1600 ms after the request on an idle
        machine, scaled with the load: _cutoff_ms); the deadline block states all four, and a pair's
        deadline is the minimum of its own (the server's maximum), finish_by_ms - the skew and,
        pooled, the cutoff. On one pool thread a_kleb, from which the query is absent, completes
        first, under the finish_by_ms clamp alone, as a single-graph request with that
        finish_by_ms; b_big resolves the query (about a second of work) after it: it stops by
        time before the cutoff and answers the resolve of a prefix, both caps stated as clamps of
        its budget (finish_by_ms's, then the cutoff's, the budget it ran under). The request's
        request_id (E4): GET /requests/{request_id} states it ended, each pair with its state,
        outcome and final usage, the response written (status, bytes); a repeat of the id
        starts no work (409, the entry with an error); an id never sent, or past the retention
        count (--request-retention 1), is 404."""
        server = self._cutoff_server('finish', '--traverse-clock-skew-ms', '500',
                                     '--request-retention', '1')
        try:
            # effective = finish_by - 500 (skew), the cutoff 250 ms (reserve) before it
            request = self._slow_resolve()
            window_ms = self._cutoff_ms(request) + 750
            finish_by = int(time.time() * 1000) + window_ms
            body = dict(request, graphs=['a_kleb', 'b_big'], finish_by_ms=finish_by,
                        request_id='it-finish-1')
            started = time.time()
            ret = server.post('resolve', body)
            took = time.time() - started
            self.assertEqual(200, ret.status_code, ret.text[:2000])
            out = ret.json()
            # by the effective deadline (the writer refuses a later envelope, 503), not after
            # b_big's whole work
            self.assertLess(took, (window_ms - 500) / 1000 + 1.9)
            self.assertEqual(('/resolve', 'it-finish-1'), (out['route'], out['request_id']))
            self.assertEqual([('a_kleb', 'answered'), ('b_big', 'answered')],
                             [(a['graph'], a['outcome']) for a in out['answers']])
            self.assertEqual(([], []), (out['not_reached'], out['unfinished']))
            kleb, big = out['answers']
            # a_kleb: finish_by_ms - the skew caps its deadline (none requested: the server's
            # maximum), stated as a clamp; it completed (stop null)
            clamped = kleb['limits']['clamped']
            self.assertEqual(1, len(clamped), clamped)
            self.assertEqual(('bounds.time_budget_ms', None),
                             (clamped[0]['field'], clamped[0]['requested']))
            self.assertLess(clamped[0]['effective'], window_ms - 500)
            self.assertEqual(clamped[0]['effective'], kleb['limits']['time_budget_ms'])
            self.assertIsNone(kleb['stop'])
            # b_big: stopped by time before the cutoff, the resolve of a prefix; both caps
            # stated, the cutoff's last (the budget it ran under)
            self.assertEqual('time', big['stop']['reason'])
            self.assertLess(big['stop']['resolved_kmers'], big['stop']['query_kmers'])
            clamped = big['limits']['clamped']
            self.assertEqual(['bounds.time_budget_ms'] * 2, [c['field'] for c in clamped])
            self.assertEqual([None, None], [c['requested'] for c in clamped])
            self.assertLess(clamped[1]['effective'], clamped[0]['effective'])
            self.assertLess(clamped[1]['effective'], window_ms - 750)
            self.assertEqual(clamped[1]['effective'], big['limits']['time_budget_ms'])
            iso = lambda ms: time.strftime('%Y-%m-%dT%H:%M:%S', time.gmtime(ms // 1000)) \
                + f'.{ms % 1000:03d}Z'
            self.assertEqual({'not_after_ms': None, 'finish_by_ms': finish_by,
                              'effective_deadline_at': iso(finish_by - 500),
                              'cutoff_at': iso(finish_by - 750)}, out['deadline'])
            # the registry: the request ended once its workers ended
            entry = self._wait_ended(server, 'it-finish-1')
            self.assertEqual(('it-finish-1', '/resolve', 'ended', 'answered'),
                             (entry['request_id'], entry['route'], entry['state'],
                              entry['reason']))
            self.assertEqual(out['server_instance'], entry['server_instance'])
            # the bytes written: the body as sent (compressed, the client accepting gzip)
            self.assertEqual({'status': 200, 'bytes': int(ret.headers['Content-Length']),
                              'written': True}, entry['response'])
            self.assertIsNotNone(entry['retained_until'])
            self.assertEqual([('a_kleb', 'ended', 'answered'), ('b_big', 'ended', 'answered')],
                             [(p['graph'], p['state'], p['outcome']) for p in entry['pairs']])
            for p in entry['pairs']:
                self.assertEqual({'started_at', 'finished_at', 'elapsed_ms', 'cpu_ms',
                                  'memory_bytes'}, set(p['usage']))
                self.assertIsNotNone(p['usage']['finished_at'])
            # a repeat starts no work: 409, the entry and an error
            ret = server.post('resolve', body)
            self.assertEqual(409, ret.status_code, ret.text[:2000])
            repeat = ret.json()
            self.assertIn("request_id 'it-finish-1' is running, or retained", repeat['error'])
            self.assertEqual(entry, {k: v for k, v in repeat.items() if k != 'error'})
            self.assertEqual(1, server.text().count('graph name(s), request_id it-finish-1'))
            # an id never sent (F3: with the server's clock, an integer)
            before = int(time.time() * 1000)
            ret = server.get('requests/it-never-sent')
            self.assertEqual(404, ret.status_code, ret.text[:2000])
            self.assertEqual({'error', 'request_id', 'server_instance', 'server_time_ms'},
                             set(ret.json()))
            self.assertIs(type(ret.json()['server_time_ms']), int)
            self.assertGreaterEqual(ret.json()['server_time_ms'], before)
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
        time_budget_ms (1500 ms) under a finish_by_ms a minute later puts the cutoff at arrival
        + 1000 ms, which caps b_big's budget (stopped by time, answered, the cap stated), the
        completed a_kleb written beside it within the budget"""
        server = self._cutoff_server('short')
        try:
            finish_by = int(time.time() * 1000) + 60000
            started = time.time()
            out = self.post_to(server, 'pattern', dict(self.SLOW_PATTERN, time_budget_ms=1500,
                                                       finish_by_ms=finish_by,
                                                       graphs=['a_kleb', 'b_big']))
            self.assertLess(time.time() - started, 5)
            self.assertEqual([('a_kleb', 'answered'), ('b_big', 'answered')],
                             [(a['graph'], a.get('refusal', {}).get('code', a['outcome']))
                              for a in out['answers']])
            kleb, big = out['answers']
            # finish_by_ms a minute later caps nothing: a_kleb states its own budget
            self.assertEqual((1500, []), (kleb['limits']['time_budget_ms'],
                                          kleb['limits']['clamped']))
            # the cutoff (arrival + 1000 ms) caps b_big's
            self.assertEqual(['time_budget_ms'], [c['field'] for c in big['limits']['clamped']])
            self.assertLess(big['limits']['time_budget_ms'], 1000)
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

    def test_resolve_and_pattern_caps_refuse_busy_and_register_nothing(self):
        """E19 on a multi-graph server (SPEC-multigraph-envelope.md section 1.1, section 8;
        SPEC-labeled-traversal-core.md §10.3, "Handler capacity"): -p 1, --http-threads 4,
        --resolve-max-concurrent 1 and --pattern-max-concurrent 1. A /resolve of b_big's whole
        record (seconds) holds the /resolve slot, the slow pattern on b_big the /pattern slot.
        Meanwhile a second /pattern and a second /resolve, each with a request_id, are refused
        at once: 503, Retry-After 1, {code: "busy", error}, unwrapped (no envelope, no
        request_id), their ids not registered (GET /requests/{id} 404); GET /capabilities (its
        multigraph block stating the caps), GET /traverse/capabilities (each probe the same
        object), GET /pattern/capabilities and GET /requests/{id} of the running requests
        answer within a second. Once each slot is free, the refused request sent again after
        Retry-After runs, its request_id admitted."""
        server = Server([self._cutoff_csv(), '--http-threads', '4', '--resolve-max-concurrent',
                         '1', '--pattern-max-concurrent', '1', '--traverse-max-time-ms', '600000',
                         '--pattern-max-contexts', '1000000'],
                        os.path.join(self.tempdir.name, 'server_busy.log'), threads=1)
        got = {}

        def post(key, route, body):
            try:
                got[key] = server.post(route, body)
            except requests.exceptions.RequestException as e:
                got[key] = e

        def running(request_id, deadline):
            while True:
                ret = server.get(f'requests/{request_id}')
                if ret.status_code == 200 and ret.json()['state'] == 'running':
                    return
                self.assertLess(time.time(), deadline, (request_id, ret.text[:2000]))
                time.sleep(0.02)
        slow_resolve = dict(self._slow_resolve(100000000), graphs=['b_big'],
                            request_id='it-busy-r1')
        slow_pattern = dict(self.SLOW_PATTERN, time_budget_ms=60000, graphs=['b_big'],
                            request_id='it-busy-p1')
        second_resolve = dict(self._slow_resolve(1000), graphs=['a_kleb'],
                              request_id='it-busy-r2')
        second_pattern = dict(self.SLOW_PATTERN, graphs=['a_kleb'], request_id='it-busy-p2')
        clients = [threading.Thread(target=post, args=('r1', 'resolve', slow_resolve)),
                   threading.Thread(target=post, args=('p1', 'pattern', slow_pattern))]
        took = {}

        def timed(name, call, status):
            started = time.time()
            ret = call()
            took[name] = time.time() - started
            self.assertEqual(status, ret.status_code, (name, ret.text[:2000]))
            return ret

        def busy(ret, route):
            self.assertEqual('1', ret.headers.get('Retry-After'))
            self.assertEqual(['code', 'error', 'server_instance'], sorted(ret.json()))
            self.assertEqual('busy', ret.json()['code'])
            self.assertEqual(server.get('capabilities').json()['server_instance'],
                             ret.json()['server_instance'])
            self.assertIn(f'at most 1 {route} request(s) at once (--{route[1:]}-max-concurrent)',
                          ret.json()['error'])
        try:
            clients[0].start()
            running('it-busy-r1', time.time() + 60)
            clients[1].start()
            running('it-busy-p1', time.time() + 60)
            # both slots held: the second of each route is refused, nothing of it registered
            p2 = timed('POST /pattern (busy)', lambda: server.post('pattern', second_pattern),
                       503)
            busy(p2, '/pattern')
            r2 = timed('POST /resolve (busy)', lambda: server.post('resolve', second_resolve),
                       503)
            busy(r2, '/resolve')
            for request_id in ('it-busy-p2', 'it-busy-r2'):
                self.assertEqual(404, server.get(f'requests/{request_id}').status_code)
            caps = timed('GET /capabilities', lambda: server.get('capabilities'), 200).json()
            stated = {'/pattern': 1, '/resolve': 1, '/traverse': 0}
            self.assertEqual(stated, caps['multigraph']['max_concurrent'])
            probes = timed('GET /traverse/capabilities',
                           lambda: server.get('traverse/capabilities'), 200).json()
            self.assertEqual([stated] * 3, [e['max_concurrent'] for e in probes['answers']])
            timed('GET /pattern/capabilities', lambda: server.get('pattern/capabilities'), 200)
            # the slow resolve (seconds) still holds its slot; the slow pattern (over a second)
            # held its own when the second /pattern was refused
            for request_id in ('it-busy-r1', 'it-busy-p1'):
                entry = timed(f'GET /requests/{request_id}',
                              lambda: server.get(f'requests/{request_id}'), 200).json()
                if request_id == 'it-busy-r1':
                    self.assertEqual('running', entry['state'], (request_id, took))
            for name, seconds in took.items():
                self.assertLess(seconds, 1.0, (name, took))
            # each slot freed: the refused request, sent again after Retry-After, runs
            for key, client, route, body, refused in ((
                    'p1', clients[1], 'pattern', second_pattern, p2), (
                    'r1', clients[0], 'resolve', second_resolve, r2)):
                client.join(timeout=300)
                self.assertEqual(200, got[key].status_code, got[key])
                time.sleep(int(refused.headers['Retry-After']))
                ret = server.post(route, body)
                self.assertEqual(200, ret.status_code, ret.text[:2000])
                self.assertEqual(body['request_id'], ret.json()['request_id'])
                self.assertEqual('ended', self._wait_ended(server, body['request_id'])['state'])
            text = server.text()
            self.assertIn('Handler threads: 4 (--http-threads); requests run at once at most: '
                          '/traverse 0, /resolve 1, /pattern 1', text)
            self.assertIn(': busy, 1 of 1 running (--pattern-max-concurrent); 503', text)
            self.assertIn(': busy, 1 of 1 running (--resolve-max-concurrent); 503', text)
            print('\n[E19 multi-graph] beside full /resolve and /pattern caps (s): '
                  + ', '.join(f'{k} {v:.3f}' for k, v in took.items()), flush=True)
        finally:
            server.stop()
            for client in clients:
                if client.ident is not None:
                    client.join(timeout=60)

    def test_single_graph_resolve_and_pattern_caps_refuse_busy(self):
        """E19 on a single-graph server (the mini index, -i/-a): --http-threads 4,
        --resolve-max-concurrent 1 and --pattern-max-concurrent 1. While a /resolve of a whole
        record (seconds) and the slow pattern (over a second) hold the slots, a second /resolve and
        a second /pattern, each with a request_id, are refused as on a multi-graph server: 503,
        Retry-After 1, {code: "busy", error, server_instance}, their ids not registered."""
        mini = (os.path.join(MINI_DIR, 'graph_k31.dbg'),
                os.path.join(MINI_DIR, 'annotation.relaxed.relabeled.row_diff_brwt_coord.annodbg'))
        if not all(os.path.exists(p) for p in mini):
            if REQUIRE_GUARDS:
                self.fail(f'$METAGRAPH_REQUIRE_GUARDS=1: the mini index is missing: {mini}')
            self.skipTest(f'the mini index is not built: {mini}')
        server = Server(['-i', mini[0], '-a', mini[1], '--http-threads', '4',
                         '--resolve-max-concurrent', '1', '--pattern-max-concurrent', '1',
                         '--traverse-max-time-ms', '600000', '--pattern-max-contexts', '1000000'],
                        os.path.join(self.tempdir.name, 'server_busy_single.log'), threads=1)
        got = {}

        def post(key, route, body):
            try:
                got[key] = server.post(route, body)
            except requests.exceptions.RequestException as e:
                got[key] = e

        def running(request_id):
            deadline = time.time() + 60
            while True:
                ret = server.get(f'requests/{request_id}')
                if ret.status_code == 200 and ret.json()['state'] == 'running':
                    return
                self.assertLess(time.time(), deadline, (request_id, ret.text[:2000]))
                time.sleep(0.02)
        clients = [threading.Thread(target=post, args=(
                       'r1', 'resolve', dict(self._slow_resolve(100000000), request_id='it-sb-r1'))),
                   threading.Thread(target=post, args=(
                       'p1', 'pattern', dict(self.SLOW_PATTERN, time_budget_ms=60000,
                                             request_id='it-sb-p1')))]
        try:
            instance = server.get('capabilities').json()['server_instance']
            for client, request_id, route, body in (
                    (clients[0], 'it-sb-r1', 'resolve',
                     dict(self._slow_resolve(1000), request_id='it-sb-r2')),
                    (clients[1], 'it-sb-p1', 'pattern',
                     dict(self.SLOW_PATTERN, request_id='it-sb-p2'))):
                client.start()
                running(request_id)
                ret = server.post(route, body)
                self.assertEqual(503, ret.status_code, ret.text[:2000])
                self.assertEqual('1', ret.headers.get('Retry-After'))
                self.assertEqual({'code': 'busy', 'server_instance': instance},
                                 {k: v for k, v in ret.json().items() if k != 'error'})
                self.assertIn(f'(--{route}-max-concurrent)', ret.json()['error'])
                self.assertEqual(404, server.get('requests/' + body['request_id']).status_code)
            for client in clients:
                client.join(timeout=300)
            self.assertEqual([200, 200], [got['r1'].status_code, got['p1'].status_code], got)
        finally:
            server.stop()
            for client in clients:
                if client.ident is not None:
                    client.join(timeout=60)

    def test_two_long_searches_leave_the_light_routes_a_thread(self):
        """E19, the plan's finding 10 (SPEC-labeled-traversal-core.md §10.3, "Handler
        capacity"): -p sets the graphs' pool and, without --http-threads, the handler threads
        too, and a /search handler waits for its graphs on that pool, so with -p 1 one long
        search holds the one handler thread and a second waits behind it: a capabilities call
        and a /resolve sent meanwhile wait until both are done. With --http-threads 4 (-p 1:
        the same pool) the two searches still run one after the other on the pool, and
        meanwhile GET /capabilities, GET /requests/{id} and a named /resolve answer within a
        second. /search has no cap: both searches are answered."""
        query = self.records['pseudo'][0][:3000000]
        search = {'FASTA': '>q\n' + query, 'graphs': ['b_big'], 'top_labels': 5}
        resolve = dict(self._slow_resolve(1000), graphs=['a_kleb'], request_id='it-search-r1')

        def run(server, flags):
            got = {}

            def post(key, body):
                try:
                    got[key] = server.post('search', body)
                except requests.exceptions.RequestException as e:
                    got[key] = e
            # how long one search takes here (a second distinct one: no cache shared)
            started = time.time()
            first = server.post('search', dict(search, FASTA='>w\n' + query[1:]))
            alone = time.time() - started
            self.assertEqual(200, first.status_code, first.text[:2000])
            clients = [threading.Thread(target=post, args=(i, dict(search, FASTA=f'>q{i}\n'
                                                                    + query[i + 2:])))
                       for i in range(2)]
            for c in clients:
                c.start()
            # both handlers are in: the first searching, the second waiting for the pool
            time.sleep(min(0.3, alone / 4))
            took = {}
            for name, call in (('GET /capabilities', lambda: server.get('capabilities')),
                               ('POST /resolve', lambda: server.post('resolve', resolve)),
                               ('GET /requests/{request_id}',
                                lambda: server.get('requests/' + resolve['request_id']))):
                started = time.time()
                ret = call()
                took[name] = time.time() - started
                self.assertEqual(200, ret.status_code, (flags, name, ret.text[:2000]))
            busy = [c.is_alive() for c in clients]
            for c in clients:
                c.join(timeout=600)
            for i in range(2):
                self.assertEqual(200, got[i].status_code, got[i])
            return alone, took, busy

        log = os.path.join(self.tempdir.name, 'server_searches_{}.log')
        server = Server([self._cutoff_csv(), '--cache-size', '0'], log.format('p1'), threads=1)
        try:
            alone, coupled, _ = run(server, '-p 1')
        finally:
            server.stop()
        if alone < 1.0:
            self.skipTest(f'one search takes {alone:.2f} s here: too fast to hold a thread')
        # without --http-threads the capabilities call waited for both searches
        self.assertGreater(coupled['GET /capabilities'], alone / 2, coupled)
        resolve['request_id'] = 'it-search-r2'
        server = Server([self._cutoff_csv(), '--cache-size', '0', '--http-threads', '4'],
                        log.format('http4'), threads=1)
        try:
            alone, took, busy = run(server, '-p 1 --http-threads 4')
        finally:
            server.stop()
        self.assertEqual([True, True], busy, (alone, took))
        for name, seconds in took.items():
            self.assertLess(seconds, 1.0, (name, took, alone))
        print(f'\n[E19 searches] one search alone {alone:.2f} s; beside two: '
              + ', '.join(f'{k} {v:.3f}' for k, v in took.items())
              + '; with -p 1 alone: ' + ', '.join(f'{k} {v:.3f}' for k, v in coupled.items()),
              flush=True)

    def test_not_after_ms_admits_resolve_and_pattern(self):
        """not_after_ms (the latest START, E20) on /resolve and /pattern: passed at the
        request's admission, 409 {error, state: expired, not_after_ms, server_time_ms} and nothing
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

    def test_a_404_past_not_after_ms_proves_a_lost_request_is_never_admitted(self):
        """F3 (SPEC-multigraph-envelope.md section 5): a request the service lost (never
        received here) is looked up by its id once its not_after_ms has passed: GET /requests
        answers 404 with server_time_ms > not_after_ms, and the request sent after that, with
        the same id, is refused 409 expired at a server time at least as late (nothing
        registered: 404 again). The same holds for a /pattern request_id and a /traverse
        attempt_id. With not_after_ms ahead, the 404 proves nothing (server_time_ms <=
        not_after_ms) and the request then runs."""
        kleb = self.records['kleb'][0]
        seed = kleb[3000:3100]
        bodies = {
            'resolve': {'sequence': seed, 'graph': 'kleb', 'discover': {'max_labels': 2}},
            'pattern': {'patterns': [{'dna': seed[:20]}], 'mode': 'count',
                        'graphs': ['kleb', 'pseudo']},
            'traverse': {'seeds': [{'sequence': seed}], 'graph': 'kleb',
                         'strategy': {'direction': 'right', 'bounds': {'max_extension_bp': 50},
                                      'output': {'detail': 'summary'}}},
        }
        handle = {'resolve': 'request_id', 'pattern': 'request_id', 'traverse': 'attempt_id'}
        instance = self.server.get('capabilities').json()['attempts']['server_instance']
        for route, body in bodies.items():
            with self.subTest(route=route):
                lost = f'it-lost-{route}'
                not_after = int(time.time() * 1000) - 1000
                ret = self.server.get(f'requests/{lost}')
                self.assertEqual(404, ret.status_code, ret.text[:2000])
                gone = ret.json()
                self.assertEqual((lost, instance),
                                 (gone['request_id'], gone['server_instance']))
                self.assertIs(type(gone['server_time_ms']), int)
                self.assertGreater(gone['server_time_ms'], not_after)
                # the request arrives after all (a late copy, or the re-send): refused
                ret = self.server.post(route, dict(body, not_after_ms=not_after,
                                                   **{handle[route]: lost}))
                self.assertEqual(409, ret.status_code, ret.text[:2000])
                out = ret.json()
                self.assertEqual(('expired', not_after), (out['state'], out['not_after_ms']))
                self.assertGreaterEqual(out['server_time_ms'], gone['server_time_ms'])
                ret = self.server.get(f'requests/{lost}')
                self.assertEqual(404, ret.status_code, ret.text[:2000])
                self.assertGreaterEqual(ret.json()['server_time_ms'], out['server_time_ms'])
                # not_after_ms ahead: no proof, and the request runs
                ahead = f'it-ahead-{route}'
                not_after = int(time.time() * 1000) + 60000
                ret = self.server.get(f'requests/{ahead}')
                self.assertEqual(404, ret.status_code, ret.text[:2000])
                self.assertLessEqual(ret.json()['server_time_ms'], not_after)
                ret = self.server.post(route, dict(body, not_after_ms=not_after,
                                                   **{handle[route]: ahead}))
                self.assertEqual(200, ret.status_code, ret.text[:2000])

    def test_a_single_graph_404_proves_a_lost_request_is_never_admitted(self):
        """The same on a single-graph server (one registry per process as well): GET
        /requests answers 404 {error, request_id, server_instance, server_time_ms} with
        server_time_ms > the lost request's not_after_ms, and the copy sent after it is refused
        409 expired at a server time at least as late, on /resolve, /pattern and /traverse."""
        kleb = self.records['kleb'][0]
        seed = kleb[3000:3100]
        bodies = {
            'resolve': {'sequence': seed[:80], 'discover': {'max_labels': 2}},
            'pattern': {'patterns': [{'dna': seed[:20]}], 'mode': 'count'},
            'traverse': {'seeds': [{'sequence': seed}],
                         'strategy': {'direction': 'right', 'bounds': {'max_extension_bp': 50},
                                      'output': {'detail': 'summary'}}},
        }
        handle = {'resolve': 'request_id', 'pattern': 'request_id', 'traverse': 'attempt_id'}
        server = self._single_server('lost')
        try:
            instance = server.get('traverse/capabilities').json()['attempts']['server_instance']
            for route, body in bodies.items():
                with self.subTest(route=route):
                    lost = f'sg-lost-{route}'
                    not_after = int(time.time() * 1000) - 1000
                    ret = server.get(f'requests/{lost}')
                    self.assertEqual(404, ret.status_code, ret.text[:2000])
                    gone = ret.json()
                    self.assertEqual({'error', 'request_id', 'server_instance', 'server_time_ms'},
                                     set(gone))
                    self.assertEqual((lost, instance),
                                     (gone['request_id'], gone['server_instance']))
                    self.assertIs(type(gone['server_time_ms']), int)
                    self.assertGreater(gone['server_time_ms'], not_after)
                    ret = server.post(route, dict(body, not_after_ms=not_after,
                                                  **{handle[route]: lost}))
                    self.assertEqual(409, ret.status_code, ret.text[:2000])
                    out = ret.json()
                    self.assertEqual(('expired', not_after), (out['state'], out['not_after_ms']))
                    self.assertGreaterEqual(out['server_time_ms'], gone['server_time_ms'])
                    self.assertEqual(404, server.get(f'requests/{lost}').status_code)
        finally:
            server.stop()

    def _expect_server_instance_panel(self, server, bodies, multi):
        """F3b on one server: each route of |bodies| (/resolve, /pattern) with request_id and
        expect_server_instance (SPEC-multigraph-envelope.md section 5)"""
        caps = server.get('capabilities').json()
        instance = caps['server_instance']
        routes = ['/traverse', '/resolve', '/pattern']
        # the capabilities key: in the multigraph block, or at the top level of a single-graph
        # server beside routes.requests (one reader for both modes)
        if multi:
            self.assertEqual(routes, caps['multigraph']['expect_server_instance'])
            self.assertNotIn('expect_server_instance', caps)
        else:
            self.assertNotIn('multigraph', caps)
            self.assertEqual(routes, caps['expect_server_instance'])
            self.assertEqual('GET /requests/{request_id}', caps['routes']['requests'])
        other = '0000000000000000'
        self.assertNotEqual(other, instance)
        for route, body in bodies.items():
            with self.subTest(route=route, multi=multi):
                rid = f'f3b-{"multi" if multi else "single"}-{route}'
                not_after = int(time.time() * 1000) + 120000
                request = dict(body, request_id=rid, not_after_ms=not_after,
                               expect_server_instance=other)
                # another instance: 409 before the admission, request-level (no envelope)
                ret = server.post(route, request)
                self.assertEqual(409, ret.status_code, ret.text[:2000])
                out = ret.json()
                self.assertEqual({'error', 'state', 'expect_server_instance', 'server_instance',
                                  'request_id', 'not_after_ms'}, set(out))
                self.assertEqual(('instance_mismatch', other, instance, rid, not_after),
                                 (out['state'], out['expect_server_instance'],
                                  out['server_instance'], out['request_id'],
                                  out['not_after_ms']))
                self.assertEqual(f"expect_server_instance '{other}' is not this server's instance "
                                 f"('{instance}'; the process restarted, or the request reached "
                                 "another server): it was not started", out['error'])
                # nothing registered: never admitted
                ret = server.get(f'requests/{rid}')
                self.assertEqual(404, ret.status_code, ret.text[:2000])
                self.assertIn('never admitted (an HTTP-level 400, the expired 409 or the '
                              'instance_mismatch 409; or never received)', ret.json()['error'])
                self.assertIs(type(ret.json()['server_time_ms']), int)
                # this server's instance: runs, and the entry ends on it
                ret = server.post(route, dict(request, expect_server_instance=instance))
                self.assertEqual(200, ret.status_code, ret.text[:2000])
                if multi:
                    self.assertEqual((rid, instance), (ret.json()['request_id'],
                                                       ret.json()['server_instance']))
                for _ in range(100):
                    entry = server.get(f'requests/{rid}')
                    if entry.status_code == 200 and entry.json()['state'] == 'ended':
                        break
                    time.sleep(0.05)
                self.assertEqual(200, entry.status_code, entry.text[:2000])
                self.assertEqual(('ended', instance, 'answered'),
                                 (entry.json()['state'], entry.json()['server_instance'],
                                  entry.json()['reason']))
                # ahead of the repeated id (retained) and of an expired not_after_ms
                ret = server.post(route, request)
                self.assertEqual((409, 'instance_mismatch'),
                                 (ret.status_code, ret.json().get('state')), ret.text[:2000])
                ret = server.post(route, dict(request, request_id=rid + '-late', not_after_ms=1))
                self.assertEqual((409, 'instance_mismatch'),
                                 (ret.status_code, ret.json().get('state')), ret.text[:2000])
                self.assertEqual(404, server.get(f'requests/{rid}-late').status_code)
                # the field's 400s: without request_id, or not a token
                for bad in ({'expect_server_instance': instance},
                            dict(request, expect_server_instance=''),
                            dict(request, expect_server_instance=None),
                            dict(request, expect_server_instance='a b'),
                            # one of the request's own 400s: ahead of an expired not_after_ms
                            dict(request, expect_server_instance='a b', not_after_ms=1)):
                    ret = server.post(route, dict(body, **bad))
                    self.assertEqual(400, ret.status_code, ret.text[:2000])
                    self.assertIn('request.expect_server_instance', ret.json()['error'])

    def test_expect_server_instance_pins_resolve_and_pattern(self):
        """F3b (SPEC-multigraph-envelope.md sections 5 and 8), on a multi-graph and a
        single-graph server: a /resolve or /pattern with request_id and expect_server_instance
        naming another instance is refused 409 instance_mismatch {error, state,
        expect_server_instance, server_instance, request_id, not_after_ms} before its admission
        (no pair starts, no envelope, GET /requests answers 404), ahead of a repeated id and of
        an expired not_after_ms; pinned to this server's instance it runs (200, the entry ended
        on it); without request_id, or malformed, a 400 naming the field. GET /capabilities
        lists the routes (the multigraph block; a single-graph server's top level)."""
        kleb = self.records['kleb'][0]
        seed = kleb[3000:3100]
        self._expect_server_instance_panel(self.server, {
            'resolve': {'sequence': seed, 'graph': 'kleb', 'discover': {'max_labels': 2}},
            'pattern': {'patterns': [{'dna': seed[:20]}], 'mode': 'count',
                        'graphs': ['kleb', 'pseudo']},
        }, True)
        server = self._single_server('f3b')
        try:
            self._expect_server_instance_panel(server, {
                'resolve': {'sequence': seed[:80], 'discover': {'max_labels': 2}},
                'pattern': {'patterns': [{'dna': seed[:20]}], 'mode': 'count'},
            }, False)
        finally:
            server.stop()

    def test_the_cli_refuses_expect_server_instance(self):
        """F3b: the CLI's /resolve (traverse --resolve) and /pattern, whose instance is their own
        and random, refuse any expect_server_instance given with request_id: the 409 body on
        stdout, exit status 1, as the CLI's /traverse; without the field the request runs"""
        graph, anno = self.pairs['kleb']
        seed = self.records['kleb'][0][3000:3100]
        cases = {
            ('traverse', '--resolve'): {'sequence': seed[:80], 'discover': {'max_labels': 2}},
            ('pattern',): {'patterns': [{'dna': seed[:20]}], 'mode': 'count'},
        }
        for command, body in cases.items():
            with self.subTest(command=command):
                path = os.path.join(self.tempdir.name, f'f3b-cli-{command[0]}.json')
                for pinned in (True, False):
                    request = dict(body, request_id='f3b-cli')
                    if pinned:
                        request['expect_server_instance'] = '0000000000000000'
                    with open(path, 'w') as f:
                        json.dump(request, f)
                    res = subprocess.run(shlex.split(METAGRAPH) + list(command)
                                         + ['--json', '-i', graph, '-a', anno, path],
                                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    if not pinned:
                        self.assertEqual(0, res.returncode, res.stderr.decode()[-2000:])
                        continue
                    self.assertEqual(1, res.returncode, res.stdout.decode()[:2000])
                    out = json.loads(res.stdout.decode())
                    self.assertEqual({'error', 'state', 'expect_server_instance',
                                      'server_instance', 'request_id'}, set(out))
                    self.assertEqual(('instance_mismatch', '0000000000000000', 'f3b-cli'),
                                     (out['state'], out['expect_server_instance'],
                                      out['request_id']))
                    self.assertNotEqual('0000000000000000', out['server_instance'])

    def _single_server(self, name, *args, binary=None):
        """A single-graph server on the kleb pair"""
        graph, anno = self.pairs['kleb']
        return Server(['-i', graph, '-a', anno] + list(args),
                      os.path.join(self.tempdir.name, f'server_single_{name}.log'), binary=binary)

    def test_single_graph_resolve_400s_register_no_request_id(self):
        """A single-graph /resolve with request_id refused 400 for its own fields (a budget
        not above the finalisation reserve, the query longer than bounds.max_query_bp or the
        server's limit, a label the index does not hold, a query shorter than k) registers
        nothing (SPEC-multigraph-envelope.md section 5, notice B4): GET /requests answers 404
        and the corrected request with the same id runs, as on a multi-graph server."""
        server = self._single_server('ids', '--resolve-max-query-bp', '5000')
        try:
            kleb = self.records['kleb'][0]
            ok = {'sequence': kleb[3000:3080], 'discover': {'max_labels': 2}}
            for n, (bad, text) in enumerate((
                    (dict(ok, bounds={'time_budget_ms': 100}),
                     'request.bounds.time_budget_ms: expected more than the finalisation '
                     'reserve of 250 ms'),
                    (dict(ok, sequence=kleb[3000:3100], bounds={'max_query_bp': 10}),
                     'request.sequence: longer than bounds.max_query_bp'),
                    (dict(ok, sequence=kleb[:6000], bounds={'max_query_bp': 10000}),
                     'request.sequence: longer than the server limit of 5000 bp'),
                    ({'sequence': ok['sequence'], 'labels': ['no-such-label']},
                     "Label not found in the annotation: 'no-such-label'"),
                    (dict(ok, sequence=kleb[:20]), 'Query shorter than k'))):
                request_id = f'sg-400-{n}'
                with self.subTest(error=text):
                    ret = server.post('resolve', dict(bad, request_id=request_id))
                    self.assertEqual(400, ret.status_code, ret.text[:2000])
                    self.assertTrue(ret.json()['error'].startswith(text), ret.json())
                    ret = server.get(f'requests/{request_id}')
                    self.assertEqual(404, ret.status_code, ret.text[:2000])
                    ret = server.post('resolve', dict(ok, request_id=request_id))
                    self.assertEqual(200, ret.status_code, ret.text[:2000])
                    self.assertEqual(('ended', 'answered'),
                                     tuple(self._wait_ended(server, request_id)[k]
                                           for k in ('state', 'reason')))
            # the same refusals with an expired not_after_ms: 400, not the expired 409
            past = int(time.time() * 1000) - 1000
            for bad in (dict(ok, bounds={'max_query_bp': 10}), {'labels': ['x']},
                        {'sequence': 5}):
                with self.subTest(expired=bad):
                    for extra in ({}, {'request_id': 'sg-expired'}):
                        ret = server.post('resolve', dict(bad, not_after_ms=past, **extra))
                        self.assertEqual(400, ret.status_code, ret.text[:2000])
                    self.assertEqual(404, server.get('requests/sg-expired').status_code)
        finally:
            server.stop()

    def test_single_graph_resolve_seed_400s_hold_no_request_id(self):
        """A 400 never holds a request_id (SPEC-multigraph-envelope.md section 5, notice B4). A
        single-graph /resolve's explicit seeds faulty by the request and k alone are refused
        before its admission, so also when its not_after_ms has passed; a seed only the work
        finds faulty is refused after it and frees the id when it answers: GET /requests
        answers 404 and the corrected request with the same id runs. The texts, with two
        faults too, are the base binary's (METAGRAPH_BASE_BINARY)."""
        server = self._single_server('seeds')
        base = None
        try:
            kleb = self.records['kleb'][0]
            label = server.post('resolve', {'sequence': kleb[3000:3080],
                                            'discover': {'max_labels': 1}}).json()
            label = label['labels'][0]['label']
            # 40 bases of the record, then bases no graph holds: a query of 110 k-mers
            absent = 'ACGT' * 25
            query = kleb[3000:3040] + absent
            kmers = len(query) - K + 1
            def explicit(*seeds):
                return {'sequence': query, 'labels': [label],
                        'select': {'policy': 'explicit',
                                   'seeds': [{'kmer_interval': list(iv), 'labels': labels}
                                             for iv, labels in seeds]}}
            good = explicit(((0, 5), [label]))
            cases = (
                ('range', explicit(((0, 1000), [label])), 'Explicit seed interval out of range'),
                ('empty', explicit(((5, 5), [label])), 'Explicit seed interval out of range'),
                ('edge', explicit(((0, kmers + 1), [label])),
                 'Explicit seed interval out of range'),
                ('label', explicit(((0, 5), ['nolabel'])),
                 "Explicit seed label was not profiled: 'nolabel'"),
                ('second', explicit(((0, 5), [label]), ((5, 5), [label])),
                 'Explicit seed interval out of range'),
                # two faults: the order of the selection, seed by seed, the interval, its
                # presence in the graph, its labels
                ('two-a', explicit(((5, 5), ['nolabel'])), 'Explicit seed interval out of range'),
                ('two-b', explicit(((20, 40), ['nolabel'])),
                 'Explicit seed interval is not fully present in the graph'),
                ('two-c', explicit(((20, 40), [label]), ((5, 5), [label])),
                 'Explicit seed interval is not fully present in the graph'),
                ('two-d', explicit(((0, 5), ['nolabel']), ((20, 40), [label])),
                 "Explicit seed label was not profiled: 'nolabel'"),
            )
            for name, bad, text in cases:
                request_id = f'seed-{name}'
                with self.subTest(case=name):
                    ret = server.post('resolve', dict(bad, request_id=request_id))
                    self.assertEqual(400, ret.status_code, ret.text[:2000])
                    self.assertEqual(text, ret.json()['error'])
                    self.assertEqual(404, server.get(f'requests/{request_id}').status_code)
                    ret = server.post('resolve', dict(good, request_id=request_id))
                    self.assertEqual(200, ret.status_code, ret.text[:2000])
            # found by the work, after the admission: the id is free once it is answered
            bad = explicit(((20, 40), [label]))
            for n in range(2):
                ret = server.post('resolve', dict(bad, request_id='seed-work'))
                self.assertEqual(400, ret.status_code, ret.text[:2000])
                self.assertEqual('Explicit seed interval is not fully present in the graph',
                                 ret.json()['error'])
                self.assertEqual(404, server.get('requests/seed-work').status_code)
            ret = server.post('resolve', dict(good, request_id='seed-work'))
            self.assertEqual(200, ret.status_code, ret.text[:2000])
            self.assertEqual('answered', self._wait_ended(server, 'seed-work')['reason'])
            # a discovery's seed label is known after the work
            ret = server.post('resolve', {'sequence': query, 'discover': {'max_labels': 2},
                                          'request_id': 'seed-discover',
                                          'select': {'policy': 'explicit', 'seeds': [
                                              {'kmer_interval': [0, 5], 'labels': ['nolabel']}]}})
            self.assertEqual(400, ret.status_code, ret.text[:2000])
            self.assertEqual(404, server.get('requests/seed-discover').status_code)
            # and it is checked before the seeds after it: their faults are not stated first
            # (recheck4: a later seed's interval was stated where the base states the label)
            discover_two = {'sequence': query, 'discover': {'max_labels': 2},
                            'select': {'policy': 'explicit', 'seeds': [
                                {'kmer_interval': [0, 5], 'labels': ['nolabel']},
                                {'kmer_interval': [0, 1000], 'labels': []}]}}
            not_profiled = "Explicit seed label was not profiled: 'nolabel'"
            ret = server.post('resolve', dict(discover_two, request_id='seed-discover-two'))
            self.assertEqual(400, ret.status_code, ret.text[:2000])
            self.assertEqual(not_profiled, ret.json()['error'])
            self.assertEqual(404, server.get('requests/seed-discover-two').status_code)
            # before the admission: a 400, not the expired 409
            past = int(time.time() * 1000) - 1000
            for name, bad, text in cases:
                with self.subTest(expired=name):
                    ret = server.post('resolve', dict(bad, not_after_ms=past))
                    self.assertEqual(400, ret.status_code, ret.text[:2000])
                    self.assertEqual(text, ret.json()['error'])
            # the texts of the base binary (which knows no request_id)
            if BASE_BINARY:
                base = self._single_server('seeds-base', binary=BASE_BINARY)
                for name, bad, text in cases:
                    with self.subTest(base=name):
                        ret = base.post('resolve', bad)
                        self.assertEqual(400, ret.status_code, ret.text[:2000])
                        self.assertEqual(text, ret.json()['error'])
                ret = base.post('resolve', discover_two)
                self.assertEqual(400, ret.status_code, ret.text[:2000])
                self.assertEqual(not_profiled, ret.json()['error'])
        finally:
            server.stop()
            if base:
                base.stop()

    def test_single_graph_pattern_400s_come_before_the_expired_409(self):
        """A single-graph /pattern refused for its own fields is refused 400 (invalid_request)
        also when its not_after_ms has passed, with or without a request_id, and registers
        nothing (SPEC-multigraph-envelope.md section 5, notice B4: the request's own 400s come
        before the expired 409); a valid request with an expired not_after_ms gets the 409."""
        server = self._single_server('pattern')
        try:
            kleb = self.records['kleb'][0]
            past = int(time.time() * 1000) - 1000
            for patterns, text in (('bad', 'request.patterns: expected a list of 1 to '),
                                   ([5], 'request.patterns[0]: expected an object')):
                for extra in ({}, {'request_id': 'sg-pattern'}):
                    with self.subTest(patterns=patterns, **extra):
                        ret = server.post('pattern', dict({'patterns': patterns,
                                                           'not_after_ms': past}, **extra))
                        self.assertEqual(400, ret.status_code, ret.text[:2000])
                        self.assertEqual('invalid_request', ret.json()['code'])
                        self.assertTrue(ret.json()['error'].startswith(text), ret.json())
                        self.assertEqual(404, server.get('requests/sg-pattern').status_code)
            ret = server.post('pattern', {'patterns': [{'dna': kleb[3000:3020]}],
                                          'not_after_ms': past, 'request_id': 'sg-pattern'})
            self.assertEqual(409, ret.status_code, ret.text[:2000])
            self.assertEqual('expired', ret.json()['state'])
            # a 400 only the index decides (no budget-aware decode on this annotation), after
            # the admission: it frees the id when it is answered, as any 400
            unbudgeted = {'patterns': [{'dna': kleb[3000:3020]}], 'mode': 'all_or_count',
                          'output': {'labels': 'all'}, 'request_id': 'sg-unbudgeted'}
            ret = server.post('pattern', unbudgeted)
            self.assertEqual((400, 'annotation_unbudgeted'),
                             (ret.status_code, ret.json()['code']), ret.text[:2000])
            self.assertEqual(404, server.get('requests/sg-unbudgeted').status_code)
            ret = server.post('pattern', {'patterns': [{'dna': kleb[3000:3020]}],
                                          'request_id': 'sg-unbudgeted'})
            self.assertEqual(200, ret.status_code, ret.text[:2000])
            self.assertEqual('answered', self._wait_ended(server, 'sg-unbudgeted')['reason'])
        finally:
            server.stop()

    def test_an_admitted_request_runs_past_not_after_ms(self):
        """not_after_ms is admission only (E20): a request admitted before it runs past it and
        finishes -- b_big's resolve of 3 Mbp, over a second (1.25 s on an idle machine),
        admitted 300 ms before not_after_ms, alone (inline) and beside c_tiny; c_tiny, queued
        behind b_big on the one pool thread, was not started by not_after_ms: not_reached,
        reason not_after_ms (decided when not_after_ms passed, not at the cutoff)."""
        server = self._cutoff_server('admission')
        try:
            request = self._slow_resolve(3000000)
            for graphs in (['b_big'], ['b_big', 'c_tiny']):
                with self.subTest(graphs=graphs):
                    not_after = int(time.time() * 1000) + 300
                    out = self.post_to(server, 'resolve', dict(request, graphs=graphs,
                                                               not_after_ms=not_after))
                    # (it ran on well past not_after_ms; the margin holds on an idle machine,
                    # where the whole resolve takes about 1.25 s)
                    self.assertGreater(time.time() * 1000, not_after + 500)
                    big = out['answers'][0]
                    self.assertEqual(('b_big', 'answered'), (big['graph'], big['outcome']))
                    self.assertGreater(big['timing']['elapsed_ms'], 800)
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

    def test_request_ids_and_attempt_ids_are_one_namespace(self):
        """GET /requests/{id} reads the request registry and the attempt ledger, so an id is
        admitted by one of them only while the other holds nothing with it (GPT A3): a
        /traverse attempt_id equal to a retained request_id is refused 409 {error, attempt}
        and a /resolve or /pattern request_id equal to a retained attempt_id 409 with the
        attempt's request object plus error, in both orders, nothing run; GET answers the id's
        holder; a repeat on one route answers as before; past the retention the id runs on
        the other route."""
        server = Server([self.csv, '--request-retention-s', '2',
                         '--traverse-attempt-retention-s', '2'],
                        os.path.join(self.tempdir.name, 'server_ids.log'))
        try:
            ecoli = self.records['ecoli'][0]
            seed = ecoli[20000:20100]
            t = {'graphs': ['ecoli'], 'seeds': [{'sequence': seed}],
                 'strategy': {'direction': 'right', 'bounds': {'max_extension_bp': 50},
                              'output': {'detail': 'summary'}}}
            r = {'graphs': ['ecoli'], 'sequence': seed, 'discover': {'max_labels': 2}}
            p = {'graphs': ['ecoli'], 'patterns': [{'dna': seed[:20]}], 'mode': 'count'}

            # a retained /resolve, then the attempt
            self.post_to(server, 'resolve', dict(r, request_id='one'))
            self._wait_ended(server, 'one')
            ret = server.post('traverse', dict(t, attempt_id='one'))
            self.assertEqual(409, ret.status_code, ret.text[:2000])
            self.assertEqual({'error', 'attempt'}, set(ret.json()))
            self.assertIn("attempt_id 'one' is running, or retained", ret.json()['error'])
            self.assertEqual(('one', 'finished', '/resolve'),
                             tuple(ret.json()['attempt'][k]
                                   for k in ('attempt_id', 'state', 'route')))
            self.assertEqual('/resolve', server.get('requests/one').json()['route'])
            self.assertEqual(404, server.get('traverse/attempt/one').status_code)

            # a retained attempt, then /resolve and /pattern
            self.post_to(server, 'traverse', dict(t, attempt_id='two'))
            self._wait_ended(server, 'two')
            for route, request in (('resolve', r), ('pattern', p)):
                with self.subTest(route=route):
                    ret = server.post(route, dict(request, request_id='two'))
                    self.assertEqual(409, ret.status_code, ret.text[:2000])
                    body = ret.json()
                    self.assertIn("request_id 'two' is running, or retained", body['error'])
                    self.assertEqual({'error'} | set(server.get('requests/two').json()),
                                     set(body))
                    self.assertEqual(('two', '/traverse', 'ended'),
                                     tuple(body[k] for k in ('request_id', 'route', 'state')))
            self.assertEqual('/traverse', server.get('requests/two').json()['route'])

            # a repeat on one route, as before
            ret = server.post('resolve', dict(r, request_id='one'))
            self.assertEqual(409, ret.status_code, ret.text[:2000])
            self.assertEqual(('/resolve', 'ended'), (ret.json()['route'], ret.json()['state']))
            self.assertNotIn('attempt', ret.json())
            ret = server.post('traverse', dict(t, attempt_id='two'))
            self.assertEqual(409, ret.status_code, ret.text[:2000])
            self.assertEqual({'error', 'attempt'}, set(ret.json()))
            self.assertNotIn('route', ret.json()['attempt'])

            # past the retention: each id runs on the other route
            time.sleep(2.5)
            self.assertEqual(404, server.get('requests/one').status_code)
            self.assertEqual(404, server.get('requests/two').status_code)
            self.post_to(server, 'traverse', dict(t, attempt_id='one'))
            self.post_to(server, 'resolve', dict(r, request_id='two'))
            self.assertEqual('/traverse', self._wait_ended(server, 'one')['route'])
            self.assertEqual('/resolve', self._wait_ended(server, 'two')['route'])
        finally:
            server.stop()

    # -------------------------------------------------------------------- per-task figures

    def test_traverse_pair_tags_live_as_long_as_the_attempt(self):
        """GET /requests/{attempt_id} names the pair a /traverse attempt ran on for as long as
        the attempt is answered (fanout F5 of the release-10 review): with the request
        registry's retention at 0 the tags are still there after the attempt finished."""
        server = Server([self.csv, '--request-retention-s', '0'],
                        os.path.join(self.tempdir.name, 'server_retention0.log'))
        try:
            seed = self.records['ecoli'][0][20000:20100]
            t = {'graphs': ['ecoli'], 'seeds': [{'sequence': seed}], 'attempt_id': 'it-ret-1',
                 'strategy': {'direction': 'right', 'bounds': {'max_extension_bp': 50},
                              'output': {'detail': 'summary'}}}
            self.post_to(server, 'traverse', t)
            entry = self._wait_ended(server, 'it-ret-1')
            (pair,) = entry['pairs']
            self.assertEqual(('ecoli', self.pairs['ecoli'][0], self.pairs['ecoli'][1]),
                             (pair['graph'], pair['graph_path'], pair['annotation_path']))
            self.assertEqual('answered', pair['outcome'])
        finally:
            server.stop()

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
        resident = single(self.post('traverse', t), '/traverse')
        loaded = single(self.post('traverse', dict(t, in_ram=True, attempt_id='in-ram-1')))
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

    def test_in_ram_is_not_loaded_past_finish_by_ms(self):
        """SPEC-multigraph-envelope.md section 4 (fanout F2 of the release-10 review): the
        in_ram memory wait is not entered once finish_by_ms has passed, also by a request of
        one pair (inline) whose memory is free: nothing is loaded, and the pair is not reached,
        on /resolve and /pattern alike."""
        kleb = self.records['kleb'][0]
        past = int(time.time() * 1000) - 60000
        for route, body in (
                ('resolve', {'graphs': ['kleb'], 'sequence': kleb[3000:3100],
                             'discover': {'max_labels': 2}}),
                ('pattern', {'graphs': ['kleb'], 'patterns': [{'dna': kleb[3000:3020]}],
                             'mode': 'count'})):
            with self.subTest(route=route):
                loads = self.server.text().count('loaded into RAM for it')
                out = self.post(route, dict(body, in_ram=True, finish_by_ms=past))
                self.assertEqual((0, 1), (out['answered'], out['refused']))
                (entry,) = out['answers']
                self.assertEqual({'http_status': 503, 'code': 'not_reached', 'reason': 'cutoff'},
                                 {k: entry['refusal'][k]
                                  for k in ('http_status', 'code', 'reason')})
                self.assertEqual([self._ref(entry)], out['not_reached'])
                self.assertEqual(loads, self.server.text().count('loaded into RAM for it'))

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
        resident = single(self.post('traverse', t))
        loaded = single(self.post('traverse', dict(t, in_ram=True)))
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
                        resident, loaded = single(resident)['results'], single(loaded)['results']
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
            out = single(via_graphs.json(), '/traverse')
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
            # graph, the older spelling: the same answer, the spelling not echoed (the request
            # was a 400 on its first absent seed before round D); the envelopes differ only in
            # their time
            via_graph = self.server.post('traverse', dict(request, graph=chunk))
            self.assertEqual(200, via_graph.status_code, via_graph.text[:2000])
            self.assertEqual(untimed(via_graphs.json()), untimed(via_graph.json()), chunk)
        # every case was met: walked, partly present, absent
        self.assertEqual({'walked', 'partly', 'absent'}, seen)
        # the seed of a chunk walks there as in a request that names the graph with graph
        one = dict(request, seeds=[{'seed_id': 'e', 'sequence': seeds['e']}])
        self.assertEqual(single(self.post('traverse', dict(one, graphs=['ecoli'])))['results'],
                         single(self.post('traverse', dict(one, graph='ecoli')))['results'])
        # a name over two graphs: graphs selects as graph does (graph_path picks the graph); a
        # traverse runs on one pair, the name's are the candidates
        ret = self.server.post('traverse', dict(one, graphs=['enterics']))
        self.assertEqual(400, ret.status_code)
        self.assertEqual('traverse_single_pair', ret.json()['code'])
        self.assertIn("Bad request: index 'enterics' spans several graphs", ret.json()['error'])
        self.assertEqual([('enterics', self.pairs['ecoli'][0]), ('enterics', self.pairs['kleb'][0])],
                         [(c['graph'], c['graph_path']) for c in ret.json()['candidates']])
        picked = single(self.post('traverse', dict(one, graphs=['enterics'],
                                                   graph_path=self.pairs['ecoli'][0])))
        self.assertEqual(single(self.post('traverse', dict(one, graph='ecoli')))['results'],
                         picked['results'])
        # the pair the name and graph_path select is the entry's
        self.assertEqual(('enterics',) + self.pairs['ecoli'],
                         (picked['graph'], picked['graph_path'], picked['annotation_path']))
        for payload, words, code in ((dict(one, graphs=['ecoli'], graph='ecoli'), 'not both',
                                      'selection_conflict'),
                                     (dict(one, graphs=['ecoli', 'kleb']), 'expected [name]',
                                      'traverse_single_pair')):
            ret = self.server.post('traverse', payload)
            self.assertEqual(400, ret.status_code)
            self.assertIn(words, ret.json()['error'])
            self.assertEqual(code, ret.json()['code'])

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

    def test_resolve_entry_states_the_boundary_tie(self):
        """A pair's discovery cut passes through the envelope unchanged (SPEC-multigraph-
        envelope.md section 7, E18): labels_truncated with its boundary tie members, as the
        single-graph answer on the pair's files states it."""
        fastas = self.columns['pseudo']
        # 100 bp of one file's record and 160 bp of the other's: two labels of unequal support
        seq = (max(records_of(fastas[:1]), key=len)[5000:5100]
               + max(records_of(fastas[1:]), key=len)[5000:5160])
        r = {'sequence': seq, 'discover': {'max_labels': 1}}
        (entry,) = self.post('resolve', dict(r, graph='pseudo'))['answers']
        self.assertEqual('answered', entry['outcome'])
        cut = entry['labels_truncated']
        # the members of the reference sample 04
        self.assertEqual({'total', 'kept', 'min_kept_kmers', 'max_dropped_kmers',
                          'dropped_full_length', 'tied_at_boundary', 'tied_kept'}, set(cut))
        self.assertEqual((2, 1, 1, 1), (cut['total'], cut['kept'], cut['tied_at_boundary'],
                                        cut['tied_kept']), cut)
        # the single-graph answer, with the unknown_labels of a multi-graph server's pair
        # (section 2: present, possibly empty)
        self.assertEqual(untimed(dict(self.cli_resolve('pseudo', r), unknown_labels=[])),
                         untimed({k: v for k, v in entry.items() if k not in TAGS}))

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
        # the single-graph answer, with the unknown_labels of a multi-graph server's pair
        # (SPEC-multigraph-envelope.md section 2: present, possibly empty)
        self.assertEqual(untimed(dict(self.cli_resolve('kleb', r), unknown_labels=[])),
                         untimed(tagless))
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


@unittest.skipIf(PROTEIN_MODE, "the fixture is DNA")
class TestHeadersUnique(TestingBase):
    """GET /capabilities' headers_unique (SPEC-multigraph-envelope.md, "headers_unique (E22)"):
    whether a sequence header names one record over the server's pairs. Five pairs of small
    random records, each with its record mapping (.seqs) except one: `a` (two files) and `b`
    distinct, `c` repeating a header of `a` (across pairs), `e` repeating a header in its two
    files (across the columns of one pair), `d` without a .seqs."""

    maxDiff = None

    # pair -> its files -> the headers of their records
    PAIRS = {
        'a': {'a1.fa': ['r1', 'r2'], 'a2.fa': ['r3']},
        'b': {'b1.fa': ['r4', 'r5']},
        'c': {'c1.fa': ['r6', 'r1']},
        'd': {'d1.fa': ['r1', 'r7']},
        'e': {'e1.fa': ['r8', 'r9'], 'e2.fa': ['r9']},
    }

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        rng = random.Random(20261031)
        d = cls.tempdir.name
        cls.pairs = {}
        cls.files = {}
        for name, files in cls.PAIRS.items():
            out = os.path.join(d, name)
            os.makedirs(out)
            fastas = []
            for fname, headers in files.items():
                path = os.path.join(out, fname)
                with open(path, 'w') as f:
                    for h in headers:
                        f.write(f'>{h}\n{"".join(rng.choice("ACGT") for _ in range(80))}\n')
                fastas.append(path)
            cls.files[name] = fastas
            graph = os.path.join(out, 'graph.dbg')
            cls._build_graph(fastas, graph, 31, 'succinct', mode='basic')
            # both passes on one thread: each lists its columns in the order of the files, so
            # the .seqs's columns are the annotation's (several threads add the files' columns
            # in the order they finish)
            cls._annotate_graph(fastas, graph, os.path.join(out, 'anno'), 'column_coord',
                                anno_type='filename', num_threads=1)
            if name != 'd':
                cls._run_command(f'{METAGRAPH} annotate -p 1 -i {graph} --anno-filename '
                                 f'--index-header-coords -o {out}/anno {" ".join(fastas)}',
                                 'annotate --index-header-coords')
                assert os.path.exists(os.path.join(out, 'anno.seqs'))
            cls.pairs[name] = (graph, os.path.join(out, 'anno.column_coord.annodbg'))
            assert all(os.path.exists(p) for p in cls.pairs[name]), cls.pairs[name]

    def serve(self, names):
        """A multi-graph server of |names|' pairs, each under its own name"""
        csv = os.path.join(self.tempdir.name, '_'.join(names) + '.csv')
        with open(csv, 'w') as f:
            for name in names:
                f.write(f'{name},{self.pairs[name][0]},{self.pairs[name][1]}\n')
        return Server([csv], os.path.join(self.tempdir.name, '_'.join(names) + '.log'))

    def certificate(self, names):
        server = self.serve(names)
        try:
            caps = server.get('capabilities').json()
            # a multi-graph server states it in graph_summary only
            self.assertNotIn('headers_unique', caps)
            return caps['graph_summary']['headers_unique'], server.text()
        finally:
            server.stop()

    def occurrence(self, name, file):
        return {'graph': name, 'graph_path': self.pairs[name][0],
                'annotation_path': self.pairs[name][1], 'column': file}

    def test_distinct_headers_are_certified(self):
        certificate, log = self.certificate(['a', 'b'])
        self.assertEqual({'value': True, 'collisions': 0, 'examples': []}, certificate)
        self.assertIn('headers_unique: true', log)

    def test_repeats_across_pairs_and_across_columns(self):
        certificate, log = self.certificate(['a', 'b', 'c', 'e'])
        a1, a2 = self.files['a']
        c1, = self.files['c']
        e1, e2 = self.files['e']
        self.assertEqual({'value': False, 'collisions': 2, 'examples': [
            {'header': 'r1', 'occurrences': 2,
             'pairs': [self.occurrence('a', a1), self.occurrence('c', c1)]},
            {'header': 'r9', 'occurrences': 2,
             'pairs': [self.occurrence('e', e1), self.occurrence('e', e2)]},
        ]}, certificate)
        self.assertIn('headers_unique: false, 2 headers occur more than once', log)

    def test_single_graph_server_states_null_while_it_loads(self):
        # e's annotation and .seqs (a repeated header: false once loaded) beside a graph that
        # never loads, a FIFO nobody writes to (its loader blocks in open() for good, as
        # test_traverse.py's test_sigterm_while_the_index_loads): GET /capabilities answers at
        # once, ready false and headers_unique null, the member present
        fifo = os.path.join(self.tempdir.name, 'never_loads.dbg')
        if not os.path.exists(fifo):
            os.mkfifo(fifo)
        port = free_port()
        log_path = os.path.join(self.tempdir.name, 'single_loading.log')
        with open(log_path, 'w') as log:
            process = subprocess.Popen(
                shlex.split(METAGRAPH) + ['server_query', '-i', fifo, '-a', self.pairs['e'][1],
                                          '--port', str(port), '--address', '127.0.0.1',
                                          '-p', '2'],
                stdout=log, stderr=subprocess.STDOUT)
            try:
                caps = None
                for _ in range(300):
                    self.assertIsNone(process.poll(), 'the server exited while loading')
                    try:
                        ret = requests.get(f'http://127.0.0.1:{port}/capabilities', timeout=2)
                        if ret.ok:
                            caps = ret.json()
                            break
                    except requests.exceptions.RequestException:
                        pass
                    time.sleep(0.1)
                self.assertIsNotNone(caps, 'the server did not listen')
                self.assertFalse(caps['ready'])
                self.assertIn('headers_unique', caps)
                self.assertIsNone(caps['headers_unique'])
                self.assertIsNone(caps['graph_summary'])
                # and its pair's profile (E17), null until the index is loaded
                self.assertIn('pair_summary', caps)
                self.assertIsNone(caps['pair_summary'])
            finally:
                process.kill()
                process.wait()

    def test_the_reference_sample_states_it_as_a_server_does(self):
        # sample 12 of the contract (tests/data/multigraph_contract) states
        # graph_summary.headers_unique in a server's shape: false with one repeat across two
        # pairs, the same members at every level, the member in graph_summary's sorted order
        # beside the summary's other members, and no top-level member
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tests', 'data',
                            'multigraph_contract', '12_capabilities_multigraph.json')
        with open(path) as f:
            body = json.load(f)['body']
        self.assertNotIn('headers_unique', body)
        sample = body['graph_summary']
        self.assertEqual(sorted(sample), list(sample))
        server = self.serve(['a', 'b', 'c'])
        try:
            summary = server.get('capabilities').json()['graph_summary']
        finally:
            server.stop()
        self.assertEqual(sorted(summary), sorted(sample))
        live, stated = summary['headers_unique'], sample['headers_unique']
        self.assertEqual((False, 1, 1), (live['value'], live['collisions'], len(live['examples'])))
        self.assertEqual((False, 1, 1),
                         (stated['value'], stated['collisions'], len(stated['examples'])))
        self.assertEqual(list(live), list(stated))
        example, stated_example = live['examples'][0], stated['examples'][0]
        self.assertEqual(list(example), list(stated_example))
        self.assertEqual((2, 2), (example['occurrences'], len(example['pairs'])))
        self.assertEqual((2, 2), (stated_example['occurrences'], len(stated_example['pairs'])))
        for occurrence, stated_occurrence in zip(example['pairs'], stated_example['pairs']):
            self.assertEqual(list(occurrence), list(stated_occurrence))
        # the sample's occurrences name pairs of its graph_summary
        pairs = {(p['graph'], p['graph_path'], p['annotation_path']) for p in sample['pairs']}
        for occurrence in stated_example['pairs']:
            self.assertIn((occurrence['graph'], occurrence['graph_path'],
                           occurrence['annotation_path']), pairs)

    def test_a_pair_listed_under_two_names_is_counted_once(self):
        # a's pair also under a second name: its headers are no repeats of themselves, and its
        # occurrences name the first name in byte order ('a', not 'z_alias')
        csv = os.path.join(self.tempdir.name, 'alias.csv')
        with open(csv, 'w') as f:
            for name, pair in (('z_alias', 'a'), ('a', 'a'), ('c', 'c')):
                f.write(f'{name},{self.pairs[pair][0]},{self.pairs[pair][1]}\n')
        server = Server([csv], os.path.join(self.tempdir.name, 'alias.log'))
        try:
            certificate = server.get('capabilities').json()['graph_summary']['headers_unique']
        finally:
            server.stop()
        a1, _ = self.files['a']
        c1, = self.files['c']
        self.assertEqual({'value': False, 'collisions': 1, 'examples': [
            {'header': 'r1', 'occurrences': 2,
             'pairs': [self.occurrence('a', a1), self.occurrence('c', c1)]}]}, certificate)

    def test_a_pair_without_record_mapping_makes_it_unknown(self):
        # d's r1 is a column label's record, not a header the server knows
        certificate, log = self.certificate(['a', 'd'])
        self.assertEqual({'value': None, 'collisions': 0, 'examples': []}, certificate)
        self.assertIn('headers_unique: null', log)
        # the pairs with a .seqs are still compared: c repeats a's r1
        certificate, _ = self.certificate(['a', 'c', 'd'])
        self.assertIsNone(certificate['value'])
        self.assertEqual(1, certificate['collisions'])
        self.assertEqual('r1', certificate['examples'][0]['header'])

    def test_single_graph_server_states_it_at_the_top(self):
        for name, value in (('a', True), ('e', False), ('d', None)):
            graph, anno = self.pairs[name]
            server = Server(['-i', graph, '-a', anno],
                            os.path.join(self.tempdir.name, f'single_{name}.log'))
            try:
                caps = server.get('capabilities').json()
            finally:
                server.stop()
            with self.subTest(pair=name):
                self.assertTrue(caps['ready'])
                self.assertIsNone(caps['graph_summary'])
                certificate = caps['headers_unique']
                self.assertEqual(value, certificate['value'])
                if name == 'e':
                    e1, e2 = self.files['e']
                    # a single-graph server's pair tags are null, as its /search cuts'
                    untagged = {'graph': None, 'graph_path': None, 'annotation_path': None}
                    self.assertEqual({'value': False, 'collisions': 1, 'examples': [
                        {'header': 'r9', 'occurrences': 2,
                         'pairs': [dict(untagged, column=e1), dict(untagged, column=e2)]}]},
                        certificate)
                else:
                    self.assertEqual((0, []), (certificate['collisions'],
                                               certificate['examples']))



@unittest.skipIf(PROTEIN_MODE, "the fixture is DNA")
class TestCapabilityProfile(TestingBase):
    """The capability profile of each pair (E17, SPEC-multigraph-envelope.md section 8.1) on one
    multi-graph server of every annotation form the integration tests build (test_query.py's),
    each over a copy of one graph of three random records (k 31, basic, no mask), annotated by
    file name, with its record mapping (.seqs) on the coordinate forms; and a second chunk
    (other records) of the column form, as two chunks of one build. The annotation form is the
    form built, decode_charged the oracle's (what a walk under a memory budget states:
    memory_bound_soft says the reads are charged before they are held exactly where it is
    true), counts on the count forms, and support_fp recomputed from each entry: one per form,
    the two chunks' one. A single-graph server states the same profile as pair_summary."""

    maxDiff = None

    FORMS = ('column', 'column_coord', 'brwt_coord', 'row_diff_coord', 'row_diff_brwt_coord',
             'row_diff_disk_coord', 'row', 'row_diff', 'row_sparse', 'row_diff_brwt',
             'row_diff_disk', 'row_diff_flat', 'row_diff_sparse', 'rb_brwt', 'brwt', 'int_brwt',
             'row_diff_int_brwt', 'row_diff_int_disk', 'rbfish', 'flat')
    # the forms whose reads a traversal charges under a memory budget (E14; the unit test
    # MultiGraphProfile.AnnotationFormAndDecodeChargedOnEveryForm covers every type)
    CHARGED = {'column', 'brwt', 'row_diff', 'row_diff_coord', 'row_diff_brwt',
               'row_diff_brwt_coord', 'row_diff_disk', 'row_diff_sparse', 'row_diff_flat',
               'row_diff_int_brwt'}
    COUNTS = {'int_brwt', 'row_diff_int_brwt', 'row_diff_int_disk'}
    # the forms whose .seqs is written as before H4 (no column labels), as r237's
    COUNT_ONLY = {'row_diff_brwt_coord'}

    @classmethod
    def write_records(cls, path, headers, seed):
        rng = random.Random(seed)
        with open(path, 'w') as f:
            for h in headers:
                f.write(f'>{h}\n{"".join(rng.choice("ACGT") for _ in range(400))}\n')

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        d = cls.tempdir.name
        cls.fasta = os.path.join(d, 'records.fa')
        cls.write_records(cls.fasta, ['r1', 'r2', 'r3'], 20261108)
        cls._build_graph([cls.fasta], os.path.join(d, 'graph.dbg'), 31, 'succinct', mode='basic')
        cls.pairs = {}
        for form in cls.FORMS:
            out = os.path.join(d, form)
            os.makedirs(out)
            graph = os.path.join(out, 'graph.dbg')
            shutil.copyfile(os.path.join(d, 'graph.dbg'), graph)
            cls._annotate_graph([cls.fasta], graph, os.path.join(out, 'anno'), form,
                                anno_type='filename', num_threads=1)
            if form.endswith('_coord'):
                cls._run_command(f'{METAGRAPH} annotate -p 1 -i {graph} --anno-filename '
                                 f'--index-header-coords -o {out}/anno {cls.fasta}',
                                 'annotate --index-header-coords')
                # r237's profile: its .seqs written before H4, without column labels
                # (seqs_order "count_only"); the other forms' state theirs ("verified")
                if form in cls.COUNT_ONLY:
                    strip_column_labels(os.path.join(out, 'anno.seqs'))
            cls.pairs[form] = (graph, os.path.join(out, f'anno.{form}.annodbg'))
            assert all(os.path.exists(p) for p in cls.pairs[form]), cls.pairs[form]
        # a second chunk of the column build: other records, other columns
        out = os.path.join(d, 'chunk2')
        os.makedirs(out)
        fasta2 = os.path.join(out, 'records2.fa')
        cls.write_records(fasta2, ['s1', 's2'], 7)
        fasta3 = os.path.join(out, 'records3.fa')
        cls.write_records(fasta3, ['s3'], 8)
        cls._build_graph([fasta2, fasta3], os.path.join(out, 'graph.dbg'), 31, 'succinct',
                         mode='basic')
        cls._annotate_graph([fasta2, fasta3], os.path.join(out, 'graph.dbg'),
                            os.path.join(out, 'anno'), 'column', anno_type='filename',
                            num_threads=1)
        cls.pairs['chunk2'] = (os.path.join(out, 'graph.dbg'),
                               os.path.join(out, 'anno.column.annodbg'))
        cls.csv = os.path.join(d, 'forms.csv')
        with open(cls.csv, 'w') as f:
            for name, (graph, anno) in cls.pairs.items():
                f.write(f'{name},{graph},{anno}\n')
        cls.server = Server([cls.csv], os.path.join(d, 'server.log'))
        cls.caps = cls.server.get('capabilities').json()
        cls.entries = {p['graph']: p for p in cls.caps['graph_summary']['pairs']}

    @classmethod
    def tearDownClass(cls):
        cls.server.stop()
        super().tearDownClass()

    def test_every_pair_states_its_profile(self):
        self.assertEqual(set(self.pairs), set(self.entries))
        for form in self.FORMS + ('chunk2',):
            p = self.entries[form]
            built = 'column' if form == 'chunk2' else form
            with self.subTest(form=form):
                self.assertEqual(self.pairs[form], (p['graph_path'], p['annotation_path']))
                self.assertEqual(built, p['annotation_form'])
                self.assertEqual(built in self.CHARGED, p['decode_charged'])
                self.assertEqual('succinct', p['graph_representation'])
                self.assertEqual('presence_only' if built in self.COUNTS else None,
                                 p['traversal']['counts'])
                coordinates = built.endswith('_coord')
                self.assertEqual((coordinates, coordinates, coordinates),
                                 (p['traversal']['has_coordinates'],
                                  p['traversal']['has_coord_to_header'],
                                  p['traversal']['supports_trace']))
                # what the load verified of the .seqs's column order (H4, section 9.1)
                self.assertEqual(None if not coordinates
                                 else 'count_only' if form in self.COUNT_ONLY else 'verified',
                                 p['seqs_order'])
                self.assertEqual(('basic', 31, True, 'absent', 'upper_bound'),
                                 (p['traversal']['regime'], p['k'], p['available'], p['mask'],
                                  p['counting']))
                # what the pair's probe states
                probe = single(self.server.get(f'traverse/capabilities?graph={form}').json())
                self.assertEqual((probe['index_meta_fp'], probe['alphabet'], probe.get('counts')),
                                 (p['index_meta_fp'], p['alphabet'], p['traversal']['counts']))
                self.assertEqual(support_fp(p), p['support_fp'])
        # one profile per form; the two chunks of the column build share theirs, their
        # identity and size apart
        fps = {form: self.entries[form]['support_fp'] for form in self.FORMS}
        self.assertEqual(len(self.FORMS), len(set(fps.values())))
        column, chunk2 = self.entries['column'], self.entries['chunk2']
        self.assertEqual(column['support_fp'], chunk2['support_fp'])
        self.assertNotEqual(column['index_meta_fp'], chunk2['index_meta_fp'])
        self.assertEqual((1, 2), (column['traversal']['num_labels'],
                                  chunk2['traversal']['num_labels']))
        log = self.server.text()
        self.assertRegex(log, r'index_meta_fp of the %d pairs \(%d columns\) computed in '
                              r'\d+\.\d+ s' % (len(self.pairs), len(self.FORMS) + 2))
        self.assertIn(f'{len(self.FORMS)} capability profile(s) (support_fp)', log)

    def test_decode_charged_is_what_a_walk_under_a_memory_budget_does(self):
        """A walk under a memory budget states memory_bound_soft: its reads charged before they
        are held where the oracle charges them (SPEC-labeled-traversal-core.md §6.8), observed
        after they are held elsewhere -- on every form, as decode_charged states it"""
        seed = self.records()[0][100:160]
        for form in self.FORMS:
            with self.subTest(form=form):
                ret = self.server.post('traverse', {
                    'graph': form, 'seeds': [{'sequence': seed}],
                    'strategy': {'labels': {'mode': 'annotate'},
                                 'bounds': {'max_memory_mb': 64, 'max_extension_bp': 20},
                                 'output': {'detail': 'summary', 'timing': False}}})
                self.assertEqual(200, ret.status_code, ret.text[:2000])
                (result,) = single(ret.json(), '/traverse')['results']
                (soft,) = [l for l in result['limitations'] if l['kind'] == 'memory_bound_soft']
                self.assertEqual(self.entries[form]['decode_charged'], CHARGED in soft['effect'],
                                 soft['effect'])

    def records(self):
        return [seq for _, seq in named_records_of([self.fasta])]

    def test_a_single_graph_server_states_its_pairs_profile(self):
        """pair_summary: the graph_summary entry of the server's one pair, its tags null, the
        same profile and support_fp (a count form and a coordinate form)"""
        for form in ('row_diff_int_brwt', 'row_diff_brwt_coord'):
            graph, anno = self.pairs[form]
            server = Server(['-i', graph, '-a', anno],
                            os.path.join(self.tempdir.name, f'single_{form}.log'))
            try:
                caps = server.get('capabilities').json()
                log = server.text()
            finally:
                server.stop()
            with self.subTest(form=form):
                self.assertTrue(caps['ready'])
                self.assertIsNone(caps['graph_summary'])
                expected = dict(self.entries[form], graph=None, graph_path=None,
                                annotation_path=None)
                self.assertEqual(expected, caps['pair_summary'])
                self.assertEqual(support_fp(caps['pair_summary']), caps['pair_summary']['support_fp'])
                self.assertIn(f'Capability profile: succinct {form} annotation, support_fp '
                              f'{expected["support_fp"]}, index_meta_fp '
                              f'{expected["index_meta_fp"]}', log)

    def test_the_reference_sample_states_the_profile_as_a_server_does(self):
        """Sample 12 (tests/data/multigraph_contract) states each pair's profile in a server's
        shape: the members of a live entry at every level, each support_fp recomputed from its
        entry, and the value of a live pair of r237's profile (succinct, k 31, basic,
        row_diff_brwt_coord with its .seqs, no mask): this server's row_diff_brwt_coord pair"""
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tests', 'data',
                            'multigraph_contract', '12_capabilities_multigraph.json')
        with open(path) as f:
            body = json.load(f)['body']
        self.assertNotIn('pair_summary', body)
        live = self.entries['row_diff_brwt_coord']
        for pair in body['graph_summary']['pairs']:
            with self.subTest(pair=pair['graph']):
                self.assertEqual(sorted(live), list(pair))
                self.assertEqual(sorted(live['traversal']), list(pair['traversal']))
                self.assertEqual(support_fp(pair), pair['support_fp'])
                self.assertEqual(live['support_fp'], pair['support_fp'])
                self.assertRegex(pair['index_meta_fp'], r'^[0-9a-f]{16}$')


def run_to_exit(args, log_path, timeout=300):
    """server_query |args| on a free port, run until it exits by itself: (exit status, log)"""
    with open(log_path, 'wb') as log:
        res = subprocess.run(shlex.split(METAGRAPH) + ['server_query'] + args
                             + ['--port', str(free_port()), '--address', '127.0.0.1', '-p', '2'],
                             stdout=log, stderr=subprocess.STDOUT, timeout=timeout)
    with open(log_path) as f:
        return res.returncode, f.read()


@unittest.skipIf(PROTEIN_MODE, "the fixture is DNA")
class TestRecordMappingAtLoad(TestingBase):
    """A record mapping (.seqs) that is not its annotation's is refused when the index loads
    (SPEC-multigraph-envelope.md section 10, H4), never by a request: another number of
    columns; stated column labels (written by annotate --index-header-coords since H4) that are
    not the annotation's, in its order; on column_coord, a column holding more coordinates than
    the mapping states k-mers for it. The fixture: one graph of big.fa (two records of 300 bp)
    and small.fa (one of 60 bp), annotated by file name on one thread (columns big.fa,
    small.fa), as column_coord and row_diff_brwt_coord; each variant of the .seqs beside a copy
    of the annotation. A .seqs without labels (as earlier binaries write it) is served with
    seqs_order "count_only"; on the row-diff form, where only the count is checked, even when
    its columns are swapped."""

    maxDiff = None
    REFUSAL = "does not belong to the annotation"

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        d = cls.d = cls.tempdir.name
        rng = random.Random(20261108)
        records = {'big.fa': [('b1', 300), ('b2', 300)], 'small.fa': [('s1', 60)]}
        for fname, recs in records.items():
            with open(os.path.join(d, fname), 'w') as f:
                for header, length in recs:
                    f.write(f'>{header}\n{"".join(rng.choice("ACGT") for _ in range(length))}\n')

        def run(cmd):
            res = subprocess.run(shlex.split(METAGRAPH) + cmd, cwd=d, stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT)
            assert res.returncode == 0, (cmd, res.stdout.decode()[-3000:])

        run(['build', '-p', '1', '--mode', 'basic', '--graph', 'succinct', '-k', '31',
             '-o', 'graph', 'big.fa', 'small.fa'])
        run(['annotate', '-p', '1', '--anno-filename', '-i', 'graph.dbg', '--anno-type', 'column',
             '--coordinates', '-o', 'anno', 'big.fa', 'small.fa'])
        run(['transform_anno', '-p', '1', '--anno-type', 'column_coord', '--coordinates',
             '-o', 'anno', 'anno.column.annodbg'])
        rd = ['transform_anno', '-p', '1', '--anno-type', 'row_diff', '--coordinates', '-o',
              'anno', '-i', 'graph.dbg', 'anno.column.annodbg']
        run(rd)
        run(rd + ['--row-diff-stage', '1'])
        run(rd + ['--row-diff-stage', '2'])
        run(['transform_anno', '-p', '1', '--anno-type', 'row_diff_brwt_coord', '--greedy',
             '-o', 'anno', '-i', 'graph.dbg', 'anno.column.annodbg'])
        run(['relax_brwt', '-p', '1', '-o', 'anno', 'anno.row_diff_brwt_coord.annodbg'])
        # an annotation of three columns over the same graph (small.fa's record again, in
        # third.fa), for an annotation replaced under a running server
        with open(os.path.join(d, 'small.fa')) as f, open(os.path.join(d, 'third.fa'), 'w') as g:
            g.write(f.read().replace('>s1', '>t1'))
        run(['annotate', '-p', '1', '--anno-filename', '-i', 'graph.dbg', '--anno-type', 'column',
             '--coordinates', '-o', 'three', 'big.fa', 'small.fa', 'third.fa'])
        run(['transform_anno', '-p', '1', '--anno-type', 'column_coord', '--coordinates',
             '-o', 'three', 'three.column.annodbg'])
        # a rebuild of the same files in the other order, with its own .seqs: consistent, and
        # of the same column count as the annotation it replaces under a running server
        run(['annotate', '-p', '1', '--anno-filename', '-i', 'graph.dbg', '--anno-type', 'column',
             '--coordinates', '-o', 'rev', 'small.fa', 'big.fa'])
        run(['transform_anno', '-p', '1', '--anno-type', 'column_coord', '--coordinates',
             '-o', 'rev', 'rev.column.annodbg'])
        run(['annotate', '-p', '1', '-i', 'graph.dbg', '--anno-filename', '--index-header-coords',
             '-o', 'rev', 'small.fa', 'big.fa'])
        # the mappings: of the files in column order, in the other order, and of one file
        for name, files in (('right', ['big.fa', 'small.fa']), ('swapped', ['small.fa', 'big.fa']),
                            ('subset', ['big.fa'])):
            run(['annotate', '-p', '1', '-i', 'graph.dbg', '--anno-filename',
                 '--index-header-coords', '-o', name, *files])
        cls.graph = os.path.join(d, 'graph.dbg')

    def pair(self, form, mapping, labels=True):
        """A copy of the |form| annotation beside a copy of the |mapping| .seqs (without its
        labels unless |labels|): (graph, annotation)"""
        out = os.path.join(self.d, f'{form}_{mapping}_{"labels" if labels else "old"}')
        if not os.path.isdir(out):
            os.makedirs(out)
            shutil.copyfile(os.path.join(self.d, f'anno.{form}.annodbg'),
                            os.path.join(out, f'anno.{form}.annodbg'))
            shutil.copyfile(os.path.join(self.d, f'{mapping}.seqs'), os.path.join(out, 'anno.seqs'))
            if not labels:
                strip_column_labels(os.path.join(out, 'anno.seqs'))
        return self.graph, os.path.join(out, f'anno.{form}.annodbg')

    def test_the_mapping_states_the_files_as_passed(self):
        # the columns, then the marker and the labels: the files' names as passed
        with open(os.path.join(self.d, 'swapped.seqs'), 'rb') as f:
            data = f.read()
        at = data.rfind(b'mtg-cols')
        self.assertGreater(at, 0)
        self.assertIn(b'small.fa', data[at:])
        self.assertLess(data[at:].index(b'small.fa'), data[at:].index(b'big.fa'))
        # stats states whether a .seqs states its labels
        stats = subprocess.run(shlex.split(METAGRAPH) + ['stats', os.path.join(self.d, 'right.seqs')],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(0, stats.returncode, stats.stderr.decode())
        self.assertIn('columns: 2', stats.stdout.decode())
        self.assertIn('column labels: stated', stats.stdout.decode())
        old = os.path.join(self.d, 'old.seqs')
        shutil.copyfile(os.path.join(self.d, 'right.seqs'), old)
        strip_column_labels(old)
        stats = subprocess.run(shlex.split(METAGRAPH) + ['stats', old],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(0, stats.returncode, stats.stderr.decode())
        self.assertIn('column labels: none', stats.stdout.decode())

    def test_a_mapping_of_the_annotation_is_served(self):
        for form in ('column_coord', 'row_diff_brwt_coord'):
            for labels, order in ((True, 'verified'), (False, 'count_only')):
                graph, anno = self.pair(form, 'right', labels)
                server = Server(['-i', graph, '-a', anno],
                                os.path.join(self.d, f'served_{form}_{labels}.log'))
                try:
                    caps = server.get('capabilities').json()
                finally:
                    server.stop()
                with self.subTest(form=form, labels=labels):
                    self.assertTrue(caps['ready'])
                    self.assertEqual(order, caps['pair_summary']['seqs_order'])

    def test_corrupt_coordinate_offsets_are_refused_at_load(self):
        # a .seqs whose first column's offsets state a length past their encoding is refused at
        # load, naming the file: a server exits 1 before it is ready, the CLI exits 1 (2^20 and
        # more segfaulted, exit 139, before the train-3 verification, V2-S3)
        graph, anno = self.pair('row_diff_brwt_coord', 'right', True)
        with open(anno[:-len('.row_diff_brwt_coord.annodbg')] + '.seqs', 'rb') as f:
            data = bytearray(f.read())
        # the number of columns, big.fa's headers b1 and b2, then its offsets' length
        at = 8 + 8 + sum(1 + len(h) for h in ('b1', 'b2'))
        self.assertEqual(2 * (300 - 31 + 1), int.from_bytes(data[at:at + 8], 'little'))
        out = os.path.join(self.d, 'corrupt_length')
        os.makedirs(out)
        corrupt = os.path.join(out, 'anno.row_diff_brwt_coord.annodbg')
        shutil.copyfile(anno, corrupt)
        seqs = os.path.join(out, 'anno.seqs')
        for length in (1 << 20, (1 << 64) - 1):
            data[at:at + 8] = length.to_bytes(8, 'little')
            with open(seqs, 'wb') as f:
                f.write(data)
            line = (f"Cannot load CoordToHeader from '{seqs}': could not read contents (possibly "
                    f"corrupted or incompatible) (caught: the coordinate offsets of column 0 state "
                    f"a length of {length}, more than their encoding holds")
            with self.subTest(length=length):
                status, log = run_to_exit(['-i', graph, '-a', corrupt],
                                          os.path.join(self.d, f'corrupt_length_{length}.log'))
                self.assertEqual(1, status, log)
                self.assertIn(line, log)
                self.assertNotIn('Annotated graph loaded', log)
                res = subprocess.run(shlex.split(METAGRAPH) + [
                        'query', '--query-mode', 'coords', '-i', graph, '-a', corrupt,
                        os.path.join(self.d, 'small.fa')],
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                self.assertEqual(1, res.returncode, res.stderr.decode())
                self.assertIn(line, res.stderr.decode())

    def test_another_mapping_is_refused_at_start_up(self):
        cases = [
            ('column_coord', 'subset', True, "it has 1 columns, the annotation 2 labels"),
            ('row_diff_brwt_coord', 'subset', False, "it has 1 columns, the annotation 2 labels"),
            ('column_coord', 'swapped', True,
             "its column 0 is 'small.fa', the annotation's column 0 'big.fa' (the files were "
             "passed in another order or under other names)"),
            ('row_diff_brwt_coord', 'swapped', True,
             "its column 0 is 'small.fa', the annotation's column 0 'big.fa'"),
            # without labels, column_coord's coordinates do not fit the swapped k-mers: big.fa's
            # 2 * (300 - 31 + 1) coordinates against small.fa's 60 - 31 + 1 k-mers
            ('column_coord', 'swapped', False,
             "column 0 ('big.fa') of the annotation holds 540 coordinates, more than the 30 "
             "k-mers the mapping states for it (the files were passed in another order)"),
        ]
        for form, mapping, labels, reason in cases:
            graph, anno = self.pair(form, mapping, labels)
            seqs = anno[:-len(f'.{form}.annodbg')] + '.seqs'
            line = (f"The record mapping '{seqs}' does not belong to the annotation '{anno}': "
                    f"{reason}")
            with self.subTest(form=form, mapping=mapping, labels=labels):
                # a single-graph server: exits 1 before it is ready
                status, log = run_to_exit(['-i', graph, '-a', anno],
                                          os.path.join(self.d, f'single_{form}_{mapping}_{labels}.log'))
                self.assertEqual(1, status, log)
                self.assertIn(line, log)
                self.assertNotIn('Annotated graph loaded', log)
                # a multi-graph server, the pair beside a good one: exits 1 before it is ready
                good = self.pair(form, 'right', True)
                csv = os.path.join(self.d, f'list_{form}_{mapping}_{labels}.csv')
                with open(csv, 'w') as f:
                    f.write(f'good,{good[0]},{good[1]}\nbad,{graph},{anno}\n')
                status, log = run_to_exit([csv], csv + '.log')
                self.assertEqual(1, status, log)
                self.assertIn(line, log)
                self.assertNotIn('Ready to serve queries', log)
                # the CLI, at its load
                res = subprocess.run(shlex.split(METAGRAPH) + [
                        'query', '--query-mode', 'coords', '-i', graph, '-a', anno,
                        os.path.join(self.d, 'small.fa')],
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                self.assertEqual(1, res.returncode, res.stderr.decode())
                self.assertIn(line, res.stderr.decode())
                # and none with --no-coord-mapping
                res = subprocess.run(shlex.split(METAGRAPH) + [
                        'query', '--query-mode', 'coords', '--no-coord-mapping', '-i', graph,
                        '-a', anno, os.path.join(self.d, 'small.fa')],
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                self.assertEqual(0, res.returncode, res.stderr.decode())

    def test_annotate_with_the_annotation_checks_the_files(self):
        # -a: each file is spot-checked against its column of the annotation (the column holds
        # the coordinates of its first and last k-mers in the graph), and the annotation's labels
        # are stated, so that files at another path pass the load's check; without -a they are
        # refused there
        def annotate(out, *args):
            return subprocess.run(shlex.split(METAGRAPH) + [
                    'annotate', '-p', '1', '-i', 'graph.dbg', '--anno-filename',
                    '--index-header-coords', '-o', out, *args],
                    cwd=self.d, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        absolute = [os.path.join(self.d, f) for f in ('big.fa', 'small.fa')]
        for form in ('column_coord', 'row_diff_brwt_coord'):
            anno = f'anno.{form}.annodbg'
            with self.subTest(form=form):
                res = annotate(f'checked_{form}', '-a', anno, *absolute)
                self.assertEqual(0, res.returncode, res.stderr.decode())
                with open(os.path.join(self.d, f'checked_{form}.seqs'), 'rb') as f:
                    data = f.read()
                trailer = data[data.rfind(b'mtg-cols'):]
                self.assertIn(b'big.fa', trailer)
                self.assertNotIn(self.d.encode(), trailer)
                res = annotate(f'swapped_{form}', '-a', anno, 'small.fa', 'big.fa')
                self.assertEqual(1, res.returncode)
                self.assertIn(f"The file 'small.fa' is not column 0 ('big.fa') of the annotation "
                              f"'{anno}': the column lacks the file's k-mer at coordinate",
                              res.stderr.decode())
                self.assertFalse(os.path.exists(os.path.join(self.d, f'swapped_{form}.seqs')))
                res = annotate(f'subset_{form}', '-a', anno, 'big.fa')
                self.assertEqual(1, res.returncode)
                self.assertIn(f"The annotation '{anno}' has 2 columns, and 1 files were passed",
                              res.stderr.decode())
        # served: with -a's labels; without, the absolute names are refused at load
        out = os.path.join(self.d, 'moved')
        os.makedirs(out)
        shutil.copyfile(os.path.join(self.d, 'anno.column_coord.annodbg'),
                        os.path.join(out, 'anno.column_coord.annodbg'))
        shutil.copyfile(os.path.join(self.d, 'checked_column_coord.seqs'),
                        os.path.join(out, 'anno.seqs'))
        server = Server(['-i', self.graph, '-a', os.path.join(out, 'anno.column_coord.annodbg')],
                        os.path.join(self.d, 'moved.log'))
        try:
            self.assertEqual('verified',
                             server.get('capabilities').json()['pair_summary']['seqs_order'])
        finally:
            server.stop()
        res = annotate(os.path.join(out, 'anno'), *absolute)
        self.assertEqual(0, res.returncode, res.stderr.decode())
        status, log = run_to_exit(['-i', self.graph, '-a',
                                   os.path.join(out, 'anno.column_coord.annodbg')],
                                  os.path.join(self.d, 'moved_refused.log'))
        self.assertEqual(1, status, log)
        self.assertIn(f"its column 0 is '{absolute[0]}', the annotation's column 0 'big.fa'", log)

    def test_annotate_with_the_annotation_is_a_spot_check(self):
        # -a's check is a spot check (SPEC-multigraph-envelope.md sections 9.1 and 10; GPT review
        # 6, finding 12): a file's probes are its first stretch with a k-mer the graph has and its
        # last one in the end of its last record, past a gap of Ns longer than a probe. A file
        # without a k-mer of the graph is refused, nothing written: one with k-mers none of which
        # the graph has (before GPT review 6 nothing was checked and the pair was served
        # "verified"), and one without k-mers, empty or of records shorter than k (so too before
        # the train-3 review, R4 S1), which no column is made from. The same sequences under
        # other headers pass: what the records are rests on the files passed, which no spot
        # check can tell
        d = os.path.join(self.d, 'spot')
        os.makedirs(d)
        rng = random.Random(20261109)
        seq = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
        gap = 'N' * 2500
        g1, g2, p1 = seq(300), seq(300), seq(300)
        files = {
            'gapped.fa': [('g1', gap + g1), ('g2', g2 + gap)],
            'plain.fa': [('p1', p1)],
            # the gaps, and plain.fa's record between them: in the graph, not in column 0
            'other.fa': [('o1', gap + p1), ('o2', p1 + gap)],
            'foreign.fa': [('f1', seq(300)), ('f2', seq(300))],
            'renamed.fa': [('r1', gap + g1), ('r2', g2 + gap)],
            # gapped.fa's records and one of Ns alone, longer than the end the tail is looked in
            'trailing.fa': [('g1', gap + g1), ('g2', g2 + gap), ('g3', 'N' * 70000)],
            # no k-mers: an empty file (an interrupted copy), records shorter than k
            'empty.fa': [],
            'short.fa': [('s1', seq(20)), ('s2', seq(25))],
        }
        for fname, recs in files.items():
            with open(os.path.join(d, fname), 'w') as f:
                f.writelines(f'>{h}\n{s}\n' for h, s in recs)

        def run(*cmd):
            return subprocess.run(shlex.split(METAGRAPH) + list(cmd), cwd=d,
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        for cmd in (['build', '-p', '1', '--mode', 'basic', '--graph', 'succinct', '-k', '31',
                     '-o', 'graph', 'gapped.fa', 'plain.fa'],
                    ['annotate', '-p', '1', '--anno-filename', '-i', 'graph.dbg', '--anno-type',
                     'column', '--coordinates', '-o', 'anno', 'gapped.fa', 'plain.fa'],
                    ['transform_anno', '-p', '1', '--anno-type', 'column_coord', '--coordinates',
                     '-o', 'anno', 'anno.column.annodbg']):
            res = run(*cmd)
            self.assertEqual(0, res.returncode, res.stderr.decode()[-3000:])
        anno = 'anno.column_coord.annodbg'

        def annotate(out, *inputs):
            return run('annotate', '-p', '1', '-i', 'graph.dbg', '--anno-filename',
                       '--index-header-coords', '-a', anno, '-o', out, *inputs)
        seqs = lambda out: os.path.join(d, f'{out}.seqs')
        # the probes pass over the gaps, and are checked
        res = annotate('gapped', 'gapped.fa', 'plain.fa')
        self.assertEqual(0, res.returncode, res.stderr.decode())
        self.assertNotIn('has no k-mer of the graph', res.stderr.decode())
        res = annotate('other', 'other.fa', 'plain.fa')
        self.assertEqual(1, res.returncode, res.stderr.decode())
        self.assertIn(f"The file 'other.fa' is not column 0 ('gapped.fa') of the annotation "
                      f"'{anno}': the column lacks the file's k-mer at coordinate 2500",
                      res.stderr.decode())
        self.assertFalse(os.path.exists(seqs('other')))
        # no k-mers: no column is made from such a file (each column holds a k-mer of its file),
        # so refusing it refuses no legitimate file
        res = run('annotate', '-p', '1', '--anno-filename', '-i', 'graph.dbg', '--anno-type',
                  'column', '--coordinates', '-o', 'four', 'gapped.fa', 'empty.fa', 'short.fa',
                  'plain.fa')
        self.assertEqual(0, res.returncode, res.stderr.decode()[-3000:])
        res = run('stats', '-a', 'four.column.annodbg', '--print-col-names')
        self.assertIn('Number of columns: 2\ngapped.fa\nplain.fa\n', res.stdout.decode())
        for name in ('empty.fa', 'short.fa'):
            out = name[:-len('.fa')]
            res = annotate(out, 'gapped.fa', name)
            self.assertEqual(1, res.returncode, res.stderr.decode())
            self.assertIn(f"The file '{name}' has no k-mers (it is empty, or each of its records is "
                          f"shorter than k = 31): no annotation column is made from such a file, so "
                          f"it cannot be column 1 ('plain.fa') of the annotation '{anno}'. Pass the "
                          f"files the annotation was made from, in its column order",
                          res.stderr.decode())
            self.assertFalse(os.path.exists(seqs(out)))
        # no k-mer of the graph: refused
        res = annotate('foreign', 'gapped.fa', 'foreign.fa')
        self.assertEqual(1, res.returncode, res.stderr.decode())
        self.assertIn(f"The file 'foreign.fa' has no k-mer of the graph: it cannot be checked "
                      f"against column 1 ('plain.fa') of the annotation '{anno}'. Pass the files "
                      f"the annotation was made from, in its column order", res.stderr.decode())
        self.assertFalse(os.path.exists(seqs('foreign')))
        # a last record without a k-mer of the graph in its end: the head alone, with a warning
        res = annotate('trailing', 'trailing.fa', 'plain.fa')
        self.assertEqual(0, res.returncode, res.stderr.decode())
        self.assertIn("The last record of 'trailing.fa' has no k-mer of the graph in its last "
                      "65536 characters: the file is checked by its first k-mers alone",
                      res.stderr.decode())
        # the limit: the same sequences under other headers pass, the mapping naming r1 and r2
        res = annotate('renamed', 'renamed.fa', 'plain.fa')
        self.assertEqual(0, res.returncode, res.stderr.decode())
        with open(seqs('renamed'), 'rb') as f:
            data = f.read()
        self.assertIn(b'r1', data)
        self.assertNotIn(b'g1', data)

    def test_a_search_load_into_ram_never_exits_the_server(self):
        # /search with in_ram on a --mmap server loads the pair with the resident pair's
        # mapping, checked at start-up, as /traverse, /resolve and /pattern do, while the
        # annotation on disk is the file loaded then: a .seqs replaced alone is not read again;
        # an annotation replaced has its .seqs read and checked. A load that fails (here the
        # annotation replaced by one of three columns, and the .seqs on disk of two columns
        # refused) fails that request alone, the route's 400 with the load's
        # message as error (as /traverse, /resolve and /pattern answer a load that fails), and
        # the server runs on (H4 review S2: it exited)
        out = os.path.join(self.d, 'live')
        os.makedirs(out)
        anno = os.path.join(out, 'anno.column_coord.annodbg')
        seqs = os.path.join(out, 'anno.seqs')
        shutil.copyfile(os.path.join(self.d, 'anno.column_coord.annodbg'), anno)
        shutil.copyfile(os.path.join(self.d, 'right.seqs'), seqs)

        def replace(source, target):
            # a new file renamed over the old: the resident's mapped file keeps its inode
            shutil.copyfile(source, target + '.new')
            os.replace(target + '.new', target)

        csv = os.path.join(self.d, 'live.csv')
        with open(csv, 'w') as f:
            f.write(f'live,{self.graph},{anno}\n')
        with open(os.path.join(self.d, 'small.fa')) as f:
            fasta = f.read()
        query = {'FASTA': fasta, 'graphs': ['live'], 'discovery_fraction': 0.5}
        server = Server([csv, '--mmap', '--mem-cap-gb', '1', '--cache-size', '0'],
                        os.path.join(self.d, 'live.log'))
        try:
            resident = server.post('search', query)
            self.assertEqual(200, resident.status_code, resident.text)
            self.assertTrue(resident.json()[0]['results'], resident.text)
            loaded = server.post('search', dict(query, in_ram=True))
            self.assertEqual((200, resident.json()), (loaded.status_code, loaded.json()))
            # the .seqs replaced (the files in another order): the resident's mapping serves
            replace(os.path.join(self.d, 'swapped.seqs'), seqs)
            loaded = server.post('search', dict(query, in_ram=True))
            self.assertEqual((200, resident.json()), (loaded.status_code, loaded.json()))
            # the annotation and its .seqs replaced by a rebuild of the files in the other
            # order (the same column count): the annotation is not the one loaded at start-up,
            # so the copy loads its .seqs from disk with the load's checks, and answers the
            # same records (with the resident's mapping it answered big.fa's for small.fa's:
            # H4's verification, finding 2)
            replace(os.path.join(self.d, 'rev.column_coord.annodbg'), anno)
            replace(os.path.join(self.d, 'rev.seqs'), seqs)
            loaded = server.post('search', dict(query, in_ram=True))
            self.assertEqual((200, resident.json()), (loaded.status_code, loaded.json()))
            self.assertEqual(['s1'], [r['sample'] for r in loaded.json()[0]['results']])
            # the annotation replaced: the request fails, the route's 400 with the load's message
            replace(os.path.join(self.d, 'three.column_coord.annodbg'), anno)
            failed = server.post('search', dict(query, in_ram=True))
            self.assertEqual(400, failed.status_code, failed.text)
            self.assertEqual(['error'], list(failed.json()))
            self.assertTrue(failed.json()['error'].startswith(
                    f"The record mapping '{seqs}' does not belong to the annotation '{anno}': "
                    f"it has 2 columns, the annotation 3 labels. Rebuild it"), failed.text)
            # the server runs on, its resident pair answering as before
            self.assertIsNone(server.process.poll())
            again = server.post('search', query)
            self.assertEqual((200, resident.json()), (again.status_code, again.json()))
        finally:
            server.stop()
        # a replaced annotation is never offered the resident's mapping: its own .seqs is read
        self.assertNotIn("The record mapping offered for " + anno, server.text())

    def test_the_row_diff_form_checks_an_old_mapping_by_its_count_only(self):
        # the documented limit: the row-diff forms state no coordinates per column, so an old
        # .seqs of the right count but swapped columns is served, its order stated unverified
        graph, anno = self.pair('row_diff_brwt_coord', 'swapped', False)
        # a query of b2's tail: its k-mers sit at coordinates 510-539 of big.fa's column, which
        # the swapped mapping states 30 k-mers for (small.fa's), so mapping them to records
        # throws. The request fails with the route's 400 naming the coordinate, and the server
        # runs on (it aborted: the query's pool ended the process on a task that threw, H4's
        # verification, finding 1); the CLI exits 1 with the message
        records = dict(named_records_of([os.path.join(self.d, 'big.fa')]))
        tail = records['b2'][-60:]
        server = Server(['-i', graph, '-a', anno], os.path.join(self.d, 'rd_swapped.log'))
        try:
            caps = server.get('capabilities').json()
            for query in ({'FASTA': f'>q\n{tail}'},
                          {'FASTA': f'>q\n{tail}', 'query_coords': True}):
                with self.subTest(query=query):
                    ret = server.post('search', query)
                    self.assertEqual(400, ret.status_code, ret.text)
                    self.assertIn('Coordinate', ret.json()['error'])
                    self.assertIn('out of range', ret.json()['error'])
                    self.assertIsNone(server.process.poll())
            self.assertEqual(200, server.get('capabilities').status_code)
        finally:
            server.stop()
        self.assertEqual('count_only', caps['pair_summary']['seqs_order'])
        tail_fa = os.path.join(self.d, 'tail.fa')
        with open(tail_fa, 'w') as f:
            f.write(f'>q\n{tail}\n')
        res = subprocess.run(shlex.split(METAGRAPH) + ['query', '--query-mode', 'coords', '-i', graph,
                                                       '-a', anno, tail_fa],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(1, res.returncode, res.stderr.decode())
        self.assertIn(f"Cannot query '{tail_fa}': Coordinate", res.stderr.decode())


@unittest.skipIf(PROTEIN_MODE, "the fixture is DNA")
class TestLabelsDisjoint(TestingBase):
    """GET /capabilities' labels_disjoint (SPEC-multigraph-envelope.md section 9.1, H4): whether no
    sequence header equals a column name of any pair, the second certificate a header label's
    record selection needs beside headers_unique. Two pairs of small random records, annotated
    by file name in their directories (a column's label is its file name), each with its .seqs
    (stating its labels: seqs_order "verified"): left x1.fa [acc1, acc2], x2.fa [acc3]; right
    y1.fa [acc4, x1.fa], y2.fa [y1.fa]: the header x1.fa is a column of left, the header y1.fa a
    column of right itself; every header is distinct (headers_unique true)."""

    maxDiff = None
    PAIRS = {
        'left': {'x1.fa': ['acc1', 'acc2'], 'x2.fa': ['acc3']},
        'right': {'y1.fa': ['acc4', 'x1.fa'], 'y2.fa': ['y1.fa']},
    }

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        rng = random.Random(20261108)
        cls.pairs = {}
        for name, files in cls.PAIRS.items():
            d = os.path.join(cls.tempdir.name, name)
            os.makedirs(d)
            for fname, headers in files.items():
                with open(os.path.join(d, fname), 'w') as f:
                    for h in headers:
                        f.write(f'>{h}\n{"".join(rng.choice("ACGT") for _ in range(80))}\n')
            fastas = list(files)
            for cmd in (['build', '-p', '1', '--mode', 'basic', '--graph', 'succinct', '-k', '31',
                         '-o', 'graph', *fastas],
                        ['annotate', '-p', '1', '--anno-filename', '-i', 'graph.dbg',
                         '--anno-type', 'column', '--coordinates', '-o', 'anno', *fastas],
                        ['transform_anno', '-p', '1', '--anno-type', 'column_coord',
                         '--coordinates', '-o', 'anno', 'anno.column.annodbg'],
                        ['annotate', '-p', '1', '-i', 'graph.dbg', '--anno-filename',
                         '--index-header-coords', '-o', 'anno', *fastas]):
                res = subprocess.run(shlex.split(METAGRAPH) + cmd, cwd=d,
                                     stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                assert res.returncode == 0, (cmd, res.stdout.decode()[-3000:])
            cls.pairs[name] = (os.path.join(d, 'graph.dbg'),
                               os.path.join(d, 'anno.column_coord.annodbg'))

    def tags(self, name):
        return {'graph': name, 'graph_path': self.pairs[name][0],
                'annotation_path': self.pairs[name][1]}

    def serve(self, names):
        csv = os.path.join(self.tempdir.name, '_'.join(names) + '.csv')
        with open(csv, 'w') as f:
            for name in names:
                f.write(f'{name},{self.pairs[name][0]},{self.pairs[name][1]}\n')
        server = Server([csv], csv + '.log')
        try:
            caps = server.get('capabilities').json()
        finally:
            server.stop()
        self.assertNotIn('labels_disjoint', caps)
        return caps, server.text()

    def test_a_header_equal_to_a_column_name_is_a_collision(self):
        caps, log = self.serve(['left', 'right'])
        summary = caps['graph_summary']
        self.assertEqual({'value': True, 'collisions': 0, 'examples': []},
                         summary['headers_unique'])
        self.assertEqual({'value': False, 'collisions': 2, 'examples': [
            {'header': 'x1.fa', 'records': [dict(self.tags('right'), column='y1.fa')],
             'columns': [self.tags('left')]},
            {'header': 'y1.fa', 'records': [dict(self.tags('right'), column='y2.fa')],
             'columns': [self.tags('right')]},
        ]}, summary['labels_disjoint'])
        self.assertIn("labels_disjoint: false, 2 headers are also column names (one of them: "
                      "'x1.fa'); a label with such a name selects the column where one has it",
                      log)
        self.assertEqual(['verified', 'verified'], [p['seqs_order'] for p in summary['pairs']])

    def test_disjoint_names_are_certified(self):
        caps, log = self.serve(['left'])
        self.assertEqual({'value': True, 'collisions': 0, 'examples': []},
                         caps['graph_summary']['labels_disjoint'])
        self.assertIn('labels_disjoint: true, computed in', log)

    def test_single_graph_server_states_it_at_the_top(self):
        graph, anno = self.pairs['right']
        server = Server(['-i', graph, '-a', anno], os.path.join(self.tempdir.name, 'single.log'))
        try:
            caps = server.get('capabilities').json()
        finally:
            server.stop()
        untagged = {'graph': None, 'graph_path': None, 'annotation_path': None}
        # right alone: its own column y1.fa shadows its record y1.fa
        self.assertEqual({'value': False, 'collisions': 1, 'examples': [
            {'header': 'y1.fa', 'records': [dict(untagged, column='y2.fa')],
             'columns': [untagged]}]}, caps['labels_disjoint'])
        self.assertTrue(caps['headers_unique']['value'])
        self.assertEqual('verified', caps['pair_summary']['seqs_order'])


if __name__ == '__main__':
    unittest.main()
