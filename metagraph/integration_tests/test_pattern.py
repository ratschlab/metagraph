import filecmp
import glob
import json
import math
import os
import random
import re
import shlex
import socket
import subprocess
import tempfile
import time
import unittest

import requests

from base import PROTEIN_MODE, TestingBase, METAGRAPH, TEST_DATA_DIR

"""
End-to-end tests of POST /pattern and `metagraph pattern` (docs/DESIGN-pattern-search.md,
increments 0-2): the counts of the graph contexts of exact DNA and IUPAC patterns, and their
label-free extraction (output.labels "none": k-mer, instance, offset, strand, node and row ids,
no annotation read), checked against PYTHON ORACLES over the very records each index was built
from — never against the engine itself.

A graph context is (strand, k-mer, offset): a k-mer of the graph that contains the oriented
pattern at that offset (§3). On a BASIC graph the k-mers are those of the records as deposited
(each record split at every symbol a DNA4 build does not keep); P is searched on '+', rc(P) on
'-', a palindromic P once as '='. On a CANONICAL graph, and on a PRIMARY one wrapped in
CanonicalDBG, they are those k-mers and their reverse complements, with orientations instead of
strands. The oracle finds every occurrence of the oriented pattern in those sequences (an
overlapping regex) and adds every k-mer that contains it, at each offset.

Fixtures:
  - TestPatternMini: the mini index (build/mini_refseq, scripts/traversal/build_mini_refseq.sh;
    or $METAGRAPH_MINI_REFSEQ), a BASIC succinct graph at k = 31 with a row_diff_brwt_coord
    annotation and its .seqs. It was built without the dummy-edge mask the route requires, so
    the test rebuilds the graph from the same FASTA with --mask-dummy in a temporary directory;
    the rebuilt .dbg is byte-identical, so the original annotation serves it as is. Skipped when
    the mini index is not built.
  - TestPatternSynthetic: random records, BASIC, CANONICAL and PRIMARY graphs at k = 15 (and a
    multi-graph server), always run.
  - TestPatternRegression: /search and /align answer byte for byte as the base binary's
    ($METAGRAPH_BASE_BINARY, e.g. the build of 804731aa) on test_api.py's fixture and requests;
    skipped without it.
"""

MINI_DIR = os.environ.get('METAGRAPH_MINI_REFSEQ', os.path.join(os.getcwd(), 'mini_refseq'))
MINI_K = 31
MINI_GRAPH = 'graph_k31.dbg'
MINI_ANNO = 'annotation.relaxed.relabeled.row_diff_brwt_coord.annodbg'
MINI_SEQS = 'annotation.relaxed.relabeled.seqs'
BASE_BINARY = os.environ.get('METAGRAPH_BASE_BINARY', '')

# the server's defaults (--pattern-* flags; §5.3, the owner's budgets of 2026-10-07): each cap
# is its field's maximum and default, but time_budget_ms defaults to 60 s under a 600 s cap
DEFAULT_CAPS = {'max_contexts': 10000, 'max_anchors': 1000, 'max_steps': 100000000,
                'time_budget_ms': 600000, 'min_information_bits': 24, 'max_patterns': 16}
DEFAULT_TIME_MS = 60000
DEFAULT_FINALIZE_MS = 250

IUPAC = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT',
         'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}
COMPLEMENT = str.maketrans('ACGTRYSWKMBDHVN', 'TGCAYRSWMKVHDBN')


def revcomp(s):
    return s.upper().translate(COMPLEMENT)[::-1]


def information_bits(pattern):
    return sum(math.log2(4 / len(IUPAC[c])) for c in pattern.upper())


def pattern_regex(pattern):
    # a lookahead: every occurrence, overlapping ones included
    return re.compile('(?=(' + ''.join('[' + IUPAC[c] + ']' for c in pattern.upper()) + '))')


def _supports_pattern():
    res = subprocess.run(shlex.split(METAGRAPH) + ['pattern'], stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    out = (res.stdout + res.stderr).decode()
    return 'Usage' in out and 'pattern' in out


def read_fasta(path):
    name, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(seq)
                name, seq = line[1:].split()[0] if len(line) > 1 else '', []
            elif line:
                seq.append(line)
    if name is not None:
        yield name, ''.join(seq)


class Records:
    """
    The records an index was built from, as (source, name, island): each record split at every
    symbol other than A, C, G, T (a DNA4 build keeps no k-mer across one). The oracles below
    scan the islands (BASIC) or the islands and their reverse complements (CANONICAL, PRIMARY).
    """

    def __init__(self, entries):
        self.islands = []
        for source, name, seq in entries:
            for piece in re.split('[^ACGT]+', seq.upper()):
                if piece:
                    self.islands.append((source, name, piece))
        self._cache = {}

    @staticmethod
    def from_fasta(paths):
        return Records([(path, name, seq) for path in paths for name, seq in read_fasta(path)])

    def sequences(self, canonical):
        for _, _, island in self.islands:
            yield island
            if canonical:
                yield revcomp(island)

    @staticmethod
    def oriented(pattern, strands, strand_stated):
        """The oriented patterns searched (§3, "Strand"), with their JSON strand or orientation."""
        p = pattern.upper()
        rc = revcomp(p)
        if rc == p:
            return [('=' if strand_stated else 'palindromic', p)]
        out = []
        if strands in ('both', 'forward'):
            out.append(('+' if strand_stated else 'forward', p))
        if strands in ('both', 'reverse'):
            out.append(('-' if strand_stated else 'reverse', rc))
        return out

    def contexts(self, pattern, k, scope='any_offset', strands='both', canonical=False):
        """{(strand or orientation, k-mer, offset)}: the graph contexts of |pattern| (L <= k)."""
        key = ('contexts', pattern, k, scope, strands, canonical)
        if key in self._cache:
            return self._cache[key]
        L = len(pattern)
        assert L <= k
        out = set()
        for strand, q in self.oriented(pattern, strands, not canonical):
            rx = pattern_regex(q)
            for seq in self.sequences(canonical):
                for m in rx.finditer(seq):
                    i = m.start()
                    # every k-mer [s, s + k) of the island that contains [i, i + L)
                    for s in range(max(0, i + L - k), min(i, len(seq) - k) + 1):
                        p = i - s
                        if scope == 'suffix' and p != k - L:
                            continue
                        out.add((strand, seq[s:s + k], p))
        self._cache[key] = out
        return out

    def anchors(self, pattern, k, strands='both', canonical=False):
        """{(strand, k-mer)}: the k-mers instantiating positions [0, k) of each oriented pattern
        (L > k, §4.2): rc(P)'s anchors are its own first k positions, not P's."""
        assert len(pattern) > k
        out = set()
        for strand, q in self.oriented(pattern, strands, not canonical):
            rx = pattern_regex(q[:k])
            for seq in self.sequences(canonical):
                for m in rx.finditer(seq):
                    if m.start() + k <= len(seq):
                        out.add((strand, seq[m.start():m.start() + k]))
        return out

    def names_with_kmer(self, kmer):
        return {name for _, name, island in self.islands if kmer in island}

    def sources_with_kmer(self, kmer):
        return {source for source, _, island in self.islands if kmer in island}


def free_port():
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(('127.0.0.1', 0))
        return s.getsockname()[1]


class Server:
    """A server_query process on a free port, ready once GET /stats answers 200."""

    def __init__(self, binary, args, log_path, timeout=180):
        os.environ['NO_PROXY'] = '127.0.0.1'
        for _ in range(20):
            self.port = free_port()
            self.log = open(log_path, 'ab')
            self.process = subprocess.Popen(
                shlex.split(binary) + ['server_query'] + args
                + ['--port', str(self.port), '--address', '127.0.0.1', '-p', '2'],
                stdout=self.log, stderr=subprocess.STDOUT)
            if self._wait_ready(timeout):
                return
            self.stop()
        raise AssertionError(f"Couldn't start the server: {binary} {args} (log: {log_path})")

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

    def post(self, route, payload, raw=False):
        return requests.post(url=self.url(route), data=payload if raw else json.dumps(payload),
                             timeout=120)

    def get(self, route):
        return requests.get(url=self.url(route), timeout=60)

    def stop(self):
        if self.process.poll() is None:
            self.process.kill()
            self.process.wait()
        self.log.close()


class PatternChecks:
    """Assertions shared by the fixtures (a mixin of unittest.TestCase)."""

    def pattern(self, server, payload, status=200):
        ret = server.post('pattern', payload)
        self.assertEqual(status, ret.status_code, ret.text)
        return ret.json()

    def assertCount(self, count, value, relation='exact', unit='graph_contexts'):
        self.assertEqual({'value': value, 'relation': relation, 'unit': unit}, count)

    def assertUnknownLabels(self, entry):
        self.assertCount(entry['counts']['labels'], None, 'unknown', 'labels')
        self.assertCount(entry['counts']['occurrences'], None, 'unknown', 'placed_occurrences')

    def assertContextCounts(self, entry, expected, k, length, scope, strand_stated=True):
        """counts.contexts against the oracle's contexts: the total, the suffix offset, every
        offset of the scope (zeros included) and every strand or orientation searched."""
        c = entry['counts']['contexts']
        self.assertCount({x: c[x] for x in ('value', 'relation', 'unit')}, len(expected))
        self.assertCount(c['suffix'], sum(1 for _, _, p in expected if p == k - length))
        offsets = range(k - length + 1) if scope == 'any_offset' else [k - length]
        self.assertEqual({str(p) for p in offsets}, set(c['by_offset']))
        for p in offsets:
            self.assertCount(c['by_offset'][str(p)], sum(1 for _, _, q in expected if q == p))
        by = c['by_strand' if strand_stated else 'by_orientation']
        self.assertNotIn('by_orientation' if strand_stated else 'by_strand', c)
        searched = entry['strands']
        keys = {('both' if s == '=' else s) for s in searched}
        self.assertEqual(keys, set(by))
        for s in searched:
            self.assertCount(by['both' if s == '=' else s],
                             sum(1 for t, _, _ in expected if t == s))
        self.assertUnknownLabels(entry)

    # the answer order's last key (§5.5): an IUPAC pattern and its reverse complement can
    # both match one instance (NA and TN both match TA), two contexts at one (node, offset)
    ORIENTATION_RANK = {'+': 0, '-': 1, '=': 2, 'forward': 0, 'reverse': 1, 'palindromic': 2}

    def assertResults(self, entry, expected, length, strand_stated=True):
        """results: exactly the oracle's contexts, each a real context, in the answer order
        (node, offset, orientation)."""
        results = entry['results']
        key = 'strand' if strand_stated else 'orientation'
        got = [(r[key], r['kmer'], r['offset']) for r in results]
        self.assertEqual(len(got), len(set(got)), 'a context returned twice')
        self.assertEqual(expected, set(got))
        for r in results:
            self.assertEqual({'kmer', 'instance', 'offset', key, 'node', 'row'}, set(r))
            self.assertEqual(r['kmer'][r['offset']:r['offset'] + length], r['instance'])
        order = [(r['node'], r['offset'], self.ORIENTATION_RANK[r[key]]) for r in results]
        self.assertEqual(sorted(order), order)
        self.assertEqual(len(set(order)), len(order), '(node, offset, orientation) is a context')
        self.assertEqual(len(results), entry['returned'])

    def assertCompleteRetrieval(self, entry):
        self.assertTrue(entry['retrieval_complete'])
        self.assertIsNone(entry['withheld'])
        self.assertIsNone(entry['cut'])
        self.assertIsNone(entry['stop'])
        self.assertEqual('full', entry['determinism'])


@unittest.skipIf(PROTEIN_MODE, "pattern search is DNA only")
@unittest.skipUnless(_supports_pattern(), "`metagraph pattern` is not available in this build")
@unittest.skipUnless(os.path.isfile(os.path.join(MINI_DIR, MINI_GRAPH))
                     and os.path.isdir(os.path.join(MINI_DIR, 'fasta')),
                     "the mini index is not built (scripts/traversal/build_mini_refseq.sh)")
class TestPatternMini(PatternChecks, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tempdir = tempfile.TemporaryDirectory()
        d = cls.tempdir.name
        fastas = sorted(glob.glob(os.path.join(MINI_DIR, 'fasta', '*.fa')))
        cls.records = Records.from_fasta(fastas)
        cls.k = MINI_K

        # the mini index carries no .edgemask: rebuilt with it from the same records and the
        # same flags (commands.log), single-threaded so that nothing depends on the order of
        # threads
        cls.graph = os.path.join(d, MINI_GRAPH)
        TestingBase._run_command(
            f'{METAGRAPH} build -p 1 --mode basic --graph succinct --state stat -k {MINI_K} '
            f'--index-ranges 12 --mask-dummy --in-ram -o {cls.graph[:-len(".dbg")]} '
            + ' '.join(fastas), 'Build the masked mini graph')
        assert os.path.isfile(cls.graph[:-len('.dbg')] + '.edgemask')
        if filecmp.cmp(cls.graph, os.path.join(MINI_DIR, MINI_GRAPH), shallow=False):
            # the same graph: its annotation (and record mapping) serves the masked copy
            for f in (MINI_ANNO, MINI_SEQS):
                os.symlink(os.path.join(MINI_DIR, f), os.path.join(d, f))
            cls.anno = os.path.join(d, MINI_ANNO)
            cls.label_of_kmer = cls.records.names_with_kmer
        else:
            TestingBase._run_command(
                f'{METAGRAPH} annotate -p 4 -i {cls.graph} --anno-filename '
                f'-o {d}/annotation ' + ' '.join(fastas), 'Annotate the masked mini graph')
            cls.anno = os.path.join(d, 'annotation.column.annodbg')
            cls.label_of_kmer = cls.records.sources_with_kmer

        cls.server = Server(METAGRAPH, ['-i', cls.graph, '-a', cls.anno],
                            os.path.join(d, 'server.log'))
        cls.stats = cls.server.get('stats').json()
        cls._choose_patterns()

    @classmethod
    def tearDownClass(cls):
        cls.server.stop()
        cls.tempdir.cleanup()

    @classmethod
    def _choose_patterns(cls):
        """Patterns taken from the records at fixed positions (so present), and absent ones."""
        islands = sorted((island for _, _, island in cls.records.islands),
                         key=lambda s: (-len(s), s[:64]))
        a, b = islands[0], islands[1]
        cls.p16 = a[1000:1016]
        cls.p14 = b[5000:5014]
        cls.p11 = a[3000:3011]           # 22 bits: below the default floor in any_offset
        cls.p40 = a[4000:4040]           # longer than k: anchors counted
        # an IUPAC pattern: a 16-mer of the records with two positions widened to two bases
        # and one to N (28 bits)
        w = list(a[2000:2016])
        w[3] = 'R' if w[3] in 'AG' else 'Y'
        w[8] = 'W' if w[8] in 'AT' else 'S'
        w[12] = 'N'
        cls.iupac16 = ''.join(w)
        # the first palindromic 12-mer (24 bits) of the longest island
        cls.pal12 = next(a[i:i + 12] for i in range(len(a) - 12)
                         if revcomp(a[i:i + 12]) == a[i:i + 12])
        rng = random.Random(4711)

        def absent(length):
            while True:
                s = ''.join(rng.choice('ACGT') for _ in range(length))
                if any(s in x or revcomp(s) in x for x in cls.records.sequences(False)):
                    continue
                # a long one without an anchor either: no path can start
                if length > cls.k and cls.records.anchors(s, cls.k):
                    continue
                return s
        cls.absent16 = absent(16)
        cls.absent40 = absent(40)

    def contexts(self, pattern, **kw):
        return self.records.contexts(pattern, self.k, **kw)

    # ------------------------------------------------------------ count mode

    def test_count_exact_dna_any_offset(self):
        out = self.pattern(self.server, {'patterns': [{'id': 'a', 'dna': self.p16},
                                                      {'id': 'b', 'dna': self.p14.lower()},
                                                      {'id': 'z', 'dna': self.absent16}],
                                         'mode': 'count'})
        self.assertEqual(1, out['pattern_contract_version'])
        self.assertEqual('count', out['mode'])
        self.assertIsNone(out['output'])
        self.assertEqual({'k': self.k, 'graph_mode': 'basic', 'alphabet': '$ACGT',
                          'strand_stated': True},
                         {x: out['index'][x] for x in ('k', 'graph_mode', 'alphabet',
                                                       'strand_stated')})
        for entry, p in zip(out['patterns'], (self.p16, self.p14, self.absent16)):
            self.assertEqual(p, entry['pattern'])
            self.assertEqual('dna', entry['kind'])
            self.assertEqual(len(p), entry['length'])
            self.assertEqual(2.0 * len(p), entry['information_bits'])
            self.assertEqual(['+', '-'], entry['strands'])
            self.assertFalse(entry['palindromic'])
            self.assertEqual('any_offset', entry['scope'])
            self.assertEqual('any_offset', entry['absence_scope'])
            self.assertFalse(entry['retrieval_complete'])
            self.assertNotIn('results', entry)
            self.assertIsNone(entry['stop'])
            self.assertEqual('full', entry['determinism'])
            self.assertContextCounts(entry, self.contexts(p), self.k, len(p), 'any_offset')
            self.assertGreaterEqual(entry['work']['steps'], entry['work']['ranges_visited'])
        self.assertEqual(['a', 'b', 'z'], [e['id'] for e in out['patterns']])
        self.assertGreater(out['patterns'][0]['counts']['contexts']['value'], 0)
        self.assertEqual(0, out['patterns'][2]['counts']['contexts']['value'])

    def test_count_suffix_scope(self):
        out = self.pattern(self.server, {'patterns': [{'dna': self.p16}, {'dna': self.p14}],
                                         'mode': 'count', 'scope': 'suffix'})
        for entry, p in zip(out['patterns'], (self.p16, self.p14)):
            self.assertEqual('suffix', entry['scope'])
            self.assertEqual('suffix_only', entry['absence_scope'])
            self.assertIsNone(entry['id'])
            expected = self.contexts(p, scope='suffix')
            self.assertEqual(expected, {c for c in self.contexts(p) if c[2] == self.k - len(p)})
            self.assertContextCounts(entry, expected, self.k, len(p), 'suffix')

    def test_count_iupac_and_palindrome(self):
        out = self.pattern(self.server, {'patterns': [{'iupac': self.iupac16},
                                                      {'dna': self.pal12}],
                                         'mode': 'count'})
        iupac, pal = out['patterns']
        self.assertEqual('iupac', iupac['kind'])
        self.assertEqual(information_bits(self.iupac16), iupac['information_bits'])
        self.assertContextCounts(iupac, self.contexts(self.iupac16), self.k, 16, 'any_offset')
        self.assertTrue(pal['palindromic'])
        self.assertEqual(['='], pal['strands'])
        self.assertEqual({'both'}, set(pal['counts']['contexts']['by_strand']))
        self.assertContextCounts(pal, self.contexts(self.pal12), self.k, 12, 'any_offset')
        self.assertGreater(pal['counts']['contexts']['value'], 0)

    def test_count_strands(self):
        for strands, symbol in (('forward', '+'), ('reverse', '-')):
            out = self.pattern(self.server, {'patterns': [{'dna': self.p14}], 'mode': 'count',
                                             'strands': strands})
            entry = out['patterns'][0]
            self.assertEqual([symbol], entry['strands'])
            expected = self.contexts(self.p14, strands=strands)
            self.assertEqual(expected, {c for c in self.contexts(self.p14) if c[0] == symbol})
            self.assertContextCounts(entry, expected, self.k, 14, 'any_offset')

    def test_information_floor(self):
        out = self.pattern(self.server, {'patterns': [{'id': 'low', 'dna': self.p11},
                                                      {'id': 'ok', 'dna': self.p16},
                                                      {'id': 'n', 'iupac': 'NNNNNNNNACGTACGTAC'}],
                                         'mode': 'count'})
        low, ok, n = out['patterns']
        for entry, bits in ((low, 22.0), (n, 20.0)):
            self.assertEqual('information_below_floor', entry['error']['code'])
            self.assertEqual(bits, entry['information_bits'])
            self.assertNotIn('counts', entry)
            self.assertNotIn('results', entry)
        self.assertEqual(self.p11, low['pattern'])
        self.assertEqual(11, low['length'])
        self.assertIn('counts', ok)
        # an exact pattern in suffix scope is admitted however short (§5.3)
        out = self.pattern(self.server, {'patterns': [{'dna': self.p11}], 'mode': 'count',
                                         'scope': 'suffix'})
        entry = out['patterns'][0]
        self.assertNotIn('error', entry)
        self.assertContextCounts(entry, self.contexts(self.p11, scope='suffix'), self.k, 11,
                                 'suffix')

    def test_bad_alphabet(self):
        out = self.pattern(self.server, {'patterns': [{'id': 'u', 'dna': 'ACGUACGTACGTACGT'},
                                                      {'id': 'r', 'dna': 'ACGRACGTACGTACGT'},
                                                      {'id': 'x', 'iupac': 'ACGTX'},
                                                      {'id': 'e', 'dna': ''},
                                                      {'id': 'ok', 'dna': self.p16}],
                                         'mode': 'count'})
        for entry in out['patterns'][:4]:
            self.assertEqual('bad_alphabet', entry['error']['code'], entry)
            self.assertNotIn('pattern', entry)
            self.assertNotIn('counts', entry)
        self.assertEqual(['dna', 'dna', 'iupac', 'dna'], [e['kind'] for e in out['patterns'][:4]])
        self.assertContextCounts(out['patterns'][4], self.contexts(self.p16), self.k, 16,
                                 'any_offset')

    def test_refusals(self):
        p = [{'dna': self.p16}]
        cases = [
            ({'patterns': p, 'bogus': 1}, 'invalid_request', "unknown field 'bogus'"),
            ({'patterns': [{'dna': self.p16, 'name': 'x'}]}, 'invalid_request',
             "unknown field 'name'"),
            ({'patterns': p, 'output': {'labels': 'none', 'format': 'x'}}, 'invalid_request',
             "unknown field 'format'"),
            ({'patterns': [{'dna': self.p16, 'iupac': self.p16}]}, 'invalid_request', 'exactly one'),
            ({'patterns': [{'id': 1, 'dna': self.p16}]}, 'invalid_request', 'id'),
            ({'patterns': []}, 'invalid_request', 'patterns'),
            ({'patterns': p * 17}, 'invalid_request', 'patterns'),
            ({}, 'invalid_request', 'patterns'),
            ([], 'invalid_request', 'object'),
            ({'patterns': p, 'mode': 'all'}, 'invalid_request', 'mode'),
            ({'patterns': p, 'scope': 'long'}, 'invalid_request', 'scope'),
            ({'patterns': p, 'strands': '+'}, 'invalid_request', 'strands'),
            ({'patterns': p, 'max_steps': 0}, 'invalid_request', 'max_steps'),
            ({'patterns': p, 'max_contexts': -1}, 'invalid_request', 'max_contexts'),
            ({'patterns': p, 'time_budget_ms': DEFAULT_FINALIZE_MS}, 'invalid_request',
             'time_budget_ms'),
            ({'patterns': p, 'stop_at_threshold': 'yes'}, 'invalid_request', 'stop_at_threshold'),
            ({'patterns': p, 'output': {'labels': 'all'}}, 'later_increment', 'labels'),
            ({'patterns': p, 'mode': 'count', 'output': {'labels': 'predicate_only'}},
             'later_increment', 'labels'),
            ({'patterns': p, 'output': {'occurrences': True}}, 'later_increment', 'occurrences'),
            ({'patterns': p, 'output': {'paths': True}}, 'later_increment', 'paths'),
            ({'patterns': p, 'predicate': {'any': ['562']}}, 'later_increment', 'predicate'),
            ({'patterns': p, 'max_labels': None}, 'later_increment', 'max_labels'),
            ({'patterns': p, 'graphs': ['x']}, 'later_increment', 'graphs'),
            ({'patterns': [{'protein': 'MKV'}]}, 'later_increment', 'protein'),
            ({'patterns': p, 'in_ram': False}, 'resident_only', 'in_ram'),
        ]
        for payload, code, words in cases:
            ret = self.server.post('pattern', payload)
            self.assertEqual(400, ret.status_code, (payload, ret.text))
            body = ret.json()
            self.assertEqual({'error', 'code'}, set(body), payload)
            self.assertEqual(code, body['code'], (payload, body))
            self.assertIn(words, body['error'], payload)
        ret = self.server.post('pattern', '{"patterns": [', raw=True)
        self.assertEqual(400, ret.status_code)
        self.assertEqual('invalid_request', ret.json()['code'])
        # false projections and occurrences are accepted and change nothing
        out = self.pattern(self.server, {'patterns': p, 'mode': 'count',
                                         'output': {'labels': 'none', 'occurrences': False,
                                                    'paths': False}})
        self.assertContextCounts(out['patterns'][0], self.contexts(self.p16), self.k, 16,
                                 'any_offset')

    def test_max_steps_and_clamps(self):
        out = self.pattern(self.server, {'patterns': [{'dna': self.p16}, {'dna': self.p14}],
                                         'mode': 'count', 'max_steps': 1})
        first, second = out['patterns']
        c = first['counts']['contexts']
        self.assertEqual('at_least', c['relation'])
        self.assertLessEqual(c['value'], len(self.contexts(self.p16)))
        self.assertEqual('max_steps', first['stop']['reason'])
        self.assertIn(first['stop']['phase'], ('discovery', 'mask_scan'))
        self.assertLessEqual(first['work']['steps'], 1)
        # the budget is the request's: the next pattern finds it spent
        self.assertEqual({'phase': 'discovery', 'reason': 'max_steps'}, second['stop'])
        self.assertCount({x: second['counts']['contexts'][x]
                          for x in ('value', 'relation', 'unit')}, None, 'unknown')
        self.assertEqual(0, second['work']['steps'])
        self.assertEqual('full', first['determinism'])
        self.assertEqual(1, out['limits']['max_steps'])
        self.assertEqual([], out['limits']['clamped'])
        # the defaults: every cap but the time budget, which defaults below its cap
        self.assertEqual(DEFAULT_TIME_MS, out['limits']['time_budget_ms'])
        self.assertEqual(DEFAULT_CAPS['max_contexts'], out['limits']['max_contexts'])

        out = self.pattern(self.server, {'patterns': [{'dna': self.p16}], 'mode': 'count',
                                         'max_steps': 10 ** 9, 'max_contexts': 10 ** 6,
                                         'time_budget_ms': 10 ** 6})
        limits = out['limits']
        self.assertEqual(DEFAULT_CAPS['max_steps'], limits['max_steps'])
        self.assertEqual(DEFAULT_CAPS['max_contexts'], limits['max_contexts'])
        self.assertEqual(DEFAULT_CAPS['time_budget_ms'], limits['time_budget_ms'])
        self.assertEqual(DEFAULT_FINALIZE_MS, limits['finalize_reserve_ms'])
        self.assertEqual(
            sorted([{'field': 'max_steps', 'requested': 10 ** 9,
                     'effective': DEFAULT_CAPS['max_steps']},
                    {'field': 'max_contexts', 'requested': 10 ** 6,
                     'effective': DEFAULT_CAPS['max_contexts']},
                    {'field': 'time_budget_ms', 'requested': 10 ** 6,
                     'effective': DEFAULT_CAPS['time_budget_ms']}],
                   key=lambda x: x['field']),
            sorted(limits['clamped'], key=lambda x: x['field']))
        self.assertEqual('exact', out['patterns'][0]['counts']['contexts']['relation'])
        # a budget under the cap is kept as asked
        out = self.pattern(self.server, {'patterns': [{'dna': self.p16}], 'mode': 'count',
                                         'time_budget_ms': 1234.5})
        self.assertEqual(1234.5, out['limits']['time_budget_ms'])
        self.assertEqual([], out['limits']['clamped'])

    def test_capabilities(self):
        caps = self.server.get('capabilities').json()
        self.assertIn('pattern', caps['features'])
        self.assertEqual('POST /pattern', caps['routes']['pattern'])
        p = caps['pattern']
        expected = {
            'pattern_contract_version': 1, 'available': True, 'unavailable_reason': None,
            'modes': ['count', 'all_or_count', 'partial'], 'default_mode': 'all_or_count',
            'projections': ['none'], 'default_projection': 'none',
            'projections_later_increment': ['all', 'predicate_only'],
            'kinds': ['dna', 'iupac'], 'scopes': ['suffix', 'any_offset'],
            'default_scope': 'any_offset', 'long_patterns': 'anchors_counted',
            'strands': ['both', 'forward', 'reverse'], 'graph_mode': 'basic', 'k': self.k,
            'alphabet': '$ACGT', 'strand_stated': True, 'mask': 'file',
            'graph_cleaned': 'unknown', 'records_shorter_than_k': 'not_indexed',
            'resident_only': True, 'caps': DEFAULT_CAPS,
            'default_time_budget_ms': DEFAULT_TIME_MS,
            'finalize_reserve_ms': DEFAULT_FINALIZE_MS,
        }
        self.assertEqual(expected, {x: p[x] for x in expected})
        self.assertEqual(['any_offset'], p['scopes_by_graph_mode']['primary'])
        # the same block on the probe a service reads (§7.3)
        probe = self.server.get('traverse/capabilities').json()
        self.assertEqual(p, probe['pattern'])
        if self.anno.endswith(MINI_ANNO):
            # coordinates and record mapping, a budgeted row-diff annotation
            self.assertEqual(('record', 'record_verified', 'budgeted'),
                             (p['placement'], p['support'], p['annotation']))

    # ------------------------------------------------------------ the label-free retrieval

    def test_retrieval_all_or_count(self):
        patterns = [self.p16, self.iupac16, self.pal12, self.absent16]
        # mode and projection omitted: all_or_count with labels "none" in this increment
        request = {'patterns': [{'iupac' if 'N' in p else 'dna': p} for p in patterns]}
        out = self.pattern(self.server, request)
        self.assertEqual('all_or_count', out['mode'])
        self.assertEqual({'labels': 'none'}, out['output'])
        for entry, p in zip(out['patterns'], patterns):
            expected = self.contexts(p)
            self.assertContextCounts(entry, expected, self.k, len(p), 'any_offset')
            self.assertCompleteRetrieval(entry)
            self.assertResults(entry, expected, len(p))
            self.assertEqual(len(expected), entry['returned'])
            for r in entry['results']:
                # the instance is P on +/=, rc(P) on -
                q = revcomp(p) if r['strand'] == '-' else p
                self.assertTrue(pattern_regex(q).match(r['instance']), (p, r))
                self.assertEqual(r['node'] - 1, r['row'])
        self.assertEqual([], out['patterns'][3]['results'])
        # deterministic: the same request answers the same contexts in the same order
        again = self.pattern(self.server, request)
        for a, b in zip(out['patterns'], again['patterns']):
            self.assertEqual(a['results'], b['results'])
            self.assertEqual(a['counts'], b['counts'])

    def test_node_and_row_ids(self):
        """Node ids are BOSS edge indices and rows are node - 1 (BASIC): checked without the
        engine's word — no route exposes node ids — by (a) the BOSS order of the returned
        k-mers, (b) the k-mer looked up again as an exact pattern of length k, (c) /search of
        the k-mer, which finds it in the records the oracle names, (d) the annotation's size."""
        out = self.pattern(self.server, {'patterns': [{'dna': self.p14}],
                                         'output': {'labels': 'none'}})
        results = out['patterns'][0]['results']
        self.assertGreater(len(results), 1)
        # (a) BOSS sorts edges by their source node's (k-1)-mer read backwards, then by the
        # last symbol: node order must be that order of the k-mers
        def boss_key(kmer):
            return (kmer[:-1][::-1], kmer[-1])
        by_node = {}
        for r in results:
            self.assertEqual(by_node.setdefault(r['node'], r['kmer']), r['kmer'])
        nodes = sorted(by_node)
        keys = [boss_key(by_node[n]) for n in nodes]
        self.assertEqual(sorted(keys), keys)
        self.assertEqual(len(set(keys)), len(keys))
        # (d) every id inside the index
        objects = self.stats['annotation']['objects']
        for r in results:
            self.assertTrue(1 <= r['node'] <= objects)
            self.assertEqual(r['node'] - 1, r['row'])
        for node in nodes[:6]:
            kmer = by_node[node]
            # (b) the k-mer as an exact pattern of length k: one context, the same node
            again = self.pattern(self.server, {'patterns': [{'dna': kmer}], 'scope': 'suffix',
                                               'strands': 'forward'})['patterns'][0]
            self.assertEqual([{'kmer': kmer, 'instance': kmer, 'offset': 0, 'strand': '+',
                               'node': node, 'row': node - 1}], again['results'])
            # (c) /search finds the k-mer where the records have it
            ret = self.server.post('search', {'FASTA': f'>q\n{kmer}\n',
                                              'discovery_fraction': 0.0, 'top_labels': 1000})
            self.assertEqual(200, ret.status_code, ret.text)
            labels = {x['sample'] for x in ret.json()[0]['results']}
            self.assertEqual(self.label_of_kmer(kmer), labels)
            self.assertTrue(labels)

    def test_all_or_count_threshold(self):
        expected = self.contexts(self.p14)
        n = len(expected)
        self.assertGreater(n, 2)
        out = self.pattern(self.server, {'patterns': [{'dna': self.p14}],
                                         'max_contexts': n - 1})
        entry = out['patterns'][0]
        # the count is exact; the contexts are withheld, all or nothing
        self.assertContextCounts(entry, expected, self.k, 14, 'any_offset')
        self.assertEqual({'reason': 'count_above_threshold'}, entry['withheld'])
        self.assertEqual([], entry['results'])
        self.assertEqual(0, entry['returned'])
        self.assertFalse(entry['retrieval_complete'])
        self.assertIsNone(entry['cut'])
        self.assertIsNone(entry['stop'])
        # at the threshold itself everything is returned
        out = self.pattern(self.server, {'patterns': [{'dna': self.p14}], 'max_contexts': n})
        self.assertCompleteRetrieval(out['patterns'][0])
        self.assertResults(out['patterns'][0], expected, 14)

    def test_partial(self):
        expected = self.contexts(self.p14)
        full = self.pattern(self.server, {'patterns': [{'dna': self.p14}]})['patterns'][0]
        self.assertResults(full, expected, 14)
        n = min(5, len(expected) - 1)
        out = self.pattern(self.server, {'patterns': [{'dna': self.p14}], 'mode': 'partial',
                                         'max_contexts': n})
        entry = out['patterns'][0]
        self.assertEqual('partial', entry['mode'])
        self.assertContextCounts(entry, expected, self.k, 14, 'any_offset')
        # the first n contexts in the deterministic order, the cut stated
        self.assertEqual(full['results'][:n], entry['results'])
        self.assertEqual({'reason': 'max_contexts'}, entry['cut'])
        self.assertIsNone(entry['withheld'])
        self.assertFalse(entry['retrieval_complete'])
        self.assertEqual(n, entry['returned'])
        # a cap above the count: complete
        out = self.pattern(self.server, {'patterns': [{'dna': self.p14}], 'mode': 'partial'})
        self.assertCompleteRetrieval(out['patterns'][0])
        self.assertEqual(full['results'], out['patterns'][0]['results'])

    def test_stops_in_retrieval(self):
        expected = self.contexts(self.p14)
        # stop_at_threshold: the count stops past max_contexts (at_least), nothing released
        out = self.pattern(self.server, {'patterns': [{'dna': self.p14}], 'max_contexts': 2,
                                         'stop_at_threshold': True})
        entry = out['patterns'][0]
        self.assertEqual('at_least', entry['counts']['contexts']['relation'])
        self.assertEqual({'phase': 'discovery', 'reason': 'max_contexts'}, entry['stop'])
        self.assertEqual({'reason': 'threshold_crossed'}, entry['withheld'])
        self.assertEqual([], entry['results'])
        self.assertTrue(out['limits']['stop_at_threshold'])
        # max_steps in all_or_count: withheld; in partial: what was built, the cut stated
        out = self.pattern(self.server, {'patterns': [{'dna': self.p14}], 'max_steps': 1})
        entry = out['patterns'][0]
        self.assertEqual({'reason': 'discovery_budget'}, entry['withheld'])
        self.assertEqual([], entry['results'])
        self.assertFalse(entry['retrieval_complete'])
        out = self.pattern(self.server, {'patterns': [{'dna': self.p14}], 'max_steps': 1,
                                         'mode': 'partial'})
        entry = out['patterns'][0]
        self.assertEqual({'reason': 'max_steps'}, entry['cut'])
        self.assertIsNone(entry['withheld'])
        self.assertFalse(entry['retrieval_complete'])
        got = {(r['strand'], r['kmer'], r['offset']) for r in entry['results']}
        self.assertLessEqual(got, expected)
        self.assertLessEqual(len(got), entry['counts']['contexts']['value'])

    def test_long_pattern(self):
        out = self.pattern(self.server, {'patterns': [{'dna': self.p40},
                                                      {'dna': self.absent40}]})
        entry, absent = out['patterns']
        self.assertEqual('long', entry['scope'])
        self.assertEqual('long', entry['absence_scope'])
        self.assertEqual(80.0, entry['information_bits'])
        self.assertEqual(62.0, entry['anchor_information_bits'])
        self.assertIn('paths_later_increment', entry['notes'])
        self.assertNotIn('contexts', entry['counts'])
        anchors = self.records.anchors(self.p40, self.k)
        a = entry['counts']['anchors']
        self.assertCount({x: a[x] for x in ('value', 'relation', 'unit')}, len(anchors),
                         unit='anchors')
        for s in ('+', '-'):
            self.assertCount(a['by_strand'][s], sum(1 for t, _ in anchors if t == s),
                             unit='anchors')
        self.assertGreater(len(anchors), 0)
        self.assertCount(entry['counts']['paths'], None, 'unknown', 'paths')
        self.assertEqual({'reason': 'paths_later_increment'}, entry['withheld'])
        self.assertEqual([], entry['results'])
        self.assertFalse(entry['retrieval_complete'])
        self.assertUnknownLabels(entry)
        # no anchor, no path: exact zeros, and that absence is complete
        self.assertCount({x: absent['counts']['anchors'][x] for x in ('value', 'relation', 'unit')},
                         0, unit='anchors')
        self.assertCount(absent['counts']['paths'], 0, unit='paths')
        self.assertTrue(absent['retrieval_complete'])
        self.assertIsNone(absent['withheld'])

    # ------------------------------------------------------------ the CLI

    def test_cli_answers_as_the_server(self):
        request = {'patterns': [{'id': 'a', 'dna': self.p16}, {'iupac': self.iupac16},
                                {'dna': self.p11}, {'dna': self.p40}],
                   'mode': 'partial', 'max_contexts': 7}
        server_out = self.pattern(self.server, request)
        path = os.path.join(self.tempdir.name, 'request.json')
        with open(path, 'w') as f:
            json.dump(request, f)
        res = subprocess.run(shlex.split(METAGRAPH) + ['pattern', '--json', '-i', self.graph,
                                                       '-a', self.anno, path],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(0, res.returncode, res.stderr.decode())
        cli_out = json.loads(res.stdout)

        def strip(out):
            out = json.loads(json.dumps(out))
            out.pop('timing')
            for e in out['patterns']:
                e.pop('timing', None)
            return out
        self.assertEqual(strip(server_out), strip(cli_out))

        # a refusal: its body, exit status 1
        with open(path, 'w') as f:
            json.dump({'patterns': [{'dna': self.p16}], 'bogus': True}, f)
        res = subprocess.run(shlex.split(METAGRAPH) + ['pattern', '--json', '-i', self.graph,
                                                       '-a', self.anno, path],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(1, res.returncode)
        self.assertEqual('invalid_request', json.loads(res.stdout)['code'])

    def test_cli_mask_required(self):
        """The mini index as built (no .edgemask): refused, since every dummy edge would count."""
        path = os.path.join(self.tempdir.name, 'request_mask.json')
        with open(path, 'w') as f:
            json.dump({'patterns': [{'dna': self.p16}], 'mode': 'count'}, f)
        res = subprocess.run(shlex.split(METAGRAPH) + [
                                 'pattern', '--json', '-i', os.path.join(MINI_DIR, MINI_GRAPH),
                                 '-a', os.path.join(MINI_DIR, MINI_ANNO), path],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(1, res.returncode, res.stderr.decode())
        self.assertEqual('mask_required', json.loads(res.stdout)['code'])


@unittest.skipIf(PROTEIN_MODE, "pattern search is DNA only")
@unittest.skipUnless(_supports_pattern(), "`metagraph pattern` is not available in this build")
class TestPatternSynthetic(PatternChecks, TestingBase):
    """Random records at k = 15, served as BASIC, CANONICAL and PRIMARY graphs, with the
    information floor lowered to 8 bits so that short motifs with many contexts (and k-mers
    holding a motif twice) are compared context by context."""

    K = 15
    FLOOR = 8
    # exact, palindromic (GAATTC, ACGT), IUPAC, and an IUPAC palindrome (WSSWWSSW) some of
    # whose instances are palindromes and others not; every one >= 8 bits
    PATTERNS = ['ACGAC', 'GAATTC', 'TTAGGA', 'RCGNYA', 'ACGT', 'WSSWWSSW']

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        d = cls.tempdir.name
        rng = random.Random(15)
        records = []
        for i in range(5):
            seq = ''.join(rng.choice('ACGT') for _ in range(3000))
            if i == 2:
                # an N splits the record: no k-mer spans it
                seq = seq[:1500] + 'NN' + seq[1502:]
            records.append((f'rec{i}', seq))
        # a k-mer holding a motif twice, and a palindrome
        records.append(('rec_motifs', 'TTTACGACGACTTTGAATTCGAATTCAAA' * 3))
        cls.fasta = os.path.join(d, 'records.fa')
        with open(cls.fasta, 'w') as f:
            for name, seq in records:
                f.write(f'>{name}\n{seq}\n')
        cls.records = Records.from_fasta([cls.fasta])

        cls.servers = {}
        for mode in ('basic', 'canonical', 'primary'):
            graph = os.path.join(d, f'graph_{mode}.dbg')
            cls._build_graph(cls.fasta, graph, cls.K, 'succinct', mode=mode,
                             extra_params='--mask-dummy --in-ram')
            cls._annotate_graph(cls.fasta, graph, os.path.join(d, f'anno_{mode}'), 'column')
            anno = os.path.join(d, f'anno_{mode}.column.annodbg')
            cls.servers[mode] = Server(METAGRAPH, ['-i', graph, '-a', anno,
                                                   '--pattern-min-information-bits',
                                                   str(cls.FLOOR)],
                                       os.path.join(d, f'server_{mode}.log'))
        cls.graph_basic = os.path.join(d, 'graph_basic.dbg')
        cls.anno_basic = os.path.join(d, 'anno_basic.column.annodbg')

    @classmethod
    def tearDownClass(cls):
        for server in cls.servers.values():
            server.stop()
        super().tearDownClass()

    def test_counts_and_contexts_in_every_graph_mode(self):
        for mode, server in self.servers.items():
            stated = mode == 'basic'
            for scope in ('any_offset', 'suffix'):
                request = {'patterns': [{'iupac': p} for p in self.PATTERNS], 'scope': scope,
                           'output': {'labels': 'none'}}
                out = self.pattern(server, request)
                self.assertEqual(mode, out['index']['graph_mode'])
                self.assertEqual(stated, out['index']['strand_stated'])
                for entry, p in zip(out['patterns'], self.PATTERNS):
                    with self.subTest(mode=mode, scope=scope, pattern=p):
                        if mode == 'primary' and scope == 'suffix':
                            # a virtual suffix is a stored prefix (§4.1): refused per pattern
                            self.assertEqual('scope_unsupported', entry['error']['code'])
                            continue
                        expected = self.records.contexts(p, self.K, scope=scope,
                                                         canonical=not stated)
                        self.assertGreater(len(expected), 0)
                        self.assertContextCounts(entry, expected, self.K, len(p), scope,
                                                 strand_stated=stated)
                        self.assertCompleteRetrieval(entry)
                        self.assertResults(entry, expected, len(p), strand_stated=stated)
                        if not stated:
                            self.assertIn('strand_unknown_canonical', entry['notes'])
                            by = entry['counts']['contexts']['by_orientation']
                            if 'forward' in by and scope == 'any_offset':
                                # both orientations of every k-mer are in the graph: (forward,
                                # x, p) has its (reverse, rc(x), k - L - p), an offset that
                                # only any_offset scope covers
                                self.assertEqual(by['forward'], by['reverse'])

    def test_rows(self):
        """The row is the annotation row of the context's k-mer: node - 1 on BASIC; on the
        other modes a k-mer and its reverse complement share one row (one annotated key)."""
        for mode, server in self.servers.items():
            out = self.pattern(server, {'patterns': [{'dna': 'GAATTC'}, {'dna': 'ACGAC'}]})
            rows = {}
            for entry in out['patterns']:
                for r in entry['results']:
                    self.assertIsNotNone(r['row'])
                    if mode == 'basic':
                        self.assertEqual(r['node'] - 1, r['row'])
                    rows.setdefault(r['kmer'], set()).add(r['row'])
            for kmer, row in rows.items():
                self.assertEqual(1, len(row))
                if mode != 'basic' and revcomp(kmer) in rows:
                    self.assertEqual(row, rows[revcomp(kmer)], (mode, kmer))

    def test_partial_order(self):
        server = self.servers['basic']
        full = self.pattern(server, {'patterns': [{'dna': 'ACGAC'}]})['patterns'][0]
        n = len(full['results'])
        self.assertGreater(n, 7)
        part = self.pattern(server, {'patterns': [{'dna': 'ACGAC'}], 'mode': 'partial',
                                     'max_contexts': 7})['patterns'][0]
        self.assertEqual(full['results'][:7], part['results'])
        self.assertEqual({'reason': 'max_contexts'}, part['cut'])

    def test_capabilities_per_graph_mode(self):
        for mode, server in self.servers.items():
            p = server.get('capabilities').json()['pattern']
            self.assertTrue(p['available'], (mode, p))
            self.assertEqual(mode, p['graph_mode'])
            self.assertEqual(self.K, p['k'])
            self.assertEqual(self.FLOOR, p['caps']['min_information_bits'])
            self.assertEqual(['any_offset'] if mode == 'primary' else ['suffix', 'any_offset'],
                             p['scopes'])
            self.assertEqual(mode == 'basic', p['strand_stated'])
            self.assertEqual('none' if mode == 'basic' else 'none_canonical', p['placement'])

    def test_multi_graph_server(self):
        d = self.tempdir.name
        csv = os.path.join(d, 'graphs.csv')
        with open(csv, 'w') as f:
            f.write(f'G1,{self.graph_basic},{self.anno_basic}\n')
        server = Server(METAGRAPH, [csv], os.path.join(d, 'server_multi.log'))
        try:
            ret = server.post('pattern', {'patterns': [{'dna': 'ACGAC'}], 'mode': 'count'})
            self.assertEqual(400, ret.status_code)
            self.assertEqual({'error': 'pattern: multi-graph servers in a later increment',
                              'code': 'later_increment'}, ret.json())
            caps = server.get('capabilities').json()
            self.assertNotIn('pattern', caps['features'])
            self.assertNotIn('pattern', caps['routes'])
            multi_block = {'pattern_contract_version': 1, 'available': False,
                           'unavailable_reason': 'multi_graph_later_increment'}
            self.assertEqual(multi_block, caps['pattern'])
            probe = server.get('traverse/capabilities?graph=G1').json()
            self.assertEqual(multi_block, probe['pattern'])
        finally:
            server.stop()

    def test_cli(self):
        request = {'patterns': [{'iupac': p} for p in self.PATTERNS], 'mode': 'count'}
        server_out = self.pattern(self.servers['basic'], request)
        path = os.path.join(self.tempdir.name, 'request.json')
        with open(path, 'w') as f:
            json.dump(request, f)
        res = subprocess.run(shlex.split(METAGRAPH) + [
                                 'pattern', '--json', '--pattern-min-information-bits',
                                 str(self.FLOOR), '-i', self.graph_basic, '-a', self.anno_basic,
                                 path],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(0, res.returncode, res.stderr.decode())
        cli_out = json.loads(res.stdout)
        for a, b in zip(server_out['patterns'], cli_out['patterns']):
            self.assertEqual(a['counts'], b['counts'])
            self.assertEqual(a['work'], b['work'])


@unittest.skipIf(PROTEIN_MODE, "the request panel is DNA")
@unittest.skipUnless(BASE_BINARY and os.path.isfile(BASE_BINARY),
                     "$METAGRAPH_BASE_BINARY (the binary before the pattern search) is not set")
class TestPatternRegression(TestingBase):
    """/search and /align answer byte for byte as before the pattern search (§11), on
    test_api.py's fixture (transcripts_100.fa, k = 6) and requests; /capabilities only gains."""

    SEQ = 'CCTCTGTGGAATCCAATCTGTCTTCCATCCTGCGTGGCCGAGGG'
    SEQ2 = 'AATAAAGGTGTGAGATAACCCCAGCGGTGCCAGGATCCGTGCA'

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        d = cls.tempdir.name
        fasta = TEST_DATA_DIR + '/transcripts_100.fa'
        graph = os.path.join(d, 'graph.dbg')
        cls._build_graph(fasta, graph, 6, 'succinct', mode='basic')
        cls._annotate_graph(fasta, graph, os.path.join(d, 'annotation'), 'column')
        args = ['-i', graph, '-a', os.path.join(d, 'annotation.column.annodbg')]
        cls.new = Server(METAGRAPH, args, os.path.join(d, 'server_new.log'))
        cls.base = Server(BASE_BINARY, args, os.path.join(d, 'server_base.log'))

    @classmethod
    def tearDownClass(cls):
        cls.new.stop()
        cls.base.stop()
        super().tearDownClass()

    def test_search_and_align_unchanged(self):
        two = f'>q1\n{self.SEQ2}\n>q2\n{self.SEQ}'
        panel = [
            ('search', {'FASTA': f'>query\n{self.SEQ}', 'top_labels': 5, 'min_exact_match': 0.1}),
            ('search', {'FASTA': two, 'discovery_fraction': 0.1}),
            ('search', {'FASTA': two, 'discovery_fraction': 0.0, 'with_signature': True}),
            ('search', {'FASTA': two, 'align': True, 'discovery_fraction': 0.0}),
            ('search', {'FASTA': '>q\nSEQUENCE_NOT_IN_GRAPH', 'discovery_fraction': 0.01,
                        'top_labels': 1}),
            ('search', {'FASTA': f'>query\n{self.SEQ}', 'query_coords': True}),
            ('search', {'FASTA': f'>query\n{self.SEQ}', 'abundance_sum': True}),
            ('search', '{"FASTA": ">query\\nAATAAAGG", "discovery_fraction": 0.1,'),
            ('align', {'FASTA': '>query0\nTCGATCGA\n>query1\nTCGATCGA', 'min_exact_match': 0}),
            ('align', {'FASTA': '>query\nNNNNNNNN', 'min_exact_match': 0}),
            ('align', {'FASTA': f'>q\n{self.SEQ}', 'max_alternative_alignments': 3}),
            ('align', {'FASTA': '>\nTCGATCGA', 'min_exact_match': 0}),
        ]
        for route, payload in panel:
            raw = isinstance(payload, str)
            a = self.base.post(route, payload, raw=raw)
            b = self.new.post(route, payload, raw=raw)
            self.assertEqual((a.status_code, a.content), (b.status_code, b.content),
                             (route, payload))
        for route in ('column_labels', 'stats'):
            a, b = self.base.get(route), self.new.get(route)
            self.assertEqual((a.status_code, a.content), (b.status_code, b.content), route)

    def test_capabilities_only_gain(self):
        def without_instance(value):
            # each process names itself (server_instance, also inside attempts)
            if isinstance(value, dict):
                return {k: without_instance(v) for k, v in value.items()
                        if k != 'server_instance'}
            if isinstance(value, list):
                return [without_instance(v) for v in value]
            return value
        a = without_instance(self.base.get('capabilities').json())
        b = without_instance(self.new.get('capabilities').json())
        self.assertIn('pattern', b)
        self.assertEqual(a['features'] + ['pattern'], b['features'])
        self.assertEqual(dict(a['routes'], pattern='POST /pattern'), b['routes'])
        for key in a:
            if key not in ('features', 'routes'):
                self.assertEqual(a[key], b[key], key)
        self.assertEqual(set(a) | {'pattern'}, set(b))
        # the probe gains the same block and nothing else
        a = without_instance(self.base.get('traverse/capabilities').json())
        b = without_instance(self.new.get('traverse/capabilities').json())
        self.assertEqual(b['pattern'], self.new.get('capabilities').json()['pattern'])
        b.pop('pattern')
        self.assertEqual(a, b)


if __name__ == '__main__':
    unittest.main()
