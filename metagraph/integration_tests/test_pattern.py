import filecmp
import glob
import json
import math
import os
import random
import re
import shlex
import shutil
import socket
import subprocess
import sys
import tempfile
import time
import unittest

import requests

from base import PROTEIN_MODE, TestingBase, METAGRAPH, TEST_DATA_DIR

"""
End-to-end tests of POST /pattern and `metagraph pattern` (docs/DESIGN-pattern-search.md,
increments 0-3): the counts of the graph contexts of exact DNA and IUPAC patterns, their
label-free extraction (output.labels "none": k-mer, instance, offset, strand, node and row ids,
no annotation read), and their labels (output.labels "all", increment 3: the columns carrying
each context and, on the mini index with its .seqs, every placed occurrence: record, seq_id,
1-based position, strand), checked against PYTHON ORACLES over the very records each index was
built from — never against the engine itself.

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
    annotation and its .seqs. It was built without the dummy-edge mask the route requires (§4),
    so a COPY of its graph in a temporary directory is given the mask by `metagraph transform
    --mask-dummy` (the one-time step for staging) and served with the original annotation, which
    the mask leaves valid; build/mini_refseq itself is never written to. The masks made the other
    ways are checked against it: a rebuild with --mask-dummy (the same .edgemask, the same
    answers) and a server started with --pattern-build-mask on the unmasked graph (the same
    answers, mask: built_at_load). Skipped when the mini index is not built.
  - TestPatternSynthetic: random records, BASIC, CANONICAL and PRIMARY graphs at k = 15 (and a
    multi-graph server), always run.
  - TestPatternFixtureBodies: the frozen fixture bodies (api/python/tests/data/traverse/pattern)
    against the SPEC (api/python/tests/test_pattern_fixtures.py), and the blanking of
    pattern_fixtures.py --check; no server, no index: always run.
  - TestPatternFixtures: pattern_fixtures.py --check, the fixture bodies against what this
    binary answers on the mini index; skipped without the mini index.
  - TestPatternRegression: /search and /align answer byte for byte as the base binary's
    ($METAGRAPH_BASE_BINARY, e.g. the build of 804731aa) on test_api.py's fixture and requests;
    skipped without it.

CI has neither the mini index nor a base binary: there TestPatternMini, TestPatternFixtures and
TestPatternRegression are skipped, which unittest counts as passed (review of 2026-10-07,
T2-01, X-TESTS-05). With $METAGRAPH_REQUIRE_GUARDS=1 a missing input fails those classes instead
of skipping them, so that a run meant to check the guarantees cannot pass without them.
"""

MINI_DIR = os.environ.get('METAGRAPH_MINI_REFSEQ', os.path.join(os.getcwd(), 'mini_refseq'))
# $METAGRAPH_REQUIRE_GUARDS=1: an input a guard class needs (the mini index, the base binary)
# missing fails the class instead of skipping it
REQUIRE_GUARDS = os.environ.get('METAGRAPH_REQUIRE_GUARDS', '') == '1'
FIXTURE_VALIDATOR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'api',
                                 'python', 'tests', 'test_pattern_fixtures.py')
MINI_K = 31
MINI_GRAPH = 'graph_k31.dbg'
MINI_ANNO = 'annotation.relaxed.relabeled.row_diff_brwt_coord.annodbg'
MINI_SEQS = 'annotation.relaxed.relabeled.seqs'
BASE_BINARY = os.environ.get('METAGRAPH_BASE_BINARY', '')
# the generator and checker of the frozen fixture bodies (SPEC-pattern-search.md §11)
FIXTURES_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts',
                               'traversal', 'pattern_fixtures.py')

# the server's defaults (--pattern-* flags; §5.3, the owner's budgets of 2026-10-07): each cap
# is its field's maximum and default, but time_budget_ms defaults to 60 s under a 600 s cap
DEFAULT_CAPS = {'max_contexts': 10000, 'max_anchors': 1000, 'max_steps': 100000000,
                'time_budget_ms': 600000, 'min_information_bits': 24, 'max_patterns': 16,
                # output.labels "all" (increment 3)
                'max_labels_per_anchor': 64, 'max_annotation_work': 100000000,
                'max_memory_mb': 256, 'max_labels': 1000, 'max_occurrences_per_label': 16}
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


# the refusal of a graph without its dummy-edge mask (§4): both remedies named
# a 10-mer whose suffix-scope count on the mini needs a mask scan (it starts a record's first
# k-mer): discovery takes 20 range steps, the scan one more, so max_steps 20 stops in the scan
# with bounds (review of 2026-10-07, X-TESTS-03)
SCAN_10 = 'ATGCCGGTGA'
SCAN_STEPS = 20
# a 31-mer of 12 bases around a run of 19 Ns (24 bits): about 2.5e7 range steps on the mini, a
# discovery a work deadline of a second stops
HEAVY31 = 'GG' + 'N' * 19 + 'TTGGCGATCT'

MASK_REQUIRED_MESSAGE = (
    'pattern: the graph was loaded without its dummy-edge mask (.edgemask): without it every '
    'dummy edge would count as a k-mer and no count would be right; give the graph its mask '
    'once with `metagraph transform --mask-dummy <graph>.dbg` (writes the .edgemask beside the '
    'graph; node ids and annotation unchanged), or pass --pattern-build-mask to server_query '
    'or pattern (builds it in memory at load); the mask is read when the graph is loaded: '
    'restart the server once the .edgemask exists')


def guard(condition, reason):
    """unittest.skipUnless(condition, reason); with $METAGRAPH_REQUIRE_GUARDS=1 a class whose
    input is missing fails instead (review of 2026-10-07, T2-01)."""
    if condition or not REQUIRE_GUARDS:
        return unittest.skipUnless(condition, reason)

    def decorate(cls):
        def fail(klass):
            raise AssertionError(f'$METAGRAPH_REQUIRE_GUARDS=1: {cls.__name__} cannot run: '
                                 f'{reason}')
        cls.setUpClass = classmethod(fail)
        return cls
    return decorate


def untimed(value):
    """An answer without its timings (they differ between any two runs)."""
    if isinstance(value, dict):
        return {k: untimed(v) for k, v in value.items() if k not in ('timing', 'elapsed_ms')}
    if isinstance(value, list):
        return [untimed(v) for v in value]
    return value


def free_port():
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(('127.0.0.1', 0))
        return s.getsockname()[1]


def half_closed_request(port, route, body, timeout=60):
    """The bytes a server sends for a complete POST /route (body: a dict, sent as JSON, or the
    raw text) whose client then half-closed its connection: shutdown(SHUT_WR) right after the
    request, as in the reviewer's repro (review GPT-2 of 2026-10-08, finding 4). The client
    still reads until the server closes. The request is held back until the FIN goes with it
    where the platform can do so (TCP_CORK on Linux, TCP_NOPUSH on macOS), so that the server
    sees the FIN as soon as it has read the request, not after a race."""
    payload = (body if isinstance(body, str) else json.dumps(body)).encode()
    head = (f'POST /{route} HTTP/1.1\r\nHost: localhost\r\nContent-Type: application/json\r\n'
            f'Content-Length: {len(payload)}\r\nConnection: close\r\n\r\n').encode()
    with socket.create_connection(('127.0.0.1', port), timeout=timeout) as s:
        hold = getattr(socket, 'TCP_CORK', None)
        if hold is None and sys.platform == 'darwin':
            hold = 4    # TCP_NOPUSH of <netinet/tcp.h>, which Python does not export
        if hold is not None:
            try:
                s.setsockopt(socket.IPPROTO_TCP, hold, 1)
            except OSError:
                pass
        s.sendall(head + payload)
        s.shutdown(socket.SHUT_WR)
        received = b''
        while True:
            part = s.recv(65536)
            if not part:
                return received
            received += part


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
@guard(os.path.isfile(os.path.join(MINI_DIR, MINI_GRAPH))
       and os.path.isdir(os.path.join(MINI_DIR, 'fasta')),
       "the mini index is not built (scripts/traversal/build_mini_refseq.sh; not in CI)")
class TestPatternMini(PatternChecks, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tempdir = tempfile.TemporaryDirectory()
        d = cls.tempdir.name
        fastas = sorted(glob.glob(os.path.join(MINI_DIR, 'fasta', '*.fa')))
        cls.records = Records.from_fasta(fastas)
        cls.k = MINI_K
        # the annotation's columns: one per FASTA file (annotate --anno-filename, renamed to
        # the file's stem, the taxid), its records in file order — their seq_id in the .seqs
        cls.columns = {os.path.basename(f)[:-len('.fa')]: list(read_fasta(f)) for f in fastas}
        cls.column_records = {c: Records([(c, name, seq) for name, seq in recs])
                              for c, recs in cls.columns.items()}

        cls.fastas = fastas
        # whatever build/mini_refseq holds, it must hold it still when the tests are done
        cls.mini_files = sorted(os.listdir(MINI_DIR))

        # the mini index carries no .edgemask: a copy of its graph gets it from `transform
        # --mask-dummy` (DESIGN §4), which writes only <graph>.edgemask beside the copy; the
        # annotation and its record mapping are the index's own (the mask changes no node id)
        cls.graph = cls._mini_copy('masked', copy_graph=True)
        cls.anno = os.path.join(os.path.dirname(cls.graph), MINI_ANNO)
        cls.label_of_kmer = cls.records.names_with_kmer
        res = TestingBase._run_command(f'{METAGRAPH} transform --mask-dummy -p 2 {cls.graph}',
                                       'Mask the copy of the mini graph')
        cls.transform_log = (res.stdout + res.stderr).decode()
        assert os.path.isfile(cls.graph[:-len('.dbg')] + '.edgemask')
        assert filecmp.cmp(cls.graph, os.path.join(MINI_DIR, MINI_GRAPH), shallow=False)
        # the mini graph as built, without a mask whatever build/mini_refseq holds: a symlink in
        # a directory of its own (the loader reads the mask beside the path it is given)
        cls.unmasked_graph = cls._mini_copy('unmasked', copy_graph=False)
        cls.unmasked_anno = os.path.join(os.path.dirname(cls.unmasked_graph), MINI_ANNO)

        cls.server = Server(METAGRAPH, ['-i', cls.graph, '-a', cls.anno],
                            os.path.join(d, 'server.log'))
        cls.stats = cls.server.get('stats').json()
        cls._choose_patterns()

    @classmethod
    def tearDownClass(cls):
        cls.server.stop()
        cls.tempdir.cleanup()

    @classmethod
    def _mini_copy(cls, name, copy_graph):
        """A directory holding the mini graph (a copy, or a symlink to it) beside symlinks to
        its annotation and record mapping; returns the graph's path there."""
        d = os.path.join(cls.tempdir.name, name)
        os.makedirs(d)
        graph = os.path.join(d, MINI_GRAPH)
        if copy_graph:
            shutil.copyfile(os.path.join(MINI_DIR, MINI_GRAPH), graph)
        else:
            os.symlink(os.path.join(MINI_DIR, MINI_GRAPH), graph)
        for f in (MINI_ANNO, MINI_SEQS):
            os.symlink(os.path.join(MINI_DIR, f), os.path.join(d, f))
        return graph

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
            ({'patterns': p, 'output': {'labels': 'predicate_only'}}, 'later_increment',
             'labels'),
            ({'patterns': p, 'mode': 'count', 'output': {'labels': 'predicate_only'}},
             'later_increment', 'labels'),
            # occurrences are placed per label (increment 3): they need labels "all"
            ({'patterns': p, 'output': {'occurrences': True}}, 'invalid_request', 'occurrences'),
            ({'patterns': p, 'output': {'paths': True}}, 'later_increment', 'paths'),
            ({'patterns': p, 'predicate': {'any': ['562']}}, 'later_increment', 'predicate'),
            ({'patterns': p, 'max_paths': None}, 'later_increment', 'max_paths'),
            ({'patterns': p, 'max_labels': None}, 'invalid_request', 'max_labels'),
            ({'patterns': p, 'max_labels_per_anchor': 0}, 'invalid_request',
             'max_labels_per_anchor'),
            ({'patterns': p, 'allow_unbudgeted_annotation': 'yes'}, 'invalid_request',
             'allow_unbudgeted_annotation'),
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
        # one RFC 8259 JSON text with unique member names (review of 2026-10-07, R1-02: each
        # was answered 200), nested at most 1,000 deep (R1-03, R2-01: a 400 without a code)
        one = '{"patterns": [{"dna": "%s"}]' % self.p16
        for raw in (one + '} GARBAGE', one + ',}', one + '} ' + one + '}',
                    one + ', /* x */ "mode": "count"}', one + ', "mode": "count", "mode": "partial"}',
                    one + ', "max_steps": 1, "max_steps": 100000}',
                    '{"patterns": [{"dna": "%s", "dna": "%s"}]}' % (self.p16, self.p14),
                    '[' * 1200 + ']' * 1200,
                    one + ', "x": ' + '[' * 1500 + ']' * 1500 + '}'):
            ret = self.server.post('pattern', raw, raw=True)
            self.assertEqual(400, ret.status_code, raw[:100])
            self.assertEqual({'error', 'code'}, set(ret.json()), raw[:100])
            self.assertEqual('invalid_request', ret.json()['code'], raw[:100])
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
            'projections': ['none', 'all'], 'default_projection': 'none',
            'projections_later_increment': ['predicate_only'], 'default_occurrences': True,
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

    def test_bounds_in_a_mask_scan(self):
        """max_steps reached in a mask scan (review of 2026-10-07, X-TESTS-03: no test answer
        had relation bounds or phase mask_scan): the count is bounds, value = lower <= the
        oracle's count <= upper, the parts summed as SPEC §7.4 says, and the stop names the
        scan."""
        truth = self.contexts(SCAN_10, scope='suffix')
        full = self.pattern(self.server, {'patterns': [{'dna': SCAN_10}], 'mode': 'count',
                                          'scope': 'suffix'})['patterns'][0]
        self.assertContextCounts(full, truth, self.k, 10, 'suffix')
        self.assertGreaterEqual(full['work']['mask_scans'], 1)
        self.assertGreater(full['work']['steps'], SCAN_STEPS)
        for mode in ('count', 'all_or_count'):
            out = self.pattern(self.server, {'patterns': [{'dna': SCAN_10}], 'mode': mode,
                                             'scope': 'suffix', 'max_steps': SCAN_STEPS})
            e = out['patterns'][0]
            self.assertEqual({'phase': 'mask_scan', 'reason': 'max_steps'}, e['stop'])
            c = e['counts']['contexts']
            self.assertEqual('bounds', c['relation'])
            self.assertEqual(c['lower'], c['value'])
            self.assertLessEqual(c['lower'], len(truth))
            self.assertLessEqual(len(truth), c['upper'])
            parts = list(c['by_strand'].values())
            self.assertEqual(c['lower'], sum(p.get('lower', p['value']) for p in parts))
            self.assertEqual(c['upper'], sum(p.get('upper', p['value']) for p in parts))
            for strand, p in c['by_strand'].items():
                n = sum(1 for t, _, _ in truth if t == strand)
                if p['relation'] == 'exact':
                    self.assertEqual(n, p['value'], strand)
                else:
                    self.assertEqual('bounds', p['relation'], strand)
                    self.assertLessEqual(p['lower'], n)
                    self.assertLessEqual(n, p['upper'])
            self.assertEqual('full', e['determinism'])
            if mode == 'all_or_count':
                self.assertEqual({'reason': 'discovery_budget'}, e['withheld'])
                self.assertEqual([], e['results'])

    def test_a_stopped_request_answers_what_it_buffered(self):
        """Review of 2026-10-07, X-EFFICIENCY-04: the work stops earlier by the estimated time
        to write what the answer holds, so that a request whose earlier patterns buffered
        many results and whose last one runs to the work deadline still answers within its
        budget, with its counts and a time stop, instead of 503 with all of its work lost.
        Fifteen patterns of 19,283 results each (a server whose max_contexts admits them all)
        and a last one whose discovery outlasts the budget: with the 250 ms reserve alone the
        work ran to 1,750 ms and the ~33 MB of results could not be written by 2,000 ms (503 in
        every run, measured on the binary before the fix)."""
        server = Server(METAGRAPH, ['-i', self.graph, '-a', self.anno,
                                    '--pattern-max-contexts', '20000'],
                        os.path.join(self.tempdir.name, 'server_volume.log'))
        try:
            request = {'patterns': [{'dna': 'ACGT'}] * 15 + [{'iupac': HEAVY31}],
                       'mode': 'partial', 'scope': 'suffix', 'strands': 'forward',
                       'time_budget_ms': 2000, 'output': {'labels': 'none'}}
            ret = requests.post(server.url('pattern'), data=json.dumps(request),
                                headers={'Accept-Encoding': 'gzip'}, timeout=120)
        finally:
            server.stop()
        self.assertEqual(200, ret.status_code, ret.text[:500])
        self.assertEqual('gzip', ret.headers.get('Content-Encoding'))
        out = ret.json()
        entries = out['patterns']
        stopped = [i for i, e in enumerate(entries)
                   if e['stop'] and e['stop']['reason'] == 'time']
        self.assertTrue(stopped, [e['stop'] for e in entries])
        first = stopped[0]
        # the patterns before the stop: complete
        expected = len(self.contexts('ACGT', scope='suffix', strands='forward'))
        for e in entries[:first]:
            self.assertCount({x: e['counts']['contexts'][x]
                              for x in ('value', 'relation', 'unit')}, expected)
            self.assertCompleteRetrieval(e)
            self.assertEqual(expected, e['returned'])
        # the stopped one, and every one after it as after any time stop, its counts kept
        self.assertEqual('time_limited', entries[first]['determinism'])
        self.assertEqual({'reason': 'time'}, entries[first]['cut'])
        self.assertLessEqual(entries[first]['returned'], expected)
        for e in entries[first + 1:]:
            self.assertEqual({'phase': 'discovery', 'reason': 'time'}, e['stop'])
            self.assertEqual('unknown', e['counts']['contexts']['relation'])
            self.assertEqual(0, e['returned'])
        self.assertLessEqual(out['timing']['elapsed_ms'], 2000)

    def test_resolve_deadline_is_503(self):
        """/resolve's 503 end to end (review of 2026-10-07, V1-04: only the exception was
        tested): a query whose k-mer mapping, which the deadline cannot interrupt, outlasts
        the whole budget is answered 503 {error, code: deadline}, uncompressed, without
        Retry-After (that is the loading 503's), and the CLI writes the same body and exits 1."""
        # the whole 7 Mbp record: its mapping takes well over the 250 ms reserve on any host
        # (2 Mbp mapped within it on a quiet machine)
        seq = ''.join(s for _, s in read_fasta(os.path.join(MINI_DIR, 'fasta', '287.fa')))
        body = {'sequence': seq, 'discover': {'max_labels': 10},
                'bounds': {'time_budget_ms': 250.001, 'max_query_bp': len(seq)}}
        expected = {'error': 'resolve: the answer could not be built and written within '
                             'bounds.time_budget_ms (250.001 ms, the finalisation reserve of '
                             '250 ms included): nothing partial is sent',
                    'code': 'deadline'}
        for headers in ({}, {'Accept-Encoding': 'gzip'}):
            ret = requests.post(self.server.url('resolve'), data=json.dumps(body),
                                headers=headers, timeout=300)
            self.assertEqual(503, ret.status_code, ret.text[:500])
            self.assertNotIn('Retry-After', ret.headers)
            self.assertNotIn('Content-Encoding', ret.headers)
            self.assertEqual(expected, ret.json())
        path = os.path.join(self.tempdir.name, 'resolve_deadline.json')
        with open(path, 'w') as f:
            json.dump(body, f)
        res = subprocess.run(shlex.split(METAGRAPH) + ['traverse', '--resolve', '--json',
                                                       '-i', self.graph, '-a', self.anno, path],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(1, res.returncode, res.stderr.decode()[-2000:])
        self.assertEqual(expected, json.loads(res.stdout))

    def test_long_pattern(self):
        out = self.pattern(self.server, {'patterns': [{'dna': self.p40},
                                                      {'dna': self.absent40}]})
        entry, absent = out['patterns']
        self.assertEqual('long', entry['scope'])
        self.assertEqual('long', entry['absence_scope'])
        self.assertEqual(80.0, entry['information_bits'])
        self.assertEqual(62.0, entry['anchor_information_bits'])
        # both windows of this 40-mer carry 62 bits (an addition of the review of 2026-10-07,
        # X-GUARANTEES-01: the least searched window, the floor's operand)
        self.assertEqual(62.0, entry['min_anchor_information_bits'])
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

    # ------------------------------------------------------------ labels (increment 3)

    NDM_F = 'GGTTTGGCGATCTGGTTTTC'   # blaNDM-1 forward primer: its k-mers carry 7 to 9 taxa

    def placed_oracle(self, pattern):
        """{(column, seq_id, 1-based start, strand)}: every occurrence of the oriented pattern
        in every record of the FASTA files (a regex scan; the mini's records hold no N)."""
        out = set()
        for column, recs in self.columns.items():
            for seq_id, (_, seq) in enumerate(recs):
                for strand, q in Records.oriented(pattern, 'both', True):
                    for m in pattern_regex(q).finditer(seq.upper()):
                        out.add((column, seq_id, m.start() + 1, strand))
        return out

    def contexts_by_column(self, pattern):
        out = {c: r.contexts(pattern, self.k) for c, r in self.column_records.items()}
        return {c: x for c, x in out.items() if x}

    def assertLabelled(self, entry, pattern):
        """output.labels "all", complete, record placement: every context with exactly the
        columns whose records hold its k-mer; every placed occurrence the FASTA scan's, its
        record's own k-mer at the occurrence being the context's (the record first, the offset
        after); by_label and the counts exact."""
        L, k = len(pattern), self.k
        by_column = self.contexts_by_column(pattern)
        expected = set().union(*by_column.values()) if by_column else set()
        self.assertCompleteRetrieval(entry)
        self.assertEqual(('record', 'budgeted'), (entry['placement'], entry['annotation']))
        self.assertEqual(([], []), (entry['rows_refused'], entry['anchors_truncated']))
        self.assertEqual(len(expected), len(entry['results']))
        rank = {b['column']: i for i, b in enumerate(entry['by_label'])}
        unions = {}
        for r in entry['results']:
            context = (r['strand'], r['kmer'], r['offset'])
            self.assertIn(context, expected)
            columns = {c for c, x in by_column.items() if context in x}
            self.assertEqual(('kmer', 'complete', len(columns)),
                             (r['support'], r['labels_status'], r['labels_total']))
            listed = [label['column'] for label in r['labels']]
            self.assertEqual(columns, set(listed), r['kmer'])
            # each result's labels in the label order of by_label
            self.assertEqual(sorted(listed, key=rank.get), listed)
            for label in r['labels']:
                self.assertEqual('kmer', label['support'])
                occ = label['occurrence_list']
                self.assertCount(label['occurrences'], len(occ), unit='placed_occurrences')
                self.assertTrue(occ)
                for o in occ:
                    name, seq = self.columns[label['column']][o['seq_id']]
                    self.assertEqual((name, len(seq), r['strand']),
                                     (o['record'], o['nt_length'], o['strand']))
                    a, b = (int(x) for x in o['nt_coords'].split('-'))
                    self.assertEqual(a + L - 1, b)
                    start = a - 1 - r['offset']
                    self.assertEqual(r['kmer'], seq[start:start + k].upper(), o)
                    unions.setdefault(label['column'], set()).add((o['seq_id'], a, o['strand']))
        placed = {(c, i, a, st) for c, x in unions.items() for i, a, st in x}
        self.assertEqual(self.placed_oracle(pattern), placed)
        # by_label: every column once, (contexts desc, column asc), exact
        self.assertEqual(set(by_column), set(rank))
        keys = [(-b['contexts']['value'], b['column']) for b in entry['by_label']]
        self.assertEqual(sorted(keys), keys)
        for b in entry['by_label']:
            c = b['column']
            self.assertIsNone(b['graph'])
            self.assertCount(b['contexts'], len(by_column[c]))
            self.assertCount(b['contexts_suffix'], sum(1 for x in by_column[c] if x[2] == k - L))
            self.assertCount(b['occurrences'], len(unions[c]), unit='placed_occurrences')
        self.assertCount(entry['counts']['labels'], len(by_column), unit='labels')
        self.assertCount(entry['counts']['occurrences'], len(placed), unit='placed_occurrences')

    def test_labels_all_against_the_fasta(self):
        patterns = [self.NDM_F, self.p16, self.p14, self.iupac16, self.pal12, self.absent16]
        request = {'patterns': [{'iupac' if set(p) - set('ACGT') else 'dna': p}
                                for p in patterns], 'output': {'labels': 'all'}}
        out = self.pattern(self.server, request)
        self.assertEqual({'labels': 'all', 'occurrences': True}, out['output'])
        for x in ('max_labels_per_anchor', 'max_annotation_work', 'max_memory_mb',
                  'max_labels', 'max_occurrences_per_label'):
            self.assertEqual(DEFAULT_CAPS[x], out['limits'][x], x)
        self.assertFalse(out['limits']['allow_unbudgeted_annotation'])
        for entry, p in zip(out['patterns'], patterns):
            with self.subTest(pattern=p):
                c = entry['counts']['contexts']
                self.assertEqual(('exact', len(self.contexts(p))), (c['relation'], c['value']))
                self.assertLabelled(entry, p)
                self.assertGreater(entry['work']['memory_bytes'], 0)
        ndm = out['patterns'][0]
        self.assertGreater(len(ndm['by_label']), 1)
        # the absent pattern: complete, exact zeros, nothing read
        absent = out['patterns'][-1]
        self.assertEqual(([], 0), (absent['by_label'], absent['work']['annotation_rows']))
        # deterministic: the same request, the same labels and occurrences
        again = self.pattern(self.server, request)
        self.assertEqual(untimed(out), untimed(again))

    def test_labels_all_withheld_truncated_and_partial(self):
        p = self.NDM_F
        n = len(self.contexts(p))
        # above the threshold: not admitted, nothing read
        entry = self.pattern(self.server, {'patterns': [{'dna': p}], 'max_contexts': n - 1,
                                           'output': {'labels': 'all'}})['patterns'][0]
        self.assertEqual({'reason': 'count_above_threshold'}, entry['withheld'])
        self.assertEqual(0, entry['work']['annotation_rows'])
        self.assertIsNone(entry['by_label'])
        self.assertUnknownLabels(entry)
        # a per-anchor cap below the rows' labels: withheld, every cut anchor stated
        by_column = self.contexts_by_column(p)
        entry = self.pattern(self.server, {'patterns': [{'dna': p}], 'max_labels_per_anchor': 2,
                                           'output': {'labels': 'all'}})['patterns'][0]
        self.assertEqual({'reason': 'anchor_labels_truncated'}, entry['withheld'])
        self.assertEqual([], entry['results'])
        totals = {}
        for c in set().union(*by_column.values()):
            totals[c[1]] = sum(1 for x in by_column.values() if c in x)
        got = {t['kmer']: (t['cap'], t['total']) for t in entry['anchors_truncated']}
        self.assertEqual({kmer: (2, t) for kmer, t in totals.items() if t > 2}, got)
        # partial: the first labels of each cut row, the counts at_least
        entry = self.pattern(self.server, {'patterns': [{'dna': p}], 'mode': 'partial',
                                           'max_labels_per_anchor': 2,
                                           'output': {'labels': 'all'}})['patterns'][0]
        self.assertIsNone(entry['withheld'])
        self.assertFalse(entry['retrieval_complete'])
        for r in entry['results']:
            self.assertEqual(('truncated', totals[r['kmer']], 2),
                             (r['labels_status'], r['labels_total'], len(r['labels'])))
            # the first two in ascending column order... of the annotation's column ids,
            # which are not the names' order: only membership is checked
            self.assertLessEqual({x['column'] for x in r['labels']},
                                 {c for c, x in by_column.items()
                                  if (r['strand'], r['kmer'], r['offset']) in x})
        self.assertEqual('at_least', entry['counts']['labels']['relation'])
        # partial's lists: labels and occurrences cut, the cuts stated, the counts whole
        full = self.pattern(self.server, {'patterns': [{'dna': p}],
                                          'output': {'labels': 'all'}})['patterns'][0]
        entry = self.pattern(self.server, {'patterns': [{'dna': p}], 'mode': 'partial',
                                           'max_labels': 3, 'max_occurrences_per_label': 2,
                                           'output': {'labels': 'all'}})['patterns'][0]
        self.assertEqual({'reason': 'max_labels', 'returned': 3}, entry['labels_cut'])
        self.assertEqual('max_occurrences_per_label', entry['occurrences_cut']['reason'])
        self.assertEqual(full['by_label'][:3], entry['by_label'])
        self.assertEqual(full['counts'], entry['counts'])
        self.assertFalse(entry['retrieval_complete'])
        listed = {}
        for r in entry['results']:
            for label in r['labels']:
                for o in label['occurrence_list']:
                    listed.setdefault(label['column'], set()).add(
                        (o['seq_id'], o['nt_coords'], o['strand']))
        self.assertEqual({b['column'] for b in full['by_label'][:3]}, set(listed))
        for column, x in listed.items():
            self.assertLessEqual(len(x), 2, column)

    def test_labels_all_count_mode_and_projection_none(self):
        """Mode count reads no annotation whatever the projection, and says so; an annotation
        field with output.labels "none" likewise; both answer the counts as without them."""
        plain = self.pattern(self.server, {'patterns': [{'dna': self.p16}], 'mode': 'count'})
        out = self.pattern(self.server, {'patterns': [{'dna': self.p16}], 'mode': 'count',
                                         'output': {'labels': 'all'}})
        self.assertIsNone(out['output'])
        self.assertEqual(['annotation_not_read'], out['patterns'][0]['notes'])
        out['patterns'][0]['notes'] = []
        self.assertEqual(untimed(plain), untimed(out))
        none = self.pattern(self.server, {'patterns': [{'dna': self.p16}]})
        out = self.pattern(self.server, {'patterns': [{'dna': self.p16}], 'max_memory_mb': 1})
        self.assertEqual(['annotation_not_read'], out['patterns'][0]['notes'])
        out['patterns'][0]['notes'] = []
        self.assertEqual(untimed(none), untimed(out))

    def test_labels_all_cli_answers_as_the_server(self):
        request = {'patterns': [{'dna': self.NDM_F}, {'iupac': self.iupac16}, {'dna': self.p40}],
                   'mode': 'partial', 'max_contexts': 9, 'max_labels': 4,
                   'output': {'labels': 'all'}}
        server_out = self.pattern(self.server, request)
        path = os.path.join(self.tempdir.name, 'request_labels.json')
        with open(path, 'w') as f:
            json.dump(request, f)
        res = subprocess.run(shlex.split(METAGRAPH) + ['pattern', '--json', '-i', self.graph,
                                                       '-a', self.anno, path],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(0, res.returncode, res.stderr.decode())
        self.assertEqual(untimed(server_out), untimed(json.loads(res.stdout)))

    def test_labels_all_without_record_mapping(self):
        """--no-coord-mapping: coordinates without the .seqs. Each label lists the column
        coordinate of the context's k-mer (kmer_coord) and the offset, placed nowhere; no
        occurrence count is claimed (note record_bounds_unknown)."""
        server = Server(METAGRAPH, ['-i', self.graph, '-a', self.anno, '--no-coord-mapping'],
                        os.path.join(self.tempdir.name, 'server_no_map.log'))
        try:
            caps = server.get('capabilities').json()['pattern']
            self.assertEqual('global', caps['placement'])
            entry = self.pattern(server, {'patterns': [{'dna': self.NDM_F}],
                                          'output': {'labels': 'all'}})['patterns'][0]
        finally:
            server.stop()
        self.assertEqual('global', entry['placement'])
        self.assertEqual(['record_bounds_unknown'], entry['notes'])
        self.assertTrue(entry['retrieval_complete'])
        self.assertCount(entry['counts']['occurrences'], None, 'unknown', 'placed_occurrences')
        by_column = self.contexts_by_column(self.NDM_F)
        self.assertCount(entry['counts']['labels'], len(by_column), unit='labels')
        # the k-mers' column coordinates: each record's k-mers numbered from its column's
        # running count, in file order
        starts = {}
        for column, recs in self.columns.items():
            s, at = [], 0
            for _, seq in recs:
                s.append(at)
                at += len(seq) - self.k + 1
            starts[column] = s
        seen = 0
        for r in entry['results']:
            for label in r['labels']:
                self.assertNotIn('occurrences', label)
                recs = self.columns[label['column']]
                for o in label['occurrence_list']:
                    seen += 1
                    self.assertEqual((r['offset'], r['strand']), (o['offset'], o['strand']))
                    c = o['kmer_coord']
                    j = max(i for i, s in enumerate(starts[label['column']]) if s <= c)
                    local = c - starts[label['column']][j]
                    self.assertEqual(r['kmer'], recs[j][1][local:local + self.k].upper())
        self.assertGreater(seen, 0)

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

        # every request file is answered, the later ones after a refused one (review of
        # 2026-10-07, R1-03, R2-01: a body nested too deep aborted the run, exit 134, and
        # left the later files unanswered)
        deep = os.path.join(self.tempdir.name, 'deep.json')
        with open(deep, 'w') as f:
            f.write('[' * 1200 + ']' * 1200)
        dup = os.path.join(self.tempdir.name, 'dup.json')
        with open(dup, 'w') as f:
            f.write('{"patterns": [{"dna": "%s"}], "mode": "count", "mode": "partial"}' % self.p16)
        ok = os.path.join(self.tempdir.name, 'ok.json')
        with open(ok, 'w') as f:
            json.dump({'patterns': [{'dna': self.p16}], 'mode': 'count'}, f)
        res = subprocess.run(shlex.split(METAGRAPH) + ['pattern', '--json', '-i', self.graph,
                                                       '-a', self.anno, deep, dup, ok],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(1, res.returncode, res.stderr.decode()[-2000:])
        lines = res.stdout.decode().strip().split('\n')
        self.assertEqual(3, len(lines), lines)
        self.assertEqual('invalid_request', json.loads(lines[0])['code'])
        self.assertEqual('invalid_request', json.loads(lines[1])['code'])
        self.assertContextCounts(json.loads(lines[2])['patterns'][0], self.contexts(self.p16),
                                 self.k, 16, 'any_offset')

    def test_cli_mask_required(self):
        """The mini index as built (no .edgemask): refused, since every dummy edge would count;
        with --pattern-build-mask the CLI answers as the server on the transformed copy."""
        request = {'patterns': [{'dna': self.p16}, {'iupac': self.iupac16}], 'mode': 'count'}
        path = os.path.join(self.tempdir.name, 'request_mask.json')
        with open(path, 'w') as f:
            json.dump(request, f)
        cli = shlex.split(METAGRAPH) + ['pattern']
        args = ['--json', '-i', self.unmasked_graph, '-a', self.unmasked_anno, path]
        res = subprocess.run(cli + args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(1, res.returncode, res.stderr.decode())
        self.assertEqual({'error': MASK_REQUIRED_MESSAGE, 'code': 'mask_required'},
                         json.loads(res.stdout))
        # the start-up note names the remedies before any request is refused
        self.assertIn('the pattern search answers mask_required. Remedies: `metagraph transform '
                      f'--mask-dummy {self.unmasked_graph}` once', res.stderr.decode())

        res = subprocess.run(cli + ['--pattern-build-mask'] + args,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(0, res.returncode, res.stderr.decode())
        self.assertEqual(untimed(self.pattern(self.server, request)),
                         untimed(json.loads(res.stdout)))

    def test_cli_unreadable_mask(self):
        """A .edgemask that is there but cannot be opened (permissions): the graph loads without
        it, as it always did, and the messages say so rather than "no mask", whose remedy
        `transform --mask-dummy` would refuse (review of 2026-10-07, M1-03)."""
        d = os.path.join(self.tempdir.name, 'unreadable')
        os.makedirs(d)
        graph = os.path.join(d, MINI_GRAPH)
        os.symlink(self.graph, graph)
        for f in (MINI_ANNO, MINI_SEQS):
            os.symlink(os.path.join(MINI_DIR, f), os.path.join(d, f))
        mask = graph[:-len('.dbg')] + '.edgemask'
        shutil.copyfile(self.graph[:-len('.dbg')] + '.edgemask', mask)
        os.chmod(mask, 0)
        try:
            if os.access(mask, os.R_OK):
                self.skipTest('the mask stays readable (running as root?)')
            request = {'patterns': [{'dna': self.p16}], 'mode': 'count'}
            path = os.path.join(d, 'request.json')
            with open(path, 'w') as f:
                json.dump(request, f)
            res = subprocess.run(shlex.split(METAGRAPH) + ['pattern', '--json', '-i', graph,
                                                           '-a', os.path.join(d, MINI_ANNO),
                                                           path],
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            err = res.stderr.decode()
            self.assertEqual(1, res.returncode, err)
            self.assertEqual('mask_required', json.loads(res.stdout)['code'])
            # the loader's warning and the start-up note both name the file and the cause
            self.assertIn(f'The dummy-edge mask {mask} exists but could not be opened', err)
            self.assertIn('--mask-dummy --force', err)
            self.assertNotIn('The graph has no dummy-edge mask', err)
        finally:
            os.chmod(mask, 0o644)

    # ------------------------------------------------------------ the mask, made three ways

    def mask_panel(self):
        """/pattern requests whose answers depend on every masked edge: counts in both scopes,
        both retrievals, IUPAC, a palindrome, a pattern longer than k, absent ones, a stop in
        discovery, and a stop inside a mask scan (its bounds count the scanned edges; review of
        2026-10-07, X-TESTS-03: the panel's only stop was in discovery)."""
        patterns = [{'dna': self.p16}, {'dna': self.p14}, {'iupac': self.iupac16},
                    {'dna': self.pal12}, {'dna': self.p40}, {'dna': self.absent16},
                    {'dna': self.absent40}]
        panel = []
        for scope in ('any_offset', 'suffix'):
            for mode in ('count', 'all_or_count', 'partial'):
                panel.append({'patterns': patterns, 'mode': mode, 'scope': scope,
                              'max_contexts': 25})
        panel.append({'patterns': [{'dna': self.p16}, {'iupac': self.iupac16}], 'mode': 'count',
                      'max_steps': 40})
        panel.append({'patterns': [{'dna': SCAN_10}], 'mode': 'count', 'scope': 'suffix',
                      'max_steps': SCAN_STEPS})
        return panel

    def assertSameAnswers(self, expected_server, server):
        for request in self.mask_panel():
            a, b = (expected_server.post('pattern', request), server.post('pattern', request))
            self.assertEqual(200, a.status_code, a.text)
            self.assertEqual((a.status_code, untimed(a.json())), (b.status_code, untimed(b.json())),
                             request)

    def test_transform_left_the_index_alone(self):
        """transform wrote <copy>.edgemask and nothing else; build/mini_refseq is unchanged."""
        self.assertEqual(self.mini_files, sorted(os.listdir(MINI_DIR)))
        self.assertEqual(sorted([MINI_GRAPH, MINI_GRAPH[:-len('.dbg')] + '.edgemask',
                                 MINI_ANNO, MINI_SEQS]),
                         sorted(os.listdir(os.path.dirname(self.graph))))
        self.assertTrue(filecmp.cmp(self.graph, os.path.join(MINI_DIR, MINI_GRAPH),
                                    shallow=False))
        # the counts it logged: every edge is a k-mer, a source or a sink dummy, and the
        # masked graph states its k-mers as its nodes
        m = re.search(r'(\d+) edges, (\d+) source dummies \(the main dummy edge included\), '
                      r'(\d+) sink dummies, (\d+) k-mers', self.transform_log)
        self.assertIsNotNone(m, self.transform_log)
        edges, source, sink, kmers = (int(x) for x in m.groups())
        self.assertEqual(edges, source + sink + kmers)
        self.assertEqual(edges, self.stats['annotation']['objects'])
        self.assertEqual(kmers, self.stats['graph']['nodes'])

    def test_transform_mask_equals_build_mask(self):
        """The mask transform gave the copy is the one `build --mask-dummy` writes, and the
        answers on a masked rebuild are the answers on the transformed copy."""
        d = os.path.join(self.tempdir.name, 'rebuilt')
        os.makedirs(d)
        graph = os.path.join(d, MINI_GRAPH)
        # the same records and flags as the mini index (commands.log), single-threaded so that
        # nothing depends on the order of threads
        TestingBase._run_command(
            f'{METAGRAPH} build -p 1 --mode basic --graph succinct --state stat -k {MINI_K} '
            f'--index-ranges 12 --mask-dummy --in-ram -o {graph[:-len(".dbg")]} '
            + ' '.join(self.fastas), 'Build the masked mini graph')
        if not filecmp.cmp(graph, os.path.join(MINI_DIR, MINI_GRAPH), shallow=False):
            self.skipTest('the rebuild is not the mini graph (built by another metagraph?): '
                          'its masks cannot be compared')
        self.assertTrue(filecmp.cmp(graph[:-len('.dbg')] + '.edgemask',
                                    self.graph[:-len('.dbg')] + '.edgemask', shallow=False))
        for f in (MINI_ANNO, MINI_SEQS):
            os.symlink(os.path.join(MINI_DIR, f), os.path.join(d, f))
        server = Server(METAGRAPH, ['-i', graph, '-a', os.path.join(d, MINI_ANNO)],
                        os.path.join(d, 'server.log'))
        try:
            self.assertSameAnswers(self.server, server)
            self.assertEqual(self.stats, server.get('stats').json())
        finally:
            server.stop()

    def test_pattern_build_mask(self):
        """A server started with --pattern-build-mask on the unmasked graph builds the same mask
        in memory: the same answers on every route of the panel, mask: built_at_load."""
        d = self.tempdir.name
        log = os.path.join(d, 'server_build_mask.log')
        server = Server(METAGRAPH, ['-i', self.unmasked_graph, '-a', self.unmasked_anno,
                                    '--pattern-build-mask'], log)
        try:
            for route in ('capabilities', 'traverse/capabilities'):
                p = server.get(route).json()['pattern']
                self.assertEqual('built_at_load', p['mask'])
                self.assertTrue(p['available'])
                # the rest of the block as on the server whose mask is a file
                q = self.server.get(route).json()['pattern']
                self.assertEqual(dict(q, mask='built_at_load'), p)
            self.assertSameAnswers(self.server, server)
            # the graph is the masked graph for every route: /stats states its k-mers
            self.assertEqual(self.stats, server.get('stats').json())
        finally:
            server.stop()
        with open(log) as f:
            text = f.read()
        self.assertRegex(text, r'--pattern-build-mask: dummy-edge mask built in [\d.]+ s')
        # built in memory only
        self.assertEqual(sorted([MINI_GRAPH, MINI_ANNO, MINI_SEQS]),
                         sorted(os.listdir(os.path.dirname(self.unmasked_graph))))
        self.assertFalse(os.path.exists(self.unmasked_graph[:-len('.dbg')] + '.edgemask'))

    def test_mask_required(self):
        """Without the mask: absent in the capabilities, 400 mask_required naming both remedies."""
        server = Server(METAGRAPH, ['-i', self.unmasked_graph, '-a', self.unmasked_anno],
                        os.path.join(self.tempdir.name, 'server_unmasked.log'))
        try:
            for route in ('capabilities', 'traverse/capabilities'):
                p = server.get(route).json()['pattern']
                self.assertEqual(('absent', False, 'mask_required'),
                                 (p['mask'], p['available'], p['unavailable_reason']))
            ret = server.post('pattern', {'patterns': [{'dna': self.p16}], 'mode': 'count'})
            self.assertEqual(400, ret.status_code)
            self.assertEqual({'error': MASK_REQUIRED_MESSAGE, 'code': 'mask_required'},
                             ret.json())
        finally:
            server.stop()


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

    def test_labels_all_on_unbudgeted_annotations(self):
        """A column annotation has no budget-aware decode: output.labels "all" is refused
        (annotation_unbudgeted) unless the request allows the unbudgeted reads; then each
        context carries the records (--anno-header columns) holding its k-mer — on CANONICAL
        and PRIMARY graphs its k-mer or the reverse complement — and nothing is placed."""
        for mode, server in self.servers.items():
            stated = mode == 'basic'
            request = {'patterns': [{'dna': 'ACGAC'}, {'dna': 'GAATTC'}],
                       'output': {'labels': 'all'}}
            ret = server.post('pattern', request)
            self.assertEqual(400, ret.status_code, ret.text)
            self.assertEqual('annotation_unbudgeted', ret.json()['code'])
            self.assertIn('allow_unbudgeted_annotation', ret.json()['error'])
            out = self.pattern(server, dict(request, allow_unbudgeted_annotation=True))
            self.assertTrue(out['limits']['allow_unbudgeted_annotation'])
            for entry, p in zip(out['patterns'], ('ACGAC', 'GAATTC')):
                with self.subTest(mode=mode, pattern=p):
                    self.assertCompleteRetrieval(entry)
                    self.assertEqual('unbudgeted', entry['annotation'])
                    self.assertIn('annotation_unbudgeted', entry['notes'])
                    self.assertEqual('none' if stated else 'none_canonical', entry['placement'])
                    self.assertCount(entry['counts']['occurrences'], None, 'unknown',
                                     'placed_occurrences')
                    labels = {}
                    for r in entry['results']:
                        kmer = r['kmer']
                        names = self.records.names_with_kmer(kmer)
                        if not stated:
                            names |= self.records.names_with_kmer(revcomp(kmer))
                        got = {x['column'] for x in r['labels']}
                        self.assertEqual(names, got, kmer)
                        self.assertEqual(len(names), r['labels_total'])
                        for x in r['labels']:
                            self.assertEqual({'column', 'support'}, set(x))
                            labels[x['column']] = labels.get(x['column'], 0) + 1
                    self.assertCount(entry['counts']['labels'], len(labels), unit='labels')
                    self.assertEqual(labels, {b['column']: b['contexts']['value']
                                              for b in entry['by_label']})

    def test_a_client_that_left_is_not_answered(self):
        """SPEC §3: a request whose client has left is not answered, a half-close counted as
        gone -- its errors no more than its answer (review GPT-2 of 2026-10-08, finding 4:
        after a complete request and shutdown(SHUT_WR), a valid count request got 0 bytes, but
        `{` got 194 (400) and {"patterns":[]} 177 (400): the error path wrote without asking).
        The same requests from a client that is still there are answered as before."""
        server = self.servers['basic']
        cases = [
            ('a count', {'patterns': [{'dna': 'ACGAC'}], 'mode': 'count'}, 200, None),
            ('a retrieval', {'patterns': [{'dna': 'ACGAC'}]}, 200, None),
            ('not JSON', '{', 400, 'invalid_request'),
            ('no patterns', {'patterns': []}, 400, 'invalid_request'),
            ('a later increment', {'patterns': [{'dna': 'ACGAC'}],
                                   'output': {'labels': 'predicate_only'}}, 400,
             'later_increment'),
            # refused once the request is parsed (the column annotation has no budgeted reads)
            ('unbudgeted labels', {'patterns': [{'dna': 'ACGAC'}], 'output': {'labels': 'all'}},
             400, 'annotation_unbudgeted'),
        ]
        for name, body, status, code in cases:
            with self.subTest(case=name):
                ret = server.post('pattern', body, raw=isinstance(body, str))
                self.assertEqual(status, ret.status_code, ret.text)
                if code:
                    self.assertEqual(code, ret.json()['code'])
                for _ in range(3):
                    self.assertEqual(b'', half_closed_request(server.port, 'pattern', body))
        # the server serves on
        self.pattern(server, {'patterns': [{'dna': 'ACGAC'}], 'mode': 'count'})

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
            # nor is this refusal written to a client that left (review GPT-2, finding 4)
            self.assertEqual(b'', half_closed_request(
                server.port, 'pattern', {'patterns': [{'dna': 'ACGAC'}], 'mode': 'count'}))
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

    def test_deadline_caps_under_the_content_timeout(self):
        """A route deadline past the server's 900 s content timeout would be accepted and the
        connection closed before any answer or 503 (review of 2026-10-07, R1-05): server_query
        refuses a /pattern cap above 899,000 ms, and /resolve's opt-in deadline is capped
        there when --traverse-max-time-ms is 0 (uncapped) or above it."""
        cmd = shlex.split(METAGRAPH) + ['server_query', '-i', self.graph_basic, '-a',
                                        self.anno_basic, '--pattern-max-time-ms', '900000',
                                        '--port', str(free_port()), '--address', '127.0.0.1']
        try:
            res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=120)
        except subprocess.TimeoutExpired:
            self.fail('server_query started with --pattern-max-time-ms 900000')
        self.assertNotEqual(0, res.returncode)
        self.assertIn('--pattern-max-time-ms must be at most 899000', res.stderr.decode())

        server = Server(METAGRAPH, ['-i', self.graph_basic, '-a', self.anno_basic,
                                    '--pattern-max-time-ms', '899000',
                                    '--traverse-max-time-ms', '0'],
                        os.path.join(self.tempdir.name, 'server_caps.log'))
        try:
            caps = server.get('capabilities').json()
            self.assertEqual(899000, caps['resolve']['time_budget']['max_time_ms'])
            self.assertEqual(899000, caps['pattern']['caps']['time_budget_ms'])
            sequence = next(self.records.sequences(False))[:300]
            ret = server.post('resolve', {'sequence': sequence, 'discover': {'max_labels': 2},
                                          'bounds': {'time_budget_ms': 1000000}})
            self.assertEqual(200, ret.status_code, ret.text)
            self.assertEqual({'time_budget_ms': 899000, 'finalize_reserve_ms': 250,
                              'clamped': [{'field': 'bounds.time_budget_ms',
                                           'requested': 1000000, 'effective': 899000}]},
                             ret.json()['limits'])
        finally:
            server.stop()

    def test_usage_states_the_threads_and_the_limits(self):
        """The usage of server_query and pattern (review of 2026-10-07): the mask built at
        load takes its threads from --threads-each on a server, from -p in the CLI (M1-04); a
        server's /pattern cap stays under the content timeout (R1-05); the finalisation
        reserve is a floor the answer's estimated writing time adds to, at stated rates (R1-04,
        X-EFFICIENCY-04)."""
        for command, needles in (
                ('server_query', ['with --threads-each threads', 'at most 899000',
                                  'longer by the estimated time to write what the answer holds',
                                  '--pattern-delivery-build-mbps', '--pattern-delivery-compress-mbps']),
                ('pattern', ['with -p threads',
                             'longer by the estimated time to write what the answer holds',
                             '--pattern-delivery-build-mbps', '--pattern-delivery-compress-mbps'])):
            res = subprocess.run(shlex.split(METAGRAPH) + [command], stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
            text = (res.stdout + res.stderr).decode()
            for needle in needles:
                self.assertIn(needle, text, command)

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


class TestPatternFixtureBodies(unittest.TestCase):
    """The frozen fixture bodies against the SPEC, with no server and no index, so that CI runs
    it (review of 2026-10-07, T2-01: api/python/tests/test_pattern_fixtures.py ran in no
    workflow); and the blanking of pattern_fixtures.py --check (C2-02)."""

    @staticmethod
    def load(path, name):
        import importlib.util
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    def test_the_bodies_are_the_spec(self):
        if not os.path.isfile(FIXTURE_VALIDATOR):
            self.skipTest('api/python/tests/test_pattern_fixtures.py is not in this checkout')
        module = self.load(FIXTURE_VALIDATOR, 'pattern_fixture_bodies')
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(module.TestPatternFixtures)
        result = unittest.TestResult()
        suite.run(result)
        problems = [f'{t.id()}: {e[-1500:]}' for t, e in result.failures + result.errors]
        self.assertEqual([], problems)
        self.assertGreater(result.testsRun, 5)
        self.assertEqual([], [t.id() for t, _ in result.skipped])

    def test_check_blanks_only_where_the_clock_stopped(self):
        """--check blanks the counts and work of the first time_limited pattern only: the ones
        after it, stopped by the same budget, are compared (a regression that let them run
        or count would otherwise pass --check)."""
        if not os.path.isfile(FIXTURES_SCRIPT):
            self.skipTest('the fixture script is not in this checkout')
        script = self.load(FIXTURES_SCRIPT, 'pattern_fixtures_script')
        path = os.path.join(os.path.dirname(FIXTURE_VALIDATOR), 'data', 'traverse', 'pattern',
                            'deadline', 'answer.json')
        with open(path) as f:
            stored = json.load(f)
        self.assertEqual(['time_limited', 'time_limited'],
                         [e['determinism'] for e in stored['patterns']])
        # the first pattern got as far as the clock let it: blanked
        varied = json.loads(json.dumps(stored))
        varied['patterns'][0]['work']['steps'] += 4096
        varied['patterns'][0]['timing']['elapsed_ms'] += 1
        self.assertEqual(script.dumps(script.blanked(stored)), script.dumps(script.blanked(varied)))
        # the next one ran or counted after the stop: seen
        for mutate in (lambda e: e['work'].update(steps=50, ranges_visited=50),
                       lambda e: e['counts']['contexts'].update(relation='at_least', value=0)):
            regressed = json.loads(json.dumps(stored))
            mutate(regressed['patterns'][1])
            self.assertNotEqual(script.dumps(script.blanked(stored)),
                                script.dumps(script.blanked(regressed)))


@unittest.skipIf(PROTEIN_MODE, "pattern search is DNA only")
@unittest.skipUnless(_supports_pattern(), "`metagraph pattern` is not available in this build")
@guard(os.path.isfile(os.path.join(MINI_DIR, MINI_GRAPH)) and os.path.isfile(FIXTURES_SCRIPT),
       "the mini index (scripts/traversal/build_mini_refseq.sh; not in CI) or the fixture "
       "script is not there")
class TestPatternFixtures(unittest.TestCase):
    """The frozen fixture bodies the search service tests against
    (api/python/tests/data/traverse/pattern/, SPEC-pattern-search.md §11) are what this binary
    answers: a change of /pattern or of either capabilities route that alters one fails here,
    not first in the service's tests -- on a machine with the mini index, which CI does not
    have (review of 2026-10-07, T2-01). (Review of milestone 1b: the bodies had been generated by
    the milestone-1 binary and missed /capabilities' new `resolve` block and the new
    mask_required message, and nothing that ran a server compared them.)"""

    def test_fixtures_are_what_this_binary_answers(self):
        binary = shlex.split(METAGRAPH)
        if len(binary) != 1:
            self.skipTest('the binary runs under a prefix: the script starts it directly')
        with tempfile.TemporaryDirectory() as work:
            # --check writes nothing: it builds its copies of the mini index in |work| (the
            # mini itself is only read), starts each server, and compares every body once its
            # run-dependent values are blanked
            res = subprocess.run([sys.executable, FIXTURES_SCRIPT, '--metagraph', binary[0],
                                  '--mini', MINI_DIR, '--work', work, '--check'],
                                 stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        out = res.stdout.decode()
        self.assertEqual(0, res.returncode, out[-4000:])
        self.assertRegex(out, r'OK: \d+ fixtures match')


@unittest.skipIf(PROTEIN_MODE, "the request panel is DNA")
@guard(bool(BASE_BINARY) and os.path.isfile(BASE_BINARY),
       "$METAGRAPH_BASE_BINARY (the binary before the pattern search) is not set (not in CI)")
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

    def test_a_half_closed_client_is_answered_as_before(self):
        """Only /pattern withholds its answers from a client that half-closed (review GPT-2 of
        2026-10-08, finding 4): /search and /align answer such a client byte for byte as
        before, their errors included."""
        panel = [
            ('search', {'FASTA': f'>query\n{self.SEQ}', 'top_labels': 5,
                        'min_exact_match': 0.1}),
            ('search', '{"FASTA": ">query\\nAATAAAGG", "discovery_fraction": 0.1,'),
            ('search', {'discovery_fraction': 0.1}),
            ('align', {'FASTA': '>query0\nTCGATCGA', 'min_exact_match': 0}),
            ('align', '{'),
        ]
        for route, payload in panel:
            a = half_closed_request(self.base.port, route, payload)
            b = half_closed_request(self.new.port, route, payload)
            self.assertTrue(a.startswith(b'HTTP/1.1 '), (route, payload, a))
            self.assertEqual(a, b, (route, payload))

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
        # the two stated gains: `pattern` (this route) and `resolve`, the block of /resolve's
        # opt-in bounds.time_budget_ms (milestone 1b, SPEC-labeled-traversal-core.md §4.5): a
        # capabilities block, not a feature, since /resolve is listed already
        self.assertEqual(set(a) | {'pattern', 'resolve'}, set(b))
        self.assertTrue(b['resolve']['time_budget']['accepted'])
        self.assertIsNone(b['resolve']['time_budget']['default'])
        # the probe gains the same two blocks and nothing else
        a = without_instance(self.base.get('traverse/capabilities').json())
        b = without_instance(self.new.get('traverse/capabilities').json())
        new_caps = self.new.get('capabilities').json()
        for block in ('pattern', 'resolve'):
            self.assertEqual(new_caps[block], b.pop(block), block)
        self.assertEqual(a, b)


if __name__ == '__main__':
    unittest.main()
