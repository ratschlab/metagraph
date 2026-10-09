"""
End-to-end tests of POST /pattern and `metagraph pattern` (docs/DESIGN-pattern-search.md): the
counts of the graph contexts of exact DNA and IUPAC patterns, their label-free extraction
(output.labels "none": k-mer, instance, offset, strand, node and row ids, no annotation read),
and their labels (output.labels "all": the columns carrying each context and, on the mini index
with its .seqs, every placed occurrence: record, seq_id, 1-based position, strand), checked
against PYTHON ORACLES over the very records each index was built from — never against the
engine itself.

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
    annotation and its .seqs. It was built without the dummy-edge mask, so a COPY of its graph
    in a temporary directory is given the mask by `metagraph transform --mask-dummy` (the
    one-time step for staging) and served with the original annotation, which the mask leaves
    valid, for exact counts; build/mini_refseq itself is never written to. The masks made the
    other ways are checked against it: a rebuild with --mask-dummy (the same .edgemask, the same
    answers) and a server started with --pattern-build-mask on the unmasked graph (the same
    answers, mask: built_at_load). The graph as built is served too: its counts upper bounds
    with estimates, checked against the masked copy's exact counts and an oracle of the source
    dummies, lists exact; once with the server's default --pattern-max-checked-entries (a
    pattern with at most 50 unchecked candidates has each tested, its counts exact) and once
    with 0 (bounds for every pattern with unchecked candidates). Skipped when the mini index
    is not built.
  - TestPatternSynthetic: random records, BASIC, CANONICAL and PRIMARY graphs at k = 15 (and a
    multi-graph server), always run.
  - TestPatternFixtureBodies: the frozen fixture bodies (api/python/tests/data/traverse/pattern)
    against the SPEC (api/python/tests/test_pattern_fixtures.py), and the blanking of
    pattern_fixtures.py --check; no server, no index: always run.
  - TestPatternFixtures: pattern_fixtures.py --check, the fixture bodies against what this
    binary answers on the mini index; skipped without the mini index.
  - TestPatternRegression: /search and /align answer byte for byte as the base binary's
    ($METAGRAPH_BASE_BINARY, an older build) on test_api.py's fixture and requests; skipped
    without it.

Without the mini index or a base binary TestPatternMini, TestPatternFixtures and
TestPatternRegression are skipped, which unittest counts as passed. With
$METAGRAPH_REQUIRE_GUARDS=1 a missing input fails those classes instead of skipping them, so
that a run meant to check the guarantees cannot pass without them: CI's mini-index job builds
the index and runs TestPatternMini and TestPatternFixtures so (it has no base binary).
"""

import filecmp
import glob
import itertools
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

# the server's defaults (--pattern-* flags; §5.3): each cap is its field's maximum and default,
# but time_budget_ms defaults to 60 s under a 600 s cap
DEFAULT_CAPS = {'max_contexts': 10000, 'max_anchors': 1000, 'max_steps': 100000000,
                'time_budget_ms': 600000, 'min_information_bits': 24, 'max_patterns': 16,
                # output.labels "all"
                'max_labels_per_anchor': 64, 'max_annotation_work': 100000000,
                'max_memory_mb': 256, 'max_labels': 1000, 'max_occurrences_per_label': 16,
                # long_search "paths"
                'max_paths': 1000,
                # no request field: the unchecked candidates a pattern on a graph without its
                # mask may have for each to be tested
                'max_checked_entries': 50,
                # a predicate (SPEC §19): the raw contexts a selection may test, its work per
                # request, and the names a predicate may list (no request field)
                'max_predicate_contexts': 100000, 'max_predicate_work': 100000000,
                'max_predicate_labels': 10000}
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


# peptides: the test's own copy of NCBI's genetic codes used here (gc.prt 4.6, ncbieaa: the
# residue of each codon in TCAG order), never the server's tables
GENETIC_CODES = {
    1: 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    2: 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG',
    11: 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    # Karyorelict and Blastocrithidia nuclear: their stop codons code a residue unless in
    # context (TAA, TAG Q and TGA W in 27; TAA, TAG E in 31): no unconditional stop
    27: 'FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    31: 'FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
}
CODONS_TCAG = [a + b + c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG']


def translate(seq, table):
    """The residues of |seq| read from its first base (an incomplete last codon dropped);
    '?' for a codon with a base other than A, C, G, T."""
    code = dict(zip(CODONS_TCAG, GENETIC_CODES[table]))
    return ''.join(code.get(seq[i:i + 3], '?') for i in range(0, len(seq) - 2, 3))


def residue_admits(residue, aa):
    """A peptide's residue (X, B, Z, J the ambiguity codes) admits the amino acid |aa|; a stop
    is admitted by the stop '*' only, which admits nothing else."""
    if aa == '?':
        return False
    if residue == '*' or aa == '*':
        return residue == aa
    return {'X': True, 'B': aa in 'DN', 'Z': aa in 'EQ', 'J': aa in 'IL'}.get(residue,
                                                                         aa == residue)


def peptide_codons(residue, table):
    return {c for c, aa in zip(CODONS_TCAG, GENETIC_CODES[table]) if residue_admits(residue, aa)}


def peptide_prefix(residues, table, s, reverse):
    """|s| is the prefix of an instance of the oriented peptide: each complete codon one of its
    residue's, a partial last codon the prefix of one (rc(P): codon i the reverse complement
    of a codon of residue m - 1 - i)."""
    m = len(residues)
    if len(s) > 3 * m:
        return False
    for i in range(0, len(s), 3):
        part = s[i:i + 3]
        r = residues[m - 1 - i // 3] if reverse else residues[i // 3]
        codons = peptide_codons(r, table)
        if not any((revcomp(c) if reverse else c).startswith(part) for c in codons):
            return False
    return True


def six_frames_reference(seq, residues, table):
    """{(0-based start, strand)}: every match of the peptide in the six frames of |seq| (its
    three frames and the three of its reverse complement), the start on |seq|'s + strand. The
    plain definition, residue by residue: six_frames must answer as it does."""
    m, L = len(residues), 3 * len(residues)
    out = set()
    for strand, t in (('+', seq), ('-', revcomp(seq))):
        for frame in range(3):
            aa = translate(t[frame:], table)
            for j in range(len(aa) - m + 1):
                if all(residue_admits(residues[x], aa[j + x]) for x in range(m)):
                    start = frame + 3 * j
                    out.add((start if strand == '+' else len(t) - start - L, strand))
    return out


# the six translations of a record, per (record, table): the mini's oracles scan the same 9 Mbp
# for every peptide, and translating them again for each would dominate the module's runtime
_FRAMES = {}


def translated_frames(seq, table):
    """[(strand, frame, residues)]: the translations of |seq|'s three frames and of the three of
    its reverse complement, as translate() reads them."""
    key = (seq, table)
    if key not in _FRAMES:
        code = dict(zip(CODONS_TCAG, GENETIC_CODES[table]))
        frames = []
        for strand, t in (('+', seq), ('-', revcomp(seq))):
            for frame in range(3):
                codons = re.findall('...', t[frame:])
                frames.append((strand, frame,
                               ''.join(map(code.get, codons, itertools.repeat('?')))))
        _FRAMES[key] = frames
    return _FRAMES[key]


def peptide_regex(residues):
    """A lookahead matching each start of the peptide in a translation, a residue the class of
    the amino acids residue_admits: '*' only the stop, X anything but a stop or a '?'."""
    classes = {'*': r'\*', 'X': '[^*?]', 'B': '[DN]', 'Z': '[EQ]', 'J': '[IL]'}
    return re.compile('(?=' + ''.join(classes.get(r, re.escape(r)) for r in residues) + ')')


def six_frames(seq, residues, table):
    """six_frames_reference(seq, residues, table), from the cached translations."""
    L = 3 * len(residues)
    regex = peptide_regex(residues)
    out = set()
    for strand, frame, aa in translated_frames(seq, table):
        for match in regex.finditer(aa):
            start = frame + 3 * match.start()
            out.add((start if strand == '+' else len(seq) - start - L, strand))
    return out


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

    def has_kmer(self, kmer, canonical=False):
        """Whether a k-mer of the graph built from these records spells |kmer|: one of the
        islands holds it (a DNA4 graph keeps every k-mer of an island and nothing else), or its
        reverse complement on a graph holding both orientations."""
        key = ('kmer', kmer, canonical)
        if key not in self._cache:
            if not hasattr(self, '_joined'):
                self._joined = '|'.join(island for _, _, island in self.islands)
            self._cache[key] = kmer in self._joined or (canonical and revcomp(kmer)
                                                        in self._joined)
        return self._cache[key]

    def paths(self, pattern, k, strands='both', canonical=False):
        """{(strand or orientation, sequence)}: the graph paths of a pattern longer than k
        (§4.2) — every string of its length that instantiates the oriented pattern and whose
        every k-window is a k-mer of the graph (has_kmer), found by extending each instance of
        its first k positions one base at a time (the graph-walk oracle over the records'
        k-mers: a path need not lie in one record)."""
        assert len(pattern) > k
        out = set()
        for strand, q in self.oriented(pattern, strands, not canonical):
            def extend(s):
                if len(s) == len(q):
                    out.add((strand, s))
                    return
                for b in IUPAC[q[len(s)]]:
                    if self.has_kmer(s[len(s) - k + 1:] + b, canonical):
                        extend(s + b)
            # the instances of the anchor window that are k-mers of the graph
            starts = {''}
            for c in q[:k]:
                starts = {x + b for x in starts for b in IUPAC[c]}
            for x in sorted(starts):
                if self.has_kmer(x, canonical):
                    extend(x)
        return out

    def branchings(self, pattern, k, strands='both', canonical=False):
        """The walks of k to L - 1 bases of the graph-walk oracle of paths (anchors included)
        that two or more k-mers of the graph extend at the next position: where the extension
        branches (work.extension_branches)."""
        assert len(pattern) > k
        count = 0
        for _, q in self.oriented(pattern, strands, not canonical):
            def extend(s):
                nonlocal count
                if len(s) == len(q):
                    return
                nxt = [b for b in IUPAC[q[len(s)]]
                       if self.has_kmer(s[len(s) - k + 1:] + b, canonical)]
                count += len(nxt) > 1
                for b in nxt:
                    extend(s + b)
            starts = {''}
            for c in q[:k]:
                starts = {x + b for x in starts for b in IUPAC[c]}
            for x in sorted(starts):
                if self.has_kmer(x, canonical):
                    extend(x)
        return count

    def source_dummies(self, k, canonical=False):
        """The source dummy k-mers of the BOSS graph of these records (what a graph without its
        dummy-edge mask counts besides its k-mers): '$' * j + x[:k - j] for 1 <= j < k and
        every node x (a (k - 1)-mer) that starts a k-mer and that no k-mer enters, i.e. the
        first k - 1 bases of a sequence (island, or its reverse complement on a graph holding
        both orientations) holding a k-mer, found at no later position of any sequence. The
        main dummy edge ($ * k: its W is $) and the sink dummies (W = $) are not among them: no
        pattern base matches $."""
        key = ('dummies', k, canonical)
        if key not in self._cache:
            seqs = [s for s in self.sequences(canonical)]
            starts = {s[:k - 1] for s in seqs if len(s) >= k}
            entered = {x for x in starts if any(s.find(x, 1) != -1 for s in seqs)}
            self._cache[key] = {'$' * j + x[:k - j] for x in starts - entered
                                for j in range(1, k)}
        return self._cache[key]

    def dummy_contexts(self, pattern, k, scope='any_offset', strands='both', canonical=False):
        """{(strand or orientation, dummy k-mer, offset)}: the contexts a graph without its mask
        counts in source dummies (the oriented pattern on bases only, never on a $)."""
        L = len(pattern)
        out = set()
        for strand, q in self.oriented(pattern, strands, not canonical):
            rx = pattern_regex(q)
            for d in self.source_dummies(k, canonical):
                for m in rx.finditer(d):
                    if scope == 'suffix' and m.start() != k - L:
                        continue
                    out.add((strand, d, m.start()))
        return out

    def names_with_kmer(self, kmer):
        return {name for _, name, island in self.islands if kmer in island}

    def sources_with_kmer(self, kmer):
        return {source for source, _, island in self.islands if kmer in island}


# a 10-mer whose suffix-scope count on the mini needs a mask scan (it starts a record's first
# k-mer): discovery takes 20 range steps, the scan one more, so max_steps 20 stops in the scan
# with bounds
SCAN_10 = 'ATGCCGGTGA'
SCAN_STEPS = 20
# a 31-mer of 12 bases around a run of 19 Ns (24 bits): about 2.5e7 range steps on the mini, a
# discovery a work deadline of a second stops
HEAVY31 = 'GG' + 'N' * 19 + 'TTGGCGATCT'

# the note of an entry whose counts carry an estimate (a graph without its dummy-edge mask)
NOTE_ESTIMATE = 'estimate_sampled_dummy_fraction'


def guard(condition, reason):
    """unittest.skipUnless(condition, reason); with $METAGRAPH_REQUIRE_GUARDS=1 a class whose
    input is missing fails instead."""
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
    request. The client still reads until the server closes. The request is held back until
    the FIN goes with it where the platform can do so (TCP_CORK on Linux, TCP_NOPUSH on
    macOS), so that the server sees the FIN as soon as it has read the request, not after a
    race."""
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
       "the mini index is not built (scripts/traversal/build_mini_refseq.sh)")
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
        TestingBase._run_command(f'{METAGRAPH} transform --mask-dummy -p 2 {cls.graph}',
                                 'Mask the copy of the mini graph')
        assert os.path.isfile(cls.graph[:-len('.dbg')] + '.edgemask')
        assert filecmp.cmp(cls.graph, os.path.join(MINI_DIR, MINI_GRAPH), shallow=False)
        # the mini graph as built, without a mask whatever build/mini_refseq holds: a symlink in
        # a directory of its own (the loader reads the mask beside the path it is given)
        cls.unmasked_graph = cls._mini_copy('unmasked', copy_graph=False)
        cls.unmasked_anno = os.path.join(os.path.dirname(cls.unmasked_graph), MINI_ANNO)

        cls.server = Server(METAGRAPH, ['-i', cls.graph, '-a', cls.anno],
                            os.path.join(d, 'server.log'))
        cls.stats = cls.server.get('stats').json()
        # the mini graph as built, without its mask: served, its counts upper bounds with
        # estimates, its lists exact
        cls.unmasked_server = Server(METAGRAPH, ['-i', cls.unmasked_graph,
                                                 '-a', cls.unmasked_anno],
                                     os.path.join(d, 'server_unmasked.log'))
        # the same with the candidate check off: every count with unchecked candidates is
        # bounds, as any count above the limit (the tests of a graph without its mask)
        cls.unchecked_server = Server(METAGRAPH, ['-i', cls.unmasked_graph,
                                                  '-a', cls.unmasked_anno,
                                                  '--pattern-max-checked-entries', '0'],
                                      os.path.join(d, 'server_unchecked.log'))
        cls._choose_patterns()

    @classmethod
    def tearDownClass(cls):
        cls.server.stop()
        cls.unmasked_server.stop()
        cls.unchecked_server.stop()
        if cls._no_map is not None:
            cls._no_map.stop()
        cls.tempdir.cleanup()

    _no_map = None

    @classmethod
    def no_map_server(cls):
        """The masked copy served with --no-coord-mapping (coordinates without the record
        mapping), started once for the tests that read it."""
        if cls._no_map is None:
            cls._no_map = Server(METAGRAPH, ['-i', cls.graph, '-a', cls.anno,
                                             '--no-coord-mapping'],
                                 os.path.join(cls.tempdir.name, 'server_no_map.log'))
        return cls._no_map

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
        # patterns at the start of an island no k-mer enters: its source dummies hold them
        # ($^j x[:k - j]), which a graph without its mask counts in its upper bounds
        dummies = cls.records.source_dummies(cls.k)
        starts = sorted(d[1:] for d in dummies if len(d.lstrip('$')) == cls.k - 1)
        cls.start16 = starts[0][:16]
        cls.start14 = starts[1][:14]
        cls.start40 = next(island[:40] for island in islands if island[:cls.k - 1] == starts[0])

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
        # a budget under the cap is kept as requested
        out = self.pattern(self.server, {'patterns': [{'dna': self.p16}], 'mode': 'count',
                                         'time_budget_ms': 1234.5})
        self.assertEqual(1234.5, out['limits']['time_budget_ms'])
        self.assertEqual([], out['limits']['clamped'])

    def test_capabilities(self):
        caps = self.server.get('capabilities').json()
        self.assertIn('pattern', caps['features'])
        self.assertEqual('POST /pattern', caps['routes']['pattern'])
        self.assertEqual('GET /pattern/capabilities', caps['routes']['pattern_capabilities'])
        p = caps['pattern']
        expected = {
            'pattern_contract_version': 1, 'available': True, 'unavailable_reason': None,
            'modes': ['count', 'all_or_count', 'partial'], 'default_mode': 'all_or_count',
            'projections': ['none', 'all', 'predicate_only'], 'default_projection': 'none',
            'projections_later_increment': [], 'default_occurrences': True,
            'kinds': ['dna', 'iupac', 'protein'], 'kinds_later_increment': [],
            # peptides: the residues (the stop '*' among them), the genetic codes and the
            # default
            'protein_residues': list('ACDEFGHIKLMNPQRSTVWYXBZJ*'),
            'genetic_codes': [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25,
                              26, 27, 28, 29, 30, 31, 32, 33],
            'default_genetic_code': 1,
            'scopes': ['suffix', 'any_offset'],
            'default_scope': 'any_offset', 'long_patterns': 'anchors_counted',
            # the paths of long patterns, opt-in
            'long_search': ['anchors', 'paths'], 'default_long_search': 'anchors',
            'strands': ['both', 'forward', 'reverse'], 'graph_mode': 'basic', 'k': self.k,
            'alphabet': '$ACGT', 'strand_stated': True, 'mask': 'file',
            # the masked graph counts exactly, no dummy fraction
            'counting': 'exact', 'dummy_fraction': None,
            'graph_cleaned': 'unknown', 'records_shorter_than_k': 'not_indexed',
            # in_ram accepted (as /search's): the budgets start after a load
            'resident_only': False, 'in_ram': 'accepted', 'caps': DEFAULT_CAPS,
            'default_time_budget_ms': DEFAULT_TIME_MS,
            'finalize_reserve_ms': DEFAULT_FINALIZE_MS,
            # the delivery rates as numbers (MB/s), the prose rules SPEC references
            'delivery_mbps': {'build': 10, 'compress': 50},
            # predicates: the operators, the strands and the access of a budget-aware index
            'predicate': {'operators': ['any', 'all', 'none', 'at_least', 'and', 'or', 'not'],
                          'strands': ['either', 'context'], 'access': 'rows'},
        }
        self.assertEqual(expected, {x: p[x] for x in expected})
        self.assertEqual(['any_offset'], p['scopes_by_graph_mode']['primary'])
        # references to the sections that state the rules (§4.5 classifies every cap)
        self.assertEqual('SPEC-pattern-search.md section 4.5', p['caps_rule'])
        self.assertEqual('SPEC-pattern-search.md section 12.2', p['protein_rule'])
        # the same block on the probe a service reads, with the route of the full block, and
        # the full block on its own route (SPEC §23), compressed when asked for
        probe = self.server.get('traverse/capabilities').json()
        self.assertEqual(dict(p, details='GET /pattern/capabilities'), probe['pattern'])
        ret = requests.get(url=self.server.url('pattern/capabilities'),
                           headers={'Accept-Encoding': 'gzip'}, timeout=60)
        self.assertEqual((200, 'gzip'), (ret.status_code, ret.headers.get('Content-Encoding')))
        self.assertEqual(p, ret.json())
        # a parameter that would select a graph is refused on a single-graph server; any other
        # is ignored
        for query in ('?graph=G1', '?graph_path=x', '?other=1&graph=G1'):
            ret = self.server.get('pattern/capabilities' + query)
            self.assertEqual(400, ret.status_code, query)
            self.assertEqual({'error': "Bad request: this server hosts a single graph; remove "
                                       "the 'graph' / 'graph_path' parameter"}, ret.json())
        self.assertEqual(p, self.server.get('pattern/capabilities?other=1').json())
        # no other method
        self.assertEqual(404, self.server.post('pattern/capabilities', {}).status_code)
        if self.anno.endswith(MINI_ANNO):
            # coordinates and record mapping, a budgeted row-diff annotation
            self.assertEqual(('record', 'record_verified', 'budgeted'),
                             (p['placement'], p['support'], p['annotation']))

    # ------------------------------------------------------------ the label-free retrieval

    def test_retrieval_all_or_count(self):
        patterns = [self.p16, self.iupac16, self.pal12, self.absent16]
        # mode and projection omitted: all_or_count with labels "none"
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
        """max_steps reached in a mask scan: the count is bounds, value = lower <= the oracle's
        count <= upper, the parts summed as SPEC §7.4 says, and the stop names the scan."""
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
        """The work stops earlier by the estimated time to write what the answer holds, so
        that a request whose earlier patterns buffered many results and whose last one runs to
        the work deadline still answers within its budget, with its counts and a time stop,
        instead of 503 with all of its work lost. Fifteen patterns of 19,283 results each (a
        server whose max_contexts admits them all) and a last one whose discovery outlasts the
        budget: with the 250 ms reserve alone the work would run to 1,750 ms and the ~33 MB of
        results could not be written by 2,000 ms (503 in every run)."""
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
        """/resolve's 503 end to end, on the mini server because it needs a query whose k-mer
        mapping, which the deadline cannot interrupt, outlasts the whole budget: 503 {error,
        code: deadline}, uncompressed although the client accepts gzip, without Retry-After
        (that is the loading 503's), and the CLI writes the same body and exits 1."""
        # the whole 7 Mbp record: its mapping takes well over the 250 ms reserve on any host
        # (2 Mbp mapped within it on a quiet machine)
        seq = ''.join(s for _, s in read_fasta(os.path.join(MINI_DIR, 'fasta', '287.fa')))
        body = {'sequence': seq, 'discover': {'max_labels': 10},
                'bounds': {'time_budget_ms': 250.001, 'max_query_bp': len(seq)}}
        expected = {'error': 'resolve: the answer could not be built and written within '
                             'bounds.time_budget_ms (250.001 ms, the finalisation reserve of '
                             '250 ms included): nothing partial is sent',
                    'code': 'deadline'}
        ret = requests.post(self.server.url('resolve'), data=json.dumps(body),
                            headers={'Accept-Encoding': 'gzip'}, timeout=300)
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
        # both windows of this 40-mer carry 62 bits (min_anchor_information_bits: the least
        # searched window, the floor's operand)
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

    # ------------------------------------------------------------ paths

    def long_patterns(self):
        """Patterns of 35-50 bases taken from the mini's records (so present): the blaNDM-1
        region around NDM-F (in several taxa), windows of the two longest islands, an IUPAC
        40-mer with two ambiguous positions; and an absent 40-mer (no anchor)."""
        ndm = next(seq for recs in self.columns.values() for _, seq in recs
                   if self.NDM_F in seq.upper()).upper()
        i = ndm.find(self.NDM_F)
        islands = sorted((island for _, _, island in self.records.islands),
                         key=lambda s: (-len(s), s[:64]))
        w = list(islands[0][3000:3040])
        w[2] = 'R' if w[2] in 'AG' else 'Y'
        w[36] = 'N'
        return [ndm[i - 10:i + 35], self.p40, islands[1][7000:7050], islands[0][6000:6035],
                ''.join(w), self.absent40]

    def path_columns(self, path):
        """The record-scan oracle of a path (§4.3): the columns carrying it (every k-mer of it
        in one of their records) and, for those whose records hold it whole (record_verified),
        {column: {(seq_id, 1-based start)}}."""
        k = self.k
        carriers = {c for c, r in self.column_records.items()
                    if all(r.has_kmer(path[i:i + k]) for i in range(len(path) - k + 1))}
        verified = {}
        for c in carriers:
            occ = {(i, m.start() + 1) for i, (_, seq) in enumerate(self.columns[c])
                   for m in re.finditer('(?=' + path + ')', seq.upper())}
            if occ:
                verified[c] = occ
        return carriers, verified

    def assertPaths(self, entry, pattern, labels=True, require_verified=False):
        """An answer of long_search "paths" against both oracles: the graph-walk oracle's paths
        (Records.paths) with their counts, each path result with its sequence, anchor k-mer and
        node and row ids; with labels, each path's columns those carrying it (or verifying it,
        with require_support), each with its support and its occurrences of the whole path, the
        summaries exact, and the placed occurrences the FASTA scan's (placed_oracle)."""
        k, L = self.k, len(pattern)
        n = L - k + 1
        expected = self.records.paths(pattern, k)
        anchors = self.records.anchors(pattern, k)
        self.assertEqual('long', entry['scope'])
        self.assertNotIn('paths_later_increment', entry['notes'])
        self.assertCount({x: entry['counts']['anchors'][x] for x in ('value', 'relation', 'unit')},
                         len(anchors), unit='anchors')
        c = entry['counts']['paths']
        self.assertCount({x: c[x] for x in ('value', 'relation', 'unit')}, len(expected),
                         unit='paths')
        self.assertEqual('completed' if anchors else 'no_anchors', c['extension'])
        for s in entry['strands']:
            self.assertCount(c['by_strand']['both' if s == '=' else s],
                             sum(1 for t, _ in expected if t == s), unit='paths')
        # every path enters its n - 1 prefixes beyond the anchor (shared prefixes once)
        if expected:
            self.assertGreaterEqual(c['candidates_examined'], max(len(expected), n - 1))
        self.assertIn('extension_edges', entry['work'])
        self.assertIn('extension_ms', entry['timing'])
        # the anchors whose extension began (every anchor: it completed) and the walks where
        # it branched, against the graph-walk oracle
        self.assertEqual((len(anchors), self.records.branchings(pattern, k)),
                         (entry['work']['extension_anchors'],
                          entry['work']['extension_branches']))
        labelled = labels and 'results' in entry
        for field in ('annotation_rows_distinct', 'verification_steps'):
            self.assertEqual(labelled, field in entry['work'], field)
        for field in ('label_intersection_ms', 'verification_ms'):
            self.assertEqual(labelled, field in entry['timing'], field)
        if 'results' not in entry:
            return
        self.assertCompleteRetrieval(entry)
        got = {(r['strand'], r['sequence']) for r in entry['results']}
        self.assertEqual(expected, got)
        self.assertEqual(len(expected), entry['returned'])
        order = []
        for r in entry['results']:
            keys = {'sequence', 'anchor_kmer', 'instance', 'offset', 'strand', 'nodes', 'rows'}
            if labels:
                keys |= {'support', 'labels_status', 'labels_total', 'labels'}
                if require_verified:
                    keys.add('labels_excluded_unverified')
            self.assertEqual(keys, set(r))
            self.assertEqual((r['sequence'], r['sequence'][:k], 0),
                             (r['instance'], r['anchor_kmer'], r['offset']))
            self.assertEqual(n, len(r['nodes']))
            self.assertEqual([x - 1 for x in r['nodes']], r['rows'])
            order.append((r['nodes'][0], self.ORIENTATION_RANK[r['strand']], r['sequence']))
        self.assertEqual(sorted(order), order)
        if not labels:
            return

        # the rows read: each distinct k-mer of the paths once (BASIC: a k-mer is its row)
        self.assertEqual(len({s[i:i + k] for _, s in expected for i in range(n)}),
                         entry['work']['annotation_rows_distinct'])
        unions, paths_of, verified_of = {}, {}, {}
        carried = 0
        for r in entry['results']:
            s, strand = r['sequence'], r['strand']
            carriers, verified = self.path_columns(s)
            carried += len(carriers)
            listed = set(verified) if require_verified else carriers
            self.assertEqual(('complete', len(carriers)), (r['labels_status'], r['labels_total']))
            if require_verified:
                self.assertEqual(len(carriers) - len(verified), r['labels_excluded_unverified'])
            self.assertEqual(listed, {x['column'] for x in r['labels']}, s)
            support = (None if not listed else 'record_verified' if listed <= set(verified)
                       else 'label_intersection' if not set(verified) & listed else 'mixed')
            self.assertEqual(support, r['support'])
            for label in r['labels']:
                col = label['column']
                occ = verified.get(col, set())
                self.assertEqual('record_verified' if occ else 'label_intersection',
                                 label['support'], (s, col))
                self.assertCount(label['occurrences'], len(occ), unit='placed_occurrences')
                got_occ = set()
                for o in label['occurrence_list']:
                    name, seq = self.columns[col][o['seq_id']]
                    self.assertEqual((name, len(seq), strand),
                                     (o['record'], o['nt_length'], o['strand']))
                    a, b = (int(x) for x in o['nt_coords'].split('-'))
                    self.assertEqual(a + L - 1, b)
                    self.assertEqual(s, seq[a - 1:b].upper())
                    got_occ.add((o['seq_id'], a))
                    unions.setdefault(col, set()).add((o['seq_id'], a, strand))
                self.assertEqual(occ, got_occ, (s, col))
                paths_of[col] = paths_of.get(col, 0) + 1
                verified_of[col] = verified_of.get(col, 0) + bool(occ)
        # by_label: (paths desc, column asc), exact
        keys = [(-b['paths']['value'], b['column']) for b in entry['by_label']]
        self.assertEqual(sorted(keys), keys)
        self.assertEqual(set(paths_of), {b['column'] for b in entry['by_label']})
        for b in entry['by_label']:
            col = b['column']
            self.assertEqual({'graph', 'column', 'paths', 'paths_record_verified', 'occurrences'},
                             set(b))
            self.assertCount(b['paths'], paths_of[col], unit='paths')
            self.assertCount(b['paths_record_verified'], verified_of[col], unit='paths')
            self.assertCount(b['occurrences'], len(unions.get(col, ())), unit='placed_occurrences')
        labels_count = entry['counts']['labels']
        self.assertCount({x: labels_count[x] for x in ('value', 'relation', 'unit')},
                         len(paths_of), unit='labels')
        rv = sum(1 for col in paths_of if verified_of[col])
        self.assertCount(labels_count['by_support']['record_verified'], rv, unit='labels')
        self.assertCount(labels_count['by_support']['label_intersection'], len(paths_of) - rv,
                         unit='labels')
        placed = {(col, i, a, st) for col, x in unions.items() for i, a, st in x}
        self.assertCount(entry['counts']['occurrences'], len(placed), unit='placed_occurrences')
        # every occurrence of the oriented pattern in the FASTA records, each once
        self.assertEqual(self.placed_oracle(pattern), placed)
        # the verification (record placement) looks up the rows of the n k-mers of every path
        # for every label carrying it, one unit each, before its joins and record steps
        self.assertGreaterEqual(entry['work']['verification_steps'], n * carried)

    def paths_request(self, patterns, **kw):
        return dict({'patterns': [{'iupac' if set(p) - set('ACGT') else 'dna': p}
                                  for p in patterns], 'long_search': 'paths'}, **kw)

    def test_long_paths_against_the_fasta(self):
        """Patterns of 35-50 bases (long_search "paths"): every path with every column carrying
        it and, where one record of the column holds the whole path, its support
        record_verified and its occurrences; the occurrences the FASTA scan's."""
        patterns = self.long_patterns()
        out = self.pattern(self.server, self.paths_request(patterns, output={'labels': 'all'}))
        self.assertEqual(('paths', 1000, 'label_intersection'),
                         (out['limits']['long_search'], out['limits']['max_paths'],
                          out['limits']['require_support']))
        for entry, p in zip(out['patterns'], patterns):
            with self.subTest(pattern=p):
                self.assertPaths(entry, p)
                self.assertEqual(('record', 'budgeted'), (entry['placement'], entry['annotation']))
        ndm = out['patterns'][0]
        self.assertEqual(2, ndm['counts']['paths']['value'])
        self.assertGreater(len(ndm['by_label']), 1)
        self.assertGreater(sum(len(self.records.paths(p, self.k)) for p in patterns), 5)
        # the absent one: no anchor, no path, the empty answer complete
        absent = out['patterns'][-1]
        self.assertEqual(('no_anchors', 0, []), (absent['counts']['paths']['extension'],
                                                 absent['counts']['paths']['value'],
                                                 absent['results']))
        # require_support "record_verified": the verified columns only
        req = self.paths_request(patterns, output={'labels': 'all'},
                                 require_support='record_verified')
        verified = self.pattern(self.server, req)
        for entry, p in zip(verified['patterns'], patterns):
            with self.subTest(pattern=p, require_support='record_verified'):
                self.assertPaths(entry, p, require_verified=True)
                self.assertEqual('exact', entry['labels_excluded_unverified']['relation'])
        # deterministic
        self.assertEqual(untimed(verified), untimed(self.pattern(self.server, req)))

    def test_long_paths_label_free_and_count(self):
        patterns = self.long_patterns()
        count = self.pattern(self.server, self.paths_request(patterns, mode='count'))
        none = self.pattern(self.server, self.paths_request(patterns, output={'labels': 'none'}))
        partial = self.pattern(self.server, self.paths_request(patterns, mode='partial'))
        for c, e, q, p in zip(count['patterns'], none['patterns'], partial['patterns'], patterns):
            with self.subTest(pattern=p):
                self.assertPaths(c, p, labels=False)
                self.assertPaths(e, p, labels=False)
                self.assertPaths(q, p, labels=False)
                self.assertNotIn('results', c)
                self.assertFalse(c['retrieval_complete'])
                self.assertEqual(e['counts'], c['counts'])
                self.assertEqual(e['results'], q['results'])
                self.assertUnknownLabels(e)
        # the option is opt-in: without it, or with "anchors", the answer by anchors
        plain = self.pattern(self.server, {'patterns': [{'dna': self.p40}]})
        anchors = self.pattern(self.server, {'patterns': [{'dna': self.p40}],
                                             'long_search': 'anchors'})
        self.assertEqual(untimed(plain), untimed(anchors))
        self.assertIn('paths_later_increment', plain['patterns'][0]['notes'])
        self.assertNotIn('long_search', plain['limits'])

    def test_long_paths_thresholds(self):
        """max_anchors admits the extension, max_paths the release (§4.2, §5.2), on the
        blaNDM-1 region: in records of both orientations, a path on each strand."""
        ndm = self.long_patterns()[0]
        paths = self.records.paths(ndm, self.k)
        self.assertGreater(len(paths), 1)
        entry = self.pattern(self.server, self.paths_request([ndm], max_anchors=0))['patterns'][0]
        self.assertEqual({'reason': 'anchors_above_threshold'}, entry['withheld'])
        self.assertCount({x: entry['counts']['paths'][x] for x in ('value', 'relation', 'unit')},
                         None, 'unknown', 'paths')
        self.assertEqual('not_admitted', entry['counts']['paths']['extension'])
        entry = self.pattern(self.server, self.paths_request([ndm], max_paths=1))['patterns'][0]
        self.assertEqual({'reason': 'count_above_threshold'}, entry['withheld'])
        self.assertEqual(('exact', len(paths)), (entry['counts']['paths']['relation'],
                                                 entry['counts']['paths']['value']))
        full = self.pattern(self.server, self.paths_request([ndm]))['patterns'][0]
        entry = self.pattern(self.server, self.paths_request([ndm], max_paths=1,
                                                             mode='partial'))['patterns'][0]
        self.assertEqual(({'reason': 'max_paths'}, 1), (entry['cut'], entry['returned']))
        self.assertEqual(full['results'][:1], entry['results'])
        self.assertFalse(entry['retrieval_complete'])
        # above the server's cap: lowered and listed
        out = self.pattern(self.server, self.paths_request([ndm], max_paths=10 ** 6))
        self.assertEqual([{'field': 'max_paths', 'requested': 10 ** 6,
                           'effective': DEFAULT_CAPS['max_paths']}], out['limits']['clamped'])

    def test_long_paths_without_record_mapping(self):
        """--no-coord-mapping: the coordinates without the .seqs. A path's labels are the
        columns carrying it, none record_verified (record bounds unknown: a chain of
        consecutive column coordinates may cross from one record into the next); each lists
        its chains (kmer_coord of the first k-mer, offset 0). require_support "record_verified"
        is refused there."""
        server = self.no_map_server()
        p = self.long_patterns()[0]
        entry = self.pattern(server, self.paths_request([p], output={'labels': 'all'})
                             )['patterns'][0]
        ret = server.post('pattern', self.paths_request([p], output={'labels': 'all'},
                                                        require_support='record_verified'))
        self.assertEqual((400, 'support_unavailable'), (ret.status_code, ret.json()['code']))
        self.assertEqual('global', entry['placement'])
        self.assertIn('record_bounds_unknown', entry['notes'])
        self.assertCompleteRetrieval(entry)
        self.assertCount(entry['counts']['occurrences'], None, 'unknown', 'placed_occurrences')
        starts = {}
        for column, recs in self.columns.items():
            s, at = [], 0
            for _, seq in recs:
                s.append(at)
                at += len(seq) - self.k + 1
            starts[column] = s
        for r in entry['results']:
            carriers, verified = self.path_columns(r['sequence'])
            self.assertEqual(carriers, {x['column'] for x in r['labels']})
            for label in r['labels']:
                self.assertEqual('label_intersection', label['support'])
                self.assertNotIn('occurrences', label)
                recs = self.columns[label['column']]
                chains = set()
                for o in label['occurrence_list']:
                    self.assertEqual((0, r['strand']), (o['offset'], o['strand']))
                    c = o['kmer_coord']
                    j = max(i for i, s in enumerate(starts[label['column']]) if s <= c)
                    chains.add((j, c - starts[label['column']][j]))
                # every whole occurrence in a record is a chain (a chain may also cross)
                self.assertLessEqual({(i, a - 1) for i, a in verified.get(label['column'], ())},
                                     chains)

    # ------------------------------------------------------------ peptides

    def peptide_contexts(self, residues, table):
        """{(strand, k-mer, offset)}: the graph contexts of a peptide of at most k bases from the
        six-frame translation of every record (the mini's records hold no N: a record is an
        island, its k-mers the graph's): every k-mer of a record that holds a match."""
        k, L = self.k, 3 * len(residues)
        out = set()
        for recs in self.columns.values():
            for _, seq in recs:
                seq = seq.upper()
                for i, strand in six_frames(seq, residues, table):
                    for s in range(max(0, i + L - k), min(i, len(seq) - k) + 1):
                        out.add((strand, seq[s:s + k], i - s))
        return out

    def peptide_placed(self, residues, table):
        """{(column, seq_id, 1-based start, strand)}: the six-frame matches in every record."""
        return {(column, seq_id, i + 1, strand)
                for column, recs in self.columns.items()
                for seq_id, (_, seq) in enumerate(recs)
                for i, strand in six_frames(seq.upper(), residues, table)}

    def peptide_paths(self, residues, table):
        """{(strand, sequence)}: the graph paths of a peptide longer than k (the graph-walk
        oracle): each anchor -- a k-mer of the records whose bases begin an instance of the
        oriented peptide, found where the six-frame translation matches the residues its whole
        codons cover -- extended one base at a time while the bases stay a prefix of an instance
        and every k-window is a k-mer of the records (Records.has_kmer)."""
        k, m = self.k, len(residues)
        whole = k // 3
        out = set()
        for strand, reverse in (('+', False), ('-', True)):
            # the residues the anchor's whole codons spell, in the record's reading direction
            head = residues[m - whole:] if reverse else residues[:whole]
            anchors = set()
            for recs in self.columns.values():
                for _, seq in recs:
                    seq = seq.upper()
                    for i, st in six_frames(seq, head, table):
                        # forward: P's first codons match on '+' where the anchor starts;
                        # reverse: rc(P)'s first bases are the reverse complement of P's last
                        # codons, which match on '-' where the anchor starts too
                        if st == ('-' if reverse else '+') and i + k <= len(seq):
                            anchors.add(seq[i:i + k])
            anchors = {x for x in anchors if peptide_prefix(residues, table, x, reverse)}

            def extend(s):
                if len(s) == 3 * m:
                    out.add((strand, s))
                    return
                for b in 'ACGT':
                    t = s + b
                    if peptide_prefix(residues, table, t, reverse) \
                            and self.records.has_kmer(t[len(t) - k:]):
                        extend(t)
            for x in sorted(anchors):
                extend(x)
        return out

    def peptides(self):
        """Peptides of the mini's records: NDM-1's first residues (blaNDM-1 on both strands of
        several taxa) within one k-mer (10 residues, 30 bases) with an ambiguity code, and
        windows of the records' translations in their bacterial code (11)."""
        rng = random.Random(1508)
        out = [('MELPNIMHPV', 11), ('MZLPBJMHPV', 1), ('MELPNXMHPV', 11)]
        recs = [seq.upper() for r in self.columns.values() for _, seq in r]
        while len(out) < 9:
            seq = rng.choice(recs)
            i = rng.randrange(len(seq) - 40)
            residues = translate(seq[i:i + 3 * 9], 11)
            if '*' in residues or '?' in residues:
                continue
            out.append((residues, 11))
        return out

    def assertPeptideEntry(self, entry, residues, table):
        self.assertEqual(('protein', residues, 3 * len(residues), len(residues), table),
                         (entry['kind'], entry['pattern'], entry['length'], entry['residues'],
                          entry['genetic_code']))
        # a residue without a codon (a stop '*' in a table without one) counts as one codon
        bits = sum(math.log2(64 / max(1, len(peptide_codons(r, table)))) for r in residues)
        self.assertAlmostEqual(bits, entry['information_bits'], places=9)
        self.assertFalse(entry['palindromic'])

    def test_peptides_against_the_six_frames(self):
        """Peptides within one k-mer: their contexts the six-frame translation's (every k-mer
        of a record holding a match), each with the columns whose records hold its k-mer, and
        the placed occurrences every match of the translation; in count mode the same
        counts."""
        peptides = self.peptides()
        for residues, table in peptides:
            with self.subTest(peptide=residues, table=table):
                request = {'patterns': [{'protein': residues}], 'genetic_code': table}
                count = self.pattern(self.server, dict(request, mode='count'))['patterns'][0]
                entry = self.pattern(self.server, dict(request, output={'labels': 'all'})
                                     )['patterns'][0]
                expected = self.peptide_contexts(residues, table)
                self.assertGreater(len(expected), 0)
                for e in (count, entry):
                    self.assertPeptideEntry(e, residues, table)
                    self.assertCount({x: e['counts']['contexts'][x]
                                      for x in ('value', 'relation', 'unit')}, len(expected))
                self.assertCompleteRetrieval(entry)
                got = {(r['strand'], r['kmer'], r['offset']) for r in entry['results']}
                self.assertEqual(expected, got)
                placed = set()
                for r in entry['results']:
                    self.assertEqual(r['kmer'][r['offset']:r['offset'] + 3 * len(residues)],
                                     r['instance'])
                    columns = {c for c, recs in self.column_records.items()
                               if recs.has_kmer(r['kmer'])}
                    self.assertEqual(columns, {x['column'] for x in r['labels']}, r['kmer'])
                    for label in r['labels']:
                        for o in label['occurrence_list']:
                            placed.add((label['column'], o['seq_id'],
                                        int(o['nt_coords'].split('-')[0]), o['strand']))
                self.assertEqual(self.peptide_placed(residues, table), placed)
                self.assertCount(entry['counts']['occurrences'], len(placed),
                                 unit='placed_occurrences')

    def test_peptide_paths_against_the_six_frames(self):
        """Peptides longer than k (14-16 residues) with long_search "paths": the paths the
        graph-walk oracle's (anchored on the six-frame translation, extended through the codon
        automaton over the records' k-mers), each label of a path record_verified exactly where
        a record holds the whole path, and the placed occurrences the six-frame matches of the
        whole peptide; without long_search the anchors only."""
        rng = random.Random(77)
        recs = [seq.upper() for r in self.columns.values() for _, seq in r]
        peptides = [('MELPNIMHPVAKLS', 11)]
        while len(peptides) < 4:
            seq = rng.choice(recs)
            i = rng.randrange(len(seq) - 60)
            residues = translate(seq[i:i + 3 * 15], 11)
            if '*' not in residues and '?' not in residues:
                peptides.append((residues, 11))
        for residues, table in peptides:
            with self.subTest(peptide=residues):
                request = {'patterns': [{'protein': residues}], 'genetic_code': table,
                           'long_search': 'paths', 'output': {'labels': 'all'}}
                entry = self.pattern(self.server, request)['patterns'][0]
                self.assertPeptideEntry(entry, residues, table)
                expected = self.peptide_paths(residues, table)
                self.assertGreater(len(expected), 0)
                c = entry['counts']['paths']
                self.assertEqual(('exact', len(expected), 'completed'),
                                 (c['relation'], c['value'], c['extension']))
                self.assertCompleteRetrieval(entry)
                got = {(r['strand'], r['sequence']) for r in entry['results']}
                self.assertEqual(expected, got)
                placed = set()
                for r in entry['results']:
                    s = r['sequence']
                    self.assertEqual(s[:self.k], r['anchor_kmer'])
                    carriers, verified = self.path_columns(s)
                    self.assertEqual(carriers, {x['column'] for x in r['labels']})
                    for label in r['labels']:
                        occ = verified.get(label['column'], set())
                        self.assertEqual('record_verified' if occ else 'label_intersection',
                                         label['support'])
                        for o in label['occurrence_list']:
                            placed.add((label['column'], o['seq_id'],
                                        int(o['nt_coords'].split('-')[0]), o['strand']))
                self.assertEqual(self.peptide_placed(residues, table), placed)
                # the option is opt-in: the anchors only without it
                plain = self.pattern(self.server, {'patterns': [{'protein': residues}],
                                                   'genetic_code': table})['patterns'][0]
                self.assertEqual({'reason': 'paths_later_increment'}, plain['withheld'])
                self.assertEqual('unknown', plain['counts']['paths']['relation'])

    def test_peptide_slots_and_genetic_codes(self):
        """The slot error of a peptide (bad_alphabet; the stop '*' is a residue), an unknown
        genetic code refused, and the genetic code read: a peptide in the vertebrate
        mitochondrial code (2) answered as the six-frame translation in that code finds it."""
        out = self.pattern(self.server, {'patterns': [{'protein': 'MELPUIMHPV'},
                                                      {'protein': 'MELPNIMHPV*'},
                                                      {'protein': 'melpnimhpv'}],
                                         'mode': 'count'})
        self.assertEqual('bad_alphabet', out['patterns'][0]['error']['code'])
        self.assertEqual({'id', 'kind', 'error'}, set(out['patterns'][0]))
        # 11 residues, 33 bases: a pattern longer than k, its anchors counted
        self.assertPeptideEntry(out['patterns'][1], 'MELPNIMHPV*', 1)
        self.assertEqual('exact', out['patterns'][1]['counts']['anchors']['relation'])
        self.assertPeptideEntry(out['patterns'][2], 'MELPNIMHPV', 1)
        ret = self.server.post('pattern', {'patterns': [{'protein': 'MELPNIMHPV'}],
                                           'genetic_code': 7})
        self.assertEqual((400, 'genetic_code_unknown'), (ret.status_code, ret.json()['code']))
        ret = self.server.post('pattern', {'patterns': [{'protein': 'MELPNIMHPV'}],
                                           'genetic_code': '11'})
        self.assertEqual((400, 'invalid_request'), (ret.status_code, ret.json()['code']))
        for residues in ('MELPNIMHPV', 'AKLSTALAAA'):
            entry = self.pattern(self.server, {'patterns': [{'protein': residues}],
                                               'genetic_code': 2, 'mode': 'count'})['patterns'][0]
            self.assertPeptideEntry(entry, residues, 2)
            self.assertCount({x: entry['counts']['contexts'][x]
                              for x in ('value', 'relation', 'unit')},
                             len(self.peptide_contexts(residues, 2)))

    def stop_peptides(self):
        """Peptides of 10 residues (30 bases, one k-mer) read from the records' translations
        with a stop codon among their middle residues, each with the table it was read in (1,
        2, 11: their stops differ, table 2 stopping at AGA and AGG and reading TGA as W)."""
        rng = random.Random(1919)
        recs = [seq.upper() for r in self.columns.values() for _, seq in r]
        out = []
        for table in (1, 2, 11, 11):
            while True:
                seq = rng.choice(recs)
                i = rng.randrange(len(seq) - 40)
                residues = translate(seq[i:i + 30], table)
                if residues.count('*') == 1 and 2 <= residues.index('*') <= 7 \
                        and '?' not in residues and (residues, table) not in out:
                    out.append((residues, table))
                    break
        return out

    def test_peptide_stops_against_the_six_frames(self):
        """The stop '*' is a stop codon of the request's table at that position, against the
        six-frame oracle in tables 1, 2 and 11 (contexts, labels and placed occurrences); X
        never a stop; in table 27 or 31, whose stop codons code a residue unless in context,
        '*' matches nothing (exact 0, the note no_stop_codon) and those codons match as their
        residue (Q at TAA and TAG in 27, E in 31)."""
        for residues, table in self.stop_peptides():
            with self.subTest(peptide=residues, table=table):
                request = {'patterns': [{'protein': residues}], 'genetic_code': table,
                           'output': {'labels': 'all'}}
                entry = self.pattern(self.server, request)['patterns'][0]
                self.assertPeptideEntry(entry, residues, table)
                expected = self.peptide_contexts(residues, table)
                self.assertGreater(len(expected), 0)
                self.assertCount({x: entry['counts']['contexts'][x]
                                  for x in ('value', 'relation', 'unit')}, len(expected))
                self.assertCompleteRetrieval(entry)
                self.assertEqual(expected, {(r['strand'], r['kmer'], r['offset'])
                                            for r in entry['results']})
                placed = {(label['column'], o['seq_id'], int(o['nt_coords'].split('-')[0]),
                           o['strand'])
                          for r in entry['results'] for label in r['labels']
                          for o in label['occurrence_list']}
                self.assertEqual(self.peptide_placed(residues, table), placed)
                self.assertNotIn('no_stop_codon', entry['notes'])
                # X in the stop's place: never a stop, so not these hits
                x = residues.replace('*', 'X')
                entry = self.pattern(self.server, dict(request, patterns=[{'protein': x}],
                                                       mode='count'))['patterns'][0]
                self.assertEqual(len(self.peptide_contexts(x, table)),
                                 entry['counts']['contexts']['value'])
                self.assertTrue(self.peptide_contexts(x, table).isdisjoint(expected))
                # tables without a stop: '*' matches nothing, and the entry says why; the
                # codon read as a stop in |table| matches as its residue there when it codes one
                i = residues.index('*')
                for other in (27, 31):
                    entry = self.pattern(self.server, dict(request, genetic_code=other)
                                         )['patterns'][0]
                    self.assertPeptideEntry(entry, residues, other)
                    self.assertEqual({'value': 0, 'relation': 'exact', 'unit': 'graph_contexts'},
                                     {x: entry['counts']['contexts'][x]
                                      for x in ('value', 'relation', 'unit')})
                    self.assertEqual([], entry['results'])
                    self.assertTrue(entry['retrieval_complete'])
                    self.assertIn('no_stop_codon', entry['notes'])
                    for aa in 'QEW':
                        r = residues[:i] + aa + residues[i + 1:]
                        if '*' in r:
                            continue
                        entry = self.pattern(self.server, {'patterns': [{'protein': r}],
                                                           'genetic_code': other,
                                                           'mode': 'count'})['patterns'][0]
                        self.assertEqual(len(self.peptide_contexts(r, other)),
                                         entry['counts']['contexts']['value'], (r, other))
                        self.assertNotIn('no_stop_codon', entry['notes'])

    # ------------------------------------------------------------ labels

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
        # the rows read, each distinct k-mer once (BASIC: a k-mer is its row); the paths'
        # counters are not a context's
        self.assertEqual(len({kmer for _, kmer, _ in expected}),
                         entry['work']['annotation_rows_distinct'])
        self.assertNotIn('verification_steps', entry['work'])
        self.assertNotIn('verification_ms', entry['timing'])
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

    def test_labels_all_without_record_mapping(self):
        """--no-coord-mapping: coordinates without the .seqs. Each label lists the column
        coordinate of the context's k-mer (kmer_coord) and the offset, placed nowhere; no
        occurrence count is claimed (note record_bounds_unknown)."""
        server = self.no_map_server()
        caps = server.get('capabilities').json()['pattern']
        self.assertEqual('global', caps['placement'])
        entry = self.pattern(server, {'patterns': [{'dna': self.NDM_F}],
                                      'output': {'labels': 'all'}})['patterns'][0]
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

    # ------------------------------------------------------------ predicates

    PREDICATE_PATTERNS = [('dna', 'GGTTTGGCGATCTGGTTTTC'), ('dna', 'CGGAATGGCTCATCACGATC'),
                          ('dna', 'GCGGCGGCGGCG'), ('iupac', 'GTGYCAGCMGCCGCGGTAA'),
                          ('iupac', 'GGACTACNVGGGTWTCTAAT')]

    def kmer_columns(self, kmers):
        """{k-mer: the columns whose records hold it as deposited}, by a scan of the FASTA files
        (each record's every k-window; the mini's records hold no N)."""
        if not hasattr(self.__class__, '_kmer_columns'):
            self.__class__._kmer_columns = {}
        cache = self.__class__._kmer_columns
        todo = set(kmers) - set(cache)
        if todo:
            found = {x: set() for x in todo}
            k = self.k
            for column, recs in self.columns.items():
                for _, seq in recs:
                    seq = seq.upper()
                    for i in range(len(seq) - k + 1):
                        w = seq[i:i + k]
                        if w in found:
                            found[w].add(column)
            cache.update(found)
        return {x: cache[x] for x in kmers}

    @staticmethod
    def random_predicate(rng, names, depth=0):
        def names_list():
            return rng.sample(names, rng.randint(1, 3))
        op = rng.randrange(4 if depth >= 2 else 7)
        if op == 0:
            return {'any': names_list()}
        if op == 1:
            return {'all': names_list()}
        if op == 2:
            return {'none': names_list()}
        if op == 3:
            labels = names_list()
            return {'at_least': {'n': rng.randint(1, len(labels)), 'labels': labels}}
        if op == 6:
            return {'not': TestPatternMini.random_predicate(rng, names, depth + 1)}
        return {('and', 'or')[op - 4]: [TestPatternMini.random_predicate(rng, names, depth + 1)
                                        for _ in range(rng.randint(1, 3))]}

    @staticmethod
    def holds(p, labels):
        (op, v), = p.items()
        if op in ('any', 'all', 'none'):
            hits = sum(n in labels for n in v)
            return hits > 0 if op == 'any' else hits == len(v) if op == 'all' else hits == 0
        if op == 'at_least':
            return sum(n in labels for n in v['labels']) >= v['n']
        if op == 'and':
            return all(TestPatternMini.holds(q, labels) for q in v)
        if op == 'or':
            return any(TestPatternMini.holds(q, labels) for q in v)
        return not TestPatternMini.holds(v, labels)

    def test_predicate_against_the_fasta(self):
        """SPEC §19 (TESTS §8): for the blaNDM-1 primers, the GCG repeat and two IUPAC primers,
        20 random predicates over the nine columns, both strand settings: the selected contexts
        are those of a scan of the mini's FASTA whose k-mer's columns (with "either" also its
        reverse complement's) satisfy the predicate, all of them returned, counted exactly;
        with "predicate_only" every result lists the predicate's columns of its own k-mer, its
        selection_labels the set it was evaluated on, and per label in selection_strands the
        orientation whose records carry it: "context" (its k-mer only), "reverse_complement"
        (with "either", its reverse complement only) or "both"."""
        rng = random.Random(5)
        names = sorted(self.columns)
        contexts = {p: self.contexts(p) for _, p in self.PREDICATE_PATTERNS}
        kmers = {kmer for c in contexts.values() for _, kmer, _ in c}
        columns = self.kmer_columns(kmers | {revcomp(x) for x in kmers})
        for t in range(20):
            pred = self.random_predicate(rng, names)
            for strands in ('either', 'context'):
                labels = 'predicate_only' if t % 4 == 0 else 'none'
                request = {'patterns': [{kind: p} for kind, p in self.PREDICATE_PATTERNS],
                           'predicate': pred, 'predicate_strands': strands,
                           'output': {'labels': labels}}
                out = self.pattern(self.server, request)
                with self.subTest(predicate=pred, strands=strands):
                    # every name a column: nothing folded away but single operands
                    self.assertEqual([], out['predicate']['unknown_labels'])
                    self.assertEqual(len(self.holds_names(pred)), out['predicate']['known'])
                    self.assertEqual(strands, out['predicate']['strands'])
                    for (_, p), e in zip(self.PREDICATE_PATTERNS, out['patterns']):
                        def evaluated(kmer):
                            s = set(columns[kmer])
                            if strands == 'either':
                                s |= columns[revcomp(kmer)]
                            return s
                        want = {c for c in contexts[p] if self.holds(pred, evaluated(c[1]))}
                        got = {(r['strand'], r['kmer'], r['offset']) for r in e['results']}
                        self.assertEqual(want, got, p)
                        self.assertEqual('completed', e['selection']['pass'], p)
                        self.assertCount(e['counts']['tested'], len(contexts[p]))
                        self.assertCount(e['counts']['selected'], len(want))
                        self.assertTrue(e['retrieval_complete'], p)
                        self.assertEqual('predicate', e['absence_filter'])
                        if labels == 'none':
                            continue
                        keep = {n for n in self.holds_names(pred)}
                        for r in e['results']:
                            own = {label['column'] for label in r['labels']}
                            self.assertEqual(columns[r['kmer']] & keep, own, r['kmer'])
                            self.assertEqual(evaluated(r['kmer']) & keep,
                                             set(r['selection_labels']), r['kmer'])
                            self.assertEqual(len(r['selection_labels']),
                                             len(r['selection_strands']), r['kmer'])
                            for n, strand in zip(r['selection_labels'], r['selection_strands']):
                                x = n in columns[r['kmer']]
                                y = strands == 'either' and n in columns[revcomp(r['kmer'])]
                                want_strand = 'both' if x and y else \
                                    'reverse_complement' if y else 'context'
                                self.assertEqual(want_strand, strand, (r['kmer'], n))

    @staticmethod
    def holds_names(p):
        (op, v), = p.items()
        if op in ('any', 'all', 'none'):
            return set(v)
        if op == 'at_least':
            return set(v['labels'])
        if op == 'not':
            return TestPatternMini.holds_names(v)
        return set().union(*(TestPatternMini.holds_names(q) for q in v))

    def test_predicate_strands_on_the_ndm_primer(self):
        """Strands on the mini (SPEC §19.5): none(546) on the blaNDM-1 forward primer selects
        no context with "either" and its 12 - contexts with "context"; a typo is reported
        unknown and folded to a constant (no row read); the compute admission and the work
        budget stop as stated; the CLI answers as the server."""
        p = [{'dna': self.NDM_F}]
        out = self.pattern(self.server, {'patterns': p, 'mode': 'count',
                                         'predicate': {'none': ['546']}})
        self.assertCount(out['patterns'][0]['counts']['selected'], 0)
        self.assertTrue(out['predicate']['vacuous'])
        out = self.pattern(self.server, {'patterns': p, 'predicate': {'none': ['546']},
                                         'predicate_strands': 'context'})
        e = out['patterns'][0]
        self.assertCount(e['counts']['selected'], 12)
        self.assertEqual({'-'}, {r['strand'] for r in e['results']})
        out = self.pattern(self.server, {'patterns': p, 'predicate': {'any': ['5622']}})
        e = out['patterns'][0]
        self.assertEqual((['5622'], False), (out['predicate']['unknown_labels'],
                                             out['predicate']['normal_form']))
        self.assertEqual(('constant', 0, ['predicate_constant']),
                         (e['selection']['pass'], e['work']['predicate_rows'], e['notes']))
        self.assertTrue(e['retrieval_complete'])
        out = self.pattern(self.server, {'patterns': [{'dna': 'GCGGCGGCGGCG'}],
                                         'predicate': {'any': ['562']},
                                         'max_predicate_contexts': 100})
        e = out['patterns'][0]
        self.assertEqual(({'reason': 'predicate_above_threshold'}, 'not_admitted', 0),
                         (e['withheld'], e['selection']['pass'], e['work']['predicate_rows']))
        out = self.pattern(self.server, {'patterns': [{'dna': 'GCGGCGGCGGCG'}], 'mode': 'partial',
                                         'predicate': {'any': ['562']},
                                         'predicate_strands': 'context', 'max_predicate_work': 1})
        e = out['patterns'][0]
        self.assertEqual(({'phase': 'selection', 'reason': 'max_predicate_work'},
                          {'reason': 'max_predicate_work'}, 1),
                         (e['stop'], e['cut'], e['counts']['tested']['value']))
        self.assertEqual('bounds', e['counts']['selected']['relation'])
        # the CLI answers as the server
        request = {'patterns': [{'dna': self.NDM_F}, {'iupac': self.iupac16}, {'dna': self.p40}],
                   'mode': 'partial', 'max_contexts': 9, 'predicate': {'any': ['562', '573']},
                   'output': {'labels': 'predicate_only'}}
        server_out = self.pattern(self.server, request)
        path = os.path.join(self.tempdir.name, 'request_predicate.json')
        with open(path, 'w') as f:
            json.dump(request, f)
        res = subprocess.run(shlex.split(METAGRAPH) + ['pattern', '--json', '-i', self.graph,
                                                       '-a', self.anno, path],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(0, res.returncode, res.stderr.decode())
        self.assertEqual(untimed(server_out), untimed(json.loads(res.stdout)))

    def test_cli_answers_as_the_server(self):
        """`metagraph pattern` answers each request file as the server: one process over a panel
        (counts and a partial retrieval, paths with labels and record verification, peptides,
        labels "all" with caps), a refusal's body with exit status 1, and every file answered,
        the later ones after a refused one (a body nested too deep must not abort the run and
        leave the later files unanswered)."""
        panel = [
            {'patterns': [{'id': 'a', 'dna': self.p16}, {'iupac': self.iupac16},
                          {'dna': self.p11}, {'dna': self.p40}],
             'mode': 'partial', 'max_contexts': 7},
            self.paths_request(self.long_patterns()[:3], mode='partial', max_paths=1,
                               output={'labels': 'all'}, require_support='record_verified'),
            {'patterns': [{'protein': 'MELPNIMHPV'}, {'protein': 'MELPNIMHPVAKLS'}],
             'genetic_code': 11, 'long_search': 'paths', 'output': {'labels': 'all'}},
            {'patterns': [{'dna': self.NDM_F}, {'iupac': self.iupac16}, {'dna': self.p40}],
             'mode': 'partial', 'max_contexts': 9, 'max_labels': 4,
             'output': {'labels': 'all'}},
        ]
        files = []
        for i, request in enumerate(panel):
            files.append(os.path.join(self.tempdir.name, f'cli_request_{i}.json'))
            with open(files[-1], 'w') as f:
                json.dump(request, f)
        refused = {
            'bogus.json': json.dumps({'patterns': [{'dna': self.p16}], 'bogus': True}),
            'deep.json': '[' * 1200 + ']' * 1200,
            'dup.json': '{"patterns": [{"dna": "%s"}], "mode": "count", "mode": "partial"}'
                        % self.p16,
        }
        for name, text in refused.items():
            files.append(os.path.join(self.tempdir.name, name))
            with open(files[-1], 'w') as f:
                f.write(text)
        files.append(os.path.join(self.tempdir.name, 'cli_ok.json'))
        with open(files[-1], 'w') as f:
            json.dump({'patterns': [{'dna': self.p16}], 'mode': 'count'}, f)
        res = subprocess.run(shlex.split(METAGRAPH) + ['pattern', '--json', '-i', self.graph,
                                                       '-a', self.anno] + files,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(1, res.returncode, res.stderr.decode()[-2000:])
        lines = res.stdout.decode().strip().split('\n')
        self.assertEqual(len(files), len(lines), lines)
        for request, line in zip(panel, lines):
            self.assertEqual(untimed(self.pattern(self.server, request)),
                             untimed(json.loads(line)), request)
        for line in lines[len(panel):len(panel) + len(refused)]:
            self.assertEqual('invalid_request', json.loads(line)['code'])
        self.assertContextCounts(json.loads(lines[-1])['patterns'][0], self.contexts(self.p16),
                                 self.k, 16, 'any_offset')

    def test_cli_unmasked(self):
        """The mini index as built (no .edgemask) is served by the CLI too: its answer states
        the counting and the dummy fraction, as the server's on the same graph (the fraction
        sampled in another process: the same draws); the start-up note names the remedies for
        exact counts; with --pattern-build-mask the CLI answers as the server on the
        transformed copy."""
        request = {'patterns': [{'dna': self.p16}, {'iupac': self.iupac16},
                                {'dna': self.start16}], 'mode': 'count'}
        path = os.path.join(self.tempdir.name, 'request_mask.json')
        with open(path, 'w') as f:
            json.dump(request, f)
        cli = shlex.split(METAGRAPH) + ['pattern']
        args = ['--json', '-i', self.unmasked_graph, '-a', self.unmasked_anno, path]
        res = subprocess.run(cli + args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(0, res.returncode, res.stderr.decode())
        out = json.loads(res.stdout)
        self.assertEqual('upper_bound', out['index']['counting'])
        caps = self.unmasked_server.get('capabilities').json()['pattern']
        self.assertEqual(caps['dummy_fraction'], out['index']['dummy_fraction'])
        self.assertEqual(untimed(self.pattern(self.unmasked_server, request)), untimed(out))
        # the candidate check off: as the server started so
        res0 = subprocess.run(cli + ['--pattern-max-checked-entries', '0'] + args,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(0, res0.returncode, res0.stderr.decode())
        self.assertEqual(untimed(self.pattern(self.unchecked_server, request)),
                         untimed(json.loads(res0.stdout)))
        # the start-up note names the remedies for exact counts
        self.assertIn('the pattern search counts upper bounds with estimates (counting: '
                      'upper_bound). For exact counts: `metagraph transform --mask-dummy '
                      f'{self.unmasked_graph}` once', res.stderr.decode())

        res = subprocess.run(cli + ['--pattern-build-mask'] + args,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(0, res.returncode, res.stderr.decode())
        self.assertEqual(untimed(self.pattern(self.server, request)),
                         untimed(json.loads(res.stdout)))

    def test_cli_unreadable_mask(self):
        """A .edgemask that is there but cannot be opened (permissions): the graph loads without
        it, and the messages say so rather than "no mask", whose remedy `transform
        --mask-dummy` would refuse."""
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
            # served without it: upper bounds, as on the unmasked graph
            self.assertEqual(0, res.returncode, err)
            out = json.loads(res.stdout)
            self.assertEqual('upper_bound', out['index']['counting'])
            self.assertEqual(untimed(self.pattern(self.unmasked_server, request)), untimed(out))
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
        discovery, and a stop inside a mask scan (its bounds count the scanned edges)."""
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
        # and the mini index itself is as it was
        self.assertEqual(self.mini_files, sorted(os.listdir(MINI_DIR)))

    # ------------------------------------------------------------ without the mask (#16)

    def exact_dummy_fraction(self):
        """The mini graph's exact fraction of real k-mers among its edges with W != $, from
        `metagraph stats --count-dummy` (the main dummy edge, W = $, among the source
        dummies)."""
        res = TestingBase._run_command(f'{METAGRAPH} stats --count-dummy {self.unmasked_graph}',
                                       'Count the dummies of the mini graph')
        text = res.stdout.decode()
        edges = int(re.search(r'edges \( k \): (\d+)', text).group(1))
        source = int(re.search(r'dummy source edges: (\d+)', text).group(1))
        sink = int(re.search(r'dummy sink edges: (\d+)', text).group(1))
        real = int(re.search(r'real edges: (\d+)', text).group(1))
        self.assertEqual(edges, source + sink + real)
        # the source dummies: the oracle's, and the main dummy edge
        self.assertEqual(source - 1, len(self.records.source_dummies(self.k)))
        return real / (edges - sink - 1)

    def test_unmasked_capabilities(self):
        """Without the mask: served, mask absent, counting upper_bound, and the dummy fraction
        sampled from 10,000 entries, its 95% interval holding the exact fraction (stats
        --count-dummy); the rest of the block as on the masked graph, which counts exactly and
        has no fraction."""
        f = self.exact_dummy_fraction()
        for route in ('capabilities', 'traverse/capabilities'):
            p = self.unmasked_server.get(route).json()['pattern']
            self.assertEqual(('absent', True, None, 'upper_bound'),
                             (p['mask'], p['available'], p['unavailable_reason'], p['counting']))
            frac = p['dummy_fraction']
            self.assertEqual({'value', 'interval', 'samples', 'source'}, set(frac))
            self.assertEqual((10000, 'sampled'), (frac['samples'], frac['source']))
            lo, hi = frac['interval']
            self.assertLessEqual(lo, f)
            self.assertLessEqual(f, hi)
            self.assertLessEqual(lo, frac['value'])
            self.assertLessEqual(frac['value'], hi)
            q = self.server.get(route).json()['pattern']
            self.assertEqual(('file', 'exact', None), (q['mask'], q['counting'],
                                                       q['dummy_fraction']))
            self.assertEqual(dict(q, mask='absent', counting='upper_bound', dummy_fraction=frac),
                             p)
        # (the probe's size under the MCP tool's ceiling: the fixture validator's
        # test_capabilities_documents_keep_a_kibibyte, on the stored documents --check compares)
        # the log names it once, at start-up
        with open(os.path.join(self.tempdir.name, 'server_unmasked.log')) as log:
            text = log.read()
        self.assertRegex(text, r'Dummy fraction sampled for the pattern search in [\d.]+ s')
        self.assertIn('The graph has no dummy-edge mask (.edgemask): the pattern search counts '
                      'upper bounds with estimates', text)

    def assertUnmaskedCount(self, count, exact, upper, fraction, where=''):
        """A count of the unmasked graph against the masked graph's exact count and the oracle's
        upper bound (the exact contexts and those of the source dummies): exact only as the
        exact count, exact 0 where nothing is a candidate, else bounds holding it, the upper
        bound the oracle's, the estimate upper x f rounded into the bounds."""
        msg = f'{where}: {count}'
        # the count's own fields (a total also carries suffix, by_offset, by_strand)
        own = {k for k in count if k in ('value', 'relation', 'unit', 'lower', 'upper',
                                         'estimate')}
        if upper == 0:
            self.assertEqual({'value': 0, 'relation': 'exact', 'unit': count['unit']},
                             {k: count[k] for k in own}, msg)
            return
        if count['relation'] == 'exact':
            self.assertEqual(exact, count['value'], msg)
            self.assertEqual({'value', 'relation', 'unit'}, own, msg)
            return
        self.assertEqual('bounds', count['relation'], msg)
        self.assertEqual({'value', 'relation', 'unit', 'lower', 'upper', 'estimate'}, own, msg)
        self.assertEqual(count['lower'], count['value'], msg)
        self.assertLessEqual(count['lower'], exact, msg)
        self.assertEqual(upper, count['upper'], msg)
        e = min(upper, max(count['lower'], round(upper * fraction)))
        self.assertEqual(e, count['estimate'], msg)

    def test_unmasked_against_the_masked(self):
        """A graph without its dummy-edge mask, on the mini index, whose rules the unit tests
        check on small graphs in every graph mode (PatternRoute.Unmasked*, PatternUnmasked.*):
        served without its mask, a count is the bounds [lower, U], U the oracle's (the contexts
        of the k-mers and of the source dummies) and the estimate U x f rounded into them, an
        absent pattern exact 0; the lists, labels and placed occurrences included, are the
        masked graph's. At the server's default --pattern-max-checked-entries a pattern with
        few unchecked candidates has each tested: its counts the masked graph's, k - 1 steps
        per candidate."""
        fraction = self.unchecked_server.get('capabilities').json()['pattern'][
            'dummy_fraction']['value']
        patterns = [self.p16, self.iupac16, self.absent16, self.start16]
        request = {'patterns': [{'iupac' if set(p) - set('ACGT') else 'dna': p}
                                for p in patterns], 'output': {'labels': 'all'}}
        out = self.pattern(self.unchecked_server, dict(request, mode='count'))
        masked = self.pattern(self.server, dict(request, mode='count'))
        self.assertEqual('upper_bound', out['index']['counting'])
        self.assertNotIn('counting', masked['index'])
        for p, e, m in zip(patterns, out['patterns'], masked['patterns']):
            exact = self.contexts(p)
            dummies = self.records.dummy_contexts(p, self.k)
            self.assertEqual(len(exact), m['counts']['contexts']['value'])
            self.assertUnmaskedCount(e['counts']['contexts'], len(exact),
                                     len(exact) + len(dummies), fraction, p)
        # start16 begins an island: its source dummies are candidates the count cannot rule out
        start = out['patterns'][3]
        self.assertGreater(len(self.records.dummy_contexts(self.start16, self.k)), 0)
        self.assertEqual('bounds', start['counts']['contexts']['relation'])
        self.assertIn(NOTE_ESTIMATE, start['notes'])
        out = self.pattern(self.unchecked_server, request)
        masked = self.pattern(self.server, request)
        for e, m in zip(out['patterns'], masked['patterns']):
            for field in ('results', 'returned', 'withheld', 'cut', 'by_label',
                          'retrieval_complete'):
                self.assertEqual(m.get(field), e.get(field), (m['pattern'], field))
        # the candidate check at the default limit: start16's few unchecked candidates are tested
        count = {'patterns': [{'dna': self.start16}], 'mode': 'count'}
        c0 = start['counts']['contexts']
        width = c0['upper'] - c0['lower']
        limit = self.unmasked_server.get('capabilities').json()['pattern']['caps'][
            'max_checked_entries']
        self.assertTrue(0 < width <= limit, (width, limit))
        e = self.pattern(self.unmasked_server, count)['patterns'][0]
        e0 = self.pattern(self.unchecked_server, count)['patterns'][0]
        self.assertEqual(self.pattern(self.server, count)['patterns'][0]['counts'], e['counts'])
        self.assertNotIn(NOTE_ESTIMATE, e['notes'])
        self.assertEqual(e0['work']['steps'] + width * (self.k - 1), e['work']['steps'])

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
        cls.unmasked_servers = {}
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
            # the same graph without its mask: the .dbg alone in a directory of its own (the
            # loader reads a mask beside the path it is given)
            u = os.path.join(d, f'unmasked_{mode}')
            os.makedirs(u)
            os.symlink(graph, os.path.join(u, os.path.basename(graph)))
            cls.unmasked_servers[mode] = Server(
                METAGRAPH, ['-i', os.path.join(u, os.path.basename(graph)), '-a', anno,
                            '--pattern-min-information-bits', str(cls.FLOOR)],
                os.path.join(d, f'server_unmasked_{mode}.log'))
        cls.graph_basic = os.path.join(d, 'graph_basic.dbg')
        cls.anno_basic = os.path.join(d, 'anno_basic.column.annodbg')

    @classmethod
    def tearDownClass(cls):
        for server in list(cls.servers.values()) + list(cls.unmasked_servers.values()):
            server.stop()
        super().tearDownClass()

    def test_refusals(self):
        """A refusal over HTTP: 400 {error, code}, the message naming the field, for one case of
        each code a request can draw and for bodies that are no single RFC 8259 JSON text
        (duplicate members, a nesting past jsoncpp's limit). The full table, every field and
        the order of the checks, is PatternRoute.Refusals and RefusalOrder."""
        server = self.servers['basic']
        p = [{'dna': 'TTAGGACGACTTTG'}]
        cases = [
            ({'patterns': p, 'bogus': 1}, 'invalid_request', "unknown field 'bogus'"),
            ({'patterns': p, 'max_steps': 0}, 'invalid_request', 'max_steps'),
            ({'patterns': p, 'predicate': {'any': [562]}}, 'invalid_request',
             'write a taxid as "562"'),
            ({'patterns': p, 'predicate': {'any': [str(i) for i in range(10001)]}},
             'predicate_too_large', '10001 names'),
            ({'patterns': p, 'graphs': ['x']}, 'invalid_request', 'single graph'),
            ({'patterns': p, 'budget_split': 1}, 'later_increment', 'budget_split'),
            ({'patterns': [{'protein': 'MKV'}], 'genetic_code': 7}, 'genetic_code_unknown',
             'genetic_code'),
            ({'patterns': p, 'in_ram': 'yes'}, 'invalid_request', 'in_ram'),
        ]
        for payload, code, words in cases:
            ret = server.post('pattern', payload)
            self.assertEqual(400, ret.status_code, (payload, ret.text))
            body = ret.json()
            self.assertEqual({'error', 'code'}, set(body), payload)
            self.assertEqual(code, body['code'], (payload, body))
            self.assertIn(words, body['error'], payload)
        one = '{"patterns": [{"dna": "TTAGGACGACTTTG"}]'
        for raw in ('{"patterns": [', one + '} GARBAGE',
                    one + ', "mode": "count", "mode": "partial"}',
                    one + ', "x": ' + '[' * 1500 + ']' * 1500 + '}'):
            ret = server.post('pattern', raw, raw=True)
            self.assertEqual(400, ret.status_code, raw[:100])
            self.assertEqual({'error', 'code'}, set(ret.json()), raw[:100])
            self.assertEqual('invalid_request', ret.json()['code'], raw[:100])
        # false projections and occurrences are accepted and change nothing; so are
        # output.paths true and long_search for a pattern of at most k bases
        out = self.pattern(server, {'patterns': p, 'mode': 'count'})
        for extra in ({'output': {'labels': 'none', 'occurrences': False, 'paths': False}},
                      {'output': {'labels': 'none', 'paths': True}},
                      {'long_search': 'paths'}, {'long_search': 'anchors'}):
            again = self.pattern(server, dict({'patterns': p, 'mode': 'count'}, **extra))
            self.assertEqual(untimed(out['patterns']), untimed(again['patterns']), extra)

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

    def test_paths_in_every_graph_mode(self):
        """long_search "paths" on BASIC, CANONICAL and PRIMARY graphs: exactly the graph-walk
        oracle's paths (over both orientations of every k-mer on the last two), and with labels
        (a column per record, unbudgeted) each path's columns the records holding every k-mer
        of it (either orientation there), support label_intersection: no coordinates, nothing
        verified."""
        rng = random.Random(23)
        islands = [island for _, _, island in self.records.islands if len(island) > 100]
        patterns = []
        for _ in range(4):
            island = rng.choice(islands)
            L = rng.randint(20, 30)
            i = rng.randrange(len(island) - L)
            patterns.append(island[i:i + L])
        w = list(patterns[0])
        w[17] = 'N'
        patterns.append(''.join(w))
        # rec_motifs: a repeat, its paths along the copies
        patterns.append('TTTACGACGACTTTGAATTCGAA')
        by_name = {}
        for _, name, island in self.records.islands:
            by_name.setdefault(name, []).append(island)
        for mode, server in self.servers.items():
            stated = mode == 'basic'
            key = 'strand' if stated else 'orientation'
            request = {'patterns': [{'iupac': p} for p in patterns], 'long_search': 'paths',
                       'output': {'labels': 'all'}, 'allow_unbudgeted_annotation': True}
            out = self.pattern(server, request)
            for entry, p in zip(out['patterns'], patterns):
                with self.subTest(mode=mode, pattern=p):
                    expected = self.records.paths(p, self.K, canonical=not stated)
                    self.assertGreater(len(expected), 0)
                    c = entry['counts']['paths']
                    self.assertEqual(('exact', len(expected), 'completed'),
                                     (c['relation'], c['value'], c['extension']))
                    self.assertCompleteRetrieval(entry)
                    self.assertEqual(expected, {(r[key], r['sequence'])
                                                for r in entry['results']})
                    self.assertIn('label_intersection_only', entry['notes'])
                    for r in entry['results']:
                        s = r['sequence']
                        self.assertEqual(s[:self.K], r['anchor_kmer'])
                        self.assertNotIn('kmer', r)
                        windows = [s[i:i + self.K] for i in range(len(s) - self.K + 1)]
                        carriers = {name for name, isl in by_name.items()
                                    if all(any(x in y or (not stated and revcomp(x) in y)
                                               for y in isl) for x in windows)}
                        self.assertEqual(carriers, {x['column'] for x in r['labels']}, s)
                        for label in r['labels']:
                            self.assertEqual('label_intersection', label['support'])

    def test_unmasked_in_every_graph_mode(self):
        """Graphs without their mask on BASIC, CANONICAL and PRIMARY: the counts of short
        motifs (many of them in source dummies at k = 15) against the record oracle and, on
        BASIC and CANONICAL, the source-dummy oracle (the upper bound exactly the contexts of
        the k-mers and of the dummies); on PRIMARY, whose wrapper searches the stored graph
        twice, an upper bound at least the exact count; every list the masked graph's, every
        complete release exact; the paths of long patterns the masked graph's."""
        for mode in ('basic', 'canonical', 'primary'):
            stated = mode == 'basic'
            caps = self.unmasked_servers[mode].get('capabilities').json()['pattern']
            self.assertEqual(('absent', 'upper_bound', True),
                             (caps['mask'], caps['counting'], caps['available']))
            fraction = caps['dummy_fraction']['value']
            for scope in ('any_offset', 'suffix'):
                if mode == 'primary' and scope == 'suffix':
                    continue
                request = {'patterns': [{'iupac': p} for p in self.PATTERNS], 'scope': scope,
                           'mode': 'count'}
                out = self.pattern(self.unmasked_servers[mode], request)
                for entry, p in zip(out['patterns'], self.PATTERNS):
                    with self.subTest(mode=mode, scope=scope, pattern=p):
                        exact = self.records.contexts(p, self.K, scope=scope,
                                                      canonical=not stated)
                        c = entry['counts']['contexts']
                        if c['relation'] == 'exact':
                            self.assertEqual(len(exact), c['value'])
                            continue
                        self.assertEqual('bounds', c['relation'])
                        self.assertLessEqual(c['lower'], len(exact))
                        self.assertIn(NOTE_ESTIMATE, entry['notes'])
                        self.assertEqual(min(c['upper'], max(c['lower'],
                                                             round(c['upper'] * fraction))),
                                         c['estimate'])
                        if mode == 'primary':
                            self.assertGreaterEqual(c['upper'], len(exact))
                            continue
                        dummies = self.records.dummy_contexts(p, self.K, scope=scope,
                                                              canonical=not stated)
                        self.assertEqual(len(exact) + len(dummies), c['upper'])
                # the lists, and their counts once complete
                for m in ('all_or_count', 'partial'):
                    request = {'patterns': [{'iupac': p} for p in self.PATTERNS],
                               'scope': scope, 'mode': m, 'max_contexts': 30}
                    out = self.pattern(self.unmasked_servers[mode], request)
                    masked = self.pattern(self.servers[mode], request)
                    for e, x in zip(out['patterns'], masked['patterns']):
                        with self.subTest(mode=mode, scope=scope, request=m,
                                          pattern=x['pattern']):
                            if m == 'partial' or x['withheld'] is None:
                                self.assertEqual(x['results'], e['results'])
                            if e['retrieval_complete']:
                                self.assertEqual(x['counts'], e['counts'])
                                self.assertEqual(x['results'], e['results'])
            # the paths of long patterns: the masked graph's
            island = max((i for _, _, i in self.records.islands), key=len)
            request = {'patterns': [{'dna': island[100:125]}, {'dna': island[:24]}],
                       'long_search': 'paths', 'output': {'labels': 'none'}}
            out = self.pattern(self.unmasked_servers[mode], request)
            masked = self.pattern(self.servers[mode], request)
            for e, x in zip(out['patterns'], masked['patterns']):
                with self.subTest(mode=mode, paths=x['pattern']):
                    self.assertEqual(x['results'], e['results'])
                    self.assertEqual(x['counts']['paths'], e['counts']['paths'])

    def test_rows(self):
        """The row is the annotation row of the context's k-mer: on a PRIMARY graph served
        wrapped in CanonicalDBG a k-mer and its reverse complement share one row (one annotated
        key). BASIC (row = node - 1) and CANONICAL are PatternRoute.RetrievalRows and
        CanonicalRowsAreTheAnnotationKeys."""
        out = self.pattern(self.servers['primary'],
                           {'patterns': [{'dna': 'GAATTC'}, {'dna': 'ACGAC'}]})
        rows = {}
        for entry in out['patterns']:
            for r in entry['results']:
                self.assertIsNotNone(r['row'])
                rows.setdefault(r['kmer'], set()).add(r['row'])
        shared = 0
        for kmer, row in rows.items():
            self.assertEqual(1, len(row))
            if revcomp(kmer) in rows:
                self.assertEqual(row, rows[revcomp(kmer)], kmer)
                shared += 1
        self.assertGreater(shared, 0)

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
        gone -- its errors no more than its answer (after a complete request and
        shutdown(SHUT_WR), a valid count request gets 0 bytes, and so do `{` and
        {"patterns":[]}, which must not get their 400). The same requests from a client that is
        still there are answered."""
        server = self.servers['basic']
        cases = [
            ('a count', {'patterns': [{'dna': 'ACGAC'}], 'mode': 'count'}, 200, None),
            ('a retrieval', {'patterns': [{'dna': 'ACGAC'}]}, 200, None),
            ('not JSON', '{', 400, 'invalid_request'),
            ('no patterns', {'patterns': []}, 400, 'invalid_request'),
            ('a later increment', {'patterns': [{'dna': 'ACGAC'}], 'budget_split': 1}, 400,
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
        """A multi-graph server answers /pattern as /search: the graphs a request names, each
        pair answered as on a single-graph server and tagged with its pair; the full block per
        pair on /pattern/capabilities?graph= and in the probe; the server-wide block on
        /capabilities with the per-pair summary (test_multigraph.py covers the rest)."""
        d = self.tempdir.name
        csv = os.path.join(d, 'graphs.csv')
        with open(csv, 'w') as f:
            f.write(f'G1,{self.graph_basic},{self.anno_basic}\n')
        server = Server(METAGRAPH, [csv], os.path.join(d, 'server_multi.log'))
        single = Server(METAGRAPH, ['-i', self.graph_basic, '-a', self.anno_basic],
                        os.path.join(d, 'server_multi_single.log'))
        try:
            request = {'patterns': [{'dna': 'ACGAC'}], 'mode': 'count'}
            ret = server.post('pattern', dict(request, graphs=['G1']))
            self.assertEqual(200, ret.status_code, ret.text)
            out = ret.json()
            self.assertEqual(['G1'], out['graphs'])
            self.assertEqual(1, len(out['answers']))
            answer = out['answers'][0]
            self.assertEqual(('G1', self.graph_basic, self.anno_basic, None),
                             (answer['graph'], answer['graph_path'], answer['annotation_path'],
                              answer['index_fp']))
            # the pair's answer is the single-graph server's on the same files
            alone = single.post('pattern', request).json()
            for a in (answer, alone):
                a.pop('timing')
            for key in ('graph', 'graph_path', 'annotation_path', 'index_fp'):
                answer.pop(key)
            self.assertEqual(alone, answer)
            # nor is a refusal written to a client that left
            self.assertEqual(b'', half_closed_request(
                server.port, 'pattern', dict(request, graphs=['nope'])))
            caps = server.get('capabilities').json()
            self.assertIn('pattern', caps['features'])
            self.assertEqual('POST /pattern', caps['routes']['pattern'])
            self.assertEqual('GET /pattern/capabilities?graph={name}[&graph_path={path}]',
                             caps['routes']['pattern_capabilities'])
            self.assertTrue(caps['pattern']['available'])
            self.assertIsNone(caps['pattern']['k'])
            self.assertEqual(['G1'], [p['graph'] for p in caps['graph_summary']['pairs']])
            block = server.get('pattern/capabilities?graph=G1').json()
            self.assertEqual(('G1', self.graph_basic), (block.pop('graph'), block.pop('graph_path')))
            self.assertEqual(single.get('pattern/capabilities').json(), block)
            probe = server.get('traverse/capabilities?graph=G1').json()
            self.assertEqual(dict(block, details='GET /pattern/capabilities'), probe['pattern'])
            ret = server.get('pattern/capabilities')
            self.assertEqual(400, ret.status_code)
            self.assertIn('needs ?graph=<name>', ret.json()['error'])
        finally:
            server.stop()
            single.stop()

    def test_deadline_caps_under_the_content_timeout(self):
        """A route deadline past the server's 900 s content timeout would be accepted and the
        connection closed before any answer or 503: server_query refuses a /pattern cap above
        899,000 ms, and /resolve's opt-in deadline is capped there when --traverse-max-time-ms
        is 0 (uncapped) or above it."""
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

    def test_predicate_caps_at_start_up(self):
        """Predicate caps (SPEC §4.5, §19.3): --pattern-max-predicate-labels is at most
        1,000,000 and --pattern-max-predicate-work at least 1, refused at start-up otherwise, on
        the server and the CLI alike; the CLI applies the caps as the server does (a predicate
        of more names than the cap: predicate_too_large; the caps in its limits); both usages
        name the flags."""
        for command in ('server_query', 'pattern'):
            for flag, value, says in (('--pattern-max-predicate-labels', '1000001',
                                       'must be an integer in [0, 1000000]'),
                                      ('--pattern-max-predicate-work', '0', 'at least 1')):
                cmd = shlex.split(METAGRAPH) + [command, '-i', self.graph_basic, '-a',
                                                self.anno_basic, flag, value]
                if command == 'server_query':
                    cmd += ['--port', str(free_port()), '--address', '127.0.0.1']
                else:
                    cmd += [os.path.join(self.tempdir.name, 'none.json')]
                try:
                    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                         timeout=120)
                except subprocess.TimeoutExpired:
                    self.fail(f'{command} started with {flag} {value}')
                self.assertNotEqual(0, res.returncode, (command, flag))
                self.assertIn(says, res.stderr.decode(), (command, flag))
            res = subprocess.run(shlex.split(METAGRAPH) + [command], stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
            text = (res.stdout + res.stderr).decode()
            for flag in ('--pattern-max-predicate-contexts', '--pattern-max-predicate-work',
                         '--pattern-max-predicate-labels'):
                self.assertIn(flag, text, command)
        # (names that are no column: folded to a constant, nothing read)
        columns = ['a', 'b', 'c']
        path = os.path.join(self.tempdir.name, 'request_predicate_cap.json')
        with open(path, 'w') as f:
            json.dump({'patterns': [{'iupac': self.PATTERNS[0]}], 'mode': 'count',
                       'predicate': {'any': columns}, 'allow_unbudgeted_annotation': True}, f)
        base = shlex.split(METAGRAPH) + ['pattern', '--json', '--pattern-min-information-bits',
                                         str(self.FLOOR), '-i', self.graph_basic, '-a',
                                         self.anno_basic]
        res = subprocess.run(base + ['--pattern-max-predicate-labels', '2', path],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(1, res.returncode)
        self.assertEqual('predicate_too_large', json.loads(res.stdout)['code'])
        res = subprocess.run(base + ['--pattern-max-predicate-labels', '3',
                                     '--pattern-max-predicate-contexts', '7',
                                     '--pattern-max-predicate-work', '1000', path],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(0, res.returncode, res.stderr.decode())
        limits = json.loads(res.stdout)['limits']
        self.assertEqual((7, 1000, 3), (limits['max_predicate_contexts'],
                                        limits['max_predicate_work'],
                                        limits['max_predicate_labels']))

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


class TestPeptideOracle(unittest.TestCase):
    """The six-frame oracle the mini's peptide tests use (cached translations and a residue
    regex) finds what its plain definition finds: no server, no index, always run."""

    def test_six_frames_is_the_reference(self):
        rng = random.Random(5)
        matches = 0
        for i in range(60):
            seq = ''.join(rng.choice('ACGTN' if i % 3 == 0 else 'ACGT')
                          for _ in range(rng.randrange(30, 400)))
            for _ in range(20):
                residues = ''.join(rng.choice('ACDEFGHIKLMNPQRSTVWY*XBZJ')
                                   for _ in range(rng.randrange(1, 6)))
                table = rng.choice(sorted(GENETIC_CODES))
                expected = six_frames_reference(seq, residues, table)
                self.assertEqual(expected, six_frames(seq, residues, table),
                                 (seq, residues, table))
                matches += len(expected)
        self.assertGreater(matches, 1000)


class TestPatternFixtureBodies(unittest.TestCase):
    """The frozen fixture bodies against the SPEC, with no server and no index, so that CI runs
    it; and the blanking of pattern_fixtures.py --check."""

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


# The mini index the fixtures were answered on: its manifest's index_fp. build_mini_refseq.sh
# picks the row-diff anchors anew on every run, so another build of the same records has
# another annotation file, whose row-diff reads decode other paths: what an answer states as
# that decoding's work (WORK_OF_THE_ANCHORS) differs, and nothing else may
FIXTURES_MINI_FP = 'bea44d604ffe08d60b348b630b4a35bd0ae0b428615727345ba191228c7ef5e7'
WORK_OF_THE_ANCHORS = ('annotation_units', 'predicate_units')


@unittest.skipIf(PROTEIN_MODE, "pattern search is DNA only")
@unittest.skipUnless(_supports_pattern(), "`metagraph pattern` is not available in this build")
@guard(os.path.isfile(os.path.join(MINI_DIR, MINI_GRAPH)) and os.path.isfile(FIXTURES_SCRIPT),
       "the mini index (scripts/traversal/build_mini_refseq.sh) or the fixture script is not "
       "there")
class TestPatternFixtures(unittest.TestCase):
    """The frozen fixture bodies the search service tests against
    (api/python/tests/data/traverse/pattern/, SPEC-pattern-search.md §11) are what this binary
    answers: a change of /pattern or of either capabilities route that alters one fails here,
    not first in the service's tests. On the build of the mini index they were answered on,
    every body byte for byte (pattern_fixtures.py --check); on another build of it (CI builds
    its own), every body but the decode work its row-diff anchors make (WORK_OF_THE_ANCHORS)."""

    def test_fixtures_are_what_this_binary_answers(self):
        binary = shlex.split(METAGRAPH)
        if len(binary) != 1:
            self.skipTest('the binary runs under a prefix: the script starts it directly')
        with open(os.path.join(MINI_DIR, MINI_ANNO[:-len('.row_diff_brwt_coord.annodbg')]
                               + '.manifest.json')) as f:
            reference = json.load(f)['index_fp'] == FIXTURES_MINI_FP
        stored = os.path.join(os.path.dirname(FIXTURE_VALIDATOR), 'data', 'traverse', 'pattern')
        with tempfile.TemporaryDirectory() as work:
            # the script writes nothing in the stored fixtures: it builds its copies of the mini
            # index in |work| (the mini itself is only read), starts each server, and compares
            # every body once its run-dependent values are blanked (--check), or writes what
            # differs into a copy of them
            command = [sys.executable, FIXTURES_SCRIPT, '--metagraph', binary[0],
                       '--mini', MINI_DIR, '--work', os.path.join(work, 'indexes')]
            if reference:
                res = subprocess.run(command + ['--check'], stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT)
                out = res.stdout.decode()
                self.assertEqual(0, res.returncode, out[-4000:])
                self.assertRegex(out, r'OK: \d+ fixtures match')
                return
            copy = os.path.join(work, 'fixtures')
            shutil.copytree(stored, copy)
            res = subprocess.run(command + ['--dir', copy], stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT)
            out = res.stdout.decode()
            self.assertEqual(0, res.returncode, out[-4000:])
            # the requests, the README and the index as stored: only answers may differ
            comparison = filecmp.dircmp(stored, copy)
            self.assertEqual([], comparison.diff_files + comparison.left_only
                             + comparison.right_only)
            generator = TestPatternFixtureBodies.load(FIXTURES_SCRIPT, 'pattern_fixtures')

            def without_anchor_work(answer):
                a = generator.blanked(answer)
                for e in a.get('patterns', []) if isinstance(a, dict) else []:
                    for k in WORK_OF_THE_ANCHORS:
                        if k in e.get('work', {}):
                            e['work'][k] = generator.BLANK
                return a

            written = re.findall(r'^wrote (.+)$', out, re.M)
            for path in written:
                name = os.path.basename(os.path.dirname(path))
                with self.subTest(fixture=name):
                    with open(os.path.join(stored, name, 'answer.json')) as f:
                        expected = without_anchor_work(json.load(f))
                    with open(path) as f:
                        self.assertEqual(expected, without_anchor_work(json.load(f)))
                    with open(os.path.join(stored, name, 'request.json')) as f, \
                            open(os.path.join(copy, name, 'request.json')) as g:
                        self.assertEqual(f.read(), g.read())


@unittest.skipIf(PROTEIN_MODE, "the request panel is DNA")
@guard(bool(BASE_BINARY) and os.path.isfile(BASE_BINARY),
       "$METAGRAPH_BASE_BINARY (the binary before the pattern search) is not set (not in CI)")
class TestPatternRegression(TestingBase):
    """/search and /align answer byte for byte as a build without the pattern search does
    (§11; $METAGRAPH_BASE_BINARY), on test_api.py's fixture (transcripts_100.fa, k = 6) and
    requests; /capabilities only gains."""

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
        """Only /pattern withholds its answers from a client that half-closed: /search and
        /align answer such a client byte for byte as the base binary does, their errors
        included."""
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

        def without_reworded(caps):
            # prose whose wording this build changed, the rule it states unchanged: compared as
            # present, not by its text (the rules are references to the SPEC's sections,
            # SPEC-labeled-traversal-core.md §10.3)
            for path in (('attempts', 'delivery_reserve', 'calibration'),
                         ('attempts', 'delivery_reserve', 'rule'), ('attempts', 'bound'),
                         ('attempts', 'not_after'), ('attempts', 'instance'),
                         ('attempts', 'suppression'), ('attempts', 'release_rule'),
                         ('deadline_check', 'rule'), ('decode_cache', 'rule'),
                         ('work_bound',), ('coordinates', 'rule'),
                         ('coordinates', 'output_bound')):
                block = caps
                for key in path[:-1]:
                    block = block.get(key, {})
                if path[-1] in block:
                    self.assertIsInstance(block[path[-1]], str, path)
                    block[path[-1]] = None
            return caps
        def without_poll_stride(caps):
            # the number of the deadline_check rule, a field beside its reference
            self.assertEqual(8, caps['deadline_check'].pop('poll_stride'))
            return caps
        a = without_reworded(without_instance(self.base.get('capabilities').json()))
        b = without_reworded(without_instance(self.new.get('capabilities').json()))
        b = without_poll_stride(b)
        self.assertIn('pattern', b)
        self.assertEqual(a['features'] + ['pattern'], b['features'])
        self.assertEqual(dict(a['routes'], pattern='POST /pattern',
                              pattern_capabilities='GET /pattern/capabilities'), b['routes'])
        for key in a:
            if key not in ('features', 'routes'):
                self.assertEqual(a[key], b[key], key)
        # the stated gains: `pattern` (this route) and `resolve`, the block of /resolve's
        # opt-in bounds.time_budget_ms (SPEC-labeled-traversal-core.md §4.5): a capabilities
        # block, not a feature, since /resolve is listed already; `in_ram` (accepted as /search
        # accepts it) and `graph_summary` (null on a single-graph server), SPEC-pattern-search.md
        # §24
        self.assertEqual(set(a) | {'pattern', 'resolve', 'in_ram', 'graph_summary'}, set(b))
        self.assertIsNone(b['graph_summary'])
        self.assertFalse(b['in_ram']['loads'])
        self.assertTrue(b['resolve']['time_budget']['accepted'])
        self.assertIsNone(b['resolve']['time_budget']['default'])
        # the probe gains the same two blocks and nothing else
        a = without_reworded(without_instance(self.base.get('traverse/capabilities').json()))
        b = without_reworded(without_instance(self.new.get('traverse/capabilities').json()))
        b = without_poll_stride(b)
        new_caps = self.new.get('capabilities').json()
        self.assertEqual(dict(new_caps['pattern'], details='GET /pattern/capabilities'),
                         b.pop('pattern'))
        self.assertEqual(new_caps['resolve'], b.pop('resolve'))
        self.assertEqual(new_caps['in_ram'], b.pop('in_ram'))
        self.assertEqual(a, b)


if __name__ == '__main__':
    unittest.main()
