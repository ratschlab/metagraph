import random
import shlex
import subprocess
import unittest

from parameterized import parameterized_class

from base import PROTEIN_MODE, METAGRAPH, start_server
from test_traverse import TraverseServer, GRAPH_MODES, K, random_seq

IUPAC = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT',
         'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}
COMPLEMENT = str.maketrans('ACGTRYSWKMBDHVN', 'TGCAYRSWMKVHDBN')


def rc(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def offsets(kmer, pattern):
    """The offsets at which |kmer| holds |pattern|."""
    return [o for o in range(len(kmer) - len(pattern) + 1)
            if all(kmer[o + i] in IUPAC[c] for i, c in enumerate(pattern))]


class PatternServer(TraverseServer):
    """The fixture's two records and some short ones, whose dummy k-mers a pattern can match."""
    mask = False

    @classmethod
    def _records(cls):
        rng = random.Random(7)
        cls.records = [cls.X + cls.P, cls.X + cls.Q]
        cls.records += [random_seq(rng.randint(K, K + 15), rng) for _ in range(6)]
        return ''.join(f'>{i}\n{record}\n' for i, record in enumerate(cls.records))

    @classmethod
    def _prepare(cls):
        if cls.mask:
            # the graph again, with the file graph.edgemask beside it
            fasta = cls.tempdir.name + '/records.fa'
            cls._build_graph(fasta, cls.graph, K, 'succinct',
                             mode=cls.mode, extra_params='--mask-dummy --in-ram')
            cls._annotate_graph(fasta, cls.graph, cls.tempdir.name + '/annotation', 'column')

    def _pattern(self, request):
        res = self._post('pattern', request)
        self.assertEqual(200, res.status_code, res.text)
        return res.json()

    def _stored_orientation(self, kmer, offset, strand, pattern):
        """The same context on the reverse complement of |kmer|."""
        if rc(pattern) != pattern:
            strand = '-' if strand == '+' else '+'
        return rc(kmer), K - len(pattern) - offset, strand

    def _expected(self, pattern):
        """The contexts of |pattern|, from a scan of the records."""
        kmers = {record[i:i + K] for record in self.records for i in range(len(record) - K + 1)}
        if self.mode == 'canonical':
            kmers |= {rc(kmer) for kmer in kmers}
        contexts = {(kmer, o, '+') for kmer in kmers for o in offsets(kmer, pattern)}
        if self.mode != 'canonical' and rc(pattern) != pattern:
            contexts |= {(kmer, o, '-') for kmer in kmers for o in offsets(kmer, rc(pattern))}
        if self.mode == 'primary':
            contexts = {self._normal(context, pattern) for context in contexts}
        return contexts

    def _normal(self, context, pattern):
        # a primary graph stores one orientation of a k-mer: compare on the smaller one
        other = self._stored_orientation(*context, pattern)
        return min(context, other)

    def _contexts(self, result, pattern):
        contexts = [(c['kmer'], c['offset'], c['strand']) for c in result['contexts']]
        self.assertEqual(len(contexts), len(set(contexts)))
        self.assertEqual(sorted(contexts, key=lambda c: (c[2], c[1])), contexts)
        if self.mode == 'primary':
            return {self._normal(context, pattern) for context in contexts}
        return set(contexts)


@unittest.skipIf(PROTEIN_MODE, 'The patterns of these tests are DNA')
@parameterized_class(('mode', 'mask'), input_values=[(mode, mask) for mode in GRAPH_MODES
                                                     for mask in [False, True]])
class TestPattern(PatternServer):
    def _patterns(self):
        rng = random.Random(3)
        exact = self.X[40:52]
        patterns = [exact, exact[:4] + 'NN' + exact[6:], 'RY' + exact[2:10] + 'SW',
                    self.records[2][:12],      # at the start of a record, beside dummy k-mers
                    self.records[3][-12:],
                    self.records[4][:K], self.records[0][5:5 + K], rc(self.records[0][5:5 + K]),
                    self.X[10:10 + K - 1], 'ACGT', 'ACN', 'acgtnnryacgt',
                    random_seq(12, rng)]
        return patterns + [random_seq(rng.randint(3, K), rng, 'ACGTRYSWKMBDHVN') for _ in range(3)]

    def test_contexts_and_counts(self):
        patterns = self._patterns()
        listed = self._pattern({'patterns': [{'sequence': p} for p in patterns],
                                'max_contexts': 10000})
        counted = self._pattern({'patterns': [{'sequence': p} for p in patterns],
                                 'max_contexts': 0})
        for out in [listed, counted]:
            self.assertEqual(K, out['k'])
            self.assertEqual(self.mode, out['mode'])
            self.assertEqual(self.mask, out['mask'])
            self.assertEqual(list(range(len(patterns))), [r['id'] for r in out['results']])
            self.assertEqual([len(p) for p in patterns], [r['length'] for r in out['results']])

        num_found = 0
        num_bounds = 0
        for pattern, a, b in zip(patterns, listed['results'], counted['results']):
            expected = self._expected(pattern.upper())
            num_found += len(expected)
            # a list is exact with and without a mask
            self.assertTrue(a['complete'], pattern)
            self.assertNotIn('stop', a)
            self.assertEqual(expected, self._contexts(a, pattern.upper()), pattern)
            self.assertEqual({'min': len(expected), 'max': len(expected)}, a['count'], pattern)
            # a count without a list is exact with a mask, and bounds without one
            self.assertTrue(b['complete'], pattern)
            if len(expected):
                self.assertNotIn('contexts', b)
            else:
                self.assertEqual([], b['contexts'])
            self.assertLessEqual(b['count']['min'], len(expected), pattern)
            self.assertLessEqual(len(expected), b['count']['max'], pattern)
            if self.mask:
                self.assertEqual(b['count']['min'], b['count']['max'], pattern)
            num_bounds += b['count']['min'] < b['count']['max']
        self.assertLess(100, num_found)
        if not self.mask:
            # the pattern at the start of a record matches dummy k-mers
            self.assertLess(0, num_bounds)

    def test_list_threshold(self):
        pattern = self.X[40:50]
        num_contexts = len(self._expected(pattern))
        self.assertLess(1, num_contexts)
        for max_contexts, listed in [(num_contexts, True), (num_contexts - 1, False)]:
            result = self._pattern({'patterns': [{'id': 'x', 'sequence': pattern}],
                                    'max_contexts': max_contexts})['results'][0]
            self.assertEqual('x', result['id'])
            self.assertEqual(listed, 'contexts' in result)
            self.assertTrue(result['complete'])
            self.assertLessEqual(result['count']['min'], num_contexts)
            self.assertLessEqual(num_contexts, result['count']['max'])
        # more matches than the default of 1000
        result = self._pattern({'patterns': [{'sequence': 'N'}]})['results'][0]
        self.assertNotIn('contexts', result)
        self.assertLess(1000, result['count']['max'])

    def test_time_budget(self):
        out = self._pattern({'patterns': [{'sequence': 'N'}, {'sequence': self.X[40:52]}],
                             'time_budget_ms': 0})
        for result in out['results']:
            self.assertFalse(result['complete'])
            self.assertEqual('time_budget', result['stop'])
            self.assertIsNone(result['count']['max'])
            self.assertNotIn('contexts', result)

    def test_invalid_requests(self):
        pattern = {'sequence': self.X[:12]}
        for request in [
            {},
            {'patterns': []},
            {'patterns': pattern},
            {'patterns': [pattern] * 17},
            {'patterns': [{}]},
            {'patterns': [{'sequence': ''}]},
            {'patterns': [{'sequence': 12}]},
            {'patterns': [{'sequence': self.X[:K + 1]}]},
            {'patterns': [{'sequence': 'ACGU'}]},
            {'patterns': [{'sequence': 'ACG$'}]},
            {'patterns': [{'sequence': 'ACGT', 'id': 1}]},
            {'patterns': [{'sequence': 'ACGT', 'mismatches': 1}]},
            {'patterns': [pattern], 'max_contexts': 10001},
            {'patterns': [pattern], 'max_contexts': -1},
            {'patterns': [pattern], 'time_budget_ms': 'long'},
            {'patterns': [pattern], 'mode': 'count'},
            {'patterns': [pattern], 'graphs': ['x']},  # no list of graphs here
        ]:
            self.assertEqual(400, self._post('pattern', request).status_code, request)


@unittest.skipIf(PROTEIN_MODE, 'The patterns of these tests are DNA')
class TestPatternLimits(PatternServer):
    """The server's time limit of a request."""
    mode = 'basic'
    server_flags = '--pattern-max-time-ms 0.5'

    def test_time_limit(self):
        # the request's budget is lowered to the server's, which is under a millisecond
        result = self._pattern({'patterns': [{'sequence': 'N'}],
                                'time_budget_ms': 30000})['results'][0]
        self.assertFalse(result['complete'])
        self.assertEqual('time_budget', result['stop'])

    def test_malformed_limit(self):
        for limit in ['0', '-1', 'inf', '1s']:
            command = (f'{METAGRAPH} server_query -i {self.graph} -a {self.annotation}'
                       f' --port {self.port + 1} --pattern-max-time-ms {limit}')
            res = subprocess.run(shlex.split(command), stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE, timeout=60)
            self.assertNotEqual(0, res.returncode, limit)


class TestPatternHashGraph(PatternServer):
    """A graph that is not succinct is refused."""
    mode = 'basic'

    @classmethod
    def _prepare(cls):
        fasta = cls.tempdir.name + '/records.fa'
        cls.hash_graph = cls.tempdir.name + '/hash.orhashdbg'
        cls._build_graph(fasta, cls.hash_graph, K, 'hash', mode=cls.mode)
        cls._annotate_graph(fasta, cls.hash_graph, cls.tempdir.name + '/hash', 'column')

    @classmethod
    def _start_server(cls, port):
        return start_server(cls.hash_graph, cls.tempdir.name + '/hash.column.annodbg',
                            cls.host, port)

    def test_refused(self):
        res = self._post('pattern', {'patterns': [{'sequence': 'ACGT'}]})
        self.assertEqual(400, res.status_code)
        self.assertIn('succinct', res.text)


if __name__ == '__main__':
    unittest.main()
