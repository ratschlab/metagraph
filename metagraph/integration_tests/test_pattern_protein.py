"""Pattern search on a protein graph (docs/DESIGN-pattern-protein.md, SPEC-pattern-search.md §29,
E16): a protein index of UniParc's shape (k 15, records named by their headers, a column
annotation), searched through `metagraph pattern` and the server. The scenarios of the suite's
plan row E16: an exact peptide of at most 15 residues, one with X, one with B, Z or J, a class and
a negated class, mismatches 1 to 3 on a 15-mer, an exact peptide longer than 15 (its paths), a
label predicate, and the refusals; and the fixtures of api/python/tests/data/traverse/pattern_protein.

Every expectation is a brute-force oracle over the records' 15-mers (oracle_contexts,
oracle_paths), which reads a pattern through its own table (sets_of), never through the server.
Skipped in a DNA build (the Protein build's module; test_pattern.py is the DNA build's).

The fixtures are compared without their timings; to rewrite them from this build, run with
METAGRAPH_WRITE_PROTEIN_FIXTURES=1.
"""
import json
import os
import shutil
import subprocess
import unittest
from tempfile import TemporaryDirectory

from base import PROTEIN_MODE, METAGRAPH
from test_pattern import Server, untimed


K = 15
# the 22 residues read as themselves, and the ambiguity letters a record may hold
LITERALS = 'ACDEFGHIKLMNPQRSTVWYUO'
AMBIGUOUS = 'BJZX'
# p1 and p2 share their first 20 residues; RDNPQE[LI]KYSTAM is in p1 (L) and p3 (I); p3 holds U,
# O and the stored ambiguity letters B, X, Z; p4 and p5 hold W15 with 1, 2 and 3 substitutions
W15 = 'EGHYPSTRLLKDNAQ'
RECORDS = {
    'p1': 'MKVLAAGIVGLLLAQSEWHTRDNPQELKYSTAMKVWEGHYPSTRLLKDNAQGF',
    'p2': 'MKVLAAGIVGLLLAQSEWHTKDNPQEIKYSTAMKVWEGHYPSTRAAQ',
    'p3': 'GGSTARDNPQEIKYSTAMUOPWQEBXZQRSTWLLKDNAQGF',
    'p4': 'CCCCC' + 'EGHYPSTRLIKDNAQ' + 'WWWWW' + 'EGWYPSTRLLKDNSQ' + 'CCCCC',
    'p5': 'FFFFF' + 'EGHYPATRLLRDNAE' + 'FFFFF',
}
FIXTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'api', 'python',
                        'tests', 'data', 'traverse', 'pattern_protein')
WRITE_FIXTURES = os.environ.get('METAGRAPH_WRITE_PROTEIN_FIXTURES', '') == '1'
RELEASE = 'fixtures-protein'

# the fixture request: the plan row's scenarios that one request can hold (a request has one
# max_mismatches and one predicate: those are their own tests)
ENTRY_REQUEST = {
    'patterns': [
        {'id': 'exact_11', 'protein': 'DNPQELKYSTA'},
        {'id': 'exact_15', 'protein': W15},
        {'id': 'x', 'protein': 'DNPQEXKYSTA'},
        {'id': 'bzj', 'protein': 'DBPQEJKYSTA'},
        {'id': 'class', 'protein': 'DNPQE[LI]KYSTA'},
        {'id': 'negated', 'protein': 'DNPQE[^P]KYS[^W]A'},
        {'id': 'long', 'protein': 'MKVLAAGIVGLLLAQSEWHT'},
        {'id': 'long_class', 'protein': 'RDNPQE[LI]KYSTAMKVW'},
        {'id': 'dna_kind', 'dna': 'ACGTACGTAC'},
        {'id': 'stop', 'protein': 'DNPQELKY*'},
    ],
    'long_search': 'paths',
}


def letters(c):
    c = c.upper()
    return {'X': LITERALS, 'B': 'DN', 'Z': 'EQ', 'J': 'IL'}.get(c, c if c in LITERALS else '')


def sets_of(text):
    """The oracle's reading of a protein pattern: per position the residues it admits."""
    sets, i = [], 0
    while i < len(text):
        if text[i] == '[':
            j = text.index(']', i)
            inside = text[i + 1:j]
            negated = inside.startswith('^')
            union = set(''.join(letters(c) for c in inside.lstrip('^')))
            sets.append({r for r in LITERALS if (r in union) != negated})
            i = j + 1
        else:
            sets.append(set(letters(text[i])))
            i += 1
    return sets


def kmers_of(records):
    return {r[i:i + K] for r in records.values() for i in range(len(r) - K + 1)}


def oracle_contexts(kmers, text, max_mismatches=0):
    """{(k-mer, offset): substitutions} of a pattern of at most k residues (any_offset): a stored
    ambiguity letter at a pattern position rules the offset out at every level."""
    sets = sets_of(text)
    found = {}
    for kmer in kmers:
        for p in range(K - len(sets) + 1):
            window = kmer[p:p + len(sets)]
            if any(c not in LITERALS for c in window):
                continue
            subs = sum(1 for c, s in zip(window, sets) if c not in s)
            if subs <= max_mismatches:
                found[(kmer, p)] = subs
    return found


def oracle_paths(kmers, text):
    """The anchors (the k-mers instantiating its first k positions) and the walks of a pattern
    longer than k: k-mers overlapping by k - 1, each residue in its position's set."""
    sets = sets_of(text)
    anchors = [x for x in kmers if all(c in s for c, s in zip(x, sets[:K]))]
    walks = []

    def walk(spelled):
        if len(spelled) == len(sets):
            walks.append(spelled)
            return
        for c in sorted(sets[len(spelled)]):
            if spelled[len(spelled) - K + 1:] + c in kmers:
                walk(spelled + c)

    for a in anchors:
        walk(a)
    return anchors, sorted(walks)


def labels_of(kmer):
    return sorted(name for name, r in RECORDS.items() if kmer in r)


@unittest.skipUnless(PROTEIN_MODE, 'pattern search on protein graphs: the Protein build')
class TestPatternProtein(unittest.TestCase):
    """The plan row's scenarios on the masked protein index (counts exact), through the server;
    the CLI and the unmasked graph against the same oracle."""

    @classmethod
    def setUpClass(cls):
        cls.tempdir = TemporaryDirectory()
        d = cls.tempdir.name
        fasta = os.path.join(d, 'records.fa')
        with open(fasta, 'w') as f:
            for name, seq in RECORDS.items():
                f.write(f'>{name}\n{seq}\n')
        cls.graph = os.path.join(d, 'graph.dbg')
        cls.anno = os.path.join(d, 'annotation.column.annodbg')
        cls.unmasked = os.path.join(d, 'unmasked.dbg')
        for command in (
                [METAGRAPH, 'build', '-k', str(K), '-o', os.path.join(d, 'graph'), fasta],
                [METAGRAPH, 'annotate', '-i', cls.graph, '--anno-header',
                 '-o', os.path.join(d, 'annotation'), fasta]):
            subprocess.run(command, check=True, stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL)
        shutil.copy(cls.graph, cls.unmasked)
        subprocess.run([METAGRAPH, 'transform', '--mask-dummy', cls.graph], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        assert os.path.exists(os.path.join(d, 'graph.edgemask'))
        cls.kmers = kmers_of(RECORDS)
        cls.server = Server(METAGRAPH, ['-i', cls.graph, '-a', cls.anno,
                                        '--index-release', RELEASE],
                            os.path.join(d, 'server.log'))

    @classmethod
    def tearDownClass(cls):
        cls.server.stop()
        cls.tempdir.cleanup()

    def pattern(self, payload, status=200):
        ret = self.server.post('pattern', payload)
        self.assertEqual(status, ret.status_code, ret.text)
        return ret.json()

    def cli(self, graph, payload, *flags):
        path = os.path.join(self.tempdir.name, 'request.json')
        with open(path, 'w') as f:
            json.dump(payload, f)
        out = subprocess.run([METAGRAPH, 'pattern', '-i', graph, '-a', self.anno,
                              '--index-release', RELEASE, *flags, path],
                             check=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        return json.loads(out.stdout)

    def assertContexts(self, entry, expected):
        """counts.contexts and the results against the oracle {(k-mer, offset): level}."""
        self.assertNotIn('error', entry)
        c = entry['counts']['contexts']
        self.assertEqual((len(expected), 'exact'), (c['value'], c['relation']))
        self.assertEqual(['+'], list(c['by_strand']))
        got = {(r['kmer'], r['offset']): r.get('mismatches', 0) for r in entry['results']}
        self.assertEqual(expected, got)
        for r in entry['results']:
            self.assertEqual('+', r['strand'])
            self.assertEqual(r['kmer'][r['offset']:r['offset'] + entry['length']], r['instance'])

    def test_capabilities(self):
        caps = self.server.get('pattern/capabilities').json()
        self.assertTrue(caps['available'])
        self.assertEqual('$ABCDEFGHIJKLMNOPQRSTUVWYZX', caps['alphabet'])
        self.assertEqual(['protein'], caps['kinds'])
        self.assertEqual('residues', caps['protein_reading'])
        self.assertEqual(['[...]', '[^...]'], caps['protein_classes'])
        self.assertEqual(list(LITERALS + 'XBZJ'), caps['protein_residues'])
        self.assertEqual([], caps['genetic_codes'])
        self.assertIsNone(caps['default_genetic_code'])
        self.assertEqual(['both', 'forward'], caps['strands'])
        self.assertEqual({'basic'}, set(caps['scopes_by_graph_mode']))
        self.assertEqual({'protein': 3}, caps['mismatches']['max'])
        self.assertEqual({'protein': 'residue'}, caps['mismatches']['unit'])
        self.assertEqual('residues', caps['mismatches']['protein_substitutes'])
        self.assertEqual('exact', caps['counting'])
        self.assertEqual(K, caps['k'])
        # GET /capabilities' pattern block and the gate block of GET /traverse/capabilities
        full = self.server.get('capabilities').json()['pattern']
        self.assertEqual(caps['kinds'], full['kinds'])
        self.assertEqual('residues', full['protein_reading'])
        gate = self.server.get('traverse/capabilities').json()
        gate = gate.get('pattern', gate)
        self.assertEqual(['protein'], gate['kinds'])
        self.assertEqual([], gate['genetic_codes'])

    def test_exact_sets_and_classes(self):
        a = self.pattern({'patterns': [p for p in ENTRY_REQUEST['patterns'][:6]]})
        self.assertEqual('$ABCDEFGHIJKLMNOPQRSTUVWYZX', a['index']['alphabet'])
        for spec, entry in zip(ENTRY_REQUEST['patterns'][:6], a['patterns']):
            with self.subTest(spec['id']):
                text = spec['protein']
                self.assertEqual('protein', entry['kind'])
                self.assertEqual(text.upper(), entry['pattern'])
                self.assertEqual(len(sets_of(text)), entry['length'])
                self.assertEqual(entry['length'], entry['residues'])
                self.assertIsNone(entry['genetic_code'])
                self.assertEqual(['+'], entry['strands'])
                self.assertFalse(entry['palindromic'])
                expected = oracle_contexts(self.kmers, text)
                self.assertContexts(entry, expected)
                self.assertGreater(len(expected), 0)
        # the 15-mer is one context (offset 0), the 11-mers several offsets
        self.assertEqual(1, a['patterns'][1]['counts']['contexts']['value'])

    def test_stored_ambiguity_is_matched_by_nothing(self):
        # p3 holds EBXZQ: no position matches its B, X or Z, X itself and a negation included
        a = self.pattern({'patterns': [{'protein': 'WQEBXZQRS'}, {'protein': 'WQEXXXQRS'},
                                       {'protein': 'WQE[^P][^P][^P]QRS'}]})
        for entry in a['patterns']:
            self.assertContexts(entry, {})

    def test_mismatches_1_to_3_on_a_15_mer(self):
        previous = {}
        for m in (1, 2, 3):
            with self.subTest(max_mismatches=m):
                a = self.pattern({'patterns': [{'protein': W15}], 'max_mismatches': m})
                entry = a['patterns'][0]
                expected = oracle_contexts(self.kmers, W15, m)
                self.assertContexts(entry, expected)
                self.assertEqual(m, entry['max_mismatches'])
                by = entry['counts']['contexts']['by_mismatches']
                levels = [by[str(u)]['value'] for u in range(m + 1)]
                self.assertEqual([sum(1 for t in expected.values() if t == u)
                                  for u in range(m + 1)], levels)
                for r in entry['results']:
                    subs = r['substitutions']
                    self.assertEqual(r['mismatches'], len(subs))
                    for s in subs:
                        self.assertEqual({'position', 'expected', 'found'}, set(s))
                        self.assertEqual(W15[s['position']], s['expected'])
                        self.assertEqual(r['instance'][s['position']], s['found'])
                        self.assertNotEqual(s['expected'], s['found'])
                # each level adds the variants one more substitution away
                self.assertTrue(set(previous) <= set(expected))
                self.assertGreater(len(expected), len(previous))
                previous = expected
        # a class is substituted as written
        a = self.pattern({'patterns': [{'protein': 'EGHYPSTRL[ST]KDNAQ'}], 'max_mismatches': 1})
        subs = [s for r in a['patterns'][0]['results'] for s in r['substitutions']]
        self.assertIn({'position': 9, 'expected': '[ST]', 'found': 'L'}, subs)
        # mismatches in a pattern longer than k: its slot, cap 0
        a = self.pattern({'patterns': [{'protein': W15 + 'G'}], 'max_mismatches': 1})
        self.assertEqual('mismatches_above_cap', a['patterns'][0]['error']['code'])
        self.assertEqual(0, a['patterns'][0]['max_mismatches_cap'])

    def test_long_exact_peptides_and_their_paths(self):
        for text in ('MKVLAAGIVGLLLAQSEWHT', 'RDNPQE[LI]KYSTAMKVW', 'GGSTARDNPQEIKYSTAMUOPWQ',
                     'MKVLAAGIVGLLLAQSEWHTW'):
            with self.subTest(text):
                a = self.pattern({'patterns': [{'protein': text}], 'long_search': 'paths'})
                entry = a['patterns'][0]
                anchors, walks = oracle_paths(self.kmers, text)
                self.assertEqual('long', entry['scope'])
                self.assertEqual(len(sets_of(text)), entry['length'])
                self.assertEqual(len(anchors), entry['counts']['anchors']['value'])
                self.assertEqual((len(walks), 'exact'), (entry['counts']['paths']['value'],
                                                         entry['counts']['paths']['relation']))
                self.assertEqual(walks, sorted(r['sequence'] for r in entry['results']))
                for r in entry['results']:
                    self.assertEqual(entry['length'] - K + 1, len(r['nodes']))
        # and the anchors alone, without long_search
        a = self.pattern({'patterns': [{'protein': 'MKVLAAGIVGLLLAQSEWHT'}]})
        self.assertIn('paths_later_increment', a['patterns'][0]['notes'])

    def test_supported_paths_and_a_predicate(self):
        """long_search "supported_paths": the walks some record carries along their whole length
        (label_intersection: this index has no coordinates), one strand; with a predicate the
        selected ones, its strands evaluated and stated as "context"."""
        texts = ('MKVLAAGIVGLLLAQSEWHT', 'RDNPQE[LI]KYSTAMKVW', 'DNPQE[LI]KYSTAMKVWEGHYPSTR')

        def carriers(walk):
            kmers = [walk[i:i + K] for i in range(len(walk) - K + 1)]
            return sorted(set.intersection(*(set(labels_of(x)) for x in kmers)))

        a = self.pattern({'patterns': [{'protein': t} for t in texts],
                          'long_search': 'supported_paths', 'output': {'labels': 'all'},
                          'allow_unbudgeted_annotation': True})
        for text, entry in zip(texts, a['patterns']):
            with self.subTest(text):
                _, walks = oracle_paths(self.kmers, text)
                supported = [w for w in walks if carriers(w)]
                self.assertEqual(len(walks), entry['counts']['paths']['value'])
                sp = entry['counts']['supported_paths']
                self.assertEqual((len(supported), 'exact'), (sp['value'], sp['relation']))
                self.assertEqual(['+'], list(sp['by_strand']))
                self.assertEqual(supported, sorted(r['sequence'] for r in entry['results']))
                for r in entry['results']:
                    self.assertEqual('+', r['strand'])
                    self.assertEqual(carriers(r['sequence']),
                                     sorted(l['column'] for l in r['labels']))
        a = self.pattern({'patterns': [{'protein': t} for t in texts],
                          'long_search': 'supported_paths', 'predicate': {'any': ['p2']},
                          'predicate_strands': 'either', 'allow_unbudgeted_annotation': True})
        self.assertEqual('context', a['predicate']['strands'])
        self.assertEqual('context', a['limits']['predicate_strands'])
        for text, entry in zip(texts, a['patterns']):
            with self.subTest(text, predicate='p2'):
                _, walks = oracle_paths(self.kmers, text)
                selected = [w for w in walks if 'p2' in carriers(w)]
                self.assertEqual(len(selected), entry['counts']['selected']['value'])
                self.assertEqual(selected, sorted(r['sequence'] for r in entry['results']))

    def test_label_predicate(self):
        # asked with "either": evaluated and stated as "context" (one strand)
        a = self.pattern({'patterns': [{'protein': 'DNPQE[LI]KYSTA'}, {'protein': W15}],
                          'predicate': {'any': ['p1']}, 'predicate_strands': 'either',
                          'output': {'labels': 'predicate_only'},
                          'allow_unbudgeted_annotation': True})
        self.assertEqual('context', a['predicate']['strands'])
        self.assertEqual('context', a['limits']['predicate_strands'])
        for spec, entry in zip(('DNPQE[LI]KYSTA', W15), a['patterns']):
            with self.subTest(spec):
                raw = oracle_contexts(self.kmers, spec)
                selected = {key for key in raw if 'p1' in labels_of(key[0])}
                self.assertEqual(len(raw), entry['counts']['contexts']['value'])
                self.assertEqual(len(selected), entry['counts']['selected']['value'])
                self.assertEqual(selected, {(r['kmer'], r['offset']) for r in entry['results']})
                for r in entry['results']:
                    self.assertEqual(['p1'], r['selection_labels'])
                    self.assertEqual(['context'], r['selection_strands'])
        # the labels of every context
        a = self.pattern({'patterns': [{'protein': 'DNPQE[LI]KYSTA'}],
                          'output': {'labels': 'all'}, 'allow_unbudgeted_annotation': True})
        for r in a['patterns'][0]['results']:
            self.assertEqual(labels_of(r['kmer']), sorted(l['column'] for l in r['labels']))

    def test_refusals(self):
        a = self.pattern({'patterns': [{'dna': 'ACGTACGTAC'}, {'iupac': 'ACGTN'},
                                       {'protein': 'DNPQELKY*'}, {'protein': 'DNPQ[]'},
                                       {'protein': 'DNPQ-ELKY'}, {'protein': 'DNPQELKYSTA'}]})
        for entry in a['patterns'][:5]:
            self.assertEqual('bad_alphabet', entry['error']['code'], entry)
        self.assertIn('protein graph', a['patterns'][0]['error']['message'])
        self.assertNotIn('error', a['patterns'][5])
        for body in ({'patterns': [{'protein': 'DNPQELKYSTA'}], 'genetic_code': 1},
                     {'patterns': [{'protein': 'DNPQELKYSTA'}], 'strands': 'reverse'}):
            with self.subTest(body):
                ret = self.server.post('pattern', body)
                self.assertEqual(400, ret.status_code, ret.text)
                self.assertEqual('invalid_request', ret.json()['code'])
        a = self.pattern({'patterns': [{'protein': 'DNPQELKYSTA'}], 'strands': 'forward'})
        self.assertEqual(['+'], a['patterns'][0]['strands'])

    def test_cli_answers_as_the_server(self):
        self.assertEqual(untimed(self.pattern(ENTRY_REQUEST)),
                         untimed(self.cli(self.graph, ENTRY_REQUEST)))

    def test_unmasked_graph_counts_upper_bounds(self):
        for checked in ('50', '0'):
            a = self.cli(self.unmasked, {'patterns': ENTRY_REQUEST['patterns'][:6],
                                         'mode': 'count'},
                         '--pattern-max-checked-entries', checked)
            self.assertEqual('upper_bound', a['index']['counting'])
            for spec, entry in zip(ENTRY_REQUEST['patterns'][:6], a['patterns']):
                with self.subTest(spec['id'], checked=checked):
                    truth = len(oracle_contexts(self.kmers, spec['protein']))
                    c = entry['counts']['contexts']
                    if c['relation'] == 'exact':
                        self.assertEqual(truth, c['value'])
                    else:
                        self.assertEqual('bounds', c['relation'])
                        self.assertLessEqual(c['lower'], truth)
                        self.assertGreaterEqual(c['upper'], truth)

    def test_fixtures(self):
        """api/python/tests/data/traverse/pattern_protein: the capabilities and one answer, as
        this build answers them (timings left out)."""
        documents = {
            os.path.join('capabilities', 'answer.json'):
                untimed(self.server.get('pattern/capabilities').json()),
            os.path.join('entry', 'request.json'): ENTRY_REQUEST,
            os.path.join('entry', 'answer.json'): untimed(self.pattern(ENTRY_REQUEST)),
        }
        for name, document in documents.items():
            path = os.path.join(FIXTURES, name)
            text = json.dumps(document, indent=1, sort_keys=True) + '\n'
            if WRITE_FIXTURES:
                os.makedirs(os.path.dirname(path), exist_ok=True)
                with open(path, 'w') as f:
                    f.write(text)
                continue
            with self.subTest(name):
                with open(path) as f:
                    self.assertEqual(f.read(), text)


if __name__ == '__main__':
    unittest.main()
