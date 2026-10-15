"""Offline tests of the lane O oracle's definitions (oracle/oracle.py) on hand-checked cases,
several of them the SPEC's own examples. No server, no engine."""

import ast
import math
import os
import sys
import unittest

SUITE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, SUITE)

from oracle import fixtures  # noqa: E402
from oracle import oracle as O  # noqa: E402
from oracle.fixtures import Record  # noqa: E402


def index(seqs, k, mode='basic'):
    """An index of one column per record ('c<i>', header 'r<i>') unless |seqs| holds Records."""
    recs = [s if isinstance(s, Record) else Record(f'c{i}', 0, f'r{i}', s)
            for i, s in enumerate(seqs)]
    return O.Index(recs, k, mode)


def ids(ctxs):
    return {(c.orientation, c.kmer, c.offset) for c in ctxs}


class Contexts(unittest.TestCase):
    def test_palindrome_is_searched_once(self):
        idx = index(['AACGTT'], 5)
        p = O.Pattern('dna', 'CG')
        self.assertTrue(p.palindromic)
        ctxs = O.contexts(idx, p, 'any_offset', 'forward')
        self.assertEqual({('=', 'AACGT', 2), ('=', 'ACGTT', 1)}, ids(ctxs))
        counts = O.counts(idx, p, ctxs, 'any_offset', 'forward', 0)
        self.assertEqual({'0': 0, '1': 1, '2': 1, '3': 0}, counts['by_offset'])
        self.assertEqual({'both': 2}, counts['by_orient'])
        self.assertEqual(0, counts['suffix'])

    def test_both_strands_and_answer_order(self):
        idx = index(['AACGTT'], 5)
        p = O.Pattern('dna', 'AAC')
        ctxs = O.contexts(idx, p, 'any_offset', 'both')
        # rc(AAC) = GTT at offset 2 of ACGTT; AACGT precedes ACGTT in BOSS order
        self.assertEqual([('+', 'AACGT', 0), ('-', 'ACGTT', 2)],
                         [(c.orientation, c.kmer, c.offset) for c in ctxs])
        self.assertEqual([('-', 'ACGTT', 2)],
                         [c.identity for c in O.contexts(idx, p, 'suffix', 'both')])
        self.assertEqual([], O.contexts(idx, p, 'suffix', 'forward'))

    def test_an_iupac_pattern_and_its_reverse_complement_at_one_offset(self):
        # SPEC 7.1: NA and TN both match TA
        idx = index(['GGTACC'], 5)
        ctxs = O.contexts(idx, O.Pattern('iupac', 'NA'), 'any_offset', 'both')
        self.assertEqual({('+', 'GGTAC', 2), ('-', 'GGTAC', 2), ('+', 'GTACC', 1),
                          ('-', 'GTACC', 1)}, ids(ctxs))

    def test_one_substitution_from_both_orientations(self):
        # SPEC 21.4: AACGTT is one substitution from AACGTA and from its reverse complement
        idx = index(['AACGTT'], 6)
        ctxs = O.contexts(idx, O.Pattern('dna', 'AACGTA'), 'any_offset', 'both', 1)
        self.assertEqual([('+', 'AACGTT', 0, 1, [(5, 'A', 'T')]),
                          ('-', 'AACGTT', 0, 1, [(5, 'A', 'T')])],
                         [(c.orientation, c.kmer, c.offset, c.level, c.subs) for c in ctxs])
        self.assertEqual([], O.contexts(idx, O.Pattern('dna', 'AACGTA'), 'any_offset', 'both', 0))

    def test_n_never_holds_a_substitution(self):
        idx = index(['ACGTACGTAA'], 5)
        p = O.Pattern('iupac', 'NNNNN')
        self.assertEqual([0] * len(idx.kmers),
                         [c.level for c in O.contexts(idx, p, 'any_offset', 'forward', 3)])

    def test_records_are_split_at_n_and_short_records_hold_nothing(self):
        # ACGTA | N | CGTAC: two islands of one k-mer each, none across the N
        idx = index(['ACGTANCGTAC', 'ACG'], 5)
        self.assertEqual({'ACGTA', 'CGTAC'}, idx.kmers)
        canon = index(['ACGTANCGTAC'], 5, 'canonical')
        self.assertEqual({'ACGTA', 'TACGT', 'CGTAC', 'GTACG'}, canon.kmers)

    def test_canonical_orientations(self):
        idx = index(['AACGTT'], 5, 'canonical')
        ctxs = O.contexts(idx, O.Pattern('dna', 'AAC'), 'any_offset', 'both')
        # both orientations of every k-mer: forward and reverse count alike in any_offset
        fwd = [c for c in ctxs if c.orientation == 'forward']
        rev = [c for c in ctxs if c.orientation == 'reverse']
        self.assertEqual(len(fwd), len(rev))
        self.assertIn(('forward', 'AACGT', 0), ids(ctxs))


class Peptides(unittest.TestCase):
    def test_codon_sets(self):
        self.assertEqual({'ATG'}, O.residue_codons('M', 1))
        self.assertEqual({'TGG'}, O.residue_codons('W', 1))
        self.assertEqual({'TGA', 'TGG'}, O.residue_codons('W', 4))
        self.assertEqual({'TAA', 'TAG', 'TGA'}, O.residue_codons('*', 1))
        self.assertEqual({'TAA', 'TAG', 'AGA', 'AGG'}, O.residue_codons('*', 2))
        self.assertEqual(set(), O.residue_codons('*', 27))
        self.assertEqual(61, len(O.residue_codons('X', 1)))
        self.assertEqual(O.residue_codons('D', 1) | O.residue_codons('N', 1),
                         O.residue_codons('B', 1))
        self.assertEqual(O.residue_codons('M', 1), O.residue_codons('M', 11))

    def test_contexts_on_both_strands(self):
        idx = index(['CATGC'], 5)
        ctxs = O.contexts(idx, O.Pattern('protein', 'M', 1), 'any_offset', 'both')
        self.assertEqual({('+', 'CATGC', 1), ('-', 'CATGC', 0)}, ids(ctxs))

    def test_sense_substitutes(self):
        # SPEC 21.3: a residue's codon outside its set is one substitution if it is a sense
        # codon; a stop codon never substitutes; X and * are never substituted
        idx = index(['CTGGC', 'CTGAC'], 5)
        p = O.Pattern('protein', 'M', 1)
        ctxs = O.contexts(idx, p, 'any_offset', 'forward', 1)
        got = {(c.kmer, c.offset): c.subs for c in ctxs}
        self.assertEqual([(0, 'M', 'L', 'CTG')], got[('CTGGC', 0)])
        self.assertEqual([(0, 'M', 'W', 'TGG')], got[('CTGGC', 1)])
        self.assertNotIn(('CTGAC', 1), got)              # TGA: a stop, no substitute
        x = O.contexts(idx, O.Pattern('protein', 'X', 1), 'any_offset', 'forward', 1)
        self.assertNotIn(('CTGAC', 1), {(c.kmer, c.offset) for c in x})

    def test_a_peptide_without_instances(self):
        p = O.Pattern('protein', 'M*', 27)
        self.assertTrue(p.no_instance)
        self.assertEqual([], O.contexts(index(['ATGTAA'], 6), p, 'any_offset', 'both', 1))
        # its bits count the residue without codons as one exact codon (SPEC 12.2)
        self.assertAlmostEqual(12.0, p.information_bits())

    def test_palindromic_peptides(self):
        # only X runs in tables without a stop codon are palindromic (SPEC 12.2)
        self.assertTrue(O.Pattern('protein', 'XX', 27).palindromic)
        self.assertFalse(O.Pattern('protein', 'XX', 1).palindromic)
        self.assertFalse(O.Pattern('protein', 'MW', 1).palindromic)


class Information(unittest.TestCase):
    def test_iupac_bits(self):
        self.assertAlmostEqual(8.0, O.Pattern('dna', 'ACGT').information_bits())
        self.assertAlmostEqual(1.0, O.Pattern('iupac', 'RN').information_bits())
        self.assertAlmostEqual(math.log2(4 / 3), O.Pattern('iupac', 'B').information_bits())

    def test_mismatch_bits_of_the_spec(self):
        # SPEC 21.5: an exact 31-mer, |B| = 94, 4,279, 125,644
        p = O.Pattern('dna', 'A' * 31)
        for m, ball in ((1, 94), (2, 4279), (3, 125644)):
            self.assertAlmostEqual(62 - math.log2(ball), p.mismatch_information_bits(m))
        self.assertAlmostEqual(p.information_bits(), p.mismatch_information_bits(0))
        # the decapeptide KDVKFGNDAR (table 1) loses about 7.8, 14.5 and 20.4 bits
        q = O.Pattern('protein', 'KDVKFGNDAR', 1)
        for m, loss in ((1, 7.8), (2, 14.5), (3, 20.4)):
            self.assertAlmostEqual(loss, q.information_bits() - q.mismatch_information_bits(m),
                                   delta=0.06)


class Labels(unittest.TestCase):
    def test_placement_of_the_spec(self):
        # SPEC 14.3: k = 5, a column of ACGTA, CCCCC, CCCCC; GTA at offset 2 of ACGTA is record
        # 0, 3-5
        recs = [Record('c', 0, 'r0', 'ACGTA'), Record('c', 1, 'r1', 'CCCCC'),
                Record('c', 2, 'r2', 'CCCCC')]
        idx = O.Index(recs, 5, 'basic')
        self.assertEqual([(0, 'r0', '+', '3-5', 5)], idx.occurrences('ACGTA', 'c', 2, 3, '+'))
        p = O.Pattern('dna', 'CCC')
        ctxs = O.contexts(idx, p, 'any_offset', 'forward')
        view = O.label_view(idx, p, ctxs)
        self.assertEqual([{'column': 'c', 'contexts': 3, 'contexts_suffix': 1, 'occurrences': 6}],
                         view['by_label'])
        self.assertEqual(6, view['occurrences'])

    def test_a_record_shorter_than_k_is_not_indexed(self):
        recs = [Record('c', 0, 'short', 'ACG'), Record('c', 1, 'long', 'ACGTA')]
        idx = O.Index(recs, 5, 'basic')
        self.assertEqual([(0, 'long', '+', '3-5', 5)], idx.occurrences('ACGTA', 'c', 2, 3, '+'))
        with self.assertRaises(ValueError):
            O.resolve(idx, 'ACGTA', ['short'])

    def test_label_order_is_contexts_desc_then_name(self):
        # AAAA: c_b holds AAAAA (offsets 0, 1), AAAAC, CAAAA, AAAAT (5 contexts); c_a AAAAA
        # (2) and AAAAG (3 contexts): c_b first, although c_a precedes it by name
        recs = [Record('c_b', 0, 'x', 'AAAAACAAAAT'), Record('c_a', 0, 'y', 'AAAAAG')]
        idx = O.Index(recs, 5, 'basic')
        p = O.Pattern('dna', 'AAAA')
        view = O.label_view(idx, p, O.contexts(idx, p, 'any_offset', 'forward'))
        self.assertEqual(['c_b', 'c_a'], [b['column'] for b in view['by_label']])
        # a tie (3 contexts each): by name
        recs = [Record('c_b', 0, 'x', 'AAAAAC'), Record('c_a', 0, 'y', 'AAAAAG')]
        idx = O.Index(recs, 5, 'basic')
        view = O.label_view(idx, p, O.contexts(idx, p, 'any_offset', 'forward'))
        self.assertEqual(['c_a', 'c_b'], [b['column'] for b in view['by_label']])

    def test_canonical_labels_hold_both_strands(self):
        idx = index(['AACGTA'], 5, 'canonical')
        self.assertEqual({'c0'}, idx.labels_of(O.revcomp('AACGT')))


class Resolve(unittest.TestCase):
    def test_runs_and_candidates(self):
        recs = [Record('c1', 0, 'r1', 'AAAAACCCCC'), Record('c2', 0, 'r2', 'CCCCCGGGGG')]
        idx = O.Index(recs, 5, 'basic')
        out = O.resolve(idx, 'AAAAACCCCCGGGGG', ['c2', 'c1', 'r1'])
        self.assertEqual(11, out['num_kmers'])
        self.assertEqual([[0, 11]], out['graph_runs'])
        self.assertEqual([{'label': 'c2', 'kind': 'column', 'kmers_supported': 6,
                           'runs': [[5, 11]]},
                          {'label': 'c1', 'kind': 'column', 'kmers_supported': 6,
                           'runs': [[0, 6]]},
                          {'label': 'r1', 'kind': 'header', 'kmers_supported': 6,
                           'runs': [[0, 6]]}], out['labels'])
        # longer first, then more labels, then the smaller start; labels in request order
        self.assertEqual([{'kmer_interval': [0, 6], 'bp_interval': [0, 10],
                           'labels': ['c1', 'r1']},
                          {'kmer_interval': [5, 11], 'bp_interval': [5, 15], 'labels': ['c2']}],
                         out['candidates'])
        self.assertEqual([], O.resolve(idx, 'AAAAACCCCC', ['c1'], min_block_kmers=7)['candidates'])

    def test_basic_resolves_the_forward_strand_only(self):
        idx = index(['AAAAACCCCC'], 5)
        self.assertEqual([], O.resolve(idx, O.revcomp('AAAAACCCCC'), ['c0'])['graph_runs'])
        canon = index(['AAAAACCCCC'], 5, 'canonical')
        self.assertEqual([[0, 6]],
                         O.resolve(canon, O.revcomp('AAAAACCCCC'), ['c0'])['labels'][0]['runs'])

    def test_trace_breaks_where_coordinates_jump(self):
        idx = index(['ACGTATTCGTAC'], 5)
        kmer = O.resolve(idx, 'ACGTAC', ['r0'])
        trace = O.resolve(idx, 'ACGTAC', ['r0'], support='trace')
        self.assertEqual([[0, 2]], kmer['labels'][0]['runs'])
        self.assertEqual([[0, 1], [1, 2]], trace['labels'][0]['runs'])
        self.assertEqual(2, trace['labels'][0]['kmers_supported'])
        # a repeat: the chain follows either copy (a set of live coordinates)
        rep = index(['ACGTACGTAC'], 5)
        self.assertEqual([[0, 6]],
                         O.resolve(rep, 'ACGTACGTAC', ['r0'], support='trace')['labels'][0]['runs'])

    def test_a_column_trace_crosses_a_record_boundary(self):
        # column coordinates concatenate the column's records: the last k-mer of one record and
        # the first of the next are consecutive (trace_record_boundaries)
        # a: GAAAA AAAAC (coordinates 0, 1); b: AAACG AACGT (2, 3)
        recs = [Record('c', 0, 'a', 'GAAAAC'), Record('c', 1, 'b', 'AAACGT')]
        idx = O.Index(recs, 5, 'basic')
        out = O.resolve(idx, 'GAAAACGT', ['c', 'a', 'b'], support='trace')
        self.assertEqual([[[0, 4]], [[0, 2]], [[2, 4]]], [l['runs'] for l in out['labels']])


class Order(unittest.TestCase):
    def test_colex_order(self):
        self.assertEqual(['AAAAC', 'ACGTA', 'TTTTA'],
                         sorted(['ACGTA', 'TTTTA', 'AAAAC'], key=O.colex_key))
        # the source node first, then the edge label
        self.assertEqual(['CAAAA', 'CAAAC'], sorted(['CAAAC', 'CAAAA'], key=O.colex_key))


class Fixtures(unittest.TestCase):
    def test_the_records_are_fixed_and_hold_their_features(self):
        cols = fixtures.build_records()
        self.assertEqual(fixtures.COLUMN_ORDER, list(cols))
        recs = {r.name: r for r in fixtures.all_records(cols)}
        self.assertIn('GAATTC', recs['alpha'].seq)
        self.assertIn('N', recs['delta1'].seq)
        self.assertLess(len(recs['delta3'].seq), 11)
        self.assertLess(len(recs['theta1'].seq), 11)
        self.assertIn('A' * 24, recs['gamma'].seq)
        self.assertIn(O.revcomp(recs['alpha'].seq[96:146]), recs['gamma'].seq)
        # eps: six planted substitutions against alpha[100:200]
        diff = [i for i, (a, b) in enumerate(zip(recs['alpha'].seq[100:200], recs['eps'].seq))
                if a != b]
        self.assertEqual([10, 40, 44, 70, 73, 77], diff)
        # the SNP bubble: alpha and beta differ in A's middle base only
        alpha_a = recs['alpha'].seq[:30] + recs['alpha'].seq[36:96]
        self.assertEqual(1, sum(a != b for a, b in zip(alpha_a, recs['beta'].seq[:90])))
        def flat(c):
            return [(r.column, r.seq_id, r.name, r.seq) for r in fixtures.all_records(c)]
        self.assertEqual(flat(cols), flat(fixtures.build_records()))
        # the annotation order is not the name order
        self.assertNotEqual(sorted(fixtures.COLUMN_ORDER), fixtures.COLUMN_ORDER)


class Independence(unittest.TestCase):
    def test_the_oracle_imports_only_the_standard_library(self):
        """The oracle shares no code with the engine, checks.py, seqtools.py, truth.py or the
        integration tests: oracle.py imports math alone, fixtures.py random and os."""
        allowed = {'oracle.py': {'math'}, 'fixtures.py': {'random', 'os'}}
        for name, ok in allowed.items():
            tree = ast.parse(open(os.path.join(SUITE, 'oracle', name)).read())
            mods = set()
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    mods |= {a.name.split('.')[0] for a in node.names}
                elif isinstance(node, ast.ImportFrom):
                    mods.add(node.module.split('.')[0] if node.module else '.')
            self.assertLessEqual(mods, ok, name)

    def test_no_module_of_the_package_imports_the_suite_or_the_engine(self):
        banned = ('checks', 'seqtools', 'truth', 'suite_common', 'metagraph', 'helpers', 'base')
        for name in os.listdir(os.path.join(SUITE, 'oracle')):
            if not name.endswith('.py'):
                continue
            tree = ast.parse(open(os.path.join(SUITE, 'oracle', name)).read())
            for node in ast.walk(tree):
                mods = []
                if isinstance(node, ast.Import):
                    mods = [a.name for a in node.names]
                elif isinstance(node, ast.ImportFrom) and node.module:
                    mods = [node.module]
                for m in mods:
                    self.assertNotIn(m.split('.')[0], banned, (name, m))


if __name__ == '__main__':
    unittest.main()
