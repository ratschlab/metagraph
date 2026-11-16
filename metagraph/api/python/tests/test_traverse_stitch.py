"""Stitching a continuation (DESIGN-graphlet-merge.md §3, phase 1a-S; Appendix E.2 T1-T9).

The CLI triples (a parent at radius R, its child from next_request(), the uninterrupted
retrieval at R + E) are in data/traverse/documents/stitch/, written by
scripts/traversal/graphlet_fixtures.py (STITCH_FIXTURES)."""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402,F401

from metagraph.traverse import ops  # noqa: E402


def rc(s):
    return s[::-1].translate(str.maketrans('ACGT', 'TGCA'))


class TestEdgeKey(unittest.TestCase):
    """The walk rule's edge identity (SPEC §6.6, walker.cpp edge_key and step_kmer)."""

    def test_step_orientation(self):
        # natural orientation on both arms: the k-mer then the base on the right, the base
        # then the k-mer on the left
        self.assertEqual('ACGTA', ops.step_of('right', 'ACGT', 'A'))
        self.assertEqual('AACGT', ops.step_of('left', 'ACGT', 'A'))

    def test_basic_is_the_forward_step(self):
        for s in ('AAAAC', 'GTTTT', 'ACGTA'):
            self.assertEqual((s, False), ops.edge_key('basic', s))

    def test_canonical_and_primary_take_the_smaller_orientation(self):
        for regime in ('canonical', 'primary'):
            # rc('GTTTT') = 'AAAAC' < 'GTTTT': the key is the reverse complement
            self.assertEqual(('AAAAC', True), ops.edge_key(regime, 'GTTTT'))
            self.assertEqual(('AAAAC', False), ops.edge_key(regime, 'AAAAC'))
            # one key for a step and its reverse complement, either orientation bit
            for s in ('ACGTTGCA', 'CCGATTAG', 'TTTTTTTA'):
                self.assertEqual(ops.edge_key(regime, s)[0], ops.edge_key(regime, rc(s))[0])
            # a palindrome keeps its forward orientation (the engine flips only when the
            # reverse complement is strictly smaller)
            self.assertEqual(('ACGT', False), ops.edge_key(regime, 'ACGT'))


if __name__ == '__main__':
    unittest.main()
