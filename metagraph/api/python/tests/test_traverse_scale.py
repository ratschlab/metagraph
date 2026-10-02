"""A generated comb-shaped trie (one terminating branch per split, §14's worst case for
materialising ancestor chains): canonical identity, interning and the cost of the
derivations at a few thousand segments."""

import os
import sys
import time
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib  # noqa: E402,F401

from metagraph.traverse import derive, is_canonical, parse  # noqa: E402

COUNTERS = ('steps=%d,successor_enumerations=%d,output_bp=%d,pair_evaluations=0,'
            'edge_reuse_probes=0,reminimisation_rounds=0,max_reminimisation_rounds=0,'
            'refusal_scans=0,switch_sources_cut=0')


def comb_annotate(n):
    """Annotate mode, right arm: a 1-bp root, then n splits, each into a 1-bp leaf
    ('C', {0}) and the next 1-bp spine segment ('G', {0, 1}); the last spine is a leaf.
    Segment ids follow the walker: leaf_i = 2i - 1, spine_i = 2i."""
    steps = 1 + 2 * n
    out = ['H mgt 1 5 basic $ACGT a k k 1000 10 0 walk * * 0123456789abcdef',
           'S aaaaaaaaaaaaaaaa 10 6 0 ACGTACGTAC',
           'L c 0 * 0 left_samples', 'L c 1 * 5 right_samples', 'O c c c i',
           'A r c %d p 0 0 1 2 0 %s * 0 * %d 0 %d %d 0 %d'
           % (n + 1, COUNTERS % (steps, steps, steps), 2 * n + 1, n + 1, n, steps),
           'B 0 2 2 3 1 %d %d 0 %d 0 0 0 0 .' % (steps, n, n),
           'G * 0 1 0-1 * * * %s * A' % ('0' if n else '*'), 'P 0 1 2 !']
    if not n:
        out.append('T D .')
    for i in range(1, n + 1):
        parent = 0 if i == 1 else 2 * (i - 1)
        # {0} against the entry {0}: a tie with '!', written explicitly
        out += ['G %d %d 1 0 * * * * * C' % (parent, i), 'P %d %d 1 0' % (i, i + 1), 'T D .',
                'G %d %d 1 ! * * * %s * G' % (parent, i, '0' if i < n else '*'),
                'P %d %d 2 !' % (i, i + 1)]
        if i == n:
            out.append('T D .')
    out.append('Z %d' % (len(out) + 1))
    return '\n'.join(out) + '\n'


class TestComb(unittest.TestCase):
    def test_small_comb_by_hand(self):
        g = parse(comb_annotate(2))
        a = g.arms['right']
        self.assertEqual([(), (0,), (0,), (2,), (2,)], [s.parents for s in a.segments])
        self.assertEqual([[0, 1], [0, 2, 3], [0, 2, 4]], [p.segments for p in derive.paths(a)])
        self.assertEqual(['C', 'G'], [b['char'] for b in derive.split_branches(
            a, derive.splits(a, g.mode)[0], g.cap)])

    def test_large_comb(self):
        n = 3000
        text = comb_annotate(n)
        t0 = time.perf_counter()
        g = parse(text)
        t_parse = time.perf_counter() - t0
        self.assertTrue(is_canonical(text))
        a = g.arms['right']
        self.assertEqual(n + 1, len(derive.paths(a)))
        t0 = time.perf_counter()
        top = g.walks('right', top=3)
        summary = derive.label_summary(g)
        claims = g.claims('right')
        t_derive = time.perf_counter() - t0
        self.assertEqual([n], [w.path_id for w in top][:1])        # the longest spine first
        self.assertEqual(n + 1, summary[1]['right']['direct_bp'])
        self.assertEqual(n + 1, summary[0]['right']['direct_bp'])
        # left_samples on every walk to its end; right_samples only along the spine
        self.assertEqual({(0, n + 1), (1, n + 1)}, {(c.label.id, c.to_bp) for c in claims
                                                   if c.to_bp == n + 1})
        self.assertLess(g.memory_bytes(), 8 * 1024 * 1024)
        # generous bounds: they catch a quadratic blow-up, not a slow machine
        self.assertLess(t_parse, 5.0, t_parse)
        self.assertLess(t_derive, 10.0, t_derive)


if __name__ == '__main__':
    unittest.main()
