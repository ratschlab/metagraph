"""The library defects the v8 test suite found (FORMAT-v8 §G), each reproduced offline on a
fixture of data/traverse/documents. One class per finding, named by its id there."""

import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402

from metagraph.traverse import GraphletStore, parse  # noqa: E402
from metagraph.traverse.mcp_tools import GraphletTools  # noqa: E402

MODES = ('claims', 'walks', 'labels', 'prefix_subset')
DEPTH_RULE = 'no common certified depth'


def _with_fp(text, fp='cd' * 32):
    """A CLI fixture with an index manifest digest (they are made without one)."""
    return text.replace(' walk fixtures * ', ' walk fixtures %s ' % fp, 1)


class TestG1TheReasonNamesTheRuleThatDecided(unittest.TestCase):
    """G1 (O1): no_bins (no index_fp, complete_to_bp 0 on both arms) compares 'unknown' by
    the depth rule, but its reason named the missing index_fp, the rule of the stronger
    verdict 'unverifiable'. The reason names the rule that decided the verdict; a rule it
    outranks stays stated in the notes."""

    def test_no_bins_without_index_fp(self):
        g = parse(T.doc_text('no_bins'))
        self.assertIsNone(g.index_fp)
        for mode in MODES:
            c = g.compare(g, mode=mode)
            self.assertEqual(('unknown', None, 0), (c.comparable, c.equal, c.depth_used))
            self.assertIn(DEPTH_RULE, c.reason, mode)
            self.assertNotIn('index_fp', c.reason, mode)
            # the identity is still unverifiable, and still said so
            self.assertTrue(any('index_fp' in n for n in c.notes), (mode, c.notes))

    def test_no_bins_with_index_fp(self):
        g = parse(_with_fp(T.doc_text('no_bins')))
        for mode in MODES:
            c = g.compare(g, mode=mode)
            self.assertEqual('unknown', c.comparable)
            self.assertIn(DEPTH_RULE, c.reason, mode)

    def test_the_depth_rule_decides_before_missing_bases(self):
        # with nothing certified on a side, bases would not make the pair comparable
        from test_traverse_regressions import without_bases
        g = parse(_with_fp(T.doc_text('no_bins')))
        bare = without_bases(g)
        for mode in ('claims', 'walks', 'prefix_subset'):
            c = g.compare(bare, mode=mode)
            self.assertEqual('unknown', c.comparable)
            self.assertIn(DEPTH_RULE, c.reason, mode)
            self.assertTrue(any('output.sequences' in n for n in c.notes), (mode, c.notes))

    def test_the_tool_answers_the_same_reason(self):
        with tempfile.TemporaryDirectory() as d:
            store = GraphletStore(d)
            h = store.put_graphlet(T.graphlet('no_bins'))
            tools = GraphletTools(store)
            for mode in MODES:
                out = tools.graphlet_compare(h, h, mode=mode)
                self.assertEqual('unknown', out['comparable'], out)
                self.assertIn(DEPTH_RULE, out['reason'], out)

    def test_the_other_verdicts_keep_or_name_their_rule(self):
        # unverifiable with a common depth: the index_fp rule decides
        merge = T.graphlet('merge')
        c = merge.compare(merge)
        self.assertEqual('unverifiable', c.comparable)
        self.assertIn('index_fp', c.reason)
        # no bases with a common depth: that rule decides
        from test_traverse_regressions import without_bases
        c = merge.compare(without_bases(merge), mode='claims')
        self.assertEqual('unknown', c.comparable)
        self.assertIn('output.sequences', c.reason)
        # qualified by a lower-bound side: the reason names the qualification, not the
        # identity that held
        constrain = parse(_with_fp(T.doc_text('merge')))
        annotate = parse(_with_fp(T.doc_text('annotate')))     # lists cut at 2
        c = constrain.compare(annotate, arm='right', mode='claims', labels=['a.fa'])
        self.assertEqual('qualified', c.comparable)
        self.assertIn('lower bound', c.reason)
        self.assertIn(c.reason, c.notes)
        # comparable: the identity and the seed
        c = constrain.compare(constrain)
        self.assertEqual((True, 'same index and seed'), (c.comparable, c.reason))


class TestG2ATailOfZeroIsEmptyOnBothArms(unittest.TestCase):
    """G2 (VO8): graphlet_walks(spell='tail', tail_bp=0) answered the whole flank on the
    right arm (s[-0:]) and nothing on the left. DESIGN §6 takes tail_bp (default 60, no
    lower bound but the library's >= 0): the tail is the walk's tail_bp outermost bases
    (the last on the right arm, the first on the left), so 0 is an empty tail on both."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.store = GraphletStore(self.tmp.name)
        self.tools = GraphletTools(self.store)

    def tearDown(self):
        self.tmp.cleanup()

    def rows(self, h, side, **kw):
        out = self.tools.graphlet_walks(h, side, rank='id', n=100, max_bytes=1 << 20, **kw)
        self.assertNotIn('error', out, out)
        return out['rows']

    def test_tails_on_both_arms(self):
        for name in ('fork', 'linear', 'merge'):
            g = T.graphlet(name)
            h = self.store.put_graphlet(g)
            for side in g.arms:
                full = {r['walk']: r['sequence'] for r in self.rows(h, side, spell='full')}
                self.assertTrue(all(full.values()), (name, side))
                for n in (0, 1, 4, 10 ** 6):
                    got = {r['walk']: r['sequence']
                           for r in self.rows(h, side, spell='tail', tail_bp=n)}
                    want = {w: (s[len(s) - min(n, len(s)):] if side == 'right'
                                else s[:n]) for w, s in full.items()}
                    self.assertEqual(want, got, (name, side, n))
                    if n == 0:
                        self.assertEqual({''}, set(got.values()), (name, side))

    def test_a_negative_tail_is_refused(self):
        h = self.store.put_graphlet(T.graphlet('fork'))
        out = self.tools.graphlet_walks(h, 'right', spell='tail', tail_bp=-1)
        self.assertEqual('bad_argument', out.get('error'), out)


if __name__ == '__main__':
    unittest.main()
