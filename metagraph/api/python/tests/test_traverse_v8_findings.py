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


if __name__ == '__main__':
    unittest.main()
