"""The library defects the v8 test suite found (FORMAT-v8 §G), each reproduced offline on a
fixture of data/traverse/documents. One class per finding, named by its id there."""

import gzip
import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402

from metagraph.traverse import GraphletStore, derive, from_response, parse  # noqa: E402
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


def _real(path):
    """A committed real retrieval (data/traverse/<path>.graphlet.json.gz) as a client sees it."""
    with gzip.open(os.path.join(T.DATA, path + '.graphlet.json.gz')) as f:
        resp = json.loads(f.read())
    return from_response(resp['results'][0], resp)


def _by_position(g, a, sid):
    """The label ends of segment |sid| as (at_bp, ref) in position order, ties by label id."""
    runs = sorted((a.runs[r] for r in derive.label_end_events(a).get(sid, ())),
                  key=lambda r: (r.to_bp, r.label))
    return [(r.to_bp, g.labels[r.label].ref) for r in runs]


class TestG3AWalkRowsEndsInPositionOrder(unittest.TestCase):
    """G3 (VO9): a walk row's ends came in run order (linear's left root: at_bp 40, 40, 1).
    SPEC §7.1 sorts events by position ("events at one at_bp are unordered: the list is
    sorted by position only") and DESIGN §2.5 has the walker sort them by at_bp. The
    ends are ordered by at_bp; equal positions by label id (the dictionary order, in
    which end_labels are listed, DESIGN §2.5), so which 8 of more a row shows follows the
    positions and never the run numbering, which at one position is not information."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.store = GraphletStore(self.tmp.name)
        self.tools = GraphletTools(self.store)

    def tearDown(self):
        self.tmp.cleanup()

    def walk_rows(self, h, side, pid):
        out = self.tools.graphlet_walk(h, side, pid, max_bytes=1 << 22)
        self.assertNotIn('error', out, out)
        self.assertNotIn('next_cursor', out)
        return out['rows']

    def test_linear_root_ends(self):
        g = T.graphlet('linear')
        h = self.store.put_graphlet(g)
        want = {'left': [(1, 'h:0:2'), (40, 'h:0:0'), (40, 'h:0:1')],
                'right': [(0, 'h:0:1'), (40, 'h:0:0'), (40, 'h:0:2')]}
        for side in ('left', 'right'):
            row = self.walk_rows(h, side, 0)[0]
            self.assertEqual(0, row['segment'])
            self.assertEqual(want[side], [(e['at_bp'], e['ref']) for e in row['ends']], side)

    def test_the_rows_shown_and_exported_follow_the_positions(self):
        # segments with more than 8 label ends whose first 8 in run order are not the
        # first 8 by position
        cases = (('real/sra/sra_rand50_00__limit1', 'right', 8),
                 ('real/sra/sra_rand50_00__limit1', 'left', 12),
                 ('coords/rep0_kmer_null', 'right', 0),
                 ('coords/rep0_cut1', 'right', 7))
        for path, side, sid in cases:
            g = _real(path)
            a = g.arms[side]
            want = _by_position(g, a, sid)
            runs = [a.runs[r] for r in derive.label_end_events(a)[sid]]
            self.assertGreater(len(want), 8, path)
            self.assertNotEqual(set(want[:8]),
                                {(r.to_bp, g.labels[r.label].ref) for r in runs[:8]}, path)
            pid = next(i for i, p in enumerate(derive.paths(a)) if sid in p.segments)
            h = self.store.put_graphlet(g)
            row, = [r for r in self.walk_rows(h, side, pid) if r['segment'] == sid]
            self.assertEqual(want[:8], [(e['at_bp'], e['ref']) for e in row['ends']], path)
            # the export of the walk lists them all, in the same order
            res = self.tools.graphlet_export(h, format='json', what='walk', arm=side,
                                             walks=[pid], path='w.json')
            self.assertNotIn('error', res, res)
            with open(res['path'], encoding='utf-8') as f:
                doc = json.load(f)
            seg, = [x for x in doc['walks'][0]['segments'] if x['segment'] == sid]
            self.assertEqual(want, [(e['at_bp'], e['ref']) for e in seg['ends']], path)


if __name__ == '__main__':
    unittest.main()
