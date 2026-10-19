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

from metagraph.traverse import (  # noqa: E402
    GraphletStore, derive, from_response, ops, parse,
)
from metagraph.traverse.mcp_tools import GraphletTools, ToolLimits  # noqa: E402

MODES = ('claims', 'walks', 'labels', 'prefix_subset')
R = ops.COMPARE_REASONS
DEPTH = R['depth']


def _with_fp(text, fp='cd' * 32):
    """A CLI fixture with an index manifest digest (they are made without one)."""
    return text.replace(' walk fixtures * ', ' walk fixtures %s ' % fp, 1)


def _coords(name):
    """A committed snapshot of data/traverse/coords as a client sees it."""
    with gzip.open(os.path.join(T.DATA, 'coords', name + '.graphlet.json.gz')) as f:
        resp = json.loads(f.read())
    return from_response(resp['results'][0], resp)


def _level5(name):
    """The fixture as a level-5 server (arrival-order parents) answered it."""
    with open(os.path.join(T.DATA, 'r21', name + '_level5.graphlet.json'),
              encoding='utf-8') as f:
        resp = json.load(f)
    return from_response(resp['results'][0], resp)


class TestG1TheReasonNamesTheRuleThatDecided(unittest.TestCase):
    """G1 (O1): no_bins (no index_fp, complete_to_bp 0 on both arms) compares 'unknown' by
    the depth rule, but its reason named the missing index_fp, the rule of the stronger
    verdict 'unverifiable'. The reason is the fixed phrase (ops.COMPARE_REASONS) of the
    rule that decided the verdict; the notes state every other rule that applied, and the
    deciding one where its phrase says less, each once."""

    def test_no_bins_without_index_fp(self):
        g = parse(T.doc_text('no_bins'))
        self.assertIsNone(g.index_fp)
        for mode in MODES:
            c = g.compare(g, mode=mode)
            self.assertEqual(('unknown', None, 0), (c.comparable, c.equal, c.depth_used))
            self.assertEqual(DEPTH, c.reason, mode)
            # the identity is still unverifiable, said once in a short note; the depth
            # rule's phrase is its whole statement, not repeated in the notes
            self.assertIn('no index_fp on both', c.notes, mode)
            self.assertNotIn(DEPTH, c.notes, mode)

    def test_no_bins_with_index_fp(self):
        g = parse(_with_fp(T.doc_text('no_bins')))
        for mode in MODES:
            c = g.compare(g, mode=mode)
            self.assertEqual(('unknown', DEPTH), (c.comparable, c.reason), mode)
            self.assertNotIn(DEPTH, c.notes, mode)

    def test_the_depth_rule_decides_before_missing_bases(self):
        # with nothing certified on a side, bases would not make the pair comparable
        from test_traverse_regressions import without_bases
        g = parse(_with_fp(T.doc_text('no_bins')))
        bare = without_bases(g)
        for mode in ('claims', 'walks', 'prefix_subset'):
            c = g.compare(bare, mode=mode)
            self.assertEqual(('unknown', DEPTH), (c.comparable, c.reason), mode)
            self.assertTrue(any('output.sequences' in n for n in c.notes), (mode, c.notes))

    def test_the_tool_answers_the_same_reason(self):
        with tempfile.TemporaryDirectory() as d:
            store = GraphletStore(d)
            h = store.put_graphlet(T.graphlet('no_bins'))
            tools = GraphletTools(store)
            for mode in MODES:
                out = tools.graphlet_compare(h, h, mode=mode)
                self.assertEqual(('unknown', DEPTH), (out['comparable'], out['reason']), out)

    def test_the_other_verdicts_name_their_rule(self):
        # unverifiable with a common depth: the index_fp rule decides; its statement is
        # its phrase, so no note repeats it
        merge = T.graphlet('merge')
        c = merge.compare(merge)
        self.assertEqual(('unverifiable', R['unverifiable_both']), (c.comparable, c.reason))
        self.assertNotIn(c.reason, c.notes)
        # no bases with a common depth: that rule decides, its whole text in the notes
        from test_traverse_regressions import without_bases
        c = merge.compare(without_bases(merge), mode='claims')
        self.assertEqual(('unknown', R['bases']), (c.comparable, c.reason))
        self.assertTrue(any(n.startswith('b: no bases; mode claims') for n in c.notes),
                        c.notes)
        # qualified by a lower-bound side: the reason names the qualification, its whole
        # text in the notes only
        constrain = parse(_with_fp(T.doc_text('merge')))
        annotate = parse(_with_fp(T.doc_text('annotate')))     # lists cut at 2
        c = constrain.compare(annotate, arm='right', mode='claims', labels=['a.fa'])
        self.assertEqual(('qualified', R['lower_bound']), (c.comparable, c.reason))
        self.assertNotIn(c.reason, c.notes)
        self.assertIn('b: label_lists (labels.max_labels_per_node): the oracle is a lower '
                      'bound', c.notes)
        # comparable: the identity and the seed
        c = constrain.compare(constrain)
        self.assertEqual((True, 'same index and seed'), (c.comparable, c.reason))
        self.assertNotIn(c.reason, c.notes)


class TestVF1AQualifiedAnswerFitsAsBefore(unittest.TestCase):
    """VF1 (VERIFY-F): G1 first made a qualified answer's reason its whole first
    qualification (up to ~480 bytes) in compare's base, which the tool never cuts, so
    answers that fit before failed and default pages lost rows or their evidence. Every
    reason is a fixed phrase of at most 80 bytes, and a qualification's at most 19, the
    length of 'same index and seed', the reason a qualified answer had before."""

    # the pages at the default max_bytes (2048) of the qualified pairs before G1
    # (bf010a984; VERIFY-F probe_tool_sizes.py): (rows, next_cursor, fields_cut)
    BEFORE = {
        ('rep0_column', 'rep0_cut1', 'claims'): (1, True, None),
        ('rep0_column', 'rep0_cut1', 'walks'): (2, True, None),
        ('rep0_cut1', 'rep0_kmer_null', 'claims'): (2, True, None),
        ('rep0_cut1', 'rep0_kmer_null', 'walks'): (2, True, None),
        ('rep0_lower_bound', 'rep0_mixed', 'claims'): (1, True, None),
        ('rep0_lower_bound', 'rep0_mixed', 'walks'): (1, True, None),
        ('merge', 'merge', 'claims'): (2, True, None),
        ('merge', 'merge', 'walks'): (2, False, None),
        ('annotate', 'annotate', 'claims'): (4, False, ['evidence']),
        ('annotate', 'annotate', 'walks'): (2, False, ['evidence']),
    }

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.store = GraphletStore(self.tmp.name)
        self.tools = GraphletTools(self.store)

    def tearDown(self):
        self.tmp.cleanup()

    def pair(self, x, y):
        if x == y:
            # level 5 (arrival order) against level 6 (majority parent), both with a
            # manifest digest: qualified by the displayed parents
            a, b = _level5(x), T.graphlet(x)
            for g in (a, b):
                g.index_fp = 'f' * 64
        else:
            a, b = _coords(x), _coords(y)
        return self.store.put_graphlet(a), self.store.put_graphlet(b)

    def test_the_smallest_page_of_rep0_cut1_and_rep0_kmer_null(self):
        ha, hb = self.pair('rep0_cut1', 'rep0_kmer_null')
        out = self.tools.graphlet_compare(ha, hb, mode='walks', max_bytes=950)
        self.assertNotIn('error', out, out)
        self.assertEqual(('qualified', R['scopes'], 1),
                         (out['comparable'], out['reason'], len(out['rows'])))

    def test_qualified_pages_at_the_default_as_before(self):
        for (x, y, mode), want in sorted(self.BEFORE.items()):
            ha, hb = self.pair(x, y)
            for mb in (None, 2048):
                out = self.tools.graphlet_compare(ha, hb, mode=mode, max_bytes=mb)
                self.assertEqual('qualified', out.get('comparable'), (x, y, mode, out))
                self.assertEqual(want, (len(out['rows']), 'next_cursor' in out,
                                        out.get('fields_cut')), (x, y, mode, mb))
                self.assertNotIn('row_truncated', out, (x, y, mode))

    def test_the_phrases_are_bounded(self):
        for k, v in R.items():
            self.assertLessEqual(len(v.encode('utf-8')), 80, k)
        for k in ('lower_bound', 'overstated', 'scopes', 'parents', 'inexact'):
            self.assertLessEqual(len(R[k].encode('utf-8')), len(R['same']), k)


class TestVF4EveryRuleAnswersItsPhrase(unittest.TestCase):
    """VF4 (VERIFY-F): each qualification answers its own phrase with its whole text in
    the notes; missing bases without index_fp state both rules once; every answer over
    the committed pairs keeps its reason a bounded fixed phrase that no note repeats."""

    def assert_once(self, c):
        self.assertNotIn(c.reason, c.notes)
        self.assertEqual(len(c.notes), len(set(c.notes)), c.notes)

    def test_each_qualification(self):
        cases = (
            # (a, b, mode, phrase, a note fragment)
            (_coords('rep0_column'), _coords('rep0_cut1'), 'claims', R['overstated'],
             'a: trace_record_boundaries (labels.seed_label_kind): reported evidence may '
             'be overstated'),
            (_coords('rep0_kmer_null'), _coords('rep0_lower_bound'), 'claims',
             R['lower_bound'], 'b: greedy_losses (branching.max_label_branches): the '
             'oracle is a lower bound'),
            (_coords('rep0_cut1'), _coords('rep0_kmer_null'), 'walks', R['scopes'],
             'left arm: completeness scopes differ (per_path vs united_history)'),
        )
        for a, b, mode, phrase, note in cases:
            c = a.compare(b, mode=mode)
            self.assertEqual(('qualified', phrase), (c.comparable, c.reason))
            self.assertTrue(any(n.startswith(note) for n in c.notes), (phrase, c.notes))
            self.assert_once(c)

    def test_the_displayed_parents(self):
        old, new = _level5('merge'), T.graphlet('merge')
        for g in (old, new):
            g.index_fp = 'f' * 64
        c = old.compare(new, mode='claims')
        self.assertEqual(('qualified', R['parents']), (c.comparable, c.reason))
        self.assertTrue(any(n.startswith('the displayed parent at a merge may follow '
                                         'different rules') for n in c.notes), c.notes)
        self.assert_once(c)

    def test_label_evidence_by_the_outcome_alone(self):
        # R2-3: the outcome's lower_bound understates, its qualified may overstate (SPEC
        # §7.0); the parser admits no third value besides complete
        for value, phrase in (('lower_bound', R['lower_bound']),
                              ('qualified', R['overstated'])):
            a, b = T.graphlet('merge'), T.graphlet('merge')
            for g in (a, b):
                g.index_fp = 'f' * 64
            b.limitations = []
            b.outcome.label_evidence = value
            c = a.compare(b, mode='labels')
            self.assertEqual(('qualified', phrase), (c.comparable, c.reason), value)
            self.assertIn('b: label_evidence %s' % value, c.notes)
            self.assert_once(c)
        self.assertNotIn('incomplete', R)

    def test_an_inexact_restriction(self):
        # test_traverse_review3's construction: a merge partition without b.fa
        from test_traverse_review3 import truncated
        g = T.graphlet('merge')
        g.index_fp = 'f' * 64
        t = truncated(g, 50)
        a = g.arms['right']
        m = next(x for x in a.segments if len(x.parents) > 1 and x.from_bp == 62)
        for part in m.partition:
            for i in range(len(part) - 1, -1, -1):
                if g.labels[part[i]].name == 'b.fa':
                    del part[i]
        g.cache.clear()
        a.cache.clear()
        c = t.compare(g, mode='claims')
        self.assertEqual(('qualified', R['inexact']), (c.comparable, c.reason))
        self.assertTrue(any('could not be reconstructed exactly' in n for n in c.notes))
        self.assert_once(c)

    def test_no_bases_without_index_fp(self):
        from test_traverse_regressions import without_bases
        merge = T.graphlet('merge')
        self.assertIsNone(merge.index_fp)
        for mode in ('claims', 'walks', 'prefix_subset'):
            c = merge.compare(without_bases(merge), mode=mode)
            self.assertEqual(('unknown', R['bases']), (c.comparable, c.reason), mode)
            # both rules stated once: the identity unverifiable, then which side has no
            # bases and what to do
            self.assertEqual(1, c.notes.count('no index_fp on both'), c.notes)
            self.assertEqual(1, sum(n.startswith('b: no bases; mode %s' % mode)
                                    for n in c.notes), c.notes)
            self.assert_once(c)

    def test_every_committed_pair(self):
        import traverse_golden as TG
        fixed = set(R.values())
        items = TG.items_of(TG.golden_files() + TG.coordinate_files())
        by = {}
        for key, r, resp in items:
            by.setdefault(from_response(r, resp).seed.sequence, []).append((r, resp))
        n = 0
        for xs in by.values():
            for (ra, pa), (rb, pb) in zip(xs, xs[1:]):
                for mode in MODES:
                    a, b = from_response(ra, pa), from_response(rb, pb)
                    for x, y in ((a, b), (b, a)):
                        c = x.compare(y, mode=mode)
                        n += 1
                        self.assertIn(c.reason, fixed)
                        self.assertLessEqual(len(c.reason.encode('utf-8')), 80)
                        if c.comparable == 'qualified':
                            self.assertLessEqual(len(c.reason), len(R['same']))
                        self.assert_once(c)
        self.assertGreater(n, 100)


class TestR21NoAnswerLongerThanBeforeG1(unittest.TestCase):
    """R2-1 (VERIFY-F-R2): no compare answer may be longer than before G1. Over the
    verifier's probe set (traverse_golden.compare_page_shapes: 336 pairs, 4 modes, max_bytes
    default, 950, 1024, 2048 and 4 MiB), against the pages of bf010a984 recorded in
    golden/compare_shapes_bf010a984.json.gz (traverse_golden.py --compare-shapes --lib
    <bf010a984>): every whole answer's bytes without its rows are at most what they were,
    and no page is worse -- no new error, no fewer rows, no row newly truncated, no more
    fields cut unless that makes a truncated row whole."""

    def test_the_probe_set(self):
        import traverse_golden as TG
        import metagraph.traverse as MT
        want = TG.read_json(os.path.join(TG.GOLDEN_DIR, 'compare_shapes_bf010a984.json.gz'))
        got = TG.compare_page_shapes(MT)
        self.assertEqual(sorted(want), sorted(got))
        whole = str(TG.SHAPE_SIZES[-1])
        longer, worse = [], []
        for k in sorted(want):
            (e0, r0, _, c0, t0, b0), (e1, r1, _, c1, t1, b1) = want[k], got[k]
            if k.endswith(':' + whole):
                if not e0 and not e1 and b1 > b0:
                    longer.append((k, b0, b1))
                continue
            if e0 is None and e1:
                worse.append((k, 'error', e1))
            elif e0 is None:
                if r1 < r0 or (t1 and not t0) or (
                        len(c1 or ()) > len(c0 or ()) and not (t0 and not t1)):
                    worse.append((k, want[k], got[k]))
        self.assertEqual([], longer[:10], '%d answers longer' % len(longer))
        self.assertEqual([], worse[:10], '%d pages worse' % len(worse))


class TestR22TheArmRuleAnswersItsPhrase(unittest.TestCase):
    """R2-2 (VERIFY-F-R2): the arm rule ('no arm in common', an explicit arm b did not
    retrieve) answers fixed phrases of COMPARE_REASONS, and an unverifiable identity it
    decides over stays said in the short note."""

    def test_an_arm_b_did_not_retrieve(self):
        # fork has both arms
        for fp in (None, 'f' * 64):
            a, b = T.graphlet('fork'), T.graphlet('fork')
            for g in (a, b):
                g.index_fp = fp
            del b.arms['right']
            b.cache.clear()
            c = a.compare(b, arm='right')
            self.assertEqual((False, R['right']), (c.comparable, c.reason))
            self.assertEqual([] if fp else ['no index_fp on both'], c.notes)
            # without an arm the left one is compared
            self.assertNotEqual(False, a.compare(b).comparable)

    def test_no_arm_in_common(self):
        a, b = T.graphlet('fork'), T.graphlet('fork')
        del a.arms['right']
        del b.arms['left']
        for g in (a, b):
            g.cache.clear()
        c = a.compare(b)
        self.assertEqual((False, R['no_arm'], ['no index_fp on both']),
                         (c.comparable, c.reason, c.notes))


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

    def test_guard_a_negative_tail_stays_refused(self):
        # a guard, not the defect: _int refused it before the fix as well
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


class TestVF2TheEndsSortIsCharged(unittest.TestCase):
    """VF2 (VERIFY-F): G3 sorts a walk row's ends by position, and a budgeted walk did not
    charge the sort (DESIGN §5.4: work units count what an algorithm visits, admitted
    before the work). Each row now charges sort_work(ends) with it: sra_hairpin's left
    root has 977 ends, whose sort (17,586 units) the rest of the walk's charge does not
    cover."""

    def test_a_budgeted_walk_charges_its_sorts(self):
        from metagraph.traverse import budget as B
        g = _real('real/sra/sra_hairpin__left_only')
        a = g.arms['left']
        p = derive.paths(a)[0]
        ends = derive.label_end_events(a)
        sorts = sum(B.sort_work(len(ends.get(sid, ()))) for sid in p.segments
                    if len(ends.get(sid, ())) > 1)
        self.assertGreater(max(len(ends.get(sid, ())) for sid in p.segments), 900)
        b = B.LocalBudget()
        rows = GraphletTools._walk_rows(g, a, p, out=[], b=b)
        self.assertGreaterEqual(b.used_work, sorts)
        # the charge changes nothing in the rows
        self.assertEqual(GraphletTools._walk_rows(g, a, p), rows)


class TestVR4ASliceCarriesABaseOrIsRefused(unittest.TestCase):
    """VR4 (VERIFY-L-R2): under local limits with a small sequence ceiling (900 bytes),
    graphlet_sequence answered an empty slice whose next_from was its from, so paging it
    never advanced, instead of result_too_large. A result that cannot be made to fit is
    result_too_large (mcp_tools' contract, DESIGN §6): a cut slice carries at least one
    base, and where not one fits beside the frame the whole slice is answered when it
    fits the ceiling, else refused."""

    def tools(self, d, ceiling, limited):
        return GraphletTools(GraphletStore(d), {}, sequence_max_bytes=ceiling,
                             local_limits=ToolLimits() if limited else None)

    def test_rep0_cut1_at_900_under_local_limits(self):
        with tempfile.TemporaryDirectory() as d:
            t = self.tools(d, 900, True)
            h = t.store.put_graphlet(_coords('rep0_cut1'))
            for side in ('left', 'right'):
                out = t.graphlet_sequence(h, side, walk=0)
                # the 300 bases fit beside the local block as answered: the whole slice
                self.assertNotIn('error', out, out)
                self.assertEqual((0, 300, 300), (out['from'], out['to'], len(out['sequence'])))
                self.assertNotIn('next_from', out)
            t = self.tools(d, 700, True)
            out = t.graphlet_sequence(h, 'left', walk=0)
            self.assertEqual('result_too_large', out.get('error'), out)

    def test_every_slice_advances_or_is_refused(self):
        n = refused = 0
        for name in ('rep0_cut1', 'rep0_kmer_null', 'rep0_lower_bound'):
            g = _coords(name)
            with tempfile.TemporaryDirectory() as d:
                for limited in (False, True):
                    for ceiling in (450, 600, 700, 900, 1200):
                        t = self.tools(d, ceiling, limited)
                        h = t.store.put_graphlet(g)
                        for side in g.arms:
                            whole = g.spell(side, 0)
                            got, start = '', 0
                            while True:
                                out = t.graphlet_sequence(h, side, walk=0, start=start)
                                n += 1
                                if out.get('error'):
                                    self.assertEqual('result_too_large', out['error'], out)
                                    refused += 1
                                    break
                                self.assertEqual(start, out['from'])
                                got += out['sequence']
                                if 'next_from' not in out:
                                    self.assertEqual(whole, got, (name, side, ceiling))
                                    break
                                self.assertGreater(out['next_from'], start, out)
                                start = out['next_from']
        self.assertGreater(refused, 0)
        self.assertGreater(n, refused)


def _formula_profile(g, side, pid):
    """DESIGN §5.1's displayed evidence along walk |pid|, per base, from the library's own
    claims (one per run, with its evidence_from): a run's label supports base x where
    evidence_from <= x < to_bp and the walk's segment at x is on the run's displayed chain
    (the first-parent chain to its anchor)."""
    a = g.arms[side]
    p = derive.paths(a)[pid]
    seg_at = []
    for sid in p.segments:
        seg_at += [sid] * a.segments[sid].length_bp
    sets = [set() for _ in seg_at]
    for c in g.claims(side, strict=False):
        if c.to_bp <= c.from_bp:
            continue
        mine = set(derive.chain(a, c.segment))
        for x in range(c.evidence_from, min(c.to_bp, len(seg_at))):
            if seg_at[x] in mine:
                sets[x].add(c.label.ref)
    return sets


def _profile_sets(g, side, pid, kind='displayed'):
    out = []
    for r in g.support_profile(side, pid, kind):
        out += [{l.ref for l in r.labels}] * (r.to_bp - r.from_bp)
    return out


class TestG4TheProfileIsDisplayedEvidence(unittest.TestCase):
    """G4 (VERIFY-O-B VO-B1): graphlet_support's profile (support_profile(kind=
    'displayed')) counted a label wherever it was present on the walk's nodes. DESIGN §5.1
    lists support_profile(kind='displayed') among the users of displayed-path evidence:
    a run's displayed bases start at its evidence_from ("the displayed bases before s
    are not this lineage's"), derived at its anchor. The profile now counts each run
    from its evidence_from, on the part of its displayed chain the walk shares."""

    def test_mini_ndm1_rc_switch1_left_walk_3(self):
        g = _real('real/mini_refseq/mini_ndm1_rc__switch1')
        a = g.arms['left']
        run = a.runs[13]
        claim, = [c for c in g.claims('left', strict=False)
                  if c.segment == run.segment and c.label.id == run.label]
        self.assertEqual((1239, 3), (claim.evidence_from, claim.path_id))
        prof = g.support_profile('left', 3)
        ref = g.labels[run.label].ref
        for r in prof:
            has = ref in {l.ref for l in r.labels}
            self.assertEqual(r.from_bp >= 1239, has, (r.from_bp, r.to_bp))
        self.assertEqual(5, len(prof[0].labels))
        # the tool answers the same rows
        with tempfile.TemporaryDirectory() as d:
            store = GraphletStore(d)
            h = store.put_graphlet(g)
            rows = GraphletTools(store).graphlet_support(h, 'left', 3, max_bytes=1 << 20)['rows']
        self.assertEqual([(r.from_bp, r.to_bp, len(r.labels)) for r in prof],
                         [(r['from_bp'], r['to_bp'], r['n_labels']) for r in rows])

    def test_merge_document(self):
        # a.fa and b.fa enter walk 0 through non-first parents (evidence_from 62 and 36):
        # they join there, and never leave at the bubbles they do not display
        g = T.body('merge')
        self.assertEqual([(36, ['b.fa'], [], 'merge'), (62, ['a.fa'], [], 'merge')],
                         [(c.at_bp, [l.name for l in c.added], [l.name for l in c.removed],
                           c.reasons[0]['why']) for c in g.support_changes('right', 0)])
        sets = _profile_sets(g, 'right', 0)
        names = {l.ref: l.name for l in g.labels}
        for x, labs in enumerate(sets):
            got = {names[r] for r in labs}
            self.assertEqual('a.fa' in got, x >= 62, x)
            self.assertEqual('b.fa' in got, x >= 36, x)

    def test_every_constrain_walk_of_the_fixtures(self):
        import glob
        n = 0
        docs = [(os.path.basename(f), parse(T.read(f)))
                for f in sorted(glob.glob(os.path.join(T.DOCS, '*.mgt')))]
        for f in sorted(glob.glob(os.path.join(T.DATA, 'real', '*', '*.graphlet.json.gz'))
                        + glob.glob(os.path.join(T.DATA, 'coords', '*.graphlet.json.gz'))
                        + glob.glob(os.path.join(T.DATA, 'review3', '*.graphlet.json'))):
            if f.endswith('.gz'):
                with gzip.open(f) as fh:
                    resp = json.loads(fh.read())
            else:
                with open(f, encoding='utf-8') as fh:
                    resp = json.load(fh)
            resp = resp['response'] if isinstance(resp.get('response'), dict) else resp
            for r in resp['results']:
                if isinstance(r.get('graphlet'), str):
                    docs.append((os.path.basename(f), from_response(r, resp)))
        for name, g in docs:
            if g.mode != 'constrain':
                continue
            for side, a in g.arms.items():
                for p in derive.paths(a):
                    if p.length_bp == 0:
                        continue
                    n += 1
                    self.assertEqual(_formula_profile(g, side, p.id),
                                     _profile_sets(g, side, p.id), (name, side, p.id))
                    # support_changes reads the same profile: past base 0 a change is
                    # exactly where its label set changes, and says what changed
                    prof = g.support_profile(side, p.id)
                    want = [(r.from_bp, {l.ref for l in r.labels} - {l.ref for l in q.labels},
                             {l.ref for l in q.labels} - {l.ref for l in r.labels})
                            for q, r in zip(prof, prof[1:])]
                    got = [(c.at_bp, {l.ref for l in c.added}, {l.ref for l in c.removed})
                           for c in g.support_changes(side, p.id) if c.at_bp > 0]
                    self.assertEqual(want, got, (name, side, p.id))
        self.assertGreater(n, 400)


def _constrain_fixtures():
    """(name, graphlet) of every constrain fixture: the documents and the real, coords and
    review3 retrievals."""
    import glob
    out = [(os.path.basename(f), parse(T.read(f)))
           for f in sorted(glob.glob(os.path.join(T.DOCS, '*.mgt')))]
    for f in sorted(glob.glob(os.path.join(T.DATA, 'real', '*', '*.graphlet.json.gz'))
                    + glob.glob(os.path.join(T.DATA, 'coords', '*.graphlet.json.gz'))
                    + glob.glob(os.path.join(T.DATA, 'review3', '*.graphlet.json'))):
        if f.endswith('.gz'):
            with gzip.open(f) as fh:
                resp = json.loads(fh.read())
        else:
            with open(f, encoding='utf-8') as fh:
                resp = json.load(fh)
        resp = resp['response'] if isinstance(resp.get('response'), dict) else resp
        for i, r in enumerate(resp['results']):
            if isinstance(r.get('graphlet'), str):
                out.append(('%s#%d' % (os.path.basename(f), i), from_response(r, resp)))
    return [(n, g) for n, g in out if g.mode == 'constrain']


class TestG41TheChangesStartFromTheDisplayedProfile(unittest.TestCase):
    """G4-1 (VERIFY-F-G4): the support-change rows started from the seed's entry set, so a
    seed label never displayed on the walk was reported removed at 0 with the fallback
    reason 'ended' (406 rows). The start set is the displayed profile at base 0: such a
    label was never there, from the walk's view; a seed label with a run ending at 0 keeps
    its K-coded end."""

    def test_no_fallback_end_on_any_constrain_walk(self):
        n = 0
        for name, g in _constrain_fixtures():
            for side, a in g.arms.items():
                for p in derive.paths(a):
                    for c in g.support_changes(side, p.id):
                        for r in c.reasons:
                            n += 1
                            self.assertFalse(r['why'] == 'ended' and 'end' not in r,
                                             (name, side, p.id, c.at_bp, r))
        self.assertGreater(n, 1000)

    def test_limits_right_walk_0(self):
        # c:2 enters the merge segment through its non-first parent (evidence_from 7): it
        # displays nothing on the walk, so no row names it (it was 'removed at 0, ended')
        g = T.body('limits')
        rows = [(c.at_bp, [l.ref for l in c.added], [l.ref for l in c.removed],
                 [r['why'] for r in c.reasons]) for c in g.support_changes('right', 0)]
        self.assertEqual([(5, [], ['c:1'], ['loss_budget'])], rows)

    def test_a_boundary_end_at_0_keeps_its_end(self):
        # linear: h:0:1 ends minority at the seed boundary (a zero-length run)
        lin = T.body('linear').support_changes('right', 0)
        self.assertEqual([(0, ['h:0:1'], 'minority')],
                         [(c.at_bp, [l.ref for l in c.removed], c.reasons[0]['why'])
                          for c in lin])


def _expected_prefixes(t):
    """(side, ref) -> the walking-order prefixes w[:j] (j > 0) over which a label has
    displayed evidence from base 0 on a walk of |t|, read off its support profile."""
    out = {}
    for side, a in t.arms.items():
        for p in derive.paths(a):
            prof = t.support_profile(side, p.id)
            if not prof:
                continue
            w = t.spell(side, p.id, 'walk')
            alive, reach = None, {}
            for r in prof:
                ids = {l.id for l in r.labels}
                alive = ids if alive is None else alive & ids
                for l in alive:
                    reach[l] = r.to_bp        # a label leaving keeps the prefix it reached
            for l, j in reach.items():
                out.setdefault((side, t.labels[l].ref), set()).add(w[:j])
    return out


def _library_prefixes(g, depth):
    got = ops._supported_prefixes(g, list(g.arms), depth)
    return {k: {w for w in v if w} for k, v in got.items() if any(v)}


class TestG5PrefixSubsetReadsDisplayedEvidence(unittest.TestCase):
    """G5 (VERIFY-F-G4, the open note): compare(mode='prefix_subset') keyed b's supported
    prefixes (_supported_prefixes) by presence contiguous from base 0, while its docstring
    promises displayed support. In mini_ndm1_rc__switch1 left, labels present on [0, 186)
    display only from 1239. b's map is now displayed evidence (DESIGN §5.1) in the DAG
    restricted to the depth: equal to the support profiles of a retrieval walked to the
    depth (test_traverse_review3.truncated)."""

    def check(self, g, depth, label):
        from test_traverse_review3 import truncated
        want = _expected_prefixes(truncated(g, depth))
        want = {k: v for k, v in want.items() if v}
        self.assertEqual(want, _library_prefixes(g, depth), label)

    def test_mini_ndm1_rc_switch1(self):
        g = _real('real/mini_refseq/mini_ndm1_rc__switch1')
        top = max(s.end_bp for a in g.arms.values() for s in a.segments)
        got = ops._supported_prefixes(g, ['left'], top)
        # label 1 is present on walk 3's first 186 bases, displayed only from 1239
        ref = g.labels[1].ref
        self.assertFalse(any(0 < len(w) <= 186 for w in got.get(('left', ref), ())), ref)
        self.check(g, top, 'full')

    def test_the_constrain_fixtures_at_their_depths(self):
        n = 0
        for name, g in _constrain_fixtures():
            top = max((s.end_bp for a in g.arms.values() for s in a.segments), default=0)
            if top <= 0 or any(not ops._has_bases(g, s) for s in g.arms):
                continue
            depths = set()
            for a in g.arms.values():
                depths.update(s.from_bp + 1 for s in a.segments if len(s.parents) > 1)
            for d in sorted(depths)[:5] + [top]:
                if d > 0:
                    n += 1
                    self.check(g, d, (name, d))
        self.assertGreater(n, 50)


if __name__ == '__main__':
    unittest.main()
