"""The format's derivations (§2.3, §2.5, §5.1) on the fixture documents: trie
structure, derived events, end labels, the chronological G.end rule, label_summary in
both modes, continuation spellings on both arms, the anchored evidence, check_rules.
Expected values come from the server's own detail: full result (full.json) where it has
them, so a derivation is checked against the walker, not against a transcription."""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402

from metagraph.traverse import derive, parse  # noqa: E402


class TestStructure(unittest.TestCase):
    def test_paths_follow_first_parents_in_leaf_order(self):
        g = T.body('merge')
        a = g.arms['right']
        # the merge at 62 stores its parents majority first (R21 (4), feature level 6): the
        # G allele (segment 5: b.fa, c.fa, both.fa) arrived second and is the first parent,
        # so the displayed path passes through it; merge_ties keeps the arrival order
        self.assertEqual([(), (0,), (0,), (1, 2), (3,), (3,), (5, 4)],
                         [s.parents for s in a.segments])
        self.assertEqual([[1, 2], [3], [3], [4, 5], [6], [6], []],
                         [s.children for s in a.segments])
        self.assertEqual([6], derive.leaves(a))
        (p,) = derive.paths(a)
        self.assertEqual(([0, 1, 3, 5, 6], 100), (p.segments, p.length_bp))
        (p,) = derive.paths(T.body('merge_ties').arms['right'])
        self.assertEqual([0, 1, 3, 4, 6], p.segments)
        for name in ('fork', 'merge', 'annotate', 'caps', 'ambiguous_split'):
            full = T.doc_json(name, 'full')['results'][0]
            g = T.body(name)
            for side, arm in full['arms'].items():
                self.assertEqual([(q['id'], q['segments'], q['length_bp']) for q in arm['paths']],
                                 [(q.id, q.segments, q.length_bp)
                                  for q in derive.paths(g.arms[side])], (name, side))

    def test_splits_group_children_and_keep_the_stored_kind(self):
        g = T.body('merge')
        sp = derive.splits(g.arms['right'], g.mode)
        # both.fa is on both alleles at each bubble
        self.assertEqual([(20, 0, [1, 2], True, 4), (46, 3, [4, 5], True, 4)],
                         [(s.at_bp, s.segment, s.children, s.ambiguous, s.labels_before)
                          for s in sp])
        # disjoint child sets, yet ambiguous (after a switch): stored, not derived
        g = T.body('ambiguous_split')
        right = g.arms['right']
        (s,) = derive.splits(right, g.mode)
        self.assertTrue(s.ambiguous)
        self.assertFalse(set(right.segments[1].entry) & set(right.segments[2].entry))
        # overlapping child sets, yet a divergence (a followed hairpin)
        left = g.arms['left']
        (s,) = derive.splits(left, g.mode)
        self.assertFalse(s.ambiguous)
        self.assertEqual(list(left.segments[1].entry), list(left.segments[2].entry))
        self.assertTrue(left.segments[0].events[0].followed)

    def test_annotate_split_counts_the_last_presence_total(self):
        g = T.body('annotate')
        s, _ = derive.splits(g.arms['right'], g.mode)
        self.assertEqual(4, s.labels_before)          # the root's last P total, not its cut list
        b = derive.split_branches(g.arms['right'], s, g.cap)
        # a.fa, c.fa, both.fa take A (3 labels, cut to 2 at the cap); b.fa, both.fa take T
        self.assertEqual([('A', 3, [0, 2]), ('T', 2, [1, 3])],
                         [(x['char'], x['labels_distinct'], x['labels']) for x in b])
        full = T.doc_json('annotate', 'full')['results'][0]['arms']['right']['splits']
        self.assertEqual([x['labels_before'] for x in full],
                         [x.labels_before for x in derive.splits(g.arms['right'], g.mode)])

    def test_split_chars_without_bases(self):
        g = T.body('limits')
        a = g.arms['right']
        self.assertIsNone(a.segments[1].walk)
        (s,) = derive.splits(a, g.mode)
        self.assertEqual(['C', 'G'], [b['char'] for b in derive.split_branches(a, s, g.cap)])


class TestEvents(unittest.TestCase):
    def test_label_end_events_come_from_ended_runs_only(self):
        g = T.body('switch_chain')
        a = g.arms['right']
        ev = derive.label_end_events(a)
        # the two Lw runs end silently on segment 0; the switches document them
        self.assertNotIn(0, ev)
        self.assertEqual({1: [2], 2: [3]}, ev)
        g = T.body('merge')
        ev = derive.label_end_events(g.arms['right'])
        # runs 3, 4 (both.fa on the minority parents) are closed by merges, not ended;
        # run 5 is both.fa's lineage through the first parents
        self.assertEqual({6: [0, 1, 2, 5]}, ev)

    def test_reconverge_events_from_multi_parent_segments(self):
        g = T.body('merge')
        self.assertEqual({3: (36, [1, 2]), 6: (62, [5, 4])},
                         derive.reconverge_events(g.arms['right']))
        self.assertEqual({3: (36, [1, 2]), 6: (62, [4, 5])},
                         derive.reconverge_events(T.body('merge_ties').arms['right']))

    def test_end_labels_ascending_with_leaf_extras(self):
        g = T.body('merge')
        ends = derive.end_labels(g.arms['right'], 6)
        # a.fa came in through the second parent of the last merge (62), b.fa through the
        # second parent of the first (36); c.fa is on both first parents; both.fa branched
        # at both bubbles and its lineage through the first parents is run 5
        self.assertEqual([(0, 0, 0, 62), (1, 1, 0, 36), (2, 2, 0, 0), (3, 5, 2, 0)],
                         [(e.label, e.run, e.branches, e.route_bp) for e in ends])
        self.assertEqual([2.0], [e.loss for e in derive.end_labels(
            T.body('switch_chain').arms['right'], 2)])
        for name in ('merge', 'merge_ties', 'switch_chain', 'caps', 'linear'):
            full = T.doc_json(name, 'full')['results'][0]
            g = T.body(name)
            for side, arm in full['arms'].items():
                a = g.arms[side]
                for q, leaf in zip(arm['paths'], derive.leaves(a)):
                    self.assertEqual(q['end_labels'],
                                     [{'label': e.label, 'loss': e.loss, 'branches': e.branches,
                                       'run': e.run, 'route_bp': e.route_bp}
                                      for e in derive.end_labels(a, leaf)], (name, side))

    def test_needed_budgets(self):
        self.assertEqual([2.0], derive.needed_budgets(T.body('limits').arms['right']))
        self.assertEqual([2.0], derive.needed_budgets(T.body('caps').arms['right']))
        self.assertEqual([], derive.needed_budgets(T.body('linear').arms['right']))


class TestEndRule(unittest.TestCase):
    def test_chronological_rule_keeps_a_reentered_label(self):
        g = T.body('reentry')
        a = g.arms['right']
        seg = a.segments[0]
        self.assertIn('end', seg.stars)
        self.assertEqual([0], list(seg.end))
        # v1's set formula: (entry U switch targets) - labels of runs ending inside
        targets = {e.to_label for e in seg.events if e.type == 'switch'}
        ended = {a.runs[r].label for r in derive.runs_by_segment(a)[0]
                 if a.runs[r].to_bp < seg.end_bp}
        self.assertEqual(set(), (set(seg.entry) | targets) - ended)

    def test_ends_before_switch_ins_at_one_position(self):
        g = T.body('switch_chain')
        a = g.arms['right']
        self.assertEqual([1], derive.rule_end(g.mode, a.segments[0], a.segments[0].entry,
                                              [], derive.runs_by_segment(a)[0], a.runs))

    def test_runs_ending_at_the_last_node_stay(self):
        g = T.body('linear')
        a = g.arms['right']
        self.assertEqual([0, 2], list(a.segments[0].end))   # acc2 ended at 0 < 40

    def test_annotate_end_is_the_last_presence_set(self):
        g = T.body('annotate')
        a = g.arms['right']
        for s in a.segments:
            self.assertIn('end', s.stars)
            self.assertEqual(list(s.presence[-1].labels), list(s.end))
        full = T.doc_json('annotate', 'full')['results'][0]['arms']['right']['segments']
        self.assertEqual([x['labels_at_end'] for x in full], [list(s.end) for s in a.segments])

    def test_entry_rules(self):
        g = T.body('merge')
        a = g.arms['right']
        self.assertEqual([0, 1, 2, 3], derive.rule_entry('constrain', 4, a.segments[0]))
        self.assertEqual([0, 1, 2, 3], derive.rule_entry('constrain', 4, a.segments[3]))
        self.assertIsNone(derive.rule_entry('constrain', 4, a.segments[1]))
        self.assertIsNone(derive.rule_entry('annotate', 0, a.segments[0]))
        self.assertEqual([[], []], [list(x) for x in derive.rule_partition((1, 2))])


class TestLabelSummary(unittest.TestCase):
    def test_constrain_lineage_roots_and_reentries(self):
        s = derive.label_summary(T.body('switch_chain'))
        # A's lineage reaches the radius (80) through B and then C or D
        self.assertEqual({'direct_bp': 30, 'reach_bp': 80, 'reentries': 0, 'runs': [0]},
                         s[0]['right'])
        self.assertEqual(1, s[1]['right']['reentries'])
        self.assertEqual(1, s[2]['right']['reentries'])   # one switch event into C
        self.assertEqual(0, s[2]['right']['reach_bp'])    # credited to the root A
        for name in ('switch_chain', 'reentry', 'merge', 'caps', 'annotate', 'same_name'):
            full = T.doc_json(name, 'full')['results'][0]['label_summary']
            mine = derive.label_summary(T.body(name))
            self.assertEqual([{side: x[side] for side in x if side != 'label'} for x in full],
                             [{side: dict(v, runs=list(v['runs'])) for side, v in x.items()}
                              for x in mine], name)

    def test_a_merge_does_not_clamp_direct_bp(self):
        s = derive.label_summary(T.body('merge'))
        self.assertEqual(100, s[1]['right']['direct_bp'])  # b.fa, through two second parents
        self.assertEqual([3, 4, 5], s[3]['right']['runs'])  # both.fa: kept and two closed runs

    def test_annotate_union_pass(self):
        s = derive.label_summary(T.body('annotate'))
        self.assertEqual((70, 70), (s[0]['right']['direct_bp'], s[0]['right']['reach_bp']))
        # c.fa was never in a recorded list at the boundary (cut at 2): reached, not direct
        self.assertEqual((0, 61), (s[2]['right']['direct_bp'], s[2]['right']['reach_bp']))
        self.assertEqual(40, s[1]['left']['direct_bp'])
        self.assertTrue(all(not x[side]['runs'] for x in s for side in x))


class TestContinuations(unittest.TestCase):
    """Right: the LAST n of seed + natural(right flank); left: the FIRST n of
    natural(left flank) + seed (v2) -- contained in the flank and crossing into it."""

    def check(self, name, side, expected=None):
        """expected None: the server's own continuation sequences (full.json)"""
        g = T.body(name)
        a = g.arms[side]
        got = [derive.continuation_sequence(g, a, p.leaf) for p in derive.paths(a)]
        if expected is None:
            arm = T.doc_json(name, 'full')['results'][0]['arms'][side]
            expected = [q['continuation']['sequence'] if 'continuation' in q else None
                        for q in arm['paths']]
            self.assertTrue(any(expected), (name, side))
        self.assertEqual(expected, got)
        for p, seq in zip(derive.paths(a), got):
            if seq is None:
                continue
            c = a.segments[p.leaf].leaf.continuation
            self.assertEqual(c.n, len(seq))
            flank = derive.natural_flank(a, p.leaf) if a.segments[0].walk is not None else None
            if flank is None:
                continue
            whole = g.seed.sequence + flank if side == 'right' else flank + g.seed.sequence
            self.assertEqual(whole[-c.n:] if side == 'right' else whole[:c.n], seq)

    def test_crossing_into_the_seed(self):
        # radius 10 < k 15: every continuation takes 5 bases of the seed
        g = T.body('fork')
        for side in ('right', 'left'):
            self.check('fork', side)
            for p in derive.paths(g.arms[side]):
                seq = derive.continuation_sequence(g, g.arms[side], p.leaf)
                self.assertEqual(15, len(seq))
                self.assertEqual(g.seed.sequence[-5:] if side == 'right' else g.seed.sequence[:5],
                                 seq[:5] if side == 'right' else seq[-5:])

    def test_contained_in_the_flank(self):
        self.check('merge', 'right')
        self.check('annotate', 'right')
        for name, side in (('merge', 'right'), ('annotate', 'right')):
            g = T.body(name)
            for p in derive.paths(g.arms[side]):
                c = g.arms[side].segments[p.leaf].leaf.continuation
                if c is not None:
                    self.assertLessEqual(c.n, p.length_bp)

    def test_explicit_sequence_without_bases(self):
        self.check('limits', 'right', ['ACAAG'])

    def test_semantic_leaves_have_none(self):
        self.check('reentry', 'right', [None])                # a dead end
        self.check('caps', 'right')                           # a truncated leaf among semantic ones
        self.check('switch_chain', 'right')                   # both at the radius

    def test_left_natural_orientation(self):
        g = T.body('fork')
        a = g.arms['left']
        self.assertEqual('GCACTCCCTA', derive.walk_bases(a, 1))
        self.assertEqual('ATCCCTCACG', derive.natural_flank(a, 1))
        # the server's segments hold the left flank leaf -> root, each reversed
        segs = T.doc_json('fork', 'full')['results'][0]['arms']['left']['segments']
        self.assertEqual(segs[1]['sequence'] + segs[0]['sequence'], derive.natural_flank(a, 1))


class TestEvidence(unittest.TestCase):
    def test_anchored_endpoint_through_two_merges(self):
        g = T.body('merge_ties')
        a = g.arms['right']
        got = [derive.evidence(a, r) for r in a.runs]
        # b.fa came through the second parent at 36 and at 62 (ties keep the arrival order)
        self.assertEqual([(0, 0), (62, 62), (0, 0), (0, 0), (0, 0)], got)
        # R keeps the EARLIEST stamp, T and the evidence the latest
        self.assertEqual(36, a.runs[1].route_bp)
        # the merge fixture: a.fa through the second parent at 62, b.fa at 36 only (each
        # merge's first parent is the majority's, R21 (4))
        a = T.body('merge').arms['right']
        self.assertEqual([(62, 62), (36, 36), (0, 0), (0, 0), (0, 0), (0, 0)],
                         [derive.evidence(a, r) for r in a.runs])

    def test_switch_moves_the_start_not_the_route(self):
        g = T.body('switch_chain')
        a = g.arms['right']
        self.assertEqual([(0, 0), (0, 30), (0, 60), (0, 60)],
                         [derive.evidence(a, r) for r in a.runs])

    def test_incoming_label_before_a_switch_at_the_merge_depth(self):
        g = T.body('switch_chain')
        a = g.arms['right']
        # run 2 (D) started at 60 by a switch from B: the incoming label at 60 is B, at 30 A
        self.assertEqual(1, derive.lineage_label_at(a, a.runs[2], 60))
        self.assertEqual(0, derive.lineage_label_at(a, a.runs[2], 30))
        self.assertEqual(3, derive.lineage_label_at(a, a.runs[2], 61))


class TestCheckRules(unittest.TestCase):
    def test_noncanonical_explicit_fields_are_reported(self):
        root = 'G * 0 40 * * * * * * TTACACCACGAGATCGCCTGGGGCTCTGACAGTTAGCATA'
        t = T.doc_text('linear').replace(root, root.replace(' * * * * * * ', ' 0-2 3 * * * * '))
        found = {(f['field'], f['kind']) for f in parse(t).check_rules()}
        self.assertEqual({('entry', 'noncanonical'), ('entry_total', 'noncanonical')}, found)

    def test_star_against_the_walkers_value(self):
        g = T.graphlet('reentry')
        ref = T.doc_json('reentry', 'full')['results'][0]
        ref['arms']['right']['segments'][0]['labels_at_end'] = []   # what v1 would give
        (f,) = g.check_rules(ref)
        self.assertEqual(('end', 'star_mismatch', [0], []),
                         (f['field'], f['kind'], f['decoded'], f['reference']))

    def test_invariants(self):
        g = T.body('merge_ties')
        g.arms['right'].segments[6].leaf.extras[1] = (0.0, 0, 36)   # the earliest stamp
        g.arms['right'].cache.clear()        # derivations cache: a parsed graphlet is immutable
        self.assertTrue(any(f.get('field') == 'route_bp' for f in g.check_rules()))
        g = T.body('switch_chain')
        g.arms['right'].runs[2].prev_run = None
        g.arms['right'].runs[3].prev_run = None
        g.arms['right'].cache.clear()
        self.assertTrue(any(f.get('run') == 1 and f['kind'] == 'invariant'
                            for f in g.check_rules()))
        g = T.body('merge')
        g.arms['right'].runs[4].segment = 0       # not a parent of the merge at 36
        g.arms['right'].cache.clear()
        self.assertTrue(any(f.get('run') == 4 for f in g.check_rules()))


if __name__ == '__main__':
    unittest.main()
