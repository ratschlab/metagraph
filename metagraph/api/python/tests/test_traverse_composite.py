"""The shared machinery of graphlet combinations (composite.py) and the ops.py pieces it
reads: the index rule split from the seed rule (no behaviour change) and the linear
per-segment label-interval accessor (docs/DESIGN-graphlet-merge.md, Appendix F.2)."""

import copy
import unittest

from tests import traverse_testlib as T  # noqa: F401  (puts the library on sys.path)

from metagraph.traverse import CombineRefused, LocalBudget, LocalBudgetExceeded
from metagraph.traverse import composite as C
from metagraph.traverse import ops

FIXTURES = sorted(T.DOCUMENTS)


def _with(g, **kw):
    """A shallow copy of graphlet |g| with fields replaced (identity edits, as W5's
    envelope tests edit H)."""
    h = copy.copy(g)
    for k, v in kw.items():
        setattr(h, k, v)
    return h


def _old_identity(a, b):
    """compare()'s identity rule as it was before the split (bc9684559, ops.py:3588): the
    oracle of 'no behaviour change'."""
    R = ops.COMPARE_REASONS
    ia, ib = ops.index_identity(a), ops.index_identity(b)
    if ia['meta_fp'] and ib['meta_fp'] and ia['meta_fp'] != ib['meta_fp']:
        return False, R['meta_fp'], R['meta_fp']
    if ia['fp'] and ib['fp'] and ia['fp'] != ib['fp']:
        return False, R['fp'], R['fp']
    for k in ('ns', 'release'):
        if ia[k] and ib[k] and ia[k] != ib[k]:
            return ops._valued(False, k, 'different index %s (%r vs %r)' % (k, ia[k], ib[k]))
    if a.seed.sequence != b.seed.sequence:
        return False, R['seed'], R['seed']
    names_a = {l.ref: l.name for l in a.labels}
    for l in b.labels:
        if l.ref in names_a and names_a[l.ref] != l.name:
            return ops._valued(False, 'names', 'label %s is %r in one and %r in the other'
                               % (l.ref, names_a[l.ref], l.name))
    if not ia['fp'] or not ib['fp']:
        why = R['unverifiable_both' if not ia['fp'] and not ib['fp'] else 'unverifiable_one']
        return 'unverifiable', why, why
    return True, R['same'], R['same']


class TestIdentitySplit(unittest.TestCase):

    def _variants(self):
        base = T.graphlet('fork')
        other = T.graphlet('linear')           # another seed and other label names
        renamed = _with(base, labels=[_with(l, name=l.name + 'x') for l in base.labels])
        out = [base, other, renamed]
        for kw in ({'index_fp': 'f1'}, {'index_fp': 'f2'}, {'index_meta_fp': 'm1'},
                   {'index_meta_fp': 'm2'}, {'index_ns': 'ns1'}, {'index_ns': 'ns2'}):
            out.append(_with(base, **kw))
            out.append(_with(other, **kw))
            out.append(_with(renamed, **kw))
        rel = _with(base, envelope=dict(base.envelope, release='r1'))
        out.append(rel)
        out.append(_with(other, envelope=dict(other.envelope, release='r2')))
        return out

    def test_identity_unchanged_on_every_pair(self):
        vs = self._variants()
        for a in vs:
            for b in vs:
                self.assertEqual(ops._identity(a, b), _old_identity(a, b))

    def test_index_rule_is_the_identity_without_its_seed_rule(self):
        vs = self._variants()
        for a in vs:
            for b in vs:
                # a with b's seed: the same index identity and labels, so the seed rule
                # passes and the old identity rule is the index rule alone
                self.assertEqual(ops.index_rule(a, b), _old_identity(_with(a, seed=b.seed), b))
                self.assertEqual(ops.seed_rule(a, b) is None,
                                 a.seed.sequence == b.seed.sequence)

    def test_order_proof_then_seed_then_names(self):
        a = T.graphlet('fork')
        renamed = _with(a, labels=[_with(l, name=l.name + 'x') for l in a.labels])
        other_seed = _with(renamed, seed=_with(a.seed, sequence='A' * len(a.seed.sequence)))
        # seed and names both differ: the seed rule speaks first in compare()
        self.assertEqual(ops._identity(a, other_seed)[1], ops.COMPARE_REASONS['seed'])
        # the index rule alone names the label
        self.assertEqual(ops.index_rule_check(a, other_seed)[0], 'names')
        # a proof of two indexes comes before the seed
        proven = _with(other_seed, index_fp='f2')
        self.assertEqual(ops._identity(_with(a, index_fp='f1'), proven)[1],
                         ops.COMPARE_REASONS['fp'])

    def test_index_rule_check_keys(self):
        a = T.graphlet('fork')
        self.assertEqual(C.index_check(a, a), ('unverifiable_both', ()))
        self.assertEqual(C.index_check(_with(a, index_fp='f'), a), ('unverifiable_one', ()))
        self.assertEqual(C.index_check(_with(a, index_fp='f'), _with(a, index_fp='f')),
                         ('same', ()))
        self.assertEqual(C.index_check(_with(a, index_meta_fp='m1'),
                                       _with(a, index_meta_fp='m2')), ('meta_fp', ()))
        self.assertEqual(C.index_check(_with(a, index_ns='x'), _with(a, index_ns='y')),
                         ('ns', ('x', 'y')))
        key, values = C.index_check(
            a, _with(a, labels=[_with(a.labels[0], name='zz')] + a.labels[1:]))
        self.assertEqual(key, 'names')
        self.assertEqual(C.different_index_clause(key, values),
                         'label %s is %r in one and %r in the other'
                         % (a.labels[0].ref, a.labels[0].name, 'zz'))
        self.assertEqual(C.different_index_clause('ns', ('x', 'y')),
                         "index namespace 'x' vs 'y'")
        self.assertEqual(C.different_index_clause('fp', ()), 'index_fp differs')
        self.assertIsNone(C.different_index_clause('unverifiable_one', ()))


class TestLabelIntervals(unittest.TestCase):

    def test_equal_to_alive_at_at_every_base(self):
        for name in FIXTURES:
            g = T.body(name)
            for side, arm in g.arms.items():
                for s in arm.segments:
                    iv = ops.label_intervals(g, side, s.id)
                    for lab, f, t in iv:
                        self.assertTrue(s.from_bp <= f < t <= s.end_bp, (name, s.id))
                    for pos in range(s.from_bp, s.end_bp):
                        got = frozenset(l for l, f, t in iv if f <= pos < t)
                        self.assertEqual(got, ops.alive_at(g, side, s.id, pos),
                                         (name, side, s.id, pos))

    def test_order_and_disjoint_per_label(self):
        for name in FIXTURES:
            g = T.body(name)
            for side, arm in g.arms.items():
                for s in arm.segments:
                    iv = ops.label_intervals(g, side, s.id)
                    self.assertEqual(iv, sorted(iv, key=lambda x: (x[1], x[0], x[2])))
                    by = {}
                    for lab, f, t in iv:
                        by.setdefault(lab, []).append((f, t))
                    for spans in by.values():
                        for (f1, t1), (f2, t2) in zip(spans, spans[1:]):
                            self.assertLessEqual(t1, f2, (name, s.id))   # disjoint

    def test_reentry_keeps_order_and_extent(self):
        g = T.body('reentry')
        got = [(g.labels[l].name, f, t) for l, f, t in ops.label_intervals(g, 'right', 0)]
        self.assertEqual(got, [('pA.fa', 0, 30), ('pB.fa', 30, 60), ('pA.fa', 60, 90)])

    def test_linear_in_entry_and_changes(self):
        for name in FIXTURES:
            g = T.body(name)
            for side, arm in g.arms.items():
                for s in arm.segments:
                    n = len(ops.label_intervals(g, side, s.id))
                    if g.mode == 'constrain':
                        self.assertLessEqual(n, len(s.entry) + len(ops._segment_ops(arm, s)))
                    else:
                        self.assertLessEqual(n, sum(len(p.labels) for p in s.presence))

    def test_zero_length_segment_has_none(self):
        g = T.body('limits')
        for side, arm in g.arms.items():
            for s in arm.segments:
                if s.length_bp == 0:
                    self.assertEqual(ops.label_intervals(g, side, s.id), [])

    def test_cut_lists(self):
        g = T.body('annotate')               # cap 2: cut lists
        cut = [(side, s.id) for side, arm in g.arms.items() for s in arm.segments
               if ops.label_intervals_cut(g, side, s.id)]
        self.assertTrue(cut)
        for side, sid in cut:
            self.assertTrue(any(p.truncated for p in g.arms[side].segments[sid].presence))
        f = T.body('fork')
        self.assertFalse(any(ops.label_intervals_cut(f, side, s.id)
                             for side, arm in f.arms.items() for s in arm.segments))

    def test_budget_charged_before_and_answer_unchanged(self):
        g = T.body('reentry')
        want = ops.label_intervals(g, 'right', 0)
        b = LocalBudget(work_units=10 ** 6)
        self.assertEqual(ops.label_intervals(g, 'right', 0, budget=b), want)
        self.assertGreater(b.used_work, 0)
        with self.assertRaises(LocalBudgetExceeded):
            ops.label_intervals(g, 'right', 0, budget=LocalBudget(work_units=1))


class TestComposite(unittest.TestCase):

    def test_refusals(self):
        e = CombineRefused('duplicate_pair', 'text')
        self.assertIsInstance(e, ValueError)
        self.assertEqual((e.code, str(e)), ('duplicate_pair', 'text'))
        n = C.NotCombinable('not this')
        self.assertIsInstance(n, TypeError)
        self.assertIsInstance(n, CombineRefused)
        self.assertEqual(n.code, 'not_combinable')

    def test_preconditions(self):
        g = T.graphlet('fork')
        v = g.subgraph([{'id': 0}])
        self.assertTrue(C.is_view(v))
        with self.assertRaises(CombineRefused) as cm:
            C.require_whole(v, 'no views')
        self.assertEqual((cm.exception.code, str(cm.exception)), ('view_unsupported', 'no views'))
        C.require_whole(g, 'no views')
        from metagraph.traverse import MissingEnvelope
        with self.assertRaises(MissingEnvelope) as cm:
            C.require_envelope(T.body('fork'), 'union_pairs()')
        self.assertTrue(str(cm.exception).startswith('union_pairs() needs the response'))
        self.assertTrue(C.has_bases(g))
        self.assertFalse(C.has_bases(T.body('caps')))      # output.sequences false
        self.assertEqual(C.alphabet_key(g), 'ACGT')
        self.assertEqual(C.case_mapped('acgT', '$ACGT'), 'ACGT')
        self.assertEqual(C.case_mapped('acgT', 'acgt'), 'acgT')

    def test_seed_difference(self):
        self.assertIsNone(C.seed_difference('ACGT', 'ACGT', 'ACGT'))
        self.assertEqual(C.seed_difference('AACG', 'CGTT', 'ACGT'), ('revcomp',))
        self.assertEqual(C.seed_difference('AACG', 'CGTT', 'ACDEFGHIK'), ('first', 0))
        self.assertEqual(C.seed_difference('ACGT', 'ACGA', 'ACGT'), ('first', 3))
        self.assertEqual(C.seed_difference('ACGT', 'ACG', 'ACGT'), ('length', 4, 3, 3))
        self.assertEqual(C.seed_difference('ACG', 'ACTT', 'ACGT'), ('length', 3, 4, 2))
        self.assertEqual(C.reverse_complement('AACGn'), 'nCGTT')

    def test_first_difference(self):
        a = {'b': {'x': 1, 'y': [1, 2]}, 'a': 1}
        self.assertIsNone(C.first_difference(a, copy.deepcopy(a)))
        b = copy.deepcopy(a)
        b['b']['y'][1] = 3
        b['b']['x'] = 2
        self.assertEqual(C.first_difference(a, b), 'b.x')       # byte order of the keys
        del b['b']['x']
        self.assertEqual(C.first_difference(a, b), 'b.x')
        b['b']['x'] = 1
        self.assertEqual(C.first_difference(a, b), 'b.y[1]')
        self.assertEqual(C.first_difference({'a': 1}, {'a': True}), 'a')
        self.assertEqual(C.first_difference([1], [1, 2]), '[1]')

    def test_retrieval_identity(self):
        g = T.graphlet('fork')
        r = C.retrieval_identity(g)
        self.assertEqual(sorted(r), ['index', 'seed', 'strategy'])
        self.assertEqual(r['seed']['sequence'], g.seed.sequence)
        self.assertNotIn('seed_id', C.canonical_json(r['seed']))
        # the caller's seed_id is left out; a knob nothing reached is not
        h = copy.copy(g)
        h.seed_summary = copy.deepcopy(g.seed_summary)
        h.seed_summary['seed']['seed_id'] = 'another'
        self.assertEqual(C.retrieval_identity(h), r)
        env = copy.deepcopy(g.envelope)
        env['strategy'].setdefault('branching', {})['max_splits_per_path'] = 128
        h2 = _with(g, envelope=env)
        self.assertEqual(C.first_difference(r, C.retrieval_identity(h2)),
                         'strategy.branching.max_splits_per_path')

    def test_weakest_outcome(self):
        self.assertIsNone(C.weakest_outcome([None]))
        got = C.weakest_outcome([
            {'walks': 'complete', 'branch_diagnostics': 'cut', 'label_evidence': 'complete',
             'delivery': 'inline'},
            None,
            {'walks': 'partial', 'branch_diagnostics': 'complete',
             'label_evidence': 'lower_bound', 'delivery': 'spooled'}])
        self.assertEqual(got, {'walks': 'partial', 'branch_diagnostics': 'cut',
                               'label_evidence': 'lower_bound', 'delivery': 'spooled'})
        self.assertEqual(list(got), list(C.OUTCOME_RANKS))

    def test_combined_evidence(self):
        a, b = T.graphlet('fork'), T.graphlet('annotate')
        ev = C.combined_evidence([('a', a, True), ('b', b, False)], 'right')
        self.assertEqual(list(ev['parts']), ['a', 'b'])
        self.assertEqual(ev['parts']['a'], ops.evidence_block(a, 'right'))
        self.assertFalse(ev['exact'])              # b unverified (and its lists are cut)
        ev1 = C.combined_evidence([('a', a, True)])
        self.assertEqual(ev1['exact'], ops.evidence_block(a)['exact'])
        self.assertFalse(C.combined_evidence([('a', a, False)])['exact'])

    def test_digests_and_tokens(self):
        g = T.graphlet('fork')
        self.assertEqual(C.part_digest(g), C.part_digest(T.body('fork')))
        self.assertNotEqual(C.part_digest(g), C.part_digest(T.body('linear')))
        d = C.combined_digest([['x', 1], {'b': 1, 'a': 2}])
        self.assertEqual(d, C.combined_digest([['x', 1], {'a': 2, 'b': 1}]))
        tok = C.token(d, 'walks', ['right'], 5)
        self.assertEqual(C.resume_pos(d, 'walks', ['right'], tok, 1), (5,))
        self.assertEqual(C.resume_pos(d, 'walks', ['right'], None, 2), (0, 0))
        with self.assertRaises(ValueError):
            C.resume_pos(d, 'walks', ['left'], tok, 1)
        with self.assertRaises(ValueError):
            C.resume_pos('other', 'walks', ['right'], tok, 1)
        with self.assertRaises(ValueError):
            C.resume_pos(d, 'segments', ['right'], tok, 1)

    def test_memory_counts_a_model_once(self):
        g, h = T.body('fork'), T.body('linear')
        one = C.models_memory_bytes([g])
        self.assertEqual(C.models_memory_bytes([g, g]), one)
        self.assertGreater(C.models_memory_bytes([g, h]), one)
        self.assertEqual(C.models_memory_bytes([g], extra=10), one + 10)

    def test_spell_pieces_equals_the_parts_spelling(self):
        from metagraph.traverse import derive
        for name in FIXTURES:
            g = T.body(name)
            if not C.has_bases(g):
                continue
            for side, arm in g.arms.items():
                for w in g.walks(side):
                    chain = derive.chain(arm, ops.leaf_segment(arm, w.path_id))
                    pieces = [(g, side, x, arm.segments[x].from_bp, arm.segments[x].end_bp)
                              for x in chain]
                    for o in ('natural', 'walk'):
                        self.assertEqual(C.spell_pieces(pieces, side, o),
                                         g.spell(side, w.path_id, o), (name, side, o))
        g = T.body('fork')
        s = g.arms['right'].segments[1]
        part = C.spell_pieces([(g, 'right', 1, s.from_bp + 1, s.end_bp)], 'right', 'walk')
        self.assertEqual(part, s.walk[1:])
        with self.assertRaises(ValueError):
            C.spell_pieces([(g, 'right', 1, s.from_bp, s.end_bp + 1)], 'right')
        with self.assertRaises(ValueError):
            C.spell_pieces([(T.body('caps'), 'right', 0, 0, 1)], 'right')
        b = LocalBudget(work_units=10 ** 6)
        with b.scope('spell_pieces'):
            self.assertEqual(C.spell_pieces([(g, 'right', 1, s.from_bp, s.end_bp)], 'right',
                                            b=b), s.walk)
        self.assertGreater(b.used_work, 0)


if __name__ == '__main__':
    unittest.main()
