"""The local queries (§5, §5.1): label selectors, spelling and orientation, walks,
claims with the anchored evidence (incl. the clipped-merge case), label walks and
routes, support profiles, continuations and next_request, views, comparisons, exports."""

import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402

from metagraph.traverse import (  # noqa: E402
    AmbiguousLabel, GraphletView, IncompleteRecording, MissingEnvelope, UnknownLabel,
    load, parse,
)
from metagraph.traverse import ops  # noqa: E402


def refs(labels):
    return [l.ref for l in labels]


class TestLabels(unittest.TestCase):
    def test_tagged_and_bare_selectors(self):
        g = T.body('same_name')
        self.assertEqual(1, g.label({'ref': 'h:1:0'}).id)
        self.assertEqual(2, g.label({'name': 'ACC2'}).id)
        self.assertEqual(2, g.label('ACC2').id)          # unique as a name
        self.assertEqual(0, g.label('h:0:0').id)         # unique as a ref
        self.assertEqual(1, g.label({'id': 1}).id)
        with self.assertRaises(AmbiguousLabel):
            g.label('ACC1')                              # two headers from two columns
        with self.assertRaises(AmbiguousLabel):
            g.label({'name': 'ACC1'})
        with self.assertRaises(UnknownLabel):
            g.label('ACC9')
        with self.assertRaises(TypeError):
            g.label({'ref': 'h:0:0', 'name': 'ACC1'})

    def test_a_name_that_looks_like_a_ref(self):
        g = T.body('limits')                   # column labels c:0 .. c:3
        g.labels[2].name = 'c:0'               # a column may be NAMED like another's ref
        with self.assertRaises(AmbiguousLabel):
            g.label('c:0')
        self.assertEqual(0, g.label({'ref': 'c:0'}).id)
        self.assertEqual(2, g.label({'name': 'c:0'}).id)

    def test_results_carry_name_and_ref(self):
        g = T.graphlet('fork')
        for c in g.claims():
            self.assertEqual({'name', 'ref'}, set(c.as_dict()['label']))
        row = g.label_summary()['h:0:0']
        self.assertEqual(('acc1', 'h:0:0'), (row['name'], row['ref']))


class TestSpelling(unittest.TestCase):
    def test_orientations(self):
        g = T.body('fork')
        seed = g.seed.sequence
        self.assertEqual('ATTGCTAAGA', g.spell('right', 0))
        self.assertEqual(seed + 'ATTGCTAAGA', g.spell('right', 0, with_seed=True))
        self.assertEqual('ATCCCTCACG', g.spell('left', 0))
        self.assertEqual('ATCCCTCACG' + seed, g.spell('left', 0, with_seed=True))
        self.assertEqual('GCACTCCCTA', g.spell('left', 0, orientation='walk'))
        self.assertEqual(seed[::-1] + 'GCACTCCCTA', g.spell('left', 0, 'walk', True))
        # T38: natural(left) + seed + natural(right) is one molecule -- here acc1's
        # (L1 E R1), spelled by the walks that carry it on each arm
        whole = g.spell('left', 0) + seed + g.spell('right', 1)
        rec = T.cli_record('fork', 'acc1')
        self.assertIn(whole, rec)

    def test_left_arm_positions_index_walking_order(self):
        g = T.body('fork')
        a = g.arms['left']
        walk = g.spell('left', 0, orientation='walk')
        for i in range(len(walk)):
            seg = next(s for s in a.segments if s.from_bp <= i < s.end_bp and s.id in (0, 1))
            self.assertEqual(walk[i], seg.walk[i - seg.from_bp])

    def test_no_bases(self):
        with self.assertRaises(ValueError):
            T.body('limits').spell('right', 0)


class TestWalks(unittest.TestCase):
    def test_support_ranking_and_fields(self):
        g = T.body('fork')
        ws = g.walks('right')
        self.assertEqual([1, 0], [w.path_id for w in ws])   # acc1 and acc3 before acc2
        w = ws[0]
        self.assertEqual((10, 'max_extension_bp', 2, True, 0),
                         (w.length_bp, w.path_reason, w.n_alive, w.complete,
                          w.beyond_certified_bp))
        self.assertEqual(['h:0:0', 'h:0:2'], refs(w.labels_full))
        self.assertEqual('TTACACCACG', w.sequence)
        self.assertEqual(['alive', 'alive'], [c.kind for c in w.claims])

    def test_route_consistent_label_filter(self):
        g = T.body('merge')
        self.assertEqual([], g.walks('right', labels=['b.fa']))       # merged in at 36
        self.assertEqual([], g.walks('right', labels=['a.fa']))       # merged in at 62
        self.assertEqual([0], [w.path_id for w in g.walks('right', labels=['b.fa'],
                                                            route_consistent=False)])
        (w,) = g.walks('right')
        # c.fa and both.fa are on both first parents (the majority's)
        self.assertEqual(['c:2', 'c:3'], refs(w.labels_full))
        self.assertEqual(4, w.n_alive)
        (w,) = T.body('merge_ties').walks('right')
        self.assertEqual(['c:0', 'c:3'], refs(w.labels_full))         # arrival order

    def test_rankings_and_filters(self):
        g = T.body('caps')
        # a truncated walk of 126 bp and two dead ends of 109 and 95 bp
        self.assertEqual([0, 1, 2], [w.path_id for w in g.walks('right', by='length')])
        self.assertEqual([0, 1], [w.path_id for w in g.walks('right', by='length', min_bp=100)])
        self.assertEqual([0], [w.path_id for w in g.walks('right', by='length', top=1)])
        with self.assertRaises(ValueError):
            g.walks('right', by='nonsense')

    def test_capped_walks_are_not_complete(self):
        (w,) = T.body('limits').walks('right')
        self.assertEqual(('max_steps', False, None), (w.path_reason, w.complete, w.sequence))

    def test_annotate_walks(self):
        g = T.body('annotate')
        (w,) = g.walks('right')
        # no label is in every (cut) recorded list along the displayed walk: it passes the
        # second bubble through the G allele (b.fa, c.fa, both.fa: the parent present in the
        # most labels), whose list (cut to b.fa, c.fa) lacks a.fa, and the lists cut to
        # a.fa, b.fa at the merge nodes lack c.fa
        self.assertEqual([], refs(w.labels_full))
        self.assertEqual(2, w.n_alive)
        (w,) = g.walks('left')
        self.assertEqual(['c:0', 'c:1'], refs(w.labels_full))


class TestClaims(unittest.TestCase):
    def test_clipped_merge_expected_claims(self):
        """§5.1: provenance at the anchored endpoint, then clipped; at D = 4 lineB is
        route_only (v2 evaluated the walk at D and gave evidence_from = 0)."""
        g = T.body('clipped_merge')
        with open(os.path.join(T.DOCS, 'clipped_merge.claims.json')) as f:
            want = json.load(f)
        for key, cut in (('all', None), ('at_4', 4), ('at_16', 16)):
            got = [{'ref': c.label.ref, 'run': c.run, 'kind': c.kind, 'from_bp': c.from_bp,
                    'to_bp': c.to_bp, 'route_bp': c.route_bp,
                    'evidence_from': c.evidence_from, 'end_class': c.end_class,
                    'reason': c.reason, 'loss': c.loss, 'branches': c.branches}
                   for c in g.claims('right', at_most_bp=cut)]
            self.assertEqual(want['claims'][key], got, key)
        for ref, routes in want['routes'].items():
            self.assertEqual(routes, g.routes({'ref': ref}, 'right'))
        (route, bases), = g.routes({'ref': 'h:0:2'}, 'right', spell=True)
        self.assertEqual(want['route_bases']['h:0:2'], bases)
        self.assertEqual(want['displayed'], g.spell('right', 0))
        for kind in ('displayed', 'route'):
            got = [[r.from_bp, r.to_bp, refs(r.labels)]
                   for r in g.support_profile('right', 0, kind=kind)]
            self.assertEqual(want['support_profile_' + kind], got)

    def test_run_start_guard(self):
        g = T.body('switch_chain')
        # B enters by a switch at 30 and C, D at 60: at D = 30 neither exists yet
        at30 = g.claims('right', at_most_bp=30)
        self.assertEqual([(0, 'diverged', 30)], [(c.run, c.kind, c.to_bp) for c in at30])
        at31 = g.claims('right', at_most_bp=31)
        self.assertEqual([(0, 'diverged', 0, 30), (1, 'stretch', 30, 31)],
                         [(c.run, c.kind, c.evidence_from, c.to_bp) for c in at31])
        self.assertIsNone(at31[1].loss)        # terminal values are unknown at a cut

    def test_boundary_claim(self):
        g = T.body('linear')
        acc2 = [c for c in g.claims('right') if c.label.name == 'acc2']
        self.assertEqual([('boundary', 0, 0, 'pruned', 'branch', 'minority')],
                         [(c.kind, c.from_bp, c.to_bp, c.end_class, c.reason, c.qualifier)
                          for c in acc2])
        self.assertEqual(['boundary'], [c.kind for c in g.claims('right', labels=['acc2'],
                                                                 at_most_bp=0)])
        # a run with bases makes no claim at a cut of 0 (the guard), in both modes
        self.assertEqual([], g.claims('right', labels=['acc1'], at_most_bp=0))
        self.assertEqual([], T.body('same_name').claims('right', labels=[{'ref': 'h:0:0'}],
                                                        at_most_bp=0))

    def test_kinds(self):
        kinds = {c.run: c.kind for c in T.body('merge').claims('right')}
        # both.fa: runs 3 and 4 closed by the merges, run 5 kept
        self.assertEqual({0: 'alive', 1: 'alive', 2: 'alive', 3: 'merged', 4: 'merged',
                          5: 'alive'}, kinds)
        kinds = [c.kind for c in T.body('switch_chain').claims('right')]
        self.assertEqual(['diverged', 'diverged', 'alive', 'alive'], kinds)
        c = T.body('switch_chain').claims('right')[3]
        self.assertEqual(('switch', 'h:0:1', 1.0, 2.0, 1),
                         (c.entered_by, c.from_label.ref, c.cost, c.loss, c.branches))
        for name in ('limits', 'caps'):
            (b,) = [c for c in T.body(name).claims('right') if c.reason == 'loss_budget']
            self.assertEqual(('end', 'pruned'), (b.kind, b.end_class))
        kinds = {c.kind for c in T.body('reentry').claims('right')}
        self.assertEqual({'diverged', 'end'}, kinds)

    def test_annotate_claims_are_the_oracle(self):
        g = T.body('same_name')
        got = [(c.label.ref, c.kind, c.to_bp, c.end_class) for c in g.claims()]
        # the two ACC1 records each run to their own end; ACC2 is not on the seed
        self.assertEqual([('h:1:0', 'end', 30, 'dead_end'), ('h:0:0', 'end', 30, 'dead_end')],
                         got)
        g = T.body('annotate')
        with self.assertRaises(IncompleteRecording):
            g.claims('right')
        self.assertTrue(all(not c.exact for c in g.claims('right', strict=False)))


class TestLabelWalksAndSupport(unittest.TestCase):
    def test_label_walks(self):
        g = T.body('merge')
        lw = g.label_walks('both.fa')
        # the two runs closed by the merges at 62 and 36 and the kept one, each with its
        # route: the kept lineage is the one through the first (majority) parents
        self.assertEqual([(3, 0, 62, 'merged', 6), (4, 0, 36, 'merged', 3),
                          (5, 0, 100, 'max_extension_bp', None)],
                         [(x.run, x.from_bp, x.to_bp, x.end, x.merged_into) for x in lw])
        self.assertEqual(['A', 'T', 'A'], [x.sequence[20] for x in lw])
        self.assertEqual(['C', 'G'], [lw[0].sequence[46], lw[2].sequence[46]])
        self.assertEqual(g.spell('right', 0), lw[2].sequence)
        (b,) = g.label_walks('b.fa')
        # b.fa's own route (T at 20, G at 46); the displayed bases are its only from 36
        self.assertEqual((36, 'T', 'G'), (b.evidence_from, b.sequence[20], b.sequence[46]))
        self.assertNotEqual(g.spell('right', 0), b.sequence)
        (a,) = g.label_walks('a.fa')
        # a.fa's own route (A at 20, C at 46): displayed only from the merge at 62
        self.assertEqual((62, 'A', 'C'), (a.evidence_from, a.sequence[20], a.sequence[46]))
        self.assertEqual(g.spell('right', 0)[62:], a.sequence[62:])
        sw = T.body('switch_chain').label_walks('B')
        self.assertEqual([('switched', 30, [1, 2])],
                         [(x.end, len(x.sequence), x.leaves_below) for x in sw])

    def test_support_changes_with_reasons(self):
        g = T.body('merge')
        ch = g.support_changes('right', 0)
        # b.fa leaves at the first bubble and returns at its merge, a.fa at the second (the
        # displayed walk follows the majority's G allele there)
        self.assertEqual([(20, [], ['c:1']), (36, ['c:1'], []), (46, [], ['c:0']),
                          (62, ['c:0'], [])],
                         [(c.at_bp, refs(c.added), refs(c.removed)) for c in ch])
        self.assertEqual({'why': 'split', 'took': ['T']},
                         {k: ch[0].reasons[0][k] for k in ('why', 'took')})
        self.assertEqual(('merge', 2), (ch[1].reasons[0]['why'], ch[1].reasons[0]['via_parent']))
        self.assertEqual({'why': 'split', 'took': ['C']},
                         {k: ch[2].reasons[0][k] for k in ('why', 'took')})
        self.assertEqual(('merge', 4), (ch[3].reasons[0]['why'], ch[3].reasons[0]['via_parent']))
        sw = T.body('switch_chain').support_changes('right', 1)
        self.assertEqual(['switched', 'switch', 'switched', 'switch'],
                         [r['why'] for c in sw for r in c.reasons])
        lin = T.body('linear').support_changes('right', 0)
        self.assertEqual([(0, ['h:0:1'], 'minority')],
                         [(c.at_bp, refs(c.removed), c.reasons[0]['why']) for c in lin])

    def test_route_profile_adds_the_runs_own_stretch(self):
        g = T.body('merge')
        disp = [(r.from_bp, r.to_bp, refs(r.labels)) for r in g.support_profile('right', 0)]
        self.assertEqual([(0, 20, ['c:0', 'c:1', 'c:2', 'c:3']), (20, 36, ['c:0', 'c:2', 'c:3']),
                          (36, 46, ['c:0', 'c:1', 'c:2', 'c:3']),
                          (46, 62, ['c:1', 'c:2', 'c:3']),
                          (62, 100, ['c:0', 'c:1', 'c:2', 'c:3'])], disp)
        route = g.support_profile('right', 0, kind='route')
        self.assertEqual([(0, 100, ['c:0', 'c:1', 'c:2', 'c:3'])],
                         [(r.from_bp, r.to_bp, refs(r.labels)) for r in route])

    def test_annotate_profile_states_exactness(self):
        g = T.body('annotate')
        prof = g.support_profile('right', 0)
        # cut lists (cap 2) along the whole displayed walk: it passes the second bubble
        # through the G allele, present in 3 labels (the majority parent)
        self.assertEqual([(0, 20, False, 4), (20, 35, False, 3), (35, 46, False, 4),
                          (46, 61, False, 3), (61, 70, False, 4)],
                         [(r.from_bp, r.to_bp, r.exact, r.total) for r in prof])
        # a row is exact exactly when its list is whole: the C allele's own stretch (a.fa,
        # both.fa: 2 labels, at the cap) is recorded whole, the merge node's is cut
        self.assertEqual([(46, 61, True, 2), (61, 62, False, 4)],
                         [(f, t, n == len(s), n) for f, t, s, n in ops._presence_profile(
                             g.arms['right'], 4, None, g) if f >= 46])


class TestContinuationAndNextRequest(unittest.TestCase):
    def test_continuation(self):
        g = T.body('fork')
        c = g.continuation('right', 0)
        # 15 = k bases: the last 5 of the seed (length 40) and the 10 of the walk
        self.assertEqual(('CTATTATTGCTAAGA', ['h:0:1'], (35, 50)),
                         (c.sequence, refs(c.labels), c.seed_coord))
        self.assertEqual({'sequence': 'CTATTATTGCTAAGA', 'labels': ['acc2']}, c.as_seed())
        c = g.continuation('left', 1)
        self.assertEqual(('GTCGGACGGGCGTAC', (-10, 5)), (c.sequence, c.seed_coord))
        with self.assertRaises(ValueError):
            T.body('caps').continuation('right', 1)     # a dead end

    def test_next_request(self):
        g = T.graphlet('switch_chain')
        req = g.next_request('right', [1], bp=50, bounds={'max_steps': 99})
        self.assertEqual([{'sequence': 'TTACTCGTAGCCGGGCGTGA', 'labels': ['C']}], req['seeds'])
        st = req['strategy']
        self.assertEqual('right', st['direction'])
        self.assertEqual(1.0, st['labels']['loss_budget'])          # 3 - loss_used 2
        self.assertEqual((50, 99), (st['bounds']['max_extension_bp'], st['bounds']['max_steps']))
        self.assertEqual('graphlet', st['output']['detail'])
        self.assertNotIn('clamped', st)
        self.assertNotIn('release', req)                            # unset: not pinned
        same = g.next_request('right', [1], reduce_budget=False, release='r7', graph='g1')
        self.assertEqual((3.0, 'r7', 'g1'), (same['strategy']['labels']['loss_budget'],
                                             same['release'], same['graph']))
        ann = T.graphlet('annotate').next_request('right', [0])
        self.assertEqual([{'sequence': 'CGGGTCGACTTGGTCCGGAC'}], ann['seeds'])
        two = T.graphlet('fork').next_request('left', [0, 1])
        self.assertEqual(['ATCCCTCACGCGTAC', 'GTCGGACGGGCGTAC'],
                         [s['sequence'] for s in two['seeds']])
        with self.assertRaises(MissingEnvelope):
            T.body('fork').next_request('left', [0])


class TestViews(unittest.TestCase):
    def test_subgraph_keeps_the_selected_chain_with_original_ids(self):
        g = T.graphlet('fork')
        v = g.subgraph(['acc2'])
        self.assertIsInstance(v, GraphletView)
        self.assertEqual({'left': [0, 1], 'right': [0, 1]}, v.segments)
        self.assertEqual({'left': [0], 'right': [0]}, v.path_ids)
        self.assertEqual([0], [w.path_id for w in v.walks('right')])
        self.assertEqual({'h:0:1'}, {c.label.ref for c in v.claims()})
        self.assertEqual('for the selected labels', v.completeness()['right']['qualified'])
        self.assertEqual(2, v.to_fasta().count('>'))
        only_right = g.subgraph(['acc2'], arm='right', mode='all')
        self.assertEqual(['right'], only_right.arms())

    def test_save_writes_the_unchanged_body_and_load_restores_the_view(self):
        g = T.graphlet('fork')
        v = g.subgraph([{'ref': 'h:0:1'}])
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'view.mgt')
            v.save(path)
            text = T.read(path)
            j = json.loads(text.split('\n')[1][2:])
            self.assertEqual({'selectors': [{'ref': 'h:0:1'}], 'mode': 'any', 'arm': None},
                             {k: j['view'][k] for k in ('selectors', 'mode', 'arm')})
            self.assertTrue(j['view']['of'].startswith('sha256:'))
            body = [l for i, l in enumerate(text.split('\n')) if i != 1]
            body[-2] = 'Z %d' % (len(body) - 1)
            self.assertEqual(T.doc_text('fork'), '\n'.join(body))
            back = load(path)
            self.assertIsInstance(back, GraphletView)
            self.assertEqual(v.path_ids, back.path_ids)
            self.assertTrue(back.backing.has_envelope)
            # a view of a body-only graphlet keeps its selector too, and no envelope
            bare = T.body('fork').subgraph(['acc2'], arm='left')
            bare.save(path)
            back = load(path)
            self.assertEqual(({'left': [0]}, 'left'), (back.path_ids, back.arm_side))
            self.assertFalse(back.backing.has_envelope)
            self.assertEqual(T.doc_text('fork'), back.backing.dump())
            self.assertEqual(T.doc_text('fork'), bare.dump(envelope=False))


class TestCompare(unittest.TestCase):
    """T39 on a hand-made locus: the constrained trie over {acc1, acc3} equals the
    annotate oracle filtered to them; a quorum-tuned run is a prefix subset of the
    exhaustive trie, the omitted walk carrying the recorded reason."""

    @classmethod
    def setUpClass(cls):
        cls.docs = {n: T.compare_text(n) for n in ('oracle_annotate', 'oracle_constrain',
                                                   'oracle_exhaustive', 'oracle_tuned')}

    def g(self, name, edit=None):
        t = self.docs[name]
        if edit:
            t = t.replace(*edit)
        return parse(t)

    def test_constrain_equals_the_annotate_oracle(self):
        c = self.g('oracle_constrain').compare(self.g('oracle_annotate'),
                                              labels=['acc1', 'acc3'])
        self.assertEqual((True, True, [], [], []),
                         (c.comparable, c.equal, c.only_in_a, c.only_in_b, c.differ))
        self.assertEqual(20, c.depth_used)
        for mode in ('claims', 'walks', 'labels'):
            c = self.g('oracle_exhaustive').compare(self.g('oracle_annotate'), mode=mode)
            self.assertTrue(c.equal, (mode, c.only_in_a, c.only_in_b, c.differ))

    def test_a_missing_claim_is_reported(self):
        c = self.g('oracle_constrain').compare(self.g('oracle_annotate'))
        self.assertFalse(c.equal)
        self.assertEqual([{'name': 'acc2', 'ref': 'h:0:2'}], [x['label'] for x in c.only_in_b])

    def test_prefix_subset_with_the_omission_reason(self):
        c = self.g('oracle_tuned').compare(self.g('oracle_exhaustive'), mode='prefix_subset')
        self.assertTrue(c.equal)
        self.assertEqual([{'arm': 'right', 'walk': 'TCCGGAAT', 'ref': 'h:0:2', 'name': 'acc2',
                           'at_bp': 0, 'reason': 'minority'}], c.only_in_b)
        rev = self.g('oracle_exhaustive').compare(self.g('oracle_tuned'), mode='prefix_subset')
        self.assertFalse(rev.equal)
        self.assertEqual(['h:0:2'], [x['ref'] for x in rev.only_in_a])

    def test_index_identity_and_seed(self):
        fp = 'a3f5b7c9d1e3f5a7b9c1d3e5f7a9b1c3d5e7f9a1b3c5d7e9f1a3b5c7d9e1f3a5'
        a = self.g('oracle_exhaustive')
        other_fp = self.g('oracle_annotate', (fp, 'b' * 64))
        c = a.compare(other_fp)
        self.assertEqual((False, None), (c.comparable, c.equal))
        no_fp = self.g('oracle_annotate', (' ' + fp + ' ', ' * '))
        c = a.compare(no_fp)
        self.assertEqual(('unverifiable', None), (c.comparable, c.equal))
        other_meta = self.g('oracle_annotate', ('3f2a9c1b0d4e5f60', '0000000000000000'))
        self.assertFalse(a.compare(other_meta).comparable)
        other_seed = self.g('oracle_annotate', ('CCATGGTACA', 'CCATGGTACC'))
        self.assertEqual('different oriented seed sequences',
                         a.compare(other_seed).reason)
        renamed = self.g('oracle_annotate', ('L h 0 3 3 3', 'L h 0 3 3 9'))
        self.assertFalse(a.compare(renamed).comparable)
        # validated_seed_id differs between the modes and does not matter
        self.assertNotEqual(a.seed.validated_seed_id,
                            self.g('oracle_annotate').seed.validated_seed_id)

    def test_scopes_qualify(self):
        merged = self.g('oracle_tuned', (' c k k 64 ', ' c k m 64 ')).dump()
        merged = parse(merged.replace('A r c 20 p ', 'A r c 20 u '))
        c = merged.compare(self.g('oracle_exhaustive'), mode='prefix_subset')
        self.assertEqual(('qualified', None), (c.comparable, c.equal))
        self.assertTrue(any('scopes differ' in n for n in c.notes))


def _with_fp(text, fp='cd' * 32):
    """A CLI fixture with an index manifest digest (they are made without one, so every
    join of two of them is 'unverifiable')."""
    return text.replace(' walk fixtures * ', ' walk fixtures %s ' % fp, 1)


def _mirrored(text):
    """The same retrieval on the LEFT arm of the mirrored index (every record reversed):
    bases are stored in walking order on both arms, so only the arm codes and the seed's
    orientation change."""
    out = []
    for line in text.split('\n'):
        if line.startswith('S '):
            f = line.split(' ')
            f[5] = f[5][::-1]
            line = ' '.join(f)
        elif line.startswith(('A r ', 'K r ')):
            line = line[:2] + 'l' + line[3:]
        out.append(line)
    return '\n'.join(out)


class TestCompareCli(unittest.TestCase):
    """compare() on CLI retrievals of one seed (documents/compare/cli_*.mgt): prefix_subset
    orientation, interior ends and recorded refusals, claims ending at the comparison depth,
    lower-bound label evidence and an empty window."""

    def g(self, name, left=False):
        t = _with_fp(T.compare_text(name))
        return parse(_mirrored(t) if left else t)

    def test_prefix_subset_reads_walks_outward_on_both_arms(self):
        for left in (False, True):
            tuned = self.g('cli_quorum_tuned', left)
            exhaustive = self.g('cli_quorum_exhaustive', left)
            names = {l.ref: l.name for l in tuned.labels}
            c = tuned.compare(exhaustive, mode='prefix_subset')
            # x.fa and z.fa end label_lost at 40 inside a walk that y.fa carries on: the
            # exhaustive trie holds the same ends, so they are no violations
            self.assertEqual((True, True, []), (c.comparable, c.equal, c.only_in_a), left)
            # the omitted walks carry the refusals the tuned body recorded in V, on the
            # left arm as on the right
            self.assertEqual([('x.fa', 'minority', 10), ('y.fa', 'minority', 40)],
                             sorted((names[x['ref']], x['reason'], x.get('at_bp'))
                                    for x in c.only_in_b), left)
            side = 'left' if left else 'right'
            self.assertEqual({side}, {x['arm'] for x in c.only_in_b})
            for x in c.only_in_b:    # rows spell walks in natural orientation
                n = len(x['walk'])
                self.assertIn(x['walk'], [w.sequence[:n] if side == 'right' else
                                          w.sequence[len(w.sequence) - n:]
                                          for w in exhaustive.walks(side)])
            claims = tuned.compare(exhaustive, mode='claims', labels=['x.fa', 'z.fa'])
            self.assertEqual(([], []), (claims.only_in_a, claims.differ))
            rev = exhaustive.compare(tuned, mode='prefix_subset')
            self.assertFalse(rev.equal)
            self.assertEqual(['x.fa', 'y.fa'], sorted(names[x['ref']] for x in rev.only_in_a))

    def test_claims_ending_at_the_depth_are_open_on_both_sides(self):
        r100 = self.g('cli_bubbles_r100')
        for name, depth in (('cli_bubbles_r60', 60), ('cli_bubbles_steps150', 65)):
            for a, b in ((self.g(name), r100), (r100, self.g(name))):
                c = a.compare(b, mode='claims')
                self.assertEqual((True, True, depth, [], [], []),
                                 (c.comparable, c.equal, c.depth_used, c.only_in_a,
                                  c.only_in_b, c.differ), name)

    def test_a_lower_bound_side_is_never_equal(self):
        constrain = parse(_with_fp(T.doc_text('merge')))
        annotate = parse(_with_fp(T.doc_text('annotate')))     # lists cut at 2
        self.assertEqual('lower_bound', annotate.outcome.label_evidence)
        c = constrain.compare(annotate, arm='right', mode='claims', labels=['a.fa'])
        self.assertEqual(('qualified', None), (c.comparable, c.equal))
        self.assertIn('b: label_lists (labels.max_labels_per_node): the oracle is a lower '
                      'bound', c.notes)
        # a qualification never upgrades a weaker verdict
        bare = parse(T.doc_text('annotate'))
        self.assertEqual('unverifiable', constrain.compare(bare, arm='right').comparable)

    def test_no_common_certified_depth_is_unknown(self):
        g = parse(_with_fp(T.doc_text('no_bins')))      # complete_to_bp 0 on both arms
        for mode in ('claims', 'walks', 'labels', 'prefix_subset'):
            c = g.compare(g, mode=mode)
            self.assertEqual(('unknown', None, 0), (c.comparable, c.equal, c.depth_used))
            self.assertIn('no common certified depth (complete_to_bp 0 on a or b)', c.notes)


class TestExports(unittest.TestCase):
    def test_fasta(self):
        g = T.graphlet('fork')
        text = g.to_fasta()
        recs = dict(zip(text.split('\n')[0::2], text.split('\n')[1::2]))
        self.assertEqual(4, len(recs))
        seed = g.seed.sequence
        self.assertEqual('ATCCCTCACG' + seed,
                         recs[[h for h in recs if h.startswith('>left_0')][0]])
        self.assertEqual(seed + 'TTACACCACG',
                         recs[[h for h in recs if h.startswith('>right_1')][0]])
        whole = seed + 'TTACACCACG'
        wrapped = g.to_fasta('right', [1], width=7)
        self.assertEqual([whole[i:i + 7] for i in range(0, len(whole), 7)],
                         wrapped.split('\n')[1:-1])
        self.assertEqual('ATTGCTAAGA\n',
                         g.to_fasta('right', [0], with_seed=False).split('\n', 1)[1])

    def test_gfa_parses_back(self):
        for name in ('fork', 'merge', 'clipped_merge'):
            g = T.graphlet(name)
            lines = [l.split('\t') for l in g.to_gfa().splitlines()]
            seqs = {l[1]: l[2] for l in lines if l[0] == 'S'}
            links = [l for l in lines if l[0] == 'L']
            paths = [l for l in lines if l[0] == 'P']
            k1 = g.k - 1
            self.assertEqual(1 + sum(len(a.segments) for a in g.arms.values()), len(seqs))
            self.assertEqual(sum(len(p) or 1 for a in g.arms.values()
                                 for p in (s.parents for s in a.segments)), len(links))
            for l in links:
                self.assertEqual('%dM' % k1, l[5])
                self.assertEqual(seqs[l[1]][-k1:], seqs[l[3]][:k1], (name, l))
            for p in paths:
                names = [x[:-1] for x in p[2].split(',')]
                spelled = seqs[names[0]] + ''.join(seqs[n][k1:] for n in names[1:])
                side, pid = p[1].split('_')
                expect = g.spell(side, int(pid), with_seed=True)
                self.assertEqual(expect, spelled, (name, p[1]))

    def test_frames_records(self):
        from metagraph.traverse.frames import records
        r = records(T.body('merge'))
        self.assertEqual(7, len(r['segments']))
        self.assertEqual(6, len(r['claims']))
        self.assertEqual({'c:0', 'c:1', 'c:2', 'c:3'}, {x['ref'] for x in r['runs']})


class TestSizeAndSummary(unittest.TestCase):
    def test_memory_and_summary(self):
        for name in T.DOCUMENTS:
            g = T.graphlet(name)
            self.assertGreater(g.memory_bytes(), 1000)
            s = g.summary()
            self.assertLessEqual(len(json.dumps(s, separators=(',', ':')).encode()), 2048)
            self.assertEqual(set(g.arms), set(s['arms']))
            for lab in s['arms'][next(iter(g.arms))]['top_labels']:
                self.assertEqual({'name', 'ref', 'direct_bp', 'reach_bp'}, set(lab))
        s = T.graphlet('limits').summary()
        self.assertIn('resource_stop', {c['kind'] for c in s['caveats']})
        self.assertFalse(s['index']['verifiable'])


if __name__ == '__main__':
    unittest.main()
