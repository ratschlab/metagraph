"""Stitching a continuation (DESIGN-graphlet-merge.md §3, phase 1a-S; Appendix E.2 T1-T9).

The CLI triples (a parent at radius R, its child from next_request(), the uninterrupted
retrieval at R + E) are in data/traverse/documents/stitch/, written by
scripts/traversal/graphlet_fixtures.py (STITCH_FIXTURES)."""

import copy
import gc
import json
import os
import re
import sys
import tempfile
import tracemalloc
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402

import metagraph.traverse as MT  # noqa: E402
from metagraph.traverse import (BadSelector, CombineRefused, Label, LocalBudget,  # noqa: E402
                                LocalBudgetExceeded, MissingEnvelope, StitchedGraphlet,
                                from_response, ops, parse, stitch)
from metagraph.traverse import composite as K  # noqa: E402
from metagraph.traverse import derive  # noqa: E402


def rc(s):
    return s[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def load(name, role, edit=None):
    """A graphlet of triple |name| from its response, |edit|(resp) applied first (the
    envelope edits of W5's tests: H's index_fp, the release, the pair)."""
    with open(os.path.join(T.STITCH, '%s.%s.json' % (name, role)), encoding='utf-8') as f:
        resp = json.load(f)
    if edit is not None:
        edit(resp)
    return from_response(resp['results'][0], resp)


def with_fp(fp):
    """The index manifest digest |fp| in the body's H line (and the capabilities)."""
    def edit(resp):
        r = resp['results'][0]
        lines = r['graphlet'].split('\n')
        h = lines[0].split(' ')
        h[14] = fp
        lines[0] = ' '.join(h)
        r['graphlet'] = '\n'.join(lines)
        r.pop('graphlet_bytes', None)
        r.pop('graphlet_lines', None)
        resp['capabilities']['index_fp'] = fp
    return edit


def with_pair(graph, path):
    def edit(resp):
        resp['graph'] = graph
        resp['graph_path'] = path
    return edit


def triple(name):
    req = T.stitch_requests(name)
    return (req['fixture_index']['arm'], load(name, 'parent'), load(name, 'child'))


def refused(test, code, text, fn, *args, **kw):
    with test.assertRaises(CombineRefused) as cm:
        fn(*args, **kw)
    test.assertEqual(code, cm.exception.code)
    test.assertEqual(text, str(cm.exception))
    return cm.exception


RANK = K.OUTCOME_RANKS


def no_stronger(test, weak, strong, what):
    """Every dimension of outcome |weak| is no stronger than |strong|'s."""
    for dim, ranks in RANK.items():
        if strong is None or strong.get(dim) is None:
            continue
        test.assertGreaterEqual(ranks.get(weak.get(dim), len(ranks)), ranks[strong[dim]],
                                '%s: %s %r stronger than %r' % (what, dim, weak.get(dim),
                                                               strong[dim]))


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


class TestStitchFixtures(unittest.TestCase):
    """The CLI triples: each child ran the request parent.next_request() builds (the
    script builds it without the library, stitch_child_request(); these pin the two
    equal), on the parent's index, seeded with the continued walk's continuation."""

    def test_child_requests_are_next_request(self):
        fixtures = T.fixture_script().STITCH_FIXTURES
        self.assertTrue(fixtures)
        for name, fx in fixtures.items():
            with self.subTest(name):
                req = T.stitch_requests(name)
                parent = T.stitch_graphlet(name, 'parent')
                walk = req['fixture_index']['walk']
                nr = parent.next_request(fx['arm'], [walk], fx['bp'],
                                         reset_branches=fx.get('reset_branches', False),
                                         **(fx.get('overrides') or {}))
                self.assertEqual(req['child'], dict(nr), T.json_diff(req['child'], dict(nr)))
                if 'sibling' in req:
                    nr = parent.next_request(fx['arm'], [fx['sibling']], fx['bp'])
                    self.assertEqual(req['sibling'], dict(nr))
                if 'grandchild' in req:
                    child = T.stitch_graphlet(name, 'child')
                    gc_ = fx['grandchild']
                    nr = child.next_request(fx['arm'], [gc_['walk']], gc_['bp'])
                    self.assertEqual(req['grandchild'], dict(nr))

    def test_children_are_seeded_with_the_continuation(self):
        for name in T.stitch_names():
            with self.subTest(name):
                req = T.stitch_requests(name)
                arm, walk = req['fixture_index']['arm'], req['fixture_index']['walk']
                parent = T.stitch_graphlet(name, 'parent')
                child = T.stitch_graphlet(name, 'child')
                self.assertEqual(parent.continuation(arm, walk).sequence, child.seed.sequence)
                self.assertEqual([arm], list(child.arms))


# ---------------------------------------------------------------------------------- T1

class TestRefusals(unittest.TestCase):
    """T1: one test per code of §3.1, in the order checked, with Appendix B.1's text."""

    def setUp(self):
        self.arm, self.P, self.C = triple('keep_forbid')

    def test_not_combinable(self):
        with self.assertRaises(TypeError) as cm:
            stitch(self.P, 'GATTACA')
        self.assertEqual('not_combinable', cm.exception.code)
        self.assertEqual('stitch() combines graphlets (Graphlet or StitchedGraphlet), not str',
                         str(cm.exception))
        with self.assertRaises(TypeError):
            stitch(None)

    def test_view_unsupported(self):
        refused(self, 'view_unsupported', 'stitch() reads whole retrievals, not views: stitch '
                'view.backing instead', stitch, self.P.subgraph([self.P.labels[0].ref]),
                self.C)
        refused(self, 'view_unsupported', 'stitch() reads whole retrievals, not views: stitch '
                'view.backing instead', stitch, self.P, self.C.subgraph([self.C.labels[0].ref]))

    def test_missing_envelope(self):
        text = '\n'.join(l for l in self.C.dump().splitlines() if not l.startswith('J')) + '\n'
        body = parse(text)
        with self.assertRaises(MissingEnvelope) as cm:
            stitch(self.P, body)
        self.assertEqual("stitch() needs the response envelope (seed_id, the normalized "
                         "strategy, annotation counters): build the graphlet with "
                         "Graphlet.from_response() or load a saved .mgt file with its J line; "
                         "a parsed body alone does not carry it", str(cm.exception))

    def _index(self, clause, c):
        refused(self, 'different_index', 'the child was retrieved from another index than its '
                'parent (%s): a continuation runs on its parent\'s index, so a child of another '
                'index continues none of its walks; union_pairs() combines one seed\'s '
                'graphlets across pairs' % clause, stitch, self.P2, c)

    def test_different_index(self):
        self.P2 = load('keep_forbid', 'parent', with_fp('a' * 16))
        self._index('index_fp differs', load('keep_forbid', 'child', with_fp('b' * 16)))

        def rel(name):
            def edit(resp):
                resp['release'] = name
                resp['capabilities']['release'] = name
            return edit
        self.P2 = load('keep_forbid', 'parent', rel('r1'))
        self._index("release 'r1' vs 'r2'", load('keep_forbid', 'child', rel('r2')))
        self.P2 = self.P
        # a ref named otherwise (compare()'s names rule), also against a StitchedGraphlet
        c = copy.copy(self.C)
        lab = c.labels[0]
        c.labels = [Label(lab.id, lab.kind, lab.column, lab.seq_id, 'other')] + \
            list(c.labels[1:])
        self._index("label %s is %r in one and 'other' in the other" % (lab.ref, lab.name), c)
        self.P2 = stitch(self.P)
        self._index("label %s is %r in one and 'other' in the other" % (lab.ref, lab.name), c)

    def test_index_unverifiable(self):
        text = ('the child states %s index manifest digest (index_fp) and its parent %s: '
                'whether they are the same index cannot be verified (not proven different), '
                'so the child is not stitched automatically; pass allow_unverified_index=True '
                'to stitch it anyway (the result then states the identity unverified)')
        p = load('keep_forbid', 'parent', with_fp('a' * 16))
        refused(self, 'index_unverifiable', text % ('no', 'one'), stitch, p, self.C)
        c = load('keep_forbid', 'child', with_fp('a' * 16))
        refused(self, 'index_unverifiable', text % ('an', 'none'), stitch, self.P, c)
        sg = stitch(p, self.C, allow_unverified_index=True)
        self.assertEqual(False, sg.identity['verified'])
        self.assertEqual('allow_unverified_index', sg.identity['accepted'])
        self.assertEqual(sg.identity, sg.evidence()['identity'])
        # both stated and equal: verified; neither: unverified, not refused (OQ6)
        sg = stitch(p, c)
        self.assertEqual({'verified': True, 'index_fp': 'a' * 16}, sg.identity)
        self.assertEqual(False, stitch(self.P, self.C).identity['verified'])
        with self.assertRaises(ValueError):
            stitch(self.P, self.C, allow_unverified_index=1)

    def test_different_pair(self):
        p = load('keep_forbid', 'parent', with_pair('g1', '/p/g1.dbg'))
        c = load('keep_forbid', 'child', with_pair('g2', '/p/g2.dbg'))
        refused(self, 'different_pair', 'the child was answered on pair '
                '{"graph":"g2","graph_path":"/p/g2.dbg","annotation_path":null}, its parent on '
                'pair {"graph":"g1","graph_path":"/p/g1.dbg","annotation_path":null}: a '
                'continuation runs on its parent\'s pair', stitch, p, c)
        # one pair stated, or the same pair: not refused
        stitch(p, self.C)
        stitch(p, load('keep_forbid', 'child', with_pair('g1', '/p/g1.dbg')))
        # with a verified identity the pair is not checked (the fingerprint decided)
        stitch(load('keep_forbid', 'parent', lambda r: (with_fp('a' * 16)(r),
                                                         with_pair('g1', None)(r))),
               load('keep_forbid', 'child', lambda r: (with_fp('a' * 16)(r),
                                                        with_pair('g2', None)(r))))

    def test_different_graph_and_mode(self):
        c = copy.copy(self.C)
        c.k = 31
        refused(self, 'different_graph', 'the child\'s graph is not its parent\'s (k 31 vs 15): '
                'its walks cannot continue the parent\'s', stitch, self.P, c)
        c = copy.copy(self.C)
        c.regime = 'canonical'
        refused(self, 'different_graph', 'the child\'s graph is not its parent\'s (regime '
                '\'canonical\' vs \'basic\'): its walks cannot continue the parent\'s',
                stitch, self.P, c)
        c = copy.copy(self.C)
        c.support = 'trace'
        refused(self, 'different_mode', 'the child was walked in constrain mode with trace '
                'support, its parent in constrain mode with kmer support: its walks answer '
                'another question', stitch, self.P, c)

    def test_arms(self):
        _, fp, fc = triple('fork_right')
        refused(self, 'not_one_arm', 'the child walked both arms: a continuation walks only the '
                'arm it continues (next_request() sets direction to it)', stitch, fc, fp)
        refused(self, 'wrong_arm', 'the child walked the right arm, not the left arm named by '
                'arm=', stitch, fp, fc, arm='l')
        with self.assertRaises(ValueError):
            stitch(fp, fc, arm='up')

    def test_continuation_codes(self):
        _, mp, mc = triple('merge_only')
        refused(self, 'no_continuation', 'walk 0 of the right arm has no continuation (it ended '
                'for a semantic reason, not at the radius or a cap)', stitch, mp, mc, walk=0)
        sib = load('keep_forbid', 'sibling')
        refused(self, 'not_a_continuation', 'the child\'s seed is not the continuation of walk 1 '
                'of the right arm: a stitched child is seeded with that walk\'s continuation '
                '(continuation(), next_request())', stitch, self.P, self.C, walk=1)
        refused(self, 'not_a_continuation', 'the child\'s seed is the continuation of no walk of '
                'the right arm: a stitched child is seeded with a walk\'s continuation '
                '(continuation(), next_request())', stitch, self.C, sib)
        seed = self.C.seed.sequence
        with mock.patch.object(derive, 'continuation_sequence', lambda g, a, leaf, w=None: seed):
            refused(self, 'ambiguous_walk', 'the child\'s seed is the continuation of walks 0, 1 '
                    'of the right arm: name the one it continues with walk=', stitch,
                    self.P, self.C)
        with self.assertRaises(BadSelector):
            stitch(self.P, self.C, walk='0')
        with self.assertRaises(IndexError):
            stitch(self.P, self.C, walk=7)

    def test_already_continued(self):
        _, lp, lc = triple('loop_branch')
        cut0 = load('loop_branch_cut0', 'child')
        sg = stitch(lp, lc)
        refused(self, 'already_continued', 'walk 1 of the right arm is already continued by '
                'another child (part right:1, digest %s): one child per walk'
                % K.part_digest(lc), stitch, sg, cut0)

        def later(resp):
            resp['strategy']['bounds']['time_budget_ms'] = 12345
        again = load('loop_branch', 'child', later)
        self.assertEqual(K.part_digest(lc), K.part_digest(again))
        refused(self, 'already_continued', 'walk 1 of the right arm is continued by this body '
                'from another request (strategy.bounds.time_budget_ms differs): equal bodies '
                'are not equal retrievals, and the two would continue differently',
                stitch, sg, again)
        # the identical retrieval again: a no-op, the stitched graphlet itself
        self.assertIs(sg, stitch(sg, load('loop_branch', 'child')))

    def test_foreign_seed_labels(self):
        c = copy.copy(self.C)
        n = c.seed.num_seed_labels
        self.assertEqual(1, n)
        c.labels = [Label(0, 'header', 0, 99, 'zz')] + list(c.labels[1:])
        refused(self, 'foreign_seed_labels', 'the child seeds {"name": "zz", "ref": "h:0:99"}, '
                'not alive at the end of walk 0 of the right arm: their runs would start new '
                'lineages at 20 bp, not continue the walk\'s', stitch, self.P, c)

    def test_not_offered(self):
        sg = stitch(self.P, self.C)
        texts = {
            'dump': 'a stitched graphlet is not an MGT document: save() writes its parts '
                    'unchanged in a bundle, and each part dumps on its own '
                    '(parts[i].graphlet.dump())',
            'to_json': 'to_json() is detail "full" of one server result, and a stitched '
                       'graphlet has none: export its parts',
            'compare': 'compare() reads one retrieval per side: compare the parts'}
        for name in ('dump', 'to_json', 'compare', 'subgraph', 'to_gfa', 'label_walks',
                     'routes', 'label_summary', 'splits', 'claims', 'support_profile',
                     'support_changes', 'save'):
            with self.subTest(name):
                refused(self, 'not_offered',
                        texts.get(name, '%s() is not offered on a stitched graphlet' % name),
                        getattr(sg, name), 'x')

    def test_label_selectors(self):
        sg = stitch(self.P, self.C)
        with self.assertRaises(BadSelector) as cm:
            sg.label({'id': 0})
        self.assertEqual('label ids are local to each part of a stitched graphlet: select by '
                         'ref or name', str(cm.exception))
        self.assertEqual('h:0:0', sg.label('acc1').ref)
        self.assertEqual('acc1', sg.label({'ref': 'h:0:0'}).name)
        self.assertEqual([l.ref for l in self.P.labels], [l.ref for l in sg.labels])


# ---------------------------------------------------------------------------------- T2

ROLES = ('parent', 'child', 'reference', 'sibling', 'grandchild')
_SUFFIX = re.compile(r' walk=\d+ junctions=\. verdict=\. certified_to=\d+$', re.M)


def every_graphlet():
    for name in T.stitch_names():
        req = T.stitch_requests(name)
        for role in ROLES:
            if role in req or role in ('parent', 'child'):
                yield '%s.%s' % (name, role), load(name, role)


class TestIdentityElement(unittest.TestCase):
    """T2: stitch(g) answers every offered query as g, with keys (w,) for path ids w."""

    def test_stitch_of_one_answers_as_it(self):
        for what, g in every_graphlet():
            with self.subTest(what):
                self._same(g, stitch(g))

    def _same(self, g, sg):
        self.assertEqual((), sg.junctions)
        self.assertEqual(1, len(sg.parts))
        self.assertEqual([l.ref for l in g.labels], [l.ref for l in sg.labels])
        self.assertEqual([s for s in ('left', 'right') if s in g.arms], sg.arms())
        has_bases = K.has_bases(g)
        for side in sg.arms():
            for by in ('support', 'length', 'loss', 'id'):
                gw = g.walks(side, by=by)
                sw = sg.walks(side, by=by)
                self.assertEqual([(w.path_id,) for w in gw], [w.key for w in sw])
                for a, b in zip(gw, sw):
                    self.assertEqual((a.length_bp, a.sequence, a.path_reason, a.end_reasons,
                                      [l.ref for l in a.labels_full], a.n_alive, a.complete,
                                      a.beyond_certified_bp),
                                     (b.length_bp, b.sequence, b.path_reason, b.end_reasons,
                                      [l.ref for l in b.labels_full], b.n_alive, b.complete,
                                      b.beyond_certified_bp))
                    self.assertIsNone(b.verdict)
                    self.assertEqual((), b.causes)
                    self.assertEqual((), b.junctions)
                    self.assertEqual(b.length_bp - b.beyond_certified_bp, b.certified_to_bp)
            for w in g.walks(side, by='id'):
                if has_bases:
                    for o in ('natural', 'walk'):
                        for ws in (False, True):
                            self.assertEqual(g.spell(side, w.path_id, o, ws),
                                             sg.spell(side, (w.path_id,), o, ws))
                leaf = derive.paths(g.arms[side])[w.path_id].leaf
                if g.arms[side].segments[leaf].leaf.continuation is not None and has_bases:
                    a, b = g.continuation(side, w.path_id), sg.continuation(side, w.path_id)
                    self.assertEqual((a.sequence, [l.ref for l in a.labels], a.seed_coord),
                                     (b.sequence, [l.ref for l in b.labels], b.seed_coord))
                    self.assertEqual(dict(g.next_request(side, [w.path_id])),
                                     dict(sg.next_request(side, [(w.path_id,)])))
            comp = sg.completeness()[side]
            self.assertEqual((g.arms[side].complete_to_bp, g.arms[side].scope, []),
                             (comp['complete_to_bp'], comp['scope'], comp['extensions']))
        if has_bases:
            for o in ('natural', 'walk'):
                self.assertEqual(g.to_fasta(orientation=o, width=7, coordinates=False),
                                 _SUFFIX.sub('', sg.to_fasta(orientation=o, width=7)))
        ev = sg.evidence()
        blk = ops.evidence_block(g, None)
        self.assertEqual({'root': blk}, ev['parts'])
        self.assertEqual((blk['outcome'], bool(blk['exact'])), (ev['outcome'], ev['exact']))
        self.assertEqual([], ev['junctions'])
        self.assertEqual(sg.digest, stitch(g).digest)


# ---------------------------------------------------------------------------------- T3

class TestAlgebra(unittest.TestCase):
    """T3 (Appendix C.3)."""

    def test_commutative_across_junctions(self):
        arm, p, c = triple('keep_forbid')
        s = load('keep_forbid', 'sibling')
        a = stitch(stitch(p, c), s)
        b = stitch(stitch(p, s), c)
        self.assertEqual(a.digest, b.digest)
        self.assertEqual([(), (('right', 0),), (('right', 1),)], [x.key for x in a.parts])
        self.assertEqual([w.key for w in a.walks(arm, by='id')],
                         [w.key for w in b.walks(arm, by='id')])
        self.assertEqual([(0, 0), (1, 0)], [w.key for w in a.walks(arm, by='id')])
        self.assertEqual(a.to_fasta(), b.to_fasta())
        # inference over continued and uncontinued candidates: one answer in every order
        self.assertEqual(a.digest, stitch(stitch(p, c, walk=0), s, walk=1).digest)
        self.assertEqual(a.digest, stitch(stitch(p, s, walk=(1,)), c, walk=(0,)).digest)

    def test_associative_along_chains(self):
        arm, p, c = triple('chain3')
        g3 = load('chain3', 'grandchild')
        a = stitch(stitch(p, c, walk=0), g3, walk=(0, 0))
        b = stitch(p, stitch(c, g3, walk=0), walk=0)
        self.assertEqual(a.digest, b.digest)
        self.assertEqual(a.digest, stitch(stitch(p, c), g3).digest)
        self.assertEqual([(), (('right', 0),), (('right', 0), ('right', 0))],
                         [x.key for x in a.parts])
        self.assertEqual([0, 15, 30], [x.offset_bp for x in a.parts])
        w = a.walks(arm, by='id')
        self.assertEqual([(0, 0, 0)], [x.key for x in w])
        self.assertEqual(45, w[0].length_bp)
        # equal to the uninterrupted retrieval's walk (the chain is gated twice)
        ref = load('chain3', 'reference')
        self.assertEqual(ref.spell(arm, 0), a.spell(arm, (0, 0, 0)))
        self.assertEqual(p.spell(arm, 0) + c.spell(arm, 0) + g3.spell(arm, 0),
                         a.spell(arm, (0, 0, 0)))
        self.assertEqual(2, len(a.junctions))
        self.assertEqual([15, 30], [j.at_bp for j in a.junctions])
        self.assertEqual([15, 15], [j.overlap_bp for j in a.junctions])
        row = a._prospective(arm, (0, 0, 0))
        self.assertEqual((45, True, 'exact'), (row['certified_to_bp'], row['complete'],
                                                row['verdict']))
        self.assertEqual((15, False, 'qualified'), (w[0].certified_to_bp, w[0].complete,
                                                     w[0].verdict))
        # next_request() and continuation() of a continued walk inside the chain, from
        # the part that owns it, the seed coordinates mapped to the root seed
        self.assertEqual(dict(c.next_request(arm, [0], 30)),
                         dict(a.next_request(arm, [(0, 0)], 30)))
        cont = a.continuation(arm, (0, 0))
        n = len(p.seed.sequence)
        self.assertEqual(g3.seed.sequence, cont.sequence)
        lo, hi = cont.seed_coord
        self.assertEqual(cont.sequence, (p.seed.sequence + a.spell(arm, (0, 0, 0)))[lo:hi])
        self.assertEqual(n + 30, hi)

    def test_parent_and_child_do_not_commute(self):
        _, fp, fc = triple('fork_right')
        with self.assertRaises(CombineRefused) as cm:
            stitch(fc, fp)
        self.assertEqual('not_one_arm', cm.exception.code)
        _, p, c = triple('keep_forbid')
        with self.assertRaises(CombineRefused) as cm:
            stitch(c, p)
        self.assertEqual('not_a_continuation', cm.exception.code)

    def test_circular_both_orders(self):
        arm, p, c = triple('circular')
        a, b = stitch(p, c), stitch(c, p)
        self.assertNotEqual(a.digest, b.digest)
        for sg in (a, b):
            self.assertEqual(['qualified'], [j.verdict for j in sg.junctions])
            for w in sg.walks(arm, by='id'):
                self.assertEqual('qualified', w.verdict)
                self.assertFalse(w.complete)
        # W8 qualifies the walk that re-enters the root seed
        self.assertEqual([('rule_violation', 49)],
                         [(x[0], x[2]) for x in a.walks(arm, by='id')[0].causes])

    def test_idempotent_and_canonical(self):
        arm, p, c = triple('keep_forbid')
        sg = stitch(p, c)
        self.assertIs(sg, stitch(sg, c))
        self.assertIs(sg, sg.stitch(c))
        self.assertEqual(sg.digest, stitch(p, load('keep_forbid', 'child')).digest)
        self.assertNotEqual(sg.digest, stitch(p).digest)
        self.assertIsInstance(stitch(sg), StitchedGraphlet)
        self.assertIs(sg, stitch(sg))
        self.assertIsNone(sg._cache.get('never'))


# ---------------------------------------------------------------------------------- T4

# per triple: the junction's causes; per stitched walk (key): (length_bp, certified_to_bp,
# complete, the walk-scope causes (code, at_bp), the prospective certified_to_bp and
# complete); and the prospective evidence's outcome (walks, branch_diagnostics,
# label_evidence, delivery). The public outcome is the parts' weakest with walks partial
# and label_evidence qualified.
EXPECT = {
    'keep_forbid': (('exact_gated',), {(0, 0): (40, 20, False, (), 40, True),
                                       (1,): (20, 20, True, (), 20, True)},
                    ('complete', 'complete', 'complete', 'inline')),
    'merge_only': (('exact_gated',), {(0,): (10, 10, True, (), 10, True),
                                      (1, 0): (60, 30, False, (), 60, True)},
                   ('complete', 'complete', 'complete', 'inline')),
    'unlimited': (('exact_gated',), {(0,): (10, 10, True, (), 10, True),
                                     (1, 0): (60, 30, False, (), 60, True)},
                  ('complete', 'complete', 'complete', 'inline')),
    'equal_finite': (('exact_gated',), {(0,): (10, 10, True, (), 10, True),
                                        (1, 0): (60, 30, False, (), 60, True)},
                     ('complete', 'complete', 'complete', 'inline')),
    'chain3': (('exact_gated',), {(0, 0): (30, 15, False, (), 30, True)},
               ('complete', 'complete', 'complete', 'inline')),
    'fork_right': (('sibling_merge_unknown',), {(0, 0): (40, 10, False, (), 10, False),
                                                (1,): (10, 10, True, (), 10, True)},
                   ('partial', 'complete', 'qualified', 'inline')),
    'fork_left': (('sibling_merge_unknown',), {(0, 0): (40, 10, False, (), 10, False),
                                               (1,): (10, 10, True, (), 10, True)},
                  ('partial', 'complete', 'qualified', 'inline')),
    'sibling_merge': (('sibling_merge_unknown',), {(0, 0): (100, 30, False, (), 30, False),
                                                   (1,): (30, 30, True, (), 30, True)},
                      ('partial', 'complete', 'qualified', 'inline')),
    'trace': (('trace_restart',), {(0, 0): (40, 20, False, (), 20, False),
                                   (1,): (20, 20, True, (), 20, True)},
              ('partial', 'complete', 'qualified', 'inline')),
    'quorum_changed': (('strategy_differs', 'label_count_knob'),
                       {(0, 0): (20, 20, False, (), 20, False),
                        (1,): (20, 20, True, (), 20, True)},
                       ('partial', 'complete', 'qualified', 'inline')),
    'reset_branches': (('branches_unequal',), {(0,): (10, 10, True, (), 10, True),
                                               (1, 0): (60, 30, False, (), 30, False)},
                       ('partial', 'complete', 'qualified', 'inline')),
    'switching': (('switching',), {(0, 0): (90, 50, False, (), 50, False),
                                   (0, 1): (90, 50, False, (), 50, False)},
                  ('partial', 'complete', 'qualified', 'inline')),
    'extra_foreign': (('extra_foreign', 'switching'), {(0, 0): (90, 50, False, (), 50, False),
                                                       (0, 1): (90, 50, False, (), 50, False)},
                      ('partial', 'complete', 'qualified', 'inline')),
    'not_continued': (('switching', 'not_continued'), {(0, 0): (60, 25, False, (), 25, False)},
                      ('partial', 'complete', 'qualified', 'inline')),
    # Appendix G.4: the child stops at 0 bp, the reference goes on to 145 bp
    'loop_branch': (('exact_gated',), {(0,): (85, 85, True, (), 85, True),
                                       (1, 0): (115, 115, False, (('candidate_changed', 115),),
                                                115, False)},
                    ('partial', 'complete', 'qualified', 'inline')),
    'loop_branch_cut0': (('exact_gated',), {(0,): (85, 85, True, (), 85, True),
                                            (1, 0): (115, 115, False,
                                                     (('decisions_unrecorded', 115),),
                                                     115, False)},
                         ('partial', 'cut', 'qualified', 'inline')),
    'loop_same_stop': (('exact_gated',), {(0,): (85, 85, True, (), 85, True),
                                          (1, 0): (145, 115, False, (('seed_reentry', 115),),
                                                   115, False)},
                       ('partial', 'complete', 'qualified', 'inline')),
    # the changed decision inside the child: the prospective result stops there (115), the
    # public one at the junction (90)
    'loop_interior': (('exact_gated',), {(0,): (85, 85, True, (), 85, True),
                                         (1, 0): (145, 90, False, (('candidate_changed', 115),),
                                                  115, False),
                                         (1, 1): (149, 90, False, (('rule_violation', 115),
                                                                   ('seed_reentry', 149)),
                                                  115, False)},
                      ('partial', 'complete', 'qualified', 'inline')),
    'earlier_stop': (('exact_gated',), {(0,): (66, 66, True, (), 66, True),
                                        (1, 0): (100, 70, False, (('seed_reentry', 100),),
                                                 100, False)},
                     ('partial', 'complete', 'qualified', 'inline')),
    # Appendix G.5: the reference stops at 50 bp with split_limit
    'splits45': (('exact_gated',), {(0,): (18, 18, True, (), 18, True),
                                    (1, 0): (85, 45, False, (('split_limit_restarted', 50),),
                                             50, False),
                                    (1, 1): (85, 45, False, (('split_limit_restarted', 50),),
                                             50, False)},
                 ('partial', 'complete', 'qualified', 'inline')),
    'events_cut': (('exact_gated',), {(0, 0): (50, 10, False, (), 50, True),
                                      (0, 1): (90, 10, False, (('decisions_unrecorded', 60),),
                                               60, False),
                                      (0, 2): (100, 10, False, (('decisions_unrecorded', 60),),
                                               60, False)},
                   ('partial', 'cut', 'qualified', 'inline')),
    'events_at_end': (('exact_gated',), {(0, 0): (60, 10, False,
                                                  (('decisions_unrecorded', 60),), 60, False),
                                         (0, 1): (50, 10, False, (), 50, True)},
                      ('partial', 'cut', 'qualified', 'inline')),
    'circular': (('exact_gated',), {(0, 0): (50, 25, False, (('rule_violation', 49),),
                                             49, False)},
                 ('partial', 'complete', 'qualified', 'inline')),
    'no_bases': (('no_bases',), {(0, 0): (40, 20, False, (), 20, False),
                                 (1,): (20, 20, True, (), 20, True)},
                 ('partial', 'complete', 'qualified', 'inline')),
}

SEED_REENTRY = {'loop_same_stop': 'same_stop', 'earlier_stop': 'earlier_stop',
                'loop_interior': 'same_stop'}


class TestOracleTriples(unittest.TestCase):
    """T4: per CLI triple, the junction's verdict and causes, every stitched walk's
    certified_to_bp, complete and causes, and the evidence, public and prospective."""

    def test_every_triple_is_pinned(self):
        self.assertEqual(sorted(T.stitch_names()), sorted(EXPECT))

    def test_triples(self):
        for name, (jcauses, walks, prospective_outcome) in EXPECT.items():
            with self.subTest(name):
                arm, p, c = triple(name)
                sg = stitch(p, c)
                self.assertEqual(1, len(sg.junctions))
                j = sg.junctions[0]
                self.assertEqual(('qualified', jcauses), (j.verdict, tuple(x[0] for x in j.causes)))
                self.assertEqual(len(c.seed.sequence), j.overlap_bp)
                got = {}
                for w in sg.walks(arm, by='id'):
                    pr = sg._prospective(arm, w.key)
                    own = tuple((x[0], x[2]) for x in w.causes if x[0] != 'exact_gated')
                    got[w.key] = (w.length_bp, w.certified_to_bp, w.complete, own,
                                  pr['certified_to_bp'], pr['complete'])
                    self.assertEqual(w.length_bp - w.certified_to_bp, w.beyond_certified_bp)
                    if len(w.key) == 1:
                        self.assertIsNone(w.verdict)
                        continue
                    # every verdict qualified; exact_gated on a walk only where its junction
                    # and the walk itself have nothing else to state
                    self.assertEqual('qualified', w.verdict)
                    self.assertFalse(w.complete)
                    gated = jcauses == ('exact_gated',) and not own
                    self.assertEqual(gated, [x[0] for x in w.causes] == ['exact_gated'])
                    self.assertEqual('exact' if gated else 'qualified', pr['verdict'])
                    self.assertLessEqual(w.certified_to_bp, j.at_bp)
                    for code, sentence, at in w.causes:
                        if code == 'seed_reentry':
                            self.assertTrue(sentence.startswith(SEED_REENTRY[name] + ': '),
                                            sentence)
                        if at is not None:
                            self.assertRegex(sentence, r'\b%d bp' % at)
                self.assertEqual(walks, got)
                # the public evidence: a qualified junction is read, so walks partial and
                # label_evidence qualified, and never exact
                ev = sg.evidence(arm)
                parts = [ops.evidence_block(x.graphlet, arm) for x in sg.parts]
                weakest = K.weakest_outcome([b['outcome'] for b in parts])
                weakest.update(walks=weakest['walks'] if weakest['walks'] == 'failed'
                               else 'partial', label_evidence='qualified')
                self.assertEqual(weakest, ev['outcome'])
                self.assertFalse(ev['exact'])
                pe = sg._prospective(arm)['evidence']
                self.assertEqual(prospective_outcome, tuple(pe['outcome'][d] for d in RANK))
                # the CLI fixtures carry no index_fp: never exact where a junction is read
                self.assertFalse(pe['exact'])
                comp = sg.completeness()[arm]['extensions'][0]
                pcomp = sg._prospective(arm)['completeness']['extensions'][0]
                self.assertEqual(j.at_bp, comp['certified_to_bp'])
                self.assertEqual(list(jcauses), comp['causes'])
                self.assertGreaterEqual(pcomp['certified_to_bp'], comp['certified_to_bp'])
                if jcauses != ('exact_gated',):
                    self.assertEqual(j.at_bp, pcomp['certified_to_bp'])

    def test_differences_lie_below_a_stated_cause(self):
        """Every stitched walk the rules would call exact (prospective) is a walk of the
        uninterrupted retrieval, base for base."""
        for name in T.stitch_names():
            req = T.stitch_requests(name)
            if 'reference' not in req:
                continue
            arm, p, c = triple(name)
            if not K.has_bases(c):
                continue
            with self.subTest(name):
                ref = load(name, 'reference')
                spelled = [ref.spell(arm, w.path_id) for w in ref.walks(arm, by='id')]
                sg = stitch(p, c)
                for w in sg.walks(arm, by='id'):
                    if sg._prospective(arm, w.key)['verdict'] != 'exact':
                        continue
                    # a walk ended by the child's radius goes on in the reference (its
                    # radius is larger); any other ends where the reference's ends
                    if w.path_reason == 'max_extension_bp':
                        self.assertTrue(any(x.startswith(w.sequence) for x in spelled), w.key)
                    else:
                        self.assertIn(w.sequence, spelled, w.key)

    def test_index_fp_and_the_root_walks(self):
        """With an index_fp on both sides: the prospective evidence is exact where the parts
        are, the public one is not; a root walk no child continues keeps the root's."""
        arm = 'right'
        p = load('keep_forbid', 'parent', with_fp('a' * 16))
        c = load('keep_forbid', 'child', with_fp('a' * 16))
        sg = stitch(p, c)
        parts_exact = all(ops.evidence_block(x.graphlet, arm)['exact'] for x in sg.parts)
        self.assertTrue(parts_exact)
        self.assertTrue(sg._prospective(arm)['evidence']['exact'])
        self.assertFalse(sg.evidence(arm)['exact'])
        root = sg.evidence(arm, walks=[(1,)])
        blk = ops.evidence_block(p, arm)
        self.assertEqual({'root': blk}, root['parts'])
        self.assertEqual((blk['outcome'], True), (root['outcome'], root['exact']))
        self.assertEqual([], root['junctions'])
        one = sg.evidence(arm, walks=[(0, 0)])
        self.assertEqual(['root', 'right:0'], list(one['parts']))
        self.assertFalse(one['exact'])
        self.assertEqual([[c[0] for c in j['causes']] for j in one['junctions']],
                         [['exact_gated']])


# ---------------------------------------------------------------------------------- T5

def every_stitch():
    for name in T.stitch_names():
        arm, p, c = triple(name)
        yield name, arm, stitch(p, c)
        req = T.stitch_requests(name)
        if 'sibling' in req:
            yield name + '+sibling', arm, stitch(stitch(p, c), load(name, 'sibling'))
        if 'grandchild' in req:
            yield name + '+grandchild', arm, stitch(stitch(p, c), load(name, 'grandchild'))


class TestNoLaundering(unittest.TestCase):
    """T5: no answer is stronger than a part or junction it reads; the public result only
    reduces the prospective one."""

    def test_no_laundering(self):
        for name, arm, sg in every_stitch():
            with self.subTest(name):
                ev = sg.evidence(arm)
                for x in sg.parts:
                    blk = ops.evidence_block(x.graphlet, arm)
                    key = 'root' if not x.key else '%s:%s' % (
                        x.key[0][0], '.'.join(str(w) for _, w in x.key))
                    self.assertEqual(blk, ev['parts'][key])
                    no_stronger(self, ev['outcome'], blk['outcome'], key)
                self.assertFalse(ev['exact'])
                self.assertEqual(len(sg.junctions), len(ev['junctions']))
                pe = sg._prospective(arm)['evidence']
                no_stronger(self, ev['outcome'], pe['outcome'], 'prospective')
                for w in sg.walks(arm, by='id'):
                    pr = sg._prospective(arm, w.key)
                    self.assertLessEqual(w.certified_to_bp, pr['certified_to_bp'])
                    self.assertTrue(pr['complete'] or not w.complete)
                    if w.complete:
                        self.assertEqual(1, len(w.pieces))
                    if len(w.pieces) > 1:
                        self.assertEqual('qualified', w.verdict)
                        self.assertFalse(w.complete)
                        self.assertLessEqual(w.certified_to_bp,
                                             sg.junctions[w.junctions[0]].at_bp)
                        wev = sg.evidence(arm, walks=[w.key])
                        self.assertFalse(wev['exact'])
                        self.assertEqual(('partial', 'qualified'),
                                         (wev['outcome']['walks'],
                                          wev['outcome']['label_evidence']))
                    else:
                        root = sg.evidence(arm, walks=[w.key])
                        blk = ops.evidence_block(sg.root, arm)
                        self.assertEqual((blk['outcome'], bool(blk['exact'])),
                                         (root['outcome'], root['exact']))
                    for code, sentence, at in w.causes:
                        if code != 'exact_gated':
                            self.assertIn({'walk': w.key, 'code': code, 'at_bp': at,
                                           'sentence': sentence}, ev['walk_causes'])

    def test_summary(self):
        for name, arm, sg in every_stitch():
            with self.subTest(name):
                s = sg.summary()
                self.assertTrue(s['experimental'])
                self.assertFalse(s['exact'])
                self.assertEqual(len(sg.junctions), len(s['junctions']))
                self.assertLessEqual(len(json.dumps(s, separators=(',', ':'))), 2048)
                small = sg.summary(max_bytes=300)
                self.assertEqual([], small['junctions'])
                self.assertTrue(small['junctions_cut'])


# ---------------------------------------------------------------------------------- T6

class TestRoundTrip(unittest.TestCase):
    """T6 (1a): the parts saved and loaded stitch to an equal digest; provenance."""

    def test_saved_parts_stitch_again(self):
        arm, p, c = triple('chain3')
        g3 = load('chain3', 'grandchild')
        sg = stitch(stitch(p, c), g3)
        with tempfile.TemporaryDirectory() as d:
            loaded = []
            for i, x in enumerate(sg.parts):
                path = os.path.join(d, '%d.mgt' % i)
                x.graphlet.save(path)
                loaded.append(MT.load(path))
            again = stitch(stitch(loaded[0], loaded[1]), loaded[2])
        self.assertEqual(sg.digest, again.digest)
        prov = sg.provenance()
        self.assertEqual(['root', 'child', 'child'], [x['role'] for x in prov])
        self.assertEqual([None, {'arm': 'right', 'walk': (0,), 'at_bp': 15, 'overlap_bp': 15},
                          {'arm': 'right', 'walk': (0, 0), 'at_bp': 30, 'overlap_bp': 15}],
                         [x['junction'] for x in prov])
        self.assertEqual([x.digest for x in sg.parts], [x['digest'] for x in prov])
        req = T.stitch_requests('chain3')
        self.assertEqual(req['child']['seeds'][0]['sequence'],
                         prov[1]['request']['seeds'][0]['sequence'])


# ---------------------------------------------------------------------------------- T7

class TestExports(unittest.TestCase):
    """T7: stitched FASTA is the parts' spellings joined at L without the overlap."""

    def test_fasta_both_arms(self):
        for name in ('fork_right', 'fork_left', 'keep_forbid', 'switching'):
            arm, p, c = triple(name)
            walk = T.stitch_requests(name)['fixture_index']['walk']
            sg = stitch(p, c)
            with self.subTest(name):
                j = sg.junctions[0]
                if name.startswith('fork'):
                    self.assertLess(j.at_bp, j.overlap_bp)       # L < n
                seed = p.seed.sequence
                for w in c.walks(arm, by='id'):
                    key = (walk, w.path_id)
                    pw, cw = p.spell(arm, walk, 'walk'), c.spell(arm, w.path_id, 'walk')
                    self.assertEqual(pw + cw, sg.spell(arm, key, 'walk'))
                    nat = (p.spell(arm, walk) + c.spell(arm, w.path_id) if arm == 'right'
                           else c.spell(arm, w.path_id) + p.spell(arm, walk))
                    self.assertEqual(nat, sg.spell(arm, key))
                    self.assertEqual(seed + nat if arm == 'right' else nat + seed,
                                     sg.spell(arm, key, with_seed=True))
                    fa = sg.to_fasta(arm, walks=[key])
                    head, body = fa.splitlines()
                    self.assertTrue(head.startswith('>%s_%d.%d length_bp=%d ' % (
                        arm, walk, w.path_id, j.at_bp + w.length_bp)), head)
                    self.assertTrue(head.endswith(' walk=%d.%d junctions=%d verdict=qualified '
                                                  'certified_to=%d' % (walk, w.path_id, j.at_bp,
                                                                       j.at_bp)), head)
                    self.assertEqual(seed + nat if arm == 'right' else nat + seed, body)
                    fw = sg.to_fasta(arm, walks=[key], with_seed=False, orientation='walk')
                    self.assertEqual(pw + cw, fw.splitlines()[1])

    def test_fasta_refusals(self):
        arm, p, c = triple('no_bases')
        sg = stitch(p, c)
        with self.assertRaises(ValueError) as cm:
            sg.to_fasta()
        self.assertEqual(derive.NO_BASES, str(cm.exception))
        with self.assertRaises(ValueError):
            p.to_fasta()
        arm, p, c = triple('keep_forbid')
        with self.assertRaises(ValueError):
            stitch(p, c).to_fasta(coordinates=True)


# ---------------------------------------------------------------------------------- T8

class TestBudgets(unittest.TestCase):
    """T8: every operation under budget= completes or stops at the same row; walks(by='id')
    resumes; the account bounds the traced peak."""

    def _sg(self):
        arm, p, c = triple('keep_forbid')
        return arm, stitch(stitch(p, c), load('keep_forbid', 'sibling'))

    def test_walks_stop_and_resume(self):
        arm, sg = self._sg()
        whole = [w.key for w in sg.walks(arm, by='id')]
        full = LocalBudget()
        sg.walks(arm, by='id', budget=full)
        used = full.used_work
        stops = []
        for limit in range(1, used + 1, max(1, used // 40)):
            arm, sg = self._sg()
            try:
                sg.walks(arm, by='id', budget=LocalBudget(work_units=limit))
                continue
            except LocalBudgetExceeded as e:
                part = e.partial
            stops.append((limit, [w.key for w in part.rows]))
            if part.rows:
                rest = sg.walks(arm, by='id', resume=part.resume)
                self.assertEqual(whole, [w.key for w in part.rows] + [w.key for w in rest])
        self.assertTrue(stops)
        # the same stop in every run
        for limit, rows in stops[:5]:
            arm, sg = self._sg()
            with self.assertRaises(LocalBudgetExceeded) as cm:
                sg.walks(arm, by='id', budget=LocalBudget(work_units=limit))
            self.assertEqual(rows, [w.key for w in cm.exception.partial.rows])
        with self.assertRaises(ValueError):
            sg.walks(arm, by='support', resume=('walks', 'x', 0))
        with self.assertRaises(ValueError):
            sg.walks(arm, by='id', resume=('walks', 'x', 0))

    def test_every_operation_under_a_budget(self):
        arm, p, c = triple('keep_forbid')
        ops_ = [lambda b: stitch(p, c, budget=b),
                lambda b: stitch(p, c).spell(arm, (0, 0), budget=b),
                lambda b: stitch(p, c).to_fasta(budget=b),
                lambda b: stitch(p, c).summary(budget=b),
                lambda b: stitch(p, c).walk(arm, (0, 0), budget=b)]
        for i, op in enumerate(ops_):
            with self.subTest(i):
                op(LocalBudget())
                with self.assertRaises(LocalBudgetExceeded):
                    op(LocalBudget(work_units=1))

    def test_account_bounds_the_peak(self):
        arm, p, c = triple('chain3')
        g3 = load('chain3', 'grandchild')

        def make():
            return stitch(stitch(p, c), g3)
        make().walks(arm, by='id')
        gc.collect()
        sg = make()
        b = LocalBudget()
        tracemalloc.start()
        try:
            base = tracemalloc.get_traced_memory()[0]
            sg.walks(arm, by='id', budget=b)
            peak = tracemalloc.get_traced_memory()[1] - base
        finally:
            tracemalloc.stop()
        self.assertLessEqual(peak, b.peak_bytes)
        self.assertGreaterEqual(sg.memory_bytes(),
                                K.models_memory_bytes([x.graphlet for x in sg.parts]))


if __name__ == '__main__':
    unittest.main()
