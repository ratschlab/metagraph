"""Stitching a continuation (DESIGN-graphlet-merge.md §3, phase 1a-S; Appendix E.2 T1-T9).

The CLI triples (a parent at radius R, its child from next_request(), the uninterrupted
retrieval at R + E) are in data/traverse/documents/stitch/, written by
scripts/traversal/graphlet_fixtures.py (STITCH_FIXTURES)."""

import copy
import dataclasses
import gc
import json
import os
import random
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


def _long_seed(g, n=100000):
    """A copy of graphlet |g| whose seed is |n| more bases long (A's before it): the arms and
    their bases unchanged, as an engine answer for a long seed with short flanks."""
    g = copy.deepcopy(g)
    s = 'A' * n + g.seed.sequence
    g.seed = dataclasses.replace(g.seed, sequence=s, length_bp=len(s),
                                 num_kmers=len(s) - g.k + 1)
    return g


def _peak(fn):
    """tracemalloc's peak of fn()'s second run (the first fills the interpreter's free lists),
    the collector held off; a budget stop ends a run."""
    gc.collect()
    gc.disable()
    try:
        for i in range(2):
            if i:
                tracemalloc.start()
                base = tracemalloc.get_traced_memory()[0]
            try:
                fn()
            except LocalBudgetExceeded:
                pass
        return tracemalloc.get_traced_memory()[1] - base
    finally:
        tracemalloc.stop()
        gc.enable()


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


def with_dropped(name):
    """An X line in the body: the child dropped seed label |name| (seed_unsupported)."""
    def edit(resp):
        r = resp['results'][0]
        lines = r['graphlet'].split('\n')
        i = next(i for i, l in enumerate(lines) if l.startswith('S '))
        lines.insert(i + 1, 'X seed_unsupported 0-0 ' + name)
        z = next(i for i, l in enumerate(lines) if l.startswith('Z '))
        lines[z] = 'Z %d' % (int(lines[z].split()[1]) + 1)
        r['graphlet'] = '\n'.join(lines)
        r.pop('graphlet_bytes', None)
        r.pop('graphlet_lines', None)
    return edit


def only_clamped(resp):
    """An envelope whose normalized strategy states nothing but clamped."""
    resp['strategy'] = {'clamped': []}


def without(block, name):
    def edit(resp):
        del resp['strategy'][block][name]
    return edit


def quorum(resp):
    resp['strategy']['branching']['min_successor_labels'] = 2


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
                    nr = child.next_request(fx['arm'], [req['fixture_index']['grandchild']],
                                            fx['grandchild']['bp'])
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

    def test_identity_over_every_part(self):
        """M4: the index proof and the pair check run against every part on both sides,
        and allow_unverified_index waives only 'unverifiable': one code in both folds
        (Appendix C.3's associativity)."""
        def chain(fps, edits=(None, None, None)):
            out = []
            for role, fp, e in zip(('parent', 'child', 'grandchild'), fps, edits):
                def edit(resp, fp=fp, e=e):
                    if fp:
                        with_fp(fp)(resp)
                    if e:
                        e(resp)
                out.append(load('chain3', role, edit))
            return out

        def folds(p, c, g, **kw):
            got = []
            for fold in (lambda: stitch(stitch(p, c, **kw), g, **kw),
                         lambda: stitch(p, stitch(c, g, **kw), **kw)):
                try:
                    got.append(fold().digest)
                except CombineRefused as e:
                    got.append(e.code)
            return got
        a, b = 'a' * 16, 'b' * 16
        allow = {'allow_unverified_index': True}
        # different digests anywhere in the chain: different_index in both folds, also
        # when a part between them states none (the review's exp12)
        for fps in ((a, None, b), (None, a, b), (a, b, None)):
            with self.subTest(fps):
                self.assertEqual(['different_index'] * 2, folds(*chain(fps), **allow))
        refused(self, 'different_index', 'the child was retrieved from another index than its '
                'parent (index_fp differs): a continuation runs on its parent\'s index, so a '
                'child of another index continues none of its walks; union_pairs() combines '
                'one seed\'s graphlets across pairs', stitch,
                stitch(*chain((a, None, b))[:2], **allow), chain((a, None, b))[2], **allow)
        # one digest on some parts: unverifiable without the opt-in, in both folds; with it
        # both folds give one digest and the identity unverified
        p, c, g = chain((a, None, a))
        self.assertEqual(['index_unverifiable'] * 2, folds(p, c, g))
        x, y = folds(p, c, g, **allow)
        self.assertEqual(x, y)
        sg = stitch(stitch(p, c, **allow), g, **allow)
        self.assertEqual((False, 'allow_unverified_index'),
                         (sg.identity['verified'], sg.identity['accepted']))
        refused(self, 'index_unverifiable', 'the parts state an index manifest digest '
                '(index_fp) unevenly (the parent\'s on 1 of 2 parts, the child\'s on 1 of 1): '
                'whether they are the same index cannot be verified (not proven different), so '
                'the child is not stitched automatically; pass allow_unverified_index=True to '
                'stitch it anyway (the result then states the identity unverified)',
                stitch, stitch(p, c, **allow), g)
        # every part states it: verified in both folds
        p, c, g = chain((a, a, a))
        x, y = folds(p, c, g)
        self.assertEqual(x, y)
        self.assertEqual({'verified': True, 'index_fp': a}, stitch(stitch(p, c), g).identity)
        # the pair check over every part too: the grandchild's pair against the root's
        p, c, g = chain((None, None, None), (with_pair('g1', '/p/g1.dbg'), None,
                                             with_pair('g2', '/p/g2.dbg')))
        self.assertEqual(['different_pair'] * 2, folds(p, c, g))

    def test_different_pair(self):
        p = load('keep_forbid', 'parent', with_pair('g1', '/p/g1.dbg'))
        c = load('keep_forbid', 'child', with_pair('g2', '/p/g2.dbg'))
        # (a pair stated without its annotation path has no n1 id: it is shown as stated;
        # an engine states both paths or the id, test_traverse_pair_id)
        refused(self, 'different_pair', 'the child was answered on pair '
                '{"graph":"g2","graph_path":"/p/g2.dbg"}, its parent on '
                'pair {"graph":"g1","graph_path":"/p/g1.dbg"}: a '
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
        # walk= on a part without the arm
        refused(self, 'not_a_continuation', 'the child\'s seed is not the continuation of walk '
                '0 of the left arm: a stitched child is seeded with that walk\'s continuation '
                '(continuation(), next_request())', stitch, self.P, load('fork_left', 'child'),
                walk=0)
        with self.assertRaises(BadSelector):
            stitch(self.P, self.C, walk='0')
        # a walk= the arm does not have: CombineRefused, not IndexError
        refused(self, 'not_a_continuation', 'walk= names no walk of the right arm (the right '
                'arm of part root has 2 walks, no walk 7): a stitched child is seeded with a '
                'walk\'s continuation (continuation(), next_request())', stitch, self.P, self.C,
                walk=7)
        _, p3, c3 = triple('chain3')
        g3 = load('chain3', 'grandchild')
        refused(self, 'not_a_continuation', 'walk= names no walk of the right arm (no walk 3.0 '
                'on the right arm: walk 3 is not continued here): a stitched child is seeded '
                'with a walk\'s continuation (continuation(), next_request())', stitch,
                stitch(p3, c3), g3, walk=(3, 0))

    def test_already_continued(self):
        _, lp, lc = triple('loop_branch')
        cut0 = load('loop_branch_cut0', 'child')
        sg = stitch(lp, lc)
        # the digest named is the subtree's there (stitch() of it alone)
        refused(self, 'already_continued', 'walk 1 of the right arm is already continued by '
                'another child (part right:1, digest %s): one child per walk'
                % stitch(lc).digest, stitch, sg, cut0)

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

    def test_retrieval_identity_is_not_cached_for_a_child(self):
        """The child's retrieval identity is computed at each stitch, not cached by id():
        the same object (so the same address) from another request is refused (the review
        of the assembled train, M2)."""
        p, c = load('keep_forbid', 'parent'), load('keep_forbid', 'child')
        sg = stitch(p, c)
        x = copy.copy(c)
        self.assertIs(sg, stitch(sg, x))
        x.envelope = copy.deepcopy(c.envelope)
        x.envelope['strategy']['branching']['min_successor_labels'] = 7
        refused(self, 'already_continued', 'walk 0 of the right arm is continued by this body '
                'from another request (strategy.branching.min_successor_labels differs): equal '
                'bodies are not equal retrievals, and the two would continue differently',
                stitch, sg, x)
        self.assertNotIn(('rid', id(x)), sg._jcore)

    def test_already_continued_compares_the_subtree(self):
        """The no-op rule compares the whole subtree, digests and retrieval identities (the
        review's exp5); the text names the existing subtree's digest."""
        p, c, g = load('chain3', 'parent'), load('chain3', 'child'), load('chain3', 'grandchild')
        sga = stitch(p, stitch(c, g))
        self.assertIs(sga, stitch(sga, stitch(c, g)))
        # a plain graphlet is compared with the part at the walk only, whatever continues it
        # (the review of the assembled train, minor 5)
        sg3 = stitch(stitch(p, c), g)
        self.assertIs(sg3, stitch(sg3, c))
        self.assertIs(sg3, stitch(sg3, g))
        self.assertIs(sg3, stitch(sg3, load('chain3', 'child')))
        refused(self, 'already_continued', 'walk 0 of the right arm is continued by this body '
                'from another request (strategy.branching.min_successor_labels differs): equal '
                'bodies are not equal retrievals, and the two would continue differently',
                stitch, sg3, load('chain3', 'child', quorum))
        # the grandchild is the same body from another request
        refused(self, 'already_continued', 'walk 0 of the right arm is continued by this body '
                'from another request (strategy.branching.min_successor_labels of part '
                'right:0.0 differs): equal bodies are not equal retrievals, and the two would '
                'continue differently', stitch, sga, stitch(c, load('chain3', 'grandchild',
                                                                     quorum)))
        # a child without the grandchild the subtree holds
        refused(self, 'already_continued', 'walk 0 of the right arm is already continued by '
                'another child (part right:0, digest %s): one child per walk'
                % stitch(c, g).digest, stitch, sga, stitch(c))
        refused(self, 'already_continued', 'walk 0 of the right arm is already continued by '
                'another child (part right:0, digest %s): one child per walk'
                % stitch(c).digest, stitch, stitch(p, c), stitch(c, g))

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

    def test_seeded_random_orders(self):
        """Commutativity and associativity as digest equality over seeded random
        construction orders: any two items of the pool that stitch (either way round) are
        replaced by their stitch until one is left."""
        def build(parts, rng):
            pool = list(parts)
            while len(pool) > 1:
                pairs = [(i, j) for i in range(len(pool)) for j in range(len(pool)) if i != j]
                rng.shuffle(pairs)
                for i, j in pairs:
                    try:
                        sg = stitch(pool[i], pool[j])
                    except CombineRefused:
                        continue
                    pool = [x for k, x in enumerate(pool) if k not in (i, j)] + [sg]
                    break
                else:
                    self.fail('no two items of the pool stitch')
            return pool[0]
        sets = {'keep_forbid': ('parent', 'child', 'sibling'),
                'chain3': ('parent', 'child', 'grandchild'),
                'chain_loop': ('parent', 'child', 'grandchild'),
                'splits_chain': ('parent', 'child', 'grandchild'),
                'comb_chain': ('parent', 'child', 'grandchild')}
        for name, roles in sets.items():
            with self.subTest(name):
                parts = [load(name, r) for r in roles]
                want = stitch(stitch(parts[0], parts[1]), parts[2])
                for seed in range(12):
                    rng = random.Random(1000 * len(name) + seed)
                    order = list(parts)
                    rng.shuffle(order)
                    got = build(order, rng)
                    self.assertEqual(want.digest, got.digest, (name, seed))
                    self.assertEqual([x.key for x in want.parts], [x.key for x in got.parts])

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
        # the digest covers each part's retrieval identity, as the union's does: the same
        # child body from another request gives another digest
        other = stitch(p, load('keep_forbid', 'child', quorum))
        self.assertEqual(K.part_digest(c), K.part_digest(other.parts[1].graphlet))
        self.assertNotEqual(sg.digest, other.digest)


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
    # the Opus code review's cases: the parent without bases, the child with them
    'no_bases_parent': (('no_bases',), {(0, 0): (30, 15, False, (), 15, False)},
                        ('partial', 'complete', 'qualified', 'inline')),
    # canonical: the step into rc(R) reuses an edge on the other strand (edge_reuse_rc)
    'invrep': (('exact_gated',), {(0,): (96, 96, True, (), 96, True),
                                  (1, 0): (142, 96, False, (('rule_violation', 111),
                                                            ('candidate_changed', 112)),
                                           111, False),
                                  (1, 1): (182, 96, False, (('rule_violation', 111),),
                                           111, False)},
               ('partial', 'complete', 'qualified', 'inline')),
    # canonical: the step into rc(S), the root seed's nodes on the other strand
    'seed_rc': (('exact_gated',), {(0, 0): (70, 20, False, (('rule_violation', 54),), 54, False),
                                   (1,): (20, 20, True, (), 20, True)},
                ('partial', 'complete', 'qualified', 'inline')),
    'chain_loop': (('exact_gated',), {(0,): (85, 85, True, (), 85, True),
                                      (1, 0): (115, 100, False, (), 115, True)},
                   ('complete', 'complete', 'complete', 'inline')),
    'splits_chain': (('exact_gated',), {(0,): (18, 18, True, (), 18, True),
                                        (1, 0): (49, 45, False, (), 49, True)},
                     ('partial', 'complete', 'complete', 'inline')),
    # W6 by length: the sibling ended at 60 without a continuation
    'beam_sibling': (('label_count_knob', 'sibling_merge_unknown'),
                     {(0, 0): (10, 1, False, (), 1, False), (1,): (60, 1, False, (), 1, False)},
                     ('partial', 'complete', 'qualified', 'inline')),
    # annotate mode: the split's other child is a step the reference blocks
    'annotate_loop': (('exact_gated',), {(0,): (85, 85, True, (), 85, True),
                                         (1, 0): (145, 90, False, (('candidate_changed', 115),),
                                                  115, False),
                                         (1, 1): (149, 90, False, (('rule_violation', 115),
                                                                   ('seed_reentry', 149)),
                                                  115, False)},
                      ('partial', 'complete', 'qualified', 'inline')),
    'annotate_merge': (('strategy_differs',), {(0, 0): (40, 20, False, (), 20, False),
                                               (1,): (20, 20, True, (), 20, True)},
                       ('partial', 'complete', 'qualified', 'inline')),
    # a followed hairpin needs no cause: the prospective walks are the reference's
    'hairpin': (('exact_gated',), {(0, 0): (56, 5, False, (), 56, True),
                                   (0, 1): (56, 5, False, (), 56, True)},
                ('complete', 'complete', 'complete', 'inline')),
}


def _comb_rows(depth):
    """The comb chain's stitched walks: the root's 20 tips, then per part below the root
    19 tips and the next part's backbone (the last part's two radius-ended walks), all
    prospective exact; a walk that reaches into a child is certified to the first junction."""
    rows = {(i,): (i + 1, i + 1, True, (), i + 1, True) for i in range(20)}
    pre = (20,)
    for d in range(1, depth + 1):
        last = d == depth
        for i in range(21 if last else 19):
            n = 20 * d + 1 + min(i, 19)
            rows[pre + (i,)] = (n, 20, False, (), n, True)
        if not last:
            rows[pre + (20,)] = (20 * (d + 1), 20, False, (), 20 * (d + 1), True)
            pre = pre + (19,)
    return rows


EXPECT['comb_chain'] = (('exact_gated',), _comb_rows(1), ('complete', 'complete', 'complete',
                                                          'inline'))

# the chains (parent, child, grandchild): per junction its causes, the walks as in EXPECT
EXPECT_CHAIN = {
    'chain3': ((('exact_gated',), ('exact_gated',)),
               {(0, 0, 0): (45, 15, False, (), 45, True)},
               ('complete', 'complete', 'complete', 'inline')),
    # M1: the root carries no bases, so neither junction runs the rule check
    'no_bases_parent': ((('no_bases',), ('no_bases',)),
                        {(0, 0, 0): (45, 15, False, (), 15, False)},
                        ('partial', 'complete', 'qualified', 'inline')),
    # the root part's edges: the grandchild takes the loop's first edge again
    'chain_loop': ((('exact_gated',), ('exact_gated',)),
                   {(0,): (85, 85, True, (), 85, True),
                    (1, 0, 0): (145, 100, False, (('candidate_changed', 115),), 115, False),
                    (1, 0, 1): (174, 100, False, (('rule_violation', 115), ('seed_reentry', 174)),
                                115, False)},
                   ('partial', 'complete', 'qualified', 'inline')),
    # the root's two splits carried through the child into the grandchild's
    'splits_chain': ((('exact_gated',), ('exact_gated',)),
                     {(0,): (18, 18, True, (), 18, True),
                      (1, 0, 0): (89, 45, False, (('split_limit_restarted', 50),), 50, False),
                      (1, 0, 1): (89, 45, False, (('split_limit_restarted', 50),), 50, False)},
                     ('partial', 'complete', 'qualified', 'inline')),
    'comb_chain': ((('exact_gated',), ('exact_gated',)), _comb_rows(2),
                   ('complete', 'complete', 'complete', 'inline')),
}

SEED_REENTRY = {'loop_same_stop': 'same_stop', 'earlier_stop': 'earlier_stop',
                'loop_interior': 'same_stop', 'annotate_loop': 'same_stop'}


class TestOracleTriples(unittest.TestCase):
    """T4: per CLI triple, the junction's verdict and causes, every stitched walk's
    certified_to_bp, complete and causes, and the evidence, public and prospective."""

    def test_every_triple_is_pinned(self):
        self.assertEqual(sorted(T.stitch_names()), sorted(EXPECT))
        self.assertEqual(sorted(n for n in T.stitch_names()
                                if 'grandchild' in T.stitch_requests(n)), sorted(EXPECT_CHAIN))

    def _rows(self, sg, arm):
        got = {}
        for w in sg.walks(arm, by='id'):
            pr = sg._prospective(arm, w.key)
            own = tuple((x[0], x[2]) for x in w.causes if x[0] != 'exact_gated')
            got[w.key] = (w.length_bp, w.certified_to_bp, w.complete, own,
                          pr['certified_to_bp'], pr['complete'])
        return got

    def test_chains(self):
        """The chains of three parts: each junction's causes, every walk, the prospective
        evidence; the completeness extensions state the walk-scope causes below each."""
        for name, (jcauses, walks, prospective_outcome) in EXPECT_CHAIN.items():
            with self.subTest(name):
                arm, p, c = triple(name)
                sg = stitch(stitch(p, c), load(name, 'grandchild'))
                self.assertEqual(jcauses, tuple(tuple(x[0] for x in j.causes)
                                                for j in sg.junctions))
                self.assertEqual(walks, self._rows(sg, arm))
                pe = sg._prospective(arm)['evidence']
                self.assertEqual(prospective_outcome, tuple(pe['outcome'][d] for d in RANK))
                self.assertFalse(sg.evidence(arm)['exact'])
                for ext, j in zip(sg.completeness()[arm]['extensions'], sg.junctions):
                    below = sorted({x[0] for w in sg.walks(arm, by='id')
                                    if w.key[:len(j.walk)] == j.walk and len(w.key) > len(j.walk)
                                    for x in w.causes if x[0] != 'exact_gated'})
                    self.assertEqual(below, sorted(ext['walk_causes']))

    def test_strategy_unknown(self):
        """M2: a strategy field W2-W5 or the split check reads, missing on either part, is
        strategy_unknown; never exact, also not in the prospective result (the review's
        exp6 and exp7)."""
        def sentence(fields):
            return ('the envelopes do not state %s of the normalized strategy: the rules that '
                    'read them (W2-W5, the split check of W8) cannot be judged' % fields)
        w5 = ('branching.min_successor_labels (both), branching.min_successor_fraction (both), '
              'bounds.min_live_labels (both), frontier.on_overflow (both)')
        constrain = w5 + (', labels.change_cost (both), labels.extra (both), '
                          'branching.max_label_branches (both), branching.max_splits_per_path '
                          '(both)')
        # then the 14 other fields W2 compares, which neither states
        constrain += ' and 14 more'
        w5 += (', exhaustive (both), schema_version (both), support (both), '
               'bounds.max_live_paths (both) and 12 more')
        cases = {'splits45': (('strategy_unknown',), constrain),
                 'quorum_changed': (('strategy_unknown',), constrain),
                 'trace': (('trace_restart', 'strategy_unknown'), constrain),
                 'annotate_loop': (('strategy_unknown',), w5),
                 # W6's keep-or-merge from the H records, not the envelopes
                 'annotate_merge': (('strategy_differs', 'strategy_unknown'), w5)}
        for name, (codes, fields) in cases.items():
            with self.subTest(name):
                arm = T.stitch_requests(name)['fixture_index']['arm']
                sg = stitch(load(name, 'parent', only_clamped), load(name, 'child', only_clamped))
                j = sg.junctions[0]
                self.assertEqual(codes, tuple(x[0] for x in j.causes))
                self.assertEqual(sentence(fields), dict((x[0], x[1]) for x in j.causes)[
                    'strategy_unknown'])
                self.assertFalse(sg._prospective(arm)['evidence']['exact'])
                for w in sg.walks(arm, by='id'):
                    if len(w.key) > 1:
                        self.assertEqual('qualified', sg._prospective(arm, w.key)['verdict'])
                        self.assertNotIn('exact_gated', [x[0] for x in w.causes])
        j = stitch(load('annotate_merge', 'parent', only_clamped),
                   load('annotate_merge', 'child', only_clamped)).junctions[0]
        self.assertEqual('the child\'s normalized strategy differs from the parent\'s at '
                         'branching.on_reconverge (keep in the parent, merge in the child, by '
                         'their H records), a field a continuation keeps (W2, W6)', j.causes[0][1])
        # a missing max_splits_per_path is unknown, not unlimited; a field one part states
        # and the other does not is unknown too
        for name, who, edit in (('splits45', 'parent', without('branching', 'max_splits_per_path')),
                                ('switching', 'child', without('labels', 'extra'))):
            with self.subTest(name):
                p = load(name, 'parent', edit if who == 'parent' else None)
                c = load(name, 'child', edit if who == 'child' else None)
                j = stitch(p, c).junctions[0]
                self.assertIn('strategy_unknown', [x[0] for x in j.causes])
                self.assertNotIn('exact_gated', [x[0] for x in j.causes])
        j = stitch(load('splits45', 'parent', without('branching', 'max_splits_per_path')),
                   load('splits45', 'child')).junctions[0]
        self.assertEqual([('strategy_unknown', sentence('branching.max_splits_per_path (the '
                                                        'parent)'))],
                         [(x[0], x[1]) for x in j.causes])
        # a null field is not stated: unknown, not a difference
        def null_splits(resp):
            resp['strategy']['branching']['max_splits_per_path'] = None
        j = stitch(load('splits45', 'parent', null_splits), load('splits45', 'child')).junctions[0]
        self.assertEqual([('strategy_unknown', sentence('branching.max_splits_per_path (the '
                                                        'parent)'))],
                         [(x[0], x[1]) for x in j.causes])
        # a field W2 compares that both envelopes omit is unknown too (the review of the
        # assembled train, minor 6)
        for block, name in (('branching', 'hairpins'), ('frontier', 'order')):
            with self.subTest(name):
                j = stitch(load('keep_forbid', 'parent', without(block, name)),
                           load('keep_forbid', 'child', without(block, name))).junctions[0]
                self.assertEqual([('strategy_unknown', sentence('%s.%s (both)' % (block, name)))],
                                 [(x[0], x[1]) for x in j.causes])

    def test_no_bases_names_the_parts(self):
        """no_bases names the parts without bases on the route from the root."""
        text = ('on the route from the root to this child, %s no bases (output.sequences '
                'false): the rule check (W8) cannot run below this junction')
        arm, p, c = triple('no_bases')
        self.assertEqual([text % 'parts root, right:0 carry'],
                         [x[1] for x in stitch(p, c).junctions[0].causes])
        arm, p, c = triple('no_bases_parent')
        sg = stitch(stitch(p, c), load('no_bases_parent', 'grandchild'))
        self.assertEqual([[text % 'part root carries']] * 2,
                         [[x[1] for x in j.causes] for j in sg.junctions])

    def test_seed_dropped(self):
        """W7: a label alive at the end of the walk that the child named and dropped (X) is
        seed_dropped, not not_continued. The CLI fixtures cannot reach it (the library's
        labels alive but not covering the tail are ones the server accepts when named), so
        the child's body states the drop."""
        arm, p, _ = triple('not_continued')
        c = load('not_continued', 'child', with_dropped('y.fa'))
        j = stitch(p, c).junctions[0]
        self.assertEqual(('switching', 'seed_dropped'), tuple(x[0] for x in j.causes))
        self.assertEqual('the child named {"name": "y.fa", "ref": "c:2"} but dropped it (X): '
                         'y.fa seed_unsupported (W7)', j.causes[1][1])
        self.assertEqual(['c:2'], [l.ref for l in j.not_continued])

    def test_whitelist_subtrees_equal_the_reference(self):
        """T4(a) and its converse: on each whitelist fixture the stitched walks below w are
        the uninterrupted retrieval's walks with w's prefix, per path."""
        for name in ('keep_forbid', 'merge_only', 'unlimited', 'equal_finite', 'chain3',
                     'comb_chain', 'hairpin'):
            with self.subTest(name):
                arm, p, c = triple(name)
                w = T.stitch_requests(name)['fixture_index']['walk']
                sg = stitch(p, c)
                if 'grandchild' in T.stitch_requests(name):
                    sg = stitch(sg, load(name, 'grandchild'))
                ref = load(name, 'reference')
                pre = p.spell(arm, w)
                spelled = [ref.spell(arm, x.path_id) for x in ref.walks(arm, by='id')]
                below = (sorted(x for x in spelled if x.startswith(pre)) if arm == 'right'
                         else sorted(x for x in spelled if x.endswith(pre)))
                mine = sorted(x.sequence for x in sg.walks(arm, by='id') if x.key[0] == w)
                self.assertEqual(below, mine)
                for x in sg.walks(arm, by='id'):
                    if x.key[0] == w:
                        self.assertEqual('exact', sg._prospective(arm, x.key)['verdict'])

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
        uninterrupted retrieval, base for base (the chains against theirs)."""
        for name in T.stitch_names():
            req = T.stitch_requests(name)
            if 'reference' not in req:
                continue
            arm, p, c = triple(name)
            if not (K.has_bases(p) and K.has_bases(c)):
                continue
            with self.subTest(name):
                ref = load(name, 'reference')
                spelled = [ref.spell(arm, w.path_id) for w in ref.walks(arm, by='id')]
                if arm == 'left':
                    spelled = [x[::-1] for x in spelled]
                sg = stitch(p, c)
                if 'grandchild' in req:
                    sg = stitch(sg, load(name, 'grandchild'))
                for w in sg.walks(arm, by='id'):
                    if sg._prospective(arm, w.key)['verdict'] != 'exact':
                        continue
                    # a walk ended by the child's radius goes on in the reference (its
                    # radius is larger); any other ends where the reference's ends
                    # (outward: reversed on the left arm)
                    seq = w.sequence[::-1] if arm == 'left' else w.sequence
                    if w.path_reason == 'max_extension_bp':
                        self.assertTrue(any(x.startswith(seq) for x in spelled), w.key)
                    else:
                        self.assertIn(seq, spelled, w.key)

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

    def test_either_arm_first(self):
        """The rule check's windows are kept per part and arm: the root walks both arms,
        and its left tails must never stand in for its right ones. The right child retakes
        the parent's seed edges (its walk edited, same length), so its walk breaks the rule
        at 25 bp whichever arm is queried first (GPT review M1)."""
        def retaken(resp):
            r = resp['results'][0]
            rows = r['graphlet'].split('\n')
            for i, line in enumerate(rows):
                if line.startswith('G '):
                    f = line.split(' ')
                    f[-1] = 'TGCAGGTTCTATTATT' + 'ACGTACGTACGTAC'
                    rows[i] = ' '.join(f)
            r['graphlet'] = '\n'.join(rows)
            r.pop('graphlet_bytes', None)
            r.pop('graphlet_lines', None)

        def make():
            return stitch(stitch(load('fork_right', 'parent'),
                                 load('fork_right', 'child', retaken)),
                          load('fork_left', 'child'))
        answers = []
        for order in (('right', 'left'), ('left', 'right')):
            sg = make()
            got = {side: [(w.key, w.causes) for w in sg.walks(side, by='id')]
                   for side in order}
            got['evidence'] = sg.evidence()
            got['summary'] = sg.summary()
            answers.append(got)
        self.assertEqual(answers[0], answers[1])
        causes = dict(answers[0]['right'])[(0, 0)]
        self.assertEqual([('rule_violation', 25)], [(c[0], c[2]) for c in causes])
        self.assertEqual((), dict(answers[0]['left'])[(0, 0)])


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
                # per junction, the walk-scope cause codes on the walks below it
                for js, j in zip(s['junctions'], sg.junctions):
                    below = sorted({x[0] for w in sg.walks(arm, by='id')
                                    if w.key[:len(j.walk)] == j.walk and len(w.key) > len(j.walk)
                                    for x in w.causes if x[0] != 'exact_gated'})
                    self.assertEqual(below, sorted(js['walk_causes']))
                # junctions are cut from the end; the statements never are
                floor = len(json.dumps(dict(s, junctions=[], junctions_cut=True),
                                       separators=(',', ':'), ensure_ascii=False).encode())
                small = sg.summary(max_bytes=floor)
                self.assertEqual([], small['junctions'])
                self.assertTrue(small['junctions_cut'])
                refused(self, 'summary_too_small', 'the summary\'s statements alone take %d '
                        'bytes, more than max_bytes=%d: they are never cut; pass max_bytes >= '
                        '%d' % (floor, floor - 1, floor), sg.summary, max_bytes=floor - 1)

    def test_summary_too_small_names_its_least_size(self):
        """summary_too_small carries the least answer's size as least_bytes, the number its
        text states: summary(max_bytes=least_bytes) answers, and a byte less is refused with
        the same size (as the union's refusal)."""
        for name, arm, sg in every_stitch():
            with self.subTest(name):
                with self.assertRaises(CombineRefused) as cm:
                    sg.summary(max_bytes=10)
                e = cm.exception
                self.assertEqual('summary_too_small', e.code)
                m = re.search(r'take (\d+) bytes, more than max_bytes=10: .*pass max_bytes >= '
                              r'(\d+)$', str(e))
                n = e.least_bytes
                self.assertEqual((n, n), (int(m.group(1)), int(m.group(2))))
                got = sg.summary(max_bytes=n)
                self.assertLessEqual(len(json.dumps(got, separators=(',', ':'),
                                                    ensure_ascii=False).encode()), n)
                with self.assertRaises(CombineRefused) as cm:
                    sg.summary(max_bytes=n - 1)
                self.assertEqual(n, cm.exception.least_bytes)
                self.assertTrue(str(cm.exception).endswith('pass max_bytes >= %d' % n))

    def test_walk_causes_below_a_gated_junction(self):
        """M3: a junction whose rules found nothing says W1-W7 hold there and W8 is judged
        per walk; summary() and completeness() state the walk-scope causes below it."""
        arm, p, c = triple('loop_branch')
        sg = stitch(p, c)
        j = sg.junctions[0]
        self.assertEqual([('exact_gated', 'W1-W7 hold at this junction, and W8 is judged per '
                           'walk below it (the walks\' causes), but exact is not switched on: '
                           'the junction rules found no difference from one uninterrupted walk '
                           'and have not yet passed the randomized comparison against '
                           'uninterrupted traversals')], [(x[0], x[1]) for x in j.causes])
        self.assertEqual([{'arm': 'right', 'walk': '1', 'at_bp': 115, 'verdict': 'qualified',
                           'causes': ['exact_gated'], 'walk_causes': ['candidate_changed']}],
                         sg.summary()['junctions'])
        self.assertEqual([{'walk': (1,), 'at_bp': 115, 'extended_to_bp': 1000,
                           'certified_to_bp': 115, 'verdict': 'qualified',
                           'causes': ['exact_gated'], 'walk_causes': ['candidate_changed']}],
                         sg.completeness()[arm]['extensions'])


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

    def test_walk_keys(self):
        """walks= is one key (an int, or a tuple of ints) or an iterable of keys."""
        arm, p, c = triple('keep_forbid')
        sg = stitch(stitch(p, c), load('keep_forbid', 'sibling'))
        self.assertEqual(sg.to_fasta(arm, walks=[(0, 0)]), sg.to_fasta(arm, walks=(0, 0)))
        self.assertEqual(sg.evidence(arm, walks=[(0, 0)]), sg.evidence(arm, walks=(0, 0)))
        self.assertEqual(sg.to_fasta(arm, walks=[(0, 0), (1, 0)]),
                         sg.to_fasta(arm, walks=iter([(0, 0), (1, 0)])))
        _, lp, lc = triple('loop_branch')
        one = stitch(lp, lc)
        self.assertEqual(one.to_fasta(arm, walks=[(0,)]), one.to_fasta(arm, walks=0))
        self.assertEqual(one.to_fasta(arm, walks=[(0,)]), one.to_fasta(arm, walks=[0]))

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
                lambda b: stitch(p, c).walk(arm, (0, 0), budget=b),
                lambda b: stitch(p, c).walks(arm, budget=b),
                lambda b: stitch(p, c).evidence(budget=b),
                lambda b: stitch(p, c).completeness(budget=b),
                lambda b: stitch(p, c).continuation(arm, (1,), budget=b),
                lambda b: stitch(p, c).next_request(arm, (1,), budget=b),
                lambda b: stitch(p, c).next_requests(arm, [(1,)], budget=b)]
        for i, op in enumerate(ops_):
            with self.subTest(i):
                op(LocalBudget())
                with self.assertRaises(LocalBudgetExceeded):
                    op(LocalBudget(work_units=1))
                # without budget=, the ambient one
                with MT.local_budget() as b:
                    op(None)
                self.assertGreater(b.used_work, 0)

    # objects built from freshly parsed parts (their own caches cold too): a chain of comb
    # continuations (61 walks, two junctions in a row), its first junction alone, and a
    # root with a child on each arm (tails per arm)

    @staticmethod
    def _comb():
        _, p, c = triple('comb_chain')
        return stitch(stitch(p, c), load('comb_chain', 'grandchild'))

    @staticmethod
    def _comb_one():
        _, p, c = triple('comb_chain')
        return stitch(p, c)

    @staticmethod
    def _fork():
        return stitch(stitch(load('fork_right', 'parent'), load('fork_right', 'child')),
                      load('fork_left', 'child'))

    def _budgeted(self):
        """(name, make, op): every public StitchedGraphlet method that takes budget=, each on
        an object make() builds."""
        arm = 'right'
        comb, comb_one, fork = self._comb, self._comb_one, self._fork
        g3 = load('comb_chain', 'grandchild')
        deep = (20, 19, 19)             # through both junctions, with a continuation
        return [
            ('stitch()', comb_one, lambda sg, b: stitch(sg, g3, budget=b)),
            ('.stitch()', comb_one, lambda sg, b: sg.stitch(g3, budget=b)),
            ('walks id', comb, lambda sg, b: sg.walks(arm, by='id', budget=b)),
            ('walks support', comb, lambda sg, b: sg.walks(arm, budget=b)),
            ('walks length', comb, lambda sg, b: sg.walks(arm, by='length', min_bp=60,
                                                          budget=b)),
            ('walk', comb, lambda sg, b: sg.walk(arm, deep, budget=b)),
            ('spell', comb, lambda sg, b: sg.spell(arm, deep, 'walk', True, budget=b)),
            ('continuation', comb, lambda sg, b: sg.continuation(arm, deep, budget=b)),
            ('next_request', comb, lambda sg, b: sg.next_request(arm, deep, budget=b)),
            ('next_requests', comb, lambda sg, b: sg.next_requests(arm, [deep], budget=b)),
            ('to_fasta', comb, lambda sg, b: sg.to_fasta(budget=b)),
            ('to_fasta walks', comb, lambda sg, b: sg.to_fasta(arm, [deep, (19,)], budget=b)),
            ('summary', comb, lambda sg, b: sg.summary(budget=b)),
            ('evidence', comb, lambda sg, b: sg.evidence(budget=b)),
            ('evidence walks', comb, lambda sg, b: sg.evidence(arm, [deep], budget=b)),
            ('completeness', comb, lambda sg, b: sg.completeness(budget=b)),
            ('fork walk left', fork, lambda sg, b: sg.walk('left', (0, 0), budget=b)),
            ('fork walk right', fork, lambda sg, b: sg.walk('right', (0, 0), budget=b)),
            ('fork summary', fork, lambda sg, b: sg.summary(budget=b)),
            ('fork evidence', fork, lambda sg, b: sg.evidence(budget=b)),
            ('fork completeness', fork, lambda sg, b: sg.completeness(budget=b)),
            ('fork to_fasta', fork, lambda sg, b: sg.to_fasta(budget=b))]

    def test_walk_charges_its_row_and_spelling(self):
        """walk() charges the walk it builds as walks() charges each walk it returns (its
        row, pieces and spelling), beside the rows and causes it reads (GPT review M3: a
        warmed walk() returned its sequence charging nothing)."""
        class Log(LocalBudget):
            def __init__(self):
                super().__init__()
                self.log = []

            def charge(self, work, nbytes=0):
                self.log.append((work, nbytes))
                return super().charge(work, nbytes)
        for make, arm, keys in ((self._comb, 'right', [(19,), (20, 19, 19)]),
                                (self._fork, 'left', [(0, 0), (1,)])):
            sg = self._warm(make())
            listed = Log()
            sg.walks(arm, by='id', budget=listed)
            for key in keys:
                one = Log()
                w = sg.walk(arm, key, budget=one)
                price = StitchedGraphlet._walk_price(w.length_bp, sg._chain(arm, key))
                self.assertIn(price, one.log, key)
                self.assertIn(price, listed.log, key)
                self.assertGreaterEqual(price[1], 3 * len(w.sequence))

    def test_summary_charges_the_evidence_it_reads(self):
        """summary() states evidence()'s outcome and exact and reads it under its own budget
        (not unbudgeted, not charged twice): it charges at least what evidence() does."""
        for make in (self._comb, self._comb_one, self._fork):
            sg = make()
            ev, sm = LocalBudget(), LocalBudget()
            sg.evidence(budget=ev)
            sg.summary(budget=sm)
            self.assertGreaterEqual(sm.used_work, ev.used_work)
            self.assertGreaterEqual(sm.peak_bytes, ev.peak_bytes)

    @staticmethod
    def _warm(sg):
        """Every query of |sg| once, unbudgeted, both arms: every cache it keeps filled."""
        for side in sg.arms():
            for w in sg.walks(side, by='id'):
                sg.walk(side, w.key)
                sg.spell(side, w.key)
            sg.walks(side)
            sg.evidence(side)
        sg.evidence()
        sg.completeness()
        sg.summary()
        sg.to_fasta()
        return sg

    @staticmethod
    def _did(op, sg, b):
        """What a budgeted call did: its charge, or its stop (with the partial's rows)."""
        try:
            op(sg, b)
            return ('done', b.used_work, b.peak_bytes)
        except LocalBudgetExceeded as e:
            s = e.stop
            rows = None if e.partial is None else [w.key for w in e.partial.rows]
            return ('stop', s.resource, s.op, s.phase, s.used, s.needed_at_least, rows)

    def test_cold_and_warm_charge_the_same(self):
        """A call charges the same whether its derivations are cached or not, a hit paying
        the cold price, and so stops at the same point (GPT review M3): every public method
        that takes budget=, on a fresh object and on one every query has warmed, at its whole
        charge and at work and memory limits below it."""
        for name, make, op in self._budgeted():
            with self.subTest(name):
                cold = LocalBudget()
                op(make(), cold)
                warm = LocalBudget()
                op(self._warm(make()), warm)
                self.assertEqual((cold.used_work, cold.peak_bytes),
                                 (warm.used_work, warm.peak_bytes))
                again = self._warm(make())
                op(again, LocalBudget())
                twice = LocalBudget()
                op(again, twice)
                self.assertEqual((cold.used_work, cold.peak_bytes),
                                 (twice.used_work, twice.peak_bytes))
                work, peak = cold.used_work, cold.peak_bytes
                limits = sorted({1, work // 7, work // 3, work // 2, (3 * work) // 4,
                                 work - 1, work} - {0})
                for lim in limits:
                    stops = [self._did(op, sg, LocalBudget(work_units=lim))
                             for sg in (make(), self._warm(make()))]
                    self.assertEqual(stops[0], stops[1], (name, 'work', lim))
                    self.assertEqual(lim >= work, stops[0][0] == 'done', (name, lim))
                for frac in (0.25, 0.5, 0.9):
                    mb = peak * frac / 2 ** 20
                    stops = [self._did(op, sg, LocalBudget(memory_mb=mb))
                             for sg in (make(), self._warm(make()))]
                    self.assertEqual(stops[0], stops[1], (name, 'memory', frac))
                    self.assertEqual('stop', stops[0][0], (name, frac))

    def test_account_bounds_the_peak(self):
        """On a chain of comb continuations (each part splits at every node)."""
        arm, p, c = triple('comb_chain')
        g3 = load('comb_chain', 'grandchild')

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

    def test_stitching_admits_what_it_builds(self):
        """GPT final review M3: stitch() builds the wrapper of a plain parent inside its
        budget's scope, and admits every graphlet it builds -- its labels and its digest, the
        canonical JSON of each part's strategy and seed -- before building it. With a
        1,000,000-character strategy label: a zero budget stops stitch(g) and stitch(g, c)
        before the wrapper is built (the traced peak before the stop stays far below the
        label), and an unlimited budget's account bounds the traced peak of stitch(g),
        stitch(p, c) with the label on either side, and of a 100,000-base root seed (warm
        free lists, as test_traverse_union does); the answers are the same with a budget."""
        arm, p0, c0 = triple('keep_forbid')

        def strategy(g):
            g = copy.deepcopy(g)
            g.envelope['strategy']['labels']['extra'] = ['A' * 10 ** 6]
            return g
        cases = {'stitch(g)': lambda: (strategy(p0),),
                 'stitch(g, long seed)': lambda: (_long_seed(p0),),
                 'stitch(p, c), p long': lambda: (strategy(p0), copy.deepcopy(c0)),
                 'stitch(p, c), c long': lambda: (copy.deepcopy(p0), strategy(c0)),
                 'stitch(p, c), long seed': lambda: (_long_seed(p0), copy.deepcopy(c0))}
        for name, make in cases.items():
            with self.subTest(name):
                plain = stitch(*make())
                b = LocalBudget()
                got = stitch(*make(), budget=b)
                self.assertEqual((got.digest, got.labels, got.identity),
                                 (plain.digest, plain.labels, plain.identity))
                args = [make(), make()]
                self.assertLessEqual(_peak(lambda: stitch(*args.pop())), b.peak_bytes)
                zero = LocalBudget(work_units=0, memory_mb=0)
                with self.assertRaises(LocalBudgetExceeded):
                    stitch(*make(), budget=zero)
                self.assertEqual(zero.call_usage(), {'work_units': 0, 'memory_bytes': 0})
                args = [make(), make()]
                self.assertLess(_peak(lambda: stitch(*args.pop(), budget=LocalBudget(
                    work_units=0, memory_mb=0))), 100000)
        # a stitched parent with no child is returned as it is: nothing is built or charged
        sg = stitch(p0, c0)
        zero = LocalBudget(work_units=0, memory_mb=0)
        self.assertIs(stitch(sg, budget=zero), sg)

    def test_spell_charges_its_seed(self):
        """spell(with_seed=True) joins the root seed to the flank (natural orientation) or its
        reversal before it (walking order, left arm) and charges them before it spells:
        GPT final review M2's sibling, found by the lane's sweep. With a 100,000-base root
        seed, on both arms of the root alone and of a junction, in both orientations: the
        account bounds the traced peak (warm free lists), and the text is unchanged."""
        arm, p0, c0 = triple('keep_forbid')
        for name, make in (('fork', lambda: stitch(_long_seed(T.graphlet('fork')))),
                           ('keep_forbid', lambda: stitch(_long_seed(p0), copy.deepcopy(c0)))):
            sg = make()
            for side in sg.arms():
                key = sg.walks(side, by='id')[0].key
                for orientation in ('natural', 'walk'):
                    with self.subTest(name, side=side, orientation=orientation):
                        b = LocalBudget()
                        text = sg.spell(side, key, orientation, True, budget=b)
                        self.assertEqual(text, make().spell(side, key, orientation, True))
                        self.assertIn(sg.seed.sequence if orientation == 'natural' or
                                      side == 'right' else sg.seed.sequence[::-1], text)
                        objs = [make(), make()]
                        self.assertLessEqual(_peak(lambda: objs.pop().spell(
                            side, key, orientation, True)), b.peak_bytes)

    def test_memory_counts_the_caches(self):
        arm, p, c = triple('comb_chain')
        sg = stitch(stitch(p, c), load('comb_chain', 'grandchild'))
        before = sg.memory_bytes()
        sg.walks(arm, by='id')
        sg.evidence(arm)
        sg.summary()
        sg.to_fasta()
        after = sg.memory_bytes()
        self.assertGreater(after, before)
        self.assertGreaterEqual(after, K.models_memory_bytes([x.graphlet for x in sg.parts]))
        self.assertGreaterEqual(after, ops._deep_bytes([sg._cache, sg._jcore]))

    def test_junctions_are_computed_once(self):
        """Each junction's core is computed once per pair of parts and carried into every
        later stitch: adding J children one at a time costs J computations, charged."""
        calls = []
        core = StitchedGraphlet._core

        def counting(self, *a, **kw):
            calls.append(a[2:4])
            return core(self, *a, **kw)
        arm, p, c = triple('keep_forbid')
        s = load('keep_forbid', 'sibling')
        with mock.patch.object(StitchedGraphlet, '_core', counting):
            sg = stitch(stitch(p, c), s)
            self.assertEqual(2, len(calls))
            stitch(sg, load('keep_forbid', 'child'))       # a no-op
            self.assertEqual(2, len(calls))
            g3 = load('comb_chain', 'grandchild')
            _, cp, cc = triple('comb_chain')
            del calls[:]
            stitch(stitch(cp, cc), g3)
            self.assertEqual(2, len(calls))
        b = LocalBudget()
        stitch(stitch(p, c), s, budget=b)
        one = LocalBudget()
        stitch(p, c, budget=one)
        self.assertGreater(b.used_work, one.used_work)


if __name__ == '__main__':
    unittest.main()
