"""The cross-pair union of one seed's graphlets (docs/DESIGN-graphlet-merge.md §2, phase
1a-U): Appendix E.1's P1-P12 on the committed fixtures, the split-graph fixtures under
documents/union/ (scripts/traversal/graphlet_fixtures.py UNION_FIXTURES) and hand-built
absent records in the shapes of SPEC-labeled-traversal-core.md §6.1."""

import copy
import json
import os
import random
import unittest

from tests import traverse_testlib as T  # noqa: F401  (puts the library on sys.path)

from metagraph.traverse import (BadSelector, CombineRefused, LocalBudget, LocalBudgetExceeded,
                                MissingEnvelope, TraverseResponse, union_pairs)
from metagraph.traverse import export as X
from metagraph.traverse import ops
from metagraph.traverse import union as U
from metagraph.traverse._codec import GraphletFormatError

UNION = os.path.join(T.DOCS, 'union')
WITH_BASES = [n for n in sorted(T.DOCUMENTS) if n not in ('caps', 'limits')]
# the fixtures that share the merge locus's seed (annotate, merge, merge_ties)
LOCUS = ('annotate', 'merge', 'merge_ties')


def union_doc(name):
    """A union fixture (documents/union/<name>.graphlet.json) as a client sees it."""
    from metagraph.traverse import from_response
    with open(os.path.join(UNION, name + '.graphlet.json'), encoding='utf-8') as f:
        resp = json.load(f)
    return from_response(resp['results'][0], resp)


def stated(**parts):
    """{pair, graphlet} for each tag=fixture name (a committed fixture, or 'union/<name>')."""
    out = []
    for tag, name in parts.items():
        g = union_doc(name[6:]) if name.startswith('union/') else T.graphlet(name)
        out.append({'pair': tag, 'graphlet': g})
    return out


def locus(**kw):
    return union_pairs(stated(gA='annotate', gM='merge', gT='merge_ties'), **kw)


def absent_result(kind, seed_len=30, **extra):
    """A seed result without a graphlet in the shapes of SPEC §6.1 steps 2-4."""
    base = {'seed': {'seed_id': 'x', 'length_bp': seed_len, 'labels_from_seed': True},
            'limitations': []}
    if kind == 'not_in_graph':
        base.update({'error': 'Seed is not fully present in the graph; graph runs: 0-3',
                     'not_in_graph': {'kmers': 16, 'kmers_present': extra.get('present', 0)},
                     'outcome': {'walks': 'not_in_graph', 'branch_diagnostics': 'complete',
                                 'label_evidence': 'complete', 'delivery': 'inline'}})
    elif kind == 'labels_not_in_graph':
        base.update({'unknown_labels': [{'name': 'x.fa', 'reason': 'not_in_graph'}],
                     'outcome': {'walks': 'labels_not_in_graph'}})
    elif kind == 'failed':
        base.update({'error': 'no label carries the whole seed',
                     'limitations': [{'kind': 'derivation', 'cause': 'no_carrier',
                                      'knob': 'seeds[].sequence'}],
                     'outcome': {'walks': 'failed'}})
    return base


def admitted_paths(gfa):
    """Every maximal path a GFA file admits (oriented walks along its L lines from any
    segment with no incoming link to a segment with no outgoing one), spelled."""
    seqs, out_links, has_in = {}, {}, set()
    for line in gfa.split('\n'):
        f = line.split('\t')
        if f[0] == 'S':
            seqs[f[1]] = f[2]
        elif f[0] == 'L':
            out_links.setdefault(f[1], []).append(f[3])
            has_in.add(f[3])
    paths = []

    def walk(name, acc):
        nxt = out_links.get(name, [])
        if not nxt:
            paths.append(acc)
            return
        for n in nxt:
            walk(n, acc + seqs[n])

    for name in seqs:
        if name not in has_in:
            walk(name, seqs[name])
    return paths


class TestOnePart(unittest.TestCase):
    """P2: the union of one part is that part."""

    def test_walks_are_the_parts_walks(self):
        for name in WITH_BASES:
            g = T.graphlet(name)
            u = union_pairs([g])
            for side in g.arms:
                parts = {g.spell(side, w.path_id): w for w in g.walks(side, by='id')}
                kinds = [w.kind for w in u.walks(side, by='id', stubs=True)]
                self.assertEqual(kinds.count('walk'), len(parts), (name, side))
                for w in u.walks(side, by='id', stubs=True):
                    sp = u.carriers(side, w.id)[0]
                    if w.kind == 'walk':
                        pw = parts[u.spell(side, w.id)]
                        self.assertEqual((sp.how, sp.walk, sp.to_bp), ('whole', pw.path_id,
                                                                       pw.length_bp))
                        self.assertEqual(sp.end_reasons['labels'], pw.end_reasons)
                        self.assertEqual(sp.end_reasons['path_reason'], pw.path_reason)
                        self.assertEqual(sp.complete, pw.complete)
                        self.assertEqual(sp.certified_to_bp,
                                         pw.length_bp - pw.beyond_certified_bp)
                        self.assertEqual({x.label.ref for x in sp.labels if x.from_bp == 0},
                                         {l.ref for l in pw.labels_full}, (name, side, w.id))
                        self.assertEqual(w.whole_in, (u.parts[0].tag,))
                    else:
                        self.assertEqual(sp.how, 'merges')
                        self.assertEqual(sp.route.how, 'held')
                        self.assertEqual(sp.route.to_bp, w.length_bp)
                        self.assertTrue(sp.route.goes_on)
                        self.assertEqual(w.held_in, (u.parts[0].tag,))

    def test_coverage_is_the_arm_certificate(self):
        for name in WITH_BASES:
            g = T.graphlet(name)
            u = union_pairs([g])
            for side, arm in g.arms.items():
                c = u.coverage(side, at_bp=arm.complete_to_bp)
                self.assertIn(c.per_pair[u.parts[0].tag]['category'],
                              ('complete', 'complete_qualified'))
                self.assertEqual(c.per_pair[u.parts[0].tag]['certified_to_bp'],
                                 arm.complete_to_bp)

    def test_no_bases_refused(self):
        with self.assertRaises(CombineRefused) as cm:
            union_pairs([{'pair': 'g', 'graphlet': T.graphlet('caps')}])
        self.assertEqual(cm.exception.code, 'no_bases')
        self.assertEqual(str(cm.exception), U.TEXTS['no_bases'] % 'g')

    def test_spell_equals_the_parts_spelling(self):
        for name in WITH_BASES:
            g = T.graphlet(name)
            u = union_pairs([g])
            for side in g.arms:
                for w in u.walks(side, by='id'):
                    sp = u.carriers(side, w.id)[0]
                    for o in ('natural', 'walk'):
                        for s in (False, True):
                            self.assertEqual(u.spell(side, w.id, o, s),
                                             g.spell(side, sp.walk, o, s))


class TestTrie(unittest.TestCase):
    """P5 on the merge locus: three pairs (annotate, merge, merge_ties) of one seed."""

    def test_segments_break_where_pairs_change_and_children_in_base_order(self):
        u = locus()
        for side in u.arms():
            segs = {s.id: s for s in u.segments(side)}
            t = u._arms[side]
            for us in t.useg:
                firsts = [t.useg[c].first for c in us.children]
                self.assertEqual(firsts, sorted(firsts, key='ACGT'.index))
                s = segs[us.id]
                for c in s.children:
                    self.assertEqual((segs[c].parent, segs[c].from_bp), (s.id, s.to_bp))
                if len(us.children) == 1:
                    c = t.useg[us.children[0]]
                    self.assertTrue(set(c.cover) != set(us.cover) or us.end in t.walk_end
                                    or us.end in t.stub_end or us.end in t.merge_at
                                    or us.end in t.zero_at, (side, us.id))

    def test_every_union_walk_is_whole_or_merges_somewhere(self):
        u = locus()
        for side in u.arms():
            for w in u.walks(side, by='id', stubs=True):
                hows = {s.how for s in u.carriers(side, w.id)}
                self.assertTrue(hows & {'whole', 'merges'}, (side, w.id, hows))

    def test_merge_span_labels_on_the_bubble(self):
        # merge.mgt: at 40 bp of walk 0, a.fa carries [0, 40) and b.fa from 36; the stub
        # [0, 2] has b.fa in continues
        g = T.graphlet('merge')
        u = union_pairs([g])
        walk = next(w for w in u.walks('right', by='id') if w.kind == 'walk')
        c = u.carriers('right', walk.id, at_bp=40)
        got = {x.label.name: x.from_bp for x in c[0].labels}
        self.assertEqual(got['a.fa'], 0)
        self.assertEqual(got['b.fa'], 36)
        stub = next(w for w in u.walks('right', by='id', stubs=True)
                    if w.kind == 'stub' and w.length_bp == 36)
        sp = u.carriers('right', stub.id)[0]
        self.assertEqual(sp.how, 'merges')
        self.assertIn('b.fa', [x.name for x in sp.continues])

    def test_merges_not_diverges_and_route_held(self):
        u = locus()
        # the C allele at 62: merge_ties displays it, merge holds it through its merge
        walks = {w.id: w for w in u.walks('right', by='id', stubs=True)}
        for w in walks.values():
            for s in u.carriers('right', w.id):
                if s.how == 'merges':
                    self.assertIsNone(s.absence)
                    self.assertIsNotNone(s.route)
                    if s.route.how == 'held':
                        self.assertIn(s.pair, w.held_in)
                        self.assertNotIn(s.pair, w.whole_in)

    def test_terminal_marker_keeps_a_shorter_walk(self):
        # gA (annotate, radius 70) ends at 70 inside the 100 bp walks of the others: its
        # walk is kept as a union walk by its terminal marker
        u = locus()
        w70 = [w for w in u.walks('right', by='id') if w.length_bp == 70]
        self.assertEqual(len(w70), 1)
        self.assertEqual(w70[0].whole_in, ('gA',))
        self.assertIn(w70[0].id, [w.id for w in u.walks('right', pairs=['gA'])])
        self.assertEqual(len(u.walks('right', pairs=['gA'])),
                         len(union_pairs(stated(gA='annotate')).walks('right')))

    def test_ends_with_its_reason_and_absence(self):
        u = locus()
        w100 = [w for w in u.walks('right', by='id') if w.length_bp == 100]
        for w in w100:
            s = {x.pair: x for x in u.carriers('right', w.id)}['gA']
            self.assertIn(s.how, ('ends', 'merges'))
            if s.how == 'ends':
                self.assertEqual((s.to_bp, s.end_reasons['path_reason'], s.absence),
                                 (70, 'max_extension_bp', 'uncertified'))
            else:
                self.assertEqual((s.route.how, s.route.to_bp, s.route.absence),
                                 ('ends', 70, 'uncertified'))

    def test_arm_not_walked_is_none(self):
        u = locus()
        left = u.walks('left', by='id')
        for s in u.carriers('left', left[0].id):
            if s.pair in ('gM', 'gT'):
                self.assertEqual((s.how, s.to_bp, s.absence), ('none', 0, 'under_strategy'))

    def test_stub_walks_left_out_unless_asked(self):
        u = locus()
        all_ = u.walks('right', by='id', stubs=True)
        self.assertEqual([w for w in all_ if w.kind == 'walk'], u.walks('right', by='id'))
        fa = u.to_fasta('right')
        self.assertNotIn('kind=stub', fa)
        self.assertIn('kind=stub', u.to_fasta('right', stubs=True))

    def test_stub_end_on_two_merges_refused(self):
        g = copy.deepcopy(T.graphlet('merge'))
        a = g.arms['right']
        # make the stub end (segment 2, a non-first parent of merge 3) a non-first parent of
        # merge 6 too
        a.segments[2].children.append(6)
        a.segments[6].parents = a.segments[6].parents + (2,)
        with self.assertRaises(GraphletFormatError):
            union_pairs([g])


class TestExports(unittest.TestCase):
    """P9 (source-supported mode)."""

    def test_fasta_spellings_and_tokens(self):
        u = locus(expected_pairs=['gA', 'gM', 'gT', 'gX'])
        fa = u.to_fasta()
        recs = fa.strip().split('\n')
        for head, seq in zip(recs[::2], recs[1::2]):
            arm, wid = head[1:].split(' ')[0].split('.u')
            side = 'left' if arm == 'left' else 'right'
            self.assertEqual(seq, u.spell(side, int(wid), with_seed=True))
            self.assertIn('scope=within_one_graph', head)
            self.assertIn('not_included=gX', head)
        self.assertIn('modes=gA:annotate,gM:constrain,gT:constrain', fa)
        self.assertIn('ends=gA@70:max_extension_bp/uncertified', fa)
        self.assertIn('labels=gA:2+/', fa)          # a cut list: a lower bound

    def test_escaping_round_trips(self):
        for name in ('café', 'a,b;c=d@e/f:g>h~i+j|k[l]m$n%o p', 'tab\tx'):
            esc = X.union_escape(name)
            self.assertTrue(all('!' <= ch <= '~' for ch in esc), esc)
            for ch in ',;=@/:>~+|[]$ ':
                self.assertNotIn(ch, esc)
            self.assertEqual(X.union_unescape(esc), name)
        self.assertEqual(X.union_escape('café'), 'caf%C3%A9')

    def test_gfa_admits_only_paths_one_pair_spells(self):
        u = locus()
        gfa = u.to_gfa()
        for line in gfa.split('\n'):
            self.assertTrue(all(' ' <= ch <= '~' or ch == '\t' for ch in line), line)
        spellings = []
        for p in u.parts:
            g = p.graphlet
            for side in g.arms:
                for w in g.walks(side):
                    spellings.append((side, g.spell(side, w.path_id, 'natural', True)))
                for w in u.walks(side, stubs=True, by='id') if side in u.arms() else []:
                    if p.tag in w.held_in or any(s.pair == p.tag and s.how in ('merges',)
                                                 for s in w.spans):
                        spellings.append((side, u.spell(side, w.id, 'natural', True)))
        seed = u.seed.sequence
        for path in admitted_paths(gfa):
            ok = any(path in s for _, s in spellings)
            self.assertTrue(ok, path)
            # no path joins a left flank to a right flank through the seed
            self.assertFalse(path.count(seed) > 1)

    def test_gfa_seed_once_per_arm_and_stub_tags(self):
        u = locus()
        gfa = u.to_gfa()
        self.assertIn('S\tseed_l\t', gfa)
        self.assertIn('S\tseed_r\t', gfa)
        self.assertIn('\tmg:Z:', gfa)
        self.assertIn('\tmt:Z:s', gfa)
        self.assertTrue(gfa.startswith('H\tVN:Z:1.0\n# union_pairs: seed 30 bp'))
        self.assertIn('# scope: ' + U.SCOPE, gfa)


if __name__ == '__main__':
    unittest.main()
