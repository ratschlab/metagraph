"""The cross-pair union of one seed's graphlets (docs/DESIGN-graphlet-merge.md §2, phase
1a-U): Appendix E.1's P1-P12 on the committed fixtures, the split-graph fixtures under
documents/union/ (scripts/traversal/graphlet_fixtures.py UNION_FIXTURES) and hand-built
absent records in the shapes of SPEC-labeled-traversal-core.md §6.1."""

import copy
import dataclasses
import json
import os
import random
import re
import unittest

from tests import traverse_testlib as T  # noqa: F401  (puts the library on sys.path)

from metagraph.traverse import (BadSelector, CombineRefused, LocalBudget, LocalBudgetExceeded,
                                MissingEnvelope, TraverseResponse, union_pairs)
from metagraph.traverse import export as X
from metagraph.traverse import derive, ops
from metagraph.traverse import union as U
from metagraph.traverse.union import _core_sentence as _core
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


def long_seed(g, n):
    """A copy of graphlet |g| whose seed is |n| more bases long (A's before it): the arms and
    their bases unchanged, as an engine answer for a long seed with short flanks."""
    g = copy.deepcopy(g)
    s = 'A' * n + g.seed.sequence
    g.seed = dataclasses.replace(g.seed, sequence=s, length_bp=len(s),
                                 num_kmers=len(s) - g.k + 1)
    return g


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


def admitted_paths(gfa, names=False):
    """Every maximal path a GFA file admits (oriented walks along its L lines from any
    segment with no incoming link to a segment with no outgoing one), spelled; with
    |names|, (spelling, [segment names]) pairs."""
    seqs, out_links, has_in = {}, {}, set()
    for line in gfa.split('\n'):
        f = line.split('\t')
        if f[0] == 'S':
            seqs[f[1]] = f[2]
        elif f[0] == 'L':
            out_links.setdefault(f[1], []).append(f[3])
            has_in.add(f[3])
    paths = []

    def walk(name, acc, via):
        nxt = out_links.get(name, [])
        if not nxt:
            paths.append((acc, via) if names else acc)
            return
        for n in nxt:
            walk(n, acc + seqs[n], via + [n])

    for name in seqs:
        if name not in has_in:
            walk(name, seqs[name], [name])
    return paths


def gfa_segment_pairs(gfa):
    """segment name -> (its bases, the pairs of its pp tag)."""
    out = {}
    for line in gfa.split('\n'):
        f = line.split('\t')
        if f[0] == 'S':
            pp = next((x[5:] for x in f[3:] if x.startswith('pp:Z:')), '')
            out[f[1]] = (f[2], [X.union_unescape(t) for t in pp.split(',') if t])
    return out


def check_gfa_paths(test, u, gfa):
    """P9, by the paths the file admits: every admitted path, and every part of it from the
    seed out to one of its segments, is spelled whole by each pair spelling that segment:
    a root-to-leaf spelling of that pair's own segment DAG (its walks and its routes
    through merges, the seed included) in that arm starts with it (right arm) or ends with
    it (left arm, written on the seed strand). No path holds both arms' seed."""
    segs = gfa_segment_pairs(gfa)
    seed = u.seed.sequence
    refs = {}
    for p in u.parts:
        for side in p.graphlet.arms:
            refs[p.tag, side] = [seed + w if side == 'right' else w[::-1] + seed
                                 for w in dag_spellings(p.graphlet, side)]
    paths = admitted_paths(gfa, names=True)
    test.assertTrue(paths)
    for spelled, via in paths:
        test.assertFalse('seed_l' in via and 'seed_r' in via, via)
        side = 'right' if 'seed_r' in via else 'left'
        for k, name in enumerate(via):
            if name in ('seed_l', 'seed_r'):
                continue
            part = ''.join(segs[x][0] for x in (via[:k + 1] if side == 'right' else via[k:]))
            for tag in segs[name][1]:
                ok = any(r.startswith(part) if side == 'right' else r.endswith(part)
                         for r in refs[tag, side])
                test.assertTrue(ok, (tag, side, via[:k + 1] if side == 'right' else via[k:]))


def limits_with_bases():
    """limits.mgt (no bases) with bases filled in: a zero-length merge leaf that max_steps
    stopped, entered from two parents."""
    g = copy.deepcopy(T.graphlet('limits'))
    for seg, w in zip(g.arms['right'].segments, ['ACG', 'CAAA', 'GAAA', '']):
        assert len(w) == seg.length_bp
        seg.walk = w
    return g


def one_part_graphlets():
    """(name, graphlet, union input) of each committed fixture with bases, limits.mgt's
    with bases filled in (it names no graph, so it is passed with a pair name)."""
    out = [(name, T.graphlet(name)) for name in WITH_BASES]
    out = [(name, g, g) for name, g in out]
    g = limits_with_bases()
    return out + [('limits+bases', g, {'pair': 'L', 'graphlet': g})]


class TestOnePart(unittest.TestCase):
    """P2: the union of one part is that part."""

    def test_walks_are_the_parts_walks(self):
        stubs_to_a_leaf = 0
        for name, g, given in one_part_graphlets():
            u = union_pairs([given])
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
                        self.assertEqual(w.held_in, (u.parts[0].tag,))
                        segs = g.arms[side].segments
                        at_leaf = [pw for pw in g.walks(side, by='id')
                                   if pw.leaf in segs[sp.stub_end].children
                                   and not segs[pw.leaf].length_bp]
                        if sp.route.goes_on:
                            self.assertEqual(at_leaf, [], (name, side, w.id))
                            continue
                        # a stub into a zero-length merge leaf: the route ends at the part
                        # leaf and takes the part walk's own end and certificate there
                        stubs_to_a_leaf += 1
                        pw = at_leaf[0]
                        self.assertEqual(sp.route.end_reasons,
                                         {'path_reason': pw.path_reason,
                                          'labels': pw.end_reasons}, (name, side, w.id))
                        self.assertEqual((sp.route.complete, sp.route.certified_to_bp),
                                         (pw.complete, pw.length_bp - pw.beyond_certified_bp))
                        self.assertEqual((sp.complete, sp.certified_to_bp),
                                         (pw.complete, pw.length_bp - pw.beyond_certified_bp))
        self.assertGreater(stubs_to_a_leaf, 0)

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

    def test_siblings_with_one_first_base_refused(self):
        g = copy.deepcopy(T.graphlet('merge'))
        segs = g.arms['right'].segments
        segs[2].walk = segs[1].walk[0] + segs[2].walk[1:]   # 1 and 2 both start with 'A'
        with self.assertRaises(GraphletFormatError) as cm:
            union_pairs([{'pair': 'g', 'graphlet': g}])
        self.assertIn("segments 1 and 2 of the right arm both start with 'A' at 20 bp",
                      str(cm.exception))

    def test_no_bases_refused(self):
        with self.assertRaises(CombineRefused) as cm:
            union_pairs([{'pair': 'g', 'graphlet': T.graphlet('caps')}])
        self.assertEqual(cm.exception.code, 'no_bases')
        self.assertEqual(str(cm.exception), U.TEXTS['no_bases'] % 'g')

    def test_spell_equals_the_parts_spelling(self):
        for name, g, given in one_part_graphlets():
            u = union_pairs([given])
            for side in g.arms:
                for w in u.walks(side, by='id'):
                    sp = u.carriers(side, w.id)[0]
                    for o in ('natural', 'walk'):
                        for s in (False, True):
                            self.assertEqual(u.spell(side, w.id, o, s),
                                             g.spell(side, sp.walk, o, s))


class TestRoutes(unittest.TestCase):
    """§2.4 'Routes through merges': a route that ends at a part leaf takes the part walk's
    own complete and certificate there."""

    def test_a_route_states_its_own_label_cut(self):
        # GPT review M2: a route through a merge reads the labels where it ends, and its
        # cut flag is that list's, not the stub's: carriers, FASTA labels=
        a = T.graphlet('annotate')
        for seg in a.arms['right'].segments[:6]:
            for p in seg.presence:
                p.total = len(p.labels)     # uncut before the merged segment 6
        b = copy.deepcopy(a)
        b.arms['right'].segments[6].parents = (4, 5)    # b displays the other parent
        b.arms['right'].segments[6].partition.reverse()
        self.assertTrue(ops.label_intervals_cut(a, 'right', 6))
        u = union_pairs([{'pair': 'A', 'graphlet': a}, {'pair': 'B', 'graphlet': b}])
        (w,) = [w for w in u.walks('right') if w.length_bp == 70 and w.whole_in == ('B',)]
        self.assertIn('labels 2+, 1 of them from 0 bp', str(u.carriers('right', w.id)).split(
            '\n')[0])
        self.assertIn('labels=A:2+/1,B:2+/1', u.to_fasta('right', walks=[w.id]))
        sa, sb = u.carriers('right', w.id)
        self.assertEqual((sa.how, sa.route.how, sa.route.to_bp, sa.labels_cut,
                          getattr(sa.route, 'labels_cut', None), len(sa.route.labels)),
                         ('merges', 'held', 70, False, True, 2))
        # GFA states each union segment's own cut: each pair's merged segment (62-70 bp)
        gfa = u.to_gfa().split('\n')
        merged = [x for x in u.segments('right') if (x.from_bp, x.to_bp) == (62, 70)]
        self.assertEqual(sorted(sup.pair for x in merged for sup in x.support), ['A', 'B'])
        for x in merged:
            (sup,) = x.support
            line = next(y for y in gfa if y.startswith('S\ts%d\t' % x.id))
            self.assertIn('pl:Z:%s=2+' % sup.pair, line)
        # constrain: the merged segment's entry list is cut, the stub's is not
        parts = stated(gA='annotate', gT='merge_ties')
        g = parts[1]['graphlet']
        g.arms['right'].segments[6].entry_total += 1
        u = union_pairs(parts)
        (w,) = [w for w in u.walks('right') if w.length_bp == 70]
        self.assertIn('gT:3+/0', u.to_fasta('right', walks=[w.id]))
        s = next(s for s in u.carriers('right', w.id) if s.pair == 'gT')
        self.assertEqual((s.how, s.route.how, s.labels_cut,
                          getattr(s.route, 'labels_cut', None)),
                         ('merges', 'held', False, True))
        # at a depth before the merged segment the route's labels are the stub's, uncut
        s = next(s for s in u.carriers('right', w.id, at_bp=40) if s.pair == 'gT')
        self.assertIs(getattr(s.route, 'labels_cut', None), False)

    def test_labels_before_a_zero_length_leaf(self):
        # GPT review m1: a walk ending at a zero-length merge leaf (7 bp); its labels at a
        # depth before the leaf are those a leaf with a base there gives (the leaf lies
        # beyond the cut either way), never the zero-length leaf's empty intervals
        from tests.test_traverse_union_stitched import _zero_leaf_root

        def at_depths(z):
            g = _zero_leaf_root(z_bases=z)
            u = union_pairs([{'pair': 'A', 'graphlet': g}])
            (w,) = [w for w in u.walks('right') if u.spell('right', w.id) == 'ACGCAAC' + z]
            return g, u, w, [[(x.label.name, x.from_bp) for x in
                              u.carriers('right', w.id, at_bp=at)[0].labels]
                             for at in range(1, 8)]

        g, u, w, zero = at_depths('')
        self.assertEqual(zero, at_depths('T')[3])
        self.assertTrue(all(zero))
        # at 1 bp: the labels alive on the first base
        self.assertEqual(sorted(n for n, _ in zero[0]),
                         sorted(g.labels[i].name for i in ops.alive_at(g, 'right', 0, 0)))
        fa = u.to_fasta('right', walks=[w.id], at_bp=1).split('\n')[0]
        self.assertIn('labels=A:3/3', fa)

    def test_a_route_held_to_a_leaf_a_resource_stopped(self):
        m = copy.deepcopy(T.graphlet('merge'))
        (leaf,) = [s for s in m.arms['right'].segments if s.leaf is not None]
        leaf.leaf.path_reason = 'M'                     # time_budget: a resource stop
        (own,) = m.walks('right')
        self.assertEqual((own.length_bp, own.path_reason, own.complete),
                         (100, 'time_budget', False))
        u = union_pairs([{'pair': 'gM', 'graphlet': m},
                         {'pair': 'gT', 'graphlet': T.graphlet('merge_ties')}])
        seen = 0
        for w in u.walks('right', by='id', stubs=True):
            s = {x.pair: x for x in u.carriers('right', w.id)}['gM']
            if s.how != 'merges' or s.route.goes_on:
                continue
            # the route is held to gM's leaf: gM's own walk there is incomplete
            seen += 1
            self.assertEqual((s.route.how, s.route.to_bp, s.route.end_reasons['path_reason'],
                              s.route.complete, s.complete),
                             ('held', 100, 'time_budget', False, False))
        self.assertEqual(seen, 1)


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

    def test_fasta_labels_at_a_depth(self):
        # at_bp= gives the labels= counts at that depth, as carriers(at_bp=) lists them
        u = locus()
        for at in (None, 40, 10):
            fa = u.to_fasta('right', at_bp=at)
            for head in fa.split('\n')[::2]:
                if not head:
                    continue
                wid = int(head.split(' ')[0].split('.u')[1])
                tok = next(t for t in head.split(' ') if t.startswith('labels='))
                want = []
                for s in u.carriers('right', wid, at_bp=at):
                    w = u.walks('right', stubs=True, by='id')[wid]
                    if s.pair not in w.held_in:
                        continue
                    labs = s.route.labels if s.how == 'merges' else s.labels
                    want.append('%s:%d%s/%d' % (s.pair, len(labs), '+' if s.labels_cut else '',
                                                sum(1 for x in labs if x.from_bp == 0)))
                self.assertEqual(tok, 'labels=' + ','.join(want), (at, wid))
        self.assertNotEqual(u.to_fasta('right', at_bp=10), u.to_fasta('right'))

    def test_escaping_round_trips(self):
        for name in ('café', 'a,b;c=d@e/f:g>h~i+j|k[l]m$n%o p', 'tab\tx'):
            esc = X.union_escape(name)
            self.assertTrue(all('!' <= ch <= '~' for ch in esc), esc)
            for ch in ',;=@/:>~+|[]$ ':
                self.assertNotIn(ch, esc)
            self.assertEqual(X.union_unescape(esc), name)
        self.assertEqual(X.union_escape('café'), 'caf%C3%A9')

    def test_gfa_admits_only_paths_one_pair_spells(self):
        for u in (locus(), split_union('annotate', 'keep')[3], split_union('constrain',
                                                                            'merge')[3]):
            gfa = u.to_gfa()
            for line in gfa.split('\n'):
                self.assertTrue(all(' ' <= ch <= '~' or ch == '\t' for ch in line), line)
            check_gfa_paths(self, u, gfa)

    def test_g4_must_2_no_cross_source_path(self):
        """Appendix G.4, Must 2: A spells AACGT and, through a merge, TACGT; B spells AACGA.
        The file admits no TACGA. Here A is merge.mgt (its stub at 20 bp, 'T...', merges at
        36 bp into the segment its walk 'A...' reaches) and B the same graph with that stub
        starting 'C' and the walk's last base changed: a stub linked to the merge's union
        segment would admit A's stub followed by B's last base."""
        a = T.graphlet('merge')
        b = copy.deepcopy(a)
        segs = b.arms['right'].segments
        self.assertEqual((segs[2].walk[0], segs[3].parents), ('T', (1, 2)))
        segs[2].walk = 'C' + segs[2].walk[1:]
        last = segs[6].walk[-1]
        segs[6].walk = segs[6].walk[:-1] + next(c for c in 'ACGT' if c != last)
        u = union_pairs([{'pair': 'A', 'graphlet': a}, {'pair': 'B', 'graphlet': b}])
        gfa = u.to_gfa()
        check_gfa_paths(self, u, gfa)
        sa = a.arms['right'].segments
        tacga = u.seed.sequence + sa[0].walk + sa[2].walk + sa[3].walk + sa[5].walk + \
            segs[6].walk
        self.assertFalse(any(tacga.startswith(p) and len(p) > len(u.seed.sequence) + 36
                             for p in admitted_paths(gfa)))
        self.assertNotIn(tacga, admitted_paths(gfa))
        # A's stub is a tag on its own union segment, never a link
        (stub_line,) = [x for x in gfa.split('\n') if x.startswith('S\t') and '\tmg:Z:A' in x
                        and x.split('\t')[2].startswith(sa[2].walk[:4])]
        name = stub_line.split('\t')[1]
        self.assertFalse(any(x.startswith('L\t%s\t' % name) for x in gfa.split('\n')))

    def test_tokens_write_trie_nodes_and_reasons_as_they_are(self):
        # A.5: a merges token names the trie node as 's<id>' or 's<id>$', and mt:Z the same
        g = limits_with_bases()
        u = union_pairs([{'pair': 'L', 'graphlet': g}])
        fa, gfa = u.to_fasta(stubs=True), u.to_gfa()
        self.assertIn('merges=L@7>s1$', fa)
        self.assertIn('\tmt:Z:s1$', gfa)
        self.assertNotIn('%24', fa + gfa)
        # constrain mode ends= joins its label reasons with a bare '|'
        s = X._span_token(U.PairSpan('A', 'constrain', 'ends', 76, 0, None, None, None, (),
                                     {'path_reason': None,
                                      'labels': {'dead_end': 1, 'max_steps': 2}},
                                     (), False, False, 76, 'graph'))
        self.assertEqual(s, 'A@76:dead_end|max_steps/graph')

    def test_gfa_seed_once_per_arm_and_stub_tags(self):
        u = locus()
        gfa = u.to_gfa()
        self.assertIn('S\tseed_l\t', gfa)
        self.assertIn('S\tseed_r\t', gfa)
        self.assertIn('\tmg:Z:', gfa)
        self.assertIn('\tmt:Z:s', gfa)
        self.assertTrue(gfa.startswith('H\tVN:Z:1.0\n# union_pairs: seed 30 bp'))
        self.assertIn('# scope: ' + U.SCOPE, gfa)

    def test_gfa_full_labels_split_back_into_tag_and_names(self):
        # labels='full': each entry is the tag once, '=', then each name escaped on its
        # own, so a release-qualified tag (with '/') and a name holding the separators
        # decode back to what segments() states
        a, m = T.graphlet('annotate'), T.graphlet('merge')
        a.labels[0].name = 'x,y;z=w/v%'
        u = union_pairs([{'pair': {'graph': 'gA', 'release': 'r1'}, 'graphlet': a},
                         {'pair': {'graph': 'gM', 'release': 'r2'}, 'graphlet': m}])
        self.assertEqual(u.pairs, ('r1/gA', 'r2/gM'))
        want = {}
        for side in u.arms():
            for s in u.segments(side):
                want['s%d' % s.id] = {sup.pair: sorted(li.label.name for li in sup.labels)
                                      for sup in s.support}
        seen = 0
        for line in u.to_gfa(labels='full').split('\n'):
            f = line.split('\t')
            if f[0] != 'S' or f[1].startswith('seed_'):
                continue
            (pl,) = [x[len('pl:Z:'):] for x in f if x.startswith('pl:Z:')]
            got = {}
            for entry in pl.split(';'):
                tag, names = entry.split('=')
                got[X.union_unescape(tag)] = sorted(
                    X.union_unescape(x) for x in names.split(',')) if names else []
            self.assertEqual(got, want[f[1]], line)
            seen += 'x,y;z=w/v%' in got.get('r1/gA', ())
        self.assertGreater(seen, 0)


def response_of(*names, extra_results=(), seeds_sent=None):
    """A TraverseResponse as the client builds it from the fixtures' graphlet responses
    (their results concatenated) and extra seed results without a graphlet."""
    from metagraph.traverse import TraverseClient
    raw = None
    results = []
    for n in names:
        d = T.doc_json(n, 'graphlet')
        if raw is None:
            raw = {k: v for k, v in d.items() if k != 'results'}
        results.extend(d['results'])
    if raw is None:
        raw = {k: v for k, v in T.doc_json('merge', 'graphlet').items() if k != 'results'}
    results.extend(extra_results)
    raw['results'] = results
    r = TraverseClient.response(raw)
    r.seeds_sent = seeds_sent
    return r


def with_env(g, **changes):
    """A shallow copy of |g| whose envelope has |changes| (None removes a key)."""
    h = copy.copy(g)
    env = dict(g.envelope)
    for k, v in changes.items():
        if v is None:
            env.pop(k, None)
        else:
            env[k] = v
    h.envelope = env
    return h


def with_fields(g, **changes):
    h = copy.copy(g)
    for k, v in changes.items():
        setattr(h, k, v)
    return h


def refused(test, code, fn, text=None):
    with test.assertRaises(CombineRefused) as cm:
        fn()
    e = cm.exception
    test.assertEqual(e.code, code, str(e))
    if text is not None:
        test.assertEqual(str(e), text)
    # .least_bytes: the size a retry needs, where the text states one (summary_too_small:
    # the number it states twice), else None (the service's min_max_bytes)
    m = re.fullmatch(r'.* take (\d+) bytes, .* pass max_bytes >= (\d+)', str(e))
    if code == 'summary_too_small':
        test.assertEqual((getattr(e, 'least_bytes', None), int(m.group(2))),
                         (int(m.group(1)), int(m.group(1))))
    else:
        test.assertIsNone(m)
        test.assertIsNone(getattr(e, 'least_bytes', None))
    return str(e)


class TestRefusals(unittest.TestCase):
    """P1: one test per code, with the exact text; check-major order."""

    def test_arguments(self):
        g = T.graphlet('merge')
        refused(self, 'bad_expected', lambda: union_pairs([g], expected_pairs='g1'),
                U.TEXTS['bad_expected_str'])
        refused(self, 'bad_expected', lambda: union_pairs([g], expected_pairs=5),
                U.TEXTS['bad_expected'] % '5')
        refused(self, 'bad_expected', lambda: union_pairs([g], expected_pairs=['g', 3]),
                U.TEXTS['bad_expected'] % 'expected_pairs[1] (3)')
        refused(self, 'bad_seqs_order', lambda: union_pairs([g], seqs_order={'g': 'maybe'}),
                U.TEXTS['bad_seqs_order'] % "seqs_order['g'] = 'maybe'")
        refused(self, 'bad_seqs_order', lambda: union_pairs([g], seqs_order=['g']))
        with self.assertRaises(TypeError):
            union_pairs([g], seed=5)
        # arguments first: a string expected_pairs before a bad element
        refused(self, 'bad_expected', lambda: union_pairs([42], expected_pairs='g'))

    def test_not_combinable_is_a_type_error_and_comes_first(self):
        with self.assertRaises(TypeError) as cm:
            union_pairs([{'pair': 'g', 'graphlet': T.graphlet('caps')}, 42])
        self.assertEqual(cm.exception.code, 'not_combinable')
        self.assertEqual(str(cm.exception), U.TEXTS['not_combinable'] % 'result 1 (a int)')
        for bad in ({'pair': 'g'}, {'graphlet': T.graphlet('merge')},
                    {'pair': 'g', 'graphlet': T.graphlet('merge'), 'error': 'x'},
                    {'pair': 'g', 'result': {}, 'other': 1}, {'pair': 3, 'error': 'x'},
                    {'response': 1, 'index': 0}, [1]):
            refused(self, 'not_combinable', lambda: union_pairs([bad]))

    def test_view_unsupported(self):
        g = T.graphlet('merge')
        refused(self, 'view_unsupported', lambda: union_pairs([g.subgraph([{'id': 0}])]),
                U.TEXTS['view_unsupported'])

    def test_not_one_seed(self):
        two = response_of('merge', 'merge')
        refused(self, 'not_one_seed', lambda: union_pairs([two]),
                U.TEXTS['not_one_seed'] % (0, 2))
        zero = response_of()
        refused(self, 'not_one_seed', lambda: union_pairs([zero]),
                U.TEXTS['not_one_seed_zero'] % 0)
        refused(self, 'not_one_seed', lambda: union_pairs([{'response': two, 'index': 2}]),
                U.TEXTS['not_one_seed_range'] % (0, 2, 2))

    def test_result_has_graphlet_and_no_graphlet(self):
        r = T.doc_json('merge', 'graphlet')['results'][0]
        refused(self, 'result_has_graphlet', lambda: union_pairs([{'pair': 'g', 'result': r}]),
                U.TEXTS['result_has_graphlet'] % 0)
        walked = {k: v for k, v in r.items() if k != 'graphlet'}
        refused(self, 'no_graphlet', lambda: union_pairs([{'pair': 'g', 'result': walked}]),
                U.TEXTS['no_graphlet'] % (0, r['outcome']['walks']))

    def test_no_seed(self):
        a = {'pair': 'g2', 'result': absent_result('not_in_graph')}
        refused(self, 'no_seed', lambda: union_pairs([a]), U.TEXTS['no_seed'])
        u = union_pairs([a], seed=T.graphlet('merge').seed.sequence)
        self.assertEqual(u.arms(), [])

    def test_missing_envelope(self):
        with self.assertRaises(MissingEnvelope) as cm:
            union_pairs([T.body('merge')])
        self.assertTrue(str(cm.exception).startswith('union_pairs() needs the response '
                                                     'envelope'))

    def test_pair_mismatch_and_no_pair(self):
        g = with_env(T.graphlet('merge'), graph='g1')
        refused(self, 'pair_mismatch', lambda: union_pairs([{'pair': 'g2', 'graphlet': g}]),
                U.TEXTS['pair_mismatch'] % (0, U._show(U._env_desc(g.envelope, g)), 'g2'))
        union_pairs([{'pair': 'g1', 'graphlet': g}])          # the pair it was answered on
        nameless = with_fields(T.graphlet('merge'), index_ns=None)
        refused(self, 'no_pair', lambda: union_pairs([nameless]), U.TEXTS['no_pair'] % 0)

    def test_ambiguous_pair(self):
        a = with_fields(with_env(T.graphlet('merge'), graph='g2'), index_fp='f1')
        b = with_fields(with_env(T.graphlet('merge_ties'), graph='g2'), index_fp='f2')
        u = union_pairs([a, b])                   # two pairs, told apart by index_fp
        self.assertEqual(len(u.parts), 2)
        self.assertEqual(len({p.tag for p in u.parts}), 2)
        text = refused(self, 'ambiguous_pair', lambda: union_pairs([a, b],
                                                                   expected_pairs=['g2']))
        self.assertTrue(text.startswith('pair g2 matches both '))

    def test_different_alphabet(self):
        a = T.graphlet('merge')
        b = with_fields(T.graphlet('merge_ties'), alphabet='$ACDEFGHIK')
        refused(self, 'different_alphabet',
                lambda: union_pairs([{'pair': 'a', 'graphlet': a}, {'pair': 'b', 'graphlet': b}]),
                U.TEXTS['different_alphabet'] % ('b', 'ACDEFGHIK', 'a', 'ACGT'))
        # '$' does not count: succinct and hash graphs are one DNA alphabet
        c = with_fields(T.graphlet('merge_ties'), alphabet='ACGT')
        union_pairs([{'pair': 'a', 'graphlet': a}, {'pair': 'c', 'graphlet': c}])

    def test_different_seed_variants(self):
        a = T.graphlet('merge')
        seq = a.seed.sequence
        pa = {'pair': 'a', 'graphlet': a}
        # unequal lengths
        refused(self, 'different_seed', lambda: union_pairs([pa, {'pair': 'b', 'graphlet':
                                                                  T.graphlet('fork')}]),
                U.TEXTS['different_seed_length'] % ("part b's seed", 40, "part a's", 30,
                                                    sum(1 for _ in [0]) - 1 + _common(
                                                        seq, T.graphlet('fork').seed.sequence)))
        # the first difference
        other = T.graphlet('ambiguous_split')
        i = _common(seq, other.seed.sequence)
        refused(self, 'different_seed', lambda: union_pairs([pa, {'pair': 'b',
                                                                  'graphlet': other}]),
                U.TEXTS['different_seed'] % ("part b's seed", "part a's", i))
        # the reverse complement
        rc = with_fields(T.graphlet('merge_ties'), seed=with_fields(
            a.seed, sequence=X_rc(seq)))
        refused(self, 'different_seed', lambda: union_pairs([pa, {'pair': 'b', 'graphlet': rc}]),
                U.TEXTS['different_seed_revcomp'] % ("part b's seed", "part a's", 'part b',
                                                     'part a'))
        # seed=, an absent record's length_bp and a stated sequence
        bad = seq[:5] + ('A' if seq[5] != 'A' else 'C') + seq[6:]
        refused(self, 'different_seed', lambda: union_pairs([pa], seed=bad),
                U.TEXTS['different_seed'] % ('seed=', "part a's", 5))
        union_pairs([pa], seed=seq.lower())            # case mapped
        refused(self, 'different_seed', lambda: union_pairs(
            [pa, {'pair': 'g2', 'result': absent_result('not_in_graph', seed_len=31)}]),
            U.TEXTS['different_seed_length_absent'] % (
                'the seed of the absent record of pair g2', 31, "part a's", 30))
        self.assertNotIn('equal for the first', U.TEXTS['different_seed_length_absent'])
        refused(self, 'different_seed', lambda: union_pairs(
            [pa, {'pair': 'g2', 'result': absent_result('not_in_graph'), 'sequence': bad}]),
            U.TEXTS['different_seed'] % ('the seed of the absent record of pair g2',
                                         "part a's", 5))
        # without parts, seed= is the reference
        refused(self, 'different_seed', lambda: union_pairs(
            [{'pair': 'g2', 'result': absent_result('not_in_graph'), 'sequence': bad}],
            seed=seq), U.TEXTS['different_seed'] % ('the seed of the absent record of pair g2',
                                                    'seed=', 5))

    def test_duplicate_pair(self):
        a, b = T.graphlet('merge'), T.graphlet('merge_ties')
        refused(self, 'duplicate_pair', lambda: union_pairs(
            [{'pair': 'g', 'graphlet': a}, {'pair': 'g', 'graphlet': b}]),
            U.TEXTS['duplicate_pair_unknown'] % (0, 1, 'g'))
        fa, fb = with_fields(a, index_fp='f'), with_fields(b, index_fp='f')
        refused(self, 'duplicate_pair', lambda: union_pairs(
            [{'pair': 'g', 'graphlet': fa}, {'pair': 'g', 'graphlet': fb}]),
            U.TEXTS['duplicate_pair'] % (0, 1, 'g'))
        # one body from two requests (a knob nothing reached)
        env = copy.deepcopy(a.envelope)
        env['strategy']['branching']['max_splits_per_path'] = 128
        a2 = with_fields(a, envelope=env)
        refused(self, 'duplicate_pair', lambda: union_pairs(
            [{'pair': 'g', 'graphlet': a}, {'pair': 'g', 'graphlet': a2}]),
            U.TEXTS['duplicate_pair_request'] % (0, 1, 'g',
                                                 'strategy.branching.max_splits_per_path'))
        # absent records: the remedy names {pair, result}
        refused(self, 'duplicate_pair', lambda: union_pairs(
            [{'pair': 'g', 'result': absent_result('not_in_graph')},
             {'pair': 'g', 'result': absent_result('not_in_graph', present=3)}], seed=a.seed.sequence),
            U.TEXTS['duplicate_pair_unknown'] % (0, 1, 'g'))

    def test_identical_retrieval_is_a_noop_and_envelope_form_wins(self):
        a = with_env(T.graphlet('merge'), graph='g2')
        stated_ = {'pair': 'g2', 'graphlet': a}
        for order in ([a, stated_], [stated_, a]):
            u = union_pairs(order)
            self.assertEqual(len(u.parts), 1)
            self.assertEqual(u.parts[0].pair_stated_by, 'envelope')
        self.assertEqual(union_pairs([a, stated_]).digest, union_pairs([stated_, a]).digest)
        # the same body on two stated pairs: two parts, no model shared or merged
        a = T.graphlet('merge')
        u = union_pairs([{'pair': 'g1', 'graphlet': a}, {'pair': 'g2', 'graphlet': a}])
        self.assertEqual(len(u.parts), 2)
        self.assertIs(u.parts[0].graphlet, u.parts[1].graphlet)
        one = ops.memory_bytes(a)
        self.assertLess(u.memory_bytes(), 2 * one)    # the model given twice counted once

    def test_not_expected(self):
        refused(self, 'not_expected', lambda: locus(expected_pairs=['gA', 'gM']),
                U.TEXTS['not_expected'] % 'gT')

    def test_query_refusals(self):
        u = locus()
        refused(self, 'ambiguous_name', lambda: u.label({'name': 'a.fa'}),
                U.TEXTS['ambiguous_name'] % ('a.fa', 'gA, gM, gT'))
        with self.assertRaises(BadSelector) as cm:
            u.label({'id': 0})
        self.assertEqual(str(cm.exception), U.TEXTS['label_id'])
        full = u.summary(max_bytes=10 ** 6)
        no_arms = ['arms left out (max_bytes)']
        size = min(len(U.C.canonical_json(o)) for o in (
            {'statements': full['statements'], 'pairs': full['pairs'], 'notes': no_arms},
            {'statements': full['statements'], 'notes': ['pairs left out (max_bytes)'] +
             no_arms}))
        refused(self, 'summary_too_small', lambda: u.summary(max_bytes=100),
                U.TEXTS['summary_too_small'] % (size, 100, size))
        refused(self, 'not_offered', u.dump, U.TEXTS['not_offered_dump'])
        refused(self, 'not_offered', u.to_json, U.TEXTS['not_offered_json'])
        refused(self, 'not_offered', u.compare, U.TEXTS['not_offered'] % 'compare')
        refused(self, 'not_offered', lambda: u.walks('right', paths='inferred'),
                U.TEXTS['not_offered_inferred'])
        refused(self, 'not_offered', lambda: u.walks('right', max_transitions=2))
        refused(self, 'not_offered', lambda: u.spell('right', 'i:1.2'))


class TestSummaryBound(unittest.TestCase):
    """summary(max_bytes=) answers at most max_bytes of canonical JSON, its notes included,
    or refuses summary_too_small (GPT review M6; the service turns a larger answer into
    result_too_large)."""

    def sweep(self, u, caps):
        got, least = {}, set()
        for mb in caps:
            try:
                out = u.summary(max_bytes=mb)
            except CombineRefused as e:
                self.assertEqual(e.code, 'summary_too_small', mb)
                self.assertGreater(e.least_bytes, mb)
                least.add(e.least_bytes)
                got[mb] = None
                continue
            n = len(U.C.canonical_json(out).encode('utf-8'))
            self.assertLessEqual(n, mb, (mb, out.get('notes')))
            self.assertIn('statements', out)
            got[mb] = n
        # every refusal names one least size, and summary(max_bytes=) of it answers
        self.assertLessEqual(len(least), 1)
        for n in least:
            self.assertLessEqual(len(U.C.canonical_json(u.summary(max_bytes=n))), n)
        return got

    def boundary(self, u, upto):
        """Every ceiling from 1 to |upto|: refused below the least answer, which is answered
        at exactly the size the refusal names."""
        got = self.sweep(u, range(1, upto))
        least = min(mb for mb, n in got.items() if n is not None)
        self.assertTrue(all(got[mb] is None for mb in range(1, least)))
        self.assertEqual(got[least], least)
        refused(self, 'summary_too_small', lambda: u.summary(max_bytes=least - 1),
                U.TEXTS['summary_too_small'] % (least, least - 1, least))
        with self.assertRaises(CombineRefused) as cm:
            u.summary(max_bytes=1)
        self.assertEqual(cm.exception.least_bytes, least)
        return least

    def test_every_ceiling_on_the_locus(self):
        u = locus()
        self.boundary(u, 2400)
        # the exact boundary of the whole answer: at its size it is whole, one byte less
        # leaves something out and says so
        full = u.summary(max_bytes=10 ** 6)
        self.assertNotIn('notes', full)
        n = len(U.C.canonical_json(full))
        self.assertEqual(u.summary(max_bytes=n), full)
        less = u.summary(max_bytes=n - 1)
        self.assertLessEqual(len(U.C.canonical_json(less)), n - 1)
        self.assertEqual(less['notes'], ['top walks cut (max_bytes)'])

    def test_a_union_without_arms(self):
        # absent records only: the least answer carries the empty arms, not the note that
        # they are left out (Opus fix review N7)
        u = union_pairs([{'pair': 'g%d' % i, 'result': absent_result('not_in_graph')}
                         for i in range(3)], seed='A' * 30)
        self.assertEqual(u.arms(), [])
        least = self.boundary(u, 1500)
        self.assertEqual(u.summary(max_bytes=least).get('arms'), {})

    def test_the_services_sweep_on_r237(self):
        # the service's repro: r237 (2 parts, 9 absent, 11 expected), 2404..6000 bp in steps
        # of 7; at 68e29b5e0 199 of the 514 answered over the limit, by up to 56 bytes
        _, u = r237_union()
        got = self.sweep(u, range(2404, 6001, 7))
        self.assertEqual(len(got), 514)
        self.assertGreater(sum(1 for n in got.values() if n is not None), 400)


def _common(a, b):
    i = 0
    while i < min(len(a), len(b)) and a[i] == b[i]:
        i += 1
    return i


def X_rc(s):
    return s[::-1].translate(str.maketrans('ACGT', 'TGCA'))


class TestInputs(unittest.TestCase):
    """P4: every row of the decision table; pairs, tags and seed checks."""

    def test_forms(self):
        seq = T.graphlet('merge').seed.sequence
        r1 = response_of('merge')
        u = union_pairs([r1])
        self.assertEqual(len(u.parts), 1)
        self.assertEqual(u.parts[0].pair_stated_by, 'envelope')
        # a one-seed response with an error result; {response, index}; {response, sequence}
        err = response_of(extra_results=[absent_result('not_in_graph')],
                          seeds_sent=[{'sequence': seq}])
        u = union_pairs([err], seed=seq)
        self.assertEqual((u.absent[0].category, u.absent[0].seed_check),
                         ('seed_not_in_graph', 'sent'))
        two = response_of('merge', extra_results=[absent_result('labels_not_in_graph')])
        u = union_pairs([{'response': two, 'index': 0}])
        self.assertEqual(len(u.parts), 1)
        u = union_pairs([{'response': two, 'index': 1}], seed=seq)
        self.assertEqual((u.absent[0].category, u.absent[0].seed_check, u.absent[0].pair_stated_by),
                         ('labels_not_in_graph', 'unchecked', 'envelope'))
        one_err = response_of(extra_results=[absent_result('failed')])
        u = union_pairs([{'response': one_err, 'sequence': seq.lower()}], seed=seq)
        self.assertEqual((u.absent[0].category, u.absent[0].seed_check), ('not_walked', 'stated'))
        # {pair, result} with and without sequence, and with the error beside it
        u = union_pairs([{'pair': 'g2', 'result': absent_result('not_in_graph', present=4),
                          'error': 'Seed is not fully present'},
                         {'pair': 'g3', 'result': absent_result('not_in_graph'), 'sequence': seq}],
                        seed=seq)
        self.assertEqual([a.seed_check for a in u.absent], ['unchecked', 'stated'])
        # {pair, error}: a string and an exception
        u = union_pairs([T.graphlet('merge'), {'pair': 'g4', 'error': 'HTTP 503'},
                         {'pair': 'g5', 'error': ValueError('deadline exceeded')}])
        self.assertEqual([(a.tag, a.category, a.detail['error']) for a in u.absent],
                         [('g4', 'failed', 'HTTP 503'), ('g5', 'failed', 'deadline exceeded')])

    def test_absent_only_union_is_armless(self):
        seq = T.graphlet('merge').seed.sequence
        u = union_pairs([{'pair': 'g2', 'result': absent_result('not_in_graph')},
                         {'pair': 'g3', 'result': absent_result('labels_not_in_graph')}],
                        seed=seq, expected_pairs=['g2', 'g3'])
        self.assertEqual((u.arms(), u.walks('right')), ([], []))
        c = u.coverage()
        self.assertTrue(c.sentence.startswith('both arms: seed not in graph: g2 (0 of 16 '
                                              'k-mers); labels not in graph: g3.'))
        self.assertIsNone(c.every_pair_certified_to_bp)

    def test_bare_name_and_envelope_pair_are_one(self):
        g = with_env(T.graphlet('merge'), graph='g2')
        u = union_pairs([g], expected_pairs=['g2', 'g9'])
        self.assertEqual(sorted(u.expected), ['g2', 'g9'])
        self.assertEqual(u._not_included, ('g9',))
        self.assertEqual(u.parts[0].tag, 'g2')

    def test_colliding_tags_get_suffixes(self):
        a = with_env(T.graphlet('merge'), graph='g2', graph_path='/x/a.dbg')
        b = with_env(T.graphlet('merge_ties'), graph='g2', graph_path='/y/b.dbg')
        u = union_pairs([a, b])
        tags = [p.tag for p in u.parts]
        self.assertEqual(len(set(tags)), 2)
        for t in tags:
            self.assertRegex(t, r'^g2~[0-9a-f]{6,}$')
        self.assertEqual({s.pair for s in u.carriers('right', 0)}, set(tags))

    def test_failure_then_answer_keeps_the_answer(self):
        g = with_env(T.graphlet('merge'), graph='g2')
        for order in ([{'pair': 'g2', 'error': 'timeout'}, g],
                      [g, {'pair': 'g2', 'error': 'timeout'}]):
            u = union_pairs(order)
            self.assertEqual((len(u.parts), len(u.absent)), (1, 0))
            prov = u.provenance()
            self.assertIn({'tag': 'g2', 'kind': 'earlier_failure', 'error': 'timeout'}, prov)

    def test_identical_absent_records_are_one(self):
        seq = T.graphlet('merge').seed.sequence
        a = {'pair': 'g2', 'result': absent_result('not_in_graph')}
        b = {'pair': 'g2', 'result': dict(absent_result('not_in_graph'),
                                          seed={'seed_id': 'other', 'length_bp': 30,
                                                'labels_from_seed': True})}
        u = union_pairs([a, b], seed=seq)
        self.assertEqual(len(u.absent), 1)

    def test_identical_absent_records_keep_the_stronger_check_in_either_order(self):
        """Of two forms of one absent record the stronger seed check is kept, whatever the
        order, and the check enters the digest (found by the r237 round trip: a stated and
        an unchecked union had one digest and different sentences)."""
        seq = T.graphlet('merge').seed.sequence
        bare = {'pair': 'g2', 'result': absent_result('not_in_graph')}
        told = dict(bare, sequence=seq.lower())
        got = []
        for order in ([bare, told], [told, bare]):
            u = union_pairs(order, seed=seq)
            got.append((u.digest, u.absent[0].seed_check, u.coverage().sentence))
        self.assertEqual(got[0], got[1])
        self.assertEqual(got[0][1], 'stated')
        self.assertNotEqual(union_pairs([bare], seed=seq).digest, got[0][0])

    def test_failures_are_listed_in_one_order(self):
        g = with_env(T.graphlet('merge'), graph='g2')
        fa, fb = {'pair': 'g2', 'error': 'timeout'}, {'pair': 'g2', 'error': 'HTTP 503'}
        got = [(u.digest, u.provenance()) for u in (union_pairs([fa, g, fb]),
                                                     union_pairs([fb, fa, g]))]
        self.assertEqual(got[0], got[1])
        self.assertEqual([e['error'] for e in got[0][1] if e['kind'] == 'earlier_failure'],
                         ['HTTP 503', 'timeout'])
        self.assertNotEqual(union_pairs([g]).digest, got[0][0])

    def test_a_single_graph_pair_is_named_by_its_tag_everywhere(self):
        # a single-graph server's answer names no graph: its tag is its index_ns, and one
        # rule (A.4) matches that name in expected_pairs=, seqs_order= and against records
        g = T.graphlet('merge')
        ns = g.index_ns
        self.assertIsNone(g.envelope.get('graph'))
        u = union_pairs([g], expected_pairs=[ns, 'g9'], seqs_order={ns: 'verified'})
        self.assertEqual((u.parts[0].tag, u.expected, u._not_included, u.parts[0].seqs_order),
                         (ns, ('g9', ns), ('g9',), 'verified'))
        refused(self, 'not_expected', lambda: union_pairs([g], expected_pairs=['g9']))
        refused(self, 'duplicate_pair', lambda: union_pairs(
            [g, {'pair': ns, 'result': absent_result('not_in_graph')}]))

    def test_two_seqs_order_values_for_one_pair_come_before_the_results(self):
        a, b = T.graphlet('merge'), T.graphlet('merge_ties')
        so = {'g': 'verified', ('g', None, None): 'count_only'}
        # an argument error: reported before the two answers on one pair are
        text = refused(self, 'bad_seqs_order', lambda: union_pairs(
            [{'pair': 'g', 'graphlet': a}, {'pair': 'g', 'graphlet': b}], seqs_order=so))
        self.assertIn("seqs_order['g'] and seqs_order[('g', None, None)] (two values for one "
                      "pair)", text)
        # two keys that only a result's pair joins: found once the pairs are known
        self.assertEqual(union_pairs([{'pair': 'g', 'graphlet': a}],
                                     seqs_order={'g': 'verified'}).parts[0].seqs_order,
                         'verified')
        text = refused(self, 'bad_seqs_order', lambda: union_pairs(
            [{'pair': 'g', 'graphlet': a}], seqs_order={
                ('g', '/x/a.dbg'): 'verified', ('g', '/y/b.dbg'): 'count_only'}))
        self.assertIn('two values for pair g', text)

    def test_the_mapping_status_is_in_the_digest(self):
        # GPT review m2: seqs_order changes the warnings and evidence, so it changes the
        # digest, and a resume token does not cross it
        parts = [{'pair': 'g', 'graphlet': T.graphlet('same_name')}]
        absent = [{'pair': 'h', 'result': absent_result('not_in_graph')}]
        us = {v: union_pairs(parts + absent, seqs_order=None if v is None else {'g': v})
              for v in (None, 'verified', 'count_only')}
        self.assertEqual(len({u.digest for u in us.values()}), 3)
        self.assertEqual(us['verified'].digest,
                         union_pairs(parts + absent, seqs_order={'g': 'verified'}).digest)
        self.assertNotEqual(union_pairs(parts + absent, seqs_order={'h': 'verified'}).digest,
                            us[None].digest)
        wid = us['verified'].walks('right')[0].id
        with self.assertRaises(LocalBudgetExceeded) as cm:
            us['verified'].carriers('right', wid, budget=LocalBudget(work_units=1))
        token = cm.exception.partial.resume
        for v in (None, 'count_only'):
            with self.assertRaises(ValueError):
                us[v].carriers('right', wid, resume=token)
        again = union_pairs(parts + absent, seqs_order={'g': 'verified'})
        self.assertEqual([str(s) for s in cm.exception.partial.rows] +
                         [str(s) for s in again.carriers('right', wid, resume=token)],
                         [str(s) for s in again.carriers('right', wid)])

    def test_seed_check_ignores_a_matching_seed_id(self):
        seq = T.graphlet('merge').seed.sequence
        r = absent_result('not_in_graph')
        r['seed']['seed_id'] = seq          # a caller's text proves nothing
        u = union_pairs([{'pair': 'g2', 'result': r}], seed=seq)
        self.assertEqual(u.absent[0].seed_check, 'unchecked')


def seed_id_variant(g, seed_id):
    """A shallow copy of |g| asked with another seed_id (one retrieval: seed_id is not part
    of the retrieval identity)."""
    h = copy.copy(g)
    h.seed_summary = copy.deepcopy(g.seed_summary)
    h.seed_summary['seed']['seed_id'] = seed_id
    return h


class TestOrder(unittest.TestCase):
    """P3: any order of the results gives the same digest and byte-equal answers, also with
    duplicate forms of one retrieval (envelope-stated and caller-stated, seed_id variants),
    absent records that differ only in seed_id or in the seed check, two failures of one
    pair and parts whose alphabets are spelled differently."""

    def answers(self, u):
        return (u.digest, u.alphabet, u.seed.sequence, str(u), u.to_fasta(stubs=True),
                u.to_gfa(), u.to_gfa(labels='full'), json.dumps(u.summary(), sort_keys=True),
                repr(u.walks('right', stubs=True)), str(u.carriers('right', 1)),
                json.dumps(u.provenance(), sort_keys=True, default=repr),
                json.dumps(u.evidence('right'), sort_keys=True, default=repr),
                repr(u.coverage('right')))

    def test_permutations(self):
        seq = T.graphlet('merge').seed.sequence
        e = with_env(T.graphlet('merge'), graph='gE')
        x = absent_result('not_in_graph')
        y = copy.deepcopy(x)
        y['seed']['seed_id'] = 'other'
        items = stated(gA='annotate') + [
            # one retrieval in four forms: envelope-stated and caller-stated, two seed_ids
            e, seed_id_variant(e, 'zz'), {'pair': 'gE', 'graphlet': e},
            {'pair': 'gE', 'graphlet': seed_id_variant(e, 'aa')},
            # an alphabet spelled without the sentinel (the same key)
            {'pair': 'gT', 'graphlet': with_fields(T.graphlet('merge_ties'),
                                                   alphabet='ACGT')},
            {'pair': 'g2', 'result': x}, {'pair': 'g2', 'result': y},
            {'pair': 'g3', 'result': absent_result('labels_not_in_graph')},
            {'pair': 'g3', 'result': absent_result('labels_not_in_graph'), 'sequence': seq},
            {'pair': 'g5', 'error': 'HTTP 503'}, {'pair': 'g5', 'error': 'timeout'}]
        expected = ['gA', 'gE', 'gT', 'g2', 'g3', 'g5', 'g9']
        ref = union_pairs(items, expected_pairs=expected)
        self.assertEqual((len(ref.parts), len(ref.absent)), (3, 3))
        want = self.answers(ref)
        rng = random.Random(7)
        for _ in range(12):
            perm = items[:]
            rng.shuffle(perm)
            exp = expected[:]
            rng.shuffle(exp)
            self.assertEqual(self.answers(union_pairs(perm, expected_pairs=exp)), want)
        # what is kept: the envelope-stated form with the first seed summary, the record
        # sorting first, the stronger seed check, the alphabet of the canonical first part
        prov = {p['tag']: p for p in ref.provenance() if p['kind'] == 'part'}
        self.assertEqual(prov['gE']['pair_stated_by'], 'envelope')
        self.assertEqual(prov['gE']['seed_id'], 'merge')
        absent = {a.tag: a for a in ref.absent}
        self.assertEqual(absent['g2'].detail['seed']['seed_id'], 'other')
        self.assertEqual(absent['g3'].seed_check, 'stated')
        self.assertEqual(ref.alphabet, ref.parts[0].graphlet.alphabet)
        self.assertEqual(seq, ref.seed.sequence)


class TestWalksByLabels(unittest.TestCase):
    """walks(labels=): the walks those labels carry whole in some pair -- a span that holds
    the walk, its own labels or its route's, never a stub's labels before its merge nor a
    span that ends or turns off (review B1)."""

    def check(self, u, side, pair, label, excluded):
        got = u.walks(side, by='id', stubs=True, labels=[{'pair': pair, 'ref': label.ref}])
        self.assertNotIn(excluded, [w.id for w in got])
        for w in got:
            self.assertIn(pair, w.held_in)
            (s,) = [s for s in u.carriers(side, w.id) if s.pair == pair]
            carried = s.labels if s.how in ('whole', 'extends') else s.route.labels
            self.assertIn((label, 0), [(x.label, x.from_bp) for x in carried], w.id)
        return got

    def test_a_stub_label_before_its_merge_carries_nothing(self):
        u = union_pairs(stated(gM='merge', gT='merge_ties'))
        seen = False
        for w in u.walks('right', by='id', stubs=True):
            for s in u.carriers('right', w.id):
                if s.pair == 'gM' and s.how == 'merges' and s.route.how == 'held':
                    stub0 = {x.label.name for x in s.labels if x.from_bp == 0}
                    route0 = {x.label.name for x in s.route.labels if x.from_bp == 0}
                    for lab in stub0 - route0:          # both.fa: from 62 bp on the route
                        seen = True
                        self.check(u, 'right', 'gM', u.label({'pair': 'gM', 'name': lab}), w.id)
                    for lab in route0:                  # a.fa carries it whole
                        got = u.walks('right', by='id', stubs=True,
                                      labels=[{'pair': 'gM', 'name': lab}])
                        self.assertIn(w.id, [x.id for x in got])
        self.assertTrue(seen)

    def test_r237_a_pair_that_ends_carries_nothing(self):
        parts, u = r237_union()
        pair, seen = 'bact_pseudomonadota', 0
        for w in u.walks('right', by='id', stubs=True):
            (s,) = [s for s in u.carriers('right', w.id) if s.pair == pair]
            if s.how in ('ends', 'diverges'):
                whole = [x.label for x in s.labels if x.from_bp == 0]
                if whole:
                    seen += 1
                    self.check(u, 'right', pair, whole[0], w.id)
        self.assertGreaterEqual(seen, 3)                # walks 4, 6 and 10 at 56d9a56fc


class TestLabels(unittest.TestCase):
    """P7: labels per pair, shared names, the header-name warnings (OQ20)."""

    def test_labels_with_their_pair(self):
        u = locus()
        labs = u.labels()
        self.assertTrue(all(str(l) == '%s/%s' % (l.pair, l.name) for l in labs))
        self.assertEqual(u.labels('gM'), [l for l in labs if l.pair == 'gM'])
        self.assertEqual(u.label({'pair': 'gM', 'name': 'a.fa'}).pair, 'gM')
        ref = u.labels('gT')[0].ref
        self.assertEqual(u.label({'pair': 'gT', 'ref': ref}).ref, ref)
        self.assertEqual(u.shared_names()['a.fa'], ('gA', 'gM', 'gT'))
        self.assertIn('names occur in more than one pair', u.statements('right').labels)
        g = union_pairs(stated(gM='merge', gF='fork' if False else 'merge_ties'))
        self.assertTrue(g.shared_names())

    def test_same_name_in_one_pair(self):
        from metagraph.traverse import AmbiguousLabel
        u = union_pairs(stated(gS='same_name'))
        with self.assertRaises(AmbiguousLabel):
            u.label({'pair': 'gS', 'name': 'ACC1'})
        with self.assertRaises(AmbiguousLabel):
            u.label({'name': 'ACC1'})
        refs = [l.ref for l in u.labels('gS') if l.name == 'ACC1']
        self.assertEqual(u.label({'pair': 'gS', 'ref': refs[1]}).ref, refs[1])

    def test_header_name_warnings(self):
        # same_name has header labels; merge has column labels only
        parts = stated(gS='same_name', gM='merge') if False else None  # noqa: F841
        hs = union_pairs(stated(gS='same_name'))
        self.assertEqual(hs._mapping_warnings(), [('gS', 'unknown')])
        self.assertIn('header names of every pair rest on an unchecked record mapping (unknown)',
                      hs.statements('right').labels)
        self.assertIn('mapping=unknown', hs.to_fasta())
        self.assertIn('unchecked record mapping', hs.to_gfa())
        self.assertIn('unchecked record mapping', str(hs.carriers('right', 0)))
        for v, warned in (('count_only', True), ('verified', False), (None, False)):
            u = union_pairs(stated(gS='same_name'), seqs_order={'gS': v})
            self.assertEqual(bool(u._mapping_warnings()), warned, v)
            self.assertEqual('mapping=' in u.to_fasta(), warned)
        # column labels carry no warning
        self.assertEqual(locus()._mapping_warnings(), [])


class TestCoverage(unittest.TestCase):
    """P8: every category of §2.6, the sentences, the default depth."""

    def test_every_category(self):
        seq = T.graphlet('merge').seed.sequence
        m = T.graphlet('merge')
        lim = copy.copy(m.limitations)
        from metagraph.traverse.model import Limitation
        seedlim = Limitation(None, 'seed_labels', 'labels.max_seed_labels', 1000, 4312, None,
                             [], 'carriers of the seed were cut')
        q = with_fields(T.graphlet('merge_ties'), limitations=list(lim) + [seedlim])
        items = [{'pair': 'g1', 'graphlet': m}, {'pair': 'g6', 'graphlet': q},
                 {'pair': 'g7', 'graphlet': T.graphlet('annotate')},
                 {'pair': 'g2', 'result': absent_result('not_in_graph')},
                 {'pair': 'g3', 'result': absent_result('not_in_graph', present=4)},
                 {'pair': 'g4', 'result': absent_result('labels_not_in_graph')},
                 {'pair': 'g5', 'result': absent_result('failed'), 'sequence': seq},
                 {'pair': 'g10', 'error': 'deadline exceeded'}]
        u = union_pairs(items, expected_pairs=['g%d' % i for i in range(1, 12)])
        c = u.coverage('right', at_bp=80)
        self.assertEqual(c.by_category, {
            'complete_qualified': ['g1', 'g6'], 'incomplete': ['g7'],
            'seed_not_in_graph': ['g2', 'g3'], 'labels_not_in_graph': ['g4'],
            'not_walked': ['g5'], 'failed': ['g10'], 'not_included': ['g11', 'g8', 'g9']})
        self.assertIn('complete for 1000 of 4312 seed labels', c.sentence)
        self.assertIn('incomplete in g7 [annotate, kmer, k 15] (certified to 70 bp: radius)',
                      c.sentence)
        self.assertIn('seed not in graph: g2 (0 of 16 k-mers), g3 (4 of 16 k-mers)',
                      c.sentence)
        self.assertIn('not walked: g5 (seed in graph; derivation: no label carries the whole '
                      'seed)', c.sentence)
        self.assertIn('failed: g10 (deadline exceeded)', c.sentence)
        self.assertIn('Seed of the absent records: stated by the caller for g5, not checked '
                      'for g2, g3, g4.', c.sentence)
        self.assertTrue(c.sentence.endswith(U.TURN_OFFS + ' ' + U.SCOPE))
        self.assertIsNone(c.every_pair_certified_to_bp)
        left = u.coverage('left', at_bp=50)
        self.assertEqual(left.by_category['arm_not_walked'], ['g1', 'g6'])

    def test_evidence_states_the_five_statements(self):
        u = locus()
        ev = u.evidence()
        self.assertEqual(u.arms(), ['left', 'right'])
        self.assertEqual(ev['statements'], [
            _core(u.coverage('left')), _core(u.coverage('right')),
            'walks: ' + ', '.join(u.statements(s).provenance[len('walks: '):]
                                  for s in u.arms()),
            U.TURN_OFFS, u._labels_sentence(), U.SCOPE])
        self.assertEqual(u.evidence('right')['statements'],
                         u.statements('right').sentences())
        # an arm-less union: its own five
        v = union_pairs([{'pair': 'g2', 'result': absent_result('not_in_graph')}],
                        seed=T.graphlet('merge').seed.sequence)
        self.assertEqual(v.evidence()['statements'], v.statements().sentences())
        self.assertEqual(len(v.evidence()['statements']), 5)
    def test_depths_and_every_pair_certified(self):
        m = T.graphlet('merge')               # complete_to_bp 100
        u = union_pairs([{'pair': 'g1', 'graphlet': m}], expected_pairs=['g1'])
        self.assertEqual(u.coverage('right').at_bp, 100)          # the largest radius
        # merge's scope limitation makes it complete_qualified: every_pair None
        self.assertIsNone(u.coverage('right').every_pair_certified_to_bp)
        plain = with_fields(m, limitations=[l for l in m.limitations if l.kind != 'scope'])
        u = union_pairs([{'pair': 'g1', 'graphlet': plain},
                         {'pair': 'g2', 'result': absent_result('not_in_graph')}],
                        expected_pairs=['g1', 'g2'])
        self.assertEqual(u.coverage('right', at_bp=60).by_category['complete'], ['g1'])
        self.assertEqual(u.coverage('right', at_bp=60).every_pair_certified_to_bp, 100)
        self.assertEqual(u.coverage('right', at_bp=150).by_category['incomplete'], ['g1'])
        u2 = union_pairs([{'pair': 'g1', 'graphlet': plain}])
        self.assertIsNone(u2.coverage('right', at_bp=60).every_pair_certified_to_bp)
        self.assertIn(U.OVER_GIVEN, u2.coverage('right').sentence)

    def test_mixed_strategies_are_stated(self):
        # P11 on the committed fixtures: constrain and annotate parts
        u = locus()
        s = u.coverage('right').sentence
        self.assertIn('gA [annotate, kmer, k 15]', s)
        self.assertIn('gM [constrain, kmer, k 15]', s)
        for w in u.walks('right', by='id', stubs=True):
            for sp in u.carriers('right', w.id):
                self.assertEqual(sp.mode, u._part_of[sp.pair].graphlet.mode)
                self.assertEqual(sp.k, 15)


def _first_label(u):
    """{pair, ref} of the first part's first label (read off the model, unbudgeted and
    nothing allocated beyond the ref)."""
    p = u.parts[0]
    return {'pair': p.tag, 'ref': p.graphlet.labels[0].ref}


class TestBudgets(unittest.TestCase):
    """P12: every operation under budget= completes or stops at the same row every run."""

    def _ops(self, u):
        return {
            'walks': lambda b, r=None: u.walks('right', budget=b, resume=r),
            'segments': lambda b, r=None: u.segments('right', budget=b, resume=r),
            'carriers': lambda b, r=None: u.carriers('right', 1, budget=b, resume=r),
            'fasta': lambda b, r=None: u.to_fasta(budget=b),
            'gfa': lambda b, r=None: u.to_gfa(budget=b),
            'summary': lambda b, r=None: u.summary(budget=b),
            'spell': lambda b, r=None: u.spell('right', 1, budget=b),
        }

    def test_stops_are_deterministic_and_resume(self):
        for limit in (50, 400, 2000, 10 ** 7):
            outcome = []
            for _ in range(2):
                u = locus()
                got = {}
                for name, fn in self._ops(u).items():
                    try:
                        fn(LocalBudget(work_units=limit))
                        got[name] = 'done'
                    except LocalBudgetExceeded as e:
                        p = e.partial
                        got[name] = ('stop', None if p is None else len(p.rows))
                        if p is not None and p.resume is not None and name != 'carriers':
                            rest = fn(None, p.resume)
                            whole = fn(None)
                            self.assertEqual(repr(p.rows + rest), repr(whole), name)
                outcome.append(got)
            self.assertEqual(outcome[0], outcome[1], limit)
        u = locus()
        for name, fn in self._ops(u).items():
            fn(LocalBudget(work_units=10 ** 9))

    # every public PairUnion method that takes budget= (GPT review M3), as (name, call);
    # arms(), provenance(), memory_bytes() and the stored attributes are exempt
    BUDGETED = (
        ('walks', lambda u, **k: u.walks('right', stubs=True, **k)),
        ('walks(top=3)', lambda u, **k: u.walks('right', top=3, **k)),
        ('walks(labels=)', lambda u, **k: u.walks('right', labels=[_first_label(u)], **k)),
        ('segments', lambda u, **k: u.segments('right', **k)),
        ('carriers', lambda u, **k: u.carriers('right', 1, **k)),
        ('spell', lambda u, **k: u.spell('right', 1, **k)),
        ('statements', lambda u, **k: u.statements('right', **k)),
        ('coverage', lambda u, **k: u.coverage('right', **k)),
        ('labels', lambda u, **k: u.labels(**k)),
        ('label', lambda u, **k: u.label(_first_label(u), **k)),
        ('label(name)', lambda u, **k: u.label({'pair': u.parts[0].tag, 'name': 'a.fa'}, **k)),
        ('shared_names', lambda u, **k: u.shared_names(**k)),
        ('evidence', lambda u, **k: u.evidence(**k)),
        ('evidence(arm)', lambda u, **k: u.evidence('right', **k)),
        ('summary', lambda u, **k: u.summary(**k)),
        ('to_fasta', lambda u, **k: u.to_fasta(stubs=True, **k)),
        ('to_gfa', lambda u, **k: u.to_gfa(labels='full', **k)))

    def test_cold_and_warm_stop_alike(self):
        """GPT review M3: each budgeted method charges the same cold (a fresh union) and
        warm (after an unbudgeted call filled its caches), so a budget one unit or one byte
        short of a cold call's charge stops the warm call too, at the same row, and the
        cold call's charge lets both finish. (At d958961a8 the five methods without
        budget= fail in their own subtest on the signature; the rest pass there.)"""
        def run(u, fn, **lim):
            try:
                fn(u, budget=LocalBudget(**lim))
                return 'done'
            except LocalBudgetExceeded as e:
                p, s = e.partial, e.stop
                return ('stop', s.resource, s.op, s.phase, s.used, s.needed_at_least,
                        None if p is None else repr(p.rows))

        for name, fn in self.BUDGETED:
            with self.subTest(name):
                b = LocalBudget()
                fn(locus(), budget=b)
                work, peak = b.used_work, b.peak_bytes
                warm = locus()
                fn(warm)
                again = LocalBudget()
                fn(warm, budget=again)
                self.assertEqual((again.used_work, again.peak_bytes), (work, peak))
                for lim in ({'work_units': work - 1}, {'memory_mb': (peak - 1) / (1 << 20)}):
                    cold = run(locus(), fn, **lim)
                    self.assertNotEqual(cold, 'done', lim)
                    self.assertEqual(run(warm, fn, **lim), cold, lim)
                for lim in ({'work_units': work}, {'memory_mb': peak / (1 << 20)}):
                    self.assertEqual((run(locus(), fn, **lim), run(warm, fn, **lim)),
                                     ('done', 'done'), lim)

    def test_the_ambient_budget_alone_and_around_another(self):
        """Opus fix review N1: under an ambient budget alone (local_budget(), no budget=)
        every budgeted method charges the same cold and warm, as with budget=; and a call
        given budget= inside local_budget() charges only its own budget, the ambient one
        nothing (a part's segment labels were built under the ambient budget, outside the
        price the union charges for them). Locus, and a union with a stitched part."""
        from metagraph.traverse import local_budget
        from tests.test_traverse_union_stitched import pairs, role, stitched
        cases = (('locus', locus, self.BUDGETED + (('str', None),)),
                 ('comb_chain', lambda: union_pairs(pairs(stitched('comb_chain'),
                                                          role('comb_chain', 'reference'))),
                  [x for x in self.BUDGETED if x[0] != 'label(name)']))
        for case, make, ops_ in cases:
            for name, fn in ops_:
                call = fn or (lambda u: str(u))
                with self.subTest(case=case, op=name):
                    got = []
                    for warm in (False, True):
                        u = make()
                        if warm:
                            call(u)
                        with local_budget() as a:
                            call(u)
                        got.append((a.used_work, a.peak_bytes))
                    self.assertEqual(got[0], got[1])
                    self.assertGreater(got[0][0], 0)
                    if fn is None:
                        continue
                    alone = LocalBudget()
                    fn(make(), budget=alone)
                    inner, u = LocalBudget(), make()
                    with local_budget() as a:
                        fn(u, budget=inner)
                    self.assertEqual((a.used_work, a.peak_bytes), (0, 0))
                    self.assertEqual((inner.used_work, inner.peak_bytes),
                                     (alone.used_work, alone.peak_bytes))
                    self.assertEqual((alone.used_work, alone.peak_bytes), got[0])

    def test_a_stitched_parts_evidence_is_the_union_calls(self):
        """Decision 2 (train 7): a stitched part's own evidence() runs under the union
        call's budget, cold as warm, the ambient budget untouched; a stop inside it is the
        union call's -- op union_evidence, the union's levers (Opus fix review N5) -- and its
        nested scope leaves neither its name nor its levers to a later stop."""
        from metagraph.traverse import StitchedGraphlet, local_budget
        from tests.test_traverse_union_stitched import pairs, role, stitched
        make = lambda: union_pairs(pairs(stitched('comb_chain'), role('comb_chain',
                                                                       'reference')))
        own, sg = LocalBudget(), stitched('comb_chain')
        with local_budget(own):
            sg.evidence('right')
        cold = LocalBudget()
        make().evidence('right', budget=cold)
        warm = make()
        warm.evidence('right')
        for p in warm.parts:
            if isinstance(p.graphlet, StitchedGraphlet):
                p.graphlet.evidence('right')
        again, inner, u = LocalBudget(), LocalBudget(), make()
        warm.evidence('right', budget=again)
        with local_budget() as a:
            u.evidence('right', budget=inner)
        self.assertEqual((again.used_work, again.peak_bytes), (cold.used_work, cold.peak_bytes))
        self.assertEqual((inner.used_work, inner.peak_bytes), (cold.used_work, cold.peak_bytes))
        self.assertEqual((a.used_work, a.peak_bytes), (0, 0))
        self.assertGreaterEqual(cold.used_work, own.used_work)
        # every stop of a sweep of work limits names the union's call
        for lim in range(1, cold.used_work, max(1, cold.used_work // 25)):
            try:
                make().evidence('right', budget=LocalBudget(work_units=lim))
            except LocalBudgetExceeded as e:
                self.assertEqual((e.stop.op, e.stop.actions),
                                 ('union_evidence', ['raise_local_budget', 'process_locally']),
                                 lim)
                self.assertIsNone(e.partial)
        # a part whose evidence() opens its own scope with levers the union has not: a stop
        # inside it, and one of the union's after it, both state the union's call
        from metagraph.traverse import budget as B

        def scoped(self, arm=None, *, budget=None):
            def body(b):
                b.charge(4321, 1000)
                return real(self, arm)
            return B.run_scoped(budget, 'evidence', ('select_walks', 'narrow_arm'), body)
        real = StitchedGraphlet.evidence
        StitchedGraphlet.evidence = scoped
        try:
            full = LocalBudget()
            make().evidence('right', budget=full)
            stops = []
            for lim in range(1, full.used_work, max(1, full.used_work // 40)):
                try:
                    make().evidence('right', budget=LocalBudget(work_units=lim))
                except LocalBudgetExceeded as e:
                    s = e.stop
                    stops.append((lim, s.needed_at_least - s.used == 4321))
                    self.assertEqual((s.op, s.actions), ('union_evidence', [
                        'raise_local_budget', 'process_locally']), lim)
                    self.assertNotIn('select_walks', str(e))
        finally:
            StitchedGraphlet.evidence = real
        # stops inside the part's evidence, and stops of the union's own after it
        inside = [lim for lim, x in stops if x]
        self.assertTrue(inside)
        self.assertTrue(any(lim > max(inside) for lim, x in stops if not x))

    def test_the_ambient_budget(self):
        """A budget set by local_budget() is every budgeted method's (print(u) included,
        charged as statements()); the exempt accessors return under it."""
        from metagraph.traverse import local_budget
        u = locus()
        for name, fn in self.BUDGETED + (('str', lambda u: str(u)),):
            with local_budget(work_units=0) as b:
                with self.assertRaises(LocalBudgetExceeded, msg=name):
                    fn(u)
            with local_budget(work_units=10 ** 9) as b:
                fn(u)
            self.assertGreater(b.used_work, 0, name)
        with local_budget(work_units=0) as b:
            u.arms(), u.provenance(), u.memory_bytes(), u.digest, u.parts
        self.assertEqual(b.used_work, 0)

    def test_memory_bytes_bounds_the_heap_on_11_parts(self):
        """P12: memory_bytes() bounds the traced heap the union holds on 11 parts, after
        construction and after every query (their answers dropped): the figures and label
        intervals the queries keep are counted. Warm (the pattern of test_traverse_memory):
        free lists filled by a first run, the collector held off while tracing. The models
        are loaded before tracing; what the union adds to their caches is in both sides.
        memory_bytes() is a getsizeof deep size and tracemalloc counts allocations: they
        differ by a few percent, so the bound after the queries has 10% slack."""
        import gc
        import tracemalloc
        from metagraph.traverse import composite as C
        names = [n for n in UNION_DOCS if not n.startswith('repeat')][:11]
        self.assertEqual(len(names), 11)

        def queries(u):
            out = []
            for side in u.arms():
                for w in u.walks(side, stubs=True):
                    out.append(u.carriers(side, w.id))
                out += [u.segments(side), u.statements(side), u.coverage(side)]
            out += [str(u), u.summary(), u.to_fasta(stubs=True), u.to_gfa(labels='full'),
                    u.evidence(), u.provenance()]
            return len(out)

        warm = union_pairs([{'pair': n, 'graphlet': union_doc(n)} for n in names])
        queries(warm)
        warm.memory_bytes()
        del warm
        gs = [union_doc(n) for n in names]
        before = C.models_memory_bytes(gs)
        gc.collect()
        gc.disable()
        tracemalloc.start()
        try:
            base = tracemalloc.get_traced_memory()[0]
            u = union_pairs([{'pair': n, 'graphlet': g} for n, g in zip(names, gs)])
            gc.collect()
            built = tracemalloc.get_traced_memory()[0] - base
            counted_built = u.memory_bytes() - before
            queries(u)
            gc.collect()
            held = tracemalloc.get_traced_memory()[0] - base
            counted = u.memory_bytes() - before
        finally:
            tracemalloc.stop()
            gc.enable()
        self.assertLessEqual(built, 1.1 * counted_built, (built, counted_built))
        self.assertLessEqual(held, 1.1 * counted, (held, counted))
        # the queries keep their figures and label intervals, and memory_bytes() says so
        self.assertGreater(counted - counted_built, 0.8 * (held - built))

    def test_the_account_bounds_the_peak(self):
        """budget.py: "the account bounds the call's peak", for union_pairs() and every union
        query: the budget's peak_bytes against tracemalloc's peak of the same call on fresh
        models and a fresh union (every derivation and cache cold), on the split fixtures'
        11 parts, r237's export and a comb. The interpreter's free lists are filled first by
        the same call on other fresh copies, the collector held off (the pattern of
        test_traverse_memory: a collection empties the free lists, and freed tuples parked
        there stay traced; warm, the peak is the library's own). Each charge is the cold
        price: a second call on the warm union charges the same. GPT review M4: also on
        long metadata, charged before it is serialized -- a 100,000-character ASCII tag, one
        of escaped and astral characters, a 1 MB error and an absent record whose error is
        100,000 quotes (escaped twice in the digest)."""
        import gc
        import tracemalloc

        def split11():
            return [{'pair': n, 'graphlet': union_doc(n)}
                    for n in [n for n in UNION_DOCS if not n.startswith('repeat')][:11]]

        def comb():
            return [{'pair': 'C', 'graphlet': comb_part(200)},
                    {'pair': 'D', 'graphlet': comb_part(3)}]

        def long_tags():
            return [{'pair': 'X' * 100000, 'graphlet': T.graphlet('merge')},
                    {'pair': ',;= \u00e9\U0001F600' * 5000, 'graphlet': T.graphlet('merge_ties')}]

        def long_errors():
            r = dict(absent_result('not_in_graph'), error='"' * 100000)
            return [{'pair': 'gM', 'graphlet': T.graphlet('merge')},
                    {'pair': 'g', 'error': 'X' * 10 ** 6}, {'pair': 'h', 'result': r}]

        def escaped_tag():
            # Opus fix review N2: an absent record's tag of 100,000 token characters, each
            # escaped to three bytes in every header and GFA line naming it
            return [{'pair': 'gM', 'graphlet': T.graphlet('merge')},
                    {'pair': '/:' * 50000, 'error': 'boom'}]

        def long_names():
            # Opus fix review N4/N11: label names of 100,000 characters and more, equal in
            # the two parts where the label ids are (shared names), read by shared_names(),
            # the statements, evidence, summary and the exports, and as a walks(labels=)
            # selector
            out = []
            for tag, name in (('gM', 'merge'), ('gT', 'merge_ties')):
                g = copy.deepcopy(T.graphlet(name))
                g.labels = [dataclasses.replace(l, name=('N%d' % i) * 50000)
                            for i, l in enumerate(g.labels)]
                out.append({'pair': tag, 'graphlet': g})
            return out

        def peak(fn, warm):
            gc.collect()
            gc.disable()
            try:
                warm()
                tracemalloc.start()
                base = tracemalloc.get_traced_memory()[0]
                fn()
                return tracemalloc.get_traced_memory()[1] - base
            finally:
                tracemalloc.stop()
                gc.enable()

        ops_ = {'walks': lambda u, **k: u.walks('right', stubs=True, **k),
                'walks(top=3)': lambda u, **k: u.walks('right', top=3, **k),
                'walks(spell)': lambda u, **k: u.walks('right', stubs=True, spell=True, **k),
                'segments': lambda u, **k: u.segments('right', **k),
                'carriers': lambda u, **k: u.carriers('right', 1, **k),
                'statements': lambda u, **k: u.statements('right', **k),
                'summary': lambda u, **k: u.summary(ceiling[0], **k),
                'to_fasta': lambda u, **k: u.to_fasta(stubs=True, **k),
                'to_gfa': lambda u, **k: u.to_gfa(labels='full', **k),
                'spell': lambda u, **k: u.spell('right', 1, **k),
                'coverage': lambda u, **k: u.coverage('right', **k),
                'labels': lambda u, **k: u.labels(**k),
                'label': lambda u, **k: u.label(_first_label(u), **k),
                'shared_names': lambda u, **k: u.shared_names(**k),
                'evidence': lambda u, **k: u.evidence(**k)}
        by_name = {'walks(labels=name)': lambda u, **k: u.walks('right', labels=[
            {'pair': u.parts[0].tag, 'name': u.parts[0].graphlet.labels[0].name}], **k)}
        ceiling = [4096]
        for name, make in (('split11', split11), ('r237', lambda: r237_union()[0]),
                           ('comb', comb), ('long tags', long_tags),
                           ('long errors', long_errors), ('escaped tag', escaped_tag),
                           ('long names', long_names)):
            # the long statements never fit 4096 bytes: summary() is measured whole there
            ceiling[0] = 4096 if name in ('split11', 'r237', 'comb') else 1 << 21
            b = LocalBudget()
            union_pairs(make(), budget=b)
            parts, other = make(), make()
            got = peak(lambda: union_pairs(parts), lambda: union_pairs(other))
            self.assertLessEqual(got, b.peak_bytes, (name, 'union_pairs'))
            for op, fn in dict(ops_, **(by_name if name == 'long names' else {})).items():
                b = LocalBudget()
                v = union_pairs(make())
                fn(v, budget=b)
                # the cold price: the same charge again on the warm union
                again = LocalBudget()
                fn(v, budget=again)
                self.assertEqual((again.used_work, again.peak_bytes),
                                 (b.used_work, b.peak_bytes), (name, op))
                u, w = union_pairs(make()), union_pairs(make())
                got = peak(lambda: fn(u), lambda: fn(w))
                self.assertLessEqual(got, b.peak_bytes, (name, op))

    def test_the_account_bounds_a_long_seed(self):
        """Opus fix review N2: seed= of 1 Mbp over absent records only: its case-mapped copy
        and the digest's JSON and bytes of it are charged before they are made; the queries
        of the arm-less union are bounded too (warm free lists, as above)."""
        import gc
        import tracemalloc

        def make():
            return [{'pair': 'g', 'error': 'x'},
                    {'pair': 'h', 'result': absent_result('not_in_graph', seed_len=10 ** 6)}]
        seed = 'acgt' * 250000

        def peak(fn):
            gc.collect()
            gc.disable()
            try:
                fn()
                tracemalloc.start()
                base = tracemalloc.get_traced_memory()[0]
                fn()
                return tracemalloc.get_traced_memory()[1] - base
            finally:
                tracemalloc.stop()
                gc.enable()
        b = LocalBudget()
        union_pairs(make(), seed=seed, budget=b)
        parts = make()
        self.assertLessEqual(peak(lambda: union_pairs(parts, seed=seed)), b.peak_bytes)
        u = union_pairs(make(), seed=seed)
        for name, fn in (('summary', lambda **k: u.summary(1 << 21, **k)),
                         ('to_gfa', lambda **k: u.to_gfa(**k)),
                         ('evidence', lambda **k: u.evidence(**k))):
            b = LocalBudget()
            fn(budget=b)
            self.assertLessEqual(peak(fn), b.peak_bytes, name)

    def test_the_preparation_is_admitted_first(self):
        """GPT final review M4: union_pairs() admits its preparation -- classifying the
        inputs, their descriptions, the matching passes quadratic in their number -- before
        it runs. Allocation is measured before the refusal, not only the refusal: under the
        reviewer's zero work and 1 KB budget, 1,000 failed records, 1,000 expected pairs,
        1,000 seqs_order keys and 120 failed records of one pair (whose ambiguous_pair check
        is cubic) stop with nothing charged and a traced peak below 64 KB (it was 1.26 MB and
        2.33 s for the first); unbudgeted and budgeted answers agree, and the account bounds
        the traced peak of the whole call (warm free lists, as above)."""
        import gc
        import tracemalloc
        forms = {
            'failed': lambda n: ([{'pair': 'g%d' % i, 'error': 'request failed'}
                                  for i in range(n)], {}),
            'one pair': lambda n: ([{'pair': 'g', 'error': 'request failed %d' % i}
                                    for i in range(n)], {}),
            'expected': lambda n: ([{'pair': 'g0', 'error': 'x'}],
                                   {'expected_pairs': ['g%d' % i for i in range(n)]}),
            'seqs_order': lambda n: ([{'pair': 'g0', 'error': 'x'}],
                                     {'seqs_order': {'g%d' % i: 'verified' for i in range(n)}}),
            'parts': lambda n: ([{'pair': 'g%d' % i, 'graphlet': T.graphlet('merge')}
                                 for i in range(n)], {})}

        def peak(fn):
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
        for name, make in forms.items():
            with self.subTest(name):
                n = {'one pair': 120, 'parts': 100}.get(name, 1000)
                args = [make(n), make(n)]
                zero = []

                def refused():
                    r, kw = args.pop()
                    b = LocalBudget(work_units=0, memory_mb=0.001)
                    zero.append(b)
                    union_pairs(r, seed='ACGT' if name != 'parts' else None, budget=b, **kw)
                self.assertLess(peak(refused), 64 * 1024)
                self.assertEqual(zero[-1].call_usage(), {'work_units': 0, 'memory_bytes': 0})
                self.assertEqual(zero[-1].stops[-1].resource, 'work')
                # the whole call: the same union, its account above its traced peak
                n = {'one pair': 40, 'parts': 20}.get(name, 200)
                r, kw = make(n)
                seed = 'ACGT' if name != 'parts' else None
                b = LocalBudget()
                self.assertEqual(union_pairs(r, seed=seed, budget=b, **kw).digest,
                                 union_pairs(r, seed=seed, **kw).digest)
                args = [make(n), make(n)]

                def whole():
                    r, kw = args.pop()
                    union_pairs(r, seed=seed, **kw)
                self.assertLessEqual(peak(whole), b.peak_bytes)

    def test_a_spelled_walk_charges_its_seed(self):
        """GPT final review M2: walks(spell=True) returns the seed with each flank (a
        100,000-base seed, 40-base flanks), and charges the seed, the flank and the
        temporary copies before it spells them: the account bounds the traced peak on both
        arms (warm free lists, as above), the same walks come back with and without a
        budget, and under the reviewer's 26,214-byte budget the call stops before it
        allocates the seed's copy (the traced peak before the stop stays below it)."""
        import gc
        import tracemalloc
        g = long_seed(T.graphlet('fork'), 100000)
        n = g.seed.length_bp

        def peak(fn):
            gc.collect()
            gc.disable()
            try:
                try:
                    fn()
                except LocalBudgetExceeded:
                    pass
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
        for side in ('left', 'right'):
            u = union_pairs([{'pair': 'F', 'graphlet': g}])
            plain = u.walks(side, spell=True, stubs=True)
            self.assertTrue(plain)
            for w in plain:
                self.assertEqual(len(w.sequence), n + w.length_bp)
                self.assertEqual(w.sequence, u.spell(side, w.id, with_seed=True))
            b = LocalBudget()
            got = u.walks(side, spell=True, stubs=True, budget=b)
            self.assertEqual([(w.id, w.sequence) for w in got],
                             [(w.id, w.sequence) for w in plain])
            self.assertLessEqual(peak(lambda: u.walks(side, spell=True, stubs=True)),
                                 b.peak_bytes, side)
            self.assertGreaterEqual(b.peak_bytes, 3 * n)
            tight = LocalBudget(memory_mb=0.025)
            with self.assertRaises(LocalBudgetExceeded) as cm:
                u.walks(side, spell=True, stubs=True, budget=tight)
            self.assertEqual(cm.exception.stop.resource, 'memory')
            self.assertGreater(cm.exception.stop.needed_at_least, 3 * n)
            self.assertLess(peak(lambda: u.walks(side, spell=True, stubs=True,
                                                 budget=LocalBudget(memory_mb=0.025))), n,
                            side)

    def test_the_account_bounds_a_gzip_load(self):
        """GPT review M4: load() of a gzip file charges the decompressor's working buffers
        before it opens (fixed: about 155 KB on CPython 3.11, 500 KB on 3.14), so the account
        bounds the traced peak for small files as for the r237 parts; warm free lists, as
        above."""
        import gc
        import gzip
        import tempfile
        import tracemalloc
        from metagraph.traverse import load
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        paths = []
        for name in ('merge', 'annotate'):
            paths.append(os.path.join(tmp.name, name + '.mgt.gz'))
            with open(paths[-1], 'wb') as f:
                f.write(gzip.compress(T.doc_text(name).encode('utf-8')))
        paths += [os.path.join(R237, p + '.mgt.gz') for p in sorted(R237_BODIES)]
        for path in paths:
            b = LocalBudget()
            load(path, budget=b)
            gc.collect()
            gc.disable()
            try:
                load(path, budget=LocalBudget())
                tracemalloc.start()
                base = tracemalloc.get_traced_memory()[0]
                load(path, budget=LocalBudget())
                got = tracemalloc.get_traced_memory()[1] - base
            finally:
                tracemalloc.stop()
                gc.enable()
            self.assertLessEqual(got, b.peak_bytes, path)

    def test_construction_charged(self):
        b = LocalBudget(work_units=10 ** 9)
        union_pairs(stated(gA='annotate', gM='merge'), budget=b)
        self.assertGreater(b.used_work, 0)
        with self.assertRaises(LocalBudgetExceeded):
            union_pairs(stated(gA='annotate', gM='merge'), budget=LocalBudget(work_units=10))

    def test_resume_token_bound_to_the_union(self):
        u, v = locus(), union_pairs(stated(gA='annotate'))
        try:
            u.walks('right', budget=LocalBudget(work_units=400))
        except LocalBudgetExceeded as e:
            if e.partial is not None and e.partial.resume is not None:
                with self.assertRaises(ValueError):
                    v.walks('right', resume=e.partial.resume)


SPLIT_MODES = [(m, r) for m in ('annotate', 'constrain') for r in ('merge', 'keep')]
UNION_DOCS = sorted(f[:-len('.graphlet.json')] for f in os.listdir(UNION)
                    if f.endswith('.graphlet.json'))


def split_union(mode, rec):
    """Graph A's and graph B's parts of the split index, the unsplit index's retrieval, and
    the union of A and B (both expected)."""
    a, b, al = (union_doc('%s_%s_%s' % (x, mode, rec)) for x in ('a', 'b', 'all'))
    return a, b, al, union_pairs([{'pair': 'A', 'graphlet': a}, {'pair': 'B', 'graphlet': b}],
                                 expected_pairs=['A', 'B'])


def dag_spellings(g, side):
    """Every root-to-leaf spelling of an arm's segment DAG (through every parent of a merge),
    in walking order: what the retrieval walked, displayed or not."""
    segs = g.arms[side].segments
    kids = {s.id: [] for s in segs}
    for s in segs:
        for p in s.parents:
            kids[p].append(s.id)
    out = []

    def go(i, acc):
        acc += segs[i].walk
        if not kids[i]:
            out.append(acc)
        for c in kids[i]:
            go(c, acc)

    for s in segs:
        if not s.parents:
            go(s.id, '')
    return out


def split_record(header):
    return T.cli_record('split_all', header)


class TestSplitGraphs(unittest.TestCase):
    """P6, the split-graph oracle and the non-claim: one tiny index split by record into two
    graphs (A: a1 and the plasmid pl; B: b1, b2 and pl) and the unsplit index, keep and merge,
    annotate and constrain."""

    def test_every_union_spelling_is_a_prefix_of_the_unsplit_retrieval(self):
        for mode, rec in SPLIT_MODES:
            a, b, al, u = split_union(mode, rec)
            for side in u.arms():
                dag = dag_spellings(al, side)
                for w in u.walks(side, by='id', stubs=True):
                    sp = u.spell(side, w.id, 'walk')
                    self.assertTrue(any(d.startswith(sp) for d in dag), (mode, rec, side, w.id))

    def test_the_shared_record_is_carried_by_both_pairs(self):
        flank = split_record('pl')[30:]             # pl after the 30 bp seed
        for mode, rec in SPLIT_MODES:
            a, b, al, u = split_union(mode, rec)
            (w,) = [w for w in u.walks('right', by='id') if u.spell('right', w.id) == flank]
            self.assertEqual(w.held_in, ('A', 'B'), (mode, rec))
            if mode == 'constrain' and rec == 'keep':
                # pl's own label carries it whole in both pairs
                for s in u.carriers('right', w.id):
                    self.assertIn(('pl', 0), [(x.label.name, x.from_bp) for x in s.labels])

    def test_a_walk_no_graph_holds_whole_is_not_in_the_union(self):
        a1, b1 = split_record('a1'), split_record('b1')
        cross = a1[30:126] + b1[50:]                # S M1 C M2 U J (a1), then T (b1)
        for mode, rec in SPLIT_MODES:
            a, b, al, u = split_union(mode, rec)
            walked = cross in dag_spellings(al, 'right')
            # the unsplit annotate retrieval walks it; no label carries it in constrain mode
            self.assertEqual(walked, mode == 'annotate', (mode, rec))
            spellings = [u.spell('right', w.id, 'walk')
                         for w in u.walks('right', by='id', stubs=True)]
            self.assertNotIn(cross, spellings)
            self.assertFalse(any(s.startswith(cross) for s in spellings))
            # the statements, FASTA and GFA state the non-claim
            self.assertIn(U.SCOPE, u.coverage('right').sentence)
            self.assertEqual(u.statements('right').scope, U.SCOPE)
            heads = [h for h in u.to_fasta().split('\n') if h.startswith('>')]
            self.assertTrue(heads and all(h.endswith(' scope=within_one_graph') for h in heads))
            self.assertIn('# scope: ' + U.SCOPE, u.to_gfa().split('\n'))

    def kmers(self, headers):
        out = set()
        for h in headers:
            r = split_record(h)
            out.update(r[i:i + 15] for i in range(len(r) - 14))
        return out

    def test_annotate_every_walk_of_one_graph_is_in_the_union(self):
        graphs = {'A': self.kmers(['pl', 'a1']), 'B': self.kmers(['pl', 'b1', 'b2'])}
        for rec in ('merge', 'keep'):
            a, b, al, u = split_union('annotate', rec)
            seed = al.seed.sequence
            spelled = {u.spell('right', w.id): w for w in u.walks('right', by='id')}
            one_graph = 0
            for w in al.walks('right', by='id'):
                sp = al.spell('right', w.path_id)
                whole = seed + sp
                ks = {whole[i:i + 15] for i in range(len(whole) - 14)}
                for g, kmers in graphs.items():
                    if ks <= kmers:
                        one_graph += 1
                        self.assertIn(sp, spelled, (rec, w.path_id))
                        self.assertIn(g, spelled[sp].held_in)
            self.assertGreaterEqual(one_graph, 2)

    def test_constrain_walks_of_a_graphs_record_are_in_its_part(self):
        records = {'A': {'pl', 'a1'}, 'B': {'pl', 'b1', 'b2'}}
        for rec in ('merge', 'keep'):
            a, b, al, u = split_union('constrain', rec)
            parts = {'A': a, 'B': b}
            for w in al.walks('right', by='id'):
                sp = al.spell('right', w.path_id)
                for g, names in records.items():
                    if {l.name for l in w.labels_full} & names:
                        part = parts[g]
                        self.assertIn(sp, [part.spell('right', x.path_id)
                                           for x in part.walks('right', by='id')], (rec, g))


class TestAbsence(unittest.TestCase):
    """P10: each absence value of §2.4 on the union fixtures."""

    def spans(self, u, side, spelling):
        (w,) = [w for w in u.walks(side, by='id', stubs=True) if u.spell(side, w.id) == spelling]
        return {s.pair: s for s in u.carriers(side, w.id)}

    def test_graph_under_annotate_per_path(self):
        a, b, al, u = split_union('annotate', 'keep')
        # a1's allele C is not in graph B, nor U after M2 on pl's allele A
        s = self.spans(u, 'right', split_record('a1')[30:])['B']
        self.assertEqual((s.how, s.to_bp, s.absence), ('diverges', 25, 'graph'))
        s = self.spans(u, 'right', split_record('pl')[30:81] + split_record('a1')[81:])['B']
        self.assertEqual((s.how, s.to_bp, s.absence), ('diverges', 51, 'graph'))

    def test_under_strategy_in_constrain_and_under_united_histories(self):
        for mode, rec in (('constrain', 'keep'), ('constrain', 'merge'), ('annotate', 'merge')):
            a, b, al, u = split_union(mode, rec)
            # group A's 126 bp walk, which graph B does not follow
            (w,) = [w for w in u.walks('right', by='id') if w.length_bp == 126]
            s = {x.pair: x for x in u.carriers('right', w.id)}['B']
            self.assertEqual((s.how, s.absence), ('diverges', 'under_strategy'), (mode, rec))

    def test_dead_end_is_graph_also_at_the_certificate_edge(self):
        b = union_doc('b_annotate_keep')
        z = split_record('b2')[30:]
        z = split_record('pl')[30:] + z[len(split_record('pl')) - 30:]   # pl's allele, then Z
        a = union_doc('a_annotate_keep')
        r76 = union_doc('a_annotate_keep_r76')
        u = union_pairs([{'pair': 'A', 'graphlet': a}, {'pair': 'B', 'graphlet': b},
                         {'pair': 'R76', 'graphlet': r76}])
        s = self.spans(u, 'right', z)
        # A's graph ends where pl ends, B's goes on with Z: a dead end inside A's certificate
        self.assertEqual((s['A'].how, s['A'].to_bp, s['A'].end_reasons['path_reason'],
                          s['A'].absence), ('ends', 76, 'dead_end', 'graph'))
        # cut at the radius there: at the certificate's edge, uncertified
        self.assertEqual((s['R76'].how, s['R76'].to_bp, s['R76'].end_reasons['path_reason'],
                          s['R76'].absence), ('ends', 76, 'max_extension_bp', 'uncertified'))
        # a dead end at the certificate's edge is still the graph's (it comes first)
        edge = copy.deepcopy(a)
        edge.arms['right'].complete_to_bp = 76
        u = union_pairs([{'pair': 'A', 'graphlet': edge}, {'pair': 'B', 'graphlet': b}])
        s = self.spans(u, 'right', z)['A']
        self.assertEqual((s.to_bp, s.absence), (76, 'graph'))

    def test_constrain_dead_end_is_graph_and_a_hairpin_end_is_not(self):
        # constrain mode ends every label at a structural dead end with D and no path reason
        a, b = union_doc('a_constrain_keep'), union_doc('b_annotate_keep')
        z = split_record('pl')[30:] + split_record('b2')[len(split_record('pl')):]
        u = union_pairs([{'pair': 'A', 'graphlet': a}, {'pair': 'B', 'graphlet': b}])
        s = self.spans(u, 'right', z)['A']
        self.assertEqual((s.how, s.to_bp, s.end_reasons, s.absence),
                         ('ends', 76, {'path_reason': None, 'labels': {'dead_end': 1}}, 'graph'))
        self.assertIn('ends=A@76:dead_end/graph', u.to_fasta('right'))
        # the same leaf ended by a hairpin dead end (Dh): a hairpin successor exists there
        h = copy.deepcopy(a)
        arm = h.arms['right']
        (leaf,) = [x for x in arm.segments if x.leaf is not None and x.end_bp == 76]
        for e in derive.end_labels(arm, leaf.id):
            run = arm.runs[e.run]
            self.assertEqual(run.end, 'D')
            run.end, run.qualifier = 'Dh', 'hairpin'
        u = union_pairs([{'pair': 'A', 'graphlet': h}, {'pair': 'B', 'graphlet': b}])
        s = self.spans(u, 'right', z)['A']
        self.assertEqual((s.how, s.to_bp, s.absence), ('ends', 76, 'under_strategy'))
        # G.5's shape: that hairpin dead end skipped the very base the other pair spells
        from metagraph.traverse.model import Event
        leaf.events = list(leaf.events) + [Event(76, 'hairpin', char=z[76], followed=False)]
        u = union_pairs([{'pair': 'A', 'graphlet': h}, {'pair': 'B', 'graphlet': b}])
        s = self.spans(u, 'right', z)['A']
        self.assertEqual((s.how, s.to_bp, s.absence), ('ends', 76, 'walk_rule'))
        self.assertIn('ends=A@76:dead_end/walk_rule', u.to_fasta('right'))

    def test_walk_rule_from_a_skipped_hairpin_against_follow(self):
        # ambiguous_split (canonical, k 15) followed a hairpin at 25 bp of its left arm into
        # the child spelling 'A'; the same graph under hairpins: skip records the hairpin
        # (followed false) and has no such child (here its other child spells 'G' there)
        follow = T.graphlet('ambiguous_split')
        skip = copy.deepcopy(follow)
        segs = skip.arms['left'].segments
        (ev,) = [e for e in segs[0].events if e.type == 'hairpin']
        self.assertEqual((ev.at_bp, ev.char, ev.followed, segs[1].walk[0]), (25, 'A', True, 'A'))
        ev.followed = False
        segs[1].walk = 'G' + segs[1].walk[1:]
        u = union_pairs([{'pair': 'follow', 'graphlet': follow},
                         {'pair': 'skip', 'graphlet': skip}])
        whole = follow.spell('left', 0, 'walk')
        self.assertTrue(whole.startswith(segs[0].walk + 'A'))
        s = self.spans(u, 'left', whole[::-1])
        self.assertEqual((s['follow'].how, s['skip'].how, s['skip'].to_bp, s['skip'].absence),
                         ('whole', 'diverges', 25, 'walk_rule'))
        self.assertEqual((s['skip'].k, s['skip'].regime), (15, 'canonical'))
        # followed, the same event decides nothing: the pair has the child
        ev.followed = True
        u = union_pairs([{'pair': 'follow', 'graphlet': follow},
                         {'pair': 'skip', 'graphlet': skip}])
        self.assertNotEqual(self.spans(u, 'left', whole[::-1])['skip'].absence, 'walk_rule')

    def test_walk_rule_against_another_k(self):
        r11, r21 = union_doc('repeat11'), union_doc('repeat21')
        u = union_pairs([{'pair': 'k11', 'graphlet': r11}, {'pair': 'k21', 'graphlet': r21}])
        whole = r21.spell('right', 0)               # the record's 96 bp after the seed
        s = self.spans(u, 'right', whole)
        # k 11 meets the repeat's (k+1)-mers again and stops by edge reuse; k 21 spells on
        self.assertEqual((s['k11'].how, s['k11'].to_bp, s['k11'].end_reasons['path_reason'],
                          s['k11'].absence), ('ends', 67, 'edge_reuse', 'walk_rule'))
        self.assertEqual((s['k11'].k, s['k21'].k, s['k11'].regime), (11, 21, 'basic'))
        self.assertIn('ends=k11@67:edge_reuse/walk_rule', u.to_fasta('right'))

    def test_an_event_decides_only_its_own_character(self):
        r11, r21 = union_doc('repeat11'), union_doc('repeat21')
        whole = r21.spell('right', 0)
        other = copy.deepcopy(r11)
        for seg in other.arms['right'].segments:
            for ev in seg.events:
                if ev.type == 'blocked':
                    self.assertEqual(ev.char, whole[ev.at_bp])
                    # the event blocked another character: it decides nothing for this one
                    ev.char = next(c for c in 'ACGT' if c != ev.char)
        u = union_pairs([{'pair': 'k11', 'graphlet': other}, {'pair': 'k21', 'graphlet': r21}])
        s = self.spans(u, 'right', whole)['k11']
        self.assertEqual((s.to_bp, s.absence), (67, 'graph'))


class TestMixedFixtures(unittest.TestCase):
    """P11 on the union fixtures: parts at two k, and an arm one pair did not walk."""

    def test_two_k(self):
        u = union_pairs([{'pair': 'A15', 'graphlet': union_doc('a_annotate_merge')},
                         {'pair': 'A21', 'graphlet': union_doc('a21_annotate_merge')}])
        s = u.coverage('right').sentence
        self.assertTrue(s.startswith('right arm, to 200 bp: complete, qualified, in A15 '
                                     '[annotate, kmer, k 15] ('), s)
        self.assertIn('A21 [annotate, kmer, k 21] (', s)
        for w in u.walks('right', by='id', stubs=True):
            self.assertEqual({sp.pair: sp.k for sp in u.carriers('right', w.id)},
                             {'A15': 15, 'A21': 21})
        # one record at two k: the same walks, each whole in both pairs
        for w in u.walks('right', by='id'):
            self.assertEqual(w.whole_in, ('A15', 'A21'))

    def test_arm_not_walked(self):
        u = union_pairs([{'pair': 'A', 'graphlet': union_doc('a_annotate_merge')},
                         {'pair': 'B', 'graphlet': union_doc('b_annotate_merge_right')}])
        c = u.coverage('left')
        self.assertEqual(c.per_pair['B'], {'category': 'arm_not_walked', 'certified_to_bp': 0,
                                           'note': 'direction right'})
        for w in u.walks('left', by='id'):
            s = {x.pair: x for x in u.carriers('left', w.id)}['B']
            self.assertEqual((s.how, s.to_bp, s.absence), ('none', 0, 'under_strategy'))

    def test_union_of_each_fixture_is_that_fixture(self):
        # P2 on the union fixtures
        for name in UNION_DOCS:
            g = union_doc(name)
            u = union_pairs([g])
            for side in g.arms:
                parts = sorted(g.spell(side, w.path_id) for w in g.walks(side, by='id'))
                self.assertEqual(sorted(u.spell(side, w.id) for w in u.walks(side, by='id')),
                                 parts, (name, side))


R237 = os.path.join(UNION, 'r237')
# the body sha256 of each part as the service's record states it (the README of r237/)
R237_BODIES = {'bact_other': 'cda902d0df71afab3d9fe1bcf33a04630468c22a1748777941686867fcc31eea',
               'bact_pseudomonadota':
                   '6d48f657aa3af50bd79c3ed39126d42c4c8c1e1bc431f8f9d6d5155cc732d522'}


def comb_part(n):
    """test_traverse_scale's comb (n splits, each into a 1-bp leaf and the next spine
    segment) with the merge fixture's envelope and seed, so it unions."""
    from metagraph.traverse import parse
    from tests.test_traverse_scale import comb_annotate
    ref = T.graphlet('merge')
    g = parse(comb_annotate(n))
    g.envelope = copy.deepcopy(ref.envelope)
    g.seed_summary = copy.deepcopy(ref.seed_summary)
    return g


def spine_part(n):
    """The comb's flank in one segment: comb_part(n)'s spine 'A' + 'G' * n as a single
    walk (one segment, one P run), the shape of the GPT review's M5 (another part segmenting
    the same flank base by base)."""
    from metagraph.traverse import parse
    from tests.test_traverse_scale import comb_annotate
    ref = T.graphlet('merge')
    lines = comb_annotate(0).splitlines()
    for i, line in enumerate(lines):
        f = line.split()
        if f[0] == 'A':
            f[3] = f[-1] = str(n + 1)
        elif f[0] == 'G':
            f[3], f[-1] = str(n + 1), 'A' + 'G' * n
        elif f[0] == 'P':
            f[2] = str(n + 1)
        else:
            continue
        lines[i] = ' '.join(f)
    g = parse('\n'.join(lines) + '\n')
    g.envelope = copy.deepcopy(ref.envelope)
    g.seed_summary = copy.deepcopy(ref.seed_summary)
    return g


class TestScale(unittest.TestCase):
    """P5: construction, walks(top=3) and print() linear on a comb, n = 250, 1000 and 2000,
    checked as ratios to n = 250, never as seconds: a linear cost grows 4 and 8 times, a
    quadratic one 16 and 64 times; the bound is three times linear. Timed in process CPU time
    (time.process_time()): wall time on a loaded machine counts the time others run. A second part that
    stops after 2 bp makes every walk's figures carry a stop below the walk's end."""

    SIZES = (250, 1000, 2000)

    @staticmethod
    def best(fn, reps=5):
        out = None
        for _ in range(reps):
            dt = fn()
            out = dt if out is None else min(out, dt)
        return out

    def test_comb_is_linear(self):
        import gc
        import time
        short = comb_part(1)
        cost = {}
        for n in self.SIZES:
            parts = [{'pair': 'C', 'graphlet': comb_part(n)}, {'pair': 'D', 'graphlet': short}]

            def timed(query):
                u = union_pairs(parts)
                t0 = time.process_time()
                query(u)
                return time.process_time() - t0

            def build():
                t0 = time.process_time()
                union_pairs(parts)
                return time.process_time() - t0

            gc.disable()
            try:
                cost[n] = (self.best(build),
                           self.best(lambda: timed(lambda u: u.walks('right', top=3))),
                           self.best(lambda: timed(str)))
            finally:
                gc.enable()
            u = union_pairs(parts)
            top = u.walks('right', top=3)
            # D's two walks first (held by both pairs, one whole in both), then C's longest
            self.assertEqual([(w.length_bp, w.held_in, w.whole_in) for w in top],
                             [(2, ('C', 'D'), ('C', 'D')), (2, ('C', 'D'), ('D',)),
                              (n + 1, ('C',), ('C',))])
            self.assertIn('walks: right %d (2 held by two or more pairs)' % (n + 2), str(u))
        base = cost[self.SIZES[0]]
        for n in self.SIZES[1:]:
            bound = 3.0 * n / self.SIZES[0]
            for what, t, t0 in zip(('construction', 'walks(top=3)', 'print'), cost[n], base):
                self.assertLess(t / t0, bound, '%s at n=%d: %.4f s against %.4f s at n=%d'
                                % (what, n, t, t0, self.SIZES[0]))

    def test_one_flank_segmented_two_ways_is_linear(self):
        """GPT review M5: one part holds the comb's spine as one segment, the other is the
        comb (a segment per base of it), 752 / 1,502 / 3,002 / 6,002 bases. The insertion
        compares each segment with the edges it meets in place: the bases it reads (the
        trie's count) and the bases _bases() slices stay within the input's, and the CPU
        time grows linearly (ratios to n = 250 below three times linear). At these sizes a
        quadratic copy costs little time, so the counts are the test; the time a guard."""
        import gc
        import time
        sizes = (250, 500, 1000, 2000)
        cost, sliced = {}, {}
        orig = U._Trie._bases

        def counted(trie, e):
            text = orig(trie, e)
            sliced[n] = sliced.get(n, 0) + len(text)
            return text

        for n in sizes:
            parts = [{'pair': 'A', 'graphlet': spine_part(n)},
                     {'pair': 'B', 'graphlet': comb_part(n)}]
            bases = sum(s.length_bp for p in parts for s in p['graphlet'].arms['right'].segments)
            self.assertEqual(bases, 3 * n + 2)
            U._Trie._bases = counted
            try:
                u = union_pairs(parts)
            finally:
                U._Trie._bases = orig
            self.assertLessEqual(sliced.get(n, 0), bases, n)
            self.assertLessEqual(u._arms['right'].bases_read, 2 * bases, n)
            (w,) = [w for w in u.walks('right') if w.held_in == ('A', 'B')]
            self.assertEqual((w.length_bp, w.whole_in, u.spell('right', w.id)),
                             (n + 1, ('A', 'B'), 'A' + 'G' * n))

            def build():
                t0 = time.process_time()
                union_pairs(parts)
                return time.process_time() - t0

            gc.disable()
            try:
                cost[n] = self.best(build)
            finally:
                gc.enable()
        for n in sizes[1:]:
            self.assertLess(cost[n] / cost[sizes[0]], 3.0 * n / sizes[0],
                            'n=%d: %.4f s against %.4f s at n=%d'
                            % (n, cost[n], cost[sizes[0]], sizes[0]))


class _FanOutClient:
    """A TraverseClient whose server is a dict: GET /capabilities, the index's own
    capabilities, and per graph what traverse() answers or raises."""

    @staticmethod
    def make(caps, answers, ns=None, own=(None, None, None)):
        from metagraph.traverse import TraverseClient

        class Client(TraverseClient):
            sent = []

            def server_capabilities(self):
                return caps

            def capabilities(self, *a, **kw):
                return {'index_ns': ns}

            def selection(self, *a, **kw):
                return own

            def traverse(self, seeds, strategy=None, **kw):
                self.sent.append((seeds, kw['graph']))
                got = answers[kw['graph']]
                if isinstance(got, BaseException):
                    raise got
                return got
        return Client('localhost', 1)


class TestTraversePairs(unittest.TestCase):
    """TraverseClient.traverse_pairs() (§5.2) against a fake server: the pairs from GET
    /capabilities, what becomes {pair, error} and what is raised (review minor 6)."""

    CAPS = {'graph_summary': {'pairs': [{'graph': 'g1', 'seqs_order': 'verified'},
                                        {'graph': 'g2'}, {'graph': 'g3'}]}}

    def test_every_listed_pair_and_its_failures(self):
        from metagraph.traverse import TraverseError
        refusal = TraverseError(409, 'the attempt is busy')
        cut = ConnectionResetError('connection reset by peer')
        c = _FanOutClient.make(self.CAPS, {'g1': 'answer-g1', 'g2': refusal, 'g3': cut})
        r = c.traverse_pairs('ACGT')
        self.assertEqual(r.results, ['answer-g1', {'pair': {'graph': 'g2'}, 'error': refusal},
                                     {'pair': {'graph': 'g3'}, 'error': cut}])
        self.assertEqual((r.seed, r.expected_pairs, r.seqs_order),
                         ('ACGT', ({'graph': 'g1'}, {'graph': 'g2'}, {'graph': 'g3'}),
                          {('g1', None, None): 'verified'}))
        self.assertEqual([g for _, g in c.sent], ['g1', 'g2', 'g3'])

    def test_a_bad_argument_or_a_bug_is_raised(self):
        for exc in (ValueError('strategy.mode is not a mode'), KeyError('oops'),
                    TypeError('a library bug')):
            c = _FanOutClient.make(self.CAPS, {'g1': 'answer-g1', 'g2': exc, 'g3': 'x'})
            with self.assertRaises(type(exc)):
                c.traverse_pairs('ACGT')

    def test_a_single_graph_server(self):
        caps = {'graph_summary': {}, 'pair_summary': {'seqs_order': 'count_only'}}
        c = _FanOutClient.make(caps, {None: 'answer'}, ns='refseq_ns')
        r = c.traverse_pairs('ACGT')
        self.assertEqual((r.results, r.expected_pairs, r.seqs_order),
                         (['answer'], ({'index_ns': 'refseq_ns'},), {'refseq_ns': 'count_only'}))
        c = _FanOutClient.make(caps, {None: OSError('timed out')}, ns='refseq_ns')
        r = c.traverse_pairs('ACGT')
        self.assertEqual(r.results[0]['pair'], {'index_ns': 'refseq_ns'})

    def test_an_empty_selection_sends_nothing(self):
        # GPT review m3: pairs=[] selects no pair; only pairs=None discovers them
        for caps in (self.CAPS, {'graph_summary': {}}):
            c = _FanOutClient.make(caps, {None: 'unexpected request'}, ns='index')
            asked = []
            c.server_capabilities = lambda: asked.append(1) or caps
            r = c.traverse_pairs('ACGT', pairs=[])
            self.assertEqual((c.sent, asked), ([], []))
            self.assertEqual((r.results, r.seed, r.expected_pairs, r.seqs_order),
                             ([], 'ACGT', (), {}))
        with self.assertRaises(ValueError):
            _FanOutClient.make(self.CAPS, {}).traverse_pairs('ACGT', pairs='g1')

    def test_a_malformed_answer_is_that_pairs_error(self):
        # GPT review m4: an answer that is not a response object (a JSON array, results
        # that are not objects) is that pair's ResponseFormatError; the other pairs are
        # still requested, and a programming error is still raised. Opus fix review N3: a
        # graphlet that is not text and a refusal that is not an object too
        import metagraph.traverse as MT
        from metagraph.traverse import TraverseClient, TraverseResponse
        # (the base had no ResponseFormatError: the test then fails on the fan-out itself)
        ResponseFormatError = getattr(MT, 'ResponseFormatError', MT.TraverseError)
        caps = self.CAPS

        class Client(TraverseClient):
            def __init__(self, answers):
                super().__init__('localhost', 1)
                self.answers, self.sent = answers, []

            def server_capabilities(self):
                return caps

            def traverse_raw(self, request):
                self.sent.append(request['graph'])
                got = self.answers[request['graph']]
                if isinstance(got, BaseException):
                    raise got
                return got

        c = Client({'g1': [], 'g2': {'results': ['not a result']}, 'g3': {'results': []}})
        r = c.traverse_pairs('ACGT', coordinates=None)
        self.assertEqual(c.sent, ['g1', 'g2', 'g3'])
        (e1, e2), ok = [x['error'] for x in r.results[:2]], r.results[2]
        self.assertEqual([x['pair'] for x in r.results[:2]], [{'graph': 'g1'}, {'graph': 'g2'}])
        self.assertIsInstance(e1, ResponseFormatError)
        self.assertIsInstance(e2, ResponseFormatError)
        self.assertEqual(str(e1), 'not a /traverse response: the answer is a JSON array, not '
                         'an object')
        self.assertEqual(e1.body, [])
        self.assertIsInstance(ok, TraverseResponse)
        c = Client({'g1': {'results': []}, 'g2': KeyError('a library bug'), 'g3': []})
        with self.assertRaises(KeyError):
            c.traverse_pairs('ACGT', coordinates=None)
        shapes = {
            'results[0].graphlet is a JSON number, not a string': {'results': [{'graphlet': 5}]},
            'results[0].graphlet is a JSON null, not a string': {'results': [{'graphlet': None}]},
            "the refused entry's refusal is a JSON string, not an object":
                {'envelope_version': 1, 'answers': [{'outcome': 'refused', 'refusal': 'no'}]},
            "the refused entry's refusal is a JSON array, not an object":
                {'envelope_version': 1, 'answers': [{'outcome': 'refused', 'refusal': ['x']}]}}
        for text, bad in shapes.items():
            c = Client({'g1': bad, 'g2': {'results': []}, 'g3': {'results': []}})
            r = c.traverse_pairs('ACGT', coordinates=None)
            self.assertEqual(c.sent, ['g1', 'g2', 'g3'], text)
            e = r.results[0]['error']
            self.assertIsInstance(e, ResponseFormatError, text)
            self.assertEqual(str(e), 'not a /traverse response: ' + text)
            self.assertEqual(e.body, bad)


class TestLoadGzip(unittest.TestCase):
    """load() reads a gzip file (the service's .mgt.gz export) as the text it holds,
    inflated as a stream within the local limits; a corrupt or truncated one is refused."""

    def setUp(self):
        import tempfile
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def write(self, name, data):
        path = os.path.join(self.tmp.name, name)
        with open(path, 'wb') as f:
            f.write(data)
        return path

    def test_gzip_equals_the_text(self):
        import gzip
        from metagraph.traverse import load, parse
        for name in ('merge', 'annotate'):
            text = T.doc_text(name)
            plain = self.write(name + '.mgt', text.encode('utf-8'))
            gz = self.write(name + '.mgt.gz', gzip.compress(text.encode('utf-8')))
            want = load(plain).dump(envelope=True)
            self.assertEqual(load(gz).dump(envelope=True), want)
            self.assertEqual(load(gz, budget=LocalBudget(work_units=10 ** 9)).dump(envelope=True),
                             want)
        for part in R237_BODIES:
            path = os.path.join(R237, part + '.mgt.gz')
            with open(path, 'rb') as f:
                text = gzip.decompress(f.read()).decode('utf-8')
            self.assertEqual(load(path).dump(envelope=True), parse(text).dump(envelope=True))

    def test_save_writes_plain_text(self):
        from metagraph.traverse import load, save
        g = load(os.path.join(R237, 'bact_other.mgt.gz'))
        out = os.path.join(self.tmp.name, 'out.mgt')
        save(g, out)
        with open(out, 'rb') as f:
            self.assertEqual(f.read(2), b'H ')

    def test_a_high_ratio_gzip_stops_at_the_default_bound(self):
        import gzip
        from metagraph.traverse import load
        from metagraph.traverse import parser as P
        text = T.doc_text('merge')
        # 64 MiB of one byte in a comment-free document: about 64 KiB of gzip, so the
        # default bound is GZIP_RATIO_FLOOR (16 MiB), far below GZIP_TEXT_MAX: not patched
        bomb = self.write('bomb.mgt.gz', gzip.compress(text.encode('utf-8') + b'A' * (64 << 20)))
        size = os.path.getsize(bomb)
        self.assertLess(size * P.GZIP_RATIO_MAX, P.GZIP_RATIO_FLOOR)
        self.assertEqual((P.GZIP_TEXT_MAX, P.GZIP_RATIO_MAX, P.GZIP_RATIO_FLOOR),
                         (512 << 20, 64, 16 << 20))
        with self.assertRaises(GraphletFormatError) as cm:
            load(bomb)
        self.assertEqual(str(cm.exception),
                         'the gzip file bomb.mgt.gz (%d bytes) inflates to more than 16777216 '
                         'bytes, the most a load() without budget= inflates (GZIP_TEXT_MAX, or '
                         'GZIP_RATIO_MAX times the file beyond GZIP_RATIO_FLOOR): refused, not '
                         'inflated further; pass budget= to bound it by its memory account '
                         'instead' % size)
        # GZIP_TEXT_MAX caps whatever the ratio allows (lowered here to see it act)
        old = P.GZIP_TEXT_MAX
        P.GZIP_TEXT_MAX = 1 << 20
        try:
            with self.assertRaises(GraphletFormatError) as cm:
                load(bomb)
        finally:
            P.GZIP_TEXT_MAX = old
        self.assertIn('inflates to more than 1048576 bytes', str(cm.exception))
        # with budget=, the memory account stops it as it stops a plain text that size
        plain = self.write('big.mgt', text.encode('utf-8') + b'A' * (64 << 20))
        stops = []
        for path in (bomb, plain):
            with self.assertRaises(LocalBudgetExceeded) as cm:
                load(path, budget=LocalBudget(memory_mb=16))
            stops.append((cm.exception.stop.resource, cm.exception.stop.op,
                          cm.exception.stop.limit))
        self.assertEqual(stops[0], stops[1])
        self.assertEqual(stops[0][0], 'memory')

    def test_real_exports_load_far_below_the_ratio(self):
        import gzip
        from metagraph.traverse import parser as P
        for name in ('bact_other', 'bact_pseudomonadota'):
            path = os.path.join(R237, name + '.mgt.gz')
            with gzip.open(path, 'rb') as f:
                ratio = len(f.read()) / os.path.getsize(path)
            self.assertLess(ratio, P.GZIP_RATIO_MAX / 8, name)

    def test_corrupt_or_truncated_gzip_refused(self):
        import gzip
        from metagraph.traverse import load
        data = gzip.compress(T.doc_text('merge').encode('utf-8'))
        bad_crc = bytearray(data)
        bad_crc[-6] ^= 0xff
        for name, blob in (('cut.mgt.gz', data[:-20]), ('crc.mgt.gz', bytes(bad_crc)),
                           ('tail.mgt.gz', data + b'junk')):
            with self.subTest(name=name):
                with self.assertRaises(GraphletFormatError) as cm:
                    load(self.write(name, blob))
                self.assertTrue(str(cm.exception).startswith(
                    'the gzip file %s is corrupt or truncated (' % name), str(cm.exception))
        # a gzip of text that is not UTF-8: the parser's refusal
        with self.assertRaises(GraphletFormatError) as cm:
            load(self.write('latin.mgt.gz', gzip.compress(b'H \xe9\n')))
        self.assertIn('not UTF-8', str(cm.exception))


def r237_union(**kw):
    """The r237 export (documents/union/r237/): the two parts as delivered, the nine absent
    records, seed= from seed.fa and the 11 pairs expected."""
    from metagraph.traverse import load
    parts = [load(os.path.join(R237, p + '.mgt.gz')) for p in sorted(R237_BODIES)]
    with open(os.path.join(R237, 'absent.json'), encoding='utf-8') as f:
        absent = json.load(f)
    with open(os.path.join(R237, 'seed.fa'), encoding='utf-8') as f:
        seed = ''.join(line.strip() for line in f if not line.startswith('>'))
    expected = [{k: g.envelope[k] for k in ('graph', 'graph_path', 'annotation_path')}
                for g in parts] + [dict(a['pair']) for a in absent]
    return parts, union_pairs(parts + absent, seed=seed, expected_pairs=expected, **kw)


class TestR237RoundTrip(unittest.TestCase):
    """OQ21: the service's r237 export of one 16S seed over the 11 graphs (2 answered, 9
    not_in_graph) loads as delivered and unites as §2.3 and §2.6 expect."""

    @classmethod
    def setUpClass(cls):
        cls.parts, cls.u = r237_union()

    def test_a_merge_that_splits_at_once(self):
        # bact_pseudomonadota's stubs at 35 bp (left) and 52 bp (right) enter zero-length
        # merge segments with two children: the trie node is the end of the union segment
        # there, 's<id>$', in the span, the merges token and mt:Z
        u = self.u
        got = []
        for side in u.arms():
            for w in u.walks(side, stubs=True, by='id'):
                for s in u.carriers(side, w.id):
                    if s.how == 'merges' and s.into.endswith('$'):
                        a = u._part_of[s.pair].graphlet.arms[side]
                        (m,) = [c for c in a.segments[s.stub_end].children
                                if a.segments[c].parents[0] != s.stub_end]
                        got.append((side, s.pair, s.to_bp, s.into, a.segments[m].length_bp,
                                    len(a.segments[m].children), s.route.how,
                                    s.route.goes_on))
        self.assertEqual(got, [('left', 'bact_pseudomonadota', 35, 's12$', 0, 2, 'held', True),
                               ('right', 'bact_pseudomonadota', 52, 's61$', 0, 2, 'held',
                                True)])
        fa, gfa = u.to_fasta(stubs=True), u.to_gfa()
        self.assertIn('merges=bact_pseudomonadota@35>s12$', fa)
        self.assertIn('merges=bact_pseudomonadota@52>s61$', fa)
        self.assertIn('\tmt:Z:s12$', gfa)
        self.assertIn('\tmt:Z:s61$', gfa)
        check_gfa_paths(self, u, gfa)

    def test_parts_load_as_delivered(self):
        import hashlib
        for g in self.parts:
            name = g.envelope['graph']
            self.assertEqual(hashlib.sha256(g.dump(envelope=False).encode('utf-8')).hexdigest(),
                             R237_BODIES[name])
            self.assertEqual((g.mode, g.k, sorted(g.arms)), ('constrain', 31, ['left', 'right']))
            self.assertEqual(g.seed.sequence, self.u.seed.sequence)

    def test_union_states_the_pairs(self):
        u = self.u
        absent = ['inv_insecta', 'inv_other', 'plant', 'vm_laurasiatheria', 'vm_other',
                  'vo_actinopterygii', 'vo_amphibia_lobefins', 'vo_other', 'vo_sauropsida']
        present = {'inv_insecta': 23, 'plant': 15, 'vm_laurasiatheria': 31, 'vo_sauropsida': 13}
        self.assertEqual(([p.tag for p in u.parts], [a.tag for a in u.absent], u._not_included),
                         (['bact_other', 'bact_pseudomonadota'], absent, ()))
        self.assertEqual(str(u).split('\n')[0], 'union of one seed (200 bp) over 11 expected '
                                                'pairs: 2 parts, 9 absent records')
        seeds = ', '.join('%s (%d of 170 k-mers)' % (t, present.get(t, 0)) for t in absent)
        self.assertEqual(
            U._core_sentence(u.coverage('left')),
            'left arm, to 5000 bp (constrain, kmer, k 31, in every pair): incomplete in '
            'bact_other (certified to 342 bp: time_budget), bact_pseudomonadota (certified to '
            '55 bp: resource_limit); seed not in graph: %s. Pairs stated by the caller: %s. '
            'Seed of the absent records: not checked for %s.'
            % (seeds, ', '.join(absent), ', '.join(absent)))
        self.assertEqual(u.statements('left').provenance,
                         'walks: left 11 (2 held by two or more pairs; 9 stub walks not listed)')
        self.assertEqual(u.statements('right').provenance,
                         'walks: right 9 (3 held by two or more pairs; 6 stub walks not listed)')
        self.assertTrue(u.statements('left').labels.endswith(
            'header names of every pair rest on an unchecked record mapping (unknown)'))
        for side in ('left', 'right'):
            c = u.coverage(side, at_bp=50)
            self.assertEqual(c.per_pair['bact_other'],
                             {'category': 'complete_qualified',
                              'certified_to_bp': {'left': 342, 'right': 341}[side],
                              'note': 'a merge united histories: complete for the '
                                      'united-history rule, not per path'})
            self.assertEqual(c.per_pair['bact_pseudomonadota'],
                             {'category': 'complete_qualified', 'certified_to_bp': 55,
                              'note': 'complete for 1000 of 12313 seed labels; a merge united '
                                      'histories: complete for the united-history rule, not '
                                      'per path'})
            self.assertEqual(c.by_category['seed_not_in_graph'], absent)
            self.assertIsNone(c.every_pair_certified_to_bp)

    def test_bact_other_extends_the_walks_both_hold(self):
        """The walks both pairs hold are bact_pseudomonadota's whole 55 bp walks (its
        certified depth), which bact_other spells and goes on beyond."""
        for side, ids in (('left', [3, 5]), ('right', [3, 5, 9])):
            both = self.u.walks(side, pairs=['bact_other', 'bact_pseudomonadota'])
            self.assertEqual([w.id for w in both], ids)
            for w in both:
                self.assertEqual((w.length_bp, w.whole_in), (55, ('bact_pseudomonadota',)))
                self.assertEqual({s.pair: (s.how, s.to_bp) for s in w.spans},
                                 {'bact_other': ('extends', 55),
                                  'bact_pseudomonadota': ('whole', 55)})

    def test_the_union_of_each_part_is_that_part(self):
        for g in self.parts:
            u = union_pairs([g])
            for side in g.arms:
                walks = {g.spell(side, w.path_id): w for w in g.walks(side, by='id')}
                mine = u.walks(side, by='id')
                self.assertEqual(sorted(walks), sorted(u.spell(side, w.id) for w in mine))
                for w in mine:
                    (sp,) = u.carriers(side, w.id)
                    pw = walks[u.spell(side, w.id)]
                    self.assertEqual((sp.how, sp.walk, sp.to_bp),
                                     ('whole', pw.path_id, pw.length_bp))
                    self.assertEqual({x.label.ref for x in sp.labels if x.from_bp == 0},
                                     {l.ref for l in pw.labels_full})

    def test_exports_read_the_union(self):
        fasta = self.u.to_fasta()
        heads = [line for line in fasta.split('\n') if line.startswith('>')]
        self.assertEqual(len(heads), 20)
        self.assertTrue(all('absent=inv_insecta:seed_not_in_graph,' in h and
                            h.endswith(' mapping=unknown scope=within_one_graph') for h in heads))
        self.assertTrue(self.u.to_gfa().startswith('H\tVN:Z:1.0\n# union_pairs: seed 200 bp; '
                                                   'parts bact_other bact_pseudomonadota; '))


if __name__ == '__main__':
    unittest.main()
