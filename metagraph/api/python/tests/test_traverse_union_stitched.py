"""A StitchedGraphlet as a union part (docs/DESIGN-graphlet-merge.md §2.3, train 7): the
stitching triples of documents/stitch/ (scripts/traversal/graphlet_fixtures.py
STITCH_FIXTURES), the stitched graphlet as one pair's part beside its reference as another's.
Where no junction is read the union answers as with the root alone; past a junction nothing
reads complete or exact."""

import gc
import json
import re
import tracemalloc
import unittest

from tests import traverse_testlib as T  # noqa: F401  (puts the library on sys.path)

from metagraph.traverse import (CombineRefused, LocalBudget, SpanJunction, StitchedGraphlet,
                                derive, stitch, union_pairs)
from tests.test_traverse_union import admitted_paths, dag_spellings, gfa_segment_pairs

# the triples a union reads: every one but those without bases (refused no_bases)
NO_BASES = ('no_bases', 'no_bases_parent')
NAMES = [n for n in T.stitch_names() if n not in NO_BASES]


def role(name, r):
    try:
        return T.stitch_graphlet(name, r)
    except FileNotFoundError:
        return None


def stitched(name):
    """The triple's stitched graphlet: the parent and its child, and the grandchild or the
    sibling where the triple has one (fresh models every call)."""
    sg = stitch(role(name, 'parent'), role(name, 'child'))
    for extra in ('grandchild', 'sibling'):
        g = role(name, extra)
        if g is not None:
            sg = sg.stitch(g)
    return sg


def pairs(a, ref):
    return [{'pair': 'A', 'graphlet': a}] + ([] if ref is None else
                                              [{'pair': 'B', 'graphlet': ref}])


def answers(u, names=None):
    """Every answer of union |u| as text, by query. |names|: a function renaming the trie
    nodes a span names ('s31', 's30$'), for unions whose segment ids differ."""
    out = {'str': str(u), 'summary': json.dumps(u.summary(), sort_keys=True),
           'fasta': u.to_fasta(stubs=True), 'gfa': u.to_gfa(labels='full'),
           'gfa_counts': u.to_gfa(), 'evidence': repr(u.evidence()),
           'labels': repr(u.labels()), 'shared': repr(u.shared_names()),
           'provenance': repr([{k: v for k, v in e.items() if k != 'digest'}
                               for e in u.provenance()])}
    for side in u.arms():
        out[side + ' walks'] = repr(u.walks(side, stubs=True, spell=True))
        out[side + ' segments'] = repr(u.segments(side))
        out[side + ' coverage'] = repr(u.coverage(side))
        out[side + ' statements'] = repr(u.statements(side))
        for w in u.walks(side, stubs=True):
            out['%s carriers %d' % (side, w.id)] = repr(u.carriers(side, w.id))
            out['%s carriers %d at 10' % (side, w.id)] = repr(u.carriers(side, w.id, at_bp=10))
    return out


def node_names(u, side):
    """'s<id>' / 's<id>$' -> the spelling of that trie node (the union segment's start or
    end), so that spans of two unions with other segment ids compare."""
    t = u._arms[side]
    out = {'seed': ''}
    for us in t.useg:
        end = u._spell_nodes(t, t.path_usegs(us.end), 'walk', False)
        out['s%d$' % us.id] = end
        out['s%d' % us.id] = end[:len(end) - (us.to_bp - us.from_bp)]
    return out


def span_key(sp, names):
    """A span as text with its trie nodes by spelling."""
    d = {f: getattr(sp, f) for f in sp.__slots__}
    if d['into'] is not None:
        d['into'] = names[d['into']]
    if d['route'] is not None:
        r = {f: getattr(sp.route, f) for f in sp.route.__slots__}
        r['merges'] = tuple((at, names[n]) for at, n in r['merges'])
        d['route'] = r
    return repr(d)


def stitched_spellings(sg, side):
    """dag_spellings() of a stitched graphlet: a root-to-leaf spelling of a part that ends
    at the walk a child continues goes on with each of the child's (its seed, the
    continuation, not spelled twice)."""
    def of(part):
        a = part.graphlet.arms.get(side)
        if a is None:
            return []
        segs = a.segments
        kids = {s.id: [] for s in segs}
        for s in segs:
            for p in s.parents:
                kids[p].append(s.id)
        cont = {derive.paths(a)[q.key[-1][1]].leaf: q for q in sg.parts
                if q.key and q.key[:-1] == part.key and q.key[-1][0] == side}
        out = []

        def go(i, acc):
            acc += segs[i].walk
            if not kids[i]:
                q = cont.get(i)
                out.extend([acc] if q is None else [acc + c for c in of(q)])
            for c in kids[i]:
                go(c, acc)

        for s in segs:
            if not s.parents:
                go(s.id, '')
        return out
    return of(sg.parts[0])


def check_gfa(test, u, gfa):
    """test_traverse_union.check_gfa_paths() with a stitched part's spellings
    (stitched_spellings()): every path the file admits, and every part of it from the seed
    out to one of its segments, is spelled whole by each pair spelling that segment."""
    segs = gfa_segment_pairs(gfa)
    seed = u.seed.sequence
    refs = {}
    for p in u.parts:
        g = p.graphlet
        stitched_ = isinstance(g, StitchedGraphlet)
        for side in (g.arms() if stitched_ else g.arms):
            ws = stitched_spellings(g, side) if stitched_ else dag_spellings(g, side)
            refs[p.tag, side] = [seed + w if side == 'right' else w[::-1] + seed for w in ws]
    paths = admitted_paths(gfa, names=True)
    test.assertTrue(paths)
    for spelled, via in paths:
        side = 'right' if 'seed_r' in via else 'left'
        for k, name in enumerate(via):
            if name in ('seed_l', 'seed_r'):
                continue
            part = ''.join(segs[x][0] for x in (via[:k + 1] if side == 'right' else via[k:]))
            for tag in segs[name][1]:
                test.assertTrue(any(r.startswith(part) if side == 'right' else r.endswith(part)
                                    for r in refs[tag, side]), (tag, side, via))


def walks_by_spelling(u, side):
    return {(w.kind, u.spell(side, w.id, with_seed=True)): w
            for w in u.walks(side, stubs=True)}


class TestStitchedPart(unittest.TestCase):

    def test_stitch_of_one_graphlet_answers_as_it(self):
        """stitch(g), no junction, answers every union query as g, beside a plain part
        and alone; only its digest and retrieval identity are its own (the coordinator's
        Q5)."""
        for name in NAMES:
            with self.subTest(name):
                ref = role(name, 'reference')
                for with_ref in (ref, None):
                    p = role(name, 'parent')
                    sg = stitch(p)
                    a = union_pairs(pairs(p, with_ref))
                    b = union_pairs(pairs(sg, with_ref))
                    got, want = answers(b), answers(a)
                    for k in want:
                        self.assertEqual(got[k], want[k], k)
                    self.assertEqual(b.parts[0].digest, sg.digest)
                    self.assertNotEqual(b.digest, a.digest)
                    self.assertEqual(b.parts[0].identity, a.parts[0].identity)

    def test_the_union_walks_are_the_stitched_walks(self):
        """A union of the stitched part alone: one union walk per stitched walk, of its
        length, bases, certified_to_bp, complete and junctions, the walk named by its
        key."""
        for name in NAMES:
            with self.subTest(name):
                sg = stitched(name)
                u = union_pairs(pairs(sg, None))
                for side in sg.arms():
                    sws = {w.key: w for w in sg.walks(side, by='id')}
                    got = {}
                    for w in u.walks(side):
                        (sp,) = u.carriers(side, w.id)
                        self.assertEqual(sp.how, 'whole')
                        key = sp.walk if isinstance(sp.walk, tuple) else (sp.walk,)
                        got[key] = (w.length_bp, u.spell(side, w.id), sp.certified_to_bp,
                                    sp.complete, tuple(j.at_bp for j in sp.junctions))
                    self.assertEqual(got, {k: (w.length_bp, w.sequence, w.certified_to_bp,
                                               w.complete, tuple(sg.junctions[i].at_bp
                                                                 for i in w.junctions))
                                           for k, w in sws.items()})

    def test_nothing_reads_complete_or_exact_past_a_junction(self):
        """Every span that reads a junction names it, is not complete, is certified at most
        to the first junction, states absence 'graph' or 'uncertified' only, and its labels
        past the junction start there (not joined across it). The statements, evidence,
        coverage and exports state the junctions."""
        seen = 0
        for name in NAMES:
            with self.subTest(name):
                sg = stitched(name)
                u = union_pairs(pairs(sg, role(name, 'reference')))
                fasta = u.to_fasta(stubs=True)
                heads = {line.split(' ')[0][1:]: line for line in fasta.split('\n')
                         if line.startswith('>')}
                for side in u.arms():
                    junctions = [j for j in sg.junctions if j.arm == side]
                    for w in u.walks(side, stubs=True):
                        for sp in u.carriers(side, w.id):
                            if sp.pair != 'A' or not sp.junctions:
                                continue
                            seen += 1
                            first = sp.junctions[0].at_bp
                            self.assertTrue(all(isinstance(j, SpanJunction)
                                                and j.verdict == 'qualified'
                                                and j.causes for j in sp.junctions))
                            self.assertIn(first, [j.at_bp for j in junctions])
                            self.assertFalse(sp.complete)
                            self.assertLessEqual(sp.certified_to_bp, first)
                            r = sp.route
                            absence = r.absence if r is not None else sp.absence
                            self.assertIn(absence, (None, 'graph', 'uncertified'))
                            if r is not None:
                                self.assertFalse(r.complete)
                                self.assertLessEqual(r.certified_to_bp, first)
                            last = sp.junctions[-1].at_bp
                            labs = r.labels if r is not None else sp.labels
                            end = r.to_bp if r is not None else sp.to_bp
                            if end > last:
                                self.assertTrue(all(x.from_bp >= last for x in labs), labs)
                            self.assertIn(' junctions=A@%d:qualified' % first,
                                          heads['%s.u%d' % (side, w.id)])
                    if junctions:
                        cov = u.coverage(side)
                        plain = union_pairs(pairs(role(name, 'parent'),
                                                  role(name, 'reference'))).coverage(side)
                        self.assertEqual(cov.per_pair['A']['category'],
                                         plain.per_pair['A']['category'])
                        self.assertIn('stitched (experimental): junction',
                                      cov.per_pair['A']['note'])
                        ev = u.evidence(side)
                        self.assertFalse(ev['exact'])
                        self.assertEqual([j['verdict'] for j in ev['parts']['A']['junctions']],
                                         ['qualified'] * len(junctions))
                self.assertIn('# stitched (experimental): A joins ', u.to_gfa())
                self.assertEqual(len(u.summary()['pairs']['A']['junctions']),
                                 len(sg.junctions))
                self.assertEqual(len(u.provenance()[0]['stitched']), len(sg.parts) - 1)
        self.assertGreater(seen, 40)

    def test_off_junction_answers_are_the_roots(self):
        """Beside the reference: every union walk on which the stitched part reads no
        junction is a walk of the union with the root instead, and every pair's span on it
        is the same (trie nodes compared by spelling); the reference's spans are the same
        on every walk both unions have, and so is its coverage."""
        for name in NAMES:
            with self.subTest(name):
                ref = role(name, 'reference')
                us = union_pairs(pairs(stitched(name), ref))
                up = union_pairs(pairs(role(name, 'parent'), role(name, 'reference')))
                for side in us.arms():
                    ns, np_ = node_names(us, side), node_names(up, side)
                    ws, wp = walks_by_spelling(us, side), walks_by_spelling(up, side)
                    for key, w in ws.items():
                        cs = us.carriers(side, w.id)
                        a = next(sp for sp in cs if sp.pair == 'A')
                        other = wp.get(key)
                        if not a.junctions:
                            self.assertIsNotNone(other, key)
                        if other is None:
                            continue
                        cp = up.carriers(side, other.id)
                        for x, y in zip(cs, cp):
                            if x.pair == 'B' or not a.junctions:
                                self.assertEqual(span_key(x, ns), span_key(y, np_),
                                                 (side, key[0], len(key[1]), x.pair))
                    if ref is not None:
                        self.assertEqual(us.coverage(side).per_pair['B'],
                                         up.coverage(side).per_pair['B'])

    def test_a_route_crossing_a_junction(self):
        """A route through the root's merge that reaches the walk a child continues goes on
        in the child: splits45's merge segment is the continued walk, and a reference whose
        first parent there is the stub spells the stub through the junction into the
        child."""
        p, c = role('splits45', 'parent'), role('splits45', 'child')
        ref = role('splits45', 'reference')
        merge = ref.arms['right'].segments[5]
        self.assertEqual(merge.parents, (2, 4))
        merge.parents = (4, 2)
        merge.partition = list(reversed(merge.partition))
        got = {}
        for label, a in (('stitched', stitch(p, c)), ('plain', p)):
            u = union_pairs(pairs(a, ref))
            (w,) = [w for w in u.walks('right') if w.length_bp == 50]
            got[label] = next(sp for sp in u.carriers('right', w.id) if sp.pair == 'A')
        sp = got['stitched']
        r = sp.route
        self.assertEqual((sp.how, sp.to_bp, r.how, r.to_bp, r.goes_on, r.complete,
                          r.certified_to_bp, sp.complete, sp.certified_to_bp),
                         ('merges', 31, 'held', 50, True, False, 45, False, 45))
        self.assertEqual(sp.junctions, (SpanJunction(45, 'qualified', (('right', 1),),
                                                     ('exact_gated',)),))
        self.assertTrue(r.labels and all(x.from_bp == 45 for x in r.labels))
        # the root alone: the route ends at its walk's end, uncertified there
        r = got['plain'].route
        self.assertEqual((r.how, r.to_bp, r.absence, got['plain'].junctions),
                         ('ends', 45, 'uncertified', ()))

    def test_segments_and_labels_name_the_sub_part(self):
        for name in ('chain3', 'sibling_merge', 'fork_left', 'comb_chain'):
            with self.subTest(name):
                sg = stitched(name)
                u = union_pairs(pairs(sg, role(name, 'reference')))
                keys = {p.key: p for p in sg.parts}
                for side in u.arms():
                    for s in u.segments(side):
                        for sup in s.support:
                            if sup.pair == 'A':
                                self.assertIn(sup.part, keys)
                            else:
                                self.assertEqual(sup.part, ())
                            off = keys[sup.part].offset_bp if sup.pair == 'A' else 0
                            self.assertEqual(sup.pieces[0][1], s.from_bp)
                            self.assertEqual(sup.pieces[-1][2], s.to_bp)
                            self.assertGreaterEqual(s.from_bp, off)
                            for iv in sup.labels:
                                self.assertTrue(s.from_bp <= iv.from_bp < iv.to_bp <= s.to_bp)
                    check_gfa(self, u, u.to_gfa())
                # labels: the stitched part's, by ref or name, local ids refused
                names = {l.name for l in sg.labels}
                self.assertEqual({l.name for l in u.labels('A')}, names)
                for n in sorted(names):
                    self.assertEqual(u.label({'pair': 'A', 'name': n}).name, n)
                gfa = u.to_gfa(labels='full')
                self.assertTrue(re.search(r'pl:Z:[^\t]*A=', gfa))

    def test_refusals(self):
        p = role('keep_forbid', 'parent')
        sg = stitch(p, role('keep_forbid', 'child'))
        with self.assertRaises(CombineRefused) as cm:
            union_pairs([{'pair': 'A', 'graphlet': sg}, {'pair': 'A', 'graphlet': p}])
        self.assertEqual(cm.exception.code, 'duplicate_pair')
        with self.assertRaises(CombineRefused) as cm:
            union_pairs([{'pair': 'A', 'graphlet': stitch(p)}, {'pair': 'A', 'graphlet': p}])
        self.assertEqual(cm.exception.code, 'duplicate_pair')
        # the same stitched graphlet twice: a no-op
        u = union_pairs([{'pair': 'A', 'graphlet': sg}, {'pair': 'A', 'graphlet': sg}])
        self.assertEqual(len(u.parts), 1)
        for name in NO_BASES:
            with self.assertRaises(CombineRefused) as cm:
                union_pairs(pairs(stitched(name), None))
            self.assertEqual(cm.exception.code, 'no_bases')
        with self.assertRaises(CombineRefused) as cm:
            union_pairs([42])
        self.assertEqual(cm.exception.code, 'not_combinable')
        self.assertTrue(str(cm.exception).startswith(
            'union_pairs() takes results: a Graphlet, a StitchedGraphlet, a TraverseResponse,'))
        # a bare stitched graphlet: its pair is its root's, as for the root itself
        outcome = []
        for a in (sg, p):
            try:
                union_pairs([a])
                outcome.append(None)
            except CombineRefused as e:
                outcome.append(e.code)
        self.assertEqual(outcome[0], outcome[1])

    def test_the_account_bounds_the_peak(self):
        """budget.py: "the account bounds the call's peak" with a stitched part: the
        pattern of test_traverse_union's test (fresh models and a fresh union, the free
        lists filled by the same call on other fresh copies, the collector held off), on
        comb_chain's three parts (61 walks), sibling_merge's routes in the child and a route
        crossing a junction; each charge the cold price."""

        def crossing():
            ref = role('splits45', 'reference')
            merge = ref.arms['right'].segments[5]
            merge.parents = (4, 2)
            merge.partition = list(reversed(merge.partition))
            return pairs(stitch(role('splits45', 'parent'), role('splits45', 'child')), ref)

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
                'segments': lambda u, **k: u.segments('right', **k),
                'carriers': lambda u, **k: u.carriers('right', 1, **k),
                'statements': lambda u, **k: u.statements('right', **k),
                'summary': lambda u, **k: u.summary(**k),
                'to_fasta': lambda u, **k: u.to_fasta(stubs=True, **k),
                'to_gfa': lambda u, **k: u.to_gfa(labels='full', **k),
                'spell': lambda u, **k: u.spell('right', 1, **k)}
        for name, make in (('comb_chain', lambda: pairs(stitched('comb_chain'),
                                                         role('comb_chain', 'reference'))),
                           ('sibling_merge', lambda: pairs(stitched('sibling_merge'),
                                                            role('sibling_merge', 'reference'))),
                           ('crossing', crossing)):
            with self.subTest(name):
                b = LocalBudget()
                union_pairs(make(), budget=b)
                parts, other = make(), make()
                got = peak(lambda: union_pairs(parts), lambda: union_pairs(other))
                self.assertLessEqual(got, b.peak_bytes, (name, 'union_pairs'))
                for op, fn in ops_.items():
                    b = LocalBudget()
                    v = union_pairs(make())
                    fn(v, budget=b)
                    again = LocalBudget()
                    fn(v, budget=again)
                    self.assertEqual((again.used_work, again.peak_bytes),
                                     (b.used_work, b.peak_bytes), (name, op))
                    u, w = union_pairs(make()), union_pairs(make())
                    got = peak(lambda: fn(u), lambda: fn(w))
                    self.assertLessEqual(got, b.peak_bytes, (name, op))

    def test_memory_bytes_counts_the_stitched_part(self):
        from metagraph.traverse import composite as K
        sg = stitched('comb_chain')
        u = union_pairs(pairs(sg, None))
        self.assertGreaterEqual(u.memory_bytes(),
                                K.models_memory_bytes([p.graphlet for p in sg.parts]))


if __name__ == '__main__':
    unittest.main()
