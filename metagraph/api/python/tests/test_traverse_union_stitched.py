"""A StitchedGraphlet as a union part (docs/DESIGN-graphlet-merge.md §2.3, train 7): the
stitching triples of documents/stitch/ (scripts/traversal/graphlet_fixtures.py
STITCH_FIXTURES), the stitched graphlet as one pair's part beside its reference as another's.
Where no junction is read the union answers as with the root alone; past a junction nothing
reads complete or exact."""

import copy
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


# ------------------------------------------------------------------ the review of 7cbee3b03

def _edited_splits45(mirrored=False):
    """splits45 with the reference's first parent at its merge the stub (a union walk that
    spells the stub through the junction into the child), on the right arm or mirrored onto
    the left (each body's arm renamed, the child seeded with the mirrored continuation)."""
    p, c, ref = (role('splits45', r) for r in ('parent', 'child', 'reference'))
    merge = ref.arms['right'].segments[5]
    merge.parents = (4, 2)
    merge.partition = list(reversed(merge.partition))
    if mirrored:
        def mirror(g):
            a = g.arms.pop('right')
            a.side = 'left'
            g.arms['left'] = a
            a.cache.clear()
            return g
        p, c, ref = mirror(p), mirror(c), mirror(ref)
        c.seed.sequence = derive.continuation_sequence(p, p.arms['left'], 5)
    return stitch(p, c), p, ref


def _zero_leaf_root(first_parent_of_m=1, z_bases=''):
    """A root on limits' model: ACG -> {C, G} -> M 'AA' -> {C, G} -> Z, Z a merge leaf of
    |z_bases| (zero-length: a merge at the radius, max_steps, with a continuation)."""
    from array import array
    from metagraph.traverse.model import Run, Segment

    def seg(i, parents, f, w, entry, end, children, partition=(), leaf=None, split=None):
        return Segment(i, tuple(parents), f, len(w), w, w[0] if w else None,
                       array('I', entry), len(entry), array('I', end),
                       [array('I', x) for x in partition], split, [], [], leaf, list(children),
                       0)
    g = copy.deepcopy(T.graphlet('limits'))
    a = g.arms['right']
    old_leaf = a.segments[3].leaf
    mp = (1, 2) if first_parent_of_m == 1 else (2, 1)
    segs = [seg(0, (), 0, 'ACG', [0, 1, 2], [0, 1, 2], [1, 2], split=False),
            seg(1, (0,), 3, 'C', [0, 1], [0], [3]),
            seg(2, (0,), 3, 'G', [2], [2], [3]),
            seg(3, mp, 4, 'AA', [0, 2], [0, 2], [4, 5],
                partition=[[0], [2]] if mp == (1, 2) else [[2], [0]], split=False),
            seg(4, (3,), 6, 'C', [0, 2], [0, 2], [6]),
            seg(5, (3,), 6, 'G', [0, 2], [0, 2], [6]),
            seg(6, (4, 5), 7, z_bases, [0, 2], [0, 2], [], partition=[[0, 2], [0, 2]],
                leaf=copy.deepcopy(old_leaf))]
    for s in segs:
        s.depth = 0 if not s.parents else segs[s.parents[0]].depth + 1
    a.segments = segs
    n = 7 + len(z_bases)
    a.runs = [Run(0, 6, 0, 0, n, 'S', 'max_steps', None, False, False, 0, None, None, None,
                  0, 0, 0.0, None),
              Run(1, 1, 1, 0, 4, 'B', 'loss_budget', None, False, False, 0, None, None, None,
                  2, 0, 0.0, 2.0),
              Run(2, 6, 2, 0, n, 'S', 'max_steps', None, False, False, 7, None, None, None,
                  0, 0, 0.0, None)]
    a.cap_trigger.segment = 6
    a.complete_to_bp = 7
    a.cache.clear()
    return g


class TestReviewFixes(unittest.TestCase):
    """The review of the assembly at 7cbee3b03 (REVIEW-opus-py-7cbee3b03.md): M1, minors 1-4,
    the nits, and the claims it found untested."""

    def test_a_child_end_the_graph_holds_is_never_graph(self):
        """M1: annotate_loop's child walk (1, 0) ends with path reason D; the reference goes
        on with A. With an unfollowed hairpin for A at that end (annotate records a skipped
        hairpin so) the graph holds A: uncertified. Without it, a real dead end: graph."""
        from metagraph.traverse.model import Event
        name, side = 'annotate_loop', 'right'
        p, c, ref = (role(name, r) for r in ('parent', 'child', 'reference'))
        sg = stitch(p, c)
        sw = next(w for w in sg.walks(side, by='id') if len(w.key) == 2)
        rp = next(q for q in derive.paths(ref.arms[side]) if ref.spell(side, q.id) == sw.sequence)
        ref.arms[side].segments[rp.leaf].walk += 'A'
        ref.arms[side].segments[rp.leaf].length_bp += 1
        ref.arms[side].cache.clear()
        got = {}
        for hairpin in (True, False):
            child = copy.deepcopy(c)
            a = child.arms[side]
            s = a.segments[derive.paths(a)[sw.key[1]].leaf]
            s.leaf.path_reason = 'D'
            if hairpin:
                s.events.append(Event(at_bp=s.end_bp, type='hairpin', char='A', total=0,
                                      followed=False))
            a.cache.clear()
            u = union_pairs(pairs(stitch(p, child), ref))
            (w,) = [w for w in u.walks(side, stubs=True)
                    if u.spell(side, w.id) == sw.sequence + 'A']
            sp = next(x for x in u.carriers(side, w.id) if x.pair == 'A')
            self.assertEqual((sp.how, sp.to_bp, [j.at_bp for j in sp.junctions]),
                             ('ends', sw.length_bp, [sg.junctions[0].at_bp]))
            got[hairpin] = sp.absence
        self.assertEqual(got, {True: 'uncertified', False: 'graph'})

    def test_a_route_crossing_at_a_zero_length_leaf(self):
        """Minor 1: the walk a child continues ends at a zero-length merge leaf Z that the
        route reaches through the merge M: it crosses into the child as at a leaf with bases."""
        p = _zero_leaf_root()
        child = T.graphlet('limits')
        child = copy.deepcopy(child)
        for s, w in zip(child.arms['right'].segments, ['ACG', 'CAAA', 'GAAA', '']):
            s.walk = w
        child.seed.sequence = 'ACAAG'
        child.seed.num_seed_labels = 1
        child.arms['right'].cache.clear()
        sg = stitch(p, child)
        self.assertEqual([j.at_bp for j in sg.junctions], [7])
        ref = _zero_leaf_root(first_parent_of_m=2, z_bases='ACG')
        got = {}
        for label, a in (('stitched', sg), ('root', p)):
            u = union_pairs(pairs(a, ref))
            (w,) = [w for w in u.walks('right', stubs=True)
                    if u.spell('right', w.id) == 'ACGGAACACG']
            sp = next(x for x in u.carriers('right', w.id) if x.pair == 'A')
            r = sp.route
            got[label] = (sp.how, r.how, r.to_bp, r.goes_on, r.complete, r.certified_to_bp,
                          [j.at_bp for j in sp.junctions], w.held_in)
        self.assertEqual(got['stitched'], ('merges', 'held', 10, True, False, 7, [7],
                                           ('A', 'B')))
        self.assertEqual(got['root'], ('merges', 'ends', 7, False, False, 7, [], ('B',)))

    def test_walk_causes_reach_every_statement(self):
        """Minor 2: a whole span past a junction states its stitched walk's walk-scope
        causes in its junctions; summary() states each junction's walk_causes as the stitched
        graphlet's summary does; FASTA, GFA and the coverage note carry them (loop_branch:
        candidate_changed)."""
        from metagraph.traverse._stitch import WALK_CAUSES
        for name in NAMES:
            with self.subTest(name):
                sg = stitched(name)
                u = union_pairs(pairs(sg, None))
                for side in sg.arms():
                    sws = {w.key: w for w in sg.walks(side, by='id')}
                    for w in u.walks(side):
                        (sp,) = u.carriers(side, w.id)
                        if not sp.junctions:
                            continue
                        mine = {c for j in sp.junctions for c in j.causes
                                if c in WALK_CAUSES and c != 'exact_gated'}
                        theirs = {c[0] for c in sws[sp.walk].causes if c[0] != 'exact_gated'}
                        self.assertEqual(mine, theirs, sp.walk)
                want = {(j['at_bp'], j['walk']): j['walk_causes']
                        for j in sg.summary(max_bytes=10 ** 6)['junctions']}
                from metagraph.traverse._stitch import walk_name
                got = {(j['at_bp'], walk_name(x.walk)): j['walk_causes'] for j, x in zip(
                    u.summary(max_bytes=10 ** 6)['pairs']['A']['junctions'], sg.junctions)}
                self.assertEqual(got, want)
        u = union_pairs(pairs(stitched('loop_branch'), role('loop_branch', 'reference')))
        self.assertIn(' junctions=A@115:qualified:exact_gated|candidate_changed ', u.to_fasta())
        self.assertIn('(qualified: exact_gated; walk causes below: candidate_changed)',
                      u.to_gfa())
        self.assertIn('walk causes below: candidate_changed',
                      u.coverage('right').per_pair['A']['note'])

    def test_the_form_kept_does_not_depend_on_the_order(self):
        """Minor 3: two forms of one stitched graphlet (equal digests) whose children differ
        in seed_id only: the union keeps the same form in either order (§2.10)."""
        p, c = role('keep_forbid', 'parent'), role('keep_forbid', 'child')
        c2 = copy.deepcopy(c)
        c2.seed_summary = copy.deepcopy(c2.seed_summary)
        c2.seed_summary['seed']['seed_id'] = 'zzz-other-id'
        s1, s2 = stitch(p, c), stitch(p, c2)
        self.assertEqual(s1.digest, s2.digest)
        ref = role('keep_forbid', 'reference')
        outs = []
        for order in ((s1, s2), (s2, s1)):
            u = union_pairs([{'pair': 'A', 'graphlet': x} for x in order]
                            + [{'pair': 'B', 'graphlet': ref}])
            self.assertEqual(len(u.parts), 2)
            outs.append((u.digest, answers(u), repr(u.provenance())))
        self.assertEqual(outs[0], outs[1])

    def test_an_arm_without_a_junction_is_judged_by_the_root(self):
        """Minor 4: fork_right's junction is on the right arm; the root states index_fp and
        the child none (allow_unverified_index). evidence('left') reads no junction: exact
        as the root alone."""
        p, c = role('fork_right', 'parent'), role('fork_right', 'child')
        p.index_fp = 'f' * 64
        sg = stitch(p, c, allow_unverified_index=True)
        self.assertFalse(sg.identity['verified'])
        alone = union_pairs([{'pair': 'A', 'graphlet': p}])
        u = union_pairs([{'pair': 'A', 'graphlet': sg}])
        self.assertTrue(alone.evidence('left')['exact'])
        self.assertEqual(u.evidence('left')['exact'], alone.evidence('left')['exact'])
        self.assertEqual(u.evidence('left')['exact'], sg.evidence('left')['exact'])
        self.assertFalse(u.evidence('right')['exact'])

    def test_traverse_pairs_reports_a_body_that_does_not_parse(self):
        """Nit: a pair's body that does not parse is that pair's {pair, error}."""
        from metagraph.traverse._codec import GraphletFormatError
        from tests.test_traverse_union import TestTraversePairs, _FanOutClient
        bad = GraphletFormatError(3, 'not MGT')
        cl = _FanOutClient.make(TestTraversePairs.CAPS, {'g1': 'answer', 'g2': bad, 'g3': 'x'})
        r = cl.traverse_pairs('ACGT')
        self.assertEqual(r.results, ['answer', {'pair': {'graph': 'g2'}, 'error': bad}, 'x'])

    def test_a_junction_certified_below_its_depth(self):
        """A root certified below the junction (a lowered complete_to_bp on keep_forbid, and
        beam_sibling, whose walk stopped at a cap at 1 bp): every span past the junction is
        certified to the root's certificate, the stitched walk's certified_to_bp."""
        p, c = role('keep_forbid', 'parent'), role('keep_forbid', 'child')
        p.arms['right'].complete_to_bp = 12
        cases = [(stitch(p, c), role('keep_forbid', 'reference'), 12),
                 (stitched('beam_sibling'), role('beam_sibling', 'reference'), 1)]
        for sg, ref, cert in cases:
            u = union_pairs(pairs(sg, ref))
            sws = {w.key: w for w in sg.walks('right', by='id')}
            seen = 0
            for w in u.walks('right', stubs=True):
                for sp in u.carriers('right', w.id):
                    if sp.pair == 'A' and sp.junctions:
                        seen += 1
                        self.assertEqual(sp.certified_to_bp, cert)
                        if sp.how == 'whole':
                            self.assertEqual(sws[sp.walk].certified_to_bp, cert)
            self.assertGreater(seen, 0)

    def test_a_walk_ending_at_a_graft_node(self):
        """A union walk that ends where the child continues (keep_forbid's parent beside the
        stitched graphlet): the stitched part extends there into the child, reading the
        junction, with the parent's labels (no child base is read)."""
        p, c = role('keep_forbid', 'parent'), role('keep_forbid', 'child')
        sg = stitch(p, c)
        (j,) = sg.junctions
        u = union_pairs([{'pair': 'A', 'graphlet': sg},
                         {'pair': 'B', 'graphlet': role('keep_forbid', 'parent')}])
        (w,) = [w for w in u.walks('right') if w.length_bp == j.at_bp
                and next(x for x in u.carriers('right', w.id) if x.pair == 'A').junctions]
        a, b = u.carriers('right', w.id)
        self.assertEqual((a.how, a.to_bp, a.walk, a.complete, a.certified_to_bp),
                         ('extends', j.at_bp, j.walk + (0,), False, j.at_bp))
        self.assertEqual((b.how, b.walk), ('whole', j.walk[0]))
        self.assertEqual([(x.label.ref, x.from_bp) for x in a.labels],
                         [(x.label.ref, x.from_bp) for x in b.labels])

    def test_labels_before_a_junction_are_the_parents(self):
        """carriers(at_bp=) at a depth before the junction on a span past it: the parent's
        labels along the continued walk, as the root alone gives them there."""
        for name in ('keep_forbid', 'chain3', 'sibling_merge'):
            with self.subTest(name):
                sg = stitched(name)
                p = role(name, 'parent')
                u, up = union_pairs(pairs(sg, None)), union_pairs(pairs(p, None))
                plain = walks_by_spelling(up, 'right')
                checked = 0
                for w in u.walks('right'):
                    (sp,) = u.carriers('right', w.id)
                    if not sp.junctions:
                        continue
                    at = sp.junctions[0].at_bp
                    spelled = u.spell('right', w.id, with_seed=True)
                    pw = plain.get(('walk', spelled[:len(spelled) - w.length_bp + at]))
                    if pw is None or at < 2:
                        continue
                    (mine,) = u.carriers('right', w.id, at_bp=at - 1)
                    (root,) = up.carriers('right', pw.id, at_bp=at - 1)
                    self.assertEqual(mine.labels, root.labels)
                    checked += 1
                self.assertGreater(checked, 0)

    def test_walks_by_labels_never_select_past_a_junction(self):
        """walks(labels=) keeps a walk a label carries whole (from 0 bp); a stitched part's
        labels past a junction start there, so no walk with bases past a junction is
        selected by them (graphlets.rst says so). The stitched graphlet itself selects
        keep_forbid's continued walk by acc1; the union does not."""
        for name in NAMES:
            sg = stitched(name)
            u = union_pairs(pairs(sg, None))
            for side in u.arms():
                for lab in sg.labels:
                    for w in u.walks(side, labels=[{'pair': 'A', 'ref': lab.ref}]):
                        (sp,) = u.carriers(side, w.id)
                        self.assertTrue(not sp.junctions
                                        or sp.junctions[0].at_bp >= w.length_bp,
                                        (name, lab.name, w.id))
        sg = stitched('keep_forbid')
        acc1 = next(x for x in sg.labels if x.name == 'acc1')
        self.assertEqual([w.key for w in sg.walks('right', by='id', labels=['acc1'])],
                         [(1, 0)])
        u = union_pairs(pairs(sg, None))
        self.assertEqual(u.walks('right', labels=[{'pair': 'A', 'ref': acc1.ref}]), [])

    def test_a_left_arm_crossing(self):
        """splits45's crossing mirrored onto the left arm: the same facts as on the right."""
        got = []
        for mirrored, side in ((False, 'right'), (True, 'left')):
            sg, _, ref = _edited_splits45(mirrored)
            u = union_pairs(pairs(sg, ref))
            (w,) = [w for w in u.walks(side) if w.length_bp == 50]
            sp = next(x for x in u.carriers(side, w.id) if x.pair == 'A')
            r = sp.route
            got.append((sp.how, sp.to_bp, r.how, r.to_bp, r.goes_on, r.complete,
                        r.certified_to_bp, [j.at_bp for j in sp.junctions]))
        self.assertEqual(got[0], ('merges', 31, 'held', 50, True, False, 45, [45]))
        self.assertEqual(got[0], got[1])


if __name__ == '__main__':
    unittest.main()
