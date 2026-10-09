"""The library speed-ups found by profiling (P1), each checked to give exactly the answers
of the code it replaced:

  * merge partitions read as sets (derive.partition_sets): membership in an array('I') of
    up to every label id was a linear scan per merge passed (evidence, routes, the
    comparison's route reconstruction); Segment.partition keeps its public type;
  * runs indexed by label (derive.runs_by_label) for routes() and label_walks(), which
    scanned every run of the arm per label;
  * a name/ref -> label index for label() (it scanned every label per lookup);
  * chains and spellings of many walks built together with their common prefixes shared
    (derive.walk_batch: walks(), to_fasta(), to_json(), to_gfa's paths, the comparison's
    keyers), and to_gfa's k-1 bases of context read back only that far;
  * walks(top=N) builds the claims of the N walks it returns, not of the whole arm.

The gate: the digests of every public operation and of the MCP tools over the committed
real and review3 retrievals equal those recorded with the library before these changes
(traverse_golden.py, data/traverse/golden/committed.json.gz).
"""

import glob
import os
import random
import sys
import time
import tracemalloc
import unittest
from array import array

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402
import traverse_golden as TG  # noqa: E402
from test_traverse_scale import comb_annotate  # noqa: E402

import metagraph.traverse as MT  # noqa: E402
from metagraph.traverse import UnknownLabel, AmbiguousLabel, derive, export, ops, parse  # noqa: E402

GOLDEN = os.path.join(TG.GOLDEN_DIR, 'committed.json.gz')


def _all_graphlets():
    """Every committed retrieval with its envelope (documents, review3, real)."""
    files = sorted(glob.glob(os.path.join(T.DOCS, '**', '*.graphlet.json'), recursive=True))
    out = [(k, MT.from_response(r, resp))
           for k, r, resp in TG.items_of(files + TG.golden_files())]
    return out


GRAPHLETS = _all_graphlets()


class TestGoldenUnchanged(unittest.TestCase):
    """Every public operation, compare() and the tools over the committed real and review3
    retrievals: the same digests as the library before the speed-ups and stage L."""

    @classmethod
    def setUpClass(cls):
        cls.want = TG.read_json(GOLDEN)
        items = TG.items_of(TG.golden_files(), base=TG.DATA)
        cls.got = TG.all_digests(MT, items)

    def test_the_same_keys(self):
        self.assertEqual(sorted(self.want), sorted(self.got))

    def test_every_digest_unchanged(self):
        bad = sorted(k for k in self.want if self.got.get(k) != self.want[k])
        self.assertEqual([], bad[:20], '%d of %d digests differ' % (len(bad), len(self.want)))


class TestPartitionSets(unittest.TestCase):
    def test_sets_equal_the_arrays_and_the_field_keeps_its_type(self):
        n = 0
        for key, g in GRAPHLETS:
            for a in g.arms.values():
                ps = derive.partition_sets(a)
                merges = [s for s in a.segments if len(s.parents) > 1]
                self.assertEqual(sorted(s.id for s in merges), sorted(ps), key)
                for s in merges:
                    n += 1
                    self.assertEqual(tuple(frozenset(p) for p in s.partition), ps[s.id])
                    for p in s.partition:
                        self.assertIsInstance(p, array)       # the public field: unchanged
        self.assertGreater(n, 10)

    def test_equal_arrays_share_one_set(self):
        g = T.graphlet('merge')
        a = g.arms['right']
        ps = derive.partition_sets(a)
        by_id = {}
        for s in a.segments:
            for p, fs in zip(s.partition, ps.get(s.id, ())):
                by_id.setdefault(id(p), set()).add(id(fs))
        self.assertTrue(all(len(v) == 1 for v in by_id.values()))

    def test_a_cleared_cache_follows_an_edited_partition(self):
        g = T.graphlet('merge')
        a = g.arms['right']
        m = next(s for s in a.segments if len(s.parents) > 1 and len(s.partition[0]))
        before = derive.partition_sets(a)[m.id][0]
        lab = m.partition[0][0]
        self.assertIn(lab, before)
        del m.partition[0][0]
        a.cache.clear()
        self.assertNotIn(lab, derive.partition_sets(a)[m.id][0])


def _routes_reference(g, sel, arm, spell=False):
    """routes() as it was: every run of the arm scanned, partitions read as arrays."""
    a = g.arm(arm)
    lab = ops.label(g, sel)
    out = []
    if g.mode != 'constrain':
        return ops.routes(g, sel, arm, spell)
    segs = a.segments
    for run in a.runs:
        if run.label != lab.id:
            continue
        route = [run.segment]
        s = run.segment
        while segs[s].parents:
            seg = segs[s]
            nxt = seg.parents[0]
            if len(seg.parents) > 1:
                at = derive.lineage_label_at(a, run, seg.from_bp)
                for j, part in enumerate(seg.partition):
                    if at in part:
                        nxt = seg.parents[j]
                        break
            route.append(nxt)
            s = nxt
        route.reverse()
        if spell:
            bases = ''.join(segs[x].walk or '' for x in route)[:run.to_bp]
            if segs[route[0]].walk is None:
                bases = None
            elif a.side == 'left':
                bases = bases[::-1]
            out.append((route, bases))
        else:
            out.append(route)
    return out


class TestRunsByLabel(unittest.TestCase):
    def test_the_index_is_the_scan(self):
        for key, g in GRAPHLETS:
            for a in g.arms.values():
                idx = derive.runs_by_label(a)
                for l in range(len(g.labels)):
                    self.assertEqual([r.id for r in a.runs if r.label == l],
                                     idx.get(l, []), key)

    def test_routes_and_label_walks_as_before(self):
        n = 0
        for key, g in GRAPHLETS:
            if g.mode != 'constrain':
                continue
            for side in g.arms:
                for l in range(min(len(g.labels), 12)):
                    for spell in (False, True):
                        self.assertEqual(_routes_reference(g, {'id': l}, side, spell),
                                         g.routes({'id': l}, side, spell=spell), key)
                    n += 1
        self.assertGreater(n, 50)


class TestLabelIndex(unittest.TestCase):
    def reference(self, g, sel):
        if isinstance(sel, str):
            return [l for l in g.labels if l.ref == sel or l.name == sel]
        (k, v), = sel.items()
        return [l for l in g.labels if getattr(l, k) == v]

    def check(self, g, sel):
        want = self.reference(g, sel)
        if len(want) == 1:
            self.assertIs(want[0], g.label(sel))
        elif not want:
            with self.assertRaises(UnknownLabel):
                g.label(sel)
        else:
            with self.assertRaises(AmbiguousLabel) as cm:
                g.label(sel)
            self.assertNotIsInstance(cm.exception, UnknownLabel)
            for l in want:
                self.assertIn(l.ref, str(cm.exception))

    def test_a_name_like_a_ref_a_shared_name_and_non_strings(self):
        g = T.body('limits')
        g.labels[2].name = 'c:0'
        for sel in ('c:0', {'name': 'c:0'}, {'ref': 'c:0'}, 'nothing'):
            self.check(g, sel)
        g = T.body('same_name')
        for sel in ('ACC1', {'name': 'ACC1'}):
            self.check(g, sel)
        for bad in ({'name': 5}, {'ref': None}, {'name': ['ACC1']}):
            with self.assertRaises(UnknownLabel):
                g.label(bad)


class TestWalkBatch(unittest.TestCase):
    def test_the_kept_prefixes_are_bounded_by_the_output(self):
        # a comb: every walk shares the spine; spelling them one by one walks the spine
        # per walk, keeping a prefix per spine segment would be quadratic in memory
        g = parse(comb_annotate(3000))
        a = g.arms['right']
        leaves = derive.leaves(a)
        out_bytes = sum(a.segments[x].end_bp for x in leaves)
        tracemalloc.start()
        t0 = time.perf_counter()
        got = derive.walk_batch(a, leaves)
        dt = time.perf_counter() - t0
        peak = tracemalloc.get_traced_memory()[1]
        tracemalloc.stop()
        self.assertEqual(len(leaves), len(got))
        # the spellings themselves, their chains (8 B per id) and the kept prefixes
        self.assertLess(peak, 4 * out_bytes + 16 * sum(len(c) for c, _ in got.values()))
        self.assertLess(dt, 10.0)
        one = derive.walk_batch(a, [leaves[-1]])
        self.assertEqual(derive.walk_bases(a, leaves[-1]), one[leaves[-1]][1])


def _walks_reference_claims(g, a, ranked_leaves):
    by_leaf = {}
    leaves_ = set(ranked_leaves)
    for c in ops.claims(g, a, strict=False):
        s = a.segments[c.segment]
        if c.segment in leaves_ and c.to_bp == s.end_bp:
            by_leaf.setdefault(c.segment, []).append(c)
    return by_leaf


class TestClaimsOfTheReturnedWalksOnly(unittest.TestCase):
    def test_the_claims_of_each_walk_as_before(self):
        n = 0
        for key, g in GRAPHLETS:
            for side, a in g.arms.items():
                for top in (1, 3, None):
                    ws = g.walks(side, top=top)
                    ref = _walks_reference_claims(g, a, [w.leaf for w in ws])
                    for w in ws:
                        n += 1
                        self.assertEqual([c.as_dict() for c in ref.get(w.leaf, [])],
                                         [c.as_dict() for c in w.claims], key)
        self.assertGreater(n, 100)

    def test_a_walk_asked_twice_gets_its_own_lists(self):
        g = T.graphlet('fork')
        w1, w2 = ops.walks_at(g, 'right', [0, 0, 1][:2] + [0], with_claims=True)[:2]
        self.assertEqual(w1.segments, w2.segments)
        self.assertIsNot(w1.segments, w2.segments)


def _gfa_reference_context(g, arm, s, k1):
    seed = g.seed.sequence
    before = derive.walk_bases(arm, s.parents[0]) if s.parents else ''
    outward = seed if arm.side == 'right' else seed[::-1]
    return (outward + before)[-k1:] if k1 > 0 else ''


class TestGfaContext(unittest.TestCase):
    def test_the_context_read_back_only_k_minus_1_bases(self):
        n = 0
        for key, g in GRAPHLETS:
            for a in g.arms.values():
                if not a.segments or a.segments[0].walk is None:
                    continue
                seed = g.seed.sequence
                outward = seed if a.side == 'right' else seed[::-1]
                for k1 in (0, 1, 4, g.k - 1, 200):
                    for s in a.segments:
                        before = export._chain_tail(a.segments, s.parents[0], k1) \
                            if s.parents and k1 > 0 else ''
                        got = (outward + before)[-k1:] if k1 > 0 else ''
                        self.assertEqual(_gfa_reference_context(g, a, s, k1), got)
                        n += 1
        self.assertGreater(n, 1000)


if __name__ == '__main__':
    unittest.main()
