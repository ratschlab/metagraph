"""Regression tests for the product bugs the real-index suite found (api/python/tests/real),
each reproduced offline: on a hand-made document, a fixture of data/traverse/documents,
or a stand-in backend. One class per bug, named by its id in the triage."""

import copy
import gc
import json
import os
import re
import sys
import tempfile
import time
import tracemalloc
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402

from metagraph.traverse import (  # noqa: E402
    GraphletFormatError, GraphletStore, TraverseClient, derive, from_response, is_canonical,
    ops, parse,
)
from metagraph.traverse import _codec  # noqa: E402
from metagraph.traverse.mcp_tools import (  # noqa: E402
    MIN_MAX_BYTES, GraphletTools,
)


def size(obj):
    return len(json.dumps(obj, separators=(',', ':'), ensure_ascii=False).encode('utf-8'))


# ------------------------------------------------------------------ hand-made documents

# annotate + on_reconverge merge, right arm, k 5. The root [0, 10) records a.fa, b.fa and
# c.fa; it splits into segment 1 ('A', b.fa only) and segment 2 ('C', a.fa and b.fa); the
# two merge into segment 3 [20, 30) (parents 1, 2: segment 1 is the FIRST parent), which
# records a.fa and b.fa and is the only leaf. The displayed walk is 0-1-3; segment 2 is
# on no displayed walk. a.fa: on the displayed walk [0, 10) only; along the route 0-2-3
# on every node to 30 (its direct_bp). c.fa: on the root only (no child records it).
ANNOT_MERGE = '''\
H mgt 1 5 basic $ACGT a k m 8 0 0 walk fixtures * *
S * 10 6 0 TTGCAACGTA
L c 0 * 0 a.fa
L c 1 * 0 b.fa
L c 2 * 0 c.fa
O c c c i
A r c 30 u 0 0 1 3 0 . * 0 * 4 0 1 1 1 40
B 0 2 3 5 1 40 1 0 1 1 1 0 0 .
G * 0 10 0-2 * * * 0 * ACGTACGTAC
P 0 10 3 !
G 0 10 10 1 * * * * * ATTTTTTTTT
P 10 20 1 1
G 0 10 10 !2 * * * * * CTTTTTTTTT
P 10 20 2 !
G 1,2 20 10 ! * * * * * GGGGGGGGGG
P 20 30 2 !
T D .
Z 18
'''


def annot_merge():
    return parse(ANNOT_MERGE)


def without_bases(g):
    """The same retrieval as output.sequences false would have written it: no G bases,
    a single-parent segment's first base, the C sequences explicit (canonical text)."""
    g = copy.deepcopy(g)
    for a in g.arms.values():
        a.cache.clear()
        for s in a.segments:
            if s.leaf is not None and s.leaf.continuation is not None:
                c = s.leaf.continuation
                if c.sequence is None and c.n:
                    c.sequence = derive.continuation_sequence(g, a, s.id)
        for s in a.segments:
            if len(s.parents) == 1 and s.walk:
                s.first_base = s.walk[0]
            s.walk = None
    g.cache.clear()
    text = g.dump(envelope=False)
    out = parse(text)
    out.envelope, out.seed_summary = g.envelope, g.seed_summary
    return out


class TestAnnotateMergeDocument(unittest.TestCase):
    def test_the_document_is_valid_and_canonical(self):
        self.assertTrue(is_canonical(ANNOT_MERGE))
        self.assertTrue(is_canonical(ANNOT_KEEP))
        g = annot_merge()
        a = g.arms['right']
        self.assertEqual([0, 1, 3], derive.paths(a)[0].segments)
        self.assertIsNone(ops._first_paths(a)[2])          # segment 2: no displayed walk
        self.assertEqual(30, g.label_summary()['c:0']['right']['direct_bp'])


# ------------------------------------------------------------------ BUG-ANNOT-CLAIMS

class TestAnnotateClaimsAreMaximalOnDisplayedWalks(unittest.TestCase):
    """A drop at a child's first base is not an end when another child carries the label
    on. Superseded in part by the GPT review's finding 3 (test_traverse_review): the
    claims are now the maximal ends of the §5.1 union rule, so a sibling on no displayed
    walk that carries the label to a merge makes ONE claim through the merge (displayed
    from the merge on), not a claim at the drop."""

    def test_a_carrier_on_no_displayed_walk_does_not_suppress_the_claim(self):
        g = annot_merge()
        got = {(c.label.name, c.from_bp, c.to_bp, c.segment, c.path_id, c.kind,
                c.evidence_from) for c in g.claims('right')}
        self.assertEqual({('a.fa', 0, 30, 3, 0, 'end', 20), ('b.fa', 0, 30, 3, 0, 'end', 0),
                          ('c.fa', 0, 10, 0, 0, 'end', 0)}, got)

    def test_a_displayed_carrier_still_suppresses_it(self):
        # no merge: segment 2 runs on to 30 (a leaf, displayed) and carries a.fa on, so
        # the drop at segment 1's first base is not maximal; segment 2's end is
        g = parse(ANNOT_KEEP)
        a = g.arms['right']
        self.assertEqual([[0, 1], [0, 2]], [p.segments for p in derive.paths(a)])
        claims = {(c.label.name, c.to_bp, c.segment) for c in g.claims('right')}
        self.assertNotIn(('a.fa', 10, 0), claims)
        self.assertIn(('a.fa', 30, 2), claims)


# the same locus without the merge (on_reconverge keep would keep both): segment 2 runs
# on to 30 as a leaf of its own
ANNOT_KEEP = '''\
H mgt 1 5 basic $ACGT a k k 8 0 0 walk fixtures * *
S * 10 6 0 TTGCAACGTA
L c 0 * 0 a.fa
L c 1 * 0 b.fa
L c 2 * 0 c.fa
O c c c i
A r c 30 p 0 0 1 3 0 . * 0 * 3 0 2 1 0 40
B 0 2 3 5 1 40 1 0 1 0 0 0 0 .
G * 0 10 0-2 * * * 0 * ACGTACGTAC
P 0 10 3 !
G 0 10 10 1 * * * * * ATTTTTTTTT
P 10 20 1 1
T D .
G 0 10 20 !2 * * * * * CTTTTTTTTTGGGGGGGGGG
P 10 30 2 !
T D .
Z 17
'''


# ------------------------------------------------------------------ BUG-ANNOT-ROUTES

class TestAnnotateRoutesReachDirectBp(unittest.TestCase):
    """routes() / label_walks() in annotate mode: one witness route per maximal end of the
    §5.1 union rule (the rule of label_summary's direct_bp), passing every merge through
    its first parent (in parents order) that records the label."""

    def test_the_direct_bp_route_goes_through_the_non_first_parent(self):
        g = annot_merge()
        self.assertEqual([([0, 2, 3], 'ACGTACGTACCTTTTTTTTTGGGGGGGGGG')],
                         g.routes('a.fa', 'right', spell=True))
        # b.fa can follow the displayed chain, so it does
        self.assertEqual([[0, 1, 3]], g.routes('b.fa', 'right'))
        self.assertEqual([([0], 'ACGTACGTAC')], g.routes('c.fa', 'right', spell=True))
        summ = g.label_summary()
        for lab in g.labels:
            d = summ[lab.ref]['right']['direct_bp']
            self.assertIn(d, [len(b) for _, b in g.routes(lab, 'right', spell=True)], lab)

    def test_label_walks_say_what_no_walk_displays(self):
        g = annot_merge()
        (lw,) = g.label_walks('a.fa')
        # the route meets the displayed chain at the merge (20): the walk through its
        # anchor (leaf 3, displayed 0-1-3) shows it from there on, not before
        self.assertEqual((0, 30, 20, 'dead_end', [3]),
                         (lw.from_bp, lw.to_bp, lw.evidence_from, lw.end, lw.leaves_below))
        self.assertEqual('ACGTACGTACCTTTTTTTTTGGGGGGGGGG', lw.sequence)
        (lb,) = g.label_walks('b.fa')
        self.assertEqual((30, 0, [3]), (lb.to_bp, lb.evidence_from, lb.leaves_below))
        (lc,) = g.label_walks('c.fa')
        self.assertEqual((10, 'label_lost', [3]), (lc.to_bp, lc.end, lc.leaves_below))

    def test_keep_mode_routes_are_the_displayed_stretches(self):
        # a tree: every route is a first-parent chain, and the routes are the claims
        for ga in (parse(T.compare_text('oracle_annotate')), parse(ANNOT_KEEP)):
            for side in ga.arms:
                for lab in ga.labels:
                    stretches = sorted((c.segment, c.to_bp) for c in
                                       ga.claims(side, labels=[lab], strict=False))
                    lws = ga.label_walks(lab, side)
                    rts = ga.routes(lab, side)
                    self.assertEqual(stretches, sorted((r[-1], lw.to_bp)
                                                       for r, lw in zip(rts, lws)))
                    for r, lw in zip(rts, lws):
                        self.assertEqual(derive.chain(ga.arms[side], r[-1]), r)
                        self.assertEqual(0, lw.evidence_from)
                        self.assertTrue(lw.leaves_below)

    def test_the_tool_shows_the_route(self):
        with tempfile.TemporaryDirectory() as d:
            store = GraphletStore(d)
            h = store.put_graphlet(annot_merge())
            tools = GraphletTools(store)
            out = tools.graphlet_labels(h, name='a.fa')
            self.assertNotIn('error', out, out)
            self.assertEqual([{'arm': 'right', 'from_bp': 0, 'to_bp': 30, 'evidence_from': 20,
                               'end': 'dead_end', 'walks_below': 1, 'merged_into': None,
                               'route': [0, 2, 3]}], out['rows'])


# ------------------------------------------------------------------ SPLIT-REASON-EMPTY

class TestSplitReasonsNameATakenBranch(unittest.TestCase):
    def test_a_label_no_branch_took_is_absent_not_split(self):
        g = annot_merge()
        (ch, back), = [g.support_changes('right', 0)]
        self.assertEqual(10, ch.at_bp)
        why = {r['label']['name']: r for r in ch.reasons}
        # a.fa went on along segment 2 ('C'); c.fa is on no child's first node
        self.assertEqual(('split', ['C']), (why['a.fa']['why'], why['a.fa']['took']))
        self.assertEqual('absent', why['c.fa']['why'])
        self.assertNotIn('took', why['c.fa'])
        self.assertEqual((20, ['a.fa']), (back.at_bp, [l.name for l in back.added]))


# ------------------------------------------------------------------ SRV-ANNOT-CONT-LABELS

def _tail_labels(g, side, pid, nodes):
    """The labels recorded on every one of the walk's last |nodes| nodes (flank nodes
    only: seed nodes record nothing); None when the range holds no flank node."""
    a = g.arms[side]
    L = a.segments[derive.paths(a)[pid].leaf].end_bp
    lo = max(0, L - nodes)
    have = None
    for r in g.support_profile(side, pid):
        if r.to_bp <= lo:
            continue
        ids = {l.id for l in r.labels}
        have = ids if have is None else have & ids
    return None if have is None else sorted(have)


class TestAnnotateContinuationLabelsAreTheTailKmers(unittest.TestCase):
    """Server fix, seen through the CLI fixtures (regenerated by graphlet_fixtures.py
    --from-cli): an annotate continuation of n bases names the labels recorded on every
    node whose k-mer lies in its tail -- the walk's last n - k + 1 nodes -- not on the
    last n nodes, whose first k - 1 begin before the tail (§7.1: the labels covering that
    whole tail). The C++ test Walker.AnnotateContinuationLabelsCoverTheTailKmers pins the
    walker; this pins the documents."""

    def test_the_fixtures_name_the_labels_of_the_tail_kmers(self):
        checked = changed = 0
        for name in ('annotate', 'same_name', 'no_bins'):
            g = T.graphlet(name)
            for side, a in g.arms.items():
                for p in derive.paths(a):
                    if a.segments[p.leaf].leaf.continuation is None:
                        continue
                    cont = g.continuation(side, p.id)
                    n = len(cont.sequence)
                    got = [l.id for l in cont.labels]
                    want = _tail_labels(g, side, p.id, max(1, n - g.k + 1))
                    with self.subTest(fixture=name, arm=side, walk=p.id):
                        if want is not None:
                            self.assertEqual(want, got)
                            checked += 1
                        old = _tail_labels(g, side, p.id, n)
                        changed += old is not None and old != got
        self.assertGreaterEqual(checked, 2)
        # the old rule (the last n nodes) names fewer labels on both CLI fixtures
        self.assertEqual(2, changed)
        self.assertEqual(['a.fa', 'b.fa'],
                         [l.name for l in T.graphlet('annotate').continuation('right', 0).labels])
        self.assertEqual(['ACC2'],
                         [l.name for l in T.graphlet('same_name').continuation('right', 1).labels])


# ------------------------------------------------------------------ COMPARE-NO-BASES

class TestCompareWithoutBases(unittest.TestCase):
    """A retrieval with output.sequences false has no displayed bases to key claims and
    walks by: claims, walks and prefix_subset are 'unknown' with the reason naming
    output.sequences (the GPT review's finding 5; the first fix matched claims as a
    multiset, 'qualified'), never every claim one-sided with comparable True; labels
    mode needs no bases."""

    def pairs(self):
        for name in ('merge', 'switch_chain', 'annotate', 'fork'):
            g = T.graphlet(name)
            yield name, g, without_bases(g)
        # one index identity on both sides: comparable True where bases exist
        for name in ('oracle_exhaustive', 'oracle_constrain'):
            g = parse(T.compare_text(name))
            yield name, g, without_bases(g)

    def test_the_same_trie_without_bases(self):
        verified = 0
        for name, g, nb in self.pairs():
            self.assertTrue(all(s.walk is None for a in nb.arms.values() for s in a.segments))
            for x, y in ((g, nb), (nb, g), (nb, nb)):
                with self.subTest(fixture=name, a=x is g):
                    for mode in ('claims', 'walks', 'prefix_subset'):
                        cmp = x.compare(y, mode=mode)
                        self.assertEqual(('unknown', None, [], [], []),
                                         (cmp.comparable, cmp.equal, cmp.only_in_a,
                                          cmp.only_in_b, cmp.differ))
                        self.assertIn('output.sequences', cmp.reason)
                        self.assertFalse(any('no difference found' in n for n in cmp.notes))
                    # labels need no bases: equal wherever the identities make it comparable
                    cmp = x.compare(y, mode='labels')
                    self.assertEqual([], cmp.only_in_a)
                    if cmp.comparable is True:
                        self.assertTrue(cmp.equal)
                        verified += 1
        self.assertGreater(verified, 0)

    def test_with_bases_the_verdict_is_unchanged(self):
        g = parse(T.compare_text('oracle_exhaustive'))
        cmp = g.compare(g)
        self.assertEqual((True, True), (cmp.comparable, cmp.equal))

    def test_different_tries_are_not_compared_either(self):
        # no difference is reported from labels and intervals alone: claims of different
        # walks share them, so neither equality nor a difference is claimed
        a = without_bases(parse(T.compare_text('cli_bubbles_r100')))
        b = without_bases(parse(T.compare_text('cli_bubbles_r60')))
        for mode in ('claims', 'walks', 'prefix_subset'):
            cmp = a.compare(b, mode=mode)
            self.assertEqual(('unknown', None, [], [], []),
                             (cmp.comparable, cmp.equal, cmp.only_in_a, cmp.only_in_b,
                              cmp.differ))

    def test_the_tool_answers(self):
        with tempfile.TemporaryDirectory() as d:
            store = GraphletStore(d)
            g = T.graphlet('merge')
            ha, hb = store.put_graphlet(g), store.put_graphlet(without_bases(g))
            tools = GraphletTools(store)
            for mode in ('claims', 'walks', 'labels', 'prefix_subset'):
                out = tools.graphlet_compare(ha, hb, mode=mode)
                self.assertNotIn('error', out, (mode, out))
                self.assertEqual(0, out['total'], (mode, out))
                if mode != 'labels':
                    self.assertEqual(('unknown', None), (out['comparable'], out['equal']), out)
                    self.assertIn('output.sequences', out['reason'])


# ------------------------------------------------------------------ COMPARE-SELECTORS-A-ONLY

class TestCompareSelectorsResolveOnEitherSide(unittest.TestCase):
    def test_a_ref_only_b_knows(self):
        a = parse(T.compare_text('oracle_constrain'))      # acc1, acc3
        b = parse(T.compare_text('oracle_exhaustive'))     # acc1, acc2, acc3
        self.assertNotIn('h:0:2', {l.ref for l in a.labels})
        for sel in ({'ref': 'h:0:2'}, 'h:0:2', {'name': b.label({'ref': 'h:0:2'}).name}):
            ab = a.compare(b, labels=[sel])
            ba = b.compare(a, labels=[sel])
            self.assertEqual([], ab.only_in_a)
            self.assertTrue(ab.only_in_b)
            self.assertEqual({'h:0:2'}, {r['label']['ref'] for r in ab.only_in_b})
            self.assertEqual(ab.only_in_b, ba.only_in_a)

    def test_neither_side_knows_it(self):
        a = parse(T.compare_text('oracle_constrain'))
        b = parse(T.compare_text('oracle_exhaustive'))
        for sel in ({'ref': 'h:9:9'}, 'nope', {'name': 'nope'}):
            with self.assertRaises(ops.UnknownLabel):
                a.compare(b, labels=[sel])

    def test_an_id_is_a_ids(self):
        a = parse(T.compare_text('oracle_constrain'))
        b = parse(T.compare_text('oracle_exhaustive'))
        # id 1 is acc3 in a (acc2 in b): ids are per graphlet, resolved against a
        cmp = a.compare(b, labels=[{'id': 1}])
        refs = {r['label']['ref'] for r in cmp.only_in_a + cmp.only_in_b}
        self.assertLessEqual(refs, {a.labels[1].ref})
        with self.assertRaises(ops.UnknownLabel):
            a.compare(b, labels=[{'id': 2}])


# ------------------------------------------------------------------ B1-COMPARE-WALKS-COST

def _reference_profile(g, a, leaf):
    """The whole chain's displayed support profile root -> |leaf| (any segment): the
    definition support_profile() applies to a walk's leaf."""
    if g.mode == 'constrain':
        return [(f, t, set(labs)) for f, t, labs in ops._alive_profile(g, a, leaf)]
    return [(f, t, set(labs)) for f, t, labs, _ in ops._presence_profile(a, leaf)]


def _reference_cut_labels(g, a, leaf, m):
    """The pre-fix definition (the whole chain's support profile, cut at m)."""
    if m == 0:
        return {g.labels[l].ref for l in a.segments[0].entry}
    have = None
    for f, t, ids in _reference_profile(g, a, leaf):
        if f >= m:
            break
        have = ids if have is None else have & ids
    return {g.labels[l].ref for l in (have or ())}


def _reference_leaves(a, depth):
    """The walks of the DAG restricted to [0, depth) (third review, finding 6): every
    segment starting before the depth none of whose children does -- written out here
    independently of ops.restricted_leaves. At depth 0 every walk is cut to nothing."""
    segs = a.segments
    if depth <= 0:
        return [p.leaf for p in derive.paths(a)]
    return [s.id for s in segs if s.from_bp < depth
            and not any(segs[c].from_bp < depth for c in s.children)]


def _reference_keys(g, side, depth):
    a = g.arms[side]
    walks, pairs, supported = {}, set(), {}
    for leaf in _reference_leaves(a, depth):
        m = min(a.segments[leaf].end_bp, depth)
        w = derive.walk_bases(a, leaf)[:m]
        refs = _reference_cut_labels(g, a, leaf, m)
        walks.setdefault((side, w if side == 'right' else w[::-1]), set()).update(refs)
        for ref in refs:
            pairs.add((side, w, ref))
        have = {l: 0 for l in a.segments[0].entry}
        cur = None
        for f, t, ids in _reference_profile(g, a, leaf):
            if f >= m:
                break
            cur = ids if cur is None else cur & ids
            for l in cur:
                have[l] = min(t, m)
        for l, j in have.items():
            supported.setdefault((side, g.labels[l].ref), set()).add(w[:j])
    return walks, pairs, supported


class TestCompareKeysAreDepthBounded(unittest.TestCase):
    """compare(mode='walks'|'prefix_subset') keys each walk by its cut [0, m): computed
    once per segment holding base m - 1, reading the chain only up to m. The keys must
    equal the whole-chain definition, at every depth -- over the walks of the DAG
    restricted to the depth (third review, finding 6: a segment entering a merge beyond
    the depth through a non-first parent is a walk of its own there)."""

    def docs(self):
        for name in ('merge', 'switch_chain', 'reentry', 'fork', 'annotate', 'quorum',
                     'ambiguous_split', 'linear'):
            yield name, T.graphlet(name)
        for name in ('cli_bubbles_r100', 'cli_quorum_exhaustive', 'oracle_annotate'):
            yield name, parse(T.compare_text(name))
        yield 'annot_merge', annot_merge()

    def test_the_keys_equal_the_whole_chain_definition(self):
        n = 0
        for name, g in self.docs():
            for side, a in g.arms.items():
                top = max((p.length_bp for p in derive.paths(a)), default=0)
                for depth in sorted({0, 1, 2, 5, 9, 10, 11, 20, 36, 37, 62, top, top + 5}):
                    memo = {}
                    want_w, want_p, want_s = _reference_keys(g, side, depth)
                    with self.subTest(doc=name, arm=side, depth=depth):
                        self.assertEqual(want_w, ops._walk_keys(g, [side], depth, None, memo))
                        self.assertEqual(want_p, ops._route_pairs(g, [side], depth, memo))
                        self.assertEqual(want_s, ops._supported_prefixes(g, [side], depth,
                                                                         memo))
                    n += 1
        self.assertGreater(n, 50)

    def test_a_shallow_comparison_does_not_read_the_whole_trie(self):
        from test_traverse_scale import comb_annotate
        big, small = parse(comb_annotate(3000)), parse(comb_annotate(4))
        t0 = time.perf_counter()
        for mode in ('walks', 'prefix_subset'):
            cmp = big.compare(small, mode=mode)
            self.assertEqual(5, cmp.depth_used)
        dt = time.perf_counter() - t0
        # the whole-chain definition reads 3001 chains of up to 3001 segments here (~4.5 M
        # profile pieces, tens of seconds); the cut reads the root's few segments
        self.assertLess(dt, 2.0, dt)


# ------------------------------------------------------------------ BUG-1, FORMAT-ERROR-ECHO

class TestHugeTokens(unittest.TestCase):
    """A digit string longer than sys.get_int_max_str_digits() reached int(), whose
    ValueError escaped parse(); and every error echoed the whole token."""

    def mutate(self, text, pattern, repl):
        out, n = re.subn(pattern, repl, text, count=1, flags=re.M)
        self.assertEqual(1, n, pattern)
        line = text[:re.search(pattern, text, flags=re.M).start()].count('\n') + 1
        return out, line

    def test_huge_integers_are_format_errors_with_a_line(self):
        text = T.doc_text('merge')
        cases = [(r'^Z \d+$', 'Z ' + '9' * 5000),
                 (r'^(A r c )\d+', r'\g<1>' + '7' * 5000),
                 (r'^(G \* 0 )\d+', r'\g<1>' + '1' * 4400),
                 (r'^(R \d+ )\d+', r'\g<1>' + '3' * 6000)]
        for pattern, repl in cases:
            bad, line = self.mutate(text, pattern, repl)
            with self.subTest(pattern=pattern):
                with self.assertRaises(GraphletFormatError) as cm:
                    parse(bad)
                self.assertEqual(line, cm.exception.line_no)
                self.assertLess(len(str(cm.exception)), 200, str(cm.exception)[:300])

    def test_25_digits_still_does_not_fit_64_bits(self):
        for digits in (21, 25, 5000):
            bad, line = self.mutate(T.doc_text('merge'), r'^Z \d+$', 'Z ' + '9' * digits)
            with self.assertRaisesRegex(GraphletFormatError, 'does not fit 64 bits'):
                parse(bad)
        with self.assertRaisesRegex(_codec.CodecError, 'does not fit 64 bits'):
            _codec.decode_ranges('1-' + '9' * 30)
        self.assertEqual(18446744073709551615, _codec.parse_int('18446744073709551615'))
        with self.assertRaises(_codec.CodecError):
            _codec.parse_int('18446744073709551616')

    def test_the_codecs_bound_their_tokens(self):
        for call, token in ((_codec.parse_int, '9' * 5000),
                            (_codec.decode_ranges, '1-' + '9' * 5000),
                            (_codec.decode_ranges, '9' * 5000),
                            (_codec.parse_num, '9' * 5000),
                            (_codec.decode_kvalue, 'i:' + '9' * 5000),
                            (_codec.decode_kvalue, 's:' + '%' * 5000),
                            (_codec.pct_unescape, '%zz' + 'x' * 5000)):
            with self.subTest(call=call.__name__, n=len(token)):
                with self.assertRaises(_codec.CodecError) as cm:
                    call(token)
                self.assertLess(len(str(cm.exception)), 120, str(cm.exception)[:200])
        self.assertEqual("'abc'", _codec.tok('abc'))
        self.assertEqual("'xxxxx'... (100 characters)", _codec.tok('x' * 100, 5))
        self.assertEqual('an integer of 16610 bits', _codec.tok(10 ** 5000))

    def test_a_long_record_is_not_echoed(self):
        text = T.doc_text('switch_chain')
        bad, line = self.mutate(text, r'^(E \d+) [a-z]', r'\1 q' + 'q' * 5000)
        with self.assertRaises(GraphletFormatError) as cm:
            parse(bad)
        self.assertEqual(line, cm.exception.line_no)
        self.assertLess(len(str(cm.exception)), 200)


# ------------------------------------------------------------------ BUG-2-LONE-SURROGATE

class TestLoneSurrogates(unittest.TestCase):
    def body(self, name='merge'):
        text = T.doc_text(name)
        lines = text.split('\n')
        i = next(k for k, l in enumerate(lines) if l.startswith('L '))
        lines[i] += '\udc80'
        return '\n'.join(lines), i + 1

    def test_parse_names_the_line(self):
        text, line = self.body()
        with self.assertRaises(GraphletFormatError) as cm:
            parse(text)
        self.assertEqual(line, cm.exception.line_no)
        self.assertIn('not UTF-8', str(cm.exception))

    def test_from_response_and_the_store(self):
        text, line = self.body()
        resp = T.doc_json('merge', 'graphlet')
        r = resp['results'][0]
        r['graphlet'] = text
        r['graphlet_bytes'] = len(text.encode('utf-8', 'surrogatepass'))
        with self.assertRaises(GraphletFormatError):
            from_response(r, resp)
        with tempfile.TemporaryDirectory() as d:
            store = GraphletStore(d)
            with self.assertRaises(GraphletFormatError):
                store.put(resp, {'seeds': []})
            self.assertEqual([], store.list())


# ------------------------------------------------------------------ BUG-5-RUN-ERROR-LINE

class TestRunErrorsNameTheirLine(unittest.TestCase):
    def r_lines(self, text):
        return [i for i, l in enumerate(text.split('\n')) if l.startswith('R ')]

    def test_a_run_outside_its_anchor(self):
        text = T.doc_text('merge')
        lines = text.split('\n')
        i = self.r_lines(text)[0]
        f = lines[i].split(' ')
        for to_bp, what in (('99999', 'outside its anchor'), ):
            g = list(f)
            g[4] = to_bp
            bad = '\n'.join(lines[:i] + [' '.join(g)] + lines[i + 1:])
            with self.assertRaises(GraphletFormatError) as cm:
                parse(bad)
            self.assertEqual(i + 1, cm.exception.line_no)
            self.assertRegex(str(cm.exception), r'^line %d: right arm: run 0 ends at 99999 '
                             r'outside its anchor' % (i + 1))

    def test_a_run_that_ends_before_it_starts(self):
        text = T.doc_text('merge')
        lines = text.split('\n')
        for i in self.r_lines(text):
            f = lines[i].split(' ')
            if int(f[4]) > 0:
                break
        f[3] = str(int(f[4]) + 1)
        bad = '\n'.join(lines[:i] + [' '.join(f)] + lines[i + 1:])
        with self.assertRaises(GraphletFormatError) as cm:
            parse(bad)
        self.assertEqual(i + 1, cm.exception.line_no)
        self.assertIn('arm: run', str(cm.exception))
        self.assertIn('ends before it starts', str(cm.exception))

    def test_runs_in_annotate_mode(self):
        lines = ANNOT_MERGE.split('\n')
        at = lines.index('T D .') + 1
        lines.insert(at, 'R 0 0 0 5 D 0 * * 0 0 0')
        lines[-2] = 'Z 19'
        with self.assertRaises(GraphletFormatError) as cm:
            parse('\n'.join(lines))
        self.assertEqual(at + 1, cm.exception.line_no)
        self.assertIn('right arm: R records in annotate mode', str(cm.exception))


# ------------------------------------------------------------------ the tool layer

class Canned(TraverseClient):
    """A backend answering every /traverse with one (possibly corrupt) response."""

    def __init__(self, response):
        super().__init__('canned', 0)
        self.response = response

    def traverse_raw(self, request):
        return copy.deepcopy(self.response)

    def capabilities(self):
        return {'release': '', 'index_fp': None}

    def resolve(self, sequence, **kw):
        return {'k': 5, 'num_kmers': 6, 'labels': [], 'candidates': []}


def corrupt_response(name, edit):
    resp = copy.deepcopy(T.doc_json(name, 'graphlet'))
    r = resp['results'][0]
    r['graphlet'] = edit(r['graphlet'])
    r.pop('graphlet_bytes', None)
    r.pop('graphlet_lines', None)
    return resp


class ToolCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.store = GraphletStore(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def tools(self, response=None, **kw):
        clients = {} if response is None else {'canned': Canned(response)}
        return GraphletTools(self.store, clients, **kw)


class TestCorruptBodiesThroughTheTools(ToolCase):
    """BUG-1, BUG-2 and FORMAT-ERROR-ECHO through traverse_fetch: format_error (never
    bad_argument, never result_too_large), short, nothing stored."""

    def fetch(self, resp, **kw):
        out = self.tools(resp).traverse_fetch('canned', {'sequence': 'ACGT'}, **kw)
        self.assertLessEqual(size(out), kw.get('max_bytes', 2048))
        self.assertEqual([], self.store.list())
        return out

    def test_huge_tokens(self):
        for junk in ('9' * 5000, 'x' * 3000):
            out = self.fetch(corrupt_response(
                'merge', lambda t: re.sub(r'^Z \d+$', 'Z ' + junk, t, flags=re.M)))
            self.assertEqual('format_error', out['error'], out)
            self.assertLess(len(out['message']), 200)

    def test_a_lone_surrogate(self):
        out = self.fetch(corrupt_response('merge', lambda t: t.replace('\nL ', '\nL \udc80',
                                                                       1)))
        self.assertEqual('format_error', out['error'], out)
        resp = corrupt_response('merge', lambda t: t + '\udc80')
        resp['results'][0]['graphlet_bytes'] = 10
        self.assertEqual('format_error', self.fetch(resp)['error'])

    def test_an_error_keeps_its_code_under_a_small_ceiling(self):
        resp = corrupt_response('merge', lambda t: re.sub(r'^Z \d+$', 'Z ' + 'x' * 3000, t,
                                                          flags=re.M))
        for mb in (MIN_MAX_BYTES, 80, 120, 300):
            out = self.fetch(resp, max_bytes=mb)
            self.assertEqual('format_error', out['error'], (mb, out))
            self.assertIsInstance(out['message'], str)


class TestContinuationDryRun(ToolCase):
    """traverse_continue(execute=False) returns the request (the full normalized one,
    equal to next_request) under the 16 KB sequence ceiling it states; max_bytes is
    accepted; a result_too_large hint names what the tool takes."""

    def long_names(self):
        resp = copy.deepcopy(T.doc_json('switch_chain', 'graphlet'))
        r = resp['results'][0]
        lines = r['graphlet'].split('\n')
        for i, line in enumerate(lines):
            if line.startswith('L '):
                f = line.split(' ', 5)
                lines[i] = ' '.join(f[:4] + ['0', '/data/refseq/release-224/' * 60 + f[5]])
        r['graphlet'] = '\n'.join(lines)
        r['graphlet_bytes'] = len(r['graphlet'].encode('utf-8'))
        return resp

    def test_a_request_over_2_kb_is_delivered(self):
        resp = self.long_names()
        tools = self.tools(resp)
        h = tools.traverse_fetch('canned', {'sequence': 'ACGT'})['handle']
        g = self.store.graphlet(h)
        pid = next(p.id for p in derive.paths(g.arms['right'])
                   if g.arms['right'].segments[p.leaf].leaf.continuation is not None)
        want = g.next_request('right', [pid])
        self.assertGreater(size(want), 2048)
        out = tools.traverse_continue(h, 'right', pid, execute=False)
        self.assertEqual(want, out['request'])
        self.assertEqual(16 * 1024, out['ceiling_bytes'])
        self.assertLessEqual(size(out), out['ceiling_bytes'])
        # a smaller ceiling: result_too_large, and the hint is one this tool can act on
        small = tools.traverse_continue(h, 'right', pid, execute=False, max_bytes=1024)
        self.assertEqual('result_too_large', small['error'])
        self.assertIn('raise max_bytes', small['message'])
        big = tools.traverse_continue(h, 'right', pid, execute=False, max_bytes=1 << 16)
        self.assertEqual(want, big['request'])

    def test_max_bytes_on_the_backend_tools(self):
        tools = self.tools(T.doc_json('merge', 'graphlet'))
        self.assertNotIn('error', tools.traverse_capabilities('canned', max_bytes=4096))
        self.assertNotIn('error', tools.traverse_resolve('canned', 'ACGTACGT', max_bytes=4096))
        default = tools.traverse_fetch('canned', {'sequence': 'ACGT'})
        out = tools.traverse_fetch('canned', {'sequence': 'ACGT'}, max_bytes=1000)
        self.assertNotIn('error', out, out)
        self.assertLessEqual(size(out), 1000)
        self.assertLess(size(out['summary']), size(default['summary']))

    def test_an_argument_a_tool_does_not_take_is_a_result(self):
        tools = self.tools(T.doc_json('merge', 'graphlet'))
        h = tools.traverse_fetch('canned', {'sequence': 'ACGT'})['handle']
        for out in (tools.graphlet_export(h, max_bytes=4096), tools.graphlet_free(h, x=1),
                    tools.traverse_continue(h, 'right', 0, bogus=True)):
            self.assertEqual('bad_argument', out['error'], out)
            self.assertIsInstance(out['message'], str)

    def test_the_hint_of_a_tool_without_max_bytes(self):
        from metagraph.traverse.mcp_tools import _too_large
        self.assertIn('raise max_bytes', _too_large(5000, 2048, True)['message'])
        self.assertNotIn('raise max_bytes', _too_large(5000, 2048, False)['message'])


class TestTinyMaxBytes(ToolCase):
    """BUG-4: every return is <= max(max_bytes, MIN_MAX_BYTES); below the floor a short
    bad_argument."""

    def test_every_return_fits(self):
        tools = self.tools(T.doc_json('annotate', 'graphlet'))
        h = tools.traverse_fetch('canned', {'sequence': 'ACGT'},
                                 strategy={'labels': {'mode': 'annotate'}})['handle']
        calls = [lambda mb: tools.graphlet_walks(h, 'left', spell='full', max_bytes=mb),
                 lambda mb: tools.graphlet_claims(h, 'right', max_bytes=mb),
                 lambda mb: tools.graphlet_summary(h, max_bytes=mb),
                 lambda mb: tools.graphlet_labels(h, max_bytes=mb),
                 lambda mb: tools.graphlet_walk(h, 'right', 0, max_bytes=mb),
                 lambda mb: tools.graphlet_support(h, 'right', 0, max_bytes=mb),
                 lambda mb: tools.graphlet_splits(h, 'right', max_bytes=mb),
                 lambda mb: tools.graphlet_compare(h, h, max_bytes=mb),
                 lambda mb: tools.graphlet_list(max_bytes=mb),
                 lambda mb: tools.traverse_capabilities('canned', max_bytes=mb),
                 lambda mb: tools.traverse_fetch('canned', {'sequence': 'ACGT'},
                                                 keep=False, max_bytes=mb)]
        for call in calls:
            for mb in (0, 1, 10, 45, 56, MIN_MAX_BYTES - 1, MIN_MAX_BYTES, 65, 80, 100,
                       128, 182, 200, 300):
                out = call(mb)
                with self.subTest(mb=mb, out=out):
                    self.assertLessEqual(size(out), max(mb, MIN_MAX_BYTES))
                    if mb < MIN_MAX_BYTES:
                        self.assertEqual('bad_argument', out['error'])
                    if 'error' in out:
                        self.assertIsInstance(out['message'], str)
        for mb in (True, 'x', 1.5):
            self.assertEqual('bad_argument',
                             tools.graphlet_summary(h, max_bytes=mb)['error'])


class TestWalkPages(ToolCase):
    """B3: a page of graphlet_walks ranks the arm once (cached per handle and arguments)
    and builds only its own rows."""

    def test_pages_equal_the_ranked_walks_and_build_only_their_rows(self):
        tools = self.tools(T.doc_json('fork', 'graphlet'))
        h = tools.traverse_fetch('canned', {'sequence': 'ACGT'})['handle']
        g = self.store.graphlet(h)
        built = []
        real = ops.walks_at

        def counting(g_, arm, ids, with_claims=False):
            built.extend(ids)
            return real(g_, arm, ids, with_claims)
        for side in g.arms:
            for rank in ('support', 'length', 'loss', 'id'):
                want = [w.path_id for w in g.walks(side, by=rank)]
                ids, cursor, pages = [], None, 0
                ops.walks_at = counting
                try:
                    while True:
                        del built[:]
                        out = tools.graphlet_walks(h, side, rank=rank, n=1, cursor=cursor)
                        self.assertNotIn('error', out, out)
                        # the page's row, and at most one more looked at for the cursor
                        self.assertLessEqual(len(built), 2)
                        ids.extend(r['walk'] for r in out['rows'])
                        pages += 1
                        cursor = out.get('next_cursor')
                        if cursor is None:
                            break
                finally:
                    ops.walks_at = real
                self.assertEqual(want, ids, (side, rank))
                self.assertEqual(len(want), pages)

    def test_the_ranking_is_cached_and_the_filter_counted(self):
        tools = self.tools(T.doc_json('merge', 'graphlet'))
        h = tools.traverse_fetch('canned', {'sequence': 'ACGT'})['handle']
        w = tools.graphlet_walks(h, 'right', label={'ref': 'c:1'})
        self.assertEqual((0, {'merge_entered': 1}), (w['total'], w['filtered']))
        self.assertEqual(1, len(tools._ranked))
        calls = []
        real = ops.rank_walks
        ops.rank_walks = lambda *a, **k: calls.append(1) or real(*a, **k)
        try:
            tools.graphlet_walks(h, 'right', label={'ref': 'c:1'})
        finally:
            ops.rank_walks = real
        self.assertEqual([], calls)


# ------------------------------------------------------------------ B2-MEMORY-BYTES

class TestMemoryAccounting(unittest.TestCase):
    def test_star_sets_are_shared(self):
        g = T.graphlet('merge')
        for a in g.arms.values():
            self.assertLessEqual(len({id(r.stars) for r in a.runs}), 2)
            self.assertLessEqual(len({id(s.stars) for s in a.segments}),
                                 len({s.stars for s in a.segments}))

    def test_the_store_charges_what_the_queries_built(self):
        with tempfile.TemporaryDirectory() as d:
            store = GraphletStore(d)
            h = store.put(T.doc_json('merge', 'graphlet'), None)
            before = store._ram_bytes
            g = store.graphlet(h)
            g.walks('right')
            g.claims('right')
            g.label_summary()
            store.graphlet(h)                   # the next access re-measures
            self.assertGreater(store._ram_bytes, before)
            # the caches measured on their own count what they share with the model
            # again: a slight overcharge, never an undercharge
            self.assertGreaterEqual(store._ram[h][1], g.memory_bytes())
            self.assertLess(store._ram[h][1], 1.05 * g.memory_bytes())

    def test_grown_caches_are_evicted_under_the_budget(self):
        with tempfile.TemporaryDirectory() as d:
            g0 = T.graphlet('merge')
            one = g0.memory_bytes()
            g0.walks('right')
            g0.claims('right')
            g0.label_summary()
            grown = g0.memory_bytes()
            self.assertGreater(grown, one)
            # room for two cold models, not for two with their caches
            store = GraphletStore(d, max_ram_mb=(2 * one + (grown - one) // 2) / (1 << 20))
            ha = store.put(T.doc_json('merge', 'graphlet'), None)
            hb = store.put(T.doc_json('merge', 'graphlet'), None)
            self.assertTrue(store.in_ram(ha) and store.in_ram(hb))
            ga = store.graphlet(ha)
            store.graphlet(hb)                  # b is now the most recent
            ga.walks('right')
            ga.claims('right')
            ga.label_summary()
            store.refresh()
            self.assertLessEqual(store._ram_bytes, store.max_ram_bytes)
            self.assertFalse(store.in_ram(ha) and store.in_ram(hb))

    def test_tool_calls_refresh_the_accounting(self):
        with tempfile.TemporaryDirectory() as d:
            store = GraphletStore(d)
            tools = GraphletTools(store, {'canned': Canned(T.doc_json('merge', 'graphlet'))})
            h = tools.traverse_fetch('canned', {'sequence': 'ACGT'})['handle']
            before = store._ram_bytes
            tools.graphlet_claims(h, 'right')
            tools.graphlet_labels(h)
            self.assertGreater(store._ram_bytes, before)
            real = store.graphlet(h).memory_bytes()
            self.assertGreaterEqual(store._ram[h][1], real)
            self.assertLess(store._ram[h][1], 1.05 * real)


if __name__ == '__main__':
    unittest.main()
