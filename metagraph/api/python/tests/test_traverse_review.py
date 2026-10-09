"""Regression tests for the external (GPT) implementation review of the graphlet library
(commit cb9eac2d), each reproduced offline: the review's reproducers turned into tests on
hand-made documents, CLI fixtures or the CLI documents the reproducers produced (embedded
below), and a stand-in backend. One class per finding, named by its number in the review.
Finding 6 (per-arm attribution of pair evaluations) was fixed by stage 2 and is covered
there; findings 1 (the server side), 8 and 9 are C++ (tests/cli/test_graphlet_codec.cpp,
integration_tests/test_traverse.py)."""

import ast
import copy
import gc
import json
import os
import sys
import tempfile
import tracemalloc
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402

from metagraph.traverse import (  # noqa: E402
    GraphletStore, Label, derive, from_response, ops, parse,
)
from metagraph.traverse import _codec  # noqa: E402
from metagraph.traverse.mcp_tools import GraphletTools  # noqa: E402
from metagraph.traverse.model import (  # noqa: E402
    Branch, GrowthBin, Limitation, SplitPoint, UnverifiableLabelName,
)
from test_traverse_mcp import CONTINUATION, SEED, FakeClient  # noqa: E402


# ------------------------------------------------------------------ the reviewer's documents
#
# /tmp/graphlet-evidence-review/repro.py, re-run with build_debug/metagraph (the same bytes
# as the reviewer's): the CLI fixture indexes 'chain' and 'bubbles' of
# scripts/traversal/graphlet_fixtures.py, with an index manifest.

# 'chain', annotate, keep, radius 10: only A is recorded within [0, 10)
CHAIN_10 = '''\
H mgt 1 15 basic $ACGT a k k 1000 1000 0 walk evidence-review-chain 778b7819ce419a3f3fdd54fc23d3e272e011476b4cb63620fa4f1465e8f4f427 d99c945a0538b4eb
S cf89642ae331e875 30 16 0 ACCAGTCTTAGAATCTGATGACCGCCCGAC
L h 0 0 0 A
O c c c i
A r c 10 p 0 0 1 1 0 steps=10,successor_enumerations=10,output_bp=10,pair_evaluations=0,edge_reuse_probes=0,reminimisation_rounds=0,max_reminimisation_rounds=0,refusal_scans=0,switch_sources_cut=0 * 0 * 1 0 1 0 0 10
B 0 1 1 1 1 10 0 0 0 0 0 0 0 .
G * 0 10 0 * * * * * GCAGGCGATC
P 0 10 1 0
T X .
C 15 0 0 0
Z 11
'''

# the same, radius 80: B first recorded at 24, C and D at 54
CHAIN_80 = '''\
H mgt 1 15 basic $ACGT a k k 1000 1000 0 walk evidence-review-chain 778b7819ce419a3f3fdd54fc23d3e272e011476b4cb63620fa4f1465e8f4f427 d99c945a0538b4eb
S cf89642ae331e875 30 16 0 ACCAGTCTTAGAATCTGATGACCGCCCGAC
L h 0 0 0 A
L h 0 1 0 B
L h 0 2 0 C
L h 0 3 0 D
O c c c i
A r c 80 p 0 0 1 3 0 steps=100,successor_enumerations=99,output_bp=100,pair_evaluations=0,edge_reuse_probes=0,reminimisation_rounds=0,max_reminimisation_rounds=0,refusal_scans=0,switch_sources_cut=0 * 0 * 3 0 2 1 0 100
B 0 2 3 3 1 100 1 0 1 0 0 0 0 .
G * 0 60 0 * * * 0 * GCAGGCGATCGCACTACTAGCCATAAGCCGCAAGAAGCCTATCCAATTCCTGTGATTTAG
P 0 24 1 0
P 24 30 2 0-1
P 30 54 1 1
P 54 60 3 1-3
G 0 60 20 3 * * * * * AAGGATTCAGCGCTGACTGG
P 60 80 1 3
T X .
C 80 0 0 .
G 0 60 20 2 * * * * * TTACTCGTAGCCGGGCGTGA
P 60 80 1 2
T X .
C 80 0 0 .
Z 23
'''


# 'bubbles' (the merge fixture's locus), annotate, on_reconverge merge, radius 100, lists
# uncapped (max_labels_per_node 64): label_evidence complete. b.fa takes the second allele
# at both bubbles (segments 2 and 5), c.fa at the second (5): both reach 100 along their
# own routes, entering the displayed chain 0-1-3-4-6 through the non-first parent of the
# merge at 62
BUBBLES_ANNOTATE_MERGE = '''\
H mgt 1 15 basic $ACGT a k m 64 20 0 walk evidence-review-bubbles 836342bdffd8dc906035e41a6e5c37a7c23054cb4662ba58be5d12ffb1239bbe fbfd396e439f0d22
S 58b906f959ef3993 30 16 0 ATGAGCGTGCCTCTAACGTTCATACCAAAT
L c 0 * 0 a.fa
L c 1 * 0 b.fa
L c 2 * 0 c.fa
L c 3 * 0 both.fa
O p c c i
K r scope branching.on_reconverge s:merge i:2 * . completeness holds for the united-history rule, not per path; use "keep" for the per-path guarantee
A r c 100 u 0 0 1 4 0 steps=132,successor_enumerations=130,output_bp=132,pair_evaluations=0,edge_reuse_probes=0,reminimisation_rounds=0,max_reminimisation_rounds=0,refusal_scans=0,switch_sources_cut=0 * 0 * 7 0 1 2 2 132
B 0 2 4 5 1 70 2 0 2 1 0 0 0 .
B 50 2 4 5 1 62 0 0 0 1 0 0 0 .
B 100 1 4 4 1 0 0 0 0 0 0 0 0 .
G * 0 20 0-3 * * * 0 * TTTCCTATTTAGCCTCTGTC
P 0 20 4 !
G 0 20 16 !1 * * * * * ATTACGTTTGACAATG
P 20 35 3 !
P 35 36 4 0-3
G 0 20 16 1,3 * * * * * TTTACGTTTGACAATG
P 20 35 2 !
P 35 36 4 0-3
G 1,2 36 10 ! * * * 0 * ACCCAGCCCT
P 36 46 4 !
G 3 46 16 0,3 * * * * * CCGGCGGGTCGACTTG
P 46 61 2 !
P 61 62 4 0-3
G 3 46 16 !0 * * * * * GCGGCGGGTCGACTTG
P 46 61 3 !
P 61 62 4 0-3
G 4,5 62 38 ! * * * * * GTCCGGACGATAGCACTTAGTTCCTCACTTCACAATAG
P 62 100 4 !
T X .
C 20 0 0 !
Z 33
'''


def claim_rows(claims):
    return sorted((c.label.name, c.from_bp, c.to_bp, c.route_bp, c.evidence_from, c.kind,
                   c.end_class) for c in claims)


# ------------------------------------------------------------------ finding 1 (library)

class TestFinding1NamesAreNotResolvedUnverified(unittest.TestCase):
    """/traverse resolves label NAMES. A replaced name ('bad\\ufffd', what servers wrote
    for 'bad\\xff') was another label's name in the review's probe, and the continuation
    went on under that label. The server now refuses such a seed (C++ and integration
    tests); the library refuses to build a continuation or next_request naming labels it
    cannot verify -- a name shared by two labels of the retrieval, or one with U+FFFD --
    and never resolves one silently."""

    def renamed(self, name, i, to):
        g = T.graphlet(name)
        g.labels[i].name = to
        return g

    def test_a_shared_name_is_refused(self):
        twin = self.renamed('fork', 2, 'acc1')        # the docs' next-request-ambiguous
        with self.assertRaises(UnverifiableLabelName) as cm:
            twin.next_request('right', [1])
        self.assertIn('h:0:0', str(cm.exception))
        self.assertIsInstance(cm.exception, ValueError)
        with self.assertRaises(UnverifiableLabelName):
            twin.continuation('right', 1).as_seed()

    def test_a_replacement_character_is_refused(self):
        sw = self.renamed('switch_chain', 2, 'C�')    # walk 1 continues under C
        self.assertEqual(['C�'], [l.name for l in sw.continuation('right', 1).labels])
        with self.assertRaises(UnverifiableLabelName) as cm:
            sw.next_request('right', [1])
        self.assertIn('U+FFFD', str(cm.exception))
        # a name elsewhere in the dictionary does not stop another walk's request
        self.assertNotIn(2, [l.id for l in sw.continuation('right', 0).labels])
        self.assertTrue(sw.next_request('right', [0])['seeds'])

    def test_verified_names_are_unchanged(self):
        sw = T.graphlet('switch_chain')
        self.assertEqual([{'sequence': CONTINUATION, 'labels': ['C']}],
                         sw.next_request('right', [1])['seeds'])
        self.assertEqual(['C'], ops.resubmittable_names(sw, ['C']))

    def test_annotate_requests_name_no_labels(self):
        # same_name records two headers named ACC1 (different columns): an annotate
        # continuation names no labels, so it is built
        g = T.graphlet('same_name')
        for p in derive.paths(g.arms['right']):
            if g.arms['right'].segments[p.leaf].leaf.continuation is not None:
                req = g.next_request('right', [p.id])
                self.assertNotIn('labels', req['seeds'][0])

    def test_the_tools_refuse_without_calling_the_backend(self):
        with tempfile.TemporaryDirectory() as d:
            store = GraphletStore(d)
            client = FakeClient()
            tools = GraphletTools(store, {'mini': client})
            sw = self.renamed('switch_chain', 2, 'C�')
            resp = copy.deepcopy(T.doc_json('switch_chain', 'graphlet'))
            h = store.put_parsed(sw, resp['results'][0]['graphlet'],
                                 {'seeds': [{'sequence': SEED['switch_chain']}]},
                                 source='mini')
            for execute in (False, True):
                out = tools.traverse_continue(h, 'right', 1, execute=execute)
                self.assertEqual('unverifiable_label_name', out['error'], out)
            self.assertEqual([], client.requests)

    def test_a_loaded_file_is_not_replayed_by_unverified_names(self):
        sw = self.renamed('switch_chain', 0, 'A�')       # a seed label
        self.assertIsNone(GraphletStore.request_of(sw))
        self.assertIsNotNone(GraphletStore.request_of(T.graphlet('switch_chain')))


# ------------------------------------------------------------------ finding 2

def _traced(text, queries=False):
    """(memory_bytes, the heap tracemalloc sees the parse -- and the queries -- retain)."""
    gc.collect()
    tracemalloc.start()
    try:
        base = tracemalloc.get_traced_memory()[0]
        g = parse(text)
        if queries:
            for side in g.arms:
                g.walks(side)
                g.claims(side, strict=False)
                for p in derive.paths(g.arms[side])[:20]:
                    g.support_profile(side, p.id)
            g.label_summary()
        gc.collect()
        held = tracemalloc.get_traced_memory()[0] - base
    finally:
        tracemalloc.stop()
    return g.memory_bytes(), held


def growth_bins_document(n=20000):
    """The review's case: a model whose retained size is its growth bins (one bin per
    base of a 20,000 bp walk there; here 20,000 empty bins added to a fixture, which
    keeps every B sum the body is checked against)."""
    g = parse(T.doc_text('linear'))
    a = g.arms['right']
    for i in range(n):
        a.growth.append(GrowthBin(10 ** 6 + i, 0, 0, 0, True, 0, 0, 0, 0, 0, 0, 0, 0, {}))
    return g.dump(envelope=False)


class TestFinding2TheStoreChargesWhatIsRetained(unittest.TestCase):
    """memory_bytes() omitted growth bins (and every other non-segment structure): a
    one-segment retrieval with 20,001 bins retained 4.9 MB but reported 21 KB and stayed
    resident under a 0.1 MB store. It is now a deep, deduplicating count of the model and
    its caches (shared interned sets once), without the program's own constants."""

    def test_the_growth_bin_case_is_not_kept_under_the_budget(self):
        text = growth_bins_document()
        m, held = _traced(text)
        self.assertGreater(held, 4 * 10 ** 6)
        self.assertLess(abs(m - held), 0.25 * held, (m, held))
        g = parse(text)
        with tempfile.TemporaryDirectory() as d:
            store = GraphletStore(d, max_ram_mb=0.1)           # 104,857 bytes
            self.assertEqual(104857, store.max_ram_bytes)
            h = store.put_parsed(g, text, None)
            self.assertFalse(store.in_ram(h))
            self.assertLessEqual(store._ram_bytes, store.max_ram_bytes)
            # still answered: parsed again on demand, and not kept
            self.assertEqual(20000 + 3, len(store.graphlet(h).arms['right'].growth))
            self.assertFalse(store.in_ram(h))

    def test_memory_bytes_is_within_a_quarter_of_the_traced_heap(self):
        docs = [(name, T.doc_text(name)) for name in (
            'merge', 'annotate', 'caps', 'limits', 'quorum', 'linear', 'fork',
            'switch_chain', 'same_name', 'no_bins', 'reentry', 'ambiguous_split')]
        docs.append(('growth_bins', growth_bins_document(5000)))
        docs.append(('bubbles_annotate_merge', BUBBLES_ANNOTATE_MERGE))
        for name, text in docs:
            for queries in (False, True):
                m, held = _traced(text, queries)
                with self.subTest(doc=name, queries=queries, reported=m, traced=held):
                    self.assertLess(abs(m - held), 0.25 * held)


# ------------------------------------------------------------------ finding 3

class TestFinding3AnnotateRoutesThroughMerges(unittest.TestCase):
    """Annotate + merge: claims() and routes() followed first parents only, so b.fa and
    c.fa (direct_bp = reach_bp = 100) had no claim or route beyond the first bubble. ROUTE
    evidence is now derived through every supported parent (the union pass of the
    normative annotate label_summary) and kept apart from DISPLAYED evidence: the claim
    has route support [0, 100) and displayed support from the merge at 62."""

    def setUp(self):
        self.g = parse(BUBBLES_ANNOTATE_MERGE)

    def test_route_claims_reach_direct_bp(self):
        g = self.g
        self.assertEqual('complete', g.outcome.label_evidence)
        got = claim_rows(g.claims('right'))
        self.assertEqual([('a.fa', 0, 100, 0, 0, 'alive', 'radius'),
                          ('b.fa', 0, 100, 62, 62, 'alive', 'radius'),
                          ('both.fa', 0, 100, 0, 0, 'alive', 'radius'),
                          ('c.fa', 0, 100, 62, 62, 'alive', 'radius')], got)
        summary = g.label_summary()
        for lab in g.labels:
            longest = max(c.to_bp for c in g.claims('right', labels=[lab]))
            self.assertEqual(summary[lab.ref]['right']['direct_bp'], longest, lab)

    def test_the_constrain_retrieval_of_the_locus_agrees(self):
        # the merge fixture: the same locus and radius in constrain mode (lineages, G
        # partitions): its run ends have the same route support ...
        merge = T.graphlet('merge')
        want = [r for r in claim_rows(merge.claims('right')) if r[5] != 'merged']
        route = [(n, f, t, k, e) for n, f, t, _, _, k, e in want]
        self.assertEqual(route, [(n, f, t, k, e) for n, f, t, _, _, k, e
                                 in claim_rows(self.g.claims('right'))])
        # ... and the same displayed support once both display the same parent at each
        # merge: the fixture (feature level 6) stores the parent carried by the most labels
        # first (R21 (4)) -- at 62 the G allele (b.fa, c.fa, both.fa), which arrived second
        # and comes second in this older document
        self.assertNotEqual(want, claim_rows(self.g.claims('right')))
        level6 = parse(BUBBLES_ANNOTATE_MERGE.replace('G 4,5 62 38 !', 'G 5,4 62 38 !'))
        self.assertEqual(want, claim_rows(level6.claims('right')))

    def test_displayed_evidence_stays_displayed(self):
        g = self.g
        (w,) = g.walks('right')
        self.assertEqual(['a.fa', 'both.fa'], [l.name for l in w.labels_full])
        self.assertEqual([], g.walks('right', labels=['b.fa']))
        self.assertEqual([0], [x.path_id for x in g.walks('right', labels=['b.fa'],
                                                          route_consistent=False)])
        self.assertEqual({'b.fa', 'c.fa'}, {c.label.name for c in w.claims if c.route_bp})

    def test_cut_claims_are_route_only_without_displayed_support(self):
        g = self.g
        at50 = {c.label.name: c for c in g.claims('right', at_most_bp=50)}
        self.assertEqual(('route_only', None, 50, 'open'),
                         (at50['b.fa'].kind, at50['b.fa'].evidence_from, at50['b.fa'].to_bp,
                          at50['b.fa'].end_class))
        self.assertEqual(('stretch', 0), (at50['a.fa'].kind, at50['a.fa'].evidence_from))
        at70 = {c.label.name: c for c in g.claims('right', at_most_bp=70)}
        self.assertEqual(('stretch', 62), (at70['c.fa'].kind, at70['c.fa'].evidence_from))

    def test_routes_and_label_walks(self):
        g = self.g
        self.assertEqual([[0, 2, 3, 5, 6]], g.routes('b.fa', 'right'))
        self.assertEqual([[0, 1, 3, 5, 6]], g.routes('c.fa', 'right'))
        self.assertEqual([[0, 1, 3, 4, 6]], g.routes('a.fa', 'right'))
        (lw,) = g.label_walks('b.fa')
        self.assertEqual((0, 100, 62, 'max_extension_bp', [6]),
                         (lw.from_bp, lw.to_bp, lw.evidence_from, lw.end, lw.leaves_below))
        (route, bases), = g.routes('b.fa', 'right', spell=True)
        self.assertEqual(100, len(bases))
        self.assertNotEqual(g.spell('right', 0), bases)
        self.assertEqual(g.spell('right', 0)[62:], bases[62:])

    def test_the_tool_counts_them_as_merge_entered(self):
        with tempfile.TemporaryDirectory() as d:
            store = GraphletStore(d)
            h = store.put_graphlet(parse(BUBBLES_ANNOTATE_MERGE))
            tools = GraphletTools(store)
            c = tools.graphlet_claims(h, 'right', max_bytes=8192)
            self.assertEqual({'merge_entered': 2}, c['filtered'])
            self.assertEqual(['a.fa', 'both.fa'], [r['label']['name'] for r in c['rows']])
            every = tools.graphlet_claims(h, 'right', route_consistent=False, max_bytes=8192)
            self.assertEqual(4, every['total'])
            rows = tools.graphlet_labels(h, name='b.fa', max_bytes=8192)['rows']
            self.assertEqual([(100, 62, [0, 2, 3, 5, 6])],
                             [(r['to_bp'], r['evidence_from'], r['route']) for r in rows])

    def test_a_tree_is_unchanged(self):
        # keep mode: every route is the displayed chain, so the claims are the displayed
        # stretches (evidence from 0) -- the oracle of §6.9
        for doc in (parse(T.compare_text('oracle_annotate')), parse(T.doc_text('same_name'))):
            for side in doc.arms:
                for c in doc.claims(side, strict=False):
                    self.assertEqual(0, c.route_bp)
                    self.assertIn(c.evidence_from, (0, None))


# ------------------------------------------------------------------ finding 4

class TestFinding4LabelComparisonsAreClippedFirst(unittest.TestCase):
    """Comparing the chain locus at radii 10 and 80 used depth 10 but reported B, C and D
    (first recorded at 24 and 54) in the deeper retrieval only: the completed summaries
    were clamped. Presence intervals and runs are now clipped to the depth BEFORE the
    summaries are derived."""

    def test_the_reviewers_case(self):
        shallow, deep = parse(CHAIN_10), parse(CHAIN_80)
        cmp = shallow.compare(deep, mode='labels')
        self.assertEqual((True, True, 10, [], [], []),
                         (cmp.comparable, cmp.equal, cmp.depth_used, cmp.only_in_a,
                          cmp.only_in_b, cmp.differ))
        self.assertEqual({('right', 'h:0:0'): {'direct_bp': 10, 'reach_bp': 10}},
                         ops._label_keys(deep, ['right'], 10, None))
        # deeper, B appears with its own interval (24..)
        keys = ops._label_keys(deep, ['right'], 30, None)
        self.assertEqual({'direct_bp': 0, 'reach_bp': 30}, keys[('right', 'h:0:1')])
        self.assertNotIn(('right', 'h:0:2'), keys)

    def test_unclipped_it_is_label_summary(self):
        docs = [T.body(n) for n in ('merge', 'switch_chain', 'reentry', 'annotate', 'caps',
                                    'quorum', 'linear', 'fork', 'same_name', 'no_bins')]
        docs += [parse(CHAIN_80), parse(BUBBLES_ANNOTATE_MERGE)]
        for g in docs:
            rows = derive.label_summary(g)
            for side, a in g.arms.items():
                for lid, v in ops._clipped_summary(g, a, 10 ** 9).items():
                    r = rows[lid][side]
                    self.assertEqual((r['direct_bp'], r['reach_bp']),
                                     (v['direct_bp'], v['reach_bp']))

    def test_a_run_starting_beyond_the_depth_does_not_exist_there(self):
        sw = T.body('switch_chain')        # A -> B at 30, B -> C | D at 60
        names = {l.ref: l.name for l in sw.labels}
        for depth, want in ((30, {'A'}), (31, {'A', 'B'}), (61, {'A', 'B', 'C', 'D'})):
            keys = ops._label_keys(sw, ['right'], depth, None)
            self.assertEqual(want, {names[k[1]] for k in keys}, depth)
        # B's run [30, 60) is credited to the lineage root A (reach_bp), cut at the depth
        self.assertEqual({('right', 'h:0:0'): {'direct_bp': 30, 'reach_bp': 45},
                          ('right', 'h:0:1'): {'direct_bp': 0, 'reach_bp': 0}},
                         ops._label_keys(sw, ['right'], 45, None))


# ------------------------------------------------------------------ finding 7

class Caps(FakeClient):
    """The stand-in backend with an index identity of its own."""

    def __init__(self, fp=None, ns=None, release=''):
        super().__init__()
        self.fp = fp
        self.caps = {'release': release, 'index_fp': fp, 'index_ns': ns}

    def capabilities(self):
        return dict(self.caps)


class TestFinding7ContinuationChecksTheIndex(unittest.TestCase):
    """traverse_continue() called the backend without checking the parent's identity:
    after the backend's fingerprint changed, replay returned index_mismatch but the
    continuation succeeded and recorded a new-index result as derived from the old-index
    parent. The parent's identity is now checked against the capabilities before the
    request is sent and against the returned graphlet before the link is stored."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.store = GraphletStore(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def fetch(self, client):
        tools = GraphletTools(self.store, {'mini': client})
        out = tools.traverse_fetch(index='mini', seed={'sequence': SEED['switch_chain']})
        self.assertIn('handle', out, out)
        return tools, out['handle']

    def test_the_reviewers_case(self):
        client = Caps(fp='a' * 64)
        tools, parent = self.fetch(client)
        n = len(self.store.list())
        client.fp = client.caps['index_fp'] = 'b' * 64
        replay = tools.traverse_fetch(replay=parent)
        cont = tools.traverse_continue(parent, 'right', 1)
        self.assertEqual(('index_mismatch', 'index_mismatch'),
                         (replay.get('error'), cont.get('error')), (replay, cont))
        self.assertEqual(n, len(self.store.list()))           # nothing stored
        self.assertEqual(1, len(client.requests))              # nothing sent

    def test_the_same_index_is_verified(self):
        client = Caps(fp='a' * 64)
        tools, parent = self.fetch(client)
        cont = tools.traverse_continue(parent, 'right', 1)
        self.assertEqual({'verified': True, 'index_fp': 'a' * 64}, cont['identity'])
        self.assertEqual(parent, self.store.get(cont['handle']).parent['handle'])

    def test_unverifiable_on_both_is_stated_not_refused(self):
        client = Caps(fp=None)
        tools, parent = self.fetch(client)
        cont = tools.traverse_continue(parent, 'right', 1)
        self.assertIn('handle', cont, cont)
        self.assertFalse(cont['identity']['verified'])
        self.assertIn('cannot be verified', cont['identity']['note'])

    def test_a_manifest_gained_is_unverifiable(self):
        # refused as UNVERIFIABLE, not as proven different; the parent with a manifest and
        # the explicit opt-in: test_traverse_review3.py, TestD4OneSidedManifestIsUnverifiable
        client = Caps(fp=None)
        tools, parent = self.fetch(client)
        client.fp = client.caps['index_fp'] = 'a' * 64
        self.assertEqual('index_unverifiable',
                         tools.traverse_continue(parent, 'right', 1)['error'])

    def test_what_answered_is_checked_too(self):
        # the capabilities state the parent's digest, the response another one: the
        # result is not stored as the parent's continuation
        client = Caps(fp='a' * 64)
        tools, parent = self.fetch(client)
        n = len(self.store.list())
        client.fp = 'c' * 64                                   # caps still say 'a' * 64
        out = tools.traverse_continue(parent, 'right', 1)
        self.assertEqual('index_mismatch', out['error'], out)
        self.assertEqual(n, len(self.store.list()))

    def test_a_namespace_or_release_change_is_refused(self):
        client = Caps(fp=None)
        tools, parent = self.fetch(client)                    # ns 'fixtures' (H)
        client.caps['index_ns'] = 'another'
        self.assertEqual('index_mismatch', tools.traverse_continue(parent, 'right', 1)['error'])
        client.caps['index_ns'] = None
        self.store.get(parent).index['release'] = 'r1'
        client.caps['release'] = 'r2'
        self.assertEqual('release_mismatch',
                         tools.traverse_continue(parent, 'right', 1)['error'])


# ------------------------------------------------------------------ finding 8 (parity)

class TestFinding8RangesAtTheTopOfTheIdSpace(unittest.TestCase):
    """The C++ codec accepted '18446744073709551615,1' and encoded {MAX, 0} as 'MAX-0'
    (tests/cli: GraphletCodec.RangesAtTheTopOfTheIdSpace); the Python codec, which the
    review found right, keeps rejecting both."""

    def test_python_rejects_what_cpp_now_rejects(self):
        mx = 2 ** 64 - 1
        for tok in ('%d,1' % mx, '%d,0' % mx, '%d,%d' % (mx, mx), '%d-%d,0' % (mx - 1, mx),
                    '%d,%d' % (mx - 1, mx), '%d' % (mx + 1)):
            with self.assertRaises(_codec.CodecError, msg=tok):
                _codec.decode_ranges(tok)
        self.assertEqual([mx], _codec.decode_ranges('%d' % mx))
        self.assertEqual([0, mx], _codec.decode_ranges('0,%d' % mx))
        for ids in ([mx, 0], [mx - 1, mx, 0], [mx, mx]):
            with self.assertRaises(_codec.CodecError):
                _codec.encode_ranges(ids)
        self.assertEqual('%d-%d' % (mx - 1, mx), _codec.encode_ranges([mx - 1, mx]))


# ------------------------------------------------------------------ finding 10

class TestFinding10WalksComparedUnderLabels(unittest.TestCase):
    """compare(mode='walks', labels=...) kept the walks the filter emptied (with no
    labels), so the constrained trie over {acc1, acc3} differed from the label-free oracle
    by the oracle's acc2 walks although claims mode called them equal."""

    def test_the_constrained_trie_is_the_filtered_oracle(self):
        constrained = parse(T.compare_text('oracle_constrain'))
        oracle = parse(T.compare_text('oracle_annotate'))
        claims = constrained.compare(oracle, labels=['acc1', 'acc3'])
        walks = constrained.compare(oracle, labels=['acc1', 'acc3'], mode='walks')
        self.assertEqual((True, True), (claims.comparable, claims.equal))
        self.assertEqual((True, True, [], [], []),
                         (walks.comparable, walks.equal, walks.only_in_a, walks.only_in_b,
                          walks.differ))
        # the filter does drop something here: unfiltered, the acc2 walks differ
        every = constrained.compare(oracle, mode='walks')
        self.assertFalse(every.equal)
        keys = ops._walk_keys(oracle, ['right'], walks.depth_used, {'h:0:1'})
        self.assertTrue(all(v for v in keys.values()))


# ------------------------------------------------------------------ finding 12

def with_arm_limitation(name, side, kind='greedy_losses', knob='branching.max_label_branches'):
    """A two-arm fixture with a label-class limitation stated for one arm only, written
    and read back (the O record says lower_bound, as the rule requires)."""
    g = T.graphlet(name)
    g.limitations.append(Limitation(side, kind, knob, 1, 1, None, [], 'losses may be '
                                    'overestimated'))
    g.outcome.label_evidence = 'lower_bound'
    out = parse(g.dump(envelope=False))
    out.envelope, out.seed_summary = g.envelope, g.seed_summary
    return out


class TestFinding12ExactIsPerArm(unittest.TestCase):
    """evidence_block's exact was the seed's (outcome.label_evidence), so an arm read with
    the other arm's limitation, and with none of its own. It is per arm now (and so is
    Claim.exact); the outcome beside it stays the seed's."""

    def test_one_arm_limited(self):
        g = with_arm_limitation('fork', 'left')
        self.assertEqual([], g.check_rules())
        self.assertFalse(ops.evidence_block(g, 'left')['exact'])
        self.assertTrue(ops.evidence_block(g, 'right')['exact'])
        self.assertFalse(ops.evidence_block(g)['exact'])
        self.assertEqual('lower_bound', ops.evidence_block(g, 'right')['outcome']['label_evidence'])
        self.assertEqual({False}, {c.exact for c in g.claims('left')})
        self.assertEqual({True}, {c.exact for c in g.claims('right')})

    def test_a_seed_limitation_applies_to_both(self):
        g = with_arm_limitation('fork', None, 'seed_labels', 'labels.max_seed_labels')
        self.assertFalse(ops.evidence_block(g, 'left')['exact'])
        self.assertFalse(ops.evidence_block(g, 'right')['exact'])

    def test_an_unexplained_outcome_is_exact_nowhere(self):
        # a body whose K record was lost: the O record still says lower_bound
        g = T.graphlet('fork')
        g.outcome.label_evidence = 'lower_bound'
        for side in (None, 'left', 'right'):
            self.assertFalse(ops.evidence_block(g, side)['exact'])

    def test_cut_lists_are_per_arm(self):
        g = T.graphlet('annotate')             # right arm only, lists cut at 2
        self.assertFalse(ops.evidence_block(g, 'right')['exact'])
        self.assertEqual({False}, {c.exact for c in g.claims('right', strict=False)})


# ------------------------------------------------------------------ finding 13

class TestFinding13Splits(unittest.TestCase):
    """Graphlet.splits(): the split records of an arm, oriented like the other accessors
    (outward at_bp, labels as {name, ref}, each branch's first base and first walk)."""

    def test_merge(self):
        g = T.graphlet('merge')
        sps = g.splits('right')
        self.assertTrue(all(isinstance(sp, SplitPoint) for sp in sps))
        self.assertEqual([(20, 0, 'ambiguous', 4, [('A', 3), ('T', 2)]),
                          (46, 3, 'ambiguous', 4, [('C', 2), ('G', 3)])],
                         [(sp.at_bp, sp.segment, sp.kind, sp.labels_before,
                           [(b.char, b.labels_total) for b in sp.branches]) for sp in sps])
        for sp in sps:
            self.assertEqual('right', sp.arm)
            self.assertTrue(sp.ambiguous)
            for b in sp.branches:
                self.assertIsInstance(b, Branch)
                self.assertTrue(all(isinstance(l, Label) for l in b.labels))
        self.assertEqual(derive.splits(g.arms['right'], g.mode),
                         [type(derive.splits(g.arms['right'], g.mode)[0])(
                             sp.at_bp, sp.segment, [b.segment for b in sp.branches],
                             sp.ambiguous, sp.labels_before) for sp in sps])

    def test_orientation_on_both_arms(self):
        for name in ('fork', 'linear', 'ambiguous_split', 'quorum'):
            g = T.graphlet(name)
            for side in g.arms:
                for sp in g.splits(side):
                    for b in sp.branches:
                        if b.walk is None or b.char is None:
                            continue
                        walking = g.spell(side, b.walk, orientation='walk')
                        natural = g.spell(side, b.walk)
                        with self.subTest(doc=name, arm=side, at=sp.at_bp):
                            self.assertEqual(b.char, walking[sp.at_bp])
                            # natural: the left flank reads inward, the same base
                            self.assertEqual(b.char, natural[sp.at_bp] if side == 'right'
                                             else natural[len(natural) - 1 - sp.at_bp])
                            self.assertIn(b.segment,
                                          derive.paths(g.arms[side])[b.walk].segments)

    def test_filter_and_labels(self):
        g = T.graphlet('fork')
        for side in g.arms:
            every = g.splits(side)
            self.assertEqual([sp for sp in every if sp.labels_before >= 2],
                             g.splits(side, min_labels_before=2))
            for sp in every:
                for b in sp.branches:
                    child = g.arms[side].segments[b.segment]
                    self.assertEqual([g.labels[l].ref for l in child.entry],
                                     [l.ref for l in b.labels])


# ------------------------------------------------------------------ finding 14

class TestFinding14PackageMetadata(unittest.TestCase):
    def test_python_requires(self):
        path = os.path.join(T.API, 'setup.py')
        with open(path) as f:
            tree = ast.parse(f.read())
        kw = {k.arg: k.value for n in ast.walk(tree) if isinstance(n, ast.Call)
              for k in n.keywords}
        self.assertEqual('>=3.10', ast.literal_eval(kw['python_requires']))
        classifiers = ast.literal_eval(kw['classifiers'])
        self.assertNotIn('Programming Language :: Python :: 3.6', classifiers)
        self.assertIn('Programming Language :: Python :: 3.10', classifiers)
        self.assertGreaterEqual(sys.version_info[:2], (3, 10))


if __name__ == '__main__':
    unittest.main()
