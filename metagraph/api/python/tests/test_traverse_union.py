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
    test.assertEqual(cm.exception.code, code, str(cm.exception))
    if text is not None:
        test.assertEqual(str(cm.exception), text)
    return str(cm.exception)


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
            U.TEXTS['different_seed_length'] % ('the seed of the absent record of pair g2',
                                                31, "part a's", 30, 0))
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
        size = len(U.C.canonical_json({'statements': u.summary(max_bytes=10 ** 6)['statements']}))
        refused(self, 'summary_too_small', lambda: u.summary(max_bytes=100),
                U.TEXTS['summary_too_small'] % (size, 100, size))
        refused(self, 'not_offered', u.dump, U.TEXTS['not_offered_dump'])
        refused(self, 'not_offered', u.to_json, U.TEXTS['not_offered_json'])
        refused(self, 'not_offered', u.compare, U.TEXTS['not_offered'] % 'compare')
        refused(self, 'not_offered', lambda: u.walks('right', paths='inferred'),
                U.TEXTS['not_offered_inferred'])
        refused(self, 'not_offered', lambda: u.walks('right', max_transitions=2))
        refused(self, 'not_offered', lambda: u.spell('right', 'i:1.2'))


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

    def test_seed_check_ignores_a_matching_seed_id(self):
        seq = T.graphlet('merge').seed.sequence
        r = absent_result('not_in_graph')
        r['seed']['seed_id'] = seq          # a caller's text proves nothing
        u = union_pairs([{'pair': 'g2', 'result': r}], seed=seq)
        self.assertEqual(u.absent[0].seed_check, 'unchecked')


class TestOrder(unittest.TestCase):
    """P3: any order of the results gives the same digest and byte-equal answers."""

    def test_permutations(self):
        seq = T.graphlet('merge').seed.sequence
        items = stated(gA='annotate', gM='merge', gT='merge_ties') + [
            {'pair': 'g2', 'result': absent_result('not_in_graph')},
            {'pair': 'g5', 'error': 'HTTP 503'}]
        ref = union_pairs(items, expected_pairs=['gA', 'gM', 'gT', 'g2', 'g5', 'g9'])
        want = (ref.digest, str(ref), ref.to_fasta(stubs=True), ref.to_gfa(),
                json.dumps(ref.summary(), sort_keys=True), repr(ref.walks('right')),
                str(ref.carriers('right', 1)))
        rng = random.Random(7)
        for _ in range(6):
            perm = items[:]
            rng.shuffle(perm)
            u = union_pairs(perm, expected_pairs=['g9', 'g5', 'g2', 'gT', 'gM', 'gA'])
            self.assertEqual((u.digest, str(u), u.to_fasta(stubs=True), u.to_gfa(),
                              json.dumps(u.summary(), sort_keys=True), repr(u.walks('right')),
                              str(u.carriers('right', 1))), want)
        self.assertEqual(seq, ref.seed.sequence)


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


R237 = os.path.join(UNION, 'r237')
# the body sha256 of each part as the service's record states it (the README of r237/)
R237_BODIES = {'bact_other': 'cda902d0df71afab3d9fe1bcf33a04630468c22a1748777941686867fcc31eea',
               'bact_pseudomonadota':
                   '6d48f657aa3af50bd79c3ed39126d42c4c8c1e1bc431f8f9d6d5155cc732d522'}


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

    def test_a_high_ratio_gzip_stops_at_the_cap(self):
        import gzip
        from metagraph.traverse import load
        from metagraph.traverse import parser as P
        text = T.doc_text('merge')
        # 64 MiB of one byte in a comment-free document: about 64 KiB of gzip
        bomb = self.write('bomb.mgt.gz', gzip.compress(text.encode('utf-8') + b'A' * (64 << 20)))
        self.assertLess(os.path.getsize(bomb), 1 << 20)
        old = P.GZIP_TEXT_MAX
        P.GZIP_TEXT_MAX = 1 << 20
        try:
            with self.assertRaises(GraphletFormatError) as cm:
                load(bomb)
        finally:
            P.GZIP_TEXT_MAX = old
        self.assertEqual(str(cm.exception),
                         'the gzip file bomb.mgt.gz inflates to more than 1048576 bytes, the most '
                         'a load() without budget= inflates (GZIP_TEXT_MAX): refused, not '
                         'inflated further; pass budget= to bound it by its memory account '
                         'instead')
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
