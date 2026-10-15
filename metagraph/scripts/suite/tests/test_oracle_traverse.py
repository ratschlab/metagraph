"""Offline tests of the lane O /traverse oracle (ORACLE-DESIGN.md section 10): the walk rule and the
reference walker on the traversal fixture's cases, whose outcomes are known analytically (several
are the SPEC's own statements); the two independent computations of section 6.9's E against each
other; and the mutation triplets on recorded engine answers, offline (their build:
tests/data/oracle_traverse_answers.json, field engine)."""

import copy
import json
import os
import random
import sys
import unittest

SUITE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, SUITE)

from oracle import compare as C  # noqa: E402
from oracle import fixtures  # noqa: E402
from oracle import mutations as M  # noqa: E402
from oracle import oracle as O  # noqa: E402
from oracle import traverse_compare as TC  # noqa: E402
from oracle import traverse_fixtures as TF  # noqa: E402
from oracle import traverse_oracle as T  # noqa: E402

DATA = os.path.join(SUITE, 'tests', 'data', 'oracle_traverse_answers.json')
COLUMNS, CASES = TF.build()
_WALKERS = {}


def walker(mode='basic', fixture='trav', k=11):
    key = (mode, fixture, k)
    if key not in _WALKERS:
        cols = COLUMNS if fixture == 'trav' else fixtures.build_records()
        recs = [r for c in cols for r in cols[c]]
        _WALKERS[key] = T.Walker(O.Index(recs, k, mode))
    return _WALKERS[key]


def config_walker(config):
    """The walker of a run config's index (indexes.CONFIGS: mode, k, fixture set)."""
    from oracle import indexes
    mode, k, _, _, fset = indexes.CONFIGS[config]
    return walker(mode, fset, k)


def run(case, strategy, arm='right', mode='basic'):
    w = walker(mode)
    seed, labels, extra = CASES[case]
    if mode != 'basic':
        labels = list(dict.fromkeys(next(r.column for r in w.index.records if r.name == n)
                                    if n not in w.index.columns else n for n in labels))
        extra = list(dict.fromkeys(next(r.column for r in w.index.records if r.name == n)
                                   for n in extra))
    st = T.Strategy(strategy)
    labs = [w.label(n) for n in labels]
    return T.constrain_walk(w, seed, arm, labs, st, [w.label(n) for n in extra])


def ends(res, name=None):
    return {(len(c[0]), c[1], c[2], c[4]) for c in res.claims if name is None or c[1] == name}


class Rules(unittest.TestCase):
    """Each rule permits one case and forbids another (GPT's P1 traversal row)."""

    def test_divergence_and_the_split_limit(self):
        res = run('divergence', {'branching': {'on_reconverge': 'keep'}, 'bounds':
                                 {'max_extension_bp': 45}})
        self.assertEqual({(30, 'div1', 'dead_end', 0.0), (30, 'div2', 'dead_end', 0.0)}, ends(res))
        res = run('divergence', {'branching': {'on_reconverge': 'keep', 'max_splits_per_path': 0},
                                 'bounds': {'max_extension_bp': 45}})
        self.assertEqual({(0, 'div1', 'split_limit', 0.0), (0, 'div2', 'split_limit', 0.0)},
                         ends(res))

    def test_the_branch_allowance(self):
        st = {'branching': {'on_reconverge': 'keep', 'max_label_branches': 0},
              'bounds': {'max_extension_bp': 45}}
        self.assertEqual({(0, 't_amb', 'branch', 0.0)}, ends(run('ambiguous', st)))
        st['branching']['max_label_branches'] = 1
        res = run('ambiguous', st)
        self.assertEqual({(30, 't_amb', 'dead_end', 0.0)}, {e for e in ends(res)})
        self.assertEqual({1}, {c[3] for c in res.claims})            # branches + 1 on both
        self.assertEqual(2, len(res.claims))

    def test_a_cycle_is_unrolled_up_to_edge_reuse(self):
        st = {'branching': {'on_reconverge': 'keep', 'max_label_branches': 0},
              'bounds': {'max_extension_bp': 60}}
        self.assertEqual({'branch'}, {c[2] for c in run('cycle', st).claims})
        res = run('cycle', {'exhaustive': True, 'bounds': {'max_extension_bp': 60}})
        self.assertIn('edge_reuse', {c[2] for c in res.claims})

    def test_label_transitions_under_a_constant_cost(self):
        base = {'branching': {'on_reconverge': 'keep'}, 'bounds': {'max_extension_bp': 110}}
        self.assertEqual({(40, 'sw1', 'label_lost', 0.0)}, ends(run('transition', base)))
        one = dict(base, labels={'change_cost': {'model': 'constant', 'value': 1},
                                 'loss_budget': 1})
        res = run('transition', one)
        # sw1 switches into sw2 at 40 (a switched source ends silently, label_lost served);
        # sw2 is lost at 70, where entering sw3 needs a loss of 2 > 1
        self.assertEqual({(40, 'sw1', 'switched', 0.0), (70, 'sw2', 'loss_budget', 1.0)},
                         ends(res))
        self.assertEqual([2.0], res.needed)
        two = dict(base, labels={'change_cost': {'model': 'constant', 'value': 1},
                                 'loss_budget': 2})
        self.assertIn((100, 'sw3', 'dead_end', 2.0), ends(run('transition', two)))

    def test_a_switch_into_two_successors_is_an_ambiguous_branch(self):
        st = {'branching': {'on_reconverge': 'keep', 'max_label_branches': 0},
              'labels': {'change_cost': {'model': 'constant', 'value': 1}, 'loss_budget': 2},
              'bounds': {'max_extension_bp': 45}}
        self.assertIn((40, 'fa', 'branch', 0.0), ends(run('fork', st)))
        st['branching']['max_label_branches'] = 1
        res = run('fork', st)
        self.assertEqual({'fb', 'fc'}, {c[1] for c in res.claims if c[6] == 'switch'})

    def test_the_quorum_and_min_live_labels(self):
        base = {'branching': {'on_reconverge': 'keep'}, 'bounds': {'max_extension_bp': 45}}
        self.assertNotIn('minority', {c[2] for c in run('quorum', base).claims})
        for knob in ({'min_successor_labels': 2}, {'min_successor_fraction': 0.5}):
            st = dict(base, branching=dict(base['branching'], **knob))
            self.assertIn((0, 'qd', 'minority', 0.0), ends(run('quorum', st)))
        st = dict(base, bounds={'max_extension_bp': 45, 'min_live_labels': 2})
        self.assertIn((0, 'qd', 'below_min_labels', 0.0), ends(run('quorum', st)))

    def test_a_circle_re_enters_the_seed(self):
        res = run('circle', {'exhaustive': True, 'bounds': {'max_extension_bp': 60}})
        self.assertEqual({(40, 'circ', 'rejoined_seed', 0.0)}, ends(res))
        self.assertEqual({'rejoined_seed'}, {e[3] for e in res.events})

    def test_hairpins_on_canonical_and_primary(self):
        for mode in ('canonical', 'primary'):
            res = run('hairpin', {'exhaustive': True, 'bounds': {'max_extension_bp': 45}},
                      mode=mode)
            self.assertEqual({'hairpin'}, {c[2] for c in res.claims})
            # served as dead_end with the qualifier (SPEC 13, implementation note I3-2)
            self.assertEqual({'dead_end'}, {c[2] for c in TC.served_claims(res.claims)})
            res = run('hairpin', {'exhaustive': True, 'branching': {'hairpins': 'follow'},
                                  'bounds': {'max_extension_bp': 45}}, mode=mode)
            # followed, the walk turns back on itself: the reverse use of an edge it took is
            # blocked (edge_reuse_rc), and the label goes on to the record's end
            self.assertIn('followed', {e[3] for e in res.events})
            self.assertIn(('blocked', 'edge_reuse_rc'), {(e[2], e[3]) for e in res.events})
        # BASIC holds one strand: no hairpin
        self.assertNotIn('hairpin', {c[2] for c in run('hairpin', {'exhaustive': True}).claims})

    def test_a_followed_hairpin_never_counts_toward_ambiguity(self):
        """SPEC 13, note 11 (final breadth run seed 4, T-primary_k10-0036): at the palindromic
        10-mer ATATATATAT both followed hairpin steps carry c_zeta; under the allowance 0 the
        label goes on along both, with branches 0, instead of ending with branch."""
        w = walker('primary', 'records', 10)
        seed = 'CGGCGGCGGCGCACTAGAAGTAAATCTGCCTC'
        st = T.Strategy({'branching': {'hairpins': 'follow', 'on_reconverge': 'keep'},
                         'bounds': {'max_extension_bp': 40}})
        res = T.constrain_walk(w, seed, 'left', [w.label('c_zeta')], st)
        self.assertEqual({('GTGCAGCTTTTACACTCTGATTATATATATAGAATGGGGT', 'max_extension_bp', 0),
                          ('GTGCAGCTTTTACACTCTGATTATATATATAT', 'edge_reuse_rc', 0)},
                         {(c[0], c[2], c[3]) for c in res.claims})

    def test_an_ambiguous_label_counts_the_branch_on_its_hairpin_step_too(self):
        """Final breadth run seed 5 (T-primary_k11-0042): after TA the label holds a followed
        hairpin (A) and two plain successors (G, T); it is ambiguous over G and T, and every
        continuation, the hairpin's included, has branches 1."""
        w = walker('primary', 'records', 11)
        st = T.Strategy({'exhaustive': True, 'branching': {'hairpins': 'follow'},
                         'bounds': {'max_extension_bp': 40}})
        res = T.constrain_walk(w, 'ATATATATTAGTCTCACA', 'left', [w.label('c_zeta')], st)
        self.assertIn(('TAA', 'edge_reuse_rc', 1), {(c[0], c[2], c[3]) for c in res.claims})

    def test_trace_follows_one_coordinate_chain(self):
        kmer = run('trace', {'branching': {'on_reconverge': 'keep'},
                             'bounds': {'max_extension_bp': 45}})
        self.assertEqual({(0, 'tr', 'branch', 0.0)}, ends(kmer))
        trace = run('trace', {'branching': {'on_reconverge': 'keep'}, 'support': 'trace',
                              'bounds': {'max_extension_bp': 45}})
        self.assertNotIn('branch', {c[2] for c in trace.claims})
        tb = run('trace_break', {'exhaustive': True, 'support': 'trace',
                                 'bounds': {'max_extension_bp': 45}})
        self.assertEqual({(0, 'tb', 'trace_break', 0.0)}, ends(tb))

    def test_the_homopolymer_consequence_of_the_spec(self):
        """SPEC 6.6: for a record X A^(k+3) Z the trie holds A^(k-1) Z, A^k Z and A^(k+1) Z and
        never the record's own A^(k+3) Z."""
        w = walker()
        seed = CASES['homopolymer'][0]
        z0 = next(r for r in w.index.records if r.name == 'hom').seq[30 + 14]
        trie = T.annotate_trie(w, seed, 'right', 30, [])
        runs = {len(x) - len(x.lstrip('A')) for x in trie.nodes if x.lstrip('A')[:1] == z0}
        self.assertEqual({10, 11, 12}, runs)


class SpecContract(unittest.TestCase):
    """SPEC 6.9: under exhaustive forbid the walker's claims are E, the union of single-label
    maximal walks; two computations that share no label logic must agree."""

    def check(self, w, seed, labels, arm, radius):
        st = T.Strategy({'exhaustive': True, 'bounds': {'max_extension_bp': radius}})
        res = T.constrain_walk(w, seed, arm, labels, st)
        e = T.label_walks(w, seed, arm, labels, radius)
        self.assertEqual(e, {(c[0], c[1], c[2], c[3]) for c in res.claims}, (seed, arm))

    def test_on_every_case_of_the_traversal_fixture(self):
        for mode in ('basic', 'canonical', 'primary'):
            w = walker(mode)
            for case, (seed, names, _) in CASES.items():
                labels = [l for l in (w.label(n) for n in names) if l is not None]
                if mode != 'basic':
                    labels = [w.label(r.column) for r in w.index.records if r.name in names]
                    labels = list({l.name: l for l in labels}.values())
                labels = [l for l in labels if T.supports_seed(w, l, seed, 'kmer')]
                for arm in ('left', 'right'):
                    self.check(w, seed, labels, arm, 50)

    def test_on_random_seeds_of_the_phase_one_records(self):
        rng = random.Random(7)
        for mode, k in (('basic', 11), ('canonical', 11), ('canonical', 10), ('primary', 10)):
            w = walker(mode, 'records', k)
            for _ in range(25):
                r = rng.choice([r for r in w.index.indexed if 'N' not in r.seq])
                s = rng.randrange(0, len(r.seq) - k)
                seed = r.seq[s:s + rng.randint(k, 25)]
                kind = 'header' if mode == 'basic' else 'column'
                labels = [l for l in w.labels_of_kind(kind)
                          if T.supports_seed(w, l, seed, 'kmer')]
                for arm in ('left', 'right'):
                    self.check(w, seed, labels, arm, 20)


class RecordedTriplets(unittest.TestCase):
    """Members 1-3 on recorded engine answers: each answer passes (or, truncated or merged,
    is inconclusive on completeness with its invariants passing), and each corruption fails the
    check meant to catch it."""

    @classmethod
    def setUpClass(cls):
        with open(DATA) as f:
            cls.answers = json.load(f)['answers']

    def walker_for(self, config):
        return config_walker(config)

    def test_every_recorded_answer_and_its_corruptions(self):
        kinds = {'complete': 0, 'incomplete': 0, 'mutations': 0}
        for item in self.answers:
            w = self.walker_for(item['config'])
            req, ans = item['request'], item['answer']
            results = TC.compare_traverse(w, req, 200, ans)
            for pos, checks in results:
                self.assertNotEqual(C.FAIL, checks.outcome(), (item['id'], checks.failed()))
                kinds['complete' if checks.outcome() == C.PASS else 'incomplete'] += 1
                for name, fam, mutated in M.mutations_for_traverse(req, ans, pos, w):
                    mc = dict(TC.compare_traverse(w, req, 200, mutated))
                    self.assertTrue(M.caught(mc[pos], fam), (item['id'], name))
                    kinds['mutations'] += 1
        self.assertGreater(kinds['complete'], 30)
        self.assertGreater(kinds['incomplete'], 3)
        self.assertGreater(kinds['mutations'], 300)

    def test_a_truncated_answer_claiming_completeness_fails(self):
        item = next(i for i in self.answers
                    if any(a.get('status') == 'truncated'
                           for r in i['answer']['results'] for a in r.get('arms', {}).values()))
        w = self.walker_for(item['config'])
        ans = json.loads(json.dumps(item['answer']))
        radius = item['request']['strategy']['bounds']['max_extension_bp']
        for a in ans['results'][0]['arms'].values():
            a.update({'status': 'complete', 'complete_to_bp': radius})
        (pos, checks), = TC.compare_traverse(w, item['request'], 200, ans)
        self.assertEqual(C.FAIL, checks.outcome())


if __name__ == '__main__':
    unittest.main()


class VerifierInjections(unittest.TestCase):
    """The verification's injected engine bugs (verify/scratch/inject.py of 2026-10-15), applied to
    the recorded answers: each must be flagged. Before the fixes they named, each was missed."""

    @classmethod
    def setUpClass(cls):
        with open(DATA) as f:
            cls.answers = json.load(f)['answers']

    def flagged(self, inject, min_applicable=3):
        applicable = 0
        for item in self.answers:
            w = config_walker(item['config'])
            base = TC.compare_traverse(w, item['request'], 200, item['answer'])
            if any(c.outcome() == C.FAIL for _, c in base):
                continue
            mutated = inject(item)
            if mutated is None:
                continue
            applicable += 1
            out = TC.compare_traverse(w, item['request'], 200, mutated)
            self.assertTrue(any(c.outcome() == C.FAIL for _, c in out), item['id'])
        self.assertGreaterEqual(applicable, min_applicable)

    @staticmethod
    def complete_keep_arms(answer):
        for res in answer['results']:
            for aj in (res.get('arms') or {}).values():
                if aj.get('status') == 'complete' and aj.get('completeness_scope') == 'per_path':
                    yield aj

    def test_g11_a_run_listed_twice(self):
        def inject(item):
            if (item['request']['strategy'].get('labels') or {}).get('mode') == 'annotate':
                return None
            a = copy.deepcopy(item['answer'])
            for aj in self.complete_keep_arms(a):
                if aj.get('runs'):
                    aj['runs'].append(copy.deepcopy(aj['runs'][-1]))
                    return a
            return None
        self.flagged(inject)

    def test_g12_a_leaf_walk_listed_twice(self):
        def inject(item):
            a = copy.deepcopy(item['answer'])
            for aj in self.complete_keep_arms(a):
                segs, paths = aj['segments'], aj['paths']
                cand = [p for p in paths if len(p['segments']) >= 2]
                if not cand:
                    continue
                p = cand[-1]
                last = next(s for s in segs if s['id'] == p['segments'][-1])
                clone = copy.deepcopy(last)
                clone['id'] = max(s['id'] for s in segs) + 1
                clone['events'] = [e for e in clone['events'] if e['type'] != 'label_end']
                segs.append(clone)
                pc = copy.deepcopy(p)
                pc['id'] = max(x['id'] for x in paths) + 1
                pc['segments'] = p['segments'][:-1] + [clone['id']]
                pc['end_labels'], pc['end_reasons'] = [], {}
                paths.append(pc)
                return a
            return None
        self.flagged(inject)

    def test_g12b_a_successor_enumerated_twice(self):
        """A dummy edge beside the real one with the same label: the leaf segment cloned whole,
        with its runs, events and path."""
        def inject(item):
            a = copy.deepcopy(item['answer'])
            for aj in self.complete_keep_arms(a):
                segs, paths, runs = aj['segments'], aj['paths'], aj.get('runs') or []
                cand = [p for p in paths if len(p['segments']) >= 2]
                if not cand:
                    continue
                p = cand[-1]
                last = next(s for s in segs if s['id'] == p['segments'][-1])
                clone = copy.deepcopy(last)
                clone['id'] = max(s['id'] for s in segs) + 1
                segs.append(clone)
                remap = {}
                for i, r in enumerate(list(runs)):
                    if r.get('segment') == last['id']:
                        rc = copy.deepcopy(r)
                        rc['segment'] = clone['id']
                        remap[i] = len(runs)
                        runs.append(rc)
                pc = copy.deepcopy(p)
                pc['id'] = max(x['id'] for x in paths) + 1
                pc['segments'] = p['segments'][:-1] + [clone['id']]
                for e in pc.get('end_labels') or []:
                    if e.get('run') in remap:
                        e['run'] = remap[e['run']]
                paths.append(pc)
                return a
            return None
        self.flagged(inject)


class ConstantCostReasons(unittest.TestCase):
    """Verification H2: the verifier's probes T1-T4 (constant cost 1, keep) and their 81c620ef
    answers. SPEC 6.7 with 13 I3 note 3: a label with no target at any admissible successor is
    label_lost (T2, T4's sw2 at 30); an over-budget switch into a target a cheaper predecessor
    carries is a plain label_lost (T3); a within-budget one is superseded (T4); extras
    unreachable within the budget are refused (T1, SPEC 5)."""

    @classmethod
    def setUpClass(cls):
        with open(DATA) as f:
            cls.probes = {p['id']: p for p in json.load(f)['probes']}

    def check(self, name):
        p = self.probes[name]
        out = TC.compare_traverse(walker('basic'), p['request'], p['status'], p['answer'])
        self.assertEqual([], [c for _, ch in out for c in ch.failed()], name)
        return out

    def test_t1_unreachable_extras_are_refused(self):
        self.check('T1_no_target_over_budget')

    def test_t2_no_target_within_budget_is_label_lost(self):
        self.check('T2_no_target_within_budget')

    def test_t3_an_over_budget_switch_into_a_carried_target_is_label_lost(self):
        self.check('T3_target_carried_over_budget')

    def test_t4_a_within_budget_switch_into_a_carried_target_is_superseded(self):
        self.check('T4_target_carried_within_budget')

    def test_the_re_minimisation_keeps_the_stay_and_switch_on(self):
        """SPEC 6.3, 6.4: after a branch exclusion a removed target is re-minimised over the
        remaining sources, its own stay included, and under switch_on loss never from a source
        that continues at the successor. The verification found the loop ignoring both."""
        w = walker('basic')
        lab = {n: w.label(n) for n in ('sw1', 'sw2', 'sw3')}
        by_name = {n: l for n, l in lab.items()}
        # v: a k-mer of sw2's W2 part, which sw2 holds and sw1 does not
        rec = next(r for r in w.index.records if r.name == 'sw2')
        v = rec.seq[-11:]
        self.assertTrue(w.holds(lab['sw2'], v) and not w.holds(lab['sw1'], v))
        sigma = {'sw1': T.Entry(0.0, 0, 'sw1', ('sw1', 0, 'seed')),
                 'sw2': T.Entry(1.0, 0, 'sw2', ('sw2', 0, 'switch'))}
        st = T.Strategy({'labels': {'change_cost': {'model': 'constant', 'value': 0.5},
                                    'loss_budget': 2}})
        best = T.derive_entry(w, st, by_name, sigma, lab['sw2'], v, 5)
        self.assertEqual(('sw1', 0.5), (best.pred, best.loss))       # the switch beats the stay
        again = T.derive_entry(w, st, by_name, sigma, lab['sw2'], v, 5, exclude={'sw1'})
        self.assertEqual(('sw2', 1.0), (again.pred, again.loss))     # the stay, once sw1 is out
        # a source that continues at v is no source under switch_on loss
        v3 = next(r for r in w.index.records if r.name == 'sw3').seq[:11]
        self.assertTrue(w.holds(lab['sw2'], v3) and w.holds(lab['sw3'], v3))
        sigma = {'sw2': T.Entry(0.0, 0, 'sw2', ('sw2', 0, 'seed'))}
        self.assertIsNone(T.derive_entry(w, st, by_name, sigma, lab['sw3'], v3, 5))
        anyst = T.Strategy({'labels': {'change_cost': {'model': 'constant', 'value': 0.5},
                                       'loss_budget': 2, 'switch_on': 'any'}})
        self.assertEqual('sw2', T.derive_entry(w, anyst, by_name, sigma, lab['sw3'], v3, 5).pred)


class MergeArms(unittest.TestCase):
    """Verification H5: a merge arm that merged nothing is held strictly to the keep walker; one
    that merged keeps its invariants and states a reconverge event at each merge."""

    @classmethod
    def setUpClass(cls):
        with open(DATA) as f:
            cls.answers = json.load(f)['answers']

    def arms(self, merged):
        for item in self.answers:
            for res in item['answer']['results']:
                for arm, aj in (res.get('arms') or {}).items():
                    if aj.get('completeness_scope') != 'united_history':
                        continue
                    did = any(len(s.get('parents') or []) > 1 for s in aj['segments'])
                    if did == merged and aj.get('status') == 'complete':
                        yield item, arm

    def test_g8_a_wrong_reason_in_a_merge_arm_that_merged_nothing_fails(self):
        n = 0
        for item, arm in self.arms(merged=False):
            a = copy.deepcopy(item['answer'])
            runs = [r for res in a['results'] for r in res['arms'][arm]['runs']
                    if r['end_reason'] in ('label_lost', 'branch', 'max_extension_bp')]
            if not runs:
                continue
            runs[0]['end_reason'] = 'dead_end'
            w = config_walker(item['config'])
            out = TC.compare_traverse(w, item['request'], 200, a)
            self.assertTrue(any(c.outcome() == C.FAIL for _, c in out), item['id'])
            n += 1
        self.assertGreater(n, 0)

    def test_a_merge_without_its_reconverge_event_fails(self):
        n = 0
        for item, arm in self.arms(merged=True):
            w = config_walker(item['config'])
            self.assertFalse(any(c.outcome() == C.FAIL for _, c in
                                 TC.compare_traverse(w, item['request'], 200, item['answer'])))
            a = copy.deepcopy(item['answer'])
            for res in a['results']:
                for s in res['arms'][arm]['segments']:
                    s['events'] = [e for e in s['events'] if e['type'] != 'reconverge']
            out = TC.compare_traverse(w, item['request'], 200, a)
            self.assertIn(f'traverse.{arm}.reconverge_events',
                          {x['check'] for _, c in out for x in c.failed()})
            n += 1
        self.assertGreater(n, 0)


MERGE_DATA = os.path.join(SUITE, 'tests', 'data', 'oracle_merge_answers.json')


class MergedArmAnswers(unittest.TestCase):
    """The four merged-arm answers of the verification round's final runs that the oracle failed;
    all four oracle errors (VERIFY-ORACLE.md Round 2, R2.1; their build: the data file's engine)."""

    @classmethod
    def setUpClass(cls):
        with open(MERGE_DATA) as f:
            cls.items = {x['id']: x for x in json.load(f)['answers']}

    def failed(self, item, answer=None):
        w = config_walker(item['config'])
        out = TC.compare_traverse(w, item['request'], item['status'],
                                  item['answer'] if answer is None else answer)
        return {x['check'] for _, c in out for x in c.failed()}

    def test_r2b_a_switched_in_run_is_supported_from_its_first_step(self):
        """T-canonical_k10-0076 right: runs [22, 23) switched in at 22 and closed by the merge at
        23 hold their label on node 23 only; node 22 is the source's (SPEC 6.2)."""
        item = self.items['T-canonical_k10-0076']
        self.assertNotIn('traverse.right.displayed_support', self.failed(item))
        # the check keeps its teeth: the run handed to a label that lacks the node it entered fails
        w = config_walker(item['config'])
        res = item['answer']['results'][0]
        view = TC.ArmView(res['arms']['right'], 'right',
                          {j: d['name'] for j, d in enumerate(res['label_dict'])})
        i, (claim, sw) = next((i, x) for i, x in enumerate(view.runs)
                              if x[1]['entered_by'] == 'switch')
        seed = item['request']['seeds'][0]['sequence']
        node = w.walk_node(seed, 'right', claim[0], sw['to_bp'])
        other = next(j for j, d in enumerate(res['label_dict'])
                     if not w.holds(w.label(d['name']), node))
        a = copy.deepcopy(item['answer'])
        a['results'][0]['arms']['right']['runs'][i]['label'] = other
        self.assertIn('traverse.right.displayed_support', self.failed(item, a))
