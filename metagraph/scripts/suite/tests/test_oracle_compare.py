"""Offline tests of the lane O comparison (oracle/compare.py) by mutation triplets (ORACLE-DESIGN.md
section 6): for each assertion a valid answer passes, a one-property corruption fails the check
meant to catch it, and a legitimately incomplete answer is inconclusive on completeness while
its observable invariants still pass, and its own corruptions fail. The answers are rendered
from the oracle (oracle/render.py); no server."""

import copy
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oracle import compare as C  # noqa: E402
from oracle import fixtures, mutations as M  # noqa: E402
from oracle import oracle as O  # noqa: E402
from oracle import render as R  # noqa: E402
from oracle import requests_gen  # noqa: E402

CONFIGS = {
    'basic': C.Config('basic', 'basic', 11),
    'canonical': C.Config('canonical', 'canonical', 11, annotation='unbudgeted'),
    'primary': C.Config('primary', 'primary', 11, annotation='unbudgeted'),
}
_INDEXES = {}


def fixture_index(mode, k=11):
    if (mode, k) not in _INDEXES:
        _INDEXES[(mode, k)] = O.Index(fixtures.all_records(fixtures.build_records()), k, mode)
    return _INDEXES[(mode, k)]


def outcome(results):
    outs = [c.outcome() for _, c, _ in results]
    return C.FAIL if C.FAIL in outs else (C.INCONCLUSIVE if C.INCONCLUSIVE in outs else C.PASS)


class Triplets(unittest.TestCase):
    """Members 1 and 2 on every generated request of a fixed seed, in every graph mode."""

    def run_mode(self, mode, n=14):
        config = CONFIGS[mode]
        index = fixture_index(mode)
        reqs = requests_gen.generate(index, f'test_{mode}', 3, labels_ok=True,
                                     resolve_ok=mode == 'basic', header_labels=mode == 'basic',
                                     n_pattern=n, n_budget=0, n_resolve=6)
        applied = 0
        for item in reqs:
            req = item['request']
            if item['route'] == 'resolve':
                try:
                    answer = R.render_resolve(index, req)
                except ValueError:
                    continue
                checks, _ = C.compare_resolve(index, config, req, answer)
                self.assertEqual(C.PASS, checks.outcome(), (req, checks.failed()))
                for name, fam, mutated in M.mutations_for_resolve(answer):
                    mc, _ = C.compare_resolve(index, config, req, mutated)
                    self.assertTrue(M.caught(mc, fam), (name, req))
                    applied += 1
                continue
            mc_eff = min(req.get('max_contexts', 100000), 100000)
            answer = R.render_answer(index, config, req, mc_eff)
            results = C.compare_pattern(index, config, req, answer, mc_eff)
            for pos, checks, exp in results:
                self.assertEqual(C.PASS, checks.outcome(),
                                 (mode, req['patterns'][pos], checks.failed()))
                entry = answer['patterns'][pos]
                for name, fam, mutated in M.mutations_for_entry(index, config, exp, entry):
                    mc, _ = C.compare_entry(index, config, req, pos, mutated, None, mc_eff, exp)
                    self.assertTrue(M.caught(mc, fam),
                                    (mode, name, req['patterns'][pos], req.get('mode')))
                    applied += 1
        self.assertGreater(applied, 50)

    def test_basic(self):
        self.run_mode('basic')

    def test_canonical(self):
        self.run_mode('canonical')

    def test_primary(self):
        self.run_mode('primary')


class IncompleteAnswers(unittest.TestCase):
    """Member 3: a stated max_steps stop. Completeness is not judged; the invariants are."""

    def setUp(self):
        self.index = fixture_index('basic')
        self.config = CONFIGS['basic']

    def check(self, req, stop_at):
        answer = R.render_answer(self.index, self.config, req, 100000, stop_at)
        results = C.compare_pattern(self.index, self.config, req, answer, 100000)
        self.assertEqual(C.INCONCLUSIVE, outcome(results),
                         [c.failed() for _, c, _ in results])
        for pos, checks, exp in results:
            if pos < stop_at[0]:
                continue
            self.assertIn('completeness', [c['check'] for c in checks.items
                                           if c['outcome'] == C.INCONCLUSIVE])
            entry = answer['patterns'][pos]
            n = 0
            for name, fam, mutated in M.mutations_for_entry(self.index, self.config, exp, entry,
                                                             incomplete=True):
                mc, _ = C.compare_entry(self.index, self.config, req, pos, mutated,
                                        C.sticky_stop(answer['patterns'][stop_at[0]]['stop'], req)
                                        if pos > stop_at[0] else None, 100000, exp)
                self.assertTrue(M.caught(mc, fam), (name, pos))
                n += 1
            self.assertGreater(n, 0)
        return answer

    def test_partial_stopped_in_discovery(self):
        req = {'mode': 'partial', 'max_steps': 10 ** 6,
               'patterns': [{'iupac': 'GNNNA'}, {'dna': 'ACGTA'}, {'iupac': 'TTNNA'}]}
        answer = self.check(req, (1, 0))
        self.assertEqual({'reason': 'max_steps'}, answer['patterns'][1]['cut'])

    def test_all_or_count_stopped(self):
        req = {'mode': 'all_or_count', 'max_steps': 10 ** 6,
               'patterns': [{'iupac': 'RNNA'}, {'iupac': 'GNNC'}]}
        self.check(req, (0, 3))

    def test_count_with_mismatches_stopped(self):
        req = {'mode': 'count', 'max_steps': 10 ** 6, 'max_mismatches': 2,
               'patterns': [{'dna': 'ACCGTTAACC'}, {'dna': 'GAATTCAAG'}]}
        self.check(req, (0, 1))

    def test_a_stopped_answer_that_claims_completeness_fails(self):
        req = {'mode': 'partial', 'max_steps': 10 ** 6, 'patterns': [{'iupac': 'GNNNA'}]}
        answer = R.render_answer(self.index, self.config, req, 100000, (0, 2))
        answer['patterns'][0]['retrieval_complete'] = True
        results = C.compare_pattern(self.index, self.config, req, answer, 100000)
        self.assertEqual(C.FAIL, outcome(results))

    def test_a_later_pattern_must_carry_the_sticky_stop(self):
        req = {'mode': 'count', 'max_steps': 10 ** 6,
               'patterns': [{'iupac': 'GNNNA'}, {'dna': 'ACGTA'}]}
        answer = R.render_answer(self.index, self.config, req, 100000, (0, 2))
        answer['patterns'][1]['stop'] = None
        results = C.compare_pattern(self.index, self.config, req, answer, 100000)
        self.assertIn('stop.sticky', [c['check'] for c in results[1][1].failed()])


class Probes(unittest.TestCase):
    """The review's own probes, as regressions of the comparison."""

    def setUp(self):
        self.index = fixture_index('basic')
        self.config = CONFIGS['basic']

    def compare(self, req, answer):
        return C.compare_pattern(self.index, self.config, req, answer, 100000)

    def test_doubling_an_exact_count_and_all_its_partitions_fails(self):
        # GPT's P-exact-0003 probe: double the total and every offset/strand/suffix partition
        alpha = self.index.records[1].seq
        req = {'mode': 'count', 'patterns': [{'dna': alpha[40:47]}]}
        answer = R.render_answer(self.index, self.config, req)
        self.assertEqual(C.PASS, outcome(self.compare(req, answer)))
        doubled = copy.deepcopy(answer)
        for count in M._count_paths(doubled['patterns'][0]):
            M._double(count)
        failed = {c['check'].split('[')[0] for c in self.compare(req, doubled)[0][1].failed()}
        self.assertTrue({'counts.total', 'counts.suffix', 'counts.by_offset',
                         'counts.by_strand'} <= failed, failed)

    def test_a_duplicate_that_replaces_a_missing_result_fails(self):
        # review finding 5: one duplicated row and one dropped row keep the length
        req = {'mode': 'all_or_count', 'patterns': [{'iupac': 'GNNNA'}]}
        answer = R.render_answer(self.index, self.config, req)
        rs = answer['patterns'][0]['results']
        self.assertGreater(len(rs), 3)
        rs[2] = copy.deepcopy(rs[1])
        failed = {c['check'] for c in self.compare(req, answer)[0][1].failed()}
        self.assertIn('results.unique', failed)
        self.assertIn('results.complete', failed)

    def test_a_missing_by_offset_key_fails(self):
        req = {'mode': 'count', 'patterns': [{'dna': 'GAATTC'}]}
        answer = R.render_answer(self.index, self.config, req)
        del answer['patterns'][0]['counts']['contexts']['by_offset']['0']
        failed = {c['check'] for c in self.compare(req, answer)[0][1].failed()}
        self.assertIn('counts.keys.by_offset', failed)

    def test_an_unmasked_bounds_count_holding_the_truth_passes(self):
        config = C.Config('unmasked', 'basic', 11, masked=False)
        req = {'mode': 'count', 'patterns': [{'dna': 'GAATTC'}]}
        answer = R.render_answer(self.index, config, req)
        c = answer['patterns'][0]['counts']['contexts']
        truth = c['value']
        dummies = len(O.dummy_contexts(self.index, O.Pattern('dna', 'GAATTC')))
        # U: the masked count plus the source dummies holding the pattern (SPEC 7.4)
        c.update({'relation': 'bounds', 'lower': truth - 1, 'upper': truth + dummies,
                  'value': truth - 1, 'estimate': truth})
        res = C.compare_pattern(self.index, config, req, answer, 100000)
        self.assertEqual(C.PASS, outcome(res), res[0][1].failed())
        for wrong in ({'lower': truth + 1, 'value': truth + 1},       # the truth outside
                      {'upper': truth + dummies + 1},                 # an upper that is not U
                      {'lower': truth, 'value': truth, 'upper': 2 * (truth + dummies)}):
            bad = copy.deepcopy(answer)
            bad['patterns'][0]['counts']['contexts'].update(wrong)
            res = C.compare_pattern(self.index, config, req, bad, 100000)
            self.assertIn('counts.total', [x['check'] for x in res[0][1].failed()], wrong)
        # on a masked graph a completed discovery must be exact
        res = C.compare_pattern(self.index, self.config, req, answer, 100000)
        self.assertEqual(C.FAIL, outcome(res))

    def test_a_stop_leaves_the_complete_orientation_bounds_at_u(self):
        # breadth run seed 4 (B-basic_k11_unchecked-0140): a max_steps stop leaves the total
        # at_least and the orientation searched first complete; without the mask that one is
        # bounds, its upper U = its truth plus its own source dummies (SPEC 7.4, 7.9)
        config = C.Config('unmasked', 'basic', 11, masked=False)
        pat = O.Pattern('iupac', 'ACGKAC')
        req = {'mode': 'count', 'max_steps': 57, 'patterns': [{'iupac': 'ACGKAC'}]}
        answer = R.render_answer(self.index, config, req, 100000, (0, 3))
        c = answer['patterns'][0]['counts']['contexts']
        truth = len([x for x in O.contexts(self.index, pat, 'any_offset', 'forward')])
        dummies = len(O.dummy_contexts(self.index, pat, 'any_offset', 'forward'))
        self.assertGreater(dummies, 0)
        c['by_strand']['+'] = {'value': 2, 'relation': 'bounds', 'unit': 'graph_contexts',
                               'lower': 2, 'upper': truth + dummies, 'estimate': truth}
        res = C.compare_pattern(self.index, config, req, answer, 100000)
        self.assertEqual(C.INCONCLUSIVE, outcome(res), res[0][1].failed())
        c['by_strand']['+']['upper'] = truth
        res = C.compare_pattern(self.index, config, req, answer, 100000)
        self.assertIn('counts.by_strand[+]', [x['check'] for x in res[0][1].failed()])

    def test_partial_cut_lists_the_first_contexts_in_answer_order(self):
        req = {'mode': 'partial', 'max_contexts': 5, 'patterns': [{'iupac': 'GNNNA'}]}
        answer = R.render_answer(self.index, self.config, req, 5)
        self.assertEqual({'reason': 'max_contexts'}, answer['patterns'][0]['cut'])
        self.assertEqual(C.PASS, outcome(C.compare_pattern(self.index, self.config, req,
                                                           answer, 5)))
        full = R.render_answer(self.index, self.config, dict(req, max_contexts=100000))
        # the sixth context in place of the fifth: still five sound, sorted results
        answer['patterns'][0]['results'][4] = full['patterns'][0]['results'][5]
        failed = {c['check'] for c in C.compare_pattern(self.index, self.config, req, answer,
                                                        5)[0][1].failed()}
        self.assertEqual({'results.complete'}, failed)

    def test_a_primary_cut_must_not_skip_an_earlier_context_of_a_listed_kmer(self):
        # final run seed 1 (P-primary_k11-0030): AAAAAA within 3 substitutions, max_contexts 1;
        # a result moved to the next offset of its k-mer is still a context, but no longer the
        # first in answer order
        index, config = fixture_index('primary'), CONFIGS['primary']
        req = {'mode': 'partial', 'max_contexts': 1, 'max_mismatches': 3, 'strands': 'forward',
               'patterns': [{'dna': 'AAAAAA'}]}
        answer = R.render_answer(index, config, req, 1)
        res = C.compare_pattern(index, config, req, answer, 1)
        self.assertEqual(C.PASS, outcome(res), res[0][1].failed())
        r = answer['patterns'][0]['results'][0]
        exp = res[0][2]
        later = next(x for x in exp.ctxs if x.kmer == r['kmer'] and x.offset > r['offset']
                     and x.level == r['mismatches'])
        r.update({'offset': later.offset, 'instance': later.kmer[later.offset:later.offset + 6],
                  'substitutions': R._subs(exp.pattern, later.subs)})
        failed = {x['check'] for x in C.compare_pattern(index, config, req, answer,
                                                        1)[0][1].failed()}
        self.assertEqual({'results.complete'}, failed)

    def test_a_withheld_answer_on_u_claims_no_completeness(self):
        # final breadth runs (P-basic_k11_unchecked-0043 and others): all_or_count withheld on
        # U while the truth fits max_contexts must still say retrieval_complete false
        config = C.Config('unchecked', 'basic', 11, masked=False)
        req = {'mode': 'all_or_count', 'max_contexts': 20, 'patterns': [{'protein': 'TK'}]}
        answer = R.render_answer(self.index, config, req, 20)
        e = answer['patterns'][0]
        truth = e['counts']['contexts']['value']
        self.assertLessEqual(truth, 20)
        e['counts']['contexts'].update({'relation': 'bounds', 'lower': truth,
                                        'upper': 10 ** 6, 'estimate': truth})
        e.update({'withheld': {'reason': 'count_above_threshold'}, 'results': [],
                  'returned': 0, 'retrieval_complete': True,
                  'notes': e['notes'] + ['threshold_upper_bound']})
        failed = {x['check'] for x in C.compare_pattern(self.index, config, req, answer,
                                                        20)[0][1].failed()}
        self.assertIn('release.all_or_count', failed)

    def test_slots(self):
        req = {'mode': 'count', 'patterns': [{'dna': 'ACGU'}, {'protein': 'MOK'}, {'dna': ''}]}
        answer = R.render_answer(self.index, self.config, req)
        self.assertEqual(C.PASS, outcome(self.compare(req, answer)))
        answer['patterns'][0] = R.render_answer(self.index, self.config,
                                                {'patterns': [{'dna': 'ACGT'}]})['patterns'][0]
        self.assertEqual(C.FAIL, outcome(self.compare(req, answer)))
        primary = CONFIGS['primary']
        req = {'mode': 'count', 'scope': 'suffix', 'patterns': [{'dna': 'GAATTC'}]}
        answer = R.render_answer(fixture_index('primary'), primary, req)
        self.assertEqual('scope_unsupported', answer['patterns'][0]['error']['code'])
        res = C.compare_pattern(fixture_index('primary'), primary, req, answer, 100000)
        self.assertEqual(C.PASS, outcome(res))

    def test_a_resolve_label_of_a_record_shorter_than_k_is_refused(self):
        req = {'sequence': 'ACGTTGCAACGTTGCA', 'labels': ['delta3']}
        body = {'error': "Label not found in the annotation: 'delta3'"}
        checks, _ = C.compare_resolve(self.index, self.config, req, body, 400)
        self.assertEqual(C.PASS, checks.outcome())
        checks, _ = C.compare_resolve(self.index, self.config, req, {'labels': []}, 200)
        self.assertEqual(C.FAIL, checks.outcome())


if __name__ == '__main__':
    unittest.main()


class RowsAndTwins(unittest.TestCase):
    """Verification H4: rows on CANONICAL and PRIMARY (SPEC 7.10), and the PRIMARY partial cut
    held to its ample twin (ORACLE-DESIGN 4.5, verification M5)."""

    def test_g7_a_canonical_row_off_by_one_fails(self):
        index, config = fixture_index('canonical'), CONFIGS['canonical']
        req = {'mode': 'all_or_count', 'patterns': [{'iupac': 'GNNNA'}]}
        answer = R.render_answer(index, config, req)
        self.assertEqual(C.PASS, outcome(C.compare_pattern(index, config, req, answer)))
        for r in answer['patterns'][0]['results']:
            r['row'] += 1
        failed = {x['check'] for x in C.compare_pattern(index, config, req, answer)[0][1].failed()}
        self.assertIn('results.row', failed)

    def test_a_primary_row_of_the_wrong_pair_fails(self):
        index, config = fixture_index('primary'), CONFIGS['primary']
        req = {'mode': 'all_or_count', 'patterns': [{'iupac': 'GNNNA'}]}
        answer = R.render_answer(index, config, req)
        self.assertEqual(C.PASS, outcome(C.compare_pattern(index, config, req, answer)))
        rs = answer['patterns'][0]['results']
        rs[0]['row'] = rs[-1]['row'] if rs[-1]['row'] != rs[0]['row'] else rs[0]['row'] + 7
        failed = {x['check'] for x in C.compare_pattern(index, config, req, answer)[0][1].failed()}
        self.assertIn('results.row', failed)

    def test_g6_a_primary_cut_listing_a_later_kmer_fails_against_its_twin(self):
        index, config = fixture_index('primary'), CONFIGS['primary']
        req = {'mode': 'partial', 'max_contexts': 5, 'patterns': [{'iupac': 'GNNNA'}]}
        answer = R.render_answer(index, config, req, 5)
        twin = R.render_answer(index, config, dict(req, max_contexts=100000))
        twins = {0: twin['patterns'][0]}
        self.assertEqual(C.PASS, outcome(C.compare_pattern(index, config, req, answer, 5,
                                                           twins)))
        # the last listed context replaced by a context the complete list has later
        later = twin['patterns'][0]['results'][7]
        answer['patterns'][0]['results'][4] = later
        failed = {x['check'] for x in C.compare_pattern(index, config, req, answer, 5,
                                                        twins)[0][1].failed()}
        self.assertIn('results.complete', failed)


class UnmaskedRelations(unittest.TestCase):
    """Verification H3 (G3-G5) and L1-L2: on a graph without its mask a completed count is exact
    or bounds; a withhold on U needs U > max_contexts, with the note threshold_upper_bound exactly
    when lower <= max_contexts; a sticky later pattern answers unknown counts with zero work."""

    def setUp(self):
        self.index = fixture_index('basic')
        self.config = C.Config('unmasked', 'basic', 11, masked=False)

    def failed(self, req, answer, mc=100000):
        return {x['check'].split('[')[0] for _, ch, _ in
                C.compare_pattern(self.index, self.config, req, answer, mc) for x in ch.failed()}

    def test_g3_g4_a_completed_count_is_never_unknown_or_at_least(self):
        req = {'mode': 'count', 'patterns': [{'dna': 'GAATTC'}]}
        for relation in ('unknown', 'at_least'):
            answer = R.render_answer(self.index, self.config, req)
            c = answer['patterns'][0]['counts']['contexts']
            c.update({'relation': relation, 'value': None if relation == 'unknown' else 0})
            self.assertIn('counts.total', self.failed(req, answer), relation)

    def record_start_pattern(self):
        """A record's first bases with at least two contexts, some source dummy holding them
        (U above the truth)."""
        for r in self.index.records:
            for L in (4, 5, 6, 3):
                p = O.Pattern('dna', r.seq[:L])
                if p.slot or len(O.contexts(self.index, p)) < 2:
                    continue
                if O.dummy_contexts(self.index, p):
                    return r.seq[:L]
        self.fail('no record-start pattern with two contexts')

    def bounded_withhold(self, pattern, mc, note, slack=0):
        req = {'mode': 'all_or_count', 'max_contexts': mc, 'patterns': [{'dna': pattern}]}
        exp = C.Expectation(self.index, req, req['patterns'][0])
        truth = exp.counts['total']
        dummies = exp.with_upper(self.index)['0']
        answer = R.render_answer(self.index, self.config, dict(req, max_contexts=100000))
        e = answer['patterns'][0]
        c = e['counts']['contexts']

        def b(v, extra):
            return {'relation': 'bounds', 'unit': 'graph_contexts', 'value': v, 'lower': v,
                    'upper': v + extra, 'estimate': v}
        c.update(b(truth - slack, dummies['total'] + slack))
        c['suffix'] = b(exp.counts['suffix'], dummies['suffix'])
        for p in c['by_offset']:
            c['by_offset'][p] = b(exp.counts['by_offset'][p], dummies['by_offset'][p])
        for o in c['by_strand']:
            c['by_strand'][o] = b(exp.counts['by_orient'][o], dummies['by_orient'][o])
        e.update({'withheld': {'reason': 'count_above_threshold'}, 'results': [],
                  'returned': 0, 'retrieval_complete': False,
                  'notes': (['threshold_upper_bound'] if note else [])
                  + ['estimate_sampled_dummy_fraction']})
        return req, answer, truth, truth + dummies['total']

    def test_g5_a_withhold_on_u_needs_u_above_max_contexts(self):
        # a start of a record: source dummies hold it, U > truth
        pattern = self.record_start_pattern()
        req, answer, truth, U = self.bounded_withhold(pattern, 10 ** 4, True)
        self.assertGreater(U, truth)
        self.assertIn('release.all_or_count', self.failed(req, answer, 10 ** 4))
        req, answer, truth, U = self.bounded_withhold(pattern, U - 1, True)
        self.assertNotIn('release.all_or_count', self.failed(req, answer, U - 1))

    def test_l2_the_note_exactly_when_lower_is_within_max_contexts(self):
        """The truth above max_contexts but lower within it: the note is due (SPEC 7.5: when
        its lower <= max_contexts), which the comparison required only when the truth fit."""
        pattern = self.record_start_pattern()
        exp_truth = len(O.contexts(self.index, O.Pattern('dna', pattern)))
        self.assertGreater(exp_truth, 1)
        mc = exp_truth - 1
        req, answer, truth, U = self.bounded_withhold(pattern, mc, False, slack=1)
        self.assertIn('release.all_or_count', self.failed(req, answer, mc))
        req, answer, truth, U = self.bounded_withhold(pattern, mc, True, slack=1)
        self.assertNotIn('release.all_or_count', self.failed(req, answer, mc))

    def test_l1_a_sticky_later_pattern_states_no_work(self):
        index, config = fixture_index('basic'), CONFIGS['basic']
        req = {'mode': 'count', 'max_steps': 10 ** 6,
               'patterns': [{'iupac': 'GNNNA'}, {'dna': 'ACGTA'}]}
        answer = R.render_answer(index, config, req, 100000, (0, 2))
        self.assertEqual(C.INCONCLUSIVE, outcome(C.compare_pattern(index, config, req, answer)))
        answer['patterns'][1]['work']['steps'] = 3
        failed = {x['check'] for x in C.compare_pattern(index, config, req, answer)[1][1].failed()}
        self.assertIn('stop.sticky_unknown', failed)
