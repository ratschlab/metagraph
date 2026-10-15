"""Offline tests of the lane O run's plumbing: the request generator, the report and the command's
refusals. No server, no engine."""

import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oracle import fixtures, report, requests_gen  # noqa: E402
from oracle import oracle as O  # noqa: E402
from oracle import run_oracle  # noqa: E402


def fixture_index(mode='basic', k=11):
    return O.Index(fixtures.all_records(fixtures.build_records()), k, mode)


class Generator(unittest.TestCase):
    def test_deterministic_from_the_seed(self):
        idx = fixture_index()
        a = requests_gen.generate(idx, 'basic_k11', 1, True, True, True, 10, 4, 5)
        b = requests_gen.generate(idx, 'basic_k11', 1, True, True, True, 10, 4, 5)
        c = requests_gen.generate(idx, 'basic_k11', 2, True, True, True, 10, 4, 5)
        self.assertEqual(json.dumps(a, sort_keys=True), json.dumps(b, sort_keys=True))
        self.assertNotEqual(json.dumps(a, sort_keys=True), json.dumps(c, sort_keys=True))

    def test_requests_stay_within_phase_one(self):
        for k in (11, 15):
            idx = fixture_index(k=k)
            for item in requests_gen.generate(idx, f'basic_k{k}', 5, True, True, True, 40, 8, 10):
                req = item['request']
                if item['route'] == 'resolve':
                    self.assertGreaterEqual(len(req['sequence']), 1)
                    continue
                self.assertLessEqual(len(req['patterns']), 16)
                m = req.get('max_mismatches', 0)
                self.assertIn(m, (0, 1, 2, 3))
                if m:
                    # substitutions with labels are a later increment (SPEC 21.16)
                    self.assertNotEqual('all', (req.get('output') or {}).get('labels'))
                for p in req['patterns']:
                    kind = next(x for x in ('dna', 'iupac', 'protein') if x in p)
                    pat = O.Pattern(kind, p[kind], req.get('genetic_code', 1))
                    if not pat.slot:
                        self.assertLessEqual(pat.length, k, p)

    def test_resolve_requests_name_their_sources(self):
        idx = fixture_index()
        reqs = requests_gen.generate(idx, 'basic_k11', 1, True, True, True, 0, 0, 30)
        supported = 0
        for item in reqs:
            try:
                out = O.resolve(idx, item['request']['sequence'], item['request']['labels'])
            except ValueError:
                continue
            supported += sum(1 for l in out['labels'] if l['kmers_supported'])
        self.assertGreater(supported, 15)


class Breadth(unittest.TestCase):
    """Verification M2, M3: every seed draws /resolve trace requests that tell trace from k-mer
    presence, and varied constant-cost /traverse requests."""

    def test_resolve_trace_requests_discriminate_in_every_seed(self):
        idx = fixture_index()
        for seed in range(1, 7):
            reqs = requests_gen.generate(idx, 'basic_k11', seed, True, True, True, 0, 0, 30)
            telling = 0
            for item in reqs:
                req = item['request']
                if req.get('support') != 'trace':
                    continue
                try:
                    trace = O.resolve(idx, req['sequence'], req['labels'], support='trace')
                    kmer = O.resolve(idx, req['sequence'], req['labels'])
                except ValueError:
                    continue            # a refused label (a record shorter than k)
                telling += trace['labels'] != kmer['labels']
            self.assertGreaterEqual(telling, 3, seed)

    def test_constant_cost_requests_vary_in_every_seed(self):
        from oracle import traverse_gen, traverse_oracle
        idx = fixture_index()
        w = traverse_oracle.Walker(idx)
        settings = set()
        for seed in range(1, 7):
            reqs = traverse_gen.generate(idx, w, 'basic_k11', seed, None, 80, 0)
            costs = [(r['request']['strategy']['labels']['change_cost']['value'],
                      r['request']['strategy']['labels']['loss_budget'],
                      r['request']['strategy']['labels'].get('switch_on', 'loss'))
                     for r in reqs
                     if ((r['request']['strategy'].get('labels') or {}).get('change_cost') or {})
                     .get('model') == 'constant']
            self.assertGreaterEqual(len(costs), 8, seed)
            settings |= set(costs)
        self.assertGreaterEqual(len(settings), 8)


class Report(unittest.TestCase):
    def test_markdown_states_scope_disagreements_and_missed_mutations(self):
        rep = {'engine': {'path': '/x/metagraph', 'version': 'v'}, 'seed': 1,
               'started_utc': 'a', 'ended_utc': 'b',
               'summary': {'verdict': 'FAIL', 'entries': 2, 'requests': 1, 'pass': 1, 'fail': 1,
                           'inconclusive': 0, 'exploratory': 0},
               'configs': {'basic_k11': {'graph_mode': 'basic', 'k': 11, 'counting': 'exact',
                                         'requests': 1, 'entries': 2,
                                         'outcomes': {'pass': 1, 'fail': 1}}},
               'disagreements': [{'config': 'basic_k11', 'request_id': 'P-1', 'position': 0,
                                  'request': {'patterns': [{'dna': 'ACGT'}]},
                                  'failed': [{'check': 'counts.total', 'expected': 3, 'got': 4}],
                                  'oracle': {'total': 3}}],
               'exploratory': [],
               'mutations': {'applied': 2, 'caught': 1, 'by_mutation': {'drop_result':
                                                                        {'applied': 2,
                                                                         'caught': 1}},
                             'missed': [{'request_id': 'P-1', 'position': 0,
                                         'mutation': 'drop_result',
                                         'families': ['results.complete']}]}}
        md = report.markdown(rep)
        for words in ('**Verdict: FAIL.**', '## Scope', 'Not covered', '/traverse (phase 2)',
                      '### basic_k11: counts.total', '`P-1` pattern 0', 'drop_result'):
            self.assertIn(words, md)


class Provenance(unittest.TestCase):
    """Verification M1: the report names the binary by its content and its build tree's HEAD."""

    def test_a_binary_in_a_git_tree(self):
        here = os.path.abspath(__file__)
        p = report.engine_provenance(here)
        self.assertEqual(64, len(p['sha256']))
        self.assertEqual(os.path.getsize(here), p['size'])
        self.assertRegex(p['build_head'] or '', '^[0-9a-f]{40}$')
        self.assertIsInstance(p['build_src_uncommitted'], int)
        self.assertIn(p['sha256'], report.markdown(self.minimal(p)))

    def test_a_binary_outside_a_git_tree(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'metagraph')
            with open(path, 'wb') as f:
                f.write(b'not an engine')
            p = report.engine_provenance(path)
            self.assertIsNone(p['build_head'])
            self.assertEqual(13, p['size'])

    def minimal(self, engine):
        return {'engine': engine, 'seed': 1, 'started_utc': 'a', 'ended_utc': 'b',
                'summary': {'verdict': 'PASS', 'entries': 0, 'requests': 0, 'pass': 0,
                            'fail': 0, 'inconclusive': 0, 'exploratory': 0},
                'configs': {}, 'disagreements': [], 'exploratory': [],
                'mutations': {'applied': 0, 'caught': 0, 'by_mutation': {}, 'missed': []}}


class Command(unittest.TestCase):
    def test_a_non_empty_out_dir_is_refused(self):
        with tempfile.TemporaryDirectory() as d:
            open(os.path.join(d, 'keep'), 'w').close()
            self.assertEqual(2, run_oracle.main(['--metagraph', '/nonexistent', '--out', d]))
            self.assertEqual(['keep'], os.listdir(d))

    def test_an_unknown_config_is_refused(self):
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, 'run')
            self.assertEqual(2, run_oracle.main(['--metagraph', '/nonexistent', '--out', out,
                                                 '--configs', 'basic_k11,nope']))
            self.assertFalse(os.path.exists(out))


if __name__ == '__main__':
    unittest.main()
