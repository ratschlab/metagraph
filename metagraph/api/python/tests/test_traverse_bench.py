"""The benchmark scripts' comparisons (scripts/traversal/bench_traverse.py, bench_library.py),
as the external review of pass 5 asked for them: a deadline-limited request is compared by
its CERTIFIED REACH (complete_to_bp, steps, complete walks, the stop) and its CONSUMED
WORK (work units, else the annotation rows requested), not by its elapsed time, which is
about the budget on both sides; a completed request's latency only on IDENTICAL COMPLETED
WORK (the same result and logical work counters), a pair whose work differs flagged and
left out; the physical path-cache reuse and peak bytes recorded where a server reports
them and tolerated where it does not; results written by the first script version (the
staging baseline, mini_before) still read. Synthetic responses: no server, no network.
"""

import gzip
import importlib.util
import json
import os
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.normpath(os.path.join(HERE, '..', '..', '..', 'scripts', 'traversal'))


def _load(name):
    spec = importlib.util.spec_from_file_location(name + '_under_test', os.path.join(SCRIPTS, name + '.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


BT = _load('bench_traverse')
BLIB = _load('bench_library')


def _arm(bp, leaves_by_reason, steps, cap=None, work_units=None, max_bp=None):
    counters = {'steps': steps, 'successor_enumerations': 2 * steps}
    if work_units is not None:
        counters['work_units'] = work_units
    a = {'status': 'truncated' if cap else 'complete', 'complete_to_bp': bp,
         'counts': {'max_bp': max_bp if max_bp is not None else bp, 'bases': bp,
                    'leaves': sum(leaves_by_reason.values()), 'segments': 3,
                    'leaves_by_reason': leaves_by_reason},
         'counters': counters}
    if cap:
        a['cap_trigger'] = {'reason': cap}
    return a


def _seed(text, elapsed, arms, rows_requested=10, timing_extra=None, outcome='complete'):
    timing = {'elapsed_ms': elapsed, 'annotation_fetch_ms': elapsed / 2, 'rows_fetched': 20,
              'tuple_rows_fetched': 0, 'coords_mapped': 0, 'cache_hits': 1}
    timing.update(timing_extra or {})
    return {'seed': {'seed_id': 's'}, 'arms': arms, 'outcome': {'walks': outcome},
            'annotation': {'rows_requested': rows_requested, 'keys_mapped': 2 * rows_requested},
            'timing': timing, 'graphlet': text, 'limitations': []}


def _record(n, rid, suite, resp, route='/traverse', request_sha='r'):
    meta = {'n': n, 'id': rid, 'suite': suite, 'route': route, 'method': 'POST', 'params': {},
            'windows': [], 'wall_ms': 1.0, 'http': 200, 'request_sha': request_sha + rid}
    return meta, {'body': rid}, resp


def _traverse(seed, budget_ms=5000, usage=None):
    resp = {'strategy': {'bounds': {'time_budget_ms': budget_ms}}, 'results': [seed],
            'timing': {'elapsed_ms': seed['timing']['elapsed_ms'], 'serialize_ms': 0.0}}
    if usage is not None:
        resp['usage'] = usage
    return resp


def _write_run(d, label, entries, raw=False, strip_v2=False):
    """A run directory: results.json from the entries (as the script assembles them) and,
    with |raw|, its raw/ responses; |strip_v2| writes the record as script version 1 did."""
    os.makedirs(d, exist_ok=True)
    records = [BT.analyze(m, req, resp) for m, req, resp in entries]
    if strip_v2:
        for rec in records:
            for s in rec.get('results') or []:
                for k in ('reach', 'work', 'physical'):
                    s.pop(k, None)
                for a in s['arms'].values():
                    a.pop('leaves_by_reason', None)
                    a.pop('counters', None)
    meta = {'label': label, 'plan_version': 1, 'panel': {'name': 'p', 'sha': 'x'}}
    res = BT.assemble(meta, records)
    BT.write_json_atomic(os.path.join(d, 'results.json'), res)
    if raw:
        os.makedirs(os.path.join(d, 'raw'))
        for m, req, resp in entries:
            with gzip.open(os.path.join(d, 'raw', '%03d_%s.json.gz' % (m['n'], m['id'])), 'wt') as f:
                json.dump({'meta': m, 'request': req, 'response': resp}, f)
    return os.path.join(d, 'results.json')


TEXT = 'H mgt 1\nS seed ACGT\nX a\n'
TEXT2 = 'H mgt 1\nS seed ACGT\nX b\n'


def _runs():
    """Run A and run B over the same five requests."""
    A, B = [], []
    # completed, identical work: B twice as fast
    A.append(_record(1, 'walks.x.named.warm', 'walks',
                     _traverse(_seed(TEXT, 100.0, {'right': _arm(50, {'dead_end': 2}, 50)}))))
    B.append(_record(1, 'walks.x.named.warm', 'walks',
                     _traverse(_seed(TEXT, 50.0, {'right': _arm(50, {'dead_end': 2}, 50)},
                                     timing_extra={'seed_phase_ms': 1.0,
                                                   'path_cache': {'hits': 7, 'stored_rows_read': 3,
                                                                  'rows_kept': 9, 'bytes_kept': 900,
                                                                  'peak_bytes': 4096}}))))
    # completed, the same result but other logical work (the counting changed)
    A.append(_record(2, 'walks.y.named.warm', 'walks',
                     _traverse(_seed(TEXT, 100.0, {'right': _arm(50, {'dead_end': 1}, 50)}, rows_requested=10))))
    B.append(_record(2, 'walks.y.named.warm', 'walks',
                     _traverse(_seed(TEXT, 10.0, {'right': _arm(50, {'dead_end': 1}, 50)}, rows_requested=12))))
    # completed, a different result
    A.append(_record(3, 'walks.z.named.warm', 'walks',
                     _traverse(_seed(TEXT, 100.0, {'right': _arm(50, {'dead_end': 1}, 50)}))))
    B.append(_record(3, 'walks.z.named.warm', 'walks',
                     _traverse(_seed(TEXT2, 10.0, {'right': _arm(50, {'dead_end': 1}, 50)}))))
    # deadline-limited: about the budget on both sides; B reached three times as far with
    # twice the work units, and one of its walks completed
    A.append(_record(4, 'walks.w.named.first', 'walks',
                     _traverse(_seed(TEXT, 5001.0, {'right': _arm(100, {'time_budget': 2}, 100, 'time_budget',
                                                                  work_units=1000)}, outcome='partial'))))
    B.append(_record(4, 'walks.w.named.first', 'walks',
                     _traverse(_seed(TEXT2, 5003.0, {'right': _arm(300, {'dead_end': 1, 'time_budget': 1}, 300,
                                                                   'time_budget', work_units=2000)},
                                     outcome='partial'))))
    # /resolve, completed
    for run, ms in ((A, 8.0), (B, 6.0)):
        run.append(_record(5, 'resolve.x.explicit_warm', 'resolve',
                           {'labels': [], 'num_kmers': 30, 'timing': {'elapsed_ms': ms}}, route='/resolve'))
    return A, B


class TestBenchTraverseAnalysis(unittest.TestCase):
    def test_reach_work_and_physical_fields(self):
        A, B = _runs()
        rec = BT.analyze(*B[0])
        s = rec['results'][0]
        self.assertEqual({'bp': 50, 'bp_left': None, 'bp_right': 50, 'walked_bp': 50, 'steps': 50,
                          'segments': 3, 'leaves': 2, 'complete_walks': 2, 'stopped_walks': 0,
                          'stop': 'complete'}, s['reach'])
        self.assertEqual(10, s['work']['rows_requested'])
        self.assertEqual(100, s['work']['successor_enumerations'])
        self.assertIsNone(s['work']['work_units'])
        self.assertEqual((7, 3, 4096), (s['physical']['reuse'], s['physical']['stored_rows_read'],
                                        s['physical']['peak_cache_bytes']))
        self.assertEqual(7, s['timing']['path_cache.hits'])
        # an older server: no path-cache fields, nothing invented
        s = BT.analyze(*A[0])['results'][0]
        self.assertEqual((None, None), (s['physical']['reuse'], s['physical']['peak_cache_bytes']))
        # a deadline-limited seed: the stop and the budget's work units
        s = BT.analyze(*B[3])['results'][0]
        self.assertEqual(('time_budget', 1, 1, 2000), (s['reach']['stop'], s['reach']['complete_walks'],
                                                       s['reach']['stopped_walks'], s['work']['work_units']))
        self.assertTrue(s['time_bound'])

    def test_a_server_of_the_same_release_without_the_object_reused_nothing(self):
        m, req, resp = _record(1, 'w', 'walks', _traverse(_seed(TEXT, 5.0, {'right': _arm(5, {'dead_end': 1}, 5)},
                                                                timing_extra={'seed_phase_ms': 0.5})))
        s = BT.analyze(m, req, resp)['results'][0]
        self.assertEqual((0, 0), (s['physical']['reuse'], s['physical']['peak_cache_bytes']))

    def test_an_attempts_usage_is_the_seeds_work(self):
        usage = {'work_units': 1254, 'reason': 'complete', 'elapsed_ms': 12.0,
                 'memory': {'peak_admitted_bytes': 4096}}
        m, req, resp = _record(1, 'w', 'walks', _traverse(_seed(TEXT, 5.0, {'right': _arm(5, {'dead_end': 1}, 5)}),
                                                          usage=usage))
        rec = BT.analyze(m, req, resp)
        self.assertEqual(1254, rec['usage']['work_units'])
        self.assertEqual(4096, rec['usage']['peak_admitted_bytes'])
        self.assertEqual(1254, rec['results'][0]['work']['usage_work_units'])
        self.assertEqual((1254, 'wu'), BT.work_amount(rec['results'][0]['work']))

    def test_the_summary_has_the_reach_and_work_table(self):
        A, _ = _runs()
        res = BT.assemble({'label': 'a'}, [BT.analyze(*e) for e in A])
        text = BT.summary_md(res)
        self.assertIn('## certified reach and consumed work (per seed)', text)
        self.assertIn('| walks.w.named.first | time_budget | -+100 |', text)


class TestBenchTraverseCompare(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def compare(self, **kw):
        A, B = _runs()
        pa = _write_run(os.path.join(self.tmp.name, 'a'), 'A', A, **kw)
        pb = _write_run(os.path.join(self.tmp.name, 'b'), 'B', B)
        return BT.compare(pa, pb)

    def row(self, text, rid):
        return next(l for l in text.split('\n') if l.startswith('| %s |' % rid))

    def test_completed_latency_only_on_identical_work(self):
        text, flags = self.compare()
        # the suite's latency: only the identical pair (0.5); the other two are excluded
        self.assertIn('| walks | 1 | 0.500 | 2 | 1 | 3.000 | 1/0/0 | 2.000 |', text)
        self.assertIn('| resolve | 1 | 0.750 | 0 | 0 |', text)
        self.assertIn('| 100.0 | 50.0 | 0.50 |', self.row(text, 'walks.x.named.warm'))
        self.assertIn('| 7 | - | 4096 |', self.row(text, 'walks.x.named.warm'))
        row = self.row(text, 'walks.y.named.warm')
        self.assertIn('(0.10)', row)
        self.assertIn('rows_requested 10->12', row)
        self.assertIn('Completed with equal results but different logical work: walks.y.named.warm', text)
        self.assertIn('(0.10)', self.row(text, 'walks.z.named.warm'))
        self.assertEqual(['walks.z.named.warm'], flags)

    def test_deadline_limited_by_reach_and_work(self):
        text, _ = self.compare()
        row = [c.strip() for c in self.row(text, 'walks.w.named.first').split('|')[1:-1]]
        self.assertEqual(['walks.w.named.first', '5001', '5003', 'time_budget', 'time_budget', '100', '300',
                          '3.00', '100', '300', '0/2', '1/2', '1000 wu', '2000 wu', '2.00', '-', '-', '-', '-',
                          'time-bound'], row)

    def test_a_first_version_results_json_is_read_from_its_raw_responses(self):
        text, flags = self.compare(raw=True, strip_v2=True)
        self.assertIn('analysed again from its raw/ responses', text)
        self.assertIn('| walks | 1 | 0.500 | 2 | 1 | 3.000 | 1/0/0 | 2.000 |', text)

    def test_a_first_version_results_json_without_raw_uses_what_it_recorded(self):
        text, flags = self.compare(strip_v2=True)
        row = [c.strip() for c in self.row(text, 'walks.w.named.first').split('|')[1:-1]]
        # bp and the rows requested were recorded; steps and complete walks were not
        self.assertEqual(['100', '300', '3.00', '-', '300', '-/2', '1/2', '10 rows', '2000 wu', '-'],
                         row[5:15])
        self.assertEqual(['walks.z.named.warm'], flags)


class TestBenchTraverseCoordinates(unittest.TestCase):
    """--coordinates (feature level 6): support-trace walks ask for record coordinates, a run
    with them pairs with an opt-out run of the same windows (the same walks: the pair's
    latency is the coordinates' cost), the summary and the comparison state the block's
    size; --coordinates-estimate states the block an opt-out run's trace walks would
    carry, from its recorded responses alone."""

    BLOCK = {'arms': {'right': [{'from_bp': 0, 'label': 0, 'occurrences': [[1000000, 1000050],
                                                                           [2000000, 2000050]],
                                 'run': 0, 'to_bp': 50}]},
             'complete': True, 'k': 31, 'kind': 'record', 'max_occurrences': 16,
             'seed': [{'label': 0, 'occurrences': [[999960, 1000000], [1999960, 2000000]]}]}
    BODY = ('H mgt 1 31 basic $ACGT c t k 64 1000 0 walk x * *\nS s 40 10 1 ' + 'A' * 40
            + '\nL h 0 0 0 acc\nO c c c i\nA r c 50 p 0 0 1 1 0 . * 1 * 1 1 0 0 0 50\n'
            'G . 0 50 ' + 'A' * 50 + ' * * * . * .\nR 0 0 0 50 X * * * *\nZ 8\n')

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def test_the_plan_asks_only_trace_walks(self):
        panel = BT.Panel(os.path.join(SCRIPTS, 'bench_panel_mini_refseq.json'))
        plain = BT.Planner(panel, 0).build(['walks'])
        coords = BT.Planner(panel, 0, 16).build(['walks'])
        self.assertEqual([r['id'] for r in plain], [r['id'] for r in coords])
        for a, b in zip(plain, coords):
            out = b['body']['strategy']['output']
            if a['params']['mode'] == 'trace':
                self.assertEqual((True, 16), (out['coordinates'], out['max_coordinate_occurrences']))
                self.assertNotEqual(BT.sha(a['body']), BT.sha(b['body']))
            else:
                self.assertNotIn('coordinates', out)
            # apart from the coordinates, the same request
            self.assertEqual(BT.sha(a['body']), BT.sha(BT.without_coordinates(b['body'])))

    def test_facts_and_pairing(self):
        def run(name, with_block):
            seed = _seed(TEXT, 10.0, {'right': _arm(50, {'max_extension_bp': 1}, 50)})
            if with_block:
                seed['coordinates'] = self.BLOCK
            meta, req, resp = _record(1, 'walks.x.trace.warm', 'walks', _traverse(seed))
            body = {'strategy': {'output': {'detail': 'graphlet'}}}
            if with_block:
                body['strategy']['output'].update(coordinates=True, max_coordinate_occurrences=16)
            meta['request_sha'] = BT.sha(body)
            meta['request_sha_nocoords'] = BT.sha(BT.without_coordinates(body))
            return _write_run(os.path.join(self.tmp.name, name), name, [(meta, body, resp)])
        pa, pb = run('a', False), run('b', True)
        rec = BT.load_results(pb)['requests'][0]['results'][0]
        c = rec['coordinates']
        self.assertEqual(('record', 2, 4, 4, 0, 1), (c['kind'], c['max_list'], c['occurrences'],
                                                     c['occurrences_total'], c['lists_cut'],
                                                     c['runs']))
        self.assertIn('## record coordinates', BT.summary_md(BT.load_results(pb)))
        text, flags = BT.compare(pa, pb)
        self.assertEqual([], flags)
        # the same request apart from the coordinates: a completed pair, its latency compared
        self.assertIn('| walks | 1 |', text)
        self.assertIn('## record coordinates', text)
        self.assertNotIn('request differs', text.split('## completed')[1].split('\n\n')[2])

    def test_the_estimate_of_a_recorded_run(self):
        d = os.path.join(self.tmp.name, 'run')
        os.makedirs(os.path.join(d, 'raw'))
        resp = {'results': [{'seed': {'num_seed_labels': 1, 'length_bp': 40}, 'graphlet': self.BODY}]}
        req = {'strategy': {'support': 'trace'}}
        with gzip.open(os.path.join(d, 'raw', '001_walks.x.trace.first.json.gz'), 'wt') as f:
            json.dump({'meta': {'id': 'walks.x.trace.first'}, 'request': req, 'response': resp}, f)
        lines = []
        rows = BT.coordinates_estimate([d], 7, out=lines.append)
        self.assertEqual(1, len(rows))
        r = rows[0]
        self.assertEqual((1, 1, 'record'), (r['seed_labels'], r['runs'], r['kind']))
        # exactly the server's shape with one occurrence per entry and 7-digit positions
        block = {'arms': {'right': [{'from_bp': 0, 'label': 0, 'occurrences': [[1000000, 1000050]],
                                     'run': 0, 'to_bp': 50}]},
                 'complete': True, 'k': 31, 'kind': 'record', 'max_occurrences': 16,
                 'seed': [{'label': 0, 'occurrences': [[1000000, 1000040]]}]}
        self.assertEqual(len(json.dumps(block, separators=(',', ':'), sort_keys=True)),
                         r['block_bytes_at_least'])
        self.assertTrue(os.path.exists(os.path.join(d, 'coordinates_estimate.json')))


class TestBenchLibraryCompare(unittest.TestCase):
    def test_time_only_on_the_same_completed_work(self):
        def run(label, rows):
            return {'meta': {'label': label, 'library': {'path': label, 'source_sha': 'x', 'git_head': 'h'}},
                    'items': [{'key': 'a'}, {'key': 'b'}], 'rows': rows, 'compare': []}
        A = run('A', [{'key': 'a', 'op': 'to_fasta', 't_min_ms': 1.0, 'fp': 'f1', 'usage': {'lwu': 10}},
                      {'key': 'b', 'op': 'to_fasta', 't_min_ms': 1.0, 'fp': 'f2', 'usage': {'lwu': 10}},
                      {'key': 'a', 'op': 'summary', 't_min_ms': 2.0, 'fp': 'f3', 'usage': {'lwu': 5}}])
        B = run('B', [{'key': 'a', 'op': 'to_fasta', 't_min_ms': 2.0, 'fp': 'f1', 'usage': {'lwu': 10}},
                      {'key': 'b', 'op': 'to_fasta', 't_min_ms': 0.1, 'fp': 'CHANGED', 'usage': {'lwu': 10}},
                      {'key': 'a', 'op': 'summary', 't_min_ms': 1.0, 'fp': 'f3', 'usage': {'lwu': 6}}])
        with tempfile.TemporaryDirectory() as d:
            pa, pb = os.path.join(d, 'a.json'), os.path.join(d, 'b.json')
            for p, x in ((pa, A), (pb, B)):
                with open(p, 'w') as f:
                    json.dump(x, f)
            text, mism = BLIB.compare_runs(pa, pb)
        row = next(l for l in text.split('\n') if l.startswith('| to_fasta |'))
        # the changed item is left out of the ratio: 2.0 / 1.0 alone
        self.assertTrue(row.startswith('| to_fasta | 1 | 2.000 | 2.000 |'), row)
        self.assertIn('| 1 | 10 -> 10 |', row)
        self.assertEqual(1, len(mism))
        self.assertIn('summary a: 5 -> 6 lwu', text)



class TestBenchLibraryInterleaved(unittest.TestCase):
    def test_the_summary_leaves_changed_results_out_and_states_the_noise(self):
        lib = {'path': 'x', 'source_sha': 's', 'git_head': 'h'}
        rows = [
            {'key': 'a', 'op': 'routes_one', 'err': None, 'a_ms': 1.0, 'b_ms': 1.1, 'a1_ms': 1.0, 'a2_ms': 1.05,
             'fp_a': 'f', 'fp_b': 'f'},
            {'key': 'b', 'op': 'routes_one', 'err': None, 'a_ms': 2.0, 'b_ms': 1.8, 'a1_ms': 2.0, 'a2_ms': 2.0,
             'fp_a': 'g', 'fp_b': 'g'},
            {'key': 'c', 'op': 'routes_one', 'err': None, 'a_ms': 1.0, 'b_ms': 0.1, 'a1_ms': 1.0, 'a2_ms': 1.0,
             'fp_a': 'h', 'fp_b': 'CHANGED'},
            {'key': 'a', 'op': 'claims', 'err': 'A: ValueError', 'a_ms': None, 'b_ms': None, 'a1_ms': None,
             'a2_ms': None, 'fp_a': None, 'fp_b': None}]
        res = {'meta': {'label': 't', 'library_a': lib, 'library_b': lib, 'method': 'm', 'rounds': 2,
                        'python': '3', 'started': 'now', 'seconds': 1, 'ops': ['routes_one', 'claims']},
               'rows': rows}
        text, changed = BLIB.interleave_md(res)
        self.assertEqual(['c'], [r['key'] for r in changed])
        row = next(l for l in text.split('\n') if l.startswith('| routes_one |'))
        # sqrt(1.1 * 0.9) and the A2/A1 noise of the two compared items
        self.assertTrue(row.startswith('| routes_one | 2 | 0.995 |'), row)
        self.assertIn('| 1.025 | 1.05 |', row)
        self.assertIn('1 results changed', text)
        self.assertIn('not timed (an op that does not apply): claims', text)

    def test_two_workers_of_one_library_agree(self):
        lib = os.path.normpath(os.path.join(HERE, '..'))
        with tempfile.TemporaryDirectory() as d:
            import contextlib
            import io
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                code = BLIB.main(['--interleave', lib, lib, '--max-items', '2', '--ops', 'summary,to_fasta',
                                  '--rounds', '1', '--out', d, '--label', 'same'])
            self.assertEqual(0, code, buf.getvalue())
            with open(os.path.join(d, 'ab_results.json')) as f:
                res = json.load(f)
        self.assertEqual(4, len(res['rows']))
        self.assertTrue(all(r['fp_a'] == r['fp_b'] for r in res['rows'] if not r['err']))


if __name__ == '__main__':
    unittest.main()
