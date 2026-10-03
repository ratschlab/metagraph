"""The MCP tool functions (§6) against an in-process fake backend: handles, summaries,
paging bound to handle+tool+arguments, per-tool ceilings, oversized rows, {name, ref}
everywhere, the evidence block, the oversize policy, replay of an expired handle."""

import copy
import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402

from metagraph.traverse import GraphletStore, GraphletView, TraverseClient  # noqa: E402
from metagraph.traverse.mcp_tools import GraphletTools  # noqa: E402

def seed_of(name):
    """The seed sequence of a fixture (its S record)."""
    return T.doc_text(name).split('\n')[1].split(' ')[5]


# (seed sequence, label mode) -> the fixture whose detail: graphlet response the fake
# server returns (merge and annotate share their seed; linear is never fetched here)
BY_SEED = {(seed_of(n), 'annotate' if n in ('annotate', 'same_name') else 'constrain'): n
           for n in ('fork', 'merge', 'switch_chain', 'annotate', 'same_name', 'reentry')}
SEED = {n: seed_of(n) for n in ('fork', 'merge', 'switch_chain', 'annotate', 'same_name')}
# switch_chain's walk 1 continuation, answered by a stand-in document (the fake has no index)
CONTINUATION = 'TTACTCGTAGCCGGGCGTGA'
BY_SEED[(CONTINUATION, 'constrain')] = 'reentry'
# ... which stands in for the parent's index: it carries the parent's index_meta_fp (a
# continuation is checked to come from its parent's index, §3.1)
STAND_IN = {CONTINUATION: 'switch_chain'}


def meta_fp_of(name):
    return T.doc_text(name).split('\n')[0].split(' ')[-1]


class FakeClient(TraverseClient):
    def __init__(self):
        super().__init__('fake', 0)
        self.requests = []
        # an index manifest digest to stamp into H (the CLI fixtures were made without one)
        self.fp = None

    def traverse_raw(self, request):
        self.requests.append(copy.deepcopy(request))
        seq = request['seeds'][0]['sequence']
        mode = ((request.get('strategy') or {}).get('labels') or {}).get('mode', 'constrain')
        if (seq, mode) not in BY_SEED:
            return {'release': '', 'results': [{'seed': {'seed_id': ''}, 'error': 'no carrier',
                                                'outcome': {'walks': 'failed'},
                                                'limitations': []}]}
        resp = copy.deepcopy(T.doc_json(BY_SEED[(seq, mode)], 'graphlet'))
        if seq in STAND_IN:
            r = resp['results'][0]
            own, parent = meta_fp_of(BY_SEED[(seq, mode)]), meta_fp_of(STAND_IN[seq])
            r['graphlet'] = r['graphlet'].replace(' walk fixtures * %s\n' % own,
                                                  ' walk fixtures * %s\n' % parent, 1)
            r['graphlet_bytes'] = len(r['graphlet'].encode('utf-8'))
            resp['capabilities']['index_meta_fp'] = parent
        if self.fp:
            r = resp['results'][0]
            r['graphlet'] = r['graphlet'].replace(' walk fixtures * ', ' walk fixtures %s '
                                                  % self.fp, 1)
            r['graphlet_bytes'] = len(r['graphlet'].encode('utf-8'))
            resp['capabilities']['index_fp'] = self.fp
        return resp

    def capabilities(self):
        return {'release': '', 'index_fp': None}

    def resolve(self, sequence, **kw):
        return {'k': 5, 'num_kmers': 6, 'labels': [{'label': 'acc%d' % i, 'kind': 'header',
                                                    'kmers_supported': 6} for i in range(30)],
                'candidates': [{}]}


class Clock:
    t = 0.0

    def __call__(self):
        return self.t


def walk_dicts(obj):
    if isinstance(obj, dict):
        yield obj
        for v in obj.values():
            yield from walk_dicts(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from walk_dicts(v)


class ToolsCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.clock = Clock()
        self.store = GraphletStore(self.tmp.name, clock=self.clock, ttl_disk_s=100)
        self.client = FakeClient()
        self.tools = GraphletTools(self.store, {'mini': self.client})

    def tearDown(self):
        self.tmp.cleanup()

    def fetch(self, seq, **kw):
        out = self.tools.traverse_fetch(seed={'sequence': seq}, **kw)
        self.assertNotIn('error', out, out)
        return out

    def assert_label_contract(self, out):
        """Labels leave the server as {name, ref}, never with an id."""
        for d in walk_dicts(out):
            if 'ref' in d and isinstance(d.get('ref'), str) and ':' in d['ref']:
                self.assertIn('name', d)
                self.assertNotIn('id', d)

    @staticmethod
    def size(obj):
        return len(json.dumps(obj, separators=(',', ':')))


class TestTools(ToolsCase):

    def test_fetch_summary_and_evidence(self):
        out = self.fetch(SEED['merge'])
        self.assertRegex(out['handle'], r'^g_[0-9a-f]{12}$')
        self.assertEqual('inline', out['delivery'])
        self.assertEqual({'complete_to_bp', 'exact', 'support', 'reconverge', 'scope',
                          'outcome', 'limitations'}, set(out['evidence']))
        self.assertEqual('graphlet', self.client.requests[0]['strategy']['output']['detail'])
        s = self.tools.graphlet_summary(out['handle'], 'right')
        self.assertEqual(100, s['summary']['arms']['right']['complete_to_bp'])
        self.assert_label_contract(s)
        self.assertNotIn('handle', self.fetch(SEED['merge'], keep=False))

    def test_oversize_is_spooled_complete(self):
        out = self.fetch(SEED['merge'], keep=False, max_graphlet_mb=0.0001)
        self.assertEqual('spooled', out['delivery'])
        self.assertIn('handle', out)                    # stored complete behind the handle
        self.assertNotIn(T.doc_text('merge')[:40], json.dumps(out))   # never the body
        # one delivery value: the tool's, the summary's outcome, the evidence, the entry
        self.assertEqual('spooled', out['summary']['outcome']['delivery'])
        self.assertEqual('spooled', out['evidence']['outcome']['delivery'])
        # validated, then not kept resident
        self.assertFalse(self.store.in_ram(out['handle']))
        g = self.store.graphlet(out['handle'])
        self.assertEqual('spooled', g.outcome.delivery)
        self.assertEqual('spooled', self.tools.graphlet_summary(
            out['handle'])['summary']['outcome']['delivery'])
        # the spool holds the server's body, complete and unchanged
        bodies = os.path.join(self.tmp.name, 'bodies')
        self.assertEqual([T.doc_text('merge')], [T.read(os.path.join(bodies, f))
                                                 for f in os.listdir(bodies)])
        # the reverse: a server-stated delivery is what the tool reports when it does not
        # spool itself
        self.assertEqual('inline', self.fetch(SEED['merge'])['delivery'])

    def test_failed_seed(self):
        out = self.tools.traverse_fetch(seed={'sequence': 'GGGGGGGGGG'})
        self.assertEqual('seed_failed', out['error'])

    def test_walks_paging_is_bound_to_the_arguments(self):
        h = self.fetch(SEED['fork'])['handle']
        both = self.tools.graphlet_walks(h, 'right', max_bytes=8192)
        self.assertEqual(2, len(both['rows']))
        budget = self.size(both) - 1         # one row and a cursor fit, two rows do not
        first = self.tools.graphlet_walks(h, 'right', max_bytes=budget)
        self.assertEqual(2, first['total'])
        self.assertEqual(1, len(first['rows']))
        self.assertLessEqual(self.size(first), budget)
        nxt = self.tools.graphlet_walks(h, 'right', max_bytes=budget,
                                        cursor=first['next_cursor'])
        self.assertEqual([0], [r['walk'] for r in nxt['rows']])   # ranked after walk 1
        self.assertNotIn('next_cursor', nxt)
        bad = self.tools.graphlet_walks(h, 'right', rank='length', cursor=first['next_cursor'])
        self.assertEqual('bad_cursor', bad['error'])
        forged = first['next_cursor'][:-2] + ('00' if first['next_cursor'][-2:] != '00' else '11')
        self.assertEqual('bad_cursor', self.tools.graphlet_walks(h, 'right',
                                                                 cursor=forged)['error'])
        h2 = self.fetch(SEED['fork'])['handle']
        self.assertEqual('bad_cursor', self.tools.graphlet_walks(
            h2, 'right', max_bytes=budget, cursor=first['next_cursor'])['error'])
        self.assert_label_contract(first)

    def test_an_oversized_row_comes_alone_and_says_what_was_cut(self):
        h = self.fetch(SEED['merge'])['handle']
        full = self.tools.graphlet_walks(h, 'right', spell='full', max_bytes=8192)
        # the walk's row does not fit alone (merge has one walk: no cursor)
        self.assertEqual(1, full['total'])
        budget = self.size(full) - len(full['rows'][0]['sequence']) // 2
        out = self.tools.graphlet_walks(h, 'right', spell='full', max_bytes=budget)
        self.assertTrue(out['row_truncated'])
        self.assertEqual(1, len(out['rows']))
        # the field contributing the most bytes is cut first, and named
        self.assertTrue(set(out['cut_fields']) <= {'sequence', 'labels_full'})
        self.assertTrue(out['cut_fields'])
        self.assertLessEqual(self.size(out), budget)
        walk = self.store.graphlet(h).spell('right', 0)
        self.assertEqual(walk, full['rows'][0]['sequence'])
        tail = self.tools.graphlet_walks(h, 'right', spell='tail', tail_bp=4)
        self.assertEqual(walk[-4:], tail['rows'][0]['sequence'])

    def test_walk_support_labels_splits_claims(self):
        h = self.fetch(SEED['merge'])['handle']
        w = self.tools.graphlet_walk(h, 'right', 0, max_bytes=4096)
        self.assertEqual([0, 1, 3, 4, 6], [r['segment'] for r in w['rows']])
        self.assertEqual({'length_bp': 38, 'loss_used': 0.0, 'labels': 4}, w['continuation'])
        s = self.tools.graphlet_support(h, 'right', 0, max_bytes=4096)
        self.assertEqual(['split', 'merge', 'split', 'merge'],
                         [r['why'][0]['why'] for r in s['rows'] if 'why' in r])
        lab = self.tools.graphlet_labels(h, max_bytes=4096)
        self.assertEqual(4, lab['total'])
        one = self.tools.graphlet_labels(h, name='both.fa', max_bytes=4096)
        self.assertEqual(('both.fa', 'c:3', 3), (one['name'], one['ref'], one['total']))
        sp = self.tools.graphlet_splits(h, 'right', max_bytes=4096)
        self.assertEqual(['ambiguous', 'ambiguous'], [r['kind'] for r in sp['rows']])
        cl = self.tools.graphlet_claims(h, 'right', max_bytes=4096)
        self.assertEqual({'c:0', 'c:3'}, {r['label']['ref'] for r in cl['rows']})
        cl = self.tools.graphlet_claims(h, 'right', route_consistent=False, max_bytes=4096)
        self.assertEqual(6, cl['total'])
        for out in (w, s, lab, one, sp, cl):
            self.assertIn('evidence', out)
            self.assert_label_contract(out)

    def test_ambiguous_names_need_a_ref(self):
        h = self.tools.traverse_fetch(seed={'sequence': SEED['same_name']},
                                      strategy={'labels': {'mode': 'annotate'}})['handle']
        out = self.tools.graphlet_labels(h, name='ACC1')
        self.assertEqual('ambiguous_label', out['error'])
        out = self.tools.graphlet_labels(h, name={'ref': 'h:1:0'}, max_bytes=4096)
        self.assertEqual('h:1:0', out['ref'])

    def test_sequence_ceiling(self):
        h = self.fetch(SEED['merge'])['handle']
        out = self.tools.graphlet_sequence(h, 'right', walk=0, with_seed=True)
        g = self.store.graphlet(h)
        self.assertEqual(SEED['merge'] + g.spell('right', 0), out['sequence'])
        self.assertEqual(16384, out['ceiling_bytes'])
        # a long walk (the RAM copy of the entry, edited) so that a cut leaves bases
        g.arms['right'].segments[6].walk = 'ACGT' * 100
        out = self.tools.graphlet_sequence(h, 'right', walk=0, with_seed=True)
        self.assertEqual(492, len(out['sequence']))   # seed 30 + 20 + 16 + 10 + 16 + 400
        ceiling = self.size(out) - 200
        self.tools.sequence_max_bytes = ceiling
        part = self.tools.graphlet_sequence(h, 'right', walk=0, with_seed=True)
        self.assertGreater(len(part['sequence']), 100)
        self.assertTrue(part['truncated'])
        self.assertLessEqual(self.size(part), ceiling)
        got, page = part['sequence'], part
        while page.get('truncated'):           # every page under the ceiling, in order
            page = self.tools.graphlet_sequence(h, 'right', walk=0, with_seed=True,
                                                start=page['next_from'])
            self.assertLessEqual(self.size(page), ceiling)
            got += page['sequence']
        self.assertEqual(out['sequence'], got)
        self.tools.sequence_max_bytes = 16384
        seg = self.tools.graphlet_sequence(h, 'right', segment=2)
        self.assertEqual('TTTACGTTTGACAATG', seg['sequence'])

    def test_exports(self):
        h = self.fetch(SEED['fork'])['handle']
        exports = os.path.realpath(os.path.join(self.tmp.name, 'exports'))
        out = self.tools.graphlet_export(h, 'fasta', path='x.fa')
        self.assertEqual((4, os.path.join(exports, 'x.fa')), (out['records'], out['path']))
        out = self.tools.graphlet_export(h, 'gfa', path='x.gfa')
        self.assertGreater(out['records'], 10)
        out = self.tools.graphlet_export(h, 'json', path='x.json')
        with open(out['path']) as f:
            self.assertIn('arms', json.load(f))
        out = self.tools.graphlet_export(h, 'mgt', path=os.path.join(exports, 'x.mgt'))
        self.assertTrue(T.read(out['path']).split('\n')[1].startswith('J {'))
        out = self.tools.graphlet_export(h, 'json', what='compare', other=h, path='c.json')
        self.assertEqual(0, out['records'])
        # what=walk: the complete chains of the named walks (no label list cut at 8)
        out = self.tools.graphlet_export(h, 'json', what='walk', arm='right', walks=[1],
                                         path='w.json')
        with open(out['path']) as f:
            got = json.load(f)
        self.assertEqual((1, [1]), (out['records'], [w['walk'] for w in got['walks']]))
        self.assertEqual(self.store.graphlet(h).spell('right', 1), got['walks'][0]['sequence'])
        self.assertEqual('bad_argument', self.tools.graphlet_export(
            h, 'json', what='everything', path='z.json')['error'])
        self.assertEqual('bad_argument', self.tools.graphlet_export(
            h, 'json', what='walk', path='z.json')['error'])
        self.assertFalse(os.path.exists(os.path.join(exports, 'z.json')))

    def test_compare_subtrie_and_lifecycle(self):
        # the fixture index has no manifest: a join of two retrievals is unverifiable
        a = self.fetch(SEED['fork'])['handle']
        b = self.fetch(SEED['fork'])['handle']
        cmp = self.tools.graphlet_compare(a, b)
        self.assertEqual(('unverifiable', None), (cmp['comparable'], cmp['equal']))
        self.tools.graphlet_free(a)
        self.tools.graphlet_free(b)
        self.client.fp = 'ab' * 32               # ... and with a manifest digest, comparable
        a = self.fetch(SEED['fork'])['handle']
        b = self.fetch(SEED['fork'])['handle']
        cmp = self.tools.graphlet_compare(a, b)
        self.assertEqual((True, True, 0), (cmp['comparable'], cmp['equal'], cmp['total']))
        sub = self.tools.graphlet_subtrie(a, [{'ref': 'h:0:1'}])
        self.assertEqual(a, sub['of'])
        view = self.store.view(sub['handle'])
        self.assertIsInstance(view, GraphletView)
        self.assertEqual({'left': [0], 'right': [0]}, view.path_ids)
        self.assertEqual(1, len(os.listdir(os.path.join(self.tmp.name, 'bodies'))))
        listed = {e['handle'] for e in self.tools.graphlet_list(max_bytes=8192)['rows']}
        self.assertEqual({a, b, sub['handle']}, listed)
        path = os.path.realpath(os.path.join(self.tmp.name, 'exports', 'saved.mgt'))
        self.assertEqual(path, self.tools.graphlet_save(a, 'saved.mgt')['path'])
        loaded = self.tools.graphlet_load('saved.mgt')
        self.assertIn('summary', loaded)
        self.assertEqual({'freed': b}, self.tools.graphlet_free(b))
        self.assertEqual('unknown_handle', self.tools.graphlet_summary(b)['error'])

    def test_continue(self):
        h = self.fetch(SEED['switch_chain'])['handle']
        out = self.tools.traverse_continue(h, 'right', 1, execute=False)
        self.assertEqual([{'sequence': CONTINUATION, 'labels': ['C']}], out['request']['seeds'])
        self.assertEqual({'handle': h, 'arm': 'right', 'walk': 1, 'overlap_bp': 20},
                         out['parent'])
        done = self.tools.traverse_continue(h, 'right', 1)
        self.assertRegex(done['handle'], r'^g_')
        self.assertEqual(CONTINUATION, self.client.requests[-1]['seeds'][0]['sequence'])

    def test_replay_after_expiry(self):
        h = self.fetch(SEED['merge'])['handle']
        self.clock.t += 1000
        self.store.sweep()
        out = self.tools.graphlet_summary(h)
        self.assertEqual('unknown_handle', out['error'])
        self.assertIn('traverse_fetch(replay=', out['hint'])
        again = self.tools.traverse_fetch(replay=h)        # the tombstone names the index
        self.assertIn('handle', again)
        self.assertEqual(self.client.requests[0], self.client.requests[-1])

    def test_replay_refuses_another_release(self):
        resp = T.doc_json('merge', 'graphlet')
        resp['release'] = 'r1'
        h = self.store.put(resp, {'seeds': [{'sequence': SEED['merge']}], 'strategy': {}},
                           source='mini')
        self.client.capabilities = lambda: {'release': 'r2'}
        out = self.tools.traverse_fetch(replay=h)
        self.assertEqual('release_mismatch', out['error'])

    def test_bad_arguments_are_results(self):
        h = self.fetch(SEED['merge'])['handle']
        self.assertEqual('bad_argument', self.tools.graphlet_walk(h, 'right', 7)['error'])
        self.assertEqual('bad_arm', self.tools.graphlet_walks(h, 'left')['error'])
        self.assertEqual('bad_argument',
                         self.tools.graphlet_sequence(h, 'right', walk=0, segment=1)['error'])

    def test_resolve(self):
        out = self.tools.traverse_resolve('mini', 'ACGTACGTAC', discover={'max_labels': 5})
        self.assertEqual((30, 10), (out['labels_total'], len(out['labels'])))
        self.assertEqual('unknown_index', self.tools.traverse_capabilities('nope')['error'])


class Raising(FakeClient):
    """A backend that fails every call with |exc|."""

    def __init__(self, exc):
        super().__init__()
        self.exc = exc

    def traverse_raw(self, request):
        raise self.exc

    def capabilities(self):
        raise self.exc

    def resolve(self, sequence, **kw):
        raise self.exc


class LongNames(FakeClient):
    """Responses whose label names are long paths (column names of a real deployment)."""

    def traverse_raw(self, request):
        resp = super().traverse_raw(request)
        r = resp['results'][0]
        if 'graphlet' in r:
            lines = r['graphlet'].split('\n')
            for i, line in enumerate(lines):
                if line.startswith('L '):
                    f = line.split(' ', 5)
                    lines[i] = ' '.join(f[:4] + ['0', '/data/refseq/release-224/' * 6
                                                 + f[5]])
            r['graphlet'] = '\n'.join(lines)
            r['graphlet_bytes'] = len(r['graphlet'].encode('utf-8'))
        return resp


class TestReviewFindings(ToolsCase):
    """The tool-layer findings of the graphlet review, one test each."""

    def handle(self, name, **kw):
        strategy = {'labels': {'mode': 'annotate'}} if name in ('annotate', 'same_name') \
            else None
        out = self.tools.traverse_fetch(seed={'sequence': seed_of(name)}, strategy=strategy,
                                        **kw)
        self.assertNotIn('error', out, out)
        return out['handle']

    def test_labels_at_a_position_are_the_labels_alive_there(self):
        h = self.handle('switch_chain')
        alive = {at: [r['name'] for r in self.tools.graphlet_labels(
            h, arm='right', at_bp=at, max_bytes=8192)['rows']] for at in (10, 50)}
        # A ends at 30 where B takes over: a segment's entry/end sets are not the answer
        self.assertEqual({10: ['A'], 50: ['B']}, alive)
        h = self.handle('reentry')
        self.assertEqual(['pB.fa'], [r['name'] for r in self.tools.graphlet_labels(
            h, arm='right', at_bp=45, max_bytes=8192)['rows']])

    def test_each_run_of_a_label_has_its_own_route(self):
        h = self.handle('merge')
        rows = self.tools.graphlet_labels(h, arm='right', name={'ref': 'c:3'},
                                          max_bytes=8192)['rows']
        self.assertEqual([(100, [0, 1, 3, 4, 6]), (36, [0, 2]), (62, [0, 1, 3, 5])],
                         [(r['to_bp'], r['route']) for r in rows])

    def test_a_derived_handle_answers_for_its_view(self):
        a = self.handle('fork')
        sub = self.tools.graphlet_subtrie(a, [{'ref': 'h:0:1'}])['handle']
        walks = self.tools.graphlet_walks(sub, 'right', max_bytes=8192)
        self.assertEqual((1, [0]), (walks['total'], [r['walk'] for r in walks['rows']]))
        self.assertEqual('for the selected labels', walks['evidence']['qualified'])
        self.assertEqual(a, walks['evidence']['view']['of'])
        s = self.tools.graphlet_summary(sub)['summary']
        self.assertEqual({'left': 1, 'right': 1}, {k: v['walks'] for k, v in s['arms'].items()})
        self.assertEqual('for the selected labels', s['arms']['right']['qualified'])
        claims = self.tools.graphlet_claims(sub, 'right', route_consistent=False,
                                            max_bytes=8192)
        self.assertEqual({'h:0:1'}, {r['label']['ref'] for r in claims['rows']})
        labels = self.tools.graphlet_labels(sub, max_bytes=8192)
        self.assertEqual({'h:0:1'}, {r['ref'] for r in labels['rows']})
        fasta = self.tools.graphlet_export(sub, 'fasta', path='v.fa')
        self.assertEqual(2, fasta['records'])           # one walk per arm, not four
        self.assertEqual('not_in_view', self.tools.graphlet_walk(sub, 'right', 1)['error'])
        for out in (self.tools.graphlet_export(sub, 'gfa', path='v.gfa'),
                    self.tools.graphlet_compare(sub, a),
                    self.tools.graphlet_subtrie(sub, ['acc1'])):
            self.assertEqual(('view_unsupported', a), (out['error'], out['of']))
        mgt = self.tools.graphlet_export(sub, 'mgt', path='v.mgt')
        self.assertIn('"view"', T.read(mgt['path']).split('\n')[1])

    def test_evidence_states_outcome_and_limitations(self):
        h = self.handle('annotate')                     # lists cut at 2: a lower bound
        g = self.store.graphlet(h)
        self.assertEqual('lower_bound', g.outcome.label_evidence)
        outs = [self.tools.graphlet_claims(h, 'right', max_bytes=8192),
                self.tools.graphlet_walks(h, 'right', max_bytes=8192),
                self.tools.graphlet_labels(h, max_bytes=8192),
                self.tools.graphlet_splits(h, 'right', max_bytes=8192),
                self.tools.graphlet_support(h, 'right', 0, max_bytes=8192),
                self.tools.graphlet_summary(h)]
        for out in outs:
            ev = out['evidence']
            self.assertFalse(ev['exact'])
            self.assertEqual('lower_bound', ev['outcome']['label_evidence'])
            self.assertIn('label_lists', ev['limitations']['right'])
        # a complete retrieval is exact, and an informational scope statement (no merge)
        # is no limitation
        ev = self.tools.graphlet_summary(self.handle('fork'))['evidence']
        self.assertEqual((True, {}), (ev['exact'], ev['limitations']))
        self.client.fp = 'ab' * 32
        a, b = self.handle('fork'), self.handle('fork')
        cmp = self.tools.graphlet_compare(a, b)
        self.assertEqual({'a', 'b'}, set(cmp['evidence']))
        self.assertEqual('complete', cmp['evidence']['b']['outcome']['label_evidence'])

    def test_failures_are_results(self):
        from urllib.error import URLError
        from metagraph.traverse import ServerInitializing, TraverseError
        cases = [(TraverseError(400, 'bad strategy'), {'error': 'backend_error', 'status': 400}),
                 (ServerInitializing(503, 'loading', retry_after=30),
                  {'error': 'backend_error', 'status': 503, 'retry_after_s': 30}),
                 (ConnectionRefusedError(61, 'Connection refused'),
                  {'error': 'backend_unreachable'}),
                 (URLError('refused'), {'error': 'backend_unreachable'})]
        h = self.handle('switch_chain')
        for exc, want in cases:
            tools = GraphletTools(self.store, {'mini': Raising(exc)})
            for out in (tools.traverse_fetch(seed={'sequence': SEED['fork']}),
                        tools.traverse_continue(h, 'right', 1),
                        tools.traverse_capabilities(),
                        tools.traverse_resolve(sequence='ACGTACGTAC')):
                self.assertEqual(want, {k: out.get(k) for k in want}, (exc, out))
        # malformed arguments
        bad = [self.tools.graphlet_walk(h, 'right', '0'),
               self.tools.graphlet_labels(h, name={'ref': 'c:0', 'name': 'A'}),
               self.tools.graphlet_labels(h, name={'tag': 'A'}),
               self.tools.graphlet_walks(h, 'right', n=0),
               self.tools.graphlet_summary(123),
               self.tools.graphlet_sequence(h, 'right', walk=0, start=-5),
               self.tools.graphlet_sequence(h, 'right', walk=0, start=5, end=2),
               self.tools.graphlet_sequence(h, 'right', walk=0, orientation='rc'),
               self.tools.graphlet_sequence(h, 'right', segment=99)]
        for out in bad:
            self.assertEqual('bad_argument', out['error'], out)
        self.assertEqual('bad_cursor', self.tools.graphlet_walks(h, 'right', cursor=123)['error'])
        out = self.tools.graphlet_labels(h, name='nope.fa')
        self.assertEqual('unknown_label', out['error'])
        self.assertIn('graphlet_labels', out['hint'])
        self.assertEqual('io_error', self.tools.graphlet_load('nonexistent.mgt')['error'])
        # a body-only file has no request to replay
        exports = os.path.join(self.tmp.name, 'exports')
        os.makedirs(exports, exist_ok=True)
        with open(os.path.join(exports, 'body.mgt'), 'w') as f:
            f.write(T.doc_text('fork'))
        loaded = self.tools.graphlet_load('body.mgt')['handle']
        self.assertEqual('not_replayable', self.tools.traverse_fetch(replay=loaded)['error'])

    def test_paths_are_confined(self):
        h = self.handle('fork')
        outside = os.path.join(self.tmp.name, 'outside')
        os.makedirs(outside)
        victim = os.path.join(outside, 'victim.txt')
        with open(victim, 'w') as f:
            f.write('keep me')
        exports = os.path.join(self.tmp.name, 'exports')
        os.makedirs(exports, exist_ok=True)
        os.symlink(outside, os.path.join(exports, 'link'))
        for path in ('../outside/victim.txt', victim, 'link/victim.txt', '/etc/x'):
            self.assertEqual('path_not_allowed',
                             self.tools.graphlet_export(h, 'fasta', path=path)['error'], path)
            self.assertEqual('path_not_allowed', self.tools.graphlet_save(h, path)['error'])
        self.assertEqual('keep me', T.read(victim))
        self.assertEqual('path_not_allowed', self.tools.graphlet_load('/etc/hosts')['error'])
        # the spool is readable: a body stored there loads
        body = os.path.join(self.tmp.name, 'bodies', os.listdir(
            os.path.join(self.tmp.name, 'bodies'))[0])
        self.assertIn('handle', self.tools.graphlet_load(body))
        # a parse error names the line, never what the file holds
        with open(os.path.join(exports, 'notes.mgt'), 'w') as f:
            f.write('## top secret line\n')
        out = self.tools.graphlet_load('notes.mgt')
        self.assertEqual('format_error', out['error'])
        self.assertNotIn('secret', out['message'])
        self.assertNotIn('##', out['message'])

    def assert_pages_fit(self, call, max_bytes):
        """Every page of a list tool is <= max_bytes -- a real page, not an error -- and
        the pages hold every row."""
        rows, out = [], call(max_bytes=max_bytes)
        while True:
            self.assertNotIn('error', out)
            self.assertLessEqual(self.size(out), max_bytes, out)
            rows.extend(out['rows'])
            if 'next_cursor' not in out:
                break
            out = call(max_bytes=max_bytes, cursor=out['next_cursor'])
        self.assertEqual(out['total'], len(rows))
        return rows

    def frame(self, call):
        """The least max_bytes a page can honour: the base, the longest cursor and the
        largest row with every cuttable field cut (scalars are never cut)."""
        from metagraph.traverse.mcp_tools import _shrink
        out = call(max_bytes=1 << 20)
        cursor = self.tools._cursor('x' * 20, 'g_' + '0' * 12, {}, out['total'])
        floor = max((_shrink(r, 0)[0] for r in out['rows']), key=self.size, default={})
        fields = sorted({k for r in out['rows'] for k, v in r.items()
                         if isinstance(v, (str, list))})
        return self.size(dict(out, rows=[floor], row_truncated=True, cut_fields=fields,
                              next_cursor=cursor))

    def test_every_return_fits_max_bytes(self):
        hs = {n: self.handle(n) for n in ('merge', 'annotate', 'switch_chain', 'fork')}
        h, ha = hs['merge'], hs['annotate']
        calls = [lambda **kw: self.tools.graphlet_walks(h, 'right', spell='full', **kw),
                 lambda **kw: self.tools.graphlet_walks(ha, 'left', spell='full', **kw),
                 lambda **kw: self.tools.graphlet_walk(ha, 'right', 0, **kw),
                 lambda **kw: self.tools.graphlet_support(ha, 'right', 0, **kw),
                 lambda **kw: self.tools.graphlet_labels(ha, **kw),
                 lambda **kw: self.tools.graphlet_labels(h, name='both.fa', **kw),
                 lambda **kw: self.tools.graphlet_claims(h, 'right', route_consistent=False,
                                                         **kw),
                 lambda **kw: self.tools.graphlet_splits(h, 'right', min_labels_before=0,
                                                         **kw),
                 lambda **kw: self.tools.graphlet_summary(ha, detail='limitations', **kw),
                 lambda **kw: self.tools.graphlet_compare(h, ha, arm='right', **kw)]
        for call in calls:
            # every budget from a page's fixed frame up: pages, never result_too_large
            # (the cursor's room was a guessed 120 B, 13-18 B short of real cursors)
            for budget in range(self.frame(call) + 40, 2600, 7):
                self.assert_pages_fit(call, budget)
        # graphlet_list pages like the other lists
        for _ in range(12):
            self.handle('fork')
        self.assertEqual(16, len(self.assert_pages_fit(self.tools.graphlet_list, 2048)))
        # the summary wrappers with long label names
        tools = GraphletTools(self.store, {'mini': LongNames()})
        for name in ('merge', 'annotate'):
            strategy = {'labels': {'mode': 'annotate'}} if name == 'annotate' else None
            out = tools.traverse_fetch(seed={'sequence': SEED[name]}, strategy=strategy)
            self.assertLessEqual(self.size(out), 2048)
            self.assertIn('summary', out)
            self.assertLessEqual(self.size(tools.graphlet_summary(out['handle'])), 2048)
            self.assertLessEqual(self.size(tools.graphlet_load(
                os.path.relpath(tools.graphlet_save(out['handle'], name + '.mgt')['path'],
                                os.path.join(self.tmp.name, 'exports')))), 2048)

    def test_filters_say_what_they_removed(self):
        h = self.handle('merge')
        # b.fa reaches 100 bp only along its own route through a merge
        w = self.tools.graphlet_walks(h, 'right', label={'ref': 'c:1'})
        self.assertEqual((0, {'route_only': 1}), (w['total'], w['filtered']))
        self.assertIn('route_consistent=false', w['hint'])
        w = self.tools.graphlet_walks(h, 'right', label={'ref': 'c:1'}, route_consistent=False)
        self.assertEqual((1, None), (w['total'], w.get('filtered')))
        c = self.tools.graphlet_claims(h, 'right', max_bytes=8192)
        # merge-entered claims (displayed from 62), not route_only ones (GPT review,
        # finding 11)
        self.assertEqual({'merge_entered': 2}, c['filtered'])
        every = self.tools.graphlet_claims(h, 'right', route_consistent=False, max_bytes=8192)
        self.assertEqual(c['total'] + 2, every['total'])

    def test_limitations_carry_what_the_agent_needs(self):
        resp = T.doc_json('limits', 'graphlet')
        h = self.store.put(resp, None)
        caveats = self.tools.graphlet_summary(h, max_bytes=8192)['summary']['caveats']
        walk = [c for c in caveats if c['kind'] == 'walk_domain'][0]
        self.assertEqual(('bounds.max_steps', 11, 12, 7),
                         (walk['knob'], walk['limit'], walk['observed'], walk['complete_to_bp']))
        self.assertIn('effect', walk)
        stop = [c for c in caveats if c['kind'] == 'resource_stop'][0]
        self.assertEqual(['narrow_seed', 'constrain_labels'], stop['actions'])
        full = self.tools.graphlet_summary(h, detail='limitations', max_bytes=8192)
        self.assertEqual(len(self.store.graphlet(h).limitations), full['total'])
        self.assertEqual(['narrow_seed', 'constrain_labels'], full['resource_stop']['actions'])
        fork = self.tools.graphlet_summary(self.handle('fork'))['summary']['caveats']
        self.assertEqual([True, True], [c.get('informational') for c in fork
                                        if c['kind'] == 'scope'])

    def test_cut_lists_report_their_true_totals(self):
        h = self.handle('annotate')
        rows = self.tools.graphlet_walk(h, 'right', 0, max_bytes=8192)['rows']
        self.assertEqual((2, 4, 2, 4, False),
                         tuple(rows[0][k] for k in ('labels_in', 'labels_in_total', 'labels_out',
                                                    'labels_out_total', 'exact')))
        rows = self.tools.graphlet_support(h, 'right', 0, max_bytes=8192)['rows']
        self.assertEqual((2, 4, False), (rows[0]['n_labels'], rows[0]['n_labels_total'],
                                         rows[0]['exact']))

    def test_cursors_and_parents_survive_a_restart(self):
        h = self.handle('annotate')
        first = self.tools.graphlet_labels(h, n=1)
        sc = self.handle('switch_chain')
        new = self.tools.traverse_continue(sc, 'right', 1)
        parent = {'handle': sc, 'arm': 'right', 'walk': 1, 'overlap_bp': 20}
        self.assertEqual(parent, new['parent'])
        restarted = GraphletTools(GraphletStore(self.tmp.name, clock=self.clock,
                                                ttl_disk_s=100), {'mini': self.client})
        second = restarted.graphlet_labels(h, n=1, cursor=first['next_cursor'])
        self.assertNotIn('error', second)
        self.assertNotEqual(first['rows'], second['rows'])
        self.assertEqual(parent, restarted.graphlet_summary(new['handle'])['parent'])
        listed = [e for e in restarted.graphlet_list(max_bytes=8192)['rows']
                  if e['handle'] == new['handle']][0]
        self.assertEqual((parent, sc), (listed['parent'], listed['derived_from']))


if __name__ == '__main__':
    unittest.main()
