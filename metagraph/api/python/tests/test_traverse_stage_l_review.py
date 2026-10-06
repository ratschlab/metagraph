"""Stage L, the review's confirmed findings (one regression test each; the stage-L tests
themselves are test_traverse_stage_l.py):

  R1  graphlet_walks charges a cached ranking as the ranking itself charges (L1: the same
      call charges the same units cold and warm, and a later page as on a fresh server);
      every budgeted tool charges the same twice in a row
  R2  graphlet_sequence keeps room for the local block: a long walk is a truncated slice
      with next_from under local limits too, never result_too_large
  R3  next_request(budget={...}) is the keyword override it was before stage L, also
      through traverse_continue; overrides of the tool's own arguments are bad_argument
  R4  a resume token binds the graphlet's body: refused on another body that shares the
      seed and the counts, accepted on the same body parsed again
  R5  graphlet_export(format=json) under local limits never makes the whole text
  R6  a parse on demand runs under the parse class raised by the call's budget argument;
      the store's own parse limits, which no call raises, are not offered as a lever
  R7  a body-level copy of an unparsed entry says in its receipt that it is unvalidated
  R8  label_walks resolves the label before the arm (UnknownLabel with a bad arm too)
  R9  the first label lookups on a fresh model scan the labels (ops._LABEL_SCANS of
      them); the index is built only after them, and both answer as a scan
  R10 (found auditing R2's other ceiling-fitted results) a budgeted receipt is checked
      with its local block before the operation runs: never a handle made and the answer
      result_too_large
"""

import collections
import copy
import itertools
import json
import os
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402
from test_traverse_scale import comb_annotate  # noqa: E402
from test_traverse_stage_l import BUDGETED_CALLS, ITEMS, _real, _tools  # noqa: E402

import metagraph.traverse as MT  # noqa: E402
from metagraph.traverse import (  # noqa: E402
    AmbiguousLabel, GraphletStore, LocalBudget, LocalBudgetExceeded, LocalLimits,
    UnknownLabel, derive, local_budget, ops,
)
from metagraph.traverse import mcp_tools  # noqa: E402
from metagraph.traverse.mcp_tools import (  # noqa: E402
    MIN_MAX_BYTES, GraphletTools, ToolLimits, _size,
)

SECRET = b'k' * 32


def fresh(item):
    return MT.from_response(item[1], item[2])


def put(st, name):
    item = _real(name)
    g = fresh(item)
    return st.put_parsed(g, item[1]['graphlet'], {'seeds': []}), g


def _usage(out):
    return out['local']['usage']


# ======================================================================= R1

class TestR1CachedRanking(unittest.TestCase):
    NAMES = ('sra_hub__strict', 'mini_ndm1__annotate_exh', 'switch_chain', 'merge')

    def test_the_same_units_cold_and_warm(self):
        for name in self.NAMES:
            with tempfile.TemporaryDirectory() as d:
                st, t = _tools(d)
                h, g = put(st, name)
                for s in g.arms:
                    calls = [lambda: t.graphlet_walks(h, s, n=3)]
                    if g.labels:
                        # a label filter ranks twice (route_consistent and every walk)
                        # and counts the walks it removed
                        calls.append(lambda: t.graphlet_walks(
                            h, s, n=2, label={'ref': g.labels[0].ref}))
                    for call in calls:
                        got = [call() for _ in range(3)]
                        self.assertEqual(_usage(got[0]), _usage(got[1]), (name, s))
                        self.assertEqual(_usage(got[0]), _usage(got[2]), (name, s))
                        # the cold usage as the budget: complete on every call
                        u = _usage(got[0])
                        lim = {'work_units': u['work_units'],
                               'memory_mb': u['memory_bytes'] / (1 << 20)}
                        for _ in range(2):
                            out = t.graphlet_walks(h, s, n=3, budget=lim) \
                                if call is calls[0] else t.graphlet_walks(
                                    h, s, n=2, label={'ref': g.labels[0].ref}, budget=lim)
                            self.assertNotIn('error', out, (name, s, out))
                            self.assertTrue(out['local']['complete'], (name, s))

    def test_a_stop_cold_and_warm_is_the_same_answer(self):
        with tempfile.TemporaryDirectory() as d:
            st, t = _tools(d)
            h, g = put(st, 'sra_hub__strict')
            s = list(g.arms)[0]
            u = _usage(t.graphlet_walks(h, s, n=3))['work_units']
            for limit in (u - 1, u // 2, u // 5):
                cold = GraphletTools(st, secret=SECRET, export_dir=os.path.join(d, 'x'),
                                     local_limits=ToolLimits())
                first = cold.graphlet_walks(h, s, n=3, budget={'work_units': limit})
                again = cold.graphlet_walks(h, s, n=3, budget={'work_units': limit})
                warm = t.graphlet_walks(h, s, n=3, budget={'work_units': limit})
                self.assertEqual(first, again, limit)
                self.assertEqual(first, warm, limit)

    def test_a_later_page_as_on_a_fresh_server(self):
        with tempfile.TemporaryDirectory() as d:
            st, t = _tools(d)
            h, g = put(st, 'sra_hub__strict')
            s = 'right'
            p1 = t.graphlet_walks(h, s, n=2)
            self.assertIn('next_cursor', p1)
            warm = t.graphlet_walks(h, s, n=2, cursor=p1['next_cursor'])
            other = GraphletTools(st, secret=SECRET, export_dir=os.path.join(d, 'x'),
                                  local_limits=ToolLimits())
            cold = other.graphlet_walks(h, s, n=2, cursor=p1['next_cursor'])
            self.assertEqual(warm, cold)

    def test_every_budgeted_tool_charges_the_same_twice(self):
        for name in ('sra_hub__strict', 'mini_ndm1__annotate_exh'):
            with tempfile.TemporaryDirectory() as d:
                st, t = _tools(d)
                h, g = put(st, name)
                s = list(g.arms)[0]
                for tool, call in BUDGETED_CALLS.items():
                    a, b = call(t, h, s), call(t, h, s)
                    self.assertNotIn('error', a, (name, tool, a))
                    self.assertEqual(_usage(a), _usage(b), (name, tool))


# ======================================================================= R2

class TestR2SequenceSlice(unittest.TestCase):
    def test_a_long_walk_under_local_limits(self):
        text = comb_annotate(20000)
        res = {'graphlet': text, 'seed': {'seed_id': 'comb'}, 'annotation': {}}
        resp = {'strategy': {'direction': 'right', 'labels': {'mode': 'annotate'}},
                'results': [res]}
        with tempfile.TemporaryDirectory() as d:
            st = GraphletStore(os.path.join(d, 'sp'))
            g = MT.from_response(res, resp)
            h = st.put_parsed(g, text, {'seeds': []})
            spine = max(derive.paths(g.arms['right']), key=lambda p: p.length_bp)
            whole = ops.spell(g, 'right', spine.id)
            self.assertGreater(len(whole), 16384)
            plain = GraphletTools(st, secret=SECRET)
            budgeted = GraphletTools(st, secret=SECRET, local_limits=ToolLimits())
            for t in (plain, budgeted):
                got, start = '', 0
                while True:
                    out = t.graphlet_sequence(h, 'right', walk=spine.id, start=start)
                    self.assertNotIn('error', out, out)
                    self.assertLessEqual(_size(out), t.sequence_max_bytes)
                    self.assertEqual(t is budgeted, 'local' in out)
                    got += out['sequence']
                    if not out.get('truncated'):
                        break
                    self.assertEqual(out['to'], out['next_from'])
                    self.assertGreater(out['next_from'], start)
                    start = out['next_from']
                self.assertEqual(whole, got)

    def test_a_small_ceiling_on_a_real_walk(self):
        item = _real('mini_ndm1__trace')
        with tempfile.TemporaryDirectory() as d:
            st = GraphletStore(os.path.join(d, 'sp'))
            g = fresh(item)
            h = st.put_parsed(g, item[1]['graphlet'], {'seeds': []})
            side, p = max(((s, p) for s, a in g.arms.items() for p in derive.paths(a)),
                          key=lambda x: x[1].length_bp)
            t = GraphletTools(st, secret=SECRET, sequence_max_bytes=4096,
                              local_limits=ToolLimits())
            out = t.graphlet_sequence(h, side, walk=p.id)
            self.assertNotIn('error', out, out)
            self.assertTrue(out['truncated'])
            self.assertLessEqual(_size(out), 4096)
            self.assertTrue(out['local']['complete'])


# ======================================================================= R3

class TestR3BudgetOverride(unittest.TestCase):
    def setUp(self):
        self.item = _real('mini_ndm1__annotate_exh')
        self.g = fresh(self.item)
        self.side = 'left'
        a = self.g.arms[self.side]
        self.pid = next(p.id for p in derive.paths(a)
                        if a.segments[p.leaf].leaf.continuation is not None)
        self.override = {'max_work_units': 5}

    def expected(self, value):
        base = self.g.next_request(self.side, [self.pid])
        strategy = copy.deepcopy(base['strategy'])
        ops._merge_overrides(strategy, {'budget': value})
        return strategy

    def test_the_library(self):
        g, s, pid = self.g, self.side, self.pid
        want = self.expected(self.override)
        self.assertEqual(want, g.next_request(s, [pid], budget=self.override)['strategy'])
        self.assertEqual(want, ops.next_request(g, s, [pid], budget=self.override)['strategy'])
        self.assertEqual(want, g.next_requests(s, [pid], budget=self.override)[0]['strategy'])
        # None given explicitly is the override too, as it always was
        self.assertIsNone(g.next_request(s, [pid], budget=None)['strategy']['budget'])
        # a LocalBudget charges the call and overrides nothing
        b = LocalBudget()
        self.assertEqual(dict(g.next_request(s, [pid])), dict(g.next_request(s, [pid],
                                                                            budget=b)))
        self.assertGreater(b.used_work, 0)
        b2 = LocalBudget()
        g.next_requests(s, [pid], budget=b2)
        self.assertGreater(b2.used_work, 0)
        # both at once: the local budget ambient
        with local_budget() as lb:
            self.assertEqual(want, g.next_request(s, [pid], budget=self.override)['strategy'])
        self.assertGreater(lb.used_work, 0)
        # LocalLimits is meant as a local budget: refused as by every other operation
        with self.assertRaises(TypeError):
            g.next_request(s, [pid], budget=LocalLimits())

    def test_traverse_continue(self):
        want = self.expected(self.override)
        for limits in (None, ToolLimits()):
            with tempfile.TemporaryDirectory() as d:
                st = GraphletStore(os.path.join(d, 'sp'))
                h = st.put_parsed(fresh(self.item), self.item[1]['graphlet'], {'seeds': []})
                t = GraphletTools(st, secret=SECRET, local_limits=limits)
                out = t.traverse_continue(h, self.side, self.pid, execute=False,
                                          overrides={'budget': self.override})
                self.assertNotIn('error', out, out)
                self.assertEqual(want, out['request']['strategy'])
                self.assertEqual(limits is not None, 'local' in out)
                out = t.traverse_continue(h, self.side, self.pid, execute=False,
                                          overrides={'bp': 77})
                self.assertEqual(77, out['request']['strategy']['bounds']['max_extension_bp'])
                for key in ('reset_branches', 'arm', 'leaves', 'self', 'g'):
                    out = t.traverse_continue(h, self.side, self.pid, execute=False,
                                              overrides={key: 1})
                    self.assertEqual('bad_argument', out.get('error'), (key, out))
                    self.assertIn(key, out['message'])


# ======================================================================= R4

def _old_fingerprint(g):
    """What a resume token was bound to before: the seed's first bases and the counts."""
    return json.dumps([g.seed.sequence[:64], len(g.labels),
                       [(s, a.counts.segments, a.counts.runs) for s, a in g.arms.items()]])


def _token(g):
    """The resume token of a stopped claims(strict=False) of |g| (None: none found)."""
    b = LocalBudget()
    g.claims(strict=False, budget=b)
    for f in (0.9, 0.7, 0.5, 0.3, 0.1):
        try:
            g.claims(strict=False, budget=LocalBudget(work_units=int(b.used_work * f)))
        except LocalBudgetExceeded as e:
            if e.partial is not None and e.partial.resume is not None:
                return e.partial
    return None


class TestR4TokenBinding(unittest.TestCase):
    def test_bodies_that_share_the_old_fingerprint(self):
        groups = collections.defaultdict(list)
        for item in ITEMS:
            g = fresh(item)
            groups[_old_fingerprint(g)].append((item[0], g))
        pairs = [(a, b) for xs in groups.values() for a, b in itertools.combinations(xs, 2)
                 if a[1].dump(envelope=False) != b[1].dump(envelope=False)]
        # the review found five pairs, two of them among the committed real and review3
        # fixtures (which no generator rewrites)
        self.assertGreaterEqual(len(pairs), 2)
        tried = 0
        for (ka, A), (kb, B) in pairs:
            self.assertNotEqual(A.body_digest, B.body_digest, (ka, kb))
            part = _token(A)
            if part is None:
                continue
            tried += 1
            with self.assertRaises(ValueError, msg=(ka, kb)):
                B.claims(strict=False, resume=part.resume)
            # the same body parsed again takes it, and resuming gives the whole list
            again = MT.parse(A.dump(envelope=True))
            self.assertEqual(A.body_digest, again.body_digest)
            self.assertEqual(A.claims(strict=False),
                             part.rows + again.claims(strict=False, resume=part.resume))
        self.assertGreater(tried, 0)

    def test_the_digest_is_the_body_and_survives_save_and_load(self):
        for item in ITEMS[::3]:
            g = fresh(item)
            self.assertEqual(MT.parse(item[1]['graphlet']).body_digest, g.body_digest)
            with tempfile.TemporaryDirectory() as d:
                p = os.path.join(d, 'g.mgt')
                g.save(p)
                self.assertEqual(g.body_digest, MT.load(p).body_digest, item[0])
                st = GraphletStore(os.path.join(d, 'sp'))
                h = st.put_parsed(g, item[1]['graphlet'], {'seeds': []}, resident=False)
                self.assertEqual(g.body_digest, st.graphlet(h).body_digest, item[0])


# ======================================================================= R5

class TestR5JsonExport(unittest.TestCase):
    def test_no_whole_text_under_local_limits(self):
        real_dumps = json.dumps
        made = []

        def spy(obj, *a, **k):
            s = real_dumps(obj, *a, **k)
            made.append(len(s))
            return s

        for name in ('sra_hub__strict', 'mini_ndm1__annotate_exh'):
            with tempfile.TemporaryDirectory() as d:
                st, t = _tools(d)
                h, g = put(st, name)
                plain = GraphletTools(st, secret=SECRET, export_dir=os.path.join(d, 'x'))
                made.clear()
                with mock.patch.object(mcp_tools.json, 'dumps', spy):
                    ref = plain.graphlet_export(h, format='json', path='plain.json')
                n = os.path.getsize(ref['path'])
                self.assertIn(n, made)                  # the spy sees the whole text
                made.clear()
                with mock.patch.object(mcp_tools.json, 'dumps', spy):
                    out = t.graphlet_export(h, format='json', path='budgeted.json')
                self.assertTrue(out['local']['complete'])
                self.assertLess(max(made, default=0), n // 2, name)
                with open(ref['path'], 'rb') as f1, open(out['path'], 'rb') as f2:
                    self.assertEqual(f1.read(), f2.read())


# ======================================================================= R6

class TestR6ParseOnDemand(unittest.TestCase):
    def stored(self, d, store_limits=None):
        item = _real('sra_hub__strict')
        # never resident (no RAM for a model): parsed again on every call
        st = GraphletStore(os.path.join(d, 'sp'), max_ram_mb=0, parse_limits=store_limits)
        g = fresh(item)
        h = st.put_parsed(g, item[1]['graphlet'], {'seeds': []}, resident=False)
        return st, h

    def test_the_call_budget_raises_the_parse_class(self):
        with tempfile.TemporaryDirectory() as d:
            st, h = self.stored(d)
            t = GraphletTools(st, secret=SECRET,
                              local_limits=ToolLimits(parse=LocalLimits(50_000, 1024)))
            out = t.graphlet_walks(h, 'right', n=3)
            self.assertEqual('local_budget_exceeded', out['error'])
            self.assertEqual(('parse', 50_000), (out['stop']['op'], out['stop']['limit']))
            self.assertIn('raise_local_budget', out['stop']['actions'])
            self.assertEqual(50_000, out['local']['parsed']['limits']['work_units'])
            # the lever works: the call's budget raises the parse too
            out = t.graphlet_walks(h, 'right', n=3, budget={'work_units': 10 ** 8})
            self.assertNotIn('error', out, out)
            self.assertEqual(10 ** 8, out['local']['parsed']['limits']['work_units'])
            self.assertEqual(1024, out['local']['parsed']['limits']['memory_mb'])
            # held to the parse class's ceiling, and never lowered below its default
            out = t.graphlet_walks(h, 'right', n=3, budget={'work_units': 10 ** 13})
            self.assertEqual(ToolLimits().parse_ceiling.work_units,
                             out['local']['parsed']['limits']['work_units'])
            t2 = GraphletTools(st, secret=SECRET, local_limits=ToolLimits())
            out = t2.graphlet_walks(h, 'right', n=3, budget={'work_units': 300})
            self.assertEqual('local_budget_exceeded', out['error'])
            self.assertNotEqual('parse', out['stop']['op'])
            self.assertEqual(ToolLimits().parse.work_units,
                             out['local']['parsed']['limits']['work_units'])

    def test_the_store_limits_are_no_lever(self):
        with tempfile.TemporaryDirectory() as d:
            st, h = self.stored(d, LocalLimits(work_units=50_000))
            t = GraphletTools(st, secret=SECRET, local_limits=ToolLimits())
            for budget in (None, {'work_units': 10 ** 8}):
                out = t.graphlet_walks(h, 'right', n=3, budget=budget)
                self.assertEqual('local_budget_exceeded', out['error'])
                self.assertEqual('parse', out['stop']['op'])
                self.assertNotIn('raise_local_budget', out['stop']['actions'])
                self.assertIn('export_mgt', out['stop']['actions'])
                self.assertIn("store's parse limits", out['message'])
                self.assertEqual(50_000, out['local']['parsed']['limits']['work_units'])

    def test_a_stop_selecting_a_stored_view_is_the_call_s_own(self):
        with tempfile.TemporaryDirectory() as d:
            st, t = _tools(d)
            h, g = put(st, 'sra_hub__strict')
            v = t.graphlet_subtrie(h, [{'id': 0}])['handle']
            out = t.graphlet_summary(v, budget={'work_units': 1})
            self.assertEqual('local_budget_exceeded', out['error'])
            self.assertNotEqual('parse', out['stop']['op'])
            self.assertNotIn('parsed on demand', out['message'])
            self.assertIn('raise_local_budget', out['stop']['actions'])


# ======================================================================= R7

class TestR7UnvalidatedCopies(unittest.TestCase):
    def test_receipts_of_an_unparsed_entry(self):
        k, r, resp = _real('sra_hub__strict')
        lines = r['graphlet'].split('\n')
        # one R record's end token made unknown: the frame (line count, H, Z) intact
        i = next(j for j, line in enumerate(lines) if line.startswith('R '))
        f = lines[i].split(' ')
        f[5] = '#bogus#'
        lines[i] = ' '.join(f)
        r2 = dict(r, graphlet='\n'.join(lines))
        r2.pop('graphlet_bytes', None)
        resp2 = dict(resp, results=[r2])
        with self.assertRaises(MT.GraphletFormatError):
            MT.from_response(r2, resp2)
        with tempfile.TemporaryDirectory() as d:
            st, t = _tools(d, store_limits=LocalLimits(work_units=1000))
            with self.assertRaises(LocalBudgetExceeded) as cm:
                st.put(resp2, {'seeds': []})
            bad = cm.exception.handle
            st2 = GraphletStore(os.path.join(d, 'sp2'))
            good = st2.put_parsed(fresh((k, r, resp)), r['graphlet'], {'seeds': []})
            t2 = GraphletTools(st2, secret=SECRET, export_dir=os.path.join(d, 'x'),
                               local_limits=ToolLimits())
            for out in (t.graphlet_export(bad, format='mgt', path='e.mgt'),
                        t.graphlet_save(bad, 's.mgt')):
                self.assertTrue(out['local']['complete'])
                self.assertIs(False, out['parsed'])
                self.assertIn('validates it', out['validation'])
            for out in (t2.graphlet_export(good, format='mgt', path='g.mgt'),
                        t2.graphlet_save(good, 'gs.mgt')):
                self.assertNotIn('parsed', out)
                self.assertNotIn('validation', out)
            # under a ceiling too small for the note, `parsed` stays and the note is cut
            small = GraphletTools(st, secret=SECRET, export_dir=os.path.join(d, 'x'),
                                  local_limits=ToolLimits(), max_bytes=256)
            out = small.graphlet_save(bad, 's2.mgt')
            self.assertIs(False, out['parsed'])
            self.assertIn('validation', out['fields_cut'])
            self.assertLessEqual(_size(out), 256)


# ======================================================================= R8

class TestR8LabelBeforeArm(unittest.TestCase):
    def test_unknown_label_with_a_bad_arm(self):
        for name in ('sra_hub__strict', 'switch_chain'):
            g = fresh(_real(name))
            other = 'left' if 'left' not in g.arms else 'bogus'
            for sel in ('no-such-label', {'id': 10 ** 6}, {'name': 'nothing'}):
                for arm in ('bogus', other):
                    for fn in (g.label_walks, lambda *a: ops.label_walks_routes(g, *a)):
                        with self.assertRaises(UnknownLabel):
                            fn(sel, arm)
                        with self.assertRaises(UnknownLabel):
                            with local_budget():
                                fn(sel, arm)
            # a good label with a bad arm is still the arm's error
            with self.assertRaises((KeyError, ValueError)):
                g.label_walks({'id': 0}, 'bogus')


# ======================================================================= R9

class TestR9FirstLookupScans(unittest.TestCase):
    @staticmethod
    def reference(g, sel):
        if isinstance(sel, str):
            return [l for l in g.labels if l.ref == sel or l.name == sel]
        (k, v), = sel.items()
        return [l for l in g.labels if isinstance(v, str) and getattr(l, k) == v]

    def answer(self, g, sel):
        try:
            return ('ok', ops.label(g, sel))
        except AmbiguousLabel as e:
            if isinstance(e, UnknownLabel):
                return ('unknown', str(e))
            return ('ambiguous', str(e))

    def test_the_first_lookups_build_no_index(self):
        g = fresh(_real('sra_hub__strict'))
        self.assertGreaterEqual(ops._LABEL_SCANS, 1)
        for i in range(ops._LABEL_SCANS):
            ops.label(g, g.labels[i].name if i % 2 else {'ref': g.labels[i].ref})
            self.assertNotIn('label_index', g.cache)
        ops.label(g, g.labels[4].ref)
        self.assertIn('label_index', g.cache)

    def test_scan_and_index_answer_alike(self):
        n = 0
        for item in ITEMS[::2]:
            g = fresh(item)
            sels = ['no-such', {'name': 5}, {'ref': None}, {'name': 'no-such'}]
            for l in g.labels[:30]:
                sels += [l.name, l.ref, {'name': l.name}, {'ref': l.ref}]
            for sel in sels:
                g.cache.clear()
                scanned = self.answer(g, sel)
                self.assertNotIn('label_index', g.cache)
                g.cache['label_scans'] = ops._LABEL_SCANS
                indexed = self.answer(g, sel)
                self.assertIn('label_index', g.cache)
                self.assertEqual(scanned, indexed, (item[0], sel))
                want = self.reference(g, sel)
                if len(want) == 1:
                    self.assertEqual(('ok', want[0]), scanned)
                else:
                    self.assertEqual('unknown' if not want else 'ambiguous', scanned[0])
                n += 1
        self.assertGreater(n, 500)

    def test_a_shared_name(self):
        g = T.body('same_name')
        for sel in ('ACC1', {'name': 'ACC1'}):
            g.cache.clear()
            first = self.answer(g, sel)
            self.assertEqual('ambiguous', first[0])
            g.cache['label_scans'] = ops._LABEL_SCANS
            self.assertEqual(first, self.answer(g, sel))
            self.assertIn('label_index', g.cache)


# ======================================================================= R10

class TestR10ReceiptsUnderTheCeiling(unittest.TestCase):
    def test_every_ceiling_near_the_least_receipt(self):
        k, r, resp = _real('merge')
        with tempfile.TemporaryDirectory() as d:
            st = GraphletStore(os.path.join(d, 'sp'))
            h = st.put_parsed(fresh((k, r, resp)), r['graphlet'], {'seeds': []})
            made = 0
            for mb in range(MIN_MAX_BYTES, 360):
                t = GraphletTools(st, secret=SECRET, export_dir=os.path.join(d, 'x'),
                                  local_limits=ToolLimits(), max_bytes=mb)
                for name, fn in (
                        ('save', lambda: t.graphlet_save(h, 's.mgt')),
                        ('mgt', lambda: t.graphlet_export(h, format='mgt', path='e.mgt')),
                        ('fasta', lambda: t.graphlet_export(h, format='fasta', path='e.fa')),
                        ('subtrie', lambda: t.graphlet_subtrie(h, [{'id': 0}]))):
                    n = len(st.list())
                    out = fn()
                    self.assertNotEqual('result_too_large', out.get('error'), (name, mb))
                    self.assertLessEqual(_size(out), mb, (name, mb))
                    if 'error' in out:
                        self.assertEqual('receipt_too_large', out['error'], (name, mb))
                        self.assertEqual(n, len(st.list()), (name, mb))
                    else:
                        made += 1
                        self.assertTrue('path' in out or 'handle' in out, (name, mb))
            self.assertGreater(made, 100)

    def test_an_interrupted_receipt(self):
        key, r, resp = _real('sra_hub__strict')
        for mb in range(MIN_MAX_BYTES, 420, 3):
            with tempfile.TemporaryDirectory() as d:
                st, t = _tools(d, ToolLimits(parse=LocalLimits(2000, 1024)), resp=resp)
                out = t.traverse_fetch('idx', seed={'sequence': 'ACGT'}, max_bytes=mb)
                self.assertNotEqual('result_too_large', out.get('error'), mb)
                self.assertLessEqual(_size(out), mb, mb)
                if 'handle' in out:
                    self.assertFalse(out['local']['complete'], mb)
                    self.assertFalse(st.get(out['handle']).parsed)
                else:
                    self.assertEqual([], st.list(), (mb, out))


if __name__ == '__main__':
    unittest.main()
