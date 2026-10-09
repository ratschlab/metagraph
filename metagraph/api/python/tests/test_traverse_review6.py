"""The external review of levels 4-5 (pinned at dcc0cebd), the library's part, and the open
items of the level-5 library batch (bf8dcee4):

  * finding 2: a budgeted operation is admitted before the derived caches it reads are
    built -- spell(), support_profile(), support_changes(), walks_at(), continuation
    and the ranking resolved a walk (building the arm's paths, a Path per walk) before
    their budget was charged, so a refused zero-budget call on comb_annotate(20000)
    left 20,001 Path objects (2.7 MB) behind and reported 0 work / 0 memory (the
    reviewer's /tmp/review_l5_admission_probe.py). Stated more widely: stopped at ANY
    charge point, a call leaves only derivations it charged (an audit of every operation
    also found path_of_leaf built uncharged under annotate routes, segment_ops under a
    constrain comparison's cuts, and the name index under prefix_subset);
  * finding 4: compare_cost()'s at_least is never above what the completed comparison
    charges -- annotate claims were priced from the whole DAG's ends, not the ends of
    the DAG restricted to the common certified depth (comb 1 vs comb 100: 9,197 lwu
    advertised, 7,963 charged; /tmp/review_l5_cost_probe.py), and labels mode charged
    its keys at their upper bound into at_least (above a comparison with labels=);
  * finding 5: graphlet_load returns the handle of the entry it stored -- it parsed the
    stored entry again under the store's parse limits, which refused it, and answered
    local_budget_exceeded without the handle (GraphletStore(max_ram_mb=0,
    parse_limits=LocalLimits(work_units=1)), /tmp/review_l5_load_receipt_probe.py); a
    stop on the store's parse limits offers no raise_local_budget (no call's budget
    raises them);
  * the batch's memory gaps: the stage-L account of compare() is at least its traced
    peak on every pair of retrievals of one seed, both orders, every mode (it was 0.82x
    on a cross pair under prefix_subset: the name index uncharged, the index itself
    2.5-2.9x the price it was charged at, the cut's segment pieces read whole).

The review of these fixes (the working tree on bf8dcee4), findings 1-3 (4 and 5 are in
test_traverse_vendored_review.py, beside VMD-04 and VMD-02):

  * finding 1: the J-line charge of standalone_text() (graphlet_export(format=mgt),
    graphlet_save) is the same on every run -- it was sized by the entry file, whose
    created and accessed times print to a length that follows the clock;
  * finding 2: to_fasta(leaves=[]) resolves no walk, so a refused call builds no paths;
  * finding 3: graphlet_load states a stop of the caller's ambient budget as that
    budget's (raise_local_budget offered), not as one on the store's parse limits; the
    same classification for graphlet_save, which parses the stored body on demand.
"""

import collections
import gc
import glob
import os
import sys
import tempfile
import tracemalloc
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402
import traverse_golden as TG  # noqa: E402
from test_traverse_scale import comb_annotate  # noqa: E402

import metagraph.traverse as MT  # noqa: E402
from metagraph.traverse import (  # noqa: E402
    GraphletStore, LocalBudget, LocalBudgetExceeded, LocalLimits, derive, ops, parse,
)
from metagraph.traverse import budget as B  # noqa: E402
from metagraph.traverse.mcp_tools import GraphletTools, ToolLimits  # noqa: E402
from metagraph.traverse.model import Label, Path  # noqa: E402


def _items():
    files = sorted(glob.glob(os.path.join(T.DOCS, '**', '*.graphlet.json'), recursive=True))
    return TG.items_of(files + TG.golden_files())


ITEMS = _items()


def _fresh(item):
    return MT.from_response(item[1], item[2])


def _real(name):
    return next(x for x in ITEMS if os.path.basename(x[0]).startswith(name + '.'))


def _caches(g):
    out = {('g', k) for k in g.cache}
    for side, a in g.arms.items():
        out |= {(side, k) for k in a.cache}
    return out


# the derivations with a price (derive._derivation_prices) and the graphlet-level ones
PRICED = set(derive._derivation_prices(collections.defaultdict(lambda: 1), 1, 'constrain'))
G_PRICED = {'label_index', 'name_counts', 'label_summary'}


def _ops(g):
    """name -> fn(budget=...) of the budgeted operations on |g| (each a fresh call)."""
    t = {'to_json': lambda **k: g.to_json(**k), 'to_fasta': lambda **k: g.to_fasta(**k),
         'to_fasta_walk': lambda **k: g.to_fasta(orientation='walk', with_seed=False,
                                                 width=7, **k),
         'to_gfa': lambda **k: g.to_gfa(**k), 'dump': lambda **k: g.dump(**k),
         'summary': lambda **k: g.summary(**k),
         'label_summary': lambda **k: g.label_summary(**k),
         'claims': lambda **k: g.claims(strict=False, **k)}
    if g.labels:
        t['subgraph'] = lambda **k: g.subgraph([{'id': 0}], **k)
        t['claims_label'] = lambda **k: g.claims(labels=[{'id': 0}], strict=False, **k)
    for s in g.arms:
        t[s + ':walks'] = lambda s=s, **k: g.walks(s, **k)
        t[s + ':walks_id'] = lambda s=s, **k: g.walks(s, by='id', **k)
        t[s + ':rank'] = lambda s=s, **k: ops.rank_walks(g, s, **k)
        t[s + ':walks_at'] = lambda s=s, **k: ops.walks_at(g, s, [0], with_claims=True, **k)
        t[s + ':splits'] = lambda s=s, **k: g.splits(s, **k)
        t[s + ':spell'] = lambda s=s, **k: g.spell(s, 0, **k)
        t[s + ':spell_walk'] = lambda s=s, **k: g.spell(s, 0, 'walk', True, **k)
        t[s + ':profile'] = lambda s=s, **k: g.support_profile(s, 0, **k)
        t[s + ':profile_route'] = lambda s=s, **k: g.support_profile(s, 0, 'route', **k)
        t[s + ':changes'] = lambda s=s, **k: g.support_changes(s, 0, **k)
        t[s + ':continuation'] = lambda s=s, **k: g.continuation(s, 0, **k)
        t[s + ':next_request'] = lambda s=s, **k: g.next_request(s, [0], **k)
        t[s + ':claims_cut'] = lambda s=s, **k: g.claims(s, at_most_bp=3, strict=False, **k)
        t[s + ':fasta_none'] = lambda s=s, **k: g.to_fasta(s, leaves=[], **k)
        if g.labels:
            t[s + ':routes'] = lambda s=s, **k: g.routes({'id': 0}, s, spell=True, **k)
            t[s + ':label_walks'] = lambda s=s, **k: g.label_walks({'id': 0}, s, **k)
            t[s + ':label_walks_name'] = lambda s=s, **k: g.label_walks(
                g.labels[0].name, s, **k)
            t[s + ':walks_label'] = lambda s=s, **k: g.walks(s, labels=[{'id': 0}], **k)
    return t


def _charge_points(fn):
    """The work used before each charge of an unlimited call: a budget of exactly that
    stops the call at that charge."""
    pts = []
    orig = B.LocalBudget.charge

    def rec(self, work, nbytes=0):
        pts.append(self.used_work)
        return orig(self, work, nbytes)
    B.LocalBudget.charge = rec
    try:
        try:
            fn(budget=LocalBudget())
        except Exception:
            pass
    finally:
        B.LocalBudget.charge = orig
    return sorted(set(pts))


def _charged(b):
    """The derivation names a budget's current call charged (its once-per-call keys)."""
    return {k[1] for k in b._used if isinstance(k, tuple) and len(k) >= 2}


def _uncharged(made, b):
    return {(w, n) for w, n in made if (n in PRICED or n in G_PRICED)
            and n not in _charged(b)}


class TestAdmissionBeforeDerivedCaches(unittest.TestCase):
    """Finding 2."""

    def test_the_reviewers_probe(self):
        text = comb_annotate(5000)
        calls = [('spell', lambda g, b: g.spell('right', 0, budget=b)),
                 ('support_profile', lambda g, b: g.support_profile('right', 0, budget=b)),
                 ('support_changes', lambda g, b: g.support_changes('right', 0, budget=b)),
                 ('walks_at', lambda g, b: ops.walks_at(g, 'right', [0], budget=b)),
                 ('rank_walks', lambda g, b: ops.rank_walks(g, 'right', budget=b)),
                 ('walks', lambda g, b: g.walks('right', top=1, budget=b)),
                 ('to_fasta', lambda g, b: g.to_fasta(budget=b))]
        for name, fn in calls:
            with self.subTest(op=name):
                g = parse(text)
                a = g.arms['right']
                before = _caches(g)
                b = LocalBudget(work_units=0, memory_mb=0)
                gc.collect()
                tracemalloc.start()
                try:
                    with self.assertRaises(LocalBudgetExceeded) as e:
                        fn(g, b)
                    peak = tracemalloc.get_traced_memory()[1]
                finally:
                    tracemalloc.stop()
                self.assertEqual(set(), _caches(g) - before)
                self.assertNotIn('paths', a.cache)
                self.assertEqual(0, sum(1 for o in gc.get_objects() if type(o) is Path
                                        and getattr(o, 'arm', None) is a))
                self.assertEqual(0, b.usage()['work_units'])
                self.assertEqual(0, b.usage()['memory_bytes'])
                self.assertIn(e.exception.stop.phase, ('setup', 'rank'))
                # a few KB of frames, not the 0.7 MB of the arm's paths
                self.assertLess(peak, 100_000)

    def test_a_refused_call_builds_no_derived_cache(self):
        shapes = [('comb', lambda: parse(comb_annotate(60)))]
        shapes += [(x[0], lambda x=x: _fresh(x)) for x in ITEMS[::3]]
        bad = []
        for name, make in shapes:
            for op in _ops(make()):
                g = make()
                fn = _ops(g)[op]
                before = _caches(g)
                try:
                    fn(budget=LocalBudget(work_units=0, memory_mb=0))
                except Exception:
                    pass
                made = _caches(g) - before
                if made:
                    bad.append((name, op, sorted(made)))
        self.assertEqual([], bad)

    def test_a_stop_at_any_charge_leaves_only_charged_derivations(self):
        shapes = [('comb', lambda: parse(comb_annotate(30)))]
        shapes += [(x[0], lambda x=x: _fresh(x)) for x in ITEMS[::6]]
        bad = []
        n = 0
        for name, make in shapes:
            for op in _ops(make()):
                pts = _charge_points(_ops(make())[op])
                for w in pts[::max(1, len(pts) // 10)]:
                    g = make()
                    fn = _ops(g)[op]
                    before = _caches(g)
                    b = LocalBudget(work_units=w)
                    try:
                        fn(budget=b)
                    except LocalBudgetExceeded:
                        pass
                    except Exception:
                        continue
                    n += 1
                    for x in _uncharged(_caches(g) - before, b):
                        bad.append((name, op, w, x))
        self.assertGreater(n, 2000)
        self.assertEqual([], bad)

    def test_a_stopped_comparison_leaves_only_charged_derivations(self):
        bad = []
        n = 0
        for (ka, ma), (kb, mb) in _pairs(both_orders=False):
            for mode in ('claims', 'walks', 'labels', 'prefix_subset'):
                pts = _charge_points(lambda budget: ma().compare(mb(), mode=mode,
                                                                 budget=budget))
                for w in [0] + pts[::max(1, len(pts) // 6)]:
                    a, b = ma(), mb()
                    ca, cb_ = _caches(a), _caches(b)
                    bud = LocalBudget(work_units=w)
                    a.compare(b, mode=mode, budget=bud)
                    n += 1
                    for g, before in ((a, ca), (b, cb_)):
                        for x in _uncharged(_caches(g) - before, bud):
                            bad.append((ka, kb, mode, w, x))
        self.assertGreater(n, 1000)
        self.assertEqual([], bad)

    def test_a_refused_tool_builds_no_derived_cache(self):
        # the tools resolve a walk argument (graphlet_walk, graphlet_sequence, an export
        # of walks) or a view's label names before the operation charges: admitted first
        zero = ToolLimits(view=LocalLimits(work_units=0), heavy=LocalLimits(work_units=0),
                          parse=LocalLimits(work_units=0))
        bad = []
        n = 0
        with tempfile.TemporaryDirectory() as d:
            for x in ITEMS[1::5]:
                names = [c[0] for c in _tool_calls(None, None, _fresh(x), None)]
                for i, name in enumerate(names):
                    store = GraphletStore(os.path.join(d, 'sp%d' % n))
                    tools = GraphletTools(store, secret=b'k' * 32, export_dir=d,
                                          local_limits=zero)
                    g = _fresh(x)
                    h = store.put_parsed(g, x[1]['graphlet'], {'seeds': []})
                    view = None
                    if name == 'view_claims':
                        # (the view's own construction reads the paths: made first)
                        view = GraphletTools(store, secret=b'k' * 32).graphlet_subtrie(
                            h, [{'ref': g.labels[0].ref}])['handle']
                    before = _caches(g)
                    name, call = _tool_calls(tools, h, g, view)[i]
                    out = call()
                    n += 1
                    made = {c for c in _caches(g) - before if c[1] in PRICED | G_PRICED}
                    if made:
                        bad.append((x[0], name, out.get('error'), sorted(made)))
        self.assertGreater(n, 150)
        self.assertEqual([], bad)

    def test_cold_and_warm_charge_the_same(self):
        # admission moved, the charges did not: a warm call (paths cached) charges what a
        # cold one does
        for name, fn in (('spell', lambda g, b: g.spell('right', 3, budget=b)),
                         ('changes', lambda g, b: g.support_changes('right', 3, budget=b)),
                         ('walks_at', lambda g, b: ops.walks_at(g, 'right', [1, 3],
                                                                with_claims=True,
                                                                budget=b))):
            g = parse(comb_annotate(50))
            cold, warm = LocalBudget(), LocalBudget()
            out = fn(g, cold)
            self.assertEqual(out, fn(g, warm), name)
            self.assertEqual(cold.usage(), warm.usage(), name)
            self.assertEqual(out, fn(parse(comb_annotate(50)), None), name)


def _tool_calls(t, h, g, view):
    out = []
    for s in g.arms:
        out += [(s + ':walks', lambda s=s: t.graphlet_walks(h, s)),
                (s + ':walks_id', lambda s=s: t.graphlet_walks(h, s, rank='id')),
                (s + ':walk', lambda s=s: t.graphlet_walk(h, s, 0)),
                (s + ':support', lambda s=s: t.graphlet_support(h, s, 0)),
                (s + ':splits', lambda s=s: t.graphlet_splits(h, s)),
                (s + ':claims', lambda s=s: t.graphlet_claims(h, s)),
                (s + ':sequence', lambda s=s: t.graphlet_sequence(h, s, walk=0)),
                (s + ':continue', lambda s=s: t.traverse_continue(h, s, 0, execute=False)),
                (s + ':export_walk', lambda s=s: t.graphlet_export(
                    h, format='json', what='walk', arm=s, walks=[0], path='w.json')),
                (s + ':export_none', lambda s=s: t.graphlet_export(
                    h, format='fasta', arm=s, walks=[], path='n.fa'))]
    out += [('summary', lambda: t.graphlet_summary(h)), ('labels', lambda: t.graphlet_labels(h))]
    if g.labels:
        out += [('labels_name', lambda: t.graphlet_labels(h, name={'ref': g.labels[0].ref})),
                ('subtrie', lambda: t.graphlet_subtrie(h, [{'ref': g.labels[0].ref}])),
                ('view_claims', lambda: t.graphlet_claims(view, list(g.arms)[0],
                                                          labels=[g.labels[0].name]))]
    return out


def _pairs(both_orders=True):
    """[((key, make a), (key, make b))]: the retrievals of one seed, pairwise."""
    by = collections.defaultdict(list)
    for x in ITEMS:
        by[_fresh(x).seed.sequence].append(x)
    out = []
    for xs in by.values():
        for i, x in enumerate(xs):
            for j, y in enumerate(xs):
                if i != j and (both_orders or i < j):
                    out.append(((x[0], lambda x=x: _fresh(x)), (y[0], lambda y=y: _fresh(y))))
    return out


def _comb_pairs():
    out = []
    for n, m in ((0, 0), (0, 1), (1, 100), (0, 100), (10, 1000), (100, 100), (100, 10),
                 (50, 7), (3, 3)):
        out.append((('comb%d' % n, lambda n=n: parse(comb_annotate(n))),
                    ('comb%d' % m, lambda m=m: parse(comb_annotate(m)))))
    return out


class TestCompareCostIsALowerBound(unittest.TestCase):
    """Finding 4: at_least <= what the completed comparison charges, every mode; at most
    the estimate, under the library's work model, and all it charges where it is exact."""

    def test_the_reviewers_probe(self):
        a, b = parse(comb_annotate(1)), parse(comb_annotate(100))
        cost = ops.compare_cost(a, b, mode='claims')
        bud = LocalBudget()
        parse(comb_annotate(1)).compare(parse(comb_annotate(100)), mode='claims', budget=bud)
        self.assertLessEqual(cost['work_units']['at_least'], bud.used_work)
        self.assertLessEqual(cost['memory_bytes']['at_least'], bud.peak_bytes)
        # the restricted DAG's ends: comb 1's two walks at the depth, not comb 100's 101
        self.assertLess(cost['work_units']['at_least'], 9197)

    def test_at_least_is_below_the_charge_everywhere(self):
        bad = []
        n = 0
        for (ka, ma), (kb, mb) in _pairs(both_orders=False) + _comb_pairs():
            a0 = ma()
            kws = [{}] + ([{'labels': [{'ref': a0.labels[-1].ref}]}] if a0.labels else [])
            for mode in ('claims', 'walks', 'labels', 'prefix_subset'):
                for kw in kws:
                    cost = ops.compare_cost(ma(), mb(), mode=mode, **kw)
                    bud = LocalBudget()
                    ma().compare(mb(), mode=mode, budget=bud, **kw)
                    n += 1
                    w = cost['work_units']
                    self.assertLessEqual(w['at_least'], w['estimate'], (ka, kb, mode))
                    self.assertEqual(B.WORK_MODEL, cost['work_model'])
                    if cost['exact'] and w['at_least'] != bud.used_work:
                        # nothing keyed (no bases in a mode keyed by them): all it charges
                        bad.append((ka, kb, mode, kw, 'exact', w['at_least'], bud.used_work))
                    if w['at_least'] > bud.used_work \
                            or cost['memory_bytes']['at_least'] > bud.peak_bytes:
                        bad.append((ka, kb, mode, kw, w['at_least'],
                                    bud.used_work, cost['memory_bytes']['at_least'],
                                    bud.peak_bytes))
        self.assertGreater(n, 500)
        self.assertEqual([], bad)

    def test_incomparable_is_exact(self):
        # nothing is keyed: the estimate is all the comparison charges
        a, b = T.graphlet('fork'), T.graphlet('merge')
        est = ops.compare_cost(a, b)
        bud = LocalBudget()
        a.compare(b, budget=bud)
        self.assertTrue(est['exact'])
        self.assertEqual(bud.used_work, est['work_units']['at_least'])


class TestLoadKeepsItsReceipt(unittest.TestCase):
    """Finding 5."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = self.tmp.name
        self.path = os.path.join(self.d, 'saved.mgt')
        _fresh(_real('fork')).save(self.path)

    def tearDown(self):
        self.tmp.cleanup()

    def tools(self, spool='sp', **kw):
        store = GraphletStore(os.path.join(self.d, spool), **kw)
        return store, GraphletTools(store, secret=b'k' * 32, export_dir=self.d,
                                    local_limits=ToolLimits())

    def test_the_reviewers_probe(self):
        store, tools = self.tools(max_ram_mb=0, parse_limits=LocalLimits(work_units=1))
        out = tools.graphlet_load(self.path)
        self.assertNotIn('error', out, out)
        self.assertEqual([out['handle']], [e['handle'] for e in store.list()])
        self.assertIsNotNone(out['summary'])
        self.assertTrue(out['local']['complete'])
        # the answer of a store that keeps the model resident, handle aside
        _, resident = self.tools('resident')
        same = resident.graphlet_load(self.path)
        self.assertEqual(dict(same, handle=None), dict(out, handle=None))

    def test_a_stored_entry_is_always_answered_with_its_handle(self):
        # every budget from nothing to enough: a load that stops before the entry is
        # stored answers local_budget_exceeded and stores nothing; once the entry is
        # stored, the answer carries its handle -- without the summary when the summary
        # stopped (local states the stop)
        _, probe = self.tools('probe', max_ram_mb=0)
        full = probe.graphlet_load(self.path)['local']['usage']['work_units']
        seen = set()
        for w in sorted(set(range(0, full + 2, max(1, full // 60))) | {full - 1, full}):
            store, tools = self.tools('sp%d' % w, max_ram_mb=0)
            out = tools.graphlet_load(self.path, budget={'work_units': w})
            stored = [e['handle'] for e in store.list()]
            if 'error' in out:
                self.assertEqual('local_budget_exceeded', out['error'], out)
                self.assertEqual([], stored, w)
                seen.add('refused')
                continue
            self.assertEqual([out['handle']], stored, w)
            if out['summary'] is None:
                self.assertFalse(out['local']['complete'], out)
                self.assertEqual('work', out['local']['stop']['resource'])
                seen.add('no summary')
            else:
                seen.add('complete')
        self.assertEqual({'refused', 'no summary', 'complete'}, seen)

    def test_a_stop_on_the_stores_parse_limits_offers_no_raise(self):
        store = GraphletStore(os.path.join(self.d, 'sp'),
                              parse_limits=LocalLimits(work_units=1))
        tools = GraphletTools(store, secret=b'k' * 32, export_dir=self.d)
        out = tools.graphlet_load(self.path)
        self.assertEqual('local_budget_exceeded', out['error'])
        self.assertNotIn('raise_local_budget', out['stop']['actions'])
        self.assertIn('store\'s parse limits', out['message'])
        self.assertEqual([], store.list())
        # a stored body parsed on demand: the same classification without local limits
        st2 = GraphletStore(os.path.join(self.d, 'sp2'), max_ram_mb=0)
        h = st2.load(self.path)
        st2.parse_limits = LocalLimits(work_units=1)
        t2 = GraphletTools(st2, secret=b'k' * 32, export_dir=self.d)
        out = t2.graphlet_summary(h)
        self.assertEqual('local_budget_exceeded', out['error'])
        self.assertNotIn('raise_local_budget', out['stop']['actions'])

    def test_the_store_returns_the_model_it_parsed(self):
        store = GraphletStore(os.path.join(self.d, 'sp'), max_ram_mb=0)
        h, g = store.load_graphlet(self.path)
        self.assertFalse(store.in_ram(h))
        self.assertEqual(MT.dump(g, envelope=True), MT.dump(store.graphlet(h), envelope=True))
        self.assertEqual(h, store.list()[0]['handle'])


def _peak(fn):
    gc.collect()
    tracemalloc.start()
    base = tracemalloc.get_traced_memory()[0]
    try:
        fn()
        return tracemalloc.get_traced_memory()[1] - base
    finally:
        tracemalloc.stop()


class TestCompareAccountBoundsThePeak(unittest.TestCase):
    """The batch's memory gaps: the account of compare() >= its traced peak, on every
    pair of retrievals of one seed and on the combs, every mode (the cross pair it failed
    on, in both orders: test_the_reported_cross_pair)."""

    def test_every_pair(self):
        low = []
        n = 0
        for (ka, ma), (kb, mb) in _pairs(both_orders=False) + _comb_pairs():
            for mode in ('claims', 'walks', 'labels', 'prefix_subset'):
                b = LocalBudget()
                ma().compare(mb(), mode=mode, budget=b)
                x, y = ma(), mb()
                peak = _peak(lambda: x.compare(y, mode=mode))
                n += 1
                if b.peak_bytes < peak:
                    low.append((ka, kb, mode, b.peak_bytes, peak))
        self.assertGreater(n, 200)
        self.assertEqual([], low)

    def test_the_reported_cross_pair(self):
        # 1,354,074 charged against a traced 1,645,185 (0.82x) at bf8dcee4
        a, b = _real('sra_hairpin__left_only_follow'), _real('sra_hairpin__left_only')
        for mode in ('prefix_subset', 'walks'):
            bud = LocalBudget()
            _fresh(a).compare(_fresh(b), mode=mode, budget=bud)
            x, y = _fresh(a), _fresh(b)
            self.assertGreaterEqual(bud.peak_bytes, _peak(lambda: x.compare(y, mode=mode)),
                                    mode)

    def test_the_name_index_is_charged_at_its_size(self):
        for n in (1, 10, 977, 5000):
            g = parse(comb_annotate(1))
            g.labels = [Label(i, 'header', 7, 1000 + i, 'name%d' % (i % max(1, n // 2)))
                        for i in range(n)]
            g.cache.clear()
            b = LocalBudget()
            with b.scope('probe'):
                ops._uses_g(b, g, 'label_index')
            peak = _peak(lambda: ops._label_index(g))
            g.cache.clear()
            # with the call's base (its frames: a few hundred bytes here), and by the
            # per-label price alone once the labels dominate
            self.assertGreaterEqual(b.peak_bytes, peak, n)
            if n >= 977:
                self.assertGreaterEqual(b.peak_bytes - B.CALL_BASE, peak, n)

    def test_the_index_is_the_same(self):
        # the index built without a list per key: the same keys, order and ids
        for x in ITEMS[::4]:
            g = _fresh(x)
            by_ref, by_name = {}, {}
            for l in g.labels:
                by_ref.setdefault(l.ref, []).append(l.id)
                by_name.setdefault(l.name, []).append(l.id)
            want = ({k: tuple(v) for k, v in by_ref.items()},
                    {k: tuple(v) for k, v in by_name.items()})
            got = ops._label_index(g)
            self.assertEqual(want, got)
            self.assertEqual([list(d) for d in want], [list(d) for d in got])


# =========================================================================================
# The review of these fixes (the working tree on bf8dcee4): findings 1-3 on the library.

class _Clock:
    def __init__(self, t):
        self.t = t

    def __call__(self):
        return self.t


def _j_line_bytes(text):
    second = text.split('\n', 2)[1]
    return len(second[2:].encode('utf-8')) if second.startswith('J ') else 0


class TestTheJLineChargeIsDeterministic(unittest.TestCase):
    """Review finding 1: standalone_text() sized its J-line charge by the entry file,
    which holds the created and accessed times -- their printed length follows the clock,
    so the same export of the same entry was charged differently from run to run (the
    reviewer's det_probe.py: 52,493 against 52,541 bytes, and a memory limit of 52,493
    stopped one run while the other answered)."""

    CLOCKS = (1791204329.5, 1791204329.1234567)

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = self.tmp.name
        self.item = _real('fork')

    def tearDown(self):
        self.tmp.cleanup()

    def put(self, spool, clock, **kw):
        store = GraphletStore(os.path.join(self.d, spool), clock=clock, **kw)
        h = store.put_parsed(_fresh(self.item), self.item[1]['graphlet'], {'seeds': []})
        tools = GraphletTools(store, secret=b'k' * 32, export_dir=self.d,
                              local_limits=ToolLimits())
        return store, h, tools

    @staticmethod
    def usages(store, h, tools):
        b = LocalBudget()
        store.standalone_text(h, budget=b)
        lib = (b.usage()['work_units'], b.usage()['memory_bytes'])
        exp = tools.graphlet_export(h, format='mgt', path='e.mgt', budget={})
        sav = tools.graphlet_save(h, 's.mgt', budget={})
        return lib, exp['local']['usage'], sav['local']['usage']

    def test_two_clocks_charge_the_same(self):
        got = []
        sizes = set()
        for i, t0 in enumerate(self.CLOCKS):
            store, h, tools = self.put('sp%d' % i, _Clock(t0))
            sizes.add(os.path.getsize(store._entry_path(h)))
            got.append(self.usages(store, h, tools))
        self.assertEqual(2, len(sizes))          # the entry files differ in length
        self.assertEqual(got[0], got[1])
        # a memory limit of exactly that usage answers on both clocks
        m = got[0][1]['memory_bytes']
        for i, t0 in enumerate(self.CLOCKS):
            store, h, tools = self.put('lim%d' % i, _Clock(t0))
            for out in (tools.graphlet_export(h, format='mgt', path='l.mgt',
                                              budget={'memory_mb': m / float(1 << 20)}),
                        tools.graphlet_save(h, 'l2.mgt',
                                            budget={'memory_mb': m / float(1 << 20)})):
                self.assertNotIn('error', out, (t0, out))

    def test_a_rewritten_entry_and_a_restart_charge_the_same(self):
        clock = _Clock(1791204329.5)
        store, h, tools = self.put('sp', clock, ttl_disk_s=1000)
        sizes, got = set(), set()
        # each use lags by more than a tenth of the TTL: get() rewrites the entry file
        for t in (1791204329.5, 1791204429.6234567, 1791204529.75):
            clock.t = t
            got.add(repr(self.usages(store, h, tools)))
            sizes.add(os.path.getsize(store._entry_path(h)))
        self.assertGreater(len(sizes), 1)
        reopened = GraphletStore(os.path.join(self.d, 'sp'), clock=clock, ttl_disk_s=1000)
        t2 = GraphletTools(reopened, secret=b'k' * 32, export_dir=self.d,
                           local_limits=ToolLimits())
        got.add(repr(self.usages(reopened, h, t2)))
        self.assertEqual(1, len(got), got)
        # a completed parse rewrites the entry without its parsed flag: the same bound
        e = store.get(h)
        bound = e.j_bound
        e.parsed = False
        store._write_entry(e)
        self.assertEqual(bound, e.j_bound)

    def test_the_bound_holds_the_j_line(self):
        # every retrieval and a view of it, and an envelope with text beyond ASCII
        low = []
        store = GraphletStore(os.path.join(self.d, 'sp'))
        tools = GraphletTools(store, secret=b'k' * 32, export_dir=self.d,
                              local_limits=ToolLimits())
        for x in ITEMS:
            g = _fresh(x)
            hs = [store.put_parsed(g, x[1]['graphlet'], {'seeds': []})]
            if g.labels:
                hs.append(tools.graphlet_subtrie(hs[0], [{'id': g.labels[0].id}])['handle'])
            for h in hs:
                n = _j_line_bytes(store.standalone_text(h))
                if n > store.get(h).j_bound:
                    low.append((x[0], h, n, store.get(h).j_bound))
        g = _fresh(self.item)
        g.envelope = dict(g.envelope, note='\u00e9t\u00e9 \u2014 \u6f22\u5b57 \U0001f600')
        h = store.put_parsed(g, self.item[1]['graphlet'], {'seeds': []})
        n = _j_line_bytes(store.standalone_text(h))
        self.assertIn('\u6f22', store.standalone_text(h))
        self.assertLessEqual(n, store.get(h).j_bound)
        self.assertEqual([], low)
        # an entry without an envelope, a view or derived_from writes no J line: no charge
        body = MT.dump(T.graphlet('linear'), envelope=False)
        bare = store.put_parsed(parse(body), body, None)
        self.assertEqual(0, _j_line_bytes(store.standalone_text(bare)))
        self.assertEqual(0, store.get(bare).j_bound)


class TestAnEmptySelectionResolvesNoWalk(unittest.TestCase):
    """Review finding 2: to_fasta(leaves=[]) built the arm's paths (a Path per walk)
    without charging them, and a refused zero-budget call (stopped at its join) left them
    behind -- 20,001 Paths, usage 0/0 (the reviewer's fasta_empty.py and
    fasta_empty_tool.py)."""

    def setUp(self):
        # the first call imports the export module: its code objects are no allocation
        # of the call measured below
        parse(comb_annotate(2)).to_fasta('right', leaves=[])

    def test_the_library_call(self):
        g = parse(comb_annotate(20000))
        a = g.arms['right']
        before = _caches(g)
        b = LocalBudget(work_units=0, memory_mb=0)
        gc.collect()
        tracemalloc.start()
        try:
            with self.assertRaises(LocalBudgetExceeded):
                g.to_fasta('right', leaves=[], budget=b)
            peak = tracemalloc.get_traced_memory()[1]
        finally:
            tracemalloc.stop()
        self.assertEqual(set(), _caches(g) - before)
        self.assertNotIn('paths', a.cache)
        self.assertEqual(0, b.usage()['work_units'])
        self.assertEqual(0, b.usage()['memory_bytes'])
        self.assertLess(peak, 100_000)
        # completed: no text and no paths; the charge (the call's base and the join, as
        # before) does not follow the arm's size
        b = LocalBudget()
        self.assertEqual('', g.to_fasta('right', leaves=[], budget=b))
        small = LocalBudget()
        parse(comb_annotate(3)).to_fasta('right', leaves=[], budget=small)
        self.assertEqual(small.usage(), b.usage())
        self.assertEqual('', g.to_fasta(leaves=[]))
        self.assertNotIn('paths', a.cache)

    def test_the_tool(self):
        with tempfile.TemporaryDirectory() as d:
            store = GraphletStore(os.path.join(d, 'sp'))
            r, resp = _comb_response(20000)
            g = MT.from_response(r, resp)
            h = store.put_parsed(g, r['graphlet'], {'seeds': []})
            tools = GraphletTools(store, secret=b'k' * 32, export_dir=d,
                                  local_limits=ToolLimits())
            a = g.arms['right']
            gc.collect()
            tracemalloc.start()
            try:
                out = tools.graphlet_export(h, format='fasta', walks=[], path='x.fa',
                                            budget={'work_units': 0, 'memory_mb': 0})
                peak = tracemalloc.get_traced_memory()[1]
            finally:
                tracemalloc.stop()
            self.assertEqual('local_budget_exceeded', out['error'])
            self.assertEqual({'work_units': 0, 'memory_bytes': 0}, out['local']['usage'])
            self.assertNotIn('paths', a.cache)
            self.assertLess(peak, 100_000)


def _comb_response(n):
    import test_traverse_memory as TM
    return TM._comb(n)


class TestALoadStoppedByTheCallersBudget(unittest.TestCase):
    """Review finding 3: graphlet_load (tools without local limits) stated every stop as
    one on the store's parse limits -- also the stop of the caller's ambient budget on
    the dump that stores the loaded graphlet: raise_local_budget dropped and a store's
    limits named that need not exist (the reviewer's amb_load.py)."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = self.tmp.name
        self.path = os.path.join(self.d, 'f.mgt')
        _fresh(_real('fork')).save(self.path)

    def tearDown(self):
        self.tmp.cleanup()

    def test_the_reviewers_probe(self):
        for w in (0, 5, 50, 200):
            store = GraphletStore(os.path.join(self.d, 'sp%d' % w))
            tools = GraphletTools(store, secret=b'k' * 32, export_dir=self.d)
            with MT.local_budget(LocalLimits(work_units=w)):
                out = tools.graphlet_load(self.path)
            self.assertEqual('local_budget_exceeded', out['error'], w)
            self.assertIn('raise_local_budget', out['stop']['actions'])
            self.assertNotIn('parse limits', out['message'])
            self.assertEqual([], store.list())
        store = GraphletStore(os.path.join(self.d, 'enough'))
        tools = GraphletTools(store, secret=b'k' * 32, export_dir=self.d)
        with MT.local_budget(LocalLimits(work_units=10 ** 6)):
            out = tools.graphlet_load(self.path)
        self.assertEqual([out['handle']], [e['handle'] for e in store.list()])

    def test_the_stores_parse_limits_are_still_stated(self):
        # the same load under the store's parse limits, ambient budget or none; and a save
        # without local limits, which parses the stored body on demand (it offered
        # raise_local_budget for the store's limits)
        for ambient in (None, LocalLimits(work_units=10 ** 6)):
            store = GraphletStore(os.path.join(self.d, 'pl%s' % bool(ambient)),
                                  parse_limits=LocalLimits(work_units=1))
            tools = GraphletTools(store, secret=b'k' * 32, export_dir=self.d)
            if ambient is None:
                out = tools.graphlet_load(self.path)
            else:
                with MT.local_budget(ambient):
                    out = tools.graphlet_load(self.path)
            self.assertEqual('local_budget_exceeded', out['error'])
            self.assertNotIn('raise_local_budget', out['stop']['actions'])
            self.assertIn("store's parse limits", out['message'])
        store = GraphletStore(os.path.join(self.d, 'save'), max_ram_mb=0)
        h = store.load(self.path)
        store.parse_limits = LocalLimits(work_units=1)
        tools = GraphletTools(store, secret=b'k' * 32, export_dir=self.d)
        out = tools.graphlet_save(h, 's.mgt')
        self.assertEqual('local_budget_exceeded', out['error'])
        self.assertNotIn('raise_local_budget', out['stop']['actions'])
        self.assertIn("store's parse limits", out['message'])
        # the caller's own budget on that save (no parse limits): its lever is offered
        store.parse_limits = None
        with MT.local_budget(LocalLimits(work_units=0)):
            out = tools.graphlet_save(h, 's.mgt')
        self.assertEqual('local_budget_exceeded', out['error'])
        self.assertIn('raise_local_budget', out['stop']['actions'])
        self.assertNotIn('parse limits', out['message'])


if __name__ == '__main__':
    unittest.main()
