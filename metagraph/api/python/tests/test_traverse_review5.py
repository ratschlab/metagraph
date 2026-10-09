"""The external review of pass 5 (HEAD 5801aea1), the library's part, and the audit of the
efficiency batch's library speed-ups that the review asked for.

  * finding 7: a continuation's change_cost override is checked whole -- entry shape,
    string endpoints, finite costs >= 0, the default, a constant's value, the model's
    type -- before any endpoint enters a set or any cost the arithmetic: an entry such as
    [["C"], "A", 0.5] raised TypeError from a set membership test in _switch_reach(),
    which escaped next_request() and traverse_continue(execute=False) instead of their
    field-qualified ValueError / bounded bad_argument (the reviewer's
    /tmp/review_pass5_library_probes.py, turned into the tests below, plus a fuzz over
    malformed shapes); and, from the review of those fixes, an entries field that is not a
    list under any model (a budgeted call charged len(entries) whatever the model: TypeError
    from next_request(budget=...) and a ToolLimits tool), and an infinite or NaN
    branching.max_label_branches override (int(inf) raised OverflowError);
  * the speed-ups (P8): ambiguity and iteration order kept, the canonical arrays kept for
    serialisation, the indexes the stage L memory model counts, no cache of every prefix
    spelling;
  * the small-call regression (P10): an index built on a graphlet's first single-label
    call made routes()/label_walks() and small exports on 1,000-label graphlets several
    times slower than the scan it replaced; the indexes are now built only once they pay
    off, with outputs identical either way.
"""

import collections
import copy
import glob
import math
import os
import random
import sys
import tempfile
import unittest
from array import array

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402
import traverse_golden as TG  # noqa: E402
from test_traverse_mcp import FakeClient  # noqa: E402
from test_traverse_speedups import _routes_reference  # noqa: E402

import metagraph.traverse as MT  # noqa: E402
from metagraph.traverse import GraphletStore, LocalBudget, derive, ops, parser  # noqa: E402
from metagraph.traverse.mcp_tools import (MIN_MAX_BYTES, GraphletTools, ToolLimits,  # noqa: E402
                                         _size)


def _reference_reach(cost, sources, targets, budget, sinks=()):
    """An independent reachability (Bellman-Ford over every name, the reviewer's) for a
    valid cost: {target: cheapest chain loss} within |budget|; a sink ends a chain."""
    names = list(dict.fromkeys(list(sources) + list(targets)))
    model = cost.get('model', 'forbid')
    if model == 'forbid' or not sources:
        return {}
    if model == 'constant':
        c = float(cost['value'])
        return {t: c for t in targets if t not in sources and c <= budget}
    table = {}
    for u, v, c in cost.get('entries') or ():
        table[(u, v)] = float(c)
    d = cost.get('default', 'forbid')
    default = math.inf if d == 'forbid' else float(d)
    dist = {x: 0.0 for x in sources}
    for _ in range(len(names)):
        changed = False
        for u in names:
            if u not in dist or u in sinks:
                continue
            for v in names:
                if u == v:
                    continue
                x = dist[u] + table.get((u, v), default)
                if x <= budget and x < dist.get(v, math.inf):
                    dist[v] = x
                    changed = True
        if not changed:
            break
    return {t: dist[t] for t in targets if t in dist and t not in sources}


# malformed entries: each is refused, whatever else the table holds
BAD_ENTRIES = [
    [['C'], 'A', 0.5],                 # the review's: a list as an endpoint
    ['C', ['A'], 0.5],
    [['C', 'A'], 'B', 1],
    [1, 'A', 0.5],                     # numbers as endpoints
    ['C', 2, 0.5],
    [None, 'A', 0.5],
    [{'name': 'C'}, 'A', 0.5],
    [True, 'A', 0.5],
    ['C', 'A', float('nan')],          # costs that are not finite numbers >= 0
    ['C', 'A', float('inf')],
    ['C', 'A', -float('inf')],
    ['C', 'A', -0.5],
    ['C', 'A', -1],
    ['C', 'A', True],
    ['C', 'A', '0.5'],
    ['C', 'A', None],
    ['C', 'A', [0.5]],
    ['C', 'A', 10 ** 400],             # beyond a double
    [],                                # wrong arity
    ['C'],
    ['C', 'A'],
    ['C', 'A', 0.5, 1],
    'C->A',                            # not an entry at all
    5,
    None,
    {'from': 'C', 'to': 'A', 'cost': 0.5},
]

BAD_COSTS = [
    5, 'table', ['table'],                                       # not an object
    {'model': ['table']}, {'model': 7}, {'model': None},         # the model's type
    {'model': 'constant', 'value': float('nan')},
    {'model': 'constant', 'value': float('inf')},
    {'model': 'constant', 'value': -1},
    {'model': 'constant', 'value': [1]},
    {'model': 'constant', 'value': True},
    {'model': 'table', 'default': float('nan'), 'entries': []},
    {'model': 'table', 'default': float('inf'), 'entries': []},
    {'model': 'table', 'default': -0.1, 'entries': []},
    {'model': 'table', 'default': ['forbid'], 'entries': []},
    {'model': 'table', 'default': 'never', 'entries': []},
    {'model': 'table', 'entries': {'C': 'A'}},
    {'model': 'table', 'entries': 'C,A,1'},
    {'model': 'table', 'entries': 3},
    {'model': 'table', 'default': 1},                # no entries (the server needs them)
    {'model': 'table', 'default': 1, 'entries': None},
] + [{'model': 'table', 'default': 'forbid', 'entries': [['A', 'B', 1], e]}
     for e in BAD_ENTRIES]

# entries that are not a list, under every model and none: a budgeted call charges their
# count whatever the model, and len(5) raised TypeError there (the review of the pass-5
# fixes: repro_entries_len.py, and its budgeted fuzz's example {'entries': -0.5})
NON_LIST_ENTRIES = [5, True, 1.5, -0.5, 'C,A,1', {'C': 'A'}]
BAD_COSTS += [dict(base, entries=e)
              for base in ({'model': 'constant', 'value': 1}, {'model': 'forbid'}, {},
                           {'model': 'mystery'})
              for e in NON_LIST_ENTRIES]

# branching.max_label_branches overrides that are no count: int(inf) raised OverflowError,
# which escaped the tools' bad_argument, and int(nan) a ValueError naming no field (the
# review of the pass-5 fixes: repro_mlb.py; Python's json reads Infinity and 1e999)
BAD_MAX_LABEL_BRANCHES = [float('inf'), -float('inf'), float('nan'), 1e999, -1, 1.5, -0.5,
                          True, None, 'x', '2', [2]]


class TestFinding7MalformedContinuationCosts(unittest.TestCase):
    """The reviewer's probe: entries [[["C"], "A", 0.5]] raised TypeError from
    next_request() and from traverse_continue(execute=False)."""

    def setUp(self):
        self.g = T.graphlet('switch_chain')
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.store = GraphletStore(self.tmp.name)
        self.h = self.store.put_graphlet(self.g, source='mini')
        self.tools = GraphletTools(self.store, {'mini': FakeClient()})

    def overrides(self, cost):
        return {'labels': {'change_cost': copy.deepcopy(cost)}}

    def test_the_reviewers_probe(self):
        invalid = self.overrides({'model': 'table', 'default': 'forbid',
                                  'entries': [[['C'], 'A', 0.5]]})
        with self.assertRaises(ValueError) as e:
            self.g.next_request('right', [1], **copy.deepcopy(invalid))
        self.assertIn('labels.change_cost.entries[0]', str(e.exception))
        out = self.tools.traverse_continue(self.h, 'right', 1, overrides=invalid,
                                           execute=False)
        self.assertEqual('bad_argument', out.get('error'), out)
        self.assertIn('labels.change_cost.entries[0]', out['message'])

    def test_every_malformed_cost_is_a_field_qualified_value_error(self):
        for cost in BAD_COSTS:
            with self.subTest(cost=cost):
                with self.assertRaises(ValueError) as e:
                    ops._switch_reach(cost, ['A'], ['B', 'C'], 2)
                self.assertTrue(str(e.exception).startswith('labels.change_cost'),
                                str(e.exception))
                with self.assertRaises(ValueError) as e:
                    self.g.next_request('right', [1], **self.overrides(cost))
                self.assertTrue(str(e.exception).startswith('labels.change_cost')
                                or str(e.exception).startswith('strategy.labels.change_cost'),
                                str(e.exception))

    def test_a_constant_without_its_value(self):
        # (an override {'model': 'constant'} merges with the retrieval's value: valid)
        with self.assertRaises(ValueError) as e:
            ops._switch_reach({'model': 'constant'}, ['A'], ['B'], 2)
        self.assertIn('labels.change_cost.value', str(e.exception))

    def test_the_index_of_the_bad_entry_is_named(self):
        for i in (0, 1, 3):
            entries = [['A', 'B', 1]] * 4
            entries = entries[:i] + [['A', 7, 1]] + entries[i + 1:]
            with self.assertRaises(ValueError) as e:
                ops._switch_reach({'model': 'table', 'entries': entries}, ['A'], ['B'], 2)
            self.assertIn('entries[%d]' % i, str(e.exception))

    def test_the_tools_answer_bad_argument_within_the_ceiling(self):
        for cost in BAD_COSTS:
            for execute in (False, True):
                out = self.tools.traverse_continue(self.h, 'right', 1,
                                                   overrides=self.overrides(cost),
                                                   execute=execute, max_bytes=MIN_MAX_BYTES)
                with self.subTest(cost=cost, execute=execute):
                    self.assertEqual('bad_argument', out.get('error'), out)
                    self.assertLessEqual(_size(out), MIN_MAX_BYTES)

    def test_a_malformed_cost_is_refused_even_when_the_caller_gives_extra(self):
        # no switch search reads the cost then, but the request would carry it
        ov = self.overrides({'model': 'table', 'entries': [['A', 'B', -1]]})
        ov['labels']['extra'] = ['A']
        with self.assertRaises(ValueError):
            self.g.next_request('right', [1], **ov)

    def test_a_loss_budget_that_is_not_finite_and_non_negative_is_refused(self):
        for b in (float('nan'), float('inf'), -1, -0.5, 10 ** 400):
            with self.subTest(loss_budget=b):
                with self.assertRaises(ValueError) as e:
                    self.g.next_request('right', [1], labels={'loss_budget': b})
                self.assertIn('labels.loss_budget', str(e.exception))
                out = self.tools.traverse_continue(self.h, 'right', 1,
                                                   overrides={'labels': {'loss_budget': b}},
                                                   execute=False)
                self.assertEqual('bad_argument', out.get('error'), out)

    def test_valid_costs_still_build_the_request(self):
        for cost in ({'model': 'table', 'default': 'forbid', 'entries': [['C', 'A', 0.5]]},
                     {'model': 'table', 'default': 0, 'entries': []},
                     {'model': 'table', 'default': 1, 'entries': (('C', 'A', 0),)},
                     {'model': 'constant', 'value': 0},
                     {'model': 'forbid'}):
            with self.subTest(cost=cost):
                req = self.g.next_request('right', [1], **self.overrides(cost))
                self.assertEqual([], T.request_violations(req))
                out = self.tools.traverse_continue(self.h, 'right', 1,
                                                   overrides=self.overrides(cost),
                                                   execute=False)
                self.assertNotIn('error', out)

    def test_next_requests_and_a_budgeted_request_refuse_it_too(self):
        bad = self.overrides({'model': 'table', 'entries': [[['C'], 'A', 0.5]]})
        with self.assertRaises(ValueError) as e:
            self.g.next_requests('right', [1], **copy.deepcopy(bad))
        self.assertIn('labels.change_cost.entries[0]', str(e.exception))
        with self.assertRaises(ValueError) as e:
            self.g.next_request('right', [1], budget=LocalBudget(), **copy.deepcopy(bad))
        self.assertIn('labels.change_cost.entries[0]', str(e.exception))

    def test_tuples_price_as_lists(self):
        lists = {'model': 'table', 'default': 'forbid',
                 'entries': [['A', 'B', 1], ['B', 'C', 0.5], ['A', 'B', 0.25]]}
        tuples = dict(lists, entries=tuple(tuple(e) for e in lists['entries']))
        self.assertEqual(ops._switch_reach(lists, ['A'], ['B', 'C'], 2),
                         ops._switch_reach(tuples, ['A'], ['B', 'C'], 2))
        self.assertEqual({'B': 0.25, 'C': 0.75}, ops._switch_reach(tuples, ['A'], ['B', 'C'], 2))


class TestFinding7Fuzz(unittest.TestCase):
    """Random tables: one with a malformed entry anywhere is a ValueError (never another
    exception), whatever its valid entries; every valid one prices as the independent
    reachability (the reviewer's 5,000 valid tables)."""

    def test_malformed_entries_anywhere(self):
        rng = random.Random(7105)
        names = ['A', 'B', 'C', 'D', 'E']
        for cell in range(2000):
            entries = [[rng.choice(names), rng.choice(names), rng.choice([0, 0.5, 1, 2])]
                       for _ in range(rng.randrange(0, 6))]
            bad = copy.deepcopy(rng.choice(BAD_ENTRIES))
            entries.insert(rng.randrange(len(entries) + 1), bad)
            cost = {'model': 'table', 'default': rng.choice(['forbid', 0, 1, 3]),
                    'entries': entries}
            sources = rng.sample(names, rng.randrange(1, 3))
            targets = [x for x in names if x not in sources]
            try:
                ops._switch_reach(cost, sources, targets, rng.choice([0, 1, 2, 5]),
                                  sinks=set(rng.sample(targets, rng.randrange(len(targets)))))
            except ValueError as e:
                self.assertIn('labels.change_cost.entries[%d]' % entries.index(bad), str(e))
            else:
                self.fail('cell %d: %r was accepted' % (cell, cost))

    def test_valid_tables_match_an_independent_reachability(self):
        rng = random.Random(963)
        vals = [0, 0.1, 0.3, 1, 2, 5]
        for cell in range(5000):
            n = rng.randrange(2, 13)
            names = [str(i) for i in range(n)]
            sources = rng.sample(names, rng.randrange(1, n))
            targets = [x for x in names if x not in sources]
            default = rng.choice(['forbid'] + vals)
            entries = [[u, v, rng.choice(vals)] for u in names for v in names
                       if rng.random() < .4]
            cost = {'model': 'table', 'default': default, 'entries': entries}
            budget = rng.choice(vals)
            sinks = set(rng.sample(targets, rng.randrange(len(targets) + 1)))
            want = _reference_reach(cost, sources, targets, budget, sinks)
            got = ops._switch_reach(cost, sources, targets, budget, sinks=sinks)
            self.assertEqual(want, got, (cell, cost, sources, targets, sinks, budget))


class TestFinding7UnderLocalBudgets(unittest.TestCase):
    """The review of the pass-5 fixes: under a local budget (next_request(budget=...), a
    GraphletTools with local_limits) a malformed cost escaped as TypeError -- the switch
    search charged len(entries) of a model that is not a table, whose entries were never
    checked. Every malformed cost is the same field-qualified ValueError / bounded
    bad_argument with a budget as without one."""

    def setUp(self):
        self.g = T.graphlet('switch_chain')
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.store = GraphletStore(self.tmp.name)
        self.h = self.store.put_graphlet(self.g, source='mini')
        self.tools = GraphletTools(self.store, {'mini': FakeClient()},
                                   local_limits=ToolLimits())

    @staticmethod
    def overrides(cost):
        return {'labels': {'change_cost': copy.deepcopy(cost)}}

    def assertFieldError(self, e, field='labels.change_cost'):
        self.assertTrue(str(e).startswith(field) or str(e).startswith('strategy.' + field),
                        str(e))

    def test_the_reviewers_repro(self):
        for cost in ({'model': 'constant', 'value': 1, 'entries': 5},
                     {'model': 'forbid', 'entries': True},
                     {'model': 'constant', 'value': 1, 'entries': 1.5},
                     {'entries': -0.5}):
            with self.subTest(cost=cost):
                for budget in (None, LocalBudget()):
                    kw = {} if budget is None else {'budget': budget}
                    with self.assertRaises(ValueError) as e:
                        self.g.next_request('right', [1], **kw, **self.overrides(cost))
                    self.assertIn('labels.change_cost.entries', str(e.exception))
                for execute in (False, True):
                    out = self.tools.traverse_continue(self.h, 'right', 1,
                                                       overrides=self.overrides(cost),
                                                       execute=execute)
                    self.assertEqual('bad_argument', out.get('error'), out)
                    self.assertIn('labels.change_cost.entries', out['message'])

    def test_every_malformed_cost_under_a_local_budget(self):
        for cost in BAD_COSTS:
            with self.subTest(cost=cost):
                with self.assertRaises(ValueError) as e:
                    self.g.next_request('right', [1], budget=LocalBudget(),
                                        **self.overrides(cost))
                self.assertFieldError(e.exception)
                with self.assertRaises(ValueError) as e:
                    self.g.next_requests('right', [1], budget=LocalBudget(),
                                         **self.overrides(cost))
                self.assertFieldError(e.exception)
                if isinstance(cost, dict):
                    with self.assertRaises(ValueError) as e:
                        ops._switch_reach(cost, ['A'], ['B', 'C'], 2, lb=LocalBudget())
                    self.assertFieldError(e.exception)
                for execute in (False, True):
                    out = self.tools.traverse_continue(self.h, 'right', 1,
                                                       overrides=self.overrides(cost),
                                                       execute=execute, max_bytes=MIN_MAX_BYTES)
                    self.assertEqual('bad_argument', out.get('error'), (execute, out))
                    self.assertLessEqual(_size(out), MIN_MAX_BYTES)

    def test_a_list_or_no_entries_is_still_charged_by_its_count(self):
        # entries are charged by their count before the check reads them (2 lwu per entry
        # in work model 1; 2 W_ELEM, the check and the build, in work model 2). Under a model
        # that does not read them they are refused, as the server's Strict parse refuses
        # them (L8, the review of 2026-10-06): a constant override deep-merged onto a table's
        # strategy kept that table's list, which passed here and was refused there; an
        # override that names another model now replaces the cost whole (TestL8 in
        # test_traverse_review_p3.py)
        def usage(cost):
            b = LocalBudget()
            try:
                got = ops._switch_reach(cost, ['A'], ['B', 'C'], 2, lb=b)
            except ValueError as e:
                got = str(e)
            return got, b.usage()['work_units']
        base = {'model': 'constant', 'value': 1}
        got0, w0 = usage(base)
        self.assertEqual({'B': 1.0, 'C': 1.0}, got0)
        entries = [['A', 'B', 0.5], ['B', 'C', 0.5], ['X', 'Y', 1]]
        per = 2 * ops.W_ELEM
        for given, k in ((None, 0), ([], 0), (entries, 3), (tuple(map(tuple, entries)), 3)):
            got, w = usage(dict(base, entries=given))
            self.assertIn('labels.change_cost.entries is not a field of model', got)
            self.assertEqual(w0 + per * k, w)
        for model in ({'model': 'forbid'}, {}):
            got, w = usage(dict(model, entries=entries))
            self.assertIn('labels.change_cost.entries is not a field of model', got)
        with self.assertRaises(ValueError) as e:
            self.g.next_request('right', [1], budget=LocalBudget(),
                                **self.overrides(dict(base, entries=entries)))
        self.assertIn('labels.change_cost.entries', str(e.exception))

    def test_max_label_branches_that_is_no_count(self):
        plain = GraphletTools(self.store, {'mini': FakeClient()})
        for v in BAD_MAX_LABEL_BRANCHES:
            ov = {'branching': {'max_label_branches': v}}
            with self.subTest(max_label_branches=v):
                for budget in (None, LocalBudget()):
                    kw = {} if budget is None else {'budget': budget}
                    with self.assertRaises(ValueError) as e:
                        self.g.next_request('right', [1], **kw, **copy.deepcopy(ov))
                    self.assertIn('branching.max_label_branches', str(e.exception))
                    with self.assertRaises(ValueError) as e:
                        self.g.next_requests('right', [1], **kw, **copy.deepcopy(ov))
                    self.assertIn('branching.max_label_branches', str(e.exception))
                for tools in (plain, self.tools):
                    for execute in (False, True):
                        out = tools.traverse_continue(self.h, 'right', 1,
                                                      overrides=copy.deepcopy(ov),
                                                      execute=execute, max_bytes=MIN_MAX_BYTES)
                        self.assertEqual('bad_argument', out.get('error'), out)
                        # (fitted to the ceiling: the field's name leads the message)
                        self.assertTrue(out['message'].startswith('branching.max_label_'),
                                        out['message'])
                        self.assertLessEqual(_size(out), MIN_MAX_BYTES)
        # the counts it takes are unchanged: integers >= 0 (an integral float too), and
        # "unlimited"
        for v in (0, 1, 2, 2.0, 10 ** 30, 'unlimited'):
            with self.subTest(max_label_branches=v):
                req = self.g.next_request('right', [1],
                                          branching={'max_label_branches': v})
                self.assertEqual(v, req['strategy']['branching']['max_label_branches'])


class TestFinding7BudgetedFuzz(unittest.TestCase):
    """The reviewer's budgeted fuzz (fuzz_p7_budget.py: 65 TypeError escapes in 10,000
    cases, all from len(entries)), here at a fixed seed: random overrides -- change_cost,
    loss_budget, extra, max_label_branches -- through next_request / next_requests under a
    LocalBudget and through a ToolLimits tool with execute False and True. Only a
    ValueError, or a bad_argument result, may come of a malformed one."""

    ATOMS = [None, True, False, 0, 1, -1, 0.5, -0.5, float('nan'), float('inf'),
             -float('inf'), 10 ** 400, '', 'A', 'forbid', 'table', 'constant', [], {}, ['A'],
             {'a': 1}, [['A', 'B', 1]], ('A', 'B', 1), 2 ** 64, 1e308, 'nan']

    def rand_val(self, rng, d=0):
        r = rng.random()
        if d > 2 or r < 0.6:
            return rng.choice(self.ATOMS)
        if r < 0.8:
            return [self.rand_val(rng, d + 1) for _ in range(rng.randrange(4))]
        return {rng.choice(['model', 'value', 'default', 'entries', 'x']):
                self.rand_val(rng, d + 1) for _ in range(rng.randrange(3))}

    def rand_cost(self, rng, names):
        c = {}
        m = rng.choice(['table', 'constant', 'forbid', None, self.rand_val(rng)])
        if m is not None:
            c['model'] = m
        if rng.random() < 0.5:
            c['value'] = rng.choice(self.ATOMS)
        if rng.random() < 0.5:
            c['default'] = rng.choice(self.ATOMS)
        if rng.random() < 0.8:
            if rng.random() < 0.8:
                c['entries'] = [
                    [rng.choice(names + self.ATOMS[:12]), rng.choice(names + self.ATOMS[:12]),
                     rng.choice(self.ATOMS[:14])] if rng.random() < 0.5
                    else self.rand_val(rng) for _ in range(rng.randrange(4))]
            else:
                c['entries'] = self.rand_val(rng)
        return c if rng.random() < 0.95 else self.rand_val(rng)

    def rand_overrides(self, rng, names):
        # labels stays an object: the library leaves an annotate retrieval's labels section
        # to the server, and the fake backend of execute=True reads it as one
        lab = {}
        if rng.random() < 0.9:
            lab['change_cost'] = self.rand_cost(rng, names)
        if rng.random() < 0.3:
            lab['loss_budget'] = rng.choice(self.ATOMS)
        if rng.random() < 0.2:
            lab['extra'] = rng.choice([names[:2], self.rand_val(rng)])
        o = {'labels': lab}
        if rng.random() < 0.3:
            o['branching'] = {'max_label_branches': rng.choice(
                [1, 2, 'unlimited', -1, 'x', float('inf'), float('nan'), 2.0, 2.5])}
        return o

    def test_no_exception_but_value_error_escapes(self):
        rng = random.Random(20261004)
        graphs = [(n, T.graphlet(n)) for n in ('switch_chain', 'merge', 'reentry', 'quorum',
                                                'same_name', 'fork', 'linear')]
        for f in ('merge_constrain_r100', 'budgetchain_r80', 'wide_forbid_r40'):
            (_, r, resp), = TG.items_of([os.path.join(T.DATA, 'review3', f + '.graphlet.json')])
            graphs.append((f, MT.from_response(r, resp)))
        with tempfile.TemporaryDirectory() as d:
            store = GraphletStore(d)
            tools = GraphletTools(store, {'mini': FakeClient()}, local_limits=ToolLimits())
            handles = {n: store.put_graphlet(g, source='mini') for n, g in graphs}
            kinds = collections.Counter()
            for case in range(1500):
                name, g = rng.choice(graphs)
                names = [l.name for l in g.labels] or ['A']
                side = rng.choice(sorted(g.arms))
                npaths = len(derive.paths(g.arms[side]))
                w = rng.randrange(npaths)
                ov = self.rand_overrides(rng, names)
                which = rng.choice(['next_request', 'next_requests', 'tool', 'tool_execute'])
                try:
                    if which == 'next_request':
                        g.next_request(side, [w], budget=LocalBudget(), **copy.deepcopy(ov))
                    elif which == 'next_requests':
                        g.next_requests(side, list(range(min(npaths, 3))),
                                        budget=LocalBudget(), **copy.deepcopy(ov))
                    else:
                        out = tools.traverse_continue(handles[name], side, w,
                                                      overrides=copy.deepcopy(ov),
                                                      execute=which == 'tool_execute')
                        kinds[out.get('error')] += 1
                except ValueError:
                    kinds['ValueError'] += 1
                except Exception as e:      # noqa: BLE001 - what the test is about
                    self.fail('case %d (%s on %s, %s %d, overrides %r): %s: %s'
                              % (case, which, name, side, w, ov, type(e).__name__, e))
            # the fuzz reached both outcomes, refusals and requests
            self.assertGreater(kinds['ValueError'] + kinds['bad_argument'], 300)
            self.assertGreater(kinds[None], 20)


# ------------------------------------------------------------------ P8 / P10

def _items():
    """(key, result, response) of every committed retrieval with a graphlet: fresh models
    are parsed from these (the derived caches start empty)."""
    files = sorted(glob.glob(os.path.join(T.DOCS, '**', '*.graphlet.json'), recursive=True))
    return TG.items_of(files + TG.golden_files())


ITEMS = _items()


def _fresh(item):
    return MT.from_response(item[1], item[2])


def _with_merges(min_labels=1):
    """Fresh constrain models with merges, the most labels first."""
    out = []
    for it in ITEMS:
        g = _fresh(it)
        if g.mode == 'constrain' and len(g.labels) >= min_labels and any(
                len(s.parents) > 1 for a in g.arms.values() for s in a.segments):
            out.append((it, g))
    out.sort(key=lambda x: -len(x[1].labels))
    return out


LAZY_KEYS = ('merge_scanned', 'merge_parts', 'partition_shared', 'run_scans')


class TestP8MergeSetsAndRunIndex(unittest.TestCase):
    """The indexes answer as the arrays and the scans they replace, in the same order,
    at every stage of their lazy construction (P10), and keep one copy."""

    def test_merge_parts_answers_as_the_arrays_at_every_stage(self):
        n = 0
        for it, g in _with_merges():
            labels = list(range(len(g.labels))) + [len(g.labels), 10 ** 6]
            for a in g.arms.values():
                merges = [s for s in a.segments if len(s.parents) > 1]
                for m in merges:
                    want = [[l in p for p in m.partition] for l in labels]
                    for stage in range(derive._MERGE_SCANS + 2):
                        parts = derive.merge_parts(a, m.id)
                        self.assertEqual(len(m.parents), len(parts))
                        self.assertEqual(want, [[l in p for p in parts] for l in labels],
                                         (it[0], stage))
                        if stage < derive._MERGE_SCANS:
                            self.assertIs(parts, m.partition)       # scanned: the arrays
                        else:
                            self.assertTrue(all(isinstance(p, frozenset) for p in parts))
                    n += 1
                built = dict(a.cache.get('merge_parts') or {})
                ps = derive.partition_sets(a)
                for k, v in built.items():
                    self.assertIs(v, ps[k])     # the lazily built sets are reused: one copy
                for k in LAZY_KEYS[:3]:
                    self.assertNotIn(k, a.cache)
                for m in merges:
                    self.assertIs(ps[m.id], derive.merge_parts(a, m.id))
        self.assertGreater(n, 10)

    def test_label_runs_scans_then_indexes_in_run_order(self):
        for it in ITEMS:
            g = _fresh(it)
            for a in g.arms.values():
                for l in list(range(len(g.labels))) + [len(g.labels)]:
                    self.assertEqual([r.id for r in a.runs if r.label == l],
                                     list(derive.label_runs(a, l)), it[0])
                if len(g.labels) > derive._RUN_SCANS:
                    self.assertIn('runs_by_label', a.cache)
                    self.assertNotIn('run_scans', a.cache)

    def test_the_store_sees_the_lazy_state_grow(self):
        # the store re-measures a model only when cache_signature() changes: a step that
        # adds memory must change it (the lazy state was a fixed-size tuple at first)
        it, g = _with_merges()[0]
        a = next(x for x in g.arms.values() if any(len(s.parents) > 1 for s in x.segments))
        m = next(s for s in a.segments if len(s.parents) > 1)
        sigs = [ops.cache_signature(g)]
        for _ in range(derive._MERGE_SCANS + 1):
            derive.merge_parts(a, m.id)
            sigs.append(ops.cache_signature(g))
        self.assertNotEqual(sigs[0], sigs[1])                  # the first scan's entry
        self.assertNotEqual(sigs[-2], sigs[-1])                # the sets built
        with tempfile.TemporaryDirectory() as d:
            store = GraphletStore(d)
            h = store.put_graphlet(_fresh(it), source='mini')
            g2 = store.graphlet(h)
            for side in g2.arms:
                for l in range(len(g2.labels)):
                    g2.routes({'id': l}, side)
            store.refresh(h)
            self.assertEqual(store._meter[h][0] + ops.cache_bytes(g2), store._ram[h][1])
            self.assertTrue(any(x.cache.get('merge_parts') or x.cache.get('partition_sets')
                                for x in g2.arms.values()))

    def test_the_store_charges_the_run_index_on_a_single_label_arm(self):
        # the review of the pass-5 fixes (sig_probe.py, store_probe.py): on an arm with one
        # label, runs_by_label (one entry) replaced label_runs()'s count (an int) and left
        # every length equal; cache_signature() did not change, and the store never charged
        # the index (504 B unaccounted on documents/ambiguous_split). It reads the keys too
        g = _fresh(ITEMS[0])
        a = next(iter(g.arms.values()))
        a.cache['x'] = 1
        before = ops.cache_signature(g)
        del a.cache['x']
        a.cache['y'] = 1
        self.assertNotEqual(before, ops.cache_signature(g))      # a swap of equal lengths
        n = 0
        for it in ITEMS:
            g = _fresh(it)
            sides = [side for side, x in g.arms.items() if len({r.label for r in x.runs}) == 1]
            if not sides:
                continue
            with tempfile.TemporaryDirectory() as d:
                store = GraphletStore(d)
                h = store.put_graphlet(g, source='mini')
                g2 = store.graphlet(h)
                for side in sides:
                    a = g2.arms[side]
                    label = a.runs[0].label
                    g2.summary()
                    g2.walks(side)
                    g2.claims(side, strict=False)
                    sigs = []
                    for _ in range(derive._RUN_SCANS + 1):
                        sigs.append(ops.cache_signature(g2))
                        g2.routes({'id': label}, side)
                        store.refresh(h)
                        self.assertEqual(store._meter[h][0] + ops.cache_bytes(g2),
                                         store._ram[h][1], (it[0], side))
                    self.assertIn('runs_by_label', a.cache)
                    self.assertNotIn('run_scans', a.cache)
                    self.assertNotEqual(sigs[-1], ops.cache_signature(g2), (it[0], side))
                    n += 1
        self.assertGreaterEqual(n, 3)       # the committed retrievals' single-label arms

    def test_a_budgeted_call_charges_the_same_whatever_the_lazy_state(self):
        # charges are cold prices (L1): a model whose unbudgeted calls left scans counted
        # and some merges' sets built charges a budgeted call exactly as a fresh one, and
        # the call builds the whole priced indexes (the memory the model counts)
        n = 0
        for it, _ in _with_merges(2)[:6]:
            for side in ('left', 'right'):
                fresh = _fresh(it)
                if side not in fresh.arms:
                    continue
                b1 = LocalBudget()
                r1 = fresh.routes({'id': 0}, side, spell=True, budget=b1)
                warm = _fresh(it)
                for l in range(min(len(warm.labels), 4)):
                    warm.routes({'id': l}, side, spell=True)
                    warm.label_walks({'id': l})
                b2 = LocalBudget()
                r2 = warm.routes({'id': 0}, side, spell=True, budget=b2)
                self.assertEqual(r1, r2)
                self.assertEqual(b1.usage(), b2.usage(), (it[0], side))
                for g in (fresh, warm):
                    a = g.arms[side]
                    self.assertIn('partition_sets', a.cache)
                    self.assertIn('runs_by_label', a.cache)
                    for k in LAZY_KEYS:
                        self.assertNotIn(k, a.cache)
                n += 1
        self.assertGreater(n, 3)

    def test_the_indexes_leave_dump_and_the_public_arrays_unchanged(self):
        for it in ITEMS[:40]:
            g = _fresh(it)
            before = parser.dump(g)
            for side in g.arms:
                for l in range(min(len(g.labels), 8)):
                    g.routes({'id': l}, side, spell=True)
                    g.label_walks({'id': l})
            g.claims(strict=False)
            for a in g.arms.values():
                derive.partition_sets(a)
                for s in a.segments:
                    for p in s.partition:
                        self.assertIsInstance(p, array)
            self.assertEqual(before, parser.dump(g), it[0])


def _strings(x, seen, out):
    if id(x) in seen:
        return
    seen.add(id(x))
    if isinstance(x, str):
        out.append(x)
    elif isinstance(x, dict):
        for k, v in x.items():
            _strings(k, seen, out)
            _strings(v, seen, out)
    elif isinstance(x, (list, tuple, set, frozenset)):
        for v in x:
            _strings(v, seen, out)


class TestP8NoSpellingCached(unittest.TestCase):
    def test_exports_and_walks_cache_no_spelling(self):
        # walk_batch() shares prefixes within one call; no spelling of a walk (bases of
        # two segments or more) outlives it in the model's caches
        n = 0
        for it in ITEMS:
            g = _fresh(it)
            if not any(a.segments and a.segments[0].walk for a in g.arms.values()):
                continue
            g.to_fasta()
            g.to_gfa()
            g.to_json()
            for side in g.arms:
                g.walks(side, top=5)
            longest = max(len(s.walk or '') for a in g.arms.values() for s in a.segments)
            got = []
            _strings([g.cache] + [a.cache for a in g.arms.values()], set(), got)
            spelled = [x for x in got if len(x) > longest and set(x) <= set('ACGTN')]
            self.assertEqual([], spelled, it[0])
            n += 1
        self.assertGreater(n, 10)


class TestP10WalkStrategies(unittest.TestCase):
    """walk_batch() spells target by target below its bound, as a batch above it and one
    target directly: the same chains and bases every way (P10)."""

    def _all(self, a, targets):
        return {t: (derive.chain(a, t),
                    derive.walk_bases(a, t) if a.segments[0].walk is not None else None)
                for t in targets}

    def test_every_strategy_gives_the_same_chains_and_bases(self):
        rng = random.Random(1010)
        saved = derive._EACH_STEPS
        n = 0
        try:
            for it in ITEMS:
                g = _fresh(it)
                for a in g.arms.values():
                    if not a.segments:
                        continue
                    ids = [s.id for s in a.segments]
                    leaves = derive.leaves(a)
                    sets = [leaves, leaves[:1], leaves[:5], rng.sample(ids, min(3, len(ids))), ids]
                    for steps in (0, 1, saved, 10 ** 9):
                        derive._EACH_STEPS = steps
                        for tg in sets:
                            want = self._all(a, tg)
                            for form in (list(tg), tuple(tg), set(tg), iter(tg)):
                                got = derive.walk_batch(a, form)
                                self.assertEqual(want, got, (it[0], steps))
                                self.assertEqual(sorted(set(tg)), list(got))
                                self.assertEqual(len(got), len({id(c) for c, _ in got.values()}))
                            # the bases alone, for the callers that read only them, and
                            # the chains alone, a list of its own per target
                            self.assertEqual({t: (None, w) for t, (_, w) in want.items()},
                                             derive.walk_batch(a, tg, chains=False),
                                             (it[0], steps))
                            chains = derive.walk_batch(a, tg, spell=False)
                            self.assertEqual({t: (c, None) for t, (c, _) in want.items()},
                                             chains, (it[0], steps))
                            self.assertEqual(len(chains),
                                             len({id(c) for c, _ in chains.values()}))
                            n += 1
        finally:
            derive._EACH_STEPS = saved
        self.assertGreater(n, 200)

    def test_bases_on_some_segments_only_is_a_value_error_every_way(self):
        saved = derive._EACH_STEPS
        try:
            for steps in (0, saved, 10 ** 9):
                derive._EACH_STEPS = steps
                g = T.graphlet('merge')
                a = g.arms['right']
                a.segments[3].walk = None
                below = [x for x in derive.leaves(a) if 3 in derive.chain(a, x)]
                self.assertTrue(below)
                for tg in (below[:1], below, derive.leaves(a)):
                    for chains in (True, False):
                        with self.assertRaises(ValueError):
                            derive.walk_batch(a, tg, chains=chains)
        finally:
            derive._EACH_STEPS = saved

    def test_a_depth_the_model_does_not_keep_changes_no_answer(self):
        g = T.graphlet('merge')
        for a in g.arms.values():
            leaves = derive.leaves(a)
            want = self._all(a, leaves)
            for s in a.segments:
                s.depth = 0                 # a hand-built model's default
            self.assertEqual(want, derive.walk_batch(a, leaves))
            self.assertEqual(want, derive.walk_batch(a, list(reversed(leaves))))
            self.assertEqual({t: (None, w) for t, (_, w) in want.items()},
                             derive.walk_batch(a, leaves, chains=False))


class TestP10OneLabelBuildsNoIndex(unittest.TestCase):
    def test_one_label_scans_and_a_listing_builds_the_indexes(self):
        it, g = _with_merges(100)[0]
        for side in g.arms:
            self.assertEqual(_routes_reference(g, {'id': 0}, side, True),
                             g.routes({'id': 0}, side, spell=True))
        g.label_walks({'id': 0})
        for a in g.arms.values():
            self.assertNotIn('runs_by_label', a.cache)
            self.assertNotIn('partition_sets', a.cache)
            self.assertFalse(a.cache.get('merge_parts'))   # no merge's sets built
        for side in g.arms:
            for l in range(len(g.labels)):
                self.assertEqual(_routes_reference(g, {'id': l}, side, True),
                                 g.routes({'id': l}, side, spell=True))
        self.assertTrue(all('runs_by_label' in a.cache for a in g.arms.values()))
        self.assertTrue(any(a.cache.get('merge_parts') for a in g.arms.values()))


if __name__ == '__main__':
    unittest.main()
