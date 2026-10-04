"""Stage L (DESIGN §21, owner decisions L1-L9): work and allocation budgets for the
library's local operations, the store's parses and the MCP tools.

  T-L1/L20  off by default: no budget, no change -- and an unlimited budget changes no
            answer either (the golden digests of every operation)
  T-L2      cold price: the same units cold, warm, and after other calls
  T-L3      determinism: the same stop in this process and under other hash seeds
  T-L4      the account and the work never exceed their limits; needed_at_least is above
  T-L5      every operation stops at one unit (no uncharged path), and completes at its
            usage
  T-L7      caches after a stop: the unbudgeted answer afterwards is a fresh model's
  T-L8      partial lists: whole rows, a prefix, and resuming concatenates to the list
  T-L9      compare: comparable 'unknown', equal None, no lists, local_stop
  T-L10     exports and saves: no text, no file
  T-L11     parse: admission, a local failure that keeps the text, adversarial bodies
  T-L12     the account is at least the traced peak (a model, checked on the fixtures)
  T-L13     deadline and cancel: non-deterministic, stated so
  T-L14     the ambient budget
  T-L16     the store's parse limits: unparsed entries and body-level copies
  T-L17-21  the tools: local blocks, local_budget_exceeded within max_bytes, partial pages
            and resume cursors, clamps and hooks, receipts
  T-L22     the real offline retrievals under the default ToolLimits
"""

import gc
import glob
import json
import os
import subprocess
import sys
import tempfile
import threading
import time
import tracemalloc
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402
import traverse_golden as TG  # noqa: E402
from test_traverse_scale import comb_annotate  # noqa: E402

import metagraph.traverse as MT  # noqa: E402
from metagraph.traverse import (  # noqa: E402
    GraphletFormatError, GraphletStore, LocalBudget, LocalBudgetExceeded, LocalLimits,
    derive, local_budget, ops, parse,
)
from metagraph.traverse import budget as B  # noqa: E402
from metagraph.traverse.mcp_tools import (  # noqa: E402
    MIN_MAX_BYTES, TOOL_CLASS, GraphletTools, ToolLimits, _size,
)

HERE = os.path.dirname(os.path.abspath(__file__))
API = os.path.dirname(HERE)


def _items():
    files = sorted(glob.glob(os.path.join(T.DOCS, '**', '*.graphlet.json'), recursive=True))
    return TG.items_of(files + TG.golden_files())


ITEMS = _items()
REAL = [x for x in ITEMS if '/real/' in x[0] or x[0].startswith('real') or 'sra_' in x[0]
        or 'uhgg_' in x[0] or 'mini_' in x[0]]
# the heaviest checks run on every other retrieval: every operation still runs on both
# modes, both arms and every kind of fixture (CLI documents, review3, real cells)
SAMPLE = ITEMS[::2]


def fresh(item):
    return MT.from_response(item[1], item[2])


def op_table(g):
    """name -> fn(budget=..., **kw) for every budgeted library operation on |g|."""
    t = {}
    sides = list(g.arms)
    t['to_json'] = lambda **k: g.to_json(**k)
    t['to_fasta'] = lambda **k: g.to_fasta(**k)
    t['to_gfa'] = lambda **k: g.to_gfa(**k)
    t['dump'] = lambda **k: g.dump(**k)
    t['summary'] = lambda **k: g.summary(**k)
    t['label_summary'] = lambda **k: g.label_summary(**k)
    t['claims'] = lambda **k: g.claims(strict=False, **k)
    if g.labels:
        t['subgraph'] = lambda **k: g.subgraph([{'id': 0}], **k)
    for s in sides:
        t[s + ':walks'] = lambda s=s, **k: g.walks(s, **k)
        t[s + ':walks_id'] = lambda s=s, **k: g.walks(s, by='id', **k)
        t[s + ':rank'] = lambda s=s, **k: ops.rank_walks(g, s, **k)
        t[s + ':splits'] = lambda s=s, **k: g.splits(s, **k)
        t[s + ':profile'] = lambda s=s, **k: g.support_profile(s, 0, **k)
        t[s + ':changes'] = lambda s=s, **k: g.support_changes(s, 0, **k)
        t[s + ':spell'] = lambda s=s, **k: g.spell(s, 0, **k)
        if g.labels:
            t[s + ':routes'] = lambda s=s, **k: g.routes({'id': 0}, s, spell=True, **k)
            t[s + ':label_walks'] = lambda s=s, **k: g.label_walks({'id': 0}, s, **k)
    return t


def run(fn, b):
    try:
        return fn(budget=b), None
    except LocalBudgetExceeded as e:
        return None, e


def applicable(g, name):
    """Whether the operation answers on |g| unbudgeted (not every one applies)."""
    try:
        op_table(g)[name]()
        return True
    except Exception:
        return False


# ======================================================================= the library

class TestOffByDefault(unittest.TestCase):
    """T-L1: no budget anywhere is the old library (the golden digests in
    test_traverse_speedups.py); an unlimited budget changes no answer either."""

    def test_an_unlimited_budget_changes_no_answer(self):
        for key, r, resp in ITEMS:
            plain = TG.op_digests(MT, MT.from_response(r, resp))
            budgeted = TG.op_digests(MT, MT.from_response(r, resp),
                                     lambda: {'budget': LocalBudget()})
            bad = [k for k in plain if plain[k] != budgeted.get(k)]
            self.assertEqual([], bad, key)

    def test_the_golden_digests_hold_under_an_unlimited_ambient_budget(self):
        want = TG.read_json(os.path.join(TG.GOLDEN_DIR, 'committed.json.gz'))
        items = TG.items_of(TG.golden_files(), base=TG.DATA)
        with local_budget() as b:
            got = TG.all_digests(MT, items, tools=False)
        self.assertGreater(b.used_work, 0)
        bad = sorted(k for k in got if want.get(k) != got[k])
        self.assertEqual([], bad[:10])


class TestColdPrice(unittest.TestCase):
    """T-L2 / L1: a derivation is charged at its structural price whether a cache holds it
    or not, so the same call charges the same units cold, warm, and after other calls
    built other caches."""

    def test_the_same_units_cold_warm_and_after_other_calls(self):
        for item in SAMPLE:
            g = fresh(item)
            names = [n for n in op_table(g) if applicable(fresh(item), n)]
            cold = {}
            for n in names:
                b = LocalBudget()
                op_table(fresh(item))[n](budget=b)
                cold[n] = (b.used_work, b.peak_bytes)
                self.assertGreater(b.used_work, 0, (item[0], n))
            warm = fresh(item)
            for n in names:                      # every cache built by every other call
                op_table(warm)[n]()
            for n in reversed(names):
                b = LocalBudget()
                op_table(warm)[n](budget=b)
                self.assertEqual(cold[n], (b.used_work, b.peak_bytes), (item[0], n))


class TestStops(unittest.TestCase):
    """T-L3 to T-L5 on every operation of every committed retrieval: one unit stops it, its
    usage completes it, usage - 1 stops it with the numbers of a stop, half of it stops it
    the same way twice; the memory account likewise."""

    def test_every_operation(self):
        n = 0
        for item in SAMPLE:
            for name in op_table(fresh(item)):
                if not applicable(fresh(item), name):
                    continue
                n += 1
                b = LocalBudget()
                op_table(fresh(item))[name](budget=b)
                use, peak = b.used_work, b.peak_bytes
                fn = op_table(fresh(item))[name]
                self.assertIsNone(run(fn, LocalBudget(work_units=use))[1], (item[0], name))
                _, e = run(fn, LocalBudget(work_units=use - 1))
                self.assertIsNotNone(e, (item[0], name))
                st = e.stop
                self.assertEqual(('work', True, 'lwu', B.WORK_MODEL),
                                 (st.resource, st.deterministic, st.unit, st.work_model))
                self.assertTrue(st.used <= st.limit < st.needed_at_least, st.as_dict())
                self.assertIn('raise_local_budget', st.actions)
                self.assertEqual('process_locally', st.actions[-1])
                if use > 1:
                    _, e1 = run(fn, LocalBudget(work_units=1))
                    self.assertIsNotNone(e1, (item[0], name))
                h = use // 2
                _, ea = run(fn, LocalBudget(work_units=h))
                _, eb = run(op_table(fresh(item))[name], LocalBudget(work_units=h))
                if h:
                    self.assertEqual(ea.stop.as_dict(), eb.stop.as_dict())
                self.assertIsNone(run(fn, LocalBudget(memory_mb=peak / (1 << 20)))[1])
                _, em = run(fn, LocalBudget(memory_mb=(peak - 1) / (1 << 20)))
                self.assertEqual('memory', em.stop.resource)
                self.assertTrue(em.stop.used <= em.stop.limit < em.stop.needed_at_least)
        self.assertGreater(n, 200)

    def test_a_stop_is_no_bad_argument_and_no_incomplete_recording(self):
        self.assertFalse(issubclass(LocalBudgetExceeded, ValueError))
        self.assertFalse(issubclass(LocalBudgetExceeded, RuntimeError))
        self.assertTrue(issubclass(LocalBudgetExceeded, Exception))

    def test_the_same_stop_under_other_hash_seeds(self):
        # T-L3: a set of strings iterates in another order under another seed; no charge
        # may depend on it
        script = r'''
import json, sys
sys.path.insert(0, %r); sys.path.insert(0, %r)
import traverse_testlib as T
from test_traverse_stage_l import op_table, run, fresh, ITEMS
from metagraph.traverse import LocalBudget
out = {}
for item in ITEMS[:12]:
    for name, fn in op_table(fresh(item)).items():
        b = LocalBudget()
        try:
            fn(budget=b)
        except Exception:
            continue
        _, e = run(op_table(fresh(item))[name], LocalBudget(work_units=b.used_work // 3))
        out['%%s:%%s' %% (item[0], name)] = [b.used_work, b.peak_bytes,
                                             e.stop.as_dict() if e else None]
print(json.dumps(out, sort_keys=True))
''' % (API, HERE)
        got = []
        for seed in ('0', '1', '12345'):
            env = dict(os.environ, PYTHONHASHSEED=seed, PYTHONDONTWRITEBYTECODE='1')
            p = subprocess.run([sys.executable, '-c', script], capture_output=True,
                               text=True, env=env, timeout=600)
            self.assertEqual(0, p.returncode, p.stderr[-2000:])
            got.append(json.loads(p.stdout))
        self.assertGreater(len(got[0]), 50)
        self.assertEqual(got[0], got[1])
        self.assertEqual(got[0], got[2])


class TestCachesAfterAStop(unittest.TestCase):
    """T-L7: a stop leaves only complete cache entries: the unbudgeted answer afterwards is
    a fresh model's, and the rules still hold."""

    def test_stopped_calls_then_the_unbudgeted_answer(self):
        for item in ITEMS[::3]:
            g = fresh(item)
            names = [n for n in op_table(g) if applicable(fresh(item), n)]
            for n in names:
                b = LocalBudget()
                op_table(fresh(item))[n](budget=b)
                for frac in (0.1, 0.4, 0.7, 0.95):
                    run(op_table(g)[n], LocalBudget(work_units=int(b.used_work * frac)))
            want = TG.op_digests(MT, fresh(item))
            got = TG.op_digests(MT, g)
            self.assertEqual(want, got, item[0])
            self.assertEqual(fresh(item).check_rules(), g.check_rules())


RESUMABLE = ('claims', 'walks_id', 'splits', 'routes', 'label_walks', 'walks_at')


def lists(g):
    out = {'claims': lambda **k: g.claims(strict=False, **k)}
    for s in g.arms:
        n = len(derive.paths(g.arms[s]))
        out[s + ':walks_id'] = lambda s=s, **k: g.walks(s, by='id', **k)
        out[s + ':walks_at'] = lambda s=s, n=n, **k: ops.walks_at(
            g, s, list(range(n))[::-1], with_claims=True, **k)
        out[s + ':splits'] = lambda s=s, **k: g.splits(s, **k)
        for l in range(min(2, len(g.labels))):
            out['%s:routes%d' % (s, l)] = lambda s=s, l=l, **k: g.routes(
                {'id': l}, s, spell=True, **k)
            out['%s:label_walks%d' % (s, l)] = lambda s=s, l=l, **k: g.label_walks(
                {'id': l}, s, **k)
    if g.labels:
        out['label_walks_both'] = lambda **k: g.label_walks({'id': 0}, **k)
    return out


def resume_all(fn, limit):
    """The list made under budgets of |limit| (doubled whenever not one row fits),
    resuming each stop with its token -> (rows, stops, every partial was a prefix)."""
    rows, token, stops = [], None, 0
    while True:
        try:
            return rows + fn(budget=LocalBudget(work_units=limit), resume=token), stops
        except LocalBudgetExceeded as e:
            stops += 1
            assert e.partial is not None
            rows += e.partial.rows
            if e.partial.resume is not None:
                token = e.partial.resume
            if not e.partial.rows:
                limit *= 2
            if stops > 5000:
                raise AssertionError('no progress')


class TestPartialLists(unittest.TestCase):
    """T-L8: a stopped list returns whole rows, a prefix in its documented order, and a
    token that resumes after the last; resuming until done is the whole list."""

    def test_resuming_concatenates_to_the_list(self):
        n = 0
        for item in ITEMS:
            g = fresh(item)
            for name, fn in lists(g).items():
                try:
                    full = fn()
                except Exception:
                    continue
                b = LocalBudget()
                fn(budget=b)
                want = TG.digest(full)
                for frac in (0.07, 0.35, 0.8):
                    got, stops = resume_all(fn, max(1, int(b.used_work * frac)))
                    self.assertEqual(want, TG.digest(got), (item[0], name, frac))
                    n += stops
        self.assertGreater(n, 300)

    def test_partial_rows_are_a_prefix(self):
        g = T.graphlet('merge')
        full = g.claims('right', strict=False)
        b = LocalBudget()
        g.claims('right', strict=False, budget=b)
        for k in range(1, b.used_work, max(1, b.used_work // 40)):
            try:
                g.claims('right', strict=False, budget=LocalBudget(work_units=k))
            except LocalBudgetExceeded as e:
                rows = e.partial.rows
                self.assertEqual([c.as_dict() for c in full[:len(rows)]],
                                 [c.as_dict() for c in rows])
                self.assertEqual(len(rows), e.stop.done['rows'])

    def test_a_token_of_other_arguments_is_refused(self):
        g = T.graphlet('merge')
        b = LocalBudget()
        g.claims('right', strict=False, budget=b)
        token = None
        for k in range(b.used_work // 2, b.used_work):
            try:
                g.claims('right', strict=False, budget=LocalBudget(work_units=k))
            except LocalBudgetExceeded as e:
                if e.partial.rows:
                    token = e.partial.resume
                    break
        self.assertIsNotNone(token)
        with self.assertRaises(ValueError):
            g.claims('right', labels=[{'id': 0}], strict=False, resume=token)
        with self.assertRaises(ValueError):
            g.splits('right', resume=token)
        with self.assertRaises(ValueError):
            T.graphlet('fork').claims('right', strict=False, resume=token)
        self.assertEqual(g.claims('right', strict=False)[-1].as_dict(),
                         g.claims('right', strict=False, resume=token)[-1].as_dict())

    def test_ranked_lists_and_summaries_have_no_partial(self):
        g = T.graphlet('merge')
        for fn in (lambda b: g.walks('right', budget=b),
                   lambda b: ops.rank_walks(g, 'right', budget=b),
                   lambda b: g.label_summary(budget=b), lambda b: g.summary(budget=b),
                   lambda b: g.support_profile('right', 0, budget=b),
                   lambda b: g.support_changes('right', 0, budget=b)):
            with self.assertRaises(LocalBudgetExceeded) as cm:
                fn(LocalBudget(work_units=5))
            self.assertIsNone(cm.exception.partial)
        with self.assertRaises(ValueError):
            g.walks('right', by='support', resume=('walks', 'x', 0))


class TestCompare(unittest.TestCase):
    """T-L9: an interrupted comparison asserts neither equality nor difference."""

    def pairs(self):
        by = {}
        for item in ITEMS:
            by.setdefault(fresh(item).seed.sequence, []).append(item)
        return [(x, y) for xs in by.values() for x, y in zip(xs, xs[1:])]

    def test_stopped_in_every_mode_and_phase(self):
        phases = set()
        n = 0
        for x, y in self.pairs():
            for mode in ('claims', 'walks', 'labels', 'prefix_subset'):
                full = fresh(x).compare(fresh(y), mode=mode)
                b = LocalBudget()
                c = fresh(x).compare(fresh(y), mode=mode, budget=b)
                self.assertEqual(TG.digest(full), TG.digest(c))
                self.assertNotIn('local_stop', c.as_dict())
                for frac in (0, 0.2, 0.5, 0.8, 0.99):
                    lim = int(b.used_work * frac)
                    if lim >= b.used_work:
                        continue
                    c = fresh(x).compare(fresh(y), mode=mode, budget=LocalBudget(work_units=lim))
                    n += 1
                    self.assertEqual(('unknown', None, [], [], []),
                                     (c.comparable, c.equal, c.only_in_a, c.only_in_b,
                                      c.differ))
                    st = c.local_stop
                    self.assertEqual('work', st['resource'])
                    self.assertEqual(st, c.as_dict()['local_stop'])
                    self.assertIn('neither equality nor a difference', c.reason)
                    phases.add(st['phase'])
        self.assertGreater(n, 50)
        self.assertTrue({'keys:a', 'keys:b'} <= phases, phases)

    def test_compare_cost(self):
        for x, y in self.pairs():
            for mode in ('claims', 'walks', 'labels', 'prefix_subset'):
                est = ops.compare_cost(fresh(x), fresh(y), mode=mode)
                b = LocalBudget()
                fresh(x).compare(fresh(y), mode=mode, budget=b)
                w = est['work_units']
                self.assertLessEqual(w['at_least'], b.used_work, (x[0], mode))
                self.assertLessEqual(w['at_least'], w['estimate'])
                self.assertEqual(B.WORK_MODEL, est['work_model'])
                if est['exact']:
                    # nothing keyed (no bases in a mode keyed by them): all it charges
                    self.assertEqual(b.used_work, w['at_least'])
        # incomparable: nothing is keyed, and the estimate is exact
        a, b_ = T.graphlet('fork'), T.graphlet('merge')
        est = ops.compare_cost(a, b_)
        bb = LocalBudget()
        a.compare(b_, budget=bb)
        self.assertTrue(est['exact'])
        self.assertEqual(bb.used_work, est['work_units']['at_least'])


class TestExports(unittest.TestCase):
    """T-L10 / L4: a stopped export returns no text, a stopped save writes no file."""

    def test_no_text_and_no_file(self):
        for item in ITEMS[::2]:
            for name in ('to_json', 'to_fasta', 'to_gfa', 'dump'):
                g = fresh(item)
                if not applicable(g, name):
                    continue
                b = LocalBudget()
                want = op_table(g)[name]()
                self.assertEqual(want, op_table(g)[name](budget=b))
                out, e = run(op_table(g)[name], LocalBudget(work_units=b.used_work // 2))
                self.assertIsNone(out)
                self.assertIn('records', e.stop.done)
            with tempfile.TemporaryDirectory() as d:
                p = os.path.join(d, 'x.mgt')
                b = LocalBudget()
                g.save(p, budget=b)
                with open(p, encoding='utf-8', newline='') as f:
                    self.assertEqual(g.dump(envelope=True), f.read())
                os.unlink(p)
                with self.assertRaises(LocalBudgetExceeded):
                    g.save(p, budget=LocalBudget(work_units=b.used_work - 1))
                self.assertEqual([], os.listdir(d))           # no file, no temporary


def _body(lines_of_g=None):
    return T.doc_text('merge')


class TestParse(unittest.TestCase):
    """T-L11 / §14: an interrupted parse is a local failure: no model, the text kept."""

    def test_admission_before_any_record(self):
        text = T.doc_text('merge')
        n = text.count('\n')
        with self.assertRaises(LocalBudgetExceeded) as cm:
            parse(text, budget=LocalBudget(work_units=B.W_LINE * n - 1))
        st = cm.exception.stop
        self.assertEqual(('parse', 'admission', 0), (st.op, st.phase, st.done['lines']))
        self.assertEqual(n, st.done['of'])
        self.assertGreaterEqual(st.needed_at_least, B.W_LINE * n)
        with self.assertRaises(LocalBudgetExceeded) as cm:
            # the call's base fits; the line list (the text and a header per line) does not
            parse(text, budget=LocalBudget(memory_mb=(B.CALL_BASE + len(text)) / (1 << 20)))
        self.assertEqual(('memory', 'admission'), (cm.exception.stop.resource,
                                                   cm.exception.stop.phase))

    def test_a_stopped_parse_is_no_format_error_and_keeps_the_text(self):
        text = T.doc_text('merge')
        copy_ = str(text)
        b = LocalBudget()
        g = parse(text, budget=b)
        self.assertEqual(text, g.dump())
        for k in range(1, b.used_work, max(1, b.used_work // 25)):
            try:
                parse(text, budget=LocalBudget(work_units=k))
            except LocalBudgetExceeded as e:
                self.assertNotIsInstance(e, GraphletFormatError)
                self.assertLessEqual(e.stop.done['lines'], e.stop.done['of'])
            self.assertEqual(copy_, text)
        self.assertEqual(text, parse(text, budget=LocalBudget(work_units=b.used_work)).dump())

    def test_from_response_and_load(self):
        resp = T.doc_json('merge', 'graphlet')
        r = resp['results'][0]
        b = LocalBudget()
        MT.from_response(r, resp, budget=b)
        with self.assertRaises(LocalBudgetExceeded):
            MT.from_response(r, resp, budget=LocalBudget(work_units=b.used_work // 2))
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, 'g.mgt')
            g = MT.from_response(r, resp)
            g.save(p)
            before = _read_bytes(p)
            b = LocalBudget()
            self.assertEqual(g.dump(envelope=True), MT.load(p, budget=b).dump(envelope=True))
            with self.assertRaises(LocalBudgetExceeded):
                MT.load(p, budget=LocalBudget(work_units=b.used_work - 1))
            self.assertEqual(before, _read_bytes(p))

    def test_adversarial_bodies_stop_within_their_budget(self):
        # a comb of 50,000 segments, and one-token deltas over a wide base (each copies it)
        bodies = [comb_annotate(25000), _wide(400, 3000)]
        for text in bodies:
            parse(text)                                   # valid documents
        for text in bodies:
            for lim in (200_000, 2_000_000):
                t0 = time.perf_counter()
                with self.assertRaises(LocalBudgetExceeded) as cm:
                    parse(text, budget=LocalBudget(work_units=lim))
                dt = time.perf_counter() - t0
                self.assertLessEqual(cm.exception.stop.used, lim)
                # a stop well before the whole parse: 1 lwu is about 0.1 us nominal
                self.assertLess(dt, 2.0 + lim * 1e-6)

    def test_many_labels_charge_their_ids(self):
        text = _wide(20_000, 3)
        b = LocalBudget()
        parse(text, budget=b)
        # 20,000 L records, and each of the 14 sets copies 20,000 ids
        self.assertGreater(b.used_work, 20_000 * B.W_LINE + 14 * 20_000 // 2)


def _wide(n_labels, n_splits):
    """An annotate comb (test_traverse_scale.comb_annotate) whose every set is wide: the
    root records |n_labels| labels, and every child's entry and every P run are the
    one-token delta '!' -- a copy of the whole base set per record."""
    steps = 1 + 2 * n_splits
    counters = ('steps=%d,successor_enumerations=%d,output_bp=%d,pair_evaluations=0,'
                'edge_reuse_probes=0,reminimisation_rounds=0,max_reminimisation_rounds=0,'
                'refusal_scans=0,switch_sources_cut=0') % (steps, steps, steps)
    out = ['H mgt 1 5 basic $ACGT a k k %d 10 0 walk * * 0123456789abcdef' % n_labels,
           'S aaaaaaaaaaaaaaaa 10 6 0 ACGTACGTAC']
    prev = ''
    for i in range(n_labels):
        name = 'l%07d' % i
        k = 0
        while k < min(len(prev), len(name)) and prev[k] == name[k]:
            k += 1
        out.append('L c %d * %d %s' % (i, k, name[k:]))
        prev = name
    out += ['O c c c i',
            'A r c %d p 0 0 1 %d 0 %s * 0 * %d 0 %d %d 0 %d'
            % (n_splits + 1, n_labels, counters, 2 * n_splits + 1, n_splits + 1, n_splits,
               steps),
            'B 0 2 %d %d 1 %d %d 0 %d 0 0 0 0 .' % (n_labels, 2 * n_labels, steps, n_splits,
                                                    n_splits),
            'G * 0 1 0-%d * * * 0 * A' % (n_labels - 1), 'P 0 1 %d !' % n_labels]
    for i in range(1, n_splits + 1):
        parent = 0 if i == 1 else 2 * (i - 1)
        out += ['G %d %d 1 ! * * * * * C' % (parent, i), 'P %d %d %d !' % (i, i + 1, n_labels),
                'T D .',
                'G %d %d 1 ! * * * %s * G' % (parent, i, '0' if i < n_splits else '*'),
                'P %d %d %d !' % (i, i + 1, n_labels)]
        if i == n_splits:
            out.append('T D .')
    out.append('Z %d' % (len(out) + 1))
    return '\n'.join(out) + '\n'


class TestDeadlineAndCancel(unittest.TestCase):
    """T-L13: both stop at a charge point, are non-deterministic, and say so."""

    def test_deadline_zero_stops_at_the_first_clock_check(self):
        g = T.graphlet('merge')
        with self.assertRaises(LocalBudgetExceeded) as cm:
            g.claims(budget=LocalBudget(deadline_s=0))
        st = cm.exception.stop
        self.assertEqual(('deadline', False, 's'), (st.resource, st.deterministic, st.unit))
        self.assertEqual(['retry_later', 'raise_local_budget'], st.actions)

    def test_cancel_from_another_thread(self):
        g = parse(comb_annotate(4000))
        b = LocalBudget()
        out = {}

        def work():
            try:
                g.to_fasta(budget=b)
                out['done'] = True
            except LocalBudgetExceeded as e:
                out['stop'] = e.stop
        t = threading.Thread(target=work)
        b.cancel()                                       # before it starts: the first charge
        t.start()
        t.join(60)
        self.assertEqual(('cancelled', False), (out['stop'].resource, out['stop'].deterministic))
        self.assertEqual(parse(comb_annotate(4000)).to_fasta(), g.to_fasta())


class TestAmbient(unittest.TestCase):
    """T-L14: local_budget() budgets the calls that pass none; an explicit one wins;
    nested blocks restore the outer; a thread starts without one."""

    def test_ambient_explicit_nested_and_threads(self):
        g = T.graphlet('merge')
        with local_budget(work_units=10 ** 9) as outer:
            g.claims()
            used = outer.used_work
            self.assertGreater(used, 0)
            mine = LocalBudget()
            g.claims(budget=mine)
            self.assertEqual(used, outer.used_work)       # the explicit one was charged
            self.assertEqual(used, mine.used_work)
            with local_budget(work_units=1) as inner:
                self.assertIs(inner, B.current())
                with self.assertRaises(LocalBudgetExceeded):
                    g.claims()
            self.assertIs(outer, B.current())
            seen = []
            t = threading.Thread(target=lambda: seen.append(B.current()))
            t.start()
            t.join()
            self.assertEqual([None], seen)
        self.assertIsNone(B.current())
        with B.unbudgeted():
            self.assertIsNone(B.current())

    def test_work_accumulates_and_memory_holds_per_call(self):
        g = T.graphlet('merge')
        one = LocalBudget()
        g.claims(budget=one)
        b = LocalBudget()
        g.claims(budget=b)
        g.claims(budget=b)
        self.assertEqual(2 * one.used_work, b.used_work)
        self.assertEqual(one.peak_bytes, b.peak_bytes)    # the account starts anew per call
        shared = LocalBudget(work_units=one.used_work + one.used_work // 2)
        g.claims(budget=shared)
        with self.assertRaises(LocalBudgetExceeded):
            g.claims(budget=shared)


class TestLimits(unittest.TestCase):
    def test_validation_and_clamps(self):
        for bad in ({'work_units': -1}, {'work_units': 1.5}, {'memory_mb': -2},
                    {'deadline_s': float('inf')}, {'work_units': True}):
            with self.assertRaises(ValueError):
                LocalLimits(**bad)
        top = LocalLimits(100, 10, None)
        self.assertEqual((LocalLimits(50, 10, 3), ['memory_mb']),
                         LocalLimits(50, None, 3).clamped(top))
        self.assertEqual((LocalLimits(100, 10, None), ['work_units', 'memory_mb']),
                         LocalLimits(200, 20).clamped(top))
        self.assertEqual((LocalLimits(1, 1, 1), []), LocalLimits(1, 1, 1).clamped(None))


class TestChargeDensity(unittest.TestCase):
    """T-L6: no heavy loop is left uncharged. Measured without a clock, so that it holds
    on any machine and under any load: the Python lines an operation executes (a trace
    hook counts them) per unit it charges stay under a gate. One lwu is about 0.1 us, a
    few lines of Python; on the committed retrievals every operation runs 0.1-6.3 lines
    per lwu (median 0.3-2.6): a loop that does work it does not charge multiplies its
    lines, not its units."""

    GATE = 16.0

    @staticmethod
    def lines_of(fn):
        n = [0]

        def tracer(frame, event, arg):
            if event == 'line':
                n[0] += 1
            return tracer
        sys.settrace(tracer)
        try:
            fn()
        finally:
            sys.settrace(None)
        return n[0]

    def test_lines_per_unit_under_the_gate(self):
        bad = []
        n = 0
        for item in SAMPLE:
            for name in op_table(fresh(item)):
                if not applicable(fresh(item), name):
                    continue
                b = LocalBudget()
                op_table(fresh(item))[name](budget=b)
                lines = self.lines_of(op_table(fresh(item))[name])
                n += 1
                if lines > self.GATE * b.used_work:
                    bad.append((item[0], name, lines, b.used_work))
        self.assertGreater(n, 200)
        self.assertEqual([], bad)


class TestMemoryModel(unittest.TestCase):
    """T-L12 / L9: the account is a model; on the committed retrievals it is at least the
    heap tracemalloc sees each operation take (cold, a fresh model)."""

    def test_the_account_is_at_least_the_traced_peak(self):
        low = []
        for item in REAL:
            if len(item[1]['graphlet']) < 8000:
                continue                       # every peak there is under the 50 KB looked at
            for name in op_table(fresh(item)):
                if not applicable(fresh(item), name):
                    continue
                b = LocalBudget()
                op_table(fresh(item))[name](budget=b)
                g = fresh(item)
                fn = op_table(g)[name]
                gc.collect()
                tracemalloc.start()
                base = tracemalloc.get_traced_memory()[0]
                res = fn()
                peak = tracemalloc.get_traced_memory()[1] - base
                tracemalloc.stop()
                del res
                if peak > 50_000 and b.peak_bytes < peak:
                    low.append((item[0], name, b.peak_bytes, peak))
        self.assertEqual([], low)


# ======================================================================= the store

class _Client:
    """A backend double serving one recorded response."""

    def __init__(self, resp):
        self.resp = resp

    def build_request(self, seeds, strategy, detail='graphlet'):
        return {'seeds': seeds, 'strategy': strategy}

    def traverse_raw(self, req):
        return json.loads(json.dumps(self.resp))

    def capabilities(self):
        return {}


def _real(name):
    return next(x for x in ITEMS if name in x[0])


def _read_bytes(path):
    with open(path, 'rb') as f:
        return f.read()


class TestStoreParseLimits(unittest.TestCase):
    """T-L16 / L3: a parse that stops keeps the body, checked without a parse, as an
    unparsed entry; it can always be copied out without one."""

    def test_unparsed_entries(self):
        key, r, resp = _real('sra_hub__strict')
        with tempfile.TemporaryDirectory() as d:
            st = GraphletStore(os.path.join(d, 'sp'), parse_limits=LocalLimits(work_units=1000))
            with self.assertRaises(LocalBudgetExceeded) as cm:
                st.put(resp, {'seeds': []})
            h = cm.exception.handle
            e = st.get(h)
            self.assertFalse(e.parsed)
            self.assertEqual({'parsed': False}, {k: v for k, v in e.meta().items()
                                                 if k == 'parsed'})
            self.assertEqual(r['graphlet_bytes'], e.bytes)
            # no model is made from it while parses stop
            with self.assertRaises(LocalBudgetExceeded):
                st.graphlet(h)
            # the escape hatch: the standalone file without a parse, which parses to the
            # model of the response itself
            p = os.path.join(d, 'x.mgt')
            st.save_body(h, p)
            direct = MT.from_response(r, resp)
            self.assertEqual(direct.dump(envelope=True), MT.load(p).dump(envelope=True))
            # a parse that completes (a larger budget) makes it a graphlet like any other
            g = st.graphlet(h, parse_budget=LocalBudget())
            self.assertEqual(direct.dump(), g.dump())
            self.assertTrue(st.get(h).parsed)
            self.assertNotIn('parsed', st.get(h).meta())
            # entries reload with their flag
            st2 = GraphletStore(os.path.join(d, 'sp'))
            self.assertTrue(st2.get(h).parsed)

    def test_a_truncated_body_is_still_refused(self):
        key, r, resp = _real('sra_hub__strict')
        with tempfile.TemporaryDirectory() as d:
            st = GraphletStore(os.path.join(d, 'sp'), parse_limits=LocalLimits(work_units=1000))
            for cut in (lambda b: b[:-30], lambda b: b.replace('Z ', 'Z 9', 1)):
                bad = json.loads(json.dumps(resp))
                bad['results'][0]['graphlet'] = cut(r['graphlet'])
                bad['results'][0].pop('graphlet_bytes', None)
                bad['results'][0].pop('graphlet_lines', None)
                with self.assertRaises(GraphletFormatError):
                    st.put(bad, {'seeds': []})
            bad = json.loads(json.dumps(resp))
            bad['results'][0]['graphlet'] = r['graphlet'][:-30] + '\n'
            with self.assertRaises(GraphletFormatError):
                st.put(bad, {'seeds': []})
            self.assertEqual([], st.list())

    def test_unbudgeted_store_parses_ignore_an_ambient_budget(self):
        key, r, resp = _real('sra_hub__strict')
        with tempfile.TemporaryDirectory() as d:
            st = GraphletStore(os.path.join(d, 'sp'))
            with local_budget(work_units=10) as b:
                h = st.put(resp, {'seeds': []})
                st._drop_ram(h)
                st.graphlet(h)
            self.assertEqual(0, b.used_work)


# ======================================================================= the tools

def _tools(d, limits=ToolLimits(), store_limits=None, resp=None):
    st = GraphletStore(os.path.join(d, 'sp'), parse_limits=store_limits)
    clients = {'idx': _Client(resp)} if resp is not None else {}
    return st, GraphletTools(st, clients, secret=b'k' * 32,
                             export_dir=os.path.join(d, 'x'), local_limits=limits)


BUDGETED_CALLS = {
    'graphlet_summary': lambda t, h, s: t.graphlet_summary(h),
    'graphlet_walks': lambda t, h, s: t.graphlet_walks(h, s, spell='tail'),
    'graphlet_walk': lambda t, h, s: t.graphlet_walk(h, s, 0),
    'graphlet_support': lambda t, h, s: t.graphlet_support(h, s, 0),
    'graphlet_labels': lambda t, h, s: t.graphlet_labels(h),
    'graphlet_splits': lambda t, h, s: t.graphlet_splits(h, s, min_labels_before=0),
    'graphlet_claims': lambda t, h, s: t.graphlet_claims(h, s),
    'graphlet_sequence': lambda t, h, s: t.graphlet_sequence(h, s, walk=0),
    'graphlet_compare': lambda t, h, s: t.graphlet_compare(h, h),
    'graphlet_export': lambda t, h, s: t.graphlet_export(h, format='json', path='e.json'),
    'graphlet_subtrie': lambda t, h, s: t.graphlet_subtrie(h, [{'id': 0}]),
    'graphlet_save': lambda t, h, s: t.graphlet_save(h, 's.mgt'),
}


class TestTools(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = self.tmp.name

    def tearDown(self):
        self.tmp.cleanup()

    def put(self, st, name='sra_hub__strict'):
        key, r, resp = _real(name)
        g = MT.from_response(r, resp)
        return st.put_parsed(g, r['graphlet'], {'seeds': []}), list(g.arms)[0], g

    def test_without_local_limits_no_budget_is_taken(self):
        st = GraphletStore(os.path.join(self.d, 'sp'))
        t = GraphletTools(st, secret=b'k' * 32, export_dir=os.path.join(self.d, 'x'))
        h, s, g = self.put(st)
        out = t.graphlet_claims(h, s, budget={'work_units': 5})
        self.assertEqual('bad_argument', out['error'])
        self.assertNotIn('local', t.graphlet_claims(h, s))

    def test_every_tool_reports_its_local_block(self):
        st, t = _tools(self.d)
        h, s, g = self.put(st)
        for name, call in BUDGETED_CALLS.items():
            out = call(t, h, s)
            self.assertNotIn('error', out, (name, out))
            loc = out['local']
            self.assertTrue(loc['complete'], (name, loc))
            self.assertGreater(loc['usage']['work_units'], 0, name)
            self.assertEqual((1, 'model'), (loc['work_model'], loc['memory_bound']))
            cls = TOOL_CLASS[name]
            self.assertEqual(getattr(ToolLimits(), cls).as_dict(), loc['limits'])
        for name in ('traverse_capabilities', 'graphlet_list', 'graphlet_free'):
            self.assertNotIn(name, TOOL_CLASS)

    def test_tiny_budgets_answer_structured_within_max_bytes(self):
        # T-L17: every budgeted tool, one unit, several ceilings
        for max_bytes in (MIN_MAX_BYTES, 300, 2048):
            st, t = _tools(os.path.join(self.d, str(max_bytes)),
                           ToolLimits(view=LocalLimits(1), heavy=LocalLimits(1)))
            t.max_bytes = max_bytes
            h, s, g = self.put(st)
            for name, call in BUDGETED_CALLS.items():
                out = call(t, h, s)
                ceiling = t.sequence_max_bytes if name == 'graphlet_sequence' else max_bytes
                self.assertLessEqual(_size(out), ceiling, (name, out))
                if name == 'graphlet_compare':
                    if 'error' not in out:
                        self.assertEqual(('unknown', None, None, []),
                                         (out['comparable'], out['equal'], out['counts'],
                                          out['rows']))
                        self.assertFalse(out['local']['complete'])
                    continue
                if out.get('error') == 'receipt_too_large':
                    # a file name or a handle that could not be returned: refused first
                    self.assertIn(name, ('graphlet_export', 'graphlet_save',
                                         'graphlet_subtrie'))
                    continue
                self.assertEqual('local_budget_exceeded', out['error'], (name, out))
                if max_bytes >= 300:
                    self.assertEqual('work', out['stop']['resource'])
                if max_bytes >= 2048:
                    self.assertIn('actions', out['stop'])
                    self.assertNotIn('stop', out['local'])      # stated once
            # nothing was made
            self.assertEqual(1, len(st.list()))
            self.assertEqual([], os.listdir(t._export_root()))

    def test_paging_under_a_budget_that_stops_every_page(self):
        # T-L18: the rows of a stopped list resume after the last, and paging to the end
        # gives exactly the unbudgeted rows
        st0 = GraphletStore(os.path.join(self.d, 'plain'))
        plain = GraphletTools(st0, secret=b'k' * 32, export_dir=os.path.join(self.d, 'x0'))
        for fixture in ('sra_hub__strict', 'mini_ndm1__annotate_exh'):
            h0, s, g = self.put(st0, fixture)
            label = {'ref': g.labels[0].ref}
            cases = {
                'graphlet_claims': lambda tl, h, c, b=None: tl.graphlet_claims(
                    h, s, route_consistent=False, cursor=c, max_bytes=65536, budget=b),
                'graphlet_labels': lambda tl, h, c, b=None: tl.graphlet_labels(
                    h, name=label, cursor=c, max_bytes=65536, budget=b),
                'graphlet_splits': lambda tl, h, c, b=None: tl.graphlet_splits(
                    h, s, min_labels_before=0, cursor=c, max_bytes=65536, budget=b),
                'graphlet_walk': lambda tl, h, c, b=None: tl.graphlet_walk(
                    h, s, 0, cursor=c, max_bytes=65536, budget=b),
                'graphlet_walks': lambda tl, h, c, b=None: tl.graphlet_walks(
                    h, s, rank='id', n=50, cursor=c, max_bytes=65536, budget=b),
            }
            for name, call in cases.items():
                want = self.rows(lambda c, bud: call(plain, h0, c))
                full = call(plain, h0, None)
                # the whole list's usage, then a budget of a fraction of it per page
                d = os.path.join(self.d, '%s-%s' % (fixture, name))
                st, t = _tools(d)
                h = st.put_parsed(MT.from_response(*_real(fixture)[1:]),
                                  _real(fixture)[1]['graphlet'], {'seeds': []})
                use = call(t, h, None)['local']['usage']['work_units']
                for frac in (0.3, 0.6):
                    pages = []
                    got = self.rows(lambda c, bud: call(t, h, c, bud), pages,
                                    max(1, int(use * frac)))
                    self.assertEqual(want, got, (fixture, name, frac))
                    last = pages[-1]
                    self.assertIn('total', last)
                    self.assertTrue(last['local']['complete'] or 'next_cursor' not in last)
                    interrupted = [p for p in pages if p.get('complete') is False]
                    for p in interrupted:
                        self.assertNotIn('total', p) if name != 'graphlet_walks' else None
                        self.assertFalse(p['local']['complete'])
                        self.assertIn('local', p['evidence']['limitations'])
                    if name != 'graphlet_walks':
                        self.assertEqual(full['total'], last['total'])

    @staticmethod
    def rows(page, pages=None, units=None):
        """Every page to the end; each call under budget {work_units: units}, doubled
        while not one row fits (the tool's stop is then the answer)."""
        out, cur = [], None
        for _ in range(5000):
            res = page(cur, None if units is None else {'work_units': units})
            if res.get('error') == 'local_budget_exceeded' and units is not None:
                units *= 2
                continue
            assert 'error' not in res, res
            if pages is not None:
                pages.append(res)
            out += res['rows']
            cur = res.get('next_cursor')
            if not cur:
                return out
        raise AssertionError('no end')

    def test_a_resume_cursor_with_other_arguments(self):
        st0, t0 = _tools(os.path.join(self.d, 'probe'))
        h0, s, g = self.put(st0)
        use = t0.graphlet_claims(h0, s, route_consistent=False,
                                 max_bytes=65536)['local']['usage']['work_units']
        st, t = _tools(self.d, ToolLimits(view=LocalLimits(use * 2 // 3, 1024)))
        h, s, g = self.put(st)
        cur = None
        for _ in range(100):
            out = t.graphlet_claims(h, s, route_consistent=False, cursor=cur, max_bytes=65536)
            cur = out['next_cursor']
            if '"r"' in _cursor_body(cur):
                break                               # a cursor that resumes the claims
        self.assertFalse(out['complete'])
        self.assertEqual('bad_cursor', t.graphlet_claims(h, s, cursor=cur)['error'])
        self.assertEqual('bad_cursor', t.graphlet_splits(h, s, cursor=cur)['error'])
        plain = GraphletTools(st, secret=b'k' * 32)
        self.assertEqual('bad_cursor', plain.graphlet_claims(
            h, s, route_consistent=False, cursor=cur)['error'])

    def test_clamps_and_hooks(self):
        # T-L19
        seen = []
        lim = ToolLimits(on_usage=lambda tool, u: seen.append((tool, u)))
        st, t = _tools(self.d, lim)
        h, s, g = self.put(st)
        out = t.graphlet_claims(h, s, budget={'work_units': 10 ** 12, 'memory_mb': 1})
        self.assertEqual(['work_units'], out['local']['clamped'])
        self.assertEqual(ToolLimits().view_ceiling.work_units,
                         out['local']['limits']['work_units'])
        self.assertEqual(1, out['local']['limits']['memory_mb'])
        out = t.graphlet_claims(h, s, budget={'work_units': 5})
        self.assertEqual('local_budget_exceeded', out['error'])
        self.assertEqual(['graphlet_claims', 'graphlet_claims'], [x[0] for x in seen])
        self.assertTrue(seen[0][1]['complete'])
        self.assertFalse(seen[1][1]['complete'])
        self.assertEqual('work', seen[1][1]['stop']['resource'])
        self.assertEqual('bad_argument', t.graphlet_claims(h, s, budget={'cpu': 1})['error'])
        self.assertEqual('bad_argument', t.graphlet_claims(h, s, budget={'work_units': -1})
                         ['error'])
        # the service's own budget object: one allowance across two calls
        allowance = LocalBudget(work_units=out_usage(t, h, s) * 3 // 2)
        st2, t2 = _tools(os.path.join(self.d, 'b'), ToolLimits(
            budget_for=lambda tool, req: allowance))
        h2 = st2.put_parsed(MT.from_response(*_real('sra_hub__strict')[1:]),
                            _real('sra_hub__strict')[1]['graphlet'], {'seeds': []})
        self.assertNotIn('error', t2.graphlet_labels(h2))
        self.assertEqual('local_budget_exceeded', t2.graphlet_labels(h2)['error'])

    def test_receipts(self):
        # T-L21: a stopped export / subtrie / save makes nothing; a completed receipt
        # keeps its essentials, the local block cut first
        st, t = _tools(self.d)
        h, s, g = self.put(st)
        for name in ('graphlet_export', 'graphlet_subtrie', 'graphlet_save'):
            out = BUDGETED_CALLS[name](t, h, s)
            self.assertIn('local', out)
        st2, t2 = _tools(os.path.join(self.d, 'b'))
        h2, s2, g2 = self.put(st2)
        t2.max_bytes = 200
        out = t2.graphlet_save(h2, 'small.mgt')
        self.assertIn('path', out)
        self.assertEqual('local', out['fields_cut'][0])
        for name, args in (('graphlet_export', dict(format='fasta', path='f.fa')),
                           ('graphlet_subtrie', dict(labels=[{'id': 0}]))):
            st3, t3 = _tools(os.path.join(self.d, name),
                             ToolLimits(heavy=LocalLimits(50, 1024)))
            h3, s3, g3 = self.put(st3)
            out = getattr(t3, name)(h3, **args)
            self.assertEqual('local_budget_exceeded', out['error'])
            self.assertEqual(1, len(st3.list()))
            self.assertEqual([], os.listdir(t3._export_root()))

    def test_body_level_copies_need_no_parse(self):
        st, t = _tools(self.d)
        h, s, g = self.put(st)
        t0 = GraphletTools(GraphletStore(os.path.join(self.d, 'p')), secret=b'k' * 32,
                           export_dir=os.path.join(self.d, 'x0'))
        h0 = t0.store.put_parsed(MT.from_response(*_real('sra_hub__strict')[1:]),
                                 _real('sra_hub__strict')[1]['graphlet'], {'seeds': []})
        for tools, hh in ((t, h), (t0, h0)):
            tools.graphlet_save(hh, 'same.mgt')
            tools.graphlet_export(hh, format='mgt', path='same2.mgt')
        for f in ('same.mgt', 'same2.mgt'):
            with open(os.path.join(t._export_root(), f), 'rb') as a, \
                    open(os.path.join(t0._export_root(), f), 'rb') as b:
                self.assertEqual(a.read(), b.read())

    def test_a_fetch_whose_parse_stops_keeps_the_body(self):
        key, r, resp = _real('sra_hub__strict')
        st, t = _tools(self.d, ToolLimits(parse=LocalLimits(2000, 1024)), resp=resp)
        out = t.traverse_fetch('idx', seed={'sequence': 'ACGT'})
        self.assertFalse(out['parsed'])
        self.assertEqual('local_budget_exceeded'.split('_')[0], 'local')
        self.assertFalse(out['local']['complete'])
        self.assertEqual('parse', out['local']['stop']['op'])
        self.assertIn('export_mgt', out['local']['stop']['actions'])
        h = out['handle']
        if 'server_summary' in out:
            self.assertEqual(r['outcome'], out['server_summary']['outcome'])
        else:
            self.assertIn('server_summary', out['fields_cut'])
        self.assertFalse(st.get(h).parsed)
        # every local tool on it needs a parse first: it stops under the parse class
        got = t.graphlet_summary(h)
        self.assertEqual('local_budget_exceeded', got['error'])
        self.assertEqual('parse', got['stop']['op'])
        # the escape hatch: the body copied out without a parse
        got = t.graphlet_export(h, format='mgt', path='kept.mgt')
        self.assertTrue(got['local']['complete'], got)
        g = MT.load(os.path.join(t._export_root(), 'kept.mgt'))
        self.assertEqual(MT.from_response(r, resp).dump(envelope=True), g.dump())
        # with a parse budget that fits it, the entry is a graphlet like any other
        t.local_limits = ToolLimits()
        got = t.graphlet_summary(h)
        self.assertNotIn('error', got)
        self.assertIn('parsed', got['local'])
        self.assertTrue(st.get(h).parsed)

    def test_a_load_whose_parse_stops_stores_nothing(self):
        st, t = _tools(self.d, ToolLimits(parse=LocalLimits(500, 1024)))
        root = t._export_root()
        p = os.path.join(root, 'in.mgt')
        T.graphlet('merge').save(p)
        before = _read_bytes(p)
        out = t.graphlet_load('in.mgt')
        self.assertEqual('local_budget_exceeded', out['error'])
        self.assertEqual([], st.list())
        self.assertEqual(before, _read_bytes(p))

    def test_a_stopped_comparison_is_unknown(self):
        st, t = _tools(self.d, ToolLimits(heavy=LocalLimits(3000, 1024)))
        h, s, g = self.put(st)
        out = t.graphlet_compare(h, h)
        self.assertEqual(('unknown', None, None, []),
                         (out['comparable'], out['equal'], out['counts'], out['rows']))
        self.assertFalse(out['local']['complete'])
        for tag in ('a', 'b'):
            self.assertEqual(['local_work'], out['evidence'][tag]['limitations']['local'])


def _cursor_body(cur):
    import base64
    b64 = cur.rsplit('.', 1)[0]
    return base64.urlsafe_b64decode(b64 + '=' * (-len(b64) % 4)).decode()


def out_usage(t, h, s):
    return t.graphlet_labels(h)['local']['usage']['work_units']


class TestRealOfflineUnderDefaults(unittest.TestCase):
    """T-L22: every committed real retrieval completes every view tool under the default
    ToolLimits (the defaults catch what is pathological, not what is real)."""

    def test_every_cell_every_view_tool(self):
        with tempfile.TemporaryDirectory() as d:
            st, t = _tools(d)
            for key, r, resp in REAL:
                g = MT.from_response(r, resp)
                h = st.put_parsed(g, r['graphlet'], {'seeds': []})
                for s in g.arms:
                    for name, call in BUDGETED_CALLS.items():
                        if TOOL_CLASS[name] != 'view':
                            continue
                        out = call(t, h, s)
                        if 'error' in out:
                            self.assertNotEqual('local_budget_exceeded', out['error'],
                                                (key, name, out))
                        else:
                            self.assertTrue(out['local']['complete'], (key, name))


if __name__ == '__main__':
    unittest.main()
