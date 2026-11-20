"""The memory of the unbudgeted exports and of compare() on comb-shaped graphlets (where
holding every spelling or a copy of the text multiplies the peak: to_fasta() at about 10x
its text, compare() several times higher than a stream), with their answers unchanged:

  * derive.walk_iter() streams walk_batch(): the same targets, chains and bases in the same
    order, a kept prefix dropped after its last dependent, no chain kept without chains=
    (8 bytes per segment, 8x the bases on a comb of 1-bp segments), a missing base raised
    before the first target;
  * to_fasta(), to_gfa(), to_json() and dump() hold their text and one walk at a time, and
    join their final line feed instead of copying the text once more;
  * compare()'s keys are cut from a stream, and prefix_subset's divergence and recorded
    refusals follow a sequence down the chains' tree instead of keeping every spelling --
    the same answers as spelling each walk on its own (restated here);
  * the local-limits account still bounds the traced peak of each on the combs;
  * walk_iter() reads the targets' chains only, never the whole arm, for a missing base
    (50x the time of 4 shallow walks on an 8,000-segment comb otherwise), and save() and the
    MCP tools' writes stay within their accounts (a buffered file's 128 KiB buffer,
    uncharged, would put save() above its account on bodies of 10-40 KB with dump()'s join
    charged at 1x).

The digests of every public operation over the committed retrievals are checked by the
golden gate (test_traverse_speedups.py); these tests add the shapes where the memory is
at stake and the paths no committed comparison reaches.
"""

import gc
import glob
import os
import random
import sys
import tempfile
import tracemalloc
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402
import traverse_golden as TG  # noqa: E402
from test_traverse_scale import comb_annotate  # noqa: E402

import metagraph.traverse as MT  # noqa: E402
from metagraph.traverse import GraphletStore, LocalBudget, derive, ops, parse  # noqa: E402
from metagraph.traverse._codec import REASON  # noqa: E402
from metagraph.traverse.mcp_tools import GraphletTools, ToolLimits  # noqa: E402


def _items():
    files = sorted(glob.glob(os.path.join(T.DOCS, '**', '*.graphlet.json'), recursive=True))
    return TG.items_of(files + TG.golden_files())


ITEMS = _items()


def _fresh(item):
    return MT.from_response(item[1], item[2])


def _comb(n):
    """comb_annotate(n) with the annotate fixture's envelope (to_json() and compare() need
    one): n splits of 1-bp segments, the walks' text quadratic in n."""
    resp = T.doc_json('annotate', 'graphlet')
    r = dict(resp['results'][0])
    r.pop('graphlet_bytes', None)
    r.pop('graphlet_lines', None)
    r['graphlet'] = comb_annotate(n)
    return r, resp


def _peak(fn):
    gc.collect()
    tracemalloc.start()
    base = tracemalloc.get_traced_memory()[0]
    try:
        out = fn()
        return tracemalloc.get_traced_memory()[1] - base, out
    finally:
        tracemalloc.stop()


def _warm_peak(warm, fn):
    """_peak() of fn after warm() ran once, the collector held off from before warm() to the
    end. CPython keeps freed small objects on per-size free lists, memory tracemalloc still
    counts as allocated where they were first made, and a full collection (_peak()'s
    gc.collect(), or an automatic one) empties the lists: a cold peak includes refilling
    them, and what that costs depends on the interpreter, not on the library. On Python
    3.11 a keyword call of the slots dataclass Claim (20 fields, the tp_call path) packs two
    20-item tuples that are parked there: 1,604 tuples of 200 bytes on the comb of 800
    splits, 1.0x its leaves' bases, which 3.12 and later do not allocate (their calls go
    through vectorcall). Warm, the peak is the library's own, the same on every run."""
    gc.collect()
    gc.disable()
    try:
        warm()
        tracemalloc.start()
        base = tracemalloc.get_traced_memory()[0]
        out = fn()
        return tracemalloc.get_traced_memory()[1] - base, out
    finally:
        tracemalloc.stop()
        gc.enable()


class TestWalkIter(unittest.TestCase):
    def test_the_stream_is_walk_batch(self):
        rng = random.Random(55)
        n = 0
        for item in ITEMS:
            g = _fresh(item)
            for a in g.arms.values():
                ids = [s.id for s in a.segments]
                if not ids:
                    continue
                for targets in (derive.leaves(a), rng.sample(ids, min(5, len(ids))), ids):
                    for spell in (True, False):
                        for chains in (True, False):
                            want = derive.walk_batch(a, targets, spell, chains)
                            got = list(derive.walk_iter(a, targets, spell, chains))
                            self.assertEqual([(t, c, w) for t, (c, w) in want.items()], got)
                            self.assertEqual(sorted(set(targets)), [t for t, _, _ in got])
                            if chains:
                                self.assertEqual(len(got), len({id(c) for _, c, _ in got}))
                    n += 1
        self.assertGreater(n, 50)

    def test_a_missing_base_is_raised_before_the_first_target(self):
        saved = derive._EACH_STEPS
        try:
            # target by target, as a batch, and target by target falling back to the batch
            for steps in (0, saved, 10 ** 9):
                derive._EACH_STEPS = steps
                a = parse(comb_annotate(40)).arms['right']
                leaves = derive.leaves(a)
                x = leaves[20]                      # a leaf in the middle of the comb
                fine = [l for l in leaves if l != x]
                want = [(l, None, derive.walk_bases(a, l)) for l in fine]
                a.segments[x].walk = None
                for targets in ([x], leaves, [s.id for s in a.segments]):
                    # a few targets raise from the call itself, more from the first next()
                    with self.assertRaises(ValueError):
                        next(derive.walk_iter(a, targets, chains=False))
                # the targets whose chains have their bases are spelled
                self.assertEqual(want, list(derive.walk_iter(a, fine, chains=False)))
        finally:
            derive._EACH_STEPS = saved

    def test_without_check_first_the_error_comes_at_its_target(self):
        saved = derive._EACH_STEPS
        try:
            a = parse(comb_annotate(40)).arms['right']
            leaves = derive.leaves(a)
            x = leaves[20]
            a.segments[x].walk = None
            some = leaves[15:25]
            # target by target: the targets before |x| are given, then its error
            derive._EACH_STEPS = 10 ** 9
            got = []
            with self.assertRaises(ValueError):
                for t, _, w in derive.walk_iter(a, some, chains=False, check_first=False):
                    got.append((t, w))
            self.assertEqual([(t, derive.walk_bases(a, t)) for t in some if t < x], got)
            with self.assertRaises(ValueError):
                next(derive.walk_iter(a, some, chains=False))
            # a few targets: from the call, built at once
            with self.assertRaises(ValueError):
                derive.walk_iter(a, [some[0], x], chains=False, check_first=False)
            self.assertEqual([(some[0], None, derive.walk_bases(a, some[0]))],
                             list(derive.walk_iter(a, [some[0]], chains=False)))
            # as a batch the union shows it before the first target either way
            derive._EACH_STEPS = 0
            with self.assertRaises(ValueError):
                next(derive.walk_iter(a, some, chains=False, check_first=False))
        finally:
            derive._EACH_STEPS = saved

    def test_the_chains_are_read_not_the_arm(self):
        # a scan of every segment for a missing base would cost 4 shallow walks of an
        # 8,000-segment comb 140-170 us against 2.7 us; only the targets' chains are read
        # (the arm is never iterated)
        class Unread(list):
            def __iter__(self):
                raise AssertionError('the whole arm was read')

        a = parse(comb_annotate(400)).arms['right']
        leaves = derive.leaves(a)
        cases = (leaves[:4], leaves[:6], leaves[-4:], leaves)
        want = [[(t, derive.walk_bases(a, t)) for t in sorted(c)] for c in cases]
        a.segments = Unread(a.segments)
        for c, w in zip(cases, want):
            for check_first in (True, False):
                self.assertEqual(w, [(t, b) for t, _, b in derive.walk_iter(
                    a, c, chains=False, check_first=check_first)])
            self.assertEqual(dict(w), {t: b for t, (_, b) in derive.walk_batch(a, c).items()})

    def test_the_stream_holds_a_few_prefixes_on_a_comb(self):
        g = parse(comb_annotate(3000))
        a = g.arms['right']
        leaves = derive.leaves(a)
        total = sum(a.segments[x].end_bp for x in leaves)

        def consume():
            n = 0
            for _, _, w in derive.walk_iter(a, leaves, chains=False):
                n += len(w)
            return n
        peak, n = _peak(consume)
        self.assertEqual(total, n)
        # the union's set and counts (a few hundred KB) and a few walks of 3 kb, against
        # 4.5 MB of bases (and 36 MB of chains kept at every spine segment before)
        self.assertLess(peak, total // 4)
        peak_chains, _ = _peak(
            lambda: sum(1 for _ in derive.walk_iter(a, leaves, spell=False)))
        self.assertLess(peak_chains, total // 2)


class TestExportPeaks(unittest.TestCase):
    """On a comb of 800 splits (0.3 MB of FASTA): each export's peak is about twice its
    text -- the records and the joined text -- where holding every walk peaked at 10x and
    adding the final line feed by copying the text at 3.1x."""

    @classmethod
    def setUpClass(cls):
        cls.r, cls.resp = _comb(800)

    def g(self):
        return MT.from_response(self.r, self.resp)

    def test_to_fasta(self):
        for kw in ({}, {'width': 60}, {'orientation': 'walk', 'with_seed': False}):
            g = self.g()
            peak, text = _peak(lambda: g.to_fasta(**kw))
            self.assertGreater(len(text), 250_000)
            self.assertLess(peak, 2.6 * len(text), kw)
        # walks in another order and repeated: each record in its slot
        g = self.g()
        n = len(derive.paths(g.arms['right']))
        order = list(range(n))[::-1][:200] + [0, 1, 2, 0]
        want = ''.join(g.to_fasta('right', leaves=[x]) for x in order)
        self.assertEqual(want, self.g().to_fasta('right', leaves=order))
        self.assertEqual(g.to_fasta('right', leaves=[0, 5, 9]),
                         ''.join(g.to_fasta('right', leaves=[x]) for x in (0, 5, 9)))

    def test_to_gfa_and_dump(self):
        g = self.g()
        peak, text = _peak(g.to_gfa)
        self.assertLess(peak, 2.6 * len(text))
        # dump(): its lines are short (22 bytes on average here), so each line's header
        # (49 bytes) and slot weigh beside the text
        g = self.g()
        peak, text = _peak(g.dump)
        self.assertLess(peak, 2 * len(text) + 80 * text.count('\n'))
        self.assertTrue(text.endswith('\n') and not text.endswith('\n\n'))

    def test_to_json(self):
        # the output holds every path's chain (8 bytes per segment): the bound is the
        # chains' slots and the text the JSON would have
        g = self.g()
        steps = sum(len(p.segments) for p in derive.paths(g.arms['right']))
        g = self.g()
        peak, out = _peak(g.to_json)
        self.assertEqual(steps, sum(len(p['segments']) for p in out['arms']['right']['paths']))
        self.assertLess(peak, 2 * 9 * steps + 4_000_000)

    def test_compare(self):
        a = self.g()
        bases = sum(s.end_bp for s in a.arms['right'].segments if s.leaf is not None)
        # claims: the keys' prefixes (each claim's bases, in the keys of both sides) and
        # the claims behind them -- 4x the leaves' bases here, 14x with every spelling held
        # and a chain kept at every spine segment; prefix_subset: the walks cut and every
        # label's supported prefixes besides -- 8-9x here, 26x with every spelling held.
        # Measured warm (_warm_peak()): 5.60x and 10.30x on Python 3.11, 5.30x and 9.71x on
        # 3.13, 5.31x and 9.90x on 3.14 (cold, as before, 3.11 counted its free lists'
        # refill too, 6.85x and 11.20x: over the bound)
        for mode, bound in (('claims', 6), ('prefix_subset', 12)):
            x, y = self.g(), self.g()
            peak, c = _warm_peak(lambda: self.g().compare(self.g(), mode=mode),
                                 lambda: x.compare(y, mode=mode))
            self.assertTrue(c.equal is not False, mode)
            self.assertLess(peak, bound * bases, mode)

    def test_the_account_bounds_the_peak(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'x.mgt')
            ops_ = {'to_fasta': lambda g, h, **k: g.to_fasta(**k),
                    'to_gfa': lambda g, h, **k: g.to_gfa(**k),
                    'to_json': lambda g, h, **k: g.to_json(**k),
                    'dump': lambda g, h, **k: g.dump(**k),
                    'save': lambda g, h, **k: g.save(path, **k),
                    'compare': lambda g, h, **k: g.compare(h, **k),
                    'prefix_subset': lambda g, h, **k: g.compare(h, mode='prefix_subset',
                                                                 **k)}
            for name, fn in ops_.items():
                b = LocalBudget()
                fn(self.g(), self.g(), budget=b)
                g, h = self.g(), self.g()
                peak, _ = _peak(lambda: fn(g, h))
                self.assertGreaterEqual(b.peak_bytes, peak, name)


class TestSeedAndSelections(unittest.TestCase):
    """The lane's sweep of the GPT final review (M2, M4): an export that writes the seed
    charges it, and a call given a list of walks admits its copies of the list before it
    makes them -- allocation measured before the stop, not only the stop."""

    @staticmethod
    def _long_seed(name, n=10 ** 6):
        """The fixture with |n| A's inserted in the middle of its seed: its first and last
        15 bases (more than k - 1) unchanged."""
        import copy
        import dataclasses
        g = copy.deepcopy(T.graphlet(name))
        h = len(g.seed.sequence) // 2
        s = g.seed.sequence[:h] + 'A' * n + g.seed.sequence[h:]
        g.seed = dataclasses.replace(g.seed, sequence=s, length_bp=len(s),
                                     num_kmers=len(s) - g.k + 1)
        return g

    def test_a_long_seed_is_charged(self):
        """to_gfa() writes the seed's S line (uncharged: 26,912 B charged against 3.0 MB
        traced for a 1 Mbp seed on the fork fixture) and to_json() spells each continuation
        from seed and flank (the whole seed copied per leaf, 1.0 MB against 60 KB); now the
        S line is charged with its copy in the joined text, and a continuation copies only
        the bases it takes from the seed. The account bounds the warm traced peak; the
        output is the one of the same graphlet's short seed where the seed is not written."""
        import copy
        for name in ('fork', 'merge', 'annotate'):
            g0 = self._long_seed(name)
            for op, fn in (('to_gfa', lambda g, **k: g.to_gfa(**k)),
                           ('to_json', lambda g, **k: g.to_json(**k)),
                           ('to_fasta', lambda g, **k: g.to_fasta(**k))):
                with self.subTest(name, op=op):
                    b = LocalBudget()
                    got = fn(copy.deepcopy(g0), budget=b)
                    self.assertEqual(got, fn(copy.deepcopy(g0)))
                    xs = [copy.deepcopy(g0), copy.deepcopy(g0)]
                    peak, _ = _warm_peak(lambda: fn(xs.pop()), lambda: fn(xs.pop()))
                    self.assertLessEqual(peak, b.peak_bytes)
            # a segment's context reads at most k - 1 bases of the seed, its last (its first
            # on the left arm): the same as on the graphlet's own seed, whose ends the long
            # one keeps
            short, long_ = T.graphlet(name), self._long_seed(name)
            self.assertLess(short.k - 1, len(short.seed.sequence) // 2)
            gfa = [x for x in long_.to_gfa().split('\n') if not x.startswith('S\tseed')]
            self.assertEqual(gfa, [x for x in short.to_gfa().split('\n')
                                   if not x.startswith('S\tseed')])

    def test_a_list_of_walks_is_admitted_before_it_is_copied(self):
        """to_fasta(leaves=<list>) resolved the list, its set of leaves and price lists, and
        ops.walks_at() copied its ids (the list, the token's arguments, the slice) before
        the first charge: 1,000,000 ids allocated 33.8 MB and 31 MB before a 1 MiB budget
        stopped them. Now the copies are admitted first: the traced peak before the stop
        stays below 64 KB under 1 MiB, and under 2,000 lwu."""
        import copy
        ids = [0] * 10 ** 6
        g = T.graphlet('fork')
        for lim in ({'memory_mb': 1}, {'work_units': 2000}):
            for name, fn in (
                    ('to_fasta', lambda h: h.to_fasta('right', leaves=ids,
                                                      budget=LocalBudget(**lim))),
                    ('walks_at', lambda h: ops.walks_at(h, 'right', ids,
                                                        budget=LocalBudget(**lim)))):
                with self.subTest(name, **lim):
                    xs = [copy.deepcopy(g), copy.deepcopy(g)]

                    def call():
                        try:
                            fn(xs.pop())
                        except MT.LocalBudgetExceeded:
                            return 'stop'
                    peak, out = _warm_peak(call, call)
                    self.assertEqual(out, 'stop')
                    self.assertLess(peak, 64 * 1024)
        # within budget: the same answers
        sel = [1, 0, 1, 0]
        b = LocalBudget()
        self.assertEqual(g.to_fasta('right', leaves=sel, budget=b),
                         g.to_fasta('right', leaves=sel))
        self.assertEqual([w.path_id for w in ops.walks_at(g, 'right', sel, budget=b)],
                         [w.path_id for w in ops.walks_at(g, 'right', sel)])
        self.assertEqual([w.path_id for w in ops.walks_at(g, 'right', iter(sel), budget=b)],
                         sel)
        # a stop while the ids are admitted is still a list's stop: no rows, and the call's
        # own resume= (None for a fresh call, the token it was given when resumed); the
        # walks are made in chunks of 32, so a stop with rows needs more ids than that
        sel = [1, 0] * 40
        want = [w.path_id for w in ops.walks_at(g, 'right', sel)]
        full = LocalBudget()
        ops.walks_at(g, 'right', sel, budget=full)
        token = None
        for k in range(1, full.used_work):
            try:
                ops.walks_at(g, 'right', sel, budget=LocalBudget(work_units=k))
            except MT.LocalBudgetExceeded as e:
                if e.partial.rows:
                    first, token = [w.path_id for w in e.partial.rows], e.partial.resume
                    break
        self.assertIsNotNone(token)
        for head, resume in (([], None), (first, token)):
            with self.assertRaises(MT.LocalBudgetExceeded) as cm:
                ops.walks_at(g, 'right', sel, budget=LocalBudget(work_units=0), resume=resume)
            self.assertEqual(cm.exception.partial.rows, [])
            self.assertEqual(cm.exception.partial.resume, resume)
            self.assertEqual(head + [w.path_id for w in ops.walks_at(g, 'right', sel,
                                                                     resume=resume)], want)


# ------------------------------------------------------------------ the comparison's scans
# restated without streaming: each walk spelled on its own

def _divergence_reference(arm, seq, depth):
    best = 0
    for leaf in ops.restricted_leaves(arm, depth):
        try:
            w = derive.walk_bases(arm, leaf)[:depth]
        except ValueError:
            return None
        n = 0
        for x, y in zip(w, seq):
            if x != y:
                break
            n += 1
        best = max(best, n)
    return best


def _refusal_reference(g, arm, seq, ref):
    lid = next((l.id for l in g.labels if l.ref == ref), None)
    if lid is None:
        return None

    def on_chain(seg, at):
        if at >= len(seq):
            return False
        try:
            w = derive.walk_bases(arm, seg)
        except ValueError:
            return False
        return len(w) >= at and w[:at] == seq[:at]

    for be in arm.branch_events:
        if not on_chain(be.segment, be.at_bp):
            continue
        for r in be.refused:
            if r.char == seq[be.at_bp] and lid in r.labels:
                return be.at_bp, r.cause
    for s in arm.segments:
        for ev in s.events:
            if ev.type == 'blocked' or (ev.type == 'hairpin' and not ev.followed):
                if on_chain(s.id, ev.at_bp) and ev.char == seq[ev.at_bp] and lid in ev.labels:
                    return ev.at_bp, REASON[ev.reason] if ev.type == 'blocked' else 'hairpin'
    return None


def _probes(g, a, rng):
    """Sequences in walking order that leave the arm at its recorded refusals (the refused
    base after the chain's bases), that follow its walks and part from them, and random
    ones; with the labels to ask about."""
    segs = a.segments
    out = []
    events = [(be.segment, be.at_bp, r.char, r.labels)
              for be in a.branch_events for r in be.refused]
    events += [(s.id, ev.at_bp, ev.char, ev.labels) for s in segs for ev in s.events
               if ev.type == 'blocked' or (ev.type == 'hairpin' and not ev.followed)]
    for seg, at, char, labels in events:
        w = derive.walk_bases(a, seg)[:at]
        tail = ''.join(rng.choice('ACGT') for _ in range(rng.randint(0, 20)))
        refs = [g.labels[l].ref for l in list(labels)[:3]] + [g.labels[0].ref]
        out.append((w + char + tail, refs))
        if w:
            i = rng.randrange(len(w))
            out.append((w[:i] + 'ACGT'[('ACGT'.index(w[i]) + 1) % 4] + w[i + 1:] + char,
                        refs))
    for leaf in derive.leaves(a):
        w = derive.walk_bases(a, leaf)
        cut = rng.randint(0, len(w))
        out.append((w[:cut] + ''.join(rng.choice('ACGT') for _ in range(5)),
                    [g.labels[rng.randrange(len(g.labels))].ref]))
        out.append((w + 'A', [g.labels[0].ref]))
    out.append(('', [g.labels[0].ref]))
    out.append((''.join(rng.choice('ACGT') for _ in range(50)), [g.labels[0].ref]))
    return out


class TestComparisonScans(unittest.TestCase):
    def test_divergence_and_refusals_as_spelled_walk_by_walk(self):
        rng = random.Random(2024)
        hits = n = 0
        for item in ITEMS:
            g = _fresh(item)
            if not g.labels:
                continue
            for a in g.arms.values():
                if not a.segments or a.segments[0].walk is None:
                    continue
                memo = {}
                depths = sorted({0, 1, a.complete_to_bp, max(0, a.complete_to_bp // 2),
                                 a.complete_to_bp + 7})
                for seq, refs in _probes(g, a, rng):
                    for depth in depths:
                        want = _divergence_reference(a, seq, depth)
                        self.assertEqual(want, ops._divergence(a, seq, depth, memo),
                                         (item[0], a.side, depth))
                        self.assertEqual(want, ops._divergence(a, seq, depth), item[0])
                    for ref in refs:
                        want = _refusal_reference(g, a, seq, ref)
                        self.assertEqual(want, ops._recorded_refusal(g, a, seq, ref, memo),
                                         (item[0], a.side, seq, ref))
                        self.assertEqual(want, ops._recorded_refusal(g, a, seq, ref))
                        hits += want is not None
                        n += 1
        self.assertGreater(n, 500)
        self.assertGreater(hits, 20)

    def test_without_bases_on_a_chain(self):
        g = T.graphlet('quorum')
        a = g.arms['right']
        segs = {be.segment for be in a.branch_events}
        self.assertTrue(segs)
        seq, refs = _probes(g, a, random.Random(1))[0]
        self.assertIsNotNone(ops._recorded_refusal(g, a, seq, refs[0]))
        for s in a.segments:
            s.walk = None
        self.assertIsNone(_refusal_reference(g, a, seq, refs[0]))
        self.assertIsNone(ops._recorded_refusal(g, a, seq, refs[0]))
        self.assertIsNone(ops._divergence(a, seq, a.complete_to_bp))
        self.assertEqual(_divergence_reference(a, seq, 0), ops._divergence(a, seq, 0))


class TestCallerSizedMetadata(unittest.TestCase):
    """G4P's sweep left these (T7/G4/VERIFY-python.md G5-G14): text and JSON values whose
    size their sender sets -- a seed summary's work block, a stop's actions and message,
    label names, a seed id, a namespace or release, a strategy's lists, an envelope -- were
    copied, serialized, escaped or encoded beyond what was charged (up to 900x the account,
    an astral name compared). Each is now charged before it is made: the account bounds the
    warm traced peak, the answer is the unbudgeted one, and a store's parse stops before it
    reads a body it cannot hold."""

    N = 200_000

    @staticmethod
    def _named(name, names):
        import dataclasses
        g = T.graphlet(name)
        g.labels = [dataclasses.replace(l, name=names[i]) if i < len(names) else l
                    for i, l in enumerate(g.labels)]
        g.cache.clear()
        for a in g.arms.values():
            a.cache.clear()
        return g

    def _bounded(self, make, fn, same=True):
        b = LocalBudget()
        got = fn(make(), b)
        if same:
            self.assertEqual(got, fn(make(), None))
        xs = [make(), make()]
        peak, _ = _warm_peak(lambda: fn(xs.pop(), None), lambda: fn(xs.pop(), None))
        self.assertLessEqual(peak, b.peak_bytes)

    def test_to_json_and_summary(self):
        import copy
        import dataclasses
        n = self.N

        def work():
            g = T.graphlet('merge')
            g.seed_summary = dict(g.seed_summary, work={'x': [[i] for i in range(n // 10)]})
            return g

        def stop(message='m', k=10_000):
            def make():
                g = T.graphlet('limits')
                g.resource_stop = dataclasses.replace(
                    g.resource_stop, message=message, actions=['a%05d' % i for i in range(k)])
                return g
            return make

        def seed_id():
            g = T.graphlet('merge')
            g.seed_summary = copy.deepcopy(g.seed_summary)
            g.seed_summary['seed']['seed_id'] = 'q' * n
            return g

        self.assertIsNotNone(T.graphlet('limits').resource_stop)
        cases = [('to_json work', work, lambda g, b: g.to_json(budget=b)),
                 ('to_json actions', stop(), lambda g, b: g.to_json(budget=b))]
        for name, make in (
                ('names', lambda: self._named('merge', ['N%d' % i * (n // 8)
                                                        for i in range(4)])),
                ('astral names', lambda: self._named('merge', ['\U0001F600' * (n // 4)] * 4)),
                ('seed_id', seed_id), ('stop', stop('m' * n, 2_000))):
            for cap in (1 << 24, 2048):
                cases.append(('summary %s %d' % (name, cap), make,
                              lambda g, b, cap=cap: g.summary(max_bytes=cap, budget=b)))
        for name, make, fn in cases:
            with self.subTest(name):
                self._bounded(make, fn)

    def test_next_request(self):
        import copy
        big = ['lab%07d' % i + 'y' * 90 for i in range(10_000)]

        def kw(b):
            # (budget=None would be a keyword override of the request's loss budget)
            return {} if b is None else {'budget': b}
        for fx in ('merge', 'annotate'):
            def env_list(fx=fx):
                g = T.graphlet(fx)
                g.envelope = copy.deepcopy(g.envelope)
                g.envelope['strategy']['labels']['extra'] = list(big)
                return g
            for name, make, fn in (
                    ('extra override', lambda fx=fx: T.graphlet(fx),
                     lambda g, b: g.next_request('right', [0], labels={'extra': big}, **kw(b))),
                    ('note override', lambda fx=fx: T.graphlet(fx),
                     lambda g, b: g.next_request('right', [0], note=big, **kw(b))),
                    ('envelope list', env_list,
                     lambda g, b: g.next_request('right', [0], **kw(b))),
                    ('envelope list, per walk', env_list,
                     lambda g, b: g.next_requests('right', [0], **kw(b)))):
                with self.subTest(fx, case=name):
                    self._bounded(make, fn)

    def test_compare_and_subgraph(self):
        n = self.N
        for name, x, y in (('names', 'A' * n, 'B' * n),
                           ('astral names', '\U0001F600' * n, '\U0001F601' * n)):
            with self.subTest(name):
                self._bounded(lambda: (self._named('merge', [x]), self._named('merge', [y])),
                              lambda gg, b: gg[0].compare(gg[1], budget=b).notes)
        for name, x, y in (('ns', 'n' * n, 'm' * n), ('release', 'r' * n, 's' * n)):
            def pair(x=x, y=y, name=name):
                out = []
                for v in (x, y):
                    g = T.graphlet('merge')
                    if name == 'ns':
                        g.index_ns = v
                    else:
                        g.envelope = dict(g.envelope, release=v)
                    out.append(g)
                return out
            with self.subTest(name):
                self._bounded(pair, lambda gg, b: gg[0].compare(gg[1], budget=b).notes)

        from metagraph.traverse.model import Label

        def many():
            g = T.graphlet('merge')
            k = len(g.labels)
            g.labels = g.labels + [Label(k + i, 'column', 1000 + i, None, 'L%d' % i)
                                   for i in range(20_000)]
            return g
        for sel in ([{'id': i} for i in range(20_000)],
                    [{'name': 'L%d' % i} for i in range(20_000)]):
            with self.subTest('subgraph', by=list(sel[0])[0]):
                self._bounded(many, lambda g, b: ops.subgraph(g, sel, budget=b, of='x').selectors)

    def test_the_writer_and_the_reader(self):
        n = self.N

        def env(value):
            def make():
                g = T.graphlet('merge')
                g.envelope = dict(g.envelope, algorithm_version=value)
                return g
            return make
        cases = [('J ints', env(list(range(300, 300 + n // 4))),
                  lambda g, b: g.dump(envelope=True, budget=b)),
                 ('J nested', env([{'k': [i]} for i in range(n // 20)]),
                  lambda g, b: g.dump(envelope=True, budget=b))]
        for cname, name in (('ascii', 'a' * n), ('ascii with %', 'a%' * (n // 2)),
                            ('bmp', '一' * n), ('astral', '\U0001F600' * n)):
            cases.append(('dump %s' % cname, lambda name=name: self._named('merge', [name]),
                          lambda g, b: g.dump(budget=b)))
            text = self._named('merge', [name]).dump(envelope=False)
            cases.append(('parse %s' % cname, lambda text=text: text,
                          lambda t, b: parse(t, budget=b).labels[0].name))

        def explicit_c():
            # a continuation whose explicit sequence differs from the one the walk spells
            # (G17): the writer spells it to compare, here across a seed of n bases
            g = TestSeedAndSelections._long_seed('merge', n)
            arm = g.arms['right']
            s = next(x for x in arm.segments if x.leaf is not None
                     and x.leaf.continuation is not None and x.walk is not None)
            c = s.leaf.continuation
            c.n = len(g.seed.sequence) + s.end_bp
            c.sequence = 'C' * c.n
            return g
        cases.append(('dump explicit C', explicit_c, lambda g, b: g.dump(budget=b)))
        for name, make, fn in cases:
            with self.subTest(name):
                self._bounded(make, fn)
        # G17's comparison spelling, pinned by its work: an explicit sequence equal to the
        # one the walk spells writes the same text as none, and is charged the n bases it
        # spells to compare (nothing was, so the peak above holds on the base too: the
        # spelling is dropped before the join)
        from metagraph.traverse import parser as P

        def leaf(g):
            arm = g.arms['right']
            return arm, next(x for x in arm.segments if x.leaf is not None
                             and x.leaf.continuation is not None and x.walk is not None)
        same, none = explicit_c(), explicit_c()
        arm, s = leaf(same)
        s.leaf.continuation.sequence = P._derived_continuation(same, arm, s)
        leaf(none)[1].leaf.continuation.sequence = None
        b1, b2 = LocalBudget(), LocalBudget()
        self.assertEqual(same.dump(budget=b1), none.dump(budget=b2))
        self.assertGreaterEqual(b1.used_work - b2.used_work,
                                (s.end_bp + s.leaf.continuation.n) >> 9)

    def test_mixed_widths(self):
        """One character beyond ASCII widens a whole joined str to four bytes a character,
        the ASCII text beside it included: summary(), compare() and the J line charge all
        of it so (Opus review of b8fb7a8d6, M2: 2.6-3.1x the account before), and the
        store's parse reads a wide body only once its text fits."""
        n = self.N
        astral = '\U0001F600'

        def seed_id_and_astral():
            import copy
            g = self._named('merge', ['x', astral])
            g.seed_summary = copy.deepcopy(g.seed_summary)
            g.seed_summary['seed']['seed_id'] = 'q' * n
            return g

        def releases(x, y):
            def make():
                out = []
                for v in (x, y):
                    g = T.graphlet('merge')
                    g.envelope = dict(g.envelope, release=v)
                    out.append(g)
                return out
            return make

        def env(value):
            def make():
                g = T.graphlet('merge')
                g.envelope = dict(g.envelope, algorithm_version=value)
                return g
            return make
        for name, make, fn in (
                ('summary seed_id and an astral name', seed_id_and_astral,
                 lambda g, b: g.summary(max_bytes=1 << 24, budget=b)),
                ('summary a long ASCII name and an astral one',
                 lambda: self._named('merge', ['a' * n, astral]),
                 lambda g, b: g.summary(max_bytes=1 << 24, budget=b)),
                ('compare ASCII and astral releases', releases('r' * n, astral),
                 lambda gg, b: gg[0].compare(gg[1], budget=b).notes),
                ('J line ASCII and astral', env(['x' * n, astral]),
                 lambda g, b: g.dump(envelope=True, budget=b)),
                ('J line ASCII only', env('x' * n),
                 lambda g, b: g.dump(envelope=True, budget=b))):
            with self.subTest(name):
                self._bounded(make, fn)
        # under a limit: stopped before the wide text is made
        for mb in (3, 5, 7):
            g = env(['x' * n * 5, astral])()

            def call():
                try:
                    g.dump(envelope=True, budget=LocalBudget(memory_mb=mb))
                except MT.LocalBudgetExceeded:
                    return 'stop'
            peak, _ = _warm_peak(call, call)
            with self.subTest('J line under a limit', mb=mb):
                self.assertLess(peak, mb << 20)
        with tempfile.TemporaryDirectory() as d:
            st = GraphletStore(os.path.join(d, 'spool'), max_ram_mb=0)
            resp = T.doc_json('merge', 'graphlet')
            r = dict(resp['results'][0])
            r['graphlet'] = self._named('merge', ['a' * 5 * n, astral]).dump(envelope=False)
            r.pop('graphlet_bytes', None)
            r.pop('graphlet_lines', None)
            h = st.put(dict(resp, results=[r]), {'seeds': []})
            for mb in (2, 6, 12):
                def call():
                    try:
                        st.graphlet(h, parse_budget=LocalBudget(memory_mb=mb))
                    except MT.LocalBudgetExceeded:
                        return 'stop'
                peak, out = _warm_peak(call, call)
                with self.subTest('store parse of a wide body', mb=mb):
                    self.assertEqual(out, 'stop')
                    self.assertLess(peak, mb << 20)
            # the body without a parse, at its width (Opus check of 0b2965d96, N2: 24 MB
            # traced against 6 MB charged for this body, and a completed call 3x over its
            # limit), and a wide J line beside an ASCII body
            p = os.path.join(d, 'out.mgt')
            rj = dict(resp, algorithm_version=['x' * n, astral])
            hj = st.put(rj, {'seeds': []})
            for name, fn in (('standalone_text', lambda b, h_=h: st.standalone_text(h_, budget=b)),
                             ('save_body', lambda b, h_=h: st.save_body(h_, p, budget=b)),
                             ('standalone_text, wide J line',
                              lambda b: st.standalone_text(hj, budget=b))):
                with self.subTest(name):
                    b = LocalBudget()
                    fn(b)
                    peak, _ = _warm_peak(lambda: fn(None), lambda: fn(None))
                    self.assertLessEqual(peak, b.peak_bytes)
                for mb in (8, 12):
                    def call():
                        try:
                            fn(LocalBudget(memory_mb=mb))
                        except MT.LocalBudgetExceeded:
                            return 'stop'
                    peak, _ = _warm_peak(call, call)
                    with self.subTest(name, mb=mb):
                        self.assertLess(peak, mb << 20)
        # parse() of UTF-8 bytes (a gzip file's load): the digest's chunk at the text's width
        # (a 1 MB body), and the split charged before it is made, so that a limit that holds
        # the decoded text and not its lines stops before the copy (a wide line of 3M
        # characters: 24 MB traced under 20 MB when charged after)
        for name, names, mbs in (('digest', ['a' * n * 5, astral], range(2, 16)),
                                 ('split', ['a' * n * 15 + astral], range(16, 24, 2))):
            data = self._named('merge', names).dump().encode('utf-8')
            for mb in mbs:
                def call():
                    try:
                        parse(data, budget=LocalBudget(memory_mb=mb))
                    except MT.LocalBudgetExceeded:
                        return 'stop'
                peak, _ = _warm_peak(call, call)
                with self.subTest('parse of wide bytes', case=name, mb=mb):
                    self.assertLess(peak, mb << 20)

    def test_cyclic_values(self):
        """A caller's value that holds itself (a cycle no JSON can hold) is sized without
        walking it twice: a budgeted call answers as the unbudgeted one -- the same
        ValueError, or the same request -- where the walk never ended (Opus review of
        b8fb7a8d6, M1). Each budget has a work limit, so that a walk that does not end
        stops (a failure) instead of growing until killed."""
        def lim():
            return LocalBudget(work_units=10**7)
        cyc = ['a']
        cyc.append(cyc)
        g = T.graphlet('merge')
        for budget in (None, lim()):
            kw = {} if budget is None else {'budget': budget}
            with self.subTest('next_request labels.extra', budgeted=budget is not None):
                with self.assertRaisesRegex(ValueError, 'labels.extra is a list of label names'):
                    g.next_request('right', [0], labels={'extra': cyc}, **kw)
        out = g.next_requests('right', [0], note=cyc,
                              budget=LocalBudget(memory_mb=64, work_units=10**7))
        self.assertEqual(repr(out), repr(g.next_requests('right', [0], note=cyc)))
        h = T.graphlet('merge')
        h.envelope = dict(h.envelope, algorithm_version=cyc)
        for budget in (None, lim()):
            with self.subTest('dump(envelope=True)', budgeted=budget is not None):
                with self.assertRaisesRegex(ValueError, 'Circular reference detected'):
                    h.dump(envelope=True, budget=budget)
        from metagraph.traverse import budget as B
        self.assertEqual(B.value_size(cyc, LocalBudget(work_units=10**6))[:4],
                         (3, 1, B.list_bytes(2), 9))

    def test_shared_values(self):
        """A caller's value that holds one container more than once, without a cycle, is
        priced as deepcopy copies it: each container once (its memo). 41 lists, each holding
        the next twice, are 2^40 paths: a walk per path never ended under a budget, where
        the unbudgeted call answers at once (Opus check of 0b2965d96, N1; a regression
        from dada7f430). Work-limited, so that a regression fails instead of hanging."""
        from metagraph.traverse import budget as B
        x = ['a']
        for _ in range(40):
            x = [x, x]
        g = T.graphlet('merge')

        def plain(out):
            # the requests without the note (whose repr() writes every path), and the
            # note's copy: as deepcopy makes it, each level one list held twice
            notes = [r['strategy'].pop('note') for r in out]
            for n in notes:
                for _ in range(40):
                    self.assertIs(n[0], n[1])
                    n = n[0]
                self.assertEqual(n, ['a'])
            return repr(out)
        want = plain(g.next_requests('right', [0], note=x))
        for kw in ({}, {'memory_mb': 64}):
            b = LocalBudget(work_units=10**7, **kw)
            with self.subTest(**kw):
                self.assertEqual(plain(g.next_requests('right', [0], note=x, budget=b)), want)
                self.assertLess(b.used_work, 10**4)
        b = LocalBudget(work_units=10**6)
        self.assertEqual(B.value_size(x, b, shared=True)[1:3],
                         (41, 40 * B.list_bytes(2) + B.list_bytes(1)))

    def test_values_shared_between_overrides(self):
        """One value given under several keyword overrides is copied once per key: the merge
        into the strategy deep-copies each override (or each entry of a section the strategy
        holds) with a memo of its own, so each is charged (Opus check of 63e7fdaa3, R1: one
        list as note0..note7 charged one copy, 7.5x under its traced peak; 3.75x for four
        entries of the annotation section)."""
        x = ['a'] * 100_000
        cases = (('8 keys', {'note%d' % i: x for i in range(8)}),
                 ('a section', {'annotation': {'n%d' % i: x for i in range(4)}}))
        for name, kw in cases:
            for op in ('next_request', 'next_requests'):
                def call(b=None, kw=kw, op=op):
                    g = T.graphlet('merge')
                    extra = {} if b is None else {'budget': b}
                    return getattr(g, op)('right', [0], **extra, **kw)
                with self.subTest(name, op=op):
                    b = LocalBudget(work_units=10**7)
                    call(b)
                    peak, _ = _warm_peak(call, call)
                    self.assertLessEqual(peak, b.peak_bytes)

    def test_many_overrides_grow_the_strategy(self):
        """The merge writes every override into the strategy (and into the merged copy the
        labels.extra rebuild reads): those dicts grow by its entries, charged with them (Opus
        check of feda889fd, F1: 200,000 scalar entries charged as none, 660x under in a
        section, 1,562x at the top level). The factor left (about 3x in a section, 7-9x at
        the top level with the call's own **kwargs) is the dict model's (F2, a follow-up:
        DICT_KEY's 16 bytes against about 38 for a dict built by insertion), so the bound
        here is a floor, not traced <= account."""
        n = 50_000
        flat = {'k%d' % i: i for i in range(n)}
        for name, kw, floor in (('a section', {'annotation': flat}, 4),
                                ('the top level', dict(flat), 10)):
            def call(b=None, kw=kw):
                g = T.graphlet('merge')
                return g.next_request('right', [0], **({} if b is None else {'budget': b}),
                                      **kw)
            with self.subTest(name):
                b = LocalBudget(work_units=10**8)
                call(b)
                peak, _ = _warm_peak(call, call)
                self.assertGreaterEqual(b.peak_bytes * floor, peak)

    def test_charges_are_the_same_on_every_interpreter(self):
        """No charge reads the interpreter's Unicode tables or object sizes (Opus review of
        b8fb7a8d6, M4): a character of Unicode 15 (U+1FAE8, unknown to 3.11's tables) and
        a short Latin-1 text are charged these figures on 3.11 and 3.14 alike. The figures
        are literals: a change of these charges moves them -- re-pin them from a run on
        both interpreters, which must agree (a figure only one of them gives is the bug
        this test is for)."""
        from metagraph.traverse import budget as B
        new = '\U0001FAE8'
        b = LocalBudget()
        self._named('merge', [new * 10000]).summary(budget=b)
        self.assertEqual((b.peak_bytes, b.used_work), (179808, 17738))
        ga, gb = T.graphlet('merge'), T.graphlet('merge')
        ga.envelope = dict(ga.envelope, release=new * 10000)
        gb.envelope = dict(gb.envelope, release='r')
        b = LocalBudget()
        ga.compare(gb, budget=b)
        self.assertEqual(b.peak_bytes, 818928)
        for k, want in ((30, 16477), (34, 16485), (38, 16493)):
            b = LocalBudget()
            try:
                parse('\u00e9' * k, budget=b)
            except MT.GraphletFormatError:
                pass
            self.assertEqual(b.peak_bytes, want, k)
        self.assertEqual(B.text_bytes(new * 1000), 4008)
        self.assertEqual(B.repr_bytes(new * 1000), 40080)
        self.assertEqual(B.value_size({'extra': [{1, 2, 3}, b'abc', object()]}),
                         (6, 2, 1099, 211, False, 4))

    def test_the_store(self):
        """standalone_text()'s J line with the encoder's pieces, and a parse on demand or a
        load that reads its body only once the parse's budget can hold it: a zero budget
        stops with under 64 KB traced, where it read and decoded the body first (2 MB)."""
        import copy
        n = self.N
        with tempfile.TemporaryDirectory() as d:
            st = GraphletStore(os.path.join(d, 'spool'), max_ram_mb=0)
            resp = T.doc_json('merge', 'graphlet')
            r = dict(resp['results'][0])
            r['graphlet'] = self._named('merge', ['a' * 5 * n]).dump(envelope=False)
            r.pop('graphlet_bytes', None)
            r.pop('graphlet_lines', None)
            h = st.put(dict(resp, results=[r]), {'seeds': []})
            path = os.path.join(d, 'x.mgt')
            parse(r['graphlet']).save(path)
            for name, fn in (('graphlet', lambda b: st.graphlet(h, parse_budget=b)),
                             ('view', lambda b: st.view(h, parse_budget=b)),
                             ('load', lambda b: st.load_graphlet(path, parse_budget=b))):
                with self.subTest(name):
                    def call():
                        try:
                            fn(LocalBudget(work_units=0, memory_mb=0))
                        except MT.LocalBudgetExceeded:
                            return 'stop'
                    peak, out = _warm_peak(call, call)
                    self.assertEqual(out, 'stop')
                    self.assertLess(peak, 64 * 1024)
            rr = copy.deepcopy(resp)
            rr['algorithm_version'] = list(range(300, 300 + n // 4))
            h2 = st.put(rr, {'seeds': []})
            b = LocalBudget()
            self.assertEqual(st.standalone_text(h2, budget=b), st.standalone_text(h2))
            peak, _ = _warm_peak(lambda: st.standalone_text(h2), lambda: st.standalone_text(h2))
            self.assertLessEqual(peak, b.peak_bytes)
            # and its work: the J line's encoder per atom (b8fb7a8d6; 9.8x under before)
            self.assertGreaterEqual(b.used_work, 2 * (n // 4))


class TestWriteAccounts(unittest.TestCase):
    """The local-limits account of save() and of the MCP tools' writes is at least their
    traced peak: the files are written without a buffer (one write), and the JSON export's
    buffered file has a stated size that is charged. A buffered file's
    io.DEFAULT_BUFFER_SIZE (128 KiB from Python 3.14) would otherwise be held uncharged."""

    def cases(self):
        out = [('comb%d' % n, lambda n=n: MT.from_response(*_comb(n))) for n in (100, 200)]
        out += [(it[0], lambda it=it: _fresh(it)) for it in ITEMS]
        return out

    def test_save(self):
        low = []
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'x.mgt')
            parse(comb_annotate(3)).save(path)        # the first call's own allocations
            n = 0
            for name, make in self.cases():
                b = LocalBudget()
                make().save(path, budget=b)
                g = make()
                peak, _ = _peak(lambda: g.save(path))
                n += 1
                if b.peak_bytes < peak:
                    low.append((name, b.peak_bytes, peak))
        self.assertGreater(n, 20)
        self.assertEqual([], low)

    def test_the_tools_writes(self):
        low = []
        with tempfile.TemporaryDirectory() as d:
            st = GraphletStore(os.path.join(d, 'sp'))
            tools = GraphletTools(st, {}, secret=b'k' * 32, export_dir=os.path.join(d, 'x'),
                                  local_limits=ToolLimits())
            n = 0
            for name, make in self.cases()[::3]:
                g = make()
                h = st.put_parsed(g, g.dump(), {'seeds': []})
                for what, call in (
                        ('mgt', lambda: tools.graphlet_export(h, format='mgt', path='e.mgt')),
                        ('json', lambda: tools.graphlet_export(h, format='json',
                                                               path='e.json')),
                        ('save', lambda: tools.graphlet_save(h, 's.mgt'))):
                    out = call()                       # the first call's own allocations
                    if 'error' in out:
                        continue                       # json without an envelope
                    peak, out = _peak(call)
                    n += 1
                    account = out['local']['usage']['memory_bytes']
                    if account < peak:
                        low.append((name, what, account, peak))
                st.free(h)
        self.assertGreater(n, 20)
        self.assertEqual([], low)


if __name__ == '__main__':
    unittest.main()
