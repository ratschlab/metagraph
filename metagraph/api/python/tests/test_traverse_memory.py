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
        # label's supported prefixes besides -- 8-9x here, 26x with every spelling held
        for mode, bound in (('claims', 6), ('prefix_subset', 12)):
            x, y = self.g(), self.g()
            peak, c = _peak(lambda: x.compare(y, mode=mode))
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
