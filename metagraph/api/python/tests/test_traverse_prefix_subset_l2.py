"""prefix_subset's lookups and charges (review of 2026-10-06, L2): a's recorded refusals
are indexed once per comparison and only those on the segments a walk reaches are tested;
the best claim of an omission is found by following the walk down a's claim chains; the
violation scan bisects b's sorted supported prefixes. Every answer is the one the scans
gave (the scans are kept here as the oracle, as the library had them), and the charges
follow the work: on the review's comb of 8,001 segments the comparison ran ~70x longer
than it charged.

The graphlets: the committed fixtures and synthetic ones with refusals injected where the
oracle's every branch is reached -- V refusals and blocked or skipped successors on
segments of a comb and a bushy tree, at the segment's start, inside it, past its end, and
before it (an event filed under an ancestor's bases) -- and claims moved below their own
segment, whose anchor is an ancestor too."""

import dataclasses
import os
import random
import sys
import time
import unittest
from array import array

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402
from test_traverse_scale import comb_annotate  # noqa: E402
from test_traverse_stage_l import _wide  # noqa: E402

import traverse_golden as TG  # noqa: E402
from metagraph.traverse import LocalBudget, derive, from_response, model, ops  # noqa: E402
from metagraph.traverse._codec import REASON  # noqa: E402
from metagraph.traverse.parser import parse  # noqa: E402


# ------------------------------------------------------------------ the scans (oracle)

def old_recorded_refusal(g, arm, seq, ref):
    """_recorded_refusal() as it scanned every refusal of the arm (6897db99), its chain
    test spelled out (the chain's bases [0, at) equal seq's)."""
    ids = ops._label_index(g)[0].get(ref)
    if not ids:
        return None
    lid = ids[0]
    segs = arm.segments

    def on_chain(seg, at):
        if at >= len(seq):
            return False
        try:
            w = derive.walk_bases(arm, seg)
        except ValueError:
            return False
        if len(w) < at:
            return False
        return w[:at] == seq[:at]

    for be in arm.branch_events:
        if not on_chain(be.segment, be.at_bp):
            continue
        for r in be.refused:
            if r.char == seq[be.at_bp] and lid in r.labels:
                return be.at_bp, r.cause
    for s in segs:
        for ev in s.events:
            if ev.type == 'blocked' or (ev.type == 'hairpin' and not ev.followed):
                if on_chain(s.id, ev.at_bp) and ev.char == seq[ev.at_bp] and lid in ev.labels:
                    return ev.at_bp, REASON[ev.reason] if ev.type == 'blocked' else 'hairpin'
    return None


def old_best_claim(arm, rows, seq, ref):
    """The best claim as the scan found it: the longest prefix, the first among equals."""
    best = None
    for i, c in enumerate(rows):
        if c.label.ref != ref:
            continue
        prefix = derive.walk_bases(arm, c.segment)[:c.to_bp]
        if seq.startswith(prefix) and (best is None or len(prefix) > best[0]):
            best = (len(prefix), i, c)
    return best


# ------------------------------------------------------------------ graphlets

def comb_with_refusals(n):
    """comb_annotate(n), label 1 on the spine, with a V refusal of 'C' for label 1 at every
    split (the review's r3)."""
    a = parse(comb_annotate(n))
    arm = a.arms['right']
    segs = arm.segments
    for i in range(1, n + 1):
        parent = 0 if i == 1 else 2 * (i - 1)
        arm.branch_events.append(model.BranchEvent(
            segs[parent].end_bp, parent, 'CG', [1, 2], None, None,
            [model.Refusal('C', 'quorum', array('I', [1]))]))
    return a


def inject(g, rng, n_events):
    """Refusals at random places of g's arms: V refusals and blocked or skipped
    successors, at a segment's start, inside it, at and past its end, and before it."""
    for arm in g.arms.values():
        segs = arm.segments
        if not segs or segs[0].walk is None:
            continue
        nl = len(g.labels)
        for _ in range(n_events):
            s = rng.randrange(len(segs))
            seg = segs[s]
            where = rng.choice(('start', 'inside', 'end', 'past', 'before'))
            at = {'start': seg.from_bp, 'inside': seg.from_bp + rng.randrange(
                max(1, seg.length_bp)), 'end': seg.end_bp, 'past': seg.end_bp + 3,
                'before': rng.randrange(max(1, seg.from_bp + 1))}[where]
            char = rng.choice('ACGT')
            labels = array('I', sorted(rng.sample(range(nl), min(nl, rng.randint(1, 3)))))
            kind = rng.choice(('V', 'blocked', 'hairpin', 'hairpin_followed'))
            if kind == 'V':
                arm.branch_events.append(model.BranchEvent(
                    at, s, char + 'A', [1, 1], None, None,
                    [model.Refusal(char, rng.choice(('quorum', 'split_limit')), labels)]))
            else:
                ev = model.Event(at, 'blocked' if kind == 'blocked' else 'hairpin')
                ev.char = char
                ev.labels = labels
                ev.total = len(labels)
                if kind == 'blocked':
                    ev.reason = next(iter(REASON))
                else:
                    ev.followed = kind == 'hairpin_followed'
                seg.events.append(ev)
    return g


def walks_of(g, side, rng, extra=8):
    """The arm's walks in walking order, their prefixes and a few mutated ones."""
    arm = g.arms[side]
    out = set()
    for p in derive.paths(arm):
        w = derive.walk_bases(arm, p.leaf)
        out.add(w)
        for j in range(0, len(w), max(1, len(w) // 5)):
            out.add(w[:j])
    base = sorted(out)
    for _ in range(extra):
        if not base:
            break
        w = rng.choice(base)
        if w:
            j = rng.randrange(len(w))
            out.add(w[:j] + rng.choice('ACGT') + w[j + 1:])
    return sorted(out)


def graphlets():
    rng = random.Random(20261006)
    out = []
    for name in sorted(T.DOCUMENTS):
        try:
            g = T.graphlet(name)
        except Exception:
            continue
        if any(a.segments and a.segments[0].walk is None for a in g.arms.values()):
            continue
        out.append(('doc:' + name, g))
        out.append(('doc+events:' + name, inject(T.graphlet(name), rng, 12)))
    out.append(('comb_refusals', comb_with_refusals(40)))
    out.append(('comb+events', inject(parse(comb_annotate(30)), rng, 40)))
    out.append(('wide+events', inject(parse(_wide(6, 12)), rng, 30)))
    # the recorded retrievals of the golden gate (real indexes), as they are and with
    # refusals injected
    for key, r, resp in TG.items_of(TG.golden_files())[:12]:
        g = from_response(r, resp)
        if any(a.segments and a.segments[0].walk is None for a in g.arms.values()):
            continue
        out.append(('real:' + key, g))
        out.append(('real+events:' + key, inject(from_response(r, resp), rng, 25)))
    return out


LABELS_PER_GRAPHLET = 30
WALKS_PER_ARM = 80


def labels_of(g):
    """The graphlet's labels, at most LABELS_PER_GRAPHLET of them, spread over its
    dictionary (the oracle spells a chain per refusal: the loop is kept small)."""
    step = max(1, len(g.labels) // LABELS_PER_GRAPHLET)
    return g.labels[::step][:LABELS_PER_GRAPHLET]


class TestTheIndexAnswersAsTheScan(unittest.TestCase):
    def test_recorded_refusals(self):
        rng = random.Random(7)
        n = hits = 0
        for name, g in graphlets():
            for side, arm in g.arms.items():
                memo = {}
                for seq in walks_of(g, side, rng)[:WALKS_PER_ARM]:
                    for lab in labels_of(g):
                        want = old_recorded_refusal(g, arm, seq, lab.ref)
                        got = ops._recorded_refusal(g, arm, seq, lab.ref, memo)
                        self.assertEqual(want, got, (name, side, seq, lab.ref))
                        n += 1
                        hits += want is not None
        self.assertGreater(n, 5000)
        self.assertGreater(hits, 100)        # the refusals are found, not only missed

    def test_the_best_claim(self):
        rng = random.Random(11)
        n = found = 0
        for name, g in graphlets():
            for side, arm in g.arms.items():
                depth = arm.complete_to_bp
                rows = [c for c in ops._restricted_claims(g, side, depth)
                        if c.kind != 'merged']
                if not rows:
                    continue
                index = ops._claim_index(arm, rows)
                for seq in walks_of(g, side, rng)[:WALKS_PER_ARM]:
                    matched = ops._follow(arm.segments, index[0], index[1], seq, len(seq))
                    for lab in labels_of(g):
                        want = old_best_claim(arm, rows, seq, lab.ref)
                        got = ops._best_claim(index, matched, lab.ref)
                        self.assertEqual(want and want[:2], got and got[:2],
                                         (name, side, seq, lab.ref))
                        n += 1
                        found += want is not None
        self.assertGreater(n, 5000)
        self.assertGreater(found, 1000)

    def test_the_best_claim_of_claims_before_their_segment(self):
        # The review of the L2 fix: a claim whose to_bp lies before its own segment (the
        # model does not exclude one) is filed under the ancestor that holds its last base,
        # found by one walk of the tree now (_chain_anchors()), not a climb per claim. Such
        # claims are made here by moving a claim down its first-parent chain: its prefix,
        # the chain's bases [0, to_bp), stays the same, and the oracle spells it
        rng = random.Random(23)
        n = moved_n = found = 0
        for name, g in graphlets():
            for side, arm in g.arms.items():
                depth = arm.complete_to_bp
                rows = [c for c in ops._restricted_claims(g, side, depth)
                        if c.kind != 'merged']
                if not rows:
                    continue
                kids = {}
                for x in arm.segments:
                    if x.parents:
                        kids.setdefault(x.parents[0], []).append(x.id)
                moved = []
                for c in rows:
                    d = c.segment
                    for _ in range(rng.randrange(4)):
                        if not kids.get(d):
                            break
                        d = rng.choice(kids[d])
                    if d != c.segment and c.to_bp > 0:
                        moved_n += 1
                    moved.append(dataclasses.replace(c, segment=d))
                index = ops._claim_index(arm, moved)
                for seq in walks_of(g, side, rng)[:WALKS_PER_ARM // 2]:
                    matched = ops._follow(arm.segments, index[0], index[1], seq, len(seq))
                    for lab in labels_of(g)[:LABELS_PER_GRAPHLET // 2]:
                        want = old_best_claim(arm, moved, seq, lab.ref)
                        got = ops._best_claim(index, matched, lab.ref)
                        self.assertEqual(want and want[:2], got and got[:2],
                                         (name, side, seq, lab.ref))
                        n += 1
                        found += want is not None
        self.assertGreater(moved_n, 300)
        self.assertGreater(n, 2000)
        self.assertGreater(found, 300)

    def test_the_anchors_are_those_of_the_climb(self):
        # _chain_anchors() against the climb it replaced, on every segment and offset of
        # the chains of a bushy and a deep arm (equal offsets included: a segment of 0
        # bases starts where its child does)
        for g in (parse(_wide(6, 12)), parse(comb_annotate(30)), T.graphlet('switch_chain')):
            for arm in g.arms.values():
                segs = arm.segments
                if not segs or segs[0].walk is None:
                    continue
                tree = ops._chain_tree(segs, {x.id for x in segs})
                kids, offs, roots = tree
                asks = [(x.id, at) for x in segs for at in range(0, offs[x.id] + 1)]

                def climb(seg, at):
                    a = seg
                    while offs[a] > at:
                        a = segs[a].parents[0]
                    return a
                self.assertEqual([climb(x, at) for x, at in asks],
                                 ops._chain_anchors(kids, roots, offs, asks))


# ------------------------------------------------------------------ charges

def only_label0(nl, ns):
    """The review's r8: a's segments carry label 0 only (its dictionary still names all nl
    labels), so every (walk, label != 0) of a shallow wide b is an omission without a
    claim, and each asked a's refusals by scanning every segment."""
    res = []
    for line in _wide(nl, ns).splitlines():
        p = line.split(' ')
        if p[0] == 'G' and p[1] == '*':
            p[4] = '0'
        if p[0] == 'P':
            p[3] = '1'
        res.append(' '.join(p))
    return '\n'.join(res) + '\n'


class TestCharges(unittest.TestCase):
    def test_a_deep_a_is_not_scanned_per_omission(self):
        # 32,001 segments and 2,000 labels: 11,994 omissions. The scan visited every segment
        # per omission (384M visits, uncharged: 47-72x the nominal time of the charge,
        # 11-17 s; now 0.6x, 0.15 s)
        a, b = parse(only_label0(2000, 16000)), parse(_wide(2000, 5))
        visits = [0]
        real = ops._follow

        def counted(*args, **kw):
            got = real(*args, **kw)
            visits[0] += len(got)
            return got
        ops._follow = counted
        try:
            bud = LocalBudget()
            t = time.process_time()
            c = ops.compare(a, b, mode='prefix_subset', budget=bud)
            cpu = time.process_time() - t
        finally:
            ops._follow = real
        self.assertGreater(len(c.only_in_b), 10000)
        # every segment a walk was followed through is charged (at least W_ELEM each)
        self.assertLessEqual(visits[0] * ops.W_ELEM, bud.used_work)
        # and the time stays within a small factor of the charge's nominal time (1 lwu ~
        # 0.1 us); generous, for a loaded machine
        self.assertLess(cpu, 5 * bud.used_work * 1e-7 + 0.2, (cpu, bud.used_work))

    def test_events_before_their_segment_are_not_climbed_per_event(self):
        # The review of the L2 fix: each refusal whose at_bp lies before its own segment
        # (the parser does not check it; a crafted or corrupt body) was filed by a climb of
        # its chain, one uncharged step per ancestor and per event: 8,000 such events on
        # the deepest segment of a 16,001-segment comb ran 8.8 s against 0.45 s nominal
        # (19.6x; 16,000 on 32,001 segments 44x, and a 2 s deadline stopped only after
        # 9.96 s). Now one walk of the tree, charged
        nl, ns, k = 20, 8000, 8000
        a = parse(only_label0(nl, ns))
        segs = a.arms['right'].segments
        deep = max(range(len(segs)), key=lambda i: segs[i].depth)
        self.assertGreater(segs[deep].from_bp, 7000)
        for _ in range(k):
            ev = model.Event(0, 'hairpin')
            ev.char = 'T'
            ev.labels = array('I', [nl - 1])
            ev.total = 1
            ev.followed = False
            segs[deep].events.append(ev)
        a = parse(a.dump())          # the library's own parser takes the body as it is
        b = parse(_wide(nl, 5))
        free = ops.compare(a, b, mode='prefix_subset')
        bud = LocalBudget()
        t = time.process_time()
        got = ops.compare(a, b, mode='prefix_subset', budget=bud)
        cpu = time.process_time() - t
        self.assertEqual((free.only_in_a, free.only_in_b), (got.only_in_a, got.only_in_b))
        self.assertLess(cpu, 5 * bud.used_work * 1e-7 + 0.2, (cpu, bud.used_work))
        # the walk is charged: a work limit below the charge stops the comparison
        stopped = ops.compare(a, b, mode='prefix_subset',
                              budget=LocalBudget(work_units=bud.used_work - 1))
        self.assertEqual('unknown', stopped.comparable)

    def test_claims_before_their_segment_are_not_climbed_per_claim(self):
        # The same for a's claims (_claim_index): a claim whose to_bp lies before its
        # segment was filed by a climb of its chain per claim. 16,000 claims of to_bp 1 on
        # the deepest segment of a 16,000-deep comb: one walk of the tree, charged
        a = parse(comb_annotate(16000))
        arm = a.arms['right']
        segs = arm.segments
        deep = max(range(len(segs)), key=lambda i: segs[i].depth)
        self.assertGreater(segs[deep].depth, 15000)
        base = [c for c in ops._restricted_claims(a, 'right', arm.complete_to_bp)
                if c.kind != 'merged'][0]
        rows = [dataclasses.replace(base, segment=deep, to_bp=1 + (i % 3))
                for i in range(16000)]
        bud = LocalBudget()
        t = time.process_time()
        with bud.scope('probe'):
            index = ops._claim_index(arm, rows, bud)
        cpu = time.process_time() - t
        # each filed under the ancestor that holds its last base, as the climb found it:
        # the first row per (anchor, to_bp)
        _, offs, _ = ops._chain_tree(segs, {deep})
        want, climbed = {}, {}
        for i, c in enumerate(rows):
            x = climbed.get(c.to_bp)
            if x is None:            # three to_bp values: the climb once each
                x = deep
                while offs[x] >= c.to_bp:
                    x = segs[x].parents[0]
                climbed[c.to_bp] = x
            want.setdefault(x, {}).setdefault(c.to_bp, i)
        _, _, by_ref, _ = index
        self.assertEqual(want, {x: {k: e[1] for k, e in zip(ks, entries)}
                                for x, (ks, entries) in by_ref[base.label.ref].items()})
        # the caller charges each row (W_ELEM x 3) and the chain tree; the anchors' walk
        # is charged here, and the time stays within a small factor of the two
        nominal = (bud.used_work + 3 * ops.W_ELEM * len(rows)
                   + 2 * ops.W_STEP * (segs[deep].depth + 1)) * 1e-7
        self.assertLess(cpu, 5 * nominal + 0.2, (cpu, bud.used_work))
        self.assertGreaterEqual(bud.used_work, ops.W_STEP * len(segs) // 2)

    def test_the_same_answer_budgeted_and_a_stop(self):
        a, b = parse(only_label0(200, 400)), parse(_wide(200, 5))
        free = ops.compare(a, b, mode='prefix_subset')
        bud = LocalBudget()
        got = ops.compare(parse(only_label0(200, 400)), b, mode='prefix_subset', budget=bud)
        self.assertEqual((free.only_in_a, free.only_in_b), (got.only_in_a, got.only_in_b))
        stopped = ops.compare(parse(only_label0(200, 400)), b, mode='prefix_subset',
                              budget=LocalBudget(work_units=bud.used_work // 2))
        self.assertEqual('unknown', stopped.comparable)
        self.assertEqual('work', stopped.local_stop['resource'])


if __name__ == '__main__':
    unittest.main()
