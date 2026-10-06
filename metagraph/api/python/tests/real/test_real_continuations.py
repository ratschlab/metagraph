"""continuation() and next_request() on real retrievals, offline and against the servers.

Offline (the retrieval cache): every continuation spells the walk's tail as the
independent speller reads it (right arm: the last n bases of seed + flank; left arm: the
first n bases of flank + seed), equals the server's detail-full continuation, sits at the
stated seed coordinates, has the stated length rule, exists exactly on the leaves the
spec gives one, and next_request() builds the request the design describes.

Live (skips per index when its server is down):
  * the generated request is ACCEPTED by the server (200, a graphlet per seed, no seed
    error, no label dropped), its seed IS the parent walk's tail at the stated overlap,
    on both arms, with continuations inside the flank and crossing into the seed;
  * the child CONTINUES the parent walk: its extensions beyond the overlap are the
    extensions a deeper retrieval of the parent shows beyond that leaf (keep
    strategies, leaves cut by the radius);
  * the continuation's labels against /resolve: constrain -- exactly the leaf's labels
    that cover the tail; annotate -- sound (each covers the tail);
  * TraverseClient.deepen() and the MCP traverse_continue() session tool agree with it.

    cd api/python/tests/real && python3 -m unittest -v test_real_continuations
"""

import collections
import copy
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import realdata as R  # noqa: E402
from test_real_ops import (  # noqa: E402  (helpers only: no TestCase is imported)
    has_bases, live_rows, op_rows, require_rows, sample)

from metagraph.traverse import (  # noqa: E402
    IncompatibleContinuations, MissingEnvelope, TraverseClient, from_response, parse)
from metagraph.traverse._codec import REASON, RESOURCE_CODES  # noqa: E402
from metagraph.traverse.mcp_tools import _size  # noqa: E402


def _env_int(name, default):
    try:
        return int(os.environ.get(name, '') or default)
    except ValueError:
        return default


CONT_CAP = _env_int('METAGRAPH_REAL_CONT_LEAF_CAP', 200)       # per arm, offline checks
LIVE_PER_ARM = _env_int('METAGRAPH_REAL_CONT_LIVE_PER_ARM', 2)  # live requests per arm
DEEP_PER_INDEX = _env_int('METAGRAPH_REAL_CONT_DEEP_PER_INDEX', 8)
BP = _env_int('METAGRAPH_REAL_CONT_BP', 150)                    # the child's radius

CONTINUED_NAMES = frozenset(REASON[c] for c in RESOURCE_CODES) | {'max_extension_bp'}

COVERAGE = collections.defaultdict(set)


def tearDownModule():
    if COVERAGE:
        allc = set().union(*COVERAGE.values())
        sys.stderr.write('\n[test_real_continuations] %d distinct cells covered: %s\n' % (
            len(allc), ', '.join('%s=%d' % (k, len(v)) for k, v in sorted(COVERAGE.items()))))


def results(rows, tag, *, full=True, sequences=None, mode=None):
    for row in rows:
        c = R.load_cell(row['index'], row)
        for i in range(c.n_results()):
            g = c.graphlet(i)
            if g is None:
                continue
            if sequences is not None and has_bases(g) != sequences:
                continue
            if mode is not None and g.mode != mode:
                continue
            COVERAGE[tag].add(c.index + '/' + c.cell)
            yield row, c, i, g, (c.full_result(i) if full else None)


def continued_paths(g, side):
    """Path ids of the leaves that carry a continuation (C record)."""
    a = g.arms[side]
    out = []
    pid = 0
    for s in a.segments:
        if s.leaf is None:
            continue
        if s.leaf.continuation is not None:
            out.append((pid, s.id))
        pid += 1
    return out


def molecule(g, side, leaf_seg):
    """The walk's whole molecule, natural orientation (independent speller), and the
    offset of seed coordinate 0 in it."""
    mol = R.spell_with_seed(g, side, leaf_seg)
    flank = len(mol) - len(g.seed.sequence)
    return mol, (0 if side == 'right' else flank), flank


def tail(g, side, leaf_seg, n):
    mol, _, _ = molecule(g, side, leaf_seg)
    return mol[len(mol) - n:] if side == 'right' else mol[:n]


def lineage_switches(g, side, leaf_seg):
    a = g.arms[side]
    out = []
    s = leaf_seg
    while True:
        out.extend(ev.at_bp for ev in a.segments[s].events if ev.type == 'switch')
        if not a.segments[s].parents:
            return out
        s = a.segments[s].parents[0]


# ================================================================== offline

class TestContinuationSpelling(unittest.TestCase):
    """Every continuation of every cached retrieval (capped per arm), on both arms."""

    @classmethod
    def setUpClass(cls):
        R.require_cache()
        cls.rows = require_rows(op_rows())

    def test_continuations_spell_the_tail_at_the_stated_coordinates(self):
        kinds = collections.Counter()
        for row, c, i, g, full in results(self.rows, 'cont.spelling'):
            bases = has_bases(g)
            seed_len = g.seed.length_bp
            for side in g.arms:
                fa = full['arms'][side]
                with self.subTest(cell=c.cell, result=i, arm=side):
                    for pid, leaf in sample(continued_paths(g, side), CONT_CAP):
                        p = fa['paths'][pid]
                        fc = p['continuation']
                        cont = g.continuation(side, pid)
                        L = p['length_bp']
                        n = len(fc['sequence'])
                        self.assertEqual(fc['sequence'], cont.sequence)
                        self.assertEqual(fc['labels'], [l.id for l in cont.labels])
                        self.assertEqual(fc['loss_used'], cont.loss_used)
                        self.assertEqual(fc['branches_used'], cont.branches_used)
                        self.assertEqual((side, pid), (cont.arm, cont.leaf))
                        if side == 'right':
                            self.assertEqual((seed_len + L - n, seed_len + L), cont.seed_coord)
                        else:
                            self.assertEqual((-L, -L + n), cont.seed_coord)
                        if bases:
                            mol, off, flank = molecule(g, side, leaf)
                            self.assertEqual(L, flank)
                            a_, b_ = cont.seed_coord
                            self.assertEqual(cont.sequence, mol[a_ + off:b_ + off])
                            self.assertEqual(tail(g, side, leaf, n), cont.sequence)
                        seed = cont.as_seed()
                        self.assertEqual(cont.sequence, seed['sequence'])
                        self.assertEqual([l.name for l in cont.labels], seed.get('labels', []))
                        kinds[(side, 'crossing' if n > L else 'inside')] += 1
        self.__class__.kinds = kinds
        self.assertTrue(kinds)

    def test_both_arms_have_continuations_inside_the_flank_and_crossing_into_the_seed(self):
        got = collections.Counter()
        for row, c, i, g, full in results(self.rows, 'cont.kinds', full=False):
            for side in g.arms:
                a = g.arms[side]
                for pid, leaf in continued_paths(g, side):
                    n = a.segments[leaf].leaf.continuation.n
                    got[(side, 'crossing' if n > a.segments[leaf].end_bp else 'inside')] += 1
        missing = [k for k in (('left', 'inside'), ('left', 'crossing'), ('right', 'inside'),
                               ('right', 'crossing')) if not got[k]]
        if missing:
            self.skipTest('the cache holds no continuation of kind %s' % missing)
        sys.stderr.write('\n[continuations] %s\n' % dict(got))

    def test_a_continuation_exactly_where_the_walk_was_cut(self):
        """A leaf carries a continuation iff it ended at the radius, a resource reason or
        the beam (spec §7.1); a semantic end raises."""
        for row, c, i, g, full in results(self.rows, 'cont.presence', full=False):
            for side, a in g.arms.items():
                pid = 0
                with self.subTest(cell=c.cell, result=i, arm=side):
                    for s in a.segments:
                        if s.leaf is None:
                            continue
                        code = s.leaf.path_reason
                        cut = code is not None and (code == 'X' or code in RESOURCE_CODES)
                        self.assertEqual(cut, s.leaf.continuation is not None, (pid, code))
                        if not cut and pid % 17 == 0:
                            with self.assertRaises(ValueError):
                                g.continuation(side, pid)
                        pid += 1

    def test_the_length_rule(self):
        """n = min(continuation_bp, max(k, bp since the tail's labels became evidence for
        the displayed walk)): since the seed boundary, the last switch (spec §7.1) or --
        under on_reconverge merge -- the latest merge the labels entered through a later
        parent (their T route_bp; the displayed bases before it are not theirs). A
        continuation crosses into the seed only to make up k bases."""
        checked = merged = 0
        for row, c, i, g, full in results(self.rows, 'cont.length'):
            cbp = g.continuation_bp
            for side, a in g.arms.items():
                fa = full['arms'][side]
                with self.subTest(cell=c.cell, result=i, arm=side):
                    for pid, leaf in continued_paths(g, side):
                        L = a.segments[leaf].end_bp
                        n = a.segments[leaf].leaf.continuation.n
                        labels = set(fa['paths'][pid]['continuation']['labels'])
                        route = max((e['route_bp'] for e in fa['paths'][pid]['end_labels']
                                     if e['label'] in labels), default=0)
                        merged += route > 0
                        sw = lineage_switches(g, side, leaf)
                        if not sw:
                            self.assertEqual(min(cbp, max(g.k, L - route)), n, (pid, L, route))
                        else:
                            self.assertLessEqual(n, min(cbp, max(g.k, L)), (pid, L, sw))
                            self.assertGreaterEqual(n, min(cbp, g.k), (pid, L, sw))
                        if n > L:
                            self.assertEqual(g.k, n, 'crossing into the seed beyond k')
                        checked += 1
        self.assertGreater(checked, 0)
        sys.stderr.write('\n[continuation length] %d checked, %d restarted at a merge\n'
                         % (checked, merged))

    def test_without_bases_the_continuations_equal_those_with_bases(self):
        """output.sequences false: C carries the sequence; the same strategy with bases
        (limit1) gives the same continuations walk by walk."""
        n = 0
        for row in op_rows(strategy='no_sequences'):
            other = R.cell_id(row['seed'], 'limit1')
            meta = R.cell_meta(row['index'], other)
            if meta is None or meta.get('status') != 'ok':
                continue
            a_cell = R.load_cell(row['index'], row)
            b_cell = R.load_cell(row['index'], other)
            COVERAGE['cont.no_sequences'].update({row['index'] + '/' + row['cell'],
                                                  row['index'] + '/' + other})
            for i in range(a_cell.n_results()):
                ga, gb = a_cell.graphlet(i), b_cell.graphlet(i)
                if ga is None:
                    continue
                self.assertFalse(has_bases(ga))
                for side in ga.arms:
                    pa, pb = continued_paths(ga, side), continued_paths(gb, side)
                    self.assertEqual([x[0] for x in pa], [x[0] for x in pb])
                    for (pid, _), _ in zip(pa, pb):
                        ca, cb = ga.continuation(side, pid), gb.continuation(side, pid)
                        with self.subTest(cell=row['cell'], arm=side, walk=pid):
                            self.assertEqual(cb.sequence, ca.sequence)
                            self.assertEqual([l.ref for l in cb.labels],
                                             [l.ref for l in ca.labels])
                            self.assertEqual(cb.seed_coord, ca.seed_coord)
                            self.assertEqual(_without_sequences_flag(gb.next_request(side, [pid])),
                                             _without_sequences_flag(ga.next_request(side, [pid])))
                        n += 1
        if n == 0:
            self.skipTest('no no_sequences / limit1 pair with continuations cached')


def _without_sequences_flag(req):
    r = copy.deepcopy(req)
    r['strategy'].get('output', {}).pop('sequences', None)
    return r


def _request_violations(req):
    """traverse_testlib.request_violations (the server's label-list checks restated)."""
    tests = os.path.dirname(HERE)
    if tests not in sys.path:
        sys.path.insert(0, tests)
    import traverse_testlib
    return traverse_testlib.request_violations(req)


def _largest_loss(cs):
    """The largest terminal loss of the continued labels (third review, finding 4: C's
    loss_used is the SMALLEST, and subtracting it let a label at a higher loss go on with
    more than its original budget left)."""
    return max((x for c in cs for x in c.losses), default=0.0)


class TestNextRequestShape(unittest.TestCase):
    """next_request() offline: one seed per walk (its continuation; constrain names its
    labels, annotate has no permitted set), the retrieval's normalized strategy with the
    arm as direction, bp as the radius, the loss budget reduced by the largest terminal
    loss of the continued labels, labels.extra rebuilt around the new seeds (valid for the
    server), overrides deep-merged, release kept; walks whose rebuilt pools differ are
    refused as one request (IncompatibleContinuations) and built one per walk; a body
    without its envelope refuses."""

    @classmethod
    def setUpClass(cls):
        R.require_cache()
        cls.rows = require_rows(op_rows())

    def test_request_fields(self):
        n = reduced = split = 0
        for row, c, i, g, full in results(self.rows, 'next_request.shape', full=False):
            base = copy.deepcopy(g.envelope['strategy'])
            base.pop('clamped', None)
            for side in g.arms:
                conts = continued_paths(g, side)
                if not conts:
                    continue
                pids = [x[0] for x in sample(conts, 3)]
                with self.subTest(cell=c.cell, result=i, arm=side):
                    try:
                        groups = [(pids, g.next_request(side, pids, bp=BP))]
                    except IncompatibleContinuations:
                        # the walks continue under different labels: one request each,
                        # each valid and with its own labels' budget
                        self.assertEqual('constrain', g.mode)
                        groups = [([p], r) for p, r in
                                  zip(pids, g.next_requests(side, pids, bp=BP))]
                        split += 1
                    for ids, req in groups:
                        self.assertEqual(len(ids), len(req['seeds']))
                        cs = [g.continuation(side, p) for p in ids]
                        for seed, cont in zip(req['seeds'], cs):
                            want = cont.as_seed()
                            if g.mode != 'constrain':
                                want.pop('labels', None)
                            self.assertEqual(want, seed)
                        self.assertEqual([], _request_violations(req))
                        st = req['strategy']
                        want = copy.deepcopy(base)
                        want['direction'] = side
                        want.setdefault('bounds', {})['max_extension_bp'] = BP
                        budget = (base.get('labels') or {}).get('loss_budget')
                        used = _largest_loss(cs)
                        if g.mode == 'constrain' and isinstance(budget, (int, float)) \
                                and budget > 0 and used > 0:
                            want['labels']['loss_budget'] = max(0.0, budget - used)
                            reduced += 1
                            try:
                                keeps = [g.next_request(side, ids, reduce_budget=False)]
                            except IncompatibleContinuations:
                                # with the whole budget more targets are reachable, and
                                # they differ per walk
                                keeps = g.next_requests(side, ids, reduce_budget=False)
                            for keep in keeps:
                                self.assertEqual(budget,
                                                 keep['strategy']['labels']['loss_budget'])
                                self.assertEqual([], _request_violations(keep))
                        if g.mode == 'constrain':
                            # rebuilt around the new seeds: the retrieval's labels minus
                            # them, in dictionary order (each kept one valid, above)
                            seeds = {x for sd in req['seeds'] for x in sd.get('labels', ())}
                            pool = [l.name for l in g.labels if l.name not in seeds]
                            got = st['labels']['extra']
                            self.assertEqual(got, [x for x in pool if x in got])
                            want['labels']['extra'] = got
                        if req.branch_budget is not None:
                            # the branch allowance reduced by the largest terminal branch
                            # count (v5.4, DESIGN §19.3): its derivation is checked in
                            # branch_budget, its value is what the request carries
                            want.setdefault('branching', {})['max_label_branches'] = \
                                req.branch_budget['effective']
                        # everything else is the retrieval's own normalized strategy
                        self.assertEqual(want, st)
                    if g.envelope.get('release'):
                        self.assertEqual(g.envelope['release'], req['release'])
                    # bp None keeps the radius; overrides deep-merge; graph is top level
                    same = g.next_request(side, pids[:1], graph='other')
                    self.assertEqual((base.get('bounds') or {}).get('max_extension_bp'),
                                     (same['strategy'].get('bounds') or {})
                                     .get('max_extension_bp'))
                    self.assertEqual('other', same['graph'])
                    ov = g.next_request(side, pids[:1], bp=BP,
                                        bounds={'max_live_paths': 7}, output={'detail': 'full'})
                    self.assertEqual(7, ov['strategy']['bounds']['max_live_paths'])
                    self.assertEqual(BP, ov['strategy']['bounds']['max_extension_bp'])
                    self.assertEqual('full', ov['strategy']['output']['detail'])
                    n += 1
        self.assertGreater(n, 0)
        sys.stderr.write('\n[next_request] %d arms, %d with a reduced loss budget, %d built '
                         'one request per walk\n' % (n, reduced, split))

    def test_a_body_alone_and_a_semantic_leaf_refuse(self):
        for row, c, i, g, full in results(sample(self.rows, 30), 'next_request.refuse',
                                          full=False):
            body = parse(c.text(i))
            for side in g.arms:
                conts = continued_paths(g, side)
                if conts:
                    with self.assertRaises(MissingEnvelope):
                        body.next_request(side, [conts[0][0]])
                ended = [p for p in range(g.arms[side].counts.leaves)
                         if p not in {x[0] for x in conts}]
                if ended:
                    with self.assertRaises(ValueError):
                        g.next_request(side, [ended[0]])


# ================================================================== live

def _client(index):
    return R.client(index)


def _pick_leaves(g, side, k):
    """Up to k continued walks of an arm: one crossing into the seed and one inside
    the flank first, when there are."""
    a = g.arms[side]
    conts = continued_paths(g, side)
    cross = [x for x in conts if a.segments[x[1]].leaf.continuation.n > a.segments[x[1]].end_bp]
    inside = [x for x in conts if x not in cross]
    out = cross[:1] + inside[:1]
    for x in sample(conts, k):
        if len(out) >= k:
            break
        if x not in out:
            out.append(x)
    return out[:k]


LIVE_STRATEGIES = ('strict', 'limit1', 'limit2_keep', 'exhaustive', 'switch1', 'trace',
                   'annotate_exh', 'beam1', 'cut_lists', 'no_sequences', 'cut_events')


class TestNextRequestLive(unittest.TestCase):
    """The generated request against the server that produced the parent."""

    @classmethod
    def setUpClass(cls):
        R.require_cache()
        cls.rows = live_rows(require_rows(op_rows(tier='core', strategy=LIVE_STRATEGIES)))

    def test_requests_are_accepted_and_overlap_the_parent_walk(self):
        kinds = collections.Counter()
        for row, c, i, g, full in results(self.rows, 'live.accepted', full=False):
            cl = _client(c.index)
            for side in g.arms:
                for pid, leaf in _pick_leaves(g, side, LIVE_PER_ARM):
                    cont = g.continuation(side, pid)
                    req = g.next_request(side, [pid], bp=BP)
                    with self.subTest(cell=c.cell, result=i, arm=side, walk=pid):
                        out = cl.traverse_raw(req)
                        self.assertEqual(1, len(out['results']))
                        res = out['results'][0]
                        self.assertNotIn('error', res, res.get('error'))
                        self.assertEqual(side, out['strategy']['direction'])
                        self.assertEqual(BP, out['strategy']['bounds']['max_extension_bp'])
                        child = from_response(res, out)
                        # the child's seed is the parent's tail, exactly as stated
                        self.assertEqual(cont.sequence, child.seed.sequence)
                        self.assertEqual(len(cont.sequence), child.seed.length_bp)
                        if has_bases(g):
                            mol, off, flank = molecule(g, side, leaf)
                            a_, b_ = cont.seed_coord
                            self.assertEqual(child.seed.sequence, mol[a_ + off:b_ + off])
                            if side == 'right':
                                self.assertTrue(mol.endswith(child.seed.sequence))
                            else:
                                self.assertTrue(mol.startswith(child.seed.sequence))
                        self.assertEqual([side], list(child.arms))
                        self.assertEqual(g.mode, child.mode)
                        self.assertEqual(g.support, child.support)
                        self.assertEqual((g.index_ns, g.index_meta_fp, g.index_fp),
                                         (child.index_ns, child.index_meta_fp, child.index_fp))
                        self.assertEqual([], child.dropped)
                        self.assertEqual(0, res['seed'].get('labels_dropped', 0))
                        if g.mode == 'constrain':
                            self.assertEqual(sorted(l.name for l in cont.labels),
                                             sorted(l.name for l in
                                                    child.labels[:child.seed.num_seed_labels]))
                            self.assertFalse(res['seed'].get('labels_from_seed'))
                        # every child walk extends the parent's molecule on the arm's side
                        if not has_bases(child):
                            # output.sequences false travels with the strategy
                            self.assertFalse(has_bases(g))
                        for w in sample(child.walks(side, by='id'), 5):
                            if w.sequence is None:
                                continue
                            cm = R.spell_with_seed(child, side, w.leaf)
                            if has_bases(g):
                                joined = mol + w.sequence if side == 'right' else \
                                    w.sequence + mol
                                n = len(cont.sequence)
                                self.assertEqual(cm, joined[len(joined) - len(cm):]
                                                 if side == 'right' else joined[:len(cm)])
                                self.assertEqual(n + w.length_bp, len(cm))
                        cross = len(cont.sequence) > flank_len(g, side, leaf)
                        kinds[(c.index, side, 'crossing' if cross else 'inside')] += 1
        self.assertTrue(kinds)
        sys.stderr.write('\n[next_request live] %s\n' % dict(kinds))
        for side in ('left', 'right'):
            for kind in ('inside', 'crossing'):
                if not any(k[1] == side and k[2] == kind for k in kinds):
                    sys.stderr.write('[next_request live] no %s continuation on the %s arm '
                                     'was submitted\n' % (kind, side))

    def test_a_reduced_loss_budget_is_accepted_and_echoed(self):
        """A walk that spent part of the loss budget (a switch) continues with the rest:
        the request carries budget minus the largest terminal loss of its labels, is
        valid (labels.extra rebuilt around the new seeds), the server echoes it, and the
        child's seed labels are the continuation's (the switched-to label included)."""
        n = 0
        for row, c, i, g, full in results([r for r in self.rows if r['strategy'] == 'switch1'],
                                          'live.budget', full=False):
            budget = (g.envelope['strategy'].get('labels') or {}).get('loss_budget') or 0
            cl = _client(c.index)
            for side in g.arms:
                for pid, leaf in continued_paths(g, side):
                    cont = g.continuation(side, pid)
                    if _largest_loss([cont]) <= 0:
                        continue
                    req = g.next_request(side, [pid], bp=BP)
                    with self.subTest(cell=c.cell, arm=side, walk=pid):
                        self.assertEqual([], _request_violations(req))
                        want = max(0.0, budget - _largest_loss([cont]))
                        self.assertEqual(want, req['strategy']['labels']['loss_budget'])
                        out = cl.traverse_raw(req)
                        res = out['results'][0]
                        self.assertNotIn('error', res, res.get('error'))
                        self.assertEqual(want, out['strategy']['labels']['loss_budget'])
                        child = from_response(res, out)
                        self.assertEqual(sorted(l.name for l in cont.labels),
                                         sorted(l.name for l in
                                                child.labels[:child.seed.num_seed_labels]))
                        self.assertEqual([], child.dropped)
                    n += 1
        if n == 0:
            self.skipTest('no continued walk that spent loss budget')
        sys.stderr.write('\n[reduced budget] %d continuations submitted\n' % n)

    def test_several_walks_in_one_request(self):
        n = 0
        for row, c, i, g, full in results(self.rows, 'live.batch', full=False):
            cl = _client(c.index)
            for side in g.arms:
                picks = _pick_leaves(g, side, 3)
                if len(picks) < 2:
                    continue
                pids = [x[0] for x in picks]
                try:
                    groups = [(pids, g.next_request(side, pids, bp=BP))]
                except IncompatibleContinuations:
                    # different labels need different labels.extra: one request each
                    groups = [([p], r) for p, r in
                              zip(pids, g.next_requests(side, pids, bp=BP))]
                with self.subTest(cell=c.cell, arm=side):
                    for ids, req in groups:
                        self.assertEqual([], _request_violations(req))
                        resp = cl.response(cl.traverse_raw(req))
                        self.assertEqual([], resp.errors)
                        self.assertEqual(len(ids), len(resp.graphlets))
                        for pid, child in zip(ids, resp.graphlets):
                            self.assertEqual(g.continuation(side, pid).sequence,
                                             child.seed.sequence)
                n += 1
                break
            if n >= 12:
                break
        if n == 0:
            self.skipTest('no arm with two continued walks')

    def test_deepen_is_next_request_plus_traverse(self):
        n = 0
        rows = [r for r in self.rows if r['index'] == 'mini_refseq'
                and r['strategy'] in ('strict', 'annotate_exh')]
        for row, c, i, g, full in results(rows, 'live.deepen', full=False):
            cl = _client(c.index)
            for side in g.arms:
                for pid, leaf in _pick_leaves(g, side, 1):
                    with self.subTest(cell=c.cell, arm=side, walk=pid):
                        via = cl.deepen(g, side, [pid], BP)
                        req = g.next_request(side, [pid], BP)
                        direct = cl.response(cl.traverse_raw(req))
                        self.assertEqual(direct.graphlets[0].dump(envelope=False),
                                         via.graphlets[0].dump(envelope=False))
                    n += 1
        if n == 0:
            self.skipTest('no mini_refseq strict / annotate_exh cell with continuations')


REPEAT_ENDS = frozenset({'edge_reuse', 'edge_reuse_rc', 'rejoined_seed'})


def _consistent(x, y):
    return x.startswith(y) or y.startswith(x)


def flank_len(g, side, leaf):
    return g.arms[side].segments[leaf].end_bp


class TestChildContinuesTheParentWalk(unittest.TestCase):
    """The child continues the parent walk: beyond the overlap, the child's walks are
    the walks a deeper retrieval of the parent (the same request, this arm only, the
    radius raised by BP) takes beyond that leaf. Keep strategies (no history is united,
    so the deeper trie displays the leaf's own continuation), complete arms, leaves cut
    by the radius. The child is a fresh seed (branch and edge-reuse state reset), so it
    may continue MORE than the deeper parent; on the exhaustive tries the two are equal,
    labels included."""

    @classmethod
    def setUpClass(cls):
        R.require_cache()
        cls.rows = live_rows(require_rows(op_rows(
            tier='core', strategy=('exhaustive', 'limit2_keep', 'trace', 'annotate_exh'))))

    def test_child_extensions_are_the_deeper_parents(self):
        per_index = collections.Counter()
        equal = repeats = 0
        for row in self.rows:
            if per_index[row['index']] >= DEEP_PER_INDEX:
                continue
            if ((row.get('graphlet') or {}).get('elapsed_s') or 0) > 2.0:
                continue
            for _, c, i, g, full in results([row], 'live.deeper', full=False, sequences=True):
                cl = _client(c.index)
                for side, a in g.arms.items():
                    if a.status != 'complete':
                        continue
                    radius = [x for x in continued_paths(g, side)
                              if a.segments[x[1]].leaf.path_reason == 'X']
                    if not radius:
                        continue
                    deep_req = copy.deepcopy(c.request)
                    st = deep_req['strategy']
                    st['direction'] = side
                    r0 = (st.get('bounds') or {}).get('max_extension_bp', 5000)
                    st.setdefault('bounds', {})['max_extension_bp'] = r0 + BP
                    deep_out = cl.traverse_raw(deep_req)
                    deep = from_response(deep_out['results'][0], deep_out)
                    if deep.arms[side].status != 'complete':
                        continue
                    deep_walks = [(R.spell_with_seed(deep, side, w.leaf), w)
                                  for w in deep.walks(side, by='id')]
                    for pid, leaf in radius[:2]:
                        pm = R.spell_with_seed(g, side, leaf)
                        cont = g.continuation(side, pid)
                        n = len(cont.sequence)
                        out = cl.traverse_raw(g.next_request(side, [pid], bp=BP))
                        child = from_response(out['results'][0], out)
                        if child.arms[side].status != 'complete':
                            continue
                        # extensions read OUTWARD (walking order) on both arms
                        child_ext, repeat = {}, False
                        for w in child.walks(side, by='id'):
                            cm = R.spell_with_seed(child, side, w.leaf)
                            ext = cm[n:] if side == 'right' else cm[:len(cm) - n][::-1]
                            child_ext.setdefault(ext, set()).update(l.name for l in w.labels_full)
                            repeat |= bool(REPEAT_ENDS & set(w.end_reasons))
                        deep_ext = {}
                        for m, w in deep_walks:
                            if side == 'right' and m.startswith(pm):
                                ext = m[len(pm):]
                            elif side == 'left' and m.endswith(pm):
                                ext = m[:len(m) - len(pm)][::-1]
                            else:
                                continue
                            deep_ext.setdefault(ext, set()).update(l.name for l in w.labels_full)
                            repeat |= bool(REPEAT_ENDS & set(w.end_reasons))
                        # a repeat is handled relative to each retrieval's own seed and
                        # edge history (edge reuse, the seed rejoined, a blocked successor):
                        # there the two may stop at different places, legitimately
                        repeat |= any(ev.type in ('blocked', 'hairpin')
                                      for x in child.arms[side].segments for ev in x.events)
                        with self.subTest(cell=c.cell, arm=side, walk=pid):
                            self.assertTrue(deep_ext, 'the deeper parent lost the leaf')
                            # the first bases agree (unless the child's very first step
                            # was blocked by the repeat)
                            first = {x[:1] for x in deep_ext if x}
                            cfirst = {y[:1] for y in child_ext if y}
                            if cfirst or not repeat:
                                self.assertEqual(first, cfirst)
                            if repeat:
                                repeats += 1
                                continue
                            self.assertLessEqual(set(deep_ext), set(child_ext))
                            if row['strategy'] == 'exhaustive':
                                self.assertEqual(deep_ext, child_ext)
                                equal += 1
                        per_index[c.index] += 1
        if not per_index:
            self.skipTest('no complete keep retrieval with a radius-cut leaf')
        sys.stderr.write('\n[deeper parent] %s, %d exhaustive leaves equal with labels, %d '
                         'repeat-affected compared on their first base only\n'
                         % (dict(per_index), equal, repeats))
        self.assertGreater(equal, 0)


class TestContinuationLabelsAgainstTheOracle(unittest.TestCase):
    """The labels of a continuation on real k-mers (/resolve on the tail sequence)."""

    @classmethod
    def setUpClass(cls):
        R.require_cache()
        cls.rows = live_rows(require_rows(op_rows(tier='core')))

    def _tails(self, mode):
        for row, c, i, g, full in results(self.rows, 'oracle.cont_labels.' + mode,
                                          full=False, mode=mode):
            if len(g.labels) > 3000:
                continue
            for side in g.arms:
                for pid, leaf in sample(continued_paths(g, side), 3):
                    cont = g.continuation(side, pid)
                    if len(cont.sequence) < g.k:
                        continue
                    yield c, i, g, side, pid, leaf, cont

    def test_constrain_continuation_labels_are_the_leaf_labels_covering_the_tail(self):
        n = 0
        from metagraph.traverse import derive
        for c, i, g, side, pid, leaf, cont in self._tails('constrain'):
            a = g.arms[side]
            alive = [g.labels[e.label] for e in derive.end_labels(a, leaf)]
            if not alive:
                continue
            k = R.oracle_num_kmers(c.index, cont.sequence)
            runs = R.oracle_label_runs(c.index, cont.sequence, alive)
            cover = sorted(l.name for l in alive if runs[l.name] == [(0, k)])
            with self.subTest(cell=c.cell, result=i, arm=side, walk=pid):
                self.assertEqual(cover, sorted(l.name for l in cont.labels))
            n += 1
        self.assertGreater(n, 0)

    def test_annotate_continuation_labels_cover_the_tail(self):
        """Sound: every label the continuation names covers every k-mer of its tail.
        And it is the walker's set: the labels recorded on every node whose k-mer lies in
        the tail -- the last n - k + 1 steps of the walk (its first k - 1 bases are no
        node's newest base inside the tail), from the seed boundary when they reach past
        it (seed nodes carry no recorded labels)."""
        n = 0
        for c, i, g, side, pid, leaf, cont in self._tails('annotate'):
            L = g.arms[side].segments[leaf].end_bp
            lo = max(0, L - (len(cont.sequence) - g.k + 1))
            have = None
            for r in g.support_profile(side, pid):
                if r.to_bp <= lo:
                    continue
                ids = {l.id for l in r.labels}
                have = ids if have is None else have & ids
            with self.subTest(cell=c.cell, result=i, arm=side, walk=pid):
                if have is not None and g.arms[side].labels_per_node.nodes_truncated == 0:
                    self.assertEqual(sorted(have), [l.id for l in cont.labels])
                if cont.labels:
                    k = R.oracle_num_kmers(c.index, cont.sequence)
                    runs = R.oracle_label_runs(c.index, cont.sequence, cont.labels)
                    for l in cont.labels:
                        self.assertEqual([(0, k)], runs[l.name], l.name)
            n += 1
        self.assertGreater(n, 0)

    def test_annotate_continuation_names_every_label_covering_the_tail(self):
        """Spec §7.1: a continuation names 'the labels covering that whole tail'.
        Fixed SRV-ANNOT-CONT-LABELS (walker.cpp, make_continuation): in annotate mode the
        walker intersected the labels of the tail's last n STEPS, whose first k-1 nodes
        begin before the tail, so a label carried by every k-mer of the tail but not by
        those k-1 straddling nodes was left out. Old repro: sra_16s_PZ326290__annotate_exh,
        left arm, continuation(left, 6): 45 bp, labels [] -- /resolve on that 45 bp
        sequence reports .../SRR14483890.fasta.gz on all 15 k-mers (it is recorded on
        steps 23-44 of the walk only). Needs the cache filled from the fixed server."""
        missing = []
        for c, i, g, side, pid, leaf, cont in self._tails('annotate'):
            if g.arms[side].labels_per_node.nodes_truncated:
                continue
            k = R.oracle_num_kmers(c.index, cont.sequence)
            runs = R.oracle_label_runs(c.index, cont.sequence, g.labels)
            cover = {name for name, r in runs.items() if r == [(0, k)]}
            got = {l.name for l in cont.labels}
            if cover - got:
                missing.append((c.cell, side, pid, len(cover - got)))
        self.assertEqual([], missing[:10], '%d annotate continuations understate their '
                         'labels' % len(missing))


class TestTraverseContinueSession(unittest.TestCase):
    """The MCP session tool: traverse_fetch -> traverse_continue keeps the parent walk,
    states the overlap and stores the child under a new handle."""

    def _sessions(self):
        R.require_server('mini_refseq')
        from metagraph.traverse import GraphletStore
        from metagraph.traverse.mcp_tools import GraphletTools
        seeds = [s for s in R.catalog_seeds('mini_refseq', tier='core') if s['kind'] != 'batch']
        if not seeds:
            self.skipTest('no mini_refseq catalog')
        hp = R.server('mini_refseq')
        with tempfile.TemporaryDirectory() as spool:
            tools = GraphletTools(GraphletStore(spool),
                                  {'mini_refseq': TraverseClient(hp[0], hp[1])})
            for seed in seeds[:3]:
                for strat in ('strict', 'annotate_exh'):
                    got = tools.traverse_fetch('mini_refseq', {'sequence': seed['sequence']},
                                               R.strategy(strat, 'mini_refseq', seed))
                    self.assertNotIn('error', got, got)
                    handle = got['handle']
                    g = tools.store.graphlet(handle)
                    for side in g.arms:
                        for pid, leaf in _pick_leaves(g, side, 1):
                            yield tools, handle, g, seed, strat, side, pid

    def test_traverse_continue_states_the_overlap(self):
        n = 0
        for tools, handle, g, seed, strat, side, pid in self._sessions():
            cont = g.continuation(side, pid)
            ov = {'bounds': {'max_extension_bp': BP}}
            res = tools.traverse_continue(handle, side, pid, ov)
            with self.subTest(seed=seed['name'], strategy=strat, arm=side):
                self.assertNotIn('error', res, res)
                self.assertEqual({'handle': handle, 'arm': side, 'walk': pid,
                                  'overlap_bp': len(cont.sequence)}, res['parent'])
                self.assertNotEqual(handle, res['handle'])
                child = tools.store.graphlet(res['handle'])
                self.assertEqual(cont.sequence, child.seed.sequence)
                self.assertEqual([side], list(child.arms))
                self.assertIn('evidence', res)
                self.assertEqual(BP, child.envelope['strategy']['bounds']['max_extension_bp'])
            n += 1
        if n == 0:
            self.skipTest('no continued walk on the mini_refseq seeds')

    def test_traverse_continue_dry_run_returns_the_request(self):
        """execute=False returns the request (design §6) under the 16 KB sequence ceiling
        it states. Fixed CONTINUE-DRYRUN-TOO-LARGE (mcp_tools.traverse_continue): the
        request of a default continuation (1000 bp of sequence, its label names and the
        full normalized strategy, ~2.1 KB) was held to the 2048-byte list-tool ceiling
        and came back as result_too_large, with no max_bytes to raise it. Old repro:
        traverse_fetch('mini_refseq', {sequence: mini_win200_00}, strict) then
        traverse_continue(handle, 'left', 0, {'bounds': {'max_extension_bp': 150}},
        execute=False) -> {'error': 'result_too_large', 'bytes': 2130, 'max_bytes': 2048}."""
        bad = []
        n = 0
        for tools, handle, g, seed, strat, side, pid in self._sessions():
            ov = {'bounds': {'max_extension_bp': BP}}
            dry = tools.traverse_continue(handle, side, pid, ov, False)
            if dry.get('request') != g.next_request(side, [pid], **copy.deepcopy(ov)):
                bad.append((seed['name'], strat, side, pid, dry.get('error'), dry.get('bytes')))
            else:
                self.assertEqual(tools.sequence_max_bytes, dry['ceiling_bytes'])
                self.assertLessEqual(_size(dry), dry['ceiling_bytes'])
            n += 1
        if n == 0:
            self.skipTest('no continued walk on the mini_refseq seeds')
        self.assertEqual([], bad)


if __name__ == '__main__':
    unittest.main()
