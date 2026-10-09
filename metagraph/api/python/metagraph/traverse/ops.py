"""The local queries over a graphlet (DESIGN-traverse-graphlet.md §5, §5.1).

Everything here is combinatorics on the parsed structure: nothing reads the graph.
Labels are resolved by tagged selectors and every result carries Label objects
({name, ref}), never a bare name. Walk arguments accept a path id (paths[].id, the
"walk" of the MCP tools), a Path or a Walk.

Evidence semantics (§5, §5.1):
  * displayed-path support of a label on a walk = the labels alive at each base of the
    walk's first-parent chain (entry sets, run ends, switch-ins), so a lineage merged
    in from a non-first parent supports the displayed bases only from the merge on --
    the same as the per-claim evidence_from of §5.1 (max of the merge-derived route
    start, evaluated at the run's ANCHORED endpoint, and the run's own start);
  * route support = the run's own [from_bp, to_bp), whatever route it took;
  * a claim cut at D: the run-start guard first (no claim for a run starting at or
    after D), then [evidence_from, to_bp) ∩ [0, D), route_only when that is empty,
    zero-length runs are boundary claims; terminal loss/branches are unknown at a cut;
  * annotate mode has no runs: a claim is a maximal end of a label's routes under the
    §5.1 union rule of direct_bp (through any parent of a merge), its displayed support
    starting at the last merge its route enters through a non-first parent -- route and
    displayed evidence kept apart exactly as for a run. Each end is shown by ONE witness
    route, chosen deterministically by the stored parent order (at every merge the first
    parent that records the label): it does not enumerate the label's other routes to
    that end, and -- like every annotate claim, which rests on per-node recorded sets --
    it does not establish a contiguous source occurrence (one indexed sequence spelling
    the route).

Resources: without a budget (the default) everything here runs locally with no work or
allocation budget -- compare(), routes() and the exports included. Their time and peak
allocation follow the graphlet's size (a comparison reads both DAGs up to the comparison
depth, routes() walks one route per run, an export spells every chosen walk), and
nothing interrupts them; an MCP tool's max_bytes bounds the bytes it RETURNS, not the
computation behind them. Every public operation takes budget= (local limits, budget.py:
work units at the cold price and a modelled memory account): it then completes or stops and
says so -- a list whose order allows it with its whole rows and a resume token, compare()
with comparable 'unknown', everything else with no partial answer.
"""

import bisect
import collections
import copy
import hashlib
import heapq
import itertools
import json
import math
import sys

from . import budget as _B
from . import coords as _C
from . import derive
from ._codec import REASON, RESOURCE_CODES, UNLIMITED, GraphletFormatError
from .budget import (
    BASE_SHIFT, DICT, DICT_KEY, ID_SHIFT_C, INT, LIST, LIST_ITEM, SET, SET_ITEM, STR, TUPLE,
    W_ELEM, W_ROW, W_STEP, LocalBudgetExceeded, Partial, dict_bytes, list_bytes, record_bytes,
    set_bytes, sort_work, str_bytes,
)
from .model import (
    AmbiguousLabel, ARM_SIDES, BadSelector, Branch, Change, Claim, Comparison,
    Continuation, CoordClaim, CoordLabelWalk, CoordWalk, IncompatibleContinuations,
    IncompleteRecording, Label, LabelWalk, NextRequest, Path, SplitPoint, SupportRun,
    UnknownLabel, UnverifiableLabelName, Walk,
)


# ------------------------------------------------------------------ local-limits plumbing
# (budget.py has the contract; these are the operations' charge points)

def _graphlet_key(g):
    """What a resume token binds besides the call's arguments: the graphlet's text, by
    the digest of its records (Graphlet.body_digest, made by the parse). A fingerprint of
    the seed and the counts can be shared by graphlets with different bodies (five pairs of
    the fixtures share one), so a token of one would be accepted on the other and resume a
    list it does not belong to; the digest is the same for the same body parsed again (a store's
    re-parse) or saved and loaded."""
    d = g.body_digest
    if d is None:
        # a model the parser did not make: its canonical text, once (not charged: no
        # model of the library is made this way)
        from .parser import body_digest, dump
        with _B.unbudgeted():
            d = g.body_digest = body_digest(dump(g, envelope=False))
    return d


def _token_digest(g, op, args):
    return hashlib.sha256(json.dumps([op, _graphlet_key(g), args], sort_keys=True,
                                     default=str).encode()).hexdigest()[:16]


def _resume_pos(g, op, args, resume, n):
    """The position a resume token (Partial.resume of a stopped call) carries: |n| ints;
    zeros without one. A token of another operation, other arguments or another graphlet
    is a ValueError (a bad argument), never a silently shifted list."""
    if resume is None:
        return (0,) * n
    if not (isinstance(resume, (tuple, list)) and len(resume) == 2 + n and resume[0] == op
            and all(isinstance(x, int) and not isinstance(x, bool) and x >= 0
                    for x in resume[2:])):
        raise ValueError('resume= takes the Partial.resume of a stopped %s() call, not %r'
                         % (op, resume))
    if resume[1] != _token_digest(g, op, args):
        raise ValueError('this resume token belongs to a %s() call with other arguments or '
                         'on another graphlet: repeat that call\'s arguments with it' % op)
    return tuple(resume[2:])


def _token(g, op, args, *pos):
    return (op, _token_digest(g, op, args)) + tuple(pos)


def _uses_g(b, g, *names):
    """Charge the graphlet-level derivations a call uses (once per call, cold price)."""
    if b is None:
        return
    n = len(g.labels)
    for name in names:
        if name == 'label_index':
            # two dicts at CPython's fill (up to 48 bytes a key after a resize), a 1-tuple
            # per key and the ref string Label.ref makes for each (measured: 200-260 bytes
            # a label; dict_bytes() models a dict at 16 bytes a key and would leave the index
            # at 0.7x its size)
            b.uses((id(g), name), 4 * n,
                   2 * (DICT + 48 * n) + n * (2 * (TUPLE + 8) + STR + 24))
        elif name == 'name_counts':
            b.uses((id(g), name), 2 * n, dict_bytes(n) + n * INT)
        elif name == 'label_summary':
            w = 6 * n * len(g.arms)
            m = len(g.arms) * n * (dict_bytes(4) + LIST) + n * (LIST + LIST_ITEM)
            for a in g.arms.values():
                derive.uses(b, g, a)
                z = derive.arm_sizes(a)
                if g.mode == 'constrain':
                    w += 6 * z['runs'] + 2 * z['switches'] + 2 * z['segments']
                    m += z['runs'] * (LIST_ITEM + 8)
                else:
                    w += W_ELEM * (3 * z['segments'] + z['presence']) \
                        + ((2 * z['presence_ids'] + 2 * z['entry_ids']) >> 2)
                    m += z['segments'] * (SET + LIST_ITEM) + z['entry_ids'] * SET_ITEM
            b.uses((id(g), name), w, m)
        else:
            raise KeyError(name)


def _label_ids_b(g, selectors, b):
    """_label_ids() under a budget: each selector resolved at a row's price, the name
    index charged once when a selector names a label by its name or ref."""
    if b is not None and selectors is not None:
        sels = [selectors] if isinstance(selectors, (str, dict, int, Label)) else selectors
        if any(isinstance(s, (str, dict, Label)) for s in sels):
            _uses_g(b, g, 'label_index')
        b.charge(W_ELEM * (1 + len(sels)), set_bytes(len(sels)))
    return _label_ids(g, selectors)


_BLOCK = 64


def _pay_rows(b, price, i, n, keep, m_row):
    """Admit rows i.. of a list (row i is one it makes) in a block of up to _BLOCK at
    their known prices -- price[t] lwu and m_row bytes each (a list: per row), plus the
    caller's row_extra; |keep|: which of them the list makes, None for all -> the index
    below which rows are paid. A block that does not fit is paid row by row, so the list
    stops at exactly the row it would stop at without blocks, having charged exactly the
    same: one charge per block instead of one per row (a charge per row of a cheap loop
    made routes() 409 % slower)."""
    j = min(n, i + _BLOCK)
    x = b.row_extra
    per_row = isinstance(m_row, list)
    if keep is None:
        k = j - i
        w = sum(price[i:j])
        m = sum(m_row[i:j]) if per_row else k * m_row
    else:
        sel = [t for t in range(i, j) if keep(t)]
        k = len(sel)
        w = sum([price[t] for t in sel])
        m = sum([m_row[t] for t in sel]) if per_row else k * m_row
    if b.admit(w + k * x[0], m + k * x[1]):
        return j
    b.charge(price[i] + x[0], (m_row[i] if per_row else m_row) + x[1])   # this row, or stop
    return i + 1


def _charge_all(b, price, m_each, block=1024):
    """Charge every element of a list whose operation has no partial answer, in blocks
    (a deadline is still noticed between them)."""
    n = len(price)
    for i in range(0, n, block):
        k = min(n, i + block) - i
        b.charge(sum(price[i:i + k]), k * m_each)


# the claim, its evidence cache entry, its ints and floats
_CLAIM_BYTES = record_bytes(Claim) + 360
FLOAT_BYTES = _B.FLOAT
_WALK_BYTES = record_bytes(Walk) + 3 * LIST + 2 * DICT_KEY + 200
_LEVERS_LIST = ('narrow_labels', 'narrow_arm')
# record coordinates (a graphlet with a coordinates block only): what a CoordClaim adds to
# a claim -- its slot and a clipped RunCoordinates with its interval tuple and ints -- and
# per occurrence a new (start, end) tuple and its slot (a clip makes new ones)
_COORD_ROW_BYTES = 8 + record_bytes(_C.RunCoordinates) + TUPLE + 3 * INT
_COORD_OCC_BYTES = _C._OCC_BYTES + 8


def _arm_coords(g, side, b=None):
    """The arm's RunCoordinates (a tuple indexed by run id), or None: the graphlet carries
    no coordinates block (or is in annotate mode, which never does). |b|: the parsed index
    charged at its cold price (once per call)."""
    if g.mode != 'constrain' or not _C.present(g):
        return None
    c = _C.of(g, b)
    return c.arms.get(side)


def _coord_price(rc):
    """(lwu, bytes) of attaching run coordinates |rc| to a row (clipped or not)."""
    n = len(rc.intervals)
    return W_ELEM + _B.W_COORD * n, _COORD_ROW_BYTES + n * _COORD_OCC_BYTES

__all__ = [
    'label', 'labels_matching', 'path_id', 'leaf_segment', 'spell', 'walks', 'claims',
    'rank_walks', 'walks_at', 'walk_filter_reason',
    'label_walks', 'routes', 'support_profile', 'support_changes', 'label_summary',
    'splits', 'continuation', 'next_request', 'next_requests', 'resubmittable_names',
    'subgraph',
    'GraphletView', 'view_from_saved', 'view_from_spec',
    'compare', 'compare_cost', 'comparability', 'index_identity', 'memory_bytes',
    'cache_bytes', 'cache_signature', 'segment_support',
    'summary', 'evidence_block',
    'END_CLASS', 'alive_at', 'limitation_dict', 'resource_stop_dict', 'informational',
]


# ------------------------------------------------------------------ labels

def label(g, selector):
    """Tagged selectors: {'ref': 'h:<column>:<seq_id>' | 'c:<column>'} |
    {'name': str} | {'id': int}. A bare string resolves only when exactly one label
    matches it as a ref OR as a name (a column may be NAMED 'c:0'; annotate mode can
    record two headers 'ACC1' from different columns); a bare int is an id."""
    if isinstance(selector, Label):
        if selector.id < len(g.labels) and g.labels[selector.id].ref == selector.ref:
            return g.labels[selector.id]
        return label(g, {'ref': selector.ref})
    if isinstance(selector, bool):
        raise BadSelector('not a label selector: %r' % (selector,))
    if isinstance(selector, int):
        selector = {'id': selector}
    if isinstance(selector, str):
        index = _label_index_after_scans(g)
        if index is None:
            hits = [l for l in g.labels if l.ref == selector or l.name == selector]
        else:
            r, n = index[0].get(selector), index[1].get(selector)
            if r is None or n is None:
                ids = r or n or ()
            else:
                ids = sorted(set(r) | set(n))     # a label may be NAMED like its own ref
            hits = [g.labels[i] for i in ids]
        if len(hits) == 1:
            return hits[0]
        if not hits:
            raise UnknownLabel('no label is named or referenced %r' % selector)
        raise AmbiguousLabel(
            '%r matches %d labels (%s): use {"ref": ...}' % (
                selector, len(hits), ', '.join('%s %r' % (l.ref, l.name) for l in hits)))
    if not isinstance(selector, dict) or len(selector) != 1:
        raise BadSelector('a label selector is {"ref": ...}, {"name": ...} or {"id": ...}: '
                          '%r' % (selector,))
    (key, value), = selector.items()
    if key == 'id':
        if isinstance(value, bool):
            # a bool is an int to Python: {'id': True} would select label 1, where a bare
            # bool is refused
            raise BadSelector('a label id is an int, not %r' % (value,))
        if not isinstance(value, int) or not 0 <= value < len(g.labels):
            raise UnknownLabel('no label id %r' % (value,))
        return g.labels[value]
    if key == 'ref' or key == 'name':
        # refs and names are strings: any other value matches no label
        index = _label_index_after_scans(g)
        if not isinstance(value, str):
            hits = []
        elif index is None:
            hits = [l for l in g.labels if l.ref == value] if key == 'ref' else \
                [l for l in g.labels if l.name == value]
        else:
            hits = [g.labels[i] for i in index[0 if key == 'ref' else 1].get(value, ())]
    else:
        raise BadSelector('unknown selector tag %r' % key)
    if len(hits) == 1:
        return hits[0]
    if not hits:
        raise UnknownLabel('no label with %s %r' % (key, value))
    raise AmbiguousLabel('%d labels named %r (%s): select by ref' % (
        len(hits), value, ', '.join(l.ref for l in hits)))


# lookups by name or ref answered by a scan of the labels before a graphlet's index is
# built: building it costs about as much as four or five scans (measured on the real
# fixtures), so a caller that looks up a few labels on a fresh model pays no more than the
# scans did, and one that looks up many pays at most about twice what the index alone would
# have cost (building the index on the first lookup would make it 2x a scan, up to 5.7x on
# large tables)
_LABEL_SCANS = 4


def _label_index_after_scans(g):
    """The name/ref index, or None while this graphlet's lookups are still answered by a
    scan (the first _LABEL_SCANS of them)."""
    got = g.cache.get('label_index')
    if got is not None:
        return got
    n = g.cache.get('label_scans', 0)
    if n < _LABEL_SCANS:
        g.cache['label_scans'] = n + 1
        return None
    return _label_index(g)


def _label_index(g):
    """(ref -> label ids, name -> label ids), ids ascending, so that a lookup by name or ref
    does not scan every label; names may repeat (annotate mode: two headers of one name),
    so each maps to a tuple."""
    got = g.cache.get('label_index')
    if got is None:
        got = g.cache['label_index'] = (_ids_by((l.ref, l.id) for l in g.labels),
                                        _ids_by((l.name, l.id) for l in g.labels))
    return got


def _ids_by(pairs):
    """{key: (ids in order)} of (key, id) pairs, keys in first-seen order. Built as the
    final tuples (a list per key only for a key seen twice): a list per key turned into
    tuples in a second dict would hold both dicts, every list and every tuple at once --
    2.5-2.9x the index, more than a lookup by name charges for it."""
    out = {}
    dup = None
    for k, i in pairs:
        got = out.get(k)
        if got is None:
            out[k] = (i,)
            continue
        if dup is None:
            dup = {}
        more = dup.get(k)
        if more is None:
            more = dup[k] = list(got)
        more.append(i)
    if dup:
        for k, v in dup.items():
            out[k] = tuple(v)
    return out


def labels_matching(g, selectors):
    if selectors is None:
        return None
    if isinstance(selectors, (str, dict, int, Label)):
        selectors = [selectors]
    return [label(g, s) for s in selectors]


def _label_ids(g, selectors):
    got = labels_matching(g, selectors)
    return None if got is None else {l.id for l in got}


# ------------------------------------------------------------------ walks by id

def path_id(arm, x):
    if isinstance(x, (Walk, Path)):
        return x.path_id if isinstance(x, Walk) else x.id
    if isinstance(x, bool) or not isinstance(x, int):
        raise BadSelector('a walk is a path id, a Path or a Walk: %r' % (x,))
    if not 0 <= x < len(derive.paths(arm)):
        raise IndexError('the %s arm has %d walks, no walk %d'
                         % (arm.side, len(derive.paths(arm)), x))
    return x


def leaf_segment(arm, x):
    return derive.paths(arm)[path_id(arm, x)].leaf


# ------------------------------------------------------------------ spelling

def spell(g, arm, leaf, orientation='natural', with_seed=False, *, budget=None):
    """The walk's bases. natural (§2.4): right = seed + flank, left = flank + seed
    (seed only with with_seed). walk: walking order (outward index i is [i]); with the
    seed, the seed is read in the walking direction too (reversed on the left arm),
    so position |seed| + i holds outward base i on both arms. budget= (local limits): the
    spelling is charged as one step at its price (the walk's chain and bases)."""
    arm = g.arm(arm)
    b = _B.resolve(budget)
    if b is None:
        seg = leaf_segment(arm, leaf) if not isinstance(leaf, Path) else leaf.leaf
    else:
        with b.scope('spell', ('select_walks',)):
            # admitted before the walk is resolved: resolving a path id builds the arm's
            # paths (a Path per walk), which a refused call must not leave behind
            derive.uses(b, g, arm, 'leaves', 'paths')
            seg = leaf_segment(arm, leaf) if not isinstance(leaf, Path) else leaf.leaf
            n = arm.segments[seg].end_bp + len(g.seed.sequence)
            d = arm.segments[seg].depth + 1
            b.charge(W_ROW + W_STEP * d + (n >> 8), list_bytes(d) + 3 * (STR + n))
            b.phase = 'spell'
    if orientation == 'natural':
        flank = derive.natural_flank(arm, seg)
        if not with_seed:
            return flank
        return g.seed.sequence + flank if arm.side == 'right' else flank + g.seed.sequence
    if orientation == 'walk':
        w = derive.walk_bases(arm, seg)
        if not with_seed:
            return w
        return (g.seed.sequence if arm.side == 'right' else g.seed.sequence[::-1]) + w
    raise ValueError("orientation is 'natural' or 'walk'")


def _chain_bases(arm, seg_id, lo, hi):
    """Outward bases [lo, hi) along seg's first-parent chain, natural orientation."""
    w = derive.walk_bases(arm, seg_id)[lo:hi]
    return w if arm.side == 'right' else w[::-1]


# ------------------------------------------------------------------ alive sets

def _segment_ops(arm, s):
    """Constrain: the label changes inside segment |s| in the order of the G.end rule
    (derive.label_changes())."""
    return derive.label_changes(s, derive.runs_by_segment(arm)[s.id], arm.runs)


def alive_at(g, arm, seg, pos):
    """The labels alive at outward base |pos| of segment |seg| (empty unless from_bp <=
    pos < from_bp + length_bp). Constrain: the entry set with every run end and switch-in
    of the segment up to |pos| applied in order (the rule of G.end and of the displayed
    support profile); annotate: the recorded set of the P run holding |pos|."""
    a = g.arm(arm)
    s = a.segments[seg]
    if not s.from_bp <= pos < s.end_bp:
        return frozenset()
    if g.mode != 'constrain':
        for p in s.presence:
            if p.from_bp <= pos < p.to_bp:
                return frozenset(p.labels)
        return frozenset()
    cur = set(s.entry)
    for at, kind, lab in _segment_ops(a, s):
        if at > pos:
            break
        if kind == 0:
            cur.discard(lab)
        else:
            cur.add(lab)
    return frozenset(cur)


def _alive_profile(g, arm, leaf_seg, b=None):
    """Constrain: the labels alive at each base of the leaf's first-parent chain, as
    maximal runs [(from, to, frozenset)]. Entry sets, then per segment the run ends
    anchored there (removed from to_bp on) and the switch-ins (added from at_bp on) in
    the chronological order of the G.end rule. |b|: each segment charged before it is
    read (its pieces and their sets)."""
    out = []
    chain = derive.chain(arm, leaf_seg)
    if b is not None:
        pw, pm = _alive_prices(b, g, arm)
        b.charge(sum([pw[x] for x in chain]), sum([pm[x] for x in chain]))
    for sid in chain:
        out.extend(_alive_pieces(arm, arm.segments[sid]))
    return _merge_runs(out)


def _alive_prices(b, g, arm):
    """(lwu, bytes) per segment of its pieces of displayed support (constrain): at most one
    piece per label change inside it and one per base, each a set (cached)."""
    def prices():
        ops_ = derive.segment_ops(arm)
        w, m = [], []
        for s in arm.segments:
            k = ops_[s.id]
            pieces = 1 + min(k, s.length_bp)
            n = len(s.entry) + k
            w.append(2 + pieces + k + ((pieces + 1) * n >> 7))
            m.append(pieces * (TUPLE + 2 * INT + set_bytes(n)) + set_bytes(n))
        return w, m
    derive.uses(b, g, arm, 'segment_ops')
    return derive.price_list(b, g, arm, 'alive_prices', prices)


def _alive_pieces(arm, s):
    """Constrain: the alive sets inside one segment as [(from, to, frozenset)] (not
    merged; nothing for a zero-length segment)."""
    return list(_iter_alive_pieces(arm, s))


def _iter_alive_pieces(arm, s):
    """_alive_pieces() one piece at a time: a comparison's cut reads them in order and
    stops at the cut, holding one piece's set instead of every piece's."""
    if s.length_bp == 0:
        return
    ops = _segment_ops(arm, s)
    cur = set(s.entry)
    pos = s.from_bp
    i = 0
    while pos < s.end_bp:
        while i < len(ops) and ops[i][0] <= pos:
            if ops[i][1] == 0:
                cur.discard(ops[i][2])
            else:
                cur.add(ops[i][2])
            i += 1
        nxt = ops[i][0] if i < len(ops) else s.end_bp
        nxt = min(max(nxt, pos + 1), s.end_bp)
        yield pos, nxt, frozenset(cur)
        pos = nxt


def _presence_profile(arm, leaf_seg, b=None, g=None):
    out = []
    chain = derive.chain(arm, leaf_seg)
    if b is not None:
        def prices():
            w, m = [], []
            for s in arm.segments:
                pr = s.presence
                n = sum(len(p.labels) for p in pr)
                w.append(W_ELEM * (1 + 2 * len(pr)) + (n >> _B.ID_SHIFT_C))
                m.append(len(pr) * (TUPLE + 3 * INT + SET) + n * SET_ITEM)
            return w, m
        pw, pm = derive.price_list(b, g, arm, 'presence_prices', prices)
        b.charge(sum([pw[x] for x in chain]), sum([pm[x] for x in chain]))
    for sid in chain:
        for p in arm.segments[sid].presence:
            out.append((p.from_bp, p.to_bp, frozenset(p.labels), p.total))
    merged = []
    for f, t, s, tot in out:
        if merged and merged[-1][2] == s and merged[-1][3] == tot and merged[-1][1] == f:
            merged[-1] = (merged[-1][0], t, s, tot)
        else:
            merged.append((f, t, s, tot))
    return merged


def _merge_runs(runs):
    out = []
    for f, t, s in runs:
        if out and out[-1][2] == s and out[-1][1] == f:
            out[-1] = (out[-1][0], t, s)
        else:
            out.append((f, t, s))
    return out


def _labels(g, ids):
    return [g.labels[i] for i in sorted(ids)]


# ------------------------------------------------------------------ claims

END_CLASS = {
    'D': 'dead_end', 'T': 'lost', 'L': 'lost', 'B': 'pruned', 'R': 'pruned', 'W': 'pruned',
    'U': 'blocked', 'V': 'blocked', 'J': 'blocked', 'X': 'radius',
    'S': 'capped', 'P': 'capped', 'N': 'capped', 'O': 'capped', 'M': 'capped', 'Y': 'capped',
}


def _first_paths(arm):
    """seg -> the smallest path id whose first-parent chain holds it (None for a
    segment reached only through non-first parents)."""
    got = arm.cache.get('first_paths')
    if got is None:
        segs = arm.segments
        got = [None] * len(segs)
        pol = derive.path_of_leaf(arm)
        for s in reversed(segs):
            best = pol.get(s.id)
            for c in s.children:
                if segs[c].parents[0] == s.id and got[c] is not None:
                    best = got[c] if best is None else min(best, got[c])
            got[s.id] = best
        arm.cache['first_paths'] = got
    return got


def _run_claim(g, arm, run, cut, exact=True, coords=None):
    """The claim of one run end, or None when the run-start guard drops it at |cut|.
    |exact|: the arm's label evidence is exact (arm_exact). |coords|: the arm's
    RunCoordinates (_arm_coords) -- a CoordClaim then carries the run's occurrences clipped
    to the cut. None (the default) makes a plain Claim: only the public paths that price
    coordinates pass them (claims(), walks()); compare() never does, so it never clips
    coordinates it would not charge."""
    route_from, ev_from = derive.evidence(arm, run)
    seg = arm.segments[run.segment]
    zero = run.from_bp == run.to_bp
    if cut is not None:
        if zero:
            if run.from_bp > cut:
                return None
        elif run.from_bp >= cut:
            return None              # the run does not exist yet at the cut
    if run.merged:
        kind, end_class = 'merged', 'open'
    elif run.silent:
        kind, end_class = 'diverged', 'lost'
    elif zero:
        kind, end_class = 'boundary', END_CLASS[run.code]
    elif (seg.leaf is not None and run.to_bp == seg.end_bp
          and seg.leaf.path_reason is not None):
        kind, end_class = 'alive', END_CLASS[run.code]
    else:
        kind, end_class = 'end', END_CLASS[run.code]
    to_bp = run.to_bp
    evidence_from = ev_from
    reason, qualifier = run.reason, run.qualifier
    loss, branches = run.loss, run.branches
    if cut is not None and to_bp > cut:
        to_bp = cut
        end_class = 'open'
        reason = qualifier = None
        loss = branches = None       # terminal values are valid at to_bp only
        if ev_from < cut:
            kind = 'stretch'
        else:
            kind, evidence_from = 'route_only', None
    elif zero:
        evidence_from = None
    if coords is not None:
        return CoordClaim(
            label=g.labels[run.label], arm=arm.side,
            path_id=_first_paths(arm)[run.segment], segment=run.segment,
            from_bp=run.from_bp, to_bp=to_bp, route_bp=route_from,
            evidence_from=evidence_from, kind=kind, reason=reason, qualifier=qualifier,
            end_class=end_class, entered_by=run.entered_by,
            from_label=None if run.from_label is None else g.labels[run.from_label],
            cost=run.cost, loss=loss, branches=branches, support=g.support, exact=exact,
            run=run.id, coordinates=_C.clip(coords[run.id], arm.side, cut))
    return Claim(
        label=g.labels[run.label], arm=arm.side, path_id=_first_paths(arm)[run.segment],
        segment=run.segment, from_bp=run.from_bp, to_bp=to_bp, route_bp=route_from,
        evidence_from=evidence_from, kind=kind, reason=reason, qualifier=qualifier,
        end_class=end_class, entered_by=run.entered_by,
        from_label=None if run.from_label is None else g.labels[run.from_label],
        cost=run.cost, loss=loss, branches=branches, support=g.support, exact=exact,
        run=run.id)


def _displayed_presence(arm):
    """Annotate mode, one pass parents-first along the displayed (first-parent) chains:
    seg -> (labels recorded on every node from the seed boundary to the segment's end,
    [(label, j)] dropped inside the segment at base j). O(recorded runs), where walking
    every walk's chain would be quadratic on a comb-shaped trie."""
    got = arm.cache.get('displayed_presence')
    if got is None:
        segs = arm.segments
        got = [None] * len(segs)
        sets = {}                            # equal sets shared (a comb repeats a few)
        for s in segs:
            alive = set(s.entry) if not s.parents else set(got[s.parents[0]][0])
            drops = []
            for pr in s.presence:
                if not alive:
                    break
                still = alive.intersection(pr.labels)
                drops.extend((l, pr.from_bp) for l in sorted(alive - still))
                alive = still
            fs = frozenset(alive)
            got[s.id] = (sets.setdefault(fs, fs), drops or ())
        arm.cache['displayed_presence'] = got
    return got


def _annotate_evidence(arm, l, anchor):
    """Annotate mode, the §5.1 merge walk for an end of label |l| anchored at |anchor|:
    the start of the merge nearest the anchor on its first-parent chain p through whose
    FIRST parent no route of |l| arrives -- the witness route of routes() enters it
    through another parent, so the displayed bases of p before it are not the label's
    (route_bp); 0 when a route of the label follows p from the seed boundary. The
    annotate reading of derive.evidence(): 'the lineage's label is in partition[0]'
    becomes 'a route of |l| reaches the end of parents[0]' (the union rule of direct_bp;
    annotate partitions are empty). Relative to p even when no walk displays p (the
    claim's path_id is then None), as in constrain mode. Cached only where p has a
    merge (elsewhere it is 0 at once)."""
    above = derive.merge_above(arm)
    m = above[anchor]
    if m is None:
        return 0
    cache = arm.cache.setdefault('annotate_evidence', {})
    got = cache.get((l, anchor))
    if got is not None:
        return got
    alive_end, _ = _annotate_route_ends(arm)
    segs = arm.segments
    got = 0
    while m is not None:
        seg = segs[m]
        if l not in alive_end[seg.parents[0]]:
            got = seg.from_bp
            break
        m = above[seg.parents[0]]
    cache[(l, anchor)] = got
    return got


def _annotate_ends(arm):
    """Annotate mode: [(label, anchor, j, route_bp)], one per maximal end of a label's
    label-consistent routes under the §5.1 union rule (the rule of direct_bp), with the
    merge-derived start of its displayed support; ordered by (first displayed walk
    through the anchor, label, j). Built per call from the cached ends (a list per
    end would cost as much as the ends again)."""
    _, ends = _annotate_route_ends(arm)
    first = _first_paths(arm)
    none = len(first)
    got = [(l, s, j, _annotate_evidence(arm, l, s)) for l, es in ends.items() for s, j in es]
    got.sort(key=lambda x: (none if first[x[1]] is None else first[x[1]], x[0], x[2], x[1]))
    return got


def _annotate_claim(g, arm, l, anchor, j, ev, cut, exact):
    """One route end [0, j) of label |l| as a claim: route support [0, j), displayed
    support [ev, j) on the anchor's first-parent chain (route_bp = ev, as for a run in
    constrain mode); None when the run-start guard drops it (a route with bases makes no
    claim at a cut of 0). At a cut D < j: a stretch when displayed support survives the
    cut, route_only when none does."""
    if cut is not None and j > 0 and cut <= 0:
        return None
    seg = arm.segments[anchor]
    code = seg.leaf.path_reason if seg.leaf is not None and j == seg.end_bp else None
    if code is None:
        reason, end_class, kind = 'label_lost', 'lost', 'end'
    else:
        reason, end_class = REASON[code], END_CLASS[code]
        kind = 'alive' if code == 'X' or code in RESOURCE_CODES else 'end'
    to_bp = j
    evidence_from = ev
    if j == 0:
        kind, evidence_from = 'boundary', None
    if cut is not None and j > cut:
        to_bp, end_class, reason = cut, 'open', None
        if ev < cut:
            kind = 'stretch'
        else:
            kind, evidence_from = 'route_only', None
    return Claim(label=g.labels[l], arm=arm.side, path_id=_first_paths(arm)[anchor],
                 segment=anchor, from_bp=0, to_bp=to_bp, route_bp=ev,
                 evidence_from=evidence_from, kind=kind, reason=reason, qualifier=None,
                 end_class=end_class, entered_by='seed', from_label=None, cost=None,
                 loss=None, branches=None, support=g.support,
                 exact=exact, run=None)


def _annotate_route_ends(arm):
    """Annotate mode, the §5.1 union rule of label_summary's direct_bp: -> (alive_end,
    ends). alive_end[seg] = the labels recorded on every node of SOME route from the
    seed boundary to the segment's end (the union over the parents' sets, then every P
    run of the segment); ends = {label: [(seg, j)]}, the maximal ends of its routes: a
    drop inside a segment (anywhere in the root), or the end of a segment where it is
    alive and no child records it at its first node (a zero-length child, which records
    nothing, carries every label on to its own end). A drop at a child's first base is
    not an end of its own: the route ends with the parent, unless another child carries
    the label on. The longest end of a label is its direct_bp."""
    got = arm.cache.get('annotate_route_ends')
    if got is not None:
        return got
    got = _route_ends_pass(arm, None)
    arm.cache['annotate_route_ends'] = got
    return got


def _route_ends_pass(arm, depth):
    """The pass of _annotate_route_ends() -> (alive_end, ends), over the whole DAG (|depth|
    None: equal alive sets shared, a comb repeats a few) or over the DAG restricted to
    [0, depth) (_annotate_route_ends_at(): the segments starting before the depth, P runs
    and children cut there, ends at most at the depth)."""
    lim = math.inf if depth is None else depth
    segs = arm.segments
    alive_end = [None] * len(segs)
    ends = {}
    sets = {} if depth is None else None
    for s in segs:
        if s.from_bp >= lim:
            continue
        if s.parents:
            alive = set()
            for p in s.parents:
                alive.update(alive_end[p])
        else:
            alive = set(s.entry)
        for pr in s.presence:
            if not alive or pr.from_bp >= lim:
                break
            still = alive.intersection(pr.labels)
            if pr.from_bp > s.from_bp or not s.parents:
                for l in alive - still:
                    ends.setdefault(l, []).append((s.id, pr.from_bp))
            alive = still
        fs = frozenset(alive)
        alive_end[s.id] = fs if sets is None else sets.setdefault(fs, fs)
    for s in segs:
        alive = alive_end[s.id]
        if not alive:
            continue
        carried = set()
        for c in s.children:
            if segs[c].from_bp >= lim:
                continue
            if not segs[c].presence:
                carried = alive
                break
            carried.update(alive.intersection(segs[c].presence[0].labels))
        for l in sorted(alive - carried):
            ends.setdefault(l, []).append((s.id, min(s.end_bp, lim)))
    return alive_end, ends


def _annotate_routes(arm, l):
    """Annotate mode: one witness route per maximal end of label |l| -> [(route = the
    segments root -> anchor, anchor, j)]. Back from the anchor, every merge is passed
    through its FIRST parent (in parents order) whose set holds the label, so a route
    that can follow the displayed chain does. Displayed routes come first (by their
    first walk), then route-only ones by length."""
    alive_end, ends = _annotate_route_ends(arm)
    segs = arm.segments
    first = _first_paths(arm)
    out = []
    for s, j in ends.get(l, ()):
        route = [s]
        while segs[route[-1]].parents:
            ps = segs[route[-1]].parents
            route.append(next((p for p in ps if l in alive_end[p]), ps[0]))
        route.reverse()
        out.append((route, s, j))
    npaths = len(derive.paths(arm))

    def key(x):
        route, s, j = x
        shown = first[s] is not None and route == derive.chain(arm, s)
        return (first[s] if shown else npaths, j, s)
    out.sort(key=key)
    return out


def _route_bases(arm, route, hi):
    """Bases [0, hi) along |route| (segment ids root -> anchor, any parents), natural
    orientation; None without bases."""
    parts = []
    for x in route:
        w = arm.segments[x].walk
        if w is None:
            return None
        parts.append(w)
    b = ''.join(parts)[:hi]
    return b if arm.side == 'right' else b[::-1]


def _displayed_leaves_below(arm, seg):
    """The leaves whose displayed (first-parent) chain passes through |seg|."""
    segs = arm.segments
    out, stack = [], [seg]
    while stack:
        x = stack.pop()
        if segs[x].leaf is not None:
            out.append(x)
        stack.extend(c for c in segs[x].children if segs[c].parents[0] == x)
    return sorted(out)


def claims(g, arm=None, labels=None, at_most_bp=None, strict=True, *, budget=None,
           resume=None):
    """The §6.9 unit: one claim per RUN END (claims, not leaves). In annotate mode the
    claims are the oracle's: per label recorded at the seed boundary its maximal
    label-consistent routes -- the §5.1 union rule of direct_bp, through ANY parent of a
    merge -- each with route support [0, to_bp) and displayed support [evidence_from,
    to_bp) on the first-parent chain of its anchor, kept apart exactly as for a run in
    constrain mode (route_bp > 0: the route joins that chain through a non-first merge
    parent; at a cut with no displayed support left it is route_only). One annotate claim
    per maximal END, with one witness route (routes()) chosen by the stored parent order:
    it does not enumerate every route to that end, nor establish a contiguous source
    occurrence. strict=True raises IncompleteRecording when a recorded list was cut (the
    oracle would be a lower bound).

    at_most_bp keeps the §5.1 cut of a claim (evidence at the run's anchored endpoint,
    then clipped): what the deeper DAG says about [0, D). compare() keys the claims of
    each DAG RESTRICTED to the depth instead, which is what a retrieval walked to D says
    (_restricted_claims).

    A retrieval with record coordinates (Graphlet.coordinates) gives CoordClaims: each
    carries its run's occurrences clipped to the claim (.coordinates, a RunCoordinates:
    where the run's own bases [from_bp, to_bp) lie in its record -- a route_only claim's
    route too; 0-based half-open on the forward strand, DESIGN §18).

    budget= (local limits, budget.py): a stop raises LocalBudgetExceeded whose .partial holds
    the whole claims made so far, in this order, and a resume token (resume=) that goes on
    after the last row; without a budget (the default) no work or allocation budget
    applies."""
    b = _B.resolve(budget)
    if b is None:
        if resume is None:
            return _claims_plain(g, arm, labels, at_most_bp, strict)
        return _claims(g, arm, labels, at_most_bp, strict, None, resume)
    with b.scope('claims', _LEVERS_LIST):
        return _claims(g, arm, labels, at_most_bp, strict, b, resume)


def _claims_plain(g, arm, labels, at_most_bp, strict):
    """claims() without a budget or a resume token: the plain loop (no position kept, no
    charge point tested per run)."""
    sides = [g.arm(arm).side] if arm is not None else [s for s in ARM_SIDES if s in g.arms]
    wanted = _label_ids(g, labels)
    out = []
    for side in sides:
        a = g.arms[side]
        exact = arm_exact(g, side)
        if g.mode == 'constrain':
            rcs = _arm_coords(g, side)
            for run in a.runs:
                if wanted is not None and run.label not in wanted:
                    continue
                c = _run_claim(g, a, run, at_most_bp, exact, rcs)
                if c is not None:
                    out.append(c)
        else:
            if strict and a.labels_per_node.nodes_truncated:
                raise IncompleteRecording(
                    'the %s arm has %d recorded label lists cut at labels.max_labels_per_node:'
                    ' claims filtered from them are lower bounds (pass strict=False to accept,'
                    ' or raise the knob)' % (side, a.labels_per_node.nodes_truncated))
            for l, anchor, j, ev in _annotate_ends(a):
                if wanted is not None and l not in wanted:
                    continue
                c = _annotate_claim(g, a, l, anchor, j, ev, at_most_bp, exact)
                if c is not None:
                    out.append(c)
    return out


def _claims(g, arm, labels, at_most_bp, strict, b, resume):
    sides = [g.arm(arm).side] if arm is not None else [s for s in ARM_SIDES if s in g.arms]
    out = []
    si = i = 0
    args = None
    try:
        wanted = _label_ids_b(g, labels, b)
        args = [sides, sorted(wanted) if wanted is not None else None, at_most_bp,
                bool(strict)]
        s0, i0 = _resume_pos(g, 'claims', args, resume, 2)
        si, i = s0, i0
        for si, side in enumerate(sides):
            if si < s0:
                continue
            a = g.arms[side]
            exact = arm_exact(g, side)
            i = start = i0 if si == s0 else 0
            if g.mode == 'constrain':
                runs = a.runs
                paid = len(runs)
                if b is not None:
                    derive.uses(b, g, a, 'leaves', 'paths', 'path_of_leaf', 'merge_above',
                                'merge_depth', 'first_paths', 'partition_sets')
                    price = _claim_prices(b, g, a)
                    m_row = _CLAIM_BYTES
                    rcs = _arm_coords(g, side, b)
                    if rcs is not None:
                        # each claim also carries its run's occurrences
                        price, m_row = _claim_coord_prices(b, g, a, price, rcs)
                    keep = None if wanted is None else (lambda t: runs[t].label in wanted)
                    b.charge(max(0, len(runs) - start))         # the scan of the runs
                    b.phase = 'rows'
                    paid = start
                else:
                    rcs = _arm_coords(g, side)
                for i in range(start, len(runs)):
                    run = runs[i]
                    if wanted is not None and run.label not in wanted:
                        continue
                    if i >= paid:
                        paid = _pay_rows(b, price, i, len(runs), keep, m_row)
                    c = _run_claim(g, a, run, at_most_bp, exact, rcs)
                    if c is not None:
                        out.append(c)
                i = len(runs)
            else:
                if strict and a.labels_per_node.nodes_truncated:
                    raise IncompleteRecording(
                        'the %s arm has %d recorded label lists cut at labels.max_labels_per_'
                        'node: claims filtered from them are lower bounds (pass strict=False '
                        'to accept, or raise the knob)'
                        % (side, a.labels_per_node.nodes_truncated))
                if b is not None:
                    _annotate_ends_price(b, g, a)
                    b.phase = 'rows'
                ends = _annotate_ends(a)
                paid = len(ends)
                if b is not None:
                    price = [W_ROW] * len(ends)
                    keep = None if wanted is None else (lambda t: ends[t][0] in wanted)
                    paid = start
                for i in range(start, len(ends)):
                    l, anchor, j, ev = ends[i]
                    if wanted is not None and l not in wanted:
                        continue
                    if i >= paid:
                        paid = _pay_rows(b, price, i, len(ends), keep, _CLAIM_BYTES)
                    c = _annotate_claim(g, a, l, anchor, j, ev, at_most_bp, exact)
                    if c is not None:
                        out.append(c)
                i = len(ends)
    except LocalBudgetExceeded as e:
        # a stop before the call knew its arguments resumes where this call began
        raise e.restate(Partial(out, resume if args is None else
                                _token(g, 'claims', args, si, i)), rows=len(out),
                        position=i)
    return out


def _claim_prices(b, g, a):
    """Per run of the arm the price of its claim: a row and the evidence walk over the
    merges above its anchor (structural, cached)."""
    md = derive.merge_depth(a)
    return derive.price_list(b, g, a, 'claim_prices',
                             lambda: [W_ROW + 2 + 3 * md[r.segment] for r in a.runs])


def _claim_coord_prices(b, g, a, price, rcs):
    """(lwu, bytes) per run of a claim that carries its run's coordinates: the claim's
    price plus the occurrences attached (a structural list, cached on the arm)."""
    def build():
        w, m = [], []
        for p, rc in zip(price, rcs):
            cw, cm = _coord_price(rc)
            w.append(p + cw)
            m.append(_CLAIM_BYTES + cm)
        return w, m
    return derive.price_list(b, g, a, 'claim_coord_prices', build)


def _restricted_claim_prices(b, g, a):
    md = derive.merge_depth(a)
    return derive.price_list(b, g, a, 'restricted_claim_prices',
                             lambda: [W_ROW + 4 + 6 * md[r.segment] for r in a.runs])


def _annotate_ends_price(b, g, a):
    """Charge _annotate_ends(): the union-rule route ends, one evidence walk per end and
    the sort, at their cold price (the ends are a fact of the model: the count is read
    from the derivation, charged first)."""
    derive.uses(b, g, a, 'leaves', 'paths', 'path_of_leaf', 'merge_above', 'merge_depth',
                'first_paths', 'annotate_route_ends')
    md = derive.merge_depth(a)
    ends = _annotate_route_ends(a)[1]
    n = sum(len(v) for v in ends.values())
    b.charge(W_ELEM * len(ends) + n)
    walk = sum(3 * md[s] for es in ends.values() for s, _ in es)
    b.uses((id(a), 'annotate_ends'), 2 * n + walk + sort_work(n),
           list_bytes(n) + n * (2 * TUPLE + 2 * INT + 120))


# ------------------------------------------------------------------ walks

def _labels_full(g, arm, p):
    """Labels with displayed support on the whole walk [0, length)."""
    if g.mode == 'constrain':
        ends = derive.end_labels(arm, p.leaf)
        return [g.labels[e.label] for e in ends
                if derive.evidence(arm, arm.runs[e.run])[1] == 0]
    return _labels(g, _displayed_presence(arm)[p.leaf][0])


def walks(g, arm, *, top=None, by='support', labels=None, route_consistent=True, min_bp=0,
          budget=None, resume=None):
    """Walk(path_id, leaf, segments, length_bp, sequence, path_reason, end_reasons,
    claims, labels_full, n_alive, complete, beyond_certified_bp). 'support' ranks by
    (-len(labels_full), -n_alive, -length_bp, path_id); 'length' by (-length_bp,
    path_id); 'loss' by (min end-label loss, -n_alive, -length_bp, path_id); 'id' keeps
    path order. Ranking uses the per-leaf records only; chains, spellings and claims
    are produced for the walks returned. A retrieval with record coordinates gives
    CoordWalks: .coordinates maps each label alive at the leaf (its ref) to the
    RunCoordinates of its run reaching the leaf, and the claims are CoordClaims.

    budget= (local limits): a stop raises LocalBudgetExceeded. With by='id' its .partial holds
    the whole walks made so far and a resume token (resume=); a ranked list never returns
    a partial (a ranking of part of the arm would be another answer)."""
    b = _B.resolve(budget)
    if resume is not None and by != 'id':
        raise ValueError("resume= continues walks(by='id') only: a ranked list is never "
                         "partial")
    if b is None and resume is None:
        a = g.arm(arm)
        ranked = _ranked_walks(g, a, by, _label_ids(g, labels), route_consistent, min_bp)
        if top is not None:
            ranked = ranked[:top]
        return _walk_objects(g, a, ranked, with_claims=True)
    if b is None:
        return _walks_b(g, arm, top, by, labels, route_consistent, min_bp, None, resume)
    with b.scope('walks', ('rank_id', 'select_walks', 'narrow_labels', 'narrow_arm')):
        return _walks_b(g, arm, top, by, labels, route_consistent, min_bp, b, resume)


def _walks_b(g, arm, top, by, labels, route_consistent, min_bp, b, resume):
    a = g.arm(arm)
    if by != 'id':
        wanted = _label_ids_b(g, labels, b)
        ranked = _ranked_walks(g, a, by, wanted, route_consistent, min_bp, b)
        if top is not None:
            ranked = ranked[:top]
        return _walk_objects(g, a, ranked, True, b)
    out = []
    args = None
    k0 = 0
    try:
        wanted = _label_ids_b(g, labels, b)
        args = [a.side, top, sorted(wanted) if wanted is not None else None,
                bool(route_consistent), min_bp]
        k0, = _resume_pos(g, 'walks', args, resume, 1)
        ranked = _ranked_walks(g, a, by, wanted, route_consistent, min_bp, b)
        if top is not None:
            ranked = ranked[:top]
        return _walk_objects(g, a, ranked[k0:], True, b, out)
    except LocalBudgetExceeded as e:
        raise e.restate(Partial(out, resume if args is None else
                                _token(g, 'walks', args, k0 + len(out))),
                        rows=len(out), position=k0 + len(out))


_WALK_KEYS = {
    'support': lambda x: (-len(x[1]), -len(x[2]), -x[0].length_bp, x[0].id),
    'length': lambda x: (-x[0].length_bp, x[0].id),
    'loss': lambda x: (x[3], -len(x[2]), -x[0].length_bp, x[0].id),
    'id': lambda x: x[0].id,
}


def _walk_keys_of(g, a, p):
    """(full, alive_ids, consistent, min_loss) of one walk, from its per-leaf records."""
    if g.mode == 'constrain':
        ends = derive.end_labels(a, p.leaf)
        alive_ids = [e.label for e in ends]
        consistent = {e.label for e in ends if e.route_bp == 0}
        min_loss = min((e.loss for e in ends), default=math.inf)
        return None, alive_ids, consistent, min_loss
    full = _labels_full(g, a, p)
    return full, list(a.segments[p.leaf].end), {l.id for l in full}, 0.0


def _rank_setup(b, g, a):
    """Charge what ranking an arm's walks reads, at its cold price; -> (lwu, bytes) per
    path id: the price of a walk's keys (structural, cached)."""
    if g.mode == 'constrain':
        derive.uses(b, g, a, 'leaves', 'paths', 'end_labels', 'merge_above', 'merge_depth',
                    'partition_sets')
        md = derive.merge_depth(a)
        rbs = derive.runs_by_segment(a)

        def prices():
            w, m = [], []
            for p in derive.paths(a):
                n = len(rbs[p.leaf])
                w.append(2 * W_ELEM + n * (4 + 3 * md[p.leaf]))
                m.append(TUPLE + 2 * LIST + 2 * _B.INT_ID * n + SET + SET_ITEM * n + 100 * n)
            return w, m
        return derive.price_list(b, g, a, 'rank_prices', prices)
    derive.uses(b, g, a, 'leaves', 'paths', 'displayed_presence')
    dp = _displayed_presence(a)
    segs = a.segments

    def prices():
        w, m = [], []
        for p in derive.paths(a):
            n, k = len(dp[p.leaf][0]), len(segs[p.leaf].end)
            w.append(2 * W_ELEM + sort_work(n) + ((n + k) >> 1))
            m.append(TUPLE + 3 * LIST + _B.INT_ID * (2 * n + k) + 2 * set_bytes(n))
        return w, m
    return derive.price_list(b, g, a, 'rank_prices_annotate', prices)


def _ranked_walks(g, a, by, wanted, route_consistent, min_bp, b=None):
    """[(path, labels_full, alive ids, min loss)] in walks() order; nothing is spelled."""
    if by not in _WALK_KEYS:
        raise ValueError("by is 'support', 'length', 'loss' or 'id'")
    ps = None
    ranked = []
    try:
        if b is not None:
            _rank_keys_charge(b, g, a, min_bp)
        # built once admitted: a refused call leaves no paths behind
        ps = derive.paths(a)
        for p in ps:
            if p.length_bp < min_bp:
                continue
            full, alive_ids, consistent, min_loss = _walk_keys_of(g, a, p)
            if wanted is not None:
                pool = consistent if route_consistent else set(alive_ids)
                if not wanted & pool:
                    continue
            if full is None:
                full = _labels_full(g, a, p)
            ranked.append((p, full, alive_ids, min_loss))
        if b is not None and by != 'id':
            b.charge(sort_work(len(ranked)))
    except LocalBudgetExceeded as e:
        raise e.restate(phase='rank', ranked=len(ranked), of=_n_walks(a))
    ranked.sort(key=_WALK_KEYS[by])
    return ranked


def _n_walks(a):
    """The arm's walk count without building its paths (a stop states it: the stop may
    come before the paths are admitted)."""
    got = a.cache.get('leaves')
    if got is not None:
        return len(got)
    return sum([1 for s in a.segments if s.leaf is not None])


def _rank_keys_charge(b, g, a, min_bp):
    """What ranking charges before it reads a key: its setup (cold price) and every walk's
    keys (a ranking has no partial answer), in blocks. The paths are read once their
    setup is admitted."""
    pw, pm = _rank_setup(b, g, a)
    ps = derive.paths(a)
    b.phase = 'rank'
    for i in range(0, len(ps), 1024):
        block = [p.id for p in ps[i:i + 1024] if p.length_bp >= min_bp]
        b.charge(sum([pw[t] for t in block]), sum([pm[t] for t in block]))


_RANK_LEVERS = ('rank_id', 'narrow_labels', 'narrow_arm')


def _rank_walks_replay(g, arm, by, labels, min_bp, n, b):
    """Charge |b| what rank_walks(g, arm, by=by, labels=labels, min_bp=min_bp) charges
    when it returns |n| ids -- at the same charge points, with the same derivations marked
    used -- without ranking: a pager that keeps a ranking (graphlet_walks) pays its cold
    price on every page. A lump of the units a ranking once charged would add the call's
    base a second time and leave its derivations unmarked, so that the walks a page then
    builds would pay for them again (1.4-1.8x on cached pages)."""
    a = g.arm(arm)
    with b.scope('rank_walks', _RANK_LEVERS):
        _label_ids_b(g, labels, b)
        got = 0
        try:
            _rank_keys_charge(b, g, a, min_bp)
            got = n
            if by != 'id':
                b.charge(sort_work(n))
        except LocalBudgetExceeded as e:
            raise e.restate(phase='rank', ranked=got, of=_n_walks(a))
        b.charge(n, list_bytes(n))


def rank_walks(g, arm, *, by='support', labels=None, route_consistent=True, min_bp=0,
               budget=None):
    """The path ids walks() returns, in its order, from the per-leaf records only: no
    chain, spelling or claim is built (a pager ranks the whole arm, then materializes
    one page with walks_at()). budget=: a stop raises, with no partial."""
    b = _B.resolve(budget)
    a = g.arm(arm)
    if b is None:
        return [x[0].id for x in _ranked_walks(g, a, by, _label_ids(g, labels),
                                               route_consistent, min_bp)]
    with b.scope('rank_walks', _RANK_LEVERS):
        ranked = _ranked_walks(g, a, by, _label_ids_b(g, labels, b), route_consistent,
                               min_bp, b)
        b.charge(len(ranked), list_bytes(len(ranked)))
        return [x[0].id for x in ranked]


def walk_filter_reason(g, arm, walk, labels):
    """Why route_consistent=True leaves out a walk that carries one of |labels| at its
    end: 'merge_entered' -- the label reaches the walk's end only along its own route,
    which joins the displayed walk through a non-first merge parent (constrain: every
    such label has T.route_bp > 0; annotate: a route under the §5.1 union rule records it
    on every node from the seed boundary to the end, the displayed chain does not) -- or,
    annotate mode only, 'not_from_seed': no route records it on every node from the seed
    boundary (it is recorded at the end but appeared, or had a gap, on the way)."""
    a = g.arm(arm)
    if g.mode == 'constrain':
        return 'merge_entered'
    p = derive.paths(a)[path_id(a, walk)]
    alive_end, _ = _annotate_route_ends(a)
    wanted = _label_ids(g, labels)
    return 'merge_entered' if wanted is None or wanted & alive_end[p.leaf] \
        else 'not_from_seed'


def walks_at(g, arm, path_ids, with_claims=False, *, budget=None, resume=None):
    """The Walk objects of |path_ids|, in that order (claims only with with_claims).
    budget=: a stop raises LocalBudgetExceeded whose .partial holds the whole walks made
    so far and a resume token (resume=)."""
    a = g.arm(arm)
    b = _B.resolve(budget)
    if b is None and resume is None:
        return _walks_at(g, a, path_ids, with_claims, None)
    path_ids = list(path_ids)
    args = [a.side, [x if isinstance(x, int) else path_id(a, x) for x in path_ids],
            bool(with_claims)]
    k0, = _resume_pos(g, 'walks_at', args, resume, 1)
    if b is None:
        return _walks_at(g, a, path_ids[k0:], with_claims, None)
    with b.scope('walks_at', ('select_walks',)):
        out = []
        try:
            return _walks_at(g, a, path_ids[k0:], with_claims, b, out)
        except LocalBudgetExceeded as e:
            raise e.restate(Partial(out, _token(g, 'walks_at', args, k0 + len(out))),
                            rows=len(out), position=k0 + len(out), of=len(path_ids))


def _walks_at(g, a, path_ids, with_claims, b, out=None):
    # the setup admitted before the paths are built (a refused call leaves none)
    price = _rank_setup(b, g, a) if b is not None else None
    ps = derive.paths(a)
    rows = []
    for pid in path_ids:
        p = ps[path_id(a, pid)]
        if price is not None:
            b.charge(price[0][p.id], price[1][p.id])
        full, alive_ids, _, min_loss = _walk_keys_of(g, a, p)
        rows.append((p, full if full is not None else _labels_full(g, a, p), alive_ids,
                     min_loss))
    return _walk_objects(g, a, rows, with_claims, b, out)


def _leaf_claims(g, a, leaves_, rcs=None):
    """leaf -> the claims ending at its last node, in claims() order: exactly what
    filtering claims(g, a, strict=False) by (segment in |leaves_|, to_bp = the leaf's end)
    gives, built for those leaves only (the claims of the whole arm cost as much as
    every walk's, for a page of a few)."""
    by_leaf = {}
    exact = arm_exact(g, a.side)
    segs = a.segments
    if g.mode == 'constrain':
        rbs = derive.runs_by_segment(a)
        for leaf in sorted(leaves_):
            end = segs[leaf].end_bp
            got = [_run_claim(g, a, a.runs[r], None, exact, rcs) for r in rbs[leaf]
                   if a.runs[r].to_bp == end]
            if got:
                by_leaf[leaf] = got
        return by_leaf
    # annotate: _annotate_ends() orders by (first walk through the anchor, label, j,
    # anchor); within one anchor that is (label, j) -- each leaf's own ends, read from the
    # arm's index of them (_annotate_leaf_ends), not from a pass over every route end of
    # the arm per call (a budgeted walks() would make that pass per chunk of 32 walks)
    ends = _annotate_leaf_ends(a)
    for s in sorted(leaves_):
        for l, j in ends.get(s, ()):
            by_leaf.setdefault(s, []).append(
                _annotate_claim(g, a, l, s, j, _annotate_evidence(a, l, s), None, exact))
    return by_leaf


_WALK_CHUNK = 32


def _annotate_leaf_ends(a):
    """Annotate mode: segment -> the route ends at its last node, [(label, j)] sorted (the
    claims of a walk ending there, in claims() order), from one pass over the arm's route
    ends (cached beside them)."""
    got = a.cache.get('annotate_leaf_ends')
    if got is None:
        segs = a.segments
        got = {}
        ends = _annotate_route_ends(a)[1]
        # the labels in ascending order: a label has at most one end at a segment's last
        # node, so each segment's list comes out sorted by (label, j)
        for l in sorted(ends):
            for s, j in ends[l]:
                if j == segs[s].end_bp:
                    got.setdefault(s, []).append((l, j))
        a.cache['annotate_leaf_ends'] = got
    return got


def _annotate_leaf_ends_price(b, g, a):
    """Charge _annotate_leaf_ends() at its cold price (after the route ends it reads): the
    sort of the labels, a visit per route end, and the index's dict, lists and pairs."""
    ends = _annotate_route_ends(a)[1]
    e = sum(len(v) for v in ends.values())
    b.uses((id(a), 'annotate_leaf_ends'), sort_work(len(ends)) + W_ELEM * e,
           list_bytes(len(ends)) + dict_bytes(e) + e * (LIST + LIST_ITEM + TUPLE + 16))


def _leaf_end_counts(a):
    """Annotate mode: leaf -> the route ends at its last node (its claims in a walk)."""
    segs = a.segments
    return {s: len(v) for s, v in _annotate_leaf_ends(a).items() if segs[s].leaf is not None}


def _walk_objects(g, a, ranked, with_claims, b=None, out=None):
    """The Walk objects of |ranked|. Under a budget, chunk by chunk: each chunk is
    charged at its price before it is built (its chains, spellings and claims), so a
    stop leaves |out| with whole walks only."""
    if b is not None:
        out = [] if out is None else out
        segs = a.segments
        md = rbs = ends_n = None
        # a retrieval with record coordinates: each walk maps its end labels to their
        # runs' occurrences, and each claim carries them (charged with the walk)
        rcs = _arm_coords(g, a.side, b)
        if rcs is not None:
            rbs_c = derive.runs_by_segment(a)
            derive.uses(b, g, a, 'end_labels')
        if with_claims:
            if g.mode == 'constrain':
                derive.uses(b, g, a, 'merge_above', 'merge_depth', 'first_paths',
                            'partition_sets', 'path_of_leaf')
                md = derive.merge_depth(a)
                rbs = derive.runs_by_segment(a)
            else:
                _annotate_ends_price(b, g, a)
                # the index of each leaf's route ends, built once: each chunk reads its own
                # leaves' (not a pass over every route end of the arm per chunk of 32 walks)
                _annotate_leaf_ends_price(b, g, a)
                ends_n = _leaf_end_counts(a)
        b.phase = 'rows'
        for k in range(0, len(ranked), _WALK_CHUNK):
            chunk = ranked[k:k + _WALK_CHUNK]
            w = 0
            m = 0
            for x in chunk:
                p = x[0]
                d = segs[p.leaf].depth + 1
                w += W_ROW + 2 * W_STEP * d + (p.length_bp >> 8)
                # the chain and the spelling, and the batch's kept prefixes (at most as
                # much again)
                m += _WALK_BYTES + 2 * list_bytes(d) + 4 * (STR + p.length_bp) \
                    + LIST_ITEM * len(x[1]) + b.row_extra[1]
                w += b.row_extra[0]
                if with_claims:
                    if rbs is not None:
                        n = len(rbs[p.leaf])
                        w += n * (W_ROW + 2 + 3 * md[p.leaf])
                    else:
                        n = ends_n.get(p.leaf, 0)
                        w += n * W_ROW
                    m += n * _CLAIM_BYTES
                if rcs is not None:
                    # the end labels' map (its entries refer to the index's records) and,
                    # with claims, every claim's occurrences at the leaf (at most its runs')
                    ends = derive.end_labels(a, p.leaf)
                    w += W_ELEM * len(ends)
                    m += _B.dict_bytes(len(ends)) + 8
                    if with_claims:
                        for r in rbs_c[p.leaf]:
                            cw, cm = _coord_price(rcs[r])
                            w += cw
                            m += cm
            b.charge(w, m)
            out.extend(_walk_objects(g, a, chunk, with_claims))
        return out
    by_leaf = {}
    rcs = _arm_coords(g, a.side)
    if ranked and with_claims:
        by_leaf = _leaf_claims(g, a, {x[0].leaf for x in ranked}, rcs)
    # chains and spellings of the walks returned, sharing their common prefixes
    batch = derive.walk_batch(a, [x[0].leaf for x in ranked]) if len(ranked) > 1 else None
    seen = set()
    out = []
    for p, full, alive_ids, _ in ranked:
        seg = a.segments[p.leaf]
        reasons = {}
        if g.mode == 'constrain':
            for e in derive.end_labels(a, p.leaf):
                name = a.runs[e.run].reason
                reasons[name] = reasons.get(name, 0) + 1
        code = seg.leaf.path_reason
        if batch is None:
            chain = p.segments
            try:
                sequence = derive.natural_flank(a, p.leaf)
            except ValueError:
                sequence = None
        else:
            chain, sequence = batch[p.leaf]
            if sequence is not None and a.side == 'left':
                sequence = sequence[::-1]
            if p.leaf in seen:
                chain = list(chain)          # a walk asked twice: each Walk owns its list
            seen.add(p.leaf)
        fields = dict(
            path_id=p.id, leaf=p.leaf, segments=chain, length_bp=p.length_bp,
            sequence=sequence, path_reason=REASON[code] if code else None,
            end_reasons=reasons, claims=by_leaf.get(p.leaf, []), labels_full=full,
            n_alive=len(alive_ids),
            complete=(code is None or code not in RESOURCE_CODES)
            and p.length_bp <= a.complete_to_bp,
            beyond_certified_bp=max(0, p.length_bp - a.complete_to_bp))
        if rcs is None:
            out.append(Walk(**fields))
        else:
            # each label alive at the leaf -> the occurrences of its run reaching the leaf
            out.append(CoordWalk(**fields, coordinates={
                g.labels[e.label].ref: rcs[e.run] for e in derive.end_labels(a, p.leaf)}))
    return out


# ------------------------------------------------------------------ per label

def routes(g, sel, arm, spell=False, *, budget=None, resume=None):
    """Per run of the label (constrain): the segments root -> anchor of the LABEL'S OWN
    route, chosen at every merge through the G partition holding the lineage's
    incoming label (not the displayed first-parent chain). spell=True adds the route's
    bases (natural orientation, from the seed boundary to the run's end).

    Annotate mode (no runs, empty partitions): one witness route per maximal end of the
    label's label-consistent routes under the §5.1 union rule (the rule of direct_bp,
    so the direct_bp route is always among them), passing every merge through its first
    parent, in parents order, that records the label -- the displayed chain where it
    can. A route through a non-first merge parent is shown by no walk. The witness is
    one route, deterministically chosen: it does not enumerate the other routes reaching
    the same end, and it establishes no contiguous source occurrence (the recorded sets
    say each node carries the label, not that one indexed sequence spells the route).

    One route per run or end, each as long as the DAG makes it. Without a budget (the
    default) no work or allocation budget applies; budget= (local limits): a stop raises
    LocalBudgetExceeded whose .partial holds the whole routes made so far, in this order,
    and a resume token (resume=)."""
    return _B.run_scoped(budget, 'routes', ('narrow_arm',),
                         lambda b: _routes(g, sel, arm, spell, b, resume))


def _route_depths(b, g, a):
    """Charge what a label's own routes read; -> seg -> the most segments a route root ->
    seg can have (the longest chain through any parents): a route's price."""
    derive.uses(b, g, a, 'runs_by_label', 'partition_sets', 'route_depth')
    return derive.route_depth(a)


def _run_route(a, run, spell, segs):
    """routes() of one run: (route, bases) with |spell|, else the route."""
    route = [run.segment]
    s = run.segment
    seg = segs[s]
    while seg.parents:
        nxt = seg.parents[0]
        if len(seg.parents) > 1:
            at = derive.lineage_label_at(a, run, seg.from_bp)
            for j, part in enumerate(derive.merge_parts(a, s)):
                if at in part:
                    nxt = seg.parents[j]
                    break
        route.append(nxt)
        s = nxt
        seg = segs[s]
    route.reverse()
    if not spell:
        return route
    # a list, not a generator: join() builds one from a generator first (the same text)
    bases = ''.join([segs[x].walk or '' for x in route])[:run.to_bp]
    if segs[route[0]].walk is None:
        bases = None
    elif a.side == 'left':
        bases = bases[::-1]
    return route, bases


def _route_price(a, segs, rd, anchor, to_bp, spell):
    d = rd[anchor]
    w = W_ROW + W_STEP * d
    m = TUPLE + list_bytes(d)
    if spell:
        n = segs[anchor].end_bp
        w += d + (n >> 8)
        m += 3 * (STR + n)
    return w, m


def _routes(g, sel, arm, spell, b, resume):
    a = g.arm(arm)
    segs = a.segments
    out = []
    args = None
    k0 = 0
    try:
        if b is not None and not isinstance(sel, (int, Label)):
            _uses_g(b, g, 'label_index')
        lab = label(g, sel)
        if b is not None or resume is not None:
            # the token's arguments: only a budgeted call stops, only a resumed one reads them
            args = [a.side, lab.id, bool(spell)]
            k0, = _resume_pos(g, 'routes', args, resume, 1)
        if g.mode != 'constrain':
            if b is not None:
                rd = _annotate_routes_price(b, g, a, lab.id)
            rows = _annotate_routes(a, lab.id)
            for k in range(k0, len(rows)):
                route, anchor, j = rows[k]
                if b is not None and spell:
                    b.charge(*_route_price(a, segs, rd, anchor, j, True))
                out.append((route, _route_bases(a, route, j)) if spell else route)
            return out
        ids = derive.label_runs(a, lab.id) if b is None else None
        if b is not None:
            rd = _route_depths(b, g, a)
            ids = derive.label_runs(a, lab.id)
            b.phase = 'rows'
        for k in range(k0, len(ids)):
            run = a.runs[ids[k]]
            if b is not None:
                b.charge(*_route_price(a, segs, rd, run.segment, run.to_bp, spell))
            out.append(_run_route(a, run, spell, segs))
        return out
    except LocalBudgetExceeded as e:
        raise e.restate(Partial(out, resume if args is None else
                                _token(g, 'routes', args, k0 + len(out))),
                        rows=len(out), position=k0 + len(out))


def _annotate_routes_price(b, g, a, l):
    """Charge _annotate_routes() of label |l| (every witness route of its ends, the sort
    by displayed walk); -> the route depths. The ends are read from the union-rule
    derivation, charged first."""
    # first_paths is built from path_of_leaf: charged with it
    derive.uses(b, g, a, 'annotate_route_ends', 'first_paths', 'leaves', 'paths',
                'path_of_leaf', 'route_depth')
    rd = derive.route_depth(a)
    es = _annotate_route_ends(a)[1].get(l, ())
    b.charge(W_ELEM * len(es) + sum(5 * W_STEP * rd[s] for s, _ in es) + sort_work(len(es)),
             list_bytes(len(es)) + sum(TUPLE + 2 * list_bytes(rd[s]) for s, _ in es))
    b.phase = 'rows'
    return rd


def label_walks_routes(g, sel, arm=None, *, budget=None, resume=None):
    """label_walks() rows with each walk's own route (the segments routes() gives for
    it), made and charged together: [(LabelWalk, route)] (a stop's Partial holds such
    pairs). What graphlet_labels(name=) pages."""
    return _B.run_scoped(budget, 'label_walks', ('narrow_arm',),
                         lambda b: _label_walks(g, sel, arm, b, resume, True))


def label_walks(g, sel, arm=None, *, budget=None, resume=None):
    """Per run of the label: (from, to, evidence_from, sequence = the label's own route
    bases [from_bp, to_bp), end = the end reason as the label_end event states it
    ('switched' for a silent Lw end, 'merged' for a run closed by a merge),
    leaves_below, merged_into).

    Annotate mode: per witness route of routes() (from the seed boundary, run None):
    leaves_below the leaves whose displayed walk passes through the route's anchor (the
    walks that show where the stretch ends, as a run's leaves_below in constrain mode; []
    when the anchor is on no displayed walk), evidence_from the claim's: where the
    route's bases become those of the anchor's first-parent chain (0 when the route is
    that chain, the start of the last merge it enters through a non-first parent
    otherwise -- §5.1, as for a run in constrain mode).

    A retrieval with record coordinates gives CoordLabelWalks: .coordinates is the run's
    RunCoordinates (where its own bases lie in its record).

    budget= (local limits): a stop raises LocalBudgetExceeded whose .partial holds the whole
    label walks made so far, in this order, and a resume token (resume=)."""
    return _B.run_scoped(budget, 'label_walks', ('narrow_arm',),
                         lambda b: _label_walks(g, sel, arm, b, resume))


_LABEL_WALK_BYTES = record_bytes(LabelWalk) + 120


def _label_walks(g, sel, arm, b, resume, with_routes=False):
    out = []
    si = k = 0
    args = None
    try:
        if b is not None and not isinstance(sel, (int, Label)):
            _uses_g(b, g, 'label_index')
        # the label before the arm: an unknown label with a bad arm is UnknownLabel
        lab = label(g, sel)
        sides = [g.arm(arm).side] if arm is not None else \
            [s for s in ARM_SIDES if s in g.arms]
        args = [sides, lab.id]
        s0, k0 = _resume_pos(g, 'label_walks', args, resume, 2)
        si, k = s0, k0
        for si, side in enumerate(sides):
            if si < s0:
                continue
            a = g.arms[side]
            segs = a.segments
            k = start = k0 if si == s0 else 0
            if g.mode != 'constrain':
                if b is not None:
                    rd = _annotate_routes_price(b, g, a, lab.id)
                    derive.uses(b, g, a, 'merge_above', 'merge_depth', 'subtree_bounds')
                    md = derive.merge_depth(a)
                    first = derive.subtree_bounds(a)[1]
                rows = _annotate_routes(a, lab.id)
                for k in range(start, len(rows)):
                    route, anchor, j = rows[k]
                    if b is not None:
                        w, m = _route_price(a, segs, rd, anchor, j, True)
                        x = b.row_extra
                        b.charge(w + 2 + 3 * md[anchor] + 2 * first[anchor] + x[0],
                                 m + _LABEL_WALK_BYTES + list_bytes(first[anchor]) + x[1])
                    seg = segs[anchor]
                    if seg.leaf is not None and j == seg.end_bp:
                        code = seg.leaf.path_reason
                        end = REASON[code] if code else None
                    else:
                        end = 'label_lost'
                    lw = LabelWalk(None, lab, side, 0, j,
                                   _annotate_evidence(a, lab.id, anchor),
                                   _route_bases(a, route, j), end,
                                   _displayed_leaves_below(a, anchor), None)
                    out.append((lw, route) if with_routes else lw)
                k = len(rows)
                continue
            if b is not None:
                rd = _route_depths(b, g, a)
                derive.uses(b, g, a, 'merge_above', 'merge_depth', 'subtree_bounds')
                md = derive.merge_depth(a)
                every = derive.subtree_bounds(a)[0]
                b.phase = 'rows'
            rcs = _arm_coords(g, side, b)
            ids = derive.label_runs(a, lab.id)
            for k in range(start, len(ids)):
                run = a.runs[ids[k]]
                if b is not None:
                    w, m = _route_price(a, segs, rd, run.segment, run.to_bp, True)
                    n = every[run.segment]
                    x = b.row_extra
                    if rcs is not None:
                        # the run's occurrences, the index's own record (one more slot)
                        w += W_ELEM + _B.W_COORD * len(rcs[run.id].intervals)
                        m += 8
                    b.charge(w + 2 + 3 * md[run.segment] + 3 * n
                             + len(segs[run.segment].children) + x[0],
                             m + _LABEL_WALK_BYTES + set_bytes(n) + list_bytes(n) + x[1])
                route, bases = _run_route(a, run, True, segs)
                if bases is not None:
                    bases = bases[run.from_bp:run.to_bp] if side == 'right' else \
                        bases[len(bases) - run.to_bp:len(bases) - run.from_bp]
                merged_into = None
                if run.merged:
                    for c in segs[run.segment].children:
                        if len(segs[c].parents) > 1:
                            merged_into = c
                end = 'merged' if run.merged else 'switched' if run.silent \
                    else run.event_reason
                if rcs is None:
                    lw = LabelWalk(run.id, lab, side, run.from_bp, run.to_bp,
                                   derive.evidence(a, run)[1], bases, end,
                                   derive.run_leaves(a, run), merged_into)
                else:
                    lw = CoordLabelWalk(run.id, lab, side, run.from_bp, run.to_bp,
                                        derive.evidence(a, run)[1], bases, end,
                                        derive.run_leaves(a, run), merged_into,
                                        coordinates=rcs[run.id])
                out.append((lw, route) if with_routes else lw)
            k = len(ids)
    except LocalBudgetExceeded as e:
        raise e.restate(Partial(out, resume if args is None else
                                _token(g, 'label_walks', args, si, k)), rows=len(out),
                        position=k)
    return out


def support_profile(g, arm, leaf, kind='displayed', *, budget=None):
    """Maximal SupportRun(from_bp, to_bp, labels, exact) over the walk [0, length).
    displayed: the labels supporting the displayed bases (§5.1); route: also every run
    anchored on the walk over its own [from_bp, to_bp) (route support, scope-free).
    Annotate mode records one set per node: both kinds are the recorded sets.

    budget= (local limits): a stop raises LocalBudgetExceeded with no partial -- the last
    maximal run of a cut profile might go on past the cut, so no prefix of runs is whole."""
    return _B.run_scoped(budget, 'support_profile', ('select_walks',),
                         lambda b: _support_profile(g, arm, leaf, kind, b))


_SUPPORT_RUN_BYTES = record_bytes(SupportRun) + LIST


def _support_out(g, b, runs, exact_of):
    """The SupportRun rows of |runs| [(from, to, set, total?)], charged before they are
    built (their labels sorted into Label lists; the profile has no partial answer)."""
    if b is not None:
        b.charge(sum([W_ELEM + sort_work(len(r[2])) + (len(r[2]) >> 1) for r in runs]),
                 sum([_SUPPORT_RUN_BYTES + LIST_ITEM * len(r[2]) for r in runs]))
    return [exact_of(r) for r in runs]


def _support_profile(g, arm, leaf, kind, b):
    a = g.arm(arm)
    if b is not None:
        # admitted before the walk is resolved (its path id builds the arm's paths)
        derive.uses(b, g, a, 'leaves', 'paths')
        if g.mode == 'constrain':
            derive.uses(b, g, a, 'segment_ops')
    seg_id = leaf_segment(a, leaf)
    if b is not None:
        b.charge(W_STEP * (a.segments[seg_id].depth + 1),
                 list_bytes(a.segments[seg_id].depth + 1))
        b.phase = 'profile'
    if g.mode != 'constrain':
        return _support_out(g, b, _presence_profile(a, seg_id, b, g),
                            lambda r: SupportRun(r[0], r[1], _labels(g, r[2]),
                                                 r[3] == len(r[2]), r[3]))
    prof = _alive_profile(g, a, seg_id, b)
    if kind == 'route':
        extra = []
        rbs = derive.runs_by_segment(a)
        chain = derive.chain(a, seg_id)
        if b is not None:
            b.charge(sum(len(rbs[sid]) for sid in chain) + len(chain))
        for sid in chain:
            for r in rbs[sid]:
                run = a.runs[r]
                if run.to_bp > run.from_bp:
                    extra.append(run)
        if extra:
            cuts = sorted({0, a.segments[seg_id].end_bp}
                          | {x for f, t, _ in prof for x in (f, t)}
                          | {x for r in extra for x in (r.from_bp, r.to_bp)})
            if b is not None:
                n = len(cuts)
                widest = max((len(x[2]) for x in prof), default=0) + len(extra)
                b.charge(sort_work(n) + n * (len(prof) + len(extra))
                         + ((n * widest) >> _B.ID_SHIFT_C),
                         n * (TUPLE + 2 * INT + set_bytes(widest)))
            pieces = []
            for f, t in zip(cuts, cuts[1:]):
                base = next((s for pf, pt, s in prof if pf <= f < pt), frozenset())
                add = {r.label for r in extra if r.from_bp <= f < r.to_bp}
                pieces.append((f, t, frozenset(base | add)))
            prof = _merge_runs(pieces)
    elif kind != 'displayed':
        raise ValueError("kind is 'displayed' or 'route'")
    return _support_out(g, b, prof,
                        lambda r: SupportRun(r[0], r[1], _labels(g, r[2]), True, len(r[2])))


def support_changes(g, arm, leaf, *, budget=None):
    """Where the displayed support changes along a walk, and why: an end (the R code,
    qualifier text included), a switch (Lw away / a switch-in), a split (labels that
    took another branch), a merge (a lineage joining through another parent).

    budget= (local limits): a stop raises LocalBudgetExceeded with no partial (the changes are
    read off the whole profile)."""
    return _B.run_scoped(budget, 'support_changes', ('select_walks',),
                         lambda b: _support_changes(g, arm, leaf, b))


_CHANGE_BYTES = record_bytes(Change) + 2 * LIST + 200
# per label that changed: its reason dict and its label's dict (with the ref's text), its
# slot in the added or removed set, the sorted list and the Change's list -- measured 521
# bytes kept and 636-695 at the peak per label on a split of 2,000-20,000 labels
_REASON_BYTES = 2 * dict_bytes(4) + SET_ITEM + 2 * LIST_ITEM + STR + 16


def _support_changes(g, arm, leaf, b):
    a = g.arm(arm)
    if b is not None:
        # admitted before the walk is resolved (its path id builds the arm's paths): the
        # profile below charges the same derivations first, once per call
        derive.uses(b, g, a, 'leaves', 'paths')
    seg_id = leaf_segment(a, leaf)
    chain = derive.chain(a, seg_id)
    segs = a.segments
    prof = _support_profile(g, a, leaf, 'displayed', b)
    prev = set(segs[chain[0]].entry)
    out = []
    on_chain = set(chain)
    rbs = derive.runs_by_segment(a)
    if b is not None:
        derive.uses(b, g, a, 'partition_sets')
        b.phase = 'changes'
        # a reason scans the walk's chain -- its segments, their runs and events -- per label
        # that changed, until it finds the label's: charged at that bound (the added-label
        # scan never stops early)
        per_label = len(chain) + sum(len(rbs[x]) + len(segs[x].events) for x in chain)
    for run in prof:
        if b is not None:
            b.charge(W_ELEM + (len(run.labels) >> 1), set_bytes(len(run.labels)))
        cur = {l.id for l in run.labels}
        at = run.from_bp
        added, removed = cur - prev, prev - cur
        if added or removed:
            if b is not None:
                n = len(added) + len(removed)
                b.charge(W_ROW + n * (3 + per_label), _CHANGE_BYTES + n * _REASON_BYTES)
            out.append(Change(at, _labels(g, added), _labels(g, removed),
                              _change_reasons(g, a, chain, on_chain, rbs, at, added, removed,
                                              b)))
        prev = cur
    # labels alive at the leaf's last node stay; nothing to report past the end
    return out


def _has_id(ids, l):
    """Whether the ascending label set |ids| (an array('I'): every set the parser makes is
    sorted, the codec refuses any other) holds |l|: a bisection, where `l in ids` scanned
    the array from its start."""
    i = bisect.bisect_left(ids, l)
    return i < len(ids) and ids[i] == l


def _split_point(segs, chain, at, b):
    """The split a label removed at |at| may have taken: the first segment of the walk's
    chain that starts there with one parent -> (the parent's end set, [(each other child's
    entry set, its first base)], the price of testing one label against them); False when
    there is none. One per Change, so the sets are not scanned again per removed label."""
    if b is not None:
        b.charge(W_STEP * len(chain))
    for sid in chain:
        s = segs[sid]
        if s.from_bp == at and len(s.parents) == 1:
            parent = segs[s.parents[0]]
            k = len(parent.children)
            if b is not None:
                b.charge(W_ELEM * k, list_bytes(k) + k * (TUPLE + STR + 8))
            sibs = [(segs[c].entry, segs[c].walk[:1] if segs[c].walk
                     else segs[c].first_base or '')
                    for c in parent.children if c != sid]
            # a bisection of each set per label: W_ELEM and a unit per halving
            price = sum(W_ELEM + len(x).bit_length() for x in [parent.end] + [e for e, _ in sibs])
            return parent.end, sibs, price
    return False


def _change_reasons(g, a, chain, on_chain, rbs, at, added, removed, b=None):
    segs = a.segments
    reasons = []
    split = None
    for l in sorted(removed):
        why = None
        if g.mode == 'constrain':
            for sid in chain:
                for r in rbs[sid]:
                    run = a.runs[r]
                    if run.label == l and run.to_bp == at and not run.merged:
                        if run.silent:
                            to = sorted({ev.to_label for s2 in chain for ev in segs[s2].events
                                         if ev.type == 'switch' and ev.at_bp == at
                                         and ev.from_label == l})
                            why = {'why': 'switched', 'to': [g.labels[x].as_dict() for x in to]}
                        else:
                            why = {'why': run.event_reason, 'end': run.end}
                        break
                if why:
                    break
        if why is None:
            if split is None:
                split = _split_point(segs, chain, at, b)
            if split:
                end, sibs, price = split
                if b is not None:
                    b.charge(price)
                took = sorted({x for entry, bases in sibs for x in bases if _has_id(entry, l)})
                if took and _has_id(end, l):
                    # a split only when another branch took the label: on the parent's last
                    # node and on no child's first node it is absent (annotate) or ended,
                    # not split
                    why = {'why': 'split', 'took': took}
        if why is None:
            why = {'why': 'absent' if g.mode != 'constrain' else 'ended'}
        reasons.append(dict(why, label=g.labels[l].as_dict(), change='removed'))
    for l in sorted(added):
        why = None
        for sid in chain:
            s = segs[sid]
            for ev in s.events:
                if ev.type == 'switch' and ev.at_bp == at and ev.to_label == l:
                    why = {'why': 'switch', 'from': g.labels[ev.from_label].as_dict(),
                           'cost': ev.cost}
            if len(s.parents) > 1 and s.from_bp == at:
                for j, part in enumerate(derive.merge_parts(a, sid)):
                    if l in part and j > 0:
                        why = {'why': 'merge', 'via_parent': s.parents[j]}
        if why is None:
            why = {'why': 'present' if g.mode != 'constrain' else 'entered'}
        reasons.append(dict(why, label=g.labels[l].as_dict(), change='added'))
    return reasons


def splits(g, arm, min_labels_before=0, *, budget=None, resume=None):
    """The split records of an arm, in the walker's order (at_bp, then first child), as
    SplitPoint(arm, at_bp, segment, kind, labels_before, branches): oriented like the
    other accessors -- at_bp is the outward distance from the seed boundary (walking
    order on both arms, as a claim's interval), labels are Label objects ({name, ref}),
    and each Branch carries its first base, the labels at its first node (the recorded
    list, cut at labels.max_labels_per_node, with the true count) and the first walk
    whose displayed chain takes it. A branch's label count is not a share of
    labels_before: a label may follow several branches (kind 'ambiguous').

    budget= (local limits): a stop raises LocalBudgetExceeded whose .partial holds the whole
    split points made so far, in this order, and a resume token (resume=)."""
    return _B.run_scoped(budget, 'splits', ('narrow_arm',),
                         lambda b: _splits(g, arm, min_labels_before, b, resume))


_SPLIT_BYTES = record_bytes(SplitPoint) + LIST
_BRANCH_BYTES = record_bytes(Branch) + LIST


def _splits(g, arm, min_labels_before, b, resume):
    a = g.arm(arm)
    args = [a.side, min_labels_before]
    k0, = _resume_pos(g, 'splits', args, resume, 1)
    out = []
    k = k0
    sps = ()
    try:
        paid = None
        if b is not None:
            derive.uses(b, g, a, 'splits', 'leaves', 'paths', 'path_of_leaf', 'first_paths')
            b.phase = 'rows'
        first = _first_paths(a)
        sps = derive.splits(a, g.mode)
        if b is not None:
            b.charge(max(0, len(sps) - k0))           # the scan of the split records
            pw, pm = _split_prices(b, g, a, sps)
            keep = lambda t: sps[t].labels_before >= min_labels_before  # noqa: E731
            paid = k0
        for k in range(k0, len(sps)):
            sp = sps[k]
            if sp.labels_before < min_labels_before:
                continue
            if paid is not None and k >= paid:
                paid = _pay_rows(b, pw, k, len(sps), keep, pm)
            branches = []
            for c in sp.children:
                child = a.segments[c]
                branches.append(Branch(c, derive.first_base(child),
                                       [g.labels[l] for l in child.entry[:g.cap]],
                                       child.entry_total, first[c]))
            out.append(SplitPoint(a.side, sp.at_bp, sp.segment,
                                  'ambiguous' if sp.ambiguous else 'divergence',
                                  sp.labels_before, branches))
    except LocalBudgetExceeded as e:
        raise e.restate(Partial(out, _token(g, 'splits', args, k)), rows=len(out),
                        position=k)
    return out


def _split_prices(b, g, a, sps):
    """(lwu, bytes) per split record of the arm's SplitPoint row (its branches and their
    labels, cut at the cap), cached."""
    def prices():
        w, m = [], []
        for sp in sps:
            n = sum(min(len(a.segments[c].entry), g.cap) for c in sp.children)
            w.append(W_ROW + W_ELEM * len(sp.children) + (n >> 1))
            m.append(_SPLIT_BYTES + len(sp.children) * _BRANCH_BYTES + LIST_ITEM * n)
        return w, m
    return derive.price_list(b, g, a, 'split_prices', prices)


def label_summary(g, *, budget=None):
    """{ref: {name, ref, <arm>: {direct_bp, reach_bp, reentries, runs}}} (both modes,
    §5.1). direct_bp is label-consistent ROUTE support, not a contiguous occurrence.
    budget= (local limits): a stop raises, with no partial (a summary of some labels is not
    the table)."""
    return _B.run_scoped(budget, 'label_summary', ('narrow_labels',),
                         lambda b: _label_summary(g, b))


def _label_summary_prices(b, g):
    """(lwu, bytes) per label of label_summary()'s row (its runs copied), cached."""
    b.uses((id(g), 'label_summary_prices'), 2 * len(g.labels))
    got = g.cache.get('label_summary_prices')
    if got is None:
        rows = derive.label_summary(g)
        w, m = [], []
        for per in rows:
            n = sum(len(v['runs']) for v in per.values())
            w.append(W_ELEM * (2 + len(per)) + (n >> 2))
            m.append(dict_bytes(2 + len(per)) + len(per) * (dict_bytes(4) + LIST)
                     + LIST_ITEM * n + 2 * STR + 24)
        got = g.cache['label_summary_prices'] = (w, m)
    return got


def _label_summary(g, b):
    if b is not None:
        _uses_g(b, g, 'label_summary')
        pw, pm = _label_summary_prices(b, g)
        b.phase = 'rows'
        # every row (the table has no partial answer), in blocks
        for i in range(0, len(pw), 1024):
            b.charge(sum(pw[i:i + 1024]), sum(pm[i:i + 1024]))
    rows = derive.label_summary(g)
    out = {}
    try:
        for lab, per in zip(g.labels, rows):
            row = lab.as_dict()
            for side, v in per.items():
                row[side] = dict(v, runs=list(v['runs']))
            out[lab.ref] = row
    except LocalBudgetExceeded as e:
        raise e.restate(labels=len(out), of=len(g.labels))
    return out


# ------------------------------------------------------------------ continuation

def continuation(g, arm, leaf, *, budget=None):
    """Continuation(sequence, labels, loss_used, branches_used, seed_coord, as_seed(),
    losses, note): the walker's label logic from C, the spelling derived (both arms).
    seed_coord is the half-open whole-molecule interval in seed coordinates (seed =
    [0, |seed|)). loss_used is C's: the SMALLEST terminal loss of the labels; losses
    holds each label's own (T), and note says when one loss budget cannot serve them
    exactly (labels that ended the walk at different losses). budget= (local limits): charged
    as one step at its price; a stop raises."""
    return _B.run_scoped(budget, 'continuation', ('select_walks',),
                         lambda b: _continuation(g, arm, leaf, b))


_CONT_BYTES = record_bytes(Continuation) + 2 * LIST + 200


def _continuation(g, arm, leaf, b):
    a = g.arm(arm)
    if b is not None:
        # admitted before the walk is resolved (its path id builds the arm's paths): a
        # refused call leaves no derived cache behind
        derive.uses(b, g, a, 'leaves', 'paths', 'end_labels')
        _uses_g(b, g, 'name_counts')
    seg_id = leaf_segment(a, leaf)
    seg = a.segments[seg_id]
    c = seg.leaf.continuation if seg.leaf else None
    if c is None:
        raise ValueError('walk %d of the %s arm has no continuation (it ended for a '
                         'semantic reason, not at the radius or a cap)'
                         % (path_id(a, leaf), a.side))
    if b is not None:
        n = seg.end_bp + len(g.seed.sequence)
        k = len(c.labels)
        b.charge(W_ROW + W_STEP * (seg.depth + 1) + (n >> 8) + 4 * k
                 + 4 * len(derive.runs_by_segment(a)[seg_id]),
                 _CONT_BYTES + 3 * (STR + n) + k * (LIST_ITEM + FLOAT_BYTES + 100))
    seq = derive.continuation_sequence(g, a, seg_id)
    L = seg.end_bp
    if a.side == 'right':
        coord = (g.seed.length_bp + L - c.n, g.seed.length_bp + L)
    else:
        coord = (-L, -L + c.n)
    labels = _labels(g, c.labels)
    losses = _terminal_losses(g, a, seg_id, labels)
    note = None
    budget = _strategy_budget(g)
    if losses and budget and len(set(losses)) > 1:
        top = max(losses)
        note = ('its labels ended the walk at different losses (%s): a continuation '
                'request carries one loss budget for all of them, so next_request() '
                'reduces labels.loss_budget %s by the largest, %s -- exact for the labels '
                'at that loss, conservative for the others (one uninterrupted walk would '
                'let each spend its own remaining %s)'
                % (_loss_list(labels, losses), _num(budget), _num(top),
                   ', '.join('%s %s' % (_shown(l.name), _num(max(0.0, budget - x)))
                             for l, x in zip(labels, losses) if x < top)))
    return Continuation(seq, labels, c.loss_used, c.branches_used, coord, a.side,
                        path_id(a, leaf), _unverifiable_names(g, labels), losses, note)


def _num(x):
    return '%g' % x


def _loss_list(labels, losses):
    return ', '.join('%s %s' % (_shown(l.name), _num(x)) for l, x in zip(labels, losses))


def _strategy_budget(g):
    """The retrieval's labels.loss_budget (constrain), None when there is none to spend."""
    if g.mode != 'constrain' or not g.envelope:
        return None
    b = ((g.envelope.get('strategy') or {}).get('labels') or {}).get('loss_budget')
    return b if isinstance(b, (int, float)) and not isinstance(b, bool) and b > 0 else None


def _terminal_losses(g, a, leaf_seg, labels):
    """Constrain: each label's terminal loss at the leaf, from its T extras (a label
    alive at the leaf with no extra has loss 0, the §2.3 rule). The C record's loss_used
    is the SMALLEST of them: subtracting it from the budget would let a continued label at
    a higher loss go on with more than its original budget left (E at 0 and C at 2 would
    keep the whole budget 3). () in annotate mode."""
    if g.mode != 'constrain':
        return ()
    ends = {e.label: e.loss for e in derive.end_labels(a, leaf_seg)}
    extras = a.segments[leaf_seg].leaf.extras
    out = []
    for l in labels:
        if l.id in ends:
            out.append(ends[l.id])
        else:
            # not alive at the leaf by its runs: the T extras still hold its loss when
            # the writer stated one; a label with neither is at loss 0 (§2.3)
            out.append(extras.get(l.id, (0.0, 0, 0))[0])
    return tuple(out)


def _is_number(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def _finite_nonneg(v):
    """A switch cost or loss budget the server accepts: a number >= 0 and finite (JSON
    carries no NaN or infinity, an int beyond a double is none either, and a NaN compares
    false with every budget, so a search over it would answer quietly wrong)."""
    if not _is_number(v):
        return False
    try:
        return math.isfinite(v) and v >= 0
    except OverflowError:
        return False


def _check_change_cost(change_cost):
    """Validates the fields of a change_cost the library prices -- the model's type, the
    entries' type (a list, under every model), a constant's value, a table's default and
    every entry of a table: [from, to, cost] with two label names (strings) and a finite
    cost >= 0, the server's rule (traverse.cpp parse_cost) --
    BEFORE any of them is used as a dict key, in a set or in arithmetic: otherwise an entry
    such as [["C"], "A", 0.5] would raise TypeError from a set membership test, which
    escapes the tool layer's bad_argument, and a malformed entry would be skipped silently
    where the server refuses the request. A field that fails is a ValueError naming it,
    and so is a field the named model does not read (forbid reads model only, constant
    model and value, table model, default and entries): the server's Strict parse refuses
    it as an unknown field. A model the
    library does not know passes: its callers answer for it. None (no change_cost) is
    forbid. -> change_cost."""
    if change_cost is None:
        return None
    if not isinstance(change_cost, dict):
        raise ValueError('labels.change_cost is an object, not %r' % (change_cost,))
    model = change_cost.get('model', 'forbid')
    if not isinstance(model, str):
        raise ValueError('labels.change_cost.model is "forbid", "constant" or "table", not %r'
                         % (model,))
    fields = _COST_FIELDS.get(model)
    if fields is not None:
        foreign = sorted(k for k in change_cost if k not in fields)
        if foreign:
            # a field of another model: what a deep merge of a model change kept (a
            # constant's value under a table, a table's entries under a constant), which
            # the server refuses: the library must not call the request valid
            raise ValueError('labels.change_cost.%s is not a field of model %r (it reads %s): '
                             'the server refuses it as an unknown field'
                             % (foreign[0], model, ', '.join(sorted(fields))))
    # entries under every model the library does not know too: a budgeted call charges
    # their count whatever the model (_switch_reach), and {'model': 'forbid', 'entries':
    # True} would raise TypeError from len() there; None reads as none
    entries = change_cost.get('entries')
    if entries is not None and not isinstance(entries, (list, tuple)):
        raise ValueError('labels.change_cost.entries is a list of [from, to, cost], not %r'
                         % (entries,))
    if model == 'constant':
        v = change_cost.get('value')
        if not _finite_nonneg(v):
            raise ValueError('labels.change_cost.value is a finite number >= 0, not %r' % (v,))
    elif model == 'table':
        d = change_cost.get('default', 'forbid')
        if not (isinstance(d, str) and d == 'forbid') and not _finite_nonneg(d):
            raise ValueError('labels.change_cost.default is "forbid" or a finite number >= 0, '
                             'not %r' % (d,))
        # required, as the server requires it: a table without entries is refused there
        if entries is None:
            raise ValueError('labels.change_cost.entries is a list of [from, to, cost], not %r'
                             % (entries,))
        for i, e in enumerate(entries):
            if not (isinstance(e, (list, tuple)) and len(e) == 3
                    and isinstance(e[0], str) and isinstance(e[1], str)):
                raise ValueError('labels.change_cost.entries[%d] is [from, to, cost] with two '
                                 'label names (strings), not %r' % (i, e))
            if not _finite_nonneg(e[2]):
                raise ValueError('labels.change_cost.entries[%d]: the cost of %r -> %r is a '
                                 'finite number >= 0, not %r' % (i, e[0], e[1], e[2]))
    return change_cost


def _check_annotate_cost(change_cost):
    """The merged labels.change_cost of an annotate-mode continuation, checked as the
    server reads it: labels are recorded, not filtered, so there is no permitted set and
    any model but "forbid" is refused there (walker.cpp: a finite change cost conflicts
    with labels.mode "annotate"), as is a field forbid does not read (the Strict parse).
    The model is checked first, so no entries of a table are read. None (absent) passes."""
    if change_cost is None:
        return
    model = change_cost.get('model', 'forbid') if isinstance(change_cost, dict) else None
    if isinstance(model, str) and model != 'forbid':
        raise ValueError('labels.change_cost.model %r does not apply in labels.mode "annotate" '
                         '(labels are recorded, not filtered: there is no permitted set); the '
                         'server refuses any model but "forbid": omit the override' % model)
    _check_change_cost(change_cost)


# the fields each change_cost model reads (traverse.cpp parse_cost: a Strict object whose
# unread members are refused)
_COST_FIELDS = {'forbid': frozenset(('model',)), 'constant': frozenset(('model', 'value')),
                'table': frozenset(('model', 'default', 'entries'))}


def _section(strategy, *path):
    """The strategy section at |path| as a dict ({} when absent). The rebuild of a
    continuation reads the merged sections as objects; one that a caller's override made
    something else is refused with a ValueError naming it (the tool layer answers
    bad_argument), as the server would refuse it (never an AttributeError, which would
    escape the tool contract)."""
    d = strategy
    for i, k in enumerate(path):
        v = d.get(k)
        if v is None:
            return {}
        if not isinstance(v, dict):
            raise ValueError('strategy.%s is an object, not %r' % ('.'.join(path[:i + 1]), v))
        d = v
    return d


# What a change_cost table costs a budgeted next_request() (work model 2: the entries'
# copies are charged, or a table of 900 entries would leave the account at the 0-entry
# figure while the traced peak is 9.8x it). Sizes as CPython makes them:
# an entry copied by copy.deepcopy() -- its list of three, built by appends (four slots),
# and its slot in the copied list: 97 B traced per entry -- and the memo's record of it
# while the copy is built (an int key, its dict slot, both of the dict's tables while it
# grows, a keep-alive slot: up to 128 B traced), dropped when deepcopy() returns
_ENTRY_COPY = LIST + 4 * LIST_ITEM + 2 * LIST_ITEM
_ENTRY_MEMO = INT + 6 * DICT_KEY + 2 * LIST_ITEM
_W_ENTRY_COPY = 3 * W_ELEM          # measured 0.6 us per entry (3.14, M5 Max)
# _switch_reach's table: a (from, to) key tuple, its float cost and its dict slot (a dict
# keeps a third of its slots free); an explicit edge (to, cost) in its source's list; a heap
# entry (loss, order, name) with its new float and int
_TABLE_PAIR = TUPLE + 2 * 8 + FLOAT_BYTES + 3 * DICT_KEY
_EDGE = TUPLE + 2 * 8 + LIST_ITEM
_HEAP_ITEM = TUPLE + 3 * 8 + FLOAT_BYTES + INT + LIST_ITEM
# the tuples of one length CPython keeps on its free list once they are freed
# (PyTuple_MAXFREELIST); a pair is TUPLE + 2 * 8 bytes, 64 on 3.14 with its cached hash
_FREE_TUPLES = 2000


def _entries_of(change_cost):
    """The entries of a change_cost section as given, () when it has none or is malformed
    (what _check_change_cost() refuses is refused there, after the charge)."""
    if isinstance(change_cost, dict):
        e = change_cost.get('entries')
        if isinstance(e, (list, tuple)):
            return e
    return ()


def _override_entries(strategy):
    """The labels.change_cost entries a strategy, or the keyword overrides of
    next_request(), carry; () where there are none or a section is not an object."""
    lab = strategy.get('labels') if isinstance(strategy, dict) else None
    return _entries_of(lab.get('change_cost')) if isinstance(lab, dict) else ()


def _charge_entry_copy(lb, k):
    """Charge a deep copy of |k| change_cost entries before it is made: its work, the
    copy (held by the caller until _release_entry_copy()) and the memo deepcopy() drops on
    return."""
    if k:
        lb.charge(_W_ENTRY_COPY * k, list_bytes(k) + k * (_ENTRY_COPY + _ENTRY_MEMO))
        lb.release(k * _ENTRY_MEMO)


def _release_entry_copy(lb, k):
    if k:
        lb.release(list_bytes(k) + k * _ENTRY_COPY)


def _switch_reach(change_cost, sources, targets, budget, sinks=(), lb=None):
    """{name: loss} for every name of |targets| that a chain of switches from a name of
    |sources| enters within |budget|, at the cheapest such chain's loss (summed left to
    right, as the walk sums a lineage's loss) -- the server's rule for an extra label
    (walker.cpp, switch_reach): the walk enforces the cumulative loss switch by switch, so a
    label is a valid switch target when SOME chain reaches it, not only one switch from a
    seed label (A -> B = 1, B -> C = 1 under a budget of 2 reaches C). A chain passes only
    through names of |sources| and |targets|: the labels the request names; a name of
    |sinks| (also a target) ends a chain but never continues one. Prices a switch as the
    server does (forbid: none reachable; constant: its value; table: a table's last entry
    for a pair wins, else its default, 'forbid' none; a label to itself costs 0; the one
    pricing of the library); None for a model the library does not know. A malformed field
    of a model it knows is a ValueError naming it, raised before any entry is used
    (_check_change_cost)."""
    if lb is not None:
        # admitted before the check reads the entries
        n = len(sources) + len(targets)
        k = len(_entries_of(change_cost))
        lb.charge(4 * n + 2 * W_ELEM * k, 4 * set_bytes(n) + dict_bytes(n) + n * FLOAT_BYTES)
    # every field is checked before an endpoint enters a set or a cost the arithmetic
    _check_change_cost(change_cost)
    sources = list(dict.fromkeys(sources))
    # the set once (rebuilt per target it would cost sources x targets steps uncharged --
    # 0.11 s of 0.14 s at 3,000 x 3,000 names)
    src = set(sources)
    targets = [t for t in dict.fromkeys(targets) if t not in src]
    model = (change_cost or {}).get('model', 'forbid')
    if not sources or model == 'forbid':
        return {}
    if model == 'constant':
        # one switch reaches every label, and a chain only costs more
        c = float(change_cost['value'])
        return {t: c for t in targets if c <= budget}
    if model != 'table':
        return None
    names = set(sources) | set(targets)
    if lb is not None:
        # the table before it is built: a pair per entry between two of the names (counted
        # first, so that a long table of other labels' pairs is not charged as held), at
        # most one per two names; what duplicates leave unbuilt is given back after
        lb.charge(W_STEP * k)
        bound = min(sum(1 for e in change_cost.get('entries') or ()
                        if e[0] in names and e[1] in names), len(names) ** 2)
        lb.charge(0, dict_bytes(bound) + bound * _TABLE_PAIR)
    table = {}
    for e in change_cost.get('entries') or ():
        if e[0] in names and e[1] in names:
            table[(e[0], e[1])] = float(e[2])     # the last entry for a pair wins
    m = len(table)
    if lb is not None:
        lb.release((bound - m) * (DICT_KEY + _TABLE_PAIR))
        # the explicit edges, a list per source name
        lb.charge(W_ELEM * m, dict_bytes(len(names)) + len(names) * LIST + m * _EDGE
                  + len(sources) * _HEAP_ITEM)
    d = change_cost.get('default', 'forbid')
    fallback = math.inf if d == 'forbid' else float(d)
    out_edges = collections.defaultdict(list)
    for (u, v), c in table.items():
        if u != v and c != math.inf:
            out_edges[u].append((v, c))
    # the default is relaxed lazily, as on the server: the names no popped name has given
    # the default yet; a popped name gives it to each of them it has no explicit entry to.
    # A list in name order, compacted per pop as walker.cpp's `kept` loop does: a set's
    # order followed PYTHONHASHSEED, and with it the push order of labels at equal loss,
    # the pop sequence and so the charges below (a set's order depends on the hash seed;
    # budget.py promises the same units in any process)
    pending = sorted(names) if fallback <= budget else []
    dist = {x: 0.0 for x in sources}
    heap = [(0.0, i, x) for i, x in enumerate(sorted(sources))]
    heapq.heapify(heap)
    order = len(heap)
    done = set()

    def relax(v, loss):
        nonlocal order
        if loss <= budget and loss < dist.get(v, math.inf):
            dist[v] = loss
            heapq.heappush(heap, (loss, order, v))
            order += 1

    while heap:
        if lb is not None:
            # a pop, its explicit edges and the names still owed the default (their copy
            # is dropped before the next pop), and the heap entries this pop may push: what
            # it did not push is given back after it, with the entry it popped
            k_pop = len(out_edges.get(heap[0][2], ())) + len(pending)
            lb.charge(4 + k_pop, list_bytes(len(pending)) + k_pop * _HEAP_ITEM)
            lb.release(list_bytes(len(pending)))
            pushed = order
        loss, _, u = heapq.heappop(heap)
        if not (u in done or loss > dist[u]):
            done.add(u)
            if u not in sinks:
                for v, c in out_edges.get(u, ()):
                    relax(v, loss + c)
                if pending:
                    kept = []
                    for v in pending:
                        if v == u or (u, v) in table:
                            kept.append(v)  # no default from u to v: a later pop gives it
                        else:
                            relax(v, loss + fallback)
                    pending = kept
        if lb is not None:
            lb.release((k_pop - (order - pushed) + 1) * _HEAP_ITEM)
    if lb is not None:
        # the table and its edges are dropped on return -- but CPython keeps up to 2,000
        # freed tuples of each length for reuse (its free list), so the keys' and edges'
        # pairs stay allocated for the rest of the call (traced: 118 KB after a 900-entry
        # table): kept in the account once per call, since a later search reuses them
        held = dict_bytes(m) + m * _TABLE_PAIR + dict_bytes(len(names)) \
            + len(names) * LIST + m * _EDGE
        if lb.first('switch_reach_free_tuples'):
            held -= min(_FREE_TUPLES, 2 * m) * (TUPLE + 3 * 8)
        lb.release(held)
    return {t: dist[t] for t in targets if dist.get(t, math.inf) <= budget}


def _rebuilt_extra(g, seed_labels, budget, strategy, lb=None):
    """labels.extra around a continuation's seed labels: the retrieval's permitted pool
    (its seed labels and its extra labels: every label of a constrain retrieval) minus the
    new seed labels -- the server refuses an extra label that duplicates a seed label
    (seed ['C'] with extra ['B', 'C', 'D'] is a 400) -- keeping each
    remaining label that the server accepts as a switch target: some chain of switches
    from a new seed label enters it within |budget| (_switch_reach; an unreachable extra
    label is refused, and under forbid, or a constant above the budget, any extra label
    is). -> (extra names in dictionary order, the labels left out as unreachable, the
    labels left out because their names cannot be verified). The retrieval's original seed
    labels become switch targets too: they were in the original pool."""
    lab = strategy.get('labels') or {}
    if lb is not None:
        # the check below reads every entry
        lb.charge(W_ELEM * len(_entries_of(lab.get('change_cost') if isinstance(lab, dict)
                                           else None)))
    cost = _check_change_cost(lab.get('change_cost')) or {'model': 'forbid'}
    seed_names = [l.name for l in seed_labels]
    seed_ids = {l.id for l in seed_labels}
    seeds = set(seed_names)
    if lb is not None:
        _uses_g(lb, g, 'name_counts')
        n = len(g.labels)
        lb.charge(8 * n, 4 * list_bytes(n) + 2 * set_bytes(n))
    # /traverse resolves an extra label by NAME too: one whose name the library cannot verify
    # to be its own is never sent (left out and stated if a chain reaches it, unreachable
    # otherwise), so no chain is counted through it -- the request does not name it
    suspect = {l.id for l in g.labels if l.id not in seed_ids and l.name not in seeds
               and _unverifiable_names(g, [l])}
    others = [l for l in g.labels if l.id not in seed_ids and l.name not in seeds]
    reach = _switch_reach(cost, seed_names, [l.name for l in others], budget,
                          sinks={l.name for l in others if l.id in suspect}, lb=lb)
    if reach is None:
        raise ValueError('change_cost model %r is not one the library can price: the '
                         'continuation\'s labels.extra cannot be validated'
                         % cost.get('model'))
    # a chain within the budget enters each of its labels within it, so the reachable labels
    # are closed: leaving out those no chain reaches changes no chain of those kept
    extra, dropped, unverifiable = [], [], []
    for l in g.labels:
        if l.id in seed_ids:
            continue
        if l.name in seeds:
            # another label under a seed label's name: it would resolve to the seed label
            unverifiable.append(l)
        elif l.name not in reach:
            dropped.append(l)
        elif l.id in suspect:
            unverifiable.append(l)
        else:
            extra.append(l.name)
    if cost.get('model', 'forbid') == 'forbid' or (cost.get('model') == 'constant'
                                                  and cost.get('value', 0) > budget):
        # the request-level rule: no extra label at all (each is unreachable)
        dropped = [l for l in g.labels if l.id not in seed_ids]
        extra, unverifiable = [], []
    return extra, dropped, unverifiable


# the left-out labels the note of next_request() names (the rest are counted)
_LEFT_OUT_SHOWN = 8


def _left_out_reach(g, cost, mine, budget, dropped, lb=None):
    """For the left-out note of next_request(): which continued labels |mine| ((label,
    terminal loss) pairs) reach each of the first _LEFT_OUT_SHOWN labels of |dropped| by a
    chain of switches through the retrieval's pool within their own remaining budget
    (budget - loss), and whether any label of |dropped| is so reached. -> ([the reaching
    labels of each shown one, in |mine| order], any).

    Exactly the note of one _switch_reach() per continued label over the whole pool (S
    searches of the pool, all S results held at once, a membership test per dropped label
    and continued label: 1.26 s for one walk at 1,000 labels), with the searches confined to
    the names the change_cost table names (its pairs between labels of the pool) plus ONE
    name it does not: a name no entry names is reached from the source by the default switch
    alone -- the first pop (the source, at loss 0) gives it 0 + default, and every other
    chain into it ends with a default switch too, so costs at least as much -- and as an
    intermediate it gives every other name the default, as each such name does; so one
    stands for all, and the names the table names get the same smallest chain sums, added in
    the same order. Under a constant model no name is named: each search is one switch."""
    model = (cost or {}).get('model', 'forbid')
    pool = list(dict.fromkeys(l.name for l in g.labels))
    named = set()
    if model == 'table':
        names = set(pool)
        for e in cost.get('entries') or ():
            if e[0] in names and e[1] in names:
                named.add(e[0])
                named.add(e[1])
    if lb is not None:
        # the names read once and the reach tests below: a reach per continued label for
        # each label shown and each left-out label the table names, one per left-out label
        n_named_dropped = sum(1 for l in dropped if l.name in named)
        lb.charge(W_ELEM * (len(pool) + len(dropped) + len(mine)
                            * (min(len(dropped), _LEFT_OUT_SHOWN) + n_named_dropped)),
                  set_bytes(len(named)) + 2 * list_bytes(len(pool)) + dict_bytes(len(pool))
                  + len(mine) * TUPLE)
    in_table = [x for x in pool if x in named]
    # the first two names no entry names (a stand-in for each source: not the source)
    others = list(itertools.islice((x for x in pool if x not in named), 2))
    reach_of = []
    for p, x in mine:
        stand_in = next((u for u in others if u != p.name), None)
        targets = [t for t in in_table if t != p.name]
        if stand_in is not None:
            targets.append(stand_in)
        r = _switch_reach(cost, [p.name], targets, budget - x, lb=lb) or {}
        reach_of.append((p, r, stand_in is not None and stand_in in r))

    def reaches(name):
        if name in named:
            return [p for p, r, _ in reach_of if p.name != name and name in r]
        return [p for p, _, rest in reach_of if p.name != name and rest]

    shown = [reaches(l.name) for l in dropped[:_LEFT_OUT_SHOWN]]
    # a boolean, not a substring test of the note's text: a label NAME holding "reachable
    # for" claimed a reach that no chain makes
    any_reach = any(shown)
    if not any_reach:
        for l in dropped[_LEFT_OUT_SHOWN:]:
            if l.name in named and reaches(l.name):
                any_reach = True
                break
        else:
            # the names no entry names: reached alike, by the continued labels that reach
            # the stand-in, each one except under its own name
            rest = {p.name for p, _, ok in reach_of if ok}
            any_reach = any(l.name not in named and (len(rest) > 1
                                                     or (rest and l.name not in rest))
                            for l in dropped[_LEFT_OUT_SHOWN:])
    if lb is not None:
        # the note's text of the shown labels: the names of the labels that reach them, each
        # written by _shown() and joined, the entries joined again into the note
        k = sum(len(v) for v in shown)
        n = sum(_shown_bytes(p.name) + 2 for v in shown for p in v)
        lb.charge(W_ELEM * k + (n >> 4), 3 * (n + STR * _LEFT_OUT_SHOWN) + k * STR)
    return shown, any_reach


def _continuation_seeds(g, a, leaves, lb=None):
    conts = [_continuation(g, a, x, lb) for x in leaves]
    seeds = []
    for c in conts:
        if not c.sequence:
            raise ValueError('walk %d has no continuation sequence (output.continuation_bp '
                             'was 0): nothing to resubmit' % c.leaf)
        if g.mode != 'constrain':
            seeds.append({'sequence': c.sequence})   # annotate mode has no permitted set
        else:
            seeds.append(c.as_seed())    # raises UnverifiableLabelName, never resolves
    return conts, seeds


# budget= of next_request() and next_requests() when the caller gives none
_NO_BUDGET_ARG = object()


def _local_or_override(budget, overrides):
    """The LocalBudget a next_request() runs under (or None). `budget` is also a keyword
    override like any other, deep-merged into the strategy under that key (the tools pass
    an agent's overrides on as given): a LocalBudget charges the call, any other value
    (None included) is that override; a LocalLimits is refused as everywhere
    else (it is meant as a local limit, never a strategy field)."""
    if budget is _NO_BUDGET_ARG:
        return _B.resolve(None)
    if isinstance(budget, _B.LocalBudget):
        return budget
    if isinstance(budget, _B.LocalLimits):
        _B.resolve(budget)                    # the TypeError of every other operation
    overrides['budget'] = budget
    return _B.resolve(None)


def next_request(g, arm, leaves, bp=None, reduce_budget=True, reset_branches=False, *,
                 budget=_NO_BUDGET_ARG, **overrides):
    """A resubmittable /traverse request continuing the given walks, as a NextRequest (a
    dict; .notes, .loss_budget and .branch_budget state what the request format cannot
    carry): one seed per walk (its continuation, natural orientation, with its labels named
    explicitly), the retrieval's normalized strategy with direction = this arm and
    max_extension_bp = bp when given.

    Constrain mode, three rules no continuation may break:
      * no route exceeds its original loss budget: with reduce_budget, labels.loss_budget
        is reduced by the LARGEST terminal loss of the continued labels over all the
        walks (T, not C's loss_used, which is the smallest). The request carries one
        budget and no per-label starting loss, so it is exact for the labels at that loss
        and conservative for the others: their continuation may stop earlier than one
        uninterrupted walk would -- stated in .notes, with each label's remaining budget
        in .loss_budget. reduce_budget=False keeps the original budget: every label
        restarts at loss 0 and may spend it again (stated);
      * no lineage branches beyond its original allowance: branching.max_label_branches
        is reduced the same way, by the LARGEST terminal branch count of the continued
        labels ("unlimited" stays unlimited), exact for the labels at that count and
        conservative for the others -- stated in .notes and .branch_budget ({original,
        effective, largest_terminal_branches, reset}). reset_branches=True keeps the
        original allowance: every lineage restarts at 0 branches and may branch again
        (stated);
      * the request is valid: labels.extra is rebuilt around the new seed labels -- the
        retrieval's permitted pool minus them, each label kept that a chain of switches
        from a seed label enters within the budget (the server's rule); every other label
        is listed in .left_out, per walk (and the first few named in .notes). A
        labels.extra the caller gives in the overrides is sent as given. An override
        section the rebuild reads that is not an object (labels, labels.change_cost,
        branching -- None included), or a malformed field of it, is a ValueError naming
        it; so is None for any section the strategy holds as an object, which the server
        refuses as null. An override's labels.change_cost that names another model
        replaces the retrieval's whole; one that names the same model (or none) updates
        its fields; a field the resulting model does not read is a ValueError naming it
        (the server's Strict parse refuses it). Other override fields are sent as given:
        the server answers for them.
    Annotate mode (labels recorded, not filtered) rebuilds nothing; the overrides merge the
    same way, and the resulting labels.change_cost is checked too: any model but "forbid"
    is a ValueError naming it (the server refuses a change cost in that mode), and so is a
    field forbid does not read.
    A label alive at a leaf that does not cover the continuation's whole tail is not among
    the continuation's labels, so it is not seeded and the continuation may lack its
    lineage: stated in .notes and in .left_out (why 'alive_not_seeded').
    Several walks share one strategy, so one labels.extra: they are built into one request
    only when the rebuilt pool is the same for every walk (their continuations carry the
    same labels, or no switch is possible); otherwise IncompatibleContinuations names
    next_requests(), one request per walk. Keyword overrides deep-merge into the
    strategy (labels.change_cost as described above); release, graph and graph_path are
    request-level. Raises MissingEnvelope on
    a body-only graphlet. budget= (local limits) a LocalBudget, not the request's loss
    budget: the continuations, the rebuilt labels.extra and the switch searches are charged;
    a stop raises, with no partial. Any other value of budget= is a keyword override like
    the others (strategy['budget']); a local budget then comes from local_budget() only."""
    lb = _local_or_override(budget, overrides)
    if lb is None:
        return _next_request(g, arm, leaves, bp, reduce_budget, reset_branches, None,
                             overrides)
    with lb.scope('next_request', ('select_walks',)):
        return _next_request(g, arm, leaves, bp, reduce_budget, reset_branches, lb,
                             overrides)


def _next_request(g, arm, leaves, bp, reduce_budget, reset_branches, lb, overrides):
    if not isinstance(reset_branches, bool):
        raise ValueError('reset_branches is true or false, not %r' % (reset_branches,))
    g.require_envelope('next_request()')
    a = g.arm(arm)
    if isinstance(leaves, (int, Path, Walk)):
        leaves = [leaves]
    # (three arguments without a budget: the tests replace this helper)
    conts, seeds = _continuation_seeds(g, a, leaves) if lb is None \
        else _continuation_seeds(g, a, leaves, lb)
    # A change_cost table is the one part of the strategy whose size the caller sets: its
    # entries are copied with the strategy (the request's copy, and the merged one the
    # labels.extra rebuild reads), so each copy is charged -- uncharged, 900 entries would
    # leave the account at the 0-entry figure while the traced peak is 9.8x it. The rest of
    # the strategy is a fixed handful of fields
    k_env = len(_override_entries(g.envelope.get('strategy'))) if lb is not None else 0
    k_ov = len(_override_entries(overrides)) if lb is not None else 0
    if lb is not None:
        _charge_entry_copy(lb, k_env)
    strategy = copy.deepcopy(g.envelope.get('strategy') or {})
    strategy.pop('clamped', None)
    strategy['direction'] = a.side
    if bp is not None:
        strategy.setdefault('bounds', {})['max_extension_bp'] = bp
    notes = []
    budget_info = None
    branch_info = None
    left_out = []
    if g.mode == 'constrain':
        labels_ = strategy.setdefault('labels', {})
        budget = labels_.get('loss_budget')
        budget = budget if isinstance(budget, (int, float)) and \
            not isinstance(budget, bool) else 0.0
        pairs = [(l, x) for c in conts for l, x in zip(c.labels, c.losses)
                 if x != math.inf]
        top = max((x for _, x in pairs), default=0.0)
        effective = budget
        if reduce_budget and budget > 0 and top > 0:
            effective = max(0.0, budget - top)
            labels_['loss_budget'] = effective
            lower = [(l, x) for l, x in pairs if x < top]
            if lower:
                notes.append(
                    'labels.loss_budget %s is the original %s minus the largest terminal loss '
                    'of the continued labels, %s (%s): exact for the labels at that loss, '
                    'conservative for %s, whose continuation may stop (loss_budget) earlier '
                    'than one uninterrupted walk would -- the request carries one loss '
                    'budget and no per-label starting loss'
                    % (_num(effective), _num(budget), _num(top),
                       ', '.join(_shown(l.name) for l, x in pairs if x == top),
                       ', '.join('%s (loss %s, could spend %s)'
                                 % (_shown(l.name), _num(x), _num(budget - x))
                                 for l, x in _unique_pairs(lower))))
        elif not reduce_budget and budget > 0 and top > 0:
            notes.append('labels.loss_budget is the original %s although the continued '
                         'labels have spent up to %s (reduce_budget=False): each restarts '
                         'at loss 0, so a route of the continuation may exceed its original '
                         'budget by up to %s' % (_num(budget), _num(top), _num(top)))
        budget_info = {
            'original': budget, 'effective': effective, 'largest_terminal_loss': top,
            'labels': [dict(l.as_dict(), loss=x, remaining=max(0.0, budget - x))
                       for l, x in _unique_pairs(pairs)]}
        # the caller's overrides are final: labels.extra is rebuilt against the budget
        # and the cost model the request will carry (a list the caller gives is theirs)
        if lb is not None:
            # both copies are held until the call returns (an override's entries replace
            # the retrieval's in |final| only once they are copied)
            _charge_entry_copy(lb, k_env)
            _charge_entry_copy(lb, k_ov)
        final = copy.deepcopy(strategy)
        _merge_overrides(final, {k: v for k, v in overrides.items()
                                 if k not in ('release', 'graph', 'graph_path')})
        flab = _section(final, 'labels')
        if lb is not None:
            lb.charge(W_ELEM * len(_entries_of(flab.get('change_cost'))))
        # the merged cost the request will carry, checked whole even when the caller gives
        # labels.extra (no switch search reads it then, and the server would refuse it)
        _check_change_cost(_section(final, 'labels', 'change_cost') or None)
        _section(final, 'branching')
        given = overrides.get('labels') or {}
        if 'loss_budget' in given and not _finite_nonneg(given['loss_budget']):
            # a NaN budget would compare false with every loss and leave out every label
            raise ValueError('labels.loss_budget is a finite number >= 0, not %r'
                             % (given['loss_budget'],))
        if given.get('extra') is not None and not (
                isinstance(given['extra'], list) and all(isinstance(x, str) for x in given['extra'])):
            raise ValueError('labels.extra is a list of label names, not %r' % (given['extra'],))
        if 'loss_budget' in given:
            effective = flab['loss_budget']
            budget_info['effective'] = effective
            if effective > budget - top:
                notes.append('labels.loss_budget %s is the caller\'s (overrides): above the '
                             'original %s minus the largest terminal loss %s, a route of the '
                             'continuation may exceed its original budget'
                             % (_num(effective), _num(budget), _num(top)))
        if 'extra' in given:
            # the caller's own list: the request carries it as given, unchecked
            extra = list(given['extra'] or [])
            per_walk = []
        else:
            pools = [_rebuilt_extra(g, c.labels, effective, final, lb) for c in conts]
            if len({(tuple(p[0]), tuple(l.id for l in p[2])) for p in pools}) > 1:
                raise IncompatibleContinuations(
                    'the walks %s continue under different labels, so the switch targets '
                    'each seed needs (labels.extra) differ and one request cannot carry them '
                    'without duplicating a seed label or dropping a switch target: build one '
                    'request per walk with next_requests()'
                    % ', '.join(str(c.leaf) for c in conts))
            extra = pools[0][0] if pools else []
            # what each walk leaves out, from ITS seed labels: the walks share labels.extra,
            # not their seed labels, so a label one walk leaves out may be another's seed
            # label (the first walk's omissions reported for all would name a label another
            # walk seeded as unreachable)
            per_walk = [(c, p[1], p[2]) for c, p in zip(conts, pools)]
        labels_['extra'] = extra
        for c, dropped, unverifiable in per_walk:
            left_out += [dict(l.as_dict(), why='unreachable', walk=c.leaf) for l in dropped] + \
                [dict(l.as_dict(), why='unverifiable_name', walk=c.leaf) for l in unverifiable]
        unverifiable = per_walk[0][2] if per_walk else []      # the same for every walk
        # A label alive at the leaf that does not cover the whole tail (it switched in
        # within the continuation's last bases) is not among C's labels, so it is no seed
        # label of the continuation: it can come back only as a switch target from a seed
        # label (under switch_on: loss, only where one is lost), and the continuation may
        # lack the lineage one uninterrupted walk keeps. Conservative -- no route above the
        # budget is accepted -- but never silent.
        alive_not_seeded = []
        for c in conts:
            seeded = {l.id for l in c.labels}
            for e in derive.end_labels(a, leaf_segment(a, c.leaf)):
                if e.label not in seeded:
                    alive_not_seeded.append((c.leaf, g.labels[e.label], e.loss))
        if alive_not_seeded:
            left_out += [dict(l.as_dict(), why='alive_not_seeded', walk=w, loss=x)
                         for w, l, x in alive_not_seeded]
            shown = ['%s (walk %d, loss %s)' % (_shown(l.name), w, _num(x))
                     for w, l, x in alive_not_seeded]
            notes.append(
                '%s %s alive at the leaf but not among the continuation\'s labels (%s not '
                'cover its whole tail): not seeded, %s can be entered again only by a switch '
                'from a seed label (under switch_on: loss, only where one is lost), so the '
                'continuation may lack %s lineage that one uninterrupted walk keeps'
                % (', '.join(shown[:8]) + (' and %d more' % (len(shown) - 8)
                                           if len(shown) > 8 else ''),
                   'is' if len(shown) == 1 else 'are',
                   'it does' if len(shown) == 1 else 'they do',
                   'it' if len(shown) == 1 else 'they',
                   'its' if len(shown) == 1 else 'their'))
        if unverifiable:
            notes.append('labels.extra leaves out %s, a switch target of the retrieval whose '
                         'name the library cannot verify to resolve back to it (%s): the '
                         'continuation cannot switch to it'
                         % (', '.join(l.ref for l in unverifiable),
                            _unverifiable_names(g, unverifiable)))
        cost = _section(final, 'labels', 'change_cost') or {'model': 'forbid'}
        for c, dropped, _ in per_walk:
            if not dropped or cost.get('model', 'forbid') == 'forbid':
                # under forbid nothing was a switch target in the retrieval either
                continue
            # a label left out is a target one uninterrupted walk may still have entered from
            # a continued label at a lower loss, whose own remaining budget is larger -- by a
            # chain through ANY label of the retrieval's pool, the kept extra labels and the
            # other continued labels included, as that walk could switch through every one
            # (otherwise a chain through a kept extra label would not be searched, and the
            # note would understate what the continuation may miss)
            mine = _unique_pairs([(l, x) for l, x in zip(c.labels, c.losses) if x != math.inf])
            via, any_reach = _left_out_reach(g, cost, mine, budget, dropped, lb)
            lost = ['%s%s' % (_shown(l.name), ' (reachable for %s within its own remaining '
                                              'budget)' % ', '.join(_shown(p.name) for p in v)
                              if v else '')
                    for l, v in zip(dropped, via)]
            notes.append(
                'labels.extra leaves out %s%s: no chain of switches from a label of the '
                'continuation enters %s within loss_budget %s, and the server refuses an '
                'unreachable extra label%s'
                % (', '.join(lost) + (' and %d more' % (len(dropped) - 8)
                                      if len(dropped) > 8 else ''),
                   ' for walk %d' % c.leaf if len(per_walk) > 1 else '',
                   'it' if len(dropped) == 1 else 'them', _num(effective),
                   '; one uninterrupted walk could still have entered %s from a label at a '
                   'lower loss' % ('it' if len(dropped) == 1 else 'them')
                   if any_reach else ''))
        # Each continued label's own branches at the leaf (T), not C's smallest. The server
        # resets branch state for a new seed, so the allowance is reduced as the loss budget
        # is: by the largest count, conservative for the labels that used fewer;
        # reset_branches keeps it, and the note says so
        branches = 0
        per_label = []
        for c in conts:
            ends = {e.label: e.branches for e in derive.end_labels(a, leaf_segment(a, c.leaf))}
            per_label += [(l, ends.get(l.id, 0)) for l in c.labels]
            branches = max([branches, c.branches_used]
                           + [ends.get(l.id, 0) for l in c.labels])
        original = _section(strategy, 'branching').get('max_label_branches')
        given_b = overrides.get('branching') or {}
        mlb = _section(final, 'branching').get('max_label_branches')
        # integral without int(): int(inf) raises OverflowError, which escapes the tool
        # layer's bad_argument, and int(nan) a ValueError that names no field (Python's
        # json reads Infinity and 1e999); is_integer() is False for both
        if 'max_label_branches' in given_b and not (
                mlb == 'unlimited' or (_is_number(mlb) and mlb >= 0 and (
                    isinstance(mlb, int) or mlb.is_integer()))):
            raise ValueError('branching.max_label_branches is a non-negative integer or '
                             '"unlimited", not %r' % (mlb,))
        if original == 'unlimited' or (isinstance(original, int)
                                       and not isinstance(original, bool) and original > 0):
            reduced = original if original == 'unlimited' or reset_branches \
                else max(0, original - branches)
            strategy.setdefault('branching', {})['max_label_branches'] = reduced
            branch_info = {'original': original, 'effective': reduced,
                           'largest_terminal_branches': branches, 'reset': reset_branches}
            if 'max_label_branches' in given_b:
                branch_info['effective'] = mlb
                if mlb == 'unlimited' and original != 'unlimited' or \
                        _is_number(mlb) and original != 'unlimited' and mlb > reduced:
                    notes.append('branching.max_label_branches %s is the caller\'s (overrides): '
                                 'above the original %d minus the largest terminal branch count '
                                 '%d, a lineage of the continuation may branch further than one '
                                 'uninterrupted walk would' % (mlb, original, branches))
            elif original != 'unlimited' and branches:
                # as the loss budget: exact for the lineages at the largest count (stated in
                # branch_budget), and a note where it is conservative for others
                fewer = [(l, b) for l, b in per_label if b < branches]
                if reset_branches:
                    notes.append('branching.max_label_branches is the original %d although the '
                                 'continued lineages had used up to %d (reset_branches=True): '
                                 'each restarts at 0 branches, so the continuation may branch '
                                 'further than one uninterrupted walk would' % (original, branches))
                elif fewer:
                    seen, shown = set(), []
                    for l, b in fewer:
                        if l.id not in seen:
                            seen.add(l.id)
                            shown.append('%s (%d, could use %d)'
                                         % (_shown(l.name), b, max(0, original - b)))
                    notes.append('branching.max_label_branches %d is the original %d minus the '
                                 'largest terminal branch count of the continued labels, %d: '
                                 'exact for the lineages at that count, conservative for %s, '
                                 'which may stop branching earlier than one uninterrupted walk '
                                 'would -- the request carries one branch allowance and no '
                                 'per-label starting count'
                                 % (reduced, original, branches,
                                    ', '.join(shown[:8]) + (' and %d more' % (len(shown) - 8)
                                                            if len(shown) > 8 else '')))
        if not (budget_info['original'] > 0 or budget_info['effective'] > 0):
            # no loss budget to spend (forbid, or a budget of 0): each label's remaining
            # budget would be 0 and say nothing, and a long list of them would cost a
            # receipt its room
            budget_info = None
    request = NextRequest({'seeds': seeds}, notes, budget_info, left_out, branch_info)
    release = g.envelope.get('release')
    if release:
        request['release'] = release
    for k in ('release', 'graph', 'graph_path'):
        if k in overrides:
            request[k] = overrides.pop(k)
    if lb is not None:
        _charge_entry_copy(lb, k_ov)          # the request's own copy of the override's
    # the same merge as |final|'s (a model change replaces change_cost; None for an object
    # section is refused), in annotate mode too
    _merge_overrides(strategy, overrides)
    if g.mode == 'annotate':
        # constrain mode checks the merged cost before its rebuild reads it (|final|);
        # annotate mode has no rebuild, so its request is checked here: unchecked, the
        # server would refuse it
        labs = strategy.get('labels')
        _check_annotate_cost(labs.get('change_cost') if isinstance(labs, dict) else None)
    # the echo carries the retrieval's cap: an override output.coordinates false (the
    # resource stop's drop_coordinates) would otherwise make every continuation, deepen()
    # and traverse_continue a 400
    _C.strip_coordinate_cap(strategy)
    request['strategy'] = strategy
    if lb is not None and g.mode == 'constrain':
        # |final| is dropped on return; the request's copies stay with the request
        _release_entry_copy(lb, k_env + k_ov)
    return request


def next_requests(g, arm, leaves, bp=None, reduce_budget=True, reset_branches=False, *,
                  budget=_NO_BUDGET_ARG, **overrides):
    """One next_request() per walk, in |leaves| order: each with its own seed labels, its
    own rebuilt labels.extra and its budgets reduced by ITS labels' largest terminal loss
    and branch count (never less conservative than one request over all the walks).
    budget= (local limits): one LocalBudget for all of them; a stop raises, with no partial.
    Any other value of budget= is a keyword override, as in next_request()."""
    if isinstance(leaves, (int, Path, Walk)):
        leaves = [leaves]
    lb = _local_or_override(budget, overrides)
    if lb is None:
        return [_next_request(g, arm, [x], bp, reduce_budget, reset_branches, None,
                              copy.deepcopy(overrides))
                for x in leaves]
    out = []
    k_ov = len(_override_entries(overrides))
    with lb.scope('next_requests', ('select_walks',)):
        for x in leaves:
            with lb.scope('next_request', ('select_walks',)):
                # each walk's own copy of the overrides, dropped once its request is built
                _charge_entry_copy(lb, k_ov)
                out.append(_next_request(g, arm, [x], bp, reduce_budget, reset_branches,
                                         lb, copy.deepcopy(overrides)))
                _release_entry_copy(lb, k_ov)
    return out


def _unique_pairs(pairs):
    """(label, loss) pairs once per label, the largest loss kept, in label order."""
    best = {}
    for l, x in pairs:
        if l.id not in best or x > best[l.id][1]:
            best[l.id] = (l, x)
    return [best[k] for k in sorted(best)]


# what a server before the refusal of names that are not UTF-8 wrote for the bytes it could
# not carry: such a name may stand for another one in the index
_REPLACEMENT = '\ufffd'


def _unverifiable_names(g, labels):
    """Why the names of |labels| cannot be resubmitted as an explicit label list (None:
    they can). /traverse resolves a NAME to a label, so a request that names labels
    holds only when each name is exactly one label's: not shared by two labels of the
    retrieval, and not a name the library cannot know to be the index's own -- one with
    U+FFFD, which servers that replaced the bytes of names that are not UTF-8 wrote, so
    that the index may hold another label of exactly that name (a server that refuses such
    a seed writes none). The library has only the retrieval to go by."""
    if not labels:
        return None
    counts = g.cache.get('name_counts')
    if counts is None:
        counts = g.cache['name_counts'] = collections.Counter(l.name for l in g.labels)
    why = []
    for l in labels:
        if counts[l.name] > 1:
            why.append('%s: %d labels of this retrieval are named %s (%s)' % (
                l.ref, counts[l.name], _shown(l.name),
                ', '.join(x.ref for x in g.labels if x.name == l.name)))
        elif _REPLACEMENT in l.name:
            why.append('%s: its name %s holds U+FFFD, which servers wrote in place of '
                       'bytes that are not UTF-8: the index may hold another label of '
                       'exactly this name' % (l.ref, _shown(l.name)))
    if not why:
        return None
    return ('the request would name labels by names that may resolve to other labels '
            '(/traverse resolves names): %s; not built -- name the labels explicitly in a '
            'request of your own once verified (e.g. with /resolve)' % '; '.join(why[:3])
            + (' and %d more' % (len(why) - 3) if len(why) > 3 else ''))


def _shown(name, limit=32):
    return repr(name if len(name) <= limit else name[:limit] + '...')


def _shown_bytes(name, limit=32):
    """A bound of the bytes _shown(name) writes, without writing them: its quotes and each
    of at most limit + 3 characters -- an ASCII one kept or escaped by a backslash, any
    other escaped as \\U........ at worst."""
    n = min(len(name), limit + 3)
    return 2 + (2 * n if name.isascii() else 10 * n)


def resubmittable_names(g, labels):
    """The names of |labels| (selectors) as a request would carry them; raises
    UnverifiableLabelName when one may resolve to another label (see continuation)."""
    labs = labels_matching(g, labels) or []
    why = _unverifiable_names(g, labs)
    if why:
        raise UnverifiableLabelName(why)
    return [l.name for l in labs]


def _merge_overrides(strategy, overrides, path='strategy'):
    """Deep-merge the keyword overrides of next_request() into |strategy| as the server will
    read the result:
      * labels.change_cost is one unit: an override that names another model replaces the
        retrieval's whole (a field-by-field merge would keep the old model's fields -- a
        constant's value under a table, a table's default and entries under a constant --
        which the server's Strict parse refuses, so no continuation could change its cost
        model);
        one that names the same model, or none, updates its fields;
      * None for a section the strategy holds as an object (labels, labels.change_cost,
        branching, bounds, ...) is a ValueError naming it: it would be read as {} or as
        forbid by the rebuild and then sent as JSON null, which the server refuses
        ("expected an object"), dropping the rebuilt labels.extra, loss_budget or
        max_label_branches."""
    for k, v in overrides.items():
        here = '%s.%s' % (path, k)
        cur = strategy.get(k)
        if v is None and isinstance(cur, dict):
            raise ValueError('%s is an object, not None: the server refuses null there (omit '
                             'the override to keep the retrieval\'s)' % here)
        if isinstance(v, dict) and isinstance(cur, dict):
            if here == 'strategy.labels.change_cost' and 'model' in v \
                    and v['model'] != cur.get('model', 'forbid'):
                strategy[k] = copy.deepcopy(v)
            else:
                _merge_overrides(cur, v, here)
        else:
            strategy[k] = copy.deepcopy(v)


# ------------------------------------------------------------------ views

class GraphletView:
    """A view over a backing graphlet (§5): every id is the backing graphlet's
    original id (no remapping, no new record type). save() writes the UNCHANGED
    backing body plus, in J, {view: {selectors, mode, arm, of}}; load() restores the
    view. Its completeness is the backing complete_to_bp, qualified 'for the selected
    labels'."""

    __slots__ = ('backing', 'selectors', 'mode', 'arm_side', 'labels', 'segments',
                 'path_ids', 'of')

    def __init__(self, backing, selectors, arm=None, mode='any', of=None, *, budget=None):
        """budget= (local limits): building the view is charged (the selected segments'
        chains, and without |of| the dump and digest of the whole body); a stop raises, with
        no partial."""
        b = _B.resolve(budget)
        if b is None:
            self._build(backing, selectors, arm, mode, of, None)
            return
        with b.scope('subgraph', ('narrow_labels', 'narrow_arm')):
            self._build(backing, selectors, arm, mode, of, b)

    def _build(self, backing, selectors, arm, mode, of, b):
        if mode not in ('any', 'all'):
            raise ValueError("mode is 'any' or 'all'")
        self.backing = backing
        if b is not None:
            # charged first: resolving names builds the name index
            _label_ids_b(backing, selectors, b)
        labs = labels_matching(backing, selectors)
        if not labs:
            # a view of no labels: None would crash with a TypeError, and [] in mode 'all'
            # would select every segment (an empty set is a subset of every one)
            raise BadSelector('a view needs at least one label selector, not %r'
                              % (selectors,))
        self.selectors = [{'ref': l.ref} for l in labs]
        self.mode = mode
        self.arm_side = None if arm is None else backing.arm(arm).side
        self.labels = labs
        if of:
            self.of = of
        else:
            from . import parser
            text = parser.dump(backing, envelope=False, budget=b)
            if b is not None:
                b.phase = 'digest'
                b.charge(len(text) >> 8, len(text) + STR)
            self.of = 'sha256:' + hashlib.sha256(text.encode('utf-8')).hexdigest()
            del text
        ids = {l.id for l in self.labels}
        self.segments = {}
        self.path_ids = {}
        for side, a in backing.arms.items():
            if self.arm_side is not None and side != self.arm_side:
                continue
            hits = []
            rbs = derive.runs_by_segment(a)
            if b is not None:
                derive.uses(b, backing, a, 'leaves', 'paths', 'end_labels')
                b.charge(W_ELEM * len(a.segments), set_bytes(len(a.segments)))
                b.phase = 'segments'
            for s in a.segments:
                if b is not None:
                    n = len(s.entry) + len(s.end) + len(rbs[s.id]) \
                        + sum(len(p.labels) for p in s.presence)
                    b.charge(W_ELEM * (2 + len(s.presence)) + (n >> _B.ID_SHIFT_PY),
                             set_bytes(n))
                present = set(s.entry) | set(s.end) | {a.runs[r].label for r in rbs[s.id]}
                for p in s.presence:
                    present.update(p.labels)
                hit = ids & present
                if b is not None:
                    b.release(set_bytes(n))
                if (mode == 'any' and hit) or (mode == 'all' and hit == ids):
                    if b is not None:
                        # (priced as a walk of the hit's chain: the pass below is cheaper,
                        # the work model charges the walk)
                        b.charge(2 * W_STEP * (s.depth + 1), list_bytes(s.depth + 1))
                        b.release(list_bytes(s.depth + 1))
                    hits.append(s.id)
            # the union of the hits' first-parent chains in one pass from the last segment
            # (a parent's id is below its child's): walking each hit's chain to the root
            # would be quadratic on a comb
            keep = set(hits)
            for s in reversed(a.segments):
                if s.id in keep and s.parents:
                    keep.add(s.parents[0])
            if b is not None:
                ps = derive.paths(a)
                b.charge(sort_work(len(keep)) + W_ELEM * len(ps),
                         list_bytes(len(keep)) + list_bytes(len(ps)))
            self.segments[side] = sorted(keep)
            self.path_ids[side] = [p.id for p in derive.paths(a) if p.leaf in keep
                                   and self._leaf_carries(a, p, ids)]

    def _leaf_carries(self, a, p, ids):
        g = self.backing
        if g.mode == 'constrain':
            have = {e.label for e in derive.end_labels(a, p.leaf)}
        else:
            have = set(a.segments[p.leaf].end)
        hit = ids & have
        return bool(hit) if self.mode == 'any' else hit == ids

    def arms(self):
        return list(self.segments)

    def walks(self, arm, *, top=None, by='support', labels=None, route_consistent=True,
              min_bp=0, budget=None, resume=None):
        """walks() of the backing graphlet restricted to the view's walks, in the same
        order: the view's walks are ranked, then cut to |top| -- cutting the backing arm's
        ranking to |top| first would return the view's walks among the arm's top ones, fewer
        than the view has, or none. Only the walks returned are built. budget=:
        as walks() -- by='id' stops with the whole walks made so far and a resume token,
        a ranked list with no partial."""
        a = self.backing.arm(arm)
        if resume is not None and by != 'id':
            raise ValueError("resume= continues walks(by='id') only: a ranked list is never "
                             "partial")
        # the view's walks, ascending path ids (built in path order), tested by bisection:
        # a set of them, built before the call's first charge, would be an allocation no
        # budget admits -- the membership test needs none
        keep = self.path_ids.get(a.side, ())
        b = _B.resolve(budget)
        if b is None:
            return walks_at(self.backing, a, self._ranked(a, keep, top, by, labels,
                                                          route_consistent, min_bp, None),
                            with_claims=True, resume=resume)
        with b.scope('walks', ('rank_id', 'select_walks', 'narrow_labels', 'narrow_arm')):
            ids = self._ranked(a, keep, top, by, labels, route_consistent, min_bp, b)
            try:
                return walks_at(self.backing, a, ids, with_claims=True, budget=b,
                                resume=resume)
            except LocalBudgetExceeded as e:
                if by != 'id':
                    e.partial = None     # a ranked list never returns a partial (walks())
                raise

    def _ranked(self, a, keep, top, by, labels, route_consistent, min_bp, b):
        n = len(keep)
        ids = [i for i in rank_walks(self.backing, a, by=by, labels=labels,
                                     route_consistent=route_consistent, min_bp=min_bp,
                                     budget=b)
               if n and keep[min(bisect.bisect_left(keep, i), n - 1)] == i]
        return ids if top is None else ids[:top]

    def claims(self, arm=None, at_most_bp=None, strict=True, *, budget=None, resume=None):
        if arm is None and self.arm_side is not None:
            arm = self.arm_side
        return self.backing.claims(arm, self.labels, at_most_bp, strict, budget=budget,
                                   resume=resume)

    def to_fasta(self, arm=None, with_seed=True, orientation='natural', width=None, *,
                 coordinates=None, budget=None):
        from . import export
        b = _B.resolve(budget)
        if b is not None:
            with b.scope('to_fasta', ('select_walks', 'narrow_arm', 'export_mgt')):
                return self._to_fasta(export, arm, with_seed, orientation, width, b,
                                      coordinates)
        return self._to_fasta(export, arm, with_seed, orientation, width, None, coordinates)

    def _to_fasta(self, export, arm, with_seed, orientation, width, b, coordinates=None):
        out = []
        for side in self.segments:
            if arm is not None and self.backing.arm(arm).side != side:
                continue
            out.append(export.to_fasta(self.backing, side, self.path_ids[side], with_seed,
                                       orientation, width, coordinates=coordinates,
                                       budget=b))
        if b is not None:
            n = sum(len(x) for x in out)
            b.charge(n >> 8, STR + n)
        return ''.join(out)

    def completeness(self):
        return {side: {'complete_to_bp': self.backing.arms[side].complete_to_bp,
                       'qualified': 'for the selected labels',
                       'scope': self.backing.arms[side].scope}
                for side in self.segments}

    def summary(self):
        return {'view': self.view_spec(), 'labels': [l.as_dict() for l in self.labels],
                'arms': {side: {'segments': len(self.segments[side]),
                                'walks': len(self.path_ids[side]),
                                **self.completeness()[side]}
                         for side in self.segments}}

    def view_spec(self):
        return {'selectors': self.selectors, 'mode': self.mode, 'arm': self.arm_side,
                'of': self.of}

    def _with_view(self):
        g = copy.copy(self.backing)
        g.view = self.view_spec()
        return g

    def dump(self, envelope=True, *, budget=None):
        """The backing body, unchanged; the view lives in J only (envelope=False: the
        body alone, which does not say it is a view)."""
        return self._with_view().dump(envelope=envelope, budget=budget)

    def save(self, path, *, budget=None):
        from . import parser
        return parser.save(self._with_view(), path, budget=budget)


def subgraph(g, selectors, arm=None, mode='any', *, budget=None, of=None):
    """GraphletView(g, selectors, arm, mode). |of|: the view's backing identity when the
    caller has one (a store handle): without it the backing body is dumped and hashed
    for its digest."""
    return GraphletView(g, selectors, arm, mode, of, budget=budget)


def view_from_spec(backing, spec, *, budget=None):
    """The GraphletView a saved J `view` spec (or a store entry's) describes over
    |backing|: the one restoration of a spec (parser.load() and GraphletStore.view() both
    call it). A spec that is no object, or names no label,
    is a GraphletFormatError -- never a view of every label or of none."""
    if not isinstance(spec, dict) or not spec.get('selectors'):
        raise GraphletFormatError(2, 'J: a saved view names no label selectors: %r'
                                  % (spec,))
    return GraphletView(backing, spec['selectors'], spec.get('arm'), spec.get('mode', 'any'),
                        spec.get('of'), budget=budget)


def view_from_saved(g, *, budget=None):
    """The view a loaded file saved (parser.load()): its graphlet becomes the backing,
    its own view cleared."""
    spec = g.view
    g.view = None
    return view_from_spec(g, spec, budget=budget)


# ------------------------------------------------------------------ comparison

def index_identity(g):
    """(§3.1) {ns, fp, meta_fp, release}: fp is the index bundle's manifest digest
    (None: no manifest, joins are unverifiable), meta_fp a NEGATIVE check only."""
    release = (g.envelope or {}).get('release') or None
    return {'ns': g.index_ns, 'fp': g.index_fp, 'meta_fp': g.index_meta_fp,
            'release': release}


def comparability(a, b):
    """Whether graphlets |a| and |b| can be compared at all -> (comparable, reason, notes):
    False when they come from different indexes (index_meta_fp, index_fp, namespace or
    release differ), have different oriented seed sequences, or name one label ref
    differently; 'unverifiable' when either side has no index manifest digest (index_fp);
    else True. |notes| is a new empty list for the caller's notes. compare() starts from
    this verdict and only ever weakens it."""
    ia, ib = index_identity(a), index_identity(b)
    notes = []
    if ia['meta_fp'] and ib['meta_fp'] and ia['meta_fp'] != ib['meta_fp']:
        return False, 'different indexes (index_meta_fp differs)', notes
    if ia['fp'] and ib['fp'] and ia['fp'] != ib['fp']:
        return False, 'different indexes (index_fp differs)', notes
    for k in ('ns', 'release'):
        if ia[k] and ib[k] and ia[k] != ib[k]:
            return False, 'different index %s (%r vs %r)' % (k, ia[k], ib[k]), notes
    if a.seed.sequence != b.seed.sequence:
        return False, 'different oriented seed sequences', notes
    # names must agree per ref (a cheap guard on top of the digest)
    names_a = {l.ref: l.name for l in a.labels}
    for l in b.labels:
        if l.ref in names_a and names_a[l.ref] != l.name:
            return False, 'label %s is %r in one and %r in the other' % (
                l.ref, names_a[l.ref], l.name), notes
    if not ia['fp'] or not ib['fp']:
        return 'unverifiable', ('no index manifest digest (index_fp) on %s: label joins '
                                'cannot be verified' % ('both' if not ia['fp'] and not ib['fp']
                                                        else 'one side')), notes
    return True, 'same index and seed', notes


# the private name, kept for callers that import it
_comparability = comparability


# how weak a comparability verdict is: a comparison only ever moves to a weaker one (an
# unverifiable index identity is not made 'qualified' by a scope note, and so on)
_COMPARABLE_RANK = {True: 0, 'qualified': 1, 'unverifiable': 2, 'unknown': 3, False: 4}

# the label-class limitations (spec §7.0, design §14): evidence may be missing or
# understated (lower bound) or overstated (qualified)
LOWER_BOUND_KINDS = ('label_lists', 'inexact_counts', 'seed_labels', 'switch_sources',
                     'greedy_losses')
# (and the derivation limitation of a walked result whose permitted set was derived from
# part of the seed, SPEC §7.0: derive.qualifies())
OVERSTATED_KINDS = ('trace_record_boundaries',)


# compare() of retrievals with record coordinates
COORDINATES_NOTE = 'coordinates are not compared'


def _displayed_parent_price(a, b, sides, depth, rules=None):
    """The lwu _displayed_parent_rules() charges to judge the merges, or None when it judges
    none (one rule on both sides, or no merge within [0, depth) of a compared arm): an
    element per merge parent and per run (constrain, the parse's run index) or presence
    run (annotate) of it, which the carried counts read. compare_cost() prices it too, so
    that its exact estimate stays exact."""
    ra, rb = rules or (derive.displayed_parent_rule(a), derive.displayed_parent_rule(b))
    if ra is not None and ra == rb:
        return None                     # one rule on both sides
    n = 0
    for g in (a, b):
        for s in sides:
            arm_ = g.arms[s]
            rbs = None
            for x in arm_.segments:
                if len(x.parents) < 2 or x.from_bp >= depth:
                    continue
                if g.mode == 'constrain' and rbs is None:
                    rbs = derive.runs_by_segment(arm_)
                for pid in x.parents:
                    n += 1 + (len(rbs[pid]) if g.mode == 'constrain'
                              else len(arm_.segments[pid].presence))
    return W_ELEM * n if n else None


def _displayed_parent_rules(a, b, sides, depth, cb=None):
    """The note of a comparison whose sides may display a merge through different parents,
    or None. From feature level 6 a merge displays the parent carried by the most labels,
    before it the first to arrive (the majority-parent rule, SPEC §7.1;
    derive.displayed_parent_rule(): a side's rule, or not known for a body without its
    envelope). The displays can differ only where a side that may follow arrival order
    shows, at a merge within [0, depth) of a compared arm, a parent the majority-parent
    rule would not (derive.majority_first() not True: a later parent carries more labels
    than the first, or the counts cannot be judged) while the other side may follow the
    majority-parent rule. A side whose merges all show their majority
    parent displays what either rule would, so two retrievals of one rule, or of any rules
    whose merges agree with the majority, are not qualified. |cb|: compare()'s budget,
    charged the carried counts where they are read (an element per merge parent and per
    run or presence run of it)."""
    ra, rb = derive.displayed_parent_rule(a), derive.displayed_parent_rule(b)
    work = _displayed_parent_price(a, b, sides, depth, (ra, rb))
    if work is None:
        return None
    if cb is not None:
        cb.charge(work)

    def shows_arrival_minority(g):
        # the side's display differs from the majority-parent rule's somewhere (or may:
        # None)
        return any(derive.majority_first(g, s, depth) is not True for s in sides)

    hit = None
    for old, rule_old, rule_new in ((a, ra, rb), (b, rb, ra)):
        if rule_old != 'majority' and rule_new != 'arrival' \
                and shows_arrival_minority(old):
            hit = 'a' if old is a else 'b'
            break
    if hit is None:
        return None
    level = derive.DISPLAYED_PARENT_LEVEL

    def stated(rule, g):
        if rule is None:
            return 'not known (no envelope)'
        lv = derive.feature_level(g)
        return 'feature level %d' % lv if lv else 'a level before %d' % level
    return ('the displayed parent at a merge may follow different rules on the two sides '
            '(%s on a, %s on b; from level %d the parent carried by the most labels, before '
            'it the first to arrive), and %s shows a merge within the depth through a '
            'parent carried by fewer labels than another (or one whose counts cannot be '
            'judged): a walk spelled through it, and the claims keyed by it, may differ '
            'where the retrieval does not; compare mode labels, or retrievals of one level'
            % (stated(ra, a), stated(rb, b), level, hit))


def _weaker(cur, new):
    return new if _COMPARABLE_RANK[new] > _COMPARABLE_RANK[cur] else cur


def _label_evidence_notes(g, sides):
    """Why |g|'s label evidence is not complete on |sides|: one note per stated
    label-class limitation, else the outcome dimension or the cut lists themselves (a
    graphlet whose K records were lost still says so through O and A)."""
    notes = []
    seen = set()
    for lim in g.limitations:
        if lim.arm is not None and lim.arm not in sides:
            continue
        if lim.kind in LOWER_BOUND_KINDS or derive.qualifies(g, lim):
            if lim.kind in seen:
                continue
            seen.add(lim.kind)
            notes.append('%s (%s): %s' % (
                lim.kind, lim.knob, 'the oracle is a lower bound' if lim.kind in
                LOWER_BOUND_KINDS else 'reported evidence may be overstated'))
    cut = sum(g.arms[s].labels_per_node.nodes_truncated for s in sides)
    if not notes and cut:
        notes.append('label_lists (labels.max_labels_per_node): %d recorded lists cut, the '
                     'oracle is a lower bound' % cut)
    if not notes and g.outcome is not None and g.outcome.label_evidence != 'complete':
        notes.append('label_evidence %s' % g.outcome.label_evidence)
    return notes


def compare(a, b, *, arm=None, labels=None, mode='claims', budget=None):
    """Comparison(comparable, reason, depth_used, equal, only_in_a, only_in_b, notes,
    differ). Keyed by LabelRef, restricted to min(complete_to_bp), claims cut there get
    end_class 'open'. Comparable only for the same index identity (§3.1) and the same
    oriented seed; a per_path vs united_history pair is 'qualified' (equal None), never
    equal, and so is a pair where either side's label evidence is a lower bound or
    qualified (cut lists cannot tell 'absent' from 'cut'); without index_fp on either side
    'unverifiable'; with no common certified depth (complete_to_bp 0) 'unknown' (equal
    None in every case but True).

    modes: claims (the §6.9 unit: run ends, keyed (arm, ref, evidence_from, displayed
    prefix)), walks (walks cut at the depth, with the labels supporting the whole cut
    prefix), labels (label_summary cut at the depth), prefix_subset (every walk of |a|
    under each label is a prefix of a walk of |b| under it; |b|'s walks |a| omits are
    listed with the reason |a| recorded for them).

    Both sides are compared over their DAGs RESTRICTED to [0, depth): the walks are the
    restricted DAG's leaves (restricted_leaves) and a claim reaching past the depth is
    anchored on its own route at the depth. Cutting the deeper DAG's displayed walks and
    claims instead would let merges beyond the depth decide what lies before it: a branch
    that enters such a merge through a non-first parent is shown by no walk of the deeper
    DAG, and its lineage would be route_only on the other branch's bases, so a retrieval and
    a deeper one of the same seed would differ (bubbles of radius 30 vs 100 at depth 30,
    merges at 36 and 62). A merge AT the depth lies outside the restricted DAG too: a
    retrieval walked to a merge position records it as a zero-length segment, and the runs
    it closes are open claims at the depth, not merged ones. Where a
    route cannot be reconstructed through a merge (its lineage in no partition) the
    comparison is 'qualified', never a difference.

    Label evidence and the per-node label limit: 'qualified' can come from
    labels.max_labels_per_node alone. A retrieval whose label lists were cut at a node (more
    labels than the limit, 64 by default: a seed carried by many labels, or a conserved
    region) states label_evidence lower_bound, and a comparison with it is qualified, never
    equal -- a cut list cannot tell 'absent' from 'cut'. compare() does not raise the limit
    for you: fetch both sides again with a larger labels.max_labels_per_node (or name fewer
    labels) when equality must be decided.

    Record coordinates are not compared: two retrievals with them compare
    exactly as without them, and the answer notes it ('coordinates are not compared').

    Displayed parents (the majority-parent rule, SPEC §7.1): from feature level 6 a merge
    displays the parent carried by the most labels, before it the first to arrive. Claims,
    walks and prefix_subset are keyed by the displayed bases, so a comparison is
    'qualified', and a note says why, when one side may follow the arrival-order rule (it
    states a level below 6, its capabilities state none, or it is a body without its
    envelope) and shows, at a merge before the depth, a parent carried by fewer labels than
    another, while the other side may follow the majority-parent rule (states 6 or more, or
    is a body without its envelope); mode labels is not affected. A side whose merges all
    show their majority parent displays what either rule would.

    The cost follows both DAGs up to the depth -- each segment read once per side (in modes
    walks and prefix_subset too: reading each walk's chain again per walk would be quadratic
    on a comb) -- and, in those modes, the walks' spellings to it, which on a deep comb grow
    as the square of the depth, as an export of its walks does. Without a budget (the
    default) no work or allocation budget applies; budget= (local limits) never raises for
    the budget: a stopped comparison RETURNS comparable 'unknown', equal None, no difference
    lists and local_stop (the stop: resource, phase, how far keying got) -- a one-sided list
    from a half-keyed side would be a false difference. compare_cost() estimates the
    charge."""
    cb = _B.resolve(budget)
    if cb is None:
        return _compare(a, b, arm, labels, mode, None, {})
    state = {}
    with cb.scope('compare', ('narrow_labels', 'narrow_arm')):
        try:
            return _compare(a, b, arm, labels, mode, cb, state)
        except LocalBudgetExceeded as e:
            return _stopped_comparison(a, b, mode, state, e)


def _stopped_comparison(a, b, mode, state, e):
    """The answer of a comparison its budget stopped: neither equality nor a difference."""
    st = e.stop
    done = {k: v for k, v in state.items() if k in ('keys_a', 'keys_b')}
    st.restate(**done)
    reason = ('the comparison stopped on the %s in phase %s (%s): neither equality nor a '
              'difference is asserted; %s'
              % ({'work': 'local work budget', 'memory': 'local memory budget',
                  'deadline': 'local deadline', 'cancelled': 'cancellation'}[st.resource],
                 st.phase, ', '.join('%s %s' % (k, v) for k, v in done.items())
                 or 'before keying', 'raise the budget, or narrow the labels or the arm'))
    return Comparison('unknown', reason, state.get('depth'), None, [], [],
                      list(state.get('notes') or []), [], mode=mode,
                      support=(a.support, b.support),
                      scopes=state.get('scopes') or ('per_path', 'per_path'),
                      strategies=((a.envelope or {}).get('strategy'),
                                  (b.envelope or {}).get('strategy')),
                      local_stop=st.as_dict())


def _compare_sides(a, b, arm):
    """The arms a comparison reads -> (sides, the reason there are none). An explicit arm
    is resolved on |a| and must be one |b| retrieved too (looked up on |b| it would raise a
    KeyError)."""
    if arm is not None:
        side = a.arm(arm).side
        if side not in b.arms:
            return [], 'the %s arm was not retrieved in b' % side
        return [side], None
    sides = [s for s in ARM_SIDES if s in a.arms and s in b.arms]
    return sides, None if sides else 'no arm in common'


_COMPARE_MODES = ('claims', 'walks', 'labels', 'prefix_subset')


def _compare(a, b, arm, labels, mode, cb, state):
    # the mode before anything is answered: otherwise an incomparable pair would echo a mode
    # that a comparable one refuses
    if mode not in _COMPARE_MODES:
        raise ValueError("mode is 'claims', 'walks', 'labels' or 'prefix_subset'")
    if cb is not None:
        cb.charge(W_ELEM * (2 * len(a.labels) + 2 * len(b.labels)),
                  dict_bytes(len(a.labels) + len(b.labels))
                  + (len(a.labels) + len(b.labels)) * (2 * STR + 24))
    comparable, reason, notes = comparability(a, b)
    strategies = ((a.envelope or {}).get('strategy'), (b.envelope or {}).get('strategy'))
    if comparable is False:
        return Comparison(False, reason, None, None, [], [], notes, mode=mode,
                          support=(a.support, b.support), strategies=strategies)
    sides, none = _compare_sides(a, b, arm)
    if not sides:
        # the support kinds and strategies of both sides, as every other answer states
        # them (not the defaults, ('kmer', 'kmer'), for any pair)
        return Comparison(False, none, None, None, [], [], notes, mode=mode,
                          support=(a.support, b.support), strategies=strategies)
    scopes = (a.arms[sides[0]].scope, b.arms[sides[0]].scope)
    depth = min(min(a.arms[s].complete_to_bp, b.arms[s].complete_to_bp) for s in sides)
    state.update(depth=depth, scopes=scopes, notes=notes)
    if cb is not None:
        cb.charge(W_ELEM * sum(len(g.arms[s].segments) for g in (a, b) for s in sides))
    if depth == 0:
        # equality over an empty window is not equality: nothing was certified on one side
        comparable = _weaker(comparable, 'unknown')
        notes.append('no common certified depth (complete_to_bp 0 on a or b)')
    for tag, g in (('a', a), ('b', b)):
        evidence = _label_evidence_notes(g, sides)
        if evidence:
            comparable = _weaker(comparable, 'qualified')
            notes.extend('%s: %s' % (tag, n) for n in evidence)
    for s in sides:
        if a.arms[s].scope != b.arms[s].scope:
            comparable = _weaker(comparable, 'qualified')
            notes.append('%s arm: completeness scopes differ (%s vs %s): a structural block '
                         'on a merged route need not be where the per-path trie ends a label'
                         % (s, a.arms[s].scope, b.arms[s].scope))
        if a.arms[s].scope == 'united_history' and any(
                len(x.parents) > 1 for x in a.arms[s].segments):
            notes.append('%s arm of a: histories were united at merges' % s)
    if a.support != b.support:
        notes.append('support kinds differ (%s vs %s)' % (a.support, b.support))
    if _C.present(a) or _C.present(b):
        # the record coordinates ride beside the claims, never in their keys, and are
        # neither clipped nor charged here
        notes.append(COORDINATES_NOTE)
    if mode in ('claims', 'walks', 'prefix_subset'):
        shown = _displayed_parent_rules(a, b, sides, depth, cb)
        if shown is not None:
            # which parent a merge displays is the server's choice, and it differs below and
            # from feature level 6 (the majority-parent rule); claims and walks are keyed by
            # the displayed bases, so the same trie can spell a walk through a merge
            # differently on the two sides. A difference there may be the display's, not the
            # retrieval's: qualified, as when the completeness scopes differ
            comparable = _weaker(comparable, 'qualified')
            notes.append(shown)
    if cb is not None and labels is not None:
        _uses_g(cb, a, 'label_index')
        _uses_g(cb, b, 'label_index')
    sel = _compare_selectors(a, b, labels)
    only_a, only_b, differ = [], [], []
    names = {l.ref: l.name for l in b.labels}
    names.update({l.ref: l.name for l in a.labels})

    def lab(ref):
        return {'name': names.get(ref), 'ref': ref}

    # claims, walks and prefix_subset are keyed by displayed bases (a claim's prefix, a
    # walk's spelling): a side retrieved with output.sequences false has none. Keying it
    # by anything else would report every claim of the identical trie as one-sided, or
    # match claims of different walks, so the comparison is not made: 'unknown'
    no_bases = [tag for tag, g in (('a', a), ('b', b))
                if any(not _has_bases(g, s) for s in sides)]
    memo = {}
    evaluated = True
    inexact = []
    if mode in ('claims', 'walks', 'prefix_subset') and not no_bases:
        beyond = {tag: sum(1 for s in sides for x in g.arms[s].segments
                           if len(x.parents) > 1 and x.from_bp >= depth)
                  for tag, g in (('a', a), ('b', b))}
        if any(beyond.values()):
            notes.append('compared over each DAG restricted to [0, %d): %s merge%s at or '
                         'beyond the depth (%s) %s not part of it, so a walk entering one '
                         'through a non-first parent is a walk of its own at the depth, as '
                         'in a retrieval walked to %d'
                         % (depth, sum(beyond.values()),
                            '' if sum(beyond.values()) == 1 else 's',
                            ', '.join('%d on %s' % (n, t) for t, n in beyond.items() if n),
                            'is' if sum(beyond.values()) == 1 else 'are', depth))
    if mode in ('claims', 'walks', 'prefix_subset') and no_bases:
        comparable = _weaker(comparable, 'unknown')
        evaluated = False
        reason = ('%s carries no bases (output.sequences false): mode %s matches %s by '
                  'their displayed bases, so it cannot compare them; retrieve both with '
                  'output.sequences true, or compare mode labels'
                  % (' and '.join(no_bases) + (' each' if len(no_bases) > 1 else ''), mode,
                     'claims' if mode == 'claims' else 'walks'))
    elif mode == 'claims':
        ka = _keys_of(cb, state, 'a', lambda: _claim_keys(a, sides, depth, sel, inexact, cb))
        kb = _keys_of(cb, state, 'b', lambda: _claim_keys(b, sides, depth, sel, inexact, cb))
        _charge_differences(cb, ka, kb, dict_bytes(6) + dict_bytes(2))
        only_a = [_key_json(k, v, lab) for k, v in sorted(ka.items()) if k not in kb]
        only_b = [_key_json(k, v, lab) for k, v in sorted(kb.items()) if k not in ka]
        differ = [{'key': _key_json(k, None, lab), 'a': ka[k], 'b': kb[k]}
                  for k in sorted(ka) if k in kb and ka[k] != kb[k]]
    elif mode == 'walks':
        wa = _keys_of(cb, state, 'a', lambda: _walk_keys(a, sides, depth, sel, memo, cb))
        wb = _keys_of(cb, state, 'b', lambda: _walk_keys(b, sides, depth, sel, memo, cb))
        _charge_differences(cb, wa, wb, dict_bytes(3) + LIST, per_item=dict_bytes(2))
        only_a = [{'arm': k[0], 'walk': k[1], 'labels': [lab(r) for r in sorted(v)]}
                  for k, v in sorted(wa.items()) if k not in wb]
        only_b = [{'arm': k[0], 'walk': k[1], 'labels': [lab(r) for r in sorted(v)]}
                  for k, v in sorted(wb.items()) if k not in wa]
        differ = [{'arm': k[0], 'walk': k[1], 'a': [lab(r) for r in sorted(wa[k])],
                   'b': [lab(r) for r in sorted(wb[k])]}
                  for k in sorted(wa) if k in wb and wa[k] != wb[k]]
    elif mode == 'labels':
        la = _keys_of(cb, state, 'a', lambda: _label_keys(a, sides, depth, sel, cb))
        lb = _keys_of(cb, state, 'b', lambda: _label_keys(b, sides, depth, sel, cb))
        _charge_differences(cb, la, lb, dict_bytes(5))
        only_a = [dict(lab(k[1]), arm=k[0]) for k in sorted(la) if k not in lb]
        only_b = [dict(lab(k[1]), arm=k[0]) for k in sorted(lb) if k not in la]
        differ = [dict(lab(k[1]), arm=k[0], a=la[k], b=lb[k])
                  for k in sorted(la) if k in lb and la[k] != lb[k]]
    else:
        only_a, only_b = _prefix_subset(a, b, sides, depth, sel, memo, inexact, cb, state)
        for row in only_a + only_b:
            row['name'] = names.get(row['ref'])
    if inexact:
        # a route that cannot be reconstructed keys a claim on the wrong branch: what
        # differs may be that, so no difference (and no equality) is asserted
        comparable = _weaker(comparable, 'qualified')
        notes.append('the DAG restricted to the depth could not be reconstructed exactly: '
                     '%s%s' % ('; '.join(inexact[:3]),
                               ' and %d more' % (len(inexact) - 3) if len(inexact) > 3
                               else ''))
    if mode == 'prefix_subset':
        holds = not only_a
    else:
        holds = not (only_a or only_b or differ)
    equal = holds if comparable is True else None
    if comparable is not True and holds and evaluated:
        notes.append('no difference found, but equality is not claimed (%s)' % comparable)
    return Comparison(comparable, reason, depth, equal, only_a, only_b, notes, differ,
                      mode=mode, support=(a.support, b.support), scopes=scopes,
                      strategies=strategies)


def _keys_of(cb, state, tag, fn):
    """One side's comparison keys, the phase named (keys:a, keys:b) and their count kept
    for the statement of a stop."""
    if cb is not None:
        cb.phase = 'keys:' + tag
    got = fn()
    state['keys_' + tag] = len(got)
    return got


def _charge_differences(cb, ka, kb, row_bytes, per_item=0):
    """The difference lists: both sides' keys sorted and looked up in the other's, and as
    many rows as keys at most (each with its labels, |per_item| apiece, where it has
    some) -- charged as that bound, before they are built."""
    if cb is None:
        return
    cb.phase = 'differences'
    n = len(ka) + len(kb)
    items = 0
    if per_item:
        items = sum(len(v) for v in ka.values()) + sum(len(v) for v in kb.values())
    cb.charge(sort_work(len(ka)) + sort_work(len(kb)) + 2 * n + W_ROW * n + 2 * items,
              2 * list_bytes(n) + n * row_bytes + items * per_item)


class _Tally:
    """The charging interface of a LocalBudget with no limits, which only adds up: what
    compare_cost() prices compare()'s structural phases with (the same expressions)."""

    def __init__(self):
        self.work = self.mem = self.peak = 0
        self._used = set()
        self.phase = None
        self.done = {}

    def charge(self, work, nbytes=0):
        self.work += work
        self.mem += nbytes
        self.peak = max(self.peak, self.mem)

    def uses(self, key, work, nbytes=0):
        if key not in self._used:
            self._used.add(key)
            self.charge(work, nbytes)

    def release(self, nbytes):
        self.mem = max(0, self.mem - nbytes)

    def first(self, key):
        if key in self._used:
            return False
        self._used.add(key)
        return True


def compare_cost(a, b, arm=None, mode='claims', *, labels=None):
    """An estimate, before running it, of what compare(a, b, arm=, labels=, mode=) charges
    a budget (local limits, the current WORK_MODEL): {mode, work_model, work_units: {at_least,
    estimate}, memory_bytes: {at_least, estimate}, exact, phases, unpriced}.

    at_least is what compare() certainly charges once it keys the two sides: the checks,
    the derivations its keyers use (cold price) and the per-run and per-cut charges, all
    structural. estimate adds the parts whose size depends on the answer -- the spellings
    of the claims' prefixes, the keys, and as many difference rows as keys (an upper
    bound for those) -- and, for prefix_subset, its scans priced on the supported labels
    of each cut, which are only known while comparing: those phases are named in
    |unpriced|. exact: compare() stops before keying (incomparable indexes, no arm in
    common, a side without bases in a mode keyed by bases), so at_least is all it
    charges. Cost: one pass over each arm's segments and runs (its structural sizes)."""
    if mode not in _COMPARE_MODES:
        raise ValueError("mode is 'claims', 'walks', 'labels' or 'prefix_subset'")
    t = _Tally()
    t.charge(_B.W_CALL + W_ELEM * (2 * len(a.labels) + 2 * len(b.labels)),
             _B.CALL_BASE + dict_bytes(len(a.labels) + len(b.labels))
             + (len(a.labels) + len(b.labels)) * (2 * STR + 24))
    phases = {}

    def result(exact, unpriced=(), extra_w=0, extra_m=0):
        phases.setdefault('setup', t.work)
        return {'mode': mode, 'work_model': _B.WORK_MODEL,
                'work_units': {'at_least': t.work, 'estimate': t.work + extra_w},
                'memory_bytes': {'at_least': t.peak, 'estimate': t.peak + extra_m},
                'exact': exact, 'phases': dict(phases), 'unpriced': list(unpriced)}

    comparable, _, _ = comparability(a, b)
    if comparable is False:
        return result(True)
    sides, _ = _compare_sides(a, b, arm)
    if not sides:
        return result(True)
    depth = min(min(a.arms[s].complete_to_bp, b.arms[s].complete_to_bp) for s in sides)
    t.charge(W_ELEM * sum(len(g.arms[s].segments) for g in (a, b) for s in sides))
    if mode in ('claims', 'walks', 'prefix_subset'):
        # the carried counts of the merges, where the two sides' display rules may differ
        # (the majority-parent rule)
        w = _displayed_parent_price(a, b, sides, depth)
        if w:
            t.charge(w)
    if labels is not None:
        _uses_g(t, a, 'label_index')
        _uses_g(t, b, 'label_index')
    phases['setup'] = t.work
    no_bases = any(not _has_bases(g, s) for g in (a, b) for s in sides)
    if mode in ('claims', 'walks', 'prefix_subset') and no_bases:
        return result(True)
    est_w = est_m = 0
    unpriced = []
    for tag, g in (('a', a), ('b', b)):
        w0 = t.work
        for side in sides:
            arm_ = g.arms[side]
            segs = arm_.segments
            if mode in ('claims',) or (mode == 'prefix_subset' and tag == 'a'):
                n_rows = _price_restricted_claims(t, g, side, depth)
                bp = sum(min(segs[r.segment].end_bp, depth) for r in arm_.runs) \
                    if g.mode == 'constrain' else n_rows * depth
                steps = sum(segs[r.segment].depth + 1 for r in arm_.runs) \
                    if g.mode == 'constrain' else n_rows * 4
                if mode == 'prefix_subset':
                    # the tree of the claims' chains and the claims filed under their
                    # anchors (_claim_index(): no prefix is spelled)
                    est_w += 2 * W_STEP * steps + (bp >> 8) + 3 * W_ELEM * n_rows \
                        + sort_work(n_rows)
                    est_m += min(steps, len(segs)) * (SET_ITEM + 2 * DICT_KEY + INT
                                                     + LIST_ITEM) \
                        + n_rows * (2 * TUPLE + 5 * 8 + 2 * LIST_ITEM + 2 * DICT_KEY
                                    + 2 * LIST + INT)
                else:
                    est_w += 2 * W_STEP * steps + (bp >> 8) + 3 * W_ELEM * n_rows \
                        + bp // 128
                    # the claims' prefixes (2x their bases, as charged) and their spelling
                    # streamed (derive.walk_iter_bytes(): the kept prefixes, at most the
                    # bases)
                    est_m += 2 * n_rows * (STR + TUPLE + LIST + 48 + 2 * STR) + 3 * bp
                if mode == 'claims':
                    est_w += sort_work(n_rows) + 2 * n_rows + W_ROW * n_rows
                    est_m += 2 * list_bytes(n_rows) + n_rows * (dict_bytes(6) + dict_bytes(2))
            if mode in ('walks', 'prefix_subset'):
                cuts = _price_cuts(t, g, side, depth)
                n0 = len(segs[0].entry) if segs else 0
                est_w += cuts * (2 * W_ELEM + W_ELEM * n0 + (depth >> 8))
                est_m += cuts * (STR + depth + TUPLE + set_bytes(n0) + n0 * (STR + 24))
                if mode == 'prefix_subset':
                    unpriced.append('prefix_scan:%s:%s' % (tag, side))
                else:
                    # as many difference rows as keys, each with its labels
                    est_w += sort_work(cuts) + (2 + W_ROW) * cuts + 2 * cuts * n0
                    est_m += list_bytes(cuts) + cuts * (dict_bytes(3) + LIST) \
                        + cuts * n0 * dict_bytes(2)
            if mode == 'labels':
                n = _clipped_summary_price(t, g, arm_)
                # the keys: compare() charges the labels the clipped summary holds, which
                # n (the most it can hold) bounds from above -- an estimate, not at_least
                # (charged into at_least, it exceeded a comparison with labels=)
                est_w += W_ELEM * 2 * n + sort_work(n) + (2 + W_ROW) * n
                est_m += n * (TUPLE + 2 * STR + 24 + DICT_KEY) + list_bytes(n) \
                    + n * dict_bytes(5)
        phases['keys:' + tag] = t.work - w0
    if mode == 'prefix_subset':
        unpriced.append('omissions')
    else:
        unpriced.append('differences')
    return result(False, unpriced, est_w, est_m)


def _price_restricted_claims(t, g, side, depth):
    """_restricted_claims()'s charges into |t|; -> the claims it prices (rows)."""
    a = g.arms[side]
    if g.mode == 'constrain':
        derive.uses(t, g, a, 'leaves', 'paths', 'path_of_leaf', 'merge_above', 'merge_depth',
                    'first_paths', 'partition_sets')
        t.uses((id(a), 'first_parent_at', depth - 1),
               *derive.price(a, 'first_parent_at', len(g.labels), g.mode))
        _charge_all(t, _restricted_claim_prices(t, g, a), _CLAIM_BYTES)
        return len(a.runs)
    derive.uses(t, g, a, 'leaves', 'paths', 'path_of_leaf', 'merge_above', 'merge_depth',
                'first_paths', 'annotate_route_ends')
    w, m = derive.price(a, 'annotate_route_ends', len(g.labels), g.mode)
    t.uses((id(a), 'route_ends_at', depth), w, m)
    # the ends of the DAG restricted to the depth, as compare() keys them: the whole DAG's
    # ends would price every end beyond the depth too, and at_least would exceed what the
    # comparison charges
    ends = _annotate_route_ends_at(a, depth)
    md = derive.merge_depth(a)
    n = sum(len(v) for v in ends.values())
    t.charge(sort_work(len(ends)) + n * (W_ROW + 2)
             + sum(3 * md[x] for v in ends.values() for x, _ in v),
             n * (_CLAIM_BYTES + TUPLE) + list_bytes(len(ends)))
    return n


def _price_cuts(t, g, side, depth):
    """_cuts()'s charges into |t| (the cuts found as _cuts() finds them, and the one pass
    that makes them priced by _cut_prices(), as _cuts() charges them); -> the distinct
    cuts."""
    a = g.arms[side]
    segs = a.segments
    derive.uses(t, g, a, 'leaves', 'paths', 'cut_cost',
                *(('segment_ops',) if g.mode == 'constrain' else ()))
    cost = derive.cut_cost(a, g.mode)
    held = _cut_transient(t, g, a)
    t.charge((W_ELEM + W_STEP) * len(segs))
    if depth <= 0:
        return 1
    plan = _cut_plan(a, depth)
    pr = _cut_prices(g, a, plan, cost, held)
    # the pass's charges in the pass's own order, so that the estimate's account rises and
    # falls as the comparison's does
    if plan['roots']:
        t.charge(0, pr['pass'])
    for step, x, m, _ in _cut_steps(plan, segs):
        if step == 'cut':
            w, keep, tr = pr['cut'](x, m)
            t.charge(w, keep + tr)
            t.release(tr)
        elif step == 'enter':
            t.charge(pr['seg'](x), held[x])
            t.release(held[x])
    if plan['roots']:
        t.release(pr['pass'])
    for _ in plan['rows']:
        t.charge(W_ELEM, TUPLE + 3 * LIST_ITEM)
    return len(plan['keys'])


def _clipped_summary_price(t, g, arm):
    derive.uses(t, g, arm)
    z = derive.arm_sizes(arm)
    n = min(len(g.labels), z['entry_ids'] + z['runs'] + z['presence_ids'])
    if g.mode == 'constrain':
        t.charge(6 * z['runs'], list_bytes(z['runs']) + n * (dict_bytes(2) + 2 * INT))
    else:
        t.charge(W_ELEM * (3 * z['segments'] + z['presence'])
                 + ((2 * z['presence_ids'] + 2 * z['entry_ids']) >> 2),
                 list_bytes(z['segments']) + z['segments'] * SET
                 + z['entry_ids'] * SET_ITEM + n * (dict_bytes(2) + 2 * INT))
    return n


def _compare_selectors(a, b, labels):
    """The LabelRefs a comparison is restricted to (None: every label). Comparisons are
    keyed by LabelRef, and a label only one side recorded (a cut list, a lower-bound
    retrieval) is exactly what they must show: a selector resolves against |a|, then
    |b|; {'ref': ...} is the key itself. {'id': ...} (and a bare int) are ids of |a|."""
    if labels is None:
        return None
    if isinstance(labels, (str, dict, int, Label)):
        labels = [labels]
    out = set()
    known = None
    for x in labels:
        if isinstance(x, dict) and len(x) == 1 and 'ref' in x and isinstance(x['ref'], str):
            if known is None:
                known = {l.ref for l in a.labels} | {l.ref for l in b.labels}
            if x['ref'] not in known:
                raise UnknownLabel('no label in either graphlet has ref %r' % x['ref'])
            out.add(x['ref'])
            continue
        if isinstance(x, bool) or isinstance(x, int) or (
                isinstance(x, dict) and len(x) == 1 and 'id' in x):
            out.add(label(a, x).ref)
            continue
        try:
            out.add(label(a, x).ref)
        except UnknownLabel:
            try:
                out.add(label(b, x).ref)
            except UnknownLabel:
                raise UnknownLabel('no label in either graphlet matches %r' % (x,)) from None
    return out


def _has_bases(g, side):
    segs = g.arms[side].segments
    return not segs or segs[0].walk is not None


def _key_json(k, v, lab):
    j = {'arm': k[0], 'label': lab(k[1]), 'evidence_from': None if k[2] < 0 else k[2],
         'prefix': k[3], 'to_bp': k[4]}
    if v is not None:
        j['end_class'] = v
    return j


def _open_at(c, depth):
    """A claim reaching the comparison depth is censored there on both sides: what ends
    AT the depth is decided by base |depth|, which the shallower retrieval never looked
    at (§5: claims cut at D get end_class open)."""
    return c.end_class == 'open' or (depth > 0 and c.to_bp >= depth)


def _claim_keys(g, sides, depth, sel, inexact=None, cb=None):
    """The claims of each arm's DAG RESTRICTED to [0, depth) (_restricted_claims) as
    comparison keys -> {(side, ref, evidence_from or -1, displayed prefix, to_bp):
    end class}."""
    out = {}
    for side in sides:
        a = g.arms[side]
        rows = [c for c in _restricted_claims(g, side, depth, inexact, cb)
                # a merged claim is not an end: the lineage goes on in the kept entry
                if c.kind != 'merged' and (sel is None or c.label.ref in sel)]
        if cb is not None:
            _charge_spellings(cb, a, {c.segment for c in rows})
            cb.charge(W_ELEM * 3 * len(rows) + sum(c.to_bp for c in rows) // 128,
                      len(rows) * (TUPLE + 48 + 2 * STR + 2 * DICT_KEY + 16)
                      + 2 * sum(c.to_bp for c in rows))
        prefixes = _claim_prefixes(a, rows, natural=True)
        for c, prefix in zip(rows, prefixes):
            key = (side, c.label.ref, c.evidence_from if c.evidence_from is not None else -1,
                   prefix, c.to_bp)
            out[key] = 'open' if _open_at(c, depth) else c.end_class
    return out


def _claim_prefixes(arm, rows, strict=False, natural=False):
    """[the walking-order bases [0, to_bp) of each row's segment chain] for the claims
    |rows|, in their order (|natural|: in natural orientation, _natural()): each segment
    spelled once as derive.walk_iter() reaches it and dropped once its rows are cut from
    it, so the prefixes are held and not every spelling with them (holding them all, and
    the batch's prefixes, made compare() peak 4x higher on a comb than spelling claim by
    claim had). Without bases (or with bases on some segments only): all None, a key
    without a prefix, or with |strict| (a side that must carry them) the ValueError
    walk_bases() raises."""
    flip = natural and arm.side == 'left'
    out = [None] * len(rows)
    at = {}
    for i, c in enumerate(rows):
        at.setdefault(c.segment, []).append(i)
    try:
        for s, _, w in derive.walk_iter(arm, set(at), chains=False, check_first=False):
            if w is None:
                if strict:
                    raise ValueError(derive.NO_BASES)
                continue
            for i in at[s]:
                out[i] = w[:rows[i].to_bp][::-1] if flip else w[:rows[i].to_bp]
    except ValueError:
        # bases on some segments only, raised at the first chain without them: the
        # prefixes cut before it are dropped with |out|, so no key carries part of a side
        if strict:
            raise
        return [None] * len(rows)
    return out


def _charge_spellings(cb, arm, segs):
    """Charge derive.walk_iter() of |segs|, its bases streamed to a consumer that keeps
    none of them: the chains' steps (work) and derive.walk_iter_bytes()."""
    sg = arm.segments
    cb.charge(len(segs))
    steps = sum(sg[x].depth + 1 for x in segs)
    bp = sum(sg[x].end_bp for x in segs)
    big = max([sg[x].end_bp for x in segs], default=0)
    cb.charge(2 * W_STEP * steps + (bp >> 8),
              derive.walk_iter_bytes(len(segs), bp, min(steps, len(sg)), big))


def _charge_chain_tree(cb, arm, segs):
    """Charge a tree of the chains of |segs| (_refusal_tree(), _restricted_tree()): the
    chains' steps, priced as their spelling is (_charge_spellings()), and per segment of
    their union its children slot, its offset and the set that collects it -- no bases."""
    sg = arm.segments
    cb.charge(len(segs))
    steps = sum(sg[x].depth + 1 for x in segs)
    bp = sum(sg[x].end_bp for x in segs)
    cb.charge(2 * W_STEP * steps + (bp >> 8),
              min(steps, len(sg)) * (SET_ITEM + 2 * DICT_KEY + INT + LIST_ITEM)
              + len(segs) * LIST_ITEM)


def segment_support(g, arm, s):
    """The displayed support inside segment |s| of arm |arm| (an Arm and a Segment of
    graphlet |g|) as (from, to, frozenset of label ids), in order and one at a time: the
    alive sets (constrain) or the recorded P sets (annotate)."""
    if g.mode == 'constrain':
        return _iter_alive_pieces(arm, s)
    return ((p.from_bp, p.to_bp, frozenset(p.labels)) for p in s.presence)


# the private name, kept for callers that import it
_segment_support = segment_support


def _cut_transient(cb, g, a):
    """seg -> the bytes a cut ending in the first-parent subtree's chain through |seg|
    holds for one segment of that chain while _cut_info() reads it (the most over the
    chain root -> seg): constrain, the segment's label changes as tuples, their sort keys
    and list, the running set and one piece's set; annotate, one P run's set. Released
    once the cut is made: a cut read a whole segment's pieces at once uncharged (on the
    hairpin pair with 977 runs ending in one segment, 2.7x the cut's charge)."""
    def prices():
        segs = a.segments
        out = []
        if g.mode == 'constrain':
            ops_ = derive.segment_ops(a)
            for s in segs:
                k = ops_[s.id]
                own = k * (2 * TUPLE + 40 + LIST_ITEM) + list_bytes(k) \
                    + 2 * set_bytes(len(s.entry) + k)
                out.append(max(out[s.parents[0]], own) if s.parents else own)
        else:
            for s in segs:
                own = TUPLE + 24 + set_bytes(max((len(p.labels) for p in s.presence),
                                                 default=0))
                out.append(max(out[s.parents[0]], own) if s.parents else own)
        return out
    return derive.price_list(cb, g, a, 'cut_transient', prices)


def _cut_info(g, arm, anchor, m):
    """One cut [0, m) of the first-parent chain root -> |anchor| (the segment holding
    base m - 1; None for m = 0) -> (the walking-order prefix, None without bases; the
    labels with displayed support on all of [0, m); {label: j} the longest [0, j),
    j <= m, each label supports from the boundary on). Only the segments up to the
    anchor are read, each only up to m."""
    segs = arm.segments
    entry = segs[0].entry if segs else ()
    have = {l: 0 for l in entry}
    if anchor is None:
        return ('' if not segs or segs[0].walk is not None else None, frozenset(entry), have)
    parts = []
    # annotate: the P runs cover the nodes entered by steps from_bp + 1.. only, so the
    # seed boundary's set (the root's entry) is intersected in first, as every other
    # annotate reading does (_displayed_presence, the route ends, label_summary): from
    # the first P run on, a label absent at the boundary would count as supporting [0, m).
    # Constrain: the first piece is the root's entry with the switch-ins at
    # its first base applied already -- seeding it would drop those.
    inter = None if g.mode == 'constrain' else frozenset(entry)
    for sid in derive.chain(arm, anchor):
        s = segs[sid]
        if s.walk is None:
            parts = None
        elif parts is not None:
            parts.append(s.walk if s.end_bp <= m else s.walk[:m - s.from_bp])
        for f, t, labs in segment_support(g, arm, s):
            if f >= m:
                break
            inter = labs if inter is None else inter & labs
            for l in inter:
                have[l] = min(t, m)
    return (None if parts is None else ''.join(parts), inter or frozenset(), have)


def restricted_leaves(arm, depth):
    """The leaves of the arm's DAG restricted to [0, depth) -- what a retrieval walked to
    that depth would end its walks at: every segment starting before the depth with no
    child starting before it (children start where their parent ends), i.e. the original
    leaves above the depth and every segment reaching or crossing it. Among the latter is
    a segment whose only child is a merge at or beyond the depth that it enters through a
    NON-first parent: no displayed walk of the deeper DAG passes through it, yet at the
    depth it is a walk of its own (bubbles merged at 36 and 62, compared at 30: comparing
    the deeper DAG's displayed walks cut at 30 would report differences that the merges
    beyond 30 cause)."""
    return [s.id for s in arm.segments
            if s.from_bp < depth and (s.leaf is not None or s.end_bp >= depth)]


def _cut_plan(a, depth):
    """The cuts of _cuts() for the DAG restricted to [0, depth) (depth > 0) and the tree one
    pass makes them in: {rows: [(leaf, k, m)] in restricted_leaves() order, m = min(its
    end, depth), k the segment holding base m - 1 (None for m = 0); keys: the distinct
    (k, m), in their first row's order; wanted: k -> its m's; marked: the segments a cut's
    chain passes through whole (each anchor's first-parent ancestors); kids: segment -> its
    first-parent children among the marked segments and the anchors, in id order; roots;
    dep: segment -> the length of its chain}. Structural: _price_cuts() reads the same
    plan."""
    segs = a.segments
    rows, keys, wanted = [], [], {}
    for leaf in restricted_leaves(a, depth):
        m = min(segs[leaf].end_bp, depth)
        if m == 0:
            k = None
        else:
            k = leaf                          # holds base m - 1 when it reaches the depth
            while segs[k].from_bp >= m:      # zero-length segments at the leaf's end
                k = segs[k].parents[0]
        rows.append((leaf, k, m))
        ms = wanted.setdefault(k, [])
        if m not in ms:
            ms.append(m)
            keys.append((k, m))
    marked = set()
    for k in wanted:
        if k is None:
            continue
        p = segs[k].parents
        x = p[0] if p else None
        while x is not None and x not in marked:
            marked.add(x)
            p = segs[x].parents
            x = p[0] if p else None
    kids, roots, dep = {}, [], {}
    # parents come first (segment ids): a node's first parent is marked, so already seen
    for x in sorted(marked.union(k for k in wanted if k is not None)):
        p = segs[x].parents
        if p:
            kids.setdefault(p[0], []).append(x)
            dep[x] = dep[p[0]] + 1
        else:
            roots.append(x)
            dep[x] = 1
    return {'rows': rows, 'keys': keys, 'wanted': wanted, 'marked': marked, 'kids': kids,
            'roots': roots, 'dep': dep}


def _cut_prices(g, a, plan, cost, held):
    """What _cuts() charges for |plan| (work model 2): each
    segment the pass reads whole once -- its own part of derive.cut_cost (its pieces of
    displayed support, their sets, the labels supporting it, its bases), with the
    transient of reading it (held, released) -- and per cut: its anchor's own part when it
    is read up to the cut only, the prefix joined from the pieces on the chain (a unit per
    8 pieces, and per 256 bases), its label set and its {label: j} built from the running
    state (a unit per 2 labels), the cut kept (its tuple, prefix, set and dict) -- and the
    pass's running state (the intersection, the dropped labels and their undo lists, the
    pieces and the stack), held for the pass -- not every cut its whole chain's cut_cost,
    which would re-read the chains per cut (quadratic on a comb)."""
    segs = a.segments
    n0 = len(segs[0].entry) if segs else 0
    nodes = list(plan['dep'])
    if g.mode == 'constrain':
        ops_ = derive.segment_ops(a)
        h_inter = max((len(segs[x].entry) + ops_[x] for x in nodes), default=0)
        h_have = n0 + h_inter
    else:
        h_inter = h_have = n0
    out_bytes = TUPLE + LIST + dict_bytes(h_have) + set_bytes(h_inter) + STR
    marked = plan['marked']
    dep = plan['dep']

    def own(x):
        p = segs[x].parents
        return cost[x] - (cost[p[0]] if p else 0)

    def cut(k, m):
        """-> (work, bytes kept, transient bytes released after) of cut (k, m)."""
        if k is None:
            return W_ELEM, out_bytes, 0
        part = 0 if k in marked else own(k)
        tr = 0 if k in marked else held[k]
        return (W_ELEM + part + (dep[k] >> 3) + (m >> 8) + ((n0 + h_have) >> 1),
                out_bytes + m, tr)
    pass_bytes = set_bytes(h_inter) + dict_bytes(h_have) + h_have * (INT + LIST_ITEM) \
        + list_bytes(len(nodes)) + len(nodes) * (LIST + TUPLE + 32)
    return {'seg': own, 'cut': cut, 'pass': pass_bytes}


def _cuts(g, side, depth, memo, cb=None):
    """[(leaf, m, cut info)] for every leaf of the arm's DAG restricted to [0, depth)
    (restricted_leaves), m = min(its end, depth). Leaves sharing the segment that holds
    base m - 1 share one cut, and all cuts are made in ONE parents-first pass over their
    chains (_cut_pass()): each segment's pieces are read once, not once per cut whose chain
    passes through it (on a comb of 4,000 walks that would re-read 16 M segment pieces and
    rewrite {label: j} for every supporting label at every piece). The cut infos are
    _cut_info()'s, value for value. Memoized per comparison in |memo|."""
    key = (id(g), side, depth)
    got = memo.get(key)
    if got is not None:
        return got
    a = g.arms[side]
    segs = a.segments
    out = []
    if cb is not None:
        # cut_cost reads segment_ops (constrain): charged with it
        derive.uses(cb, g, a, 'leaves', 'paths', 'cut_cost',
                    *(('segment_ops',) if g.mode == 'constrain' else ()))
        cost = derive.cut_cost(a, g.mode)
        held = _cut_transient(cb, g, a)
        # the restricted leaves' scan and the marking of their chains (work model 2)
        cb.charge((W_ELEM + W_STEP) * len(segs))
    if depth <= 0:
        # nothing is certified: one empty cut with the boundary labels, when there is a walk
        if derive.paths(a):
            out.append((None, 0, _cut_info(g, a, None, 0)))
        memo[key] = out
        return out
    plan = _cut_plan(a, depth)
    pr = _cut_prices(g, a, plan, cost, held) if cb is not None else None
    infos = _cut_pass(g, a, plan, cb, pr, held if cb is not None else None)
    for leaf, k, m in plan['rows']:
        if cb is not None:
            cb.charge(W_ELEM, TUPLE + 3 * LIST_ITEM)
        out.append((leaf, m, infos[(k, m)]))
    memo[key] = out
    return out


def _cut_steps(plan, segs):
    """The steps of the pass over plan's tree, in order (depth first, parents first, each
    node's children in id order): ('cut', k, m, partial) -- the cut (k, m), read up to m on
    the parent's state when |partial|, else from the state after k -- ('enter', x, None,
    None) -- segment x read whole -- and ('leave', x, None, None). _cuts() works through
    them; _price_cuts() charges them in the same order."""
    wanted, marked, kids = plan['wanted'], plan['marked'], plan['kids']
    for m in wanted.get(None, ()):
        yield 'cut', None, m, False
    stack = [(r, False) for r in reversed(plan['roots'])]
    while stack:
        x, leaving = stack.pop()
        if leaving:
            yield 'leave', x, None, None
            continue
        if x not in marked:
            # an anchor below which no cut lies: each of its cuts on its parent's state
            for m in wanted[x]:
                yield 'cut', x, m, True
            continue
        # a segment read whole: first the cuts ending inside it (none, as a rule: a
        # segment a cut's chain passes through ends before that cut's anchor, so the cut
        # of its own reaches its end), then its pieces, its cuts and the subtrees below
        end = segs[x].end_bp
        for m in wanted.get(x, ()):
            if m < end:
                yield 'cut', x, m, True
        yield 'enter', x, None, None
        for m in wanted.get(x, ()):
            if m >= end:
                yield 'cut', x, m, False
        stack.append((x, True))
        for c in reversed(kids.get(x, ())):
            stack.append((c, False))


def _cut_pass(g, a, plan, cb, pr, held):
    """The cut infos of plan['keys'] -> {(k, m): info}, in one depth-first pass over the
    first-parent tree of their chains (_cut_steps()): a running state -- the walking-order
    pieces, the intersection of the displayed support sets, the labels dropped from it
    with the end of their support -- is extended by each segment's pieces on the way down
    and restored on the way up (an undo per segment), so a segment is read once for all
    the cuts below it. A cut (k, m) is its anchor's pieces up to m on its parent's state.
    Each info is _cut_info(g, a, k, m)'s: (the prefix, None when a segment of the chain
    has no bases; the labels supporting all of [0, m); {label: j} -- 0 for a root entry
    label that never supported a piece, the end of the last piece a label supported
    before it dropped out, min(the last piece's end, m) for those supporting it all) --
    the {label: j} in another key order (its readers do not depend on it)."""
    segs = a.segments
    entry = segs[0].entry if segs else ()
    constrain = g.mode == 'constrain'
    st = {'inter': None if constrain else set(entry), 't': 0, 'none': 0}
    pieces = []
    dropped = {}

    def apply(x, m):
        """Extend the state by segment x's pieces up to m (None: all) -> its undo."""
        s = segs[x]
        t0 = st['t']
        init = False
        drops = []
        if s.walk is None:
            st['none'] += 1
        else:
            pieces.append(s.walk if m is None or s.end_bp <= m else s.walk[:m - s.from_bp])
        inter = st['inter']
        for f, t, labs in segment_support(g, a, s):
            if m is not None and f >= m:
                break
            if inter is None:
                inter = st['inter'] = set(labs)
                init = True
            else:
                d = inter - labs
                if d:
                    # each label leaves the intersection once on a chain: its support
                    # ends with the previous piece (0 before any piece)
                    for l in d:
                        dropped[l] = st['t']
                    inter -= d
                    drops.extend(d)
            st['t'] = t
        return (s, init, drops, t0)

    def undo(u):
        s, init, drops, t0 = u
        for l in drops:
            del dropped[l]
        if init:
            st['inter'] = None
        elif drops:
            st['inter'].update(drops)
        st['t'] = t0
        if s.walk is None:
            st['none'] -= 1
        else:
            pieces.pop()

    def info(m):
        inter = st['inter']
        have = {l: 0 for l in entry}
        have.update(dropped)
        if inter:
            tm = min(st['t'], m)
            for l in inter:
                have[l] = tm
        return (None if st['none'] else ''.join(pieces),
                frozenset(inter) if inter else frozenset(), have)

    out = {}
    undos = {}
    if cb is not None and plan['roots']:
        cb.charge(0, pr['pass'])
    for step, x, m, partial in _cut_steps(plan, segs):
        if step == 'cut':
            if cb is not None:
                w, keep, tr = pr['cut'](x, m)
                cb.charge(w, keep + tr)
            if x is None:
                got = _cut_info(g, a, None, m)
            elif partial:
                u = apply(x, m)
                got = info(m)
                undo(u)
            else:
                got = info(m)
            if cb is not None:
                cb.release(tr)
            out[(x, m)] = got
        elif step == 'enter':
            if cb is not None:
                cb.charge(pr['seg'](x), held[x])
            undos[x] = apply(x, None)
            if cb is not None:
                cb.release(held[x])
        else:
            undo(undos.pop(x))
    if cb is not None and plan['roots']:
        cb.release(pr['pass'])
    return out


def _route_parent(a, run, seg):
    """The parent of merge |seg| through which |run|'s lineage arrived (its incoming label
    in that parent's partition) -> (parent, exact). Not found: parents[0], not exact."""
    at = derive.lineage_label_at(a, run, seg.from_bp)
    for j, part in enumerate(derive.merge_parts(a, seg.id)):
        if at in part:
            return seg.parents[j], True
    return seg.parents[0], False


def _first_parent_at(a, pos):
    """seg -> the segment of its first-parent chain holding outward base |pos| (None for
    a segment ending at or before it), one pass parents first; the last position asked
    is kept in the arm's cache."""
    got = a.cache.get('first_parent_at')
    if got is not None and got[0] == pos:
        return got[1]
    segs = a.segments
    out = [None] * len(segs)
    for s in segs:
        if s.from_bp <= pos < s.end_bp:
            out[s.id] = s.id
        elif s.from_bp > pos and s.parents:
            out[s.id] = out[s.parents[0]]
    a.cache['first_parent_at'] = (pos, out)
    return out


def _route_segment_at(a, run, pos):
    """The segment of |run|'s OWN route (back from its anchor, through the partition
    holding the lineage at every merge, as routes() chooses) that holds outward base
    |pos| -> (segment, exact). Between merges the route is the first-parent chain, so the
    walk jumps from merge to merge (merge_above) and only the merges beyond |pos| are
    decided per lineage."""
    segs = a.segments
    above = derive.merge_above(a)
    at = _first_parent_at(a, pos)
    s = run.segment
    exact = True
    while True:
        m = above[s]
        if m is None or segs[m].from_bp <= pos:
            return (at[s] if at[s] is not None else s), exact
        s, ok = _route_parent(a, run, segs[m])
        exact = exact and ok


def _evidence_at(a, run, anchor, t):
    """derive.evidence() for |run| evaluated at |anchor| instead of its own anchor: the
    §5.1 merge walk over the first-parent chain root -> anchor, merges at or before |t|.
    For a run cut at the comparison depth, |anchor| is the segment of its own route at
    the depth, where the restricted DAG anchors it."""
    segs = a.segments
    above = derive.merge_above(a)
    route_from = 0
    m = above[anchor]
    while m is not None:
        seg = segs[m]
        if seg.from_bp <= t and derive.lineage_label_at(a, run, seg.from_bp) \
                not in derive.merge_parts(a, m)[0]:
            route_from = seg.from_bp
            break
        m = above[seg.parents[0]]
    return route_from, max(route_from, run.from_bp)


def _restricted_run_claim(g, a, run, depth, exact, inexact):
    """The claim of |run| in the DAG restricted to [0, depth): the run-start guard of a
    cut, then -- for a run reaching past the depth, closed by a merge AT the depth, or
    anchored on a segment starting there -- anchored at the segment of its own route that
    holds base depth - 1 (not at its deeper anchor, whose displayed chain may enter a later
    merge through another parent), with displayed support evaluated there. Everything at
    or beyond the depth is outside the restricted DAG: a merge at the depth is one a
    retrieval walked to the depth never makes, so a run it closes is an open claim there
    and a run anchored on it is anchored on its own route (otherwise a radius equal to a
    merge position would compare unequal with a deeper retrieval). Any other run
    ending within the depth is its uncut claim (its anchor and every merge on its chain lie
    below the depth)."""
    zero = run.from_bp == run.to_bp
    if zero:
        if run.from_bp > depth:
            return None
    elif run.from_bp >= depth:
        return None
    reanchor = run.to_bp > depth or (0 < depth == run.to_bp and (
        run.merged or a.segments[run.segment].from_bp >= depth))
    if not reanchor:
        return _run_claim(g, a, run, depth, exact)
    t, ok = _route_segment_at(a, run, depth - 1)
    if not ok and inexact is not None:
        inexact.append('%s arm, run %d (%s): its route through a merge is not in any '
                       'partition' % (a.side, run.id, g.labels[run.label].ref))
    route_from, ev_from = _evidence_at(a, run, t, depth)
    to_bp = min(run.to_bp, depth)
    if zero:
        kind, evidence_from = 'boundary', None
    elif ev_from < depth:
        kind, evidence_from = 'stretch', ev_from
    else:
        kind, evidence_from = 'route_only', None
    return Claim(
        label=g.labels[run.label], arm=a.side, path_id=_first_paths(a)[t], segment=t,
        from_bp=run.from_bp, to_bp=to_bp, route_bp=route_from, evidence_from=evidence_from,
        kind=kind, reason=None, qualifier=None, end_class='open', entered_by=run.entered_by,
        from_label=None if run.from_label is None else g.labels[run.from_label],
        cost=run.cost, loss=None, branches=None, support=g.support, exact=exact, run=run.id)


def _annotate_route_ends_at(arm, depth):
    """_annotate_route_ends() of the DAG restricted to [0, depth): the union rule over
    the segments starting before the depth, P runs cut there, and a label alive at a
    segment's restricted end with no kept child recording it ends there (at the depth for
    a segment reaching it). -> {label: [(seg, j)]}."""
    segs = arm.segments
    # the whole DAG is the restricted one only when nothing starts at or beyond the depth:
    # a merge AT the depth (zero length, ending there) is outside it
    if all(s.end_bp <= depth and (s.from_bp < depth or not s.parents) for s in segs):
        return _annotate_route_ends(arm)[1]
    return _route_ends_pass(arm, depth)[1]


def _restricted_claims(g, side, depth, inexact=None, cb=None):
    """The claims of one arm's DAG restricted to [0, depth): what a retrieval walked to
    the depth would claim, and so what a comparison at the depth must key -- not the
    deeper DAG's claims cut at the depth (claims(at_most_bp)), whose displayed chains run
    through merges beyond the depth: a lineage that enters such a merge through a
    non-first parent is route_only there, on the other branch's bases, while the shallower
    retrieval shows it as a stretch on its own branch. Constrain: _restricted_run_claim
    per run. Annotate: the maximal route ends of the restricted DAG; their displayed
    support is the claims' own (every merge on an anchor's chain lies below the depth).
    |inexact| collects what could not be reconstructed exactly."""
    a = g.arms[side]
    exact = arm_exact(g, side)
    out = []
    if g.mode == 'constrain':
        if cb is not None:
            derive.uses(cb, g, a, 'leaves', 'paths', 'path_of_leaf', 'merge_above',
                        'merge_depth', 'first_paths', 'partition_sets')
            cb.uses((id(a), 'first_parent_at', depth - 1), *derive.price(
                a, 'first_parent_at', len(g.labels), g.mode))
            # every run is keyed (a comparison has no partial answer)
            _charge_all(cb, _restricted_claim_prices(cb, g, a), _CLAIM_BYTES)
        for run in a.runs:
            c = _restricted_run_claim(g, a, run, depth, exact, inexact)
            if c is not None:
                out.append(c)
        return out
    if cb is not None:
        derive.uses(cb, g, a, 'leaves', 'paths', 'path_of_leaf', 'merge_above', 'merge_depth',
                    'first_paths', 'annotate_route_ends')
        w, m = derive.price(a, 'annotate_route_ends', len(g.labels), g.mode)
        cb.uses((id(a), 'route_ends_at', depth), w, m)
        md = derive.merge_depth(a)
    ends_at = _annotate_route_ends_at(a, depth)
    if cb is not None:
        n = sum(len(v) for v in ends_at.values())
        cb.charge(sort_work(len(ends_at)) + n * (W_ROW + 2)
                  + sum(3 * md[x] for v in ends_at.values() for x, _ in v),
                  n * (_CLAIM_BYTES + TUPLE) + list_bytes(len(ends_at)))
    for l, es in sorted(ends_at.items()):
        for anchor, j in es:
            c = _annotate_claim(g, a, l, anchor, j, _annotate_evidence(a, l, anchor), depth,
                                exact)
            if c is None:
                continue
            if j >= depth and a.segments[anchor].end_bp > depth:
                # the restricted DAG ends the route at the depth, where the deeper one
                # goes on: censored there (open), not lost
                c.reason, c.end_class = None, 'open'
                if c.evidence_from is not None and c.evidence_from < depth:
                    c.kind = 'stretch'
                else:
                    c.kind, c.evidence_from = 'route_only', None
            out.append(c)
    return out


def _walk_keys(g, sides, depth, sel, memo=None, cb=None):
    memo = {} if memo is None else memo
    out = {}
    for side in sides:
        seen = set()
        for p, m, info in _cuts(g, side, depth, memo, cb):
            if cb is not None:
                n = len(info[1])
                cb.charge(W_ELEM + (0 if id(info) in seen else
                                    W_ELEM * n + (m >> 8)),
                          0 if id(info) in seen else
                          STR + m + TUPLE + set_bytes(n) + n * (STR + 24))
            if id(info) in seen:
                # walks sharing a cut share its key and its labels: one update per cut,
                # not one per walk (a shallow cut of a wide arm is shared by thousands
                # of walks, each with every seed label)
                continue
            seen.add(id(info))
            prefix, ids, _ = info
            seq = _natural(side, _cut_prefix(prefix))
            refs = {g.labels[l].ref for l in ids}
            if sel is not None:
                refs &= sel
                if not refs:
                    # restricted to labels, the comparison is over the walks UNDER them:
                    # a walk none of them supports is in neither side's set (as for the
                    # claims of a label filter), not a walk with no labels
                    continue
            out.setdefault((side, seq), set()).update(refs)
    return out


def _label_keys(g, sides, depth, sel, cb=None):
    """(side, ref) -> {direct_bp, reach_bp} of every label the arm holds within [0, depth),
    each summary DERIVED from the records clipped to the depth (_clipped_summary), not a
    completed summary clamped to it: a label first recorded beyond the depth is in
    neither retrieval's [0, depth), and clamping its reach_bp would report it in the
    deeper one only."""
    out = {}
    for side in sides:
        got = _clipped_summary(g, g.arms[side], depth, cb)
        if cb is not None:
            cb.charge(W_ELEM * 2 * len(got), len(got) * (TUPLE + 2 * STR + 24 + DICT_KEY))
        for lid, v in got.items():
            ref = g.labels[lid].ref
            if sel is not None and ref not in sel:
                continue
            out[(side, ref)] = v
    return out


def _run_exists_at(run, depth):
    """The run-start guard of a claim at a cut: a run starting at or after the cut does
    not exist yet; a zero-length run (a boundary claim) exists at its position."""
    if run.from_bp == run.to_bp:
        return run.from_bp <= depth
    return run.from_bp < depth


def _clipped_summary(g, arm, depth, cb=None):
    """label id -> {direct_bp, reach_bp} of one arm's label_summary (§5.1, both modes)
    computed over [0, depth) only, for the labels the arm holds there. Constrain: the
    runs that exist at the cut, each cut to it (direct_bp: the seed runs from 0; reach_bp:
    credited to the lineage root). Annotate: the P runs starting before the cut, each cut
    to it, and the union pass of direct_bp stopped at the cut; a label recorded at the
    seed boundary only has both 0. At a depth beyond the arm it is label_summary itself
    (for the labels it holds)."""
    out = {}
    if cb is not None:
        derive.uses(cb, g, arm)
        z = derive.arm_sizes(arm)
        n = min(len(g.labels), z['entry_ids'] + z['runs'] + z['presence_ids'])
        if g.mode == 'constrain':
            cb.charge(6 * z['runs'], list_bytes(z['runs']) + n * (dict_bytes(2) + 2 * INT))
        else:
            cb.charge(W_ELEM * (3 * z['segments'] + z['presence'])
                      + ((2 * z['presence_ids'] + 2 * z['entry_ids']) >> 2),
                      list_bytes(z['segments']) + z['segments'] * SET
                      + z['entry_ids'] * SET_ITEM + n * (dict_bytes(2) + 2 * INT))
    if g.mode == 'constrain':
        root = [0] * len(arm.runs)
        for r in arm.runs:
            root[r.id] = root[r.prev_run] if (r.from_label is not None
                                              and r.prev_run is not None) else r.label
        for r in arm.runs:
            if not _run_exists_at(r, depth):
                continue
            to = min(r.to_bp, depth)
            own = out.setdefault(r.label, {'direct_bp': 0, 'reach_bp': 0})
            if r.from_label is None and r.from_bp == 0:
                own['direct_bp'] = max(own['direct_bp'], to)
            lin = out.setdefault(root[r.id], {'direct_bp': 0, 'reach_bp': 0})
            lin['reach_bp'] = max(lin['reach_bp'], to)
        return out
    segs = arm.segments
    if segs:
        for l in segs[0].entry:
            out.setdefault(l, {'direct_bp': 0, 'reach_bp': 0})
    alive_end = [None] * len(segs)
    for s in segs:
        for p in s.presence:
            if p.from_bp >= depth:
                break
            to = min(p.to_bp, depth)
            for l in p.labels:
                v = out.setdefault(l, {'direct_bp': 0, 'reach_bp': 0})
                if to > v['reach_bp']:
                    v['reach_bp'] = to
        if s.parents:
            alive = set()
            for p in s.parents:
                alive.update(alive_end[p])
        else:
            alive = set(s.entry)
        if s.from_bp < depth:
            for p in s.presence:
                if not alive or p.from_bp >= depth:
                    break
                still = alive.intersection(p.labels)
                for l in alive - still:
                    v = out[l]
                    if p.from_bp > v['direct_bp']:
                        v['direct_bp'] = p.from_bp
                alive = still
            end = min(s.end_bp, depth)
            for l in alive:
                v = out[l]
                if end > v['direct_bp']:
                    v['direct_bp'] = end
        else:
            alive = set()
        alive_end[s.id] = alive
    return out


def _walk_prefix(arm, seg_id, hi):
    """Outward bases [0, hi) along seg's first-parent chain in WALKING order. Prefix tests
    read walks outward on both arms: the natural spelling of the left arm is the reverse
    of the walking order, where an outward prefix is a suffix."""
    return derive.walk_bases(arm, seg_id)[:hi]


def _natural(side, w):
    return w if side == 'right' else w[::-1]


def _cut_prefix(prefix):
    if prefix is None:
        raise ValueError('this retrieval carries no bases (output.sequences: false)')
    return prefix


def _route_pairs(g, sides, depth, memo=None, cb=None):
    """(side, walk prefix cut at depth in walking order, ref) of every walk under each
    label with displayed support over the whole cut walk."""
    memo = {} if memo is None else memo
    out = set()
    for side in sides:
        seen = set()
        for p, m, info in _cuts(g, side, depth, memo, cb):
            if id(info) in seen:
                continue                     # walks sharing the cut give the same pairs
            seen.add(id(info))
            if cb is not None:
                n = len(info[1])
                cb.charge(W_ELEM * (1 + 2 * n), n * (TUPLE + 24 + STR + 16 + SET_ITEM))
            seq = _cut_prefix(info[0])
            for l in info[1]:
                out.add((side, seq, g.labels[l].ref))
    return out


def _supported_prefixes(g, sides, depth, memo=None, cb=None):
    """(side, ref) -> the walking-order prefixes w[:j] of every walk w (cut at depth) on
    which the label has displayed support over all of [0, j), j maximal per walk: a
    label that ends inside a walk (at a node where the walk goes on with other labels)
    still supports the walk up to its end, which whole-walk pairs alone never see."""
    memo = {} if memo is None else memo
    out = {}
    for side in sides:
        seen = set()
        for p, m, info in _cuts(g, side, depth, memo, cb):
            if id(info) in seen:
                continue
            seen.add(id(info))
            if cb is not None:
                # one prefix string per label the cut supports (the quadratic part of
                # prefix_subset: a label per walk per base of its support)
                n = len(info[2])
                j = sum(info[2].values())
                cb.charge(W_ELEM * (1 + 2 * n) + (j >> 8),
                          n * (STR + SET_ITEM + TUPLE + 24) + j)
            w = _cut_prefix(info[0])
            for l, j in info[2].items():
                out.setdefault((side, g.labels[l].ref), set()).add(w[:j])
    return out


def _recorded_refusal(g, arm, seq, ref, memo=None, cb=None):
    """-> (at_bp, reason) of a successor that |g| recorded as not taken by label |ref|
    where walk |seq| (walking order) leaves g's trie: a V refusal (quorum, split limit), a
    blocked successor or a skipped hairpin, on a segment whose chain spells seq up to
    the branch and for the base seq has there. None when nothing was recorded. |memo|:
    the comparison's, where the index of the arm's refusals is kept (one omission after
    another asks about the same ones). |cb| (local limits): each segment followed and each
    recorded refusal tested is charged before it is.

    The chains are not spelled: seq is followed down the tree of the refusals' chains
    (_refusal_index()), and a refusal at |at| is on seq's path when the chain matches seq
    up to |at| -- read at the segment that holds base at - 1 of the chain (its anchor), so
    only the refusals anchored on the segments seq reaches are tested, not every refusal of
    the arm on every call (an arm of 32,001 segments scanned per omission); the first
    recorded one in the arm's order (branch events, then segments' events) wins."""
    ids = _label_index(g)[0].get(ref)
    if not ids:
        return None
    lid = ids[0]
    key = ('refusal_index', id(g), arm.side)
    index = None if memo is None else memo.get(key)
    if index is None:
        index = _refusal_index(arm, cb)
        if memo is not None:
            memo[key] = index
    if index is False:
        # a refusal's chain lacks bases: no spelling matches seq (a retrieval without
        # bases is not compared by its walks)
        return None
    kids, roots, by_anchor = index
    if not by_anchor:
        return None
    matched = _follow(arm.segments, kids, roots, seq, len(seq), cb)
    best = None                      # (priority, at_bp, reason)
    for s, m in matched.items():
        for prio, at, kind, x in by_anchor.get(s, ()):
            if best is not None and prio >= best[0]:
                break                # in priority order: nothing later here can win
            if cb is not None:
                cb.charge(W_ELEM)
            if at >= len(seq) or m < at:
                continue
            c = seq[at]
            if kind == 'V':
                for r in x.refused:
                    if cb is not None:
                        cb.charge(W_ELEM + (len(r.labels) >> ID_SHIFT_C))
                    if r.char == c and lid in r.labels:
                        best = (prio, at, r.cause)
                        break
            else:
                if cb is not None:
                    cb.charge(W_ELEM + (len(x.labels) >> ID_SHIFT_C))
                if x.char == c and lid in x.labels:
                    best = (prio, at, REASON[x.reason] if x.type == 'blocked' else 'hairpin')
            if best is not None and best[0] == prio:
                break
    if cb is not None:
        cb.release(_follow_bytes(len(matched)))
    return None if best is None else best[1:]


def _refusal_index(arm, cb=None):
    """(kids, roots, {anchor segment: [(priority, at_bp, kind, event)] in priority order})
    over the refusals _recorded_refusal() reads (V refusals, blocked successors, skipped
    hairpins) -- the first-parent tree of their chains (_refusal_tree()) and each refusal
    filed under the segment of its chain that holds base at_bp - 1 (the chain's first
    segment for at_bp 0: an event before its own segment is filed under the ancestor that
    holds it) -- or False when one of those chains has a segment without bases. Priority:
    the arm's order, branch events first. A refusal past its segment's end is on no
    chain and is left out. |cb| (local limits): the search for the ancestors of events before
    their segment (_chain_anchors()), where there are any; the rest is charged by the
    caller."""
    tree = _refusal_tree(arm)
    if tree is False:
        return False
    kids, offs, roots = tree
    segs = arm.segments
    by_anchor = {}
    before = []                      # (segment, at_bp, entry) of events before their segment
    prio = 0

    def add(seg, at, kind, x):
        nonlocal prio
        prio += 1
        if offs[seg] + len(segs[seg].walk) < at:
            return
        if offs[seg] <= at:
            by_anchor.setdefault(seg, []).append((prio, at, kind, x))
        else:
            before.append((seg, at, (prio, at, kind, x)))

    for be in arm.branch_events:
        add(be.segment, be.at_bp, 'V', be)
    for s in segs:
        for ev in s.events:
            if ev.type == 'blocked' or (ev.type == 'hairpin' and not ev.followed):
                add(s.id, ev.at_bp, 'E', ev)
    if before:
        # An event whose at_bp lies before its own segment (the parser does not check it,
        # the server never writes one) is filed under an ancestor. One walk of the tree finds
        # them all, charged -- a climb of each event's chain would take O(events x depth) on a
        # deep comb -- and the lists they join are put back in the arm's order
        got = _chain_anchors(kids, roots, offs, before, cb)
        joined = set()
        for (_, _, entry), a in zip(before, got):
            by_anchor.setdefault(a, []).append(entry)
            joined.add(a)
        if cb is not None:
            cb.charge(sort_work(sum(len(by_anchor[a]) for a in joined)))
        for a in joined:
            by_anchor[a].sort()      # by priority: unique, so nothing else is compared
    return kids, roots, by_anchor


def _refusal_tree(arm):
    """({segment: its first-parent children}, {segment: the length of its chain's bases
    before it}, [roots]) over the chains of the segments _recorded_refusal() reads (branch
    events, blocked successors, skipped hairpins), or False when one of those chains has a
    segment without bases."""
    segs = arm.segments
    return _chain_tree(segs, {be.segment for be in arm.branch_events}
                       | {s.id for s in segs for ev in s.events
                          if ev.type == 'blocked' or (ev.type == 'hairpin' and not ev.followed)})


def _chain_tree(segs, targets):
    """({segment: its first-parent children}, {segment: the length of its chain's bases
    before it}, [roots]) over the first-parent chains of |targets|, or False when one of
    those chains has a segment without bases."""
    on = set()
    for x in targets:
        while x not in on:
            on.add(x)
            if segs[x].walk is None:
                return False
            p = segs[x].parents
            if not p:
                break
            x = p[0]
    kids, offs, roots = {}, {}, []
    for x in sorted(on):                 # parents first
        p = segs[x].parents
        if p:
            q = p[0]
            offs[x] = offs[q] + len(segs[q].walk)
            kids.setdefault(q, []).append(x)
        else:
            offs[x] = 0
            roots.append(x)
    return kids, offs, roots


def _chain_anchors(kids, roots, offs, asks, cb=None):
    """[the deepest segment of the chain root -> |segment| whose offset is at most |at|, for
    each (segment, at, ...) of |asks|] (at >= 0, so a root at least) -- the segment of the
    tree (kids, offs, roots: _chain_tree()) that holds base |at| of that chain, or the
    chain's last segment that starts at it. One walk of the tree down from its roots, the
    chain to the current segment kept with its offsets, and one bisection per ask: O(tree +
    asks x log depth), where a climb per ask would be O(asks x depth). |cb| (local limits):
    a step per segment of the tree and a bisection per ask, and the walk's stack, chain and
    answer, charged before the walk (the stack and the chain are dropped after it)."""
    n = len(offs)
    by_seg = {}
    for i, (seg, at) in enumerate((q[0], q[1]) for q in asks):
        by_seg.setdefault(seg, []).append((i, at))
    if cb is not None:
        held = 3 * list_bytes(n) + n * (TUPLE + 2 * 8)
        cb.charge(W_STEP * n + W_ELEM * len(asks) * max(1, n.bit_length()),
                  held + list_bytes(len(asks)) + dict_bytes(len(by_seg))
                  + len(asks) * (TUPLE + 2 * 8 + LIST_ITEM))
    out = [None] * len(asks)
    chain_off, chain_seg = [], []
    stack = [(r, 0) for r in reversed(roots)]
    while stack:
        s, d = stack.pop()
        del chain_off[d:], chain_seg[d:]
        chain_off.append(offs[s])
        chain_seg.append(s)
        for i, at in by_seg.get(s, ()):
            j = bisect.bisect_right(chain_off, at) - 1
            out[i] = chain_seg[j if j > 0 else 0]
        ks = kids.get(s)
        if ks:
            stack.extend((c, d + 1) for c in reversed(ks))
    if cb is not None:
        cb.release(held + dict_bytes(len(by_seg)) + len(asks) * (TUPLE + 2 * 8 + LIST_ITEM))
    return out


def _follow_bytes(n):
    """What _follow() holds per segment it followed (its entry in the answer)."""
    return n * (DICT_KEY + 2 * INT)


def _follow(segs, kids, roots, seq, lim, cb=None):
    """{segment: the length of seq's prefix its chain spells, its own bases included} for
    every segment of the tree (kids, roots) whose chain BEFORE it is seq's (bases [0,
    lim)): seq followed down the tree, a segment's bases compared only where its parent
    matched in full. |cb| (local limits): each segment charged before it is compared -- its
    bases (two slices of them, dropped at once) and its entry in the answer, which the
    caller releases (_follow_bytes()) once it drops the answer."""
    out = {}
    stack = [(r, 0) for r in roots]
    while stack:
        s, off = stack.pop()
        w = segs[s].walk
        if cb is not None:
            hi = max(0, min(len(w), lim - off))
            cb.charge(W_ELEM + (hi >> BASE_SHIFT), 2 * str_bytes(hi) + _follow_bytes(1))
            cb.release(2 * str_bytes(hi))
        n = _match_len(w, seq, off, lim)
        out[s] = off + n
        if n == len(w) and off + n < lim:
            stack.extend((c, off + n) for c in kids.get(s, ()))
    return out


def _match_len(w, seq, off, lim):
    """How many leading bases of |w| equal seq[off:], at most lim - off."""
    hi = max(0, min(len(w), lim - off))
    if w[:hi] == seq[off:off + hi]:
        return hi
    n = 0
    while w[n] == seq[off + n]:
        n += 1
    return n


def _divergence(arm, seq, depth, memo=None):
    """The length of the longest prefix of |seq| (walking order) that some walk of the
    arm's DAG restricted to [0, depth) spells: where seq leaves the arm's trie there.
    |memo|: the comparison's, where the restricted walks' first-parent tree is kept.

    The walks are not spelled: seq is followed down the first-parent tree of the
    restricted walks (_restricted_tree()), a segment's bases compared where its chain so
    far matches, and the deepest match is the longest common prefix with some walk cut at
    the depth. Keeping every restricted walk's spelling for the comparison (one omission
    after another asks again) would hold the whole arm's text spelled out -- on a comb as
    much as to_fasta() writes -- and spelling them per omission would be quadratic in
    time."""
    key = ('divergence_tree', id(arm), depth)
    tree = None if memo is None else memo.get(key)
    if tree is None:
        tree = _restricted_tree(arm, depth)
        if memo is not None:
            memo[key] = tree
    if tree is False:
        return None
    kids, roots = tree
    # the deepest match: every walk below a followed segment spells its matched prefix
    return max(_follow(arm.segments, kids, roots, seq, min(depth, len(seq))).values(),
               default=0)


def _restricted_tree(arm, depth):
    """({segment: its first-parent children whose subtrees hold a walk of the DAG
    restricted to [0, depth)}, [such roots]) -- the tree _divergence() follows a sequence
    down, the restricted leaves' chains and nothing else -- or False when one of those
    chains has a segment without bases (the walks cannot be spelled: no divergence)."""
    segs = arm.segments
    on = set(restricted_leaves(arm, depth))
    # parents come first (segment ids), so one pass from the last segment marks every
    # chain of a restricted leaf
    for s in reversed(segs):
        if s.id in on and s.parents:
            on.add(s.parents[0])
    kids = {}
    roots = []
    for s in segs:
        if s.id not in on:
            continue
        if s.walk is None:
            return False
        if s.parents:
            kids.setdefault(s.parents[0], []).append(s.id)
        else:
            roots.append(s.id)
    return kids, roots


def _claim_index(arm, rows, cb=None):
    """The claims |rows| of one arm of |a| as prefix_subset() looks them up: (kids, roots,
    {ref: {anchor: ([to_bp ascending], [(to_bp, row, claim) of the first row at it])}},
    {ref: (0, row, claim) of the first claim at to_bp 0}) -- the first-parent tree of the
    claims' chains (_chain_tree()) and each claim filed under the segment of its chain that
    holds its last base (to_bp - 1). A claim's prefix -- its chain's bases [0, to_bp) -- is
    a prefix of a walk exactly when the walk, followed down the tree (_follow()), matches
    that segment's chain at least to to_bp. A ValueError (derive.NO_BASES) where a chain
    lacks bases: every walk of |a| holds them here (the caller compares no side without).
    |cb| (local limits): the search for the anchors of claims whose to_bp lies before their
    segment (_chain_anchors()), where there are any; the rest is charged by the caller."""
    segs = arm.segments
    tree = _chain_tree(segs, {c.segment for c in rows})
    if tree is False:
        raise ValueError(derive.NO_BASES)
    kids, offs, roots = tree
    by_ref, zero = {}, {}
    before = []                      # (segment, to_bp - 1, row) of claims before their segment
    for i, c in enumerate(rows):
        ref = c.label.ref
        if c.to_bp <= 0:
            zero.setdefault(ref, (0, i, c))
            continue
        x = c.segment
        if offs[x] >= c.to_bp:
            # its last base before its own segment: an ancestor holds it, found below by
            # one walk of the tree (a climb per claim would be O(claims x depth), uncharged)
            before.append((x, c.to_bp - 1, i))
            continue
        at, first = by_ref.setdefault(ref, {}).setdefault(x, ([], {}))
        if c.to_bp not in first:
            first[c.to_bp] = (c.to_bp, i, c)       # rows in order: the first one wins
    if before:
        for (_, _, i), x in zip(before, _chain_anchors(kids, roots, offs, before, cb)):
            c = rows[i]
            at, first = by_ref.setdefault(c.label.ref, {}).setdefault(x, ([], {}))
            cur = first.get(c.to_bp)
            if cur is None or i < cur[1]:
                first[c.to_bp] = (c.to_bp, i, c)   # filed late: the first row still wins
    for anchors in by_ref.values():
        for x, (at, first) in anchors.items():
            ks = sorted(first)
            anchors[x] = (ks, [first[k] for k in ks])
    return kids, roots, by_ref, zero


def _best_claim(index, matched, ref, cb=None):
    """-> (to_bp, row, claim) of the claim of |ref| whose prefix is the longest prefix of
    the walk that |matched| (_follow() over the index's tree) followed -- the first in row
    order among equals -- or None: what the scan of every claim of the label found
    (seq.startswith(prefix), the longest, the first), by following the walk's chain."""
    _, _, by_ref, zero = index
    best = zero.get(ref)
    anchors = by_ref.get(ref)
    if anchors:
        # the smaller of the two sides is iterated: the label's anchors or the segments
        # the walk reached
        pairs = ((x, matched[x]) for x in anchors if x in matched) \
            if len(anchors) <= len(matched) else \
            ((x, m) for x, m in matched.items() if x in anchors)
        if cb is not None:
            cb.charge(W_ELEM * min(len(anchors), len(matched)))
        for x, m in pairs:
            ks, entries = anchors[x]
            if cb is not None:
                cb.charge(W_ELEM * max(1, len(ks).bit_length()))
            j = bisect.bisect_right(ks, m) - 1
            if j >= 0:
                e = entries[j]
                if best is None or e[0] > best[0] or (e[0] == best[0] and e[1] < best[1]):
                    best = e
    return best


def _prefix_subset(a, b, sides, depth, sel, memo=None, inexact=None, cb=None, state=None):
    memo = {} if memo is None else memo
    state = {} if state is None else state
    if cb is not None:
        cb.phase = 'keys:a'
    pa = _route_pairs(a, sides, depth, memo, cb)
    state['keys_a'] = len(pa)
    if cb is not None:
        cb.phase = 'keys:b'
    pb = _route_pairs(b, sides, depth, memo, cb)
    state['keys_b'] = len(pb)
    if sel is not None:
        pa = {x for x in pa if x[2] in sel}
        pb = {x for x in pb if x[2] in sel}
    # a's walk under a label must be a prefix of a walk of b on which b's label supports
    # at least that prefix -- not necessarily the whole walk (it may end inside it)
    supported = _supported_prefixes(b, sides, depth, memo, cb)
    if cb is not None:
        cb.phase = 'prefix_scan'
        cb.charge(sort_work(len(pa)) + sort_work(len(pb)), list_bytes(len(pa) + len(pb)))
    violations = []
    ordered = {}
    for side, seq, ref in sorted(pa):
        cand = ordered.get((side, ref))
        if cand is None:
            # the label's supported prefixes in order, once: the strings beginning with seq
            # are a run of that order starting where seq would be inserted, so one bisection
            # and one startswith decide (testing every prefix with startswith per walk would
            # make 40M calls on a wide pair)
            got = supported.get((side, ref), ())
            if cb is not None:
                cb.charge(sort_work(len(got)) * (1 + (depth >> 9)), list_bytes(len(got)))
            cand = ordered[(side, ref)] = sorted(got)
        if cb is not None:
            steps = 1 + max(1, len(cand).bit_length())
            cb.charge(W_ELEM + steps * (1 + (len(seq) >> 9)), dict_bytes(3) + STR + len(seq))
        i = bisect.bisect_left(cand, seq)
        if not (i < len(cand) and cand[i].startswith(seq)):
            violations.append({'arm': side, 'walk': _natural(side, seq), 'ref': ref})
    omissions = []
    claim_index = {}
    for side in sides:
        arm = a.arms[side]
        rows = [c for c in _restricted_claims(a, side, depth, inexact, cb)
                if c.kind != 'merged']
        if rows and cb is not None:
            # the tree of the claims' chains (no bases: the walks are followed down it, not
            # spelled) and per claim its entry under its label and anchor, the anchors'
            # sorted lists and the per-label dicts (a claim each at most)
            _charge_chain_tree(cb, arm, {c.segment for c in rows})
            cb.charge(W_ELEM * 3 * len(rows) + sort_work(len(rows)),
                      len(rows) * (TUPLE + 3 * 8 + 2 * LIST_ITEM + 2 * DICT_KEY + 2 * LIST
                                   + INT + TUPLE + 2 * 8)
                      + len({c.label.ref for c in rows}) * (DICT + DICT_KEY + TUPLE))
        claim_index[side] = _claim_index(arm, rows, cb) if rows else None
    if cb is not None:
        cb.phase = 'omissions'
        n_refusal = {s: len(a.arms[s].branch_events)
                     + sum(len(x.events) for x in a.arms[s].segments) for s in sides}
        n_leaves = {s: len(restricted_leaves(a.arms[s], depth)) for s in sides}
        cb.charge(W_ELEM * sum(len(a.arms[s].segments) for s in sides))
    followed = None                 # (side, seq, the claims' tree followed by seq)
    for side, seq, ref in sorted(pb):
        if (side, seq, ref) in pa:
            continue
        if cb is not None:
            cb.charge(W_ROW, dict_bytes(5) + STR + len(seq))
        index = claim_index.get(side)
        best = None
        if index is not None:
            if followed is None or followed[:2] != (side, seq):
                # pb is in (side, seq, ref) order: one walk is followed once for all its
                # labels, and the answer of the walk before is dropped first
                if cb is not None and followed is not None:
                    cb.release(_follow_bytes(len(followed[2])))
                followed = (side, seq, _follow(a.arms[side].segments, index[0], index[1],
                                               seq, len(seq), cb))
            best = _best_claim(index, followed[2], ref, cb)
        if best is not None and best[0] >= len(seq):
            continue                 # a covers the whole walk under this label
        entry = {'arm': side, 'walk': _natural(side, seq), 'ref': ref}
        arm = a.arms[side]
        ev_to = arm.evidence_complete_to_bp
        if best is None:
            if cb is not None:
                # the refusal is looked up by the label's id: a's name/ref index, charged
                # (the cross pair builds it)
                _uses_g(cb, a, 'label_index')
                # a's recorded refusals indexed once per comparison (their chains' tree and
                # each refusal under its anchor), then per omission the segments seq
                # reaches and the refusals anchored there (_recorded_refusal() charges them)
                if cb.first((id(arm), 'refusal_spellings')):
                    _charge_chain_tree(cb, arm, {be.segment for be in arm.branch_events}
                                       | {x.id for x in arm.segments if x.events})
                    cb.charge(W_ELEM * n_refusal[side],
                              n_refusal[side] * (TUPLE + 4 * 8 + LIST_ITEM + DICT_KEY + LIST))
            # the label went on along another successor in a: the reason, when there
            # is one, is a recorded refusal of b's successor (V, or a blocked or skipped
            # E), not a claim
            got = _recorded_refusal(a, arm, seq, ref, memo, cb)
            if got is not None:
                entry['at_bp'], entry['reason'] = got
            else:
                if cb is not None:
                    # the divergence of every restricted walk of a (a base at a time), the
                    # tree of their chains once per comparison
                    if cb.first((id(arm), 'divergence_spellings', depth)):
                        _charge_chain_tree(cb, arm, restricted_leaves(arm, depth))
                    cb.charge(n_leaves[side] * (1 + (min(depth, len(seq)) >> 1)))
                div = _divergence(arm, seq, depth, memo)
                if ev_to is not None and div is not None and div >= ev_to:
                    # beyond the evidence boundary: the refusal may have been dropped
                    entry['at_bp'], entry['reason'] = div, 'unexplained_capped'
                else:
                    entry['reason'] = 'unexplained'
        else:
            _, _, c = best
            entry['at_bp'] = c.to_bp
            entry['reason'] = c.qualifier or c.reason or c.kind
            if _open_at(c, depth):
                entry['reason'] = 'open'
            elif ev_to is not None and c.to_bp >= ev_to and c.reason is None:
                entry['reason'] = 'unexplained_capped'
        omissions.append(entry)
    if cb is not None and followed is not None:
        cb.release(_follow_bytes(len(followed[2])))
    return violations, omissions


# ------------------------------------------------------------------ size, summary

# the model's record types: a cache refers to them, it never copies them
_MODEL_TYPES = None


def _model_types():
    global _MODEL_TYPES
    if _MODEL_TYPES is None:
        from . import model as m
        _MODEL_TYPES = frozenset((
            m.Graphlet, m.Arm, m.Segment, m.Run, m.Label, m.SeedInfo, m.Dropped, m.Event,
            m.PresenceRun, m.LeafContinuation, m.Leaf, m.GrowthBin, m.Refusal,
            m.BranchEvent, m.CapTrigger, m.Frontier, m.LabelsPerNode, m.Counts,
            m.Limitation, m.Outcome, m.ResourceStop))
    return _MODEL_TYPES


_SLOTS = {}


def _slots_of(t):
    got = _SLOTS.get(t)
    if got is None:
        names = []
        for k in t.__mro__:
            for n in getattr(k, '__slots__', ()):
                if n not in ('__dict__', '__weakref__') and n not in names:
                    names.append(n)
        got = _SLOTS[t] = tuple(names)
    return got


_SHARED = None


def _shared_ids():
    """The ids of the program's own constants a model refers to without owning them: the
    strings, tuples and frozensets of the package modules' globals (code tables such as
    REASON, the parser's shared star sets) and of their functions' code constants ('right',
    'constrain', 'label_lost', cache keys, ...). They exist whether or not a model does,
    so charging them overstated a small model by a fifth (they are module-lifetime
    objects: their ids stay valid)."""
    global _SHARED
    if _SHARED is None:
        import types
        from . import _codec, derive as d, export, model, parser
        ids = set()
        stack = []

        def code(c):
            for k in c.co_consts:
                if isinstance(k, types.CodeType):
                    code(k)
                else:
                    stack.append(k)
        for m in (_codec, parser, model, d, export, sys.modules[__name__]):
            for v in vars(m).values():
                if isinstance(v, types.FunctionType):
                    code(v.__code__)
                elif isinstance(v, type) and v.__module__ == m.__name__:
                    for x in vars(v).values():
                        f = getattr(x, '__func__', x)
                        if isinstance(f, types.FunctionType):
                            code(f.__code__)
                        elif isinstance(x, property) and x.fget is not None:
                            code(x.fget.__code__)
                elif not isinstance(v, types.ModuleType):
                    stack.append(v)
        while stack:
            x = stack.pop()
            t = type(x)
            if t is str or t is tuple or t is frozenset or t is dict:
                if id(x) in ids:
                    continue
                if t is not dict:            # a module's dict is mutable: entered only
                    ids.add(id(x))
                if t is dict:
                    stack.extend(x.keys())
                    stack.extend(x.values())
                elif t is not str:
                    stack.extend(x)
        _SHARED = frozenset(ids)
    return _SHARED


def _deep_bytes(roots, stop=None):
    """sys.getsizeof over everything reachable from |roots|, each object once: what the
    objects hold in the heap. Not counted: None, booleans, small ints and strings of at
    most one character (interpreter singletons), the program's own constants
    (_shared_ids), and the objects of the types in |stop| (not entered either).
    Iterative: a deep trie does not recurse."""
    from array import array
    getsizeof = sys.getsizeof
    seen = set(_shared_ids())
    total = 0
    stack = list(roots)
    while stack:
        x = stack.pop()
        if x is None:
            continue
        i = id(x)
        if i in seen:
            continue
        seen.add(i)
        t = type(x)
        if t is bool:
            continue
        if t is int:
            if not -5 <= x <= 256:
                total += getsizeof(x)
            continue
        if t is str:
            if len(x) > 1:
                total += getsizeof(x)
            continue
        if stop is not None and t in stop:
            continue
        total += getsizeof(x)
        if t is float or t is array or t is bytes:
            continue
        if t is list or t is tuple or t is set or t is frozenset:
            stack.extend(x)
        elif t is dict:
            stack.extend(x.keys())
            stack.extend(x.values())
        elif hasattr(t, '__slots__'):
            for name in _slots_of(t):
                stack.append(getattr(x, name, None))
        elif hasattr(x, '__dict__'):
            stack.append(x.__dict__)
    return total


def memory_bytes(g):
    """The model's footprint in the heap: a deduplicating deep sizeof over the slotted
    model -- every record, label set, string, int and float it holds -- plus the caches
    the queries built (Graphlet.cache, Arm.cache), the envelope and the summary. What a
    GraphletStore charges against max_ram_mb (a shallow count would see 35-67 % of the
    traced heap)."""
    return _deep_bytes([g])


def cache_bytes(g):
    """The derived caches alone (Graphlet.cache and every Arm.cache), the model's records
    they refer to not counted: what the queries added since the model was measured. The
    store re-measures this, not the whole model, after a query (cost ~ the caches)."""
    roots = [g.cache] + [a.cache for a in g.arms.values()]
    return _deep_bytes(roots, stop=_model_types())


def cache_signature(g):
    """A cheap fingerprint of the caches: each cache's keys and the number of entries of
    each value, which changes whenever a query adds to them. The keys too, not only their
    count: a step may swap one key for another of the same length -- runs_by_label (one
    entry on a single-label arm) replacing label_runs()'s count leaves the lengths equal,
    and the store would never charge the index."""
    sig = [tuple(g.cache)]
    for a in g.arms.values():
        sig.append(tuple(a.cache))
        for v in a.cache.values():
            sig.append(len(v) if isinstance(v, (dict, list, tuple)) else 1)
    for v in g.cache.values():
        sig.append(len(v) if isinstance(v, (dict, list, tuple)) else 1)
    return tuple(sig)


def _kv(v):
    """A K value as JSON: the UNLIMITED sentinel reads 'unlimited' (as in the JSON
    summary); ints, floats and strings stay as they are."""
    return 'unlimited' if v is UNLIMITED else v


def informational(lim):
    """A stated limitation that limits nothing: the scope statement of an arm on which no
    history was united (observed 0 merges). The C++ outcome rule does not count it
    either (walks stays complete)."""
    return lim.kind == 'scope' and lim.observed == 0


def limitation_dict(lim, effect=True):
    """One K record as JSON: everything the agent needs to act on it -- the knob and its
    current limit, what was observed, how far the evidence holds, and the effect."""
    d = {'kind': lim.kind, 'knob': lim.knob, 'limit': _kv(lim.limit),
         'observed': _kv(lim.observed)}
    if lim.arm is not None:
        d['arm'] = lim.arm
    if lim.complete_to_bp is not None:
        d['complete_to_bp'] = lim.complete_to_bp
    for name, value in lim.extra:
        d.setdefault(name, _kv(value))
    if informational(lim):
        d['informational'] = True
    if effect:
        d['effect'] = lim.effect
    return d


def resource_stop_dict(q):
    if q is None:
        return None
    return {'scope': q.scope, 'resource': q.resource, 'phase': q.phase,
            'requested': _kv(q.requested), 'effective': _kv(q.effective),
            'used': _kv(q.used), 'remaining': _kv(q.remaining), 'actions': list(q.actions),
            'message': q.message}


def arm_exact(g, side):
    """Whether the label evidence of one arm is exact: no label-class limitation (a
    lower bound or an overstatement) stated for that arm or for the seed, and no recorded
    list of the arm cut. Per ARM: a limitation of the other arm does not make this one
    inexact. A graphlet whose outcome says label evidence is not complete while no K
    record says where (a body whose K records were lost) is exact on no arm."""
    if g.arms[side].labels_per_node.nodes_truncated:
        return False
    stated = False
    for lim in g.limitations:
        # a walked result's partial derivation overstates like trace_record_boundaries does
        if lim.kind in LOWER_BOUND_KINDS or derive.qualifies(g, lim):
            stated = True
            if lim.arm is None or lim.arm == side:
                return False
    return stated or g.outcome is None or g.outcome.label_evidence == 'complete'


def evidence_block(g, side=None, view=None):
    """The block every local answer carries (§6): what it is certified to and how. Beyond
    §6's {complete_to_bp, exact, support, reconverge, scope}: the four outcome dimensions
    and the kinds of the limitations that apply (seed level and per arm), so an answer
    read on its own never looks more certain than the retrieval. exact is PER ARM
    (arm_exact): it holds only when the label evidence of every arm the answer covers
    (|side|, or every arm) is complete and none of their recorded lists was cut; the
    outcome's label_evidence beside it is the seed's. |view|: a GraphletView the answer
    was restricted to (its completeness is qualified 'for the selected labels'). A
    retrieval that asked for record coordinates adds `coordinates`: {kind, complete} of
    its block, or {reason} when it has none."""
    if side is None:
        sides = [s for s in ARM_SIDES if s in g.arms]
    else:
        sides = [g.arm(side).side]
    kinds = {}
    for lim in g.limitations:
        if informational(lim) or (lim.arm is not None and lim.arm not in sides):
            continue
        got = kinds.setdefault('seed' if lim.arm is None else lim.arm, [])
        if lim.kind not in got:
            got.append(lim.kind)
    if g.resource_stop is not None:
        kinds.setdefault('seed', []).append('resource_stop')
    out = {
        'complete_to_bp': {s: g.arms[s].complete_to_bp for s in sides},
        'exact': all(arm_exact(g, s) for s in sides),
        'support': g.support, 'reconverge': g.reconverge,
        'scope': {s: g.arms[s].scope for s in sides},
        'outcome': g.outcome.as_dict() if g.outcome is not None else None,
        'limitations': kinds,
    }
    if _C.requested(g):
        # only when the request asked for coordinates: no other evidence block carries the
        # key. The block's kind and completeness (a cut list or a lower-bound run is
        # complete false), or the reason there is none
        raw = g.seed_summary['coordinates']
        out['coordinates'] = {'kind': raw.get('kind'), 'complete': raw.get('complete')} \
            if isinstance(raw, dict) else {'reason': g.seed_summary.get('coordinates_reason')}
    if view is not None:
        out['view'] = view.view_spec()
        out['qualified'] = 'for the selected labels'
    return out


def summary(g, arm=None, max_bytes=2048, *, budget=None):
    """Agent-facing, <= max_bytes of JSON: per-arm status / completeness / counts, the
    top labels by direct_bp, the stated limitations (caveats: kind, knob, limit,
    observed, complete_to_bp and effect; informational ones flagged) and the resource
    stop with its suggested actions. A seed whose permitted set was derived from part of
    it says so (seed.derivation: {partial, kmers_read, of}); a retrieval that asked
    for record coordinates states their kind and completeness (coordinates), or the reason
    it has none. Under the cap the effect texts go first, then top
    labels, then caveats (each step stated). budget= (local limits): a stop raises, with no
    partial."""
    return _B.run_scoped(budget, 'summary', ('narrow_arm',),
                         lambda b: _summary(g, arm, max_bytes, b))


def _summary(g, arm, max_bytes, b):
    g.require_envelope('summary()')
    sides = [g.arm(arm).side] if arm is not None else [s for s in ARM_SIDES if s in g.arms]
    if b is not None:
        _uses_g(b, g, 'label_summary')
        n = len(g.labels)
        for side in sides:
            derive.uses(b, g, g.arms[side], 'leaves', 'paths')
        text = 600 + sum(300 + len(l.effect) for l in g.limitations)
        if _C.requested(g) or derive.partial_derivation(g) is not None:
            # the coordinates entry and the derivation statement (at most ~200 bytes)
            text += 200
        # per arm the labels ranked: the list and a key tuple per label
        b.charge(len(sides) * (sort_work(n) + 2 * n),
                 len(sides) * (list_bytes(n) + n * (TUPLE + 3 * LIST_ITEM + 2 * INT))
                 + 8 * text)
        b.phase = 'fit'
        # each size() is one json.dumps of the answer, which only shrinks
        size_price = 32 + (text >> _B.TEXT_SHIFT)
    rows = derive.label_summary(g)
    out = {'seed': {'seed_id': (g.seed_summary.get('seed') or {}).get('seed_id', ''),
                    'length_bp': g.seed.length_bp, 'mode': g.mode,
                    'labels': len(g.labels)},
           'outcome': g.outcome.as_dict(),
           'index': {'ns': g.index_ns, 'verifiable': g.index_fp is not None},
           'arms': {}, 'caveats': []}
    pd = derive.partial_derivation(g)
    if pd is not None:
        # the partial derivation named (SPEC §7.0): the permitted set was derived from part
        # of the seed, a superset of its carriers, and the walk stopped at the seed
        out['seed']['derivation'] = {'partial': True, 'kmers_read': _kv(pd.observed),
                                     'of': g.seed.num_kmers}
    if _C.requested(g):
        raw = g.seed_summary['coordinates']
        if isinstance(raw, dict):
            got = {'kind': raw.get('kind'), 'complete': raw.get('complete'),
                   'max_occurrences': raw.get('max_occurrences')}
            if raw.get('runs_lower_bound'):
                got['runs_lower_bound'] = raw['runs_lower_bound']
            out['coordinates'] = got
        else:
            out['coordinates'] = {'reason': g.seed_summary.get('coordinates_reason')}
    for side in sides:
        a = g.arms[side]
        ranked = sorted(range(len(g.labels)),
                        key=lambda l: (-rows[l][side]['direct_bp'], -rows[l][side]['reach_bp'],
                                       l))
        top = [dict(g.labels[l].as_dict(), direct_bp=rows[l][side]['direct_bp'],
                    reach_bp=rows[l][side]['reach_bp']) for l in ranked[:5]]
        out['arms'][side] = {
            'status': a.status, 'complete_to_bp': a.complete_to_bp, 'scope': a.scope,
            'walks': a.counts.leaves, 'segments': a.counts.segments,
            'splits': a.counts.splits, 'merges': a.counts.merges, 'bases': a.counts.bases,
            'max_bp': max((p.length_bp for p in derive.paths(a)), default=0),
            'top_labels': top}
    for lim in g.limitations:
        if lim.arm is None or lim.arm in sides:
            d = limitation_dict(lim)
            d.setdefault('arm', None)
            out['caveats'].append(d)
    if g.resource_stop is not None:
        q = resource_stop_dict(g.resource_stop)
        out['caveats'].append({'kind': 'resource_stop', 'resource': q['resource'],
                               'phase': q['phase'], 'actions': q['actions'],
                               'message': q['message']})

    def size():
        if b is not None:
            b.charge(size_price, 2 * text)
            b.release(2 * text)
        return len(json.dumps(out, separators=(',', ':'), ensure_ascii=False).encode())

    # the effect texts are the longest and the most derivable part (they restate the
    # kind and the knob): they go first, and the cut is stated
    if size() > max_bytes:
        for c in out['caveats']:
            if c.pop('effect', None) is not None:
                out['caveat_effects_cut'] = True
            if c.get('kind') == 'resource_stop' and size() > max_bytes:
                c.pop('message', None)
    while size() > max_bytes:
        trimmed = False
        for side in sides:
            if out['arms'][side]['top_labels']:
                out['arms'][side]['top_labels'].pop()
                trimmed = True
        if not trimmed:
            if out['caveats']:
                out['caveats'].pop()
                out['caveats_truncated'] = True
            else:
                break
    return out
