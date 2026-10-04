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

Resources: everything here runs locally with NO work or allocation budget -- compare(),
routes() and the exports included (stage 2 of the design's resource contract bounds the
backend walk only). Their time and peak allocation follow the graphlet's size (a
comparison reads both DAGs up to the comparison depth, routes() walks one route per run,
an export spells every chosen walk), and nothing interrupts them; an MCP tool's max_bytes
bounds the bytes it RETURNS, not the computation behind them.
"""

import collections
import copy
import hashlib
import heapq
import json
import math
import sys

from . import derive
from ._codec import REASON, RESOURCE_CODES, UNLIMITED
from .model import (
    AmbiguousLabel, ARM_SIDES, BadSelector, Branch, Change, Claim, Comparison,
    Continuation, IncompatibleContinuations, IncompleteRecording, Label, LabelWalk,
    NextRequest, Path, SplitPoint, SupportRun, UnknownLabel, UnverifiableLabelName, Walk,
)

__all__ = [
    'label', 'labels_matching', 'path_id', 'leaf_segment', 'spell', 'walks', 'claims',
    'rank_walks', 'walks_at', 'walk_filter_reason',
    'label_walks', 'routes', 'support_profile', 'support_changes', 'label_summary',
    'splits', 'continuation', 'next_request', 'next_requests', 'resubmittable_names',
    'subgraph',
    'GraphletView', 'view_from_saved',
    'compare', 'index_identity', 'memory_bytes', 'cache_bytes', 'cache_signature',
    'summary', 'evidence_block',
    'END_CLASS', 'alive_at', 'limitation_dict', 'resource_stop_dict', 'informational',
]


# ------------------------------------------------------------------ labels

def label(g, selector):
    """Tagged selectors (v3): {'ref': 'h:<column>:<seq_id>' | 'c:<column>'} |
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
        hits = [l for l in g.labels if l.ref == selector or l.name == selector]
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
        if not isinstance(value, int) or not 0 <= value < len(g.labels):
            raise UnknownLabel('no label id %r' % (value,))
        return g.labels[value]
    if key == 'ref':
        hits = [l for l in g.labels if l.ref == value]
    elif key == 'name':
        hits = [l for l in g.labels if l.name == value]
    else:
        raise BadSelector('unknown selector tag %r' % key)
    if len(hits) == 1:
        return hits[0]
    if not hits:
        raise UnknownLabel('no label with %s %r' % (key, value))
    raise AmbiguousLabel('%d labels named %r (%s): select by ref' % (
        len(hits), value, ', '.join(l.ref for l in hits)))


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

def spell(g, arm, leaf, orientation='natural', with_seed=False):
    """The walk's bases. natural (§2.4): right = seed + flank, left = flank + seed
    (seed only with with_seed). walk: walking order (outward index i is [i]); with the
    seed, the seed is read in the walking direction too (reversed on the left arm),
    so position |seed| + i holds outward base i on both arms."""
    arm = g.arm(arm)
    seg = leaf_segment(arm, leaf) if not isinstance(leaf, Path) else leaf.leaf
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
    """Constrain: the label changes inside segment |s| as (position, 0 = a run end | 1 = a
    switch-in, label) in the chronological order of the G.end rule (ties: ends first). A
    run ending at to_bp is absent from base to_bp on; a switch-in is present from at_bp."""
    rbs = derive.runs_by_segment(arm)
    ops = []
    for r in rbs[s.id]:
        run = arm.runs[r]
        if run.to_bp < s.end_bp:
            ops.append((run.to_bp, 0, run.label))
    for ev in s.events:
        if ev.type == 'switch':
            ops.append((ev.at_bp, 1, ev.to_label))
    ops.sort(key=lambda o: (o[0], o[1]))
    return ops


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


def _alive_profile(g, arm, leaf_seg):
    """Constrain: the labels alive at each base of the leaf's first-parent chain, as
    maximal runs [(from, to, frozenset)]. Entry sets, then per segment the run ends
    anchored there (removed from to_bp on) and the switch-ins (added from at_bp on) in
    the chronological order of the G.end rule."""
    out = []
    for sid in derive.chain(arm, leaf_seg):
        out.extend(_alive_pieces(arm, arm.segments[sid]))
    return _merge_runs(out)


def _alive_pieces(arm, s):
    """Constrain: the alive sets inside one segment as [(from, to, frozenset)] (not
    merged; nothing for a zero-length segment)."""
    out = []
    if s.length_bp == 0:
        return out
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
        out.append((pos, nxt, frozenset(cur)))
        pos = nxt
    return out


def _presence_profile(arm, leaf_seg):
    out = []
    for sid in derive.chain(arm, leaf_seg):
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


def _run_claim(g, arm, run, cut, exact=True):
    """The claim of one run end, or None when the run-start guard drops it at |cut|.
    |exact|: the arm's label evidence is exact (arm_exact)."""
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
    every walk's chain was quadratic on a comb-shaped trie."""
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
    segs = arm.segments
    alive_end = [None] * len(segs)
    ends = {}
    sets = {}                                # equal sets shared (a comb repeats a few)
    for s in segs:
        if s.parents:
            alive = set()
            for p in s.parents:
                alive.update(alive_end[p])
        else:
            alive = set(s.entry)
        for pr in s.presence:
            if not alive:
                break
            still = alive.intersection(pr.labels)
            if pr.from_bp > s.from_bp or not s.parents:
                for l in alive - still:
                    ends.setdefault(l, []).append((s.id, pr.from_bp))
            alive = still
        fs = frozenset(alive)
        alive_end[s.id] = sets.setdefault(fs, fs)
    for s in segs:
        alive = alive_end[s.id]
        if not alive:
            continue
        carried = set()
        for c in s.children:
            if not segs[c].presence:
                carried = alive
                break
            carried.update(alive.intersection(segs[c].presence[0].labels))
        for l in sorted(alive - carried):
            ends.setdefault(l, []).append((s.id, s.end_bp))
    got = (alive_end, ends)
    arm.cache['annotate_route_ends'] = got
    return got


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


def claims(g, arm=None, labels=None, at_most_bp=None, strict=True):
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
    (_restricted_claims)."""
    sides = [g.arm(arm).side] if arm is not None else [s for s in ARM_SIDES if s in g.arms]
    wanted = _label_ids(g, labels)
    out = []
    for side in sides:
        a = g.arms[side]
        exact = arm_exact(g, side)
        if g.mode == 'constrain':
            for run in a.runs:
                if wanted is not None and run.label not in wanted:
                    continue
                c = _run_claim(g, a, run, at_most_bp, exact)
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


# ------------------------------------------------------------------ walks

def _labels_full(g, arm, p):
    """Labels with displayed support on the whole walk [0, length)."""
    if g.mode == 'constrain':
        ends = derive.end_labels(arm, p.leaf)
        return [g.labels[e.label] for e in ends
                if derive.evidence(arm, arm.runs[e.run])[1] == 0]
    return _labels(g, _displayed_presence(arm)[p.leaf][0])


def walks(g, arm, *, top=None, by='support', labels=None, route_consistent=True, min_bp=0):
    """Walk(path_id, leaf, segments, length_bp, sequence, path_reason, end_reasons,
    claims, labels_full, n_alive, complete, beyond_certified_bp). 'support' ranks by
    (-len(labels_full), -n_alive, -length_bp, path_id); 'length' by (-length_bp,
    path_id); 'loss' by (min end-label loss, -n_alive, -length_bp, path_id); 'id' keeps
    path order. Ranking uses the per-leaf records only; chains, spellings and claims
    are produced for the walks returned."""
    a = g.arm(arm)
    ranked = _ranked_walks(g, a, by, _label_ids(g, labels), route_consistent, min_bp)
    if top is not None:
        ranked = ranked[:top]
    return _walk_objects(g, a, ranked, with_claims=True)


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


def _ranked_walks(g, a, by, wanted, route_consistent, min_bp):
    """[(path, labels_full, alive ids, min loss)] in walks() order; nothing is spelled."""
    if by not in _WALK_KEYS:
        raise ValueError("by is 'support', 'length', 'loss' or 'id'")
    ranked = []
    for p in derive.paths(a):
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
    ranked.sort(key=_WALK_KEYS[by])
    return ranked


def rank_walks(g, arm, *, by='support', labels=None, route_consistent=True, min_bp=0):
    """The path ids walks() returns, in its order, from the per-leaf records only: no
    chain, spelling or claim is built (a pager ranks the whole arm, then materializes
    one page with walks_at())."""
    a = g.arm(arm)
    return [x[0].id for x in _ranked_walks(g, a, by, _label_ids(g, labels),
                                           route_consistent, min_bp)]


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


def walks_at(g, arm, path_ids, with_claims=False):
    """The Walk objects of |path_ids|, in that order (claims only with with_claims)."""
    a = g.arm(arm)
    ps = derive.paths(a)
    rows = []
    for pid in path_ids:
        p = ps[path_id(a, pid)]
        full, alive_ids, _, min_loss = _walk_keys_of(g, a, p)
        rows.append((p, full if full is not None else _labels_full(g, a, p), alive_ids,
                     min_loss))
    return _walk_objects(g, a, rows, with_claims)


def _walk_objects(g, a, ranked, with_claims):
    by_leaf = {}
    if ranked and with_claims:
        leaves_ = {x[0].leaf for x in ranked}
        for c in claims(g, a, strict=False):
            s = a.segments[c.segment]
            if c.segment in leaves_ and c.to_bp == s.end_bp:
                by_leaf.setdefault(c.segment, []).append(c)
    out = []
    for p, full, alive_ids, _ in ranked:
        seg = a.segments[p.leaf]
        reasons = {}
        if g.mode == 'constrain':
            for e in derive.end_labels(a, p.leaf):
                name = a.runs[e.run].reason
                reasons[name] = reasons.get(name, 0) + 1
        code = seg.leaf.path_reason
        try:
            sequence = derive.natural_flank(a, p.leaf)
        except ValueError:
            sequence = None
        out.append(Walk(
            path_id=p.id, leaf=p.leaf, segments=p.segments, length_bp=p.length_bp,
            sequence=sequence, path_reason=REASON[code] if code else None,
            end_reasons=reasons, claims=by_leaf.get(p.leaf, []), labels_full=full,
            n_alive=len(alive_ids),
            complete=(code is None or code not in RESOURCE_CODES)
            and p.length_bp <= a.complete_to_bp,
            beyond_certified_bp=max(0, p.length_bp - a.complete_to_bp)))
    return out


# ------------------------------------------------------------------ per label

def routes(g, sel, arm, spell=False):
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

    No work or allocation budget applies (one route per run or end, each as long as the
    DAG makes it): see the module docstring."""
    a = g.arm(arm)
    lab = label(g, sel)
    out = []
    if g.mode != 'constrain':
        for route, anchor, j in _annotate_routes(a, lab.id):
            out.append((route, _route_bases(a, route, j)) if spell else route)
        return out
    segs = a.segments
    for run in a.runs:
        if run.label != lab.id:
            continue
        route = [run.segment]
        s = run.segment
        while segs[s].parents:
            seg = segs[s]
            nxt = seg.parents[0]
            if len(seg.parents) > 1:
                at = derive.lineage_label_at(a, run, seg.from_bp)
                for j, part in enumerate(seg.partition):
                    if at in part:
                        nxt = seg.parents[j]
                        break
            route.append(nxt)
            s = nxt
        route.reverse()
        if spell:
            bases = ''.join(segs[x].walk or '' for x in route)[:run.to_bp]
            if segs[route[0]].walk is None:
                bases = None
            elif a.side == 'left':
                bases = bases[::-1]
            out.append((route, bases))
        else:
            out.append(route)
    return out


def label_walks(g, sel, arm=None):
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
    otherwise -- §5.1, as for a run in constrain mode)."""
    lab = label(g, sel)
    sides = [g.arm(arm).side] if arm is not None else [s for s in ARM_SIDES if s in g.arms]
    out = []
    for side in sides:
        a = g.arms[side]
        if g.mode != 'constrain':
            for route, anchor, j in _annotate_routes(a, lab.id):
                seg = a.segments[anchor]
                if seg.leaf is not None and j == seg.end_bp:
                    code = seg.leaf.path_reason
                    end = REASON[code] if code else None
                else:
                    end = 'label_lost'
                out.append(LabelWalk(None, lab, side, 0, j,
                                     _annotate_evidence(a, lab.id, anchor),
                                     _route_bases(a, route, j), end,
                                     _displayed_leaves_below(a, anchor), None))
            continue
        rts = routes(g, lab, a, spell=True)
        k = 0
        for run in a.runs:
            if run.label != lab.id:
                continue
            route, bases = rts[k]
            k += 1
            if bases is not None:
                bases = bases[run.from_bp:run.to_bp] if side == 'right' else \
                    bases[len(bases) - run.to_bp:len(bases) - run.from_bp]
            merged_into = None
            if run.merged:
                for c in a.segments[run.segment].children:
                    if len(a.segments[c].parents) > 1:
                        merged_into = c
            end = 'merged' if run.merged else 'switched' if run.silent else run.event_reason
            out.append(LabelWalk(run.id, lab, side, run.from_bp, run.to_bp,
                                 derive.evidence(a, run)[1], bases, end,
                                 derive.run_leaves(a, run), merged_into))
    return out


def support_profile(g, arm, leaf, kind='displayed'):
    """Maximal SupportRun(from_bp, to_bp, labels, exact) over the walk [0, length).
    displayed: the labels supporting the displayed bases (§5.1); route: also every run
    anchored on the walk over its own [from_bp, to_bp) (route support, scope-free).
    Annotate mode records one set per node: both kinds are the recorded sets."""
    a = g.arm(arm)
    seg_id = leaf_segment(a, leaf)
    if g.mode != 'constrain':
        return [SupportRun(f, t, _labels(g, s), tot == len(s), tot)
                for f, t, s, tot in _presence_profile(a, seg_id)]
    prof = _alive_profile(g, a, seg_id)
    if kind == 'route':
        extra = []
        rbs = derive.runs_by_segment(a)
        for sid in derive.chain(a, seg_id):
            for r in rbs[sid]:
                run = a.runs[r]
                if run.to_bp > run.from_bp:
                    extra.append(run)
        if extra:
            cuts = sorted({0, a.segments[seg_id].end_bp}
                          | {x for f, t, _ in prof for x in (f, t)}
                          | {x for r in extra for x in (r.from_bp, r.to_bp)})
            pieces = []
            for f, t in zip(cuts, cuts[1:]):
                base = next((s for pf, pt, s in prof if pf <= f < pt), frozenset())
                add = {r.label for r in extra if r.from_bp <= f < r.to_bp}
                pieces.append((f, t, frozenset(base | add)))
            prof = _merge_runs(pieces)
    elif kind != 'displayed':
        raise ValueError("kind is 'displayed' or 'route'")
    return [SupportRun(f, t, _labels(g, s), True, len(s)) for f, t, s in prof]


def support_changes(g, arm, leaf):
    """Where the displayed support changes along a walk, and why: an end (the R code,
    qualifier text included), a switch (Lw away / a switch-in), a split (labels that
    took another branch), a merge (a lineage joining through another parent)."""
    a = g.arm(arm)
    seg_id = leaf_segment(a, leaf)
    chain = derive.chain(a, seg_id)
    segs = a.segments
    prof = support_profile(g, a, leaf)
    prev = set(segs[chain[0]].entry)
    out = []
    on_chain = set(chain)
    rbs = derive.runs_by_segment(a)
    for run in prof:
        cur = {l.id for l in run.labels}
        at = run.from_bp
        added, removed = cur - prev, prev - cur
        if added or removed:
            out.append(Change(at, _labels(g, added), _labels(g, removed),
                              _change_reasons(g, a, chain, on_chain, rbs, at, added, removed)))
        prev = cur
    # labels alive at the leaf's last node stay; nothing to report past the end
    return out


def _change_reasons(g, a, chain, on_chain, rbs, at, added, removed):
    segs = a.segments
    reasons = []
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
            for sid in chain:
                s = segs[sid]
                if s.from_bp == at and len(s.parents) == 1:
                    parent = segs[s.parents[0]]
                    took = sorted({x for c in parent.children if c != sid
                                   for x in (segs[c].walk[:1] if segs[c].walk
                                             else segs[c].first_base or '')
                                   if l in segs[c].entry})
                    if l in parent.end and took:
                        # a split only when another branch took the label: on the
                        # parent's last node and on no child's first node it is
                        # absent (annotate) or ended, not split
                        why = {'why': 'split', 'took': took}
                    break
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
                for j, part in enumerate(s.partition):
                    if l in part and j > 0:
                        why = {'why': 'merge', 'via_parent': s.parents[j]}
        if why is None:
            why = {'why': 'present' if g.mode != 'constrain' else 'entered'}
        reasons.append(dict(why, label=g.labels[l].as_dict(), change='added'))
    return reasons


def splits(g, arm, min_labels_before=0):
    """The split records of an arm, in the walker's order (at_bp, then first child), as
    SplitPoint(arm, at_bp, segment, kind, labels_before, branches): oriented like the
    other accessors -- at_bp is the outward distance from the seed boundary (walking
    order on both arms, as a claim's interval), labels are Label objects ({name, ref}),
    and each Branch carries its first base, the labels at its first node (the recorded
    list, cut at labels.max_labels_per_node, with the true count) and the first walk
    whose displayed chain takes it. A branch's label count is not a share of
    labels_before: a label may follow several branches (kind 'ambiguous')."""
    a = g.arm(arm)
    first = _first_paths(a)
    out = []
    for sp in derive.splits(a, g.mode):
        if sp.labels_before < min_labels_before:
            continue
        branches = []
        for c in sp.children:
            child = a.segments[c]
            branches.append(Branch(c, derive.first_base(child),
                                   [g.labels[l] for l in child.entry[:g.cap]],
                                   child.entry_total, first[c]))
        out.append(SplitPoint(a.side, sp.at_bp, sp.segment,
                              'ambiguous' if sp.ambiguous else 'divergence',
                              sp.labels_before, branches))
    return out


def label_summary(g):
    """{ref: {name, ref, <arm>: {direct_bp, reach_bp, reentries, runs}}} (both modes,
    §5.1). direct_bp is label-consistent ROUTE support, not a contiguous occurrence."""
    rows = derive.label_summary(g)
    out = {}
    for lab, per in zip(g.labels, rows):
        row = lab.as_dict()
        for side, v in per.items():
            row[side] = dict(v, runs=list(v['runs']))
        out[lab.ref] = row
    return out


# ------------------------------------------------------------------ continuation

def continuation(g, arm, leaf):
    """Continuation(sequence, labels, loss_used, branches_used, seed_coord, as_seed(),
    losses, note): the walker's label logic from C, the spelling derived (both arms).
    seed_coord is the half-open whole-molecule interval in seed coordinates (seed =
    [0, |seed|)). loss_used is C's: the SMALLEST terminal loss of the labels; losses
    holds each label's own (T), and note says when one loss budget cannot serve them
    exactly (labels that ended the walk at different losses)."""
    a = g.arm(arm)
    seg_id = leaf_segment(a, leaf)
    seg = a.segments[seg_id]
    c = seg.leaf.continuation if seg.leaf else None
    if c is None:
        raise ValueError('walk %d of the %s arm has no continuation (it ended for a '
                         'semantic reason, not at the radius or a cap)'
                         % (path_id(a, leaf), a.side))
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
    is the SMALLEST of them: subtracting it from the budget let a continued label at a
    higher loss go on with more than its original budget left (the review's
    budgetchain: E at 0 and C at 2 kept the whole budget 3). () in annotate mode."""
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


def _switch_cost(change_cost, src, dst):
    """cost(src -> dst) by label NAME under the request's change_cost, as the server
    prices a switch (forbid: +inf; constant: value; table: the last entry for the pair,
    else its default, 'forbid' = +inf; a label to itself costs 0). None for a model the
    library does not know: reachability cannot be decided then. A malformed field of a
    model it knows is a ValueError naming the field (a caller's override; the tool layer
    answers bad_argument), never an exception of another kind."""
    if src == dst:
        return 0.0
    model = (change_cost or {}).get('model', 'forbid')
    if model == 'forbid':
        return math.inf
    if model == 'constant':
        v = change_cost.get('value')
        if not _is_number(v):
            raise ValueError('labels.change_cost.value is a number, not %r' % (v,))
        return float(v)
    if model == 'table':
        entries = change_cost.get('entries') or []
        if not isinstance(entries, list):
            raise ValueError('labels.change_cost.entries is a list of [from, to, cost], not %r'
                             % (entries,))
        got = None
        for e in entries:
            if isinstance(e, list) and len(e) == 3 and e[0] == src and e[1] == dst:
                if not _is_number(e[2]):
                    raise ValueError('labels.change_cost.entries: the cost of %r is a number, '
                                     'not %r' % (e[:2], e[2]))
                got = float(e[2])
        if got is not None:
            return got
        d = change_cost.get('default', 'forbid')
        if d != 'forbid' and not _is_number(d):
            raise ValueError('labels.change_cost.default is a number or "forbid", not %r' % (d,))
        return math.inf if d == 'forbid' else float(d)
    return None


def _section(strategy, *path):
    """The strategy section at |path| as a dict ({} when absent). The rebuild of a
    continuation reads the merged sections as objects; one that a caller's override made
    something else is refused with a ValueError naming it (the tool layer answers
    bad_argument), as the server would refuse it (round 3, finding B: an AttributeError
    escaped the tool contract)."""
    d = strategy
    for i, k in enumerate(path):
        v = d.get(k)
        if v is None:
            return {}
        if not isinstance(v, dict):
            raise ValueError('strategy.%s is an object, not %r' % ('.'.join(path[:i + 1]), v))
        d = v
    return d


def _switch_reach(change_cost, sources, targets, budget, sinks=()):
    """{name: loss} for every name of |targets| that a chain of switches from a name of
    |sources| enters within |budget|, at the cheapest such chain's loss (summed left to
    right, as the walk sums a lineage's loss) -- the server's rule for an extra label
    (walker.cpp, switch_reach): the walk enforces the cumulative loss switch by switch, so
    a label is a valid switch target when SOME chain reaches it, not only one switch from a
    seed label (the stage-2 recheck's design answer: A -> B = 1, B -> C = 1 under a budget
    of 2 left C out). A chain passes only through names of |sources| and |targets|: the
    labels the request names; a name of |sinks| (also a target) ends a chain but never
    continues one. Prices as _switch_cost does (a table's last entry for a pair
    wins, else its default; a label to itself costs 0); None for a model the library does
    not know. A malformed field of a model it knows is a ValueError naming it."""
    sources = list(dict.fromkeys(sources))
    targets = [t for t in dict.fromkeys(targets) if t not in set(sources)]
    model = (change_cost or {}).get('model', 'forbid')
    if not sources or model == 'forbid':
        return {}
    if model == 'constant':
        # one switch reaches every label, and a chain only costs more
        c = change_cost.get('value')
        if not _is_number(c):
            raise ValueError('labels.change_cost.value is a number, not %r' % (c,))
        return {t: float(c) for t in targets if float(c) <= budget}
    if model != 'table':
        return None
    names = set(sources) | set(targets)
    entries = change_cost.get('entries') or []
    if not isinstance(entries, list):
        raise ValueError('labels.change_cost.entries is a list of [from, to, cost], not %r'
                         % (entries,))
    table = {}
    for e in entries:
        if isinstance(e, list) and len(e) == 3 and e[0] in names and e[1] in names:
            if not _is_number(e[2]):
                raise ValueError('labels.change_cost.entries: the cost of %r is a number, '
                                 'not %r' % (e[:2], e[2]))
            table[(e[0], e[1])] = float(e[2])     # the last entry for a pair wins
    d = change_cost.get('default', 'forbid')
    if d != 'forbid' and not _is_number(d):
        raise ValueError('labels.change_cost.default is a number or "forbid", not %r' % (d,))
    fallback = math.inf if d == 'forbid' else float(d)
    out_edges = collections.defaultdict(list)
    for (u, v), c in table.items():
        if u != v and c != math.inf:
            out_edges[u].append((v, c))
    # the default is relaxed lazily, as on the server: the names no popped name has given
    # the default yet; a popped name gives it to each of them it has no explicit entry to
    pending = set(names) if fallback <= budget else set()
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
        loss, _, u = heapq.heappop(heap)
        if u in done or loss > dist[u]:
            continue
        done.add(u)
        if u in sinks:
            continue
        for v, c in out_edges.get(u, ()):
            relax(v, loss + c)
        for v in list(pending):
            if v != u and (u, v) not in table:
                relax(v, loss + fallback)
                pending.discard(v)
    return {t: dist[t] for t in targets if dist.get(t, math.inf) <= budget}


def _rebuilt_extra(g, seed_labels, budget, strategy):
    """labels.extra around a continuation's seed labels: the retrieval's permitted pool
    (its seed labels and its extra labels: every label of a constrain retrieval) minus the
    new seed labels -- the server refuses an extra label that duplicates a seed label
    (the review: seed ['C'] with extra ['B', 'C', 'D'] was a 400) -- keeping each
    remaining label that the server accepts as a switch target: some chain of switches
    from a new seed label enters it within |budget| (_switch_reach; an unreachable extra
    label is refused, and under forbid, or a constant above the budget, any extra label
    is). -> (extra names in dictionary order, the labels left out as unreachable, the
    labels left out because their names cannot be verified). The retrieval's original seed
    labels become switch targets too: they were in the original pool."""
    lab = strategy.get('labels') or {}
    cost = lab.get('change_cost') or {'model': 'forbid'}
    if cost.get('model') == 'constant' and not _is_number(cost.get('value')):
        raise ValueError('labels.change_cost.value is a number, not %r' % (cost.get('value'),))
    seed_names = [l.name for l in seed_labels]
    seed_ids = {l.id for l in seed_labels}
    seeds = set(seed_names)
    # /traverse resolves an extra label by NAME too: one whose name the library cannot verify
    # to be its own is never sent (left out and stated if a chain reaches it, unreachable
    # otherwise), so no chain is counted through it -- the request does not name it
    suspect = {l.id for l in g.labels if l.id not in seed_ids and l.name not in seeds
               and _unverifiable_names(g, [l])}
    others = [l for l in g.labels if l.id not in seed_ids and l.name not in seeds]
    reach = _switch_reach(cost, seed_names, [l.name for l in others], budget,
                          sinks={l.name for l in others if l.id in suspect})
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


def _continuation_seeds(g, a, leaves):
    conts = [continuation(g, a, x) for x in leaves]
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


def next_request(g, arm, leaves, bp=None, reduce_budget=True, reset_branches=False,
                 **overrides):
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
        branching), or a malformed field of it, is a ValueError naming it.
    A label alive at a leaf that does not cover the continuation's whole tail is not among
    the continuation's labels, so it is not seeded and the continuation may lack its
    lineage: stated in .notes and in .left_out (why 'alive_not_seeded').
    Several walks share one strategy, so one labels.extra: they are built into one request
    only when the rebuilt pool is the same for every walk (their continuations carry the
    same labels, or no switch is possible); otherwise IncompatibleContinuations names
    next_requests(), one request per walk. Keyword overrides deep-merge into the
    strategy; release, graph and graph_path are request-level. Raises MissingEnvelope on
    a body-only graphlet."""
    if not isinstance(reset_branches, bool):
        raise ValueError('reset_branches is true or false, not %r' % (reset_branches,))
    g.require_envelope('next_request()')
    a = g.arm(arm)
    if isinstance(leaves, (int, Path, Walk)):
        leaves = [leaves]
    conts, seeds = _continuation_seeds(g, a, leaves)
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
        final = copy.deepcopy(strategy)
        _deep_merge(final, {k: v for k, v in overrides.items()
                            if k not in ('release', 'graph', 'graph_path')})
        flab = _section(final, 'labels')
        _section(final, 'labels', 'change_cost')
        _section(final, 'branching')
        given = overrides.get('labels') or {}
        if 'loss_budget' in given and not _is_number(given['loss_budget']):
            raise ValueError('labels.loss_budget is a number, not %r' % (given['loss_budget'],))
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
            pools = [_rebuilt_extra(g, c.labels, effective, final) for c in conts]
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
            # label (the stage-2 recheck, P3: the first walk's omissions were reported for all,
            # naming a label another walk seeded as unreachable)
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
        # budget is accepted -- but never silent (round 3, finding D).
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
            # (the review of the stage-3 fixes, P3: a chain through a kept extra label was
            # not searched, and the note understated what the continuation may miss)
            mine = _unique_pairs([(l, x) for l, x in zip(c.labels, c.losses) if x != math.inf])
            pool = [l.name for l in g.labels]
            reach_of = [(p, _switch_reach(cost, [p.name], pool, budget - x) or {})
                        for p, x in mine]
            lost = []
            for l in dropped:
                via = [p for p, r in reach_of if l.name in r]
                lost.append('%s%s' % (_shown(l.name), ' (reachable for %s within its own '
                                      'remaining budget)' % ', '.join(_shown(p.name)
                                                                      for p in via)
                                      if via else ''))
            notes.append(
                'labels.extra leaves out %s%s: no chain of switches from a label of the '
                'continuation enters %s within loss_budget %s, and the server refuses an '
                'unreachable extra label%s'
                % (', '.join(lost[:8]) + (' and %d more' % (len(lost) - 8)
                                          if len(lost) > 8 else ''),
                   ' for walk %d' % c.leaf if len(per_walk) > 1 else '',
                   'it' if len(lost) == 1 else 'them', _num(effective),
                   '; one uninterrupted walk could still have entered %s from a label at a '
                   'lower loss' % ('it' if len(lost) == 1 else 'them')
                   if any('reachable for' in x for x in lost) else ''))
        # Each continued label's own branches at the leaf (T), not C's smallest. The server
        # resets branch state for a new seed, so the allowance is reduced as the loss budget
        # is: by the largest count, conservative for the labels that used fewer (the stage-2
        # recheck's design answer 3); reset_branches keeps it, and the note says so
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
        if 'max_label_branches' in given_b and not (
                mlb == 'unlimited' or (_is_number(mlb) and int(mlb) == mlb and mlb >= 0)):
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
            # budget would be 0 and say nothing, and a long list of them cost a receipt its
            # room (round 3, finding C)
            budget_info = None
    request = NextRequest({'seeds': seeds}, notes, budget_info, left_out, branch_info)
    release = g.envelope.get('release')
    if release:
        request['release'] = release
    for k in ('release', 'graph', 'graph_path'):
        if k in overrides:
            request[k] = overrides.pop(k)
    _deep_merge(strategy, overrides)
    request['strategy'] = strategy
    return request


def next_requests(g, arm, leaves, bp=None, reduce_budget=True, reset_branches=False,
                  **overrides):
    """One next_request() per walk, in |leaves| order: each with its own seed labels, its
    own rebuilt labels.extra and its budgets reduced by ITS labels' largest terminal loss
    and branch count (never less conservative than one request over all the walks)."""
    if isinstance(leaves, (int, Path, Walk)):
        leaves = [leaves]
    return [next_request(g, arm, [x], bp, reduce_budget, reset_branches,
                         **copy.deepcopy(overrides))
            for x in leaves]


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
    that the index may hold another label of exactly that name (a server now refuses
    such a seed instead). The library has only the retrieval to go by."""
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


def resubmittable_names(g, labels):
    """The names of |labels| (selectors) as a request would carry them; raises
    UnverifiableLabelName when one may resolve to another label (see continuation)."""
    labs = labels_matching(g, labels) or []
    why = _unverifiable_names(g, labs)
    if why:
        raise UnverifiableLabelName(why)
    return [l.name for l in labs]


def _deep_merge(dst, src):
    for k, v in src.items():
        if isinstance(v, dict) and isinstance(dst.get(k), dict):
            _deep_merge(dst[k], v)
        else:
            dst[k] = copy.deepcopy(v)


# ------------------------------------------------------------------ views

class GraphletView:
    """A view over a backing graphlet (§5, v3): every id is the backing graphlet's
    original id (no remapping, no new record type). save() writes the UNCHANGED
    backing body plus, in J, {view: {selectors, mode, arm, of}}; load() restores the
    view. Its completeness is the backing complete_to_bp, qualified 'for the selected
    labels'."""

    __slots__ = ('backing', 'selectors', 'mode', 'arm_side', 'labels', 'segments',
                 'path_ids', 'of')

    def __init__(self, backing, selectors, arm=None, mode='any', of=None):
        if mode not in ('any', 'all'):
            raise ValueError("mode is 'any' or 'all'")
        self.backing = backing
        self.selectors = _jsonable_selectors(backing, selectors)
        self.mode = mode
        self.arm_side = None if arm is None else backing.arm(arm).side
        self.labels = labels_matching(backing, selectors)
        self.of = of or ('sha256:' + hashlib.sha256(
            backing.dump(envelope=False).encode('utf-8')).hexdigest())
        ids = {l.id for l in self.labels}
        self.segments = {}
        self.path_ids = {}
        for side, a in backing.arms.items():
            if self.arm_side is not None and side != self.arm_side:
                continue
            keep = set()
            rbs = derive.runs_by_segment(a)
            for s in a.segments:
                present = set(s.entry) | set(s.end) | {a.runs[r].label for r in rbs[s.id]}
                for p in s.presence:
                    present.update(p.labels)
                hit = ids & present
                if (mode == 'any' and hit) or (mode == 'all' and hit == ids):
                    keep.update(derive.chain(a, s.id))
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

    def walks(self, arm, **kw):
        a = self.backing.arm(arm)
        keep = set(self.path_ids.get(a.side, ()))
        return [w for w in self.backing.walks(a, **kw) if w.path_id in keep]

    def claims(self, arm=None, at_most_bp=None, strict=True):
        if arm is None and self.arm_side is not None:
            arm = self.arm_side
        return self.backing.claims(arm, self.labels, at_most_bp, strict)

    def to_fasta(self, arm=None, with_seed=True, orientation='natural', width=None):
        from . import export
        out = []
        for side in self.segments:
            if arm is not None and self.backing.arm(arm).side != side:
                continue
            out.append(export.to_fasta(self.backing, side, self.path_ids[side], with_seed,
                                       orientation, width))
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

    def dump(self, envelope=True):
        """The backing body, unchanged; the view lives in J only (envelope=False: the
        body alone, which no longer says it is a view)."""
        return self._with_view().dump(envelope=envelope)

    def save(self, path):
        from . import parser
        return parser.save(self._with_view(), path)


def _jsonable_selectors(g, selectors):
    out = []
    for lab in labels_matching(g, selectors):
        out.append({'ref': lab.ref})
    return out


def subgraph(g, selectors, arm=None, mode='any'):
    return GraphletView(g, selectors, arm, mode)


def view_from_saved(g):
    spec = g.view
    backing = g
    backing.view = None
    sels = spec.get('selectors') or []
    return GraphletView(backing, sels, spec.get('arm'), spec.get('mode', 'any'),
                        spec.get('of'))


# ------------------------------------------------------------------ comparison

def index_identity(g):
    """(§3.1) {ns, fp, meta_fp, release}: fp is the index bundle's manifest digest
    (None: no manifest, joins are unverifiable), meta_fp a NEGATIVE check only."""
    release = (g.envelope or {}).get('release') or None
    return {'ns': g.index_ns, 'fp': g.index_fp, 'meta_fp': g.index_meta_fp,
            'release': release}


def _comparability(a, b):
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


# how weak a comparability verdict is: a comparison only ever moves to a weaker one (an
# unverifiable index identity is not made 'qualified' by a scope note, and so on)
_COMPARABLE_RANK = {True: 0, 'qualified': 1, 'unverifiable': 2, 'unknown': 3, False: 4}

# the label-class limitations (spec §7.0, design §14): evidence may be missing or
# understated (lower bound) or overstated (qualified)
LOWER_BOUND_KINDS = ('label_lists', 'inexact_counts', 'seed_labels', 'switch_sources',
                     'greedy_losses')
OVERSTATED_KINDS = ('trace_record_boundaries',)


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
        if lim.kind in LOWER_BOUND_KINDS or lim.kind in OVERSTATED_KINDS:
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


def compare(a, b, *, arm=None, labels=None, mode='claims'):
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
    claims instead let merges beyond the depth decide what lies before it: a branch that
    enters such a merge through a non-first parent is shown by no walk of the deeper DAG,
    and its lineage was route_only on the other branch's bases, so a retrieval and a
    deeper one of the same seed differed (the review's bubbles, radius 30 vs 100 at depth
    30, merges at 36 and 62). A merge AT the depth lies outside the restricted DAG too: a
    retrieval walked to a merge position records it as a zero-length segment, and the runs
    it closes are open claims at the depth, not merged ones (round 3, finding A). Where a
    route cannot be reconstructed through a merge (its lineage in no partition) the
    comparison is 'qualified', never a difference.

    No work or allocation budget applies: the cost follows both DAGs up to the depth (see
    the module docstring)."""
    comparable, reason, notes = _comparability(a, b)
    strategies = ((a.envelope or {}).get('strategy'), (b.envelope or {}).get('strategy'))
    if comparable is False:
        return Comparison(False, reason, None, None, [], [], notes, mode=mode,
                          support=(a.support, b.support), strategies=strategies)
    sides = [a.arm(arm).side] if arm is not None else \
        [s for s in ARM_SIDES if s in a.arms and s in b.arms]
    if not sides:
        return Comparison(False, 'no arm in common', None, None, [], [], notes, mode=mode)
    scopes = (a.arms[sides[0]].scope, b.arms[sides[0]].scope)
    depth = min(min(a.arms[s].complete_to_bp, b.arms[s].complete_to_bp) for s in sides)
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
    if mode not in ('claims', 'walks', 'labels', 'prefix_subset'):
        raise ValueError("mode is 'claims', 'walks', 'labels' or 'prefix_subset'")
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
        ka = _claim_keys(a, sides, depth, sel, inexact)
        kb = _claim_keys(b, sides, depth, sel, inexact)
        only_a = [_key_json(k, v, lab) for k, v in sorted(ka.items()) if k not in kb]
        only_b = [_key_json(k, v, lab) for k, v in sorted(kb.items()) if k not in ka]
        differ = [{'key': _key_json(k, None, lab), 'a': ka[k], 'b': kb[k]}
                  for k in sorted(ka) if k in kb and ka[k] != kb[k]]
    elif mode == 'walks':
        wa = _walk_keys(a, sides, depth, sel, memo)
        wb = _walk_keys(b, sides, depth, sel, memo)
        only_a = [{'arm': k[0], 'walk': k[1], 'labels': [lab(r) for r in sorted(v)]}
                  for k, v in sorted(wa.items()) if k not in wb]
        only_b = [{'arm': k[0], 'walk': k[1], 'labels': [lab(r) for r in sorted(v)]}
                  for k, v in sorted(wb.items()) if k not in wa]
        differ = [{'arm': k[0], 'walk': k[1], 'a': [lab(r) for r in sorted(wa[k])],
                   'b': [lab(r) for r in sorted(wb[k])]}
                  for k in sorted(wa) if k in wb and wa[k] != wb[k]]
    elif mode == 'labels':
        la, lb = _label_keys(a, sides, depth, sel), _label_keys(b, sides, depth, sel)
        only_a = [dict(lab(k[1]), arm=k[0]) for k in sorted(la) if k not in lb]
        only_b = [dict(lab(k[1]), arm=k[0]) for k in sorted(lb) if k not in la]
        differ = [dict(lab(k[1]), arm=k[0], a=la[k], b=lb[k])
                  for k in sorted(la) if k in lb and la[k] != lb[k]]
    else:
        only_a, only_b = _prefix_subset(a, b, sides, depth, sel, memo, inexact)
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


def _claim_keys(g, sides, depth, sel, inexact=None):
    """The claims of each arm's DAG RESTRICTED to [0, depth) (_restricted_claims) as
    comparison keys -> {(side, ref, evidence_from or -1, displayed prefix, to_bp):
    end class}."""
    out = {}
    for side in sides:
        a = g.arms[side]
        for c in _restricted_claims(g, side, depth, inexact):
            if c.kind == 'merged':
                continue             # not an end: the lineage goes on in the kept entry
            if sel is not None and c.label.ref not in sel:
                continue
            try:
                prefix = _chain_bases(a, c.segment, 0, c.to_bp)
            except ValueError:
                prefix = None
            key = (side, c.label.ref, c.evidence_from if c.evidence_from is not None else -1,
                   prefix, c.to_bp)
            out[key] = 'open' if _open_at(c, depth) else c.end_class
    return out


def _segment_support(g, arm, s):
    """The displayed support inside one segment as [(from, to, frozenset of ids)]:
    the alive sets (constrain) or the recorded P sets (annotate)."""
    if g.mode == 'constrain':
        return _alive_pieces(arm, s)
    return [(p.from_bp, p.to_bp, frozenset(p.labels)) for p in s.presence]


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
    inter = None
    for sid in derive.chain(arm, anchor):
        s = segs[sid]
        if s.walk is None:
            parts = None
        elif parts is not None:
            parts.append(s.walk if s.end_bp <= m else s.walk[:m - s.from_bp])
        for f, t, labs in _segment_support(g, arm, s):
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
    depth it is a walk of its own (the review's bubbles, merged at 36 and 62, compared at
    30: comparing the deeper DAG's displayed walks cut at 30 reported differences that
    the merges beyond 30 caused)."""
    return [s.id for s in arm.segments
            if s.from_bp < depth and (s.leaf is not None or s.end_bp >= depth)]


def _cuts(g, side, depth, memo):
    """[(leaf, m, cut info)] for every leaf of the arm's DAG restricted to [0, depth)
    (restricted_leaves), m = min(its end, depth). Leaves sharing the segment that holds
    base m - 1 share one cut: the work is the distinct anchors' chains up to the depth,
    not every walk's whole chain (B1). Memoized per comparison in |memo|."""
    key = (id(g), side, depth)
    got = memo.get(key)
    if got is not None:
        return got
    a = g.arms[side]
    segs = a.segments
    infos = {}
    out = []
    if depth <= 0:
        # nothing is certified: one empty cut with the boundary labels, when there is a walk
        if derive.paths(a):
            out.append((None, 0, _cut_info(g, a, None, 0)))
        memo[key] = out
        return out
    for leaf in restricted_leaves(a, depth):
        m = min(segs[leaf].end_bp, depth)
        if m == 0:
            k = None
        else:
            k = leaf                          # holds base m - 1 when it reaches the depth
            while segs[k].from_bp >= m:      # zero-length segments at the leaf's end
                k = segs[k].parents[0]
        info = infos.get((k, m))
        if info is None:
            info = infos[(k, m)] = _cut_info(g, a, k, m)
        out.append((leaf, m, info))
    memo[key] = out
    return out


def _route_parent(a, run, seg):
    """The parent of merge |seg| through which |run|'s lineage arrived (its incoming label
    in that parent's partition) -> (parent, exact). Not found: parents[0], not exact."""
    at = derive.lineage_label_at(a, run, seg.from_bp)
    for j, part in enumerate(seg.partition):
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
                not in seg.partition[0]:
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
    and a run anchored on it is anchored on its own route (round 3, finding A: a radius
    equal to a merge position compared unequal with a deeper retrieval). Any other run
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
    # a merge AT the depth (zero length, ending there) is outside it (round 3, finding A)
    if all(s.end_bp <= depth and (s.from_bp < depth or not s.parents) for s in segs):
        return _annotate_route_ends(arm)[1]
    alive_end = [None] * len(segs)
    ends = {}
    for s in segs:
        if s.from_bp >= depth:
            continue
        if s.parents:
            alive = set()
            for p in s.parents:
                alive.update(alive_end[p])
        else:
            alive = set(s.entry)
        for pr in s.presence:
            if not alive or pr.from_bp >= depth:
                break
            still = alive.intersection(pr.labels)
            if pr.from_bp > s.from_bp or not s.parents:
                for l in alive - still:
                    ends.setdefault(l, []).append((s.id, pr.from_bp))
            alive = still
        alive_end[s.id] = frozenset(alive)
    for s in segs:
        if s.from_bp >= depth or not alive_end[s.id]:
            continue
        alive = alive_end[s.id]
        carried = set()
        for c in s.children:
            if segs[c].from_bp >= depth:
                continue
            if not segs[c].presence:
                carried = alive
                break
            carried.update(alive.intersection(segs[c].presence[0].labels))
        for l in sorted(alive - carried):
            ends.setdefault(l, []).append((s.id, min(s.end_bp, depth)))
    return ends


def _restricted_claims(g, side, depth, inexact=None):
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
        for run in a.runs:
            c = _restricted_run_claim(g, a, run, depth, exact, inexact)
            if c is not None:
                out.append(c)
        return out
    for l, es in sorted(_annotate_route_ends_at(a, depth).items()):
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


def _walk_keys(g, sides, depth, sel, memo=None):
    memo = {} if memo is None else memo
    out = {}
    for side in sides:
        seen = set()
        for p, m, info in _cuts(g, side, depth, memo):
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


def _label_keys(g, sides, depth, sel):
    """(side, ref) -> {direct_bp, reach_bp} of every label the arm holds within [0, depth),
    each summary DERIVED from the records clipped to the depth (_clipped_summary), not a
    completed summary clamped to it: a label first recorded beyond the depth is in
    neither retrieval's [0, depth), and clamping its reach_bp would report it in the
    deeper one only."""
    out = {}
    for side in sides:
        for lid, v in _clipped_summary(g, g.arms[side], depth).items():
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


def _clipped_summary(g, arm, depth):
    """label id -> {direct_bp, reach_bp} of one arm's label_summary (§5.1, both modes)
    computed over [0, depth) only, for the labels the arm holds there. Constrain: the
    runs that exist at the cut, each cut to it (direct_bp: the seed runs from 0; reach_bp:
    credited to the lineage root). Annotate: the P runs starting before the cut, each cut
    to it, and the union pass of direct_bp stopped at the cut; a label recorded at the
    seed boundary only has both 0. At a depth beyond the arm it is label_summary itself
    (for the labels it holds)."""
    out = {}
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


def _route_pairs(g, sides, depth, memo=None):
    """(side, walk prefix cut at depth in walking order, ref) of every walk under each
    label with displayed support over the whole cut walk."""
    memo = {} if memo is None else memo
    out = set()
    for side in sides:
        seen = set()
        for p, m, info in _cuts(g, side, depth, memo):
            if id(info) in seen:
                continue                     # walks sharing the cut give the same pairs
            seen.add(id(info))
            seq = _cut_prefix(info[0])
            for l in info[1]:
                out.add((side, seq, g.labels[l].ref))
    return out


def _supported_prefixes(g, sides, depth, memo=None):
    """(side, ref) -> the walking-order prefixes w[:j] of every walk w (cut at depth) on
    which the label has displayed support over all of [0, j), j maximal per walk: a
    label that ends inside a walk (at a node where the walk goes on with other labels)
    still supports the walk up to its end, which whole-walk pairs alone never see."""
    memo = {} if memo is None else memo
    out = {}
    for side in sides:
        seen = set()
        for p, m, info in _cuts(g, side, depth, memo):
            if id(info) in seen:
                continue
            seen.add(id(info))
            w = _cut_prefix(info[0])
            for l, j in info[2].items():
                out.setdefault((side, g.labels[l].ref), set()).add(w[:j])
    return out


def _recorded_refusal(g, arm, seq, ref):
    """-> (at_bp, reason) of a successor that |g| recorded as not taken by label |ref|
    where walk |seq| (walking order) leaves g's trie: a V refusal (quorum, split limit), a
    blocked successor or a skipped hairpin, on a segment whose chain spells seq up to
    the branch and for the base seq has there. None when nothing was recorded."""
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


def _divergence(arm, seq, depth):
    """The length of the longest prefix of |seq| (walking order) that some walk of the
    arm's DAG restricted to [0, depth) spells: where seq leaves the arm's trie there."""
    best = 0
    for leaf in restricted_leaves(arm, depth):
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


def _prefix_subset(a, b, sides, depth, sel, memo=None, inexact=None):
    memo = {} if memo is None else memo
    pa = _route_pairs(a, sides, depth, memo)
    pb = _route_pairs(b, sides, depth, memo)
    if sel is not None:
        pa = {x for x in pa if x[2] in sel}
        pb = {x for x in pb if x[2] in sel}
    # a's walk under a label must be a prefix of a walk of b on which b's label supports
    # at least that prefix -- not necessarily the whole walk (it may end inside it)
    supported = _supported_prefixes(b, sides, depth, memo)
    violations = []
    for side, seq, ref in sorted(pa):
        if not any(w.startswith(seq) for w in supported.get((side, ref), ())):
            violations.append({'arm': side, 'walk': _natural(side, seq), 'ref': ref})
    omissions = []
    a_claims = {}
    for side in sides:
        arm = a.arms[side]
        for c in _restricted_claims(a, side, depth, inexact):
            if c.kind == 'merged':
                continue
            prefix = _walk_prefix(arm, c.segment, c.to_bp)
            a_claims.setdefault((side, c.label.ref), []).append((prefix, c, arm))
    for side, seq, ref in sorted(pb):
        if (side, seq, ref) in pa:
            continue
        best = None
        for prefix, c, arm in a_claims.get((side, ref), ()):
            if seq.startswith(prefix) and (best is None or len(prefix) > len(best[0])):
                best = (prefix, c, arm)
        if best is not None and len(best[0]) >= len(seq):
            continue                 # a covers the whole walk under this label
        entry = {'arm': side, 'walk': _natural(side, seq), 'ref': ref}
        arm = a.arms[side]
        ev_to = arm.evidence_complete_to_bp
        if best is None:
            # the label went on along another successor in a: the reason, when there
            # is one, is a recorded refusal of b's successor (V, or a blocked or skipped
            # E), not a claim
            got = _recorded_refusal(a, arm, seq, ref)
            if got is not None:
                entry['at_bp'], entry['reason'] = got
            else:
                div = _divergence(arm, seq, depth)
                if ev_to is not None and div is not None and div >= ev_to:
                    # beyond the evidence boundary: the refusal may have been dropped
                    entry['at_bp'], entry['reason'] = div, 'unexplained_capped'
                else:
                    entry['reason'] = 'unexplained'
        else:
            prefix, c, arm = best
            entry['at_bp'] = c.to_bp
            entry['reason'] = c.qualifier or c.reason or c.kind
            if _open_at(c, depth):
                entry['reason'] = 'open'
            elif ev_to is not None and c.to_bp >= ev_to and c.reason is None:
                entry['reason'] = 'unexplained_capped'
        omissions.append(entry)
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
    GraphletStore charges against max_ram_mb (the old shallow count saw 35-67 % of the
    traced heap, B2)."""
    return _deep_bytes([g])


def cache_bytes(g):
    """The derived caches alone (Graphlet.cache and every Arm.cache), the model's records
    they refer to not counted: what the queries added since the model was measured. The
    store re-measures this, not the whole model, after a query (cost ~ the caches)."""
    roots = [g.cache] + [a.cache for a in g.arms.values()]
    return _deep_bytes(roots, stop=_model_types())


def cache_signature(g):
    """A cheap fingerprint of the caches' sizes (the number of entries of each), which
    changes whenever a query adds to them."""
    sig = [len(g.cache)]
    for a in g.arms.values():
        sig.append(len(a.cache))
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
        if lim.kind in LOWER_BOUND_KINDS or lim.kind in OVERSTATED_KINDS:
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
    was restricted to (its completeness is qualified 'for the selected labels')."""
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
    if view is not None:
        out['view'] = view.view_spec()
        out['qualified'] = 'for the selected labels'
    return out


def summary(g, arm=None, max_bytes=2048):
    """Agent-facing, <= max_bytes of JSON: per-arm status / completeness / counts, the
    top labels by direct_bp, the stated limitations (caveats: kind, knob, limit,
    observed, complete_to_bp and effect; informational ones flagged) and the resource
    stop with its suggested actions. Under the cap the effect texts go first, then top
    labels, then caveats (each step stated)."""
    g.require_envelope('summary()')
    sides = [g.arm(arm).side] if arm is not None else [s for s in ARM_SIDES if s in g.arms]
    rows = derive.label_summary(g)
    out = {'seed': {'seed_id': (g.seed_summary.get('seed') or {}).get('seed_id', ''),
                    'length_bp': g.seed.length_bp, 'mode': g.mode,
                    'labels': len(g.labels)},
           'outcome': g.outcome.as_dict(),
           'index': {'ns': g.index_ns, 'verifiable': g.index_fp is not None},
           'arms': {}, 'caveats': []}
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
