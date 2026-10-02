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
    zero-length runs are boundary claims; terminal loss/branches are unknown at a cut.
"""

import copy
import hashlib
import json
import math
import sys

from . import derive
from ._codec import REASON, RESOURCE_CODES, UNLIMITED
from .model import (
    AmbiguousLabel, ARM_SIDES, BadSelector, Change, Claim, Comparison, Continuation,
    IncompleteRecording, Label, LabelWalk, Path, SupportRun, UnknownLabel, Walk,
)

__all__ = [
    'label', 'labels_matching', 'path_id', 'leaf_segment', 'spell', 'walks', 'claims',
    'label_walks', 'routes', 'support_profile', 'support_changes', 'label_summary',
    'continuation', 'next_request', 'subgraph', 'GraphletView', 'view_from_saved',
    'compare', 'index_identity', 'memory_bytes', 'summary', 'evidence_block',
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
        s = arm.segments[sid]
        if s.length_bp == 0:
            continue
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
    return _merge_runs(out)


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


def _run_claim(g, arm, run, cut):
    """The claim of one run end, or None when the run-start guard drops it at |cut|."""
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
        cost=run.cost, loss=loss, branches=branches, support=g.support, exact=True,
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
        for s in segs:
            alive = set(s.entry) if not s.parents else set(got[s.parents[0]][0])
            drops = []
            for pr in s.presence:
                if not alive:
                    break
                still = alive.intersection(pr.labels)
                drops.extend((l, pr.from_bp) for l in sorted(alive - still))
                alive = still
            got[s.id] = (frozenset(alive), drops)
        arm.cache['displayed_presence'] = got
    return got


def _annotate_stretches(g, arm):
    """Annotate mode: per walk and per label recorded at the seed boundary, the stretch
    from 0 on which the label is recorded at every node (the §6.9 oracle E, filtered
    from the recorded sets alone). -> [(label, anchor, j, path_id, leaf)], maximal per
    label. A stretch ends inside its anchor (shared by every walk through it), at a
    leaf's end, or at the first base of a child -- where it is not maximal when another
    first-parent child of the anchor carries the label on."""
    got = arm.cache.get('annotate_stretches')
    if got is not None:
        return got
    segs = arm.segments
    pres = _displayed_presence(arm)
    first = _first_paths(arm)
    pol = derive.path_of_leaf(arm)
    out = []
    for s in segs:
        alive_end, drops = pres[s.id]
        pid = first[s.id]
        if pid is None:
            continue                     # on no displayed walk (a non-first merge parent)
        for l, j in drops:
            anchor = s.id
            if j == s.from_bp and s.parents:
                # absent from the child's first base: the stretch ends with the parent,
                # maximal only if no other walk through the parent carries l on
                anchor = s.parents[0]
                if any(_carries(arm, c, l) for c in segs[anchor].children
                       if c != s.id and segs[c].parents[0] == anchor):
                    continue
            out.append((l, anchor, j, pid, derive.paths(arm)[pid].leaf))
        if s.leaf is not None:
            for l in sorted(alive_end):
                out.append((l, s.id, s.end_bp, pol[s.id], s.id))
    # one claim per (label, anchor, position): walks sharing a prefix share the stretch
    seen = set()
    uniq = []
    for x in out:
        if (x[0], x[1], x[2]) not in seen:
            seen.add((x[0], x[1], x[2]))
            uniq.append(x)
    got = sorted(uniq, key=lambda x: (x[3], x[0], x[2]))
    arm.cache['annotate_stretches'] = got
    return got


def _carries(arm, c, l):
    """Whether a walk through child |c| has |l| recorded at c's first base (through
    zero-length segments, which record nothing, to their first-parent children)."""
    seg = arm.segments[c]
    if seg.presence:
        return l in seg.presence[0].labels
    return any(_carries(arm, cc, l) for cc in seg.children
               if arm.segments[cc].parents[0] == c)


def _annotate_claim(g, arm, l, anchor, j, pid, leaf, cut):
    """One oracle stretch [0, j) as a claim; None when the run-start guard drops it
    (a stretch with bases makes no claim at a cut of 0)."""
    if cut is not None and j > 0 and cut <= 0:
        return None
    seg = arm.segments[leaf]
    length = seg.end_bp
    if j < length:
        reason, end_class, kind = 'label_lost', 'lost', 'end'
    else:
        code = seg.leaf.path_reason
        reason, end_class = REASON[code], END_CLASS[code]
        kind = 'alive' if code == 'X' or code in RESOURCE_CODES else 'end'
    to_bp = j
    ev = 0
    if j == 0:
        kind, ev = 'boundary', None
    if cut is not None and j > cut:
        to_bp, kind, end_class, reason = cut, 'stretch', 'open', None
    return Claim(label=g.labels[l], arm=arm.side, path_id=pid, segment=anchor, from_bp=0,
                 to_bp=to_bp, route_bp=0, evidence_from=ev, kind=kind, reason=reason,
                 qualifier=None, end_class=end_class, entered_by='seed', from_label=None,
                 cost=None, loss=None, branches=None, support=g.support,
                 exact=arm.labels_per_node.nodes_truncated == 0, run=None)


def claims(g, arm=None, labels=None, at_most_bp=None, strict=True):
    """The §6.9 unit: one claim per RUN END (claims, not leaves). In annotate mode the
    claims are the oracle's: per label recorded at the seed boundary its maximal
    label-consistent walks; strict=True raises IncompleteRecording when a recorded
    list was cut (the oracle would be a lower bound)."""
    sides = [g.arm(arm).side] if arm is not None else [s for s in ARM_SIDES if s in g.arms]
    wanted = _label_ids(g, labels)
    out = []
    for side in sides:
        a = g.arms[side]
        if g.mode == 'constrain':
            for run in a.runs:
                if wanted is not None and run.label not in wanted:
                    continue
                c = _run_claim(g, a, run, at_most_bp)
                if c is not None:
                    out.append(c)
        else:
            if strict and a.labels_per_node.nodes_truncated:
                raise IncompleteRecording(
                    'the %s arm has %d recorded label lists cut at labels.max_labels_per_node:'
                    ' claims filtered from them are lower bounds (pass strict=False to accept,'
                    ' or raise the knob)' % (side, a.labels_per_node.nodes_truncated))
            for l, anchor, j, pid, leaf in _annotate_stretches(g, a):
                if wanted is not None and l not in wanted:
                    continue
                c = _annotate_claim(g, a, l, anchor, j, pid, leaf, at_most_bp)
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
    wanted = _label_ids(g, labels)
    ranked = []
    for p in derive.paths(a):
        if p.length_bp < min_bp:
            continue
        seg = a.segments[p.leaf]
        if g.mode == 'constrain':
            ends = derive.end_labels(a, p.leaf)
            alive_ids = [e.label for e in ends]
            consistent = {e.label for e in ends if e.route_bp == 0}
            min_loss = min((e.loss for e in ends), default=math.inf)
            full = None
        else:
            alive_ids = list(seg.end)
            full = _labels_full(g, a, p)
            consistent = {l.id for l in full}
            min_loss = 0.0
        if wanted is not None:
            pool = consistent if route_consistent else set(alive_ids)
            if not wanted & pool:
                continue
        if full is None:
            full = _labels_full(g, a, p)
        ranked.append((p, full, alive_ids, min_loss))
    keys = {
        'support': lambda x: (-len(x[1]), -len(x[2]), -x[0].length_bp, x[0].id),
        'length': lambda x: (-x[0].length_bp, x[0].id),
        'loss': lambda x: (x[3], -len(x[2]), -x[0].length_bp, x[0].id),
        'id': lambda x: x[0].id,
    }
    if by not in keys:
        raise ValueError("by is 'support', 'length', 'loss' or 'id'")
    ranked.sort(key=keys[by])
    if top is not None:
        ranked = ranked[:top]
    by_leaf = {}
    if ranked:
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
    bases (natural orientation, from the seed boundary to the run's end)."""
    a = g.arm(arm)
    lab = label(g, sel)
    out = []
    if g.mode != 'constrain':
        for l, anchor, j, pid, leaf in _annotate_stretches(g, a):
            if l == lab.id:
                r = derive.chain(a, anchor)
                out.append((r, _chain_bases(a, anchor, 0, j)) if spell else r)
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
    leaves_below, merged_into)."""
    lab = label(g, sel)
    sides = [g.arm(arm).side] if arm is not None else [s for s in ARM_SIDES if s in g.arms]
    out = []
    for side in sides:
        a = g.arms[side]
        if g.mode != 'constrain':
            for l, anchor, j, pid, leaf in _annotate_stretches(g, a):
                if l != lab.id:
                    continue
                try:
                    seq = _chain_bases(a, anchor, 0, j)
                except ValueError:
                    seq = None
                c = _annotate_claim(g, a, l, anchor, j, pid, leaf, None)
                out.append(LabelWalk(None, lab, side, 0, j, 0, seq, c.reason,
                                     [leaf], None))
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
                    if l in parent.end:
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
    """Continuation(sequence, labels, loss_used, branches_used, seed_coord, as_seed()):
    the walker's label logic from C, the spelling derived (both arms). seed_coord is
    the half-open whole-molecule interval in seed coordinates (seed = [0, |seed|))."""
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
    return Continuation(seq, _labels(g, c.labels), c.loss_used, c.branches_used, coord,
                        a.side, path_id(a, leaf))


def _deep_merge(dst, src):
    for k, v in src.items():
        if isinstance(v, dict) and isinstance(dst.get(k), dict):
            _deep_merge(dst[k], v)
        else:
            dst[k] = copy.deepcopy(v)


def next_request(g, arm, leaves, bp=None, reduce_budget=True, **overrides):
    """A resubmittable /traverse request continuing the given walks: one seed per walk
    (its continuation, natural orientation, with its labels named explicitly), the
    retrieval's normalized strategy with direction = this arm, max_extension_bp = bp
    when given, and (reduce_budget) the loss budget reduced by the largest loss_used
    among the walks -- one strategy serves every seed, and no lineage may exceed its
    original budget. Keyword overrides deep-merge into the strategy; release, graph and
    graph_path are request-level. Raises MissingEnvelope on a body-only graphlet."""
    g.require_envelope('next_request()')
    a = g.arm(arm)
    if isinstance(leaves, (int, Path, Walk)):
        leaves = [leaves]
    conts = [continuation(g, a, x) for x in leaves]
    seeds = []
    for c in conts:
        if not c.sequence:
            raise ValueError('walk %d has no continuation sequence (output.continuation_bp '
                             'was 0): nothing to resubmit' % c.leaf)
        s = c.as_seed()
        if g.mode != 'constrain':
            s.pop('labels', None)    # annotate mode has no permitted set
        seeds.append(s)
    strategy = copy.deepcopy(g.envelope.get('strategy') or {})
    strategy.pop('clamped', None)
    strategy['direction'] = a.side
    if bp is not None:
        strategy.setdefault('bounds', {})['max_extension_bp'] = bp
    if reduce_budget and g.mode == 'constrain':
        labels_ = strategy.get('labels') or {}
        budget = labels_.get('loss_budget')
        used = max((c.loss_used for c in conts if c.loss_used != math.inf), default=0.0)
        if isinstance(budget, (int, float)) and budget > 0 and used > 0:
            labels_['loss_budget'] = max(0.0, budget - used)
            strategy['labels'] = labels_
    request = {'seeds': seeds}
    release = g.envelope.get('release')
    if release:
        request['release'] = release
    for k in ('release', 'graph', 'graph_path'):
        if k in overrides:
            request[k] = overrides.pop(k)
    _deep_merge(strategy, overrides)
    request['strategy'] = strategy
    return request


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
    listed with the reason |a| recorded for them)."""
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
    sel = None
    if labels is not None:
        sel = {label(a, x).ref for x in ([labels] if isinstance(labels, (str, dict))
                                          else labels)}
    only_a, only_b, differ = [], [], []
    names = {l.ref: l.name for l in b.labels}
    names.update({l.ref: l.name for l in a.labels})

    def lab(ref):
        return {'name': names.get(ref), 'ref': ref}

    if mode == 'claims':
        ka, kb = _claim_keys(a, sides, depth, sel), _claim_keys(b, sides, depth, sel)
        only_a = [_key_json(k, v, lab) for k, v in sorted(ka.items()) if k not in kb]
        only_b = [_key_json(k, v, lab) for k, v in sorted(kb.items()) if k not in ka]
        differ = [{'key': _key_json(k, None, lab), 'a': ka[k], 'b': kb[k]}
                  for k in sorted(ka) if k in kb and ka[k] != kb[k]]
    elif mode == 'walks':
        wa, wb = _walk_keys(a, sides, depth, sel), _walk_keys(b, sides, depth, sel)
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
    elif mode == 'prefix_subset':
        only_a, only_b = _prefix_subset(a, b, sides, depth, sel)
        for row in only_a + only_b:
            row['name'] = names.get(row['ref'])
    else:
        raise ValueError("mode is 'claims', 'walks', 'labels' or 'prefix_subset'")
    if mode == 'prefix_subset':
        holds = not only_a
    else:
        holds = not (only_a or only_b or differ)
    equal = holds if comparable is True else None
    if comparable is not True and holds:
        notes.append('no difference found, but equality is not claimed (%s)' % comparable)
    return Comparison(comparable, reason, depth, equal, only_a, only_b, notes, differ,
                      mode=mode, support=(a.support, b.support), scopes=scopes,
                      strategies=strategies)


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


def _claim_keys(g, sides, depth, sel):
    out = {}
    for side in sides:
        a = g.arms[side]
        for c in claims(g, side, at_most_bp=depth, strict=False):
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


def _cut_labels(g, a, p, m):
    """Labels with displayed support on the whole [0, m) of walk p."""
    if m == 0:
        return {g.labels[l].ref for l in a.segments[p.segments[0]].entry}
    prof = support_profile(g, a, p.id)
    have = None
    for run in prof:
        if run.from_bp >= m:
            break
        ids = {l.id for l in run.labels}
        have = ids if have is None else have & ids
    return {g.labels[l].ref for l in (have or ())}


def _walk_keys(g, sides, depth, sel):
    out = {}
    for side in sides:
        a = g.arms[side]
        for p in derive.paths(a):
            m = min(p.length_bp, depth)
            try:
                seq = _chain_bases(a, p.leaf, 0, m)
            except ValueError:
                seq = 'path:%d' % p.id
            refs = _cut_labels(g, a, p, m)
            if sel is not None:
                refs &= sel
            out.setdefault((side, seq), set()).update(refs)
    return out


def _label_keys(g, sides, depth, sel):
    rows = derive.label_summary(g)
    out = {}
    for lab, per in zip(g.labels, rows):
        if sel is not None and lab.ref not in sel:
            continue
        for side in sides:
            v = per.get(side)
            if v is None or (v['reach_bp'] == 0 and v['direct_bp'] == 0 and not v['runs']
                             and not _label_seen(g, side, lab.id)):
                continue
            out[(side, lab.ref)] = {'direct_bp': min(v['direct_bp'], depth),
                                    'reach_bp': min(v['reach_bp'], depth)}
    return out


def _label_seen(g, side, lid):
    """Recorded at the seed boundary of the arm (a label with no stretch of its own)."""
    a = g.arms[side]
    return bool(a.segments) and lid in a.segments[0].entry


def _walk_prefix(arm, seg_id, hi):
    """Outward bases [0, hi) along seg's first-parent chain in WALKING order. Prefix tests
    read walks outward on both arms: the natural spelling of the left arm is the reverse
    of the walking order, where an outward prefix is a suffix."""
    return derive.walk_bases(arm, seg_id)[:hi]


def _natural(side, w):
    return w if side == 'right' else w[::-1]


def _route_pairs(g, sides, depth):
    """(side, walk prefix cut at depth in walking order, ref) of every walk under each
    label with displayed support over the whole cut walk."""
    out = set()
    for side in sides:
        a = g.arms[side]
        for p in derive.paths(a):
            m = min(p.length_bp, depth)
            seq = _walk_prefix(a, p.leaf, m)
            for ref in _cut_labels(g, a, p, m):
                out.add((side, seq, ref))
    return out


def _supported_prefixes(g, sides, depth):
    """(side, ref) -> the walking-order prefixes w[:j] of every walk w (cut at depth) on
    which the label has displayed support over all of [0, j), j maximal per walk: a
    label that ends inside a walk (at a node where the walk goes on with other labels)
    still supports the walk up to its end, which whole-walk pairs alone never see."""
    out = {}
    for side in sides:
        a = g.arms[side]
        for p in derive.paths(a):
            m = min(p.length_bp, depth)
            w = _walk_prefix(a, p.leaf, m)
            have = {l: 0 for l in a.segments[p.segments[0]].entry} if p.segments else {}
            cur = None
            for r in support_profile(g, a, p.id):
                if r.from_bp >= m:
                    break
                ids = {l.id for l in r.labels}
                cur = ids if cur is None else cur & ids
                for l in cur:
                    have[l] = min(r.to_bp, m)
            for l, j in have.items():
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


def _divergence(arm, seq):
    """The length of the longest prefix of |seq| (walking order) that some walk of the
    arm spells: where seq leaves the arm's trie."""
    best = 0
    for p in derive.paths(arm):
        try:
            w = derive.walk_bases(arm, p.leaf)
        except ValueError:
            return None
        n = 0
        for x, y in zip(w, seq):
            if x != y:
                break
            n += 1
        best = max(best, n)
    return best


def _prefix_subset(a, b, sides, depth, sel):
    pa = _route_pairs(a, sides, depth)
    pb = _route_pairs(b, sides, depth)
    if sel is not None:
        pa = {x for x in pa if x[2] in sel}
        pb = {x for x in pb if x[2] in sel}
    # a's walk under a label must be a prefix of a walk of b on which b's label supports
    # at least that prefix -- not necessarily the whole walk (it may end inside it)
    supported = _supported_prefixes(b, sides, depth)
    violations = []
    for side, seq, ref in sorted(pa):
        if not any(w.startswith(seq) for w in supported.get((side, ref), ())):
            violations.append({'arm': side, 'walk': _natural(side, seq), 'ref': ref})
    omissions = []
    a_claims = {}
    for side in sides:
        arm = a.arms[side]
        for c in claims(a, side, at_most_bp=depth, strict=False):
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
                div = _divergence(arm, seq)
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

def memory_bytes(g):
    """An estimate of the model's footprint (interned sets counted once)."""
    seen = set()
    total = sys.getsizeof(g)

    def add(x):
        nonlocal total
        if id(x) in seen:
            return
        seen.add(id(x))
        total += sys.getsizeof(x)

    for lab in g.labels:
        add(lab)
        add(lab.name)
    add(g.seed.sequence)
    for a in g.arms.values():
        for s in a.segments:
            add(s)
            for x in (s.entry, s.end, s.walk, s.parents, s.children, s.events, s.presence):
                if x is not None:
                    add(x)
            for p in s.presence:
                add(p)
                add(p.labels)
            for e in s.events:
                add(e)
            for p in s.partition:
                add(p)
            if s.leaf is not None:
                add(s.leaf)
                add(s.leaf.extras)
        for r in a.runs:
            add(r)
        for b in a.branch_events:
            add(b)
    return total


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


def evidence_block(g, side=None, view=None):
    """The block every local answer carries (§6): what it is certified to and how. Beyond
    §6's {complete_to_bp, exact, support, reconverge, scope}: the four outcome dimensions
    and the kinds of the limitations that apply (seed level and per arm), so an answer
    read on its own never looks more certain than the retrieval. exact holds only when
    the label evidence is complete AND no recorded list was cut. |view|: a GraphletView
    the answer was restricted to (its completeness is qualified 'for the selected
    labels')."""
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
    label_evidence = g.outcome.label_evidence if g.outcome is not None else 'complete'
    out = {
        'complete_to_bp': {s: g.arms[s].complete_to_bp for s in sides},
        'exact': label_evidence == 'complete'
        and all(g.arms[s].labels_per_node.nodes_truncated == 0 for s in sides),
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
