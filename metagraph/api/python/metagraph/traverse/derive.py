"""Everything the format derives instead of storing (DESIGN-traverse-graphlet.md §2.3,
§2.5, §5.1): the '*' rules, the trie structure (children, leaves, paths, splits), the
events that are not stored (label_end from R, reconverge from G), leaf end labels,
label_summary in both modes, the displayed-path evidence of a run, continuation
spellings, and check_rules().

These are part of the format contract: the library must not need the C++ to reproduce
them. Per-arm results are cached in Arm.cache (a graphlet is immutable once parsed).
"""

from array import array

from ._codec import RESOURCE_CODES, REASON
from .model import EndLabel, Path, Split

__all__ = [
    'entry_base', 'rule_entry', 'rule_end', 'rule_partition', 'runs_by_segment',
    'leaves', 'paths', 'path_of_leaf', 'chain', 'splits', 'split_branches',
    'label_end_events', 'reconverge_events', 'end_labels', 'labels_at_end',
    'label_summary', 'needed_budgets', 'continuation_sequence', 'walk_bases',
    'natural_flank', 'evidence', 'lineage_label_at', 'merge_above', 'run_leaves',
    'check_rules',
]


def _union(sets):
    out = set()
    for s in sets:
        out.update(s)
    return sorted(out)


# ------------------------------------------------------------------ the §2.3 rules

def entry_base(segments, seg):
    """The SETEXPR base of G.entry (§2.1, v2): none for the root, the union of the
    parents' end sets for a merge, parent[0]'s end set otherwise."""
    if not seg.parents:
        return None
    if len(seg.parents) > 1:
        return _union(segments[p].end for p in seg.parents)
    return list(segments[seg.parents[0]].end)


def rule_entry(mode, num_seed_labels, seg):
    """G.entry '*' (constrain only): the root's state is the seed labels
    (walker.cpp init_arm), a merge's entry is the union of its partition. None where no
    rule applies (annotate mode; constrain split children, whose entry is primary)."""
    if mode != 'constrain':
        return None
    if not seg.parents:
        return list(range(num_seed_labels))
    if len(seg.parents) > 1:
        return _union(seg.partition)
    return None


def rule_partition(parents):
    """G.partition '*': no merge, or one empty list per parent (annotate merges, where
    labels_via_parent is empty: v3 states the arity)."""
    if len(parents) <= 1:
        return []
    return [array('I') for _ in parents]


def rule_end(mode, seg, entry, presence, runs_here, runs):
    """G.end '*' (v2, chronological). Constrain: from the entry set apply, in increasing
    position (ties: ends before switch-ins), every run end anchored here with
    to_bp < from_bp + length_bp (remove its label) and every switch event here (add its
    target). Runs ending at the segment's last node are still alive there
    (walker.cpp:2580). v1's set formula lost a label that re-entered after ending
    (A -> B -> A inside one segment). Annotate: the last P set, or the entry set."""
    if mode != 'constrain':
        return list(presence[-1].labels) if presence else list(entry)
    end_bp = seg.from_bp + seg.length_bp
    ops = []
    for r in runs_here:
        run = runs[r]
        if run.to_bp < end_bp:
            ops.append((run.to_bp, 0, run.label))
    for ev in seg.events:
        if ev.type == 'switch':
            ops.append((ev.at_bp, 1, ev.to_label))
    if not ops:
        return list(entry)
    ops.sort(key=lambda o: (o[0], o[1]))
    s = set(entry)
    for _, kind, label in ops:
        if kind == 0:
            s.discard(label)
        else:
            s.add(label)
    return sorted(s)


# ------------------------------------------------------------------ trie structure

def runs_by_segment(arm):
    got = arm.cache.get('runs_by_segment')
    if got is None:
        got = [[] for _ in arm.segments]
        for r in arm.runs:
            got[r.segment].append(r.id)
        arm.cache['runs_by_segment'] = got
    return got


def leaves(arm):
    got = arm.cache.get('leaves')
    if got is None:
        got = [s.id for s in arm.segments if s.leaf is not None]
        arm.cache['leaves'] = got
    return got


def chain(arm, seg_id):
    """The first-parent chain root -> seg (a path's segments, walker.cpp finalize)."""
    out = []
    s = seg_id
    segs = arm.segments
    while True:
        out.append(s)
        p = segs[s].parents
        if not p:
            break
        s = p[0]
    out.reverse()
    return out


def paths(arm):
    """paths[]: leaf ordinal in segment-id order, chained through parents[0], length =
    the leaf's from_bp + length_bp -- exactly finalize(), so path ids are the server's."""
    got = arm.cache.get('paths')
    if got is None:
        got = [Path(i, leaf, arm.segments[leaf].end_bp, arm)
               for i, leaf in enumerate(leaves(arm))]
        arm.cache['paths'] = got
    return got


def path_of_leaf(arm):
    got = arm.cache.get('path_of_leaf')
    if got is None:
        got = {p.leaf: p.id for p in paths(arm)}
        arm.cache['path_of_leaf'] = got
    return got


def splits(arm, mode):
    """splits[]: the single-parent children grouped by parent (a split always has >= 2
    children; nf == 1 continues the segment), ordered by (at_bp, first child id) -- the
    walker's processing order. kind ambiguous <=> the parent's stored G.split (v2)."""
    key = 'splits'
    got = arm.cache.get(key)
    if got is None:
        groups = {}
        for s in arm.segments:
            if len(s.parents) == 1:
                groups.setdefault(s.parents[0], []).append(s.id)
        got = []
        for parent, children in groups.items():
            p = arm.segments[parent]
            if mode == 'constrain':
                before = len(p.end)
            else:
                before = p.presence[-1].total if p.presence else p.entry_total
            got.append(Split(p.end_bp, parent, children, bool(p.split), before))
        got.sort(key=lambda sp: (sp.at_bp, sp.children[0]))
        arm.cache[key] = got
    return got


def first_base(seg):
    if seg.walk:
        return seg.walk[0]
    return seg.first_base


def split_branches(arm, split, cap):
    """Per branch: its first base, the labels at its first node (the child's entry, the
    true count and the list cut at the cap)."""
    out = []
    for c in split.children:
        child = arm.segments[c]
        out.append({'segment': c, 'char': first_base(child),
                    'labels_distinct': child.entry_total, 'labels': list(child.entry[:cap])})
    return out


# ------------------------------------------------------------------ derived events

def label_end_events(arm):
    """seg -> [run ids] of runs ended with a label_end event there (R end not Lw/m;
    at = to_bp)."""
    got = arm.cache.get('label_end_events')
    if got is None:
        got = {}
        for r in arm.runs:
            if not r.silent and not r.merged:
                got.setdefault(r.segment, []).append(r.id)
        arm.cache['label_end_events'] = got
    return got


def reconverge_events(arm):
    """seg -> (at_bp, parents): every G with more than one parent (walker.cpp
    merge_level: at the merged segment's from_bp, the segments joined)."""
    return {s.id: (s.from_bp, list(s.parents)) for s in arm.segments if len(s.parents) > 1}


def end_labels(arm, leaf):
    """The labels alive at a leaf, ascending (v2): the runs anchored at the leaf with
    to_bp at its end (every label alive at the leaf has an ended run there, and every
    run ending at its last node is such a label), with T's (loss, branches, route_bp)."""
    cache = arm.cache.setdefault('end_labels', {})
    got = cache.get(leaf)
    if got is None:
        seg = arm.segments[leaf]
        extras = seg.leaf.extras if seg.leaf else {}
        got = []
        for r in runs_by_segment(arm)[leaf]:
            run = arm.runs[r]
            if run.to_bp == seg.end_bp and not run.merged and not run.silent:
                loss, branches, route = extras.get(run.label, (0.0, 0, 0))
                got.append(EndLabel(run.label, loss, branches, r, route))
        got.sort(key=lambda e: e.label)
        cache[leaf] = got
    return got


def labels_at_end(arm, seg_id):
    return list(arm.segments[seg_id].end)


def needed_budgets(arm):
    """The needed_budget of every loss_budget end (a histogram: order is not
    information, §2.5); here in run order."""
    return [r.needed_budget for r in arm.runs
            if r.code == 'B' and r.needed_budget is not None]


def run_leaves(arm, run):
    """The leaves below a run's anchor (the walks that include its end)."""
    out = []
    stack = [run.segment]
    segs = arm.segments
    seen = set()
    while stack:
        s = stack.pop()
        if s in seen:
            continue
        seen.add(s)
        if segs[s].leaf is not None:
            out.append(s)
        stack.extend(segs[s].children)
    return sorted(out)


# ------------------------------------------------------------------ label summary

def label_summary(g):
    """Per label, per requested arm: {direct_bp, reach_bp, reentries, runs}, exactly
    Walker::summarize() (constrain) and Walker::summarize_annotate() (§5.1, normative)."""
    got = g.cache.get('label_summary')
    if got is not None:
        return got
    n = len(g.labels)
    out = [{} for _ in range(n)]
    for side, arm in g.arms.items():
        per = [{'direct_bp': 0, 'reach_bp': 0, 'reentries': 0, 'runs': []} for _ in range(n)]
        if g.mode == 'constrain':
            root = [0] * len(arm.runs)
            for r in arm.runs:
                # prev_run(r) < r always, so one forward pass resolves switch chains
                if r.from_label is not None and r.prev_run is not None:
                    root[r.id] = root[r.prev_run]
                else:
                    root[r.id] = r.label
            for r in arm.runs:
                own = per[r.label]
                own['runs'].append(r.id)
                if r.from_label is None and r.from_bp == 0:
                    # measured along the label's OWN route: a merge does not clamp it
                    own['direct_bp'] = max(own['direct_bp'], r.to_bp)
                lineage = per[root[r.id]]
                lineage['reach_bp'] = max(lineage['reach_bp'], r.to_bp)
            # switch EVENTS, not switched runs: split clones copy from_label and cost
            # without a second switch
            for s in arm.segments:
                for ev in s.events:
                    if ev.type == 'switch':
                        per[ev.to_label]['reentries'] += 1
        else:
            for s in arm.segments:
                for p in s.presence:
                    for l in p.labels:
                        if p.to_bp > per[l]['reach_bp']:
                            per[l]['reach_bp'] = p.to_bp
            # parents precede children (segments are created parents first), so one
            # pass in id order with the UNION of the parents' surviving sets visits
            # every segment once (walking the routes is exponential on a merged DAG)
            alive_end = [None] * len(arm.segments)
            for s in arm.segments:
                if not s.parents:
                    alive = set(s.entry)
                else:
                    alive = set()
                    for p in s.parents:
                        alive.update(alive_end[p])
                for p in s.presence:
                    if not alive:
                        break
                    still = alive.intersection(p.labels)
                    for l in alive - still:
                        if p.from_bp > per[l]['direct_bp']:
                            per[l]['direct_bp'] = p.from_bp
                    alive = still
                end = s.end_bp
                for l in alive:
                    if end > per[l]['direct_bp']:
                        per[l]['direct_bp'] = end
                alive_end[s.id] = alive
        for l in range(n):
            out[l][side] = per[l]
    g.cache['label_summary'] = out
    return out


# ------------------------------------------------------------------ spelling

def walk_bases(arm, leaf):
    """The flank in WALKING order, root -> leaf (outward index i is [i])."""
    segs = arm.segments
    parts = []
    for s in chain(arm, leaf):
        w = segs[s].walk
        if w is None:
            raise ValueError('this retrieval carries no bases (output.sequences: false)')
        parts.append(w)
    return ''.join(parts)


def natural_flank(arm, leaf):
    """Natural orientation (§2.4): right flank = walking order; left flank = its
    reverse (= concat leaf -> root of the per-segment reversals, spell_path)."""
    w = walk_bases(arm, leaf)
    return w if arm.side == 'right' else w[::-1]


def continuation_sequence(g, arm, leaf):
    """The C record's spelling (v2, both arms stated): right arm: the LAST n bases of
    seed + natural(right flank); left arm: the FIRST n bases of natural(left flank) +
    seed (n <= |seed| + flank; a continuation may cross into the seed)."""
    seg = arm.segments[leaf]
    c = seg.leaf.continuation if seg.leaf else None
    if c is None:
        return None
    if c.sequence is not None:
        return c.sequence
    if c.n == 0:
        return ''
    flank = natural_flank(arm, leaf)
    if arm.side == 'right':
        whole = g.seed.sequence + flank
        return whole[len(whole) - c.n:]
    whole = flank + g.seed.sequence
    return whole[:c.n]


# ------------------------------------------------------------------ evidence (§5.1)

def lineage_label_at(arm, run, depth):
    """The lineage's INCOMING label at |depth|: its label before any switch committed
    at that depth. A switch committed by the step out of depth d starts its run at d
    (walker.cpp commit_entries), after the merge at d, so the incoming run is the one
    with from_bp < d; prev_run walks back across switches."""
    r = run
    runs = arm.runs
    while r.from_bp >= depth and r.from_label is not None and r.prev_run is not None:
        r = runs[r.prev_run]
    return r.label


def merge_above(arm):
    """seg -> the nearest merge segment at or above it on its first-parent chain (None:
    none), so that the evidence scan visits the merges only."""
    got = arm.cache.get('merge_above')
    if got is None:
        got = []
        for seg in arm.segments:                 # parents precede children
            if len(seg.parents) > 1:
                got.append(seg.id)
            else:
                got.append(got[seg.parents[0]] if seg.parents else None)
        arm.cache['merge_above'] = got
    return got


def evidence(arm, run):
    """-> (route_from, evidence_from) of a run (§5.1, v3 normative).

    Provenance is derived at the run's ANCHORED endpoint (also for merge-closed 'm'
    runs) and never by evaluating the walk at a cut depth: with p = the first-parent
    chain root -> anchor and t = to_bp, the latest merge on p at or before t whose
    first-parent partition does not hold the lineage's incoming label is where the
    displayed bases start to be this lineage's (route_from; 0 when none). For a leaf
    label route_from equals T.route_bp. evidence_from = max(route_from, from_bp): a
    label never claims bases before its own run."""
    cache = arm.cache.setdefault('evidence', {})
    got = cache.get(run.id)
    if got is not None:
        return got
    segs = arm.segments
    t = run.to_bp
    route_from = 0
    above = merge_above(arm)
    m = above[run.segment]
    while m is not None:
        seg = segs[m]
        if seg.from_bp <= t:
            label_at = lineage_label_at(arm, run, seg.from_bp)
            if label_at not in seg.partition[0]:
                route_from = seg.from_bp
                break
        m = above[seg.parents[0]]
    got = (route_from, max(route_from, run.from_bp))
    cache[run.id] = got
    return got


# ------------------------------------------------------------------ check_rules

def check_rules(g, reference=None):
    """The §2.3 rules recomputed and compared, plus the §4/§9 invariants of the body.

    Returns a list of findings, each {arm, segment|run, field, kind, ...}; [] means
    nothing diverged. Kinds:
      'star_mismatch'  a field the source wrote as '*' decodes to a value other than
                       |reference| (a detail: full results[i]) holds -- the writer
                       emitted '*' where its rule does not reproduce the walker;
      'noncanonical'   an explicit field equal to its rule's value (the canonical form
                       writes '*'): valid, but not what a conforming writer emits;
      'invariant'      a structural invariant of §4/§9 that does not hold.
    """
    out = []
    out.extend(_outcome_findings(g))
    ref_arms = (reference or {}).get('arms') if reference else None
    for side, arm in g.arms.items():
        segs = arm.segments
        rbs = runs_by_segment(arm)
        ref = ref_arms.get(side) if ref_arms else None
        ref_segs = {s['id']: s for s in ref['segments']} if ref and 'segments' in ref else {}
        for s in segs:
            # explicit values equal to their rule: not canonical
            re_ = rule_entry(g.mode, g.seed.num_seed_labels, s)
            if 'entry' not in s.stars and re_ is not None and list(s.entry) == re_:
                out.append({'arm': side, 'segment': s.id, 'field': 'entry',
                            'kind': 'noncanonical'})
            if 'entry_total' not in s.stars and s.entry_total == len(s.entry):
                out.append({'arm': side, 'segment': s.id, 'field': 'entry_total',
                            'kind': 'noncanonical'})
            rend = rule_end(g.mode, s, s.entry, s.presence, rbs[s.id], arm.runs)
            if 'end' not in s.stars and list(s.end) == rend:
                out.append({'arm': side, 'segment': s.id, 'field': 'end',
                            'kind': 'noncanonical'})
            rs = ref_segs.get(s.id)
            if rs is not None:
                for fld, mine, theirs in (
                        ('entry', list(s.entry), rs.get('labels')),
                        ('entry_total', s.entry_total, rs.get('labels_total')),
                        ('end', list(s.end), rs.get('labels_at_end')),
                        ('partition', [list(p) for p in s.partition],
                         rs.get('labels_via_parent'))):
                    if fld in s.stars and theirs is not None and mine != theirs \
                            and not (fld == 'partition' and not mine and not theirs):
                        out.append({'arm': side, 'segment': s.id, 'field': fld,
                                    'kind': 'star_mismatch', 'decoded': mine,
                                    'reference': theirs})
        if g.mode != 'constrain':
            continue
        # §4 / T37 invariants
        n_events = sum(len(v) for v in label_end_events(arm).values())
        if n_events > len(arm.runs):
            out.append({'arm': side, 'field': 'runs', 'kind': 'invariant',
                        'msg': '#label_end events > #runs'})
        successors = {}
        for o in arm.runs:
            if o.prev_run is not None and o.from_label is not None:
                successors.setdefault(o.prev_run, []).append(o)
        merge_heads = {(p, s.from_bp) for s in segs if len(s.parents) > 1 for p in s.parents}
        for r in arm.runs:
            if r.silent:
                # a silent Lw end is the prev_run of a run entered by switch from its
                # label at its to_bp (the switch may sit on the anchor's children)
                if not any(o.from_label == r.label and o.from_bp == r.to_bp
                           for o in successors.get(r.id, ())):
                    out.append({'arm': side, 'run': r.id, 'field': 'end', 'kind': 'invariant',
                                'msg': 'Lw run is no switched run\'s prev_run at its to_bp'})
            elif r.merged:
                if (r.segment, r.to_bp) not in merge_heads:
                    out.append({'arm': side, 'run': r.id, 'field': 'segment',
                                'kind': 'invariant',
                                'msg': 'm run anchor is no parent of a merge at its to_bp'})
        for leaf in leaves(arm):
            for e in end_labels(arm, leaf):
                route_from, _ = evidence(arm, arm.runs[e.run])
                if route_from != e.route_bp:
                    out.append({'arm': side, 'segment': leaf, 'run': e.run,
                                'field': 'route_bp', 'kind': 'invariant',
                                'msg': 'merge-derived evidence %d != T.route_bp %d'
                                       % (route_from, e.route_bp)})
        if ref and 'runs' in ref:
            for r, rr in zip(arm.runs, ref['runs']):
                if 'segment' in rr and rr['segment'] != r.segment:
                    out.append({'arm': side, 'run': r.id, 'field': 'segment',
                                'kind': 'star_mismatch', 'decoded': r.segment,
                                'reference': rr['segment']})
    return out


# (v5.2) each outcome dimension is complete only when no limitation of its class applies
WALK_CLASS = ('walk_domain', 'seed_labels')
LABEL_LOWER_CLASS = ('label_lists', 'inexact_counts', 'seed_labels', 'switch_sources',
                     'greedy_losses')
LABEL_QUALIFIED_CLASS = ('trace_record_boundaries',)


def _outcome_findings(g):
    """The owner's conservative rule (§14, v5.2): a reader can never find a limitation
    whose dimension still reads complete -- and, since each other value is defined by
    the limitations that apply, none that is not complete without one (a K record lost
    from the body shows here)."""
    out = []
    o = g.outcome

    def bad(field, msg):
        out.append({'field': 'outcome.' + field, 'kind': 'invariant', 'msg': msg})

    kinds = {l.kind for l in g.limitations}
    walks = (kinds & set(WALK_CLASS)) | {
        'scope' for l in g.limitations
        if l.kind == 'scope' and isinstance(l.observed, int) and l.observed > 0}
    if walks and o.walks == 'complete':
        bad('walks', 'complete despite %s' % sorted(walks))
    if not walks and o.walks == 'partial':
        bad('walks', 'partial without a walk_domain, seed_labels or merged scope limitation')
    cut = 'branch_events' in kinds
    if cut != (o.branch_diagnostics == 'cut'):
        bad('branch_diagnostics', '%s with%s a branch_events limitation'
            % (o.branch_diagnostics, '' if cut else 'out'))
    qualified = bool(kinds & set(LABEL_QUALIFIED_CLASS))
    lower = bool(kinds & set(LABEL_LOWER_CLASS))
    want = 'qualified' if qualified else 'lower_bound' if lower else 'complete'
    if o.label_evidence != want:
        bad('label_evidence', '%s where the limitations say %s' % (o.label_evidence, want))
    return out


def is_resource_code(code):
    return code in RESOURCE_CODES


def reason_name(code):
    return REASON[code]
