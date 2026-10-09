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
    'entry_base', 'rule_entry', 'rule_end', 'label_changes', 'rule_partition',
    'runs_by_segment',
    'runs_by_label', 'label_runs', 'partition_sets', 'merge_parts',
    'leaves', 'paths', 'path_of_leaf', 'chain', 'splits', 'split_branches',
    'label_end_events', 'reconverge_events', 'end_labels', 'labels_at_end',
    'label_summary', 'needed_budgets', 'continuation_sequence', 'walk_bases',
    'walk_batch', 'walk_iter', 'walk_iter_bytes', 'natural_flank', 'evidence',
    'lineage_label_at', 'merge_above',
    'run_leaves', 'check_rules',
]

NO_BASES = 'this retrieval carries no bases (output.sequences: false)'
_EMPTY_SET = frozenset()


def _union(sets):
    out = set()
    for s in sets:
        out.update(s)
    return sorted(out)


# ------------------------------------------------------------------ the §2.3 rules

def entry_base(segments, seg):
    """The SETEXPR base of G.entry (§2.1): none for the root, the union of the
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
    labels_via_parent is empty: the format states the arity)."""
    if len(parents) <= 1:
        return []
    return [array('I') for _ in parents]


def label_changes(seg, runs_here, runs):
    """Constrain: the label changes inside segment |seg| as (position, 0 = a run end | 1 = a
    switch-in, label) in the chronological order of the G.end rule (ties: ends before
    switch-ins): every run anchored here (|runs_here|, ids into |runs|) ending before the
    segment's last node, and every switch event. A run ending at to_bp is absent from base
    to_bp on; a switch-in is present from at_bp. The one statement of the normative
    tie-break: rule_end() and the displayed support (ops) both read it."""
    end_bp = seg.from_bp + seg.length_bp
    ops = []
    for r in runs_here:
        run = runs[r]
        if run.to_bp < end_bp:
            ops.append((run.to_bp, 0, run.label))
    for ev in seg.events:
        if ev.type == 'switch':
            ops.append((ev.at_bp, 1, ev.to_label))
    ops.sort(key=lambda o: (o[0], o[1]))
    return ops


def rule_end(mode, seg, entry, presence, runs_here, runs):
    """G.end '*' (chronological). Constrain: from the entry set apply, in increasing
    position (ties: ends before switch-ins), every run end anchored here with
    to_bp < from_bp + length_bp (remove its label) and every switch event here (add its
    target). Runs ending at the segment's last node are still alive there
    (walker.cpp:2580). A set formula would lose a label that re-enters after ending
    (A -> B -> A inside one segment). Annotate: the last P set, or the entry set."""
    if mode != 'constrain':
        return list(presence[-1].labels) if presence else list(entry)
    ops = label_changes(seg, runs_here, runs)
    if not ops:
        return list(entry)
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


def runs_by_label(arm):
    """label id -> the ids of its runs, in run order. routes() and label_walks() answer
    for one label: scanning every run of the arm per label made a listing over all labels
    quadratic (labels x runs). Callers that answer for one label read it through
    label_runs(), which builds it only once it pays off."""
    got = arm.cache.get('runs_by_label')
    if got is None:
        got = {}
        for r in arm.runs:
            got.setdefault(r.label, []).append(r.id)
        arm.cache['runs_by_label'] = got
        # label_runs()'s count, which nothing reads; on a single-label arm the swap leaves
        # every length equal, which is why ops.cache_signature() reads the keys too
        arm.cache.pop('run_scans', None)
    return got


# Building runs_by_label costs about six scans of the runs (58-70 us against 10-12 us on
# the 900-1,000-label retrievals), so building it on a fresh model's first single-label
# call would make routes() and label_walks() for one label 2-5x slower than a scan. The
# first lookups of an arm therefore scan; the index is
# built by the lookup after them, so a listing over every label pays about one index more
# than building it first, and a caller that asks for a few labels pays no more than the
# scans did (the rule of ops._LABEL_SCANS for label names).
_RUN_SCANS = 6


def label_runs(arm, label):
    """The ids of |label|'s runs on |arm|, in run order (a list; () for none):
    runs_by_label() once it is built, else a scan of the runs for the arm's first
    _RUN_SCANS lookups. The same ids in the same order either way."""
    got = arm.cache.get('runs_by_label')
    if got is None:
        n = arm.cache.get('run_scans', 0)
        if n < _RUN_SCANS:
            arm.cache['run_scans'] = n + 1
            return [r.id for r in arm.runs if r.label == label]
        got = runs_by_label(arm)
    return got.get(label, ())


def _merge_sets(seg, shared):
    """One merge's partition as frozensets; equal arrays (interned: one object) share one
    set through |shared| (id(array) -> set)."""
    sets = []
    for p in seg.partition:
        fs = shared.get(id(p))
        if fs is None:
            fs = shared[id(p)] = frozenset(p) if len(p) else _EMPTY_SET
        sets.append(fs)
    return tuple(sets)


def partition_sets(arm):
    """merge segment id -> its partition (G.partition, one set per parent) as frozensets.
    Membership in the stored array('I') is a linear scan of up to every label id, and the
    evidence walk and the route reconstruction test one label per merge they pass; equal
    arrays (interned: one object) share one set. Segment.partition keeps the arrays, which
    is what dump() encodes and what the model exposes. Built whole here (what a budgeted
    call is charged for: derive.uses() builds it with its price); unbudgeted membership
    tests read merge_parts(), which builds a merge's sets only once they pay off. The
    sets merge_parts() built are reused, and its state is dropped (one copy)."""
    got = arm.cache.get('partition_sets')
    if got is None:
        built = arm.cache.pop('merge_parts', None) or {}
        shared = arm.cache.pop('partition_shared', None) or {}
        arm.cache.pop('merge_scanned', None)
        got = {}
        for s in arm.segments:
            if len(s.parents) > 1:
                x = built.get(s.id)
                got[s.id] = x if x is not None else _merge_sets(s, shared)
        arm.cache['partition_sets'] = got
    return got


# Membership in a merge's stored arrays is a scan that stops at the label; building the
# merge's sets costs about 1.3x a scan of all its arrays that finds nothing (array('I')
# boxes each element it compares: about 16 us per 1,000 ids). One label's routes() and
# label_walks() test a merge on its runs' routes at most twice (the route, then the
# evidence walk: measured over every label of the 900-1,000-label retrievals), and
# building at the second test would make label_walks() of one label 1.7x slower than
# scanning. A merge is scanned for its first _MERGE_SCANS tests and its sets are
# built at the next: a loop over many labels or runs pays at most two scans more per
# merge than building first.
_MERGE_SCANS = 2


def merge_parts(arm, seg_id):
    """The partition of merge |seg_id| (one container per parent, in parents order) for a
    membership test: its frozensets once partition_sets() is built or once this merge has
    been tested _MERGE_SCANS times; until then, the stored arrays. The same answer to
    `label in part` either way.

    Unbudgeted state, three flat dicts in Arm.cache -- merge_scanned (merge -> the tests it
    was scanned for), merge_parts (merge -> its sets) and partition_shared (id(array) ->
    set) -- at most the sets partition_sets() holds plus one entry per merge and per
    distinct array. Flat, so that every step that adds memory adds an entry: the store
    re-measures a model's caches only when ops.cache_signature() (keys and lengths) changes,
    and state kept in a fixed-size tuple would grow unseen by its max_ram_mb account."""
    c = arm.cache
    got = c.get('partition_sets')
    if got is not None:
        return got[seg_id]
    # the scans touch one dict only (a call on a fresh model costs no more than a scan)
    scanned = c.get('merge_scanned')
    if scanned is None:
        c['merge_scanned'] = {seg_id: 1}
        return arm.segments[seg_id].partition
    n = scanned.get(seg_id, 0)
    if n < _MERGE_SCANS:
        scanned[seg_id] = n + 1
        return arm.segments[seg_id].partition
    built = c.get('merge_parts')
    if built is None:
        built = c['merge_parts'] = {}
    x = built.get(seg_id)
    if x is None:
        shared = c.get('partition_shared')
        if shared is None:
            shared = c['partition_shared'] = {}
        x = built[seg_id] = _merge_sets(arm.segments[seg_id], shared)
    return x


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
    walker's processing order. kind ambiguous <=> the parent's stored G.split."""
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
    """The labels alive at a leaf, ascending: the runs anchored at the leaf with
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
    s = leaf
    while True:
        seg = segs[s]
        w = seg.walk
        if w is None:
            raise ValueError(NO_BASES)
        parts.append(w)
        p = seg.parents
        if not p:
            break
        s = p[0]
    parts.reverse()
    return ''.join(parts)


def walk_batch(arm, targets, spell=True, chains=True):
    """{seg: (chain, bases)} for every segment of |targets|: its first-parent chain root
    -> seg (a new list per target; None with chains=False, for a caller that reads only
    the bases) and, with |spell|, its bases root -> end in WALKING order (None when the
    retrieval carries no bases; a ValueError for a model with bases on some segments
    only, as walk_bases()).

    The targets share the prefixes of their common ancestors: a prefix is kept only where
    two targets' chains part (or one target goes on below another), so the work and the
    memory follow the output. Spelling each target's chain on its own re-reads every
    shared segment, which made spelling all walks quadratic in their shared depth; a
    prefix kept at every segment of a deep chain would be quadratic in memory instead.
    For one target this is chain() and walk_bases(). The dict holds every target's
    spelling at once: a caller that uses each spelling once and drops it (an export, a
    comparison's keys) iterates walk_iter() instead."""
    # built in full before it is returned: a missing base never leaves part of the dict, so
    # the stream need not read the chains for it first
    return {t: (c, w) for t, c, w in walk_iter(arm, targets, spell, chains, check_first=False)}


def walk_iter(arm, targets, spell=True, chains=True, *, check_first=True):
    """walk_batch() as a stream: an iterator of (seg, chain, bases) per target, in
    increasing segment id (walk_batch()'s order), each built when it is reached (a few
    targets are built at once, see below). Only the prefixes a target
    still to come is spelled from are held: a kept prefix is dropped once the last of the
    chains parting at its segment has been spelled from it, so on a comb (the walker's
    segment ids put each spine segment just before its two children) a few prefixes are
    held at a time, not one per spine segment -- the batch's prefixes and the dict of every
    spelling would make to_fasta() peak at 10x its text on a comb of 1-bp segments, against
    3.1x when each walk is spelled on its own. A consumer that keeps no spelling holds
    O(depth) of them.

    A model with bases on some segments only raises its ValueError before the first
    target is yielded (for a few targets: from this call), as walk_batch() does, so a
    caller that answers that error for the whole call never sees part of the targets.
    Target by target that costs a read of the
    targets' chains first (the batch finds it while it collects their union): a consumer
    that drops what it was given when the error comes passes |check_first|=False and gets
    the error at the first target without bases instead, at no extra cost (every caller in
    this library does: the exports, compare()'s keys, walk_batch())."""
    segs = arm.segments
    if isinstance(targets, (list, tuple, set, frozenset)) and len(targets) <= _FEW_WALKS:
        # a few walks (a small export, walks(top=2)): chain() and walk_bases() in one loop
        # each, without the bound and the batch's set and dict, which would make small
        # exports 25-60 % slower than walk_bases(). Target by target costs at most k times
        # the deepest chain, the batch at least three Python steps per segment of it, so for
        # k <= 3 targets it is never the slower choice. Built here and handed out as an
        # iterator over them, not by a generator: its frame would cost a one-walk to_fasta()
        # a tenth of its time (0.4 of 4 us). A missing base raises before any target is given
        if not targets:
            return iter(())
        has_bases = spell and segs[0].walk is not None
        if len(targets) == 1:
            for t in targets:
                return iter(((t,) + _walk_one(segs, t, has_bases, chains),))
        return iter([(t,) + _walk_one(segs, t, has_bases, chains)
                     for t in sorted(set(targets))])
    return _walk_stream(segs, targets, spell, chains, check_first)


def _walk_stream(segs, targets, spell, chains, check_first):
    """walk_iter() of more than a few targets, as a generator."""
    tset = set(targets)
    if not tset:
        return
    has_bases = spell and bool(segs) and segs[0].walk is not None
    # target by target while that costs less -- each target's chain walked on its own --,
    # else as a batch, where finding the shared prefixes pays off. Finding them (a set and
    # a dict over the chains' union, a sort) would make a few short walks 2-4x slower to
    # spell than walk_bases() (to_fasta(), to_gfa() and walks(top=5) on the small and the
    # 1,000-label fixtures). The total is read from the segments' depths before any
    # walk (no work is spent on a walk that is then redone); a depth the model does not
    # keep (a hand-built model's 0) only changes the choice, never the answer: with chains
    # the walk itself is counted by the chains' lengths (without chains by the depths),
    # and past _EACH_STEPS times the arm's segments the targets not yet given are spelled
    # as a batch. Written out here rather than in a generator of its own: a generator
    # level per target would cost what a small export's own work does (to_gfa() of a
    # 6-walk retrieval 10 % slower)
    budget = _EACH_STEPS * len(segs)
    if sum([segs[t].depth for t in tset]) + len(tset) > budget:
        yield from _batch_stream(segs, tset, has_bases, chains)
        return
    if has_bases and check_first:
        # bases on some segments only (a hand-built model; a parsed one has them on all or
        # none) raise before the first target: the targets' chains are read once for it,
        # never the whole arm (reading every segment would cost a few shallow walks of an
        # 8,000-segment arm 50x their own time)
        _check_bases(segs, tset)
    order = sorted(tset)
    for i, t in enumerate(order):
        c, w = _walk_one(segs, t, has_bases, chains)
        yield t, c, w
        budget -= len(c) if chains else segs[t].depth + 1
        if budget < 0:
            yield from _batch_stream(segs, set(order[i + 1:]), has_bases, chains)
            return


def walk_iter_bytes(n_targets, bp, union, big=0, chain_steps=0):
    """The local-limits memory account of walk_iter() (a bound of what it holds at once) for
    |n_targets| targets with |bp| bases in all, |big| the longest, over a union of |union|
    chain segments: the union's set and children counts; the kept prefixes -- at most |bp|,
    since the prefixes held at once wait each for its own pending segment, and no two of
    those lie on one chain, so each prefix is at most the bases of distinct targets below
    it -- each in a kept entry (at most two per target: targets and partings); and the
    spelling being given, with the pieces it is joined from. With chains, |chain_steps|
    (the targets' chain lengths) bounds the kept chains the same way. The consumer's own
    use of what it is given (a record, a key) is charged by the consumer."""
    from .budget import DICT_KEY, INT, LIST, LIST_ITEM, SET_ITEM, STR
    return (union * (SET_ITEM + DICT_KEY + INT)
            + 2 * n_targets * (LIST + 3 * LIST_ITEM + STR + DICT_KEY)
            + bp + 2 * (STR + big) + LIST_ITEM * chain_steps)


# Spelling target by target costs about one Python step per chain segment, the batch about
# three per segment of the chains' union plus C-level copies of the shared prefixes:
# measured on the 57 benchmark graphlets, target by target was 0.24-0.62x the batch's time
# up to 5 chain steps per segment and 1.5x at 14 (5.2x at 73, to_fasta() of a 2,213-segment
# beam arm), so the crossover is near 9 steps per segment of the union. The arm's segment
# count stands in for the union (at least as large, so the choice errs toward the batch
# only when the union is small, where both are cheap).
_EACH_STEPS = 8
_FEW_WALKS = 3


def _batch_stream(segs, tset, has_bases, chains=True):
    """The batch of walk_iter(): the union of the targets' chains, a prefix kept where they
    part, each kept prefix dropped after its last dependent (see walk_iter())."""
    if not tset:
        return
    # the union of the targets' chains, and per segment how many of its first-parent
    # children lie in it
    inu = set()
    kids = {}
    for t in tset:
        s = t
        while s not in inu:
            inu.add(s)
            seg = segs[s]
            # a segment without bases raises here, before the first target is yielded:
            # every segment spelled below lies in this union
            if has_bases and seg.walk is None:
                raise ValueError(NO_BASES)
            p = seg.parents
            if not p:
                break
            q = p[0]
            kids[q] = kids.get(q, 0) + 1
            s = q
    # where chains part: seg -> [chain (None without |chains|), bases, the chains still to
    # be spelled from it]. Each first-parent child of a kept segment in the union leads to
    # exactly one later segment spelled from it (the first one taken below: segments that
    # are neither targets nor partings have one child in the union), so its count is its
    # number of children in the union, and it is dropped when that reaches 0
    keep = {}
    for s in sorted(x for x in inu if x in tset or kids.get(x, 0) >= 2):
        ids = []
        x = s
        base = None
        while True:
            ids.append(x)
            p = segs[x].parents
            if not p:
                break
            x = p[0]
            base = keep.get(x)
            if base is not None:
                break
        ids.reverse()
        bases = None
        if has_bases:
            parts = [segs[i].walk for i in ids]
            bases = ''.join(parts) if base is None else base[1] + ''.join(parts)
        # the chain only for a caller that reads it: kept at every parting it is 8 bytes
        # per segment, which on a comb of 1-bp segments is 8x the spelling it comes with
        chain_ = None if not chains else ids if base is None else base[0] + ids
        if base is not None:
            base[2] -= 1
            if not base[2]:
                del keep[x]
        n = kids.get(s, 0)
        if n:
            # a parting, or a target that others go on below: kept for them
            keep[s] = [chain_, bases, n]
            if s in tset:
                # its own list: the kept one is read by the chains below it
                yield s, (list(chain_) if chains else None), bases
        elif s in tset:
            yield s, chain_, bases


def _check_bases(segs, targets):
    """walk_bases()'s ValueError when a segment on the chain of one of |targets| has no
    bases, each segment of their union read once."""
    seen = set()
    for t in targets:
        s = t
        while s not in seen:
            seen.add(s)
            seg = segs[s]
            if seg.walk is None:
                raise ValueError(NO_BASES)
            p = seg.parents
            if not p:
                break
            s = p[0]


def _walk_one(segs, t, has_bases, chains=True):
    """(chain root -> t, its bases in walking order or None): chain() and walk_bases() in
    one loop; without |chains| (None for the chain) only the bases, as walk_bases()."""
    if not chains:
        if not has_bases:
            return None, None
        parts = []
        s = t
        while True:
            seg = segs[s]
            w = seg.walk
            if w is None:
                raise ValueError(NO_BASES)
            parts.append(w)
            p = seg.parents
            if not p:
                break
            s = p[0]
        parts.reverse()
        return None, ''.join(parts)
    ids = []
    parts = []
    s = t
    while True:
        seg = segs[s]
        ids.append(s)
        parts.append(seg.walk)
        p = seg.parents
        if not p:
            break
        s = p[0]
    ids.reverse()
    if not has_bases:
        return ids, None
    if None in parts:
        raise ValueError(NO_BASES)
    parts.reverse()
    return ids, ''.join(parts)


def natural_flank(arm, leaf):
    """Natural orientation (§2.4): right flank = walking order; left flank = its
    reverse (= concat leaf -> root of the per-segment reversals, spell_path)."""
    w = walk_bases(arm, leaf)
    return w if arm.side == 'right' else w[::-1]


def continuation_sequence(g, arm, leaf, walk=None):
    """The C record's spelling (both arms stated): right arm: the LAST n bases of
    seed + natural(right flank); left arm: the FIRST n bases of natural(left flank) +
    seed (n <= |seed| + flank; a continuation may cross into the seed). |walk|: the
    leaf's bases in walking order when the caller has spelled them already."""
    seg = arm.segments[leaf]
    c = seg.leaf.continuation if seg.leaf else None
    if c is None:
        return None
    if c.sequence is not None:
        return c.sequence
    if c.n == 0:
        return ''
    if walk is None:
        flank = natural_flank(arm, leaf)
    else:
        flank = walk if arm.side == 'right' else walk[::-1]
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
    """-> (route_from, evidence_from) of a run (§5.1).

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
            if label_at not in merge_parts(arm, m)[0]:
                route_from = seg.from_bp
                break
        m = above[seg.parents[0]]
    got = (route_from, max(route_from, run.from_bp))
    cache[run.id] = got
    return got


# ------------------------------------------------------------------ local-limits prices

def arm_sizes(arm):
    """The structural sizes of an arm that the derivations' prices are made of (local
    limits, DESIGN §5.4): facts of the model, cached once computed. Callers charge the pass
    at its price first (sizes_price)."""
    got = arm.cache.get('sizes')
    if got is None:
        segs = arm.segments
        entry_ids = end_ids = presence = presence_ids = events = switches = 0
        merges = partition_ids = parts = leaves = leaf_depth = walk_bp = depth_sum = 0
        leaf_runs = 0
        rbs = runs_by_segment(arm)
        for s in segs:
            entry_ids += len(s.entry)
            end_ids += len(s.end)
            if s.presence:
                presence += len(s.presence)
                for p in s.presence:
                    presence_ids += len(p.labels)
            if s.events:
                events += len(s.events)
                for ev in s.events:
                    if ev.type == 'switch':
                        switches += 1
            if len(s.parents) > 1:
                merges += 1
                parts += len(s.partition)
                for p in s.partition:
                    partition_ids += len(p)
            d = s.depth + 1
            depth_sum += d
            if s.leaf is not None:
                leaves += 1
                leaf_depth += d
                walk_bp += s.from_bp + s.length_bp
                leaf_runs += len(rbs[s.id])
        got = {'segments': len(segs), 'runs': len(arm.runs), 'leaves': leaves,
               'merges': merges, 'parts': parts, 'partition_ids': partition_ids,
               'entry_ids': entry_ids, 'end_ids': end_ids, 'presence': presence,
               'presence_ids': presence_ids, 'events': events, 'switches': switches,
               'leaf_depth': leaf_depth, 'depth_sum': depth_sum, 'walk_bp': walk_bp,
               'leaf_runs': leaf_runs, 'splits': len(splits_of_counts(arm))}
        arm.cache['sizes'] = got
    return got


def splits_of_counts(arm):
    """The split parents (without building the Split records): for arm_sizes."""
    return {s.parents[0] for s in arm.segments if len(s.parents) == 1}


def merge_depth(arm):
    """seg -> the merges at or above it on its first-parent chain: how many merges the
    evidence walk of a run anchored there may test (the per-entry price of evidence())."""
    got = arm.cache.get('merge_depth')
    if got is None:
        got = []
        segs = arm.segments
        for s in segs:
            up = got[s.parents[0]] if s.parents else 0
            got.append(up + (1 if len(s.parents) > 1 else 0))
        arm.cache['merge_depth'] = got
    return got


def subtree_bounds(arm):
    """(seg -> an upper bound of the segments a walk over ALL its children visits, seg ->
    the segments of its first-parent subtree): the per-entry prices of run_leaves() and of
    the displayed leaves below an anchor. Children come after their parents, so the
    segments below seg are at most those after it."""
    got = arm.cache.get('subtree_bounds')
    if got is None:
        segs = arm.segments
        n = len(segs)
        every = [1] * n
        first = [1] * n
        for s in reversed(segs):
            tot = 1
            fp = 1
            for c in s.children:
                tot += every[c]
                if segs[c].parents[0] == s.id:
                    fp += first[c]
            every[s.id] = min(tot, n - s.id)
            first[s.id] = fp
        got = (every, first)
        arm.cache['subtree_bounds'] = got
    return got


def route_depth(arm):
    """seg -> the most segments a route root -> seg can have, through any parents (the
    first-parent depth bounds only the displayed chain): the price of a label's own
    route back from its anchor."""
    got = arm.cache.get('route_depth')
    if got is None:
        got = []
        for s in arm.segments:
            got.append(1 + max((got[p] for p in s.parents), default=0))
        arm.cache['route_depth'] = got
    return got


def segment_ops(arm):
    """seg -> the label changes inside it (constrain: the runs anchored there that end
    before its last node, and its switch-ins): the pieces of displayed support it splits
    into are at most one more."""
    got = arm.cache.get('segment_ops')
    if got is None:
        segs = arm.segments
        got = [0] * len(segs)
        for r in arm.runs:
            if r.to_bp < segs[r.segment].from_bp + segs[r.segment].length_bp:
                got[r.segment] += 1
        for s in segs:
            for ev in s.events:
                if ev.type == 'switch':
                    got[s.id] += 1
        arm.cache['segment_ops'] = got
    return got


def cut_cost(arm, mode):
    """seg -> the price (lwu) of a comparison cut of the first-parent chain root -> seg
    (ops._cut_info): per segment its pieces of displayed support, their sets and the
    labels kept supporting, summed along the chain. The count of a cut's work, read in
    one lookup."""
    got = arm.cache.get('cut_cost')
    if got is None:
        segs = arm.segments
        ops_ = segment_ops(arm) if mode == 'constrain' else None
        n0 = len(segs[0].entry) if segs else 0
        got = []
        for s in segs:
            if mode == 'constrain':
                pieces = 1 + min(ops_[s.id], s.length_bp)
                ids = pieces * (len(s.entry) + len(s.events))
            else:
                pieces = len(s.presence)
                ids = 0
                for p in s.presence:
                    ids += len(p.labels)
            own = 2 * (2 + pieces) + (ids >> 5) + ((pieces * n0) >> 2) + (s.length_bp >> 8)
            got.append((got[s.parents[0]] if s.parents else 0) + own)
        arm.cache['cut_cost'] = got
    return got


def _derivation_prices(z, n_labels, mode):
    """name -> (lwu, bytes) of each cached derivation at its structural (cold) price."""
    from .budget import (LIST, LIST_ITEM, SET, SET_ITEM, DICT, DICT_KEY, TUPLE, INT,
                         W_ELEM, list_bytes, dict_bytes)
    S, R, L = z['segments'], z['runs'], z['leaves']
    M, P = z['merges'], z['splits']
    ids_shift = 2
    return {
        'leaves': (S, list_bytes(L)),
        'paths': (S + 6 * L, list_bytes(L) + L * (72 + 2 * INT)),
        'path_of_leaf': (2 * L, dict_bytes(L) + L * INT),
        'splits': (3 * S + 8 * P, list_bytes(P) + P * (80 + LIST + 2 * LIST_ITEM)
                   + dict_bytes(P)),
        'merge_above': (2 * S, list_bytes(S)),
        'merge_depth': (2 * S, list_bytes(S)),
        'route_depth': (3 * S + z['parts'], list_bytes(S)),
        'cut_cost': (8 * S + z['presence'] + R, 2 * list_bytes(S) + S * INT),
        'segment_ops': (2 * S + R + z['events'], list_bytes(S)),
        'subtree_bounds': (4 * S, 2 * list_bytes(S) + 2 * S * INT),
        'first_paths': (3 * S, list_bytes(S)),
        'label_end_events': (2 * R, dict_bytes(S) + list_bytes(R)),
        'end_labels': (3 * L + 4 * z['leaf_runs'],
                       dict_bytes(L) + L * LIST + z['leaf_runs'] * (72 + LIST_ITEM + FLOAT_B)),
        'partition_sets': (2 * M + (z['partition_ids'] >> 5),
                           dict_bytes(M) + M * (TUPLE + 16) + z['parts'] * SET
                           + z['partition_ids'] * SET_ITEM),
        'runs_by_label': (2 * R, dict_bytes(n_labels) + n_labels * LIST + R * LIST_ITEM),
        'first_parent_at': (2 * S, list_bytes(S)),
        # annotate mode: the displayed sets and the union-rule route ends (shared sets)
        'displayed_presence': (W_ELEM * (2 * S + z['presence'])
                               + ((z['presence_ids'] + z['entry_ids']) >> ids_shift),
                               list_bytes(S) + S * (TUPLE + 16) + S * SET
                               + z['entry_ids'] * SET_ITEM),
        'annotate_route_ends': (W_ELEM * (3 * S + z['presence'])
                                + ((z['presence_ids'] + 2 * z['entry_ids']) >> ids_shift),
                                list_bytes(S) + S * SET + z['entry_ids'] * SET_ITEM
                                + dict_bytes(n_labels) + n_labels * LIST
                                + (z['presence'] + S) * (TUPLE + LIST_ITEM + INT)),
    }


FLOAT_B = 24


def sizes_price(arm):
    """The price of arm_sizes() over the segments and runs (admitted before the pass);
    its P runs and events are charged by uses() from the sizes it read."""
    return 3 * len(arm.segments) + len(arm.runs) // 4


def price(arm, name, n_labels, mode):
    """(lwu, bytes) of one derivation of |arm| at its cold price (arm_sizes must have
    been charged: price() reads it)."""
    table = arm.cache.get('prices')
    if table is None:
        table = arm.cache['prices'] = _derivation_prices(arm_sizes(arm), n_labels, mode)
    return table[name]


def price_list(b, g, arm, name, fn, per=1):
    """A per-element price list of |arm| (fn() -> list): a structural derivation, cached;
    its cold price (|per| lwu per element of the arm's segments and runs) charged once per
    call. What block admission reads instead of pricing each row in Python."""
    if b is not None:
        b.uses((id(arm), name), per * (len(arm.segments) + len(arm.runs)),
               8 * (len(arm.segments) + len(arm.runs)) + 56)
    got = arm.cache.get(name)
    if got is None:
        got = arm.cache[name] = fn()
    return got


def uses(b, g, arm, *names):
    """Charge budget |b| the derivations of |arm| a call uses, each at its cold price and
    once per call (whether a cache holds it or not); the sizes they are priced from
    first. No-op without a budget."""
    if b is None:
        return
    key = id(arm)
    if (key, 'sizes') not in b._used:
        b.uses((key, 'sizes'), sizes_price(arm))
        z = arm_sizes(arm)
        # the pass's P runs and events, read from what it counted (the same charge warm
        # or cold): the one part of a price charged after its work, at most a few ns each
        b.charge((z['presence'] + z['events']) >> 1)
    table = arm.cache.get('prices')
    if table is None:
        # the prices are facts of the model (its sizes): computed once
        table = arm.cache['prices'] = _derivation_prices(arm_sizes(arm), len(g.labels),
                                                         g.mode)
    used = b._used
    for n in names:
        if (key, n) not in used:
            w, m = table[n]
            b.uses((key, n), w, m)
            build = _BUILT_WHEN_CHARGED.get(n)
            if build is not None:
                # paid for whole: a budgeted call reads the whole index, as priced, never
                # the scans that label_runs() and merge_parts() answer unbudgeted calls
                # with before an index pays off (their work is not what the price models)
                build(arm)


_BUILT_WHEN_CHARGED = {'runs_by_label': runs_by_label, 'partition_sets': partition_sets}


# ------------------------------------------------------------------ check_rules

# From this feature level the server orders a merge's parents so that the first -- the one
# the displayed walk follows (paths, spellings, continuations, route_bp) -- is the head
# carrying the most labels, ties in arrival order (the majority-parent rule, SPEC §7.1).
# Retrievals of a lower feature level keep arrival order, and the library displays what
# the body stores either way.
DISPLAYED_PARENT_LEVEL = 6


def carried_labels(g, a, s, i):
    """How many labels carry parent |i| of merge segment |s| (arm |a|) into the merge, by
    the server's majority-parent rule (SPEC §7.1, Walker::merge_level): constrain -- the
    labels whose lineages the parent's head brings into the merge node: the ones the merge
    kept from it (its partition entry) and the ones it closed there (an 'm' run ending at
    the parent's last node, a label that came through another parent too); annotate -- the
    fewest labels present at a node of the parent segment, true counts (every head at the
    merge node holds the node's own labels, so the parent's own bases tell the routes apart:
    as many labels as carry it whole). A label lost on the step into the merge node is in
    neither, as it is in no head the server compares."""
    pid = s.parents[i]
    parent = a.segments[pid]
    if g.mode == 'constrain':
        closed = sum(1 for r in runs_by_segment(a)[pid]
                     if a.runs[r].merged and a.runs[r].to_bp == parent.end_bp)
        return len(s.partition[i]) + closed
    return min([parent.entry_total] + [pr.total for pr in parent.presence])


def feature_level(g):
    """The feature level the response envelope's capabilities state (0: none stated, or a
    body without its envelope)."""
    caps = (g.envelope or {}).get('capabilities')
    v = caps.get('feature_level') if isinstance(caps, dict) else None
    return v if isinstance(v, int) and not isinstance(v, bool) else 0


def displayed_parent_rule(g):
    """Which rule chose the parent g's merges display (the majority-parent rule, SPEC §7.1):
    'majority' (the envelope
    states feature level 6 or more), 'arrival' (it states a lower level, or its
    capabilities state none: a server from before the levels were stated -- every server
    of level 5 or more states its level), or None, not known (a body without its
    envelope, or an envelope without capabilities)."""
    level = feature_level(g)
    if level:
        return 'majority' if level >= DISPLAYED_PARENT_LEVEL else 'arrival'
    caps = (g.envelope or {}).get('capabilities')
    return 'arrival' if isinstance(caps, dict) else None


def majority_first(g, side, depth=None):
    """Whether every merge of arm |side| (those starting before |depth|, when given)
    displays a parent carried by the most labels -- no later parent carries more than the
    first (carried_labels(); ties keep the arrival order) -- so that its display is the one
    the majority-parent rule of feature level 6 makes, whichever rule made it: True; False
    when a merge shows a parent that rule would not (a lower level's arrival order); None
    when a merge cannot be judged (constrain mode with recorded lists cut, whose partitions
    are lower bounds, or a merge without its partition). An arm without merges is True."""
    a = g.arms[side]
    out = True
    for s in a.segments:
        if len(s.parents) < 2 or (depth is not None and s.from_bp >= depth):
            continue
        if g.mode == 'constrain' and (a.labels_per_node.nodes_truncated
                                      or len(s.partition) != len(s.parents)):
            out = None
            continue
        first = carried_labels(g, a, s, 0)
        if any(carried_labels(g, a, s, i) > first for i in range(1, len(s.parents))):
            return False
    return out


def check_rules(g, reference=None):
    """The §2.3 rules recomputed and compared, plus the §4/§9 invariants of the body.

    Returns a list of findings, each {arm, segment|run, field, kind, ...}; [] means
    nothing diverged. Kinds:
      'star_mismatch'  a field the source wrote as '*' decodes to a value other than
                       |reference| (a detail: full results[i]) holds -- the writer
                       emitted '*' where its rule does not reproduce the walker;
      'noncanonical'   an explicit field equal to its rule's value (the canonical form
                       writes '*'): valid, but not what a conforming writer emits;
      'invariant'      a structural invariant of §4/§9 that does not hold (and, for a
                       retrieval whose envelope states feature level 6 or more, the
                       majority-parent rule: at
                       every merge the first parent is carried by the most labels -- the
                       displayed walk follows it; carried_labels()).
    """
    out = []
    out.extend(_outcome_findings(g))
    level = feature_level(g)
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
        if level >= DISPLAYED_PARENT_LEVEL and not (
                g.mode == 'constrain' and arm.labels_per_node.nodes_truncated):
            # the majority-parent rule: from feature level 6 the displayed walk passes a
            # merge through the parent carried by the most labels, ties in arrival order
            # (which the body does not record: only a later parent carrying MORE is a
            # violation). In constrain mode not decidable where recorded lists were cut
            # (label_lists): the partitions are then lower bounds, and the server compares
            # whole lineage states; annotate mode compares the true counts the body states
            for s in segs:
                if len(s.parents) > 1 and len(s.partition) == len(s.parents):
                    first = carried_labels(g, arm, s, 0)
                    for i in range(1, len(s.parents)):
                        have = carried_labels(g, arm, s, i)
                        if have > first:
                            out.append({'arm': side, 'segment': s.id, 'field': 'parents',
                                        'kind': 'invariant',
                                        'msg': 'merge: the first (displayed) parent %d carries '
                                               '%d labels, parent %d carries %d (feature level '
                                               '%d: the first is the one carried by the most)'
                                               % (s.parents[0], first, s.parents[i], have,
                                                  level)})
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


# each outcome dimension is complete only when no limitation of its class applies
WALK_CLASS = ('walk_domain', 'seed_labels')
LABEL_LOWER_CLASS = ('label_lists', 'inexact_counts', 'seed_labels', 'switch_sources',
                     'greedy_losses')
LABEL_QUALIFIED_CLASS = ('trace_record_boundaries',)
# a WALKED result's derivation limitation (SPEC §7.0) is in the qualified class too: its
# permitted set was derived from part of the seed, a superset of the labels carrying the
# whole seed, so a label it reports may be one the whole seed excludes. A failed seed's
# derivation limitation says why there is no result and qualifies nothing (the server's
# outcome_of, traverse.cpp)
WALKED_QUALIFIED_CLASS = ('derivation',)


def qualifies(g, lim):
    """Whether limitation |lim| of graphlet |g| puts its label evidence in the qualified
    class (something reported may be overstated): trace_record_boundaries, and a walked
    result's derivation (a permitted set derived from part of the seed)."""
    if lim.kind in LABEL_QUALIFIED_CLASS:
        return True
    return lim.kind in WALKED_QUALIFIED_CLASS and g.outcome is not None \
        and g.outcome.walks != 'failed'


def partial_derivation(g):
    """The seed-level derivation limitation of a walked result (SPEC §7.0: the time budget
    ran out while the permitted set was derived, after observed of the seed's k-mers; the
    walk stopped at the seed), or None."""
    if g.outcome is None or g.outcome.walks == 'failed':
        return None
    return next((l for l in g.limitations if l.kind == 'derivation' and l.arm is None),
                None)


def _outcome_findings(g):
    """The conservative rule (§14): a reader can never find a limitation
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
    qualified = any(qualifies(g, l) for l in g.limitations)
    lower = bool(kinds & set(LABEL_LOWER_CLASS))
    want = 'qualified' if qualified else 'lower_bound' if lower else 'complete'
    if o.label_evidence != want:
        bad('label_evidence', '%s where the limitations say %s' % (o.label_evidence, want))
    return out


def is_resource_code(code):
    return code in RESOURCE_CODES


def reason_name(code):
    return REASON[code]
