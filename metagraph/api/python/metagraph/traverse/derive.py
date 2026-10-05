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
    'runs_by_label', 'label_runs', 'partition_sets', 'merge_parts',
    'leaves', 'paths', 'path_of_leaf', 'chain', 'splits', 'split_branches',
    'label_end_events', 'reconverge_events', 'end_labels', 'labels_at_end',
    'label_summary', 'needed_budgets', 'continuation_sequence', 'walk_bases',
    'walk_batch', 'natural_flank', 'evidence', 'lineage_label_at', 'merge_above',
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
        # label_runs()'s count, no longer read; on a single-label arm the swap leaves every
        # length equal, which is why ops.cache_signature() reads the keys too
        arm.cache.pop('run_scans', None)
    return got


# Building runs_by_label costs about six scans of the runs (58-70 us against 10-12 us on
# the 900-1,000-label retrievals), and built on a fresh model's first single-label call it
# made routes() and label_walks() for one label 2-5x slower than the scan it replaced (the
# review of the efficiency batch, P10). The first lookups of an arm scan; the index is
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
# building at the second test made label_walks() of one label 1.7x slower than the scans
# of ce949da5 (P10). A merge is scanned for its first _MERGE_SCANS tests and its sets are
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
    # the scans touch one dict only (a call on a fresh model costs what ce949da5's did)
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
    For one target this is chain() and walk_bases()."""
    segs = arm.segments
    if isinstance(targets, (list, tuple, set, frozenset)) and len(targets) <= _FEW_WALKS:
        # a few walks (a small export, walks(top=2)): chain() and walk_bases() in one loop
        # each, without the bound and the batch's set and dict, which made small exports
        # 25-60 % slower than walk_bases() had (P10). Target by target costs at most k times
        # the deepest chain, the batch at least three Python steps per segment of it, so for
        # k <= 3 targets it is never the slower choice
        if not targets:
            return {}
        has_bases = spell and segs[0].walk is not None
        if len(targets) == 1:
            for t in targets:
                return {t: _walk_one(segs, t, has_bases, chains)}
        return {t: _walk_one(segs, t, has_bases, chains) for t in sorted(set(targets))}
    tset = set(targets)
    if not tset:
        return {}
    has_bases = spell and bool(segs) and segs[0].walk is not None
    one = _walk_each(segs, tset, has_bases, chains)
    if one is not None:
        return one
    # the union of the targets' chains, and per segment how many of its first-parent
    # children lie in it
    inu = set()
    kids = {}
    for t in tset:
        s = t
        while s not in inu:
            inu.add(s)
            p = segs[s].parents
            if not p:
                break
            q = p[0]
            kids[q] = kids.get(q, 0) + 1
            s = q
    keep = {}          # where chains part: seg -> (chain, bases)
    out = {}
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
            if None in parts:
                raise ValueError(NO_BASES)
            bases = ''.join(parts) if base is None else base[1] + ''.join(parts)
        chain_ = ids if base is None else base[0] + ids
        n = kids.get(s, 0)
        target = s in tset
        if n >= 2 or (target and n):
            keep[s] = (chain_, bases)
            if target:
                out[s] = (list(chain_) if chains else None, bases)
        elif target:
            out[s] = (chain_ if chains else None, bases)
    return out


# Spelling target by target costs about one Python step per chain segment, the batch about
# three per segment of the chains' union plus C-level copies of the shared prefixes:
# measured on the 57 benchmark graphlets, target by target was 0.24-0.62x the batch's time
# up to 5 chain steps per segment and 1.5x at 14 (5.2x at 73, to_fasta() of a 2,213-segment
# beam arm), so the crossover is near 9 steps per segment of the union. The arm's segment
# count stands in for the union (at least as large, so the choice errs toward the batch
# only when the union is small, where both are cheap).
_EACH_STEPS = 8
_FEW_WALKS = 3


def _walk_each(segs, tset, has_bases, chains=True):
    """walk_batch() target by target, while that costs less: each target's chain walked on
    its own, (chain, bases) as walk_batch() gives them, in its order -- or None where the
    chains' total length exceeds _EACH_STEPS times the arm's segments, where finding the
    shared prefixes pays off. Finding them (a set and a dict over the chains' union, a
    sort) made a few short walks 2-4x slower to spell than walk_bases() had (to_fasta(),
    to_gfa() and walks(top=5) on the small and the 1,000-label fixtures, P10). The total
    is read from the segments' depths before any walk (no work is spent on a walk that is
    then redone); a depth the model does not keep (a hand-built model's 0) only changes
    the choice, never the answer, and with chains the walk itself stops at the same bound
    (counted by the chains' lengths; without chains by the depths)."""
    budget = _EACH_STEPS * len(segs)
    if sum([segs[t].depth for t in tset]) + len(tset) > budget:
        return None
    out = {}
    for t in sorted(tset):
        got = out[t] = _walk_one(segs, t, has_bases, chains)
        budget -= len(got[0]) if chains else segs[t].depth + 1
        if budget < 0:
            return None
    return out


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
    """The C record's spelling (v2, both arms stated): right arm: the LAST n bases of
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
            if label_at not in merge_parts(arm, m)[0]:
                route_from = seg.from_bp
                break
        m = above[seg.parents[0]]
    got = (route_from, max(route_from, run.from_bp))
    cache[run.id] = got
    return got


# ------------------------------------------------------------------ stage L prices

def arm_sizes(arm):
    """The structural sizes of an arm that the derivations' prices are made of (stage L,
    DESIGN §21.3): facts of the model, cached once computed. Callers charge the pass at
    its price first (sizes_price)."""
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
    once per call (L1: whether a cache holds it or not); the sizes they are priced from
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
