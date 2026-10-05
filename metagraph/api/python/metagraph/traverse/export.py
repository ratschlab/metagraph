"""Exports of a graphlet: to_json() (today's detail: full results[i], the conformance
oracle of T37), FASTA and GFA.

to_json() rebuilds every field from the body. The envelope supplies only what the
body does not carry: seed_id, seed_id_mismatch, the derived-set fields of the seed
(labels_from_seed, labels_supporting_total, labels_dropped, labels_dropped_digest),
the annotation counters, timing and duplicate. Additive fields of the design (§4, §8):
runs[].segment/branches/loss, segments[].labels_via_parent (merges), hairpin `followed`,
revisit `length_bp`/`same_distance`, label_dict[].column/seq_id.

Resources: without a budget (the default) the exports run with no work or allocation
budget. Their output is as large as the graphlet makes it (to_fasta() spells every chosen
walk's whole chain: quadratic in the walk on a comb-shaped trie); an MCP tool's max_bytes
bounds the bytes RETURNED, not this computation or its peak allocation. With budget=
(stage L, budget.py) an export either completes or raises LocalBudgetExceeded: never a
partial text (L4), and save() writes no file; the stop says how many records were done.
"""

import json

from . import budget as _B
from . import derive
from ._codec import REASON, UNLIMITED
from .budget import (DICT_KEY, FLOAT, JSON_ID, LIST, LIST_ITEM, STR, W_GFA_LINE, W_NODE,
                     W_RECORD, W_STEP, LocalBudgetExceeded, dict_bytes, list_bytes)
from .model import ARM_SIDES

_LEVERS = ('select_walks', 'narrow_arm', 'export_mgt')
_BLOCK = 64


def _blocks(bud, pw, pm, done, n=None):
    """The (start, end) blocks of _BLOCK elements of an export, each charged at the sum
    of its elements' prices before it is built (an export has no partial answer: one
    charge per block, a deadline noticed between blocks); |done| counts the records."""
    n = len(pw) if n is None else n
    for i in range(0, n, _BLOCK):
        j = min(n, i + _BLOCK)
        if bud is not None:
            bud.charge(sum(pw[i:j]), sum(pm[i:j]))
            done[0] += j - i
        yield i, j

__all__ = ['to_json', 'normalize_result', 'limitation_json', 'to_fasta', 'to_gfa',
           'value_json']


def value_json(v):
    """A typed K value as JSON: the knob's own type, "unlimited" for u."""
    return 'unlimited' if v is UNLIMITED else v


def limitation_json(lim):
    j = {'kind': lim.kind, 'knob': lim.knob, 'limit': value_json(lim.limit),
         'observed': value_json(lim.observed), 'effect': lim.effect}
    if lim.complete_to_bp is not None:
        j['complete_to_bp'] = lim.complete_to_bp
    for name, v in lim.extra:
        j[name] = value_json(v)
    return j


def _event_json(ev):
    j = {'at_bp': ev.at_bp, 'type': ev.type}
    t = ev.type
    if t == 'switch':
        j.update({'from': ev.from_label, 'to': ev.to_label, 'cost': ev.cost})
    elif t in ('blocked', 'hairpin'):
        j['char'] = ev.char
        if t == 'blocked':
            j['reason'] = REASON[ev.reason]
        j['labels'] = list(ev.labels)
        j['labels_total'] = ev.total
        j['truncated'] = ev.total > len(ev.labels)
        if t == 'hairpin':
            j['followed'] = bool(ev.followed)
    elif t == 'revisit':
        j['segments'] = [ev.segment]
        j['length_bp'] = 0 if ev.delta is None else ev.delta
        j['same_distance'] = ev.delta is None
    elif t == 'tip':
        j.update({'char': ev.char, 'length_bp': ev.length_bp})
    elif t == 'bubble':
        j.update({'length_bp': ev.length_bp, 'alleles': ev.alleles})
    return j


def _segment_price(g, arm, seg, sequences, n_ends):
    """(lwu, bytes) of one segment's JSON node: its label lists copied as ints, its events
    and P runs as nodes, its bases (reversed on the left arm)."""
    ids = len(seg.entry) + len(seg.end)
    for p in seg.presence:
        ids += len(p.labels)
    for p in seg.partition:
        ids += len(p)
    for ev in seg.events:
        if ev.labels is not None:
            ids += len(ev.labels)
    nodes = 1 + len(seg.events) + len(seg.presence) + n_ends + (len(seg.parents) > 1)
    w = W_NODE * nodes + (ids >> 1) + (seg.length_bp >> 9)
    m = dict_bytes(14) + nodes * dict_bytes(7) + ids * JSON_ID + 4 * LIST \
        + LIST_ITEM * (nodes + len(seg.partition)) + len(seg.parents) * LIST_ITEM
    if sequences:
        m += 2 * (STR + seg.length_bp)
    return w, m


def _segment_json(g, arm, seg, sequences):
    j = {'id': seg.id, 'parents': list(seg.parents), 'from_bp': seg.from_bp,
         'length_bp': seg.length_bp}
    if sequences:
        j['sequence'] = seg.walk if arm.side == 'right' else seg.walk[::-1]
    j['labels'] = list(seg.entry)
    j['labels_total'] = seg.entry_total
    j['labels_truncated'] = seg.entry_total > len(seg.entry)
    j['labels_at_end'] = list(seg.end)
    events = []
    if len(seg.parents) > 1:
        events.append({'at_bp': seg.from_bp, 'type': 'reconverge',
                       'segments': list(seg.parents)})
    events.extend(_event_json(ev) for ev in seg.events)
    for r in derive.label_end_events(arm).get(seg.id, ()):
        run = arm.runs[r]
        ev = {'at_bp': run.to_bp, 'type': 'label_end', 'label': run.label,
              'reason': run.event_reason}
        if run.code == 'B':
            ev['needed_budget'] = run.needed_budget
        ev['structural_successors'] = run.structural_successors
        events.append(ev)
    # the walker stable-sorts by at_bp only: the order among equal positions is not
    # information (§2.5)
    events.sort(key=lambda e: e['at_bp'])
    j['events'] = events
    if g.mode == 'annotate':
        j['label_sets'] = [{'from_bp': p.from_bp, 'to_bp': p.to_bp, 'labels': list(p.labels),
                            'labels_total': p.total, 'truncated': p.total > len(p.labels)}
                           for p in seg.presence]
    if len(seg.parents) > 1:
        j['labels_via_parent'] = [list(p) for p in seg.partition]
    return j


def _arm_json(g, arm, bud=None, done=None):
    sequences = bool(arm.segments) and arm.segments[0].walk is not None
    if bud is not None:
        derive.uses(bud, g, arm, 'leaves', 'paths', 'splits', 'end_labels', 'label_end_events')
        bud.phase = 'records'
        ends = derive.label_end_events(arm)
    j = {'status': arm.status,
         'frontier_remaining': {'live_paths': arm.frontier.live_paths,
                                'live_labels': arm.frontier.live_labels,
                                'exact': arm.frontier.exact},
         'complete_to_bp': arm.complete_to_bp, 'completeness_scope': arm.scope,
         'labels_per_node': {'cap': g.cap, 'max_seen': arm.labels_per_node.max_seen,
                             'nodes_truncated': arm.labels_per_node.nodes_truncated},
         'evidence': {'complete': arm.evidence_complete_to_bp is None,
                      'complete_to_bp': arm.evidence_complete_to_bp},
         'limitations': [limitation_json(l) for l in g.limitations if l.arm == arm.side]}
    c = arm.cap_trigger
    if c is not None:
        j['cap_trigger'] = {'reason': REASON[c.reason], 'at_bp': c.at_bp,
                            'segment': c.segment, 'live_paths': c.live_paths,
                            'live_labels': c.live_labels, 'exact': c.exact}
    if bud is None:
        j['segments'] = [_segment_json(g, arm, s, sequences) for s in arm.segments]
    else:
        segs = arm.segments
        pw, pm = derive.price_list(bud, g, arm, 'json_segment_prices', lambda: tuple(
            map(list, zip(*[_segment_price(g, arm, x, sequences, len(ends.get(x.id, ())))
                            for x in segs]))) if segs else ([], []))
        j['segments'] = []
        for i0, i1 in _blocks(bud, pw, pm, done):
            j['segments'].extend(_segment_json(g, arm, x, sequences) for x in segs[i0:i1])
    sps = derive.splits(arm, g.mode)
    if bud is not None:
        def split_prices():
            sw, sm = [], []
            for sp in sps:
                n = sum(min(len(arm.segments[c].entry), g.cap) for c in sp.children)
                sw.append(W_NODE * (1 + len(sp.children)) + (n >> 1))
                sm.append(dict_bytes(8) + len(sp.children) * (dict_bytes(6) + LIST
                                                              + LIST_ITEM)
                          + n * JSON_ID + LIST)
            return sw, sm
        sw, sm = derive.price_list(bud, g, arm, 'json_split_prices', split_prices)
        for _ in _blocks(bud, sw, sm, done):
            pass
    splits = []
    for sp in sps:
        branches = []
        for b in derive.split_branches(arm, sp, g.cap):
            b['labels_truncated'] = b['labels_distinct'] > len(b['labels'])
            branches.append(b)
        splits.append({'at_bp': sp.at_bp, 'prefix_bp': sp.at_bp, 'segment': sp.segment,
                       'children': list(sp.children),
                       'kind': 'ambiguous' if sp.ambiguous else 'divergence',
                       'labels_before': sp.labels_before, 'branches': branches})
    j['splits'] = splits
    paths = []
    ps = derive.paths(arm)
    # every walk's chain, and the bases of the walks whose continuation is spelled from
    # them, with their common prefixes built once
    cont_leaves = [p.leaf for p in ps
                   if arm.segments[p.leaf].leaf.continuation is not None
                   and arm.segments[p.leaf].leaf.continuation.sequence is None
                   and arm.segments[p.leaf].leaf.continuation.n > 0] if sequences else []
    if bud is not None:
        steps = sum(arm.segments[p.leaf].depth + 1 for p in ps)
        bp = sum(arm.segments[x].end_bp for x in cont_leaves)
        bud.charge(2 * W_STEP * steps + (bp >> 8),
                 2 * LIST_ITEM * steps + len(ps) * (LIST + 64) + 2 * bp)
    chains = derive.walk_batch(arm, [p.leaf for p in ps], spell=False)
    spelled = {}
    if cont_leaves:
        spelled = derive.walk_batch(arm, cont_leaves, chains=False)
    if bud is not None:
        def path_prices():
            rbs = derive.runs_by_segment(arm)
            pw, pm = [], []
            for p in ps:
                n = len(rbs[p.leaf])
                c = arm.segments[p.leaf].leaf.continuation
                k = 0 if c is None else len(c.labels)
                pw.append(W_NODE * (2 + n) + k)
                pm.append(dict_bytes(9) + n * (dict_bytes(5) + 2 * FLOAT) + k * JSON_ID
                          + (0 if c is None else dict_bytes(4) + STR + c.n + LIST))
            return pw, pm
        pw, pm = derive.price_list(bud, g, arm, 'json_path_prices', path_prices)
        for _ in _blocks(bud, pw, pm, done):
            pass
    for p in ps:
        seg = arm.segments[p.leaf]
        ends = derive.end_labels(arm, p.leaf)
        reasons = {}
        for e in ends:
            name = arm.runs[e.run].reason
            reasons[name] = reasons.get(name, 0) + 1
        pj = {'id': p.id, 'segments': chains[p.leaf][0], 'length_bp': p.length_bp,
              'end_reasons': reasons}
        if seg.leaf.path_reason is not None:
            pj['path_reason'] = REASON[seg.leaf.path_reason]
        pj['n_labels'] = len(ends)
        pj['end_labels'] = [{'label': e.label, 'loss': e.loss, 'branches': e.branches,
                             'run': e.run, 'route_bp': e.route_bp} for e in ends]
        cont = seg.leaf.continuation
        if cont is not None:
            got = spelled.get(p.leaf)
            pj['continuation'] = {
                'sequence': derive.continuation_sequence(g, arm, p.leaf,
                                                         None if got is None else got[1]),
                'labels': list(cont.labels), 'loss_used': cont.loss_used,
                'branches_used': cont.branches_used}
        paths.append(pj)
    j['paths'] = paths
    runs = []
    if bud is not None:
        n = len(arm.runs)
        bud.charge(W_NODE * n, n * (dict_bytes(10) + 2 * FLOAT + 2 * LIST_ITEM))
        done[0] += n
        nb = derive.price_list(bud, g, arm, 'branch_event_ids', lambda: [sum(
            len(be.ambiguous) + len(be.dropped) + sum(len(r.labels) for r in be.refused)
            for be in arm.branch_events)])[0]
        bud.charge(W_NODE * (len(arm.growth) + len(arm.branch_events)) + (nb >> 1),
                 len(arm.growth) * dict_bytes(14) + len(arm.branch_events) * dict_bytes(8)
                 + nb * JSON_ID)
    for r in arm.runs:
        rj = {'label': r.label, 'from_bp': r.from_bp, 'to_bp': r.to_bp,
              'entered_by': r.entered_by, 'route_bp': r.route_bp}
        if r.from_label is not None:
            rj['from'] = r.from_label
            rj['cost'] = r.cost
        if not r.merged:
            rj['end_reason'] = r.reason
        # additive (§4): the anchor, and the lineage's terminal branches and loss
        rj['segment'] = r.segment
        rj['branches'] = r.branches
        rj['loss'] = r.loss
        runs.append(rj)
    j['runs'] = runs
    j['growth'] = [{
        'from_bp': b.from_bp, 'max_live_paths': b.max_live_paths,
        'distinct_live_labels': b.distinct_live_labels, 'live_pairs': b.live_pairs,
        'exact': b.exact, 'steps': b.steps, 'divergences': b.divergences,
        'ambiguous_branches': b.ambiguous_branches, 'splits': b.splits,
        'reconvergences': b.reconvergences, 'bubbles': b.bubbles, 'tips': b.tips,
        'blocked_repeat': b.blocked_repeat,
        'label_ends': {REASON[c]: n for c, n in b.label_ends.items() if n}}
        for b in arm.growth]
    j['branch_events'] = [{
        'at_bp': be.at_bp, 'segment': be.segment, 'successors': be.chars,
        'labels_per_successor': list(be.labels_per_successor),
        'ambiguous': list(be.ambiguous), 'dropped': list(be.dropped),
        'refused': [{'char': r.char, 'cause': r.cause, 'labels': list(r.labels)}
                    for r in be.refused]} for be in arm.branch_events]
    j['branch_events_truncated'] = max(0, arm.branch_events_total - len(arm.branch_events))
    j['needed_budgets'] = derive.needed_budgets(arm)
    j['counters'] = dict(arm.counters)
    return j


def to_json(g, *, budget=None):
    """Today's detail: full results[i] (natural orientation) -- the conformance oracle:
    from_response(r, out).to_json() equals the full result after normalize_result().
    budget= (stage L): every node is charged before it is built; a stop raises
    LocalBudgetExceeded and no result is returned (done: the nodes built)."""
    bud = _B.resolve(budget)
    if bud is None:
        return _to_json(g, None)
    with bud.scope('to_json', _LEVERS):
        done = [0]
        try:
            return _to_json(g, bud, done)
        except LocalBudgetExceeded as e:
            raise e.restate(records=done[0])


def _to_json(g, bud, done=None):
    g.require_envelope('to_json()')
    meta = g.seed_summary.get('seed', {})
    nsl = g.seed.num_seed_labels
    seed = {
        'seed_id': meta.get('seed_id', ''),
        'validated_seed_id': g.seed.validated_seed_id,
        'seed_id_mismatch': meta.get('seed_id_mismatch', False),
        'length_bp': g.seed.length_bp,
        'num_kmers': g.seed.num_kmers,
        'labels': [l.name for l in g.labels[:nsl]],
        'labels_from_seed': meta.get('labels_from_seed', False),
        'labels_supporting_total': meta.get('labels_supporting_total', nsl),
        'labels_dropped': meta.get('labels_dropped', 0),
        'labels_dropped_digest': meta.get('labels_dropped_digest', ''),
        'dropped_labels': [{'label': d.name, 'reason': d.reason,
                            'runs': [[a, b] for a, b in d.runs], 'runs_kind': 'presence'}
                           for d in g.dropped],
    }
    out = {'seed': seed,
           'limitations': [limitation_json(l) for l in g.limitations if l.arm is None],
           'outcome': g.outcome.as_dict(), 'label_mode': g.mode}
    dict_ = []
    for l in g.labels:
        lj = {'name': l.name, 'kind': l.kind, 'column': l.column}
        if l.kind == 'header':
            lj['seq_id'] = l.seq_id
        dict_.append(lj)
    out['label_dict'] = dict_
    out['arms'] = {side: _arm_json(g, g.arms[side], bud, done) for side in ARM_SIDES
                   if side in g.arms}
    if bud is not None:
        from .ops import _uses_g
        _uses_g(bud, g, 'label_summary')
        n = len(g.labels)
        runs = sum(len(a.runs) for a in g.arms.values())
        bud.charge(W_NODE * n * (1 + len(g.arms)) + runs,
                 list_bytes(n) + n * (dict_bytes(3) + len(g.arms) * (dict_bytes(4) + LIST))
                 + runs * JSON_ID)
    summary = derive.label_summary(g)
    out['label_summary'] = [dict({'label': i}, **{side: dict(per[side], runs=list(per[side]['runs']))
                                                   for side in ARM_SIDES if side in per})
                            for i, per in enumerate(summary)]
    if g.resource_stop is not None:
        q = g.resource_stop
        out['resource_stop'] = {
            'scope': q.scope, 'resource': q.resource, 'phase': q.phase,
            'requested': value_json(q.requested), 'effective': value_json(q.effective),
            'used': value_json(q.used), 'remaining': value_json(q.remaining),
            'actions': list(q.actions), 'message': q.message}
    out['annotation'] = g.seed_summary.get('annotation', {})
    if 'timing' in g.seed_summary:
        out['timing'] = g.seed_summary['timing']
    if g.seed_summary.get('duplicate'):
        out['duplicate'] = True
    return out


def normalize_result(result):
    """The §2.5 normalisation for comparing a to_json() result with the server's
    detail: full result: events sorted within equal at_bp, needed_budgets sorted,
    timing removed, growth label_ends null -> {}. Returns a normalised deep copy."""
    r = json.loads(json.dumps(result))
    r.pop('timing', None)
    for arm in (r.get('arms') or {}).values():
        for s in arm.get('segments', []):
            s['events'] = sorted(s.get('events', []),
                                 key=lambda e: (e['at_bp'], json.dumps(e, sort_keys=True)))
        arm['needed_budgets'] = sorted(arm.get('needed_budgets', []))
        for b in arm.get('growth', []):
            if b.get('label_ends') is None:
                b['label_ends'] = {}
    return r


# ----------------------------------------------------------------------- FASTA

def _wrap(seq, width):
    if not width:
        return seq
    return '\n'.join(seq[i:i + width] for i in range(0, len(seq), width)) or ''


def _fasta_prices(bud, g, a):
    """(lwu, bytes) per walk (path id) of its FASTA record without the seed: the header
    from its end labels, the bases oriented and the line in the list (cached)."""
    def prices():
        rbs = derive.runs_by_segment(a)
        pw, pm = [], []
        for p in derive.paths(a):
            pw.append(W_RECORD + (p.length_bp >> 8) + len(rbs[p.leaf]))
            pm.append(3 * (STR + p.length_bp + 1) + STR + 120 + 2 * LIST_ITEM)
        return pw, pm
    return derive.price_list(bud, g, a, 'fasta_prices', prices)


def _oriented(g, side, w, orientation, with_seed):
    """ops.spell() of a walk whose walking-order bases |w| are spelled already."""
    seed = g.seed.sequence
    if orientation == 'natural':
        flank = w if side == 'right' else w[::-1]
        if not with_seed:
            return flank
        return seed + flank if side == 'right' else flank + seed
    if not with_seed:
        return w
    return (seed if side == 'right' else seed[::-1]) + w


def to_fasta(g, arm=None, leaves=None, with_seed=True, orientation='natural', width=None,
             *, budget=None):
    """One record per walk: '>{arm}_{path_id} length_bp=.. reason=.. labels=..'.
    budget= (stage L): the spellings and every record are charged before they are
    built; a stop raises LocalBudgetExceeded and no text is returned."""
    bud = _B.resolve(budget)
    if bud is None:
        return _to_fasta(g, arm, leaves, with_seed, orientation, width, None)
    with bud.scope('to_fasta', _LEVERS):
        done = [0]
        try:
            return _to_fasta(g, arm, leaves, with_seed, orientation, width, bud, done)
        except LocalBudgetExceeded as e:
            raise e.restate(records=done[0])


def _to_fasta(g, arm, leaves, with_seed, orientation, width, bud, done=None):
    from . import ops
    out = []
    total = 0
    sides = [g.arm(arm).side] if arm is not None else [s for s in ARM_SIDES if s in g.arms]
    for side in sides:
        a = g.arms[side]
        paths = derive.paths(a)
        chosen = paths if leaves is None else [paths[ops.path_id(a, x)] for x in leaves]
        if not chosen:
            continue
        if orientation not in ('natural', 'walk'):
            raise ValueError("orientation is 'natural' or 'walk'")
        # every chosen walk's bases at once, their common prefixes spelled once
        if not a.segments or a.segments[0].walk is None:
            raise ValueError(derive.NO_BASES)
        if bud is not None:
            derive.uses(bud, g, a, 'leaves', 'paths', 'end_labels')
            if leaves is None:
                # every walk: the spelling's bound is the arm's structural totals
                z = derive.arm_sizes(a)
                bud.charge(len(chosen))
                bud.charge(2 * W_STEP * z['leaf_depth'] + (z['walk_bp'] >> 8),
                           2 * len(chosen) * (STR + 40 + 56) + 2 * z['walk_bp']
                           + 2 * 9 * z['leaf_depth'])
            else:
                ops._charge_spellings(bud, a, {p.leaf for p in chosen})
            bud.phase = 'records'
            seed_n = len(g.seed.sequence) if with_seed else 0
            fw, fm = _fasta_prices(bud, g, a)
            if width:
                pw, pm = [], []
                for p in chosen:
                    lines = (p.length_bp + seed_n) // width + 1
                    pw.append(fw[p.id] + (seed_n >> 8))
                    pm.append(fm[p.id] + 3 * (seed_n + lines))
                    total += p.length_bp + seed_n + lines + 120
            elif leaves is None:
                pw, pm = fw, fm
                total += derive.arm_sizes(a)['walk_bp'] + len(chosen) * (seed_n + 121)
            else:
                pw = [fw[p.id] for p in chosen]
                pm = [fm[p.id] for p in chosen]
                total += sum([p.length_bp for p in chosen]) + len(chosen) * (seed_n + 121)
            if not width:
                # the seed in each record (its bases spelled once more per record)
                bud.charge(len(chosen) * (seed_n >> 8), 3 * len(chosen) * (seed_n + 1))
        spelled = derive.walk_batch(a, [p.leaf for p in chosen], chains=False)
        if bud is not None:
            for _ in _blocks(bud, pw, pm, done, len(chosen)):
                pass
        for p in chosen:
            seq = _oriented(g, side, spelled[p.leaf][1], orientation, with_seed)
            seg = a.segments[p.leaf]
            reason = (REASON[seg.leaf.path_reason] if seg.leaf.path_reason
                      else ','.join(sorted({a.runs[e.run].reason
                                            for e in derive.end_labels(a, p.leaf)})) or '.')
            n = len(derive.end_labels(a, p.leaf)) if g.mode == 'constrain' else len(seg.end)
            out.append('>%s_%d length_bp=%d reason=%s labels=%d seed=%s' % (
                side, p.id, p.length_bp, reason, n, 'yes' if with_seed else 'no'))
            out.append(_wrap(seq, width))
    if bud is not None:
        bud.phase = 'join'
        bud.charge(total >> 8, 2 * (STR + total))
    return '\n'.join(out) + ('\n' if out else '')


# ----------------------------------------------------------------------- GFA

def _seg_name(side, sid):
    return '%s%d' % (side[0], sid)


def _gfa_prices(bud, g, arm, k1):
    """(lwu, bytes per segment, the text) of the S and L lines of each segment, cached."""
    def prices():
        pw, pm, tt = [], [], 0
        for s in arm.segments:
            n = s.length_bp + k1
            nl = max(1, len(s.parents))
            pw.append(W_GFA_LINE * (1 + nl) + min(s.depth + 1, k1 + 1) + (n >> 8))
            pm.append(4 * (STR + n) + 40 + nl * (STR + 40 + LIST_ITEM))
            tt += n + 30 + 30 * nl
        return pw, pm, tt
    return derive.price_list(bud, g, arm, 'gfa_prices', prices)


def _gfa_path_prices(bud, g, arm):
    """(lwu, bytes per walk, the text) of the P line of each walk, cached."""
    def prices():
        segs = arm.segments
        rbs = derive.runs_by_segment(arm)
        pw, pm, tt = [], [], 0
        for p in derive.paths(arm):
            d = segs[p.leaf].depth + 2
            n = len(rbs[p.leaf]) if g.mode == 'constrain' else len(segs[p.leaf].end)
            pw.append(W_GFA_LINE + 3 * d + n)
            pm.append(d * (STR + 10 + LIST_ITEM) + 2 * (STR + 14 * d + 24 * n) + LIST)
            tt += 14 * d + 24 * n + 30
        return pw, pm, tt
    return derive.price_list(bud, g, arm, 'gfa_path_prices', prices)


def _chain_tail(segs, seg_id, n):
    """The last |n| or more bases (walking order) of the first-parent chain root -> seg:
    the whole chain when it is shorter. A suffix of derive.walk_bases(), read back only as
    far as it is needed."""
    parts = []
    have = 0
    s = seg_id
    while True:
        w = segs[s].walk
        if w is None:
            raise ValueError(derive.NO_BASES)
        parts.append(w)
        have += len(w)
        p = segs[s].parents
        if have >= n or not p:
            break
        s = p[0]
    parts.reverse()
    return ''.join(parts)


def to_gfa(g, with_seed=True, *, budget=None):
    """GFA 1.0 on the seed strand (natural orientation). Segments are the trie/DAG
    segments with k-1 bases of context, so every L line overlaps by k-1 (de Bruijn
    style): a right-arm segment carries the k-1 bases before it, a left-arm segment the
    k-1 bases after it. P per walk (seed included when with_seed); LB (end labels as
    refs) and ER (end reason) tags on P lines. Requires the bases. budget= (stage L):
    every line is charged before it is built; a stop raises LocalBudgetExceeded and no
    text is returned."""
    bud = _B.resolve(budget)
    if bud is None:
        return _to_gfa(g, with_seed, None)
    with bud.scope('to_gfa', _LEVERS):
        done = [0]
        try:
            return _to_gfa(g, with_seed, bud, done)
        except LocalBudgetExceeded as e:
            raise e.restate(records=done[0])


def _to_gfa(g, with_seed, bud, done=None):
    k1 = g.k - 1
    total = 0
    seed = g.seed.sequence
    lines = ['H\tVN:Z:1.0']
    if with_seed:
        lines.append('S\tseed\t%s\tLN:i:%d' % (seed, len(seed)))
    links = []
    paths = []
    for side in ARM_SIDES:
        arm = g.arms.get(side)
        if arm is None:
            continue
        segs = arm.segments
        if bud is not None:
            derive.uses(bud, g, arm, 'leaves', 'paths', 'end_labels')
            bud.phase = 'records'
            if segs and segs[0].walk is None:
                raise ValueError('to_gfa() needs the bases (output.sequences: true)')
            pw, pm, tt = _gfa_prices(bud, g, arm, k1)
            total += tt
            for _ in _blocks(bud, pw, pm, done):
                pass
        for s in segs:
            if s.walk is None:
                raise ValueError('to_gfa() needs the bases (output.sequences: true)')
            # the walking-order context: seed (outward) + the bases before s on its
            # first-parent chain -- only their last k-1 are used, so the chain is read
            # back only that far (spelling the whole chain per segment was quadratic in
            # the depth)
            before = _chain_tail(segs, s.parents[0], k1) if s.parents and k1 > 0 else ''
            outward_seed = seed if side == 'right' else seed[::-1]
            ctx = (outward_seed + before)[-k1:] if k1 > 0 else ''
            walk = ctx + s.walk
            text = walk if side == 'right' else walk[::-1]
            lines.append('S\t%s\t%s\tLN:i:%d' % (_seg_name(side, s.id), text, len(text)))
            if s.parents:
                for p in s.parents:
                    a, b = (_seg_name(side, p), _seg_name(side, s.id)) if side == 'right' \
                        else (_seg_name(side, s.id), _seg_name(side, p))
                    links.append('L\t%s\t+\t%s\t+\t%dM' % (a, b, k1))
            elif with_seed:
                a, b = ('seed', _seg_name(side, s.id)) if side == 'right' \
                    else (_seg_name(side, s.id), 'seed')
                links.append('L\t%s\t+\t%s\t+\t%dM' % (a, b, k1))
        ps = derive.paths(arm)
        if bud is not None:
            steps = derive.arm_sizes(arm)['leaf_depth']
            bud.charge(2 * W_STEP * steps, 2 * LIST_ITEM * steps + len(ps) * (LIST + 64))
        chains = derive.walk_batch(arm, [p.leaf for p in ps], spell=False)
        if bud is not None:
            pw, pm, tt = _gfa_path_prices(bud, g, arm)
            total += tt
            for _ in _blocks(bud, pw, pm, done):
                pass
        for p in ps:
            names = [_seg_name(side, x) + '+' for x in chains[p.leaf][0]]
            if side == 'left':
                names.reverse()
            if with_seed:
                names = names + ['seed+'] if side == 'left' else ['seed+'] + names
            leaf = segs[p.leaf]
            if g.mode == 'constrain':
                refs = [g.labels[e.label].ref for e in derive.end_labels(arm, p.leaf)]
            else:
                refs = [g.labels[l].ref for l in leaf.end]
            reason = REASON[leaf.leaf.path_reason] if leaf.leaf.path_reason else 'semantic'
            paths.append('P\t%s_%d\t%s\t%s\tLB:Z:%s\tER:Z:%s' % (
                side, p.id, ','.join(names), ','.join(['%dM' % k1] * (len(names) - 1)) or '*',
                ','.join(refs) or '.', reason))
    if bud is not None:
        bud.phase = 'join'
        bud.charge(total >> 8, 2 * (STR + total) + list_bytes(len(lines) + len(links)
                                                          + len(paths)))
    return '\n'.join(lines + links + paths) + '\n'
