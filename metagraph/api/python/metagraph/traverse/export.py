"""Exports of a graphlet: to_json() (today's detail: full results[i], the conformance
oracle of T37), FASTA and GFA.

to_json() rebuilds every field from the body. The envelope supplies only what the
body does not carry: seed_id, seed_id_mismatch, the derived-set fields of the seed
(labels_from_seed, labels_supporting_total, labels_dropped, labels_dropped_digest),
the annotation counters, timing and duplicate. Additive fields of the design (§4, §8):
runs[].segment/branches/loss, segments[].labels_via_parent (merges), hairpin `followed`,
revisit `length_bp`/`same_distance`, label_dict[].column/seq_id.

Resources: without a budget (the default) the exports run with no work or allocation
budget. Their output is as large as the graphlet makes it, and on a comb-shaped trie two of
them are quadratic in the walk: to_fasta() spells every chosen walk's whole chain, and
to_json() -- like the server's detail "full", which it reproduces -- lists every path's
whole chain of segment ids (paths[].segments): n splits make n paths of up to n segments
each, so its text grows as n^2 (the search service measured about 3.4 bytes x n^2 of JSON).
An MCP tool's max_bytes bounds the bytes RETURNED, not this computation or its peak
allocation; graphlet_export(format="mgt") is linear (the body as it is). Their peak is
about twice their text: the walks are spelled as a stream (derive.walk_iter()), each
dropped once its record is written, and the final line feed is joined with the records
rather than added to the joined text (a copy of all of it). With budget=
(stage L, budget.py) an export either completes or raises LocalBudgetExceeded: never a
partial text (L4), and save() writes no file; the stop says how many records were done.
"""

import json
from collections.abc import Sized
from operator import lt as _lt

from . import budget as _B
from . import coords as _C
from . import derive
from ._codec import REASON, UNLIMITED
from .budget import (DICT_KEY, FLOAT, ID_SHIFT_PY, JSON_ID, LIST, LIST_ITEM, STR, W_ELEM,
                     W_GFA_LINE, W_NODE, W_PAIR, W_RECORD, W_STEP, LocalBudgetExceeded,
                     dict_bytes, list_bytes)
from .model import ARM_SIDES, MissingEnvelope

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


def _limitation_price(lim):
    """(lwu, bytes) of limitation_json(lim): a node, its dict."""
    return W_NODE, dict_bytes(6 + len(lim.extra)) + LIST_ITEM


def _envelope_price(g):
    """(lwu, bytes) of the parts of to_json() built outside the arms -- the seed block
    (the seed labels' list, a dict per dropped label with a [from, to] list per run), the
    seed-level limitations, the outcome and label_dict -- and per arm side those of its
    limitations and its header dicts (frontier, labels per node, evidence, cap trigger):
    -> (lwu, bytes, {side: (lwu, bytes)}). A structural derivation of the graphlet, cached;
    to_json() charges it at its cold price (one pass over the dropped labels, limitations
    and labels) on every call. O3: these were built before the first charge, so 10,000
    dropped labels of 40 runs each (a traced peak of 34.5 MB) completed under memory_mb=1."""
    got = g.cache.get('json_envelope_price')
    if got is None:
        nsl = g.seed.num_seed_labels
        nd = len(g.dropped)
        runs = sum(len(d.runs) for d in g.dropped)
        n = len(g.labels)
        w = W_NODE * (2 + nd + n) + W_PAIR * runs
        m = dict_bytes(16) + dict_bytes(12) + dict_bytes(4) + list_bytes(nsl) \
            + list_bytes(nd) + nd * (dict_bytes(4) + LIST) \
            + runs * (LIST + 3 * LIST_ITEM) + list_bytes(n) + n * dict_bytes(4)
        arms = {}
        for lim in g.limitations:
            lw, lm = _limitation_price(lim)
            if lim.arm is None:
                w += lw
                m += lm
            else:
                aw, am = arms.get(lim.arm, (0, 0))
                arms[lim.arm] = (aw + lw, am + lm)
        m += LIST
        for side in g.arms:
            aw, am = arms.get(side, (0, 0))
            arms[side] = (aw + 5 * W_NODE, am + LIST + dict_bytes(20) + 4 * dict_bytes(6))
        got = g.cache['json_envelope_price'] = (w, m, arms)
    return got


def _uses_envelope_price(bud, g):
    bud.uses((id(g), 'json_envelope_price'),
             W_ELEM * (len(g.dropped) + len(g.limitations) + len(g.labels) + 1), 64)
    return _envelope_price(g)


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
        # the arm's header dicts and its limitations, before they are built (O3)
        bud.charge(*_uses_envelope_price(bud, g)[2].get(arm.side, (0, 0)))
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
        # the paths' chains (the output's lists) and the chains the stream keeps for them
        # (at most as many slots again), and the continuations' walks streamed
        steps = sum(arm.segments[p.leaf].depth + 1 for p in ps)
        bp = sum(arm.segments[x].end_bp for x in cont_leaves)
        big = max([arm.segments[x].end_bp for x in cont_leaves], default=0)
        bud.charge(2 * W_STEP * steps + (bp >> 8),
                   LIST_ITEM * steps + len(ps) * LIST
                   + derive.walk_iter_bytes(len(ps), 0, len(arm.segments), 0, steps)
                   + derive.walk_iter_bytes(len(cont_leaves), bp, len(arm.segments), big))
    # the continuations' sequences (their last n bases) as the stream spells their walks,
    # each walk dropped once its continuation is cut from it; the chains below are the
    # output's own lists, built as the paths are written
    conts = {leaf: derive.continuation_sequence(g, arm, leaf, w)
             for leaf, _, w in derive.walk_iter(arm, cont_leaves, chains=False,
                                                check_first=False)}
    chains = derive.walk_iter(arm, [p.leaf for p in ps], spell=False)
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
    for p, (_, chain_, _) in zip(ps, chains):
        seg = arm.segments[p.leaf]
        ends = derive.end_labels(arm, p.leaf)
        reasons = {}
        for e in ends:
            name = arm.runs[e.run].reason
            reasons[name] = reasons.get(name, 0) + 1
        pj = {'id': p.id, 'segments': chain_, 'length_bp': p.length_bp,
              'end_reasons': reasons}
        if seg.leaf.path_reason is not None:
            pj['path_reason'] = REASON[seg.leaf.path_reason]
        pj['n_labels'] = len(ends)
        pj['end_labels'] = [{'label': e.label, 'loss': e.loss, 'branches': e.branches,
                             'run': e.run, 'route_bp': e.route_bp} for e in ends]
        cont = seg.leaf.continuation
        if cont is not None:
            pj['continuation'] = {
                'sequence': conts[p.leaf] if p.leaf in conts
                else derive.continuation_sequence(g, arm, p.leaf),
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


def _json_copy(x):
    """A copy of the JSON value |x| (dicts, lists and atoms, as json.loads made them): what
    copy.deepcopy makes of it, without deepcopy's memo -- a dict entry and a kept-alive
    reference per container copied, never charged, which over a coordinates block of
    15,000 occurrences outgrew the copy itself (the traced peak passed the account by a
    third). A JSON value shares and nests nothing, so no memo is needed; a list is copied
    at its exact length (list.copy()), never over-allocated as a comprehension would."""
    if isinstance(x, dict):
        return {k: _json_copy(v) for k, v in x.items()}
    if isinstance(x, list):
        out = x.copy()
        for i, v in enumerate(out):
            if isinstance(v, (dict, list)):
                out[i] = _json_copy(v)
        return out
    return x


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
    if bud is not None:
        # the seed block, the seed-level limitations, the outcome and label_dict, charged
        # before they are built (O3: their account was missing however many labels were
        # dropped -- the server lists every explicit request label it could not use)
        w, m, _ = _uses_envelope_price(bud, g)
        bud.charge(w, m)
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
    if 'coordinates' in g.seed_summary:
        # the record coordinates are JSON beside the body (MGT v1 is frozen, §18): copied
        # verbatim -- the same block in every detail, so to_json() still equals the full
        # result -- as a copy, so that a caller's edit of the result never reaches the
        # graphlet's summary
        raw = g.seed_summary['coordinates']
        if bud is not None and isinstance(raw, dict):
            entries, occ = _C._block_counts(raw)
            bud.charge(W_NODE * entries + _B.W_COORD * occ,
                       entries * dict_bytes(8) + occ * (LIST + 2 * LIST_ITEM + 2 * _B.INT))
        out['coordinates'] = _json_copy(raw)
        if 'coordinates_reason' in g.seed_summary:
            out['coordinates_reason'] = g.seed_summary['coordinates_reason']
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
             *, coordinates=None, budget=None):
    """One record per walk: '>{arm}_{path_id} length_bp=.. reason=.. labels=..'.

    Record coordinates (§18.3): with |coordinates| None (the default) a retrieval that
    carries them adds ' coords=ACC:start-end,...' to each header -- 1-based, closed, on the
    record's forward strand, the seed included when with_seed -- for the header labels
    (accessions) alive at the walk's leaf whose run was entered by the seed and covers the
    whole walk: at most 8 (label, occurrence) pairs, ordered by label id then start, then '
    coords_more=N' for the rest, and ' coords_cut=1' when such a run's list was cut at the
    server's cap. Column labels are left out (their positions are global, C10). A
    retrieval without them gives exactly the records it always did. False: never;
    True: required (MissingEnvelope on a body alone, ValueError when the retrieval has
    none). The JSON block is 0-based half-open; FASTA headers are 1-based closed.

    budget= (stage L): the spellings and every record are charged before they are
    built; a stop raises LocalBudgetExceeded and no text is returned."""
    if coordinates is not None and not isinstance(coordinates, bool):
        # by type: 0 equals False, yet the dispatch is by identity (it wrote coords=)
        raise ValueError('coordinates is None (when present), True or False, not %r'
                         % (coordinates,))
    if coordinates is True:
        if not g.has_envelope:
            raise MissingEnvelope('to_fasta(coordinates=True) needs the response envelope: '
                                  'a parsed body alone carries no record coordinates')
        if not _C.present(g):
            raise ValueError('this retrieval carries no record coordinates (%s)'
                             % _C.reason(g))
    bud = _B.resolve(budget)
    if bud is None:
        return _to_fasta(g, arm, leaves, with_seed, orientation, width, None,
                         coordinates=coordinates)
    with bud.scope('to_fasta', _LEVERS):
        done = [0]
        try:
            return _to_fasta(g, arm, leaves, with_seed, orientation, width, bud, done,
                             coordinates=coordinates)
        except LocalBudgetExceeded as e:
            raise e.restate(records=done[0])


_COORDS_IN_HEADER = 8


def _fasta_rcs(g, side, coordinates):
    """The arm's RunCoordinates a FASTA header states (None: no coords= field)."""
    if coordinates is False or g.mode != 'constrain' or not _C.present(g):
        return None
    return _C.of(g).arms.get(side)


def _coords_header(g, a, p, rcs, with_seed):
    """' coords=ACC:s-e,...' (1-based closed) for walk |p|, '' when nothing qualifies: a
    header label alive at the leaf whose run was entered by the seed and covers the whole
    walk -- its occurrences are the walk's (with the seed: the seed's occurrence it
    continues, which precedes the walk on the record's forward strand on the right arm
    and follows it on the left)."""
    if rcs is None:
        return ''
    n_seed = g.seed.length_bp if with_seed else 0
    pairs = []
    cut = False
    for e in derive.end_labels(a, p.leaf):
        lab = g.labels[e.label]
        run = a.runs[e.run]
        if lab.kind != 'header' or run.from_label is not None or run.from_bp != 0 \
                or run.to_bp != p.length_bp:
            continue
        rc = rcs[e.run]
        cut = cut or rc.truncated
        for s, t in rc.intervals:
            if a.side == 'right':
                s -= n_seed
            else:
                t += n_seed
            if t > s:
                pairs.append((e.label, s, t))
    if not pairs:
        return ' coords_cut=1' if cut else ''
    pairs.sort()
    text = ' coords=' + ','.join('%s:%d-%d' % (g.labels[l].name, s + 1, t)
                                 for l, s, t in pairs[:_COORDS_IN_HEADER])
    if len(pairs) > _COORDS_IN_HEADER:
        text += ' coords_more=%d' % (len(pairs) - _COORDS_IN_HEADER)
    if cut:
        text += ' coords_cut=1'
    return text


def _fasta_coord_prices(bud, g, a, rcs):
    """(lwu, bytes) per walk (path id) of its coords= header field: the end labels read,
    the occurrences formatted and sorted, the field's text (cached)."""
    def prices():
        longest = max((len(l.name) for l in g.labels), default=0)
        pw, pm = [], []
        for p in derive.paths(a):
            ends = derive.end_labels(a, p.leaf)
            n = sum(len(rcs[e.run].intervals) for e in ends)
            shown = min(n, _COORDS_IN_HEADER)
            pw.append(W_RECORD * (n > 0) + len(ends) + _B.W_COORD * n + _B.sort_work(n))
            pm.append(3 * (shown * (longest + 44) + 40) + n * (_B.TUPLE + 3 * _B.INT))
        return pw, pm
    return derive.price_list(bud, g, a, 'fasta_coord_prices', prices)


def _fasta_record(g, side, a, p, w, orientation, with_seed, width, rcs=None):
    """(header, bases) of walk |p| whose walking-order bases |w| are spelled."""
    seg = a.segments[p.leaf]
    reason = (REASON[seg.leaf.path_reason] if seg.leaf.path_reason
              else ','.join(sorted({a.runs[e.run].reason
                                    for e in derive.end_labels(a, p.leaf)})) or '.')
    n = len(derive.end_labels(a, p.leaf)) if g.mode == 'constrain' else len(seg.end)
    return ('>%s_%d length_bp=%d reason=%s labels=%d seed=%s%s' % (
        side, p.id, p.length_bp, reason, n, 'yes' if with_seed else 'no',
        _coords_header(g, a, p, rcs, with_seed)),
        _wrap(_oriented(g, side, w, orientation, with_seed), width))


def _fasta_records(g, side, a, chosen, orientation, with_seed, width, out, ordered=False,
                   rcs=None):
    """Append the records of |chosen| (in its order) to |out|, each spelled as the stream
    of derive.walk_iter() reaches it and dropped once written: the text and one walk's
    bases are held, not every chosen walk's (holding them all, and the batch's prefixes,
    made to_fasta() peak at 10x its text on a comb, 3.1x before the batch). |ordered|:
    |chosen| is known to be in leaf order (every walk)."""
    order = [p.leaf for p in chosen]
    # a missing base raised part-way leaves records in |out| that the raise discards with
    # it: no read of the chains first (derive.walk_iter())
    stream = derive.walk_iter(a, order, chains=False, check_first=False)
    if ordered or all(map(_lt, order, order[1:])):
        # every walk, or leaves given in segment order: the stream's own order. The record
        # is _fasta_record() written out: a call per record made small exports 20 % slower
        segs = a.segments
        constrain = g.mode == 'constrain'
        seed = 'yes' if with_seed else 'no'
        for p, (_, _, w) in zip(chosen, stream):
            seg = segs[p.leaf]
            reason = (REASON[seg.leaf.path_reason] if seg.leaf.path_reason
                      else ','.join(sorted({a.runs[e.run].reason
                                            for e in derive.end_labels(a, p.leaf)})) or '.')
            n = len(derive.end_labels(a, p.leaf)) if constrain else len(seg.end)
            if rcs is None:
                out.append('>%s_%d length_bp=%d reason=%s labels=%d seed=%s' % (
                    side, p.id, p.length_bp, reason, n, seed))
            else:
                out.append('>%s_%d length_bp=%d reason=%s labels=%d seed=%s%s' % (
                    side, p.id, p.length_bp, reason, n, seed,
                    _coords_header(g, a, p, rcs, with_seed)))
            out.append(_wrap(_oriented(g, side, w, orientation, with_seed), width))
        return
    # leaves in another order or repeated: each record into its slot
    at = {}
    for i, leaf in enumerate(order):
        at.setdefault(leaf, []).append(i)
    first = len(out)
    out.extend([None] * (2 * len(chosen)))
    for leaf, _, w in stream:
        for i in at[leaf]:
            out[first + 2 * i], out[first + 2 * i + 1] = _fasta_record(
                g, side, a, chosen[i], w, orientation, with_seed, width, rcs)


def _to_fasta(g, arm, leaves, with_seed, orientation, width, bud, done=None,
              coordinates=None):
    from . import ops
    out = []
    total = 0
    sides = [g.arm(arm).side] if arm is not None else [s for s in ARM_SIDES if s in g.arms]
    for side in sides:
        a = g.arms[side]
        rcs = _fasta_rcs(g, side, coordinates)
        # the walks chosen on THIS arm, read once: a one-shot iterator (a generator, map(),
        # iter()) has no len(), and the test below raised TypeError on it where the
        # unbudgeted export had always iterated it (O11, the review of 2026-10-06: an
        # unbudgeted answer changed). Read per arm, as the selection always was: with
        # arm=None a one-shot iterator is exhausted by the first arm, the next reads none
        lv = leaves if leaves is None or isinstance(leaves, Sized) else list(leaves)
        if bud is not None and lv is not leaves:
            # the copy of the caller's iterator, charged once it is read: its length is
            # known only then (budget.py's stated exception for sizes known when built)
            bud.charge(len(lv) >> ID_SHIFT_PY, list_bytes(len(lv)))
        if lv is not None and not len(lv):
            # no walk chosen: nothing to resolve, so the arm's paths are not built -- they
            # were, uncharged, and a refused export (stopped at its join) left them behind
            # (the review of the level 4-5 fixes, finding 2)
            continue
        if bud is not None and a.segments:
            # admitted before the walks are resolved: the arm's paths (a Path per walk)
            # are built only once charged, so a refused export leaves none behind
            derive.uses(bud, g, a, 'leaves', 'paths', 'end_labels')
        paths = derive.paths(a)
        chosen = paths if lv is None else [paths[ops.path_id(a, x)] for x in lv]
        if not chosen:
            continue
        if orientation not in ('natural', 'walk'):
            raise ValueError("orientation is 'natural' or 'walk'")
        # the chosen walks' bases, their common prefixes spelled once (_fasta_records)
        if not a.segments or a.segments[0].walk is None:
            raise ValueError(derive.NO_BASES)
        if bud is not None:
            if leaves is None:
                # every walk: the spelling's bound is the arm's structural totals, streamed
                # (the records charge each walk's bases as they are written)
                z = derive.arm_sizes(a)
                bud.charge(len(chosen))
                bud.charge(2 * W_STEP * z['leaf_depth'] + (z['walk_bp'] >> 8),
                           derive.walk_iter_bytes(len(chosen), z['walk_bp'], z['segments']))
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
            if rcs is not None:
                # the coords= field of every header (the index charged once)
                _C.uses(bud, g)
                cw, cm = _fasta_coord_prices(bud, g, a, rcs)
                pw = [x + cw[p.id] for x, p in zip(pw, chosen)]
                pm = [x + cm[p.id] for x, p in zip(pm, chosen)]
                total += sum([cm[p.id] for p in chosen]) // 3
        if bud is not None:
            for _ in _blocks(bud, pw, pm, done, len(chosen)):
                pass
        _fasta_records(g, side, a, chosen, orientation, with_seed, width, out,
                       ordered=leaves is None, rcs=rcs)
    if bud is not None:
        bud.phase = 'join'
        bud.charge(total >> 8, STR + total)
    if out:
        # the final line feed joined with the rest: adding it to the joined text copied the
        # whole text once more at the peak
        out.append('')
    return '\n'.join(out)


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
            # the chains streamed, each dropped once its P line is built
            steps = derive.arm_sizes(arm)['leaf_depth']
            big = max([segs[p.leaf].depth + 1 for p in ps], default=0)
            bud.charge(2 * W_STEP * steps,
                       derive.walk_iter_bytes(len(ps), 0, len(segs), 0, steps)
                       + list_bytes(big))
        # each walk's chain as the stream reaches it (paths are in leaf order, the
        # stream's), dropped once its P line is written: holding every chain at once (8
        # bytes per segment of each) made to_gfa() peak 1.3x higher than chain by chain
        chains = derive.walk_iter(arm, [p.leaf for p in ps], spell=False)
        if bud is not None:
            pw, pm, tt = _gfa_path_prices(bud, g, arm)
            total += tt
            for _ in _blocks(bud, pw, pm, done):
                pass
        for p, (_, chain_, _) in zip(ps, chains):
            names = [_seg_name(side, x) + '+' for x in chain_]
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
        bud.charge(total >> 8, STR + total + list_bytes(len(lines) + len(links) + len(paths)))
    # one list and the final line feed joined with the rest (adding it to the joined text
    # copied the whole text once more at the peak)
    lines.extend(links)
    lines.extend(paths)
    lines.append('')
    return '\n'.join(lines)
