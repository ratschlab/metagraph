"""Hold a /traverse answer against the oracle (ORACLE-DESIGN.md section 10.3).

An arm is read back into walks (outward order) from its segments: each segment's walk is its first
parent's walk plus the bases it adds (a left arm's segment is written in natural orientation, so
reversed). From that come the leaves (paths), the claims (runs: the walk to the run's end, the label,
its reason, branches and loss), the labels recorded at every node (annotate mode's label_sets) and
the blocked / hairpin events. A complete arm must equal the oracle as complete sets; an arm a cap
truncated is judged on its invariants up to complete_to_bp (SPEC 6.10).
"""

import collections
import json

from .compare import Checks, INCONCLUSIVE, EXPLORATORY, PASS, FAIL  # noqa: F401
from . import traverse_oracle as T

CAP_REASONS = {'max_steps', 'time_budget', 'resource_limit', 'max_live_paths', 'max_paths',
               'max_output_bp', 'beam_pruned'}


class _Names(dict):
    def get(self, key, default=None):
        return dict.get(self, key, f'#{key}')


class ArmView:
    """An arm of a response, read back into walks."""

    def __init__(self, arm_json, arm, names):
        self.arm = arm
        self.json = arm_json
        # a label id outside label_dict reads as '#<id>', never as None (a corrupted answer)
        names = _Names(names)
        self.names = names
        self.problems = []
        seg_list = arm_json.get('segments') or []
        segs = {s['id']: s for s in seg_list}
        self.segs = segs
        # multiplicities (verification H1): sets alone would hide a listing twice
        self.duplicates = []
        if len(segs) != len(seg_list):
            self.duplicates.append(('segment id twice', len(seg_list) - len(segs)))
        self.walk = {}
        for sid in sorted(segs):
            self._walk_of(sid)
        self.leaves = {}
        children = {}
        for s in seg_list:
            for par in (s.get('parents') or [])[:1]:
                children.setdefault(par, []).append(s)
        merged = any(len(s.get('parents') or []) > 1 for s in seg_list)
        self.order_problems = []
        for par, kids in children.items():
            firsts = [self._first_base(k) for k in kids]
            if len(set(firsts)) != len(firsts):
                self.duplicates.append(('siblings with one first base', par, firsts))
        if not merged:
            seen = {}
            for sid, w in self.walk.items():
                if w in seen and self.segs[sid].get('length_bp', 0) > 0:
                    self.duplicates.append(('two segments end at one walk', seen[w], sid, w))
                seen.setdefault(w, sid)
        leaf_seen = {}
        for p in arm_json.get('paths') or []:
            chain = p.get('segments') or []
            if not chain or chain[-1] not in self.walk:
                self.problems.append(('path without segments', p.get('id')))
                continue
            w = self.walk[chain[-1]]
            if len(w) != p.get('length_bp'):
                self.problems.append(('path length', p.get('id'), len(w), p.get('length_bp')))
            if w in self.leaves:
                self.duplicates.append(('a leaf walk twice', w))
            if chain[-1] in leaf_seen:
                self.duplicates.append(('a leaf segment in two paths', chain[-1]))
            leaf_seen[chain[-1]] = p.get('id')
            self.leaves[w] = p
        self.claims = set()
        self.runs = []
        self.label_ends = set()
        self.switches = set()
        counted = collections.Counter()
        for sid, s in segs.items():
            evs = s.get('events') or []
            # SPEC 7.1: events at one at_bp are unordered, the list is sorted by position
            if [e.get('at_bp', 0) for e in evs] != sorted(e.get('at_bp', 0) for e in evs):
                self.order_problems.append(('events not sorted by at_bp', sid))
            for e in evs:
                if sid not in self.walk:
                    continue
                counted[(sid, json.dumps(e, sort_keys=True))] += 1
                if e.get('type') == 'label_end':
                    nb = e.get('needed_budget')
                    self.label_ends.add((self.walk[sid][:e['at_bp']], names.get(e.get('label')),
                                         e.get('reason'),
                                         round(float(nb), 9) if nb is not None else None))
                elif e.get('type') == 'switch':
                    key = (self.walk[sid][:e.get('at_bp', 0) + 1], names.get(e.get('from')),
                           names.get(e.get('to')), round(float(e.get('cost', 0.0)), 9))
                    if key in self.switches:
                        self.duplicates.append(('a switch twice', key))
                    self.switches.add(key)
        self.duplicates += [('an event twice', k) for k, n in counted.items() if n > 1]
        for sp in arm_json.get('splits') or []:
            chars = [b.get('char') for b in sp.get('branches') or []]
            # SPEC 6.2: successors are processed in ascending character order
            if chars != sorted(chars):
                self.order_problems.append(('split branches not in character order', chars))
        for r in arm_json.get('runs') or []:
            w = self.walk.get(r.get('segment'))
            if w is None or r.get('to_bp', 0) > len(w):
                self.problems.append(('run outside its segment', r))
                continue
            name = names.get(r.get('label'))
            claim = (w[:r['to_bp']], name, r.get('end_reason'), r.get('branches'),
                     round(float(r.get('loss', 0.0)), 9), r.get('from_bp'), r.get('entered_by'))
            if claim in self.claims and not merged:
                self.duplicates.append(('a run twice', claim))
            self.claims.add(claim)
            self.runs.append((claim, r))
        raw = collections.Counter(json.dumps(r, sort_keys=True) for r in arm_json.get('runs') or [])
        self.duplicates += [('a run record twice', k) for k, n in raw.items() if n > 1]
        self.node_labels = {}
        self.events = set()
        for sid, s in segs.items():
            base = self.walk[sid]
            start = s.get('from_bp', 0)
            for ls in s.get('label_sets') or []:
                labs = frozenset(names.get(i) for i in ls.get('labels') or [])
                for j in range(ls['from_bp'] + 1, ls['to_bp'] + 1):
                    self.node_labels[base[:j]] = labs
            if not s.get('parents'):
                self.node_labels.setdefault('', frozenset(names.get(i)
                                                          for i in s.get('labels') or []))
            for e in s.get('events') or []:
                if e.get('type') == 'blocked':
                    self.events.add((base[:e['at_bp']], e.get('char'), 'blocked', e.get('reason')))
                elif e.get('type') == 'hairpin':
                    self.events.add((base[:e['at_bp']], e.get('char'), 'hairpin',
                                     'followed' if e.get('followed') else 'hairpin'))
            del start

    def _first_base(self, seg):
        seq = seg.get('sequence') or ''
        if not seq:
            return None
        return seq[0] if self.arm == 'right' else seq[-1]

    def _walk_of(self, sid):
        if sid in self.walk:
            return self.walk[sid]
        s = self.segs[sid]
        seq = s.get('sequence') or ''
        out = seq if self.arm == 'right' else seq[::-1]
        parents = s.get('parents') or []
        prefix = self._walk_of(parents[0]) if parents and parents[0] in self.segs else ''
        w = prefix + out
        if len(w) != s.get('from_bp', 0) + s.get('length_bp', 0):
            self.problems.append(('segment length', sid, len(w), s.get('from_bp'),
                                  s.get('length_bp')))
        self.walk[sid] = w
        return w

    def all_walks(self):
        out = set()
        for w in self.walk.values():
            for j in range(len(w) + 1):
                out.add(w[:j])
        return out


def _diff(want, got, limit=10):
    want, got = set(want), set(got)
    return {'missing': sorted(want - got, key=repr)[:limit],
            'extra': sorted(got - want, key=repr)[:limit]}


def compare_arm(checks, w, seed, arm, arm_json, names, strategy, labels, extra, mode,
                label_objs=None):
    """Checks of one arm of a seed's result."""
    view = ArmView(arm_json, arm, names)
    prefix = f'traverse.{arm}'
    checks.add(f'{prefix}.shape', not view.problems, 'segments, paths and runs consistent',
               view.problems[:10])
    checks.add(f'{prefix}.unique', not view.duplicates,
               'no segment, leaf, run or event listed twice; siblings part on their first base',
               view.duplicates[:10])
    checks.add(f'{prefix}.order', not view.order_problems,
               'events sorted by at_bp, split branches in character order',
               view.order_problems[:10])
    st = T.Strategy(strategy)
    R = st.radius
    status = arm_json.get('status')
    ctb = arm_json.get('complete_to_bp')
    keep = st.on_reconverge == 'keep'
    checks.add(f'{prefix}.scope', arm_json.get('completeness_scope') ==
               ('per_path' if keep else 'united_history'),
               'per_path' if keep else 'united_history', arm_json.get('completeness_scope'))
    complete = status == 'complete'
    checks.add(f'{prefix}.status', (complete and ctb == R) or
               (status in ('truncated', 'pruned') and isinstance(ctb, int) and 0 <= ctb < R),
               f'complete with complete_to_bp {R}, or truncated below it',
               {'status': status, 'complete_to_bp': ctb})
    if not complete:
        checks.mark('completeness', INCONCLUSIVE,
                    note=f'{arm}: a cap truncated the arm at complete_to_bp {ctb}')
        limits = [l.get('kind') for l in arm_json.get('limitations') or []]
        checks.add(f'{prefix}.walk_domain', 'walk_domain' in limits
                   and isinstance(arm_json.get('cap_trigger'), dict),
                   'a walk_domain limitation and cap_trigger', limits)
        ctb = ctb if isinstance(ctb, int) else 0
    if mode == 'annotate':
        trie = T.annotate_trie(w, seed, arm, R, labels, st.hairpins)
        _compare_annotate(checks, prefix, trie, view, complete, ctb)
    else:
        res = T.constrain_walk(w, seed, arm, labels, st, extra)
        merged = any(len(s.get('parents') or []) > 1 for s in arm_json.get('segments') or [])
        if keep or not merged:
            if not keep:
                # merge with no merge done: the united history is each path's own (SPEC 6.5,
                # 6.10), so the arm is held to the keep walker strictly (verification H5)
                _scope_limitation(checks, prefix, arm_json, 0)
            _compare_constrain(checks, prefix, res, view, complete, ctb)
        else:
            _merge_invariants(checks, prefix, res, view, arm_json, w, seed, arm,
                              {l.name: l for l in list(labels) + list(extra)})
    return view


def _compare_annotate(checks, prefix, trie, view, complete, ctb):
    want_walks = set(trie.nodes)
    got_walks = view.all_walks()
    if complete:
        checks.add(f'{prefix}.leaves', view_leaves(view) == trie.leaves, 'the trie\'s leaves',
                   _leaf_diff(trie.leaves, view_leaves(view)))
        checks.add(f'{prefix}.walks', got_walks == want_walks, 'every walk of the trie',
                   _diff(want_walks, got_walks))
        got_labels = {k: v for k, v in view.node_labels.items() if k in want_walks}
        want_labels = {k: trie.nodes[k] for k in got_walks if k in trie.nodes}
        bad = [(k, sorted(want_labels[k]), sorted(got_labels.get(k, ())))
               for k in sorted(want_labels, key=len) if got_labels.get(k) != want_labels[k]]
        checks.add(f'{prefix}.node_labels', not bad, 'the labels present at every node', bad[:10])
        checks.add(f'{prefix}.events', view.events == trie.events,
                   'every blocked and hairpin step', _diff(trie.events, view.events))
        return
    # truncated: the trie up to complete_to_bp, nothing unsound beyond it
    checks.add(f'{prefix}.sound', got_walks <= want_walks, 'every walk an oracle walk',
               _diff(set(), got_walks - want_walks))
    upto = {x for x in want_walks if len(x) <= ctb}
    checks.add(f'{prefix}.complete_to_bp', upto <= got_walks,
               f'every walk of at most {ctb} bases', _diff(upto, got_walks & upto))
    early = {k: v for k, v in trie.leaves.items() if len(k) < ctb}
    got_early = {k: v for k, v in view_leaves(view).items() if len(k) < ctb}
    checks.add(f'{prefix}.leaves_before_ctb', got_early == early,
               f'every leaf ending before {ctb} with its reason', _leaf_diff(early, got_early))


def view_leaves(view):
    return {w: p.get('path_reason') for w, p in view.leaves.items()}


def _leaf_diff(want, got):
    keys = set(want) | set(got)
    return [(k, want.get(k), got.get(k)) for k in sorted(keys, key=lambda x: (len(x), x))
            if want.get(k) != got.get(k)][:10]


def _compare_constrain(checks, prefix, res, view, complete, ctb):
    got_walks = view.all_walks()
    if complete:
        checks.add(f'{prefix}.walks', got_walks == res.walks, 'every walk of the constrained trie',
                   _diff(res.walks, got_walks))
        want = served_claims(res.claims)
        checks.add(f'{prefix}.claims', view.claims == want,
                   'every label end (walk, label, reason, branches, loss, from_bp, entered_by)',
                   _diff(want, view.claims))
        want_events = {(c[0], c[1], c[2], res.needed_at.get((c[0], c[1]))
                        if c[2] == 'loss_budget' else None)
                       for c in res.claims if c[2] != 'switched'}
        checks.add(f'{prefix}.label_end_events', view.label_ends == want_events,
                   'every label end event with its SPEC 6.7 reason (and needed_budget)',
                   _diff(want_events, view.label_ends))
        got_needed = sorted(round(float(x), 9) for x in view.json.get('needed_budgets') or [])
        want_needed = sorted(round(x, 9) for x in res.needed)
        checks.add(f'{prefix}.needed_budgets', got_needed == want_needed, want_needed,
                   got_needed)
        checks.add(f'{prefix}.switches', view.switches == res.switches,
                   'every switch with its source, target and cost',
                   _diff(res.switches, view.switches))
        want_leaves = {w: ends for w, ends in res.leaves.items()}
        got_leaves = {}
        for wk, p in view.leaves.items():
            got_leaves[wk] = {view.names.get(e.get('label')): None for e in
                              p.get('end_labels') or []}
        missing = sorted(set(want_leaves) - set(got_leaves), key=lambda x: (len(x), x))
        extra = sorted(set(got_leaves) - set(want_leaves), key=lambda x: (len(x), x))
        checks.add(f'{prefix}.leaves', not missing and not extra, 'the leaves',
                   {'missing': missing[:10], 'extra': extra[:10]})
        bad = []
        for wk, p in view.leaves.items():
            if wk not in want_leaves:
                continue
            ends = want_leaves[wk]
            hist = {}
            for r in ends.values():
                hist[T.served(r)] = hist.get(T.served(r), 0) + 1
            got_names = sorted(view.names.get(e.get('label')) for e in p.get('end_labels') or [])
            if got_names != sorted(ends) or p.get('end_reasons') != hist:
                bad.append((wk, ends, got_names, p.get('end_reasons')))
        checks.add(f'{prefix}.end_labels', not bad, 'each leaf\'s end labels and reasons', bad[:10])
        checks.add(f'{prefix}.events', view.events == res.events,
                   'every blocked and hairpin step a live label would have taken',
                   _diff(res.events, view.events))
        return
    checks.add(f'{prefix}.sound', got_walks <= res.walks, 'every walk an oracle walk',
               _diff(set(), got_walks - res.walks))
    upto = {x for x in res.walks if len(x) <= ctb}
    checks.add(f'{prefix}.complete_to_bp', upto <= got_walks,
               f'every walk of at most {ctb} bases', _diff(upto, got_walks & upto))
    early = {c for c in served_claims(res.claims) if len(c[0]) < ctb}
    got_early = {c for c in view.claims if len(c[0]) < ctb}
    checks.add(f'{prefix}.claims_before_ctb', got_early == early,
               f'every label end before {ctb}', _diff(early, got_early))
    late = {c for c in view.claims if len(c[0]) >= ctb and c[2] not in CAP_REASONS}
    served_all = served_claims(res.claims)
    checks.add(f'{prefix}.claims_after_ctb', late <= served_all,
               'a label end at or beyond complete_to_bp is an oracle claim or a cap end',
               _diff(set(), late - served_all))


def served_claims(claims):
    """The oracle's claims with their reasons as a run states them (SPEC 13, note I3-2)."""
    return {(c[0], c[1], T.served(c[2])) + tuple(c[3:]) for c in claims}


def _scope_limitation(checks, prefix, arm_json, merges):
    limits = {l.get('kind'): l for l in arm_json.get('limitations') or []}
    checks.add(f'{prefix}.scope_limitation', 'scope' in limits
               and limits['scope'].get('observed') == merges,
               f'a scope limitation observing {merges} merges', limits.get('scope'))


def _merge_invariants(checks, prefix, res, view, arm_json, w, seed, arm, by_name):
    """on_reconverge merge: the united history is not modelled (SPEC 6.10); its invariants are.
    Under merge the displayed walk follows first parents (SPEC 7.1), so a run's walk is the
    label's own route only after its merges: the label must hold every node of the run's
    displayed support -- the whole run when route_bp is 0, else from the latest merge its leaf
    states (end_labels' route_bp) -- and its end node in any case."""
    checks.mark('completeness', INCONCLUSIVE, note='merge: completeness is the united history\'s')
    got_walks = view.all_walks()
    checks.add(f'{prefix}.sound', got_walks <= res.walks, 'every walk an oracle walk',
               _diff(set(), got_walks - res.walks))
    segs = arm_json.get('segments') or []
    merged_segs = [s for s in segs if len(s.get('parents') or []) > 1]
    _scope_limitation(checks, prefix, arm_json, len(merged_segs))
    # SPEC 6.5: a reconverge event {at_bp, segments} on each merged segment, at its start, naming
    # its parents; and no reconverge event elsewhere
    want = {(s['id'], s.get('from_bp'), tuple(s.get('parents'))) for s in merged_segs}
    got = {(s['id'], e.get('at_bp'), tuple(e.get('segments') or ())) for s in segs
           for e in s.get('events') or [] if e.get('type') == 'reconverge'}
    checks.add(f'{prefix}.reconverge_events', got == want,
               'a reconverge event at each merge, naming its parents', _diff(want, got))
    leaf_route = {}
    for p in arm_json.get('paths') or []:
        for e in p.get('end_labels') or []:
            leaf_route[e.get('run')] = e.get('route_bp', 0)
    bad = []
    for i, (claim, r) in enumerate(view.runs):
        walk, name, reason = claim[0], claim[1], claim[2]
        label = by_name.get(name)
        if label is None:
            bad.append(('unknown label', claim))
            continue
        start = r.get('from_bp', 0)
        rb = r.get('route_bp', 0)
        if rb:
            start = max(start, leaf_route.get(i, len(walk)))
        nodes = range(max(start, 0), len(walk) + 1) if start < len(walk) else [len(walk)]
        missing = [j for j in nodes if not w.holds(label, w.walk_node(seed, arm, walk, j))]
        if missing:
            bad.append((claim, missing[:5]))
    checks.add(f'{prefix}.displayed_support', not bad,
               'every run\'s label on the nodes of its displayed support', bad[:10])


# -- a whole /traverse answer -------------------------------------------------------------

def seed_labels(w, seed_json, req_seed, strategy, mode):
    """(validated labels, dropped names, derived) as the oracle reads the request (SPEC 6.1)."""
    st = T.Strategy(strategy)
    if mode == 'annotate':
        kind = (strategy.get('labels') or {}).get('seed_label_kind', 'auto')
        if kind == 'auto':
            kind = 'header' if w.mode == 'basic' and w.index.seq_ids else 'column'
        return w.labels_of_kind(kind), [], False
    names = req_seed.get('labels')
    if names is None:
        kind = (strategy.get('labels') or {}).get('seed_label_kind', 'auto')
        if kind == 'auto':
            kind = 'header' if w.mode == 'basic' and w.index.seq_ids else 'column'
        cands = [l for l in w.labels_of_kind(kind)
                 if T.supports_seed(w, l, req_seed['sequence'], 'kmer')]
        cands.sort(key=lambda l: l.order)
        keep = [l for l in cands if T.supports_seed(w, l, req_seed['sequence'], st.support)]
        dropped = [l.name for l in cands if l not in keep]
        return keep, dropped, True
    labels = [w.label(n) for n in names]
    keep = [l for l in labels if l is not None
            and T.supports_seed(w, l, req_seed['sequence'], st.support)]
    dropped = [l.name for l in labels if l is not None and l not in keep]
    return keep, dropped, False


def expected_refusal(request):
    """A refusal the oracle predicts (SPEC 5): extra labels are switch targets, rejected when
    unreachable -- under forbid, or when no chain of switches from a seed label enters them
    within loss_budget (a constant cost: one switch of cost c, so reachable iff c <= budget)."""
    lab = (request.get('strategy') or {}).get('labels') or {}
    if not lab.get('extra') or lab.get('mode', 'constrain') != 'constrain':
        return None
    cost = lab.get('change_cost') or {'model': 'forbid'}
    if cost.get('model', 'forbid') == 'forbid':
        return 'strategy.labels.extra'
    if cost.get('model') == 'constant' and float(cost.get('value', 0)) > float(
            lab.get('loss_budget', 0)) + 1e-12:
        return 'unreachable'
    return None


def compare_traverse(w, request, status, answer):
    """-> [(seed position, Checks)]."""
    out = []
    strategy = request.get('strategy') or {}
    mode = (strategy.get('labels') or {}).get('mode', 'constrain')
    refusal = expected_refusal(request)
    if refusal is not None:
        c = Checks()
        text = answer.get('error', '') if isinstance(answer, dict) else str(answer)
        c.add('traverse.refusal', status == 400 and refusal in text,
              f'400 naming {refusal!r}', [status, text[:300]])
        return [(None, c)]
    if status != 200 or not isinstance(answer, dict):
        c = Checks()
        c.add('traverse.status', False, 200, [status, answer if not isinstance(answer, dict)
                                              else answer.get('error')])
        return [(None, c)]
    results = answer.get('results') or []
    if len(results) != len(request['seeds']):
        c = Checks()
        c.add('traverse.results', False, len(request['seeds']), len(results))
        return [(None, c)]
    direction = strategy.get('direction', 'both')
    arms = ['left', 'right'] if direction == 'both' else [direction]
    extra_names = (strategy.get('labels') or {}).get('extra') or []
    for i, (req_seed, res) in enumerate(zip(request['seeds'], results)):
        checks = Checks()
        seed = req_seed['sequence'].upper()
        labels, dropped, derived = seed_labels(w, res.get('seed'), req_seed, strategy, mode)
        names = {j: d.get('name') for j, d in enumerate(res.get('label_dict') or [])}
        if mode == 'constrain':
            got = res.get('seed') or {}
            got_labels = got.get('labels') or []
            want = [l.name for l in labels]
            checks.add('traverse.seed.labels',
                       [x if isinstance(x, str) else x.get('name') for x in got_labels] == want
                       if derived else
                       sorted(x if isinstance(x, str) else x.get('name') for x in got_labels)
                       == sorted(want), want, got_labels)
            got_dropped = sorted(d.get('label') if isinstance(d, dict) else d
                                 for d in got.get('dropped_labels') or [])
            checks.add('traverse.seed.dropped', got_dropped == sorted(dropped)
                       and all(d.get('reason') == 'seed_unsupported'
                               for d in got.get('dropped_labels') or []),
                       sorted(dropped), got.get('dropped_labels'))
        checks.add('traverse.label_mode', res.get('label_mode') == mode, mode, res.get('label_mode'))
        extra = [w.label(n) for n in extra_names]
        for arm in arms:
            arm_json = (res.get('arms') or {}).get(arm)
            if not isinstance(arm_json, dict):
                checks.add(f'traverse.{arm}', False, 'the arm', None)
                continue
            compare_arm(checks, w, seed, arm, arm_json, names, strategy, labels,
                        [l for l in extra if l is not None], mode)
        out.append((i, checks))
    return out


# -- output representations (SPEC 7.1: tree is full without sequences, summary the per-leaf view) --

def _strip_sequences(x):
    if isinstance(x, dict):
        return {k: _strip_sequences(v) for k, v in x.items() if k not in ('sequence', 'timing')}
    if isinstance(x, list):
        return [_strip_sequences(v) for v in x]
    return x


SUMMARY_PATH_KEYS = ('id', 'length_bp', 'n_labels', 'end_reasons', 'path_reason', 'continuation')


def compare_details(full, tree, summary):
    """The same request in detail full, tree and summary: tree equals full without its sequences;
    summary keeps every field of an arm but segments, runs and splits, and of each path its
    length, label count, end reasons, path reason and continuation."""
    checks = Checks()
    for name, other in (('tree', tree), ('summary', summary)):
        if not isinstance(other, dict) or len(other.get('results') or []) != len(full['results']):
            checks.add(f'traverse.detail.{name}', False, 'one result per seed', other)
    if checks.failed():
        return checks
    for f, t, s in zip(full['results'], tree['results'], summary['results']):
        fa = _strip_sequences(f.get('arms') or {})
        ta = _strip_sequences(t.get('arms') or {})
        checks.add('traverse.detail.tree', fa == ta, 'full without its sequences',
                   [k for a in fa for k in fa[a] if fa[a][k] != (ta.get(a) or {}).get(k)][:10])
        bad = []
        for a, arm in (f.get('arms') or {}).items():
            sa = (s.get('arms') or {}).get(a) or {}
            for key, value in arm.items():
                if key in ('segments', 'runs', 'splits'):
                    if key in sa:
                        bad.append((a, key, 'present'))
                    continue
                if key == 'paths':
                    want = [{k: p[k] for k in SUMMARY_PATH_KEYS if k in p} for p in value]
                    got = [{k: p[k] for k in SUMMARY_PATH_KEYS if k in p}
                           for p in sa.get('paths') or []]
                    if _strip_timing(want) != _strip_timing(got):
                        bad.append((a, 'paths'))
                elif _strip_timing(value) != _strip_timing(sa.get(key)):
                    bad.append((a, key))
        for key in ('seed', 'label_dict', 'label_summary', 'outcome', 'limitations', 'label_mode'):
            if f.get(key) != s.get(key):
                bad.append(('seed result', key))
        checks.add('traverse.detail.summary', not bad, 'the per-leaf view of full', bad[:10])
    return checks


def _strip_timing(x):
    if isinstance(x, dict):
        return {k: _strip_timing(v) for k, v in x.items() if k != 'timing'}
    if isinstance(x, list):
        return [_strip_timing(v) for v in x]
    return x
