"""The /traverse oracle (ORACLE-DESIGN.md section 10.2): the walks of a seed computed from the
records' text, sharing no code with the engine or the suite.

SPEC-labeled-traversal-core.md section 6.9's contract, with a model built from the records:
- annotate_trie(): T, every structural walk up to the radius under the walk rule (section 6.10), each node
  with the labels present on it, each leaf with its path reason, every blocked or hairpin step;
- constrain_walk(): a reference walker of sections 6.3-6.4 under on_reconverge keep (a trie), path by path:
  the label state with its loss recurrence (forbid, constant), the branch allowance, the quorum,
  min_live_labels, the split limit and trace support; its claims are the label ends (runs) with their
  walk, reason, branches and loss;
- label_walks(): E, the union of single-label maximal walks (the direct statement of section 6.9 for
  forbid, exhaustive), a second and simpler computation of the same claims.

Walks are written in outward order (the left arm's bases from the seed outwards).
"""

from . import oracle as O

INF = float('inf')
BLOCK_RANK = {'rejoined_seed': 3, 'edge_reuse_rc': 2, 'edge_reuse': 1}
# SPEC 13, implementation note I3-2: the served EndReason has no switched, superseded, minority,
# below_min_labels, split_limit or hairpin; a run's end_reason is the value below, and the label's
# label_end event states the SPEC 6.7 reason as its qualifier (a switched source has no event)
SERVED_REASON = {'switched': 'label_lost', 'superseded': 'label_lost', 'minority': 'branch',
                 'below_min_labels': 'branch', 'split_limit': 'branch', 'hairpin': 'dead_end'}


def served(reason):
    return SERVED_REASON.get(reason, reason)


class Label:
    """A permitted label: a record's header or a column, with what it holds."""

    def __init__(self, kind, name, records, column_id, seq_id):
        self.kind = kind
        self.name = name
        self.records = records          # indexed records of the label
        self.column_id = column_id
        self.seq_id = seq_id

    @property
    def order(self):
        return (self.column_id, self.seq_id if self.seq_id is not None else -1)

    def __repr__(self):
        return f'Label({self.kind}, {self.name})'


class Walker:
    """The walk graph of an index (oracle.Index) with its labels: successors, edges, hairpins."""

    def __init__(self, index):
        self.index = index
        self.k = index.k
        self.mode = index.mode
        self.kmers = index.kmers
        self.column_ids = {c: i for i, c in enumerate(index.columns)}
        self._holds = {}
        # column coordinates of every indexed record's k-mers (SPEC 7.1: offset_i + j)
        self.coords = {}           # (kind, name) -> {coordinate: kmer}
        for r in index.indexed:
            off = index.column_offset[(r.column, r.name)]
            for p in range(len(r.seq) - self.k + 1):
                x = r.seq[p:p + self.k]
                if all(c in O.BASES for c in x):
                    self.coords.setdefault(('header', r.name), {})[p] = x
                    self.coords.setdefault(('column', r.column), {})[off + p] = x

    # -- labels ---------------------------------------------------------------------------

    def label(self, name):
        """A request label name: a column, else an indexed record's header (SPEC 4.1; a column
        first, as a name of both kinds would be)."""
        idx = self.index
        if name in idx.columns:
            return Label('column', name, [r for r in idx.indexed if r.column == name],
                         self.column_ids[name], None)
        for r in idx.indexed:
            if r.name == name:
                return Label('header', name, [r], self.column_ids[r.column],
                             idx.seq_ids[(r.column, r.name)])
        return None

    def labels_of_kind(self, kind):
        idx = self.index
        if kind == 'column':
            return [Label('column', c, [r for r in idx.indexed if r.column == c],
                          self.column_ids[c], None) for c in idx.columns]
        return [Label('header', r.name, [r], self.column_ids[r.column],
                      idx.seq_ids[(r.column, r.name)]) for r in idx.indexed]

    def holds(self, label, node):
        key = (label.kind, label.name, node)
        if key not in self._holds:
            rc = O.revcomp(node) if self.mode != 'basic' else None
            self._holds[key] = any(node in r.seq or (rc is not None and rc in r.seq)
                                   for r in label.records)
        return self._holds[key]

    def node_labels(self, labels, node):
        return frozenset(l.name for l in labels if self.holds(l, node))

    def chain_coords(self, label, node):
        """The label's coordinates of |node| (trace support, BASIC only)."""
        table = self.coords.get((label.kind, label.name), {})
        return {g for g, x in table.items() if x == node}

    # -- the walk graph -------------------------------------------------------------------

    def successors(self, node, arm):
        """[(char, next node, edge (k+1)-mer as written)] in ascending character order; '$' never."""
        out = []
        for c in 'ACGT':
            if arm == 'right':
                v, e = node[1:] + c, node + c
            else:
                v, e = c + node[:-1], c + node
            if v in self.kmers:
                out.append((c, v, e))
        return out

    def edge_key(self, e):
        """The edge's identity (SPEC 6.6): the (k+1)-mer on BASIC; canonical with its orientation
        bit on CANONICAL and PRIMARY."""
        if self.mode == 'basic':
            return (e, 0)
        rc = O.revcomp(e)
        return (e, 0) if e <= rc else (rc, 1)

    def is_hairpin(self, u, v, e):
        if self.mode == 'basic':
            return False
        if e == O.revcomp(e):
            return True
        return self.k % 2 == 0 and (u == O.revcomp(u) or v == O.revcomp(v))

    def seed_nodes(self, seed):
        k = self.k
        nodes = {seed[i:i + k] for i in range(len(seed) - k + 1)}
        if self.mode != 'basic':
            nodes |= {O.revcomp(x) for x in nodes}
        return nodes

    def block(self, u, v, e, used, seed_nodes, hairpins):
        """Why the step u -> v is not admissible, or None (SPEC 6.5, 6.6; precedence
        rejoined_seed > edge_reuse_rc > edge_reuse). A skipped hairpin comes first: the SPEC ranks
        it against no other block, and the served walker decides it on the step itself (a step
        into a seed node that is a palindromic k-mer at even k is a hairpin; the run report states
        this reading)."""
        if hairpins == 'skip' and self.is_hairpin(u, v, e):
            return 'hairpin'
        if v in seed_nodes:
            return 'rejoined_seed'
        key = self.edge_key(e)
        rc_key = (key[0], 1 - key[1])
        if self.mode != 'basic' and rc_key in used and rc_key != key:
            return 'edge_reuse_rc'
        if key in used:
            return 'edge_reuse'
        return None

    def walk_node(self, seed, arm, walk, j):
        """The node a walk has reached after j steps (j = 0: the seed boundary)."""
        k = self.k
        if arm == 'right':
            s = seed[-k:] + walk[:j]
            return s[-k:]
        s = walk[:j][::-1] + seed[:k]
        return s[:k]


def start_node(seed, arm, k):
    return seed[-k:] if arm == 'right' else seed[:k]


# -- annotate: the structural trie T --------------------------------------------------------

class Trie:
    def __init__(self):
        self.nodes = {}       # walk -> frozenset of label names at its last node ('' = boundary)
        self.leaves = {}      # walk -> path reason
        self.events = set()   # (walk, char, kind, reason)  kind: blocked | hairpin


def annotate_trie(w, seed, arm, radius, labels, hairpins='skip'):
    """T (SPEC 6.9): every walk of at most |radius| bases from the seed boundary under the walk rule,
    the labels of |labels| present on every node, each leaf's path reason."""
    t = Trie()
    seed_nodes = w.seed_nodes(seed)
    root = start_node(seed, arm, w.k)
    t.nodes[''] = w.node_labels(labels, root)

    def visit(walk, node, used):
        if len(walk) == radius:
            t.leaves[walk] = 'max_extension_bp'
            return
        nxt, blocks = [], []
        for c, v, e in w.successors(node, arm):
            b = w.block(node, v, e, used, seed_nodes, hairpins)
            if b is None:
                nxt.append((c, v, e))
                if hairpins == 'follow' and w.is_hairpin(node, v, e):
                    t.events.add((walk, c, 'hairpin', 'followed'))
            else:
                blocks.append(b)
                t.events.add((walk, c, 'hairpin' if b == 'hairpin' else 'blocked', b))
        if not nxt:
            real = [b for b in blocks if b != 'hairpin']
            t.leaves[walk] = (max(real, key=BLOCK_RANK.get) if real else 'dead_end')
            return
        for c, v, e in nxt:
            child = walk + c
            t.nodes[child] = w.node_labels(labels, v)
            visit(child, v, used | {w.edge_key(e)})

    visit('', root, frozenset())
    return t


# -- constrain: the reference walker --------------------------------------------------------

class Strategy:
    """The knobs the reference walker models, with the SPEC's defaults (section 5) or the
    exhaustive preset's (section 6.9)."""

    def __init__(self, request_strategy):
        s = request_strategy or {}
        ex = bool(s.get('exhaustive'))
        lab = s.get('labels') or {}
        br = s.get('branching') or {}
        bo = s.get('bounds') or {}
        cost = lab.get('change_cost') or {'model': 'forbid'}
        self.model = cost.get('model', 'forbid')
        self.cost_value = float(cost.get('value', 0)) if self.model == 'constant' else INF
        self.loss_budget = float(lab.get('loss_budget', 0))
        self.switch_on = lab.get('switch_on', 'loss')
        mlb = br.get('max_label_branches', 'unlimited' if ex else 0)
        self.max_branches = INF if mlb == 'unlimited' else int(mlb)
        msp = br.get('max_splits_per_path', 'unlimited' if ex else 64)
        self.max_splits = INF if msp == 'unlimited' else int(msp)
        self.min_succ_labels = int(br.get('min_successor_labels', 1))
        self.min_succ_fraction = float(br.get('min_successor_fraction', 0.0))
        self.min_live = int(bo.get('min_live_labels', 1))
        self.on_reconverge = br.get('on_reconverge', 'keep' if ex else 'merge')
        self.hairpins = br.get('hairpins', 'skip')
        self.radius = int(bo.get('max_extension_bp', 5000))
        self.support = s.get('support', 'kmer')

    def cost(self, a, b):
        if a.name == b.name:
            return 0.0
        return self.cost_value


class Entry:
    __slots__ = ('loss', 'branches', 'pred', 'run', 'coords', 'switch_cost')

    def __init__(self, loss, branches, pred, run, coords=None, switch_cost=None):
        self.loss = loss
        self.branches = branches
        self.pred = pred          # the predecessor label's name
        self.run = run            # (label name, from_bp, entered_by)
        self.coords = coords      # live chain coordinates under trace
        self.switch_cost = switch_cost   # the cost of the switch that entered it, if one did


class Result:
    def __init__(self):
        self.claims = set()       # (walk, label, reason, branches, loss, from_bp, entered_by)
        self.leaves = {}          # walk -> {label: reason} of the labels that end at the leaf
        self.walks = set()        # every walk prefix of the constrained trie
        self.needed = []          # needed budgets of loss_budget ends
        self.needed_at = {}       # (walk, label) -> needed_budget of a loss_budget end
        self.switches = set()     # (walk, from, to, cost) of every switch taken
        self.events = set()       # (walk, char, kind, reason) blocked / hairpin events


def seed_chain(w, label, seed, arm):
    """The label's coordinates at the seed boundary along a whole consecutive occurrence of the
    seed (trace), or an empty set."""
    k = w.k
    nodes = [seed[i:i + k] for i in range(len(seed) - k + 1)]
    table = w.coords.get((label.kind, label.name), {})
    starts = [g for g, x in table.items() if x == nodes[0]
              and all(table.get(g + j) == nodes[j] for j in range(len(nodes)))]
    if arm == 'right':
        return {g + len(nodes) - 1 for g in starts}
    return set(starts)


def supports_seed(w, label, seed, support):
    k = w.k
    if support == 'trace':
        return bool(seed_chain(w, label, seed, 'right'))
    return all(w.holds(label, seed[i:i + k]) for i in range(len(seed) - k + 1))


def step_chain(w, by_name, entry, v, arm):
    """The live coordinates of |entry|'s chain that continue to v (+1 right, -1 left)."""
    if entry.coords is None:
        return None
    d = 1 if arm == 'right' else -1
    lab = by_name[entry.run[0]]
    table = w.coords.get((lab.kind, lab.name), {})
    return {g + d for g in entry.coords if table.get(g + d) == v}


def label_in_a(w, by_name, sigma, tl, v, arm, trace):
    """tl in A(v): a permitted label present at v; under trace a live label only along its
    coordinate chain (SPEC 6.3)."""
    if not w.holds(tl, v):
        return False
    live = sigma.get(tl.name)
    return not (trace and live is not None and not step_chain(w, by_name, live, v, arm))


def derive_entry(w, st, by_name, sigma, tl, v, depth, arm='right', trace=False,
                 exclude=frozenset()):
    """The entry of target tl at successor v (SPEC 6.3): the cheapest of its stay and the
    switches from the sources -- under switch_on loss only a label not in A(v) -- by (loss, stay
    before switch, branches of the predecessor, column id / seq_id of the predecessor); None when
    none is within the budget. |exclude|: the sources a branch exclusion took out (6.4): neither
    their stay nor their switches count, so a removed target is re-minimised over the remaining
    sources, its own stay included."""
    best = None
    stay = sigma.get(tl.name)
    if stay is not None and tl.name not in exclude:
        coords = step_chain(w, by_name, stay, v, arm) if trace else None
        if not trace or coords:
            best = (stay.loss, 0, stay.branches, tl.order, tl.name, stay, coords)
    for sn, se in sigma.items():
        if sn == tl.name or sn in exclude:
            continue
        src = by_name[sn]
        if st.switch_on != 'any' and label_in_a(w, by_name, sigma, src, v, arm, trace):
            continue
        cost = st.cost(src, tl)
        if cost == INF:
            continue
        cand = (se.loss + cost, 1, se.branches, src.order, sn, se, None)
        if best is None or cand[:4] < best[:4]:
            best = cand
    if best is None or best[0] > st.loss_budget + 1e-12:
        return None
    loss, kind, branches, _, pred, pe, coords = best
    if kind == 0:
        return Entry(loss, branches, tl.name, pe.run, coords)
    if trace:
        coords = {g for g, x in w.coords.get((tl.kind, tl.name), {}).items() if x == v}
    return Entry(loss, branches, pred, (tl.name, depth, 'switch'), coords,
                 switch_cost=st.cost(by_name[pred], tl))


def constrain_walk(w, seed, arm, labels, strategy, extra=()):
    """The constrained trie A of |labels| (validated seed labels) and |extra| (switch targets) under
    |strategy| with on_reconverge keep: its claims, leaves, walks and events."""
    st = strategy
    res = Result()
    seed_nodes = w.seed_nodes(seed)
    root = start_node(seed, arm, w.k)
    permitted = list(labels) + list(extra)
    by_name = {l.name: l for l in permitted}
    trace = st.support == 'trace'
    sigma0 = {}
    for l in labels:
        coords = seed_chain(w, l, seed, arm) if trace else None
        sigma0[l.name] = Entry(0.0, 0, l.name, (l.name, 0, 'seed'), coords)

    def end(walk, name, entry, reason):
        res.claims.add((walk, name, reason, entry.branches, round(entry.loss, 9), entry.run[1],
                        entry.run[2]))

    def step_coords(entry, v):
        if entry.coords is None:
            return None
        d = 1 if arm == 'right' else -1
        table = w.coords.get((by_name[entry.run[0]].kind, by_name[entry.run[0]].name), {})
        return {g + d for g in entry.coords if table.get(g + d) == v}

    def carries(name, entry, v):
        """Whether the lineage would continue on v (SPEC 6.6's "would carry"): k-mer presence,
        or under trace the coordinate chain."""
        if trace and entry.coords is not None:
            return bool(step_coords(entry, v))
        return w.holds(by_name[name], v)

    def visit(walk, node, used, sigma, splits):
        res.walks.add(walk)
        depth = len(walk)
        if depth == st.radius:
            res.leaves[walk] = {n: 'max_extension_bp' for n in sigma}
            for n, e in sigma.items():
                end(walk, n, e, 'max_extension_bp')
            return
        succ = w.successors(node, arm)
        if not succ:
            res.leaves[walk] = {n: 'dead_end' for n in sigma}
            for n, e in sigma.items():
                end(walk, n, e, 'dead_end')
            return
        admissible, blocked = [], []
        for c, v, e in succ:
            b = w.block(node, v, e, used, seed_nodes, st.hairpins)
            if b is None:
                admissible.append((c, v, e))
            else:
                blocked.append((c, v, e, b))
        # the label states of the admissible successors (SPEC 6.3)
        def in_a(tl, v):
            return label_in_a(w, by_name, sigma, tl, v, arm, trace)

        def derive(tl, v, exclude=frozenset()):
            return derive_entry(w, st, by_name, sigma, tl, v, depth, arm, trace, exclude)

        states, targets = {}, {}
        trace_breaks = set()
        for c, v, e in admissible:
            sv = {}
            targets[c] = [tl for tl in permitted if in_a(tl, v)]
            for tl in permitted:
                if trace and w.holds(tl, v) and tl.name in sigma and not in_a(tl, v):
                    trace_breaks.add(tl.name)
            for tl in targets[c]:
                entry = derive(tl, v)
                if entry is not None:
                    sv[tl.name] = entry
            states[c] = sv
        vertex = {c: v for c, v, e in admissible}
        # ambiguity and the branch allowance (SPEC 6.4); a followed hairpin step never counts
        # toward ambiguity, its split is a divergence (SPEC 13, I3 note 11)
        hairpin_chars = {c for c, v, e in admissible if w.is_hairpin(node, v, e)}
        excluded = set()
        while True:
            feeds = {}
            for c, sv in states.items():
                if c in hairpin_chars:
                    continue
                for tn, te in sv.items():
                    feeds.setdefault(te.pred, set()).add(c)
            over = [s for s, cs in feeds.items() if len(cs) >= 2 and s not in excluded
                    and sigma[s].branches + 1 > st.max_branches]
            if not over:
                break
            excluded.update(over)
            # every entry derived from an excluded source leaves every successor; its target is
            # re-minimised over the remaining sources, its own stay included (6.3, 6.4)
            for c, sv in states.items():
                for tn in [tn for tn, te in sv.items() if te.pred in excluded]:
                    del sv[tn]
                    entry = derive(by_name[tn], vertex[c], frozenset(excluded))
                    if entry is not None:
                        sv[tn] = entry
        feeds = {}
        for c, sv in states.items():
            if c in hairpin_chars:
                continue
            for tn, te in sv.items():
                feeds.setdefault(te.pred, set()).add(c)
        for s, cs in feeds.items():
            if len(cs) >= 2:
                # an ambiguous source (decided over the non-hairpin successors): every entry
                # derived from it gets branches + 1, the one on a followed hairpin step too
                for c, sv in states.items():
                    for te in sv.values():
                        if te.pred == s:
                            te.branches = sigma[s].branches + 1
        # quorum and min_live_labels (SPEC 6.4)
        followed, minority, below = [], {}, {}
        for c, v, e in admissible:
            sv = states[c]
            if not sv:
                continue
            n = len(sv)
            if n < st.min_succ_labels or n < st.min_succ_fraction * len(sigma):
                for te in sv.values():
                    minority.setdefault(te.pred, True)
                continue
            if n < st.min_live:
                for te in sv.values():
                    below.setdefault(te.pred, True)
                continue
            followed.append((c, v, e))
        split_limited = len(followed) >= 2 and splits + 1 > st.max_splits
        if split_limited:
            res.leaves[walk] = {n: 'split_limit' for n in sigma}
            for n, e in sigma.items():
                end(walk, n, e, 'split_limit')
            return
        continuing = {te.pred for c, v, e in followed for te in states[c].values()}
        stays = {tn for c, v, e in followed for tn, te in states[c].items() if te.pred == tn}
        ended = {}
        for n, e in sigma.items():
            if n in stays:
                continue
            if n in excluded:
                reason = 'branch'
            elif n in continuing:
                reason = 'switched'
            elif n in minority:
                reason = 'minority'
            elif n in below:
                reason = 'below_min_labels'
            else:
                reason = _lost_reason(by_name[n], e, blocked, trace_breaks, states, targets,
                                      vertex, st, res, walk,
                                      lambda v, n=n, e=e: carries(n, e, v), in_a)
            ended[n] = reason
            end(walk, n, e, reason)
        for c, v, e, b in blocked:
            if any(carries(n, e_, v) for n, e_ in sigma.items()):
                res.events.add((walk, c, 'hairpin' if b == 'hairpin' else 'blocked', b))
        for c, v, e in followed:
            if st.hairpins == 'follow' and w.is_hairpin(node, v, e):
                res.events.add((walk, c, 'hairpin', 'followed'))
        if not followed:
            res.leaves[walk] = ended
            return
        nsplits = splits + (1 if len(followed) >= 2 else 0)
        for c, v, e in followed:
            for tn, te in states[c].items():
                if te.pred != tn:
                    # the switch event's place: the walk with the step's base (two children
                    # of a split can hold the same switch)
                    res.switches.add((walk + c, te.pred, tn, round(te.switch_cost, 9)))
        for c, v, e in followed:
            visit(walk + c, v, used | {w.edge_key(e)}, states[c], nsplits)

    visit('', root, frozenset(), sigma0, 0)
    return res


def _lost_reason(label, entry, blocked, trace_breaks, states, targets, vertex, st, res, walk,
                 carries, in_a):
    """Why a label that continues on no followed successor ends (SPEC 6.7 with 13, I3 note 3):
    a block of a successor that would carry it; trace_break; then the switches from it into a
    target -- a label of A(v) -- of an admissible successor: `superseded` when one was within the
    budget (a cheaper predecessor took the target), `loss_budget` when one existed only above
    it into a target nobody carries there (needed_budget: the cheapest such switch), otherwise
    `label_lost` (also for an over-budget switch into a target another source carries)."""
    blocks = [b for c, v, e, b in blocked if carries(v)]
    real = [b for b in blocks if b != 'hairpin']
    if real:
        return max(real, key=BLOCK_RANK.get)
    if blocks:
        return 'hairpin'
    if label.name in trace_breaks:
        return 'trace_break'
    within, needed = False, None
    for c, tls in targets.items():
        v = vertex[c]
        if st.switch_on != 'any' and in_a(label, v):
            continue            # not a source of a switch at v
        for tl in tls:
            if tl.name == label.name:
                continue
            cost = st.cost(label, tl)
            if cost == INF:
                continue
            total = entry.loss + cost
            if total <= st.loss_budget + 1e-12:
                within = True
            elif tl.name not in states[c] and (needed is None or total < needed):
                needed = total
    if within:
        return 'superseded'
    if needed is not None:
        res.needed.append(needed)
        res.needed_at[(walk, label.name)] = needed
        return 'loss_budget'
    return 'label_lost'


# -- E: the single-label maximal walks (SPEC 6.9, forbid, exhaustive) ------------------------

def label_walks(w, seed, arm, labels, radius, hairpins='skip'):
    """{(walk, label, reason, branches)}: per label, every maximal walk on which the label holds
    every node, with why it ends there; branches counts the nodes where it went on two or more
    successors. Independent of constrain_walk (no label state, no recurrence)."""
    out = set()
    seed_nodes = w.seed_nodes(seed)
    root = start_node(seed, arm, w.k)
    for l in labels:
        def visit(walk, node, used, branches):
            if len(walk) == radius:
                out.add((walk, l.name, 'max_extension_bp', branches))
                return
            succ = w.successors(node, arm)
            if not succ:
                out.add((walk, l.name, 'dead_end', branches))
                return
            go, blocks = [], []
            for c, v, e in succ:
                if not w.holds(l, v):
                    continue
                b = w.block(node, v, e, used, seed_nodes, hairpins)
                if b is None:
                    go.append((c, v, e))
                else:
                    blocks.append(b)
            if not go:
                real = [b for b in blocks if b != 'hairpin']
                reason = (max(real, key=BLOCK_RANK.get) if real
                          else 'hairpin' if blocks else 'label_lost')
                out.add((walk, l.name, reason, branches))
                return
            # a followed hairpin step never counts toward ambiguity (SPEC 13, note 11), but an
            # ambiguous label's every continuation, the hairpin's too, counts the branch
            plain = [x for x in go if not w.is_hairpin(node, x[1], x[2])]
            nb = branches + (1 if len(plain) >= 2 else 0)
            for c, v, e in go:
                visit(walk + c, v, used | {w.edge_key(e)}, nb)
        visit('', root, frozenset(), 0)
    return out
