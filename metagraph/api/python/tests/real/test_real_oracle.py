"""The search oracle: every cached traversal re-checked by independent /resolve searches
on the same server (realdata.py's oracle client; /resolve shares no code with the
walker -- it maps the query's k-mers and reads their annotation rows).

Per cached cell and seed result, on a deterministic sample of each arm:

  a_present    every sampled leaf walk, spelled with the seed in natural orientation
               (left = natural(left) + seed, right = seed + natural(right); the harness
               speller over the parsed G records, not the library's), is fully in the
               graph. Without bases (output.sequences false) the C continuations are.
  b_claims     every claim() with a non-empty displayed-support interval
               [evidence_from, to_bp): the label's own /resolve runs over the displayed
               molecule cover every k-mer of that stretch (plus the seed-boundary node
               for a stretch from 0 entered by the seed, and for zero-length boundary
               claims at 0). Annotate claims are the recorded-presence stretches.
               Under support trace (mini_refseq) a seed-entered claim is also confirmed
               by /resolve with support trace.
  b_lost       the converse at a plain label_lost end (R end L) and at a silent switch
               (Lw) inside the anchor: /resolve does NOT report the label on the node at
               outward to_bp (k-mer support only), and does on the run's last node.
  b_walks      the agent's top walks (walks(by='support')): Walk.sequence is the
               independent spelling, every labels_full label covers the whole walk
               (boundary + [0, length)), and (constrain) every label of every displayed
               support_profile() run covers that run.
  b_stretch    annotate: the §6.9 oracle E computed HERE from the recorded sets (per
               label, its longest continuous recorded presence from the boundary along a
               displayed walk) is real: /resolve confirms it.
  b_maximal    annotate, library only: claims() holds every label's maximal
               label-consistent route -- the §5.1 union rule through ANY merge parent,
               computed here independently (with keep it is E) -- and nothing beyond it;
               a claim along its anchor's first-parent chain on a displayed walk (route_bp
               0) stays within E (GPT review, finding 3: routes through non-first merge
               parents were dropped).
  c_routes     a run whose displayed support starts at a merge (evidence_from >
               from_bp): cut at D = evidence_from, claims() reports it route_only (no
               displayed bases), and routes(spell=True) spells the label's own route:
               it is a real route (parents chain to the anchor, not the displayed
               chain), it is in the graph, and the label covers its own [from_bp, to_bp)
               on it.
  d_annotate   annotate mode: the recorded label sets (P runs, the root's boundary list)
               at sampled nodes equal the labels /resolve reports for those k-mers; a
               cut list (total > |list|) is a subset with the true count as its total.
  e_direct     label_summary(): for sampled labels with direct_bp = d > 0, a route of
               length d exists on which the label is present on every k-mer (constrain:
               the run's own route from routes(spell=True), the whole molecule incl.
               the seed; annotate: a route through the recorded sets, the boundary node
               and [0, d)) -- the library's route when it has one, else a route
               reconstructed here through any parent.
  e_library    library only: routes(spell=True) spells a route of length direct_bp for
               every label with direct_bp > 0 (expected failures: BUG-ANNOT-ROUTES).

Coordinates (k-mer vs base, design §2.4, §5.1; spec §4.2). A run, a claim or a P run
[a, b) holds the nodes entered by steps a+1..b, i.e. the k-mers whose newest base is
outward base i for i in [a, b); the seed-boundary node (step 0) is the seed's last
k-mer on the right arm and its first on the left. On the molecule M (seed length S,
flank length F, natural orientation):

    right arm  M = seed + flank        outward base i = M[S + i]
               node at i = k-mer start S + i - k + 1;  boundary = S - k
               [a, b)  ->  k-mer starts [S + a - k + 1, S + b - k + 1)
    left arm   M = flank + seed        outward base i = M[F - 1 - i]
               node at i = k-mer start F - 1 - i;      boundary = F
               [a, b)  ->  k-mer starts [F - b, F - a)

/resolve runs are half-open k-mer-start intervals (spec §4.2). Many molecules go in one
/resolve call joined by a single N: every k-mer that spans the separator is absent
from the graph, so runs break exactly there and each molecule's runs are recovered by
offset (checked by test_n_separator_isolates_molecules).

The oracle answers are cached on disk by the harness (keyed by the index identity and
the payload), so a second run needs no server. Without a cache entry and without a
server, a test skips; a live server that serves another index than the cache was
filled from makes the tests skip too (no false verdicts).

    cd api/python/tests/real && python3 -m unittest -v test_real_oracle

Knobs (environment): METAGRAPH_REAL_ORACLE_LEAVES (sampled leaves per arm, default 8),
METAGRAPH_REAL_ORACLE_CLAIMS (sampled claims per arm and class, default 40),
METAGRAPH_REAL_ORACLE_CHUNK_BP (bases per /resolve call, default 60000).
"""

import collections
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import realdata as R  # noqa: E402

from metagraph.traverse import IncompleteRecording  # noqa: E402

N_LEAVES = int(os.environ.get('METAGRAPH_REAL_ORACLE_LEAVES', '8'))
N_CLAIMS = int(os.environ.get('METAGRAPH_REAL_ORACLE_CLAIMS', '40'))
N_ROUTES = 12
N_DIRECT = 12
N_ANNOT_LEAVES = 4
N_ANNOT_POS = 48
CHUNK_BP = int(os.environ.get('METAGRAPH_REAL_ORACLE_CHUNK_BP', '60000'))
CHUNK_LABELS = 800
DISCOVER_MAX = 1000000

# (cell, area) -> reason: a PRODUCT bug this suite found; the test is kept and marked
# unittest.expectedFailure. When the bug is fixed these become unexpected successes:
# remove the entries then.
#
# Fixed (2026-10-02), kept here as the record of what b_maximal and e_library found:
# BUG-ANNOT-CLAIMS (library, ops._carries): in annotate mode with on_reconverge merge,
#   a drop at a child's first base was suppressed as "not maximal" when a sibling on NO
#   displayed walk (its chains end at a non-first merge parent) carried the label, and
#   that sibling's stretch was never claimed either -- uhgg_rand100_00__annotate_merge,
#   left, c:846 had no claim; now [0, 202) on walk 93. 19 cells.
# BUG-ANNOT-ROUTES (library, ops.routes / label_walks in annotate mode): routes() spelled
#   only the first-parent stretches, so a label whose direct_bp (the §5.1 union rule) is
#   reached through a non-first merge parent had no route of length direct_bp; now one
#   witness route per maximal end of the union rule -- same cell, c:846: the 297 bp route
#   0-2-5-35-40-43-95-101-107-110-...-706. 20 cells.
KNOWN_BUGS = {}

AREAS = ('a_present', 'b_claims', 'b_lost', 'b_walks', 'b_stretch', 'b_maximal', 'c_routes',
         'd_annotate', 'e_direct', 'e_library')


# ---------------------------------------------------------------- small helpers

def pick(items, n):
    """A deterministic spread of at most n items (first and last included)."""
    items = list(items)
    if len(items) <= n:
        return items
    if n <= 1:
        return items[:n]
    idx = sorted({round(j * (len(items) - 1) / (n - 1)) for j in range(n)})
    return [items[i] for i in idx]


def covers(runs, lo, hi):
    """Some run [x, y) holds every k-mer start in [lo, hi) (lo < hi)."""
    return any(x <= lo and hi <= y for x, y in runs)


def holds(runs, j):
    return any(x <= j < y for x, y in runs)


def kmer_range(side, S, F, k, a, b, boundary=False):
    """k-mer starts of outward [a, b) on a molecule with seed length S and flank F
    (module docstring); boundary=True adds the seed-boundary node (requires a == 0)."""
    if side == 'right':
        lo, hi = S + a - k + 1, S + b - k + 1
        if boundary:
            lo = S - k
    else:
        lo, hi = F - b, F - a
        if boundary:
            hi = F + 1
    return lo, hi


def node_kmer(side, S, F, k, i):
    """The k-mer start of the node at outward base i (i = -1: the boundary node)."""
    if side == 'right':
        return S + i - k + 1
    return F - 1 - i


def chain_of(arm, sid):
    """root -> sid along parents[0] (written here, not taken from derive)."""
    out = [sid]
    while arm.segments[out[-1]].parents:
        out.append(arm.segments[out[-1]].parents[0])
    return out[::-1]


def spell_chain(g, side, route, upto=None):
    """Natural-orientation flank of a segment route (any parents), cut at |upto|."""
    a = g.arms[side]
    w = ''.join(a.segments[s].walk for s in route)
    if upto is not None:
        w = w[:upto]
    return w if side == 'right' else w[::-1]


def molecule(g, side, flank):
    return g.seed.sequence + flank if side == 'right' else flank + g.seed.sequence


def label_kind(index):
    rec = R.recorded_capabilities(index) or {}
    caps = rec.get('capabilities') or rec
    return 'header' if caps.get('has_coord_to_header') else None


# ---------------------------------------------------------------- the batched oracle

def _slice(runs, off, nk):
    out = []
    for a, b in runs or ():
        a2, b2 = max(a, off), min(b, off + nk)
        if a2 < b2:
            out.append((a2 - off, b2 - off))
    return out


class Oracle:
    """Molecules queried in N-joined /resolve chunks (deduplicated, in insertion order,
    so the payloads -- and the disk cache -- are deterministic)."""

    def __init__(self, index, k, support='kmer'):
        self.index, self.k, self.support = index, k, support
        self.want = collections.OrderedDict()      # molecule -> set of label names
        self.graph, self.runs = {}, {}
        self.done = False

    def need(self, seq, names=()):
        if len(seq) < self.k:
            raise ValueError('a molecule shorter than k')
        self.want.setdefault(seq, set()).update(names)

    def _chunks(self):
        cur, bp, labs = [], 0, set()
        for seq, names in self.want.items():
            add = len(seq) + (1 if cur else 0)
            if cur and (bp + add > CHUNK_BP or len(labs | names) > CHUNK_LABELS):
                yield cur, labs
                cur, bp, labs = [], 0, set()
                add = len(seq)
            cur.append(seq)
            bp += add
            labs |= names
        if cur:
            yield cur, labs

    def run(self, filler=None):
        if self.done:
            return
        for seqs, labs in self._chunks():
            joined = 'N'.join(seqs)
            if not labs and filler:
                labs = {filler}
            if labs:
                out = R.oracle_resolve(self.index, joined, sorted(labs), support=self.support)
            else:
                out = R.oracle_resolve(self.index, joined)
            assert out['num_kmers'] == len(joined) - self.k + 1, (out['num_kmers'], len(joined))
            per = {l['label']: l.get('runs') or [] for l in out.get('labels') or []}
            off = 0
            for seq in seqs:
                nk = len(seq) - self.k + 1
                self.graph[seq] = _slice(out.get('graph_runs'), off, nk)
                self.runs[seq] = {n: _slice(per.get(n), off, nk) for n in self.want[seq]}
                off += len(seq) + 1
        self.done = True

    def present(self, seq):
        return self.graph[seq] == [(0, len(seq) - self.k + 1)]

    def label_runs(self, seq, name):
        return self.runs[seq][name]


class Discover:
    """Annotate oracle: every label /resolve reports per k-mer of each molecule."""

    def __init__(self, index, k, kind=None):
        self.index, self.k, self.kind = index, k, kind
        self.want = collections.OrderedDict()
        self.labels = {}
        self.done = False

    def need(self, seq):
        self.want.setdefault(seq, None)

    def run(self):
        if self.done:
            return
        cur, bp = [], 0
        groups = []
        for seq in self.want:
            add = len(seq) + (1 if cur else 0)
            if cur and bp + add > CHUNK_BP:
                groups.append(cur)
                cur, bp, add = [], 0, len(seq)
            cur.append(seq)
            bp += add
        if cur:
            groups.append(cur)
        disc = {'max_labels': DISCOVER_MAX}
        if self.kind:
            disc['kind'] = self.kind
        for seqs in groups:
            joined = 'N'.join(seqs)
            out = R.oracle_resolve(self.index, joined, discover=disc)
            assert not out.get('labels_truncated'), out.get('labels_truncated')
            off = 0
            for seq in seqs:
                nk = len(seq) - self.k + 1
                got = []
                for l in out.get('labels') or []:
                    r = _slice(l.get('runs'), off, nk)
                    if r:
                        got.append((l['label'], r))
                self.labels[seq] = got
                off += len(seq) + 1
        self.done = True

    def at(self, seq, j):
        return collections.Counter(n for n, r in self.labels[seq] if holds(r, j))


# ---------------------------------------------------------------- the per-cell plan

class Plan:
    """What the oracle checks for one seed result: the sampled molecules, claims, runs,
    labels and nodes, each with its expected k-mer interval. Built from the cached
    graphlet; the oracle queries run once for all areas."""

    def __init__(self, index, row, i, g):
        self.index, self.row, self.i, self.g = index, row, i, g
        self.k = g.k
        self.S = g.seed.length_bp
        self.has_bases = all(s.walk is not None for a in g.arms.values() for s in a.segments)
        self.oracle = Oracle(index, self.k)
        self.trace = Oracle(index, self.k, support='trace') if g.support == 'trace' else None
        self.disc = Discover(index, self.k, label_kind(index)) if g.mode == 'annotate' else None
        self.dup_names = {n for n, c in collections.Counter(
            l.name for l in g.labels).items() if c > 1}
        self.filler = g.labels[0].name if g.labels else None
        self.present = []         # (what, side, molecule)
        self.claims = []          # (claim, side, molecule, lo, hi)
        self.lost = []            # (run, side, molecule, kmer)
        self.routes = []          # dicts, see _plan_routes
        self.annot = []           # (side, molecule, kmer, recorded ids, total, where)
        self.direct = []          # dicts, see _plan_direct
        self.walks = []           # (walk, side, molecule, F, [(label, lo, hi, what)])
        self.stretch = []         # (label, side, molecule, lo, hi, j, segment)
        self.notes = collections.Counter()
        # b_walks / b_stretch query through their own oracle, so that the payloads of
        # the areas above (and their disk cache) do not depend on these
        self.extra = Oracle(index, self.k)
        for side in g.arms:
            self._plan_arm(side)

    # ------------------------------------------------------------ helpers

    def _leaf_molecule(self, side, leaf):
        mol = R.spell_with_seed(self.g, side, leaf)
        return mol, self.g.arms[side].segments[leaf].end_bp

    def _need(self, mol, names=()):
        names = [n for n in names if n not in self.dup_names]
        self.oracle.need(mol, names)

    # ------------------------------------------------------------ one arm

    def _plan_arm(self, side):
        g = self.g
        a = g.arms[side]
        leaves = [s.id for s in a.segments if s.leaf is not None]
        if not self.has_bases:
            # output.sequences false: the continuations are the only bases
            for leaf in pick(leaves, N_LEAVES):
                c = a.segments[leaf].leaf.continuation
                if c is not None and c.sequence and len(c.sequence) >= self.k:
                    self.present.append(('continuation of leaf %d' % leaf, side, c.sequence))
                    self._need(c.sequence)
            return
        deepest = max(leaves, key=lambda s: (a.segments[s].end_bp, -s)) if leaves else None
        chosen = pick(leaves, N_LEAVES)
        if deepest is not None and deepest not in chosen:
            chosen.append(deepest)
        on_chain = {}                 # segment -> the chosen leaf whose chain holds it
        for leaf in chosen:
            for s in chain_of(a, leaf):
                on_chain.setdefault(s, leaf)
        mols = {}

        def mol_for_segment(sid):
            leaf = on_chain.get(sid)
            if leaf is None:
                # a segment off the sampled chains: its own molecule (up to its end)
                key = ('seg', sid)
                if key not in mols:
                    mols[key] = (R.spell_with_seed(g, side, sid), a.segments[sid].end_bp)
                return mols[key]
            key = ('leaf', leaf)
            if key not in mols:
                mols[key] = self._leaf_molecule(side, leaf)
            return mols[key]

        for leaf in chosen:
            mol, F = mol_for_segment(leaf)
            self.present.append(('leaf %d' % leaf, side, mol))
            self._need(mol)

        self._plan_claims(side, a, mol_for_segment)
        if g.mode == 'constrain':
            self._plan_lost(side, a, mol_for_segment)
            self._plan_routes(side, a)
        else:
            self._plan_annotate(side, a, chosen, mol_for_segment)
            self._plan_stretch(side, a)
        self._plan_direct(side, a)
        self._plan_walks(side, a)

    # ------------------------------------------------------------ b: claims

    def _plan_claims(self, side, a, mol_for_segment):
        g = self.g
        try:
            cl = g.claims(side, strict=True)
        except IncompleteRecording:
            # cut lists: a label recorded on every node of a stretch is still truly
            # there (recorded is a subset of the truth), so the stretches stay sound
            cl = g.claims(side, strict=False)
            self.notes['claims_from_cut_lists'] += 1
        special, normal = [], []
        for c in cl:
            if c.kind in ('diverged', 'merged', 'boundary') or c.entered_by == 'switch' \
                    or (c.evidence_from is not None and c.evidence_from > c.from_bp):
                special.append(c)
            else:
                normal.append(c)
        for c in pick(special, N_CLAIMS) + pick(normal, N_CLAIMS):
            if c.label.name in self.dup_names:
                self.notes['claim_label_name_not_unique'] += 1
                continue
            mol, F = mol_for_segment(c.segment)
            if c.kind == 'boundary':
                if c.from_bp != 0:
                    self.notes['boundary_claim_not_at_0'] += 1
                    continue
                j = node_kmer(side, self.S, F, self.k, -1)
                lo, hi = j, j + 1
            else:
                ef = c.evidence_from
                if ef is None or ef >= c.to_bp:
                    self.notes['claim_without_displayed_stretch'] += 1
                    continue
                boundary = ef == 0 and c.from_bp == 0 and c.entered_by == 'seed'
                lo, hi = kmer_range(side, self.S, F, self.k, ef, c.to_bp, boundary)
            assert 0 <= lo < hi <= len(mol) - self.k + 1, (lo, hi, len(mol))
            self.claims.append((c, side, mol, lo, hi))
            self._need(mol, [c.label.name])
            if self.trace is not None and c.entered_by == 'seed':
                self.trace.need(mol, [c.label.name])

    # ------------------------------------------------------------ b: lost

    def _plan_lost(self, side, a, mol_for_segment):
        from metagraph.traverse import derive
        g = self.g
        if g.support != 'kmer':
            return
        lost = [r for r in a.runs if r.end in ('L', 'Lw')
                and r.to_bp < a.segments[r.segment].end_bp]
        for r in pick(lost, N_CLAIMS):
            name = g.labels[r.label].name
            if name in self.dup_names:
                continue
            mol, F = mol_for_segment(r.segment)
            j = node_kmer(side, self.S, F, self.k, r.to_bp)
            assert 0 <= j <= len(mol) - self.k, (j, len(mol))
            # the run's last node carries the label on the displayed path only when its
            # displayed stretch is not empty
            ef = derive.evidence(a, r)[1]
            jb = node_kmer(side, self.S, F, self.k, r.to_bp - 1) if ef < r.to_bp else None
            self.lost.append((r, side, mol, j, jb))
            self._need(mol, [name])

    # ------------------------------------------------------------ c: routes

    def _plan_routes(self, side, a):
        from metagraph.traverse import derive
        g = self.g
        cand = []
        for r in a.runs:
            rf, ef = derive.evidence(a, r)
            if rf > 0 and ef > r.from_bp and r.to_bp > ef:
                cand.append((r, ef))
        if not cand:
            return
        rts_cache = {}
        for r, ef in pick(cand, N_ROUTES):
            name = g.labels[r.label].name
            if name in self.dup_names:
                continue
            if r.label not in rts_cache:
                rts_cache[r.label] = (g.routes({'id': r.label}, side, spell=True),
                                      [x.id for x in a.runs if x.label == r.label])
            rts, ids = rts_cache[r.label]
            route, bases = rts[ids.index(r.id)]
            item = {'run': r, 'side': side, 'cut': ef, 'route': route, 'bases': bases,
                    'name': name}
            if bases is not None and len(bases) == r.to_bp:
                mol = molecule(g, side, bases)
                boundary = r.from_bp == 0 and r.from_label is None
                item['mol'] = mol
                item['range'] = kmer_range(side, self.S, r.to_bp, self.k, r.from_bp,
                                           r.to_bp, boundary)
                self._need(mol, [name])
            self.routes.append(item)

    # ------------------------------------------------------------ d: annotate sets

    def _plan_annotate(self, side, a, chosen, mol_for_segment):
        root = a.segments[0]
        n_pos = 0
        seen = set()
        for leaf in pick(chosen, N_ANNOT_LEAVES):
            mol, F = mol_for_segment(leaf)
            self.disc.need(mol)
            if ('root',) not in seen:
                seen.add(('root',))
                j = node_kmer(side, self.S, F, self.k, -1)
                self.annot.append((side, mol, j, list(root.entry), root.entry_total,
                                   'boundary (root entry)'))
            for sid in chain_of(a, leaf):
                for p in a.segments[sid].presence:
                    for pos in sorted({p.from_bp, p.to_bp - 1}):
                        if (sid, pos) in seen or n_pos >= N_ANNOT_POS:
                            continue
                        seen.add((sid, pos))
                        n_pos += 1
                        j = node_kmer(side, self.S, F, self.k, pos)
                        self.annot.append((side, mol, j, list(p.labels), p.total,
                                           'segment %d P[%d,%d) at %d' % (
                                               sid, p.from_bp, p.to_bp, pos)))

    # ------------------------------------------------------------ e: direct_bp

    def _plan_direct(self, side, a):
        g = self.g
        summ = g.label_summary()
        rows = [(lab, summ[lab.ref][side]['direct_bp']) for lab in g.labels
                if side in summ[lab.ref] and summ[lab.ref][side]['direct_bp'] > 0
                and lab.name not in self.dup_names]
        if not rows:
            return
        top = sorted(rows, key=lambda x: (-x[1], x[0].id))[:N_DIRECT // 2]
        rest = [x for x in rows if x not in top]
        for lab, d in top + pick(rest, N_DIRECT - len(top)):
            item = {'label': lab, 'side': side, 'd': d}
            if g.mode == 'constrain':
                runs = [r for r in a.runs if r.label == lab.id and r.from_label is None
                        and r.from_bp == 0 and r.to_bp == d]
                item['found'] = bool(runs)
                if not runs:
                    self.direct.append(item)
                    continue
                r = runs[0]
                rts = g.routes(lab, side, spell=True)
                ids = [x.id for x in a.runs if x.label == lab.id]
                route, bases = rts[ids.index(r.id)]
                item.update(route=route, bases=bases, library=True)
                mol = molecule(g, side, bases)
                item['mol'] = mol
                item['range'] = (0, len(mol) - self.k + 1)     # seed label: everything
            else:
                lib = [(rt, b) for rt, b in g.routes(lab, side, spell=True)
                       if b is not None and len(b) == d]
                if lib:
                    route, bases = lib[0]
                    item.update(route=route, bases=bases, library=True)
                else:
                    route = annotate_direct_route(a, lab.id, d)
                    item['found'] = route is not None
                    if route is None:
                        self.direct.append(item)
                        continue
                    bases = spell_chain(g, side, route, d)
                    item.update(route=route, bases=bases, library=False)
                    self.notes['direct_route_not_from_library'] += 1
                item['found'] = True
                mol = molecule(g, side, bases)
                item['mol'] = mol
                item['range'] = kmer_range(side, self.S, d, self.k, 0, d, boundary=True)
            self.direct.append(item)
            self._need(item['mol'], [lab.name])

    # ------------------------------------------------------------ b: walks

    def _plan_walks(self, side, a):
        """The agent's top walks (walks(by='support')): the sequence the library spells,
        labels_full over the whole walk, and (constrain) each displayed support run."""
        g = self.g
        for w in g.walks(side, top=N_LEAVES, by='support'):
            mol, F = self._leaf_molecule(side, w.leaf)
            checks = []
            for lab in w.labels_full:
                lo, hi = kmer_range(side, self.S, F, self.k, 0, F, boundary=True)
                checks.append((lab, lo, hi, 'labels_full'))
            if g.mode == 'constrain':
                for sr in g.support_profile(side, w.path_id, kind='displayed'):
                    if sr.to_bp <= sr.from_bp:
                        continue
                    lo, hi = kmer_range(side, self.S, F, self.k, sr.from_bp, sr.to_bp)
                    for lab in sr.labels:
                        checks.append((lab, lo, hi, 'support_profile [%d,%d)'
                                       % (sr.from_bp, sr.to_bp)))
            checks = [c for c in checks if c[0].name not in self.dup_names]
            names = sorted({c[0].name for c in checks})
            if len(names) > 300:
                keep = set(pick(names, 300))
                checks = [c for c in checks if c[0].name in keep]
                self.notes['walk_labels_sampled'] += 1
            self.walks.append((w, side, mol, F, checks))
            self.extra.need(mol, {c[0].name for c in checks})

    # ------------------------------------------------------------ b: stretches (annotate)

    def _plan_stretch(self, side, a):
        """The §6.9 oracle E computed here from the recorded sets alone: per label, its
        longest continuous recorded presence from the boundary along a displayed walk."""
        best = displayed_max(a)
        rows = sorted(((l, j, sid) for l, (j, sid) in best.items() if j > 0),
                      key=lambda x: (-x[1], x[0]))
        on_walk = displayed_segments(a)
        for l, j, sid in rows[:N_CLAIMS // 2] + pick(rows[N_CLAIMS // 2:], N_CLAIMS // 2):
            lab = self.g.labels[l]
            if lab.name in self.dup_names:
                continue
            leaf = leaf_below(a, sid, on_walk)
            mol, F = self._leaf_molecule(side, leaf)
            lo, hi = kmer_range(side, self.S, F, self.k, 0, j, boundary=True)
            self.stretch.append((lab, side, mol, lo, hi, j, sid))
            self.extra.need(mol, [lab.name])

    # ------------------------------------------------------------ run

    def run(self):
        self.oracle.run(filler=self.filler)
        self.extra.run(filler=self.filler)
        if self.trace is not None:
            self.trace.run()
        if self.disc is not None:
            self.disc.run()


def annotate_direct_route(a, lid, d):
    """A segment route (through ANY parents) along which label |lid| is recorded on
    every node from the boundary to outward d: the §5.1 annotate pass with a back
    pointer per segment (written here independently of derive.label_summary)."""
    via = {}
    alive_end = {}
    for s in a.segments:
        if not s.parents:
            alive = lid in s.entry
        else:
            alive = False
            for p in s.parents:
                if alive_end.get(p):
                    alive, via[s.id] = True, p
                    break
        reach = None
        if alive:
            reach = s.end_bp
            for p in s.presence:
                if lid not in p.labels:
                    reach, alive = p.from_bp, False
                    break
            if reach is not None and reach >= d:
                route = [s.id]
                while route[-1] in via:
                    route.append(via[route[-1]])
                return route[::-1]
        alive_end[s.id] = alive
    return None


def displayed_segments(a):
    """on_walk[s]: segment s lies on some leaf's first-parent chain (a displayed walk)."""
    segs = a.segments
    on = [False] * len(segs)
    for s in segs:
        if s.leaf is None:
            continue
        x = s.id
        while not on[x]:
            on[x] = True
            if not segs[x].parents:
                break
            x = segs[x].parents[0]
    return on


def displayed_max(a):
    """The §6.9 oracle E per label, from the recorded sets alone (no ops code): label ->
    (j, segment), j = the longest continuous recorded presence from the seed boundary
    along any displayed walk (first-parent chains of leaves), segment = where it ends."""
    segs = a.segments
    on = displayed_segments(a)
    best = {}
    alive_end = [None] * len(segs)
    for s in segs:
        alive = set(s.entry) if not s.parents else set(alive_end[s.parents[0]])
        reach = {}
        for p in s.presence:
            still = alive.intersection(p.labels)
            for l in alive - still:
                reach[l] = p.from_bp
            alive = still
        for l in alive:
            reach[l] = s.end_bp
        alive_end[s.id] = alive
        if on[s.id]:
            for l, j in reach.items():
                if j > best.get(l, (-1, None))[0]:
                    best[l] = (j, s.id)
    return best


def union_max(a):
    """Per label the longest j of the §5.1 annotate pass (no ops code): the recorded sets
    from the seed boundary through ANY parent of a merge (the union of the parents'
    surviving sets) -- what label_summary's direct_bp and the annotate claims follow."""
    segs = a.segments
    best = {}
    alive_end = [None] * len(segs)
    for s in segs:
        if s.parents:
            alive = set()
            for q in s.parents:
                alive |= alive_end[q]
        else:
            alive = set(s.entry)
        for p in s.presence:
            still = alive.intersection(p.labels)
            for l in alive - still:
                best[l] = max(best.get(l, 0), p.from_bp)
            alive = still
        for l in alive:
            best[l] = max(best.get(l, 0), s.end_bp)
        alive_end[s.id] = alive
    return best


def leaf_below(a, sid, on_walk):
    """A leaf whose first-parent chain holds |sid| (descending through first-parent
    children that lie on a displayed walk; smallest id first)."""
    x = sid
    while a.segments[x].leaf is None:
        x = min(c for c in a.segments[x].children
                if a.segments[c].parents[0] == x and on_walk[c])
    return x


class _PlanCache:
    key = None
    plans = None

    @classmethod
    def get(cls, row):
        key = (row['index'], row['cell'])
        if cls.key != key:
            cls.key, cls.plans = None, None
            cell = R.load_cell(row['index'], row)
            plans = []
            for i in range(cell.n_results()):
                g = cell.graphlet(i)
                plans.append(None if g is None else Plan(row['index'], row, i, g))
            cls.key, cls.plans = key, plans
        return cls.plans


# ---------------------------------------------------------------- the checks

def _require_oracle(index):
    """Skip unless the oracle can answer for the index the cache was filled from: a
    live server must serve that index (equal index_meta_fp); a down server is fine as
    long as the answers are in the oracle cache (a miss then skips)."""
    if R.server_up(index) and not R.live_matches_cache(index):
        raise unittest.SkipTest('%s: the live server serves another index than the '
                                'cache was filled from' % index)


class _OracleBase(unittest.TestCase):
    maxDiff = 2000

    def plans(self, row):
        R.require_cache(row['index'], row['cell'])
        _require_oracle(row['index'])
        plans = [p for p in _PlanCache.get(row) if p is not None]
        for p in plans:
            p.run()
        return plans

    # ------------------------------------------------------------ a

    def check_a_present(self, row):
        n = 0
        for p in self.plans(row):
            for what, side, mol in p.present:
                n += 1
                with self.subTest(result=p.i, arm=side, what=what):
                    self.assertTrue(p.oracle.present(mol),
                                    'not fully in the graph: graph runs %s of %d k-mers'
                                    % (p.oracle.graph[mol], len(mol) - p.k + 1))
        if n == 0:
            self.skipTest('no walk with bases and no continuation >= k')

    # ------------------------------------------------------------ b

    def check_b_claims(self, row):
        n = 0
        for p in self.plans(row):
            for c, side, mol, lo, hi in p.claims:
                n += 1
                runs = p.oracle.label_runs(mol, c.label.name)
                with self.subTest(result=p.i, arm=side, label=c.label.name, kind=c.kind,
                                  run=c.run, interval=(c.evidence_from, c.to_bp)):
                    self.assertTrue(
                        covers(runs, lo, hi),
                        '%s claim [%s, %d) of %s (%s, entered by %s, segment %d): /resolve '
                        'runs %s do not cover k-mers [%d, %d)' % (
                            c.kind, c.evidence_from, c.to_bp, c.label.name, c.end_class,
                            c.entered_by, c.segment, runs[:6], lo, hi))
                    if p.trace is not None and c.entered_by == 'seed':
                        truns = p.trace.label_runs(mol, c.label.name)
                        self.assertTrue(covers(truns, lo, hi),
                                        'trace: /resolve trace runs %s do not cover [%d, %d)'
                                        % (truns[:6], lo, hi))
        if n == 0:
            self.skipTest('no claim with a displayed stretch (%s)' % dict(
                sum((p.notes for p in _PlanCache.get(row) if p), collections.Counter())))

    def check_b_lost(self, row):
        n = 0
        for p in self.plans(row):
            for r, side, mol, j, jb in p.lost:
                n += 1
                name = p.g.labels[r.label].name
                runs = p.oracle.label_runs(mol, name)
                with self.subTest(result=p.i, arm=side, label=name, run=r.id, end=r.end):
                    self.assertFalse(
                        holds(runs, j),
                        'run %d of %s ends %s at %d inside segment %d, but /resolve reports '
                        'the label on that node (k-mer %d; runs %s)' % (
                            r.id, name, r.end, r.to_bp, r.segment, j, runs[:6]))
                    # and it is there on the run's last displayed node
                    if jb is not None:
                        self.assertTrue(holds(runs, jb), ('last node', jb, runs[:6]))
        if n == 0:
            self.skipTest('no label_lost / switch end inside a segment')

    def check_b_walks(self, row):
        n = 0
        for p in self.plans(row):
            for w, side, mol, F, checks in p.walks:
                n += 1
                with self.subTest(result=p.i, arm=side, walk=w.path_id):
                    self.assertEqual(molecule(p.g, side, w.sequence), mol,
                                     'Walk.sequence vs the independent speller')
                    self.assertEqual(w.length_bp, F)
                    self.assertTrue(p.extra.present(mol))
                    bad = []
                    for lab, lo, hi, what in checks:
                        runs = p.extra.label_runs(mol, lab.name)
                        if not covers(runs, lo, hi):
                            bad.append((what, lab.name, (lo, hi), runs[:4]))
                    self.assertEqual(bad[:5], [], '%d of %d label stretches not confirmed'
                                     % (len(bad), len(checks)))
        if n == 0:
            self.skipTest('no walk')

    def check_b_stretch(self, row):
        """The maximal displayed stretches E (computed here, not by the library) are
        real: /resolve reports the label on the boundary node and on every node of
        [0, j) along the displayed walk."""
        n = 0
        for p in self.plans(row):
            for lab, side, mol, lo, hi, j, sid in p.stretch:
                n += 1
                runs = p.extra.label_runs(mol, lab.name)
                with self.subTest(result=p.i, arm=side, label=lab.name, j=j, segment=sid):
                    self.assertTrue(covers(runs, lo, hi), (runs[:6], lo, hi))
        if n == 0:
            self.skipTest('no label recorded beyond the boundary on a displayed walk')

    def check_b_maximal(self, row):
        """claims() in annotate mode holds every label's maximal label-consistent route
        (spec §6.9 with the §5.1 union rule: through ANY parent of a merge; E per label
        under keep), and no claim reaches further: per label, the longest claim == the
        union pass's j (library vs the independent pass here; no server). A claim whose
        route is its anchor's first-parent chain (route_bp 0) on a displayed walk is a
        displayed stretch: within E's j."""
        R.require_cache(row['index'], row['cell'])
        n = 0
        for p in _PlanCache.get(row):
            if p is None:
                continue
            for side, a in p.g.arms.items():
                want = {l: j for l, j in union_max(a).items() if j > 0}
                shown = {l: j for l, (j, _) in displayed_max(a).items()}
                got = {}
                over = []
                for c in p.g.claims(side, strict=False):
                    got[c.label.id] = max(got.get(c.label.id, -1), c.to_bp)
                    if c.route_bp == 0 and c.path_id is not None \
                            and c.to_bp > shown.get(c.label.id, 0):
                        over.append((c.label.ref, c.to_bp, shown.get(c.label.id)))
                n += len(want)
                with self.subTest(result=p.i, arm=side):
                    missing = sorted((p.g.labels[l].ref, j, got.get(l))
                                     for l, j in want.items() if got.get(l) != j)
                    extra = sorted(p.g.labels[l].ref for l in got if l not in want
                                   and got[l] > 0)
                    self.assertEqual(missing, [], '%d of %d labels: (ref, maximal route, '
                                     'longest claim) %s' % (len(missing), len(want),
                                                            missing[:5]))
                    self.assertEqual(extra, [], 'claims beyond any route')
                    self.assertEqual(over[:5], [], 'displayed claims beyond E')
        if n == 0:
            self.skipTest('no label recorded beyond the boundary')

    # ------------------------------------------------------------ c

    def check_c_routes(self, row):
        n = 0
        for p in self.plans(row):
            g = p.g
            for it in p.routes:
                n += 1
                r, side, cut = it['run'], it['side'], it['cut']
                a = g.arms[side]
                with self.subTest(result=p.i, arm=side, run=r.id, label=it['name'], cut=cut):
                    got = [c for c in g.claims(side, labels=[{'id': r.label}],
                                               at_most_bp=cut) if c.run == r.id]
                    self.assertEqual(len(got), 1, 'one claim for the run at the cut')
                    c = got[0]
                    self.assertEqual(c.kind, 'route_only', c)
                    self.assertIsNone(c.evidence_from)
                    self.assertEqual(c.to_bp, cut)
                    self.assertEqual(c.end_class, 'open')
                    self.assertIsNone(c.loss)
                    route = it['route']
                    self.assertEqual(route[-1], r.segment)
                    self.assertFalse(a.segments[route[0]].parents, 'starts at the root')
                    for x, y in zip(route, route[1:]):
                        self.assertIn(x, a.segments[y].parents, (x, y))
                    self.assertNotEqual(route, chain_of(a, r.segment),
                                        'a merge-derived run must not ride the displayed chain')
                    self.assertIsNotNone(it['bases'])
                    self.assertEqual(it['bases'], spell_chain(g, side, route, r.to_bp))
                    mol = it['mol']
                    self.assertTrue(p.oracle.present(mol), 'route not in the graph: %s'
                                    % p.oracle.graph[mol])
                    lo, hi = it['range']
                    runs = p.oracle.label_runs(mol, it['name'])
                    self.assertTrue(covers(runs, lo, hi),
                                    'the label does not cover its own route [%d, %d): runs %s, '
                                    'k-mers [%d, %d)' % (r.from_bp, r.to_bp, runs[:6], lo, hi))
        if n == 0:
            self.skipTest('no run whose displayed support starts at a merge')

    # ------------------------------------------------------------ d

    def check_d_annotate(self, row):
        n = 0
        for p in self.plans(row):
            names = [l.name for l in p.g.labels]
            for side, mol, j, ids, total, where in p.annot:
                n += 1
                truth = p.disc.at(mol, j)
                rec = collections.Counter(names[x] for x in ids)
                with self.subTest(result=p.i, arm=side, where=where):
                    if p.dup_names & set(rec):
                        # two labels share a name: /resolve cannot tell them apart
                        self.assertLessEqual(set(rec), set(truth))
                        continue
                    if total == len(ids):
                        self.assertEqual(set(rec), set(truth),
                                         'recorded %d labels, /resolve %d; only recorded %s, '
                                         'only /resolve %s' % (
                                             len(rec), len(truth),
                                             sorted(set(rec) - set(truth))[:5],
                                             sorted(set(truth) - set(rec))[:5]))
                    else:
                        self.assertLessEqual(set(rec), set(truth), 'cut list not a subset')
                        self.assertEqual(total, len(truth), 'cut list total')
        if n == 0:
            self.skipTest('no annotate node sampled')

    # ------------------------------------------------------------ e

    def check_e_direct(self, row):
        n = 0
        for p in self.plans(row):
            for it in p.direct:
                n += 1
                lab, side, d = it['label'], it['side'], it['d']
                with self.subTest(result=p.i, arm=side, label=lab.name, direct_bp=d):
                    self.assertTrue(it.get('found'), 'no route of length direct_bp=%d' % d)
                    mol = it['mol']
                    self.assertEqual(len(mol), p.S + d)
                    a = p.g.arms[side]
                    route = it['route']
                    self.assertFalse(a.segments[route[0]].parents)
                    for x, y in zip(route, route[1:]):
                        self.assertIn(x, a.segments[y].parents, (x, y))
                    self.assertTrue(p.oracle.present(mol))
                    lo, hi = it['range']
                    runs = p.oracle.label_runs(mol, lab.name)
                    self.assertTrue(covers(runs, lo, hi),
                                    'direct_bp %d of %s: /resolve runs %s on the route do '
                                    'not cover k-mers [%d, %d)' % (d, lab.name, runs[:6],
                                                                   lo, hi))
        if n == 0:
            self.skipTest('no label with direct_bp > 0')

    def check_e_library(self, row):
        """routes(label, arm, spell=True) spells a route of length direct_bp for every
        label with direct_bp > 0 (label_summary: 'continuous presence from the seed
        boundary along some route'; the library must be able to show that route). The
        oracle half of this is e_direct; this is the library half (no server)."""
        R.require_cache(row['index'], row['cell'])
        n = 0
        for p in _PlanCache.get(row):
            if p is None:
                continue
            g = p.g
            summ = g.label_summary()
            for side, a in g.arms.items():
                labs = [lab for lab in g.labels if side in summ[lab.ref]
                        and summ[lab.ref][side]['direct_bp'] > 0]
                if len(labs) > 2000:
                    labs = pick(labs, 2000)
                missing = []
                for lab in labs:
                    d = summ[lab.ref][side]['direct_bp']
                    n += 1
                    if not any(b is not None and len(b) == d
                               for _, b in g.routes(lab, side, spell=True)):
                        missing.append((lab.ref, d))
                with self.subTest(result=p.i, arm=side):
                    self.assertEqual(missing, [], '%d of %d labels have no route of length '
                                     'direct_bp: %s' % (len(missing), len(labs), missing[:5]))
        if n == 0:
            self.skipTest('no label with direct_bp > 0')


# ---------------------------------------------------------------- generation

def _applies(area, row):
    modes = {r.get('label_mode') for r in row.get('results') or [] if 'label_mode' in r}
    tags = set(row.get('tags') or ())
    if not modes:
        return False
    bases = 'sequences_false' not in tags
    if area == 'a_present':
        return True
    if area in ('b_claims', 'b_walks', 'e_direct', 'e_library'):
        return bases
    if area in ('b_stretch', 'b_maximal'):
        return bases and 'annotate' in modes
    if area == 'b_lost':
        return bases and 'constrain' in modes and 'trace' not in tags
    if area == 'c_routes':
        return bases and 'constrain' in modes
    if area == 'd_annotate':
        return bases and 'annotate' in modes
    return False


def _make(area, row):
    def test(self):
        getattr(self, 'check_' + area)(row)
    test.__doc__ = '%s: %s/%s' % (area, row['index'], row['cell'])
    reason = KNOWN_BUGS.get((row['cell'], area))
    if reason:
        test = unittest.expectedFailure(test)
        test.__doc__ += ' [expected failure: %s]' % reason
    return test


def _rows(index):
    try:
        return R.cells(index, heavy=None)
    except Exception:                              # noqa: BLE001 -- no cache: no tests
        return []


def _build_classes():
    out = {}
    for index in R.INDEXES:
        attrs = {}
        for row in _rows(index):
            for area in AREAS:
                if _applies(area, row):
                    attrs['test_%s__%s' % (row['cell'], area)] = _make(area, row)
        cls = type('TestOracle_%s' % index, (_OracleBase,), attrs)
        cls.__doc__ = 'The /resolve oracle on every cached %s cell' % index
        cls.__module__ = __name__
        out[cls.__name__] = cls
    return out


globals().update(_build_classes())


class TestOracleMechanics(unittest.TestCase):
    """The oracle's own assumptions, checked on the live servers."""

    def test_n_separator_isolates_molecules(self):
        """Joining molecules with one N: every k-mer across the separator is absent and
        each molecule's runs are exactly its runs when queried alone."""
        for index in R.INDEXES:
            with self.subTest(index=index):
                try:
                    R.require_catalog(index)
                    _require_oracle(index)
                except unittest.SkipTest:
                    continue
                seeds = [s for s in R.catalog_seeds(index) if s['kind'] != 'batch'][:3]
                k = R._k(index)
                labels = sorted({l for s in seeds for l in s['labels'][:2]})
                try:
                    alone = {s['name']: R.oracle_resolve(index, s['sequence'], labels)
                             for s in seeds}
                    joined_seq = 'N'.join(s['sequence'] for s in seeds)
                    joined = R.oracle_resolve(index, joined_seq, labels)
                except unittest.SkipTest:
                    continue
                off = 0
                for s in seeds:
                    nk = len(s['sequence']) - k + 1
                    self.assertEqual(_slice(joined['graph_runs'], off, nk),
                                     [tuple(x) for x in alone[s['name']]['graph_runs']])
                    per = {l['label']: l['runs'] for l in joined['labels']}
                    for l in alone[s['name']]['labels']:
                        self.assertEqual(_slice(per.get(l['label']), off, nk),
                                         [tuple(x) for x in l['runs']], (s['name'], l['label']))
                    # the k k-mers holding the separator are never in the graph
                    sep = off + len(s['sequence'])
                    if sep < len(joined_seq):
                        for a, b in joined['graph_runs']:
                            self.assertFalse(a < sep + 1 and b > sep - k + 1, (a, b, sep))
                    off += len(s['sequence']) + 1

    def test_kmer_coordinates_on_a_known_walk(self):
        """The coordinate rule on a hand-checked case: the seed's own k-mers are the
        boundary node and its predecessors, and outward base 0 is the k-mer right after
        the boundary on the right arm, right before it on the left."""
        S, F, k = 50, 10, 31
        self.assertEqual(kmer_range('right', S, F, k, 0, 1), (20, 21))
        self.assertEqual(node_kmer('right', S, F, k, -1), S - k)
        self.assertEqual(kmer_range('right', S, F, k, 0, F, boundary=True), (S - k, S + F - k + 1))
        self.assertEqual(S + F - k + 1, S + F - k + 1)             # the molecule's last k-mer
        self.assertEqual(kmer_range('left', S, F, k, 0, 1), (F - 1, F))
        self.assertEqual(node_kmer('left', S, F, k, -1), F)
        self.assertEqual(kmer_range('left', S, F, k, 0, F, boundary=True), (0, F + 1))

    def test_every_index_has_oracle_tests(self):
        R.require_cache()
        mod = sys.modules[__name__]
        for index in R.INDEXES:
            if not R.cells(index):
                continue
            cls = getattr(mod, 'TestOracle_%s' % index)
            names = [m for m in dir(cls) if m.startswith('test_')]
            for area in AREAS:
                self.assertTrue(any(m.endswith('__' + area) for m in names), (index, area))


if __name__ == '__main__':
    unittest.main()
