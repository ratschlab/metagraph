"""Stitching a continuation to the walk it continues (docs/DESIGN-graphlet-merge.md §3,
phase 1a-S). EXPERIMENTAL: the API may change until `exact` is switched on.

A continuation (next_request(), deepen(), traverse_continue) is a new traversal seeded with
a walk's continuation. stitch(parent, child) joins it to that walk, on one pair, as a
StitchedGraphlet: its parts held unchanged (by reference, never copied or modified), every
answer computed from them, positions outward from the root seed's boundary.

A continuation restarts what one uninterrupted walk would carry on: branch state, loss
budget, edge history, seed re-entry, merges with other walks and (under trace) the
followed occurrence. So each junction and each stitched walk states, as named causes, why
its walks may differ from one uninterrupted walk of the parent's request (§3.4). The rules
that would call a junction `exact` (W1-W8, the whitelist) are computed, and every failing
condition is named; but `exact` is not switched on before a seeded randomized comparison
against uninterrupted traversals passes (OQ19). Until then every verdict is `qualified`,
with the one cause `exact_gated` where the rules found no difference. Completeness,
evidence and outcomes carry every cause (§3.5): a stitched walk is certified at most to its
first junction, no walk that reaches into a child is complete, and no answer that depends
on a junction reads `exact`. An answer that depends on no junction (stitch(g) without a
child, a root walk no child continues) answers as the root does.

What is not offered (CombineRefused 'not_offered'): dump(), to_json(), compare(),
subgraph(), to_gfa(), label_walks(), routes(), label_summary(), splits(); claims(),
support_profile() and support_changes() (phase 1b), save() (phase 1c).

The local rule check compares (k+1)-mers exactly; the engine compares 128-bit hashes where a
(k+1)-mer does not fit 64 bits (DNA k >= 32, protein k > 11), so a hash collision that
blocks a step in the engine is not seen here.
"""

import json
from dataclasses import dataclass, field
from typing import Optional

from . import budget as _B
from . import composite as _K
from . import derive
from . import export as _export
from . import ops as _ops
from ._codec import REASON
from .budget import (STR, W_ELEM, W_ROW, W_STEP, LocalBudgetExceeded, Partial, record_bytes,
                     set_bytes)
from .client import answered_pair
from .composite import CombineRefused, NotCombinable
from .model import (ARM_SIDES, AmbiguousLabel, BadSelector, Continuation, Graphlet, Label,
                    UnknownLabel)

__all__ = ['stitch', 'StitchedGraphlet', 'Part', 'Junction', 'StitchedWalk', 'CAUSES',
           'JUNCTION_CAUSES', 'WALK_CAUSES']


# ------------------------------------------------------------------ refusal texts (B.1)

_T_NOT_COMBINABLE = 'stitch() combines graphlets (Graphlet or StitchedGraphlet), not %s'
_T_VIEW = 'stitch() reads whole retrievals, not views: stitch view.backing instead'
_T_DIFFERENT_INDEX = (
    'the child was retrieved from another index than its parent (%s): a continuation runs '
    'on its parent\'s index, so a child of another index continues none of its walks; '
    'union_pairs() combines one seed\'s graphlets across pairs')
_T_INDEX_UNVERIFIABLE = (
    'the child states %s index manifest digest (index_fp) and its parent %s: whether they '
    'are the same index cannot be verified (not proven different), so the child is not '
    'stitched automatically; pass allow_unverified_index=True to stitch it anyway (the '
    'result then states the identity unverified)')
_T_DIFFERENT_PAIR = ('the child was answered on pair %s, its parent on pair %s: a '
                     'continuation runs on its parent\'s pair')
_T_DIFFERENT_GRAPH = ('the child\'s graph is not its parent\'s (%s %r vs %r): its walks '
                      'cannot continue the parent\'s')
_T_DIFFERENT_MODE = ('the child was walked in %s mode with %s support, its parent in %s mode '
                     'with %s support: its walks answer another question')
_T_NOT_ONE_ARM = ('the child walked both arms: a continuation walks only the arm it '
                  'continues (next_request() sets direction to it)')
_T_WRONG_ARM = 'the child walked the %s arm, not the %s arm named by arm='
_T_NO_CONTINUATION = ('walk %s of the %s arm has no continuation (it ended for a semantic '
                      'reason, not at the radius or a cap)')
_T_NOT_A_CONTINUATION_OF = (
    'the child\'s seed is not the continuation of walk %s of the %s arm: a stitched child is '
    'seeded with that walk\'s continuation (continuation(), next_request())')
_T_NOT_A_CONTINUATION = (
    'the child\'s seed is the continuation of no walk of the %s arm: a stitched child is '
    'seeded with a walk\'s continuation (continuation(), next_request())')
_T_AMBIGUOUS_WALK = ('the child\'s seed is the continuation of walks %s of the %s arm: name '
                     'the one it continues with walk=')
_T_ALREADY = ('walk %s of the %s arm is already continued by another child (part %s, '
              'digest %s): one child per walk')
_T_ALREADY_REQUEST = (
    'walk %s of the %s arm is continued by this body from another request (%s differs): '
    'equal bodies are not equal retrievals, and the two would continue differently')
_T_FOREIGN = ('the child seeds %s, not alive at the end of walk %s of the %s arm: their runs '
              'would start new lineages at %d bp, not continue the walk\'s')
_T_NOT_OFFERED = {
    'dump': ('a stitched graphlet is not an MGT document: save() writes its parts unchanged '
             'in a bundle, and each part dumps on its own (parts[i].graphlet.dump())'),
    'to_json': ('to_json() is detail "full" of one server result, and a stitched graphlet '
                'has none: export its parts'),
    'compare': 'compare() reads one retrieval per side: compare the parts',
}
_T_NOT_OFFERED_ANY = '%s() is not offered on a stitched graphlet'
_T_LABEL_ID = ('label ids are local to each part of a stitched graphlet: select by ref or '
               'name')

_SHOWN = 8          # labels a refusal names (then "and N more"), as next_request's notes


def _label_list(labels):
    shown = [json.dumps({'name': l.name, 'ref': l.ref}, ensure_ascii=False)
             for l in labels[:_SHOWN]]
    more = len(labels) - _SHOWN
    return ', '.join(shown) + (' and %d more' % more if more > 0 else '')


# ------------------------------------------------------------------ causes (B.2)

# the order of Appendix B.2's table: a junction's and a walk's causes are listed in it
CAUSES = ('trace_restart', 'strategy_differs', 'extra_foreign', 'switching',
          'branches_unequal', 'label_count_knob', 'sibling_merge_unknown', 'not_continued',
          'seed_dropped', 'no_bases', 'rule_violation', 'candidate_changed', 'exact_gated',
          'decisions_unrecorded', 'seed_reentry', 'split_limit_restarted')
JUNCTION_CAUSES = CAUSES[:10] + ('exact_gated',)
WALK_CAUSES = ('rule_violation', 'candidate_changed', 'decisions_unrecorded', 'seed_reentry',
               'split_limit_restarted', 'exact_gated')
_RANK = {c: i for i, c in enumerate(CAUSES)}

_GATED = ('exact_gated', 'W1-W8 hold here, but exact is not switched on: the rules found no '
          'difference from one uninterrupted walk and have not yet passed the randomized '
          'comparison against uninterrupted traversals')
_GATED_WALK = ('exact_gated', 'W1-W8 hold for this walk and its junctions, but exact is not '
               'switched on: the rules found no difference from one uninterrupted walk and '
               'have not yet passed the randomized comparison', None)

# W2: the strategy fields a continuation may change (judged by W3-W5 where they matter)
_W2_ALLOWED = ('direction', 'clamped', 'output', 'annotation',
               ('bounds', 'max_extension_bp'), ('bounds', 'max_steps'), ('bounds', 'max_paths'),
               ('bounds', 'max_output_bp'), ('bounds', 'time_budget_ms'),
               ('bounds', 'max_memory_mb'), ('bounds', 'max_work_units'),
               ('labels', 'max_labels_per_node'), ('labels', 'loss_budget'),
               ('labels', 'extra'), ('branching', 'max_label_branches'))


# ------------------------------------------------------------------ data classes (C.2)

@dataclass(frozen=True)
class Part:
    """One retrieval of a stitched graphlet, held unchanged. |key|: () for the root,
    (('right', 3),) for the child continuing root walk 3, (('right', 3), ('right', 0)) for
    the child continuing that child's walk 0. |offset_bp|: where the part's seed boundary
    lies, outward from the root seed's boundary (the continued walk's end). |digest|: the
    sha-256 of its body (composite.part_digest)."""
    key: tuple
    graphlet: Graphlet = field(repr=False, compare=False)
    offset_bp: int
    digest: str


@dataclass(frozen=True)
class Junction:
    """Where a child continues a walk: |arm|, the continued |walk| (a walk key), |at_bp| the
    walk's end (outward from the root seed), |overlap_bp| the child's seed length (the
    continuation; never spelled twice), |part| the child's part key, the |verdict| and its
    |causes| ((code, sentence), in Appendix B.2's order), and the labels alive at the walk's
    end the child |continued| (seeded) and those it did not (|not_continued|; constrain)."""
    arm: str
    walk: tuple
    at_bp: int
    overlap_bp: int
    part: tuple
    verdict: str
    causes: tuple
    continued: tuple
    not_continued: tuple


@dataclass(slots=True)
class StitchedWalk:
    """A walk of a stitched graphlet: the root's walks, with every continued walk replaced by
    its child's walks. |key| (3,) is root walk 3, (3, 0) walk 0 of the child continuing it.
    |pieces|: (part key, path id, from_bp, to_bp) per part, outward from the root seed. The
    end fields (|path_reason|, |end_reasons|, |n_alive|) are the last piece's;
    |labels_full| the labels with displayed support on every piece and continued at every
    junction. |junctions|: indexes into StitchedGraphlet.junctions. |verdict|: 'qualified'
    for a walk a junction joins (every verdict in phase 1a-S); verdict is None for a walk no
    junction joins, and its evidence is the root's. |causes|: the walk's own (code,
    sentence, at_bp) of §3.4 (exact_gated: none found). |certified_to_bp|: how far
    the walk is certified (at most its first junction in phase 1a-S); |complete|: certified
    whole, its end certified by its part, no junction and no walk-scope cause on it."""
    key: tuple
    arm: str
    length_bp: int
    sequence: Optional[str]
    path_reason: Optional[str]
    end_reasons: dict
    labels_full: list
    n_alive: int
    complete: bool
    verdict: Optional[str]
    certified_to_bp: int
    beyond_certified_bp: int
    junctions: tuple
    pieces: tuple
    causes: tuple


# ------------------------------------------------------------------ helpers

def _arm_name(arm):
    if arm is None:
        return None
    s = {'l': 'left', 'r': 'right'}.get(arm, arm)
    if s not in ARM_SIDES:
        raise ValueError("arm is 'left', 'right', 'l' or 'r', not %r" % (arm,))
    return s


def walk_name(key):
    """A walk key as text: 3 for (3,), 3.0 for (3, 0)."""
    return '.'.join(str(x) for x in key)


def part_name(key):
    """A part key as text: 'root', 'right:3', 'right:3.0'."""
    if not key:
        return 'root'
    return '%s:%s' % (key[0][0], '.'.join(str(w) for _, w in key))


def _walk_prefix(key):
    """The walk-key prefix of a part's walks: (3,) for the child of root walk 3."""
    return tuple(w for _, w in key)


def _rc(s):
    return _K.reverse_complement(s)


def _kmers(seq, k, regime):
    """A seed's nodes as k-mers, natural orientation (both orientations outside basic)."""
    out = {seq[i:i + k] for i in range(len(seq) - k + 1)}
    if regime != 'basic':
        out |= {_rc(x) for x in out}
    return out


def _seed_walking(g, side):
    s = g.seed.sequence
    return s if side == 'right' else s[::-1]


def _tails(g, side):
    """Per segment the last k bases, walking order, of the seed and the chain through its
    end (any parent's: the parents of a merge end in one node)."""
    a = g.arms[side]
    k = g.k
    seed = _seed_walking(g, side)[-k:]
    out = [None] * len(a.segments)
    for s in a.segments:              # creation order: parents first
        before = seed if not s.parents else out[s.parents[0]]
        out[s.id] = (before + s.walk)[-k:]
    return out


def _before(tails, g, side, s):
    """The k bases before segment |s|'s first base (walking order)."""
    if not s.parents:
        return _seed_walking(g, side)[-g.k:]
    return tails[s.parents[0]]


def _natural(side, w):
    return w if side == 'right' else w[::-1]


def _ancestors(arm, seg):
    """|seg| and its DAG ancestors through any parent, ascending ids."""
    segs = arm.segments
    seen = {seg}
    stack = [seg]
    while stack:
        for p in segs[stack.pop()].parents:
            if p not in seen:
                seen.add(p)
                stack.append(p)
    return sorted(seen)


def _split_counts(arm):
    """The engine's split count per segment (walker.cpp: 0 at the root, one more on each
    child of a split, the maximum over the parents at a merge)."""
    segs = arm.segments
    out = [0] * len(segs)
    for s in segs:
        if s.parents:
            out[s.id] = max(out[p] + (1 if segs[p].split is not None else 0)
                            for p in s.parents)
    return out


def _strategy(g):
    return (g.envelope or {}).get('strategy') or {}


def _knob(st, block, name, default=None):
    v = (st.get(block) or {}).get(name, default)
    return v


def _without_allowed(st, other):
    """|st| without W2's allowed fields (max_live_paths too where both stop on overflow)."""
    out = json.loads(json.dumps(st))
    both_stop = _knob(st, 'frontier', 'on_overflow') == 'stop' and \
        _knob(other, 'frontier', 'on_overflow') == 'stop'
    for f in _W2_ALLOWED + ((('bounds', 'max_live_paths'),) if both_stop else ()):
        if isinstance(f, str):
            out.pop(f, None)
        elif isinstance(out.get(f[0]), dict):
            out[f[0]].pop(f[1], None)
    return out


def _is_graphlet(x):
    return isinstance(x, (Graphlet, StitchedGraphlet)) or _K.is_view(x)


def _type_name(x):
    return type(x).__name__


# ------------------------------------------------------------------ stitch()

def stitch(parent, child=None, *, arm=None, walk=None, allow_unverified_index=False,
           budget=None):
    """Join |child|, a continuation of one walk of |parent|, to that walk (§3.1) -> a new
    StitchedGraphlet (the parts held by reference, never copied or modified).

    |parent|: a Graphlet with its envelope, or a StitchedGraphlet; |child|: a Graphlet or a
    StitchedGraphlet rooted at the continuation; None gives the parent alone as a
    StitchedGraphlet with no junction (it answers every offered query as the graphlet).
    |arm|: 'left', 'right', 'l' or 'r' (default: the child's one arm). |walk|: the continued
    walk, a path id on a Graphlet parent, a walk key (tuple) on a StitchedGraphlet parent
    (an int is the root walk); default: the one walk of the arm whose continuation is the
    child's seed (several: ambiguous_walk). One junction per call.

    Refusals (CombineRefused with .code, in the order checked; Appendix B.1 has the texts):
    not_combinable (also a TypeError), view_unsupported, MissingEnvelope, different_index,
    index_unverifiable (a manifest digest on one side only; allow_unverified_index=True
    stitches it and states the identity unverified), different_pair, different_graph,
    different_mode, not_one_arm, wrong_arm, no_continuation, not_a_continuation,
    ambiguous_walk, already_continued (the same retrieval again at that walk is a no-op:
    the parent is returned), foreign_seed_labels. budget= (local limits): the junction's
    rule check is charged; a stop raises with no partial."""
    if not isinstance(allow_unverified_index, bool):
        raise ValueError('allow_unverified_index is True or False, not %r'
                         % (allow_unverified_index,))
    for x in (parent,) if child is None else (parent, child):
        if not _is_graphlet(x):
            raise NotCombinable(_T_NOT_COMBINABLE % _type_name(x))
    for x in (parent,) if child is None else (parent, child):
        _K.require_whole(x, _T_VIEW)
    for x in (parent,) if child is None else (parent, child):
        if isinstance(x, Graphlet):
            _K.require_envelope(x, 'stitch()')
    base = parent if isinstance(parent, StitchedGraphlet) else \
        StitchedGraphlet._root_of(parent, allow_unverified_index)
    if child is None:
        return base
    return _K.run_scoped(budget, 'stitch',
                         lambda b: base._stitched(child, arm, walk, allow_unverified_index,
                                                  b))


# ------------------------------------------------------------------ the graphlet

class StitchedGraphlet:
    """A continuation joined to the walk it continues (§3; experimental). Built by stitch().
    Its parts are held unchanged; every answer is computed from them, positions outward
    from the root seed's boundary, walks keyed by tuples of walk ids along the chain."""

    def __init__(self, parts, identity, allow_unverified_index):
        # parts: [Part] in canonical order (root first, then by key)
        self._parts = tuple(parts)
        self._by_key = {p.key: p for p in self._parts}
        self.identity = identity
        self._allow = allow_unverified_index
        self._cache = {}
        root = self._parts[0].graphlet
        self.k, self.regime, self.alphabet = root.k, root.regime, root.alphabet
        self.mode, self.support = root.mode, root.support
        self.seed = root.seed
        self._labels, self._label_of = self._union_labels()
        self._junctions = tuple(self._make_junction(p) for p in self._parts[1:])
        self._junctions = tuple(sorted(self._junctions, key=lambda j: (j.arm, j.walk)))
        self.digest = _K.combined_digest([[[list(s) for s in p.key], p.digest]
                                          for p in self._parts])

    # ---------------------------------------------------------------- construction

    @classmethod
    def _root_of(cls, g, allow_unverified_index=False):
        identity = _identity_alone(g)
        return cls([Part((), g, 0, _K.part_digest(g))], identity, allow_unverified_index)

    def _stitched(self, child, arm, walk, allow, b):
        """The checks of §3.1 in order, then the new graphlet (or self: a no-op)."""
        croot = child.root if isinstance(child, StitchedGraphlet) else child
        proot = self.root
        identity = self._check_index(croot, child.labels, allow)
        if not identity['verified']:
            self._check_pair(proot, croot)
        for name in ('k', 'regime'):
            if getattr(proot, name) != getattr(croot, name):
                raise CombineRefused('different_graph', _T_DIFFERENT_GRAPH % (
                    name, getattr(croot, name), getattr(proot, name)))
        if _K.alphabet_key(proot) != _K.alphabet_key(croot):
            raise CombineRefused('different_graph', _T_DIFFERENT_GRAPH % (
                'alphabet', croot.alphabet, proot.alphabet))
        if (proot.mode, proot.support) != (croot.mode, croot.support):
            raise CombineRefused('different_mode', _T_DIFFERENT_MODE % (
                croot.mode, croot.support, proot.mode, proot.support))
        if len(croot.arms) != 1:
            raise CombineRefused('not_one_arm', _T_NOT_ONE_ARM)
        side = next(iter(croot.arms))
        want = _arm_name(arm)
        if want is not None and want != side:
            raise CombineRefused('wrong_arm', _T_WRONG_ARM % (side, want))
        cseq = _K.case_mapped(croot.seed.sequence, croot.alphabet)
        part, pid, key = self._continued_walk(side, walk, cseq, b)
        pg = part.graphlet
        step = (side, pid)
        ckey = part.key + (step,)
        have = self._by_key.get(ckey)
        if have is not None:
            if isinstance(child, StitchedGraphlet):
                same = [p.digest for p in self.parts if p.key[:len(ckey)] == ckey] == \
                    [p.digest for p in child.parts]
            else:
                same = have.digest == _K.part_digest(croot)
            if same:
                ra, rb = (_K.retrieval_identity(have.graphlet),
                          _K.retrieval_identity(croot))
                diff = _K.first_difference(ra, rb)
                if diff is None:
                    return self
                raise CombineRefused('already_continued', _T_ALREADY_REQUEST % (
                    walk_name(key), side, diff))
            raise CombineRefused('already_continued', _T_ALREADY % (
                walk_name(key), side, part_name(ckey), have.digest))
        leaf = derive.paths(pg.arms[side])[pid].leaf
        L = pg.arms[side].segments[leaf].end_bp
        if croot.mode == 'constrain':
            alive = {pg.labels[l].ref for l in pg.arms[side].segments[leaf].end}
            foreign = [l for l in croot.labels[:croot.seed.num_seed_labels]
                       if l.ref not in alive]
            if foreign:
                raise CombineRefused('foreign_seed_labels', _T_FOREIGN % (
                    _label_list(foreign), walk_name(key), side, part.offset_bp + L))
        offset = part.offset_bp + L
        if isinstance(child, StitchedGraphlet):
            new = [Part(ckey + p.key, p.graphlet, offset + p.offset_bp, p.digest)
                   for p in child.parts]
        else:
            new = [Part(ckey, child, offset, _K.part_digest(child))]
        parts = sorted(self._parts + tuple(new), key=lambda p: p.key)
        identity = _weaker_identity(self.identity, identity)
        if isinstance(child, StitchedGraphlet):
            identity = _weaker_identity(identity, child.identity)
        if b is not None:
            # the junction's rule check (W8) is computed on demand; the new graphlet's own
            # index (parts, keys, junctions) is charged here
            b.charge(W_ROW * len(parts), len(parts) * (record_bytes(Part) + 200))
        return StitchedGraphlet(parts, identity, self._allow or allow)

    def _check_index(self, croot, labels, allow):
        """different_index, then index_unverifiable -> the identity statement."""
        names = self._label_of
        clause = None
        key, values = _K.index_check(self.root, croot)
        if key in ('meta_fp', 'fp', 'ns', 'release'):
            clause = _K.different_index_clause(key, values)
        else:
            for l in labels:
                if l.ref in names and names[l.ref].name != l.name:
                    clause = _K.INDEX_CLAUSES['names'] % (l.ref, names[l.ref].name, l.name)
                    break
            else:
                if key == 'names':
                    clause = _K.different_index_clause(key, values)
        if clause is not None:
            raise CombineRefused('different_index', _T_DIFFERENT_INDEX % clause)
        fp_p, fp_c = self.root.index_fp, croot.index_fp
        if fp_p and fp_c:
            return {'verified': True, 'index_fp': fp_p}
        if fp_p or fp_c:
            why = _T_INDEX_UNVERIFIABLE % ('an' if fp_c else 'no', 'one' if fp_p else 'none')
            if not allow:
                raise CombineRefused('index_unverifiable', why)
            return {'verified': False, 'accepted': 'allow_unverified_index',
                    'note': 'the child states %s index manifest digest (index_fp) and its '
                            'parent %s: whether they are the same index cannot be verified '
                            '(not proven different); stitched because the caller passed '
                            'allow_unverified_index' % ('an' if fp_c else 'no',
                                                         'one' if fp_p else 'none')}
        return {'verified': False,
                'note': 'neither the parent nor the child states an index manifest digest '
                        '(index_fp): that they are the same index cannot be verified'}

    def _check_pair(self, proot, croot):
        pp = answered_pair(proot.envelope)
        cp = answered_pair(croot.envelope)
        if pp is not None and cp is not None and pp != cp:
            def show(p):
                return json.dumps({'graph': p[0], 'graph_path': p[1], 'annotation_path': p[2]},
                                  separators=(',', ':'), ensure_ascii=False)
            raise CombineRefused('different_pair', _T_DIFFERENT_PAIR % (show(cp), show(pp)))

    def _continued_walk(self, side, walk, cseq, b):
        """-> (part, path id, walk key) of the walk the child continues (no_continuation,
        not_a_continuation, ambiguous_walk)."""
        if walk is not None:
            key = self._walk_key(walk)
            part, pid = self._part_of_walk(side, key)
            g = part.graphlet
            a = g.arms[side]
            leaf = derive.paths(a)[pid].leaf
            if a.segments[leaf].leaf.continuation is None:
                raise CombineRefused('no_continuation', _T_NO_CONTINUATION % (
                    walk_name(key), side))
            seq = derive.continuation_sequence(g, a, leaf)
            if seq is None or _K.case_mapped(seq, g.alphabet) != cseq:
                raise CombineRefused('not_a_continuation', _T_NOT_A_CONTINUATION_OF % (
                    walk_name(key), side))
            return part, pid, key
        found = []
        for part in self._parts:
            g = part.graphlet
            if side not in g.arms:
                continue
            a = g.arms[side]
            if b is not None:
                b.charge(W_ELEM * len(a.segments))
            for p in derive.paths(a):
                seg = a.segments[p.leaf]
                if seg.leaf.continuation is None:
                    continue
                seq = derive.continuation_sequence(g, a, p.leaf)
                if seq is not None and _K.case_mapped(seq, g.alphabet) == cseq:
                    found.append((part, p.id, _walk_prefix(part.key) + (p.id,)))
        if not found:
            raise CombineRefused('not_a_continuation', _T_NOT_A_CONTINUATION % side)
        if len(found) > 1:
            raise CombineRefused('ambiguous_walk', _T_AMBIGUOUS_WALK % (
                ', '.join(walk_name(k) for _, _, k in found), side))
        return found[0]

    def _walk_key(self, walk):
        if isinstance(walk, bool):
            raise BadSelector('a walk is a path id or a walk key (tuple of ints): %r'
                              % (walk,))
        if isinstance(walk, int):
            return (walk,)
        if isinstance(walk, (tuple, list)) and walk and all(
                isinstance(x, int) and not isinstance(x, bool) for x in walk):
            return tuple(walk)
        raise BadSelector('a walk is a path id or a walk key (tuple of ints): %r' % (walk,))

    def _part_of_walk(self, side, key):
        """(part, path id) of walk key |key| on arm |side|, the walk itself (continued or
        not): the part reached through key[:-1]."""
        pkey = tuple((side, w) for w in key[:-1])
        part = self._by_key.get(pkey)
        if part is None:
            raise IndexError('no walk %s on the %s arm: walk %s is not continued here'
                             % (walk_name(key), side, walk_name(key[:-1])))
        g = part.graphlet
        if side not in g.arms:
            raise KeyError('arm %r was not requested in this retrieval' % side)
        n = len(derive.paths(g.arms[side]))
        if not 0 <= key[-1] < n:
            raise IndexError('the %s arm of part %s has %d walks, no walk %d'
                             % (side, part_name(pkey), n, key[-1]))
        return part, key[-1]

    def _union_labels(self):
        """The labels joined by ref: the root's order, then labels new in later parts,
        renumbered in that order (a label id is local to each part)."""
        out, by_ref = [], {}
        for p in self._parts:
            for l in p.graphlet.labels:
                if l.ref not in by_ref:
                    nl = Label(len(out), l.kind, l.column, l.seq_id, l.name)
                    out.append(nl)
                    by_ref[l.ref] = nl
        return tuple(out), by_ref

    # ---------------------------------------------------------------- properties

    @property
    def parts(self):
        return self._parts

    @property
    def root(self):
        return self._parts[0].graphlet

    @property
    def junctions(self):
        return self._junctions

    @property
    def labels(self):
        return self._labels

    def arms(self):
        return [s for s in ARM_SIDES if s in self.root.arms]

    def arm_of(self, arm=None):
        if arm is None:
            sides = self.arms()
            if len(sides) != 1:
                raise ValueError('both arms are present: name one (left or right)')
            return sides[0]
        s = _arm_name(arm)
        if s not in self.root.arms:
            raise KeyError('arm %r was not requested in this retrieval' % arm)
        return s

    def stitch(self, child, *, arm=None, walk=None, allow_unverified_index=False,
               budget=None):
        """stitch(self, child, ...): another junction, a new StitchedGraphlet."""
        return stitch(self, child, arm=arm, walk=walk,
                      allow_unverified_index=allow_unverified_index, budget=budget)

    def label(self, selector):
        """A label of the union by ref or name: 'c:0', 'acc1', {'ref': ...}, {'name': ...}
        (AmbiguousLabel / UnknownLabel as Graphlet.label). {'id': n} is refused: ids are
        local to each part."""
        if isinstance(selector, Label):
            selector = {'ref': selector.ref}
        if isinstance(selector, dict):
            if set(selector) == {'id'}:
                raise BadSelector(_T_LABEL_ID)
            if len(selector) != 1 or not set(selector) <= {'ref', 'name'}:
                raise BadSelector('a label selector is a ref or a name: %r' % (selector,))
            (kind, value), = selector.items()
            if kind == 'ref':
                if value in self._label_of:
                    return self._label_of[value]
                raise UnknownLabel('no label with ref %r' % (value,))
            got = [l for l in self._labels if l.name == value]
        elif isinstance(selector, str):
            if selector in self._label_of:
                byname = [l for l in self._labels if l.name == selector]
                if not byname or byname == [self._label_of[selector]]:
                    return self._label_of[selector]
                raise AmbiguousLabel('%r is a ref and a name: select {"ref": ...} or '
                                     '{"name": ...}' % selector)
            got = [l for l in self._labels if l.name == selector]
        else:
            raise BadSelector('a label selector is a ref or a name: %r' % (selector,))
        if len(got) == 1:
            return got[0]
        if not got:
            raise UnknownLabel('no label named %r' % (selector,))
        raise AmbiguousLabel('%d labels are named %r: select one by ref' % (len(got), selector))

    def _union(self, g, ids):
        return [self._label_of[g.labels[i].ref] for i in ids]

    # ---------------------------------------------------------------- junctions

    def _make_junction(self, part):
        """The junction of a child part: its verdict and junction causes (W1-W7, no_bases;
        exact_gated where none fails)."""
        key = part.key
        side, pid = key[-1]
        ppart = self._by_key[key[:-1]]
        pg, cg = ppart.graphlet, part.graphlet
        a = pg.arms[side]
        leaf = derive.paths(a)[pid].leaf
        L = a.segments[leaf].end_bp
        sp, sc = _strategy(pg), _strategy(cg)
        causes = []
        if pg.support == 'trace':
            causes.append(('trace_restart', 'under support trace a lineage follows one '
                           'occurrence chain from the root seed; the child validated its seed '
                           'labels against any occurrence of the tail, so it can follow an '
                           'occurrence the parent\'s lineage never had (W1)'))
        diff = _K.first_difference(_without_allowed(sp, sc), _without_allowed(sc, sp))
        if diff is not None:
            causes.append(('strategy_differs', 'the child\'s normalized strategy differs from '
                           'the parent\'s at %s, a field a continuation keeps (W2)' % diff))
        continued, not_cont, dropped = (), (), ()
        if pg.mode == 'constrain':
            pool = {l.name for l in pg.labels[:pg.seed.num_seed_labels]} | \
                set(_knob(sp, 'labels', 'extra', []) or [])
            foreign = [n for n in (_knob(sc, 'labels', 'extra', []) or []) if n not in pool]
            if foreign:
                causes.append(('extra_foreign', 'the child\'s labels.extra names %s, outside '
                               'the parent\'s permitted pool (its seed labels and labels.extra) '
                               '(W2)' % ', '.join(json.dumps(n, ensure_ascii=False)
                                                   for n in foreign[:_SHOWN])
                               + (' and %d more' % (len(foreign) - _SHOWN)
                                  if len(foreign) > _SHOWN else '')))
            mp = (_knob(sp, 'labels', 'change_cost') or {}).get('model', 'forbid')
            mc = (_knob(sc, 'labels', 'change_cost') or {}).get('model', 'forbid')
            if mp != 'forbid' or mc != 'forbid':
                causes.append(('switching', 'labels.change_cost is %s in the parent and %s in '
                               'the child: switching restarts with the child\'s loss budget '
                               'and source lists (W3 needs forbid on both)' % (mp, mc)))
            seeds = {l.ref for l in cg.labels[:cg.seed.num_seed_labels]}
            ends = derive.end_labels(a, leaf)
            alive = [pg.labels[e.label] for e in ends]
            continued = tuple(self._label_of[l.ref] for l in alive if l.ref in seeds)
            xnames = {}
            for x in cg.dropped:
                xnames.setdefault(x.name, x.reason)
            dropped = tuple(l for l in alive if l.ref not in seeds and l.name in xnames)
            not_named = tuple(l for l in alive if l.ref not in seeds and l.name not in xnames)
            not_cont = tuple(self._label_of[l.ref] for l in alive if l.ref not in seeds)
            mlb_p = _knob(sp, 'branching', 'max_label_branches')
            mlb_c = _knob(sc, 'branching', 'max_label_branches')
            if not (mlb_p == 'unlimited' and mlb_c == 'unlimited'):
                bs = {e.branches for e in ends if pg.labels[e.label].ref in seeds}
                b = next(iter(bs)) if len(bs) == 1 else (0 if not bs else None)
                ok = b is not None and isinstance(mlb_p, int) and isinstance(mlb_c, int) \
                    and mlb_c == mlb_p - b
                if not ok:
                    causes.append(('branches_unequal', 'branching.max_label_branches is %s in '
                                   'the parent and %s in the child, and the continued labels '
                                   'ended the walk at branch counts %s: the allowance does not '
                                   'shift uniformly (W4)'
                                   % (mlb_p, mlb_c, sorted(bs) or '[]')))
        knobs = []
        for who, st in (('parent', sp), ('child', sc)):
            for block, name, bad in (
                    ('branching', 'min_successor_labels', lambda v: v > 1),
                    ('branching', 'min_successor_fraction', lambda v: v != 0),
                    ('bounds', 'min_live_labels', lambda v: v > 1),
                    ('frontier', 'on_overflow', lambda v: v != 'stop')):
                v = _knob(st, block, name)
                if v is not None and bad(v):
                    knobs.append('%s.%s %s in the %s' % (block, name, json.dumps(v), who))
        if knobs:
            causes.append(('label_count_knob', 'a knob reads other labels\' counts, which '
                           'restart at the junction: %s (W5)' % ', '.join(knobs)))
        if pg.reconverge == 'merge' and cg.reconverge == 'merge':
            sib = []
            for p in derive.paths(a):
                if p.id == pid:
                    continue
                s = a.segments[p.leaf]
                if s.leaf.continuation is not None or p.length_bp > L:
                    sib.append(p.id)
            if sib:
                pre = _walk_prefix(ppart.key)
                causes.append(('sibling_merge_unknown', 'under merge the reference joins paths '
                               'that meet beyond %d bp; walks %s of the arm went on or ended '
                               'beyond it, and the child cannot see them (W6)'
                               % (ppart.offset_bp + L,
                                  ', '.join(walk_name(pre + (w,)) for w in sib[:_SHOWN])
                                  + (' and %d more' % (len(sib) - _SHOWN)
                                     if len(sib) > _SHOWN else ''))))
        if pg.mode == 'constrain':
            if not_named:
                causes.append(('not_continued', 'the child does not seed %s, alive at the end '
                               'of the walk: their lineages do not go on in the child (W7)'
                               % _label_list(not_named)))
            if dropped:
                causes.append(('seed_dropped', 'the child named %s but dropped %s (X): %s (W7)'
                               % (_label_list(dropped), 'it' if len(dropped) == 1 else 'them',
                                  ', '.join('%s %s' % (l.name, xnames[l.name])
                                            for l in dropped[:_SHOWN]))))
        if not (_K.has_bases(pg) and _K.has_bases(cg)):
            causes.append(('no_bases', 'a part carries no bases (output.sequences false): the '
                           'rule check (W8) cannot run'))
        if not causes:
            causes.append(_GATED)
        causes.sort(key=lambda c: _RANK[c[0]])
        cseed = len(cg.seed.sequence)
        return Junction(side, _walk_prefix(ppart.key) + (pid,), ppart.offset_bp + L, cseed,
                        key, 'qualified', tuple(causes), continued, not_cont)

    def _junction_of_part(self, key):
        for i, j in enumerate(self._junctions):
            if j.part == key:
                return i, j
        raise KeyError(key)

    @staticmethod
    def _gated(j):
        """Whether junction |j| would be exact once exact is switched on (prospective)."""
        return len(j.causes) == 1 and j.causes[0][0] == 'exact_gated'

    # ---------------------------------------------------------------- the rule check (W8)

    def _context(self, key, b=None):
        """The rule check's context of the junction into part |key| (cached): the edge keys
        the reference has used on the continued walk's route (a list of sets, one per part
        of the chain), the root seed's nodes, the child seed's nodes, and the split count
        the reference carries into the child (the counts at the continued leaves, summed
        along the chain)."""
        got = self._cache.get(('ctx', key))
        if got is not None:
            return got
        side, pid = key[-1]
        ppart = self._by_key[key[:-1]]
        pg = ppart.graphlet
        a = pg.arms[side]
        leaf = derive.paths(a)[pid].leaf
        route = _ancestors(a, leaf)
        if b is not None:
            n = sum(a.segments[s].length_bp for s in route) + pg.k * len(route)
            b.charge(W_ELEM * len(route) + n, set_bytes(n) + n * (STR + pg.k + 1))
        tails = self._tails(ppart.key, side)
        edges = set()
        for s in route:
            seg = a.segments[s]
            window = _before(tails, pg, side, seg) + seg.walk
            for i in range(seg.length_bp):
                edges.add(_ops.edge_key(pg.regime, _natural(side, window[i:i + pg.k + 1]))[0])
        if ppart.key:
            up = self._context(ppart.key, b)
            chain_edges = up['edges'] + (frozenset(edges),)
            count = up['count']
        else:
            chain_edges = (frozenset(edges),)
            count = 0
        count += _split_counts(a)[leaf]
        cg = self._by_key[key].graphlet
        got = {'edges': chain_edges, 'count': count,
               'root_nodes': self._root_nodes(),
               'child_nodes': _kmers(cg.seed.sequence, cg.k, cg.regime),
               'max_splits': _knob(_strategy(pg), 'branching', 'max_splits_per_path')}
        self._cache[('ctx', key)] = got
        return got

    def _root_nodes(self):
        got = self._cache.get('root_nodes')
        if got is None:
            root = self.root
            got = self._cache['root_nodes'] = _kmers(root.seed.sequence, root.k, root.regime)
        return got

    def _tails(self, key, side):
        got = self._cache.get(('tails', key))
        if got is None:
            got = self._cache[('tails', key)] = _tails(self._by_key[key].graphlet, side)
        return got

    def _walk_causes(self, key, pid, b=None):
        """The walk-scope causes (§3.4 W8) of walk |pid| of child part |key|, in the child's
        own outward coordinates: [(code, sentence, at_bp local)], the first of each kind on
        the walk's route (the DAG ancestors of its leaf through any parent). Cached."""
        got = self._cache.get(('w8', key, pid))
        if got is not None:
            return got
        side = key[-1][0]
        part = self._by_key[key]
        cg = part.graphlet
        a = cg.arms[side]
        segs = a.segments
        k, regime = cg.k, cg.regime
        leaf = derive.paths(a)[pid].leaf
        route = _ancestors(a, leaf)
        on_route = set(route)
        if b is not None:
            n = sum(segs[s].length_bp for s in route)
            b.charge(W_ELEM * (len(route) + len(a.branch_events)) + n,
                     set_bytes(len(route)) + 4 * (STR + k + 1))
        ctx = self._context(key, b)
        edges, root_nodes, child_nodes = ctx['edges'], ctx['root_nodes'], ctx['child_nodes']
        tails = self._tails(key, side)
        found = {}

        def note(code, at, sentence):
            if code not in found or at < found[code][2]:
                found[code] = (code, sentence, at)

        def in_b(wstep):
            nat = _natural(side, wstep)
            ek = _ops.edge_key(regime, nat)[0]
            if any(ek in e for e in edges):
                return 'edge'
            target = nat[1:] if side == 'right' else nat[:k]
            if target in root_nodes and target not in child_nodes:
                return 'seed'
            return None

        def kmer_at(s, d):
            window = _before(tails, cg, side, s) + s.walk
            i = d - s.from_bp
            return window[i:i + k]

        def taken(s, d):
            """The route's next bases from the node at depth |d| of segment |s|."""
            if d < s.end_bp:
                return {s.walk[d - s.from_bp]}
            return {segs[c].walk[:1] for c in s.children if c in on_route} - {''}

        counts = _split_counts(a)
        limit = ctx['max_splits']
        for sid in route:
            s = segs[sid]
            window = _before(tails, cg, side, s) + s.walk
            for i in range(s.length_bp):
                why = in_b(window[i:i + k + 1])
                if why:
                    note('rule_violation', s.from_bp + i,
                         'the walk takes the step %s at {bp} bp, %s (a step the reference '
                         'blocks); it is not cut' % (
                             s.walk[i], 'an edge the walk used before the junction'
                             if why == 'edge' else 'into a root-seed node'))
                    break
            if s.split is not None and len(s.children) > 1:
                at = s.end_bp
                end_kmer = tails[sid]
                for c in s.children:
                    if c in on_route:
                        continue
                    ch = segs[c].walk[:1]
                    if ch and in_b(end_kmer + ch):
                        note('candidate_changed', at, 'another branch of the split at {bp} '
                             'bp, %s, is a step the reference blocks: the split may not happen '
                             'in the reference' % ch)
                if limit is not None and limit != 'unlimited' and \
                        ctx['count'] + counts[sid] + 1 > limit:
                    note('split_limit_restarted', at, 'the child split at {bp} bp with %d '
                         'splits on its route; with the %d the reference carries into the '
                         'junction that exceeds max_splits_per_path %d' % (
                             counts[sid], ctx['count'], limit))
            for ev in s.events:
                if ev.type != 'blocked' or ev.reason != 'J':
                    continue
                wstep = kmer_at(s, ev.at_bp) + ev.char
                nat = _natural(side, wstep)
                target = nat[1:] if side == 'right' else nat[:k]
                if target in root_nodes:
                    continue
                ek = _ops.edge_key(regime, nat)[0]
                sub = 'same_stop' if any(ek in e for e in edges) else 'earlier_stop'
                note('seed_reentry', ev.at_bp, '%s: the child blocked %s at {bp} bp as '
                     'rejoined_seed, a node of its seed outside the root seed; %s' % (
                         sub, ev.char,
                         'the reference blocks the same edge (edge_reuse), so only the end '
                         'reason differs' if sub == 'same_stop' else
                         'the reference takes the new edge and goes on'))
        for be in a.branch_events:
            if be.segment not in on_route:
                continue
            s = segs[be.segment]
            km = kmer_at(s, be.at_bp)
            go = taken(s, be.at_bp)
            for ch in be.chars:
                if ch in go:
                    continue
                if in_b(km + ch):
                    note('candidate_changed', be.at_bp, 'the successor %s the child admitted '
                         'at {bp} bp (a branch event) is a step the reference blocks: the '
                         'child\'s ambiguity or refusal there may not happen in the '
                         'reference' % ch)
        ecb = a.evidence_complete_to_bp
        end = segs[leaf].end_bp
        if ecb is not None and end >= ecb:
            note('decisions_unrecorded', ecb, 'the walk ends at or beyond the child\'s '
                 'evidence boundary ({bp} bp, output.max_branch_events): why it ended there '
                 'is not recorded')
        got = [found[c] for c in sorted(found, key=lambda c: (found[c][2], _RANK[c]))]
        self._cache[('w8', key, pid)] = got
        return got

    # ---------------------------------------------------------------- stitched walks

    def _rows(self, side, b=None):
        """Every stitched walk of the arm in key order, as (key, [(part, pid)]), not built."""
        got = self._cache.get(('rows', side))
        if got is not None:
            return got
        out = []

        def visit(part, chain):
            g = part.graphlet
            a = g.arms[side]
            if b is not None:
                b.charge(W_ELEM * len(a.segments))
            for p in derive.paths(a):
                here = chain + [(part, p.id)]
                kid = self._by_key.get(part.key + ((side, p.id),))
                if kid is None:
                    out.append((_walk_prefix(part.key) + (p.id,), here))
                else:
                    visit(kid, here)
        visit(self._parts[0], [])
        self._cache[('rows', side)] = out
        return out

    def _certify(self, side, chain, prospective, b=None):
        """(certified_to_bp, complete, causes in stitched coordinates, junction indexes) of a
        stitched walk (§3.5). Every walk-scope cause on the walk is stated, wherever
        certification stops. Piece i certifies its walk up to its part's complete_to_bp and,
        inside a child, its first walk-scope cause; it passes to piece i + 1 only when its
        walk ends within that depth and the junction is exact (in phase 1a-S none is;
        prospective: exact_gated read as exact)."""
        causes, jidx, first = [], [], []
        for i, (part, pid) in enumerate(chain):
            own = []
            if i > 0:
                ji, j = self._junction_of_part(part.key)
                jidx.append(ji)
                if not any(c[0] == 'no_bases' for c in j.causes):
                    own = self._walk_causes(part.key, pid, b)
            first.append(min((c[2] for c in own), default=None))
            causes.extend((c[0], c[1], part.offset_bp + c[2]) for c in own)
        cert, reached, whole = 0, False, False
        for i, (part, pid) in enumerate(chain):
            a = part.graphlet.arms[side]
            n = derive.paths(a)[pid].length_bp
            limit = a.complete_to_bp if first[i] is None else min(a.complete_to_bp, first[i])
            got = min(n, limit)
            cert = part.offset_bp + got
            whole = got >= n
            if i == len(chain) - 1:
                reached = True
                break
            if not (whole and prospective and self._gated(self._junctions[jidx[i]])):
                break
        la = chain[-1][0].graphlet.arms[side]
        end_ok = _ops._walk_end(la, derive.paths(la)[chain[-1][1]])[1]
        complete = reached and whole and end_ok and not causes
        causes = [(c, s.replace('{bp}', str(at)), at) for c, s, at in causes]
        causes.sort(key=lambda c: (c[2], _RANK[c[0]]))
        return cert, complete, causes, tuple(jidx)

    def _walk_row(self, side, key, chain, prospective=False, b=None):
        """The scalar facts of one stitched walk, as a dict (no spelling)."""
        last_part, last_pid = chain[-1]
        lg = last_part.graphlet
        la = lg.arms[side]
        p = derive.paths(la)[last_pid]
        length = last_part.offset_bp + p.length_bp
        cert, complete, causes, jidx = self._certify(side, chain, prospective, b)
        if len(chain) == 1:
            verdict = None          # no junction joins the walk: no junction verdict
        elif prospective:
            verdict = 'exact' if not causes and all(
                self._gated(self._junctions[i]) for i in jidx) else 'qualified'
        else:
            verdict = 'qualified'
        if len(chain) > 1 and not causes and all(self._gated(self._junctions[i])
                                                 for i in jidx):
            causes = [_GATED_WALK] if not prospective else []
        pieces = []
        for part, pid in chain:
            n = derive.paths(part.graphlet.arms[side])[pid].length_bp
            pieces.append((part.key, pid, part.offset_bp, part.offset_bp + n))
        return {'key': key, 'length_bp': length, 'certified_to_bp': cert,
                'complete': complete, 'verdict': verdict, 'causes': tuple(causes),
                'junctions': jidx, 'pieces': tuple(pieces), 'chain': chain}

    def _labels_full(self, side, chain):
        """Labels with displayed support on every piece (by ref), continued at every
        junction: the last piece's walk's full labels kept where each earlier piece's walk
        holds them whole too."""
        refs = None
        for part, pid in chain:
            g = part.graphlet
            a = g.arms[side]
            p = derive.paths(a)[pid]
            here = {l.ref for l in _ops._labels_full(g, a, p)}
            refs = here if refs is None else refs & here
        last = chain[-1][0].graphlet
        a = last.arms[side]
        p = derive.paths(a)[chain[-1][1]]
        return [self._label_of[l.ref] for l in _ops._labels_full(last, a, p) if l.ref in refs]

    def _walk_object(self, side, row):
        chain = row['chain']
        part, pid = chain[-1]
        g = part.graphlet
        a = g.arms[side]
        ranked = self._cache.get(('ranked', part.key, side))
        if ranked is None:
            ranked = self._cache[('ranked', part.key, side)] = {
                x[0].id: x for x in _ops._ranked_walks(g, a, 'id', None, True, 0)}
        w = _ops._walk_objects(g, a, [ranked[pid]], with_claims=False)[0]
        try:
            seq = self._spell_chain(side, chain, 'natural')
        except ValueError:
            seq = None
        full = self._labels_full(side, chain) if len(chain) > 1 else \
            [self._label_of[l.ref] for l in w.labels_full]
        return StitchedWalk(
            key=row['key'], arm=side, length_bp=row['length_bp'], sequence=seq,
            path_reason=w.path_reason, end_reasons=dict(w.end_reasons), labels_full=full,
            n_alive=w.n_alive, complete=row['complete'], verdict=row['verdict'],
            certified_to_bp=row['certified_to_bp'],
            beyond_certified_bp=row['length_bp'] - row['certified_to_bp'],
            junctions=row['junctions'], pieces=row['pieces'], causes=row['causes'])

    def _rank_keys(self, side, chain, by):
        part, pid = chain[-1]
        g = part.graphlet
        a = g.arms[side]
        p = derive.paths(a)[pid]
        full, alive_ids, consistent, min_loss = _ops._walk_keys_of(g, a, p)
        if len(chain) > 1:
            full = self._labels_full(side, chain)       # as StitchedWalk.labels_full
        elif full is None:
            full = _ops._labels_full(g, a, p)
        return full, alive_ids, consistent, min_loss

    def walks(self, arm, *, top=None, by='support', labels=None, route_consistent=True,
              min_bp=0, budget=None, resume=None):
        """[StitchedWalk] of the arm, ranked as Graphlet.walks() ranks its walks (by
        'support', 'length', 'loss' -- the loss of the part holding the walk's end; joined
        losses are phase 1b -- or 'id', the key order), ties by key; |labels| keeps the walks
        whose end carries one of them (the part's route_consistent rule), |min_bp| the walks
        of at least that stitched length. budget=: a stop raises; with by='id' its .partial
        holds the whole walks made so far and a resume token (resume=)."""
        side = self.arm_of(arm)
        if by not in ('support', 'length', 'loss', 'id'):
            raise ValueError("by is 'support', 'length', 'loss' or 'id'")
        if resume is not None and by != 'id':
            raise ValueError("resume= continues walks(by='id') only: a ranked list is never "
                             "partial")
        wanted = None
        if labels is not None:
            sel = labels if isinstance(labels, (list, tuple, set, frozenset)) else [labels]
            wanted = {self.label(s).ref for s in sel}
        args = [side, top, by, sorted(wanted) if wanted is not None else None,
                bool(route_consistent), min_bp]
        b = _B.resolve(budget)
        out = []
        k0 = 0
        try:
            if b is not None:
                scope = b.scope('walks', ('rank_id', 'select_walks', 'narrow_labels',
                                          'narrow_arm'))
                scope.__enter__()
            try:
                k0, = _K.resume_pos(self.digest, 'walks', args, resume, 1)
                rows = []
                for key, chain in self._rows(side, b):
                    part, pid = chain[-1]
                    n = part.offset_bp + derive.paths(part.graphlet.arms[side])[pid].length_bp
                    if n < min_bp:
                        continue
                    full, alive_ids, consistent, min_loss = self._rank_keys(side, chain, by)
                    if wanted is not None:
                        g = part.graphlet
                        pool = consistent if route_consistent else set(alive_ids)
                        if not wanted & {g.labels[i].ref for i in pool}:
                            continue
                    rows.append((key, chain, n, len(full), len(alive_ids), min_loss))
                if b is not None:
                    b.charge(W_ELEM * len(rows) + (_B.sort_work(len(rows)) if by != 'id'
                                                   else 0))
                order = {
                    'support': lambda r: (-r[3], -r[4], -r[2], r[0]),
                    'length': lambda r: (-r[2], r[0]),
                    'loss': lambda r: (r[5], -r[4], -r[2], r[0]),
                    'id': lambda r: r[0]}[by]
                rows.sort(key=order)
                if top is not None:
                    rows = rows[:top]
                for key, chain, n, _, _, _ in rows[k0:]:
                    if b is not None:
                        b.charge(W_ROW + W_STEP * n // 64 + (n >> 8) + 4 * len(chain),
                                 record_bytes(StitchedWalk) + 3 * (STR + n) + 400)
                    out.append(self._walk_object(side, self._walk_row(side, key, chain,
                                                                      b=b)))
                return out
            finally:
                if b is not None:
                    scope.__exit__(None, None, None)
        except LocalBudgetExceeded as e:
            if by == 'id':
                raise e.restate(Partial(out, _K.token(self.digest, 'walks', args,
                                                      k0 + len(out))),
                                rows=len(out), position=k0 + len(out))
            raise

    def _row(self, side, key, prospective=False, b=None):
        key = self._walk_key(key)
        for k, chain in self._rows(side, b):
            if k == key:
                return self._walk_row(side, k, chain, prospective, b)
        raise IndexError('no stitched walk %s on the %s arm (a continued walk is replaced by '
                         'its child\'s walks)' % (walk_name(key), side))

    def walk(self, arm, walk, *, budget=None):
        """One StitchedWalk by key."""
        side = self.arm_of(arm)
        return _K.run_scoped(budget, 'walks', lambda b: self._walk_object(
            side, self._row(side, walk, False, b)))

    # ---------------------------------------------------------------- spelling

    def _spell_chain(self, side, chain, orientation, b=None):
        pieces = []
        for part, pid in chain:
            g = part.graphlet
            a = g.arms[side]
            leaf = derive.paths(a)[pid].leaf
            for s in derive.chain(a, leaf):
                seg = a.segments[s]
                pieces.append((g, side, s, seg.from_bp, seg.end_bp))
        return _K.spell_pieces(pieces, side, orientation, b)

    def spell(self, arm, walk, orientation='natural', with_seed=False, *, budget=None):
        """The stitched walk's bases, from the parts without the overlap (§3.2): right
        flank = the parent's walk then the child's, left flank = the child's then the
        parent's; with the root seed as Graphlet.spell() adds it ('walk': walking order, the
        seed read in the walking direction)."""
        side = self.arm_of(arm)
        if orientation not in ('natural', 'walk'):
            raise ValueError("orientation is 'natural' or 'walk'")

        def run(b):
            row = self._row(side, walk, False, b)
            w = self._spell_chain(side, row['chain'], orientation, b)
            if not with_seed:
                return w
            seed = self.seed.sequence
            if orientation == 'natural':
                return seed + w if side == 'right' else w + seed
            return (seed if side == 'right' else seed[::-1]) + w
        return _K.run_scoped(budget, 'spell', run)

    # ---------------------------------------------------------------- continuation

    def continuation(self, arm, walk, *, budget=None):
        """The continuation of a stitched walk: its part's (Graphlet.continuation()), the
        seed coordinates mapped to the root seed (a child seed coordinate x is root
        coordinate x + |root seed| + offset - |child seed| on the right arm, x - offset on
        the left), the labels the union's."""
        side = self.arm_of(arm)
        key = self._walk_key(walk)
        part, pid = self._part_of_walk(side, key)
        c = part.graphlet.continuation(side, pid, budget=budget)
        g = part.graphlet
        if part.key:
            n = len(g.seed.sequence)
            if side == 'right':
                sh = len(self.seed.sequence) + part.offset_bp - n
            else:
                sh = -part.offset_bp
            coord = (c.seed_coord[0] + sh, c.seed_coord[1] + sh)
        else:
            coord = c.seed_coord
        return Continuation(c.sequence, [self._label_of[l.ref] for l in c.labels],
                            c.loss_used, c.branches_used, coord, side, c.leaf, c.unverifiable,
                            c.losses, c.note)

    def _one_part(self, side, walks):
        if isinstance(walks, (int, tuple)) and not isinstance(walks, bool) and (
                isinstance(walks, int) or all(isinstance(x, int) for x in walks)):
            walks = [walks]
        keys = [self._walk_key(w) for w in walks]
        located = [self._part_of_walk(side, k) for k in keys]
        if len({p.key for p, _ in located}) > 1:
            raise ValueError('next_request() continues walks of one part: %s lie in parts %s; '
                             'build one request per part' % (
                                 ', '.join(walk_name(k) for k in keys),
                                 ', '.join(sorted({part_name(p.key) for p, _ in located}))))
        return located[0][0], [pid for _, pid in located]

    def next_request(self, arm, walks, bp=None, reduce_budget=True, reset_branches=False,
                     **overrides):
        """The part owning the walks' next_request() (walk keys mapped to its path ids;
        walks of one part)."""
        side = self.arm_of(arm)
        part, pids = self._one_part(side, walks)
        return part.graphlet.next_request(side, pids, bp, reduce_budget, reset_branches,
                                          **overrides)

    def next_requests(self, arm, walks, bp=None, reduce_budget=True, reset_branches=False,
                      **overrides):
        """The part owning the walks' next_requests() (walks of one part)."""
        side = self.arm_of(arm)
        part, pids = self._one_part(side, walks)
        return part.graphlet.next_requests(side, pids, bp, reduce_budget, reset_branches,
                                           **overrides)

    # ---------------------------------------------------------------- exports

    def to_fasta(self, arm=None, walks=None, with_seed=True, orientation='natural', width=None,
                 *, coordinates=None, budget=None):
        """FASTA, one record per stitched walk in key order (or |walks|, keys, in that order),
        spelled from the parts without the overlap. Headers are Graphlet.to_fasta()'s
        ('>{arm}_{walk} length_bp=.. reason=.. labels=.. seed=..') plus walk=3.0,
        junctions=<at_bp,...> ('.' for none), verdict= and certified_to=. A part without
        bases is refused as its own export would be (ValueError). coordinates=True is
        refused (record coordinates across a junction are phase 1b); None and False write
        none."""
        if coordinates is True:
            raise ValueError('record coordinates across a junction are not offered on a '
                             'stitched graphlet (phase 1b): export the parts')
        if coordinates is not None and not isinstance(coordinates, bool):
            raise ValueError('coordinates is None (when present), True or False, not %r'
                             % (coordinates,))
        if orientation not in ('natural', 'walk'):
            raise ValueError("orientation is 'natural' or 'walk'")
        sides = [self.arm_of(arm)] if arm is not None else self.arms()

        def run(b):
            out = []
            for side in sides:
                rows = self._rows(side, b)
                if walks is None:
                    chosen = rows
                else:
                    by = dict(rows)
                    chosen = []
                    for w in walks:
                        k = self._walk_key(w)
                        if k not in by:
                            raise IndexError('no stitched walk %s on the %s arm'
                                             % (walk_name(k), side))
                        chosen.append((k, by[k]))
                for key, chain in chosen:
                    row = self._walk_row(side, key, chain, False, b)
                    if b is not None:
                        b.charge(W_ROW + (row['length_bp'] >> 8),
                                 3 * (STR + row['length_bp'] + len(self.seed.sequence)) + 200)
                    w = self._spell_chain(side, chain, 'walk', b)
                    part, pid = chain[-1]
                    g = part.graphlet
                    a = g.arms[side]
                    p = derive.paths(a)[pid]
                    seg = a.segments[p.leaf]
                    reason = (REASON[seg.leaf.path_reason] if seg.leaf.path_reason
                              else ','.join(sorted({a.runs[e.run].reason
                                                    for e in derive.end_labels(a, p.leaf)}))
                              or '.')
                    nlab = len(derive.end_labels(a, p.leaf)) if g.mode == 'constrain' \
                        else len(seg.end)
                    js = ','.join(str(self._junctions[i].at_bp) for i in row['junctions'])
                    out.append('>%s_%s length_bp=%d reason=%s labels=%d seed=%s walk=%s '
                               'junctions=%s verdict=%s certified_to=%d' % (
                                   side, walk_name(key), row['length_bp'], reason, nlab,
                                   'yes' if with_seed else 'no', walk_name(key), js or '.',
                                   row['verdict'] or '.', row['certified_to_bp']))
                    out.append(_export._wrap(_export._oriented(self.root, side, w, orientation,
                                                               with_seed), width))
            if out:
                out.append('')
            return '\n'.join(out)
        return _K.run_scoped(budget, 'to_fasta', run)

    # ---------------------------------------------------------------- completeness

    def completeness(self):
        """Per arm: the root's complete_to_bp and scope, unchanged (stitching never makes
        "every walk to depth D" true for a larger D), and one extension per junction on the
        arm, by (arm, walk key): {walk, at_bp, extended_to_bp (at_bp + the child's
        complete_to_bp), certified_to_bp (how far the certificate reaches through the
        junction: no junction is exact in phase 1a-S, so at most at_bp, the continued walk's
        own certified depth), verdict, causes (the junction's cause codes)}."""
        return self._completeness(False)

    def _chain_of(self, side, key):
        """[(part, path id)] of walk key |key| on arm |side|, root first."""
        return [self._part_of_walk(side, key[:i + 1]) for i in range(len(key))]

    def _completeness(self, prospective):
        out = {}
        for side in self.arms():
            a = self.root.arms[side]
            exts = []
            for j in self._junctions:
                if j.arm != side:
                    continue
                cg = self._by_key[j.part].graphlet
                ext = j.at_bp + cg.arms[side].complete_to_bp
                cert = min(j.at_bp, self._certify(side, self._chain_of(side, j.walk),
                                                  prospective)[0])
                if prospective and self._gated(j) and cert == j.at_bp:
                    # every walk below the junction certified to the extension, or as far
                    # as the first one that is not complete
                    cert = ext
                    for key, chain in self._rows(side):
                        if key[:len(j.walk)] == j.walk:
                            row = self._walk_row(side, key, chain, True)
                            if not row['complete']:
                                cert = min(cert, row['certified_to_bp'])
                exts.append({'walk': j.walk, 'at_bp': j.at_bp, 'extended_to_bp': ext,
                             'certified_to_bp': cert,
                             'verdict': ('exact' if self._gated(j) else 'qualified')
                             if prospective else j.verdict,
                             'causes': [c[0] for c in j.causes
                                        if not (prospective and c[0] == 'exact_gated')]})
            out[side] = {'complete_to_bp': a.complete_to_bp, 'scope': a.scope,
                         'extensions': exts}
        return out

    # ---------------------------------------------------------------- evidence

    def evidence(self, arm=None, walks=None):
        """The evidence of an answer (§3.5): {parts: each part's own evidence_block keyed by
        part name ('root', 'right:3', ...), junctions: those the answer reads with their
        verdicts and causes, walk_causes: the walk-scope causes on the walks it reads,
        identity, outcome: the weakest value per dimension over the parts and junctions read
        (a qualified junction or a walk-scope cause: walks partial, label_evidence
        qualified), exact}. |arm|: the arm's parts and junctions (None: all); |walks|: only
        the parts and junctions those walks (keys, on |arm|) read -- a root walk no child
        continues reads the root alone and states its evidence. exact holds only when every
        part read is exact, no junction is read (no junction is exact in phase 1a-S) and no
        walk read has a walk-scope cause."""
        return self._evidence(arm, walks, False)

    def _evidence(self, arm, walks, prospective):
        """evidence(); |prospective|: the private test hook of §3.5 (exact_gated read as
        exact; where a junction is read, exact needs the identity verified)."""
        sides = [self.arm_of(arm)] if arm is not None else self.arms()
        if walks is not None and arm is None and len(sides) > 1:
            raise ValueError('walks= names walks of one arm: give arm=')
        read_parts, read_j, wcauses = set(), set(), []
        for side in sides:
            if walks is None:
                chosen = self._rows(side)
            else:
                ws = walks if isinstance(walks, list) else [walks]
                chosen = [(r['key'], r['chain']) for r in
                          (self._row(side, w, prospective) for w in ws)]
            for key, chain in chosen:
                for part, _ in chain:
                    read_parts.add(part.key)
                row = self._walk_row(side, key, chain, prospective)
                read_j.update(row['junctions'])
                for c in row['causes']:
                    if c[0] != 'exact_gated':
                        wcauses.append({'walk': key, 'code': c[0], 'at_bp': c[2],
                                        'sentence': c[1]})
        parts = [p for p in self._parts if p.key in read_parts]
        verified = bool(self.identity.get('verified')) or not read_j
        ev = _K.combined_evidence([(part_name(p.key), p.graphlet, True) for p in parts],
                                  sides[0] if len(sides) == 1 else None)
        js = [self._junctions[i] for i in sorted(read_j)]
        outcome = ev['outcome']
        weak = any(not (prospective and self._gated(j)) for j in js) or bool(wcauses)
        if weak and outcome is not None:
            outcome = _K.weakest_outcome([outcome, dict(outcome, walks='partial',
                                                        label_evidence='qualified')])
        exact = ev['exact'] and not weak and (verified or not js)
        if not prospective:
            exact = exact and not js
        return {'parts': ev['parts'],
                'junctions': [{'arm': j.arm, 'walk': j.walk, 'at_bp': j.at_bp,
                               'verdict': ('exact' if self._gated(j) else 'qualified')
                               if prospective else j.verdict,
                               'causes': [list(c) for c in j.causes
                                          if not (prospective and c[0] == 'exact_gated')]}
                              for j in js],
                'walk_causes': wcauses,
                'identity': dict(self.identity),
                'outcome': outcome,
                'exact': exact,
                'experimental': True}

    def _prospective(self, arm, walk=None):
        """The private test hook of §3.5 (T4, T11), never part of a public answer: the
        result with exact_gated read as exact. With |walk|: that stitched walk's row
        {key, length_bp, certified_to_bp, complete, verdict, causes, junctions, pieces};
        without: {completeness: the arm's, evidence: the arm's}."""
        side = self.arm_of(arm)
        if walk is not None:
            return self._row(side, walk, True)
        return {'completeness': self._completeness(True)[side],
                'evidence': self._evidence(side, None, True)}

    # ---------------------------------------------------------------- summary, provenance

    def summary(self, max_bytes=2048, *, budget=None):
        """Agent-facing, <= max_bytes of JSON: the root seed, per arm the root's certificate
        and the stitched walks' count and longest, the junctions with their verdicts and
        cause codes, the identity, the outcome and exact of evidence(); marked
        experimental. Junctions beyond the cap are cut (stated)."""
        def run(b):
            meta = (self.root.seed_summary or {}).get('seed') or {}
            out = {'stitched': True, 'experimental': True,
                   'seed': {'seed_id': meta.get('seed_id', ''),
                            'length_bp': self.seed.length_bp, 'mode': self.mode,
                            'labels': len(self._labels)},
                   'parts': len(self._parts), 'arms': {}, 'junctions': []}
            for side in self.arms():
                a = self.root.arms[side]
                rows = self._rows(side, b)
                lengths = [chain[-1][0].offset_bp + derive.paths(
                    chain[-1][0].graphlet.arms[side])[chain[-1][1]].length_bp
                    for _, chain in rows]
                out['arms'][side] = {'complete_to_bp': a.complete_to_bp, 'scope': a.scope,
                                     'walks': len(rows), 'max_bp': max(lengths, default=0)}
            for j in self._junctions:
                out['junctions'].append({'arm': j.arm, 'walk': walk_name(j.walk),
                                         'at_bp': j.at_bp, 'verdict': j.verdict,
                                         'causes': [c[0] for c in j.causes]})
            ev = self.evidence()
            out['identity'] = {'verified': self.identity.get('verified')}
            out['outcome'] = ev['outcome']
            out['exact'] = ev['exact']

            def size():
                if b is not None:
                    b.charge(W_ROW)
                return len(json.dumps(out, separators=(',', ':'),
                                      ensure_ascii=False).encode())
            while size() > max_bytes and out['junctions']:
                out['junctions'].pop()
                out['junctions_cut'] = True
            return out
        return _K.run_scoped(budget, 'summary', run)

    def provenance(self):
        """One entry per part (§3.6): {key, role: root | child, junction: {arm, walk, at_bp,
        overlap_bp} | None, pair: the pair it was answered on, identity, seed_id, digest,
        request: the part's replay request as GraphletStore.request_of() builds it}. No
        single request reproduces a stitched graphlet: replay each part and stitch again."""
        from .store import GraphletStore
        out = []
        for p in self._parts:
            g = p.graphlet
            j = None
            if p.key:
                _, jj = self._junction_of_part(p.key)
                j = {'arm': jj.arm, 'walk': jj.walk, 'at_bp': jj.at_bp,
                     'overlap_bp': jj.overlap_bp}
            meta = (g.seed_summary or {}).get('seed') or {}
            out.append({'key': p.key, 'role': 'child' if p.key else 'root', 'junction': j,
                        'pair': answered_pair(g.envelope), 'identity': _ops.index_identity(g),
                        'seed_id': meta.get('seed_id'), 'digest': p.digest,
                        'request': GraphletStore.request_of(g)})
        return out

    def memory_bytes(self):
        """The parts' models (a model held twice counted once) and the junction indexes."""
        extra = 0
        for key, v in self._cache.items():
            if isinstance(key, tuple) and key[0] == 'ctx':
                n = sum(len(e) for e in v['edges'][-1:]) + len(v['child_nodes'])
                extra += set_bytes(n) + n * (STR + self.k + 1)
        return _K.models_memory_bytes([p.graphlet for p in self._parts], extra)

    # ---------------------------------------------------------------- not offered

    def _refuse(self, name):
        raise CombineRefused('not_offered', _T_NOT_OFFERED.get(name, _T_NOT_OFFERED_ANY % name))

    def dump(self, *args, **kw):
        self._refuse('dump')

    def to_json(self, *args, **kw):
        self._refuse('to_json')

    def compare(self, *args, **kw):
        self._refuse('compare')

    def subgraph(self, *args, **kw):
        self._refuse('subgraph')

    def to_gfa(self, *args, **kw):
        self._refuse('to_gfa')

    def label_walks(self, *args, **kw):
        self._refuse('label_walks')

    def routes(self, *args, **kw):
        self._refuse('routes')

    def label_summary(self, *args, **kw):
        self._refuse('label_summary')

    def splits(self, *args, **kw):
        self._refuse('splits')

    def claims(self, *args, **kw):
        self._refuse('claims')

    def support_profile(self, *args, **kw):
        self._refuse('support_profile')

    def support_changes(self, *args, **kw):
        self._refuse('support_changes')

    def save(self, *args, **kw):
        self._refuse('save')

    def __repr__(self):
        return '<StitchedGraphlet parts=%d junctions=%d digest=%s>' % (
            len(self._parts), len(self._junctions), self.digest[:16])


# ------------------------------------------------------------------ identity

def _identity_alone(g):
    if g.index_fp:
        return {'verified': True, 'index_fp': g.index_fp}
    return {'verified': False, 'note': 'no index manifest digest (index_fp): the identity of '
                                       'a later part cannot be verified'}


def _weaker_identity(a, b):
    """The identity statement of a graphlet built from two: verified only when both are."""
    if a.get('verified') and b.get('verified'):
        return a
    if not a.get('verified') and a.get('accepted'):
        return a
    if not b.get('verified'):
        return b
    return a


class _StitchModule(type(json)):
    """`import metagraph.traverse.stitch` binds the package attribute `stitch` to this
    module (the import system does so after the module ran); calling it calls stitch(),
    so metagraph.traverse.stitch(parent, child) works however the name was bound."""

    def __call__(self, *args, **kw):
        return stitch(*args, **kw)


import sys as _sys  # noqa: E402
_sys.modules[__name__].__class__ = _StitchModule
