"""Stitching a continuation to the walk it continues (docs/DESIGN-graphlet-merge.md §3,
phase 1a-S). EXPERIMENTAL: the API may change until `exact` is switched on. A private
module: its names are exported from metagraph.traverse, where `stitch` is the function.

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
with the one cause `exact_gated` where the rules found no difference (at a junction: W1-W7
hold there and W8 is judged per walk below it). Completeness,
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
from .budget import (STR, TUPLE, W_ELEM, W_ROW, W_STEP, LocalBudgetExceeded, Partial,
                     dict_bytes, list_bytes, record_bytes, set_bytes)
from .client import PairIdError, answered_pair
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
_T_INDEX_UNVERIFIABLE_PARTS = (
    'the parts state an index manifest digest (index_fp) unevenly (the parent\'s on %d of '
    '%d parts, the child\'s on %d of %d): whether they are the same index cannot be '
    'verified (not proven different), so the child is not stitched automatically; pass '
    'allow_unverified_index=True to stitch it anyway (the result then states the identity '
    'unverified)')
_T_NO_WALK = ('walk= names no walk of the %s arm (%s): a stitched child is seeded with a '
              'walk\'s continuation (continuation(), next_request())')
_T_DIFFERENT_PAIR = ('the child was answered on pair %s, its parent on pair %s: a '
                     'continuation runs on its parent\'s pair')
# p2_pair_id, bad_pair_id, pair_id_mismatch: the parent's or the child's pair statement the
# library does not take, with the words of client.PairIdError
_T_BAD_PAIR_STATEMENT = '%s: %s'
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
_T_SUMMARY_TOO_SMALL = ('the summary\'s statements alone take %d bytes, more than '
                        'max_bytes=%d: they are never cut; pass max_bytes >= %d')
_T_LABEL_ID = ('label ids are local to each part of a stitched graphlet: select by ref or '
               'name')

_SHOWN = 8          # labels a refusal names (then "and N more"), as next_request's notes
_EV_BLOCK = 1024    # bytes of a part's evidence block (888 deep on the stitch fixtures)


def _label_list(labels):
    shown = [json.dumps({'name': l.name, 'ref': l.ref}, ensure_ascii=False)
             for l in labels[:_SHOWN]]
    more = len(labels) - _SHOWN
    return ', '.join(shown) + (' and %d more' % more if more > 0 else '')


# ------------------------------------------------------------------ causes (B.2)

# the order of Appendix B.2's table: a junction's and a walk's causes are listed in it
CAUSES = ('trace_restart', 'strategy_differs', 'extra_foreign', 'switching',
          'branches_unequal', 'label_count_knob', 'sibling_merge_unknown', 'not_continued',
          'seed_dropped', 'no_bases', 'strategy_unknown', 'rule_violation', 'candidate_changed',
          'exact_gated', 'decisions_unrecorded', 'seed_reentry', 'split_limit_restarted')
JUNCTION_CAUSES = CAUSES[:11] + ('exact_gated',)
WALK_CAUSES = ('rule_violation', 'candidate_changed', 'decisions_unrecorded', 'seed_reentry',
               'split_limit_restarted', 'exact_gated')
_RANK = {c: i for i, c in enumerate(CAUSES)}

_GATED = ('exact_gated', 'W1-W7 hold at this junction, and W8 is judged per walk below it (the '
          'walks\' causes), but exact is not switched on: the junction rules found no '
          'difference from one uninterrupted walk and have not yet passed the randomized '
          'comparison against uninterrupted traversals')
_T_NO_BASES = ('on the route from the root to this child, %s carr%s no bases '
               '(output.sequences false): the rule check (W8) cannot run below this junction')
_GATED_WALK = ('exact_gated', 'W1-W8 hold for this walk and its junctions, but exact is not '
               'switched on: the rules found no difference from one uninterrupted walk and '
               'have not yet passed the randomized comparison', None)

# W2: the strategy fields a continuation may change (judged by W3-W5 where they matter);
# `clamped` is the server's report of the bounds it clamped, not a knob
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


def _kmers_price(seq, k, regime):
    """(lwu, bytes) of _kmers(seq, k, regime): a step per k-mer and orientation, the set and
    its strings."""
    n = max(0, len(seq) - k + 1) * (1 if regime == 'basic' else 2)
    return n, set_bytes(n) + n * (STR + k)


def _tails_price(g, side):
    """(lwu, bytes) of _tails(g, side): an element per segment, the list and its strings."""
    n = len(g.arms[side].segments)
    return W_ELEM * n, list_bytes(n) + n * (STR + g.k)


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


# the strategy fields W3-W5 and the split check read; without one of them on either part
# the junction is strategy_unknown (a missing max_splits_per_path is unknown, not
# unlimited). The split limit does not apply in annotate mode (SPEC §6.4)
_W5_KNOBS = (('branching', 'min_successor_labels', lambda v: v > 1),
             ('branching', 'min_successor_fraction', lambda v: v != 0),
             ('bounds', 'min_live_labels', lambda v: v > 1),
             ('frontier', 'on_overflow', lambda v: v != 'stop'))
# W2: the fields of the server's normalized strategy (schema_version 1) that W2 compares,
# every one the server states in either mode (the committed CLI responses) but the allowed
# ones; one absent from both envelopes is strategy_unknown too (max_live_paths not where
# both stop on overflow, when W2 allows it)
_W2_STATED = (('exhaustive',), ('schema_version',), ('support',),
              ('bounds', 'max_live_paths'), ('bounds', 'min_live_labels'),
              ('branching', 'bubble_window_bp'), ('branching', 'hairpins'),
              ('branching', 'max_splits_per_path'), ('branching', 'min_successor_fraction'),
              ('branching', 'min_successor_labels'), ('branching', 'on_reconverge'),
              ('branching', 'tip_window_bp'), ('frontier', 'on_overflow'),
              ('frontier', 'order'), ('labels', 'change_cost'), ('labels', 'max_seed_labels'),
              ('labels', 'max_switch_sources'), ('labels', 'mode'),
              ('labels', 'seed_label_kind'), ('labels', 'switch_on'))
_REQUIRED = {'annotate': tuple((b, n) for b, n, _ in _W5_KNOBS),
             'constrain': tuple((b, n) for b, n, _ in _W5_KNOBS) + (
                 ('labels', 'change_cost'), ('labels', 'extra'),
                 ('branching', 'max_label_branches'), ('branching', 'max_splits_per_path'))}


def _has(st, path):
    """Whether |st| states field |path| ((block, name) or (name,)) with a value."""
    if len(path) == 1:
        return st.get(path[0]) is not None
    block = st.get(path[0])
    return isinstance(block, dict) and block.get(path[1]) is not None


def _paths_of(d, prefix=()):
    """The leaf field paths of a JSON object (a dict's own fields; a list is a leaf)."""
    out = []
    for k, v in d.items():
        if isinstance(v, dict) and v:
            out.extend(_paths_of(v, prefix + (k,)))
        else:
            out.append(prefix + (k,))
    return out


def _common(a, b):
    """|a| restricted to the fields |b| has too (recursively), for W2's comparison of two
    strategies of which one may not state every field; a null field is not stated."""
    out = {}
    for k, v in a.items():
        if k not in b or v is None or b[k] is None:
            continue
        if isinstance(v, dict) and isinstance(b[k], dict):
            out[k] = _common(v, b[k])
        else:
            out[k] = v
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
    not_combinable (also a TypeError), view_unsupported, MissingEnvelope, p2_pair_id,
    bad_pair_id and pair_id_mismatch (a pair statement the library does not take, with the
    PairIdError's words: client._pair()), different_index,
    index_unverifiable (a manifest digest on some parts only; allow_unverified_index=True
    stitches it and states the identity unverified), different_pair (the identity checks
    run against every part on both sides), different_graph, different_mode, not_one_arm,
    wrong_arm, no_continuation, not_a_continuation, ambiguous_walk, already_continued (the
    same retrieval, or for a stitched child the same subtree, again at that walk is a
    no-op: the parent is returned), foreign_seed_labels. budget= (local limits): the identity
    checks, the new graphlet's index, labels and digest, and each junction's W1-W7 (computed
    once per pair of parts) are charged, each before it is built; a plain parent's own
    StitchedGraphlet is built inside the call and charged so (stitch(g) included); W8 runs
    per walk in the queries and is charged there. A stop raises with no partial."""
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
    # the pair each states, read by the one rule (client._pair()): a p2 or malformed
    # pair_id, an id beside one path or beside paths that are not its is refused with the
    # PairIdError's code and words (a StitchedGraphlet's parts were read when it was built)
    for x, who in ((parent, 'the parent'), (child, 'the child')):
        if isinstance(x, Graphlet):
            try:
                answered_pair(x.envelope)
            except PairIdError as e:
                raise CombineRefused(e.code, _T_BAD_PAIR_STATEMENT % (who, e)) from None
    if child is None and isinstance(parent, StitchedGraphlet):
        return parent

    def build(b):
        # a plain parent's wrapper is built inside the call's scope, its labels and digest
        # admitted first (GPT final review M3)
        base = parent if isinstance(parent, StitchedGraphlet) else \
            StitchedGraphlet._root_of(parent, b)
        if child is None:
            return base
        return base._stitched(child, arm, walk, allow_unverified_index, b)
    return _K.run_scoped(budget, 'stitch', build)


# ------------------------------------------------------------------ the graphlet

class StitchedGraphlet:
    """A continuation joined to the walk it continues (§3; experimental). Built by stitch().
    Its parts are held unchanged; every answer is computed from them, positions outward
    from the root seed's boundary, walks keyed by tuples of walk ids along the chain."""

    def __init__(self, parts, identity, inherited=(), b=None):
        # parts: [Part] in canonical order (root first, then by key); |inherited|: the
        # junction cores of the graphlets this one is built from
        self._parts = tuple(parts)
        self._by_key = {p.key: p for p in self._parts}
        self.identity = identity
        # the queries' derivations (contexts, walk causes, rows, tails, root nodes) and,
        # under the same keys, their cold prices with the derivations each read (_uses)
        self._cache = {}
        self._price = {}
        self._jcore = {}
        root = self._parts[0].graphlet
        self.k, self.regime, self.alphabet = root.k, root.regime, root.alphabet
        self.mode, self.support = root.mode, root.support
        self.seed = root.seed
        self._labels, self._label_of = self._union_labels()
        self._junctions = tuple(sorted((self._make_junction(p, inherited, b)
                                        for p in self._parts[1:]),
                                       key=lambda j: (j.arm, j.walk)))
        self._jindex = {j.part: i for i, j in enumerate(self._junctions)}
        self.digest = self._digest_of([(p.key, p) for p in self._parts])

    def _rid(self, g):
        """The retrieval identity of |g|, a graphlet this one holds (cached by id(), which
        is safe only for held objects)."""
        got = self._jcore.get(('rid', id(g)))
        if got is None:
            got = self._jcore[('rid', id(g))] = _K.retrieval_identity(g)
        return got

    def _digest_of(self, keyed):
        """The digest of a tree of parts, [(key relative to its root, Part)] in canonical
        order: the sha-256 of each part's key, body digest and retrieval identity (§3.7, as
        the union's digest states retrieval identity), so two stitched graphlets of equal
        bodies from other requests, whose junctions can read otherwise, differ."""
        return _K.combined_digest([[[list(s) for s in k], p.digest, self._rid(p.graphlet)]
                                   for k, p in keyed])

    # ---------------------------------------------------------------- construction

    @classmethod
    def _root_of(cls, g, b=None):
        if b is not None:
            w, m = _build_price([g])
            b.charge(W_ROW + w, record_bytes(Part) + 200 + m)
        return cls([Part((), g, 0, _K.part_digest(g))], _identity_of([g], [], False))

    def _stitched(self, child, arm, walk, allow, b):
        """The checks of §3.1 in order, then the new graphlet (or self: a no-op)."""
        croot = child.root if isinstance(child, StitchedGraphlet) else child
        proot = self.root
        if b is not None:
            # the names check over the child's labels, and its seed, case-mapped to find the
            # walk it continues (GPT final review M3: admitted before they run)
            b.charge(W_ELEM * len(child.labels), STR + len(croot.seed.sequence))
        cparts = child.parts if isinstance(child, StitchedGraphlet) else \
            (Part((), child, 0, _K.part_digest(child)),)
        identity = self._check_index(cparts, child.labels, allow, b)
        if not identity['verified']:
            self._check_pairs(cparts)
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
        ckey = part.key + ((side, pid),)
        if ckey in self._by_key:
            self._same_child(ckey, cparts, key, side, isinstance(child, StitchedGraphlet), b)
            return self
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
        new = tuple(Part(ckey + p.key, p.graphlet, offset + p.offset_bp, p.digest)
                    for p in cparts)
        parts = sorted(self._parts + new, key=lambda p: p.key)
        if b is not None:
            # the new graphlet's own index (parts, keys, junctions), its labels and digest
            # (_build_price); a junction's core is charged when it is computed (once per
            # pair of parts), W8 in the queries
            w, m = _build_price([p.graphlet for p in parts])
            b.charge(W_ROW * len(parts) + w, len(parts) * (record_bytes(Part) + 200) + m)
        inherited = (self._jcore,) + ((child._jcore,) if isinstance(child, StitchedGraphlet)
                                      else ())
        return StitchedGraphlet(parts, identity, inherited, b)

    def _same_child(self, ckey, cparts, key, side, subtree, b=None):
        """At a continued walk: returns when the child is already there (a no-op), else
        already_continued. A stitched child (|subtree|) must be the whole subtree there, every
        part's key, body and retrieval identity equal; a plain graphlet the part at the walk,
        whatever continues it. The child's retrieval identity is computed, never cached: the
        cache is keyed by the graphlets this one holds. |b|: the identities' case-mapped
        seeds, and the digest of the subtree a refusal names (_build_price), charged before
        they are made."""
        tree = [(p.key[len(ckey):], p) for p in self._parts if p.key[:len(ckey)] == ckey]
        have = tree if subtree else [((), self._by_key[ckey])]
        new = [(p.key, p) for p in cparts]
        if [(k, p.digest) for k, p in have] == [(k, p.digest) for k, p in new]:
            if b is not None:
                b.charge(W_ROW * len(new), sum(STR + len(p.graphlet.seed.sequence)
                                               for p in cparts))
            for (k, a), (_, c) in zip(have, new):
                diff = _K.first_difference(self._rid(a.graphlet),
                                           _K.retrieval_identity(c.graphlet))
                if diff is not None:
                    raise CombineRefused('already_continued', _T_ALREADY_REQUEST % (
                        walk_name(key), side, diff if not k else '%s of part %s' % (
                            diff, part_name(ckey + k))))
            return
        if b is not None:
            b.charge(*_build_price([p.graphlet for _, p in tree]))
        raise CombineRefused('already_continued', _T_ALREADY % (
            walk_name(key), side, part_name(ckey), self._digest_of(tree)))

    def _check_index(self, cparts, labels, allow, b=None):
        """different_index over every part on both sides (each distinct index identity of
        one side against each of the other), the names rule over the union's labels, then
        index_unverifiable (allow_unverified_index waives only that) -> the identity
        statement of the stitched graphlet."""
        mine = [p.graphlet for p in self._parts]
        theirs = [p.graphlet for p in cparts]
        done = set()
        for a in mine:
            for c in theirs:
                pair = (_index_key(a), _index_key(c))
                if pair in done:
                    continue
                done.add(pair)
                if b is not None:
                    b.charge(W_ROW)
                key, values = _K.index_check(a, c)
                if key in ('meta_fp', 'fp', 'ns', 'release'):
                    raise CombineRefused('different_index', _T_DIFFERENT_INDEX
                                         % _K.different_index_clause(key, values))
        names = self._label_of
        for l in labels:
            if l.ref in names and names[l.ref].name != l.name:
                raise CombineRefused('different_index', _T_DIFFERENT_INDEX % (
                    _K.INDEX_CLAUSES['names'] % (l.ref, names[l.ref].name, l.name)))
        return _identity_of(mine, theirs, allow)

    def _check_pairs(self, cparts):
        """different_pair: a part of either side answered on another pair than a part of the
        other (only pairs both state), compared by Pair.identity: by the n1 id where both
        state one or have both paths (an old document and its id-form twin are one pair),
        else by the paths."""
        mine = {answered_pair(p.graphlet.envelope) for p in self._parts} - {None}
        for p in cparts:
            cp = answered_pair(p.graphlet.envelope)
            if cp is None:
                continue
            for pp in sorted(mine, key=lambda x: tuple('' if v is None else v for v in x)):
                if pp.identity != cp.identity:
                    def show(x):
                        # the pair by graph and n1 id, no path (Pair.named())
                        return json.dumps(x.named(), separators=(',', ':'), ensure_ascii=False)
                    raise CombineRefused('different_pair', _T_DIFFERENT_PAIR % (show(cp),
                                                                                show(pp)))

    def _continued_walk(self, side, walk, cseq, b):
        """-> (part, path id, walk key) of the walk the child continues (no_continuation,
        not_a_continuation, ambiguous_walk)."""
        if walk is not None:
            key = self._walk_key(walk)
            owner = self._by_key.get(tuple((side, w) for w in key[:-1]))
            if owner is not None and side not in owner.graphlet.arms:
                raise CombineRefused('not_a_continuation', _T_NOT_A_CONTINUATION_OF % (
                    walk_name(key), side))
            try:
                part, pid = self._part_of_walk(side, key)
            except IndexError as e:
                raise CombineRefused('not_a_continuation', _T_NO_WALK % (side, e)) from None
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
            raise IndexError('the %s arm of part %s has %d walk%s, no walk %d'
                             % (side, part_name(pkey), n, '' if n == 1 else 's', key[-1]))
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

    def _core(self, pg, cg, side, pid, b=None):
        """The junction causes of child |cg| continuing walk |pid| of |pg|'s arm |side| that
        depend on these two parts only (W1-W7 and strategy_unknown; no_bases and the
        positions are the stitched graphlet's): {L, causes: [(code, sentence)], sib: the
        sibling walks of W6 (local ids), continued, not_continued: label refs}. Computed
        once per pair of parts and carried into every stitched graphlet holding them, so
        adding J children one at a time costs linear work; charged to |b| when computed."""
        a = pg.arms[side]
        paths = derive.paths(a)
        if b is not None:
            b.charge(W_ELEM * (len(paths) + len(a.segments)) + 2 * W_ROW,
                     list_bytes(len(paths)) + 2 * (STR + 64))
        leaf = paths[pid].leaf
        L = a.segments[leaf].end_bp
        sp, sc = _strategy(pg), _strategy(cg)
        mode = pg.mode
        causes = []
        if pg.support == 'trace':
            causes.append(('trace_restart', 'under support trace a lineage follows one '
                           'occurrence chain from the root seed; the child validated its seed '
                           'labels against any occurrence of the tail, so it can follow an '
                           'occurrence the parent\'s lineage never had (W1)'))
        # strategy_unknown: a field the rules read missing on either part, or a field one
        # part states and the other does not
        wp, wc = _without_allowed(sp, sc), _without_allowed(sc, sp)
        unknown = []
        for path in _REQUIRED.get(mode, ()):
            hp, hc = _has(sp, path), _has(sc, path)
            if not (hp and hc):
                unknown.append('%s (%s)' % ('.'.join(path), 'both' if not (hp or hc) else
                                            'the parent' if not hp else 'the child'))
        named = {u.split(' ')[0] for u in unknown}
        both_stop = _knob(sp, 'frontier', 'on_overflow') == 'stop' and \
            _knob(sc, 'frontier', 'on_overflow') == 'stop'
        for path in _W2_STATED:
            dotted = '.'.join(path)
            if dotted in named or (both_stop and path == ('bounds', 'max_live_paths')):
                continue
            if not (_has(sp, path) or _has(sc, path)):
                named.add(dotted)
                unknown.append('%s (both)' % dotted)
        for who, x, y in (('the child', wp, wc), ('the parent', wc, wp)):
            for path in _paths_of(x):
                dotted = '.'.join(path)
                mine = x
                for k in path:
                    mine = mine.get(k) if isinstance(mine, dict) else None
                if mine is None:
                    continue                # not stated here (both: an unknown W2 field)
                node = y
                for k in path:
                    node = node.get(k) if isinstance(node, dict) else None
                    if node is None:
                        break
                if node is None and dotted not in named:
                    named.add(dotted)
                    unknown.append('%s (%s)' % (dotted, who))
        diff = _K.first_difference(_common(wp, wc), _common(wc, wp))
        if diff is not None:
            causes.append(('strategy_differs', 'the child\'s normalized strategy differs from '
                           'the parent\'s at %s, a field a continuation keeps (W2)' % diff))
        elif pg.reconverge != cg.reconverge:
            causes.append(('strategy_differs', 'the child\'s normalized strategy differs from '
                           'the parent\'s at branching.on_reconverge (%s in the parent, %s in '
                           'the child, by their H records), a field a continuation keeps (W2, '
                           'W6)' % (pg.reconverge, cg.reconverge)))
        continued, not_cont, not_named, dropped, xnames = (), (), (), (), {}
        if mode == 'constrain':
            if _has(sp, ('labels', 'extra')) and _has(sc, ('labels', 'extra')):
                pool = {l.name for l in pg.labels[:pg.seed.num_seed_labels]} | \
                    set(_knob(sp, 'labels', 'extra', []) or [])
                foreign = [n for n in (_knob(sc, 'labels', 'extra', []) or [])
                           if n not in pool]
                if foreign:
                    causes.append(('extra_foreign', 'the child\'s labels.extra names %s, '
                                   'outside the parent\'s permitted pool (its seed labels and '
                                   'labels.extra) (W2)' % ', '.join(
                                       json.dumps(n, ensure_ascii=False)
                                       for n in foreign[:_SHOWN])
                                   + (' and %d more' % (len(foreign) - _SHOWN)
                                      if len(foreign) > _SHOWN else '')))
            if _has(sp, ('labels', 'change_cost')) and _has(sc, ('labels', 'change_cost')):
                mp = (_knob(sp, 'labels', 'change_cost') or {}).get('model', 'forbid')
                mc = (_knob(sc, 'labels', 'change_cost') or {}).get('model', 'forbid')
                if mp != 'forbid' or mc != 'forbid':
                    causes.append(('switching', 'labels.change_cost is %s in the parent and %s '
                                   'in the child: switching restarts with the child\'s loss '
                                   'budget and source lists (W3 needs forbid on both)'
                                   % (mp, mc)))
            seeds = {l.ref for l in cg.labels[:cg.seed.num_seed_labels]}
            ends = derive.end_labels(a, leaf)
            alive = [pg.labels[e.label] for e in ends]
            continued = tuple(l.ref for l in alive if l.ref in seeds)
            for x in cg.dropped:
                xnames.setdefault(x.name, x.reason)
            dropped = tuple(l for l in alive if l.ref not in seeds and l.name in xnames)
            not_named = tuple(l for l in alive if l.ref not in seeds and l.name not in xnames)
            not_cont = tuple(l.ref for l in alive if l.ref not in seeds)
            if _has(sp, ('branching', 'max_label_branches')) and \
                    _has(sc, ('branching', 'max_label_branches')):
                mlb_p = _knob(sp, 'branching', 'max_label_branches')
                mlb_c = _knob(sc, 'branching', 'max_label_branches')
                if not (mlb_p == 'unlimited' and mlb_c == 'unlimited'):
                    bs = {e.branches for e in ends if pg.labels[e.label].ref in seeds}
                    br = next(iter(bs)) if len(bs) == 1 else (0 if not bs else None)
                    ok = br is not None and isinstance(mlb_p, int) and \
                        isinstance(mlb_c, int) and mlb_c == mlb_p - br
                    if not ok:
                        causes.append(('branches_unequal', 'branching.max_label_branches is %s '
                                       'in the parent and %s in the child, and the continued '
                                       'labels ended the walk at branch counts %s: the '
                                       'allowance does not shift uniformly (W4)'
                                       % (mlb_p, mlb_c, sorted(bs) or '[]')))
        knobs = []
        for who, st in (('parent', sp), ('child', sc)):
            for block, name, bad in _W5_KNOBS:
                v = _knob(st, block, name)
                if v is not None and bad(v):
                    knobs.append('%s.%s %s in the %s' % (block, name, json.dumps(v), who))
        if knobs:
            causes.append(('label_count_knob', 'a knob reads other labels\' counts, which '
                           'restart at the junction: %s (W5)' % ', '.join(knobs)))
        sib = []
        if pg.reconverge == 'merge' and cg.reconverge == 'merge':
            for p in paths:
                if p.id != pid and (a.segments[p.leaf].leaf.continuation is not None
                                    or p.length_bp > L):
                    sib.append(p.id)
        if mode == 'constrain':
            if not_named:
                causes.append(('not_continued', 'the child does not seed %s, alive at the end '
                               'of the walk: their lineages do not go on in the child (W7)'
                               % _label_list(not_named)))
            if dropped:
                causes.append(('seed_dropped', 'the child named %s but dropped %s (X): %s (W7)'
                               % (_label_list(dropped), 'it' if len(dropped) == 1 else 'them',
                                  ', '.join('%s %s' % (l.name, xnames[l.name])
                                            for l in dropped[:_SHOWN]))))
        if unknown:
            causes.append(('strategy_unknown', 'the envelopes do not state %s of the normalized '
                           'strategy: the rules that read them (W2-W5, the split check of W8) '
                           'cannot be judged' % (', '.join(unknown[:_SHOWN])
                                                 + (' and %d more' % (len(unknown) - _SHOWN)
                                                    if len(unknown) > _SHOWN else ''))))
        return {'L': L, 'causes': tuple(causes), 'sib': tuple(sib), 'continued': continued,
                'not_continued': not_cont, 'overlap': len(cg.seed.sequence)}

    def _make_junction(self, part, inherited, b=None):
        """The junction of a child part: its core (computed once per pair of parts, or taken
        from |inherited|), placed in this graphlet: W6's siblings named with their keys and
        position, no_bases when a part from the root to the child carries no bases, and
        exact_gated where no cause applies."""
        key = part.key
        side, pid = key[-1]
        ppart = self._by_key[key[:-1]]
        pg, cg = ppart.graphlet, part.graphlet
        ck = (id(pg), id(cg), side, pid)
        core = None
        for d in inherited:
            core = d.get(ck)
            if core is not None:
                break
        if core is None:
            core = self._core(pg, cg, side, pid, b)
        self._jcore[ck] = core
        if b is not None:
            b.charge(W_ROW + len(key), record_bytes(Junction) + 200)
        causes = list(core['causes'])
        pre = _walk_prefix(ppart.key)
        at = ppart.offset_bp + core['L']
        sib = core['sib']
        if sib:
            causes.append(('sibling_merge_unknown', 'under merge the reference joins paths that '
                           'meet beyond %d bp; walks %s of the arm went on or ended beyond it, '
                           'and the child cannot see them (W6)'
                           % (at, ', '.join(walk_name(pre + (w,)) for w in sib[:_SHOWN])
                              + (' and %d more' % (len(sib) - _SHOWN)
                                 if len(sib) > _SHOWN else ''))))
        bare = [part_name(key[:i]) for i in range(len(key) + 1)
                if not self._bases(self._by_key[key[:i]].graphlet)]
        if bare:
            causes.append(('no_bases', _T_NO_BASES % (
                ('part ' if len(bare) == 1 else 'parts ') + ', '.join(bare),
                'ies' if len(bare) == 1 else 'y')))
        if not causes:
            causes.append(_GATED)
        causes.sort(key=lambda c: _RANK[c[0]])
        return Junction(side, pre + (pid,), at, core['overlap'], key, 'qualified',
                        tuple(causes), tuple(self._label_of[r] for r in core['continued']),
                        tuple(self._label_of[r] for r in core['not_continued']))

    def _bases(self, g):
        got = self._jcore.get(('bases', id(g)))
        if got is None:
            got = self._jcore[('bases', id(g))] = _K.has_bases(g)
        return got

    def _junction_of_part(self, key):
        i = self._jindex[key]
        return i, self._junctions[i]

    @staticmethod
    def _gated(j):
        """Whether junction |j| would be exact once exact is switched on (prospective)."""
        return len(j.causes) == 1 and j.causes[0][0] == 'exact_gated'

    # ---------------------------------------------------------------- derivations

    def _uses(self, b, ck):
        """Charge |b| this graphlet's cached derivation |ck| at its cold price, then the
        derivations it read, each once per call (budget.uses_of(), as derive.uses() charges
        a model's): a hit pays what computing it charges, in the same order, so a call
        charges the same and stops at the same point whatever earlier calls cached. No-op
        without a budget."""
        if b is None:
            return
        work, nbytes, deps = self._price[ck]
        if b.uses_of(self, 'stitched', work, nbytes, depth=ck):
            for d in deps:
                self._uses(b, d)

    def _computing(self, b, ck, price):
        """Charge |b| derivation |ck| at |price| (work, bytes) before it is computed (once
        per call)."""
        if b is not None:
            b.uses_of(self, 'stitched', price[0], price[1], depth=ck)

    def _keep(self, ck, got, price, deps=()):
        """Cache derivation |ck| with its cold price and the derivations it read, in the
        order it read them."""
        self._price[ck] = (price[0], price[1], tuple(deps))
        self._cache[ck] = got
        return got

    # ---------------------------------------------------------------- the rule check (W8)

    def _context(self, key, b=None):
        """The rule check's context of the junction into part |key| (cached): the edge keys
        the reference has used on the continued walk's route (a list of sets, one per part
        of the chain), the root seed's nodes, the child seed's nodes, and the split count
        the reference carries into the child (the counts at the continued leaves, summed
        along the chain)."""
        ck = ('ctx', key)
        got = self._cache.get(ck)
        if got is not None:
            self._uses(b, ck)
            return got
        side, pid = key[-1]
        ppart = self._by_key[key[:-1]]
        pg = ppart.graphlet
        cg = self._by_key[key].graphlet
        a = pg.arms[side]
        leaf = derive.paths(a)[pid].leaf
        route = _ancestors(a, leaf)
        n = sum(a.segments[s].length_bp for s in route) + pg.k * len(route)
        cw, cm = _kmers_price(cg.seed.sequence, cg.k, cg.regime)
        price = (W_ELEM * len(route) + n + cw, set_bytes(n) + n * (STR + pg.k + 1) + cm)
        self._computing(b, ck, price)
        deps = [('tails', ppart.key, side)]
        tails = self._tails(ppart.key, side, b)
        edges = set()
        for s in route:
            seg = a.segments[s]
            window = _before(tails, pg, side, seg) + seg.walk
            for i in range(seg.length_bp):
                edges.add(_ops.edge_key(pg.regime, _natural(side, window[i:i + pg.k + 1]))[0])
        if ppart.key:
            deps.append(('ctx', ppart.key))
            up = self._context(ppart.key, b)
            chain_edges = up['edges'] + (frozenset(edges),)
            count = up['count']
        else:
            chain_edges = (frozenset(edges),)
            count = 0
        count += _split_counts(a)[leaf]
        deps.append(('root_nodes',))
        got = {'edges': chain_edges, 'count': count,
               'root_nodes': self._root_nodes(b),
               'child_nodes': _kmers(cg.seed.sequence, cg.k, cg.regime),
               'max_splits': _knob(_strategy(pg), 'branching', 'max_splits_per_path')}
        return self._keep(ck, got, price, deps)

    def _root_nodes(self, b=None):
        """The root seed's nodes (k-mers, natural orientation), cached."""
        ck = ('root_nodes',)
        got = self._cache.get(ck)
        if got is not None:
            self._uses(b, ck)
            return got
        root = self.root
        price = _kmers_price(root.seed.sequence, root.k, root.regime)
        self._computing(b, ck, price)
        return self._keep(ck, _kmers(root.seed.sequence, root.k, root.regime), price)

    def _tails(self, key, side, b=None):
        """_tails() of part |key|'s arm |side|, cached per part and arm (the root walks
        both)."""
        ck = ('tails', key, side)
        got = self._cache.get(ck)
        if got is not None:
            self._uses(b, ck)
            return got
        g = self._by_key[key].graphlet
        price = _tails_price(g, side)
        self._computing(b, ck, price)
        return self._keep(ck, _tails(g, side), price)

    def _walk_causes(self, key, pid, b=None):
        """The walk-scope causes (§3.4 W8) of walk |pid| of child part |key|, in the child's
        own outward coordinates: [(code, sentence, at_bp local)], the first of each kind on
        the walk's route (the DAG ancestors of its leaf through any parent). Cached."""
        ck = ('w8', key, pid)
        got = self._cache.get(ck)
        if got is not None:
            self._uses(b, ck)
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
        n = sum(segs[s].length_bp for s in route)
        price = (W_ELEM * (len(route) + len(a.branch_events)) + n,
                 set_bytes(len(route)) + 4 * (STR + k + 1))
        self._computing(b, ck, price)
        ctx = self._context(key, b)
        edges, root_nodes, child_nodes = ctx['edges'], ctx['root_nodes'], ctx['child_nodes']
        tails = self._tails(key, side, b)
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
        return self._keep(ck, got, price, (('ctx', key), ('tails', key, side)))

    # ---------------------------------------------------------------- stitched walks

    def _rows(self, side, b=None):
        """Every stitched walk of the arm in key order, as (key, [(part, pid)]), not built.
        Cached; charged at its cold price (_rows_price())."""
        ck = ('rows', side)
        got = self._cache.get(ck)
        if got is not None:
            self._uses(b, ck)
            return got
        price = self._rows_price(side)
        self._computing(b, ck, price)
        out = []

        def visit(part, chain):
            for p in derive.paths(part.graphlet.arms[side]):
                here = chain + [(part, p.id)]
                kid = self._by_key.get(part.key + ((side, p.id),))
                if kid is None:
                    out.append((_walk_prefix(part.key) + (p.id,), here))
                else:
                    visit(kid, here)
        visit(self._parts[0], [])
        return self._keep(ck, out, price)

    def _rows_price(self, side):
        """(lwu, bytes) of _rows(side): W_ELEM per segment of each part on the arm, and per
        walk of each (at most one row each) its row, key, chain list and piece."""
        work = nbytes = walks = 0
        for p in self._parts:
            if p.key and p.key[0][0] != side:
                continue
            a = p.graphlet.arms[side]
            d = len(p.key) + 1
            n = _ops._n_walks(a)
            walks += n
            work += W_ELEM * len(a.segments)
            nbytes += n * (3 * TUPLE + 16 + 8 * d + list_bytes(d))
        return work, nbytes + list_bytes(walks)

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

    def _walk_length(self, side, chain):
        part, pid = chain[-1]
        return part.offset_bp + derive.paths(part.graphlet.arms[side])[pid].length_bp

    @staticmethod
    def _walk_price(n, chain):
        """(lwu, bytes) of building a StitchedWalk of |n| bases over |chain|: its row, its
        pieces and its spelling (three copies of the text at most), as walks() and walk()
        charge each walk they return."""
        return (W_ROW + W_STEP * n // 64 + (n >> 8) + 4 * len(chain),
                record_bytes(StitchedWalk) + 3 * (STR + n) + 400)

    def _walk_object(self, side, row):
        chain = row['chain']
        part, pid = chain[-1]
        g = part.graphlet
        a = g.arms[side]
        # the walk's keys as the part's walks(by='id') ranks them, for this walk alone
        p = derive.paths(a)[pid]
        full, alive_ids, _, min_loss = _ops._walk_keys_of(g, a, p)
        if full is None:
            full = _ops._labels_full(g, a, p)
        w = _ops._walk_objects(g, a, [(p, full, alive_ids, min_loss)], with_claims=False)[0]
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
                    part = chain[-1][0]
                    n = self._walk_length(side, chain)
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
                        b.charge(*self._walk_price(n, chain))
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

    def _chain(self, side, key, b=None):
        """The chain [(part, path id)] of stitched walk |key| (a walk key) on arm |side|."""
        for k, chain in self._rows(side, b):
            if k == key:
                return chain
        raise IndexError('no stitched walk %s on the %s arm (a continued walk is replaced by '
                         'its child\'s walks)' % (walk_name(key), side))

    def _row(self, side, key, prospective=False, b=None):
        key = self._walk_key(key)
        return self._walk_row(side, key, self._chain(side, key, b), prospective, b)

    def walk(self, arm, walk, *, budget=None):
        """One StitchedWalk by key. budget=: the arm's rows it is looked up in, then the walk
        as walks() charges each walk it returns (its row and spelling, its walk causes)."""
        side = self.arm_of(arm)

        def run(b):
            key = self._walk_key(walk)
            chain = self._chain(side, key, b)
            if b is not None:
                b.charge(*self._walk_price(self._walk_length(side, chain), chain))
            return self._walk_object(side, self._walk_row(side, key, chain, False, b))
        return _K.run_scoped(budget, 'walks', run)

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
            if b is not None:
                # the text beyond its pieces: the joined flank, its reversal, the seed (in
                # walking order reversed) and the joined text, at most three copies of the
                # seed and flank (the seed was not charged: GPT final review M2's sibling)
                n = row['length_bp'] + (len(self.seed.sequence) if with_seed else 0)
                b.charge(n >> 8, 3 * (STR + n))
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

    def _walk_keys(self, walks):
        """A walks= argument as walk keys: one key (a path id or a tuple of ints) or an
        iterable of keys -- (1, 0) is the one walk 1.0, ((1, 0), (2,)) two walks."""
        if isinstance(walks, int) and not isinstance(walks, bool) or \
                isinstance(walks, tuple) and walks and all(
                    isinstance(x, int) and not isinstance(x, bool) for x in walks):
            return [self._walk_key(walks)]
        return [self._walk_key(w) for w in walks]

    def _one_part(self, side, walks):
        keys = self._walk_keys(walks)
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
                    for k in self._walk_keys(walks):
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

    def completeness(self, *, budget=None):
        """Per arm: the root's complete_to_bp and scope, unchanged (stitching never makes
        "every walk to depth D" true for a larger D), and one extension per junction on the
        arm, by (arm, walk key): {walk, at_bp, extended_to_bp (at_bp + the child's
        complete_to_bp), certified_to_bp (how far the certificate reaches through the
        junction: no junction is exact in phase 1a-S, so at most at_bp, the continued walk's
        own certified depth), verdict, causes (the junction's cause codes), walk_causes (the
        walk-scope cause codes on the walks below it)}. budget=: the arm's rows, per junction
        its extension and continued walk, and the walk rule (W8) on every walk below it, the
        derivations at their cold price as summary() charges them; a stop raises with no
        partial."""
        return _K.run_scoped(budget, 'completeness', lambda b: self._completeness(False, b))

    def _causes_below(self, j, b=None):
        """The walk-scope cause codes (exact_gated left out) on the stitched walks through
        junction |j|, in Appendix B.2's order: W8 is judged per walk, so a junction whose own
        rules hold can still stand above a walk the reference walks otherwise. |b|: the scan
        of the arm's rows, each walk below (its row) and its derivations."""
        codes = set()
        rows = self._rows(j.arm, b)
        if b is not None:
            b.charge(W_ELEM * len(rows))
        for key, chain in rows:
            if key[:len(j.walk)] == j.walk and len(key) > len(j.walk):
                if b is not None:
                    b.charge(W_ROW + W_ELEM * len(chain))
                codes.update(c[0] for c in self._walk_row(j.arm, key, chain, False, b)['causes']
                             if c[0] != 'exact_gated')
        return sorted(codes, key=_RANK.get)

    def _chain_of(self, side, key):
        """[(part, path id)] of walk key |key| on arm |side|, root first."""
        return [self._part_of_walk(side, key[:i + 1]) for i in range(len(key))]

    def _completeness(self, prospective, b=None):
        out = {}
        for side in self.arms():
            a = self.root.arms[side]
            exts = []
            for j in self._junctions:
                if j.arm != side:
                    continue
                cg = self._by_key[j.part].graphlet
                ext = j.at_bp + cg.arms[side].complete_to_bp
                chain = self._chain_of(side, j.walk)
                if b is not None:
                    b.charge(W_ROW + W_ELEM * len(chain),
                             dict_bytes(7) + 2 * list_bytes(len(j.causes)))
                cert = min(j.at_bp, self._certify(side, chain, prospective, b)[0])
                if prospective and self._gated(j) and cert == j.at_bp:
                    # every walk below the junction certified to the extension, or as far
                    # as the first one that is not complete
                    cert = ext
                    for key, chain in self._rows(side, b):
                        if key[:len(j.walk)] == j.walk:
                            row = self._walk_row(side, key, chain, True, b)
                            if not row['complete']:
                                cert = min(cert, row['certified_to_bp'])
                exts.append({'walk': j.walk, 'at_bp': j.at_bp, 'extended_to_bp': ext,
                             'certified_to_bp': cert,
                             'verdict': ('exact' if self._gated(j) else 'qualified')
                             if prospective else j.verdict,
                             'causes': [c[0] for c in j.causes
                                        if not (prospective and c[0] == 'exact_gated')],
                             'walk_causes': self._causes_below(j, b)})
            out[side] = {'complete_to_bp': a.complete_to_bp, 'scope': a.scope,
                         'extensions': exts}
        return out

    # ---------------------------------------------------------------- evidence

    def evidence(self, arm=None, walks=None, *, budget=None):
        """The evidence of an answer (§3.5): {parts: each part's own evidence_block keyed by
        part name ('root', 'right:3', ...), junctions: those the answer reads with their
        verdicts and causes, walk_causes: the walk-scope causes on the walks it reads,
        identity, outcome: the weakest value per dimension over the parts and junctions read
        (a qualified junction or a walk-scope cause: walks partial, label_evidence
        qualified), exact}. |arm|: the arm's parts and junctions (None: all); |walks|: only
        the parts and junctions those walks (keys, on |arm|) read -- a root walk no child
        continues reads the root alone and states its evidence. exact holds only when every
        part read is exact, no junction is read (no junction is exact in phase 1a-S) and no
        walk read has a walk-scope cause. budget=: the rows of the walks read, the walk rule
        (W8) on each, the derivations at their cold price as summary() charges them, each
        part's evidence block and each junction and cause stated; a stop raises with no
        partial."""
        return _K.run_scoped(budget, 'evidence', lambda b: self._evidence(arm, walks, False, b),
                             ('select_walks', 'narrow_arm'))

    def _evidence(self, arm, walks, prospective, b=None):
        """evidence(); |prospective|: the private test hook of §3.5 (exact_gated read as
        exact; where a junction is read, exact needs the identity verified)."""
        sides = [self.arm_of(arm)] if arm is not None else self.arms()
        if walks is not None and arm is None and len(sides) > 1:
            raise ValueError('walks= names walks of one arm: give arm=')
        read_parts, read_j, wcauses = set(), set(), []
        for side in sides:
            if walks is None:
                chosen = self._rows(side, b)
            else:
                chosen = [(k, self._chain(side, k, b)) for k in self._walk_keys(walks)]
            for key, chain in chosen:
                for part, _ in chain:
                    read_parts.add(part.key)
                if b is not None:
                    b.charge(W_ROW + W_ELEM * len(chain))
                row = self._walk_row(side, key, chain, prospective, b)
                read_j.update(row['junctions'])
                if b is not None and row['causes']:
                    b.charge(W_ELEM * len(row['causes']), len(row['causes']) * dict_bytes(4))
                for c in row['causes']:
                    if c[0] != 'exact_gated':
                        wcauses.append({'walk': key, 'code': c[0], 'at_bp': c[2],
                                        'sentence': c[1]})
        parts = [p for p in self._parts if p.key in read_parts]
        verified = bool(self.identity.get('verified')) or not read_j
        if b is not None:
            # each part's evidence block and the limitations it reads; each junction
            # stated with its causes
            nlim = sum(len(p.graphlet.limitations) for p in parts)
            jbytes = sum(dict_bytes(5) + 2 * list_bytes(len(self._junctions[i].causes))
                         for i in read_j)
            b.charge(W_ROW * (len(parts) + len(read_j)) + W_ELEM * nlim,
                     len(parts) * _EV_BLOCK + nlim * (STR + 16) + jbytes)
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
        without: {completeness: the arm's, evidence: the arm's}. Uncharged: it takes no
        budget and ignores the ambient one (a test hook, not an operation)."""
        side = self.arm_of(arm)
        if walk is not None:
            return self._row(side, walk, True)
        return {'completeness': self._completeness(True)[side],
                'evidence': self._evidence(side, None, True)}

    # ---------------------------------------------------------------- summary, provenance

    def summary(self, max_bytes=2048, *, budget=None):
        """Agent-facing, <= max_bytes of JSON: the root seed, per arm the root's certificate
        and the stitched walks' count and longest, the junctions with their verdicts and
        cause codes and the walk-scope cause codes on the walks below each (walk_causes),
        the identity, the outcome and exact of evidence(); marked experimental. Junctions
        beyond the cap are cut from the end (junctions_cut); the rest is never cut:
        summary_too_small when it alone exceeds max_bytes, with that least size in its text
        and as .least_bytes (summary(max_bytes=e.least_bytes) answers)."""
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
                if b is not None:
                    b.charge(W_ELEM * len(rows))
                lengths = [self._walk_length(side, chain) for _, chain in rows]
                out['arms'][side] = {'complete_to_bp': a.complete_to_bp, 'scope': a.scope,
                                     'walks': len(rows), 'max_bp': max(lengths, default=0)}
            for j in self._junctions:
                out['junctions'].append({'arm': j.arm, 'walk': walk_name(j.walk),
                                         'at_bp': j.at_bp, 'verdict': j.verdict,
                                         'causes': [c[0] for c in j.causes],
                                         'walk_causes': self._causes_below(j, b)})
            ev = self._evidence(None, None, False, b)
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
            n = size()
            if n > max_bytes:
                # the least answer's size, as the union's refusal carries it (least_bytes)
                e = CombineRefused('summary_too_small', _T_SUMMARY_TOO_SMALL % (n, max_bytes, n))
                e.least_bytes = n
                raise e
            return out
        return _K.run_scoped(budget, 'summary', run)

    def provenance(self):
        """One entry per part (§3.6): {key, role: root | child, junction: {arm, walk, at_bp,
        overlap_bp} | None, pair: the pair it was answered on ({graph, pair_id}, no path;
        None on a single-graph server), identity, seed_id, digest, request: the part's replay
        request as GraphletStore.request_of() builds it, its pair named by {graph, pair_id}
        (client.request_without_paths())}. No single request reproduces a stitched graphlet:
        replay each part and stitch again."""
        from .client import request_without_paths
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
            pair = answered_pair(g.envelope)
            out.append({'key': p.key, 'role': 'child' if p.key else 'root', 'junction': j,
                        'pair': None if pair is None else pair.named(),
                        'identity': _ops.index_identity(g),
                        'seed_id': meta.get('seed_id'), 'digest': p.digest,
                        'request': request_without_paths(GraphletStore.request_of(g))})
        return out

    def memory_bytes(self):
        """What this graphlet holds: the parts' models and everything it built from them --
        its parts and labels, the junctions and their cores, and every cache of the queries
        (rule-check indexes, walk causes, rows, tails, root nodes) with its price -- each
        object once."""
        return _ops._deep_bytes([p.graphlet for p in self._parts]
                                + [self._parts, self._by_key, self._labels, self._label_of,
                                   self._junctions, self._jindex, self._jcore, self._cache,
                                   self._price, self.identity])

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

# a label of a new StitchedGraphlet: its Label, its ref (a key of the by-ref index, made
# twice while the labels are joined), its list and tuple slots and dict entry
_M_LABEL = 320


def _build_price(graphlets):
    """(lwu, bytes) of building a StitchedGraphlet over the parts |graphlets|, beyond its
    index: the labels joined (per label of each part an element and _M_LABEL) and the
    digest (GPT final review M3): per part its retrieval identity -- its case-mapped seed,
    kept -- and in the digest's canonical JSON its strategy, seed and seed labels at their
    JSON bound, four copies (the encoder's chunks, the text, its UTF-8 bytes and the next
    part's identity), and 1 lwu per 8 bytes of it."""
    from .union import _json_bound     # union.py imports this module
    nl = sum(len(g.labels) for g in graphlets)
    meta = seeds = 0
    for g in graphlets:
        n = len(g.seed.sequence)
        seeds += STR + n
        meta += 600 + n + _json_bound((g.envelope or {}).get('strategy') or {}) + 32 * min(
            g.seed.num_seed_labels, len(g.labels))
    return (W_ELEM * nl + (meta >> 3) + (seeds >> 6),
            _M_LABEL * nl + list_bytes(nl) + dict_bytes(nl) + seeds + 4 * (STR + meta))


def _index_key(g):
    i = _ops.index_identity(g)
    return (i['ns'], i['fp'], i['meta_fp'], i['release'])


def _identity_of(mine, theirs, allow):
    """The identity statement of the parts |mine| (the parent's) and |theirs| (the child's),
    no index difference between them proven: verified when every part states the index
    manifest digest (then one digest); unverified, and not refused, when none does (OQ6);
    otherwise index_unverifiable unless |allow| (allow_unverified_index), which the
    statement then records."""
    every = list(mine) + list(theirs)
    fps = [g.index_fp for g in every]
    if all(fps):
        return {'verified': True, 'index_fp': fps[0]}
    if not any(fps):
        return {'verified': False,
                'note': 'no part states an index manifest digest (index_fp): that the parts '
                        'are one index cannot be verified'}
    if not allow:
        nm = sum(1 for g in mine if g.index_fp)
        nt = sum(1 for g in theirs if g.index_fp)
        if nm in (0, len(mine)) and nt in (0, len(theirs)):
            raise CombineRefused('index_unverifiable', _T_INDEX_UNVERIFIABLE % (
                'an' if nt else 'no', 'one' if nm else 'none'))
        raise CombineRefused('index_unverifiable', _T_INDEX_UNVERIFIABLE_PARTS % (
            nm, len(mine), nt, len(theirs)))
    return {'verified': False, 'accepted': 'allow_unverified_index',
            'note': 'some parts state an index manifest digest (index_fp) and others none: '
                    'whether they are one index cannot be verified (not proven different); '
                    'stitched because the caller passed allow_unverified_index'}
