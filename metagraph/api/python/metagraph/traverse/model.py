"""The in-memory graphlet: slotted dataclasses mirroring the MGT v1 records
(DESIGN-traverse-graphlet.md §2.2, §5).

Ids are the format's ordinals and never change: label ids index Graphlet.labels, segment
and run ids index Arm.segments / Arm.runs. Label sets are interned array('I') (shared
between equal sets) and must be treated as immutable. Fields the format derives
("*") are resolved by the parser, so every field here holds a value; which ones the
source wrote as "*" is kept in Segment.stars / Run.stars for check_rules().

The operations live in derive.py, ops.py and export.py; Graphlet's methods delegate to
them (imported lazily: model must not import its users).
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

__all__ = [
    'Label', 'SeedInfo', 'Dropped', 'Event', 'PresenceRun', 'LeafContinuation', 'Leaf',
    'Segment', 'Run', 'GrowthBin', 'Refusal', 'BranchEvent', 'CapTrigger', 'Frontier',
    'LabelsPerNode', 'Counts', 'Limitation', 'Outcome', 'ResourceStop', 'Arm', 'Graphlet',
    'Walk', 'Claim', 'LabelWalk', 'SupportRun', 'Change', 'Continuation', 'Comparison',
    'CoordClaim', 'CoordWalk', 'CoordLabelWalk',
    'Split', 'SplitPoint', 'Branch', 'Path', 'EndLabel', 'MissingEnvelope',
    'AmbiguousLabel', 'UnknownLabel', 'BadSelector', 'IncompleteRecording',
    'UnverifiableLabelName', 'ARM_SIDES',
]

ARM_SIDES = ('left', 'right')


class MissingEnvelope(LookupError):
    """An operation that needs the response envelope (seed_id, the normalized strategy,
    annotation counters) on a graphlet parsed from a body alone (§5): raised instead
    of inventing defaults."""


class AmbiguousLabel(LookupError):
    """A bare string matched no label or several (a column may be NAMED 'c:0'; annotate
    mode can record two headers 'ACC1' from different columns): use a tagged selector."""


class UnknownLabel(AmbiguousLabel):
    """No label matches the selector (a subclass of AmbiguousLabel for the library's
    callers; the tool layer reports it as its own error)."""


class BadSelector(TypeError, ValueError):
    """A malformed label selector or walk argument -- not a selector that matches
    nothing. A TypeError for the library's callers, and a ValueError, so that the tool
    layer reports it as a bad argument instead of failing."""


class IncompleteRecording(RuntimeError):
    """claims(strict=True) on an annotate arm whose recorded lists were cut: the oracle
    filtered from them would be a lower bound, not the label-consistent walks."""


class UnverifiableLabelName(ValueError):
    """A request would name labels whose names the library cannot verify to resolve back
    to exactly these labels: /traverse resolves a name to a label, and two labels of the
    retrieval sharing it, or a name with U+FFFD (what a server before the refusal of
    names that are not UTF-8 wrote for bytes it could not carry, so the index's own
    name may be another), could constrain the new traversal to another label. Raised
    instead of building the request (Continuation.as_seed, next_request)."""


class IncompatibleContinuations(ValueError):
    """next_request() over several walks whose continuations carry different labels: one
    request has one strategy, so one labels.extra, and the permitted pool around each
    seed (its continuation's labels, every other label of the retrieval a switch target)
    differs per walk -- a shared list would duplicate one seed's labels or drop another
    seed's switch targets. Raised instead of building it: next_requests() builds one
    request per walk."""


# ------------------------------------------------------------------------ records

@dataclass(slots=True)
class Label:
    id: int
    kind: str                 # 'column' | 'header'
    column: int
    seq_id: Optional[int]     # header labels only
    name: str

    @property
    def ref(self):
        """The stable cross-retrieval join key (LabelRef); names are for display."""
        if self.kind == 'header':
            return 'h:%d:%d' % (self.column, self.seq_id)
        return 'c:%d' % self.column

    def as_dict(self):
        return {'name': self.name, 'ref': self.ref}


@dataclass(slots=True)
class SeedInfo:
    validated_seed_id: str
    length_bp: int
    num_kmers: int
    num_seed_labels: int
    sequence: str


@dataclass(slots=True)
class Dropped:
    reason: str
    runs: List[Tuple[int, int]]
    name: str


@dataclass(slots=True)
class Event:
    """A stored event (E record). label_end and reconverge events are derived (from R
    and from G) and are not Events here."""
    at_bp: int
    type: str                           # switch | blocked | hairpin | revisit | tip | bubble
    from_label: Optional[int] = None    # switch
    to_label: Optional[int] = None      # switch
    cost: Optional[float] = None        # switch
    char: Optional[str] = None          # blocked, hairpin, tip
    reason: Optional[str] = None        # blocked: end-reason code
    total: Optional[int] = None         # blocked, hairpin: true label count
    labels: Any = None                  # blocked, hairpin: array('I')
    followed: Optional[bool] = None     # hairpin
    segment: Optional[int] = None       # revisit: the segment first seen
    delta: Optional[int] = None         # revisit: distance delta, None = same distance
    length_bp: Optional[int] = None     # tip, bubble
    alleles: Optional[str] = None       # bubble


@dataclass(slots=True)
class PresenceRun:
    """Annotate mode: the labels present on the nodes entered by steps from_bp+1..to_bp."""
    from_bp: int
    to_bp: int
    total: int
    labels: Any                          # array('I'), cut at the cap

    @property
    def truncated(self):
        return self.total > len(self.labels)


@dataclass(slots=True)
class LeafContinuation:
    """The C record: the walker's label logic (primary); the spelling is derived from
    the bases unless the document carries none (sequence then explicit)."""
    n: int
    loss_used: float
    branches_used: int
    labels: Any                          # array('I')
    sequence: Optional[str] = None       # explicit only when G carries no bases


@dataclass(slots=True)
class Leaf:
    path_reason: Optional[str]           # end-reason code, None = semantic end
    extras: Dict[int, Tuple[float, int, int]]   # label -> (loss, branches, route_bp)
    continuation: Optional[LeafContinuation] = None


@dataclass(slots=True)
class Segment:
    id: int
    parents: Tuple[int, ...]
    from_bp: int
    length_bp: int
    walk: Optional[str]                  # bases in WALKING order; None: sequences false
    first_base: Optional[str]
    entry: Any                           # array('I'): labels at the entry node
    entry_total: int
    end: Any                             # array('I'): labels at the last node
    partition: List[Any]                 # per parent (merges); [] otherwise
    split: Optional[bool]                # Split::ambiguous of a split here; None: no split
    presence: List[PresenceRun] = field(default_factory=list)
    events: List[Event] = field(default_factory=list)
    leaf: Optional[Leaf] = None
    children: List[int] = field(default_factory=list)
    depth: int = 0                       # segments on the first-parent chain above
    stars: frozenset = field(default=frozenset(), compare=False, repr=False)

    @property
    def end_bp(self):
        return self.from_bp + self.length_bp

    @property
    def is_merge(self):
        return len(self.parents) > 1

    @property
    def is_leaf(self):
        return self.leaf is not None


@dataclass(slots=True)
class Run:
    id: int
    segment: int                         # the anchor (§4)
    label: int
    from_bp: int
    to_bp: int
    end: str                             # the R end token: 'D', 'Rm', 'Lw', 'm', ...
    reason: Optional[str]                # enum name ('label_lost'); None for 'm'
    qualifier: Optional[str]             # walker text ('minority'); None without
    silent: bool                         # Lw: switched away, no label_end event
    merged: bool                         # m: closed by a merge, not ended
    route_bp: int                        # the EARLIEST merge stamp (0: none)
    from_label: Optional[int]            # None: entered by seed
    cost: Optional[float]
    prev_run: Optional[int]
    structural_successors: Optional[int]
    branches: int                        # terminal values, valid at to_bp only
    loss: float
    needed_budget: Optional[float] = None
    stars: frozenset = field(default=frozenset(), compare=False, repr=False)

    @property
    def code(self):
        """The end-reason code letter, None for a merge-closed run."""
        return None if self.merged else self.end[0]

    @property
    def entered_by(self):
        return 'seed' if self.from_label is None else 'switch'

    @property
    def ended(self):
        return not self.merged

    @property
    def event_reason(self):
        """The label_end event's reason as today's JSON writes it: the qualifier text
        when there is one, the enum name otherwise."""
        return self.qualifier or self.reason


@dataclass(slots=True)
class GrowthBin:
    from_bp: int
    max_live_paths: int
    distinct_live_labels: int
    live_pairs: int
    exact: bool
    steps: int
    divergences: int
    ambiguous_branches: int
    splits: int
    reconvergences: int
    bubbles: int
    tips: int
    blocked_repeat: int
    label_ends: Dict[str, int]           # end-reason code -> count


@dataclass(slots=True)
class Refusal:
    char: str
    cause: str
    labels: Any


@dataclass(slots=True)
class BranchEvent:
    at_bp: int
    segment: int
    chars: str
    labels_per_successor: List[int]
    ambiguous: Any
    dropped: Any
    refused: List[Refusal]


@dataclass(slots=True)
class CapTrigger:
    reason: str                          # end-reason code
    at_bp: int
    segment: int
    live_paths: int
    live_labels: int
    exact: bool
    demand: float


@dataclass(slots=True)
class Frontier:
    live_paths: int
    live_labels: int
    exact: bool


@dataclass(slots=True)
class LabelsPerNode:
    max_seen: int
    nodes_truncated: int


@dataclass(slots=True)
class Counts:
    """The last six A fields: they validate the body (a truncated arm fails them)."""
    segments: int
    runs: int
    leaves: int
    splits: int
    merges: int
    bases: int


@dataclass(slots=True)
class Limitation:
    """One K record: a stated limitation (spec §7.0, design §14). |arm| None = seed
    level. Values are typed: int, float, str or UNLIMITED; |extra| keeps its order."""
    arm: Optional[str]
    kind: str
    knob: str
    limit: Any
    observed: Any
    complete_to_bp: Optional[int]
    extra: List[Tuple[str, Any]]
    effect: str


@dataclass(slots=True)
class Outcome:
    walks: str                           # complete | partial | failed
    branch_diagnostics: str              # complete | cut
    label_evidence: str                  # complete | lower_bound | qualified
    delivery: str                        # inline | spooled | paged

    def as_dict(self):
        return {'walks': self.walks, 'branch_diagnostics': self.branch_diagnostics,
                'label_evidence': self.label_evidence, 'delivery': self.delivery}


@dataclass(slots=True)
class ResourceStop:
    scope: str
    resource: str
    phase: str
    requested: Any
    effective: Any
    used: Any
    remaining: Any
    actions: List[str]
    message: str


@dataclass(slots=True)
class Arm:
    side: str                            # 'left' | 'right'
    status: str                          # complete | truncated | pruned
    complete_to_bp: int
    scope: str                           # per_path | united_history
    frontier: Frontier
    labels_per_node: LabelsPerNode
    cap_trigger: Optional[CapTrigger]
    counters: Dict[str, int]             # in the writer's order; unknown names kept
    growth: List[GrowthBin]
    branch_events: List[BranchEvent]
    branch_events_total: int
    evidence_complete_to_bp: Optional[int]   # None: no branch event dropped
    counts: Counts
    segments: List[Segment] = field(default_factory=list)
    runs: List[Run] = field(default_factory=list)
    cache: Dict[str, Any] = field(default_factory=dict, compare=False, repr=False)

    @property
    def code(self):
        return self.side[0]


@dataclass(slots=True)
class Graphlet:
    format: int
    k: int
    regime: str
    alphabet: str
    mode: str                            # constrain | annotate
    support: str                         # kmer | trace
    reconverge: str                      # merge | keep
    cap: int                             # labels.max_labels_per_node
    continuation_bp: int
    seed_index: int
    index_ns: Optional[str]
    index_fp: Optional[str]
    index_meta_fp: Optional[str]
    seed: SeedInfo
    labels: List[Label]
    dropped: List[Dropped]
    outcome: Outcome
    resource_stop: Optional[ResourceStop]
    limitations: List[Limitation]
    arms: Dict[str, Arm]                 # left before right, requested arms only
    # results[i] of the response without the body: the per-seed JSON summary (§3). Named
    # seed_summary because summary() is the agent-facing view of the public API (§5)
    seed_summary: Optional[dict] = None
    envelope: Optional[dict] = None      # the response without results
    view: Optional[dict] = None          # J {view: ...} of a saved GraphletView
    derived_from: Optional[str] = None
    has_j: bool = field(default=False, compare=False, repr=False)
    cache: Dict[str, Any] = field(default_factory=dict, compare=False, repr=False)
    # the sha-256 of the parsed text's records without its J line and Z record (parse()):
    # what a resume token of a stopped list binds (local limits), so that a token is refused on
    # any other graphlet and accepted on the same body parsed again or saved and loaded
    body_digest: Optional[str] = field(default=None, compare=False, repr=False)

    # ---------------------------------------------------------------- text

    # every operation takes budget= (local limits, metagraph.traverse.budget): keyword-only,
    # None (the default) for no work or allocation budget; the lists whose order allows
    # it also take resume= (the token of a stopped call's Partial)

    def dump(self, envelope=None, *, budget=None):
        """Canonical MGT text. envelope None: a J line iff the graphlet was read with
        one; True: whenever an envelope is attached; False: the body only."""
        from . import parser
        return parser.dump(self, envelope=envelope, budget=budget)

    def save(self, path, *, budget=None):
        from . import parser
        return parser.save(self, path, budget=budget)

    @classmethod
    def load(cls, path, *, budget=None):
        from . import parser
        return parser.load(path, budget=budget)

    @classmethod
    def parse(cls, text, *, budget=None):
        from . import parser
        return parser.parse(text, budget=budget)

    @classmethod
    def from_response(cls, result, response, *, budget=None):
        from . import parser
        return parser.from_response(result, response, budget=budget)

    # ---------------------------------------------------------------- lookup

    def arm(self, side=None):
        """'left' | 'right' | 'l' | 'r'; None when only one arm was requested."""
        if isinstance(side, Arm):
            return side
        if side is None:
            if len(self.arms) != 1:
                raise ValueError('both arms are present: name one (left or right)')
            return next(iter(self.arms.values()))
        s = {'l': 'left', 'r': 'right'}.get(side, side)
        try:
            return self.arms[s]
        except KeyError:
            raise KeyError('arm %r was not requested in this retrieval' % side) from None

    def label(self, selector):
        from . import ops
        return ops.label(self, selector)

    def labels_matching(self, selectors):
        from . import ops
        return ops.labels_matching(self, selectors)

    @property
    def has_envelope(self):
        return self.seed_summary is not None and self.envelope is not None

    @property
    def coordinates(self):
        """The record coordinates of the retrieval (coords.Coordinates: where each run's
        own bases lie in the indexed records, DESIGN §18), or None -- coordinates_reason
        says why. Parsed from the per-seed summary and validated eagerly against the body
        (GraphletFormatError when inconsistent), kept in the graphlet's cache: a property,
        not a slot, so a graphlet without coordinates carries no extra field."""
        from . import coords
        return coords.of(self)

    @property
    def coordinates_reason(self):
        """Why coordinates is None: 'no envelope' (a parsed body alone carries none),
        'not requested' (the request did not set strategy.output.coordinates), or the
        server's reason ('index has no coordinates', 'support kmer', 'no traversal',
        'partial derivation'); None when the retrieval has them."""
        from . import coords
        return coords.reason(self)

    def require_envelope(self, what):
        if not self.has_envelope:
            raise MissingEnvelope(
                '%s needs the response envelope (seed_id, the normalized strategy, '
                'annotation counters): build the graphlet with Graphlet.from_response() or '
                'load a saved .mgt file with its J line; a parsed body alone does not carry '
                'it' % what)

    # ---------------------------------------------------------------- operations

    def to_json(self, *, budget=None):
        from . import export
        return export.to_json(self, budget=budget)

    def summary(self, arm=None, max_bytes=2048, *, budget=None):
        """Agent-facing, <= max_bytes (2 KB): per-arm status, completeness, counts, top
        labels, caveats. Raises MissingEnvelope on a body-only graphlet."""
        from . import ops
        return ops.summary(self, arm, max_bytes, budget=budget)

    def spell(self, arm, leaf, orientation='natural', with_seed=False, *, budget=None):
        from . import ops
        return ops.spell(self, arm, leaf, orientation, with_seed, budget=budget)

    def walks(self, arm, *, top=None, by='support', labels=None, route_consistent=True,
              min_bp=0, budget=None, resume=None):
        from . import ops
        return ops.walks(self, arm, top=top, by=by, labels=labels,
                         route_consistent=route_consistent, min_bp=min_bp, budget=budget,
                         resume=resume)

    def claims(self, arm=None, labels=None, at_most_bp=None, strict=True, *, budget=None,
               resume=None):
        from . import ops
        return ops.claims(self, arm, labels, at_most_bp, strict, budget=budget,
                          resume=resume)

    def label_walks(self, label, arm=None, *, budget=None, resume=None):
        from . import ops
        return ops.label_walks(self, label, arm, budget=budget, resume=resume)

    def routes(self, label, arm, spell=False, *, budget=None, resume=None):
        from . import ops
        return ops.routes(self, label, arm, spell, budget=budget, resume=resume)

    def support_profile(self, arm, leaf, kind='displayed', *, budget=None):
        from . import ops
        return ops.support_profile(self, arm, leaf, kind, budget=budget)

    def support_changes(self, arm, leaf, *, budget=None):
        from . import ops
        return ops.support_changes(self, arm, leaf, budget=budget)

    def label_summary(self, *, budget=None):
        from . import ops
        return ops.label_summary(self, budget=budget)

    def continuation(self, arm, leaf, *, budget=None):
        from . import ops
        return ops.continuation(self, arm, leaf, budget=budget)

    def next_request(self, arm, leaves, bp=None, reduce_budget=True, reset_branches=False,
                     **overrides):
        """A NextRequest continuing |leaves| (see ops.next_request): .notes states what
        the request cannot carry exactly (a conservative loss budget or branch allowance,
        a switch target left out, an allowance the caller restarted). budget= is a
        LocalBudget (local limits), not the request's loss budget; any other value of it
        is a keyword override like the others (passed on with them)."""
        from . import ops
        return ops.next_request(self, arm, leaves, bp, reduce_budget, reset_branches,
                                **overrides)

    def next_requests(self, arm, leaves, bp=None, reduce_budget=True, reset_branches=False,
                      **overrides):
        """One NextRequest per walk (walks whose continuations carry different labels
        cannot share one request: IncompatibleContinuations). budget= as in
        next_request()."""
        from . import ops
        return ops.next_requests(self, arm, leaves, bp, reduce_budget, reset_branches,
                                 **overrides)

    def subgraph(self, selectors, arm=None, mode='any', *, budget=None, of=None):
        from . import ops
        return ops.subgraph(self, selectors, arm, mode, budget=budget, of=of)

    def to_fasta(self, arm=None, leaves=None, with_seed=True, orientation='natural',
                 width=None, *, coordinates=None, budget=None):
        """FASTA, one record per walk (export.to_fasta); |coordinates|: None states the
        record coordinates in the headers when the retrieval has them, False never,
        True requires them."""
        from . import export
        return export.to_fasta(self, arm, leaves, with_seed, orientation, width,
                               coordinates=coordinates, budget=budget)

    def to_gfa(self, with_seed=True, *, budget=None):
        from . import export
        return export.to_gfa(self, with_seed, budget=budget)

    def splits(self, arm, min_labels_before=0, *, budget=None, resume=None):
        """The arm's splits (SplitPoint: at_bp outward, kind, labels_before, branches with
        their first base, labels as {name, ref} and the first walk taking each), in the
        walker's order."""
        from . import ops
        return ops.splits(self, arm, min_labels_before, budget=budget, resume=resume)

    def compare(self, other, *, arm=None, labels=None, mode='claims', budget=None):
        from . import ops
        return ops.compare(self, other, arm=arm, labels=labels, mode=mode, budget=budget)

    def compare_cost(self, other, *, arm=None, labels=None, mode='claims'):
        """An estimate of what compare() would charge a budget (ops.compare_cost)."""
        from . import ops
        return ops.compare_cost(self, other, arm=arm, labels=labels, mode=mode)

    def memory_bytes(self):
        from . import ops
        return ops.memory_bytes(self)

    def check_rules(self, reference=None):
        from . import derive
        return derive.check_rules(self, reference)


# ------------------------------------------------------------------ query results

@dataclass(slots=True)
class EndLabel:
    label: int
    loss: float
    branches: int
    run: int
    route_bp: int


@dataclass(slots=True)
class Path:
    """paths[i]: the leaf and its length; the root -> leaf chain is produced on access
    (a comb-shaped trie would make materialising every chain quadratic)."""
    id: int
    leaf: int
    length_bp: int
    arm: Any = field(default=None, repr=False, compare=False)

    @property
    def segments(self):
        from . import derive
        return derive.chain(self.arm, self.leaf)


@dataclass(slots=True)
class Split:
    at_bp: int
    segment: int
    children: List[int]
    ambiguous: bool
    labels_before: int


@dataclass(slots=True)
class Branch:
    """One branch of a split: its first segment, the base it starts with (the base at
    outward position at_bp of the split, the same character in either orientation:
    the left arm is read reversed, not complemented), the labels at its first node (the
    recorded list, cut at labels.max_labels_per_node) with their true count, and the
    smallest walk (path id) whose displayed chain takes it (None: none does)."""
    segment: int
    char: Optional[str]
    labels: List[Label]
    labels_total: int
    walk: Optional[int]


@dataclass(slots=True)
class SplitPoint:
    """A split of an arm as Graphlet.splits() returns it: at_bp is the outward distance
    from the seed boundary (walking order, like a claim's from_bp/to_bp on both arms),
    kind 'ambiguous' (some lineage continued on more than one branch) or 'divergence',
    labels_before the labels at the parent's last node (true count)."""
    arm: str
    at_bp: int
    segment: int
    kind: str
    labels_before: int
    branches: List[Branch]

    @property
    def ambiguous(self):
        return self.kind == 'ambiguous'


@dataclass(slots=True)
class Claim:
    """The §6.9 unit: one per run end (claims, not leaves)."""
    label: Label
    arm: str
    path_id: Optional[int]
    segment: int
    from_bp: int
    to_bp: int
    route_bp: int                        # merge-derived start of displayed support
    evidence_from: Optional[int]         # None for route_only / boundary
    kind: str                            # end alive diverged merged stretch route_only boundary
    reason: Optional[str]
    qualifier: Optional[str]
    end_class: str                       # lost dead_end radius capped blocked pruned open
    entered_by: str
    from_label: Optional[Label]
    cost: Optional[float]
    loss: Optional[float]                # None when cut before to_bp (unknown there)
    branches: Optional[int]
    support: str
    exact: bool
    run: Optional[int] = None

    def as_dict(self):
        return {
            'label': self.label.as_dict(), 'arm': self.arm, 'path_id': self.path_id,
            'from_bp': self.from_bp, 'to_bp': self.to_bp, 'route_bp': self.route_bp,
            'evidence_from': self.evidence_from, 'kind': self.kind, 'reason': self.reason,
            'qualifier': self.qualifier, 'end_class': self.end_class,
            'entered_by': self.entered_by,
            'from_label': self.from_label.as_dict() if self.from_label else None,
            'cost': self.cost, 'loss': self.loss, 'branches': self.branches,
            'support': self.support, 'exact': self.exact,
        }

    # a claim of a retrieval without record coordinates has none: a class attribute, not
    # a slot, so a Claim's size (the local-limits account of every claims call) and its
    # fields (the golden digests) do not depend on coordinates; CoordClaim carries them
    coordinates = None


@dataclass(slots=True)
class CoordClaim(Claim):
    """A claim of a retrieval with record coordinates: |coordinates| is its run's
    coords.RunCoordinates clipped to the claim (the run's own bases [from_bp, to_bp), the
    cut included; a route_only claim's route too). Made only for graphlets that carry a
    coordinates block. Note: a CoordClaim never equals a Claim (dataclass equality compares
    the types), which no comparison of the library relies on."""
    coordinates: Any = None

    def as_dict(self):
        d = Claim.as_dict(self)
        if self.coordinates is not None:
            d['coordinates'] = self.coordinates.as_dict(
                'record' if self.label.kind == 'header' else 'column')
        return d


@dataclass(slots=True)
class Walk:
    path_id: int
    leaf: int
    segments: List[int]
    length_bp: int
    sequence: Optional[str]
    path_reason: Optional[str]
    end_reasons: Dict[str, int]
    claims: List[Claim]
    labels_full: List[Label]
    n_alive: int
    complete: bool
    beyond_certified_bp: int

    coordinates = None          # a class attribute (see Claim): CoordWalk carries them


@dataclass(slots=True)
class CoordWalk(Walk):
    """A walk of a retrieval with record coordinates: |coordinates| maps each label
    alive at the walk's leaf (its ref) to the coords.RunCoordinates of the run that reaches
    the leaf -- where that run's own bases [from_bp, to_bp) lie in its record (a run entered
    by the seed at 0 covers the whole walk)."""
    coordinates: Any = None


@dataclass(slots=True)
class LabelWalk:
    run: Optional[int]
    label: Label
    arm: str
    from_bp: int
    to_bp: int
    evidence_from: int
    sequence: Optional[str]
    end: str
    leaves_below: List[int]
    merged_into: Optional[int]

    coordinates = None          # a class attribute (see Claim): CoordLabelWalk carries them


@dataclass(slots=True)
class CoordLabelWalk(LabelWalk):
    """A label walk of a retrieval with record coordinates: |coordinates| is its run's
    coords.RunCoordinates (its own bases [from_bp, to_bp))."""
    coordinates: Any = None


@dataclass(slots=True)
class SupportRun:
    from_bp: int
    to_bp: int
    labels: List[Label]
    exact: bool
    total: Optional[int] = None


@dataclass(slots=True)
class Change:
    at_bp: int
    added: List[Label]
    removed: List[Label]
    reasons: List[dict]


@dataclass(slots=True)
class Continuation:
    sequence: str
    labels: List[Label]
    loss_used: float                     # the SMALLEST terminal loss of |labels| (C)
    branches_used: int
    seed_coord: Tuple[int, int]
    arm: str = 'right'
    leaf: int = 0
    # why the labels' names cannot be resubmitted (None: each is the name of exactly one
    # label of the retrieval and holds no U+FFFD)
    unverifiable: Optional[str] = None
    # constrain: each label's own terminal loss at the leaf (T), parallel to |labels|;
    # a continuation request can carry one loss budget only (see next_request)
    losses: Tuple[float, ...] = ()
    # what a continuation request cannot carry exactly for these labels (None: nothing):
    # one loss budget for labels that ended the walk at different losses
    note: Optional[str] = None

    def as_seed(self):
        """A /traverse seed: names are what the request resolves (an explicit list).
        Raises UnverifiableLabelName when a name may resolve to another label."""
        seed = {'sequence': self.sequence}
        if self.labels:
            if self.unverifiable:
                raise UnverifiableLabelName(self.unverifiable)
            seed['labels'] = [l.name for l in self.labels]
        return seed


class NextRequest(dict):
    """A /traverse request built by next_request(): to the server (and to json.dumps) a
    plain dict of the request's fields. What the request format cannot express rides
    beside it, never inside (the server rejects unknown fields): |notes| states every
    way the continuation may differ from one uninterrupted walk (a loss budget or a branch
    allowance that is conservative for some labels, or restarts when the caller asked, a
    switch target left out, a label alive at a leaf that is not seeded), |loss_budget| how
    the budget was derived ({original, effective, largest_terminal_loss, labels: [{name,
    ref, loss, remaining}]}; None in annotate mode and when there is no loss budget to
    spend), |branch_budget| how branching.max_label_branches was derived ({original,
    effective, largest_terminal_branches, reset}: the original minus the largest terminal
    branch count of the continued labels, "unlimited" kept, the original when reset; None
    in annotate mode and when the original allows no branch), and |left_out| every label
    of the retrieval's permitted pool that the continuation does not carry as one
    uninterrupted walk would ([{name, ref, why, walk}], per continued walk): neither a
    seed label of that walk nor in labels.extra (why 'unreachable' or
    'unverifiable_name'), or alive at a continued leaf without covering the
    continuation's whole tail (why 'alive_not_seeded', with its loss: not a seed label,
    though it may be a switch target in labels.extra); the notes name the first few."""

    def __init__(self, request=(), notes=(), loss_budget=None, left_out=(),
                 branch_budget=None):
        super().__init__(request)
        self.notes = list(notes)
        self.loss_budget = loss_budget
        self.left_out = list(left_out)
        self.branch_budget = branch_budget


@dataclass(slots=True)
class Comparison:
    comparable: Any                      # True | 'qualified' | 'unverifiable' | 'unknown' | False
    reason: str
    depth_used: Optional[int]
    equal: Optional[bool]                # None when equality cannot be claimed
    only_in_a: List[Any]
    only_in_b: List[Any]
    notes: List[str]
    differ: List[Any] = field(default_factory=list)
    mode: str = 'claims'
    support: Tuple[str, str] = ('kmer', 'kmer')
    scopes: Tuple[str, str] = ('per_path', 'per_path')
    strategies: Tuple[Any, Any] = (None, None)
    # local limits: the stop of a comparison its local budget interrupted (comparable 'unknown')
    local_stop: Optional[dict] = None

    def as_dict(self):
        out = {
            'comparable': self.comparable, 'reason': self.reason,
            'depth_used': self.depth_used, 'equal': self.equal, 'mode': self.mode,
            'only_in_a': list(self.only_in_a), 'only_in_b': list(self.only_in_b),
            'differ': list(self.differ), 'notes': list(self.notes),
            'support': list(self.support), 'scopes': list(self.scopes),
        }
        if self.local_stop is not None:
            out['local_stop'] = dict(self.local_stop)
        return out
