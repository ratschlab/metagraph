"""Record coordinates (DESIGN §18 and §26).

A request with `strategy.output.coordinates: true` gets, per seed result in every detail
level, either a `coordinates` block -- where each run's own bases lie in the indexed
records -- or `coordinates: null` with `coordinates_reason`. MGT v1 is frozen, so the block
is not in the body: it rides in the per-seed JSON summary (Graphlet.seed_summary), which
the J line of a saved file, the store's entries and every view keep with the body.

This module parses and validates the block (eagerly: an inconsistent block is rejected,
GraphletFormatError 'coordinates: ...'), keeps the parsed index in the graphlet's cache
(never a slot of the model: Graphlet, Claim, Walk and LabelWalk keep their layout, so the
local-limits account and the golden digests of graphlets without coordinates do not
depend on them) and clips a run's occurrences to a claim's cut.

The intervals (the frozen wire contract of level 6):
  * 0-based, half-open, on the record's forward strand. With c the chain's k-mer coordinate
    at the run's last node and L = to_bp - from_bp: right arm [c+k-L, c+k), left arm
    [c, c+L); L = 0 is an empty interval at the seed boundary. A seed occurrence is
    [c_first, c_first + |seed|).
  * Header labels (kind record): positions within the record. Column labels (kind column):
    the column's k-mer index space -- record i's k-mer j is offset_i + j, offset_{i+1} =
    offset_i + len_i - k + 1 -- so the last k-1 bases of a record share their numbers with
    the next record's first positions: a column interval cannot be attributed to one record
    without the record lengths, which the block does not carry, and two records' intervals
    are not necessarily disjoint. Kind mixed: each label's own kind says how to read its
    positions.
  * Each list keeps the first max_occurrences occurrences by start; a cut list states its
    true count (total), and the seed-level limitation `coordinates` counts the cut lists.
  * chains_ended: chains that continued on no followed path of the run's lineage (chains
    split among the children of a split are not counted; clones inherit the prefix's).
  * lower_bound: the run was entered by a switch into a label whose own lineage was live
    there; that label's chains starting at the switch node are left out, so the run's
    occurrences are a lower bound.

Stdlib only.
"""

from dataclasses import dataclass
from typing import Any, Dict, Tuple

from . import budget as _B
from . import derive
from ._codec import MAX_U64, UNLIMITED, GraphletFormatError, tok as _tok
from .budget import DICT, DICT_KEY, INT, TUPLE, W_ELEM, record_bytes
# defined in the light client (importing the client loads no model), named here too
from .client import strip_coordinate_cap  # noqa: F401

__all__ = ['Coordinates', 'RunCoordinates', 'SeedOccurrences', 'KINDS', 'REASONS',
           'present', 'requested', 'of', 'reason', 'validate', 'attach', 'clip',
           'label_kind', 'uses', 'index_price', 'strip_coordinate_cap']

KINDS = ('record', 'column', 'mixed')
# the null reasons of the contract, checked in this order by the server; a reason the
# library does not know is kept as stated (a later level may add one), never refused
REASONS = ('index has no coordinates', 'support kmer', 'no traversal', 'partial derivation')
LIMITATION_KIND = 'coordinates'
CAP_KNOB = 'output.max_coordinate_occurrences'
# the cache key of the parsed index (Graphlet.cache): set only for a graphlet that carries a
# block, so that a graphlet without one keeps its caches -- and its memory_bytes(),
# cache_signature() and every local-limits charge -- unaffected by coordinates
CACHE_KEY = 'coordinates'


@dataclass(slots=True, frozen=True)
class SeedOccurrences:
    """The seed's occurrences in one seed label: |intervals| ((start, end) pairs, ascending
    by start, each |seed| bases long), at most the cap, and |total| the true count."""
    label: int
    intervals: Tuple[Tuple[int, int], ...]
    total: int

    @property
    def truncated(self):
        return self.total > len(self.intervals)


@dataclass(slots=True, frozen=True)
class RunCoordinates:
    """One run's occurrences: where its own bases [from_bp, to_bp) lie in its label's
    record (or column). |intervals| ascend by start, each to_bp - from_bp long; |total| is
    the true count (> len(intervals) when the list was cut at the cap); |chains_ended|
    counts copies that stopped inside the run; |lower_bound| marks a switch-entered run
    whose occurrences may be understated. A claim's coordinates are its run's, clipped to
    the claim's cut (clip()): to_bp is then the cut."""
    run: int
    label: int
    from_bp: int
    to_bp: int
    intervals: Tuple[Tuple[int, int], ...]
    total: int
    chains_ended: int = 0
    lower_bound: bool = False

    @property
    def truncated(self):
        return self.total > len(self.intervals)

    def as_dict(self, kind=None):
        """The JSON form a claim, a walk or a label walk carries: intervals as [start,
        end] lists (0-based, half-open, forward strand); total, and truncated when the list
        was cut; chains_ended and lower_bound only when they say something; |kind| (the
        label's: 'record' or 'column') when given."""
        d = {'intervals': [list(x) for x in self.intervals], 'total': self.total}
        if self.truncated:
            d['truncated'] = True
        if self.chains_ended:
            d['chains_ended'] = self.chains_ended
        if self.lower_bound:
            d['lower_bound'] = True
        if kind is not None:
            d['kind'] = kind
        return d


@dataclass(slots=True, frozen=True)
class Coordinates:
    """A seed result's coordinates block, parsed and validated against its body: |kind|
    (record | column | mixed), |k|, |max_occurrences| (an int or UNLIMITED), |complete|
    (no list cut and no lower-bound run), |seed| {seed label id: SeedOccurrences},
    |arms| {side: (RunCoordinates per run, indexed by run id)}, |runs_lower_bound| the
    number of runs marked lower_bound, |lists_cut| the number of lists cut at the cap."""
    kind: str
    k: int
    max_occurrences: Any
    complete: bool
    seed: Dict[int, SeedOccurrences]
    arms: Dict[str, Tuple[RunCoordinates, ...]]
    runs_lower_bound: int = 0
    lists_cut: int = 0

    def run(self, side, run_id):
        """The RunCoordinates of run |run_id| of arm |side| (None: the arm has none)."""
        runs = self.arms.get(side)
        return None if runs is None else runs[run_id]

    def summary(self):
        """{kind, complete, max_occurrences, lists_cut?, runs_lower_bound?}: what an
        evidence block or a summary states about the block."""
        d = {'kind': self.kind, 'complete': self.complete,
             'max_occurrences': 'unlimited' if self.max_occurrences is UNLIMITED
             else self.max_occurrences}
        if self.lists_cut:
            d['lists_cut'] = self.lists_cut
        if self.runs_lower_bound:
            d['runs_lower_bound'] = self.runs_lower_bound
        return d


# ------------------------------------------------------------------ presence

def requested(g):
    """Whether the seed result answers a request for coordinates: its summary has the
    `coordinates` field (a block or null)."""
    s = g.seed_summary
    return isinstance(s, dict) and 'coordinates' in s


def present(g):
    """Whether the seed result carries a coordinates block (not null): a test of the
    summary alone, which parses nothing (compare() and the evidence block use it)."""
    s = g.seed_summary
    return isinstance(s, dict) and isinstance(s.get('coordinates'), dict)


def reason(g):
    """Why the graphlet has no coordinates: 'no envelope' (a body without its envelope,
    which never carries them), 'not requested' (the request did not ask), or the server's
    coordinates_reason; None when it has a block."""
    if not g.has_envelope:
        return 'no envelope'
    s = g.seed_summary
    if 'coordinates' not in s:
        return 'not requested'
    if isinstance(s['coordinates'], dict):
        return None
    r = s.get('coordinates_reason')
    return r if isinstance(r, str) else None


def of(g, b=None):
    """The graphlet's parsed Coordinates, or None (no block). Built and validated from the
    summary the first time (GraphletFormatError when inconsistent) and kept in g.cache.
    |b| (local limits): the index is a derivation of the summary, charged at its cold
    price once per call (uses()) whether the cache holds it or not, and before it is built
    (a stop leaves nothing in the cache)."""
    if not present(g):
        return None
    if b is not None:
        # keyed, before the build: a cold build charged through validate(g, b) is unkeyed,
        # so the call's next use (the other arm's claims) would charge it a second time,
        # and a call would cost more after the store re-parsed the model than while it was
        # resident
        uses(b, g)
    got = g.cache.get(CACHE_KEY)
    if got is None:
        got = validate(g)
        g.cache[CACHE_KEY] = got
    return got


def label_kind(g, label):
    """How a label's positions are read: 'record' (a header label: within the record) or
    'column' (global column positions)."""
    return 'record' if g.labels[label].kind == 'header' else 'column'


# ------------------------------------------------------------------ local limits

# the parsed index's modelled size: a record per entry with its interval tuple and three
# ints, and per occurrence a (start, end) tuple of two ints (positions are above 256, so
# never the interpreter's shared small ints). The ints are in fact the summary's own (a
# pair holds the objects json.loads made, never new ones), so the 2 * INT is the margin
# that pays for what validating allocates besides the index: _intervals' list beside its
# tuple (a slot and its over-allocation per occurrence, while the tuple is built)
_OCC_BYTES = TUPLE + 2 * 8 + 2 * INT


def _entry_bytes():
    return max(record_bytes(RunCoordinates), record_bytes(SeedOccurrences)) + TUPLE + 3 * INT


def _block_counts(block):
    """(entries, occurrences) of a raw block, read without validating it (the price of
    validating it): a malformed block counts what it has and is refused by validate()."""
    entries = occ = 0
    seed = block.get('seed')
    lists = [seed] if isinstance(seed, list) else []
    arms = block.get('arms')
    if isinstance(arms, dict):
        lists.extend(v for v in arms.values() if isinstance(v, list))
    for lst in lists:
        entries += len(lst)
        for e in lst:
            if isinstance(e, dict) and isinstance(e.get('occurrences'), list):
                occ += len(e['occurrences'])
    return entries, occ


def index_price(entries, occ):
    """(lwu, bytes) of validating and holding an index of |entries| entries and |occ|
    occurrences: W_ELEM per entry and W_COORD per occurrence; the records, their interval
    tuples and the dicts that index them."""
    return (W_ELEM * entries + _B.W_COORD * occ,
            record_bytes(Coordinates) + 3 * DICT + DICT_KEY * entries
            + entries * _entry_bytes() + occ * (_OCC_BYTES + 8))


def uses(b, g):
    """Charge budget |b| the graphlet's coordinate index at its cold price, once per call
    (a derivation is charged on every use, whether the cache holds it or not)."""
    if b is None or not present(g):
        return
    got = g.cache.get(CACHE_KEY)
    if got is not None:
        entries = len(got.seed) + sum(len(v) for v in got.arms.values())
        occ = sum(len(x.intervals) for x in got.seed.values()) + sum(
            len(r.intervals) for v in got.arms.values() for r in v)
    else:
        entries, occ = _block_counts(g.seed_summary['coordinates'])
    b.uses((id(g), CACHE_KEY), *index_price(entries, occ))


# ------------------------------------------------------------------ validation

def _err(msg):
    return GraphletFormatError(0, 'coordinates: ' + msg)


def _uint(v, what, minimum=0):
    """|v| as a JSON integer >= |minimum| and below 2^64 (a position or a count of the
    server's): a bool, a float, a string or a huge number is refused, never coerced."""
    if isinstance(v, bool) or not isinstance(v, int) or v < minimum or v > MAX_U64:
        raise _err('%s is not an integer >= %d below 2^64: %s' % (what, minimum, _tok(v)))
    return v


def _intervals(raw, length, what):
    """The (start, end) tuples of an occurrence list: each end - start == |length|, the
    starts strictly ascending (two chains are two coordinates)."""
    if not isinstance(raw, list):
        raise _err('%s: occurrences is not a list' % what)
    out = []
    prev = -1
    for x in raw:
        if not isinstance(x, list) or len(x) != 2:
            raise _err('%s: an occurrence is not a [start, end] pair: %s' % (what, _tok(x)))
        s = _uint(x[0], what + ' start')
        e = _uint(x[1], what + ' end')
        if e - s != length:
            raise _err('%s: [%d, %d) is %d bases, not %d' % (what, s, e, e - s, length))
        if s <= prev:
            raise _err('%s: occurrences do not ascend by start (%d after %d)' % (what, s, prev))
        prev = s
        out.append((s, e))
    return tuple(out)


def _total(entry, n, cap, what):
    """The true count of a list of |n| occurrences: occurrences_total when the list was cut
    (then n is the cap and the total above it), else n (and n within the cap)."""
    if 'occurrences_total' in entry:
        total = _uint(entry['occurrences_total'], what + ' occurrences_total')
        if cap is UNLIMITED or n != cap or total <= n:
            raise _err('%s: occurrences_total %d stated for a list of %d that was not cut at '
                       'the cap %s' % (what, total, n, 'unlimited' if cap is UNLIMITED
                                       else cap))
        return total
    if cap is not UNLIMITED and n > cap:
        raise _err('%s: %d occurrences above the cap %d' % (what, n, cap))
    return n


def _expected_kind(g):
    kinds = {l.kind for l in g.labels}
    if kinds == {'header'}:
        return 'record'
    if kinds == {'column'}:
        return 'column'
    if kinds:
        return 'mixed'
    return None


def _echo(g):
    """The request's strategy.output as the envelope echoes it (None: no echo to check)."""
    strategy = (g.envelope or {}).get('strategy')
    if not isinstance(strategy, dict):
        return None
    out = strategy.get('output')
    return out if isinstance(out, dict) else {}


def _check_null(g, s):
    r = s.get('coordinates_reason')
    if not isinstance(r, str) or not r:
        raise _err('null without a coordinates_reason')
    if r == 'support kmer' and g.support == 'trace' and g.mode == 'constrain':
        raise _err('reason "support kmer" on a support: trace retrieval in constrain mode')
    if r == 'partial derivation' and not (g.support == 'trace' and derive.partial_derivation(g) is not None):
        raise _err('reason "partial derivation" without a walked trace result whose '
                   'permitted set was derived from part of the seed (a derivation '
                   'limitation)')


def validate(g, b=None):
    """The summary's coordinates checked against the body -> Coordinates, or None for
    the null form (its reason checked). Raises GraphletFormatError('coordinates: ...') when
    anything is inconsistent: the block's presence against the request's echo; kind
    against the L records and k against H; one seed entry per seed label, ascending;
    one entry per run of every requested arm, in R order, with the R record's run, label,
    from_bp and to_bp; intervals of the right length, ascending; occurrences_total exactly
    on the lists cut at the cap; chains_ended > 0 and lower_bound true where stated (and
    lower_bound only on a switch-entered run); complete, runs_lower_bound and the
    seed-level K record `coordinates` (lists_cut, observed, limit) agreeing with the lists;
    and every seed-entered run's occurrences continuing a seed occurrence of its label
    where the seed list was not cut. |b| (local limits): the index is charged before it
    is built (a stop leaves nothing in the cache)."""
    s = g.seed_summary
    if not isinstance(s, dict):
        return None
    echo = _echo(g)
    asked = 'coordinates' in s
    if echo is not None:
        if echo.get('coordinates') is True and not asked:
            raise _err('requested (strategy.output.coordinates true) but the result carries '
                       'neither a block nor a reason')
        if asked and echo.get('coordinates') is not True:
            raise _err('the result carries coordinates the request did not ask for (no '
                       'strategy.output.coordinates true in the echo)')
    if not asked:
        return None
    block = s['coordinates']
    if block is None:
        _check_null(g, s)
        return None
    if not isinstance(block, dict):
        raise _err('neither an object nor null: %s' % _tok(block))
    if 'coordinates_reason' in s:
        raise _err('a block and a coordinates_reason')
    if g.support != 'trace' or g.mode != 'constrain':
        raise _err('a block on a %s retrieval in %s mode: coordinates exist under support '
                   'trace only (§18.1)' % (g.support, g.mode))
    if b is not None:
        b.charge(*index_price(*_block_counts(block)))
    kind = block.get('kind')
    if kind not in KINDS:
        raise _err('kind %s' % _tok(kind))
    want = _expected_kind(g)
    if want is not None and kind != want:
        raise _err('kind %s, the L records say %s' % (kind, want))
    k = _uint(block.get('k'), 'k', 1)
    if k != g.k:
        raise _err('k %d, H says %d' % (k, g.k))
    raw_cap = block.get('max_occurrences')
    cap = UNLIMITED if raw_cap == 'unlimited' else _uint(raw_cap, 'max_occurrences', 1)
    if echo is not None and 'max_coordinate_occurrences' in echo \
            and echo['max_coordinate_occurrences'] != raw_cap:
        raise _err('max_occurrences %s, the request echoes %s'
                   % (_tok(raw_cap), _tok(echo['max_coordinate_occurrences'])))
    complete = block.get('complete')
    if not isinstance(complete, bool):
        raise _err('complete is not a boolean: %s' % _tok(complete))
    cut = []                        # the true totals of the lists cut at the cap

    # the seed: one entry per seed label, ascending
    raw_seed = block.get('seed')
    if not isinstance(raw_seed, list):
        raise _err('seed is not a list')
    nsl = g.seed.num_seed_labels
    if len(raw_seed) != nsl:
        raise _err('%d seed entries for %d seed labels' % (len(raw_seed), nsl))
    seed = {}
    for i, e in enumerate(raw_seed):
        if not isinstance(e, dict):
            raise _err('seed entry %d is not an object' % i)
        label = _uint(e.get('label'), 'seed entry %d label' % i)
        if label != i:
            raise _err('seed entry %d is label %d: one entry per seed label, ascending'
                       % (i, label))
        what = 'seed label %d' % label
        iv = _intervals(e.get('occurrences'), g.seed.length_bp, what)
        total = _total(e, len(iv), cap, what)
        if total > len(iv):
            cut.append(total)
        seed[label] = SeedOccurrences(label, iv, total)

    # the arms: exactly the requested ones, one entry per run in R order
    raw_arms = block.get('arms')
    if not isinstance(raw_arms, dict) or set(raw_arms) != set(g.arms):
        raise _err('arms %s, the body has %s' % (
            _tok(sorted(raw_arms) if isinstance(raw_arms, dict) else raw_arms),
            sorted(g.arms)))
    arms = {}
    lower = 0
    for side in [x for x in ('left', 'right') if x in g.arms]:
        runs = g.arms[side].runs
        lst = raw_arms[side]
        if not isinstance(lst, list) or len(lst) != len(runs):
            raise _err('the %s arm has %s entries for %d runs' % (
                side, len(lst) if isinstance(lst, list) else _tok(lst), len(runs)))
        out = []
        for i, (e, run) in enumerate(zip(lst, runs)):
            what = '%s run %d' % (side, i)
            if not isinstance(e, dict):
                raise _err('%s is not an object' % what)
            if _uint(e.get('run'), what + ' run') != i:
                raise _err('%s names run %s: one entry per run, in R order'
                           % (what, _tok(e.get('run'))))
            for field_, have in (('label', run.label), ('from_bp', run.from_bp),
                                 ('to_bp', run.to_bp)):
                v = _uint(e.get(field_), '%s %s' % (what, field_))
                if v != have:
                    raise _err('%s: %s %d, the R record says %d' % (what, field_, v, have))
            iv = _intervals(e.get('occurrences'), run.to_bp - run.from_bp, what)
            total = _total(e, len(iv), cap, what)
            if total > len(iv):
                cut.append(total)
            ended = 0
            if 'chains_ended' in e:
                ended = _uint(e['chains_ended'], what + ' chains_ended', 1)
            lb = False
            if 'lower_bound' in e:
                if e['lower_bound'] is not True:
                    raise _err('%s: lower_bound is stated only as true' % what)
                if run.from_label is None:
                    raise _err('%s: lower_bound on a run entered by the seed (only a switch '
                               'into a live label understates)' % what)
                lb = True
                lower += 1
            out.append(RunCoordinates(i, run.label, run.from_bp, run.to_bp, iv, total, ended,
                                      lb))
        arms[side] = tuple(out)

    # the lists against the seed: a seed-entered run continues a seed occurrence of its
    # label (chains only shrink along a lineage), checkable where the seed list is whole.
    # Both lists ascend by start and the seed start a run occurrence implies grows with it
    # (right: s - from_bp - n, left: e + from_bp), so one merge pass over the two tuples
    # checks them: no set of the seed starts per run, which the index price never paid for
    # (two sets of 5,000 starts alive at once would be 1.3 MB over a block the price puts
    # at 1.8 MB, and the traced peak would pass the account)
    n = g.seed.length_bp
    for side, runs in arms.items():
        r_of = g.arms[side].runs
        for rc in runs:
            if r_of[rc.run].from_label is not None:
                continue
            so = seed.get(rc.label)
            if so is None or so.truncated:
                continue
            seed_iv = so.intervals
            j, m = 0, len(seed_iv)
            for s_, e_ in rc.intervals:
                at = s_ - rc.from_bp - n if side == 'right' else e_ + rc.from_bp
                while j < m and seed_iv[j][0] < at:
                    j += 1
                if j == m or seed_iv[j][0] != at:
                    raise _err('%s run %d: [%d, %d) continues no seed occurrence of label %d'
                               % (side, rc.run, s_, e_, rc.label))

    # runs_lower_bound, complete and the K record against the lists
    if lower:
        if _uint(block.get('runs_lower_bound'), 'runs_lower_bound', 1) != lower:
            raise _err('runs_lower_bound %s, %d runs are marked'
                       % (_tok(block.get('runs_lower_bound')), lower))
    elif 'runs_lower_bound' in block:
        raise _err('runs_lower_bound stated with no run marked lower_bound')
    if complete != (not cut and not lower):
        raise _err('complete %s with %d list(s) cut and %d lower-bound run(s)'
                   % (complete, len(cut), lower))
    lims = [l for l in g.limitations if l.kind == LIMITATION_KIND]
    if any(l.arm is not None for l in lims):
        raise _err('a coordinates limitation of an arm (it is a seed-level limitation)')
    if not cut:
        if lims:
            raise _err('a coordinates limitation (K) with no list cut')
    else:
        if len(lims) != 1:
            raise _err('%d list(s) cut and %d coordinates limitations (K): exactly one'
                       % (len(cut), len(lims)))
        lim = lims[0]
        extra = dict(lim.extra)
        if lim.knob != CAP_KNOB or lim.limit != cap or lim.observed != max(cut) \
                or extra.get('lists_cut') != len(cut):
            raise _err('the coordinates limitation (K: knob %s, limit %s, observed %s, '
                       'lists_cut %s) does not state the %d list(s) cut at %s, the largest '
                       '%d' % (_tok(lim.knob), _tok(lim.limit), _tok(lim.observed),
                               _tok(extra.get('lists_cut')), len(cut), cap, max(cut)))
    return Coordinates(kind, k, cap, complete, seed, arms, lower, len(cut))


def attach(g, b=None):
    """Validate the summary's coordinates eagerly and keep the parsed index in
    g.cache: what from_response(), a parse of a saved file's J line and the store's first
    parse of an unparsed entry do. Nothing is cached for a graphlet without a block."""
    got = validate(g, b)
    if got is not None:
        g.cache[CACHE_KEY] = got
    return got


# ------------------------------------------------------------------ requests

# strip_coordinate_cap(): imported above from the client


# ------------------------------------------------------------------ clipping

def clip(rc, side, at_most_bp=None):
    """|rc| (a run's RunCoordinates) clipped to a claim cut at |at_most_bp|: the run's own
    bases [from_bp, D') with D' = min(to_bp, at_most_bp) -- right arm [s, s + (D' - from_bp)),
    left arm [e - (D' - from_bp), e) (the left arm walks backwards along the record). A
    claim's coordinates are its run's own route: a route_only claim's too (§18.3). The
    total, chains_ended and lower_bound stay the run's (a cut changes no chain)."""
    if at_most_bp is None or at_most_bp >= rc.to_bp:
        return rc
    d = max(rc.from_bp, at_most_bp) - rc.from_bp
    if side == 'right':
        iv = tuple((s, s + d) for s, _ in rc.intervals)
    else:
        iv = tuple((e - d, e) for _, e in rc.intervals)
    return RunCoordinates(rc.run, rc.label, rc.from_bp, rc.from_bp + d, iv, rc.total,
                          rc.chains_ended, rc.lower_bound)
