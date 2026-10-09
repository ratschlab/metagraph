"""Local limits: work and allocation budgets for the library's local operations (DESIGN §5.4).

Off by default: an operation called without a budget -- no `budget=` argument and no
ambient budget (local_budget()) -- charges nothing, its output does not depend on this
module, and nothing below is consulted beyond one `is None` test per charge point.

With a budget, a local call either completes or stops and says so (DESIGN §5.4):

  * WORK is counted in local work units (lwu): a deterministic, weighted count of the
    model elements an operation's algorithm visits and of the rows and text it produces
    (work model 2, WORK_MODEL; 1 lwu is about 0.1 us of CPython 3.11 on the reference
    machine). It is not CPU time. Units are charged at the COLD PRICE: a derivation an
    operation uses (paths, splits, merge maps, evidence, label summaries, ...) is charged at
    its structural price on every call, whether a cache holds it or not, so the same call
    on the same graphlet with the same budget charges the same units and stops at the same
    row in any process and whatever earlier calls cached. Warm calls pay more than they
    spend.
  * MEMORY is a modelled account (`memory_bound: "model"`): the bytes of the result
    objects, rows, spelled sequences, output text and derived caches the call builds or
    uses, each charged at a size modelled on CPython (ASCII text at one byte per character
    plus its header, lists, sets, dicts and the slotted records at their measured sizes,
    a parse at a fit of its traced peak), released inside a call only where a loop drops a
    temporary before its next step: the account bounds the call's peak. Every call also
    carries a base (CALL_BASE bytes, W_CALL lwu) for its frames, temporaries and setup,
    charged with its first charge. It is a model, not a measurement: Python can
    allocate more (interpreter internals, fragmentation, other CPython versions, text
    beyond ASCII), so the bound is soft in-process; a hard bound needs process limits
    (the service's). A budget's memory limit applies to each call on its own; its work
    accumulates across the calls that share it (an analysis allowance).
  * Every charge is admitted before the work or allocation it pays for, with one stated
    exception: what has a size known only once it is built -- the lines of dump(), up
    to 64 at a time, a label set decoded during a parse, and the walks a one-shot iterator
    passes to to_fasta(leaves=) (a generator: copied into a list as it is read, per arm)
    -- is charged as soon as it is built, so the account may run that far ahead of the
    check. The account itself
    never exceeds the limit. Cheap loops are admitted in blocks of rows (64) at their
    known prices; a block that does not fit is paid row by row, so a list stops at
    exactly the row it would stop at with a charge per row.
  * A deadline (`deadline_s`, elapsed seconds from the budget's creation) and cancel()
    are checked at charge points: a deadline at most every DEADLINE_POLL lwu, a
    cancellation at every charge. Both are non-deterministic and say so; a derivation
    charged in one piece runs to its end first (the stop's largest_charge says how much).

A stop raises LocalBudgetExceeded (an Exception, neither a ValueError -- a bad argument --
nor a RuntimeError -- the base of IncompleteRecording), carrying the LocalStop that states
it and, for the lists whose order allows it, a Partial of whole rows with a resume token.
compare() never raises for a budget: it returns comparable 'unknown'. No `except` in the
library catches a stop, and every cache entry is stored only when complete, so a stopped
call leaves only complete caches and a later unbudgeted call answers byte for byte as if
the stop had not happened.

Not charged (stated in the docs): the store's bookkeeping (refresh, cache_signature,
cache_bytes, memory_bytes, its spool lock, and the digest check of a body it parses on
demand: about 0.8 lwu per 256 bytes, against the parse's 32 per 256), the HTTP client's
decoding of a response, cursor MACs,
check_rules(), the tables of frames() (the library calls they make are charged), the MCP
framework's serialisation of a result, and interpreter work outside the charge points.
"""

import contextvars
import math
import sys
import time
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Any, List, Optional

__all__ = ['LocalLimits', 'LocalBudget', 'LocalStop', 'LocalBudgetExceeded', 'Partial',
           'local_budget', 'unbudgeted', 'current', 'WORK_MODEL', 'DEADLINE_POLL',
           'MEMORY_BOUND']

# Work model 2 (WORK_MODEL). Where an operation's charges are not obvious from its code:
#   * next_request()/next_requests() with a change_cost table: the labels owed the default
#     switch cost are walked in name order, as walker.cpp does, so the pops and their
#     charge are the same in every process (a set's hash-seeded order would make them
#     vary); and the table is charged: each copy of its entries (the merged strategy, the
#     request's, next_requests()' per walk: 3 W_ELEM and their bytes per entry, the copy's
#     memo while it is built), each check of them (W_ELEM per entry), and per switch search
#     the check and the build (W_ELEM each per entry), the count of its pairs (W_STEP per
#     entry), its edges (W_ELEM each), and the bytes of its table, edges and heap entries
#     and of the reached labels' losses;
#   * compare(mode='prefix_subset'): b's supported prefixes are sorted once per arm and
#     label and bisected per walk; a's claims are filed in the tree of their chains; an
#     omission is charged W_ROW, each segment its walk is followed through (W_ELEM + 1 per
#     512 bases), each anchor or bisection of its label (W_ELEM) and each recorded refusal
#     tested there (W_ELEM + 1 per 32 label ids), after a's refusals are indexed once per
#     comparison (W_ELEM each); the divergence of a's walks is charged only where it is
#     computed; a refusal or claim whose position lies before its own segment (a body the
#     server never writes) is filed under its ancestor by one walk of the chains' tree
#     (W_STEP per segment, a bisection per refusal or claim, and the sort of the lists they
#     join), not by a climb per refusal (quadratic on a deep comb);
#   * GraphletStore.standalone_text() and save_body(): a call of their own, so their base
#     (W_CALL, CALL_BASE) is charged once like every other call's, and the body's digest
#     check is one more lwu per 256 bytes;
#   * to_json(): the seed block (W_NODE per dropped label, W_PAIR per [from, to] run of
#     one), the seed-level and arm-level limitations and label_dict (W_NODE each), and
#     their bytes, before they are built;
#   * support_changes(): each label that changed is charged its reason's whole scan of the
#     walk's chain (segments, runs, events) and, once per change, the split it may have
#     taken (W_STEP per chain segment, W_ELEM per sibling) and per removed label a
#     bisection of each sibling's entry set and the parent's end set (W_ELEM and a unit per
#     halving); its bytes are 729 per changed label;
#   * walks() and walks_at(with_claims=True) in annotate mode: the index of each leaf's
#     route ends is built once (W_ELEM per route end, the sorts, its bytes), not a pass over
#     all of the arm's route ends per chunk of 32 walks;
#   * next_request()'s left-out note: its switch searches run over the names the
#     change_cost table names plus one stand-in (the same answer as a search of the whole
#     pool, with smaller searches, charged as such), and the note's reach tests (W_ELEM
#     each) and text are charged;
#   * compare(mode='walks' / 'prefix_subset') and compare_cost(): each segment of the cuts'
#     chains is charged its own part of cut_cost once per comparison, and each cut its
#     anchor's part, its prefix's join and its labels (and the pass's running state, held
#     for the pass), not its whole chain's cut_cost (quadratic on a comb); the scan of the
#     restricted leaves and the marking of their chains is W_ELEM + W_STEP per segment; a
#     cut's bytes count its label set and its {label: j} at their bound (in constrain mode
#     the root's entry and the first piece's labels);
#   * to_fasta(leaves=<one-shot iterator>): the copy of the walks it gives.
# The work model changes no unbudgeted answer.
WORK_MODEL = 2
MEMORY_BOUND = 'model'
# a deadline is read from the clock at most this many lwu apart (about 6.5 ms nominal)
DEADLINE_POLL = 65536

# ------------------------------------------------------------------ the weights
# lwu per element; 1 lwu ~ 0.1 us of CPython 3.11 on the reference machine (Apple M5 Max).
# Changing a weight moves stop points: it needs a new WORK_MODEL.
W_LINE = 48            # an MGT record line parsed (+ 1 lwu per 8 bytes of it)
W_DUMP_LINE = 24       # an MGT record line written
W_CALL = 200           # what any call costs beside its elements (its setup in Python)
W_NODE = 12            # a JSON node built (segment, path, run, event, presence run)
W_RECORD = 12          # a FASTA record (header and orientation)
W_GFA_LINE = 16        # a GFA segment, link or path line
W_ROW = 12             # a result row (claim, walk, split, label walk, change)
W_STEP = 1             # a chain step of a root -> leaf or route walk
W_ELEM = 2             # a model element visited (segment, run, presence run, event)
W_SORT = 2             # a comparison of a sort, per element and level (n log2 n)
ID_SHIFT_PY = 2        # label ids iterated in Python: 1 lwu per 4
ID_SHIFT_C = 5         # label ids in a C set or array operation: 1 lwu per 32
ID_SHIFT_PARSE = 1     # label ids decoded during a parse: 1 lwu per 2
BASE_SHIFT = 9         # bases spelled or copied: 1 lwu per 512
TEXT_SHIFT = 4         # text written (JSON, sizes measured): 1 lwu per 16 bytes
# an occurrence of record coordinates validated, clipped or written (feature level 6). No
# input of an earlier level has one, so adding the weight moved no stop point of work
# model 1
W_COORD = 2
# a [from, to] pair written into a JSON list (to_json()'s runs of a dropped label, work
# model 2)
W_PAIR = 2

# ------------------------------------------------------------------ the size model
# CPython sizes, the larger of 3.10-3.14 where they differ (an ASCII str's header is 49
# bytes up to 3.11 and 41 from 3.12)
STR = 49               # + 1 per ASCII character
LIST = 56              # + 9 per item (8 for the slot, ~1/8 over-allocation)
LIST_ITEM = 9
SET = 216              # + 54 per item (at most 3/5 full, 16 bytes per slot)
SET_ITEM = 54
DICT = 232             # + 16 per key, plus the values
DICT_KEY = 16
ARRAY = 64             # array('I'): + 4 per id
INT = 28               # an int above 256 (smaller ones are shared)
FLOAT = 24
TUPLE = 40             # + 8 per item
JSON_ID = 41           # a label id copied into a JSON list (its int and its slot, measured)
INT_ID = 37            # a label id in a list of ints (an int above 256 and its slot)
# what every call holds beside what it builds: its frames, its temporaries, the allocator's
# rounding of small blocks (measured: up to ~12 KB on the smallest graphlets); charged once
# per outermost call
CALL_BASE = 16384
TEXT_LINE = 72         # a line of text in a list: its header, slot and block rounding


def str_bytes(n):
    return STR + n


def list_bytes(n):
    return LIST + LIST_ITEM * n


def set_bytes(n):
    return SET + SET_ITEM * n


def dict_bytes(k):
    return DICT + DICT_KEY * k


def sort_work(n):
    """The comparisons of sorting n elements (timsort's n log2 n, at W_SORT each)."""
    return W_SORT * n * max(1, int(math.log2(n)) if n > 1 else 1)


_RECORD_SIZES = {}


def record_bytes(cls):
    """sys.getsizeof of an instance of a slotted record type (Claim, Walk, ...), taken
    once: the record's own slots; what its fields refer to is charged by the caller."""
    got = _RECORD_SIZES.get(cls)
    if got is None:
        slots = getattr(cls, '__slots__', ())
        try:
            inst = object.__new__(cls)
            got = sys.getsizeof(inst)
        except TypeError:
            got = 56 + 8 * len(slots)
        _RECORD_SIZES[cls] = got = max(got, 56 + 8 * len(slots))
    return got


# ------------------------------------------------------------------ limits

def _check_limit(name, v, integer):
    if v is None:
        return v
    if isinstance(v, bool) or not isinstance(v, (int, float)) or (integer and not isinstance(v, int)):
        raise ValueError('%s is %s >= 0 or None (unlimited), not %r'
                         % (name, 'an integer' if integer else 'a number', v))
    if not v >= 0 or (isinstance(v, float) and math.isinf(v)):
        raise ValueError('%s is %s >= 0 or None (unlimited), not %r'
                         % (name, 'an integer' if integer else 'a finite number', v))
    return v


@dataclass(frozen=True)
class LocalLimits:
    """The limits of a local call. None: unlimited. memory_mb bounds the modelled account
    of ONE call (MiB, 2**20 bytes); work_units the units charged (accumulating across the
    calls that share a LocalBudget); deadline_s the elapsed seconds from the budget's
    creation (non-deterministic)."""
    work_units: Optional[int] = None
    memory_mb: Optional[float] = None
    deadline_s: Optional[float] = None

    def __post_init__(self):
        _check_limit('work_units', self.work_units, True)
        _check_limit('memory_mb', self.memory_mb, False)
        _check_limit('deadline_s', self.deadline_s, False)

    @property
    def memory_bytes(self):
        return None if self.memory_mb is None else int(self.memory_mb * (1 << 20))

    def clamped(self, ceiling):
        """-> (these limits held to |ceiling|, the names of the fields that were lowered).
        A field the ceiling bounds and these leave unlimited takes the ceiling's value
        (and is named): a request never escapes a ceiling by asking for 'unlimited'."""
        if ceiling is None:
            return self, []
        out, cut = {}, []
        for name in ('work_units', 'memory_mb', 'deadline_s'):
            mine, top = getattr(self, name), getattr(ceiling, name)
            if top is not None and (mine is None or mine > top):
                out[name] = top
                cut.append(name)
            else:
                out[name] = mine
        return LocalLimits(**out), cut

    def as_dict(self):
        return {'work_units': self.work_units, 'memory_mb': self.memory_mb,
                'deadline_s': self.deadline_s}


# ------------------------------------------------------------------ stops

_RESOURCE_UNIT = {'work': 'lwu', 'memory': 'bytes', 'deadline': 's', 'cancelled': None}


@dataclass
class LocalStop:
    """Why a local call stopped (modelled on the backend's resource_stop): the resource,
    the operation and its phase, the limit, what was used and -- where it is known -- a
    lower bound of what the call needs (all in |unit|: lwu, bytes or seconds), the largest
    single charge, how far the call got (|done|: {'rows': n, 'of': N} and the like),
    whether the stop is reproducible, the work model, the levers and a message."""
    resource: str                       # work | memory | deadline | cancelled
    op: str
    phase: str
    limit: Any
    used: Any
    needed_at_least: Any
    largest_charge: int
    done: dict
    deterministic: bool
    actions: List[str]
    message: str
    unit: Optional[str] = None
    work_model: int = WORK_MODEL
    scope: str = 'local'

    def as_dict(self):
        return {'scope': self.scope, 'resource': self.resource, 'op': self.op,
                'phase': self.phase, 'unit': self.unit, 'limit': self.limit,
                'used': self.used, 'needed_at_least': self.needed_at_least,
                'largest_charge': self.largest_charge, 'done': dict(self.done),
                'deterministic': self.deterministic, 'work_model': self.work_model,
                'actions': list(self.actions), 'message': self.message}

    def restate(self, phase=None, actions=(), **done):
        """How far the operation got, stated where the stop is caught (a list knows its
        rows; the charge point only its phase); the message follows."""
        if phase is not None:
            self.phase = phase
        self.done.update(done)
        for a in actions:
            if a not in self.actions:
                # before process_locally, the lever of last resort
                self.actions.insert(max(0, len(self.actions) - 1), a)
        self.message = _message(self.op, self.resource, self.phase, self.limit, self.used,
                                self.needed_at_least, self.done, self.actions)
        return self


@dataclass
class Partial:
    """What a stopped list had produced: whole rows, a prefix in the operation's
    documented order, and a token that the same operation with the same arguments takes
    as resume= to go on after the last row (None where the order allows no resumption).
    |total|: the length of the whole list when the operation knew it before emitting."""
    rows: list
    resume: Any = None
    total: Optional[int] = None


class LocalBudgetExceeded(Exception):
    """A budgeted local call stopped: .stop (LocalStop) states it, .partial (Partial or
    None) holds the whole rows of a list stopped where its order allows it."""

    def __init__(self, stop, partial=None):
        super().__init__(stop.message)
        self.stop = stop
        self.partial = partial

    def restate(self, partial=None, phase=None, actions=(), **done):
        """Complete the statement where the stop is caught: how far the call got, and
        the whole rows it had produced (a list); -> self, to be raised again."""
        self.stop.restate(phase, actions, **done)
        self.args = (self.stop.message,)
        if partial is not None:
            self.partial = partial
        return self


_ACTIONS_ALWAYS = ('raise_local_budget',)
_ACTIONS_TAIL = ('process_locally',)


def _fmt(n):
    if isinstance(n, bool) or not isinstance(n, (int, float)):
        return str(n)
    if isinstance(n, float):
        return '%.3g' % n
    return '{:,}'.format(n)


# ------------------------------------------------------------------ the budget

class _Scope:
    __slots__ = ('b',)

    def __init__(self, b):
        self.b = b

    def __enter__(self):
        return self.b

    def __exit__(self, *exc):
        self.b._depth -= 1
        return False


class LocalBudget:
    """A stateful budget: one call, or several sharing one allowance (work accumulates;
    the memory limit holds for each call). Not thread-safe, except cancel(), which another
    thread may call: the running call stops at its next charge."""

    __slots__ = ('limits', '_work_limit', '_mem_limit', '_deadline_at', 'used_work',
                 'account_bytes', 'peak_bytes', 'largest_charge', 'stops', '_cancelled',
                 '_next_poll', '_depth', 'op', 'phase', 'done', 'levers', '_call_work0',
                 '_call_peak', '_t0', '_used', 'sub_op', 'row_extra', '_wcap', '_mcap',
                 '_base', '__weakref__')

    work_model = WORK_MODEL

    def __init__(self, limits=None, *, work_units=None, memory_mb=None, deadline_s=None):
        if limits is None:
            limits = LocalLimits(work_units, memory_mb, deadline_s)
        elif not isinstance(limits, LocalLimits):
            raise TypeError('limits is a LocalLimits, not %r' % (limits,))
        elif (work_units, memory_mb, deadline_s) != (None, None, None):
            raise TypeError('give LocalLimits or the keyword limits, not both')
        self.limits = limits
        self._work_limit = limits.work_units
        self._mem_limit = limits.memory_bytes
        # the limits as numbers (inf: none), so a charge compares without a None test
        self._wcap = math.inf if self._work_limit is None else self._work_limit
        self._mcap = math.inf if self._mem_limit is None else self._mem_limit
        self._t0 = time.monotonic()
        self._deadline_at = None if limits.deadline_s is None else self._t0 + limits.deadline_s
        self._next_poll = 0 if self._deadline_at is not None else math.inf
        self.used_work = 0
        self.account_bytes = 0
        self.peak_bytes = 0
        self.largest_charge = 0
        self.stops = []
        self._cancelled = False
        self._depth = 0
        self.op = None
        self.phase = None
        self.done = {}
        self.levers = ()
        self._call_work0 = 0
        self._call_peak = 0
        self._used = set()
        # the library operation running inside a caller's own scope (an MCP tool), named
        # in a stop; and what that caller adds to each row a list makes (its own row of a
        # page, built from the library's): charged with the row, so that a list stops
        # between whole rows of the caller's too
        self.sub_op = None
        self.row_extra = (0, 0)
        self._base = None

    # ---------------------------------------------------------------- calls

    def scope(self, op, levers=()):
        """Enter an operation: the outermost one names the stop and starts a new memory
        account (nothing a previous call held is held now); a nested one (walks() calling
        claims(), an export spelling walks) only counts its depth."""
        if self._depth == 0:
            self.op = op
            self.sub_op = None
            self.phase = 'setup'
            self.done = {}
            self.levers = levers
            self.account_bytes = 0
            self._call_work0 = self.used_work
            self._call_peak = 0
            self._used = set()
            self.row_extra = (0, 0)
            # what the call costs beside its elements (its setup, its frames and
            # temporaries): charged with its first charge, which the operation makes
            # where it states its stops
            self._base = (W_CALL, CALL_BASE)
            self._depth += 1
            return _Scope(self)
        elif self._depth == 1:
            self.sub_op = op
            self.levers = tuple(dict.fromkeys(tuple(self.levers) + tuple(levers)))
        self._depth += 1
        return _Scope(self)

    def call_usage(self):
        """The usage of the current (or last) outermost call: {work_units, memory_bytes
        (the peak of its account)}."""
        return {'work_units': self.used_work - self._call_work0,
                'memory_bytes': self._call_peak}

    # ---------------------------------------------------------------- charges

    def charge(self, work, nbytes=0):
        """Admit |work| lwu and |nbytes| of the account before the step they pay for, or
        stop (raise LocalBudgetExceeded) without charging them."""
        if self._base is not None:
            work += self._base[0]
            nbytes += self._base[1]
            self._base = None
        u = self.used_work + work
        if u > self._wcap or self._cancelled:
            self._stop('cancelled' if self._cancelled else 'work', work, nbytes)
        if nbytes:
            acc = self.account_bytes + nbytes
            if acc > self._mcap:
                self._stop('memory', work, nbytes)
            self.account_bytes = acc
            if acc > self._call_peak:
                self._call_peak = acc
                if acc > self.peak_bytes:
                    self.peak_bytes = acc
        self.used_work = u
        if work > self.largest_charge:
            self.largest_charge = work
        if u >= self._next_poll:
            self._poll()

    def admit(self, work, nbytes=0):
        """charge() of a block of steps admitted together, or False -- nothing charged and
        no stop -- when the block does not fit: the caller then pays its steps one by one,
        so a stop falls on exactly the step it would fall on without blocks (and charges
        exactly what they would). A deadline or a cancellation stops as in charge()."""
        if self._cancelled:
            self._stop('cancelled', work, nbytes)
        base = self._base or (0, 0)
        u = self.used_work + work + base[0]
        acc = self.account_bytes + nbytes + base[1]
        if u > self._wcap or acc > self._mcap:
            return False
        self.charge(work, nbytes)
        return True

    def uses(self, key, work, nbytes=0):
        """Charge a derivation the call uses at its (cold) price, once per outermost call
        however often the call's parts use it: |key| names it (an arm's id and the
        derivation's name)."""
        if key not in self._used:
            self.charge(work, nbytes)
            self._used.add(key)

    def first(self, key):
        """True the first time |key| is seen in the current outermost call (a part of
        the call to charge once), and marks it."""
        if key in self._used:
            return False
        self._used.add(key)
        return True

    def require(self, work=0, nbytes=0):
        """Stop now when |work| more lwu, or an account of |nbytes|, cannot fit (an
        admission check that charges nothing: a certain lower bound refused early)."""
        base = self._base or (0, 0)
        if self.used_work + work + base[0] > self._wcap:
            self._stop('work', work + base[0], 0)
        if self.account_bytes + nbytes + base[1] > self._mcap:
            self._stop('memory', 0, nbytes + base[1])

    def refund(self, work):
        """Return the unused remainder of a block admitted in advance."""
        self.used_work -= work

    def release(self, nbytes):
        """Take back from the account what a loop dropped before its next step (a
        per-segment temporary): the account bounds the call's peak, not the sum of
        everything it ever allocated. Called only once the object is gone."""
        self.account_bytes = max(0, self.account_bytes - nbytes)

    def cancel(self):
        """Stop the running (or the next) budgeted call at its next charge."""
        self._cancelled = True

    @property
    def cancelled(self):
        return self._cancelled

    def _poll(self):
        if self._deadline_at is not None:
            if time.monotonic() > self._deadline_at:
                self._stop('deadline', 0, 0)
            self._next_poll = self.used_work + DEADLINE_POLL

    # ---------------------------------------------------------------- reports

    def remaining(self):
        return {'work_units': None if self._work_limit is None
                else max(0, self._work_limit - self.used_work),
                'memory_bytes': self._mem_limit}

    def usage(self):
        """{work_units (charged, all calls), memory_bytes (the largest account of a
        call), largest_charge, work_model}."""
        return {'work_units': self.used_work, 'memory_bytes': self.peak_bytes,
                'largest_charge': self.largest_charge, 'work_model': WORK_MODEL}

    def _stop(self, resource, work, nbytes):
        op = self.sub_op or self.op or 'local'
        if resource == 'work':
            limit, used, need = self._work_limit, self.used_work, self.used_work + work
        elif resource == 'memory':
            limit, used, need = self._mem_limit, self.account_bytes, self.account_bytes + nbytes
        elif resource == 'deadline':
            limit = self.limits.deadline_s
            used = round(time.monotonic() - self._t0, 6)
            need = None
        else:
            limit = used = need = None
        deterministic = resource in ('work', 'memory')
        if deterministic:
            actions = list(_ACTIONS_ALWAYS) + [a for a in self.levers
                                              if a not in _ACTIONS_ALWAYS] + list(_ACTIONS_TAIL)
        else:
            actions = ['retry_later'] + (['raise_local_budget'] if resource == 'deadline'
                                         else [])
        stop = LocalStop(resource=resource, op=op, phase=self.phase or 'setup', limit=limit,
                         used=used, needed_at_least=need,
                         largest_charge=max(self.largest_charge, work),
                         done=dict(self.done), deterministic=deterministic,
                         actions=actions, message=_message(op, resource, self.phase, limit,
                                                           used, need, self.done, actions),
                         unit=_RESOURCE_UNIT[resource])
        self.stops.append(stop)
        raise LocalBudgetExceeded(stop)


def _message(op, resource, phase, limit, used, need, done, actions):
    what = {'work': 'local work budget', 'memory': 'local memory budget (a modelled account)',
            'deadline': 'local deadline', 'cancelled': 'cancellation'}[resource]
    msg = '%s stopped on the %s' % (op, what)
    unit = {'work': ' lwu', 'memory': ' bytes', 'deadline': ' s'}.get(resource, '')
    if limit is not None:
        msg += ': %s%s used of %s' % (_fmt(used), unit, _fmt(limit))
        if need is not None:
            msg += ', at least %s needed' % _fmt(need)
    if done:
        of = done.get('of')
        parts = ['%s=%s' % (k, _fmt(v)) for k, v in done.items() if k != 'of']
        if parts:
            msg += '; done: %s%s' % (', '.join(parts), '' if of is None else ' of %s' % _fmt(of))
    msg += ' (phase %s)' % (phase or 'setup')
    if resource in ('work', 'memory'):
        msg += '; nothing incomplete was returned. Levers: %s' % ', '.join(actions)
    return msg


# ------------------------------------------------------------------ ambient budget

_AMBIENT = contextvars.ContextVar('metagraph_traverse_local_budget', default=None)


def current():
    """The ambient budget (local_budget()), or None."""
    return _AMBIENT.get()


def resolve(budget):
    """The budget a call runs under: an explicit one, else the ambient one, else None."""
    if budget is not None:
        if not isinstance(budget, LocalBudget):
            raise TypeError('budget is a LocalBudget, not %r' % (budget,))
        return budget
    return _AMBIENT.get()


def run_scoped(budget, op, levers, fn):
    """The prologue of a budgeted operation: fn(None) when there is no budget (|budget|
    None and no ambient one, resolve()), else fn(b) inside b.scope(op, levers)."""
    b = resolve(budget)
    if b is None:
        return fn(None)
    with b.scope(op, levers):
        return fn(b)


@contextmanager
def unbudgeted():
    """Run the enclosed code with no ambient budget: what a store's parse without parse
    limits runs under, inside a tool call that has one (the parse is not the call's
    operation)."""
    token = _AMBIENT.set(None)
    try:
        yield
    finally:
        _AMBIENT.reset(token)


@contextmanager
def local_budget(limits=None, **kw):
    """Run the enclosed code under one budget (a LocalBudget, LocalLimits or the keyword
    limits): every library call without an explicit budget= charges it. Context-local
    (contextvars): a new thread starts without it; an asyncio task inherits the budget of
    the context it was created in, so the tasks created inside one block share it."""
    b = limits if isinstance(limits, LocalBudget) else LocalBudget(limits, **kw)
    token = _AMBIENT.set(b)
    try:
        yield b
    finally:
        _AMBIENT.reset(token)
