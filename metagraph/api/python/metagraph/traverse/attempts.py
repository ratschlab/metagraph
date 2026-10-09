"""The answers about an attempt and the release rule (feature level 5, SPEC §5 and §10.3).

What a ledger reads from the server about one attempt -- a cancel's or GET
/traverse/attempt's answer (AttemptAnswer), a 409 (classify_409(): AttemptConflict,
AttemptExpired, InstanceMismatch), a response or an error that carries the attempt's usage
-- and release_verdict(), the capabilities' `release_rule` applied to one of them.

A light module: the standard library and nothing else of this package. Importing it, or
metagraph.traverse.client (which re-exports every name here), or `from metagraph.traverse
import release_verdict` (the package's names are loaded lazily), loads no parser, model or
operation, so a process that only dispatches attempts and releases their capacity (the
search service's API side) pays for none of them.

The readers take the wire as the server writes it and are typed where a release depends on
a value: an instant is an integer in [0, 2^53 - 1] -- never a bool or a float, which
Python compares equal to integers (True == 1) --, an instance id a string, a flag a JSON
bool, and a key whose value is null reads as absent (the duplicate 409 of a level-5 server
writes `"tombstone": null`). A value of another type reads as None, which no release
condition accepts.
"""

import datetime
from dataclasses import dataclass, field
from typing import Any, List, Optional, Tuple

__all__ = ['Suppression', 'AttemptSent', 'AttemptAnswer', 'TraverseError',
           'ServerInitializing', 'AttemptAtBound', 'AttemptExpired', 'AttemptConflict',
           'InstanceMismatch', 'TraverseResponse', 'ReleaseVerdict', 'release_verdict',
           'classify_409', 'WIRE_INT_MAX']

# the wire's integers (not_after_ms, suppressed_until_ms, server_time_ms): [0, 2^53 - 1]
WIRE_INT_MAX = (1 << 53) - 1


def _wire_int(v):
    """|v| when it is an integer of the wire's range, else None: a bool is no instant
    (True == 1 in Python, so `sup.not_after_ms != sent.not_after_ms` took true for 1) and a
    float is not what the server writes."""
    if isinstance(v, bool) or not isinstance(v, int) or not 0 <= v <= WIRE_INT_MAX:
        return None
    return v


def _wire_str(v):
    return v if isinstance(v, str) else None


def _wire_bool(v):
    return v if isinstance(v, bool) else None


def _obj(body, key):
    """body[key] when |body| is an object and that value is one, else None: a non-object
    `attempt` (or usage) is read as absent, never indexed."""
    v = body.get(key) if isinstance(body, dict) else None
    return v if isinstance(v, dict) else None


# the keys only a suppression statement has (not_after_ms is also the attempt's own key
# in a state, so it alone does not make one)
_SUPPRESSION_KEYS = ('suppressed_until_ms', 'covers_admission', 'covers_admission_reason')


@dataclass(frozen=True)
class Suppression:
    """What a tombstone answer states (feature level 5): |suppressed_until_ms| the hold's
    last millisecond on the server's wall clock (inclusive, Unix epoch ms), |not_after_ms|
    the one it was judged against (None: none was), |covers_admission| whether no copy of
    that request can start on this server_instance -- True or False as stated, None when the
    answer does not state it as a JSON bool --, and when not, |covers_admission_reason|
    (`no_not_after_ms` | `beyond_tombstone_max`). An instant the answer leaves out or writes
    as anything but an integer is None."""
    suppressed_until_ms: Optional[int]
    not_after_ms: Optional[int]
    covers_admission: Optional[bool]
    covers_admission_reason: Optional[str] = None

    @classmethod
    def of(cls, body):
        """The suppression a tombstone answer (or a 409's `attempt`) states; None where it
        states none: no suppressed_until_ms, covers_admission or covers_admission_reason
        (no tombstone, or a server below feature level 5). Present when any of them is,
        each field read on its own (a missing suppressed_until_ms is None, not an error)."""
        if not isinstance(body, dict) or all(body.get(k) is None for k in _SUPPRESSION_KEYS):
            return None
        return cls(_wire_int(body.get('suppressed_until_ms')),
                   _wire_int(body.get('not_after_ms')),
                   _wire_bool(body.get('covers_admission')),
                   _wire_str(body.get('covers_admission_reason')))


@dataclass(frozen=True)
class AttemptSent:
    """The attempt fields of a request as it was SENT -- what the release rule is applied to:
    a field the client did not send (expect_server_instance to a server below feature level
    5) is None here, whatever was asked for."""
    attempt_id: str
    not_after_ms: Optional[int] = None
    expect_server_instance: Optional[str] = None

    @classmethod
    def of(cls, request):
        """From a request dict (build_request()'s); None without attempt_id."""
        if isinstance(request, cls) or request is None:
            return request
        if request.get('attempt_id') is None:
            return None
        return cls(request['attempt_id'], request.get('not_after_ms'),
                   request.get('expect_server_instance'))


class AttemptAnswer(dict):
    """The JSON answer of POST /traverse/cancel or GET /traverse/attempt (a dict),
    with its HTTP |status|, the |route| (`cancel` | `attempt`) and, for a cancel, the body
    |sent| (whether not_after_ms went out with it). The typed reading: attempt_id, state
    and server_instance at the top level, else in a cancel's `attempt` object."""

    def __init__(self, body, status, route, sent=None):
        super().__init__(body if isinstance(body, dict) else {})
        self.status = status
        self.route = route
        self.sent = sent

    def _top_or_block(self, key):
        v = _wire_str(self.get(key))
        if v is None:
            v = _wire_str((_obj(self, 'attempt') or {}).get(key))
        return v

    @property
    def attempt_id(self):
        """The id the answer is about (top level, else the cancel's `attempt` object); None
        where it names none -- such an answer releases nothing (release_verdict)."""
        return self._top_or_block('attempt_id')

    @property
    def server_instance(self):
        return self._top_or_block('server_instance')

    @property
    def state(self):
        """running | stopping | finished | unknown."""
        return self._top_or_block('state')

    @property
    def cancelled(self):
        return self.get('cancelled')

    @property
    def tombstone(self):
        """True: the id is tombstoned (a 404); False: a cancel that was NOT tombstoned (a
        429); None: the answer says neither (or not as a JSON bool)."""
        return _wire_bool(self.get('tombstone'))

    @property
    def reason(self):
        """A 429's reason: `tombstones_full` (retry the cancel later) or `no_suppression` (the
        server keeps no tombstones: a request with the id arriving later runs). For an attempt
        state, why it finished."""
        return self.get('reason')

    @property
    def attempt(self):
        """The attempt's state: a cancel's `attempt` object, GET's answer itself (200)."""
        if self.route == 'attempt':
            return dict(self) if self.status == 200 else None
        return _obj(self, 'attempt')

    @property
    def usage(self):
        """What the attempt consumed, once it finished: a cancel's in its `attempt` object
        (else at the top level), GET /traverse/attempt's at the top level; None where the
        answer carries none (an attempt still running, a tombstone)."""
        if self.route == 'cancel':
            got = _obj(_obj(self, 'attempt'), 'usage')
            if got is not None:
                return got
        return _obj(self, 'usage')

    @property
    def suppression(self):
        """The tombstone's Suppression (feature level 5), or None."""
        return Suppression.of(self)

    @property
    def finished(self):
        return self.state == 'finished'

    @property
    def retryable(self):
        """A 429 that a later cancel may turn into a tombstone (tombstones_full)."""
        return self.status == 429 and self.get('reason') == 'tombstones_full'


class TraverseError(RuntimeError):
    """A non-2xx answer: |status| the HTTP status, |message| the server's `error`, |body| the
    answer's JSON. |usage|: the usage block an error to a request with attempt_id carries
    (a 400 or 500 after the request was read, a 503 at the attempt's bound), else None --
    what the attempt consumed is there, not in a TraverseResponse."""

    def __init__(self, status, message, body=None):
        self.status = status
        self.message = message
        self.body = body
        self.usage = _obj(body, 'usage')
        # traverse(): the attempt fields as sent (AttemptSent), for release_verdict()
        self.sent = None
        super().__init__('HTTP %s: %s' % (status, message))


class ServerInitializing(TraverseError):
    """503 with Retry-After: the index is still loading (nothing was registered or run);
    retry after |retry_after| seconds."""

    def __init__(self, status, message, body=None, retry_after=None):
        super().__init__(status, message, body)
        self.retry_after = retry_after


class AttemptAtBound(TraverseError):
    """503 with `usage` (reason `deadline`): the attempt reached the duration bound the
    server enforces for it while its response was built or written, so nothing of it is
    delivered. It was registered and ran, and |usage| states what it consumed. A retry under
    the same attempt_id is refused (409) only while the server keeps the id: while the
    attempt is retained (retention_s, retention_count) and, for one sent with not_after_ms,
    held through not_after_ms + clock_skew_allowance_ms (at most tombstone_max_s); after that
    -- at once with retention_s 0 -- a retry runs again (SPEC §5, §10.3). Not a loading
    server: never retry it as one."""


class AttemptExpired(TraverseError):
    """409 with `state: "expired"`: the request's not_after_ms had passed on the server's
    clock when its handler started, so it was NOT started — nothing ran, nothing was
    registered (no usage; GET /traverse/attempt answers 404 for its id). The body states
    `not_after_ms` and the server's `server_time_ms` (|not_after_ms|, |server_time_ms|,
    |server_instance|, and the |attempt_id| it echoes). The server answers this only for an
    id it does not hold: one that is running, retained, held or tombstoned is refused as such
    first (AttemptConflict), so an expired answer also says that no copy of the attempt was
    registered there when it was judged -- a release ground (release_verdict(): `expired`)."""

    @property
    def not_after_ms(self):
        return _wire_int((self.body or {}).get('not_after_ms'))

    @property
    def server_time_ms(self):
        return _wire_int((self.body or {}).get('server_time_ms'))

    @property
    def server_instance(self):
        return _wire_str((self.body or {}).get('server_instance'))

    @property
    def attempt_id(self):
        return _wire_str((self.body or {}).get('attempt_id'))


class AttemptConflict(TraverseError):
    """409 with `attempt`: the attempt_id is running, retained, held after it finished, or
    tombstoned by a cancel on this server, so the request was refused and nothing ran for it.
    |attempt| is that id's state as GET /traverse/attempt answers it ({} when the body's
    `attempt` is not an object); |attempt_id| and |state| are read from it; |tombstoned|: a
    cancel named the id first (state unknown, `tombstone: true`), and then |suppression|
    states the hold (feature level 5), judged against this request's own not_after_ms -- the
    refusal extended it. |server_instance| and |tombstoned| fall back to the body's top level
    where the object lacks them (a null counts as absent). |usage|: the usage of the attempt
    the id names -- a finished one's, what it consumed; this refused copy consumed nothing.

    For release_verdict() the object is the id's state: `finished` is a finished state,
    `running` and `stopping` hold, a tombstone holds (a 409 frees nothing early)."""

    def __init__(self, status, message, body=None):
        super().__init__(status, message, body)
        got = _obj(_obj(body, 'attempt'), 'usage')
        if got is not None:
            self.usage = got

    @property
    def attempt(self):
        return _obj(self.body, 'attempt') or {}

    def _block_or_top(self, key):
        v = self.attempt.get(key)
        if v is None and isinstance(self.body, dict):
            v = self.body.get(key)
        return v

    @property
    def attempt_id(self):
        """The refused id (the `attempt` object's, else the top level's)."""
        return _wire_str(self._block_or_top('attempt_id'))

    @property
    def state(self):
        return _wire_str(self.attempt.get('state'))

    @property
    def tombstoned(self):
        return self._block_or_top('tombstone') is True

    @property
    def suppression(self):
        return Suppression.of(self.attempt)

    @property
    def server_instance(self):
        return _wire_str(self._block_or_top('server_instance'))


class InstanceMismatch(TraverseError):
    """409 with `state: "instance_mismatch"` (feature level 5): the request named another
    |expect_server_instance| than the process's |server_instance| (it restarted, or the
    request reached another server), so it was refused before anything ran or was
    registered. The client forgets the server_instance it had read, so that
    expect_server_instance='auto' reads the new one."""

    @property
    def expect_server_instance(self):
        return _wire_str((self.body or {}).get('expect_server_instance'))

    @property
    def server_instance(self):
        return _wire_str((self.body or {}).get('server_instance'))


def classify_409(body, sent=None, message=None):
    """The exception a 409 whose decoded JSON is |body| stands for (TraverseClient raises
    it; a caller reading the wire itself gets the same reading) -> an exception instance of
    status 409, never raised here. |sent|: the request as sent (its dict, or AttemptSent);
    |message|: the error text (default: the body's `error`).

    In this order (no server writes two of these, the order states which wins):
      1. an `attempt` object that states a `state`: AttemptConflict (the id's state);
      2. a top-level `state: "expired"`: AttemptExpired;
      3. a top-level `state: "instance_mismatch"`: InstanceMismatch only when |sent| carried
         expect_server_instance (a level-5 send) -- else a plain TraverseError, a body the
         request cannot have caused, which states nothing a ledger may act on;
      4. an `attempt` object without a state: AttemptConflict;
      5. anything else (no object, an `attempt` that is not one): TraverseError."""
    b = body if isinstance(body, dict) else {}
    if message is None:
        message = b.get('error') if isinstance(b.get('error'), str) else ''
    block = _obj(b, 'attempt')
    if block is not None and isinstance(block.get('state'), str):
        return AttemptConflict(409, message, body)
    state = b.get('state')
    if state == 'expired':
        return AttemptExpired(409, message, body)
    if state == 'instance_mismatch':
        s = sent if isinstance(sent, AttemptSent) else \
            AttemptSent.of(sent) if isinstance(sent, dict) else None
        if s is not None and s.expect_server_instance is not None:
            return InstanceMismatch(409, message, body)
        return TraverseError(409, message, body)
    if block is not None:
        return AttemptConflict(409, message, body)
    return TraverseError(409, message, body)


@dataclass
class TraverseResponse:
    envelope: dict                      # the response without results
    graphlets: List[Any]                # per seed: a Graphlet, or None for an error result
    errors: List[dict] = field(default_factory=list)   # [{index, error, result}]
    raw: Optional[dict] = None
    # deepen(): what the continuation request could not carry exactly (NextRequest.notes:
    # a loss budget conservative for some labels, a switch target left out, ...)
    notes: List[str] = field(default_factory=list)
    # traverse(): the attempt fields as sent (AttemptSent; None without attempt_id)
    sent: Optional[AttemptSent] = None
    # traverse(coordinates='auto'): what the automatic rule decided ({requested, reason,
    # memory_budget?, note?}, auto_coordinates()); None where it did not apply
    coordinates_auto: Optional[dict] = None

    @property
    def usage(self):
        """The response-level usage of a request with attempt_id (what the attempt
        consumed: work units, modelled memory, seeds, elapsed ms, the bound the server
        enforced); None without attempt_id."""
        return _obj(self.envelope, 'usage')


# ------------------------------------------------------------------ the release rule

@dataclass(frozen=True)
class ReleaseVerdict:
    """release_verdict()'s answer: |release| whether the attempt's capacity may be given back
    now, |early| whether that is before the attempt finished (on a tombstone), |code| a token
    for the condition that decided (released: `tombstone`, `finished`, `expired`, `clock`;
    not: the condition not met, see release_verdict()) and |why| that condition in words.
    True as a bool exactly when |release|.

    |assumptions|: what a release rests on beyond what the server holds (empty: nothing) --
    for `finished` and `expired`, that no copy of the request runs after the release (see
    release_verdict() for the tokens; for `expired` also `server_clock_step_back`: that the
    server's clock does not step back below not_after_ms), for `clock`,
    `uninterruptible_overrun`: that the attempt's run past its bound (up to its next
    delivery check) has ended. The release is the rule's
    either way; a ledger that must not see a request run twice, or must not overlap two
    attempts on one capacity, reads them."""
    release: bool
    early: bool
    code: str
    why: str
    assumptions: Tuple[str, ...] = ()

    def __bool__(self):
        return self.release


def _hold(code, why):
    return ReleaseVerdict(False, False, code, why)


_EPOCH = datetime.datetime(1970, 1, 1, tzinfo=datetime.timezone.utc)


def _epoch_ms(stamp):
    """An attempt state's time ("2026-10-05T04:33:14.159Z") as Unix epoch ms, else None."""
    if not isinstance(stamp, str):
        return None
    try:
        t = datetime.datetime.fromisoformat(stamp.replace('Z', '+00:00'))
    except ValueError:
        return None
    if t.tzinfo is None:
        t = t.replace(tzinfo=datetime.timezone.utc)
    return (t - _EPOCH) // datetime.timedelta(milliseconds=1)


def _finish_ms(state):
    """When the attempt a state or usage block describes finished: finished_at, else
    stopped_at (the response is written after it), else received_at -- never later than the
    finish, so a hold judged from it is never longer than the server's."""
    if not isinstance(state, dict):
        return None
    for key in ('finished_at', 'stopped_at', 'received_at'):
        ms = _epoch_ms(state.get(key))
        if ms is not None:
            return ms
    return None


def _finished(why, sent, state, attempts, skew, instance=None):
    """A `finished` release, with the assumptions it rests on. The rule holds a finished
    attempt's id against a replay only for one sent with not_after_ms, until not_after_ms +
    clock_skew_allowance_ms and at most tombstone_max_s after it finished (or after a later
    refused copy), on that process; otherwise "a finished state assumes that no copy of the
    request arrives after the attempt left retention" -- and a process that keeps nothing
    (retention_s 0) or a restarted one (no expect_server_instance to refuse the copy) runs
    that copy again. |state|: the state or usage block that says it finished (its
    server_instance and not_after_ms: the copy and the process that finished; |instance|
    where the block names no process)."""
    tokens, notes = [], []
    st = state if isinstance(state, dict) else {}
    inst = _wire_str(st.get('server_instance')) or instance
    naf = _wire_int(st.get('not_after_ms'))
    if naf is None:
        naf = _wire_int((_obj(st, 'usage') or {}).get('not_after_ms'))
    # the finish must be of this copy on the process it was pinned to -- else the pinned
    # process holds nothing of it and admits a copy that reaches it, or the hold that
    # finish left is judged against another copy's not_after_ms
    if sent.expect_server_instance is not None:
        if inst is None:
            tokens.append('server_instance_not_stated')
            notes.append('the finished state names no server_instance, so it may be another '
                         "process's than the %s the attempt was pinned to: no copy reaches "
                         'that one' % sent.expect_server_instance)
        elif inst != sent.expect_server_instance:
            tokens.append('other_server_instance')
            notes.append('the attempt finished on server_instance %s, not on %s that it was '
                         'pinned to, which holds nothing of it: no copy reaches %s'
                         % (inst, sent.expect_server_instance, sent.expect_server_instance))
    if sent.not_after_ms is not None and naf != sent.not_after_ms:
        tokens.append('other_copy_not_held')
        notes.append('the finished copy was sent with not_after_ms %r, not %d: its hold is '
                     'judged against that, so no copy sent with %d arrives after the attempt '
                     'left retention' % (naf, sent.not_after_ms, sent.not_after_ms))
    if sent.not_after_ms is None:
        tokens.append('sent_without_not_after_ms')
        notes.append('sent without not_after_ms, the id is held for the retention only: '
                     'no copy of the request arrives after the attempt left retention')
    retention = (attempts or {}).get('retention_s')
    if not attempts or isinstance(retention, bool) or not isinstance(retention, int):
        # without the capabilities' attempts block nothing says how long this server
        # keeps the id (retention_s 0 keeps nothing): an empty list must mean "held"
        tokens.append('server_hold_unchecked')
        notes.append('the server\'s attempts block (retention_s, tombstone_max_s) was not '
                     'given, so its hold of the id is not checked: no copy of the request '
                     'arrives after the server dropped the id')
    elif retention == 0:
        tokens.append('no_retention')
        notes.append('this server keeps no finished attempt (retention_s 0): no copy of '
                     'the request arrives later (one would run again)')
    elif attempts.get('tombstone_max_s') is None:
        # a server below feature level 5 holds a finished id while it is retained (by
        # age and by count, retention_count), never through not_after_ms
        tokens.append('no_finished_hold')
        notes.append('the attempts block states no tombstone_max_s (a server below feature '
                     'level 5): the id is held only while the attempt is retained '
                     '(retention_s, retention_count), so no copy of the request arrives '
                     'after it left retention')
    elif sent.not_after_ms is not None:
        finish = _finish_ms(state)
        hold = None if finish is None else finish + 1000 * attempts['tombstone_max_s']
        if hold is None or sent.not_after_ms + (skew or 0) > hold:
            tokens.append('beyond_tombstone_max')
            notes.append('not_after_ms + clock_skew_allowance_ms lies beyond '
                         'tombstone_max_s after the finish%s, so the id is not held '
                         'until then: no copy of the request arrives after it was dropped'
                         % ('' if finish is not None else ' (the finish is not stated)'))
    if sent.expect_server_instance is None:
        tokens.append('sent_without_expect_server_instance')
        notes.append('sent without expect_server_instance: no copy reaches a restarted '
                     'process, which holds nothing of the attempt and would run it')
    if tokens:
        why += '; assumed: ' + '; '.join(notes)
    return ReleaseVerdict(True, False, 'finished', why, tuple(tokens))


# The server compares the bound only at its delivery checks, it is not interrupted at it:
# past the bound an attempt runs the rest of the piece it was in, the walk up to its next
# poll that reads the clock (one in poll_stride), the stopped seed's finalisation and the
# building of its result up to the first delivery check -- the attempt's bound and the
# release rule (SPEC-labeled-traversal-core.md §6.8, §10.3, which the capabilities' bound and
# release_rule refer to) state it so, and that run has no stated length, so a clock release
# assumes it has ended
_CLOCK_ASSUMED = ("the attempt's run past its bound has ended: past it the attempt runs on "
                  "until its next delivery check (then the 503) or its handler's return -- "
                  'the rest of the piece it was in when the bound passed and, when its walk '
                  'had not stopped by then, the walk up to its next poll that reads the '
                  "clock, the stopped seed's finalisation and the building of its result up "
                  'to the first delivery check --, a run of no stated length')


def release_verdict(answer, sent, *, now_ms=None, clock_skew_allowance_ms=None,
                    bound_ms=None, attempts=None, hold_clock_while_running=True):
    """The capabilities' `release_rule` (feature level 5, SPEC §10.3) applied to one answer
    about the attempt |sent| (AttemptSent, or the request dict as sent) -> ReleaseVerdict.

    A ledger may release an attempt's capacity BEFORE it finished only on a 404 with
    `tombstone: true` and `covers_admission: true` (a cancel's, or GET /traverse/attempt's)
    from the same server_instance, for an attempt it sent with exactly that not_after_ms and
    with expect_server_instance equal to that server_instance, whose suppressed_until_ms
    reaches that not_after_ms + clock_skew_allowance_ms (code `tombstone`, early; the last
    is the server's own condition for covers_admission, checked again against the ledger's
    allowance -- without one, against not_after_ms itself). It releases otherwise:

      * on a finished state (code `finished`) -- GET /traverse/attempt's or a cancel's
        (`state: finished`, 200 or 404), the `attempt` object of a 409 that refused a copy
        (AttemptConflict: the id's state), the response (a TraverseResponse), or an
        error that carries the attempt's usage (it ran);
      * on an expired 409 (AttemptExpired, code `expired`): the server checks the id
        against what it holds before not_after_ms, so no copy of the attempt was registered
        there when it was judged, and every later copy is refused there while its clock
        does not step back below not_after_ms (the server's own condition: its check is
        strict). The 409's server_time_ms - not_after_ms is the step that survives; the
        release states `server_clock_step_back` when that margin is within
        clock_skew_allowance_ms (the step every release assumes the server's clock may
        take), or cannot be judged (no server_time_ms, no skew given);
      * once its own clock passes not_after_ms + clock_skew_allowance_ms + bound_ms (code
        `clock`: give |now_ms|, the ledger's Unix epoch ms, and the two allowances -- the
        capabilities' attempts block and the attempt's bound, usage.bound_ms or at most
        hard_cap_ms; |answer| may then be None). The clock takes every copy of the attempt
        as stopped at its bound, but past it the attempt runs on until its next delivery
        check or its handler's return -- the rest of the piece it was in and, when its walk
        had not stopped, the walk up to its next poll that reads the clock, the stopped
        seed's finalisation and the building of its result up to the first delivery check
        --, a run of no stated length: its release always states
        `uninterruptible_overrun`. While the last answer says the attempt is running or
        stopping -- past its bound, so in that run -- the clock releases nothing
        (hold_clock_while_running, the default: read its state again and release on
        finished); with hold_clock_while_running=False it releases, stating the same
        assumption.

    A `finished` or `expired` release states in |assumptions| what it rests on where the
    server does not hold the id against a replay: `sent_without_not_after_ms`,
    `server_hold_unchecked` (no |attempts| block was given, or it states no retention_s),
    `no_retention` (the server's retention_s is 0), `no_finished_hold` (the block states no
    tombstone_max_s: a server below feature level 5 holds a finished id only while it is
    retained), `beyond_tombstone_max` (not_after_ms + clock_skew_allowance_ms later than
    tombstone_max_s after the finish), `other_server_instance` /
    `server_instance_not_stated` (the finish is not stated by the process the attempt was
    pinned to), `other_copy_not_held` (the finished copy was sent with another not_after_ms)
    and `sent_without_expect_server_instance` (a copy reaching a restarted process); an
    `expired` release also `server_clock_step_back` (above).
    |attempts| is the capabilities' attempts block (or the whole GET /capabilities
    document), which also gives clock_skew_allowance_ms and, as the bound, hard_cap_ms when
    they are not given. Where the clock has passed too, the `clock` release (which assumes
    only `uninterruptible_overrun`) is returned instead.

    Not released, with the first condition not met: `answer_for_other_attempt` (the answer
    names another attempt_id, or none), `running` / `stopping` (a cancel's 200 with state
    stopping, a GET or a 409 saying so, release nothing), `not_tombstoned` (429: the
    answer's reason, tombstones_full or no_suppression, is in |why|), `no_tombstone`,
    `sent_without_not_after_ms`, `sent_without_expect_server_instance` (an attempt sent
    without either is never released early on a tombstone), `no_suppression_stated` (a
    server below feature level 5), `covers_admission_unstated` (not stated as a JSON bool),
    `covers_admission_false` (its reason in |why|), `other_server_instance`,
    `other_not_after_ms`, `suppression_short` (suppressed_until_ms short of not_after_ms +
    clock_skew_allowance_ms, or not stated), `not_a_release_answer` (a 409 for a tombstoned
    id, an instance_mismatch, an expired 409 the attempt cannot have caused, or an error
    before the attempt was registered: none is in the rule's list), `no_answer`."""
    sent = AttemptSent.of(sent)
    if sent is None:
        raise ValueError('release_verdict() is about an attempt: |sent| has no attempt_id')
    if isinstance(attempts, dict) and isinstance(attempts.get('attempts'), dict):
        attempts = attempts['attempts']
    if attempts:
        if clock_skew_allowance_ms is None:
            clock_skew_allowance_ms = attempts.get('clock_skew_allowance_ms')
        if bound_ms is None:
            bound_ms = attempts.get('hard_cap_ms')
    verdict = _answer_verdict(answer, sent, attempts, clock_skew_allowance_ms)
    if (verdict.release and not verdict.assumptions) or now_ms is None:
        return verdict
    if sent.not_after_ms is None or clock_skew_allowance_ms is None or bound_ms is None:
        return verdict
    limit = sent.not_after_ms + clock_skew_allowance_ms + bound_ms
    if now_ms <= limit:
        return verdict
    running = verdict.code in ('running', 'stopping')
    if running and hold_clock_while_running:
        # an answer saying running or stopping past the clock's limit is the
        # attempt's run past its bound still going -- what the clock assumes has ended. The
        # search service's sweep holds here too
        return _hold(verdict.code,
                     "%s; the ledger's clock (%d) has passed not_after_ms + "
                     "clock_skew_allowance_ms + bound_ms (%d), but the answer says the "
                     "attempt is %s: its run past its bound, up to its next delivery check, "
                     "is still going, a run of no stated length, so the clock releases "
                     "nothing on it (hold_clock_while_running): read its state again"
                     % (verdict.why, now_ms, limit, verdict.code))
    why = ("the ledger's clock (%d) has passed not_after_ms + clock_skew_allowance_ms + "
           "bound_ms (%d): no copy of the attempt can start any more, and every copy that "
           "started is past its bound" % (now_ms, limit))
    if running:
        why += ('; the last answer said %s, released because hold_clock_while_running is '
                'false' % verdict.code)
    return ReleaseVerdict(True, False, 'clock', why + '; assumed: ' + _CLOCK_ASSUMED,
                          ('uninterruptible_overrun',))


def _other_attempt(named, sent, what):
    if named is None:
        return _hold('answer_for_other_attempt',
                     '%s names no attempt_id: it is not known to be about attempt %r'
                     % (what, sent.attempt_id))
    return _hold('answer_for_other_attempt', '%s is about attempt %r, not %r'
                 % (what, named, sent.attempt_id))


def _answer_verdict(answer, sent, attempts=None, skew=None):
    if answer is None:
        return _hold('no_answer', 'no answer: only the clock releases it')
    if isinstance(answer, TraverseResponse):
        usage = answer.usage
        if not usage or usage.get('attempt_id') != sent.attempt_id:
            return _hold('answer_for_other_attempt',
                         'the response carries no usage of attempt %r' % sent.attempt_id)
        return _finished('the response: the attempt finished', sent, usage, attempts, skew)
    if isinstance(answer, AttemptConflict):
        if answer.attempt_id != sent.attempt_id:
            return _other_attempt(answer.attempt_id, sent, 'the 409\'s attempt')
        state = answer.state
        if state == 'finished':
            # the object is the id's state as GET /traverse/attempt answers it
            return _finished('a finished state (the attempt object of the 409 that refused '
                             'a copy)', sent, answer.attempt, attempts, skew,
                             answer.server_instance)
        if state in ('running', 'stopping'):
            return _hold(state, 'the 409 says the attempt is %s: release on its finished '
                         'state' % state)
        return _hold('not_a_release_answer',
                     'a 409 for a %s id is not in the release rule: nothing ran for this '
                     'request, and a copy may run or be running; read GET '
                     '/traverse/attempt, cancel the id, or wait for the clock'
                     % ('tombstoned' if answer.tombstoned else 'held'))
    if isinstance(answer, AttemptExpired):
        return _expired(answer, sent, skew)
    if isinstance(answer, TraverseError):
        if isinstance(answer, InstanceMismatch) or answer.status == 409:
            state = answer.body.get('state') if isinstance(answer.body, dict) else None
            return _hold('not_a_release_answer',
                         'a 409 (%s) is not in the release rule: nothing ran for this '
                         'request, but an earlier copy may run or be running; read GET '
                         '/traverse/attempt, cancel the id, or wait for the clock'
                         % (state or 'unclassified'))
        usage = answer.usage
        if usage and usage.get('attempt_id') == sent.attempt_id:
            return _finished('the response (HTTP %s with the attempt\'s usage): the attempt '
                             'ran and finished' % answer.status, sent, usage, attempts, skew)
        return _hold('not_a_release_answer',
                     'HTTP %s without the attempt\'s usage is not in the release rule'
                     % answer.status)
    if not isinstance(answer, AttemptAnswer):
        raise TypeError('release_verdict() reads an AttemptAnswer (cancel(), attempt()), a '
                        'TraverseResponse or a TraverseError, not %s' % type(answer).__name__)
    if answer.attempt_id != sent.attempt_id:
        # an answer that names no id is not known to be about this attempt
        return _other_attempt(answer.attempt_id, sent, 'the answer')
    if answer.status == 429:
        return _hold('not_tombstoned',
                     'a 429 releases nothing: the id was not tombstoned (%s)'
                     % (answer.reason or 'no reason stated'))
    if answer.state == 'finished':
        state = answer.attempt
        state = state if isinstance(state, dict) else dict(answer)
        return _finished('a finished state (%s %s)' % (answer.route, answer.status), sent,
                         state, attempts, skew, answer.server_instance)
    if answer.state in ('running', 'stopping'):
        return _hold(answer.state, 'the attempt is %s: release on its finished state'
                     % answer.state)
    if answer.status != 404 or answer.tombstone is not True:
        return _hold('no_tombstone', 'no tombstone (HTTP %s, tombstone %r): nothing is '
                     'promised about a later copy' % (answer.status, answer.tombstone))
    if sent.not_after_ms is None:
        return _hold('sent_without_not_after_ms', 'the attempt was sent without not_after_ms: '
                     'never released early on a tombstone')
    if sent.expect_server_instance is None:
        return _hold('sent_without_expect_server_instance',
                     'the attempt was sent without expect_server_instance: a copy reaching a '
                     'restarted process would run there, so it is never released early on a '
                     'tombstone')
    sup = answer.suppression
    if sup is None:
        return _hold('no_suppression_stated',
                     'the tombstone states no suppression (covers_admission): a server below '
                     'feature level 5 holds it for its retention only')
    if sup.covers_admission is None:
        return _hold('covers_admission_unstated',
                     'the tombstone does not state covers_admission as a JSON bool')
    if not sup.covers_admission:
        return _hold('covers_admission_false', 'covers_admission is false (%s)'
                     % (sup.covers_admission_reason or 'no reason stated'))
    if answer.server_instance != sent.expect_server_instance:
        return _hold('other_server_instance',
                     'the tombstone is server_instance %r\'s; the attempt was sent to %r'
                     % (answer.server_instance, sent.expect_server_instance))
    if sup.not_after_ms != sent.not_after_ms:
        return _hold('other_not_after_ms',
                     'the tombstone was judged against not_after_ms %r; the attempt was sent '
                     'with %r' % (sup.not_after_ms, sent.not_after_ms))
    # covers_admission restated with the ledger's allowance -- a server that states
    # it true with a hold short of not_after_ms + skew breaks SPEC §10.3, and a copy could
    # then be admitted once the hold is gone
    need = sent.not_after_ms + (skew or 0)
    if sup.suppressed_until_ms is None or sup.suppressed_until_ms < need:
        return _hold('suppression_short',
                     'the tombstone is held until %r, short of not_after_ms + '
                     'clock_skew_allowance_ms (%d): covers_admission true does not hold'
                     % (sup.suppressed_until_ms, need))
    return ReleaseVerdict(True, True, 'tombstone',
                          'a 404 with tombstone and covers_admission from server_instance %s '
                          'for not_after_ms %d: no copy of the attempt can start there'
                          % (answer.server_instance, sent.not_after_ms))


def _expired(answer, sent, skew=None):
    """An expired 409 to the attempt as sent. |skew|: clock_skew_allowance_ms (the
    call's, else the attempts block's), the step back of the server's clock every release
    assumes it survives."""
    if answer.attempt_id != sent.attempt_id:
        return _other_attempt(answer.attempt_id, sent, 'the expired 409')
    if sent.not_after_ms is None:
        return _hold('not_a_release_answer',
                     'an expired 409 to an attempt sent without not_after_ms cannot be its '
                     'answer')
    if answer.not_after_ms != sent.not_after_ms:
        return _hold('other_not_after_ms',
                     'the expired 409 was judged against not_after_ms %r; the attempt was '
                     'sent with %r' % (answer.not_after_ms, sent.not_after_ms))
    if sent.expect_server_instance is not None and \
            answer.server_instance != sent.expect_server_instance:
        return _hold('other_server_instance',
                     'the expired 409 is server_instance %r\'s; the attempt was sent to %r'
                     % (answer.server_instance, sent.expect_server_instance))
    why = ('an expired 409: the server refuses an id it holds (running, retained, held or '
           'tombstoned) before it checks not_after_ms, so no copy of the attempt was '
           'registered on server_instance %s when this one was judged, and every later copy '
           'is refused there while its clock reads later than not_after_ms (it has passed '
           'it). It settles nothing: an '
           'earlier copy may have run, finished and left retention, its usage unknown'
           % answer.server_instance)
    tokens, notes = [], []
    # Later copies are refused only while the server's clock reads later than not_after_ms
    # (its check is strict and adds no allowance): the release_rule assumes the clock does
    # not step back below it, and the 409's server_time_ms - not_after_ms is the step it
    # survives. The step back every release assumes is at most clock_skew_allowance_ms, so
    # only a margin above that is covered; a margin within it (a request that arrived just
    # late), or one that cannot be judged, is an assumption of this release of its own
    margin = None if answer.server_time_ms is None else \
        answer.server_time_ms - sent.not_after_ms
    if margin is None or skew is None or margin <= skew:
        tokens.append('server_clock_step_back')
        notes.append("this server's clock does not step back below not_after_ms, where a later "
                     'copy would be admitted again; %s' % (
                         'the 409 states no server_time_ms, so the step it survives is unknown'
                         if margin is None else
                         'server_time_ms - not_after_ms (%d ms) is the step it survives, %s'
                         % (margin, 'and no clock_skew_allowance_ms was given to compare it '
                            'with' if skew is None else 'within the clock_skew_allowance_ms '
                            '(%d ms) the server\'s clock may step back' % skew)))
    if sent.expect_server_instance is None:
        tokens.append('sent_without_expect_server_instance')
        notes.append('sent without expect_server_instance: no copy reached a restarted '
                     'process (or another server at the address) before not_after_ms')
    if tokens:
        why += '; assumed: ' + '; '.join(notes)
    return ReleaseVerdict(True, False, 'expired', why, tuple(tokens))
