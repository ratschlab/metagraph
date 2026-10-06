"""TraverseClient: the HTTP side of "retrieve once, process locally" (§5).

The client only issues requests (/traverse, /resolve, /traverse/capabilities); every
follow-up is answered by the graphlet locally, and deepening is a new traversal whose
request the graphlet builds (Graphlet.next_request). Stdlib only (urllib): a
requests-like `session` (anything with .request(method, url, data=, headers=,
timeout=) returning .status_code/.content/.headers) can be injected instead, e.g. a
requests.Session or a test double.

The traversal routes write compact JSON and compress it on request; the client asks
explicitly for `Accept-Encoding: gzip, deflate` (gzip preferred by the server) and
decodes either itself, so the transport does not depend on the HTTP library.

A request may be an ATTEMPT of a ledger that reserved an allowance for it (DESIGN §14,
stage 4): `attempt_id` (unique per server process while it runs or is retained), with
`budget_id` and `locus_id` echoed. Every response to it carries a `usage` block
(TraverseResponse.usage), the server enforces a duration bound on it (usage.bound), and it
can be cancelled (cancel()) and its state read (attempt()) by id, also from another client.

Feature level 5 (SPEC §5, §10.3): a cancel names the cancelled request's `not_after_ms`,
so the tombstone it leaves is held until no copy of that request can start
(`covers_admission`); a request names the `server_instance` it is meant for
(`expect_server_instance`), so a copy that reaches a restarted process is refused (409
`instance_mismatch`); and the capabilities' `release_rule` says when a ledger may give an
attempt's capacity back -- release_verdict() applies it to an answer. Both fields are sent
only to a server that states `feature_level` >= 5 (read from GET /capabilities, or given to
the client): a server below it refuses them with a 400, and a 400 to a request with
attempt_id has already used that id up.

Feature level 6 (DESIGN §18, §26): record coordinates. traverse() asks for them
(strategy.output.coordinates) by itself for a support: trace strategy that does not set
them, where GET /traverse/capabilities states the coordinates block for the index with
supported true (supports_coordinates()) -- never on an index that cannot report them,
where the null form would only cost depth. Under a request memory budget
(bounds.max_memory_mb) only while the D4 gate passes (decision X-C8,
AUTO_COORDINATES_UNDER_MEMORY_BUDGET): a run's coordinate account is charged with the
walk, so a budgeted walk stops a little shallower with them; when it does not pass the
response says so (coordinates_auto, notes) and coordinates=True asks for them, and when
coordinates asked for that way share a memory stop the notes say how to drop them. It never
retries silently without them: an older server's 400 is raised as it is. A request whose
output.coordinates is not true never carries output.max_coordinate_occurrences, which the
server refuses without it (C-N6).

What GET /capabilities states is kept per server PROCESS: it is forgotten when an answer
names another server_instance than the one read (a restart at the same address), on a 409
instance_mismatch, and on a 400 that refuses a level-5 field as unknown (an older binary
now answers there: the level read is no longer the process's). A cancel refused that way
is sent again without not_after_ms, since that 400 stopped nothing. A request with
expect_server_instance that is the first to reach a replaced process cannot be saved that
way: an older binary uses its attempt_id up with the 400 (and a level-5 one refuses it,
409 instance_mismatch, without using it up).
"""

import copy
import datetime
import gzip
import json
import re
import urllib.error
import urllib.parse
import urllib.request
import zlib
from dataclasses import dataclass, field
from typing import Any, List, Optional, Tuple

from .coords import strip_coordinate_cap

__all__ = ['TraverseClient', 'TraverseResponse', 'TraverseError', 'ServerInitializing',
           'AttemptAtBound', 'AttemptExpired', 'AttemptConflict', 'InstanceMismatch',
           'UnsupportedFeature', 'AttemptAnswer', 'AttemptSent', 'Suppression',
           'ReleaseVerdict', 'release_verdict', 'ACCEPT_ENCODING', 'SUPPRESSION_LEVEL',
           'COORDINATES_LEVEL', 'AUTO_COORDINATES_UNDER_MEMORY_BUDGET', 'auto_coordinates',
           'drop_coordinates_note', 'strip_coordinate_cap']

ACCEPT_ENCODING = 'gzip, deflate'
# the feature level of the tombstone's suppression (a cancel with not_after_ms,
# covers_admission) and of expect_server_instance (SPEC §10.3)
SUPPRESSION_LEVEL = 5
# the feature level whose GET /traverse/capabilities states the coordinates block (SPEC
# §10.3): supports_coordinates() reads the block itself, not the level
COORDINATES_LEVEL = 6
# The D4 gate (decision X-C8): under a request memory budget the library asks for
# coordinates automatically only if the median depth at the budget's stop with them stays
# within 10% of the depth without -- each run's coordinate entry and occurrences are
# charged when the run is created, so a budgeted walk with them stops shallower. Measured
# at level 6 (DESIGN §26.5): on mini_refseq's 56 trace cells the median depth is 0.959-0.971
# of the depth without (a 3-4% drop), but the gate FAILS where it matters most: column labels
# with 16 or more chains per run (refseq33m's taxid columns) and seeds with many header
# labels need 2-4x the memory and can fail at depth 0 at a budget that completes without
# them (the regime test). The gate is therefore not passed, and the library does not ask
# under a memory budget: the answer says so and how to ask (_MEMORY_BUDGET_NOTE,
# coordinates=True)
AUTO_COORDINATES_UNDER_MEMORY_BUDGET = False
_MEMORY_BUDGET_NOTE = (
    'record coordinates were not requested automatically: the request has a memory budget '
    '(bounds.max_memory_mb), under which their account makes the walk stop shallower (the '
    'D4 gate, decision X-C8); pass coordinates=True to ask for them (output.coordinates), or '
    'drop the memory budget')
# where coordinates requested automatically under a memory budget share a memory stop
# (DESIGN §26.5: the single-lineage switch cells are the outliers; at 50% of their peak two
# of them fail at depth 0 with coordinates and not without)
_DROP_COORDINATES_NOTE = (
    'record coordinates were requested automatically under a memory budget '
    '(bounds.max_memory_mb), and %d seed(s) stopped on it with drop_coordinates offered: '
    'their coordinate account is part of what stopped them; send the request again with '
    'coordinates=False (output.coordinates false) to walk further, or raise the budget')
# A server refuses a field it does not know with a 400 naming it ("request: unknown field
# 'not_after_ms'", server.cpp and traverse.cpp). For the fields the client sends by the
# server's level, that 400 says the level the client kept is not the process's any more.
_UNKNOWN_FIELD = re.compile(r"unknown field '([^']*)'")
_LEVELLED_FIELDS = frozenset(('not_after_ms', 'expect_server_instance'))


# the largest max_coordinate_occurrences a server accepts (2^64 - 2: the wire contract)
_MAX_CAP = (1 << 64) - 2


def _valid_cap(v):
    if v == 'unlimited':
        return True
    return isinstance(v, int) and not isinstance(v, bool) and 1 <= v <= _MAX_CAP


def _memory_budget(strategy):
    v = ((strategy or {}).get('bounds') or {}).get('max_memory_mb')
    return isinstance(v, (int, float)) and not isinstance(v, bool) and v > 0


def auto_coordinates(client, strategy, graph=None, graph_path=None):
    """The automatic rule for record coordinates (decision C8, X-C8) -> (value, statement):
    value True (ask), or None (leave the strategy as given); statement None when the rule
    does not apply (the strategy sets output.coordinates, or is not support: trace in
    constrain mode), else {requested, reason[, memory_budget][, note]} -- reason
    'supported' (asked; memory_budget true when under a request memory budget, the D4 gate
    passing), 'unsupported' (the index states no coordinates block, or supported false) or
    'memory_budget' (a request memory budget with the gate not passing: not asked, |note|
    says how to ask). |client|: anything with supports_coordinates() (a TraverseClient);
    one without asks for none."""
    s = strategy or {}
    out = s.get('output') if isinstance(s.get('output'), dict) else {}
    if 'coordinates' in out:
        return None, None
    labels = s.get('labels') if isinstance(s.get('labels'), dict) else {}
    if s.get('support', 'kmer') != 'trace' or labels.get('mode', 'constrain') != 'constrain':
        return None, None
    probe = getattr(client, 'supports_coordinates', None)
    if probe is None or not probe(graph, graph_path):
        return None, {'requested': False, 'reason': 'unsupported'}
    if _memory_budget(s):
        if not AUTO_COORDINATES_UNDER_MEMORY_BUDGET:
            return None, {'requested': False, 'reason': 'memory_budget',
                          'note': _MEMORY_BUDGET_NOTE}
        return True, {'requested': True, 'reason': 'supported', 'memory_budget': True}
    return True, {'requested': True, 'reason': 'supported'}


def drop_coordinates_note(results, statement):
    """The note for coordinates the automatic rule asked for under a memory budget
    (|statement|: auto_coordinates()'s) where seed results of |results| (the response's
    raw `results`, graphlet or failed) stopped on memory with drop_coordinates offered,
    else None: it says how to walk further without them (DESIGN §26.5)."""
    if not statement or not statement.get('requested') or not statement.get('memory_budget'):
        return None
    n = 0
    for r in results or ():
        stop = r.get('resource_stop') if isinstance(r, dict) else None
        if isinstance(stop, dict) and 'drop_coordinates' in (stop.get('actions') or ()):
            n += 1
    return _DROP_COORDINATES_NOTE % n if n else None


def _unknown_field(status, message):
    """The field a 400's |message| refuses as unknown, else None."""
    if status != 400 or not isinstance(message, str):
        return None
    m = _UNKNOWN_FIELD.search(message)
    return m.group(1) if m else None


@dataclass(frozen=True)
class Suppression:
    """What a tombstone answer states (feature level 5): |suppressed_until_ms| the hold's
    last millisecond on the server's wall clock (inclusive, Unix epoch ms), |not_after_ms|
    the one it was judged against (None: none was), |covers_admission| whether no copy of
    that request can start on this server_instance, and when not,
    |covers_admission_reason| (`no_not_after_ms` | `beyond_tombstone_max`)."""
    suppressed_until_ms: int
    not_after_ms: Optional[int]
    covers_admission: bool
    covers_admission_reason: Optional[str] = None

    @classmethod
    def of(cls, body):
        """The suppression a tombstone answer (or a 409's `attempt`) states; None where it
        states none (no tombstone, or a server below feature level 5)."""
        if not isinstance(body, dict) or body.get('suppressed_until_ms') is None:
            return None
        return cls(body['suppressed_until_ms'], body.get('not_after_ms'),
                   body.get('covers_admission') is True, body.get('covers_admission_reason'))


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
    """The JSON answer of POST /traverse/cancel or GET /traverse/attempt (a dict, as before),
    with its HTTP |status|, the |route| (`cancel` | `attempt`) and, for a cancel, the body
    |sent| (whether not_after_ms went out with it). The typed reading:"""

    def __init__(self, body, status, route, sent=None):
        super().__init__(body)
        self.status = status
        self.route = route
        self.sent = sent

    @property
    def attempt_id(self):
        return self.get('attempt_id')

    @property
    def server_instance(self):
        return self.get('server_instance')

    @property
    def state(self):
        """running | stopping | finished | unknown."""
        return self.get('state')

    @property
    def cancelled(self):
        return self.get('cancelled')

    @property
    def tombstone(self):
        """True: the id is tombstoned (a 404); False: a cancel that was NOT tombstoned (a
        429); None: the answer says neither."""
        return self.get('tombstone')

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
        return self.get('attempt')

    @property
    def suppression(self):
        """The tombstone's Suppression (feature level 5), or None."""
        return Suppression.of(self)

    @property
    def finished(self):
        return self.get('state') == 'finished'

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
        self.usage = body.get('usage') if isinstance(body, dict) else None
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
    delivered. It was registered and ran -- its id is used up on this server (a retry under
    the same attempt_id is refused, 409) -- and |usage| states what it consumed. Not a
    loading server: never retry it as one."""


class AttemptExpired(TraverseError):
    """409 with `state: "expired"`: the request's not_after_ms had passed on the server's
    clock when its handler started, so it was NOT started — nothing ran, nothing was
    registered (no usage; GET /traverse/attempt answers 404 for its id). The body states
    `not_after_ms` and the server's `server_time_ms` (|not_after_ms|, |server_time_ms|,
    |server_instance|). Told apart from the other 409s (an attempt_id that is running,
    retained or tombstoned: AttemptConflict; another server_instance: InstanceMismatch) by its
    state."""

    @property
    def not_after_ms(self):
        return self.body.get('not_after_ms')

    @property
    def server_time_ms(self):
        return self.body.get('server_time_ms')

    @property
    def server_instance(self):
        return self.body.get('server_instance')


class AttemptConflict(TraverseError):
    """409 with `attempt`: the attempt_id is running, retained, held after it finished, or
    tombstoned by a cancel on this server, so the request was refused and nothing ran (no
    usage). |attempt| is that id's state; |tombstoned|: a cancel named the id first (state
    unknown, `tombstone: true`), and then |suppression| states the hold (feature level 5),
    judged against this request's own not_after_ms -- the refusal extended it. Not a release
    answer (release_verdict(): read the id's state, or cancel it)."""

    @property
    def attempt(self):
        return self.body.get('attempt') or {}

    @property
    def state(self):
        return self.attempt.get('state')

    @property
    def tombstoned(self):
        return self.attempt.get('tombstone') is True

    @property
    def suppression(self):
        return Suppression.of(self.attempt)

    @property
    def server_instance(self):
        return self.attempt.get('server_instance')


class InstanceMismatch(TraverseError):
    """409 with `state: "instance_mismatch"` (feature level 5): the request named another
    |expect_server_instance| than the process's |server_instance| (it restarted, or the
    request reached another server), so it was refused before anything ran or was
    registered. The client forgets the server_instance it had read, so that
    expect_server_instance='auto' reads the new one."""

    @property
    def expect_server_instance(self):
        return self.body.get('expect_server_instance')

    @property
    def server_instance(self):
        return self.body.get('server_instance')


class UnsupportedFeature(ValueError):
    """Raised before sending: the request asks for a field the server does not state (its
    |feature_level| is below |needed|). Nothing was sent -- a server below the level refuses
    the field with a 400, which for a request with attempt_id would use the id up."""

    def __init__(self, feature, feature_level, needed):
        self.feature = feature
        self.feature_level = feature_level
        self.needed = needed
        super().__init__('%s needs a server at feature_level >= %d; this one states %s'
                         % (feature, needed, feature_level))


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
        return self.envelope.get('usage')


def _decoded(data, encoding, status):
    """_decode_body(), a body its Content-Encoding does not decode (a truncated gzip or
    deflate stream) a TraverseError of the answer's status: zlib.error and EOFError are no
    error a caller of the client expects, and gzip's BadGzipFile, an OSError, read as a
    transport failure (backend_unreachable in the tools)."""
    try:
        return _decode_body(data, encoding)
    except (zlib.error, EOFError, OSError) as e:
        raise TraverseError(status, 'the response body could not be decoded (Content-Encoding '
                            '%s: %s)' % (encoding, type(e).__name__)) from None


def _decode_body(data, encoding):
    encoding = (encoding or '').strip().lower()
    if encoding == 'gzip':
        return gzip.decompress(data)
    if encoding == 'deflate':
        try:
            return zlib.decompress(data)            # zlib stream (what the server writes)
        except zlib.error:
            return zlib.decompress(data, -zlib.MAX_WBITS)   # raw deflate (RFC-lenient)
    return data


class TraverseClient:
    def __init__(self, host, port, api_path=None, *, session=None, timeout=900, release=None,
                 scheme='http', graph=None, feature_level=None):
        """|feature_level|: the server's, when the caller knows it; else it is read from GET
        /capabilities the first time a feature-level-5 field is to be sent (and never for a
        request without one, so the requests of earlier levels go out as they always did).
        Once GET /capabilities has been read (expect_server_instance='auto' reads it for the
        instance) the level it states wins over the given one, and a given level of 5 or
        more is dropped when the server refuses a level-5 field as unknown: either is the
        process's own statement."""
        self.host = host
        self.port = port
        self.base = '%s://%s:%s' % (scheme, host, port)
        if api_path:
            self.base += '/' + api_path.strip('/')
        self.session = session
        self.timeout = timeout
        self.release = release
        self.graph = graph
        self.feature_level = feature_level
        self._server_doc = None
        # (graph, graph_path) -> whether GET /traverse/capabilities states coordinates
        self._coord_support = {}

    # ---------------------------------------------------------------- transport

    def _request(self, method, path, payload=None, answers=(), with_status=False):
        """|answers|: non-2xx statuses whose JSON body is an answer, returned rather than
        raised (cancel and attempt: 404 states that nothing runs under the id).
        |with_status|: -> (HTTP status, body)."""
        url = self.base + path
        body = None if payload is None else json.dumps(payload).encode('utf-8')
        headers = {'Accept-Encoding': ACCEPT_ENCODING}
        if body is not None:
            headers['Content-Type'] = 'application/json'
        if self.session is not None:
            resp = self.session.request(method, url, data=body, headers=headers,
                                        timeout=self.timeout)
            status = resp.status_code
            hdrs = {k.lower(): v for k, v in dict(resp.headers or {}).items()}
            data = resp.content
            # requests decodes Content-Encoding itself; a raw double does not
            if getattr(resp, 'raw_encoded', False):
                data = _decoded(data, hdrs.get('content-encoding'), status)
        else:
            req = urllib.request.Request(url, data=body, headers=headers, method=method)
            try:
                with urllib.request.urlopen(req, timeout=self.timeout) as r:
                    status = r.status
                    hdrs = {k.lower(): v for k, v in r.headers.items()}
                    raw = r.read()
                data = _decoded(raw, hdrs.get('content-encoding'), status)
            except urllib.error.HTTPError as e:
                status = e.code
                hdrs = {k.lower(): v for k, v in (e.headers or {}).items()}
                try:
                    raw = e.read() or b''
                finally:
                    e.close()
                data = _decoded(raw, hdrs.get('content-encoding'), status)
        out = None
        if isinstance(data, (bytes, bytearray)):
            try:
                text = data.decode('utf-8')
            except UnicodeDecodeError:
                # not UTF-8, so no JSON answer: an error page is still classified by its
                # status below, a 2xx is 'not JSON' -- never a bare UnicodeDecodeError, which
                # no caller expects from a server's answer (VMD-07)
                text = data.decode('utf-8', 'replace')
                data = None
        else:
            text = data
        if data is not None:
            try:
                out = json.loads(text) if text else {}
            except ValueError:
                out = None
        if isinstance(out, dict) and self._server_doc and path != '/capabilities':
            self._check_instance(out)
        if status in answers and isinstance(out, dict):
            return (status, out) if with_status else out
        if not 200 <= status < 300:
            message = out.get('error') if isinstance(out, dict) else (text or '')[:500]
            refused = _unknown_field(status, message)
            if status == 400 and isinstance(message, str) and 'coordinate' in message:
                # a server that refuses the coordinates fields is not the one whose
                # capabilities were read (a rollback at the same address): read them again
                self._coord_support.clear()
            if refused in _LEVELLED_FIELDS and isinstance(payload, dict) and refused in payload:
                # a field sent by the server's level is unknown to the process answering:
                # an older binary replaced the one whose level the client kept (a rollback
                # at the same address; it cannot answer instance_mismatch). Kept, the level
                # would send the field again: a cancel that never stops anything, an 'auto'
                # request that uses its attempt_id up every time
                self._forget_level()
            if status == 503:
                # Two 503s mean opposite things to a ledger: a loading index (Retry-After, no
                # usage: nothing registered) and an attempt stopped at its bound (usage, no
                # Retry-After: registered, consumed, its id used up). Told apart by the body,
                # not by the status alone (review of the stage-4 backend, F4)
                retry = hdrs.get('retry-after')
                if isinstance(out, dict) and 'usage' in out and not retry:
                    raise AttemptAtBound(status, message, out)
                raise ServerInitializing(status, message, out,
                                         int(retry) if retry and retry.isdigit() else None)
            if status == 409 and isinstance(out, dict):
                # three refusals, none of which ran anything: told apart by the body
                if out.get('state') == 'expired':
                    raise AttemptExpired(status, message, out)
                if out.get('state') == 'instance_mismatch':
                    # the process this client knew is gone: read its successor's instance
                    # before an 'auto' request names it again
                    self._server_doc = None
                    raise InstanceMismatch(status, message, out)
                if 'attempt' in out:
                    raise AttemptConflict(status, message, out)
            raise TraverseError(status, message, out)
        if out is None:
            raise TraverseError(status, 'the response is not JSON', text[:500])
        return (status, out) if with_status else out

    # ---------------------------------------------------------------- feature level

    def _server(self, refresh=False):
        if self._server_doc is None or refresh:
            self._server_doc = self.server_capabilities()
        return self._server_doc

    def _check_instance(self, out):
        """Forget the kept GET /capabilities when an answer names another server_instance:
        the process was restarted (or replaced by another binary) at the same address, and
        its instance and maybe its level are not the ones read. Every answer about an
        attempt names its process (a cancel's, a state's, a usage block, a 409's attempt),
        so a restart is noticed before the next 'auto' request names the old instance."""
        known = (self._server_doc.get('attempts') or {}).get('server_instance')
        if known is None:
            return
        inst = out.get('server_instance')
        if inst is None:
            for key in ('usage', 'attempt'):
                sub = out.get(key)
                if isinstance(sub, dict) and sub.get('server_instance') is not None:
                    inst = sub['server_instance']
                    break
        if inst is not None and inst != known:
            self._server_doc = None

    def _forget_level(self):
        self._server_doc = None
        if self.feature_level is not None and self.feature_level >= SUPPRESSION_LEVEL:
            self.feature_level = None

    def server_feature_level(self, refresh=False):
        """The server's feature_level: GET /capabilities' once read, else the one given to
        the client, else GET /capabilities' read now (|refresh| reads it again). Kept per
        server process: forgotten when an answer names another server_instance, on an
        InstanceMismatch and when the server refuses a level-5 field as unknown (see the
        module's docstring). A server without that route (before feature level 1) is 0."""
        if not refresh:
            if self._server_doc is not None:
                return int(self._server_doc.get('feature_level') or 0)
            if self.feature_level is not None:
                return self.feature_level
        try:
            doc = self._server(refresh=True)
        except TraverseError as e:
            if e.status != 404:
                raise
            self._server_doc = doc = {}         # kept: no route is the answer
        return int(doc.get('feature_level') or 0)

    def server_instance(self, refresh=False):
        """The server process's `server_instance` (16 hex, new at every start) from GET
        /capabilities' attempts block, read once and kept for that process (forgotten as
        server_feature_level() is, or read again with |refresh|); None where the server
        states none."""
        return ((self._server(refresh) or {}).get('attempts') or {}).get('server_instance')

    def refresh_capabilities(self):
        """Read GET /capabilities again (a restarted server: a new server_instance, maybe
        another feature level)."""
        return self._server(refresh=True)

    def _level(self, needed):
        level = self.server_feature_level()
        return level, level >= needed

    # ---------------------------------------------------------------- routes

    def capabilities(self, graph=None, graph_path=None):
        """GET /traverse/capabilities: what one index supports. On a multi-graph server name
        it: |graph| (the client's own graph by default) and, when the name spans several
        graphs, |graph_path|."""
        graph = graph or self.graph
        query = []
        if graph:
            query.append('graph=' + urllib.parse.quote(graph, safe=''))
        if graph_path:
            query.append('graph_path=' + urllib.parse.quote(graph_path, safe=''))
        return self._request('GET', '/traverse/capabilities'
                             + ('?' + '&'.join(query) if query else ''))

    def supports_coordinates(self, graph=None, graph_path=None, refresh=False):
        """Whether the index (|graph|, |graph_path| on a multi-graph server) reports record
        coordinates: its GET /traverse/capabilities states the coordinates block (feature
        level 6: on every index, `supported` saying whether this one records them) with
        supported true, and supports_trace. Read once per (graph, graph_path) and kept
        (|refresh| reads it again; a 400 that refuses a coordinates field forgets it). A
        server whose capabilities route fails with an HTTP error is taken as not supporting
        them."""
        key = (graph or self.graph, graph_path)
        got = None if refresh else self._coord_support.get(key)
        if got is None:
            try:
                caps = self.capabilities() if graph is None and graph_path is None \
                    else self.capabilities(graph, graph_path)
            except TraverseError:
                return False
            block = caps.get('coordinates') if isinstance(caps, dict) else None
            got = isinstance(block, dict) and block.get('supported') is True \
                and caps.get('supports_trace') is True
            self._coord_support[key] = got
        return got

    def auto_coordinates(self, strategy, graph=None, graph_path=None):
        """auto_coordinates(self, ...): the automatic rule for record coordinates."""
        return auto_coordinates(self, strategy, graph, graph_path)

    def server_capabilities(self):
        """GET /capabilities: the server-wide document — routes and features, feature_level,
        algorithm_version, mode (single | multi) and the graph names, the attempts block and
        how deadlines are checked (deadline_check). Answered while a single index loads
        (`ready: false`)."""
        return self._request('GET', '/capabilities')

    def resolve(self, sequence, *, labels=None, discover=None, select=None, **opts):
        req = {'sequence': sequence}
        if labels is not None:
            req['labels'] = list(labels)
        if discover is not None:
            req['discover'] = discover
        if select is not None:
            req['select'] = select
        req.update(opts)
        if self.graph and 'graph' not in req:
            req['graph'] = self.graph
        return self._request('POST', '/resolve', req)

    def traverse_raw(self, request):
        return self._request('POST', '/traverse', request)

    def cancel(self, attempt_id, wait_ms=None, not_after_ms=None):
        """POST /traverse/cancel -> AttemptAnswer (the server's answer, a dict, with its HTTP
        status): `cancelled: True` and `state` stopping (or finished, within |wait_ms|) when
        the attempt was asked to stop; `cancelled: False` with `state` finished or unknown
        (HTTP 404) when nothing runs under the id. An unknown id is tombstoned (`tombstone:
        True`): a request arriving later with it is refused (409, AttemptConflict) and never
        runs there while the tombstone is held. When the server holds its maximum of
        tombstones (429, reason `tombstones_full`) or keeps none (429, `no_suppression`) the
        id is NOT tombstoned: nothing is promised.

        |not_after_ms| (feature level 5): the not_after_ms of the request being cancelled.
        The tombstone is then held until no copy of that request can start, within the
        server's tombstone_max_s, and the answer states so (`suppression`: covers_admission);
        release_verdict() says whether the attempt's capacity may be released on it. Sent
        only to a server at feature level 5 or above: to an older one the cancel goes out
        without it (the cancel still stops the attempt; |sent| shows what went out, and the
        answer states no suppression, which release_verdict() never releases on). A process
        that refuses it as unknown although the level read was 5 (an older binary now
        answers at the address) gets the cancel again without it, once."""
        payload = {'attempt_id': attempt_id}
        if wait_ms is not None:
            payload['wait_ms'] = wait_ms
        if not_after_ms is not None and self._level(SUPPRESSION_LEVEL)[1]:
            payload['not_after_ms'] = not_after_ms
        try:
            status, out = self._request('POST', '/traverse/cancel', payload,
                                        answers=(404, 429), with_status=True)
        except TraverseError as e:
            if 'not_after_ms' not in payload or \
                    _unknown_field(e.status, e.message) != 'not_after_ms':
                raise
            # the 400 refused the cancel before it was read: nothing was stopped and nothing
            # tombstoned, so sending it again is safe, and stopping the attempt is the
            # cancel's first job (the request it cancels would otherwise run). _request has
            # forgotten the level that let not_after_ms out
            payload = {k: v for k, v in payload.items() if k != 'not_after_ms'}
            status, out = self._request('POST', '/traverse/cancel', payload,
                                        answers=(404, 429), with_status=True)
        return AttemptAnswer(out, status, 'cancel', payload)

    def attempt(self, attempt_id):
        """GET /traverse/attempt/{attempt_id} -> AttemptAnswer: the attempt's state (running |
        stopping | finished, the reason, when it stopped, the response written, its usage), or
        `state: unknown` (HTTP 404) once it is no longer retained -- for a tombstoned id with
        `tombstone: True` and, at feature level 5, its suppression (reading it does not extend
        it)."""
        status, out = self._request('GET', '/traverse/attempt/' + attempt_id, answers=(404,),
                                    with_status=True)
        return AttemptAnswer(out, status, 'attempt')

    def build_request(self, seeds, strategy=None, *, detail='graphlet', timing=True,
                      graph=None, attempt_id=None, budget_id=None, locus_id=None,
                      not_after_ms=None, expect_server_instance=None, coordinates=None,
                      max_coordinate_occurrences=None):
        """The /traverse request (a dict). |expect_server_instance| (with |attempt_id| only,
        feature level 5): the server_instance the attempt is meant for -- a string, or 'auto'
        for the server's own (server_instance(), kept for that process). An explicit one to
        a server below feature level 5 raises UnsupportedFeature (nothing is sent: that
        server would refuse the field with a 400 that uses the attempt_id up); 'auto' is
        then left out. Either way AttemptSent.of(request) states what is sent.

        |coordinates| (feature level 6): True sets strategy.output.coordinates, False
        removes it (the same as false to the server, and valid to one that does not know
        the field), None leaves the strategy as given (pure: no capabilities are read here;
        traverse() applies the automatic rule). |max_coordinate_occurrences|: the cap per
        occurrence list (an int >= 1 or "unlimited"). The cap is removed whenever the
        resulting output.coordinates is not true: the server refuses it alone (C-N6)."""
        if expect_server_instance is not None and attempt_id is None:
            raise ValueError('expect_server_instance is sent with attempt_id only')
        if coordinates is not None and not isinstance(coordinates, bool):
            # by type, not by ==: 0 and 1 equal False and True, but the dispatch below is
            # by identity, so they were taken for None (coordinates=0 kept coordinates,
            # coordinates=1 dropped an explicit cap)
            raise ValueError('coordinates is True, False or None, not %r' % (coordinates,))
        if max_coordinate_occurrences is not None and not _valid_cap(max_coordinate_occurrences):
            # refused here, before anything is sent: the server's 400 would use an
            # attempt_id up (the frozen range of strategy.output.max_coordinate_occurrences)
            raise ValueError('max_coordinate_occurrences is an integer in [1, %d] or '
                             '"unlimited", not %r' % (_MAX_CAP, max_coordinate_occurrences))
        norm = []
        for s in seeds:
            norm.append({'sequence': s} if isinstance(s, str) else dict(s))
        strategy = copy.deepcopy(strategy) if strategy else {}
        output = strategy.setdefault('output', {})
        output['detail'] = detail
        output['timing'] = timing
        if coordinates is True:
            output['coordinates'] = True
        elif coordinates is False:
            output.pop('coordinates', None)
        if max_coordinate_occurrences is not None:
            output['max_coordinate_occurrences'] = max_coordinate_occurrences
        strip_coordinate_cap(strategy)
        req = {'seeds': norm, 'strategy': strategy}
        if self.release:
            req['release'] = self.release
        graph = graph or self.graph
        if graph:
            req['graph'] = graph
        for key, value in (('attempt_id', attempt_id), ('budget_id', budget_id),
                           ('locus_id', locus_id), ('not_after_ms', not_after_ms)):
            if value is not None:
                req[key] = value
        if expect_server_instance is not None:
            level, ok = self._level(SUPPRESSION_LEVEL)
            if expect_server_instance == 'auto':
                instance = self.server_instance() if ok else None
                # reading the instance read the level the process states, which wins over
                # a level given to the client
                if instance is not None and self._level(SUPPRESSION_LEVEL)[1]:
                    req['expect_server_instance'] = instance
            elif not ok:
                raise UnsupportedFeature('expect_server_instance', level, SUPPRESSION_LEVEL)
            else:
                req['expect_server_instance'] = expect_server_instance
        return req

    def traverse(self, seeds, strategy=None, *, detail='graphlet', timing=True, graph=None,
                 attempt_id=None, budget_id=None, locus_id=None, not_after_ms=None,
                 expect_server_instance=None, coordinates='auto',
                 max_coordinate_occurrences=None):
        """-> TraverseResponse(envelope, graphlets, errors). A seed whose permitted set
        could not be derived is an error result ({seed, error}), not an exception: the
        other seeds' graphlets are still returned (so is a seed an attempt never started,
        cancelled or past its bound: resource_stop phase not_started). A request error is
        raised (TraverseError; with attempt_id its usage is in the error's .usage), an
        attempt whose response was stopped at its bound raises AttemptAtBound, a request
        whose |not_after_ms| (Unix epoch ms) had passed on the server's clock raises
        AttemptExpired (not started), an attempt_id the server holds (running, retained, held
        or tombstoned) AttemptConflict, and an |expect_server_instance| (build_request())
        that is not the server's InstanceMismatch. The response's .sent, and the .sent of
        any exception raised once the request was built (a TraverseError, and a transport
        error too: a timeout, a reset connection), state the attempt fields as sent (for
        release_verdict()); an error raised while building it (UnsupportedFeature) sent
        nothing.

        |coordinates| (feature level 6): 'auto' (the default) asks for record coordinates
        for a support: trace strategy that does not set output.coordinates, where the index
        reports them (supports_coordinates()) -- under a request memory budget only while
        the D4 gate passes (X-C8, AUTO_COORDINATES_UNDER_MEMORY_BUDGET); what it decided is
        .coordinates_auto, and .notes states a decision not to ask under a memory budget,
        or, where coordinates it asked for under one share a memory stop, how to drop them
        (drop_coordinates_note()). True / False / None as in build_request(). A server that
        refuses the field (an older one: 400) raises; nothing is retried without it."""
        auto = None
        if coordinates == 'auto':
            coordinates, auto = auto_coordinates(self, strategy, graph)
        req = self.build_request(seeds, strategy, detail=detail, timing=timing, graph=graph,
                                 attempt_id=attempt_id, budget_id=budget_id, locus_id=locus_id,
                                 not_after_ms=not_after_ms,
                                 expect_server_instance=expect_server_instance,
                                 coordinates=coordinates,
                                 max_coordinate_occurrences=max_coordinate_occurrences)
        sent = AttemptSent.of(req)
        try:
            resp = self.response(self.traverse_raw(req))
        except Exception as e:
            # an unanswered attempt (the socket timed out, the connection was reset) is the
            # one a ledger must cancel and judge with release_verdict(), and with 'auto' only
            # the client knows whether an instance went out
            try:
                e.sent = sent
            except (AttributeError, TypeError):     # an exception type without attributes
                pass
            raise
        resp.sent = sent
        resp.coordinates_auto = auto
        if auto is not None and auto.get('note'):
            resp.notes.append(auto['note'])
        dropped = drop_coordinates_note((resp.raw or {}).get('results'), auto)
        if dropped is not None:
            resp.notes.append(dropped)
        return resp

    @staticmethod
    def response(out):
        from .parser import from_response
        envelope = {k: v for k, v in out.items() if k != 'results'}
        graphlets, errors = [], []
        for i, r in enumerate(out.get('results', [])):
            if 'graphlet' in r:
                graphlets.append(from_response(r, out))
            else:
                graphlets.append(None)
                errors.append({'index': i, 'error': r.get('error'), 'result': r})
        return TraverseResponse(envelope, graphlets, errors, out)

    def deepen(self, graphlet, arm, leaves, bp=None, **overrides):
        """Iterative deepening stays a backend call: the graphlet derives the
        continuations and builds the request (Graphlet.next_request), the server runs it.
        The response's notes carry the request's: how the continuation may differ from
        one uninterrupted walk (a loss budget conservative for the labels that ended at a
        lower loss, ...). Walks whose continuations carry different labels cannot share a
        request (IncompatibleContinuations): deepen them one walk at a time."""
        if self.graph and 'graph' not in overrides:
            overrides['graph'] = self.graph
        req = graphlet.next_request(arm, leaves, bp, **overrides)
        out = self.response(self.traverse_raw(req))
        out.notes = list(getattr(req, 'notes', ()))
        return out


# ------------------------------------------------------------------ the release rule

@dataclass(frozen=True)
class ReleaseVerdict:
    """release_verdict()'s answer: |release| whether the attempt's capacity may be given back
    now, |early| whether that is before the attempt finished (on a tombstone), |code| a token
    for the condition that decided (released: `tombstone`, `finished`, `clock`; not:
    the condition not met, see release_verdict()) and |why| that condition in words. True
    as a bool exactly when |release|.

    |assumptions|: for a `finished` release, the rule's assumptions it rests on instead of
    the server holding the id against a copy of the request (empty: the server holds it
    until no copy can start) -- `sent_without_not_after_ms`, `no_retention`,
    `beyond_tombstone_max`, `sent_without_expect_server_instance` (see release_verdict()).
    The release is the rule's either way; a ledger that must not see a request run twice
    reads them."""
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


def _finished(why, sent, state, attempts, skew):
    """A `finished` release, with the assumptions it rests on. The rule holds a finished
    attempt's id against a replay only for one sent with not_after_ms, until not_after_ms +
    clock_skew_allowance_ms and at most tombstone_max_s after it finished, on that process;
    otherwise "a finished state assumes that no copy of the request arrives after the
    attempt left retention" -- and a process that keeps nothing (retention_s 0) or a
    restarted one (no expect_server_instance to refuse the copy) runs that copy again."""
    tokens, notes = [], []
    if sent.not_after_ms is None:
        tokens.append('sent_without_not_after_ms')
        notes.append('sent without not_after_ms, the id is held for the retention only: '
                     'no copy of the request arrives after the attempt left retention')
    if attempts:
        if attempts.get('retention_s') == 0:
            tokens.append('no_retention')
            notes.append('this server keeps no finished attempt (retention_s 0): no copy of '
                         'the request arrives later (one would run again)')
        elif sent.not_after_ms is not None and attempts.get('tombstone_max_s') is not None:
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


def release_verdict(answer, sent, *, now_ms=None, clock_skew_allowance_ms=None,
                    bound_ms=None, attempts=None):
    """The capabilities' `release_rule` (feature level 5, SPEC §10.3) applied to one answer
    about the attempt |sent| (AttemptSent, or the request dict as sent) -> ReleaseVerdict.

    A ledger may release an attempt's capacity BEFORE it finished only on a 404 with
    `tombstone: true` and `covers_admission: true` (a cancel's, or GET /traverse/attempt's)
    from the same server_instance, for an attempt it sent with exactly that not_after_ms and
    with expect_server_instance equal to that server_instance (code `tombstone`, early). It
    releases otherwise only on a finished state -- GET /traverse/attempt's or a cancel's
    (`state: finished`, 200 or 404) or the response (a TraverseResponse, or an error that
    carries the attempt's usage: it ran) -- (code `finished`), or once its own clock passes
    not_after_ms + clock_skew_allowance_ms + bound_ms (code `clock`: give |now_ms|, the
    ledger's Unix epoch ms, and the two allowances -- the capabilities' attempts block and
    the attempt's bound, usage.bound_ms or at most hard_cap_ms; |answer| may then be None).

    A `finished` release states in |assumptions| what it rests on where the server does not
    hold the id against a replay: `sent_without_not_after_ms`, `no_retention` (the server's
    retention_s is 0), `beyond_tombstone_max` (not_after_ms + clock_skew_allowance_ms later
    than tombstone_max_s after the finish) and `sent_without_expect_server_instance` (a copy
    reaching a restarted process). The two about the server are judged only with
    |attempts|, the capabilities' attempts block (or the whole GET /capabilities document),
    which also gives clock_skew_allowance_ms and, as the bound, hard_cap_ms when they are
    not given. Where the clock has passed too, the `clock` release (which assumes none of
    them) is returned instead.

    Not released, with the first condition not met: `answer_for_other_attempt`, `running` /
    `stopping` (a cancel's 200 with state stopping releases nothing), `not_tombstoned` (429:
    the answer's reason, tombstones_full or no_suppression, is in |why|), `no_tombstone`,
    `sent_without_not_after_ms`, `sent_without_expect_server_instance` (an attempt sent
    without either is never released early on a tombstone), `no_suppression_stated` (a
    server below feature level 5), `covers_admission_false` (its reason in |why|),
    `other_server_instance`, `other_not_after_ms`, `not_a_release_answer` (a 409 --
    duplicate, tombstoned, expired or instance_mismatch --, or an error before the attempt
    was registered: none is in the rule's list), `no_answer`."""
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
    if now_ms > limit:
        return ReleaseVerdict(True, False, 'clock',
                              "the ledger's clock (%d) has passed not_after_ms + "
                              "clock_skew_allowance_ms + bound_ms (%d): the attempt cannot "
                              "start any more and has reached its bound" % (now_ms, limit))
    return verdict


def _answer_verdict(answer, sent, attempts=None, skew=None):
    if answer is None:
        return _hold('no_answer', 'no answer: only the clock releases it')
    if isinstance(answer, TraverseResponse):
        usage = answer.usage
        if not usage or usage.get('attempt_id') != sent.attempt_id:
            return _hold('answer_for_other_attempt',
                         'the response carries no usage of attempt %r' % sent.attempt_id)
        return _finished('the response: the attempt finished', sent, usage, attempts, skew)
    if isinstance(answer, TraverseError):
        if isinstance(answer, (AttemptExpired, InstanceMismatch, AttemptConflict)):
            return _hold('not_a_release_answer',
                         'a 409 (%s) is not in the release rule: nothing ran for this '
                         'request, but an earlier copy may run or be running; read GET '
                         '/traverse/attempt, cancel the id, or wait for the clock'
                         % (answer.body or {}).get('state', 'attempt_id held'))
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
    if answer.attempt_id is not None and answer.attempt_id != sent.attempt_id:
        return _hold('answer_for_other_attempt', 'the answer is about attempt %r, not %r'
                     % (answer.attempt_id, sent.attempt_id))
    if answer.status == 429:
        return _hold('not_tombstoned',
                     'a 429 releases nothing: the id was not tombstoned (%s)'
                     % (answer.reason or 'no reason stated'))
    if answer.state == 'finished':
        state = answer.get('attempt') if answer.route == 'cancel' else None
        state = state if isinstance(state, dict) else dict(answer)
        return _finished('a finished state (%s %s)' % (answer.route, answer.status), sent,
                         state, attempts, skew)
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
    return ReleaseVerdict(True, True, 'tombstone',
                          'a 404 with tombstone and covers_admission from server_instance %s '
                          'for not_after_ms %d: no copy of the attempt can start there'
                          % (answer.server_instance, sent.not_after_ms))
