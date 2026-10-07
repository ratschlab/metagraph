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
the client): a server below it refuses them with a 400 after it registered the
attempt_id, so a request with that id is refused (409) while that server retains it
(retention_s, retention_count; below level 5 nothing holds it longer) and runs again after
that.

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
way: an older binary registers its attempt_id before it answers the 400 (the id is
refused while that server retains it, and runs again after that), and a level-5 one
refuses it, 409 instance_mismatch, without registering it.
"""

import copy
import gzip
import http.client
import json
import re
import urllib.error
import urllib.parse
import urllib.request
import zlib

# the answers' readers and the release rule live in the light module attempts.py (LRG-G7):
# re-exported here, where they always were, as the same objects
from .attempts import (AttemptAnswer, AttemptAtBound, AttemptConflict, AttemptExpired,
                       AttemptSent, InstanceMismatch, ReleaseVerdict, ServerInitializing,
                       Suppression, TraverseError, TraverseResponse, classify_409,
                       release_verdict)

__all__ = ['TraverseClient', 'TraverseResponse', 'TraverseError', 'ServerInitializing',
           'AttemptAtBound', 'AttemptExpired', 'AttemptConflict', 'InstanceMismatch',
           'UnsupportedFeature', 'ResolveDeadline', 'AttemptAnswer', 'AttemptSent',
           'Suppression', 'ReleaseVerdict', 'release_verdict', 'classify_409', 'ACCEPT_ENCODING',
           'SUPPRESSION_LEVEL', 'COORDINATES_LEVEL', 'AUTO_COORDINATES_UNDER_MEMORY_BUDGET',
           'auto_coordinates', 'drop_coordinates_note', 'strip_coordinate_cap']

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


def strip_coordinate_cap(strategy):
    """|strategy| (in place) without output.max_coordinate_occurrences unless its
    output.coordinates is true: the server refuses the cap without coordinates (400,
    decision C-N6), so the library never sends it alone -- build_request(), next_request()
    (whose request carries the retrieval's echo, cap included) and the tools' requests
    (revision 2: dropping coordinates through a continuation keeps working). Defined here,
    not in coords.py (which re-exports it), so that importing the client loads no model
    (LRG-G7)."""
    out = strategy.get('output') if isinstance(strategy, dict) else None
    if isinstance(out, dict) and out.get('coordinates') is not True:
        out.pop('max_coordinate_occurrences', None)
    return strategy


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


class ResolveDeadline(TraverseError):
    """503 with `code: "deadline"` and no `usage`: a /resolve sent with bounds.time_budget_ms
    (resolve(time_budget_ms=)) whose answer could not be built and written within that budget
    -- its finalisation overran the reserve (SPEC §4.5). Nothing partial was sent and nothing
    was registered; send it again with a larger budget or a shorter sequence. Not a loading
    server (whose 503 has Retry-After and no code): never retry it as one."""


class UnsupportedFeature(ValueError):
    """Raised before sending: the request asks for a field the server does not state (its
    |feature_level| is below |needed|). Nothing was sent -- a server below the level refuses
    the field with a 400 after it registered the attempt_id of the request, whose id it
    then refuses while it retains it."""

    def __init__(self, feature, feature_level, needed):
        self.feature = feature
        self.feature_level = feature_level
        self.needed = needed
        super().__init__('%s needs a server at feature_level >= %d; this one states %s'
                         % (feature, needed, feature_level))


class _ResponseCut(ConnectionError):
    """No whole answer was received: the response was cut in transfer -- the connection
    closed before its Content-Length (http.client.IncompleteRead: a server drops a
    connection at its content timeout, a process dies with its buffers unsent, a proxy cuts
    it) -- or its status line or a header could not be read (another
    http.client.HTTPException). A ConnectionError, so an OSError like every other
    transport failure: callers that map those (the tools' backend_unreachable) map this one
    too, never a bare HTTPException that no handler expects (O35, the review of
    2026-10-06). Not a TraverseError of the status, even when one was read: a cut 2xx is
    no server answer, and a cut 409 or 503 body no longer tells its cases apart -- for an
    attempt the conservative reading is "unanswered" (cancel it, then judge it with
    release_verdict()). |reason| is the message (what the tools report)."""

    def __init__(self, message, status=None):
        super().__init__(message)
        self.reason = message
        self.status = status


def _cut(method, path, e, status=None, attempt_id=None):
    """The _ResponseCut of http.client.HTTPException |e| raised while |method| |path| was
    answered (|status|: the HTTP status, when its line was read; |attempt_id|: the
    request's, when it carried one). The message says the server was reached and may have
    run the request: the request went out before any answer was read, so a reader of
    "unreachable" alone would take a retry for free (the review of the P3 fixes)."""
    what = 'the response to %s %s' % (method, path) if status is None else \
        'the HTTP %d answer to %s %s' % (status, method, path)
    if isinstance(e, http.client.IncompleteRead):
        got = len(e.partial or b'')
        more = ', %d more expected' % e.expected if e.expected is not None else ''
        detail = 'cut in transfer (%d bytes read%s)' % (got, more)
    else:
        detail = 'unreadable (%s: %s)' % (type(e).__name__, e)
    then = 'cancel attempt %r, then judge it with release_verdict()' % attempt_id \
        if attempt_id is not None else 'a retry sends it again'
    return _ResponseCut('%s was %s: no whole answer was received, but the server was '
                        'reached and may have run the request (%s)' % (what, detail, then),
                        status)


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
            status = None
            try:
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
            except http.client.HTTPException as e:
                # a body cut short (IncompleteRead) or an unreadable status line or header:
                # neither a TraverseError nor an OSError, so it escaped every caller's
                # mapping (O35). The session= path needs no such case: requests raises its
                # ChunkedEncodingError (an OSError) for the same cut
                raise _cut(method, path, e, status,
                           payload.get('attempt_id') if isinstance(payload, dict) else None) \
                    from e
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
                # request whose attempt_id that binary registers and answers with a 400
                # every time
                self._forget_level()
            if status == 503:
                # Two 503s mean opposite things to a ledger: a loading index (Retry-After, no
                # usage: nothing registered) and an attempt stopped at its bound (usage, no
                # Retry-After: registered and consumed -- a retry under its id is refused
                # while the server holds the id, and runs again after). Told apart by the body,
                # not by the status alone (review of the stage-4 backend, F4)
                retry = hdrs.get('retry-after')
                if isinstance(out, dict) and 'usage' in out and not retry:
                    raise AttemptAtBound(status, message, out)
                if isinstance(out, dict) and out.get('code') == 'deadline' and not retry:
                    # a /resolve past its own bounds.time_budget_ms (SPEC §4.5): answered, not
                    # loading, so a ServerInitializing would invite a retry of the same budget
                    raise ResolveDeadline(status, message, out)
                raise ServerInitializing(status, message, out,
                                         int(retry) if retry and retry.isdigit() else None)
            if status == 409 and isinstance(out, dict):
                # three refusals, none of which ran anything: told apart by the body
                # (classify_409(): the attempt object's state first, an instance_mismatch
                # only to a request that named an instance)
                if out.get('state') == 'instance_mismatch':
                    # the process this client knew is gone: read its successor's instance
                    # before an 'auto' request names it again (forgetting it costs one read)
                    self._server_doc = None
                raise classify_409(out, payload if isinstance(payload, dict) else None,
                                   message)
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

    def resolve(self, sequence, *, labels=None, discover=None, select=None, time_budget_ms=None,
                **opts):
        """POST /resolve -> the answer, a dict as the server wrote it (SPEC §4).

        |time_budget_ms|: the request's deadline, sent as bounds.time_budget_ms (merged into a
        `bounds` given in |opts|; giving it there too with another value is a ValueError). It
        must exceed the server's finalisation reserve (capabilities()['resolve']['time_budget']
        ['finalize_reserve_ms'], 250 ms) and is lowered to its max_time_ms, stated in the
        answer's limits.clamped. The answer then carries `limits` and `stop`: stop None when the
        work completed, else {phase, reason: 'time', resolved_kmers, query_kmers, resolved_bp,
        query_bp, remainder_from_bp, message}, the answer being exactly the resolve of the
        sequence's first resolved_kmers k-mers -- resolve sequence[remainder_from_bp:] for the
        rest. An answer that could not be built within the budget raises ResolveDeadline (503).
        A server without the field (no `resolve` block in its capabilities) refuses it: a
        TraverseError 400 naming it. None (the default): sent exactly as before."""
        req = {'sequence': sequence}
        if labels is not None:
            req['labels'] = list(labels)
        if discover is not None:
            req['discover'] = discover
        if select is not None:
            req['select'] = select
        req.update(opts)
        if time_budget_ms is not None:
            bounds = dict(req.get('bounds') or {})
            if 'time_budget_ms' in bounds and bounds['time_budget_ms'] != time_budget_ms:
                raise ValueError('time_budget_ms=%r and bounds.time_budget_ms=%r differ'
                                 % (time_budget_ms, bounds['time_budget_ms']))
            bounds['time_budget_ms'] = time_budget_ms
            req['bounds'] = bounds
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
        server would register the attempt_id and refuse the field with a 400); 'auto' is
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
            # refused here, before anything is sent: the server's 400 would come after it
            # registered an attempt_id, refused then while the server holds it (the frozen
            # range of strategy.output.max_coordinate_occurrences)
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
