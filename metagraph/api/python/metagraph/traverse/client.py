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
"""

import copy
import gzip
import json
import urllib.error
import urllib.parse
import urllib.request
import zlib
from dataclasses import dataclass, field
from typing import Any, List, Optional

__all__ = ['TraverseClient', 'TraverseResponse', 'TraverseError', 'ServerInitializing',
           'AttemptAtBound', 'AttemptExpired', 'ACCEPT_ENCODING']

ACCEPT_ENCODING = 'gzip, deflate'


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
    `not_after_ms` and the server's `server_time_ms`. Told apart from the other 409 (an
    attempt_id that is running, retained or tombstoned: `attempt` in the body) by its state."""


@dataclass
class TraverseResponse:
    envelope: dict                      # the response without results
    graphlets: List[Any]                # per seed: a Graphlet, or None for an error result
    errors: List[dict] = field(default_factory=list)   # [{index, error, result}]
    raw: Optional[dict] = None
    # deepen(): what the continuation request could not carry exactly (NextRequest.notes:
    # a loss budget conservative for some labels, a switch target left out, ...)
    notes: List[str] = field(default_factory=list)

    @property
    def usage(self):
        """The response-level usage of a request with attempt_id (what the attempt
        consumed: work units, modelled memory, seeds, elapsed ms, the bound the server
        enforced); None without attempt_id."""
        return self.envelope.get('usage')


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
                 scheme='http', graph=None):
        self.host = host
        self.port = port
        self.base = '%s://%s:%s' % (scheme, host, port)
        if api_path:
            self.base += '/' + api_path.strip('/')
        self.session = session
        self.timeout = timeout
        self.release = release
        self.graph = graph

    # ---------------------------------------------------------------- transport

    def _request(self, method, path, payload=None, answers=()):
        """|answers|: non-2xx statuses whose JSON body is an answer, returned rather than
        raised (cancel and attempt: 404 states that nothing runs under the id)."""
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
                data = _decode_body(data, hdrs.get('content-encoding'))
        else:
            req = urllib.request.Request(url, data=body, headers=headers, method=method)
            try:
                with urllib.request.urlopen(req, timeout=self.timeout) as r:
                    status = r.status
                    hdrs = {k.lower(): v for k, v in r.headers.items()}
                    data = _decode_body(r.read(), hdrs.get('content-encoding'))
            except urllib.error.HTTPError as e:
                status = e.code
                hdrs = {k.lower(): v for k, v in (e.headers or {}).items()}
                try:
                    data = _decode_body(e.read() or b'', hdrs.get('content-encoding'))
                finally:
                    e.close()
        text = data.decode('utf-8') if isinstance(data, (bytes, bytearray)) else data
        try:
            out = json.loads(text) if text else {}
        except ValueError:
            out = None
        if status in answers and isinstance(out, dict):
            return out
        if not 200 <= status < 300:
            message = out.get('error') if isinstance(out, dict) else (text or '')[:500]
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
            if status == 409 and isinstance(out, dict) and out.get('state') == 'expired':
                raise AttemptExpired(status, message, out)
            raise TraverseError(status, message, out)
        if out is None:
            raise TraverseError(status, 'the response is not JSON', text[:500])
        return out

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

    def cancel(self, attempt_id, wait_ms=None):
        """POST /traverse/cancel: the server's answer as a dict — `cancelled: True` and
        `state` stopping (or finished, within |wait_ms|) when the attempt was asked to stop;
        `cancelled: False` with `state` finished or unknown (HTTP 404) when nothing runs
        under the id. An unknown id is tombstoned (`tombstone: True`) for the server's
        retention_s: a request arriving later with it is refused (409), so such a 404 from
        the same `server_instance` means it will not run there within that time. When the
        server holds its maximum of tombstones the id is NOT tombstoned (HTTP 429,
        `tombstone: False`): nothing is promised, and the cancel can be retried."""
        payload = {'attempt_id': attempt_id}
        if wait_ms is not None:
            payload['wait_ms'] = wait_ms
        return self._request('POST', '/traverse/cancel', payload, answers=(404, 429))

    def attempt(self, attempt_id):
        """GET /traverse/attempt/{attempt_id}: the attempt's state (running | stopping |
        finished, the reason, when it stopped, the response written, its usage), or
        `state: unknown` (HTTP 404) once it is no longer retained."""
        return self._request('GET', '/traverse/attempt/' + attempt_id, answers=(404,))

    def build_request(self, seeds, strategy=None, *, detail='graphlet', timing=True,
                      graph=None, attempt_id=None, budget_id=None, locus_id=None,
                      not_after_ms=None):
        norm = []
        for s in seeds:
            norm.append({'sequence': s} if isinstance(s, str) else dict(s))
        strategy = copy.deepcopy(strategy) if strategy else {}
        output = strategy.setdefault('output', {})
        output['detail'] = detail
        output['timing'] = timing
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
        return req

    def traverse(self, seeds, strategy=None, *, detail='graphlet', timing=True, graph=None,
                 attempt_id=None, budget_id=None, locus_id=None, not_after_ms=None):
        """-> TraverseResponse(envelope, graphlets, errors). A seed whose permitted set
        could not be derived is an error result ({seed, error}), not an exception: the
        other seeds' graphlets are still returned (so is a seed an attempt never started,
        cancelled or past its bound: resource_stop phase not_started). A request error is
        raised (TraverseError; with attempt_id its usage is in the error's .usage), an
        attempt whose response was stopped at its bound raises AttemptAtBound, and a request
        whose |not_after_ms| (Unix epoch ms) had passed on the server's clock raises
        AttemptExpired (not started)."""
        req = self.build_request(seeds, strategy, detail=detail, timing=timing, graph=graph,
                                 attempt_id=attempt_id, budget_id=budget_id, locus_id=locus_id,
                                 not_after_ms=not_after_ms)
        return self.response(self.traverse_raw(req))

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
