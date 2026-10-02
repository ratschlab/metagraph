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
"""

import copy
import gzip
import json
import urllib.error
import urllib.request
import zlib
from dataclasses import dataclass, field
from typing import Any, List, Optional

__all__ = ['TraverseClient', 'TraverseResponse', 'TraverseError', 'ServerInitializing',
           'ACCEPT_ENCODING']

ACCEPT_ENCODING = 'gzip, deflate'


class TraverseError(RuntimeError):
    """A non-2xx answer: |status| the HTTP status, |message| the server's `error`."""

    def __init__(self, status, message, body=None):
        self.status = status
        self.message = message
        self.body = body
        super().__init__('HTTP %s: %s' % (status, message))


class ServerInitializing(TraverseError):
    """503: the index is still loading; retry after |retry_after| seconds."""

    def __init__(self, status, message, body=None, retry_after=None):
        super().__init__(status, message, body)
        self.retry_after = retry_after


@dataclass
class TraverseResponse:
    envelope: dict                      # the response without results
    graphlets: List[Any]                # per seed: a Graphlet, or None for an error result
    errors: List[dict] = field(default_factory=list)   # [{index, error, result}]
    raw: Optional[dict] = None


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

    def _request(self, method, path, payload=None):
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
        if not 200 <= status < 300:
            message = out.get('error') if isinstance(out, dict) else (text or '')[:500]
            if status == 503:
                retry = hdrs.get('retry-after')
                raise ServerInitializing(status, message, out,
                                         int(retry) if retry and retry.isdigit() else None)
            raise TraverseError(status, message, out)
        if out is None:
            raise TraverseError(status, 'the response is not JSON', text[:500])
        return out

    # ---------------------------------------------------------------- routes

    def capabilities(self):
        return self._request('GET', '/traverse/capabilities')

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

    def build_request(self, seeds, strategy=None, *, detail='graphlet', timing=True,
                      graph=None):
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
        return req

    def traverse(self, seeds, strategy=None, *, detail='graphlet', timing=True, graph=None):
        """-> TraverseResponse(envelope, graphlets, errors). A seed whose permitted set
        could not be derived is an error result ({seed, error}), not an exception: the
        other seeds' graphlets are still returned."""
        req = self.build_request(seeds, strategy, detail=detail, timing=timing, graph=graph)
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
        continuations and builds the request, the server runs it."""
        if self.graph and 'graph' not in overrides:
            overrides['graph'] = self.graph
        req = graphlet.next_request(arm, leaves, bp, **overrides)
        return self.response(self.traverse_raw(req))
