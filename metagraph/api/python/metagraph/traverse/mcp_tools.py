"""Framework-agnostic MCP tool functions over graphlets (DESIGN-traverse-graphlet.md §6).

GraphletTools holds a GraphletStore and the TraverseClients of the indexes it may call;
each public method is one tool and returns a JSON-serialisable dict. The MCP server
(wherever it lives) registers them. Contract:

  * label ids never leave the server: labels appear as {name, ref} everywhere; label
    arguments are tagged selectors ({"ref": ...} | {"name": ...}) or a bare string that
    must be unique as a ref or a name;
  * every local answer carries `evidence: {complete_to_bp, exact, support, reconverge,
    scope, outcome, limitations}` -- the outcome dimensions and the kinds of the stated
    limitations travel with every answer, so none reads more certain than the retrieval;
    exact holds only when the label evidence is complete and no recorded list was cut;
  * list tools page by `cursor` and return `next_cursor` and `total`; a cursor is
    opaque and bound to the entry handle, the tool and its normalized arguments (one
    presented with other arguments is rejected), and stays valid across a restart (its
    secret lives in the spool). EVERY return is <= max(max_bytes, MIN_MAX_BYTES) (2048
    by default; graphlet_sequence and the request of traverse_continue(execute=False)
    have the 16 KB sequence ceiling, stated in their result, and traverse_capabilities a
    16 KB default ceiling, CAPABILITIES_MAX_BYTES, for the server's description of itself
    (an explicit max_bytes still holds); a max_bytes below
    MIN_MAX_BYTES = 64, the smallest error that fits, is a bad_argument): the cursor's
    room is measured, not guessed, and a result that cannot be made to fit is an
    explicit result_too_large whose hint names what this tool offers (raise max_bytes,
    or export). An error keeps its code under the ceiling (its message is cut), also
    when the call itself cannot be bound (a valid max_bytes it supplied still holds). A
    single row larger than max_bytes is returned alone with row_truncated: true and the
    fields that were cut named -- never shortened silently;
  * an operation that makes something (traverse_fetch and traverse_continue storing a
    handle, graphlet_subtrie, graphlet_export, graphlet_save, graphlet_load) always
    returns its receipt -- the new handle or the file's path: over the ceiling its
    optional fields (summary, evidence, sizes) are dropped first and named in
    fields_cut; a receipt that could not fit even then is refused BEFORE anything is
    stored or written (receipt_too_large), never answered with result_too_large after;
  * max_bytes bounds the bytes a tool returns, not the work behind them: the local tools
    (graphlet_compare, graphlet_export, the route listings) have no work or allocation
    budget yet, and their time and peak memory follow the graphlet's size;
  * a filter never hides rows silently: what a filter removed is counted (`filtered`);
  * a derived handle (graphlet_subtrie) is a view: every local tool answers for the
    selected labels and says so (`evidence.view`, `qualified`), or refuses explicitly
    (view_unsupported) where it cannot;
  * oversize: a body above max_graphlet_mb is spooled complete, not kept resident, and
    only the summary and the handle are returned (delivery "spooled", the same value in
    the summary's outcome); a hard store limit rejects the fetch explicitly; a truncated
    document is never stored or returned;
  * files: exports and saves go to the export directory, loads come from it or from the
    spool; a bare file name is placed there, any path resolving elsewhere is refused;
  * failures are results: backend errors (with the HTTP status, and retry_after_s while
    the server initialises), an unreachable backend, I/O errors and malformed arguments
    come back as {error, message, ...}; only bugs raise;
  * an unknown handle answers "replay with traverse_fetch(replay=<handle>)" when the
    store kept its request;
  * a replay or a continuation is derived from its entry only on the entry's index
    (§3.1): the backend's capabilities are checked before the request is sent and the
    returned graphlet before anything is stored -- a different release, namespace,
    index_meta_fp or index_fp is release_mismatch / index_mismatch (proven different); a
    manifest digest on one side only is index_unverifiable (not proven either way), run
    only when the caller passes allow_unverified_index=true -- never as an automatic
    fallback -- and then stated as `identity: {verified: false, accepted, note}`; with no
    digest on either side the result states `identity: {verified: false, note}` instead
    of refusing;
  * a request that would name labels by names the library cannot verify to resolve
    back to them (two labels of the retrieval share one, or it holds U+FFFD) is not
    built: unverifiable_label_name.
"""

import base64
import collections
import hashlib
import hmac
import inspect
import json
import os

from . import derive, export, ops
from ._codec import GraphletFormatError
from .client import AttemptAtBound, ServerInitializing, TraverseError
from .model import (
    ARM_SIDES, AmbiguousLabel, IncompleteRecording, MissingEnvelope, UnknownLabel,
    UnverifiableLabelName,
)
from .parser import utf8_bytes
from .store import StoreLimitExceeded, UnknownHandle

__all__ = ['GraphletTools', 'ToolError', 'tool_names', 'DEFAULT_MAX_BYTES',
           'SEQUENCE_MAX_BYTES', 'CAPABILITIES_MAX_BYTES', 'MIN_MAX_BYTES']

DEFAULT_MAX_BYTES = 2048
SEQUENCE_MAX_BYTES = 16 * 1024
# traverse_capabilities' default ceiling: the server's capabilities in one piece
CAPABILITIES_MAX_BYTES = 16 * 1024
# the smallest max_bytes a tool accepts: below it not even the shortest error fits
# ({"error":"result_too_large","max_bytes":63,"message":""} is 57 bytes)
MIN_MAX_BYTES = 64
# ranked walk lists kept for paging (graphlet_walks): (handle, arm, arguments) -> ids
_RANKED_CACHE = 16
_LABELS_IN_ROW = 8
_WALKS_HINT = {
    'merge_entered': 'walks whose selected label is alive at the end only along its own '
                     'route through a merge (it joined the displayed walk through a '
                     'non-first merge parent, so its displayed support starts there)',
    'not_from_seed': 'walks whose selected label is recorded at the last node but not on '
                     'every node from the seed boundary along any route',
}
_MERGE_ENTERED_HINT = ('claims whose label joined the displayed walk through a non-first '
                       'merge parent (route support from from_bp, displayed support only '
                       'from evidence_from) were left out: pass route_consistent=false to '
                       'include them')


def _size(obj):
    return len(json.dumps(obj, separators=(',', ':'), ensure_ascii=False).encode('utf-8'))


def _labels_json(labels, limit=None):
    out = [l.as_dict() for l in (labels if limit is None else labels[:limit])]
    return out


class ToolError(Exception):
    def __init__(self, code, message, **extra):
        super().__init__(message)
        self.code = code
        self.message = message
        self.extra = extra

    def as_dict(self):
        return dict({'error': self.code, 'message': self.message}, **self.extra)


class _Receipt(dict):
    """The result of an operation that has changed something (a file written, a handle
    created): the caller must always learn its essential fields -- the file's path or the
    new handle -- or it cannot find what was made. Over the ceiling, the optional fields
    are dropped in |optional| order and named in fields_cut; the essential fields with
    every optional one cut are checked to fit BEFORE the operation runs (_check_receipt),
    so a completed operation is never answered with result_too_large."""

    def __init__(self, fields, optional):
        super().__init__(fields)
        self.declared = list(optional)
        self.optional = [k for k in optional if k in fields]

    def minimal(self):
        """The essential fields with every declared optional field named as cut: the
        largest receipt the operation can be reduced to (checked before it runs, when
        the optional fields do not exist yet)."""
        out = {k: v for k, v in self.items() if k not in self.declared}
        if self.declared:
            out['fields_cut'] = list(self.declared)
        return out

    def fit(self, limit):
        out = dict(self)
        cut = []
        for k in self.optional:
            if _size(dict(out, fields_cut=cut) if cut else out) <= limit:
                break
            out.pop(k)
            cut.append(k)
        if cut:
            out['fields_cut'] = cut
        return out


def _check_receipt(receipt, limit, lever='raise max_bytes'):
    """Refuse an operation whose receipt could not be returned under the ceiling even with
    every optional field cut -- before it runs, so that nothing is made that the caller
    would not be told about. |lever| is what THIS tool offers (a file name only where the
    receipt carries one)."""
    n = _size(receipt.minimal())
    if n > limit:
        raise ToolError('receipt_too_large', 'the receipt of this operation (%d bytes: %s) '
                        'would not fit the %d-byte ceiling, so nothing was stored or written: '
                        '%s' % (n, ', '.join(k for k in receipt if k not in receipt.declared),
                                limit, lever), max_bytes=limit)


# ------------------------------------------------------------------ arguments

def _int(name, v, minimum=0):
    if isinstance(v, bool) or not isinstance(v, int) or v < minimum:
        raise ToolError('bad_argument', '%s is an integer >= %d, not %r' % (name, minimum, v))
    return v


def _opt_int(name, v, minimum=0):
    return None if v is None else _int(name, v, minimum)


def _bool(name, v):
    if not isinstance(v, bool):
        raise ToolError('bad_argument', '%s is true or false, not %r' % (name, v))
    return v


def _choice(name, v, choices):
    if v not in choices:
        raise ToolError('bad_argument', '%s is one of %s, not %r'
                        % (name, ', '.join(repr(c) for c in choices), v))
    return v


def _handle(h, name='handle'):
    if not isinstance(h, str):
        raise ToolError('bad_argument', '%s is a handle string (g_...), not %r' % (name, h))
    return h


def _selector(v):
    """A tagged selector ({"ref"|"name"|"id": value}) or a bare string."""
    if isinstance(v, str):
        return v
    if isinstance(v, dict) and len(v) == 1 and next(iter(v)) in ('ref', 'name', 'id'):
        return v
    raise ToolError('bad_argument', 'a label is {"ref": ...}, {"name": ...} or a string, '
                    'not %r' % (v,))


def _selectors(v):
    if v is None:
        return None
    if isinstance(v, (str, dict)):
        return [_selector(v)]
    if isinstance(v, list):
        return [_selector(x) for x in v]
    raise ToolError('bad_argument', 'labels is a label or a list of labels, not %r' % (v,))


# ------------------------------------------------------------------ the wrapper

def _too_large(n, limit, takes_max):
    """The result_too_large that replaces an answer over the ceiling, itself held to it:
    the longest of three forms that fits. The hint names what THIS tool offers."""
    if takes_max:
        hint = 'raise max_bytes, or export the complete answer to a file'
    else:
        hint = 'this tool takes no max_bytes: export the complete answer to a file'
    for out in ({'error': 'result_too_large', 'bytes': n, 'max_bytes': limit,
                 'message': 'the answer does not fit max_bytes even with its rows cut: '
                            + hint},
                {'error': 'result_too_large', 'bytes': n, 'max_bytes': limit,
                 'message': hint},
                {'error': 'result_too_large', 'max_bytes': limit, 'message': ''}):
        if _size(out) <= limit:
            return out
    return out


def _fit_error(out, limit):
    """An error over the ceiling keeps its code: its message is cut (with '...'), then
    its other fields dropped. None when not even the code fits."""
    msg = str(out.get('message', ''))
    for cand in (dict(out), {'error': out['error'], 'message': msg}):
        cand['message'] = msg
        if _size(cand) <= limit:
            return cand
        if _size(dict(cand, message='...')) > limit:
            continue
        lo, hi = 0, len(msg)
        while lo < hi:                       # the longest prefix that fits with '...'
            mid = (lo + hi + 1) // 2
            if _size(dict(cand, message=msg[:mid] + '...')) <= limit:
                lo = mid
            else:
                hi = mid - 1
        return dict(cand, message=msg[:lo] + '...')
    return None


def _valid_ceiling(v):
    return v is not None and not isinstance(v, bool) and isinstance(v, int) \
        and v >= MIN_MAX_BYTES


def _tool(fn):
    """Expected failures become {error, message, ...} results; bugs still raise. Every
    result is held to the tool's ceiling: one that cannot be made to fit is an explicit
    result_too_large, never a silently shortened answer; an error keeps its code."""
    name = fn.__name__
    sig = inspect.signature(fn)
    takes_max = 'max_bytes' in sig.parameters
    # the position of max_bytes among the positional arguments after self
    max_pos = list(sig.parameters).index('max_bytes') - 1 if takes_max else None

    def default_limit(self, execute):
        if name == 'graphlet_sequence' or (name == 'traverse_continue' and execute is False):
            return self.sequence_max_bytes
        if name == 'traverse_capabilities':
            # the server's description of itself, which grows with what it offers (2.4 KB once
            # the attempts were described): an agent's first discovery call must not fail by
            # default (review of the stage-4 backend: result_too_large at 2048)
            return max(self.max_bytes, CAPABILITIES_MAX_BYTES)
        return self.max_bytes

    def wrapped(self, *args, **kw):
        try:
            bound = sig.bind(self, *args, **kw)
        except TypeError as e:
            # an argument the tool does not take is the caller's error, not a bug. The
            # ceiling the caller supplied still holds for this error (the review's
            # graphlet_list(max_bytes=64, foo=1) answered 77 bytes under the default):
            # a valid max_bytes, by keyword or by position, is honoured even though the
            # call as a whole cannot be bound
            supplied = kw.get('max_bytes')
            if supplied is None and max_pos is not None and len(args) > max_pos:
                supplied = args[max_pos]
            limit = supplied if takes_max and _valid_ceiling(supplied) else \
                default_limit(self, kw.get('execute', True))
            out = {'error': 'bad_argument', 'message': str(e)}
            return _fit_error(out, limit) or _too_large(_size(out), limit, takes_max)
        limit = bound.arguments.get('max_bytes')
        if limit is not None and not _valid_ceiling(limit):
            return {'error': 'bad_argument', 'message': 'max_bytes: an int >= %d'
                    % MIN_MAX_BYTES}
        if limit is None:
            limit = default_limit(self, bound.arguments.get('execute', True))
        try:
            out = fn(self, *args, **kw)
        except ToolError as e:
            out = e.as_dict()
        except UnknownHandle as e:
            out = {'error': 'unknown_handle', 'message': str(e), 'replayable': e.replayable}
            if e.replayable:
                out['hint'] = 'replay with traverse_fetch(replay=%r)' % e.handle
        except UnknownLabel as e:
            # before AmbiguousLabel, its base class: nothing matched, so a ref would not help
            out = {'error': 'unknown_label', 'message': str(e),
                   'hint': 'list the labels with graphlet_labels(handle) and select one by '
                           'its {"ref": ...}'}
        except AmbiguousLabel as e:
            out = {'error': 'ambiguous_label', 'message': str(e),
                   'hint': 'select the label by {"ref": ...}'}
        except MissingEnvelope as e:
            out = {'error': 'missing_envelope', 'message': str(e)}
        except IncompleteRecording as e:
            out = {'error': 'incomplete_recording', 'message': str(e)}
        except UnverifiableLabelName as e:
            # before ValueError, its base class: never resolved by a name it cannot verify
            out = {'error': 'unverifiable_label_name', 'message': str(e),
                   'hint': 'build the request yourself with names verified against the '
                           'index (e.g. traverse_resolve)'}
        except StoreLimitExceeded as e:
            out = {'error': 'too_large', 'message': str(e)}
        except GraphletFormatError as e:
            # never stored, never returned as a graphlet (§6)
            out = {'error': 'format_error', 'message': str(e)}
        except OSError as e:
            # the backend's transport errors are mapped where the backend is called;
            # what is left is the local file system
            out = {'error': 'io_error', 'message': e.strerror or type(e).__name__}
        except (IndexError, ValueError) as e:
            out = {'error': 'bad_argument', 'message': str(e)}
        # the call may have grown the caches of the models it read: the store charges
        # them now, so that max_ram_mb holds between calls
        refresh = getattr(self.store, 'refresh', None)
        if refresh is not None:
            refresh()
        if isinstance(out, _Receipt):
            # a completed operation (a file written, a handle created) keeps its
            # receipt: its optional fields go first, and the essential ones were checked
            # to fit before the operation ran (D6)
            out = out.fit(limit)
        n = _size(out)
        if n > limit:
            fitted = None
            if isinstance(out.get('error'), str) and out['error'] != 'result_too_large':
                fitted = _fit_error(out, limit)
            out = fitted or _too_large(n, limit, takes_max)
        return out
    wrapped.__name__ = fn.__name__
    wrapped.__doc__ = fn.__doc__
    wrapped.__signature__ = sig
    return wrapped


class _LazyRows:
    """A list of rows built on first access (graphlet_walks): _page reads len() and the
    rows of one page, so only those are materialized."""

    def __init__(self, n, build):
        self.n = n
        self.build = build
        self.got = {}

    def __len__(self):
        return self.n

    def __getitem__(self, i):
        if i < 0 or i >= self.n:
            raise IndexError(i)
        row = self.got.get(i)
        if row is None:
            row = self.got[i] = self.build(i)
        return row


class GraphletTools:
    def __init__(self, store, clients=None, *, default_index=None,
                 max_bytes=DEFAULT_MAX_BYTES, sequence_max_bytes=SEQUENCE_MAX_BYTES,
                 secret=None, export_dir=None):
        self.store = store
        self.clients = dict(clients or {})
        self.default_index = default_index or (next(iter(self.clients)) if self.clients
                                               else None)
        self.max_bytes = max_bytes
        self.sequence_max_bytes = sequence_max_bytes
        # cursors outlive a restart for as long as the handles they page do: the secret
        # lives in the spool unless one is configured
        self._key = secret if secret is not None else store.cursor_secret()
        self.export_dir = export_dir
        # (handle, arm, rank, label, min_bp, route_consistent) -> (ranked path ids,
        # filtered): a page of graphlet_walks ranks on cheap keys once and spells only
        # its own rows (path ids, not model objects: a re-parsed graphlet stays valid)
        self._ranked = collections.OrderedDict()

    # ================================================================ plumbing

    def _client(self, index):
        index = index or self.default_index
        if not isinstance(index, str) or index not in self.clients:
            raise ToolError('unknown_index', 'no index %r is configured (known: %s)'
                            % (index, ', '.join(sorted(self.clients)) or 'none'))
        return index, self.clients[index]

    @staticmethod
    def _backend(call, *args):
        """A backend call with its failures as results: an HTTP error with its status
        (503: with retry_after_s), a transport failure as backend_unreachable."""
        try:
            return call(*args)
        except ServerInitializing as e:
            raise ToolError('backend_error', e.message or 'the server is initialising',
                            status=e.status, retry_after_s=e.retry_after,
                            hint='the index is still loading: retry later') from None
        except AttemptAtBound as e:
            # registered and run: its id is used up, so it is no loading server to retry
            raise ToolError('backend_error', str(e.message), status=e.status,
                            hint='the attempt reached the duration bound the server enforces '
                                 'for it; a retry needs a new attempt_id') from None
        except TraverseError as e:
            raise ToolError('backend_error', str(e.message), status=e.status) from None
        except OSError as e:
            raise ToolError('backend_unreachable', 'the index server could not be reached: '
                            '%s' % (getattr(e, 'reason', None) or e.strerror
                                    or type(e).__name__)) from None

    def _resolve(self, handle):
        """-> (the backing graphlet, the GraphletView of a derived handle or None)."""
        got = self.store.view(_handle(handle))
        if isinstance(got, ops.GraphletView):
            return got.backing, got
        return got, None

    def _graphlet(self, handle):
        return self._resolve(handle)[0]

    def _arm(self, g, arm, view=None):
        # an arm is named by a string: anything else (an object, a list) is the caller's
        # error, answered as one -- never a TypeError from using it as a key (round 3,
        # finding B)
        if arm is not None and not isinstance(arm, str):
            raise ToolError('bad_arm', 'arm is "left" or "right" ("l", "r"), not %r' % (arm,))
        try:
            a = g.arm(arm)
        except (KeyError, ValueError) as e:
            raise ToolError('bad_arm', str(e).strip('"'))
        if view is not None and a.side not in view.segments:
            raise ToolError('not_in_view', 'this derived handle is a view of the %s arm only'
                            % ', '.join(view.segments))
        return a

    @staticmethod
    def _in_view(view, side, walk=None, segment=None):
        if view is None:
            return
        if walk is not None and walk not in view.path_ids.get(side, ()):
            raise ToolError('not_in_view', 'walk %d of the %s arm is not in this view (%d '
                            'walks carry the selected labels)'
                            % (walk, side, len(view.path_ids.get(side, ()))))
        if segment is not None and segment not in view.segments.get(side, ()):
            raise ToolError('not_in_view', 'segment %d of the %s arm is not in this view'
                            % (segment, side))

    @staticmethod
    def _no_view(view, tool, handle):
        if view is not None:
            raise ToolError('view_unsupported', '%s does not answer for a view: call it on '
                            'the backing handle %r' % (tool, view.of), of=view.of)

    def _digest(self, tool, handle, args):
        norm = json.dumps({k: v for k, v in sorted(args.items())
                           if k not in ('cursor', 'max_bytes', 'n')},
                          sort_keys=True, default=str)
        return hashlib.sha256(('%s\0%s\0%s' % (tool, handle, norm)).encode()).hexdigest()[:16]

    def _cursor(self, tool, handle, args, offset):
        body = json.dumps({'h': handle, 't': tool, 'a': self._digest(tool, handle, args),
                           'o': offset}, separators=(',', ':')).encode()
        mac = hmac.new(self._key, body, hashlib.sha256).hexdigest()[:16]
        return base64.urlsafe_b64encode(body).decode().rstrip('=') + '.' + mac

    def _offset(self, tool, handle, args, cursor):
        if cursor is None or cursor == '':
            return 0
        if not isinstance(cursor, str):
            raise ToolError('bad_cursor', 'a cursor is the next_cursor string of a page')
        try:
            b64, mac = cursor.rsplit('.', 1)
            body = base64.urlsafe_b64decode(b64 + '=' * (-len(b64) % 4))
            ok = hmac.compare_digest(mac, hmac.new(self._key, body, hashlib.sha256)
                                     .hexdigest()[:16])
            j = json.loads(body)
        except (ValueError, TypeError):
            raise ToolError('bad_cursor', 'not a cursor of this server')
        if not ok or not isinstance(j, dict) or not isinstance(j.get('o'), int):
            raise ToolError('bad_cursor', 'not a cursor of this server')
        if j['h'] != handle or j['t'] != tool or j['a'] != self._digest(tool, handle, args):
            raise ToolError('bad_cursor', 'this cursor belongs to another handle, tool or '
                            'other arguments: repeat the original arguments with it')
        return j['o']

    def _page(self, tool, handle, args, rows, base, max_bytes=None, n=None):
        """Rows from the cursor's offset while the whole result, next_cursor included,
        stays <= max_bytes. The room for next_cursor is the length of the longest cursor
        this list can carry (the last offset), measured, not a fixed guess."""
        max_bytes = _opt_int('max_bytes', max_bytes, 1) or self.max_bytes
        start = self._offset(tool, handle, args, args.get('cursor'))
        if start > len(rows):
            raise ToolError('bad_cursor', 'the cursor is past the end of the list')
        out = dict(base, total=len(rows), rows=[])
        reserve = (len(',"next_cursor":')
                   + len(json.dumps(self._cursor(tool, handle, args, len(rows)))))
        i = start
        while i < len(rows) and (n is None or len(out['rows']) < n):
            out['rows'].append(rows[i])
            more = i + 1 < len(rows)
            if _size(out) + (reserve if more else 0) > max_bytes:
                out['rows'].pop()
                break
            i += 1
        if not out['rows'] and i < len(rows):
            # a single row larger than the page: alone, cut, and the cut named. The
            # frame counts every field that could be named as cut, so the final list
            # (a subset) cannot push the result over
            row = rows[i]
            more = i + 1 < len(rows)
            fields = [k for k, v in row.items() if isinstance(v, (str, list)) and v]
            frame = dict(out, rows=[], row_truncated=True, cut_fields=fields)
            budget = max_bytes - _size(frame) - (reserve if more else 0)
            row, cut = _shrink(row, budget)
            out = dict(out, rows=[row], row_truncated=True, cut_fields=cut)
            i += 1
        if i < len(rows):
            out['next_cursor'] = self._cursor(tool, handle, args, i)
        return out

    def _fit_summary(self, g, frame, side=None, max_bytes=None):
        """g.summary() sized so that |frame| with the summary in it stays <= max_bytes."""
        room = (max_bytes or self.max_bytes) - _size(dict(frame, summary={}))
        return g.summary(side, max_bytes=max(room, 0))

    # ---------------------------------------------------------------- files

    def _export_root(self):
        root = self.export_dir or os.path.join(self.store.spool_dir, 'exports')
        os.makedirs(root, exist_ok=True)
        return os.path.realpath(root)

    def _confined(self, path, roots, default=None):
        """An agent-supplied path, resolved (symlinks and '..' included), that must lie
        under one of |roots|; a relative path is taken in the first root. The tools are
        an agent's -- and any prompt injection's -- reach into the file system: nothing
        outside the export directory (and, for reading, the spool) is written or read."""
        if path is None:
            if default is None:
                raise ToolError('bad_argument', 'a path (a file name) is required')
            path = default
        if not isinstance(path, str) or not path or '\0' in path:
            raise ToolError('bad_argument', 'a path is a file name, not %r' % (path,))
        full = path if os.path.isabs(path) else os.path.join(roots[0], path)
        real = os.path.realpath(full)
        for r in roots:
            if real != r and os.path.commonpath([real, r]) == r:
                return real
        raise ToolError('path_not_allowed', 'files are read and written in the server\'s '
                        'export directory only: pass a file name')

    # ================================================================ backend tools

    @_tool
    def traverse_capabilities(self, index=None, max_bytes=None):
        name, client = self._client(index)
        return {'index': name, 'capabilities': self._backend(client.capabilities)}

    @_tool
    def traverse_resolve(self, index=None, sequence=None, labels=None, discover=None,
                         select=None, limit=10, max_bytes=None):
        name, client = self._client(index)
        if not isinstance(sequence, str) or not sequence:
            raise ToolError('bad_argument', 'sequence is a non-empty string')
        if labels is not None and not (isinstance(labels, list)
                                       and all(isinstance(l, str) for l in labels)):
            raise ToolError('bad_argument', 'labels is a list of label names, not %r' % (labels,))
        limit = _int('limit', limit, 1)
        out = self._backend(lambda: client.resolve(sequence, labels=labels,
                                                   discover=discover, select=select))
        labels_ = out.get('labels', [])
        res = {'index': name, 'k': out.get('k'), 'num_kmers': out.get('num_kmers'),
               'labels_total': len(labels_),
               'labels': [{'name': l['label'], 'kind': l.get('kind'),
                           'kmers_supported': l.get('kmers_supported')}
                          for l in labels_[:limit]],
               'candidates': len(out.get('candidates', []))}
        if 'selection' in out:
            res['seeds'] = [{'seed_id': s['seed_id'], 'sequence_bp': len(s['sequence']),
                             'labels': len(s['labels'])}
                            for s in out['selection']['seeds'][:limit]]
        while _size(res) > (max_bytes or self.max_bytes) and (res['labels']
                                                              or res.get('seeds')):
            # long names: fewer rows, and the cut stated (labels_total stays)
            res['labels_cut_to_fit'] = True
            if res.get('seeds'):
                res['seeds'].pop()
            else:
                res['labels'].pop()
        return res

    @_tool
    def traverse_fetch(self, index=None, seed=None, strategy=None, keep=True,
                       max_graphlet_mb=8, replay=None, allow_unverified_index=False,
                       max_bytes=None):
        """One retrieval -> its summary and (keep) a handle. replay=<handle> re-runs the
        stored request of an entry (also an expired one) against the same index; an index
        whose identity with the entry's cannot be verified (a manifest digest on one side
        only) is refused as index_unverifiable unless allow_unverified_index=true, and
        the result then states it (identity.verified false). A stored result's handle is
        always returned: over the ceiling its summary, evidence and graphlet_bytes are
        cut first (fields_cut)."""
        _bool('keep', keep)
        _bool('allow_unverified_index', allow_unverified_index)
        limit = max_bytes or self.max_bytes
        if isinstance(max_graphlet_mb, bool) or not isinstance(max_graphlet_mb, (int, float)) \
                or max_graphlet_mb < 0:
            raise ToolError('bad_argument', 'max_graphlet_mb is a number >= 0')
        if replay is not None:
            req, ident = self._replay_request(_handle(replay, 'replay'))
            if not req:
                raise ToolError('not_replayable', 'entry %r holds no request (it was loaded '
                                'from a body-only file, or from one whose seed label names '
                                'cannot be verified to resolve to its labels): fetch the '
                                'seed with traverse_fetch(seed=...)' % replay)
            index = index or ident.get('source')
        else:
            if not isinstance(seed, dict) or not isinstance(seed.get('sequence'), str):
                raise ToolError('bad_request', 'seed: {sequence, labels?, seed_id?} is required')
            if strategy is not None and not isinstance(strategy, dict):
                raise ToolError('bad_request', 'strategy is an object')
            req = None
        name, client = self._client(index)
        identity = None
        if req is None:
            req = client.build_request([seed], strategy or {}, detail='graphlet')
        else:
            identity = self._check_same_index(client, ident, allow_unverified_index)
        response = self._backend(client.traverse_raw, req)
        results = response.get('results') or []
        if not results:
            raise ToolError('empty_response', 'the server returned no result')
        result = results[0]
        if 'graphlet' not in result:
            return {'error': 'seed_failed', 'message': result.get('error'),
                    'outcome': result.get('outcome'), 'limitations': result.get('limitations')}
        nbytes = result.get('graphlet_bytes')
        if nbytes is None:
            nbytes = len(utf8_bytes(result['graphlet']))
        spooled = nbytes > max_graphlet_mb * (1 << 20)
        from .parser import from_response
        g = from_response(result, response)       # the whole body validated
        if identity is not None:
            # what actually answered, not only what the capabilities said before
            identity = _weakest(identity, _same_index(ident, self._identity_of(g),
                                                      'the replay', allow_unverified_index))
        # the delivery the receipt will state: a spooled body's is 'spooled', one byte longer
        # than 'inline', so the receipt is checked with it (the stage-2 recheck, P3: checked
        # with 'inline', a receipt one byte over the ceiling stored its entry and returned
        # result_too_large without the handle)
        delivery = 'spooled' if spooled else g.outcome.delivery
        if spooled or keep:
            # checked before the entry is stored: a handle the caller cannot be told
            # about is not made (a spooled body is stored even with keep false)
            receipt = {'index': name, 'delivery': delivery,
                       'handle': 'g_' + '0' * 12}
            if identity is not None:
                receipt['identity'] = identity
            _check_receipt(_Receipt(receipt, ('summary', 'evidence', 'graphlet_bytes')),
                           limit)
        if spooled:
            # max_graphlet_mb is a RAM threshold: the body goes to the spool complete, the
            # model built to validate it is not kept, and the outcome says what happened
            handle = self.store.put_parsed(g, result['graphlet'], req, source=name,
                                           resident=False, delivery='spooled')
        elif keep:
            handle = self.store.put_parsed(g, result['graphlet'], req, source=name)
        else:
            handle = None
        out = {'index': name, 'delivery': delivery, 'graphlet_bytes': nbytes}
        if handle is not None:
            out['handle'] = handle
        if identity is not None:
            out['identity'] = identity
        out['evidence'] = ops.evidence_block(g)
        out['summary'] = self._fit_summary(g, out, max_bytes=max_bytes)
        if handle is None:
            return out
        return _Receipt(out, ('summary', 'evidence', 'graphlet_bytes'))

    def _replay_request(self, handle):
        """The stored request and index identity of an entry, live or expired."""
        try:
            e = self.store.get(handle)
            return e.request, e.index
        except UnknownHandle as u:
            if u.replayable and u.request:
                return u.request, u.index
            raise

    def _check_same_index(self, client, ident, allow_unverified=False):
        """A replay or a continuation runs against the same index as its entry: the
        backend's capabilities must state the entry's identity (_same_index) before
        anything is submitted. -> the identity statement for the result."""
        caps = self._backend(client.capabilities)
        return _same_index(ident, {'ns': caps.get('index_ns') or None,
                                   'fp': caps.get('index_fp') or None,
                                   'meta_fp': caps.get('index_meta_fp') or None,
                                   'release': caps.get('release') or None},
                           'the index', allow_unverified)

    @staticmethod
    def _identity_of(g):
        return {'ns': g.index_ns, 'fp': g.index_fp, 'meta_fp': g.index_meta_fp,
                'release': (g.envelope or {}).get('release') or None}

    @_tool
    def traverse_continue(self, handle, arm, walk, overrides=None, execute=True,
                          allow_unverified_index=False, max_bytes=None, reset_branches=False):
        """A continuation is a new traversal: the library builds the request from the
        walk's continuation, the backend runs it (execute=False: the request only, under
        the 16 KB sequence ceiling -- it carries the continuation's bases and label
        names -- stated as ceiling_bytes; max_bytes raises it). The new entry keeps its
        parent walk, and is stored only when the parent's index identity holds before the
        request (capabilities) and in the result (_same_index); the result states it
        (`identity`). A manifest digest on one side only is index_unverifiable unless
        allow_unverified_index=true. Both forms carry `notes` (what the request cannot
        carry exactly: one loss budget for labels that ended at different losses is
        conservative for the lower ones, and one branch allowance likewise, a switch
        target left out, a label alive at the leaf that is not seeded, an allowance that
        restarts with reset_branches=true), where a loss budget applies `loss_budget`
        ({original, effective, largest_terminal_loss}), and where a branch allowance
        applies `branch_budget` ({original, effective, largest_terminal_branches, reset}:
        branching.max_label_branches reduced by the largest terminal branch count of the
        continued labels, "unlimited" kept; reset_branches=true keeps the original): all
        are part of the receipt, never cut. `loss_budget_labels` (each label's terminal
        loss and remaining budget) is the first optional field cut, and named in
        fields_cut, when the receipt would not fit."""
        g, view = self._resolve(handle)
        a = self._arm(g, arm, view)
        walk = _int('walk', walk)
        self._in_view(view, a.side, walk=walk)
        if overrides is not None and not isinstance(overrides, dict):
            raise ToolError('bad_argument', 'overrides is an object')
        _bool('execute', execute)
        _bool('allow_unverified_index', allow_unverified_index)
        _bool('reset_branches', reset_branches)
        req = g.next_request(a, [walk], reset_branches=reset_branches, **(overrides or {}))
        cont = ops.continuation(g, a, walk)
        parent = {'handle': handle, 'arm': a.side, 'walk': walk,
                  'overlap_bp': len(cont.sequence)}
        # The notes and the budget's derivation are essential: they state how the
        # continuation differs from one uninterrupted walk. Each label's own loss and
        # remaining budget is detail the notes already summarise, so it is the first
        # field cut (and named) when the receipt would not fit: a continuation of a walk
        # with many labels was refused under the default ceiling (round 3, finding C)
        stated = {'notes': list(req.notes)}
        per_label = None
        if req.loss_budget is not None:
            stated['loss_budget'] = {k: v for k, v in req.loss_budget.items() if k != 'labels'}
            per_label = req.loss_budget.get('labels')
        if req.branch_budget is not None:
            # how the branch allowance was derived, as essential as the loss budget's
            stated['branch_budget'] = dict(req.branch_budget)
        if not execute:
            out = dict({'request': dict(req), 'parent': parent,
                        'ceiling_bytes': max_bytes or self.sequence_max_bytes}, **stated)
            if per_label is not None:
                out['loss_budget_labels'] = per_label
            return out
        ident = self.store.get(handle).index
        name, client = self._client(ident.get('source'))
        # a continuation is derived from its parent only on the parent's index: checked
        # against the capabilities before submitting, and against what answered before
        # the parent/child link is stored
        identity = self._check_same_index(client, ident, allow_unverified_index)
        optional = ('loss_budget_labels', 'summary', 'evidence')
        _check_receipt(_Receipt(dict({'handle': 'g_' + '0' * 12, 'parent': parent,
                                      'identity': identity}, **stated), optional),
                       max_bytes or self.max_bytes)
        response = self._backend(client.traverse_raw, req)
        result = (response.get('results') or [{}])[0]
        if 'graphlet' not in result:
            return {'error': 'seed_failed', 'message': result.get('error'), 'parent': parent}
        from .parser import from_response
        g2 = from_response(result, response)
        identity = _weakest(identity, _same_index(ident, self._identity_of(g2),
                                                  'the continuation', allow_unverified_index))
        # again with what answered: a weaker identity states more
        _check_receipt(_Receipt(dict({'handle': 'g_' + '0' * 12, 'parent': parent,
                                      'identity': identity}, **stated), optional),
                       max_bytes or self.max_bytes)
        new = self.store.put_parsed(g2, result['graphlet'], dict(req), source=name,
                                    parent=parent, derived_from=handle)
        out = dict({'handle': new, 'parent': parent, 'identity': identity}, **stated)
        if per_label is not None:
            out['loss_budget_labels'] = per_label
        out['evidence'] = ops.evidence_block(g2)
        out['summary'] = self._fit_summary(g2, out, max_bytes=max_bytes)
        return _Receipt(out, optional)

    # ================================================================ local tools

    @_tool
    def graphlet_summary(self, handle, arm=None, detail='summary', cursor=None,
                         max_bytes=None):
        """detail='summary': the agent-facing summary (a view's for a derived handle; a
        continuation's parent walk). detail='limitations': every stated limitation in
        full (knob, limit, observed, complete_to_bp, effect) and the resource stop,
        paged like a list."""
        g, view = self._resolve(handle)
        side = None if arm is None else self._arm(g, arm, view).side
        _choice('detail', detail, ('summary', 'limitations'))
        e = self.store.get(handle)
        if detail == 'limitations':
            rows = [ops.limitation_dict(l) for l in g.limitations
                    if side is None or l.arm is None or l.arm == side]
            base = {'handle': handle, 'outcome': g.outcome.as_dict(),
                    'resource_stop': ops.resource_stop_dict(g.resource_stop)}
            args = dict(arm=side, detail=detail, cursor=cursor)
            return self._page('graphlet_summary', handle, args, rows, base, max_bytes)
        out = {'handle': handle, 'evidence': ops.evidence_block(g, side, view)}
        if e.parent is not None:
            out['parent'] = e.parent
        if view is not None:
            s = view.summary()
            s['outcome'] = g.outcome.as_dict()
            if side is not None:
                s['arms'] = {k: v for k, v in s['arms'].items() if k == side}
            out['summary'] = s
            return out
        out['summary'] = self._fit_summary(g, out, side, max_bytes)
        return out

    @_tool
    def graphlet_walks(self, handle, arm, rank='support', n=10, min_bp=0, label=None,
                       spell='none', tail_bp=60, route_consistent=True, cursor=None,
                       max_bytes=None):
        g, view = self._resolve(handle)
        a = self._arm(g, arm, view)
        _choice('rank', rank, ('support', 'length', 'loss', 'id'))
        _choice('spell', spell, ('none', 'tail', 'full'))
        n = _int('n', n, 1)
        min_bp = _int('min_bp', min_bp)
        tail_bp = _int('tail_bp', tail_bp)
        _bool('route_consistent', route_consistent)
        labels = None if label is None else [_selector(label)]
        keep = None if view is None else set(view.path_ids.get(a.side, ()))
        ids, filtered = self._ranked_ids(handle, g, a, rank, labels, min_bp,
                                         route_consistent, keep)

        def row_of(i):
            w, = ops.walks_at(g, a, [ids[i]])
            row = {'walk': w.path_id, 'length_bp': w.length_bp,
                   'end': w.path_reason or ','.join(sorted(w.end_reasons)) or 'semantic',
                   'n_alive': w.n_alive, 'n_full': len(w.labels_full),
                   'labels_full': _labels_json(w.labels_full, _LABELS_IN_ROW),
                   'complete': w.complete, 'beyond_certified_bp': w.beyond_certified_bp}
            if len(w.labels_full) > _LABELS_IN_ROW:
                row['labels_full_truncated'] = True
            if spell != 'none' and w.sequence is not None:
                s = w.sequence
                if spell == 'tail':
                    s = s[-tail_bp:] if a.side == 'right' else s[:tail_bp]
                row['sequence'] = s
            return row
        rows = _LazyRows(len(ids), row_of)
        args = dict(arm=a.side, rank=rank, min_bp=min_bp, label=label, spell=spell,
                    tail_bp=tail_bp, route_consistent=route_consistent, cursor=cursor)
        base = {'handle': handle, 'arm': a.side,
                'evidence': ops.evidence_block(g, a.side, view)}
        if filtered:
            # counted as what they are, as graphlet_claims counts its claims: a walk the
            # filter removes is merge-entered (route_only is a claim's kind at a cut,
            # which no uncut walk is)
            base['filtered'] = dict(sorted(filtered.items()))
            base['hint'] = ' and '.join(_WALKS_HINT[k] for k in sorted(filtered)) + \
                ' were left out: pass route_consistent=false to include them'
        return self._page('graphlet_walks', handle, args, rows, base, max_bytes, n)

    def _ranked_ids(self, handle, g, a, rank, labels, min_bp, route_consistent, keep):
        """(the ranked path ids of graphlet_walks, {why: count} of the walks the
        route_consistent filter removed), cached per handle and arguments: paging an arm
        ranks it once (B3)."""
        key = (handle, a.side, rank, json.dumps(labels, sort_keys=True, default=str),
               min_bp, route_consistent)
        got = self._ranked.get(key)
        if got is not None:
            self._ranked.move_to_end(key)
            return got
        ids = ops.rank_walks(g, a, by=rank, labels=labels, min_bp=min_bp,
                             route_consistent=route_consistent)
        if keep is not None:
            ids = [i for i in ids if i in keep]
        filtered = collections.Counter()
        if labels is not None and route_consistent:
            every = set(ops.rank_walks(g, a, by='id', labels=labels, min_bp=min_bp,
                                       route_consistent=False))
            if keep is not None:
                every &= keep
            for pid in sorted(every - set(ids)):
                filtered[ops.walk_filter_reason(g, a, pid, labels)] += 1
        got = self._ranked[key] = (ids, filtered)
        while len(self._ranked) > _RANKED_CACHE:
            self._ranked.popitem(last=False)
        return got

    @staticmethod
    def _walk_rows(g, a, p, cut=_LABELS_IN_ROW):
        """One row per segment of walk |p|: labels in and out with their TRUE totals (a
        recorded list may be cut: entry_total, the last P total), merges, label ends."""
        rows = []
        for sid in p.segments:
            s = a.segments[sid]
            out_total = (s.presence[-1].total if s.presence else s.entry_total) \
                if g.mode != 'constrain' else len(s.end)
            row = {'segment': sid, 'from_bp': s.from_bp, 'to_bp': s.end_bp,
                   'labels_in': len(s.entry), 'labels_in_total': s.entry_total,
                   'labels_out': len(s.end), 'labels_out_total': out_total,
                   'exact': s.entry_total == len(s.entry) and out_total == len(s.end),
                   'in': _labels_json([g.labels[l] for l in s.entry], cut),
                   'out': _labels_json([g.labels[l] for l in s.end], cut)}
            if len(s.parents) > 1:
                row['merge_of'] = list(s.parents)
            ends = [a.runs[r] for r in derive.label_end_events(a).get(sid, ())]
            if ends:
                row['ends'] = [dict(g.labels[r.label].as_dict(), at_bp=r.to_bp,
                                    why=r.event_reason)
                               for r in (ends if cut is None else ends[:cut])]
            rows.append(row)
        return rows

    @_tool
    def graphlet_walk(self, handle, arm, walk, cursor=None, max_bytes=None):
        """One walk's segment chain (labels in/out per segment with their true totals,
        events) and its continuation; the chain is paged like a list."""
        g, view = self._resolve(handle)
        a = self._arm(g, arm, view)
        pid = ops.path_id(a, _int('walk', walk))
        self._in_view(view, a.side, walk=pid)
        p = derive.paths(a)[pid]
        rows = self._walk_rows(g, a, p)
        base = {'handle': handle, 'arm': a.side, 'walk': pid, 'length_bp': p.length_bp,
                'evidence': ops.evidence_block(g, a.side, view)}
        if a.segments[p.leaf].leaf.continuation is not None:
            c = ops.continuation(g, a, pid)
            base['continuation'] = {'length_bp': len(c.sequence), 'loss_used': c.loss_used,
                                    'labels': len(c.labels)}
        args = dict(arm=a.side, walk=pid, cursor=cursor)
        return self._page('graphlet_walk', handle, args, rows, base, max_bytes)

    @_tool
    def graphlet_support(self, handle, arm, walk, step=None, cursor=None, max_bytes=None):
        """Support runs along a walk with who left / joined and why (an end code, a
        switch, "split: labels took C", a merge); n_labels_total is the true count where
        a recorded list was cut."""
        g, view = self._resolve(handle)
        a = self._arm(g, arm, view)
        walk = _int('walk', walk)
        step = _opt_int('step', step)
        self._in_view(view, a.side, walk=walk)
        prof = g.support_profile(a, walk)
        changes = {c.at_bp: c for c in g.support_changes(a, walk)}
        rows = []
        for r in prof:
            row = {'from_bp': r.from_bp, 'to_bp': r.to_bp, 'n_labels': len(r.labels),
                   'n_labels_total': len(r.labels) if r.total is None else r.total,
                   'exact': r.exact}
            c = changes.get(r.from_bp)
            if c is not None:
                row['joined'] = _labels_json(c.added, _LABELS_IN_ROW)
                row['left'] = _labels_json(c.removed, _LABELS_IN_ROW)
                row['why'] = [{k: v for k, v in x.items() if k != 'change'}
                              for x in c.reasons[:_LABELS_IN_ROW]]
            rows.append(row)
        if step:
            rows = [r for r in rows if r['to_bp'] - r['from_bp'] >= step or 'why' in r]
        args = dict(arm=a.side, walk=walk, step=step, cursor=cursor)
        base = {'handle': handle, 'arm': a.side, 'walk': walk,
                'evidence': ops.evidence_block(g, a.side, view)}
        return self._page('graphlet_support', handle, args, rows, base, max_bytes)

    @_tool
    def graphlet_labels(self, handle, arm=None, rank='direct_bp', n=20, min_direct_bp=0,
                        at_bp=None, name=None, cursor=None, max_bytes=None):
        g, view = self._resolve(handle)
        if arm is not None:
            sides = [self._arm(g, arm, view).side]
        else:
            sides = [s for s in ARM_SIDES if s in g.arms
                     and (view is None or s in view.segments)]
        n = _int('n', n, 1)
        in_view = None if view is None else {l.id for l in view.labels}
        if name is not None:
            lab = ops.label(g, _selector(name))
            if in_view is not None and lab.id not in in_view:
                raise ToolError('not_in_view', 'label %s is not one of this view\'s labels'
                                % lab.ref)
            walks = g.label_walks(lab, sides[0] if len(sides) == 1 else None)
            # each run's OWN route (ops.routes yields one per run of the label, in run
            # order on the arm), not the first run's; annotate mode: one witness route
            # per label walk, in label_walks order
            routes = {}
            for side in {lw.arm for lw in walks}:
                if g.mode == 'constrain':
                    ids = [r.id for r in g.arms[side].runs if r.label == lab.id]
                    routes[side] = dict(zip(ids, ops.routes(g, lab, side)))
                else:
                    routes[side] = iter(ops.routes(g, lab, side))
            rows = []
            for lw in walks:
                if g.mode == 'constrain':
                    route = routes[lw.arm][lw.run] if lw.run is not None else None
                else:
                    route = next(routes[lw.arm])
                rows.append({'arm': lw.arm, 'from_bp': lw.from_bp, 'to_bp': lw.to_bp,
                             'evidence_from': lw.evidence_from, 'end': lw.end,
                             'walks_below': len(lw.leaves_below),
                             'merged_into': lw.merged_into, 'route': route})
            base = dict(lab.as_dict(), handle=handle,
                        evidence=ops.evidence_block(g, arm, view))
            args = dict(arm=arm, name=lab.ref, cursor=cursor)
            return self._page('graphlet_labels', handle, args, rows, base, max_bytes, n)
        _choice('rank', rank, ('direct_bp', 'reach_bp'))
        min_direct_bp = _int('min_direct_bp', min_direct_bp)
        at_bp = _opt_int('at_bp', at_bp)
        alive = None if at_bp is None else {s: _alive_at(g, g.arms[s], at_bp) for s in sides}
        summary = derive.label_summary(g)
        rows = []
        for lab, per in zip(g.labels, summary):
            if in_view is not None and lab.id not in in_view:
                continue
            for side in sides:
                v = per[side]
                if v['direct_bp'] < min_direct_bp:
                    continue
                if alive is not None and lab.id not in alive[side]:
                    continue
                rows.append(dict(lab.as_dict(), arm=side, direct_bp=v['direct_bp'],
                                 reach_bp=v['reach_bp'], reentries=v['reentries'],
                                 runs=len(v['runs'])))
        rows.sort(key=lambda r: (-r[rank], r['ref'], r['arm']))
        args = dict(arm=arm, rank=rank, min_direct_bp=min_direct_bp, at_bp=at_bp,
                    cursor=cursor)
        base = {'handle': handle, 'evidence': ops.evidence_block(g, arm, view)}
        return self._page('graphlet_labels', handle, args, rows, base, max_bytes, n)

    @_tool
    def graphlet_splits(self, handle, arm, n=10, min_labels_before=2, cursor=None,
                        max_bytes=None):
        g, view = self._resolve(handle)
        a = self._arm(g, arm, view)
        n = _int('n', n, 1)
        min_labels_before = _int('min_labels_before', min_labels_before)
        keep = None if view is None else set(view.segments.get(a.side, ()))
        rows = []
        for sp in derive.splits(a, g.mode):
            if sp.labels_before < min_labels_before:
                continue
            if keep is not None and sp.segment not in keep:
                continue
            rows.append({'at_bp': sp.at_bp, 'segment': sp.segment,
                         'kind': 'ambiguous' if sp.ambiguous else 'divergence',
                         'labels_before': sp.labels_before,
                         'branches': [{'char': b['char'], 'labels': b['labels_distinct']}
                                      for b in derive.split_branches(a, sp, g.cap)]})
        args = dict(arm=a.side, min_labels_before=min_labels_before, cursor=cursor)
        base = {'handle': handle, 'arm': a.side,
                'evidence': ops.evidence_block(g, a.side, view)}
        return self._page('graphlet_splits', handle, args, rows, base, max_bytes, n)

    @_tool
    def graphlet_claims(self, handle, arm, min_bp=0, route_consistent=True, labels=None,
                        cursor=None, max_bytes=None):
        g, view = self._resolve(handle)
        a = self._arm(g, arm, view)
        min_bp = _int('min_bp', min_bp)
        _bool('route_consistent', route_consistent)
        sel = _selectors(labels)
        if view is not None:
            ids = {l.id for l in view.labels}
            if sel is not None:
                ids &= {l.id for l in g.labels_matching(sel)}
            sel = [{'id': i} for i in sorted(ids)]
        rows = []
        filtered = collections.Counter()
        for c in g.claims(a, labels=sel, strict=False):
            ev = c.evidence_from if c.evidence_from is not None else c.to_bp
            if c.to_bp - ev < min_bp:
                continue
            if route_consistent and c.route_bp:
                # a merge-entered claim (displayed support from evidence_from) is counted
                # as what it is: route_only is the kind of a claim with NO displayed
                # support left at a cut, which this list (uncut) never holds
                filtered['route_only' if c.kind == 'route_only' else 'merge_entered'] += 1
                continue
            d = c.as_dict()
            d.pop('arm')
            rows.append(d)
        args = dict(arm=a.side, min_bp=min_bp, route_consistent=route_consistent,
                    labels=labels, cursor=cursor)
        base = {'handle': handle, 'arm': a.side,
                'evidence': ops.evidence_block(g, a.side, view)}
        if filtered:
            base['filtered'] = dict(sorted(filtered.items()))
            base['hint'] = _MERGE_ENTERED_HINT
        if g.mode == 'annotate' and a.labels_per_node.nodes_truncated:
            base['caveat'] = ('recorded label lists were cut: these claims are lower bounds '
                              '(label_lists limitation)')
        return self._page('graphlet_claims', handle, args, rows, base, max_bytes)

    @_tool
    def graphlet_sequence(self, handle, arm, walk=None, segment=None, start=0, end=None,
                          orientation='natural', with_seed=False):
        """A slice [start, end) of a walk's (or a segment's) bases, at most 16 KB; the
        result states the ceiling and where to continue."""
        g, view = self._resolve(handle)
        a = self._arm(g, arm, view)
        if (walk is None) == (segment is None):
            raise ToolError('bad_argument', 'name exactly one of walk or segment')
        _choice('orientation', orientation, ('natural', 'walk'))
        _bool('with_seed', with_seed)
        start = _int('start', start)
        end = _opt_int('end', end)
        if end is not None and end < start:
            raise ToolError('bad_argument', 'end %d is before start %d' % (end, start))
        if walk is not None:
            walk = _int('walk', walk)
            self._in_view(view, a.side, walk=ops.path_id(a, walk))
        else:
            segment = _int('segment', segment)
            if segment >= len(a.segments):
                raise ToolError('bad_argument', 'the %s arm has %d segments, no segment %d'
                                % (a.side, len(a.segments), segment))
            self._in_view(view, a.side, segment=segment)
        try:
            if walk is not None:
                seq = ops.spell(g, a, walk, orientation, with_seed)
            else:
                s = a.segments[segment]
                if s.walk is None:
                    raise ValueError('this retrieval carries no bases (output.sequences: false)')
                seq = s.walk if (orientation == 'walk' or a.side == 'right') else s.walk[::-1]
        except ValueError as e:
            raise ToolError('no_bases', str(e))
        if start > len(seq):
            raise ToolError('bad_argument', 'start %d is beyond the %d bases' % (start, len(seq)))
        end = len(seq) if end is None else min(end, len(seq))
        part = seq[start:end]
        out = {'handle': handle, 'arm': a.side, 'orientation': orientation,
               'with_seed': with_seed, 'from': start, 'length': len(seq),
               'ceiling_bytes': self.sequence_max_bytes,
               'evidence': ops.evidence_block(g, a.side, view)}
        out['to'] = start + len(part)
        out['sequence'] = part
        if _size(out) > self.sequence_max_bytes:
            # the ceiling holds for the whole result, the continuation keys included
            frame = _size(dict(out, sequence='', truncated=True, to=end, next_from=end))
            part = part[:max(0, self.sequence_max_bytes - frame)]
            out.update(sequence=part, truncated=True, to=start + len(part),
                       next_from=start + len(part))
        return out

    @_tool
    def graphlet_export(self, handle, format='fasta', arm=None, walks=None, path=None,
                        what=None, other=None):
        """Write a complete export to a file in the export directory -> {path, bytes,
        records}: fasta, gfa, json (the full results[i]; what=compare with other=<handle>
        writes the complete comparison lists; what=walk with arm and walks writes the
        complete segment chains of those walks), mgt (the standalone file). |path| is a
        file name (default: <handle>.<ext>). The receipt is always returned: records,
        then bytes, are cut first (fields_cut). The export has no work or allocation
        budget: the file is as large as the graphlet makes it."""
        g, view = self._resolve(handle)
        _choice('format', format, ('fasta', 'gfa', 'json', 'mgt'))
        if format == 'json':
            _choice('what', what, (None, 'compare', 'walk'))
        elif what is not None:
            raise ToolError('bad_argument', 'what= applies to format json only')
        side = None if arm is None else self._arm(g, arm, view).side
        if walks is not None:
            if isinstance(walks, int) and not isinstance(walks, bool):
                walks = [walks]
            if not isinstance(walks, list):
                raise ToolError('bad_argument', 'walks is a walk id or a list of walk ids')
            walks = [_int('walk', w) for w in walks]
        path = self._confined(path, [self._export_root()],
                              default='%s.%s' % (handle, {'json': 'json', 'mgt': 'mgt',
                                                          'gfa': 'gfa'}.get(format, 'fa')))
        _check_receipt(_Receipt({'path': path}, ('records', 'bytes')), self.max_bytes,
                       'use a shorter file name')
        if format == 'fasta':
            if view is not None:
                parts = []
                for s in view.segments:
                    if side is not None and s != side:
                        continue
                    ids = [w for w in view.path_ids[s] if walks is None or w in walks]
                    parts.append(export.to_fasta(g, s, ids))
                text = ''.join(parts)
            else:
                text = g.to_fasta(side, walks)
            records = text.count('>')
        elif format == 'gfa':
            self._no_view(view, 'graphlet_export(format=gfa)', handle)
            text = g.to_gfa()
            records = sum(1 for l in text.splitlines() if l[:1] in 'SLP')
        elif format == 'json':
            if what == 'compare':
                self._no_view(view, 'graphlet_export(what=compare)', handle)
                if other is None:
                    raise ToolError('bad_argument', 'what=compare needs other=<handle>')
                gb, vb = self._resolve(other)
                self._no_view(vb, 'graphlet_export(what=compare)', other)
                cmp = g.compare(gb, arm=side)
                obj = cmp.as_dict()
                records = len(cmp.only_in_a) + len(cmp.only_in_b) + len(cmp.differ)
            elif what == 'walk':
                if side is None or not walks:
                    raise ToolError('bad_argument', 'what=walk needs arm and walks')
                a = g.arms[side]
                obj = {'handle': handle, 'arm': side, 'walks': []}
                for w in walks:
                    pid = ops.path_id(a, w)
                    self._in_view(view, side, walk=pid)
                    p = derive.paths(a)[pid]
                    try:
                        seq = ops.spell(g, a, pid)
                    except ValueError:
                        seq = None
                    obj['walks'].append({'walk': pid, 'length_bp': p.length_bp,
                                         'sequence': seq,
                                         'segments': self._walk_rows(g, a, p, cut=None)})
                obj['evidence'] = ops.evidence_block(g, side, view)
                records = len(obj['walks'])
            else:
                self._no_view(view, 'graphlet_export(format=json)', handle)
                obj = g.to_json()
                records = sum(len(x['paths']) for x in obj['arms'].values())
            text = json.dumps(obj, ensure_ascii=False)
        else:
            from .parser import dump
            # a derived handle's graphlet carries its view: the file is the unchanged
            # backing body with the view in J (§5)
            text = dump(g, envelope=True)
            records = text.count('\n')
        data = text.encode('utf-8')
        with open(path, 'wb') as f:
            f.write(data)
        return _Receipt({'path': path, 'bytes': len(data), 'records': records},
                        ('records', 'bytes'))

    @_tool
    def graphlet_compare(self, a, b, arm=None, mode='claims', cursor=None, max_bytes=None):
        """Graphlet.compare() of two handles over their DAGs restricted to the common
        certified depth, paged. max_bytes bounds the page returned, not the comparison:
        it runs with no work or allocation budget (its cost follows both DAGs up to the
        depth)."""
        ga, va = self._resolve(a)
        gb, vb = self._resolve(b)
        self._no_view(va, 'graphlet_compare', a)
        self._no_view(vb, 'graphlet_compare', b)
        _choice('mode', mode, ('claims', 'walks', 'labels', 'prefix_subset'))
        if arm is not None:
            self._arm(ga, arm)
        cmp = ga.compare(gb, arm=arm, mode=mode)
        rows = ([dict(side='a', **_row(x)) for x in cmp.only_in_a]
                + [dict(side='b', **_row(x)) for x in cmp.only_in_b]
                + [dict(side='differ', **_row(x)) for x in cmp.differ])
        sides = [ga.arm(arm).side] if arm is not None else \
            [s for s in ARM_SIDES if s in ga.arms and s in gb.arms]
        evidence = {}
        for tag, g in (('a', ga), ('b', gb)):
            ev = ops.evidence_block(g, sides[0] if len(sides) == 1 else None) \
                if sides else ops.evidence_block(g)
            evidence[tag] = ev
        base = {'a': a, 'b': b, 'mode': mode, 'comparable': cmp.comparable,
                'reason': cmp.reason, 'equal': cmp.equal, 'depth_used': cmp.depth_used,
                'counts': {'only_in_a': len(cmp.only_in_a), 'only_in_b': len(cmp.only_in_b),
                           'differ': len(cmp.differ)},
                'notes': cmp.notes, 'support': list(cmp.support),
                'scopes': list(cmp.scopes), 'evidence': evidence,
                'export': 'graphlet_export(format=json, what=compare) writes the complete lists'}
        args = dict(b=b, arm=arm, mode=mode, cursor=cursor)
        return self._page('graphlet_compare', a, args, rows, base, max_bytes)

    @_tool
    def graphlet_subtrie(self, handle, labels, arm=None, mode='any'):
        """A derived handle: a view with the backing graphlet's original ids; every local
        tool on it answers for the selected labels. The new handle is always returned:
        over the ceiling the summary, evidence and view are cut first (fields_cut)."""
        g, view = self._resolve(handle)
        self._no_view(view, 'graphlet_subtrie', handle)
        sel = _selectors(labels)
        if not sel:
            raise ToolError('bad_argument', 'labels: at least one label')
        _choice('mode', mode, ('any', 'all'))
        side = None if arm is None else self._arm(g, arm).side
        v = ops.subgraph(g, sel, side, mode)
        e = self.store.get(handle)
        optional = ('summary', 'evidence', 'view')
        _check_receipt(_Receipt({'handle': 'g_' + '0' * 12, 'of': handle}, optional),
                       self.max_bytes, 'this tool takes no max_bytes: the configured ceiling '
                       'is below the smallest receipt')
        import copy as _copy
        g2 = _copy.copy(g)
        g2.view = v.view_spec()
        g2.view['of'] = handle
        new = self.store.put_graphlet(g2, e.request, derived_from=handle,
                                      source=e.index.get('source'))
        v.of = handle
        return _Receipt({'handle': new, 'of': handle, 'view': g2.view, 'summary': v.summary(),
                         'evidence': ops.evidence_block(g, side, v)}, optional)

    @_tool
    def graphlet_list(self, cursor=None, max_bytes=None):
        rows = self.store.list()
        return self._page('graphlet_list', '', dict(cursor=cursor), rows, {}, max_bytes)

    @_tool
    def graphlet_free(self, handle):
        self.store.free(_handle(handle))
        return {'freed': handle}

    @_tool
    def graphlet_save(self, handle, path):
        """The standalone .mgt file of an entry, written to the export directory."""
        path = self._confined(path, [self._export_root()])
        _check_receipt(_Receipt({'path': path}, ('bytes',)), self.max_bytes,
                       'use a shorter file name')
        n = self.store.save(_handle(handle), path)
        return _Receipt({'path': path, 'bytes': n}, ('bytes',))

    @_tool
    def graphlet_load(self, path):
        """A saved .mgt file (from the export directory or the spool) -> a new handle and
        its summary (cut first, and named in fields_cut, when the answer would not fit:
        the handle is always returned)."""
        roots = [self._export_root(), os.path.realpath(self.store.spool_dir)]
        path = self._confined(path, roots)
        # a handle is 'g_' and 12 hex digits (GraphletStore._new_handle)
        _check_receipt(_Receipt({'handle': 'g_' + '0' * 12}, ('summary',)), self.max_bytes,
                       'this tool takes no max_bytes: the configured ceiling is below the '
                       'smallest receipt')
        try:
            h = self.store.load(path)
        except GraphletFormatError as e:
            # the file may be anything the agent named: say where, never what it holds
            raise ToolError('format_error', 'not a valid MGT v1 document (line %d)'
                            % e.line_no) from None
        except UnicodeDecodeError:
            raise ToolError('format_error', 'not a valid MGT v1 document (not UTF-8)') from None
        g = self.store.graphlet(h)
        return _Receipt({'handle': h, 'summary': self._fit_summary(g, {'handle': h})
                         if g.has_envelope else None}, ('summary',))


def _same_index(parent, other, what, allow_unverified=False):
    """The §3.1 identity check of a retrieval derived from (or replaying) an entry: the
    release and the namespace, where both state one, must agree, a differing
    index_meta_fp proves another index, and the manifest digest index_fp must be equal
    where both state it (index_mismatch otherwise: proven different). A digest on ONE side
    only proves nothing either way: that is index_unverifiable -- not run automatically,
    and run only when the caller opts in with allow_unverified (the result then states the
    identity unverified and that it was accepted). Neither with it (no manifest anywhere)
    is not refused: the result states the identity unverifiable. -> the identity statement
    {verified, ...}."""
    rel_a, rel_b = parent.get('release') or None, other.get('release') or None
    if rel_a and rel_b and rel_a != rel_b:
        raise ToolError('release_mismatch', 'the entry was retrieved from release %r, %s '
                        'serves %r: not run silently' % (rel_a, what, rel_b))
    ns_a, ns_b = parent.get('ns') or None, other.get('ns') or None
    if ns_a and ns_b and ns_a != ns_b:
        raise ToolError('index_mismatch', 'the entry was retrieved from index %r, %s is %r'
                        % (ns_a, what, ns_b))
    meta_a, meta_b = parent.get('meta_fp') or None, other.get('meta_fp') or None
    if meta_a and meta_b and meta_a != meta_b:
        raise ToolError('index_mismatch', 'the entry was retrieved from another index '
                        '(index_meta_fp differs from that of %s)' % what)
    fp_a, fp_b = parent.get('fp') or None, other.get('fp') or None
    if fp_a and fp_b:
        if fp_a != fp_b:
            raise ToolError('index_mismatch', 'the entry was retrieved from another index '
                            '(index_fp differs from that of %s)' % what)
        return {'verified': True, 'index_fp': fp_a}
    if fp_a or fp_b:
        # one manifest digest: equal release, namespace and index_meta_fp (where stated)
        # do not prove the same index, and nothing proves another one
        why = ('%s states %s index manifest digest (index_fp) and the entry %s: whether '
               'they are the same index cannot be verified (not proven different)'
               % (what, 'an' if fp_b else 'no', 'one' if fp_a else 'none'))
        if not allow_unverified:
            raise ToolError('index_unverifiable', why + ', so nothing is derived from the '
                            'entry automatically',
                            hint='pass allow_unverified_index=true to run it anyway; the '
                                 'result then states the identity unverified')
        return {'verified': False, 'accepted': 'allow_unverified_index',
                'note': why + '; run because the caller passed allow_unverified_index'}
    return {'verified': False,
            'note': 'neither the entry nor %s states an index manifest digest (index_fp): '
                    'that they are the same index cannot be verified (the server loaded no '
                    'index manifest)' % what}


def _weakest(first, second):
    """Of two identity statements (before the request, and of what answered), the one
    that claims less: a verified identity never hides an unverified answer."""
    if first.get('verified') and not second.get('verified'):
        return second
    return first


def _row(x):
    return x if isinstance(x, dict) else {'value': x}


def _alive_at(g, arm, at_bp):
    """The labels alive at outward base |at_bp| on any walk of the arm: per segment
    spanning the position, the set at that base (ops.alive_at: run ends and switch-ins
    up to it in constrain mode, the P run holding it in annotate mode) -- not the
    segment's entry or end set."""
    out = set()
    for s in arm.segments:
        if s.from_bp <= at_bp < s.end_bp:
            out |= ops.alive_at(g, arm, s.id, at_bp)
    return out


def _shrink(row, budget):
    """Cut one oversized row to |budget| bytes, the field contributing the most bytes
    first (a string loses its tail, a list half its items); -> (row, the names of the
    cut fields). Scalars are never cut: they are what the row is about."""
    row = dict(row)
    cut = []
    while _size(row) > budget:
        sizes = [(_size(v), k) for k, v in row.items()
                 if (isinstance(v, str) and v) or (isinstance(v, list) and v)]
        if not sizes:
            break
        _, k = max(sizes)
        v = row[k]
        if isinstance(v, str):
            over = _size(row) - budget
            row[k] = v[:max(0, len(v) - over - 1)]
        else:
            row[k] = v[:len(v) // 2]
        if k not in cut:
            cut.append(k)
    return row, cut


def tool_names():
    return [n for n in dir(GraphletTools) if n.startswith(('traverse_', 'graphlet_'))]
