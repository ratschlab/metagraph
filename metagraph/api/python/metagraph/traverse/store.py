"""GraphletStore: handles, an LRU of parsed graphlets in RAM and a disk spool (§5, §6).

ENTRY identity and BODY identity are separate (v2): an entry handle is opaque
('g_' + 12 random hex) and owns (normalized request, index identity, envelope, body
digest); bodies are deduplicated by sha256 in the spool. Two requests with different
loss budgets can produce byte-identical bodies (no switch happened) but continue
differently, so a body hash is never the handle.

Layout under spool_dir:
    bodies/<sha256>.mgt       the body as received (deduplicated)
    entries/<handle>.json     the entry (request, index identity, envelope, summary, digest)
    tombstones/<handle>.json  an expired entry's request, so that it can be replayed
    cursor.key                the MCP layer's cursor secret (0600), so that cursors stay
                              valid across a restart for as long as their handles do

A body is validated (parsed in full) before it is stored: a truncated MGT document is
never stored or returned as a graphlet. RAM holds parsed graphlets under max_ram_mb /
max_handles (least recently used first out) and drops one not used for ttl_ram_s; a
graphlet larger than the whole RAM budget is never kept resident (it is parsed on
demand), and a spooled one is not made resident when it is stored. A resident model is
charged its deep heap footprint (Graphlet.memory_bytes()), and again when the caches its
queries build have grown (refresh(): on every access, and after every MCP tool call),
so that max_ram_mb bounds the heap the resident models hold. The disk keeps an
entry for ttl_disk_s after its last use, and a tombstone for ttl_tomb_s after its
expiry. The clock is injectable.
"""

import collections
import hashlib
import json
import os
import re
import secrets
import tempfile
import time
from dataclasses import dataclass, field
from typing import Any, Optional

from .parser import dump, from_response, j_object, parse, seed_envelope, utf8_bytes

__all__ = ['GraphletStore', 'Entry', 'UnknownHandle', 'StoreLimitExceeded']


# what _new_handle() makes; anything else is unknown without touching the disk (a handle
# is joined into spool paths)
_HANDLE_RE = re.compile(r'g_[0-9a-f]{12}\Z')


class UnknownHandle(KeyError):
    """The handle is not (or no longer) stored. |replayable|: its request is kept (the
    entry expired), so traverse_fetch(replay=<handle>) can run it again."""

    def __init__(self, handle, replayable=False, request=None, index=None):
        self.handle = handle
        self.replayable = replayable
        self.request = request
        self.index = index or {}
        msg = 'unknown graphlet handle %r' % handle
        if replayable:
            msg += ': it expired; replay with traverse_fetch(replay=%r)' % handle
        super().__init__(msg)

    def __str__(self):
        return self.args[0]


class StoreLimitExceeded(ValueError):
    """A body over the store's hard limit: rejected explicitly, never truncated."""


@dataclass
class Entry:
    handle: str
    request: Optional[dict]
    index: dict
    envelope: dict
    seed_summary: dict
    digest: str
    bytes: int
    created: float
    accessed: float
    view: Optional[dict] = None
    derived_from: Optional[str] = None
    # how the body reached this store when it differs from the O record: 'spooled' when
    # the MCP layer spooled it (one value for the outcome everywhere)
    delivery: Optional[str] = None
    # a continuation's parent walk {handle, arm, walk, overlap_bp}: kept with the entry,
    # so the link survives the one-shot return value
    parent: Optional[dict] = None
    store: Any = field(default=None, repr=False, compare=False)

    @property
    def graphlet(self):
        return self.store.graphlet(self.handle)

    def meta(self):
        return {'handle': self.handle, 'index': self.index, 'digest': self.digest,
                'bytes': self.bytes, 'created': self.created, 'accessed': self.accessed,
                'view': self.view, 'derived_from': self.derived_from,
                'parent': self.parent, 'delivery': self.delivery,
                'seed_id': (self.seed_summary.get('seed') or {}).get('seed_id')}


def set_delivery(g, delivery):
    """The outcome's delivery dimension as it actually happened (the model's O record and
    the attached per-seed summary), so that summary(), the stored entry and the tool's
    own `delivery` carry one value."""
    if g.outcome is not None:
        g.outcome.delivery = delivery
    if g.seed_summary and isinstance(g.seed_summary.get('outcome'), dict):
        g.seed_summary = dict(g.seed_summary, outcome=dict(g.seed_summary['outcome'],
                                                            delivery=delivery))


def _write_atomic(path, data):
    d = os.path.dirname(path)
    fd, tmp = tempfile.mkstemp(prefix='.tmp-', dir=d)
    try:
        with os.fdopen(fd, 'wb') as f:
            f.write(data)
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


class GraphletStore:
    def __init__(self, spool_dir, max_ram_mb=512, max_handles=64, ttl_ram_s=1800,
                 ttl_disk_s=7 * 86400, *, clock=time.time, max_body_mb=None,
                 ttl_tomb_s=30 * 86400):
        self.spool_dir = spool_dir
        self.max_ram_bytes = int(max_ram_mb * (1 << 20))
        self.max_handles = max_handles
        self.ttl_ram_s = ttl_ram_s
        self.ttl_disk_s = ttl_disk_s
        self.ttl_tomb_s = ttl_tomb_s
        self.clock = clock
        self.max_body_bytes = None if max_body_mb is None else int(max_body_mb * (1 << 20))
        for sub in ('bodies', 'entries', 'tombstones'):
            os.makedirs(os.path.join(spool_dir, sub), exist_ok=True)
        self._entries = {}
        self._ram = collections.OrderedDict()     # handle -> (graphlet, bytes, last use)
        self._ram_bytes = 0
        # handle -> (bytes without the caches, the caches' signature when measured)
        self._meter = {}
        self._persisted = {}                       # handle -> accessed time on disk
        for name in os.listdir(os.path.join(spool_dir, 'entries')):
            if name.endswith('.json'):
                try:
                    e = self._read_entry(os.path.join(spool_dir, 'entries', name))
                except (OSError, ValueError):
                    continue
                self._entries[e.handle] = e

    # ---------------------------------------------------------------- paths

    def _body_path(self, digest):
        return os.path.join(self.spool_dir, 'bodies', digest + '.mgt')

    def _entry_path(self, handle):
        return os.path.join(self.spool_dir, 'entries', handle + '.json')

    def _tomb_path(self, handle):
        return os.path.join(self.spool_dir, 'tombstones', handle + '.json')

    def _read_entry(self, path):
        with open(path, 'r', encoding='utf-8') as f:
            j = json.load(f)
        return Entry(store=self, **j)

    def _write_entry(self, e):
        j = {k: getattr(e, k) for k in ('handle', 'request', 'index', 'envelope',
                                        'seed_summary', 'digest', 'bytes', 'created',
                                        'accessed', 'view', 'derived_from', 'delivery',
                                        'parent')}
        _write_atomic(self._entry_path(e.handle),
                      json.dumps(j, sort_keys=True, ensure_ascii=False).encode('utf-8'))

    # ---------------------------------------------------------------- put

    def _new_handle(self):
        while True:
            h = 'g_' + secrets.token_hex(6)
            if h not in self._entries and not os.path.exists(self._tomb_path(h)):
                return h

    def _store_body(self, body):
        data = utf8_bytes(body)
        if self.max_body_bytes is not None and len(data) > self.max_body_bytes:
            raise StoreLimitExceeded('the graphlet body has %d bytes, over the store limit of '
                                     '%d: rejected (never truncated)'
                                     % (len(data), self.max_body_bytes))
        digest = hashlib.sha256(data).hexdigest()
        path = self._body_path(digest)
        if not os.path.exists(path):
            _write_atomic(path, data)
        return digest, len(data)

    def put(self, response, request, index=0, *, source=None, resident=True, delivery=None,
            parent=None, derived_from=None):
        """Store results[index] of a detail: graphlet response -> an opaque handle. The
        body is parsed first: an invalid or truncated one raises and nothing is
        stored. |source| names the index (the client) it came from, for replays."""
        result = response['results'][index]
        g = from_response(result, response)   # validates the whole body
        return self.put_parsed(g, result['graphlet'], request, source=source,
                               resident=resident, delivery=delivery, parent=parent,
                               derived_from=derived_from)

    def put_parsed(self, g, body, request, *, source=None, resident=True, delivery=None,
                   parent=None, derived_from=None):
        """Store a graphlet already parsed from |body| (the caller validated it).
        |resident|: keep the model in RAM (False: a spooled body, parsed again on
        demand). |delivery|: the delivery that actually happened when it is not the O
        record's (the MCP layer spooled it); it replaces the outcome's delivery dimension
        in the model and in the stored summary, so one value exists."""
        if delivery is not None:
            set_delivery(g, delivery)
        if derived_from is not None:
            g.derived_from = derived_from
        return self._put_graphlet(g, body, request, source, resident=resident,
                                  delivery=delivery, parent=parent)

    def put_all(self, response, request, *, source=None):
        """-> one handle per result, None where the seed failed (no graphlet)."""
        return [self.put(response, request, i, source=source) if 'graphlet' in r else None
                for i, r in enumerate(response.get('results', []))]

    def put_graphlet(self, g, request=None, derived_from=None, source=None):
        """Store a graphlet built elsewhere (a loaded file, a view's backing body)."""
        g.derived_from = derived_from or g.derived_from
        return self._put_graphlet(g, dump(g, envelope=False), request, source)

    def _put_graphlet(self, g, body, request, source=None, *, resident=True, delivery=None,
                      parent=None):
        digest, nbytes = self._store_body(body)
        now = self.clock()
        h = self._new_handle()
        ident = self._identity(g)
        if source is not None:
            ident['source'] = source
        # the envelope of THIS seed: usage reduced to the totals and its per_seed entry, the
        # rule of a saved file's J line (parser.seed_envelope)
        e = Entry(handle=h, request=request, index=ident,
                  envelope=seed_envelope(g.envelope, g.seed_index) if g.envelope else {},
                  seed_summary=g.seed_summary or {},
                  digest=digest, bytes=nbytes, created=now, accessed=now, view=g.view,
                  derived_from=g.derived_from, delivery=delivery, parent=parent, store=self)
        self._write_entry(e)
        self._entries[h] = e
        if resident:
            self._remember(h, g, now)
        return h

    @staticmethod
    def _identity(g):
        return {'ns': g.index_ns, 'fp': g.index_fp, 'meta_fp': g.index_meta_fp,
                'release': (g.envelope or {}).get('release') or None}

    # ---------------------------------------------------------------- get

    def get(self, handle):
        e = self._entries.get(handle) if isinstance(handle, str) else None
        if e is None:
            raise self._unknown(handle)
        now = self.clock()
        if now - e.accessed > self.ttl_disk_s:
            self._expire(handle)
            raise self._unknown(handle)
        # the disk TTL counts from the last use, also across restarts: the entry file is
        # rewritten when its stored time lags by a tenth of the TTL (not on every read)
        stale = now - self._persisted.get(handle, e.accessed) > self.ttl_disk_s / 10
        e.accessed = now
        if stale:
            self._write_entry(e)
            self._persisted[handle] = now
        return e

    def _unknown(self, handle):
        if not isinstance(handle, str) or not _HANDLE_RE.match(handle):
            return UnknownHandle(handle)
        try:
            with open(self._tomb_path(handle), 'r', encoding='utf-8') as f:
                tomb = json.load(f)
            return UnknownHandle(handle, True, tomb.get('request'), tomb.get('index'))
        except (OSError, ValueError):
            return UnknownHandle(handle)

    def graphlet(self, handle):
        e = self.get(handle)
        now = e.accessed
        got = self._ram.get(handle)
        if got is not None and now - got[2] <= self.ttl_ram_s:
            self._ram.move_to_end(handle)
            self._ram[handle] = (got[0], got[1], now)
            self.refresh(handle)
            return got[0]
        if got is not None:
            self._drop_ram(handle)
        with open(self._body_path(e.digest), 'r', encoding='utf-8', newline='') as f:
            body = f.read()
        g = parse(body)
        g.envelope = e.envelope or None
        g.seed_summary = e.seed_summary or None
        g.view = e.view
        g.derived_from = e.derived_from
        if e.delivery is not None:
            g.outcome.delivery = e.delivery
        self._remember(handle, g, now)
        return g

    def view(self, handle):
        """The entry's graphlet, or the GraphletView it stores."""
        g = self.graphlet(handle)
        if g.view:
            from .ops import GraphletView
            spec = g.view
            return GraphletView(g, spec.get('selectors') or [], spec.get('arm'),
                                spec.get('mode', 'any'), spec.get('of'))
        return g

    def _remember(self, handle, g, now):
        from . import ops
        size = g.memory_bytes()
        if handle in self._ram:
            self._drop_ram(handle)
        if size > self.max_ram_bytes:
            return      # never resident above the whole budget: parsed on demand
        self._ram[handle] = (g, size, now)
        self._meter[handle] = (size - ops.cache_bytes(g), ops.cache_signature(g))
        self._ram_bytes += size
        self.refresh()

    def refresh(self, handle=None):
        """Charge the resident models (or one) again where their caches grew since they
        were measured -- a query builds derived caches (paths, evidence, label
        summaries, ...) on the model it ran on -- then evict, least recently used first,
        until the budget holds. Costs a signature per model, and the caches of the
        models whose caches changed."""
        from . import ops
        for h in ([handle] if handle is not None else list(self._ram)):
            got = self._ram.get(h)
            if got is None:
                continue
            g, size, used = got
            base, sig = self._meter.get(h, (size, None))
            now_sig = ops.cache_signature(g)
            if now_sig != sig:
                new = base + ops.cache_bytes(g)
                self._ram_bytes += new - size
                self._ram[h] = (g, new, used)
                self._meter[h] = (base, now_sig)
        while len(self._ram) > 1 and (len(self._ram) > self.max_handles
                                      or self._ram_bytes > self.max_ram_bytes):
            oldest = next(iter(self._ram))
            self._drop_ram(oldest)
        if len(self._ram) == 1 and self._ram_bytes > self.max_ram_bytes:
            # one model whose caches outgrew the whole budget: parsed again on demand
            self._drop_ram(next(iter(self._ram)))

    def _drop_ram(self, handle):
        got = self._ram.pop(handle, None)
        self._meter.pop(handle, None)
        if got is not None:
            self._ram_bytes -= got[1]

    def in_ram(self, handle):
        return handle in self._ram

    # ---------------------------------------------------------------- lifecycle

    def free(self, handle):
        """Delete the entry (and its body when no other entry holds it). Freed entries
        are not replayable: the caller asked for them to go."""
        e = self._entries.pop(handle, None) if isinstance(handle, str) else None
        if e is None:
            raise self._unknown(handle)
        self._drop_ram(handle)
        try:
            os.unlink(self._entry_path(handle))
        except OSError:
            pass
        self._gc_body(e.digest)

    def _gc_body(self, digest):
        if not any(x.digest == digest for x in self._entries.values()):
            try:
                os.unlink(self._body_path(digest))
            except OSError:
                pass

    def _expire(self, handle):
        e = self._entries.pop(handle)
        self._drop_ram(handle)
        tomb = {'handle': handle, 'request': e.request, 'index': e.index,
                'expired': self.clock()}
        _write_atomic(self._tomb_path(handle), json.dumps(tomb).encode('utf-8'))
        try:
            os.unlink(self._entry_path(handle))
        except OSError:
            pass
        self._gc_body(e.digest)

    def sweep(self):
        """Apply the TTLs now: -> {ram_dropped, expired, tombstones_removed}. A tombstone
        (an expired entry's request, for replay) goes ttl_tomb_s after the expiry, so the
        spool does not grow without bound."""
        now = self.clock()
        ram = [h for h, (_, _, used) in self._ram.items() if now - used > self.ttl_ram_s]
        for h in ram:
            self._drop_ram(h)
        expired = [h for h, e in self._entries.items() if now - e.accessed > self.ttl_disk_s]
        for h in expired:
            self._expire(h)
        removed = []
        tombs = os.path.join(self.spool_dir, 'tombstones')
        for name in sorted(os.listdir(tombs)):
            if not name.endswith('.json'):
                continue
            path = os.path.join(tombs, name)
            try:
                with open(path, 'r', encoding='utf-8') as f:
                    when = float(json.load(f).get('expired'))
            except (OSError, ValueError, TypeError, AttributeError):
                continue
            if now - when > self.ttl_tomb_s:
                try:
                    os.unlink(path)
                    removed.append(name[:-len('.json')])
                except OSError:
                    pass
        return {'ram_dropped': ram, 'expired': expired, 'tombstones_removed': removed}

    def cursor_secret(self):
        """The spool's secret for the MCP layer's cursors (created once, 0600). A cursor
        must outlive a restart for as long as its handle does: handles live in the
        spool, so the secret does too."""
        path = os.path.join(self.spool_dir, 'cursor.key')
        try:
            with open(path, 'rb') as f:
                key = f.read()
            if len(key) >= 16:
                return key
        except FileNotFoundError:
            pass
        # written aside and linked in: a concurrent reader never sees a partial key, and
        # two creators agree on whichever link landed first
        fd, tmp = tempfile.mkstemp(prefix='.key-', dir=self.spool_dir)
        try:
            os.fchmod(fd, 0o600)
            with os.fdopen(fd, 'wb') as f:
                f.write(secrets.token_bytes(32))
            try:
                os.link(tmp, path)
            except FileExistsError:
                pass
        finally:
            os.unlink(tmp)
        with open(path, 'rb') as f:
            return f.read()

    def list(self):
        return [e.meta() for e in sorted(self._entries.values(), key=lambda x: x.created)]

    def save(self, handle, path):
        """The standalone .mgt file (H, J, body) of an entry."""
        g = self.graphlet(handle)
        from .parser import save
        return save(g, path)

    def load(self, path, request=None):
        """A saved .mgt file -> a new handle (its J line restores envelope and view)."""
        with open(path, 'r', encoding='utf-8', newline='') as f:
            g = parse(f.read())
        return self.put_graphlet(g, request if request is not None else
                                 self.request_of(g))

    @staticmethod
    def request_of(g):
        """The replayable request of a graphlet with an envelope: its seed as validated
        (sequence, the seed labels by name) and the normalized strategy. None when there
        is none -- also when the seed labels' names cannot be verified to resolve back
        to them (ops.resubmittable_names): a replay must not run under other labels."""
        if not g.has_envelope:
            return None
        seed = {'sequence': g.seed.sequence}
        meta = g.seed_summary.get('seed') or {}
        if g.mode == 'constrain' and not meta.get('labels_from_seed'):
            from .model import UnverifiableLabelName
            from .ops import resubmittable_names
            try:
                seed['labels'] = resubmittable_names(
                    g, [{'id': l.id} for l in g.labels[:g.seed.num_seed_labels]])
            except UnverifiableLabelName:
                return None
        if meta.get('seed_id'):
            seed['seed_id'] = meta['seed_id']
        strategy = dict(g.envelope.get('strategy') or {})
        strategy.pop('clamped', None)
        req = {'seeds': [seed], 'strategy': strategy}
        if g.envelope.get('release'):
            req['release'] = g.envelope['release']
        return req

    def __contains__(self, handle):
        if not isinstance(handle, str):
            return False
        return handle in self._entries

    def j_of(self, handle):
        return j_object(self.graphlet(handle))
