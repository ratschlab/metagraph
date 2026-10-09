"""GraphletStore: handles, an LRU of parsed graphlets in RAM and a disk spool (§5, §6).

ENTRY identity and BODY identity are separate: an entry handle is opaque
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
never stored or returned as a graphlet. Local limits (DESIGN §5.4): with parse limits
(parse_limits, or an explicit parse budget) a parse that stops keeps the body -- stored
after the checks that need no parse (graphlet_bytes, graphlet_lines, the H record and the
Z line count) as an entry marked parsed: false, which no local answer is derived from until
a parse completes; it can always be copied out without one (save_body, body_text). RAM holds parsed graphlets under max_ram_mb /
max_handles (least recently used first out) and drops one not used for ttl_ram_s; a
graphlet larger than the whole RAM budget is never kept resident (it is parsed on
demand), and a spooled one is not made resident when it is stored. A resident model is
charged its deep heap footprint (Graphlet.memory_bytes()), and again when the caches its
queries build have grown (refresh(): on every access, and after every MCP tool call),
so that max_ram_mb bounds the heap the resident models hold. The disk keeps an
entry for ttl_disk_s after its last use, and a tombstone for ttl_tomb_s after its
expiry. The clock is injectable.

Processes sharing a spool each keep their own RAM and index of the entries; the spool
stays consistent between them: a body is deleted
only when no entry file names it any more, a handle another process stored is read from
its entry file on first use (get, free, `in`; list() lists the spool's entries), one
another process freed or expired answers as unknown (as an expired, replayable one where
its tombstone keeps a request), and an expiry re-reads the entry file's last use first. A
body missing from the spool expires its entry, on a parse and on a copy without one.

The body's lifecycle is serialized across the processes (and the threads of one store) by an
exclusive lock of the spool (an fcntl.flock of spool/.lock): a put holds it from the check
of its body to the write of its entry file, a free or an expiry from its scan of the entry
files to the unlink of the body. Without it a free in one process could delete a body
between another's check of it and the entry naming it -- the put would return a handle whose
body is gone (about half of the handles of three processes storing one body in a loop).
Limits: a platform without fcntl (Windows), a spool whose file system refuses the lock file
or flock (read-only; some network file systems: the store's _lock.unlocked names the error)
has the threads' lock only, and a process of an earlier library version sharing the spool
takes none, so the race remains against it (and the orphan pass below spares a body only by
its age).
A body that no entry file names any more is deleted by sweep() once it is older than
_ORPHAN_GRACE_S (a free that kept it for an entry another process had already removed, or a
process killed between a body and its entry, would otherwise leave it for good). Body GC
reads each entry file of another process once (its digest is cached by name and inode, so a
free or a sweep does not list and parse every such file again). Body GC
deletes nothing while it cannot tell which bodies the entries name: when the entry files
cannot be listed, or one of them is there but cannot be read (EACCES, EMFILE, EIO) -- a
free then keeps its body and a sweep its orphans, until the file can be read or is
removed. An entry file that is no JSON (damaged; no store can read it either) names no
body.

Damage: a body or entry file is written to a temporary name, flushed to the device
(fsync) and renamed, so a crash leaves the old file or the new one, never a new name on an
empty or torn file -- on Linux; macOS's fsync does not flush the drive's own cache
(F_FULLFSYNC is not used), so a power loss there can still lose the latest writes. A
stored body is checked against its digest (its name) when it is read -- a parse on demand
and a copy without one -- and before a put reuses it (once per process and file); a body
that fails is deleted and its entry expires (replayable when it kept a request), so the
next put of the same body writes it again. A damaged entry file is skipped.
"""

import collections
import dataclasses
import errno
import hashlib
import json
import os
import re
import secrets
import tempfile
import threading
import time
from dataclasses import dataclass, field
from typing import Any, Optional

try:
    import fcntl
except ImportError:          # Windows: the threads' lock only (stated in the module text)
    fcntl = None

from . import budget as _B
from . import coords
from ._codec import CodecError, GraphletFormatError, parse_int
from .budget import LocalBudget, LocalBudgetExceeded, LocalLimits
from .parser import (_check_transport, _replace_atomically, _write_all, dump, from_response, j_object, parse,
                     seed_envelope, utf8_bytes)

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
    # local limits: False when the body was stored after a parse that stopped (its integrity
    # checked without one); True once a parse of it completed
    parsed: bool = True
    store: Any = field(default=None, repr=False, compare=False)
    # the bound standalone_text() charges its J line by (_j_bound()): set whenever the
    # entry is written or read, never stored in its file
    j_bound: int = field(default=0, repr=False, compare=False)

    @property
    def graphlet(self):
        return self.store.graphlet(self.handle)

    def meta(self):
        out = {'handle': self.handle, 'index': self.index, 'digest': self.digest,
               'bytes': self.bytes, 'created': self.created, 'accessed': self.accessed,
               'view': self.view, 'derived_from': self.derived_from,
               'parent': self.parent, 'delivery': self.delivery,
               'seed_id': (self.seed_summary.get('seed') or {}).get('seed_id')}
        if not self.parsed:
            out['parsed'] = False
        return out


# the fields an entry file may carry (the Entry dataclass without its store); a file of
# another version keeps what this one knows (an unknown key, such as 'parsed' read by an
# older library, must not make the whole store fail to open)
_ENTRY_FIELDS = frozenset(f.name for f in dataclasses.fields(Entry)) - {'store', 'j_bound'}
_ENTRY_REQUIRED = frozenset(f.name for f in dataclasses.fields(Entry)
                            if f.default is dataclasses.MISSING
                            and f.default_factory is dataclasses.MISSING)

# what the J line adds to its sources besides their own text (at most): the keys it puts
# them under, with their punctuation
_J_KEYS = len(',"results":[],"view":,"derived_from":')


def _j_bound(e):
    """An upper bound, in UTF-8 bytes, of the J line standalone_text() writes for |e|
    (0 when it writes none): the line's sources -- the envelope, the seed summary, the
    view and derived_from -- as JSON with the spaces the line leaves out, plus the keys
    it puts them under. A function of what the entry holds only, so that the same export
    of the same entry is charged the same on every run (budgets are deterministic). Not
    the entry file's size, which also holds the created and accessed times: their printed
    length follows the clock (1791204329.5 against 1791204329.6234567), so a memory limit
    could stop one run of an export while an identical run answers."""
    if not (e.envelope or e.seed_summary or e.view is not None
            or e.derived_from is not None):
        return 0
    text = json.dumps([e.envelope, e.seed_summary, e.view, e.derived_from],
                      ensure_ascii=False)
    return _J_KEYS + (len(text) if text.isascii() else len(text.encode('utf-8')))


def set_delivery(g, delivery):
    """The outcome's delivery dimension as it actually happened (the model's O record and
    the attached per-seed summary), so that summary(), the stored entry and the tool's
    own `delivery` carry one value."""
    if g.outcome is not None:
        g.outcome.delivery = delivery
    if g.seed_summary and isinstance(g.seed_summary.get('outcome'), dict):
        g.seed_summary = dict(g.seed_summary, outcome=dict(g.seed_summary['outcome'],
                                                            delivery=delivery))


def _check_frame(body):
    """The checks of a body that need no parse: it ends with a line feed, its first line
    is an MGT H record of this version, its last a Z record counting its lines. -> the H
    record's fields. GraphletFormatError otherwise."""
    if not body.endswith('\n'):
        raise GraphletFormatError(0, 'truncated document: it does not end with a line feed')
    first = body.find('\n')
    head = body[:first].split(' ')
    if len(head) != 16 or head[0] != 'H' or head[1] != 'mgt' or head[2] != '1':
        raise GraphletFormatError(1, 'expected the H record of MGT v1')
    last = body.rfind('\n', 0, len(body) - 1) + 1
    z = body[last:-1].split(' ')
    n = body.count('\n')
    try:
        if len(z) != 2 or z[0] != 'Z' or parse_int(z[1]) != n:
            raise GraphletFormatError(n, 'the Z record does not count the body\'s %d lines' % n)
        parse_int(head[11])
    except CodecError as e:
        raise GraphletFormatError(n, str(e)) from None
    return head


def _unparsed_head(result, response):
    """(the body as text, its H record's fields, the index identity it states) of
    results[i], read WITHOUT a parse -- after the checks that need none (graphlet_bytes and
    graphlet_lines against the body, the H record, the Z line counting its lines): what
    put_unparsed() stores, and what an unparsed replay or continuation checks against the
    entry it derives from before storing anything."""
    body = result.get('graphlet')
    if body is None:
        raise ValueError('this result carries no graphlet')
    if isinstance(body, (bytes, bytearray)):
        body = bytes(body).decode('utf-8')
    _check_transport(body, result)
    head = _check_frame(body)
    ident = {'ns': None if head[13] == '*' else head[13],
             'fp': None if head[14] == '*' else head[14],
             'meta_fp': None if head[15] == '*' else head[15],
             'release': (response or {}).get('release') or None}
    return body, head, ident


def unparsed_identity(result, response):
    """The index identity {ns, fp, meta_fp, release} a result's H record (and the response's
    release) states, read without a parse (_unparsed_head())."""
    return _unparsed_head(result, response)[2]


def _o_with_delivery(body, start, delivery):
    """(where the O record begins in |body|, where it ends, the record with its delivery
    token replaced -- what a parsed model with the entry's delivery dumps), searched from
    |start|; (None, None, None) without an O record."""
    code = {'inline': 'i', 'spooled': 's', 'paged': 'p'}[delivery]
    i = start
    while i < len(body):
        j = body.find('\n', i)
        if j < 0:
            j = len(body)
        if body.startswith('O ', i):
            f = body[i:j].split(' ')
            f[4] = code
            return i, j, ' '.join(f)
        i = j + 1
    return None, None, None


# a body no entry file names is deleted by sweep() only once it is this old (seconds, by its
# modification time): a put of a library version without the spool lock writes its body
# before its entry, and a body it has just written must not be taken for an orphan
_ORPHAN_GRACE_S = 3600
_HASH_CHUNK = 1 << 20


def _write_atomic(path, data):
    d = os.path.dirname(path)
    fd, tmp = tempfile.mkstemp(prefix='.tmp-', dir=d)
    # one write, without the 128 KiB buffer a buffered file holds (parser._write_all()),
    # flushed before the rename: a rename to a new name is not ordered after the data on
    # every file system (ext4's delayed allocation), and a crash could leave the new name
    # on an empty file -- a body that every later put of the same body would reuse
    _replace_atomically(tmp, path, lambda: _write_all(fd, data, sync=True))


def _file_sha256(path):
    """The sha256 of the file at |path| (hex), read in pieces; None when it is gone."""
    h = hashlib.sha256()
    try:
        with open(path, 'rb', buffering=0) as f:
            while True:
                chunk = f.read(_HASH_CHUNK)
                if not chunk:
                    break
                h.update(chunk)
    except FileNotFoundError:
        return None
    return h.hexdigest()


# a lock file the spool's file system will not create (a read-only spool: nothing can be
# deleted there either) or a flock it does not support (some network file systems): the
# store goes on with the threads' lock only, as it did before the spool lock, and says so
# in unlocked (the module text states the limit)
_NO_LOCK_FILE = frozenset((errno.EROFS, errno.EACCES, errno.EPERM))
_NO_FLOCK = frozenset(x for x in (getattr(errno, n, None) for n in (
    'ENOLCK', 'EOPNOTSUPP', 'ENOTSUP', 'ENOSYS', 'EINVAL')) if x is not None)


class _SpoolLock:
    """The exclusive lock of one spool: an fcntl.flock of spool/.lock (other processes, and
    other stores of the same spool in this process: a flock belongs to its open file) and
    an RLock (the threads of this store). Re-entrant within a thread, so a locked step may
    call another (free -> GC, get -> expire). Two stores of one spool must not nest their
    locked calls in one thread: the second flock would wait for the first. |unlocked|: why
    the last hold had no process lock (None: it had one, or the platform has no fcntl)."""

    __slots__ = ('path', '_rlock', '_depth', '_fd', 'unlocked')

    def __init__(self, spool_dir):
        self.path = os.path.join(spool_dir, '.lock')
        self._rlock = threading.RLock()
        self._depth = 0
        self._fd = None
        self.unlocked = None

    def __enter__(self):
        self._rlock.acquire()
        if self._depth == 0 and fcntl is not None:
            try:
                self._fd = self._flock()
            except BaseException:
                self._rlock.release()
                raise
        self._depth += 1
        return self

    def _flock(self):
        """Open and flock the lock file -> its descriptor; None where the file system
        refuses either (_NO_LOCK_FILE, _NO_FLOCK: noted in unlocked)."""
        try:
            fd = os.open(self.path, os.O_RDWR | os.O_CREAT, 0o600)
        except OSError as e:
            if e.errno in _NO_LOCK_FILE:
                self.unlocked = errno.errorcode.get(e.errno, str(e.errno))
                return None
            raise
        try:
            fcntl.flock(fd, fcntl.LOCK_EX)
        except BaseException as e:
            os.close(fd)
            if isinstance(e, OSError) and e.errno in _NO_FLOCK:
                self.unlocked = errno.errorcode.get(e.errno, str(e.errno))
                return None
            raise
        self.unlocked = None
        return fd

    def __exit__(self, *exc):
        self._depth -= 1
        try:
            if self._depth == 0 and self._fd is not None:
                fd, self._fd = self._fd, None
                try:
                    fcntl.flock(fd, fcntl.LOCK_UN)
                finally:
                    os.close(fd)
        finally:
            self._rlock.release()
        return False


class GraphletStore:
    def __init__(self, spool_dir, max_ram_mb=512, max_handles=64, ttl_ram_s=1800,
                 ttl_disk_s=7 * 86400, *, clock=time.time, max_body_mb=None,
                 ttl_tomb_s=30 * 86400, parse_limits=None):
        self.spool_dir = spool_dir
        if parse_limits is not None and not isinstance(parse_limits, LocalLimits):
            raise TypeError('parse_limits is a LocalLimits, not %r' % (parse_limits,))
        # local limits: every parse the store performs (validation in put, a parse on demand,
        # load) runs under a fresh budget of these limits; None: unbudgeted
        self.parse_limits = parse_limits
        # the usage of the last parse the store ran under a budget ({work_units,
        # memory_bytes}, None when the last access parsed nothing)
        self.last_parse = None
        self.max_ram_bytes = int(max_ram_mb * (1 << 20))
        self.max_handles = max_handles
        self.ttl_ram_s = ttl_ram_s
        self.ttl_disk_s = ttl_disk_s
        self.ttl_tomb_s = ttl_tomb_s
        self.clock = clock
        self.max_body_bytes = None if max_body_mb is None else int(max_body_mb * (1 << 20))
        for sub in ('bodies', 'entries', 'tombstones'):
            os.makedirs(os.path.join(spool_dir, sub), exist_ok=True)
        # the body lifecycle across the processes sharing the spool
        self._lock = _SpoolLock(spool_dir)
        self._entries = {}
        # digest -> the handles of this index whose entries name it: the body GC's own test
        # (not an any() over every entry per freed body)
        self._by_digest = {}
        # entry file name -> (its inode, the digest it names), for the entry files of other
        # processes: an entry's digest never changes, and a file is only ever replaced
        # whole (a new inode), so each is read once (not again per freed body)
        self._foreign = {}
        # digest -> (inode, mtime_ns, size) of the body file last checked against it (a
        # put reuses a body once its content is known to be the digest's)
        self._verified = {}
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
                self._index(e)
                # the accessed time the entry file holds: the disk TTL's staleness test
                # reads it (unseeded, a handle used more often than every tenth of the TTL
                # would never have its use written, and would expire after a restart)
                self._persisted[e.handle] = e.accessed

    # ---------------------------------------------------------------- paths

    def _body_path(self, digest):
        return os.path.join(self.spool_dir, 'bodies', digest + '.mgt')

    def _entry_path(self, handle):
        return os.path.join(self.spool_dir, 'entries', handle + '.json')

    def _tomb_path(self, handle):
        return os.path.join(self.spool_dir, 'tombstones', handle + '.json')

    def _index(self, e):
        """Add entry |e| to this process's index (by handle and by the digest it names)."""
        old = self._entries.get(e.handle)
        if old is not None:
            self._unindex(e.handle)
        self._entries[e.handle] = e
        self._by_digest.setdefault(e.digest, set()).add(e.handle)

    def _unindex(self, handle):
        """Remove |handle| from this process's index -> its entry (None: not indexed)."""
        e = self._entries.pop(handle, None)
        if e is not None:
            hs = self._by_digest.get(e.digest)
            if hs is not None:
                hs.discard(handle)
                if not hs:
                    del self._by_digest[e.digest]
        return e

    def _read_entry(self, path):
        """An entry file -> Entry; ValueError for a file that is no entry (not an object,
        a required field missing), which the callers skip. Keys this version does not
        know are left out, not fatal."""
        with open(path, 'r', encoding='utf-8') as f:
            j = json.load(f)
        if not isinstance(j, dict):
            raise ValueError('an entry file holds an object, not %s' % type(j).__name__)
        missing = _ENTRY_REQUIRED - set(j)
        if missing:
            raise ValueError('the entry lacks %s' % ', '.join(sorted(missing)))
        e = Entry(store=self, **{k: v for k, v in j.items() if k in _ENTRY_FIELDS})
        e.j_bound = _j_bound(e)
        return e

    def _write_entry(self, e):
        j = {k: getattr(e, k) for k in ('handle', 'request', 'index', 'envelope',
                                        'seed_summary', 'digest', 'bytes', 'created',
                                        'accessed', 'view', 'derived_from', 'delivery',
                                        'parent')}
        if not e.parsed:
            j['parsed'] = False             # written only when false: older entries read True
        e.j_bound = _j_bound(e)
        _write_atomic(self._entry_path(e.handle),
                      json.dumps(j, sort_keys=True, ensure_ascii=False).encode('utf-8'))
        # what the file now says: get() rewrites it once this lags by a tenth of the TTL
        self._persisted[e.handle] = e.accessed

    # ---------------------------------------------------------------- put

    def _new_handle(self):
        while True:
            h = 'g_' + secrets.token_hex(6)
            if h not in self._entries and not os.path.exists(self._tomb_path(h)):
                return h

    def _store_body(self, body):
        """Make the spool hold |body| -> (digest, bytes). Called under the spool lock, with
        the entry that names the body written before the lock is let go."""
        data = utf8_bytes(body)
        if self.max_body_bytes is not None and len(data) > self.max_body_bytes:
            raise StoreLimitExceeded('the graphlet body has %d bytes, over the store limit of '
                                     '%d: rejected (never truncated)'
                                     % (len(data), self.max_body_bytes))
        digest = hashlib.sha256(data).hexdigest()
        path = self._body_path(digest)
        # an existing file is reused only once its content is the digest's: a crash can
        # leave it empty or torn, and every later put of the same body -- a re-fetch that
        # had just validated it -- would point its new entry at the damaged file
        if not self._body_intact(path, digest, len(data)):
            _write_atomic(path, data)
            self._note_verified(path, digest)
        return digest, len(data)

    def _body_intact(self, path, digest, n):
        """Whether the body file at |path| holds the |n| bytes of |digest|: its size, then
        its hash -- read once per process and file (inode, modification time, size), since
        a body file is only ever replaced whole."""
        try:
            st = os.stat(path)
        except FileNotFoundError:
            return False
        if st.st_size != n:
            return False
        key = (st.st_ino, st.st_mtime_ns, st.st_size)
        if self._verified.get(digest) == key:
            return True
        if _file_sha256(path) != digest:
            return False
        self._verified[digest] = key
        return True

    def _note_verified(self, path, digest):
        try:
            st = os.stat(path)
        except OSError:
            return
        self._verified[digest] = (st.st_ino, st.st_mtime_ns, st.st_size)

    def _parse_budget(self, parse_budget):
        """The budget a parse of the store runs under: the caller's, else a fresh one of
        the store's parse limits, else None (unbudgeted)."""
        if parse_budget is not None:
            return parse_budget
        if self.parse_limits is not None:
            return LocalBudget(self.parse_limits)
        return None

    def _parsing(self, fn, pb):
        """fn(pb): a parse of the store under its own budget, or unbudgeted -- never
        charged to an ambient budget the caller runs under (the parse is the store's,
        not the caller's operation; the tools report it as local.parsed)."""
        try:
            if pb is None:
                with _B.unbudgeted():
                    return fn(None)
            return fn(pb)
        finally:
            self.last_parse = None if pb is None else pb.call_usage()

    def put(self, response, request, index=0, *, source=None, resident=True, delivery=None,
            parent=None, derived_from=None, parse_budget=None):
        """Store results[index] of a detail: graphlet response -> an opaque handle. The
        body is parsed first: an invalid or truncated one raises and nothing is
        stored. |source| names the index (the client) it came from, for replays.

        Under parse limits (parse_limits, or |parse_budget|) a parse that stops raises its
        LocalBudgetExceeded AFTER storing the body as an unparsed entry (parsed: false)
        when it passes the checks that need no parse: e.handle names it."""
        result = response['results'][index]
        try:
            g = self._parsing(lambda b: from_response(result, response, budget=b),
                              self._parse_budget(parse_budget))   # validates the whole body
        except LocalBudgetExceeded as e:
            e.handle = self.put_unparsed(result, response, request, source=source,
                                         delivery=delivery, parent=parent,
                                         derived_from=derived_from)
            raise
        return self.put_parsed(g, result['graphlet'], request, source=source,
                               resident=resident, delivery=delivery, parent=parent,
                               derived_from=derived_from)

    def put_unparsed(self, result, response, request, *, source=None, delivery=None,
                     parent=None, derived_from=None):
        """Store results[i] WITHOUT parsing its body (local limits: a parse stopped on its
        budget): only after the checks that need no parse -- graphlet_bytes and
        graphlet_lines against the body, its H record, and its Z line counting its lines
        -- so a body cut in transport is still refused. The entry is marked parsed:
        false: its identity is read from the H record, and the server's per-seed summary
        is kept as it came. A local answer needs a parse that completes first."""
        body, head, ident = _unparsed_head(result, response)
        if source is not None:
            ident['source'] = source
        envelope = seed_envelope({k: v for k, v in response.items() if k != 'results'},
                                 parse_int(head[11]))
        summary = {k: v for k, v in result.items() if k != 'graphlet'}
        if delivery is not None and isinstance(summary.get('outcome'), dict):
            summary['outcome'] = dict(summary['outcome'], delivery=delivery)
        with self._lock:              # from the body's check to its entry's write
            digest, nbytes = self._store_body(body)
            now = self.clock()
            h = self._new_handle()
            e = Entry(handle=h, request=request, index=ident, envelope=envelope,
                      seed_summary=summary, digest=digest, bytes=nbytes, created=now,
                      accessed=now, derived_from=derived_from, delivery=delivery,
                      parent=parent, parsed=False, store=self)
            self._write_entry(e)
            self._index(e)
        return h

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

    def put_view(self, backing, g, request=None, derived_from=None, source=None):
        """Store the view |g| (a copy of the graphlet of entry |backing| with its view
        spec) WITHOUT dumping its body: the view entry shares the backing entry's stored
        body, which is what dump(g, envelope=False) writes for every body the store holds
        (local limits: a view costs no copy of the body). The backing body must still be in
        the spool: when it is gone, the backing entry expires (UnknownHandle)."""
        with self._lock:
            # the backing entry read, its body checked and the view's entry written under
            # one hold of the spool lock: a free of the backing entry in another process
            # could delete the body before the view's entry names it
            e = self.get(backing)
            if not os.path.exists(self._body_path(e.digest)):
                self._expire(backing)
                raise self._unknown(backing)
            g.derived_from = derived_from or g.derived_from
            now = self.clock()
            h = self._new_handle()
            ident = self._identity(g)
            if source is not None:
                ident['source'] = source
            v = Entry(handle=h, request=request, index=ident,
                      envelope=seed_envelope(g.envelope, g.seed_index) if g.envelope else {},
                      seed_summary=g.seed_summary or {}, digest=e.digest, bytes=e.bytes,
                      created=now, accessed=now, view=g.view, derived_from=g.derived_from,
                      delivery=e.delivery, store=self)
            self._write_entry(v)
            self._index(v)
        self._remember(h, g, now)
        return h

    def _put_graphlet(self, g, body, request, source=None, *, resident=True, delivery=None,
                      parent=None):
        ident = self._identity(g)
        if source is not None:
            ident['source'] = source
        # the envelope of THIS seed: usage reduced to the totals and its per_seed entry, the
        # rule of a saved file's J line (parser.seed_envelope)
        envelope = seed_envelope(g.envelope, g.seed_index) if g.envelope else {}
        with self._lock:              # from the body's check to its entry's write
            digest, nbytes = self._store_body(body)
            now = self.clock()
            h = self._new_handle()
            e = Entry(handle=h, request=request, index=ident, envelope=envelope,
                      seed_summary=g.seed_summary or {},
                      digest=digest, bytes=nbytes, created=now, accessed=now, view=g.view,
                      derived_from=g.derived_from, delivery=delivery, parent=parent,
                      store=self)
            self._write_entry(e)
            self._index(e)
        if resident:
            self._remember(h, g, now)
        return h

    @staticmethod
    def _identity(g):
        return {'ns': g.index_ns, 'fp': g.index_fp, 'meta_fp': g.index_meta_fp,
                'release': (g.envelope or {}).get('release') or None}

    # ---------------------------------------------------------------- get

    def get(self, handle):
        e = self._live(handle)
        if e is None:
            raise self._unknown(handle)
        now = self.clock()
        if now - e.accessed > self.ttl_disk_s and not self._fresher_on_disk(e, now):
            self._expire(handle)
            raise self._unknown(handle)
        # the disk TTL counts from the last use, also across restarts: the entry file is
        # rewritten when its stored time lags by a tenth of the TTL (not on every read)
        stale = now - self._persisted.setdefault(handle, e.accessed) > self.ttl_disk_s / 10
        e.accessed = now
        if stale:
            self._write_entry(e)
        return e

    def _live(self, handle):
        """The entry of |handle| as the spool holds it now, or None: one this store
        knows whose file is still there, else one another process sharing the spool
        stored after this store was opened (_adopt()). One this store knows whose file is
        gone -- freed or expired by another process -- is dropped: its tombstone (if any)
        answers, never a listed handle whose body may be gone. get(), free(), `in` and
        list() all see the spool so (a handle another process stored is known to them
        before a get() adopts it)."""
        e = self._entries.get(handle) if isinstance(handle, str) else None
        if e is None:
            return self._adopt(handle)
        if not os.path.exists(self._entry_path(handle)):
            self._forget(handle)
            return None
        return e

    def _forget(self, handle):
        """Drop |handle| from this process's index and RAM (its files are not touched)."""
        self._unindex(handle)
        self._drop_ram(handle)
        self._persisted.pop(handle, None)

    def _adopt(self, handle):
        """An entry another process sharing the spool stored after this store was opened
        (read from its file), or None."""
        if not isinstance(handle, str) or not _HANDLE_RE.match(handle):
            return None
        try:
            e = self._read_entry(self._entry_path(handle))
        except (OSError, ValueError, TypeError):
            return None
        if e.handle != handle:
            return None
        self._index(e)
        self._persisted[handle] = e.accessed
        return e

    def _fresher_on_disk(self, e, now):
        """Whether the entry file records a use within the disk TTL that this process
        has not seen (another process sharing the spool used the handle): then its time
        is taken, and the entry is not expired on this process's stale view."""
        try:
            on_disk = self._read_entry(self._entry_path(e.handle)).accessed
        except (OSError, ValueError, TypeError):
            return False
        if not isinstance(on_disk, (int, float)) or now - on_disk > self.ttl_disk_s:
            return False
        e.accessed = max(e.accessed, on_disk)
        self._persisted[e.handle] = on_disk
        return True

    def _unknown(self, handle):
        if not isinstance(handle, str) or not _HANDLE_RE.match(handle):
            return UnknownHandle(handle)
        try:
            with open(self._tomb_path(handle), 'r', encoding='utf-8') as f:
                tomb = json.load(f)
        except (OSError, ValueError):
            return UnknownHandle(handle)
        if not isinstance(tomb, dict):
            return UnknownHandle(handle)          # not a tombstone this store wrote
        # replayable only with a request to replay (an entry loaded from a file without
        # one has none: the hint would offer a replay that is then refused)
        request = tomb.get('request')
        return UnknownHandle(handle, bool(request), request, tomb.get('index'))

    def graphlet(self, handle, *, parse_budget=None):
        """The entry's parsed graphlet: the resident model, or a parse of its body (under
        the store's parse limits, or |parse_budget|: a stop raises LocalBudgetExceeded and
        keeps the entry and its body as they are; last_parse states the parse's usage)."""
        e = self.get(handle)
        now = e.accessed
        got = self._ram.get(handle)
        self.last_parse = None
        if got is not None and now - got[2] <= self.ttl_ram_s:
            self._ram.move_to_end(handle)
            self._ram[handle] = (got[0], got[1], now)
            self.refresh(handle)
            return got[0]
        if got is not None:
            self._drop_ram(handle)
        # the stored body, checked against its digest: decoded as the text file read
        # it (newline='': the same text)
        body = self._read_body(handle, e).decode('utf-8')
        # an entry stores {} for none: a graphlet had an envelope when either side holds
        # something (the rule of a J line, parser._attach_j)
        has_env = bool(e.envelope) or bool(e.seed_summary)
        # a summary with a coordinates block: every parse of the entry attaches the parsed
        # index, a re-parse of a parsed entry (its resident model expired or was evicted)
        # too -- so the model is the one from_response() made, with the same
        # memory_bytes() and cache_signature(); an entry without a block is parsed alone
        block = has_env and isinstance(e.seed_summary, dict) \
            and isinstance(e.seed_summary.get('coordinates'), dict)
        check = has_env and (not e.parsed or block)

        def parse_entry(b):
            if not check:
                return parse(body, budget=b)
            # the first whole parse of an entry stored unparsed (its parse stopped on a
            # budget), or any parse of one with a block: its record coordinates are
            # validated with it, as from_response() would have -- one parse call, so
            # a stop in either keeps the entry as it was (unparsed, or parsed and not
            # resident)
            if b is None:
                g_ = parse(body)
                g_.envelope, g_.seed_summary = e.envelope, e.seed_summary
                coords.attach(g_)
                return g_
            with b.scope('parse', ('export_mgt',)):
                g_ = parse(body, budget=b)
                g_.envelope, g_.seed_summary = e.envelope, e.seed_summary
                coords.attach(g_, b)
                return g_
        g = self._parsing(parse_entry, self._parse_budget(parse_budget))
        if not e.parsed:
            # the body now passed a whole parse: it is a graphlet like any other
            e.parsed = True
            self._write_entry(e)
        g.envelope = e.envelope if has_env else None
        g.seed_summary = e.seed_summary if has_env else None
        g.view = e.view
        g.derived_from = e.derived_from
        if e.delivery is not None:
            g.outcome.delivery = e.delivery
        self._remember(handle, g, now)
        return g

    def view(self, handle, *, parse_budget=None):
        """The entry's graphlet, or the GraphletView it stores."""
        g = self.graphlet(handle, parse_budget=parse_budget)
        if g.view is not None:
            # the same restoration and the same presence test as parser.load(): the model
            # is the shared resident one, so its view stays set
            from .ops import view_from_spec
            return view_from_spec(g, g.view)
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
        with self._lock:
            e = self._live(handle)
            if e is None:
                raise self._unknown(handle)
            self._forget(handle)
            try:
                os.unlink(self._entry_path(handle))
            except OSError:
                pass
            self._gc_bodies((e.digest,))

    def _held_here(self, digest):
        """Whether an entry of this process's index whose file is still there names body
        |digest|. One whose file another process removed (freed or expired) is dropped from
        the index here: otherwise it would keep the body for good -- a free deletes its body
        only when no entry of the index names it, stale ones included, and nothing would
        collect it after."""
        for h in sorted(self._by_digest.get(digest, ())):
            if os.path.exists(self._entry_path(h)):
                return True
            self._forget(h)
        return False

    def _named_on_disk(self):
        """The digests the spool's entry files name (each read once: _foreign), or None when
        the entries cannot be listed or one of them is there but cannot be read (then no
        body may be deleted). Under the spool lock, so no put is between its body and its
        entry. This process's entries are read from its index; one whose file is gone is
        dropped from it (_held_here's rule)."""
        entries = os.path.join(self.spool_dir, 'entries')
        named = set()
        seen = {}
        mine = set()
        try:
            with os.scandir(entries) as it:
                listing = [(d.name, d.inode()) for d in it if d.name.endswith('.json')]
        except OSError:
            return None
        for name, ino in listing:
            h = name[:-len('.json')]
            e = self._entries.get(h)
            if e is not None:
                named.add(e.digest)
                mine.add(h)
                continue
            got = self._foreign.get(name)
            if got is None or got[0] != ino:
                try:
                    with open(os.path.join(entries, name), 'r', encoding='utf-8') as f:
                        j = json.load(f)
                except FileNotFoundError:
                    continue              # removed since the listing: it names nothing now
                except OSError:
                    # an entry that is there but cannot be read now (EACCES, EMFILE, EIO):
                    # which body it names is unknown, so this scan deletes none -- read as
                    # naming nothing, the orphan pass would delete every old body such a
                    # live entry names (its handle then answered as expired). The rule of a
                    # listing that fails, per file
                    return None
                except ValueError:
                    continue              # no entry (not JSON, not UTF-8): no store reads it
                d = j.get('digest') if isinstance(j, dict) else None
                if not isinstance(d, str):
                    continue
                got = (ino, d)
            seen[name] = got
            named.add(got[1])
        # the cache keeps the files listed now only: it never outgrows the directory
        self._foreign = seen
        for h in [h for h in self._entries if h not in mine]:
            self._forget(h)
        return named

    def _gc_bodies(self, digests, *, orphans=False):
        """Delete each body of |digests| that no entry file names any more -- of this
        process's index or of another process's, read from the spool -- and with |orphans|
        every body no entry file names that is older than _ORPHAN_GRACE_S (and every
        temporary file left there as long). Under the spool lock: the scan of the entry
        files and the unlinks happen with no put between a body and its entry."""
        with self._lock:
            left = [d for d in dict.fromkeys(digests) if not self._held_here(d)]
            if not left and not orphans:
                return
            named = self._named_on_disk()
            if named is None:
                return                    # the entries cannot be read: delete nothing
            for d in left:
                if d not in named:
                    self._unlink_body(d)
            if orphans:
                self._sweep_orphans(named)

    def _sweep_orphans(self, named):
        bodies = os.path.join(self.spool_dir, 'bodies')
        cutoff = time.time() - _ORPHAN_GRACE_S        # file times are the wall clock's
        try:
            with os.scandir(bodies) as it:
                listing = [d.name for d in it]
        except OSError:
            return
        for name in listing:
            if name.endswith('.mgt') and name[:-len('.mgt')] in named:
                continue
            if not (name.endswith('.mgt') or name.startswith('.tmp-')):
                continue
            path = os.path.join(bodies, name)
            try:
                if os.stat(path).st_mtime >= cutoff:
                    continue
                os.unlink(path)
            except OSError:
                continue
            if name.endswith('.mgt'):
                self._verified.pop(name[:-len('.mgt')], None)

    def _unlink_body(self, digest):
        try:
            os.unlink(self._body_path(digest))
        except OSError:
            pass
        self._verified.pop(digest, None)

    def _read_body(self, handle, e):
        """The bytes of entry |e|'s stored body, checked against its digest: a body that is
        gone, or whose content is not the digest's (a crash left it empty or torn),
        expires the entry -- UnknownHandle, replayable when it kept a request -- and a
        damaged one is deleted, so that the next put of the same body writes it again
        (kept, a damaged body would raise a format error on every use, the entry would stay
        listed, and a replay or a fresh fetch would point its new entry at the same
        file)."""
        path = self._body_path(e.digest)
        try:
            # read at once without a buffer: a buffered file holds io.DEFAULT_BUFFER_SIZE
            # (128 KiB from Python 3.14) beside the bytes, which no account of
            # standalone_text() charges (graphlet_export(format=mgt) would peak above its
            # account on 27 of 34 retrievals)
            with open(path, 'rb', buffering=0) as f:
                raw = f.read()
        except FileNotFoundError:
            # the body is gone (removed outside this store): the entry is expired, so the
            # handle answers as an expired one -- replayable from its request -- instead
            # of failing on every use while it stays listed (a copy without a parse too)
            self._expire(handle)
            raise self._unknown(handle) from None
        if hashlib.sha256(raw).hexdigest() == e.digest:
            return raw
        with self._lock:
            # the file itself, under the lock: another put may have just replaced it with
            # the right content (then it is read again), else it is the damaged one
            if _file_sha256(path) == e.digest:
                with open(path, 'rb', buffering=0) as f:
                    raw = f.read()
                if hashlib.sha256(raw).hexdigest() == e.digest:
                    return raw
            self._unlink_body(e.digest)
            self._expire(handle)
        raise self._unknown(handle)

    def _expire(self, handle, *, gc=True):
        """Expire |handle|: a tombstone keeps its request, its entry file goes, and with |gc|
        its body when no entry names it any more (a sweep collects the bodies of all it
        expired at once)."""
        with self._lock:
            e = self._unindex(handle)
            if e is None:
                return None               # dropped meanwhile (another thread's GC)
            self._drop_ram(handle)
            self._persisted.pop(handle, None)
            tomb = {'handle': handle, 'request': e.request, 'index': e.index,
                    'expired': self.clock()}
            _write_atomic(self._tomb_path(handle), json.dumps(tomb).encode('utf-8'))
            try:
                os.unlink(self._entry_path(handle))
            except OSError:
                pass
            if gc:
                self._gc_bodies((e.digest,))
        return e

    def sweep(self):
        """Apply the TTLs now: -> {ram_dropped, expired, tombstones_removed}. A tombstone
        (an expired entry's request, for replay) goes ttl_tomb_s after the expiry, so the
        spool does not grow without bound. The bodies the expired entries named go when no
        entry file names them any more -- one scan of the entry files for all of them, not
        one per expired entry -- and so does every body no entry file names that is older
        than _ORPHAN_GRACE_S (otherwise such a body would be kept for good)."""
        now = self.clock()
        ram = [h for h, (_, _, used) in self._ram.items() if now - used > self.ttl_ram_s]
        for h in ram:
            self._drop_ram(h)
        expired = [h for h, e in self._entries.items() if now - e.accessed > self.ttl_disk_s
                   and not self._fresher_on_disk(e, now)]
        with self._lock:
            digests = []
            for h in expired:
                e = self._expire(h, gc=False)
                if e is not None:
                    digests.append(e.digest)
            self._gc_bodies(digests, orphans=True)
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
        """The spool's entries, oldest first: with those another process sharing the
        spool stored, without those another process removed (_live())."""
        try:
            names = {n[:-len('.json')] for n in os.listdir(os.path.join(self.spool_dir,
                                                                        'entries'))
                     if n.endswith('.json')}
        except OSError:
            names = None
        if names is not None:
            for h in [h for h in self._entries if h not in names]:
                self._forget(h)
            for h in sorted(names.difference(self._entries)):
                self._adopt(h)
        return [e.meta() for e in sorted(self._entries.values(), key=lambda x: x.created)]

    def save(self, handle, path):
        """The standalone .mgt file (H, J, body) of an entry."""
        g = self.graphlet(handle)
        from .parser import save
        return save(g, path)

    def body_text(self, handle):
        """The stored body, as received (no parse), checked against its digest: a body
        that is gone or damaged expires the entry (_read_body())."""
        e = self.get(handle)
        # decoded as the text file did with newline='': the same text
        return self._read_body(handle, e).decode('utf-8')

    def standalone_text(self, handle, *, budget=None):
        """The entry's standalone .mgt text -- H, the entry's J line, the body -- built
        WITHOUT a parse (local limits: the escape hatch of an entry no parse can afford): byte
        for byte what save() writes for a canonical body (every body the store holds
        from a server, a save or a view is one), the delivery the entry records written
        into the O record as a parsed model would. The text of an unparsed entry
        (Entry.parsed false) was checked for its frame only (transport, H, Z), never
        validated record by record: a parse of it may still refuse it. budget=: the copy
        is charged (a call of its own: with its base, CALL_BASE)."""
        return _B.run_scoped(budget, 'standalone_text', (),
                             lambda b: self._standalone_text(handle, b))

    def _standalone_text(self, handle, b):
        e = self.get(handle)
        has_env = bool(e.seed_summary) or bool(e.envelope)
        if b is not None:
            # the body read, the slices of it the text is joined from, and the text: three
            # bodies at once, and the body's digest check (sha256 at about 0.8 lwu per 256
            # bytes on the reference machine)
            b.charge(2 * (e.bytes >> 8), 3 * (e.bytes + 49))
            if has_env or e.view is not None or e.derived_from is not None:
                # the J line, not in e.bytes: its text, its copy in the answer and the JSON
                # encoder's pieces (a list of chunks up to Python 3.11), by the bound
                # _j_bound() takes from the entry -- uncharged, an export of a small
                # retrieval peaked above its account under 3.11
                jb = e.j_bound
                b.charge(jb >> 8, 4 * (jb + 49))
        body = self.body_text(handle)
        head = _check_frame(body)
        first = body.find('\n')
        if body.startswith('J ', first + 1):
            raise ValueError('the stored body carries a J line')
        # the O record, rewritten in place when the entry records another delivery: the
        # text is one join of slices of the body around it, never a copy of the body
        # after its H line first
        o_at = o_end = new_o = None
        if e.delivery is not None:
            o_at, o_end, new_o = _o_with_delivery(body, first + 1, e.delivery)
        last = body.rfind('\n', 0, len(body) - 1) + 1
        want_j = has_env or e.view is not None or e.derived_from is not None
        if not want_j:
            if new_o is None:
                return body
            return ''.join((body[:o_at], new_o, body[o_end:]))
        lines = parse_int(body[last:-1].split(' ')[1])
        j = seed_envelope(e.envelope, parse_int(head[11]))
        j['results'] = [{k: v for k, v in (e.seed_summary or {}).items()
                         if k not in ('graphlet', 'graphlet_bytes', 'graphlet_lines')}]
        if e.view is not None:
            j['view'] = e.view
        if e.derived_from is not None:
            j['derived_from'] = e.derived_from
        jt = json.dumps(j, sort_keys=True, separators=(',', ':'), ensure_ascii=False)
        pieces = [body[:first + 1], 'J ', jt, '\n']
        if new_o is None:
            pieces.append(body[first + 1:last])
        else:
            pieces += [body[first + 1:o_at], new_o, body[o_end:last]]
        pieces.append('Z %d\n' % (lines + 1))
        return ''.join(pieces)

    def save_body(self, handle, path, *, budget=None):
        """standalone_text() written atomically to |path| (no parse) -> bytes written."""
        b = _B.resolve(budget)
        if b is None:
            data = self._standalone_text(handle, None).encode('utf-8')
            _write_atomic(os.path.abspath(path), data)
            return len(data)
        with b.scope('save_body'):
            text = self._standalone_text(handle, b)
            # its UTF-8 bytes, made beside it, as graphlet_export(format=mgt) charges them
            # (mcp_tools._write_text_charged()): left out, graphlet_save would peak above
            # its account on 28 of 34 retrievals
            b.charge(0, 2 * (len(text) + 49))
            data = text.encode('utf-8')
            _write_atomic(os.path.abspath(path), data)
            return len(data)

    def load(self, path, request=None, *, parse_budget=None, graph=None, graph_path=None):
        """A saved .mgt file -> a new handle (its J line restores envelope and view).
        Under parse limits (or |parse_budget|) a parse that stops raises and stores
        nothing: the file is left as it is. |graph|, |graph_path|: the graph of a
        multi-graph server the graphlet was retrieved from, kept in the replay request
        (request_of())."""
        return self.load_graphlet(path, request, parse_budget=parse_budget, graph=graph,
                                  graph_path=graph_path)[0]

    def load_graphlet(self, path, request=None, *, parse_budget=None, graph=None,
                      graph_path=None):
        """load() -> (the new handle, the graphlet the load parsed). The model is the one
        a resident entry holds; returned also when the store does not keep it in RAM
        (max_ram_mb), so that a caller needs no second parse of what it just loaded -- a
        parse on demand runs under the store's parse limits, which may refuse what the
        load's own budget admitted (graphlet_load would lose the handle of an entry it had
        stored)."""
        with open(path, 'r', encoding='utf-8', newline='') as f:
            text = f.read()
        g = self._parsing(lambda b: parse(text, budget=b), self._parse_budget(parse_budget))
        h = self.put_graphlet(g, request if request is not None else
                              self.request_of(g, graph=graph, graph_path=graph_path))
        return h, g

    @staticmethod
    def request_of(g, *, graph=None, graph_path=None):
        """The replayable request of a graphlet with an envelope: its seed as validated
        (sequence, the seed labels by name) and the normalized strategy. None when there
        is none -- also when the seed labels' names cannot be verified to resolve back
        to them (ops.resubmittable_names): a replay must not run under other labels.

        |graph|, |graph_path| (else the envelope's, where the server echoed them): the
        graph a multi-graph server ran it on, kept so that a replay runs on that graph
        (without them a multi-graph server refuses the replay, or a client's other graph
        answers). A request without them is valid on a
        single-graph server only."""
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
        graph = graph or g.envelope.get('graph')
        graph_path = graph_path or g.envelope.get('graph_path')
        if graph:
            req['graph'] = graph
        if graph_path:
            req['graph_path'] = graph_path
        return req

    def __contains__(self, handle):
        return self._live(handle) is not None

    def j_of(self, handle):
        return j_object(self.graphlet(handle))
