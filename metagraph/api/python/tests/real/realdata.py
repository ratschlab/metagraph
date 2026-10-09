"""Real-index harness for the metagraph.traverse test suite (api/python/tests/real/).

Shared by the real-data test modules in this directory; never imported by the unit
tests in api/python/tests/. Stdlib + the package only. Importing it has no side
effects: no request is sent until a function asks for one.

Servers (started by the harness, left running; METAGRAPH_REAL_SERVERS overrides):

    sra          127.0.0.1:5701  primary, 42 G nodes, RowDiff<BRWT>, labels = SRA runs (columns)
    uhgg         127.0.0.1:5702  basic, 9.7 G nodes, 4,644 genome labels (columns)
    mini_refseq  127.0.0.1:5703  basic, k-mer coordinates, header (accession) labels via
                                 the .seqs sidecar, --index-manifest: index_fp is set

    METAGRAPH_REAL_SERVERS="sra=host:5701,uhgg=host:5702"   (names not given keep their
    default; NAME=off disables one; a JSON object {"sra": "host:port"} works too)
    METAGRAPH_REAL_SCRATCH=<dir>   the scratch root (catalog, cache, oracle cache)

What lives here:

  * endpoints and SKIP helpers: a test skips when its server is unreachable (or its
    cache cell is missing), so the suite is safe to run anywhere;
  * the SEED CATALOG (<scratch>/catalog.json): seeds per index chosen by real /resolve
    calls (fully in the graph; labels = the profile's full-length carriers);
  * the STRATEGY MATRIX (STRATEGIES): named strategies with per-index parameters, and
    matrix(index) = the (seed, strategy) cells;
  * the RETRIEVAL CACHE (<scratch>/cache/<index>/<cell>.{request,graphlet,full}.json.gz
    + cache/manifest.json): every cell fetched once with detail "graphlet" AND detail
    "full" (output.timing false), so tests replay without refetching;
  * the ORACLE CLIENT: oracle_present / oracle_label_runs from /resolve (the search
    path, which shares no code with the walker), cached under <scratch>/oracle_cache;
  * independent spellers: spell_with_seed(graphlet, arm, leaf) over the parsed G records
    and spell_full_with_seed(full_result, arm, path_id) over the detail-full JSON;
  * the RECORD ORACLE for record coordinates (oracle_positions): the bases at reported
    intervals and the occurrence counts of strings, read from the source FASTA of
    mini_refseq (build/mini_refseq/fasta) -- the index's own input, not the server --
    cached like the /resolve answers (so the offline fixtures replay them).

CLI (run from anywhere):

    python3 realdata.py status
    python3 realdata.py catalog [--index sra ...] [--refresh]
    python3 realdata.py fill    [--index sra ...] [--refresh] [--max-seconds 120] [--workers 2]
    python3 realdata.py report
"""

import argparse
import copy
import functools
import gzip
import hashlib
import json
import os
import random
import socket
import sys
import threading
import time
import unittest
import urllib.error
import urllib.request
import zlib
from collections import OrderedDict
from dataclasses import dataclass
from typing import Callable, Optional

HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.dirname(HERE)
API = os.path.dirname(TESTS)                       # api/python
REPO = os.path.dirname(os.path.dirname(API))       # metagraph/ (the C++ tree)
if API not in sys.path:
    sys.path.insert(0, API)

# ------------------------------------------------------------------ locations

# outside the repo (the cache is ~100 MB); kept across runs so the cells are fetched once
DEFAULT_SCRATCH = os.path.join(os.path.expanduser('~'), '.cache', 'metagraph-traverse-real')
SCRATCH = os.environ.get('METAGRAPH_REAL_SCRATCH') or DEFAULT_SCRATCH
SERVERS_DIR = os.path.join(SCRATCH, 'servers')
CAPABILITIES_PATH = os.path.join(SERVERS_DIR, 'capabilities.json')
CACHE_DIR = os.path.join(SCRATCH, 'cache')
MANIFEST_PATH = os.path.join(CACHE_DIR, 'manifest.json')
ORACLE_DIR = os.path.join(SCRATCH, 'oracle_cache')
CATALOG_PATH = os.path.join(SCRATCH, 'catalog.json')

INDEXES = ('sra', 'uhgg', 'mini_refseq')
DEFAULT_SERVERS = OrderedDict([
    ('sra', ('127.0.0.1', 5701)),
    ('uhgg', ('127.0.0.1', 5702)),
    ('mini_refseq', ('127.0.0.1', 5703)),
])

HOME = os.path.expanduser('~')
MINI_REFSEQ = os.path.join(REPO, 'build', 'mini_refseq')
# the source records of the record oracle, per index: one FASTA per column (file stem =
# the column name, records in file order), each record's header = its accession. An index
# without one has no record oracle (its coordinates are null anyway: no k-mer coordinates)
RECORDS_DIR = {'mini_refseq': os.path.join(MINI_REFSEQ, 'fasta')}
INDEX_FILES = {
    'sra': {'graph': os.path.join(HOME, 'metagraph-indexes/sra/graph_primary_small.dbg'),
            'annotation': os.path.join(
                HOME, 'metagraph-indexes/sra/annotation.relaxed.row_diff_brwt.annodbg')},
    'uhgg': {'graph': os.path.join(HOME, 'metagraph-indexes/uhgg/graph_complete_k31.dbg'),
             'annotation': os.path.join(
                 HOME, 'metagraph-indexes/uhgg/annotation.relaxed.relabeled.row_diff_brwt.annodbg')},
    'mini_refseq': {
        'graph': os.path.join(MINI_REFSEQ, 'graph_k31.dbg'),
        'annotation': os.path.join(
            MINI_REFSEQ, 'annotation.relaxed.relabeled.row_diff_brwt_coord.annodbg'),
        'manifest': os.path.join(MINI_REFSEQ, 'annotation.relaxed.relabeled.manifest.json')},
}
UHGG_DIR = os.path.join(HOME, 'metagraph-indexes', 'uhgg')
QUERY_MATERIAL = {
    'r16s': os.path.join(UHGG_DIR, 'r16s.json'),               # {accession: 16S sequence}
    'genomes': [os.path.join(UHGG_DIR, f)
                for f in ('MGYG000000041.fna', 'mgyg.fna', 'g153.fna')],
    'ndm1': os.path.join(MINI_REFSEQ, 'sanity', 'ndm1.fa'),
    'ndm1_rc': os.path.join(MINI_REFSEQ, 'sanity', 'ndm1_rc.fa'),
    'mini_records': os.path.join(MINI_REFSEQ, 'raw'),
}

PROBE_TIMEOUT = 3.0          # seconds, the reachability probe
HEAVY_FULL_BYTES = 64 << 20  # a cell whose full JSON exceeds this (or whose bodies exceed
HEAVY_BODY_BYTES = 8 << 20   # this) is flagged heavy: cells() leaves it out by default
DEFAULT_TIMEOUT = 300.0      # seconds, a live request from a test
CELL_TIMEOUT = 120.0         # seconds, one request while filling the cache


def _env_servers():
    raw = os.environ.get('METAGRAPH_REAL_SERVERS', '').strip()
    out = OrderedDict(DEFAULT_SERVERS)
    if not raw:
        return out
    if raw.startswith('{'):
        items = json.loads(raw).items()
    else:
        items = [p.split('=', 1) for p in raw.split(',') if p.strip()]
    for name, addr in items:
        name, addr = name.strip(), str(addr).strip()
        if addr.lower() in ('', 'off', 'none', 'disabled'):
            out[name] = None
            continue
        host, _, port = addr.rpartition(':')
        out[name] = (host or '127.0.0.1', int(port))
    return out


SERVERS = _env_servers()


def server(index):
    """(host, port) of |index|'s server; None when disabled by METAGRAPH_REAL_SERVERS."""
    if index not in SERVERS:
        raise KeyError('unknown index %r (known: %s)' % (index, ', '.join(SERVERS)))
    return SERVERS[index]


def endpoint(index):
    hp = server(index)
    return None if hp is None else 'http://%s:%d' % hp


# ------------------------------------------------------------------ http

class RealServerError(RuntimeError):
    """A non-2xx answer from a real server (status, message, body)."""

    def __init__(self, index, route, status, message, body=None):
        self.index, self.route, self.status, self.message, self.body = (
            index, route, status, message, body)
        super().__init__('%s %s: HTTP %s: %s' % (index, route, status, message))


class ServerUnavailable(unittest.SkipTest):
    """The index's server is down, initializing (503) or disabled: tests skip."""


@dataclass
class HttpResult:
    status: int
    json: Optional[object]       # parsed body (None when not JSON)
    body: bytes                  # decoded body bytes (what the server serialised)
    wire_bytes: int              # bytes on the wire (gzip when the server compressed)
    encoding: Optional[str]
    elapsed_s: float

    @property
    def ok(self):
        return 200 <= self.status < 300

    @property
    def error(self):
        return self.json.get('error') if isinstance(self.json, dict) else None


def http(index, method, route, payload=None, *, timeout=DEFAULT_TIMEOUT):
    """One request to |index|'s server with Accept-Encoding: gzip. Returns HttpResult
    for every HTTP status; raises ServerUnavailable when the server is disabled and the
    socket's own errors (URLError, TimeoutError) otherwise."""
    base = endpoint(index)
    if base is None:
        raise ServerUnavailable('%s: disabled by METAGRAPH_REAL_SERVERS' % index)
    url = base + '/' + route.lstrip('/')
    data = None if payload is None else json.dumps(payload, separators=(',', ':')).encode()
    headers = {'Accept-Encoding': 'gzip'}
    if data is not None:
        headers['Content-Type'] = 'application/json'
    req = urllib.request.Request(url, data=data, headers=headers, method=method)
    t0 = time.perf_counter()
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            status, enc, raw = r.status, r.headers.get('Content-Encoding'), r.read()
    except urllib.error.HTTPError as e:
        status = e.code
        enc = e.headers.get('Content-Encoding') if e.headers else None
        try:
            raw = e.read() or b''
        finally:
            e.close()
    elapsed = time.perf_counter() - t0
    enc_l = (enc or '').strip().lower()
    if enc_l == 'gzip':
        body = gzip.decompress(raw)
    elif enc_l == 'deflate':
        body = zlib.decompress(raw)
    else:
        body = raw
    try:
        obj = json.loads(body) if body else None
    except ValueError:
        obj = None
    return HttpResult(status, obj, body, len(raw), enc, elapsed)


def post(index, route, payload, *, timeout=DEFAULT_TIMEOUT):
    """POST and return the JSON; non-2xx -> RealServerError, 503 -> ServerUnavailable."""
    res = http(index, 'POST', route, payload, timeout=timeout)
    if res.status == 503:
        raise ServerUnavailable('%s: initializing (503)' % index)
    if not res.ok:
        raise RealServerError(index, route, res.status, res.error or res.body[:300], res.json)
    return res.json


def client(index, timeout=DEFAULT_TIMEOUT):
    """metagraph.traverse.TraverseClient for |index| (the code under test)."""
    from metagraph.traverse import TraverseClient
    hp = server(index)
    if hp is None:
        raise ServerUnavailable('%s: disabled by METAGRAPH_REAL_SERVERS' % index)
    return TraverseClient(hp[0], hp[1], timeout=timeout)


# ------------------------------------------------------------------ reachability / skips

_UP = {}
_UP_LOCK = threading.Lock()


def server_up(index, *, refresh=False):
    """True iff GET /traverse/capabilities answers 200 (cached per process)."""
    with _UP_LOCK:
        if not refresh and index in _UP:
            return _UP[index][0]
    caps = None
    try:
        res = http(index, 'GET', 'traverse/capabilities', timeout=PROBE_TIMEOUT)
        up = res.status == 200 and isinstance(res.json, dict)
        caps = res.json if up else None
    except (ServerUnavailable, OSError, urllib.error.URLError, socket.timeout, ValueError):
        up = False
    with _UP_LOCK:
        _UP[index] = (up, caps)
    return up


def live_capabilities(index):
    """The live server's capabilities; raises ServerUnavailable when it is down."""
    if not server_up(index):
        require_server(index)
    with _UP_LOCK:
        return copy.deepcopy(_UP[index][1])


def recorded_capabilities(index=None):
    """capabilities.json as written when the servers were started (offline-safe):
    {index: {host, port, pid, k, regime, index_ns, index_fp, index_meta_fp,
    num_labels, ..., capabilities}}; one index's entry when |index| is given."""
    try:
        with open(CAPABILITIES_PATH) as f:
            caps = json.load(f)
    except OSError:
        caps = {}
    return caps.get(index) if index is not None else caps


def identity(index):
    """{index_ns, index_fp, index_meta_fp, k, regime} -- live when the server is up,
    else as recorded."""
    src = None
    if server_up(index):
        src = live_capabilities(index)
    else:
        rec = recorded_capabilities(index)
        src = rec.get('capabilities') if rec else None
    if not src:
        return None
    return {k: src.get(k) for k in ('index_ns', 'index_fp', 'index_meta_fp', 'k', 'regime',
                                    'alphabet', 'num_labels')}


def live_identity(index):
    """How the live server's index relates to the one the cache was filled from ->
    (status, reason), status one of
      'same'          both state an index manifest digest (index_fp) and they are equal;
      'unverifiable'  index_meta_fp agrees, but at most one side states index_fp: equal
                      metadata is a negative check only (two_sided False: neither has a
                      manifest; True: one side has one and the other not);
      'different'     index_meta_fp, index_fp (where both state it) or index_ns differ;
      'down'          no recorded entry, or the server does not answer.
    Every fingerprint both sides state is compared: equal index_meta_fp alone would accept a
    live server whose manifest digest differs from the recorded one."""
    m = load_manifest().get('servers', {}).get(index)
    if not m or not server_up(index):
        return 'down', 'no recorded identity or no live server'
    live = live_capabilities(index) or {}
    for key in ('index_meta_fp', 'index_fp', 'index_ns'):
        a, b = m.get(key) or None, live.get(key) or None
        if a and b and a != b:
            return 'different', '%s differs (recorded %s, live %s)' % (key, a, b)
    fp_rec, fp_live = m.get('index_fp') or None, live.get('index_fp') or None
    if fp_rec and fp_live:
        return 'same', 'equal index manifest digests'
    if not (m.get('index_meta_fp') and live.get('index_meta_fp')):
        return 'unverifiable', ('index_meta_fp is not stated on both sides and no manifest '
                                'digest on both: nothing identifies the index')
    if fp_rec or fp_live:
        return 'unverifiable', ('the %s states an index manifest digest and the %s none: '
                                'equal index_meta_fp is a negative check only'
                                % (('cache', 'live server') if fp_rec
                                   else ('live server', 'cache')))
    return 'unverifiable', ('neither the cache nor the live server states an index '
                            'manifest digest: equal index_meta_fp is a negative check only')


def live_matches_cache(index):
    """True iff the live server may be taken to serve the index the cache was filled
    from: equal manifest digests ('same'), or -- where neither side has a manifest, as
    for sra and uhgg -- equal index_meta_fp, the only check there is (a negative one, but
    enough to catch a wrong port). A digest on one side only cannot be verified and a
    differing fingerprint proves another index: False for both (live_identity says
    which)."""
    status, _ = live_identity(index)
    if status == 'same':
        return True
    if status != 'unverifiable':
        return False
    m = load_manifest().get('servers', {}).get(index) or {}
    live = live_capabilities(index) or {}
    return not (m.get('index_fp') or live.get('index_fp')) \
        and bool(m.get('index_meta_fp')) and m.get('index_meta_fp') == live.get('index_meta_fp')


def require_server(index, testcase=None):
    """Raise unittest.SkipTest (ServerUnavailable) unless |index|'s server answers."""
    hp = server(index)
    if hp is None:
        raise ServerUnavailable('%s: disabled by METAGRAPH_REAL_SERVERS' % index)
    if not server_up(index):
        raise ServerUnavailable('%s server at %s:%d is unreachable or initializing'
                                % (index, hp[0], hp[1]))


def require_cache(index=None, cell=None):
    """Raise unittest.SkipTest unless the cache holds |cell| of |index| (or, without a
    cell, at least one ok cell of |index| / of any index)."""
    if cell is not None:
        meta = cell_meta(index, cell)
        if meta is None or meta.get('status') != 'ok':
            raise unittest.SkipTest('cache cell %s/%s missing (run realdata.py fill)'
                                    % (index, cell))
        return meta
    if not cells(index):
        raise unittest.SkipTest('no cached cells%s (run realdata.py fill)'
                                % (' for ' + index if index else ''))


def require_catalog(index=None):
    cat = load_catalog()
    if not cat.get('indexes') or (index and not cat['indexes'].get(index)):
        raise unittest.SkipTest('no seed catalog%s (run realdata.py catalog)'
                                % (' for ' + index if index else ''))
    return cat


def _skip_decorator(check):
    def deco(obj):
        if isinstance(obj, type):
            orig = obj.__dict__.get('setUpClass')

            def setUpClass(cls):
                check()
                if orig is not None:
                    orig.__func__(cls)
                else:
                    super(obj, cls).setUpClass()
            obj.setUpClass = classmethod(setUpClass)
            return obj

        @functools.wraps(obj)
        def wrapper(*a, **kw):
            check()
            return obj(*a, **kw)
        return wrapper
    return deco


def skip_unless_server(*indexes):
    """Class or method decorator: skip (at run time, not import time) unless every
    named server answers."""
    def check():
        for i in indexes:
            require_server(i)
    return _skip_decorator(check)


def skip_unless_cache(index=None, cell=None):
    """Class or method decorator: skip unless the cache holds the cell(s)."""
    return _skip_decorator(lambda: require_cache(index, cell))


def skip_unless_catalog(index=None):
    return _skip_decorator(lambda: require_catalog(index))


# ------------------------------------------------------------------ small io helpers

def _write_json_gz(path, obj_or_bytes):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    data = obj_or_bytes if isinstance(obj_or_bytes, (bytes, bytearray)) else \
        json.dumps(obj_or_bytes, separators=(',', ':')).encode()
    tmp = path + '.tmp.%d.%d' % (os.getpid(), threading.get_ident())
    with gzip.open(tmp, 'wb', compresslevel=6) as f:
        f.write(data)
    os.replace(tmp, path)


def _read_gz_bytes(path):
    with gzip.open(path, 'rb') as f:
        return f.read()


def _read_json_gz(path):
    return json.loads(_read_gz_bytes(path))


def _write_json(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + '.tmp.%d.%d' % (os.getpid(), threading.get_ident())
    with open(tmp, 'w') as f:
        json.dump(obj, f, indent=1, sort_keys=False)
    os.replace(tmp, path)


def read_fasta(path):
    """[(header, sequence)] (header without '>', sequence upper case)."""
    out, name, parts = [], None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if name is not None:
                    out.append((name, ''.join(parts).upper()))
                name, parts = line[1:].split()[0] if line[1:].split() else '', []
            elif line:
                parts.append(line)
    if name is not None:
        out.append((name, ''.join(parts).upper()))
    return out


_RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def revcomp(s):
    return s.translate(_RC)[::-1]


# ------------------------------------------------------------------ oracle client

class OracleError(RuntimeError):
    pass


def _label_name(x):
    return x.name if hasattr(x, 'name') else str(x)


def _oracle_key(index, payload):
    ident = identity(index) or {}
    blob = json.dumps({'index': index, 'meta_fp': ident.get('index_meta_fp'),
                       'fp': ident.get('index_fp'), 'payload': payload},
                      sort_keys=True, separators=(',', ':')).encode()
    return hashlib.sha256(blob).hexdigest()


def oracle_resolve(index, sequence, labels=None, *, support='kmer', discover=None):
    """POST /resolve (cached on disk under oracle_cache/<index>/). Exactly one of
    |labels| (names or Label objects) and |discover| ({max_labels, kind}); neither ->
    discover {max_labels: 1} (the cheapest way to get graph_runs). The cache key holds
    the index identity, so a different index behind the same name never hits."""
    sequence = sequence.upper()
    payload = {'sequence': sequence}
    if labels is not None:
        names = sorted({_label_name(l) for l in labels})
        if not names:
            raise ValueError('oracle_resolve: an empty label list')
        payload['labels'] = names
    else:
        payload['discover'] = dict(discover) if discover else {'max_labels': 1}
    if support != 'kmer':
        payload['support'] = support
    key = _oracle_key(index, payload)
    path = os.path.join(ORACLE_DIR, index, key[:2], key + '.json.gz')
    if os.path.exists(path):
        try:
            return _read_json_gz(path)['response']
        except (OSError, ValueError, KeyError):
            pass
    require_server(index)
    res = http(index, 'POST', 'resolve', payload)
    if res.status == 503:
        raise ServerUnavailable('%s: initializing (503)' % index)
    if not res.ok:
        raise OracleError('%s /resolve HTTP %s: %s' % (index, res.status, res.error))
    out = {k: v for k, v in res.json.items() if k not in ('capabilities', 'timing')}
    _write_json_gz(path, {'request': payload, 'response': out})
    return out


def _k(index):
    ident = identity(index)
    return int(ident['k']) if ident and ident.get('k') else 31


_RECORDS_MEMO = {}


def _source_records(index):
    """{('header', accession): [(0, seq)], ('column', name): [(offset, seq), ...]} of the
    index's source FASTA, or None without one. A column's offsets are its k-mer index
    space (DESIGN §18): record i's k-mer j is offset_i + j, offset_{i+1} = offset_i +
    len_i - k + 1, so record i's bases are numbered [offset_i, offset_i + len_i)."""
    d = RECORDS_DIR.get(index)
    if not d or not os.path.isdir(d):
        return None
    got = _RECORDS_MEMO.get(d)
    if got is None:
        k = _k(index)
        got = {}
        for fn in sorted(os.listdir(d)):
            if not fn.endswith('.fa'):
                continue
            col, off = [], 0
            for h, s in read_fasta(os.path.join(d, fn)):
                s = s.upper()
                got[('header', h)] = [(0, s)]
                col.append((off, s))
                off += max(0, len(s) - k + 1)
            got[('column', fn[:-3])] = col
        _RECORDS_MEMO[d] = got
    return got


def _count(hay, needle):
    """Occurrences of |needle| in |hay|, overlapping ones too (two coordinate chains of
    a repeat can overlap)."""
    n, i = 0, hay.find(needle)
    while i >= 0:
        n += 1
        i = hay.find(needle, i + 1)
    return n


def oracle_positions(index, kind, name, slices, needles):
    """The record oracle (cached like oracle_resolve, keyed by the index identity and
    the payload): for the label |name| of |kind| ('header' | 'column'), per [start, end)
    of |slices| the sha256 of the bases there in every source record that holds the whole
    interval (a header label: its record; a column label: its records in the column's
    k-mer index space, where a record's last k-1 bases share their numbers with the next
    record's first: up to two), and per string of |needles| its occurrence count in the
    label's record(s). Read from the source FASTA (no server); without it, a miss is the
    harness's require_server (offline: a failure, the fixtures are stale)."""
    payload = {'positions': {'kind': kind, 'name': name,
                             'slices': [list(x) for x in slices]},
               # the needles as the payload's 'sequence': the recorded answer keeps only its
               # length and digest (offline._slim), the key binds the strings
               'sequence': '\n'.join(needles)}
    key = _oracle_key(index, payload)
    path = os.path.join(ORACLE_DIR, index, key[:2], key + '.json.gz')
    if os.path.exists(path):
        try:
            return _read_json_gz(path)['response']
        except (OSError, ValueError, KeyError):
            pass
    recs = _source_records(index)
    if recs is None:
        require_server(index)
        raise unittest.SkipTest('%s: no source records for the record oracle (%s)'
                                % (index, RECORDS_DIR.get(index)))
    target = recs.get((kind, name))
    if target is None:
        raise OracleError('%s: no source record for %s label %r' % (index, kind, name))
    out = {'slices': [], 'counts': []}
    for s, e in slices:
        hits = []
        for off, seq in target:
            if off <= s and e <= off + len(seq):
                hits.append(hashlib.sha256(seq[s - off:e - off].encode()).hexdigest())
        out['slices'].append(hits)
    for needle in needles:
        out['counts'].append(sum(_count(seq, needle) for _, seq in target))
    _write_json_gz(path, {'request': payload, 'response': out})
    return out


def oracle_graph_runs(index, sequence):
    """The k-mer runs [a, b) of |sequence| that are in |index|'s graph (resolve's
    graph_runs; strand-agnostic in the primary regime, forward strand in basic)."""
    if len(sequence) < _k(index):
        raise ValueError('a sequence shorter than k=%d has no k-mers' % _k(index))
    out = oracle_resolve(index, sequence)
    return [tuple(r) for r in out.get('graph_runs') or []]


def oracle_present(index, sequence):
    """True iff /resolve reports every k-mer of |sequence| in the graph. A sequence
    shorter than k raises (no k-mer: the claim would be vacuous)."""
    if len(sequence) < _k(index):
        raise ValueError('a sequence shorter than k=%d has no k-mers' % _k(index))
    out = oracle_resolve(index, sequence)
    n = out['num_kmers']
    return [list(r) for r in out.get('graph_runs') or []] == [[0, n]]


def oracle_label_runs(index, sequence, labels, *, support='kmer', chunk=1000):
    """Per label name: the k-mer runs [(a, b), ...] of |sequence| on which /resolve
    reports the label (k-mer start intervals, half-open; base interval [a, b + k - 1)).
    |labels| are names or Label objects; a label the server reports without runs maps
    to []. support='trace' asks for coordinate-consecutive runs (header labels,
    mini_refseq). Chunked by |chunk| labels per request."""
    names = sorted({_label_name(l) for l in labels})
    out = {n: [] for n in names}
    for i in range(0, len(names), chunk):
        res = oracle_resolve(index, sequence, names[i:i + chunk], support=support)
        for l in res.get('labels') or []:
            if l['label'] in out:
                out[l['label']] = [tuple(r) for r in l.get('runs') or []]
    return out


def oracle_num_kmers(index, sequence):
    return max(0, len(sequence) - _k(index) + 1)


def oracle_label_covers(index, sequence, label, *, support='kmer'):
    """True iff |label| supports every k-mer of |sequence| (one run [0, n))."""
    runs = oracle_label_runs(index, sequence, [label], support=support)[_label_name(label)]
    return runs == [(0, oracle_num_kmers(index, sequence))]


def oracle_profile(index, sequence, *, max_labels=5000, kind=None):
    """The discovery profile: resolve with discover (every label touching the
    sequence, top max_labels by support)."""
    disc = {'max_labels': max_labels}
    if kind:
        disc['kind'] = kind
    return oracle_resolve(index, sequence, discover=disc)


def full_carriers(resolve_out):
    """Names of the labels whose runs cover every k-mer (one run [0, num_kmers))."""
    n = resolve_out['num_kmers']
    return [l['label'] for l in resolve_out.get('labels') or [] if l.get('runs') == [[0, n]]]


# ------------------------------------------------------------------ independent spellers

def _side(arm):
    s = getattr(arm, 'side', arm)
    return {'l': 'left', 'r': 'right'}.get(s, s)


def _leaf_id(leaf):
    if isinstance(leaf, bool):
        raise TypeError('a leaf is a segment id, a Walk or a Path')
    if hasattr(leaf, 'leaf'):               # Walk / Path
        return leaf.leaf
    return int(leaf)


def walk_order_flank(graphlet, arm, segment):
    """Outward bases root -> |segment| along the first-parent chain (walking order),
    from the parsed G records only (no derive/ops code)."""
    a = graphlet.arms[_side(arm)]
    sid = _leaf_id(segment)
    chain = []
    while True:
        chain.append(sid)
        parents = a.segments[sid].parents
        if not parents:
            break
        sid = parents[0]
    parts = []
    for s in reversed(chain):
        w = a.segments[s].walk
        if w is None:
            raise ValueError('this retrieval carries no bases (output.sequences: false)')
        parts.append(w)
    return ''.join(parts)


def spell_with_seed(graphlet, arm, leaf):
    """The whole molecule in natural orientation (§2.4), independent of the library's
    spelling code: right arm seed + flank, left arm reverse(walking-order flank) + seed.
    |leaf| is a leaf SEGMENT id, or a Walk / Path (their .leaf); any segment id works
    (the molecule up to that segment's end)."""
    flank = walk_order_flank(graphlet, arm, leaf)
    seed = graphlet.seed.sequence
    return seed + flank if _side(arm) == 'right' else flank[::-1] + seed


def leaf_segments(graphlet, arm):
    """Leaf segment ids in segment-id order (path id i <-> the i-th leaf, §2.5)."""
    return [s.id for s in graphlet.arms[_side(arm)].segments if s.leaf is not None]


def spell_path_with_seed(graphlet, arm, path_id):
    return spell_with_seed(graphlet, arm, leaf_segments(graphlet, arm)[path_id])


def spell_full_with_seed(full_result, arm, path_id, seed_sequence=None):
    """The same molecule from a detail-full results[i] (natural-orientation segment
    sequences; the left flank is the concat leaf -> root). |seed_sequence| defaults to
    the graphlet-free source: full_result has no seed bases, so pass the seed (e.g.
    CachedCell.seed_sequence(i))."""
    side = _side(arm)
    a = full_result['arms'][side]
    p = a['paths'][path_id]
    segs = a['segments']
    chain = p['segments']
    if side == 'right':
        flank = ''.join(segs[s]['sequence'] for s in chain)
    else:
        flank = ''.join(segs[s]['sequence'] for s in reversed(chain))
    if seed_sequence is None:
        raise ValueError('spell_full_with_seed needs the seed sequence')
    return seed_sequence + flank if side == 'right' else flank + seed_sequence


def leaf_molecules(graphlet, arm):
    """{leaf segment id: natural molecule with the seed} for every leaf of |arm|."""
    return {s: spell_with_seed(graphlet, arm, s) for s in leaf_segments(graphlet, arm)}


# ------------------------------------------------------------------ full-JSON comparison

def _canon(x):
    return json.dumps(x, sort_keys=True, separators=(',', ':'))


def normalize_full_result(result):
    """The §2.5 'not information' normalisation of a detail-full results[i], written
    independently of metagraph.traverse.export.normalize_result: a deep copy without
    'timing'; per arm, each segment's events ordered by (at_bp, canonical JSON) --
    order among equal at_bp is not information --, needed_budgets sorted, and
    growth[].label_ends null -> {}."""
    r = json.loads(_canon(result))
    r.pop('timing', None)
    for arm in (r.get('arms') or {}).values():
        for seg in arm.get('segments') or []:
            if 'events' in seg:
                seg['events'] = sorted(seg['events'], key=lambda e: (e.get('at_bp'), _canon(e)))
        if 'needed_budgets' in arm:
            arm['needed_budgets'] = sorted(arm['needed_budgets'])
        for b in arm.get('growth') or []:
            if 'label_ends' in b and b['label_ends'] is None:
                b['label_ends'] = {}
    return r


def json_diff(a, b, path='', limit=20):
    """Up to |limit| human-readable differences between two JSON values."""
    out = []

    def walk(x, y, p):
        if len(out) >= limit:
            return
        if isinstance(x, dict) and isinstance(y, dict):
            for k in sorted(set(x) | set(y), key=str):
                if k not in x:
                    out.append('%s/%s: only in b' % (p, k))
                elif k not in y:
                    out.append('%s/%s: only in a' % (p, k))
                else:
                    walk(x[k], y[k], '%s/%s' % (p, k))
        elif isinstance(x, list) and isinstance(y, list):
            if len(x) != len(y):
                out.append('%s: length %d vs %d' % (p, len(x), len(y)))
            for i, (u, v) in enumerate(zip(x, y)):
                walk(u, v, '%s/%d' % (p, i))
        elif x != y or type(x) is not type(y) and not (
                isinstance(x, (int, float)) and isinstance(y, (int, float))):
            out.append('%s: %r vs %r' % (p, x, y))
    walk(a, b, path)
    return out


# ------------------------------------------------------------------ live traversal

def live_traverse(index, seed, strategy_name=None, *, detail='graphlet', timing=False,
                  strategy_dict=None, timeout=DEFAULT_TIMEOUT):
    """POST /traverse live (skips when the server is down): |seed| is a catalog seed
    name, a seed record or a raw sequence; the strategy is a named one or
    |strategy_dict|. Returns (request, response JSON)."""
    require_server(index)
    if isinstance(seed, str) and set(seed.upper()) <= set('ACGTN') and len(seed) >= 2 \
            and not any(s['name'] == seed for s in catalog_seeds(index)):
        seed = {'name': 'adhoc', 'kind': 'adhoc', 'sequence': seed}
    elif isinstance(seed, str):
        seed = catalog_seed(index, seed)
    if strategy_dict is not None:
        st = copy.deepcopy(strategy_dict)
        st.setdefault('output', {}).update({'detail': detail, 'timing': timing})
        req = {'seeds': [{'sequence': s} for s in seed_sequences(seed)], 'strategy': st}
    else:
        req = build_request(index, seed, strategy_name or 'strict', detail, timing)
    return req, post(index, 'traverse', req, timeout=timeout)


# ------------------------------------------------------------------ strategy matrix

@dataclass
class Strategy:
    name: str
    doc: str
    build: Callable                      # (index, seed_record) -> strategy dict
    indexes: Optional[tuple] = None      # None: every index
    bulk: bool = False                   # also on 'bulk' tier seeds
    batch: bool = False                  # also on the multi-seed 'batch' pseudo-seed
    single: bool = True                  # on single seeds
    deterministic: bool = True           # False: a time budget is meant to trip
    tags: tuple = ()
    kinds: Optional[tuple] = None        # only seeds of these kinds (None: every kind)
    seed_tiers: Optional[tuple] = None   # only seeds of these tiers (None: every tier)

    def applies(self, index, seed, tiers=None):
        """The cell (index, seed, this strategy) is in the matrix. The cache holds the
        full product (the fill takes minutes); |tiers| ('core', 'bulk') restricts it:
        a 'bulk' seed then gets only the strategies marked bulk (the cheap ones). A seed of
        a NARROW_KINDS kind gets only the strategies listed for it."""
        if self.indexes is not None and index not in self.indexes:
            return False
        if seed['kind'] in NARROW_KINDS:
            return self.name in NARROW_KINDS[seed['kind']]
        if seed['kind'] == 'batch':
            return self.batch
        if not self.single:
            return False
        if self.kinds is not None and seed['kind'] not in self.kinds:
            return False
        if self.seed_tiers is not None and seed.get('tier') not in self.seed_tiers:
            return False
        if tiers is None:
            return True
        tier = seed.get('tier')
        return tier in tiers and (tier == 'core' or self.bulk)


# seed kinds whose cells are only these strategies: the repeat windows of mini_refseq are
# there for record coordinates (lists of several occurrences, cuts, switches into live
# labels), not for the whole matrix
NARROW_KINDS = {'repeat': ('trace', 'trace_coords', 'trace_coords_column')}

ANNOTATE_RADIUS = {'sra': 100, 'uhgg': 300, 'mini_refseq': 300}
BEAM20_BP = {'sra': 3000, 'uhgg': 5000, 'mini_refseq': 5000}
BEAM1_BP = 10000
TIGHT_MS = {'sra': 200, 'uhgg': 50, 'mini_refseq': 5}


def _beam(width, bp):
    return {'labels': {'mode': 'annotate'},
            'frontier': {'order': 'most_supported_first', 'on_overflow': 'beam'},
            'bounds': {'max_extension_bp': bp, 'max_live_paths': width}}


def _switch(index, seed):
    st = {'labels': {'change_cost': {'model': 'constant', 'value': 1}, 'loss_budget': 1}}
    if seed.get('switch_extra'):
        st['labels']['extra'] = list(seed['switch_extra'])
    return st


STRATEGIES = OrderedDict((s.name, s) for s in [
    Strategy('strict', 'the default strategy: constrain, derived seed labels, forbid, '
             'branch limit 0, merge, radius 5 kb', lambda i, s: {}, bulk=True, batch=True,
             tags=('constrain', 'merge')),
    Strategy('limit1', 'branch limit 1', lambda i, s: {'branching': {'max_label_branches': 1}},
             bulk=True, tags=('constrain', 'merge')),
    Strategy('limit2', 'branch limit 2, on_reconverge merge (explicit)',
             lambda i, s: {'branching': {'max_label_branches': 2, 'on_reconverge': 'merge'}},
             bulk=True, tags=('constrain', 'merge')),
    Strategy('limit2_keep', 'branch limit 2, on_reconverge keep (vs limit2: merge vs keep)',
             lambda i, s: {'branching': {'max_label_branches': 2, 'on_reconverge': 'keep'}},
             tags=('constrain', 'keep')),
    Strategy('exhaustive', 'constrain trie: exhaustive, radius 300, live-path cap 500',
             lambda i, s: {'exhaustive': True,
                           'bounds': {'max_extension_bp': 300, 'max_live_paths': 500}},
             tags=('constrain', 'keep', 'exhaustive')),
    Strategy('annotate_exh', 'label-free trie: annotate exhaustive, radius 100 (SRA) / '
             '300, live-path cap 300',
             lambda i, s: {'exhaustive': True, 'labels': {'mode': 'annotate'},
                           'bounds': {'max_extension_bp': ANNOTATE_RADIUS[i],
                                      'max_live_paths': 300}},
             bulk=True, tags=('annotate', 'keep', 'exhaustive')),
    Strategy('annotate_merge', 'annotate, NOT exhaustive, on_reconverge merge, same radius '
             'and cap as annotate_exh (vs annotate_exh: merge vs keep)',
             lambda i, s: {'labels': {'mode': 'annotate'},
                           'branching': {'on_reconverge': 'merge'},
                           'bounds': {'max_extension_bp': ANNOTATE_RADIUS[i],
                                      'max_live_paths': 300}},
             tags=('annotate', 'merge')),
    Strategy('beam1', 'annotate beam width 1, most_supported_first, 10 kb',
             lambda i, s: _beam(1, BEAM1_BP), bulk=True, batch=True,
             tags=('annotate', 'beam')),
    Strategy('beam20', 'annotate beam width 20, most_supported_first, 3 kb (SRA) / 5 kb',
             lambda i, s: _beam(20, BEAM20_BP[i]), tags=('annotate', 'beam')),
    Strategy('beam20_10k', 'the design §7 reference case: annotate beam width 20, '
             'most_supported_first, 10 kb, on the core random 50-mers of SRA (tens of MB of '
             'full JSON: heavy)', lambda i, s: _beam(20, 10000), indexes=('sra',),
             kinds=('random',), seed_tiers=('core',), tags=('annotate', 'beam', 'scale')),
    Strategy('switch1', 'one switch: change_cost constant 1, loss_budget 1 (+ labels.extra '
             'on a seed that names switch targets: the UHGG contig end)', _switch,
             tags=('constrain', 'switch')),
    Strategy('left_only', 'strict, direction left', lambda i, s: {'direction': 'left'},
             tags=('constrain', 'direction')),
    Strategy('right_only', 'strict, direction right', lambda i, s: {'direction': 'right'},
             tags=('constrain', 'direction')),
    Strategy('no_sequences', 'limit1 with output.sequences false (G without bases, C with '
             'its sequence)', lambda i, s: {'branching': {'max_label_branches': 1},
                                            'output': {'sequences': False}},
             tags=('constrain', 'sequences_false')),
    Strategy('cut_lists', 'annotate exhaustive radius 100, cap 300, max_labels_per_node 2 '
             '(cut lists: label_evidence lower_bound, claims(strict) raises)',
             lambda i, s: {'exhaustive': True,
                           'labels': {'mode': 'annotate', 'max_labels_per_node': 2},
                           'bounds': {'max_extension_bp': 100, 'max_live_paths': 300}},
             tags=('annotate', 'cut_lists')),
    Strategy('cut_events', 'constrain exhaustive radius 150, cap 500, max_branch_events 1 '
             '(cut evidence: branch_diagnostics cut)',
             lambda i, s: {'exhaustive': True,
                           'bounds': {'max_extension_bp': 150, 'max_live_paths': 500},
                           'output': {'max_branch_events': 1}},
             tags=('constrain', 'cut_events')),
    Strategy('trace', 'support trace (coordinate-consecutive header labels), on_reconverge '
             'keep (trace refuses merge); output.coordinates pinned false, so that a '
             'library that asks for record coordinates by itself (traverse_fetch, '
             'feature level 6) fetches what the cache holds',
             lambda i, s: {'support': 'trace', 'branching': {'on_reconverge': 'keep'},
                           'output': {'coordinates': False}},
             indexes=('mini_refseq',), tags=('constrain', 'trace', 'keep')),
    Strategy('trace_coords', 'trace with record coordinates (output.coordinates true, '
             'max_coordinate_occurrences "unlimited"): the positional oracle f_positions '
             'checks every occurrence against the source records (feature level 6)',
             lambda i, s: {'support': 'trace', 'branching': {'on_reconverge': 'keep'},
                           'bounds': {'max_extension_bp': 2000},
                           'output': {'coordinates': True,
                                      'max_coordinate_occurrences': 'unlimited'}},
             indexes=('mini_refseq',), tags=('constrain', 'trace', 'keep', 'coordinates')),
    Strategy('trace_coords_column', 'trace on column labels (seed_label_kind column) with '
             'record coordinates: global column positions (kind column), '
             'trace_record_boundaries',
             lambda i, s: {'support': 'trace', 'branching': {'on_reconverge': 'keep'},
                           'labels': {'seed_label_kind': 'column'},
                           'bounds': {'max_extension_bp': 2000},
                           'output': {'coordinates': True,
                                      'max_coordinate_occurrences': 'unlimited'}},
             indexes=('mini_refseq',),
             tags=('constrain', 'trace', 'keep', 'coordinates', 'column_kind')),
    Strategy('time_tight', 'exhaustive radius 2 kb, live-path cap 100000, a tight '
             'time_budget_ms (200 SRA / 50 UHGG / 5 mini) -> walks partial; NOT '
             'deterministic (graphlet and full differ)',
             lambda i, s: {'exhaustive': True,
                           'bounds': {'max_extension_bp': 2000, 'max_live_paths': 100000,
                                      'time_budget_ms': TIGHT_MS[i]}},
             deterministic=False, tags=('constrain', 'time_budget')),
    Strategy('batch_refuse', 'multi-seed batch, exhaustive radius 100, max_seed_labels 3: '
             'seeds with > 3 carriers are refused per seed ({seed, error}, no graphlet)',
             lambda i, s: {'exhaustive': True, 'labels': {'max_seed_labels': 3},
                           'bounds': {'max_extension_bp': 100, 'max_live_paths': 500}},
             batch=True, single=False, tags=('constrain', 'batch', 'derivation_error')),
])


def strategy(name, index, seed=None):
    """The strategy dict of |name| for |index| (and |seed|: switch1 adds the seed's
    switch_extra labels). A fresh deep copy; output.detail/timing are not set."""
    s = STRATEGIES[name]
    return copy.deepcopy(s.build(index, seed or {}))


# ------------------------------------------------------------------ seed catalog

def load_catalog():
    try:
        with open(CATALOG_PATH) as f:
            return json.load(f)
    except (OSError, ValueError):
        return {}


def catalog_seeds(index, kind=None, tier=None):
    """The catalog's seeds of |index| (list of dicts; see build_catalog for the keys)."""
    seeds = (load_catalog().get('indexes') or {}).get(index) or []
    return [s for s in seeds if (kind is None or s['kind'] == kind)
            and (tier is None or s.get('tier') == tier)]


def catalog_seed(index, name):
    for s in catalog_seeds(index):
        if s['name'] == name:
            return s
    raise KeyError('no seed %r in the %s catalog' % (name, index))


def seed_sequences(seed):
    """The sequences a request carries for |seed| (one, or several for 'batch')."""
    return list(seed['sequences']) if seed['kind'] == 'batch' else [seed['sequence']]



def _portable(source):
    """A seed's source with the home directory written as ~: the catalog is copied into
    the committed offline fixtures, which must not carry a local user's paths."""
    if isinstance(source, dict) and isinstance(source.get('file'), str):
        source = dict(source)
        if source['file'].startswith(HOME + os.sep):
            source['file'] = '~' + source['file'][len(HOME):]
    return source

def _seed_record(index, name, kind, tier, seq, source, res, note='', **extra):
    carriers = full_carriers(res)
    rec = OrderedDict([
        ('name', name), ('index', index), ('kind', kind), ('tier', tier),
        ('length_bp', len(seq)), ('num_kmers', res['num_kmers']),
        ('graph_runs', res.get('graph_runs')),
        ('n_labels', len(carriers)),
        ('labels_touching', len(res.get('labels') or [])),
        ('labels_truncated', res.get('labels_truncated')),
        ('label_kind', (res.get('labels') or [{}])[0].get('kind') if res.get('labels') else None),
        ('source', _portable(source)), ('note', note),
        ('labels', carriers),
        ('sequence', seq),
    ])
    rec.update(extra)
    return rec


def _covered(res):
    return res.get('graph_runs') == [[0, res['num_kmers']]]


def _windows(rng, records, width, n_max):
    """Random windows (rec_name, start, seq) over |records| [(name, seq)], ACGT only."""
    tries = 0
    usable = [(h, s) for h, s in records if len(s) >= width + 2]
    while tries < n_max and usable:
        tries += 1
        h, s = usable[rng.randrange(len(usable))]
        a = rng.randrange(0, len(s) - width)
        w = s[a:a + width]
        if set(w) <= set('ACGT'):
            yield h, a, w


def _genome_records():
    out = []
    for path in QUERY_MATERIAL['genomes']:
        for h, s in read_fasta(path):
            out.append((os.path.basename(path), h, s))
    return out


def _pick_random(index, rng, records, width, n, *, kind_arg=None, min_labels=1,
                 max_labels=None, tries=300, log=print, exclude=()):
    """n windows from |records| [(file, header, seq)] fully in the graph with
    min_labels <= full carriers (<= max_labels)."""
    got = []
    seen = set(exclude)
    flat = [(f + ':' + h, s) for f, h, s in records]
    for rec, a, w in _windows(rng, flat, width, tries):
        if w in seen:
            continue
        seen.add(w)
        res = oracle_profile(index, w, kind=kind_arg)
        c = len(full_carriers(res))
        if not _covered(res) or c < min_labels or (max_labels is not None and c > max_labels):
            continue
        f, _, h = rec.partition(':')
        got.append((f, h, a, w, res))
        if len(got) >= n:
            break
    if len(got) < n:
        log('  %s: only %d/%d windows of %d bp found' % (index, len(got), n, width))
    return got


def _best_candidate_window(index, seq, kind=None):
    """The longest /resolve candidate (labels with one identical maximal run) of |seq|
    as a window: (start_bp, end_bp, window)."""
    res = oracle_profile(index, seq, kind=kind)
    cands = res.get('candidates') or []
    if not cands:
        return None
    c = cands[0]
    a, b = c['bp_interval']
    return a, b, seq[a:b]


def _catalog_sra(rng, log):
    idx = 'sra'
    seeds = []
    r16s = json.load(open(QUERY_MATERIAL['r16s']))
    pz = r16s['PZ326290.1'][9:1405]
    res = oracle_profile(idx, pz)
    seeds.append(_seed_record(idx, 'sra_16s_PZ326290', '16s', 'core', pz,
                              {'file': QUERY_MATERIAL['r16s'], 'record': 'PZ326290.1',
                               'start': 9, 'end': 1405}, res,
                              'the 1,396 bp 16S window of PZ326290.1 (6 SRA runs carry it)'))
    ok = r16s['OK510341.1']
    a, b, w = _best_candidate_window(idx, ok)
    res = oracle_profile(idx, w)
    seeds.append(_seed_record(idx, 'sra_16s_OK510341', '16s', 'core', w,
                              {'file': QUERY_MATERIAL['r16s'], 'record': 'OK510341.1',
                               'start': a, 'end': b}, res,
                              'the longest /resolve candidate of the other 16S'))
    genomes = _genome_records()
    picks = _pick_random(idx, rng, genomes, 50, 10, log=log)
    for i, (f, h, a, w, res) in enumerate(picks):
        seeds.append(_seed_record(idx, 'sra_rand50_%02d' % i, 'random',
                                  'core' if i < 3 else 'bulk', w,
                                  {'file': f, 'record': h, 'start': a, 'end': a + 50}, res,
                                  'random 50-mer of a gut genome, fully in the graph'))
    # hub: a random window carried by >= 300 runs (<= 1000 = the default max_seed_labels,
    # so constrain and exhaustive accept the derived set). Random gut-genome 50-mers do
    # not reach 300 carriers (0 of 400 tried); ~3 % of the 16S 50-mers do, so the random
    # windows are drawn from the two 16S sequences first
    taken = {s['sequence'] for s in seeds}
    hub = []
    for pool in ([('r16s.json', k, v) for k, v in r16s.items()], genomes):
        hub = _pick_random(idx, rng, pool, 50, 1, min_labels=300, max_labels=1000,
                           tries=400, log=log, exclude=taken)
        if hub:
            break
    for f, h, a, w, res in hub:
        seeds.append(_seed_record(idx, 'sra_hub', 'hub', 'core', w,
                                  {'file': f, 'record': h, 'start': a, 'end': a + 50}, res,
                                  'a random 50-mer carried by >= 300 SRA runs'))
    return seeds


def _first_beyond_labels(index, seq, carriers, limit=8):
    """Labels recorded beyond the seed's right end by a short annotate beam (width 1,
    300 bp), minus the seed's carriers: switch targets for the switch1 cell."""
    st = _beam(1, 300)
    st['direction'] = 'right'
    st['output'] = {'detail': 'full', 'timing': False}
    try:
        out = post(index, 'traverse', {'seeds': [{'sequence': seq}], 'strategy': st},
                   timeout=CELL_TIMEOUT)
    except (RealServerError, OSError) as e:
        return [], 'probe failed: %s' % e
    r = out['results'][0]
    if 'error' in r:
        return [], 'probe failed: %s' % r['error']
    names = [l['name'] for l in r.get('label_dict') or []]
    arm = r['arms']['right']
    counts = {}
    for seg in arm.get('segments') or []:
        for ls in seg.get('label_sets') or []:
            for lid in ls.get('labels') or []:
                counts[lid] = counts.get(lid, 0) + 1
    car = set(carriers)
    beyond = [names[i] for i, _ in sorted(counts.items(), key=lambda x: (-x[1], x[0]))
              if names[i] not in car]
    maxbp = max((p['length_bp'] for p in arm.get('paths') or []), default=0)
    return beyond[:limit], 'beam1 right: %d bp, %d labels beyond the carriers' % (
        maxbp, len(beyond))


def _catalog_uhgg(rng, log):
    idx = 'uhgg'
    seeds = []
    r16s = json.load(open(QUERY_MATERIAL['r16s']))
    best = None
    for acc, s in r16s.items():
        got = _best_candidate_window(idx, s)
        if got and (best is None or got[1] - got[0] > best[2] - best[1]):
            best = (acc, got[0], got[1], got[2])
    if best:
        acc, a, b, w = best
        res = oracle_profile(idx, w)
        seeds.append(_seed_record(idx, 'uhgg_16s', '16s', 'core', w,
                                  {'file': QUERY_MATERIAL['r16s'], 'record': acc,
                                   'start': a, 'end': b}, res,
                                  'the longest /resolve candidate window of the two 16S'))
    genomes = _genome_records()
    picks = _pick_random(idx, rng, genomes, 100, 10, log=log)
    for i, (f, h, a, w, res) in enumerate(picks):
        seeds.append(_seed_record(idx, 'uhgg_rand100_%02d' % i, 'random',
                                  'core' if i < 3 else 'bulk', w,
                                  {'file': f, 'record': h, 'start': a, 'end': a + 100}, res,
                                  'random 100 bp genome window, fully in the graph'))
    # contig end: the last 100 bp of a contig; prefer one where the graph continues past
    # the end under labels that do not carry the seed (switch targets)
    chosen = None
    fallback = None
    for f, h, s in genomes:
        if len(s) < 300:
            continue
        w = s[-100:]
        if not set(w) <= set('ACGT'):
            continue
        res = oracle_profile(idx, w)
        if not _covered(res) or not full_carriers(res):
            continue
        if fallback is None:
            fallback = (f, h, len(s) - 100, w, res, [], 'no continuation probed')
        extra, note = _first_beyond_labels(idx, w, full_carriers(res))
        if extra:
            chosen = (f, h, len(s) - 100, w, res, extra, note)
            break
    chosen = chosen or fallback
    if chosen:
        f, h, a, w, res, extra, note = chosen
        seeds.append(_seed_record(idx, 'uhgg_contig_end', 'contig_end', 'core', w,
                                  {'file': f, 'record': h, 'start': a, 'end': a + 100}, res,
                                  'the last 100 bp of a contig (right arm leaves the '
                                  'genome; %s)' % note, switch_extra=extra))
    cons = _pick_random(idx, rng, genomes, 100, 1, min_labels=5, tries=400, log=log,
                        exclude={s['sequence'] for s in seeds})
    if not cons:
        cons = _pick_random(idx, rng, [('r16s.json', k, v) for k, v in r16s.items()], 100, 1,
                            min_labels=5, tries=200, log=log)
    for f, h, a, w, res in cons:
        seeds.append(_seed_record(idx, 'uhgg_conserved', 'conserved', 'core', w,
                                  {'file': f, 'record': h, 'start': a, 'end': a + 100}, res,
                                  'a 100 bp window carried by >= 5 genomes'))
    return seeds


def _catalog_mini(rng, log):
    idx = 'mini_refseq'
    seeds = []
    for name, key, note in (('mini_ndm1', 'ndm1', 'the whole blaNDM-1 gene (813 bp)'),
                            ('mini_ndm1_rc', 'ndm1_rc',
                             'blaNDM-1 reverse complement (basic graph: the other strand)')):
        h, s = read_fasta(QUERY_MATERIAL[key])[0]
        res = oracle_profile(idx, s, kind='header')
        if _covered(res) and full_carriers(res):
            seeds.append(_seed_record(idx, name, 'gene', 'core', s,
                                      {'file': QUERY_MATERIAL[key], 'record': h, 'start': 0,
                                       'end': len(s)}, res, note))
    rdir = QUERY_MATERIAL['mini_records']
    records = []
    for fn in sorted(os.listdir(rdir)):
        if fn.endswith('.fa'):
            for h, s in read_fasta(os.path.join(rdir, fn)):
                if 2000 <= len(s) <= 500000:
                    records.append((fn, h, s))
    picks = _pick_random(idx, rng, records, 200, 6, kind_arg='header', log=log)
    for i, (f, h, a, w, res) in enumerate(picks):
        seeds.append(_seed_record(idx, 'mini_win200_%02d' % i, 'record_window',
                                  'core' if i < 4 else 'bulk', w,
                                  {'file': os.path.join(rdir, f), 'record': h, 'start': a,
                                   'end': a + 200}, res,
                                  '200 bp window cut from a record'))
    seeds.extend(_repeat_seeds(idx, log))
    # a record end (the right arm leaves the record)
    for f, h, s in sorted(records, key=lambda r: len(r[2])):
        w = s[-200:]
        res = oracle_profile(idx, w, kind='header')
        if _covered(res) and full_carriers(res):
            seeds.append(_seed_record(idx, 'mini_record_end', 'record_end', 'core', w,
                                      {'file': os.path.join(rdir, f), 'record': h,
                                       'start': len(s) - 200, 'end': len(s)}, res,
                                      'the last 200 bp of a record'))
            break
    return seeds


# 150 bp windows of an insertion sequence that NC_013122.1 holds three times and
# NZ_CABHKL010000003.1 three to five times (the 6-copy repeat of the coordinates plan, seen
# from another record): seeds whose coordinate lists hold several occurrences
REPEAT_WINDOWS = (('NC_013122.1.fa', 'NC_013122.1', 59505),
                  ('NC_013122.1.fa', 'NC_013122.1', 60104),
                  ('NC_013122.1.fa', 'NC_013122.1', 59880))
REPEAT_BP = 150


def _repeat_seeds(index, log=print):
    """The repeat windows (kind 'repeat', NARROW_KINDS: trace strategies only), each
    checked fully in the graph with its carriers from /resolve."""
    out = []
    rdir = QUERY_MATERIAL['mini_records']
    for i, (fn, rec, start) in enumerate(REPEAT_WINDOWS):
        path = os.path.join(rdir, fn)
        try:
            s = dict(read_fasta(path))[rec][start:start + REPEAT_BP].upper()
        except (OSError, KeyError):
            log('  %s: no record %s in %s' % (index, rec, path))
            continue
        res = oracle_profile(index, s, kind='header')
        if _covered(res) and full_carriers(res):
            out.append(_seed_record(index, 'mini_rep_%02d' % i, 'repeat', 'core', s,
                                    {'file': path, 'record': rec, 'start': start,
                                     'end': start + REPEAT_BP}, res,
                                    'a 150 bp window of an insertion sequence repeated in '
                                    'several records'))
    return out


def extend_catalog_repeats(log=print):
    """Add the repeat windows to an existing mini_refseq catalog (the other seeds kept as
    they are, so that the cached cells stay valid)."""
    cat = load_catalog()
    seeds = (cat.get('indexes') or {}).get('mini_refseq')
    if not seeds or any(s['kind'] == 'repeat' for s in seeds):
        return cat
    require_server('mini_refseq')
    new = _repeat_seeds('mini_refseq', log)
    batch = [s for s in seeds if s['kind'] == 'batch']
    cat['indexes']['mini_refseq'] = [s for s in seeds if s['kind'] != 'batch'] + new + batch
    cat['updated'] = time.strftime('%Y-%m-%dT%H:%M:%S')
    _write_json(CATALOG_PATH, cat)
    log('mini_refseq: %d repeat seeds added' % len(new))
    return cat


def _batch_seed(index, seeds):
    """A 3-seed batch: two seeds with 1..3 carriers and one with more (so batch_refuse
    refuses exactly that one), all short."""
    few = [s for s in seeds if 1 <= s['n_labels'] <= 3 and s['length_bp'] <= 300]
    many = [s for s in seeds if 3 < s['n_labels'] <= 1000 and s['length_bp'] <= 1000]
    if len(few) < 2 or not many:
        return None
    members = [few[0], many[0], few[1]]
    return OrderedDict([
        ('name', '%s_batch3' % ('mini' if index == 'mini_refseq' else index)),
        ('index', index), ('kind', 'batch'), ('tier', 'core'),
        ('members', [m['name'] for m in members]),
        ('members_n_labels', [m['n_labels'] for m in members]),
        ('sequences', [m['sequence'] for m in members]),
        ('note', 'three seeds in one request; batch_refuse refuses the middle one '
                 '(> 3 carriers) per seed'),
    ])


CATALOG_BUILDERS = {'sra': _catalog_sra, 'uhgg': _catalog_uhgg, 'mini_refseq': _catalog_mini}
CATALOG_RNG_SEED = 20261002


def build_catalog(indexes=None, *, refresh=False, rng_seed=CATALOG_RNG_SEED, log=print):
    """Build (or extend) <scratch>/catalog.json from real /resolve calls. Per seed:
    name, index, kind (16s | random | hub | contig_end | conserved | gene |
    record_window | record_end | batch), tier (core: every strategy; bulk: the cheap
    ones), sequence, length_bp, num_kmers, graph_runs ([[0, num_kmers]]), labels (the
    full-length carriers = the derived seed labels), n_labels, labels_touching,
    label_kind, source {file, record, start, end}, note; switch_extra (contig end)."""
    cat = load_catalog() or {}
    cat.setdefault('format', 1)
    cat.setdefault('indexes', {})
    for index in indexes or INDEXES:
        if cat['indexes'].get(index) and not refresh:
            log('%s: catalog present (%d seeds), --refresh to rebuild'
                % (index, len(cat['indexes'][index])))
            if index == 'mini_refseq':
                # a catalog made before the repeat windows existed gains them
                cat = extend_catalog_repeats(log)
            continue
        require_server(index)
        t0 = time.time()
        rng = random.Random('%s:%s' % (rng_seed, index))
        seeds = CATALOG_BUILDERS[index](rng, log)
        b = _batch_seed(index, seeds)
        if b:
            seeds.append(b)
        cat['indexes'][index] = seeds
        cat.setdefault('identity', {})[index] = identity(index)
        cat['updated'] = time.strftime('%Y-%m-%dT%H:%M:%S')
        _write_json(CATALOG_PATH, cat)
        log('%s: %d seeds in %.1f s: %s' % (
            index, len(seeds), time.time() - t0,
            ', '.join('%s(%s,%s)' % (s['name'], s.get('length_bp', '-'),
                                     s.get('n_labels', s.get('members_n_labels')))
                      for s in seeds)))
    return cat


# ------------------------------------------------------------------ the matrix and the cache

def cell_id(seed_name, strategy_name):
    return '%s__%s' % (seed_name, strategy_name)


def matrix(index, tiers=None):
    """[(seed_name, strategy_name)] of |index| in catalog x STRATEGIES order (the full
    product; tiers=('core',) or ('core', 'bulk') gives the reduced matrix)."""
    out = []
    for seed in catalog_seeds(index):
        for name, st in STRATEGIES.items():
            if st.applies(index, seed, tiers):
                out.append((seed['name'], name))
    return out


def build_request(index, seed, strategy_name, detail='graphlet', timing=False):
    """The /traverse request of a cell: the seed(s) without labels (derived), the named
    strategy, output.detail and output.timing set."""
    if isinstance(seed, str):
        seed = catalog_seed(index, seed)
    st = strategy(strategy_name, index, seed)
    out = st.setdefault('output', {})
    out['detail'] = detail
    out['timing'] = timing
    return {'seeds': [{'sequence': s} for s in seed_sequences(seed)], 'strategy': st}


def with_detail(request, detail):
    """A copy of |request| with output.detail replaced."""
    r = copy.deepcopy(request)
    r['strategy'].setdefault('output', {})['detail'] = detail
    return r


def _cell_paths(index, cell):
    base = os.path.join(CACHE_DIR, index, cell)
    return {'request': base + '.request.json.gz', 'graphlet': base + '.graphlet.json.gz',
            'full': base + '.full.json.gz'}


_MANIFEST_LOCK = threading.Lock()
_MANIFEST_MEMO = {'mtime': None, 'data': None}


def load_manifest():
    try:
        mt = os.path.getmtime(MANIFEST_PATH)
    except OSError:
        return {}
    if _MANIFEST_MEMO['mtime'] != mt:
        try:
            with open(MANIFEST_PATH) as f:
                _MANIFEST_MEMO['data'] = json.load(f)
            _MANIFEST_MEMO['mtime'] = mt
        except (OSError, ValueError):
            return {}
    return _MANIFEST_MEMO['data']


def cell_meta(index, cell):
    return (load_manifest().get('cells') or {}).get('%s/%s' % (index, cell))


def cells(index=None, *, status=('ok',), strategy=None, seed=None, kind=None, tag=None,
          deterministic=None, heavy=False, max_full_bytes=None, max_body_bytes=None,
          tier=None):
    """Manifest rows (dicts) of the cached cells, filtered; default status 'ok' (both
    detail levels cached) and heavy=False (cells whose full JSON > 64 MB or bodies > 8 MB
    are left out: heavy=None includes them, heavy=True selects only them). status=None
    returns every row (skipped/error too)."""
    rows = list((load_manifest().get('cells') or {}).values())
    if isinstance(status, str):
        status = (status,)
    out = []
    for r in rows:
        if index is not None and r['index'] != index:
            continue
        if status is not None and r.get('status') not in status:
            continue
        if strategy is not None and r['strategy'] not in (
                (strategy,) if isinstance(strategy, str) else strategy):
            continue
        if seed is not None and r['seed'] != seed:
            continue
        if kind is not None and r.get('seed_kind') not in (
                (kind,) if isinstance(kind, str) else kind):
            continue
        if tag is not None and tag not in (r.get('tags') or ()):
            continue
        if tier is not None and r.get('seed_tier') not in (
                (tier,) if isinstance(tier, str) else tier):
            continue
        if deterministic is not None and bool(r.get('deterministic')) != deterministic:
            continue
        if heavy is not None and bool(r.get('heavy')) != heavy:
            continue
        fb = (r.get('full') or {}).get('json_bytes') or 0
        bb = (r.get('graphlet') or {}).get('body_bytes') or 0
        if max_full_bytes is not None and fb > max_full_bytes:
            continue
        if max_body_bytes is not None and bb > max_body_bytes:
            continue
        out.append(r)
    order = {i: n for n, i in enumerate(INDEXES)}
    sorder = {s: n for n, s in enumerate(STRATEGIES)}
    out.sort(key=lambda r: (order.get(r['index'], 99), r.get('seed_order', 0),
                            sorder.get(r['strategy'], 99)))
    return out


class CachedCell:
    """One cached (index, seed, strategy) cell. JSON is parsed on first access."""

    def __init__(self, index, cell, meta=None):
        self.index = index
        self.cell = cell
        self.meta = meta if meta is not None else cell_meta(index, cell)
        if self.meta is None:
            raise KeyError('cell %s/%s is not in the cache manifest' % (index, cell))
        self.paths = _cell_paths(index, cell)
        self._memo = {}

    def __repr__(self):
        return 'CachedCell(%s/%s, %s)' % (self.index, self.cell, self.meta.get('status'))

    @property
    def seed_name(self):
        return self.meta['seed']

    @property
    def strategy_name(self):
        return self.meta['strategy']

    @property
    def deterministic(self):
        return bool(self.meta.get('deterministic'))

    def _load(self, what):
        if what not in self._memo:
            p = self.paths[what]
            if not os.path.exists(p):
                raise FileNotFoundError(p)
            self._memo[what] = _read_gz_bytes(p)
        return self._memo[what]

    def raw(self, what):
        """The exact bytes the server sent (decoded): what in request|graphlet|full."""
        return self._load(what)

    @property
    def request(self):
        """The graphlet request as sent (detail 'graphlet'); request_for('full') is the
        full one (identical apart from output.detail)."""
        return json.loads(self._load('request'))

    def request_for(self, detail):
        return with_detail(self.request, detail)

    @property
    def graphlet_response(self):
        if 'graphlet_json' not in self._memo:
            self._memo['graphlet_json'] = json.loads(self._load('graphlet'))
        return self._memo['graphlet_json']

    @property
    def full_response(self):
        if 'full_json' not in self._memo:
            self._memo['full_json'] = json.loads(self._load('full'))
        return self._memo['full_json']

    @property
    def has_full(self):
        return os.path.exists(self.paths['full'])

    def n_results(self):
        return len(self.graphlet_response.get('results') or [])

    def graphlet_result(self, i=0):
        return self.graphlet_response['results'][i]

    def full_result(self, i=0):
        return self.full_response['results'][i]

    def seed_sequence(self, i=0):
        return self.request['seeds'][i]['sequence']

    def text(self, i=0):
        """The MGT body of result i (None for an error result)."""
        return self.graphlet_result(i).get('graphlet')

    def texts(self):
        return [r.get('graphlet') for r in self.graphlet_response.get('results') or []]

    def graphlet(self, i=0):
        """metagraph.traverse.from_response(results[i], response) -- fresh each call
        (tests may mutate it); None for an error result."""
        from metagraph.traverse import from_response
        r = self.graphlet_result(i)
        if 'graphlet' not in r:
            return None
        return from_response(r, self.graphlet_response)

    def graphlets(self):
        return [self.graphlet(i) for i in range(self.n_results())]

    @property
    def seed_record(self):
        try:
            return catalog_seed(self.index, self.seed_name)
        except KeyError:
            return None


def load_cell(index, cell):
    """CachedCell for |cell| (a cell id, or a manifest row)."""
    if isinstance(cell, dict):
        return CachedCell(cell['index'], cell['cell'], cell)
    return CachedCell(index, cell)


def iter_cells(index=None, **filters):
    for row in cells(index, **filters):
        yield CachedCell(row['index'], row['cell'], row)


def _arm_facts(arm):
    counts = arm.get('counts') or {}
    cap = arm.get('cap_trigger')
    return {'status': arm.get('status'), 'complete_to_bp': arm.get('complete_to_bp'),
            'scope': arm.get('completeness_scope'),
            'segments': counts.get('segments'), 'leaves': counts.get('leaves'),
            'runs': counts.get('runs'), 'max_bp': counts.get('max_bp'),
            'cap': cap.get('reason') if isinstance(cap, dict) else None,
            'limitations': [l.get('kind') for l in arm.get('limitations') or []]}


def _seed_facts(r):
    if 'graphlet' not in r:
        return {'error': r.get('error')}
    arms = {side: _arm_facts(a) for side, a in (r.get('arms') or {}).items()}
    return {'outcome': r.get('outcome'), 'label_mode': r.get('label_mode'),
            'num_seed_labels': (r.get('seed') or {}).get('num_seed_labels'),
            'graphlet_bytes': r.get('graphlet_bytes'), 'graphlet_lines': r.get('graphlet_lines'),
            'limitations': [l.get('kind') for l in r.get('limitations') or []],
            'arms': arms}


def _time_tripped(resp):
    for r in resp.get('results') or []:
        for a in (r.get('arms') or {}).values():
            cap = a.get('cap_trigger')
            if isinstance(cap, dict) and cap.get('reason') == 'time_budget':
                return True
            lbr = ((a.get('counts') or {}).get('leaves_by_reason') or {})
            if lbr.get('time_budget'):
                return True
    return False


def _fetch(index, request, max_seconds):
    try:
        res = http(index, 'POST', 'traverse', request, timeout=max_seconds)
    except (socket.timeout, TimeoutError) as e:
        return None, 'timeout', 'no answer within %.0f s (%s)' % (max_seconds, e)
    except urllib.error.URLError as e:
        if isinstance(e.reason, (socket.timeout, TimeoutError)):
            return None, 'timeout', 'no answer within %.0f s' % max_seconds
        return None, 'unreachable', str(e)
    except OSError as e:
        return None, 'unreachable', str(e)
    if res.elapsed_s > max_seconds:
        return res, 'slow', 'answered after %.1f s > %.0f s' % (res.elapsed_s, max_seconds)
    if not res.ok:
        return res, 'http_%d' % res.status, res.error or res.body[:300].decode('utf-8', 'replace')
    return res, 'ok', None


def fetch_cell(index, seed, strategy_name, *, max_seconds=CELL_TIMEOUT, seed_order=0):
    """Fetch one cell (graphlet, then full), write its files, return the manifest row."""
    if isinstance(seed, str):
        seed = catalog_seed(index, seed)
    cid = cell_id(seed['name'], strategy_name)
    paths = _cell_paths(index, cid)
    st = STRATEGIES[strategy_name]
    req = build_request(index, seed, strategy_name, 'graphlet', timing=False)
    row = OrderedDict([
        ('index', index), ('cell', cid), ('seed', seed['name']), ('seed_kind', seed['kind']),
        ('seed_tier', seed.get('tier')), ('seed_order', seed_order),
        ('seed_length_bp', seed.get('length_bp')), ('seed_n_labels', seed.get('n_labels')),
        ('strategy', strategy_name), ('tags', list(st.tags)),
        ('deterministic', st.deterministic),
        ('files', {k: os.path.relpath(v, CACHE_DIR) for k, v in paths.items()}),
        ('fetched_at', time.strftime('%Y-%m-%dT%H:%M:%S')),
    ])
    _write_json_gz(paths['request'], req)
    for p in (paths['graphlet'], paths['full']):
        if os.path.exists(p):
            os.remove(p)
    g, gstat, gerr = _fetch(index, req, max_seconds)
    row['graphlet'] = {'status': gstat, 'error': gerr}
    if g is not None:
        row['graphlet'].update({'http_status': g.status, 'elapsed_s': round(g.elapsed_s, 3),
                                'wire_bytes': g.wire_bytes, 'json_bytes': len(g.body)})
    if gstat != 'ok':
        row['status'] = 'skipped' if gstat in ('timeout', 'slow') else 'error'
        row['error'] = gerr
        if g is not None and g.ok:            # slow but answered: keep it
            _write_json_gz(paths['graphlet'], g.body)
        return row
    _write_json_gz(paths['graphlet'], g.body)
    resp = g.json
    row['graphlet']['body_bytes'] = sum(r.get('graphlet_bytes') or 0
                                        for r in resp.get('results') or [])
    row['graphlet']['lines'] = sum(r.get('graphlet_lines') or 0
                                   for r in resp.get('results') or [])
    row['n_results'] = len(resp.get('results') or [])
    row['n_errors'] = sum(1 for r in resp.get('results') or [] if 'graphlet' not in r)
    row['results'] = [_seed_facts(r) for r in resp.get('results') or []]
    row['time_budget_tripped'] = _time_tripped(resp)
    f, fstat, ferr = _fetch(index, with_detail(req, 'full'), max_seconds)
    row['full'] = {'status': fstat, 'error': ferr}
    if f is not None:
        row['full'].update({'http_status': f.status, 'elapsed_s': round(f.elapsed_s, 3),
                            'wire_bytes': f.wire_bytes, 'json_bytes': len(f.body)})
    if fstat == 'ok' or (f is not None and f.ok):
        _write_json_gz(paths['full'], f.body)
    row['heavy'] = (((row.get('full') or {}).get('json_bytes') or 0) > HEAVY_FULL_BYTES
                    or row['graphlet']['body_bytes'] > HEAVY_BODY_BYTES)
    if fstat == 'ok':
        row['status'] = 'ok'
    elif fstat in ('timeout', 'slow'):
        row['status'] = 'graphlet_only'
        row['error'] = 'full: ' + (ferr or fstat)
    else:
        row['status'] = 'error'
        row['error'] = 'full: ' + (ferr or fstat)
    return row


def _save_manifest(update_rows, servers_info=None, extra=None):
    with _MANIFEST_LOCK:
        m = {}
        try:
            with open(MANIFEST_PATH) as fh:
                m = json.load(fh)
        except (OSError, ValueError):
            m = {}
        m.setdefault('format', 1)
        m.setdefault('created', time.strftime('%Y-%m-%dT%H:%M:%S'))
        m['updated'] = time.strftime('%Y-%m-%dT%H:%M:%S')
        m['scratch'] = SCRATCH
        m.setdefault('cells', {})
        if servers_info:
            m.setdefault('servers', {}).update(servers_info)
        if extra:
            m.update(extra)
        for row in update_rows:
            m['cells']['%s/%s' % (row['index'], row['cell'])] = row
        _write_json(MANIFEST_PATH, m)


def fill_cache(indexes=None, *, refresh=False, max_seconds=CELL_TIMEOUT, workers=2,
               log=print, only_strategies=None, only_seeds=None):
    """Fetch every matrix cell of |indexes| once with detail graphlet and full (timing
    off). Cells already 'ok' are kept unless |refresh|. The indexes run in parallel,
    |workers| cells at a time per index. Returns the manifest."""
    indexes = [i for i in (indexes or INDEXES) if server(i) is not None]
    servers_info = {}
    for i in indexes:
        require_server(i)
        ident = identity(i)
        servers_info[i] = dict(ident, host=server(i)[0], port=server(i)[1])
    _save_manifest([], servers_info, {'max_seconds': max_seconds, 'workers_per_index': workers,
                                      'timing': False})
    t_all = time.time()

    def run_index(index):
        seeds = {s['name']: (n, s) for n, s in enumerate(catalog_seeds(index))}
        todo = []
        for seed_name, sname in matrix(index):
            if only_strategies and sname not in only_strategies:
                continue
            if only_seeds and seed_name not in only_seeds:
                continue
            cid = cell_id(seed_name, sname)
            meta = cell_meta(index, cid)
            if meta and meta.get('status') == 'ok' and not refresh:
                continue
            todo.append((seed_name, sname))
        log('%s: %d cells to fetch' % (index, len(todo)))
        lock = threading.Lock()
        it = iter(todo)
        t0 = time.time()
        done = [0]

        def worker():
            while True:
                with lock:
                    nxt = next(it, None)
                if nxt is None:
                    return
                seed_name, sname = nxt
                order, seed = seeds[seed_name]
                try:
                    row = fetch_cell(index, seed, sname, max_seconds=max_seconds,
                                     seed_order=order)
                except Exception as e:      # noqa: BLE001 -- recorded, never fatal
                    row = {'index': index, 'cell': cell_id(seed_name, sname),
                           'seed': seed_name, 'seed_kind': seed['kind'],
                           'seed_order': order, 'strategy': sname, 'status': 'error',
                           'error': '%s: %s' % (type(e).__name__, e)}
                _save_manifest([row])
                with lock:
                    done[0] += 1
                    g = row.get('graphlet') or {}
                    fl = row.get('full') or {}
                    log('  %s %3d/%d %-40s %-13s g %6.2fs %9s B | full %6.2fs %11s B%s' % (
                        index, done[0], len(todo), row['cell'], row['status'],
                        g.get('elapsed_s') or 0, g.get('body_bytes', '-'),
                        fl.get('elapsed_s') or 0, fl.get('json_bytes', '-'),
                        ('  ' + str(row.get('error'))[:120]) if row.get('error') else ''))
        ths = [threading.Thread(target=worker, daemon=True) for _ in range(max(1, workers))]
        for t in ths:
            t.start()
        for t in ths:
            t.join()
        log('%s: done in %.0f s' % (index, time.time() - t0))

    threads = [threading.Thread(target=run_index, args=(i,), daemon=True) for i in indexes]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    _save_manifest([], extra={'last_fill_s': round(time.time() - t_all, 1)})
    return load_manifest()


# ------------------------------------------------------------------ reporting

def report(log=print):
    m = load_manifest()
    rows = list((m.get('cells') or {}).values())
    log('cache: %s (%d cells, last fill %s s)' % (CACHE_DIR, len(rows), m.get('last_fill_s')))
    for index in INDEXES:
        rs = [r for r in rows if r['index'] == index]
        if not rs:
            continue
        by = {}
        for r in rs:
            by[r['status']] = by.get(r['status'], 0) + 1
        g_t = sum((r.get('graphlet') or {}).get('elapsed_s') or 0 for r in rs)
        f_t = sum((r.get('full') or {}).get('elapsed_s') or 0 for r in rs)
        log('%s: %d cells %s; request time graphlet %.0f s + full %.0f s' % (
            index, len(rs), by, g_t, f_t))
        for r in rs:
            if r['status'] != 'ok':
                log('   %-12s %s: %s' % (r['status'], r['cell'], str(r.get('error'))[:200]))
        top = sorted(rs, key=lambda r: -((r.get('full') or {}).get('json_bytes') or 0))[:5]
        for r in top:
            g = r.get('graphlet') or {}
            f = r.get('full') or {}
            log('   largest: %-40s body %9s B (%s lines), graphlet JSON %9s B, full JSON '
                '%11s B, %.1f s + %.1f s' % (
                    r['cell'], g.get('body_bytes'), g.get('lines'), g.get('json_bytes'),
                    f.get('json_bytes'), g.get('elapsed_s') or 0, f.get('elapsed_s') or 0))
        trip = [r['cell'] for r in rs if r.get('time_budget_tripped')]
        if trip:
            log('   time budget tripped: %s' % ', '.join(trip))
        errs = [(r['cell'], i, x.get('error')) for r in rs for i, x in
                enumerate(r.get('results') or []) if x.get('error')]
        for c, i, e in errs:
            log('   seed error %s[%d]: %s' % (c, i, str(e)[:160]))


def status(log=print):
    for index in SERVERS:
        hp = server(index)
        up = server_up(index) if hp else False
        ident = identity(index) or {}
        log('%-12s %-22s %s  ns=%s fp=%s meta=%s' % (
            index, '%s:%s' % hp if hp else 'disabled', 'UP  ' if up else 'DOWN',
            ident.get('index_ns'), (ident.get('index_fp') or '-')[:16],
            ident.get('index_meta_fp')))
    cat = load_catalog()
    for index in INDEXES:
        seeds = (cat.get('indexes') or {}).get(index) or []
        log('catalog %-12s %d seeds, %d matrix cells' % (index, len(seeds),
                                                        len(matrix(index))))
    log('cache   %d ok cells of %d' % (len(cells()), len(cells(status=None))))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('cmd', choices=['status', 'catalog', 'fill', 'report'])
    ap.add_argument('--index', action='append', choices=list(INDEXES))
    ap.add_argument('--refresh', action='store_true')
    ap.add_argument('--max-seconds', type=float, default=CELL_TIMEOUT)
    ap.add_argument('--workers', type=int, default=2)
    ap.add_argument('--strategy', action='append')
    ap.add_argument('--seed', action='append')
    a = ap.parse_args(argv)

    def log(*x):
        print(*x, flush=True)
    if a.cmd == 'status':
        status(log)
    elif a.cmd == 'catalog':
        build_catalog(a.index, refresh=a.refresh, log=log)
    elif a.cmd == 'fill':
        fill_cache(a.index, refresh=a.refresh, max_seconds=a.max_seconds, workers=a.workers,
                   log=log, only_strategies=a.strategy, only_seeds=a.seed)
        report(log)
    elif a.cmd == 'report':
        report(log)


if __name__ == '__main__':
    main()
