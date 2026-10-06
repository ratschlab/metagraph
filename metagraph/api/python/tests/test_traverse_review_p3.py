"""The P3/P4 library items of the review of 2026-10-06 (metagraph.traverse):

  * O21 (X-CONCURRENCY-01, U20-08): the spool's body lifecycle is serialized across the
    processes sharing it (an fcntl.flock of spool/.lock) -- a free in one process deleted
    the body of an entry another was storing, so a put returned a handle whose body was
    gone -- and a free no longer keeps a body for an in-memory entry whose file another
    process removed (a leak with no orphan sweep); sweep() deletes old orphans, and no
    body while an entry file is there but cannot be read;
  * O23 (X-CONCURRENCY-05, X-EFFICIENCY-08): body GC reads each entry file of another
    process once (cached by name and inode), and a sweep collects all its expired bodies in
    one scan -- it re-parsed every foreign entry per expired entry (M x N);
  * O17 (U20-03): a damaged body file is repaired, not reused -- a put checks an existing
    body against its digest, a read checks what it read, and a damaged body expires its
    entry (replayable) and is deleted, so the next put writes it again;
  * O11 (U19-01): to_fasta(leaves=<one-shot iterator>) answers as it did before the
    regression (a TypeError from len()), byte for byte, per arm;
  * O35 (U20-10): a response cut in transfer raises the client's ConnectionError (an
    OSError), which the tools answer as backend_unreachable -- not a bare IncompleteRead --
    saying that the server was reached and may have run the request;
  * G15 (GAP-05-05): graphlet_export without local limits writes atomically: a failed
    re-export leaves the earlier file and no temporary file;
  * L8 (U15-02): a change_cost override that names another model replaces the cost whole
    (a deep merge kept the old model's fields, which the server's strict parse refuses);
    a field the model does not read, and None for an object section, are ValueErrors --
    in annotate mode too, where any model but forbid is one (the server refuses it);
  * L7 (U15-06, X-EFFICIENCY-06): the left-out note's switch searches are confined to the
    names the change_cost table names plus one stand-in -- the same note -- and its
    "could still have entered" clause is a boolean, not a substring test of label names;
  * L3 (U15-05): support_changes' split check bisects the sorted label sets once per
    change, and its charge covers it (it scanned the arrays per removed label, uncharged);
  * L5 (U14-02): annotate walks with claims read an index of each leaf's route ends,
    built once (a pass over every route end per chunk of 32 walks, and its charge);
  * L12 (U16-03, X-EFFICIENCY-05): compare's cuts are made in one parents-first pass over
    their chains (each chain was re-read per cut), charged per segment once.
"""

import copy
import dataclasses
import gc
import hashlib
import http.client
import json
import math
import os
import random
import shutil
import socket
import subprocess
import sys
import tempfile
import threading
import time
import tracemalloc
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402
import test_traverse_mcp as M  # noqa: E402
from test_traverse_scale import COUNTERS, comb_annotate  # noqa: E402

from metagraph.traverse import (  # noqa: E402
    GraphletStore, LocalBudget, TraverseClient, UnknownHandle, export, ops, parse,
)
from metagraph.traverse import budget as B  # noqa: E402
from metagraph.traverse import mcp_tools  # noqa: E402
from metagraph.traverse import store as S  # noqa: E402
from metagraph.traverse.mcp_tools import GraphletTools, ToolLimits  # noqa: E402
from metagraph.traverse.model import Label  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
API = os.path.dirname(HERE)


def _traced(fn):
    gc.collect()
    tracemalloc.start()
    base = tracemalloc.get_traced_memory()[0]
    try:
        out = fn()
        return tracemalloc.get_traced_memory()[1] - base, out
    finally:
        tracemalloc.stop()


def _bodies(spool):
    return sorted(n for n in os.listdir(os.path.join(spool, 'bodies')) if n.endswith('.mgt'))


class _mode_0:
    """|path| unreadable (mode 0) inside the with block."""

    def __init__(self, path):
        self.path = path

    def __enter__(self):
        self.mode = os.stat(self.path).st_mode & 0o777
        os.chmod(self.path, 0)

    def __exit__(self, *exc):
        os.chmod(self.path, self.mode)


# ======================================================================= O21

class _Interleave:
    """Run |other| in a thread from inside |store|'s next entry write -- a put between
    its body and its entry -- and give it |wait| seconds: without the spool lock it
    completes there; with it, it waits for the put."""

    def __init__(self, store, other, wait=0.3):
        self.store, self.other, self.wait = store, other, wait
        self.done_inside = None
        self.thread = None
        self.error = []

    def __enter__(self):
        orig = self.store._write_entry

        def run():
            try:
                self.other()
            except BaseException as e:      # reported by the test
                self.error.append(e)

        def hooked(e):
            if self.thread is None:
                self.thread = threading.Thread(target=run)
                self.thread.start()
                self.thread.join(self.wait)
                self.done_inside = not self.thread.is_alive()
            return orig(e)
        self.store._write_entry = hooked
        return self

    def __exit__(self, *exc):
        del self.store._write_entry
        if self.thread is not None:
            self.thread.join()
        return False


class TestO21SpoolLock(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        self.a = GraphletStore(self.d)
        self.b = GraphletStore(self.d, max_ram_mb=0)
        self.hb = self.b.put(T.doc_json('fork', 'graphlet'), {'seeds': [], 'who': 'B'})

    def assert_usable(self, store, h):
        self.assertIsNotNone(store.graphlet(h))
        self.assertTrue(store.body_text(h).startswith('H mgt 1 '))
        self.assertIn(h, [m['handle'] for m in store.list()])

    def race(self, put):
        """B frees its entry of the same body while A's put is between its body and its
        entry: the free waits for the put and then keeps the body."""
        with _Interleave(self.a, lambda: self.b.free(self.hb)) as x:
            h = put()
        self.assertEqual([], x.error)
        self.assertFalse(x.done_inside, 'the free ran inside the put')
        self.assertNotIn(self.hb, self.b)
        self.assertEqual(1, len(_bodies(self.d)))
        self.assert_usable(self.a, h)
        # and the body goes with the last entry naming it
        self.a.free(h)
        self.assertEqual([], _bodies(self.d))

    def test_put(self):
        for resident in (False, True):
            with self.subTest(resident=resident):
                self.setUp()
                self.race(lambda: self.a.put(T.doc_json('fork', 'graphlet'),
                                             {'seeds': [], 'who': 'A'}, resident=resident))

    def test_load(self):
        self.race(lambda: self.a.load(os.path.join(T.DOCS, 'fork.mgt')))

    def test_put_unparsed(self):
        resp = T.doc_json('fork', 'graphlet')
        self.race(lambda: self.a.put_unparsed(resp['results'][0], resp, {'seeds': []}))

    def test_put_view(self):
        g = self.a.graphlet(self.hb)          # A adopts B's entry: a view shares its body
        gv = copy.copy(g)
        gv.view = {'labels': ['x']}
        self.race(lambda: self.a.put_view(self.hb, gv, request={'seeds': []}))

    def test_put_view_of_a_backing_whose_body_is_gone(self):
        g = self.a.graphlet(self.hb)
        os.unlink(self.a._body_path(self.a.get(self.hb).digest))
        gv = copy.copy(g)
        gv.view = {'labels': ['x']}
        with self.assertRaises(UnknownHandle):
            self.a.put_view(self.hb, gv, request={'seeds': []})
        self.assertEqual([], self.a.list())

    def test_a_stale_entry_keeps_no_body(self):
        # U20-08: A frees its entry of a body whose other entry B removed; A's in-memory
        # copy of that entry kept the body for good
        doc = os.path.join(T.DOCS, 'fork.mgt')
        h1, h2 = self.a.load(doc), self.a.load(doc)
        self.b.free(self.hb)
        self.b.free(h1)
        self.assertEqual(1, len(_bodies(self.d)))
        self.a.free(h2)
        self.assertEqual([], _bodies(self.d))
        self.assertNotIn(h1, self.a)

    def test_sweep_deletes_old_orphans_only(self):
        orphan = os.path.join(self.d, 'bodies', '0' * 64 + '.mgt')
        young = os.path.join(self.d, 'bodies', '1' * 64 + '.mgt')
        tmp = os.path.join(self.d, 'bodies', '.tmp-left')
        for p in (orphan, young, tmp):
            with open(p, 'w') as f:
                f.write('x')
        old = time.time() - S._ORPHAN_GRACE_S - 10
        os.utime(orphan, (old, old))
        os.utime(tmp, (old, old))
        named = _bodies(self.d)
        self.a.sweep()
        left = _bodies(self.d)
        self.assertNotIn('0' * 64 + '.mgt', left)
        self.assertIn('1' * 64 + '.mgt', left)            # younger than the grace
        self.assertFalse(os.path.exists(tmp))
        # the body B's entry names stays, however old
        self.assertEqual(len(named) - 1, len(left))
        self.assert_usable(self.b, self.hb)

    def test_processes(self):
        # three real processes storing, reading and freeing one body in a loop: before the
        # lock about half of the handles a put had just returned answered as expired (no
        # other entry may name the body: B's would keep it)
        self.b.free(self.hb)
        script = '''if 1:
            import os, sys, time
            sys.path.insert(0, %r); sys.path.insert(0, %r)
            import traverse_testlib as T
            from metagraph.traverse import GraphletStore, UnknownHandle
            s = GraphletStore(sys.argv[1], max_ram_mb=0)
            doc = T.doc_json('fork', 'graphlet')
            # all three start together (the imports take longer than the loop)
            open(os.path.join(sys.argv[1], 'ready.%%d' %% os.getpid()), 'w').close()
            while not os.path.exists(os.path.join(sys.argv[1], 'go')):
                time.sleep(0.01)
            lost = 0
            for i in range(int(sys.argv[2])):
                h = s.put(doc, {'seeds': [], 'i': i}, resident=False)
                try:
                    s.body_text(h)
                except UnknownHandle:
                    lost += 1
                    continue
                s.free(h)
            print(lost)
            ''' % (API, HERE)
        procs = [subprocess.Popen([sys.executable, '-c', script, self.d, '600'],
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                 for _ in range(3)]
        deadline = time.time() + 120
        while len([n for n in os.listdir(self.d) if n.startswith('ready.')]) < 3 \
                and time.time() < deadline:
            time.sleep(0.01)
        open(os.path.join(self.d, 'go'), 'w').close()
        outs = [p.communicate(timeout=600) for p in procs]
        for p, (out, err) in zip(procs, outs):
            self.assertEqual(0, p.returncode, err.decode())
            self.assertEqual('0', out.decode().strip(), err.decode())
        self.assertEqual([], _bodies(self.d))

    def test_a_file_system_without_flock(self):
        # some network file systems refuse flock: the store goes on with the threads' lock,
        # as before the spool lock, and says why (the module text states the limit)
        import errno
        if S.fcntl is None:
            self.skipTest('no fcntl on this platform')

        def refuse(fd, op):
            raise OSError(errno.ENOLCK, 'No locks available')
        with mock.patch.object(S.fcntl, 'flock', refuse):
            h = self.a.put(T.doc_json('merge', 'graphlet'), {'seeds': []})
            self.assertEqual('ENOLCK', self.a._lock.unlocked)
            self.a.free(h)
        h = self.a.put(T.doc_json('merge', 'graphlet'), {'seeds': []})
        self.assertIsNone(self.a._lock.unlocked)
        self.assert_usable(self.a, h)

    def refused_open(self, path, err):
        """A store-module open() that raises OSError |err| for |path| (an entry file that
        is there but cannot be read now: EMFILE, EACCES, EIO)."""
        real = open

        def fake(p, *a, **k):
            if isinstance(p, str) and p == path:
                raise OSError(err, os.strerror(err))
            return real(p, *a, **k)
        return mock.patch.object(S, 'open', fake, create=True)

    def test_an_unreadable_entry_file_stops_the_gc(self):
        # the review of the P3 fixes: an entry file that is there but cannot be read was
        # read as naming no body, so the orphan pass deleted every old body such a live
        # entry named (its handle then answered as expired), and a free the body it shared.
        # No body is deleted by a scan that cannot read an entry file
        import errno
        body = self.b._body_path(self.b.get(self.hb).digest)
        entry = self.b._entry_path(self.hb)
        orphan = os.path.join(self.d, 'bodies', '0' * 64 + '.mgt')
        with open(orphan, 'w') as f:
            f.write('x')
        old = time.time() - S._ORPHAN_GRACE_S - 10
        ways = [('EMFILE', lambda: self.refused_open(entry, errno.EMFILE)),
                ('EIO', lambda: self.refused_open(entry, errno.EIO))]
        if os.geteuid() != 0:                 # root reads a mode-0 file
            ways.append(('EACCES', lambda: _mode_0(entry)))
        for name, way in ways:
            with self.subTest(sweep=name):
                os.utime(body, (old, old))
                os.utime(orphan, (old, old))
                with way():
                    self.a.sweep()
                self.assertTrue(os.path.exists(body))
                self.assertTrue(os.path.exists(orphan), 'the scan deleted an orphan')
                self.assert_usable(self.b, self.hb)
            with self.subTest(free=name):
                # A frees its own entry of the body B's unreadable entry names
                h = self.a.load(os.path.join(T.DOCS, 'fork.mgt'))
                self.assertEqual(body, self.a._body_path(self.a.get(h).digest))
                with way():
                    self.a.free(h)
                self.assertNotIn(h, self.a)
                self.assertTrue(os.path.exists(body))
                self.assert_usable(self.b, self.hb)
        # readable again: the orphan goes, the body B's entry names stays however old
        os.utime(body, (old, old))
        self.a.sweep()
        self.assertFalse(os.path.exists(orphan))
        self.assertTrue(os.path.exists(body))
        self.assert_usable(self.b, self.hb)

    def test_the_gc_goes_on_past_a_removed_or_damaged_entry_file(self):
        # an entry removed between the listing and its open, or one that is no JSON (no
        # store can read it), names no body: the scan goes on and deletes what it may
        import errno
        old = time.time() - S._ORPHAN_GRACE_S - 10
        # B's entry file "removed meanwhile" (A has read no entry file yet: none cached):
        # its old body is an orphan to that scan
        body = self.b._body_path(self.b.get(self.hb).digest)
        os.utime(body, (old, old))
        with self.refused_open(self.b._entry_path(self.hb), errno.ENOENT):
            self.a.sweep()
        self.assertFalse(os.path.exists(body))
        damaged = os.path.join(self.d, 'entries', 'g_' + 'a' * 12 + '.json')
        with open(damaged, 'wb') as f:
            f.write(b'\x00\x00 not an entry')
        orphan = os.path.join(self.d, 'bodies', '0' * 64 + '.mgt')
        with open(orphan, 'w') as f:
            f.write('x')
        os.utime(orphan, (old, old))
        self.a.sweep()
        self.assertFalse(os.path.exists(orphan))

    def test_the_lock_is_re_entrant_and_released(self):
        with self.a._lock:
            with self.a._lock:
                h = self.a.put(T.doc_json('fork', 'graphlet'), {'seeds': []})
            self.a.free(h)
        done = []
        t = threading.Thread(target=lambda: done.append(self.b.free(self.hb)))
        t.start()
        t.join(10)
        self.assertEqual([None], done)


# ======================================================================= O23

class _Clock:
    def __init__(self, t=1_000_000.0):
        self.t = t

    def __call__(self):
        return self.t


class TestO23BodyGcCost(unittest.TestCase):
    """A sweep that expires M entries while N entry files of another process exist reads
    each foreign file once: N + O(M) opens, not M x N; a later free reads none of them."""

    def counting_open(self):
        calls = {'n': 0}
        real = open

        def counted(*a, **k):
            calls['n'] += 1
            return real(*a, **k)
        return calls, mock.patch.object(S, 'open', counted, create=True)

    def test_sweep_and_free(self):
        M_, N_ = 40, 300
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        clock = _Clock()
        a = GraphletStore(d, ttl_disk_s=100, clock=clock)
        doc = T.doc_json('fork', 'graphlet')
        body = doc['results'][0]['graphlet']
        mine = []
        for i in range(M_):
            # distinct bodies: a comment line is no MGT, so vary the seed id instead
            r = copy.deepcopy(doc)
            r['results'][0]['graphlet'] = body.replace('fork', 'fork%03d' % i, 1) \
                if 'fork' in body.split('\n', 1)[0] else body
            mine.append(a.put_unparsed(r['results'][0], r, {'i': i}))
        # another process stores N entries after A opened (a real one)
        script = '''if 1:
            import sys, copy
            sys.path.insert(0, %r); sys.path.insert(0, %r)
            import traverse_testlib as T
            from metagraph.traverse import GraphletStore
            s = GraphletStore(sys.argv[1], clock=lambda: 1_000_000.0 + 1000)
            doc = T.doc_json('merge', 'graphlet')
            for i in range(int(sys.argv[2])):
                s.put(doc, {'i': i}, resident=False)
            ''' % (API, HERE)
        subprocess.run([sys.executable, '-c', script, d, str(N_)], check=True, timeout=600)
        foreign = [n for n in os.listdir(os.path.join(d, 'entries'))
                   if n[:-len('.json')] not in a._entries]
        self.assertEqual(N_, len(foreign))
        clock.t += 1000                     # A's own entries are past the disk TTL
        calls, patch = self.counting_open()
        with patch:
            got = a.sweep()
        self.assertEqual(sorted(mine), sorted(got['expired']))
        # each foreign entry file once, each expired entry a few times (its TTL check,
        # its tombstone): far below M x N
        self.assertLessEqual(calls['n'], N_ + 4 * M_)
        # the foreign body is kept, A's own are gone
        names = _bodies(d)
        self.assertEqual(1, len(names))
        # a later free in A reads no foreign entry file again
        h = a.put(T.doc_json('fork', 'graphlet'), {'x': 1})
        calls, patch = self.counting_open()
        with patch:
            a.free(h)
        self.assertLessEqual(calls['n'], 2)


# ======================================================================= O17

class TestO17DamagedBody(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        self.s = GraphletStore(self.d, max_ram_mb=0)
        self.doc = T.doc_json('fork', 'graphlet')
        self.h1 = self.s.put(self.doc, {'seeds': [{'sequence': 'ACGT'}]})
        self.path = self.s._body_path(self.s.get(self.h1).digest)
        with open(self.path, 'rb') as f:
            self.good = f.read()

    def damage(self, how):
        data = b'' if how == 'empty' else (self.good[:len(self.good) // 2] if how == 'torn'
                                           else b'\0' * len(self.good))
        with open(self.path, 'wb') as f:
            f.write(data)

    def test_a_put_repairs_the_body(self):
        for how in ('empty', 'torn', 'zeros'):
            with self.subTest(how=how):
                self.damage(how)
                s2 = GraphletStore(self.d, max_ram_mb=0)     # another process, a restart
                h2 = s2.put(self.doc, {'seeds': []})
                with open(self.path, 'rb') as f:
                    self.assertEqual(self.good, f.read())
                self.assertIsNotNone(s2.graphlet(h2))
                self.assertIsNotNone(self.s.graphlet(self.h1))   # the old entry too

    def test_a_read_expires_a_damaged_body(self):
        for how, read in (('empty', 'graphlet'), ('torn', 'body_text'),
                          ('zeros', 'standalone_text'), ('zeros', 'graphlet')):
            with self.subTest(how=how, read=read):
                self.setUp()
                self.damage(how)
                with self.assertRaises(UnknownHandle) as e:
                    getattr(self.s, read)(self.h1)
                self.assertTrue(e.exception.replayable)
                self.assertFalse(os.path.exists(self.path))    # deleted, not reused
                self.assertNotIn(self.h1, self.s)
                h3 = self.s.put(self.doc, {'seeds': []})         # a re-fetch repairs it
                self.assertIsNotNone(self.s.graphlet(h3))

    def test_the_tools_answer_unknown_handle_replayable(self):
        tools = GraphletTools(self.s, {'mini': M.FakeClient()})
        self.damage('empty')
        out = tools.graphlet_summary(self.h1)
        self.assertEqual('unknown_handle', out.get('error'), out)
        self.assertTrue(out.get('replayable'), out)

    def test_writes_are_flushed_before_the_rename(self):
        synced = []
        real = os.fsync
        with mock.patch.object(os, 'fsync', lambda fd: (synced.append(fd), real(fd))[1]):
            self.s.put(T.doc_json('merge', 'graphlet'), {'seeds': []})
        self.assertGreaterEqual(len(synced), 2)          # the body and the entry

    def test_a_good_body_is_hashed_once_per_file(self):
        hashed = []
        real = S._file_sha256
        with mock.patch.object(S, '_file_sha256', lambda p: (hashed.append(p), real(p))[1]):
            s2 = GraphletStore(self.d)
            s2.put(self.doc, {'seeds': []})
            s2.put(self.doc, {'seeds': []})
        self.assertEqual(1, len(hashed))


# ======================================================================= O11

class TestO11ToFastaIterators(unittest.TestCase):
    def test_one_shot_iterators(self):
        g = parse(comb_annotate(5))
        want = g.to_fasta('right', [0, 2, 4])
        self.assertTrue(want)
        for leaves in ((i for i in (0, 2, 4)), map(int, '024'), iter([0, 2, 4])):
            with self.subTest(kind=type(leaves).__name__):
                self.assertEqual(want, g.to_fasta('right', leaves))
        self.assertEqual('', g.to_fasta('right', (i for i in ())))
        self.assertEqual(want, export.to_fasta(g, 'right', (i for i in (0, 2, 4)),
                                               budget=LocalBudget()))
        b = LocalBudget()
        export.to_fasta(g, 'right', (i for i in (0, 2, 4)), budget=b)
        b2 = LocalBudget()
        export.to_fasta(g, 'right', [0, 2, 4], budget=b2)
        # the copy of the iterator's walks is charged on top of what a list costs
        self.assertGreater(b.usage()['memory_bytes'], b2.usage()['memory_bytes'])

    def test_both_arms_read_the_iterator_in_turn(self):
        # as before the regression: a one-shot iterator is exhausted by the left arm
        g = T.graphlet('fork')
        self.assertEqual(['left', 'right'], sorted(g.arms))
        by_list = g.to_fasta(None, [0])
        by_gen = g.to_fasta(None, (x for x in [0]))
        self.assertEqual(2, by_list.count('>'))
        self.assertEqual(1, by_gen.count('>'))
        self.assertTrue(by_gen.startswith('>left_0 '))


# ======================================================================= O35

class _CuttingServer:
    """Announces Content-Length 100000, sends 37 bytes and closes cleanly (a FIN): what a
    server's content timeout or a proxy does to a response in transfer."""

    def __init__(self):
        self.status = 200
        self.sock = socket.socket()
        self.sock.bind(('127.0.0.1', 0))
        self.sock.listen(8)
        self.port = self.sock.getsockname()[1]
        threading.Thread(target=self.serve, daemon=True).start()

    def serve(self):
        while True:
            try:
                conn, _ = self.sock.accept()
            except OSError:
                return
            data = b''
            conn.settimeout(5)
            try:
                while b'\r\n\r\n' not in data:
                    data += conn.recv(65536)
                head, _, rest = data.partition(b'\r\n\r\n')
                cl = 0
                for line in head.split(b'\r\n'):
                    if line.lower().startswith(b'content-length:'):
                        cl = int(line.split(b':')[1])
                while len(rest) < cl:
                    rest += conn.recv(65536)
                reason = b'OK' if self.status == 200 else b'Bad Request'
                conn.sendall(b'HTTP/1.1 %d %s\r\nContent-Type: application/json\r\n'
                             b'Content-Length: 100000\r\n\r\n{"results": [{"graphlet": {"x"'
                             % (self.status, reason))
                conn.shutdown(socket.SHUT_RDWR)
            except OSError:
                pass
            finally:
                conn.close()

    def close(self):
        self.sock.close()


class TestO35CutResponse(unittest.TestCase):
    def setUp(self):
        self.server = _CuttingServer()
        self.addCleanup(self.server.close)
        self.client = TraverseClient('127.0.0.1', self.server.port, feature_level=4,
                                     timeout=30)

    def test_the_client_raises_a_connection_error(self):
        for status in (200, 400):
            with self.subTest(status=status):
                self.server.status = status
                with self.assertRaises(ConnectionError) as e:
                    self.client.traverse_raw({'seeds': [{'sequence': 'ACGT' * 9}]})
                self.assertIsInstance(e.exception, OSError)
                self.assertIsInstance(e.exception.__cause__, http.client.IncompleteRead)
                self.assertIn('cut in transfer', str(e.exception))
                self.assertIn('HTTP %d' % status, str(e.exception))
                # the server was reached: never "unreachable", and what a retry does (the
                # review of the P3 fixes)
                self.assertIn('the server was reached and may have run the request (a retry '
                              'sends it again)', str(e.exception))
        # an attempt: cancel it and judge it, as the docstring says
        self.server.status = 200
        with self.assertRaises(ConnectionError) as e:
            self.client.traverse_raw({'seeds': [{'sequence': 'ACGT' * 9}],
                                      'attempt_id': 'att-1'})
        self.assertIn("cancel attempt 'att-1', then judge it with release_verdict()",
                      str(e.exception))

    def test_the_tools_answer_backend_unreachable_and_report_usage(self):
        usage = []
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        for limits in (None, ToolLimits(on_usage=lambda tool, u: usage.append(tool))):
            tools = GraphletTools(GraphletStore(d), {'mini': self.client}, local_limits=limits)
            for status in (200, 400):
                self.server.status = status
                for name, call in (
                        ('traverse_fetch', lambda: tools.traverse_fetch(
                            seed={'sequence': 'ACGT' * 9})),
                        ('traverse_capabilities', tools.traverse_capabilities),
                        ('traverse_resolve', lambda: tools.traverse_resolve(
                            sequence='ACGT' * 8))):
                    with self.subTest(limits=limits is not None, status=status, tool=name):
                        out = call()
                        self.assertEqual('backend_unreachable', out.get('error'), out)
                        self.assertIn('cut in transfer', out['message'])
                        # the code of a transport failure, not its words: the server was
                        # reached and may have run the request
                        self.assertNotIn('could not be reached', out['message'])
                        self.assertIn('the server was reached and may have run the request '
                                      '(a retry sends it again)', out['message'])
        self.assertIn('traverse_fetch', usage)


# ======================================================================= G15

class TestG15AtomicExport(unittest.TestCase):
    def test_a_failed_re_export_leaves_the_earlier_file(self):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        st = GraphletStore(d)
        h = st.put(T.doc_json('switch_chain', 'graphlet'), T.doc_json('switch_chain', 'request'))
        real = mcp_tools._write_all

        def half_then_enospc(fd, data, *a, **k):
            real(fd, data[:len(data) // 2])
            raise OSError(28, 'No space left on device')
        for fmt in ('mgt', 'fasta', 'gfa', 'json'):
            for label, tools in (('unbudgeted', GraphletTools(st)),
                                 ('limits', GraphletTools(st, local_limits=ToolLimits()))):
                with self.subTest(format=fmt, tools=label):
                    first = tools.graphlet_export(h, format=fmt)
                    path = first['path']
                    with open(path, 'rb') as f:
                        full = f.read()
                    self.assertEqual(first['bytes'], len(full))
                    with mock.patch.object(mcp_tools, '_write_all', half_then_enospc):
                        again = tools.graphlet_export(h, format=fmt)
                    if fmt == 'json' and label == 'limits':
                        # written by json.dump through a buffered file, not _write_all
                        self.assertNotIn('error', again)
                    else:
                        self.assertEqual('io_error', again.get('error'), again)
                    with open(path, 'rb') as f:
                        self.assertEqual(full, f.read())
                    self.assertEqual([], [n for n in os.listdir(os.path.dirname(path))
                                          if n.startswith('.export-')])
                    os.unlink(path)

    def test_the_file_keeps_its_permissions(self):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        st = GraphletStore(d)
        h = st.put(T.doc_json('switch_chain', 'graphlet'), T.doc_json('switch_chain', 'request'))
        tools = GraphletTools(st)
        path = tools.graphlet_export(h, format='fasta')['path']
        umask = os.umask(0)
        os.umask(umask)
        self.assertEqual(0o666 & ~umask, os.stat(path).st_mode & 0o777)
        os.chmod(path, 0o640)
        tools.graphlet_export(h, format='fasta')
        self.assertEqual(0o640, os.stat(path).st_mode & 0o777)


# ======================================================================= L8

def _strict_cost(cc, path='strategy.labels.change_cost'):
    """A mirror of traverse.cpp parse_cost (a Strict object): -> None, or the server's
    error."""
    if not isinstance(cc, dict):
        return '%s: expected an object' % path
    model = cc.get('model', 'forbid')
    if model not in ('forbid', 'constant', 'table'):
        return '%s.model: unknown value' % path
    reads = {'forbid': {'model'}, 'constant': {'model', 'value'},
             'table': {'model', 'default', 'entries'}}[model]
    if model == 'constant' and not isinstance(cc.get('value'), (int, float)):
        return '%s.value is required for model constant' % path
    if model == 'table' and not isinstance(cc.get('entries'), list):
        return '%s.entries: expected an array' % path
    for k in cc:
        if k not in reads:
            return "%s: unknown field '%s'" % (path, k)
    return None


def _strict_request(req):
    st = req['strategy']
    for sec in ('labels', 'branching', 'bounds', 'output'):
        if sec in st and not isinstance(st[sec], dict):
            return 'strategy.%s: expected an object' % sec
    lab = st.get('labels') or {}
    if 'change_cost' in lab:
        why = _strict_cost(lab['change_cost'])
        if why is not None:
            return why
        # walker.cpp: a finite change cost conflicts with labels.mode "annotate"
        if lab.get('mode') == 'annotate' and lab['change_cost'].get('model', 'forbid') != 'forbid':
            return 'strategy.labels.change_cost does not apply in labels.mode "annotate"'
    return None


class TestL8Overrides(unittest.TestCase):
    TABLE = {'model': 'table', 'default': 'forbid',
             'entries': [['C', 'A', 0.5], ['C', 'B', 0.5], ['C', 'D', 0.5]]}

    def table_graphlet(self):
        g = T.graphlet('switch_chain')
        g.envelope['strategy']['labels']['change_cost'] = {
            'model': 'table', 'default': 1.0, 'entries': [['A', 'B', 1.0]]}
        return g

    def cases(self):
        yield 'constant -> table', T.graphlet('switch_chain'), {'labels': {'change_cost': self.TABLE}}
        yield 'constant -> forbid', T.graphlet('switch_chain'), \
            {'labels': {'change_cost': {'model': 'forbid'}}}
        yield 'constant -> constant', T.graphlet('switch_chain'), \
            {'labels': {'change_cost': {'value': 0.5}}}
        yield 'table -> constant', self.table_graphlet(), \
            {'labels': {'change_cost': {'model': 'constant', 'value': 0.5}}}
        yield 'table -> table', self.table_graphlet(), \
            {'labels': {'change_cost': {'default': 2.0}}}

    def test_a_model_change_replaces_the_cost(self):
        want = {
            'constant -> table': self.TABLE,
            'constant -> forbid': {'model': 'forbid'},
            'constant -> constant': {'model': 'constant', 'value': 0.5},
            'table -> constant': {'model': 'constant', 'value': 0.5},
            'table -> table': {'model': 'table', 'default': 2.0, 'entries': [['A', 'B', 1.0]]},
        }
        for name, g, ov in self.cases():
            with self.subTest(name):
                for budget in (None, LocalBudget()):
                    kw = {} if budget is None else {'budget': budget}
                    req = g.next_request('right', [1], **kw, **copy.deepcopy(ov))
                    self.assertEqual(want[name], req['strategy']['labels']['change_cost'])
                    self.assertIsNone(_strict_request(req))
                reqs = g.next_requests('right', [1], **copy.deepcopy(ov))
                self.assertIsNone(_strict_request(reqs[0]))

    def test_the_rebuild_follows_the_new_model(self):
        g = T.graphlet('switch_chain')
        req = g.next_request('right', [1], labels={'change_cost': {'model': 'forbid'}})
        self.assertEqual([], req['strategy']['labels']['extra'])

    def test_none_for_an_object_section(self):
        for ov, field in (({'labels': None}, 'strategy.labels'),
                          ({'labels': {'change_cost': None}}, 'strategy.labels.change_cost'),
                          ({'branching': None}, 'strategy.branching'),
                          ({'bounds': None}, 'strategy.bounds')):
            with self.subTest(field=field):
                with self.assertRaises(ValueError) as e:
                    T.graphlet('switch_chain').next_request('right', [1], **ov)
                self.assertTrue(str(e.exception).startswith(field + ' is an object'),
                                str(e.exception))
        # annotate mode too (no rebuild there)
        with self.assertRaises(ValueError):
            T.graphlet('annotate').next_request('right', [0], bounds=None)

    def test_a_field_the_model_does_not_read(self):
        for cc, field in (({'model': 'constant', 'value': 1, 'entries': []}, 'entries'),
                          ({'model': 'forbid', 'value': 1}, 'value'),
                          ({'model': 'table', 'entries': [], 'value': 1}, 'value')):
            with self.subTest(cc=cc):
                with self.assertRaises(ValueError) as e:
                    T.graphlet('switch_chain').next_request(
                        'right', [1], labels={'change_cost': cc})
                self.assertIn('labels.change_cost.%s is not a field of model' % field,
                              str(e.exception))

    # annotate mode has no rebuild, so its merged cost went out unchecked and the server
    # refused it (the review of the P3 fixes): -> the field the ValueError names
    ANNOTATE_REFUSED = (
        ({'model': 'forbid', 'value': 1}, 'labels.change_cost.value is not a field of model'),
        ({'model': 'constant', 'value': 1}, "labels.change_cost.model 'constant' does not "
                                            'apply in labels.mode "annotate"'),
        ({'model': 'constant'}, "labels.change_cost.model 'constant' does not apply"),
        ({'model': 'table', 'entries': [['A', 'B', 1.0]]},
         "labels.change_cost.model 'table' does not apply"),
        ({'value': 1}, 'labels.change_cost.value is not a field of model'),
        ('forbid', 'labels.change_cost is an object'),
    )

    def test_annotate_mode(self):
        g = T.graphlet('annotate')
        self.assertEqual('annotate', g.mode)
        for cc, why in self.ANNOTATE_REFUSED:
            for budget in (None, LocalBudget()):
                with self.subTest(cc=cc, budgeted=budget is not None):
                    kw = {} if budget is None else {'budget': budget}
                    with self.assertRaises(ValueError) as e:
                        g.next_request('right', [0], labels={'change_cost': copy.deepcopy(cc)},
                                       **kw)
                    self.assertIn(why, str(e.exception))
                    with self.assertRaises(ValueError):
                        g.next_requests('right', [0], labels={'change_cost': copy.deepcopy(cc)},
                                        **kw)
        # what the server accepts is still built, as before
        plain = g.next_request('right', [0])
        for ov in ({}, {'labels': {'change_cost': {'model': 'forbid'}}},
                   {'labels': {'change_cost': {}}}):
            with self.subTest(ov=ov):
                req = g.next_request('right', [0], **copy.deepcopy(ov))
                self.assertEqual(plain, req)
                self.assertIsNone(_strict_request(req))

    def test_the_tool_answers_bad_argument(self):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        st = GraphletStore(d)
        h = st.put(T.doc_json('switch_chain', 'graphlet'), T.doc_json('switch_chain', 'request'))
        tools = GraphletTools(st, {'mini': M.FakeClient()})
        out = tools.traverse_continue(h, 'right', 1, overrides={'labels': None},
                                      execute=False)
        self.assertEqual('bad_argument', out.get('error'), out)
        out = tools.traverse_continue(
            h, 'right', 1, overrides={'labels': {'change_cost': self.TABLE}}, execute=False)
        self.assertNotIn('error', out)
        ha = st.put(T.doc_json('annotate', 'graphlet'), T.doc_json('annotate', 'request'))
        out = tools.traverse_continue(
            ha, 'right', 0, overrides={'labels': {'change_cost': {'model': 'constant',
                                                                  'value': 1}}},
            execute=False)
        self.assertEqual('bad_argument', out.get('error'), out)
        self.assertIn('labels.mode "annotate"', out['message'])

    def test_the_cli_parses_every_built_request(self):
        binary = os.environ.get('METAGRAPH_CLI') or next(
            (p for p in (os.path.join(T.REPO, 'build', 'metagraph'),
                         os.path.join(T.REPO, 'build_debug', 'metagraph'))
             if os.path.exists(p)), None)
        if binary is None or not os.path.exists(binary):
            self.skipTest('needs a metagraph binary ($METAGRAPH_CLI or build/metagraph)')
        G = T.fixture_script()
        root = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, root, True)

        def run(req, index='chain'):
            if not os.path.isdir(os.path.join(root, index)):
                G.build_index(index, binary, root)
            path = os.path.join(root, 'req.json')
            with open(path, 'w') as f:
                json.dump({'seeds': req['seeds'], 'strategy': req['strategy']}, f)
            return subprocess.run([binary, 'traverse', '--index-name', G.CLI_INDEX_NAME,
                                   '-i', os.path.join(root, index, 'graph.dbg'),
                                   '-a', G.annotation_of(index, root), path],
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=300)
        for name, g, ov in self.cases():
            with self.subTest(name):
                p = run(g.next_request('right', [1], **copy.deepcopy(ov)))
                self.assertEqual(0, p.returncode, p.stderr.decode()[-400:])
        # annotate mode: what the library builds parses, and what it refuses the server
        # refuses too (built here with the check bypassed) -- the library mirrors it
        ga = T.graphlet('annotate')
        index = G.CLI_FIXTURES['annotate']['index']
        for ov in ({}, {'labels': {'change_cost': {'model': 'forbid'}}}):
            with self.subTest(annotate=ov):
                p = run(ga.next_request('right', [0], **copy.deepcopy(ov)), index)
                self.assertEqual(0, p.returncode, p.stderr.decode()[-400:])
        for cc, _ in self.ANNOTATE_REFUSED:
            with self.subTest(annotate_refused=cc):
                with mock.patch.object(ops, '_check_annotate_cost', lambda c: None):
                    req = ga.next_request('right', [0], labels={'change_cost': copy.deepcopy(cc)})
                p = run(req, index)
                self.assertNotEqual(0, p.returncode, 'the server accepted %r' % (cc,))
                self.assertIn('change_cost', p.stderr.decode())


# ======================================================================= L7

def _reference_note_reach(g, cost, mine, budget, dropped):
    """The note's reach as it was computed: one switch search of the whole pool per
    continued label -> (the reaching labels of each of the first 8 left-out labels, any
    left-out label reached)."""
    pool = [l.name for l in g.labels]
    reach_of = [(p, ops._switch_reach(cost, [p.name], pool, budget - x) or {})
                for p, x in mine]
    via = [[p for p, r in reach_of if l.name in r] for l in dropped]
    return via[:8], any(via)


class _G:
    def __init__(self, names):
        self.labels = [Label(id=i, kind='column', column=i, seq_id=None, name=n)
                       for i, n in enumerate(names)]


class TestL7LeftOutNote(unittest.TestCase):
    def test_the_reach_is_the_full_search_s(self):
        rng = random.Random(70707)
        vals = [0, 0.1, 0.25, 0.5, 1, 1.5, 2, 3]
        for cell in range(1500):
            n = rng.randrange(2, 30)
            names = [rng.choice(['L%d' % i, 'L%d' % i, 'dup', 'L%d' % (i // 2)])
                     for i in range(n)]
            g = _G(names)
            model = rng.choice(['table', 'table', 'constant'])
            if model == 'table':
                cost = {'model': 'table', 'default': rng.choice(['forbid'] + vals),
                        'entries': [[rng.choice(names), rng.choice(names), rng.choice(vals)]
                                    for _ in range(rng.randrange(0, 2 * n))]}
            else:
                cost = {'model': 'constant', 'value': rng.choice(vals)}
            mine = [(l, rng.choice(vals)) for l in rng.sample(g.labels, rng.randrange(1, n))]
            dropped = rng.sample(g.labels, rng.randrange(1, n))
            budget = rng.choice([0.5, 1, 2, 3, 5])
            want = _reference_note_reach(g, cost, mine, budget, dropped)
            got = ops._left_out_reach(g, cost, mine, budget, dropped)
            self.assertEqual(want, got, (cell, cost, mine, dropped, budget))

    def continuation(self, n_labels, n_cont, cost, losses=None, names=None):
        g = T.graphlet('switch_chain')
        base = len(g.labels)
        for k in range(n_labels - base):
            g.labels.append(Label(id=base + k, kind='column', column=1000 + k, seq_id=None,
                                  name=(names or 'lab%05d').__mod__(k)))
        g.cache.pop('name_counts', None)
        cont = g.labels[:n_cont]
        losses = losses or tuple([1.5] + [0.0] * (n_cont - 1))
        g.envelope['strategy']['labels']['loss_budget'] = 2.0
        g.envelope['strategy']['labels']['change_cost'] = cost
        orig = ops._continuation_seeds

        def stub(g_, a_, leaves, *rest):
            conts, _ = orig(g_, a_, leaves, *rest)
            c = dataclasses.replace(conts[0], labels=list(cont), losses=losses)
            return [c], [c.as_seed()]
        return g, mock.patch.object(ops, '_continuation_seeds', stub)

    def test_a_label_name_holding_the_words_claims_no_reach(self):
        # the clause was decided by a substring test of the note's text, so a label NAMED
        # "... reachable for ..." claimed a reach no chain makes
        g, patch = self.continuation(12, 2, {'model': 'constant', 'value': 3},
                                     losses=(1.5, 3), names='reachable for %d')
        with patch:
            req = g.next_request('right', [1])
        note = [x for x in req.notes if x.startswith('labels.extra leaves out')]
        self.assertEqual(1, len(note))
        self.assertIn("'reachable for 0'", note[0])
        self.assertNotIn('one uninterrupted walk could still', note[0])

    def test_the_note_names_what_reaches_a_label(self):
        g, patch = self.continuation(40, 20, {'model': 'table', 'default': 1.0, 'entries': []})
        with patch:
            req = g.next_request('right', [1])
        note = [x for x in req.notes if x.startswith('labels.extra leaves out')][0]
        self.assertIn('reachable for', note)
        self.assertIn('and 12 more', note)
        self.assertTrue(note.endswith('from a label at a lower loss'))

    def test_scale_and_charge(self):
        for cost in ({'model': 'table', 'default': 1.0, 'entries': []},
                     {'model': 'constant', 'value': 1.0}):
            with self.subTest(model=cost['model']):
                g, patch = self.continuation(2000, 1000, cost)
                b = LocalBudget()
                with patch:
                    t = time.process_time()
                    req = g.next_request('right', [1], budget=b)
                    dt = time.process_time() - t
                self.assertEqual(1000, len([x for x in req.left_out if x['why'] == 'unreachable']))
                # 1.39 s and 18,027,203 lwu (table) with a search of the whole pool per
                # continued label; well under a second and 1 M lwu now
                self.assertLess(dt, 2.0)
                self.assertLess(b.usage()['work_units'], 1_000_000)
                # the account bounds the traced peak
                g, patch = self.continuation(2000, 1000, cost)
                b = LocalBudget()
                with patch:
                    peak, _ = _traced(lambda: g.next_request('right', [1], budget=b))
                self.assertLessEqual(peak, b.peak_bytes)


# ======================================================================= L3

def _wide_split(w):
    """Annotate: a root with labels 0..w-1 split into a 1-label leaf 'C' {0} and a
    w-label sibling 'G' (repro U15-05)."""
    steps = 3
    out = ['H mgt 1 5 basic $ACGT a k k %d 10 0 walk * * 0123456789abcdef' % max(w, 1000),
           'S aaaaaaaaaaaaaaaa 10 6 0 ACGTACGTAC']
    out += ['L c %d * 0 lab%d' % (i, i) for i in range(w)]
    out += ['O c c c i',
            'A r c 2 p 0 0 1 2 0 %s * 0 * 3 0 2 1 0 %d' % (COUNTERS % (steps, steps, steps), steps),
            'B 0 2 2 3 1 %d 1 0 1 0 0 0 0 .' % steps,
            'G * 0 1 0-%d * * * 0 * A' % (w - 1), 'P 0 1 %d !' % w,
            'G 0 1 1 0 * * * * * C', 'P 1 2 1 0', 'T D .',
            'G 0 1 1 ! * * * * * G', 'P 1 2 %d !' % w, 'T D .']
    out.append('Z %d' % (len(out) + 1))
    return '\n'.join(out) + '\n'


def _reference_reasons(g, a, chain, at, removed):
    """The split check as it was: membership scans of the arrays per removed label."""
    segs = a.segments
    out = {}
    for l in sorted(removed):
        for sid in chain:
            s = segs[sid]
            if s.from_bp == at and len(s.parents) == 1:
                parent = segs[s.parents[0]]
                took = sorted({x for c in parent.children if c != sid
                               for x in (segs[c].walk[:1] if segs[c].walk
                                         else segs[c].first_base or '')
                               if l in segs[c].entry})
                if l in parent.end and took:
                    out[l] = took
                break
    return out


class TestL3SupportChanges(unittest.TestCase):
    def test_the_split_check_is_the_scan_s(self):
        for g in [T.graphlet(n) for n in T.DOCUMENTS if n not in T.HAND_MADE] + \
                [parse(_wide_split(50)), parse(comb_annotate(30))]:
            for side, a in g.arms.items():
                for p in ops.derive.paths(a):
                    chain = ops.derive.chain(a, p.leaf)
                    for ch in ops.support_changes(g, side, p.id):
                        removed = {g.labels.index(l) if False else l.id for l in ch.removed}
                        want = _reference_reasons(g, a, chain, ch.at_bp, removed)
                        got = {r['label']['ref']: r['took'] for r in ch.reasons
                               if r.get('why') == 'split'}
                        self.assertEqual({g.labels[l].ref: t for l, t in want.items()}, got)

    def test_scale_and_charge(self):
        w = 20000
        g = parse(_wide_split(w))
        t = time.process_time()
        un = ops.support_changes(g, 'right', 0)
        dt = time.process_time() - t
        self.assertLess(dt, 2.0)          # 3.9 s here (6.3 s in the review) before
        b = LocalBudget()
        self.assertEqual([c.reasons for c in un],
                         [c.reasons for c in ops.support_changes(g, 'right', 0, budget=b)])
        self.assertEqual(w - 1, sum(len(c.removed) for c in un))
        # the label ids the split check reads (one bisection of each set per removed
        # label) are charged: at least the ids / 32 the review asked for
        self.assertGreaterEqual(b.usage()['work_units'], 2 * w // 32 + w)
        g = parse(_wide_split(2000))
        b = LocalBudget()
        peak, _ = _traced(lambda: ops.support_changes(g, 'right', 0, budget=b))
        self.assertLessEqual(peak, b.peak_bytes)


# ======================================================================= L5

def _comb_wide(n, width):
    """comb_annotate(n) with |width| labels: every spine segment carries all of them,
    every leaf all but one (repro X-EFFICIENCY-05)."""
    steps = 1 + 2 * n
    out = ['H mgt 1 5 basic $ACGT a k k 1000 10 0 walk * * 0123456789abcdef',
           'S aaaaaaaaaaaaaaaa 10 6 0 ACGTACGTAC']
    out += ['L c %d * %d lab%05d' % (i, i, i) for i in range(width)]
    out += ['O c c c i',
            'A r c %d p 0 0 1 %d 0 %s * 0 * %d 0 %d %d 0 %d'
            % (n + 1, width, COUNTERS % (steps, steps, steps), 2 * n + 1, n + 1, n, steps),
            'B 0 2 2 3 1 %d %d 0 %d 0 0 0 0 .' % (steps, n, n),
            'G * 0 1 0-%d * * * %s * A' % (width - 1, '0' if n else '*'),
            'P 0 1 %d !' % width]
    for i in range(1, n + 1):
        parent = 0 if i == 1 else 2 * (i - 1)
        out += ['G %d %d 1 0-%d * * * * * C' % (parent, i, width - 2),
                'P %d %d %d !' % (i, i + 1, width - 1), 'T D .',
                'G %d %d 1 ! * * * %s * G' % (parent, i, '0' if i < n else '*'),
                'P %d %d %d !' % (i, i + 1, width)]
        if i == n:
            out.append('T D .')
    out.append('Z %d' % (len(out) + 1))
    return '\n'.join(out) + '\n'


class TestL5AnnotateWalkClaims(unittest.TestCase):
    def test_claims_and_charge(self):
        n, width = 640, 64
        g = parse(_comb_wide(n, width))
        a = g.arms['right']
        e = sum(len(v) for v in ops._annotate_route_ends(a)[1].values())
        un = ops.walks(g, 'right', by='id')
        reads = {'n': 0}
        real = ops._annotate_route_ends

        def counted(arm):
            reads['n'] += 1
            return real(arm)
        g2 = parse(_comb_wide(n, width))
        b = LocalBudget()
        with mock.patch.object(ops, '_annotate_route_ends', counted):
            bw = ops.walks(g2, 'right', by='id', budget=b)
        self.assertEqual([(w.path_id, [c.as_dict() for c in w.claims]) for w in un],
                         [(w.path_id, [c.as_dict() for c in w.claims]) for w in bw])
        chunks = -(-len(un) // ops._WALK_CHUNK)
        self.assertGreater(chunks, 10)
        # the route ends are read a fixed number of times, not once per chunk
        self.assertLess(reads['n'], chunks)
        def claims_charge(ids):
            b0, b1 = LocalBudget(), LocalBudget()
            ops.walks_at(g2, 'right', ids, with_claims=False, budget=b0)
            ops.walks_at(g2, 'right', ids, with_claims=True, budget=b1)
            got = ops.walks_at(g2, 'right', ids, with_claims=True)
            return (b1.usage()['work_units'] - b0.usage()['work_units'],
                    sum(len(w.claims) for w in got))
        ids = [w.path_id for w in un]
        w1, c1 = claims_charge(ids[:ops._WALK_CHUNK])
        w2, c2 = claims_charge(ids)
        # each further chunk adds its claims' rows only: a pass over all e route ends per
        # chunk added (chunks - 1) x e more before
        self.assertLessEqual(w2 - w1, ops.W_ROW * (c2 - c1))
        self.assertGreater((chunks - 1) * e, w2 - w1)


# ======================================================================= L12

class TestL12Cuts(unittest.TestCase):
    def graphlets(self):
        out = [T.graphlet(n) for n in T.DOCUMENTS if n not in T.HAND_MADE]
        out += [parse(T.doc_text(n)) for n in T.HAND_MADE]
        out += [parse(comb_annotate(n)) for n in (0, 1, 2, 7, 60)]
        out += [parse(_comb_wide(20, 6))]
        return out

    def test_every_cut_is_the_chain_read_s(self):
        for g in self.graphlets():
            for side, a in g.arms.items():
                depths = {0, 1}
                for s in a.segments:
                    depths.update((s.from_bp, s.end_bp, s.end_bp + 1, max(0, s.end_bp - 1)))
                for depth in sorted(depths):
                    for leaf, m, info in ops._cuts(g, side, depth, {}):
                        k = None
                        if m:
                            k = leaf
                            while a.segments[k].from_bp >= m:
                                k = a.segments[k].parents[0]
                        want = ops._cut_info(g, a, k, m)
                        self.assertEqual(want[0], info[0])
                        self.assertEqual(set(want[1]), set(info[1]))
                        self.assertEqual(want[2], info[2])

    def test_budgeted_cuts_are_the_same_and_priced_alike(self):
        for g in self.graphlets():
            for side in g.arms:
                depth = g.arms[side].complete_to_bp
                b = LocalBudget()
                with b.scope('cuts'):
                    got = ops._cuts(g, side, depth, {}, b)
                self.assertEqual([(x, m, i) for x, m, i in ops._cuts(g, side, depth, {})],
                                 [(x, m, i) for x, m, i in got])
                t = LocalBudget()
                with t.scope('cuts'):
                    ops._price_cuts(t, g, side, depth)
                self.assertEqual(b.used_work, t.used_work)
                self.assertEqual(b.peak_bytes, t.peak_bytes)

    def test_scale(self):
        for n in (2000, 4000):
            a, b = parse(comb_annotate(n)), parse(comb_annotate(n))
            t = time.process_time()
            ops.compare(a, b, mode='walks')
            dt = time.process_time() - t
            self.assertLess(dt, 2.5)     # 3.1 s / 12.9 s at 2,000 / 4,000 in the review
        a, b = parse(comb_annotate(4000)), parse(comb_annotate(4000))
        bud = LocalBudget()
        c = ops.compare(a, b, mode='walks', budget=bud)
        self.assertIsNone(c.local_stop)
        est = ops.compare_cost(a, b, mode='walks')
        self.assertLessEqual(est['work_units']['at_least'], bud.usage()['work_units'])
        # the cuts are charged per segment once: far below the chains' cut_cost per cut
        arm = a.arms['right']
        cost = ops.derive.cut_cost(arm, a.mode)
        per_cut = sum(cost[x.leaf] for x in ops.derive.paths(arm))
        self.assertLess(bud.usage()['work_units'], per_cut)


if __name__ == '__main__':
    unittest.main()
