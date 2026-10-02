"""GraphletStore: opaque entry handles separate from body identity, body dedup by
sha256, LRU by count and RAM, both TTLs on a fake clock, tombstones for replay, the
hard size limit, rejection of truncated bodies, save/load, persistence."""

import copy
import os
import re
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402

from metagraph.traverse import (  # noqa: E402
    GraphletFormatError, GraphletStore, StoreLimitExceeded, UnknownHandle,
)


class Clock:
    def __init__(self):
        self.t = 1000.0

    def __call__(self):
        return self.t


def response(name):
    return T.doc_json(name, 'graphlet')


class TestStore(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.clock = Clock()
        self.store = GraphletStore(self.tmp.name, max_handles=2, ttl_ram_s=100,
                                   ttl_disk_s=1000, clock=self.clock)

    def tearDown(self):
        self.tmp.cleanup()

    def bodies(self):
        return os.listdir(os.path.join(self.tmp.name, 'bodies'))

    def test_handles_are_opaque_and_bodies_deduplicated(self):
        req1 = {'seeds': [{'sequence': 'X'}], 'strategy': {'labels': {'loss_budget': 1}}}
        req2 = {'seeds': [{'sequence': 'X'}], 'strategy': {'labels': {'loss_budget': 2}}}
        h1 = self.store.put(response('fork'), req1)
        h2 = self.store.put(response('fork'), req2)    # byte-identical body, other request
        self.assertNotEqual(h1, h2)
        for h in (h1, h2):
            self.assertRegex(h, r'^g_[0-9a-f]{12}$')
        self.assertEqual(1, len(self.bodies()))
        self.assertEqual(req2, self.store.get(h2).request)
        self.store.free(h1)
        self.assertEqual(1, len(self.bodies()))         # still held by h2
        self.store.free(h2)
        self.assertEqual([], self.bodies())
        with self.assertRaises(UnknownHandle) as e:
            self.store.get(h1)
        self.assertFalse(e.exception.replayable)

    def test_lazy_parse_and_lru(self):
        hs = [self.store.put(response(n), {}) for n in ('fork', 'merge', 'linear')]
        self.assertFalse(self.store.in_ram(hs[0]))      # max_handles 2
        self.assertTrue(self.store.in_ram(hs[2]))
        g = self.store.graphlet(hs[0])                  # parsed again from the spool
        self.assertTrue(g.has_envelope)
        self.assertEqual(T.doc_text('fork'), g.dump())
        self.assertTrue(self.store.in_ram(hs[0]))
        self.assertFalse(self.store.in_ram(hs[1]))      # least recently used out

    def test_ram_budget(self):
        sizes = {n: T.graphlet(n).memory_bytes() for n in ('fork', 'merge')}
        # room for either graphlet, not for both: the least recently used goes
        one = (max(sizes.values()) + min(sizes.values()) // 2) / (1 << 20)
        small = GraphletStore(os.path.join(self.tmp.name, 's'), max_ram_mb=one,
                              clock=self.clock)
        a = small.put(response('fork'), {})
        b = small.put(response('merge'), {})
        self.assertFalse(small.in_ram(a))
        self.assertTrue(small.in_ram(b))
        # a graphlet larger than the whole budget is never resident, not even the
        # newest: it is parsed on demand
        tiny = GraphletStore(os.path.join(self.tmp.name, 't'), max_ram_mb=0.00001,
                             clock=self.clock)
        c = tiny.put(response('merge'), {})
        self.assertFalse(tiny.in_ram(c))
        self.assertEqual(T.doc_text('merge'), tiny.graphlet(c).dump(envelope=False))
        self.assertFalse(tiny.in_ram(c))
        self.assertEqual(0, tiny._ram_bytes)

    def test_a_spooled_body_keeps_one_delivery_value(self):
        h = self.store.put(response('merge'), {}, resident=False, delivery='spooled')
        self.assertFalse(self.store.in_ram(h))
        g = self.store.graphlet(h)
        self.assertEqual('spooled', g.outcome.delivery)
        self.assertEqual('spooled', g.seed_summary['outcome']['delivery'])
        self.assertEqual('spooled', g.summary()['outcome']['delivery'])
        again = GraphletStore(self.tmp.name, clock=self.clock)
        self.assertEqual('spooled', again.graphlet(h).outcome.delivery)
        # the body in the spool is the server's, unchanged
        with open(os.path.join(self.tmp.name, 'bodies', self.bodies()[0]),
                  encoding='utf-8', newline='') as f:
            self.assertEqual(T.doc_text('merge'), f.read())

    def test_parent_and_tombstones_and_the_cursor_secret(self):
        parent = {'handle': 'g_000000000000', 'arm': 'right', 'walk': 1, 'overlap_bp': 20}
        h = self.store.put(response('fork'), {'seeds': []}, parent=parent,
                           derived_from=parent['handle'])
        again = GraphletStore(self.tmp.name, clock=self.clock)
        self.assertEqual(parent, again.get(h).parent)
        listed = [e for e in again.list() if e['handle'] == h][0]
        self.assertEqual((parent, 'g_000000000000'), (listed['parent'], listed['derived_from']))
        # tombstones go ttl_tomb_s after the expiry
        store = GraphletStore(self.tmp.name, ttl_disk_s=1000, ttl_tomb_s=5000,
                              clock=self.clock)
        self.clock.t += 2000
        self.assertEqual([h], store.sweep()['expired'])
        self.assertTrue(store._unknown(h).replayable)
        self.clock.t += 4000
        self.assertEqual([], store.sweep()['tombstones_removed'])
        self.clock.t += 2000
        self.assertEqual([h], store.sweep()['tombstones_removed'])
        self.assertEqual([], os.listdir(os.path.join(self.tmp.name, 'tombstones')))
        self.assertFalse(store._unknown(h).replayable)
        # one cursor secret per spool, private, the same after a restart
        key = store.cursor_secret()
        self.assertGreaterEqual(len(key), 16)
        self.assertEqual(key, GraphletStore(self.tmp.name).cursor_secret())
        mode = os.stat(os.path.join(self.tmp.name, 'cursor.key')).st_mode & 0o777
        self.assertEqual(0o600, mode)
        # a handle is never joined into a spool path unless it is one
        with self.assertRaises(UnknownHandle) as e:
            store.get('../entries/x')
        self.assertFalse(e.exception.replayable)

    def test_ttls_with_a_fake_clock(self):
        req = {'seeds': [{'sequence': 'ACGT'}], 'strategy': {}}
        h = self.store.put(response('fork'), req, source='uhgg')
        self.clock.t += 150
        out = self.store.sweep()
        self.assertEqual([h], out['ram_dropped'])
        self.assertFalse(self.store.in_ram(h))
        self.store.graphlet(h)                          # back from disk, refreshes
        self.clock.t += 999
        self.assertEqual([], self.store.sweep()['expired'])
        self.clock.t += 1001
        self.assertEqual([h], self.store.sweep()['expired'])
        with self.assertRaises(UnknownHandle) as e:
            self.store.graphlet(h)
        self.assertTrue(e.exception.replayable)
        self.assertEqual(req, e.exception.request)
        self.assertIn('traverse_fetch(replay=', str(e.exception))
        self.assertEqual([], self.bodies())

    def test_use_is_remembered_across_restarts(self):
        h = self.store.put(response('fork'), {})
        self.clock.t += 600
        self.store.get(h)                               # a use after 600 s
        self.clock.t += 600                             # 1200 s after the put
        again = GraphletStore(self.tmp.name, ttl_disk_s=1000, clock=self.clock)
        self.assertEqual([], again.sweep()['expired'])  # 600 s since the last use

    def test_expiry_on_access(self):
        h = self.store.put(response('fork'), {})
        self.clock.t += 2000
        with self.assertRaises(UnknownHandle):
            self.store.get(h)

    def test_truncated_and_oversized_bodies_are_never_stored(self):
        resp = copy.deepcopy(response('merge'))
        r = resp['results'][0]
        r['graphlet'] = r['graphlet'][:r['graphlet'].rindex('R ')]
        del r['graphlet_bytes'], r['graphlet_lines']
        with self.assertRaises(GraphletFormatError):
            self.store.put(resp, {})
        self.assertEqual([], self.bodies())
        tiny = GraphletStore(os.path.join(self.tmp.name, 't'), max_body_mb=0.0001)
        with self.assertRaises(StoreLimitExceeded):
            tiny.put(response('merge'), {})

    def test_save_load_and_persistence(self):
        h = self.store.put(response('annotate'), {'seeds': [], 'strategy': {}}, source='x')
        path = os.path.join(self.tmp.name, 'a.mgt')
        self.store.save(h, path)
        h2 = self.store.load(path)
        self.assertNotEqual(h, h2)
        self.assertEqual(self.store.graphlet(h).dump(envelope=False),
                         self.store.graphlet(h2).dump(envelope=False))
        req = self.store.get(h2).request
        # the request is rebuilt from the saved envelope: the seed as the body states it
        self.assertEqual(self.store.graphlet(h).seed.sequence, req['seeds'][0]['sequence'])
        self.assertNotIn('labels', req['seeds'][0])     # annotate: no permitted set
        again = GraphletStore(self.tmp.name, clock=self.clock)
        self.assertIn(h, again)
        self.assertEqual('x', again.get(h).index['source'])
        self.assertEqual(T.doc_text('annotate'), again.graphlet(h).dump())
        self.assertEqual(2, len(again.list()))

    def test_put_all_skips_failed_seeds(self):
        resp = copy.deepcopy(response('fork'))
        resp['results'].append({'seed': {}, 'error': 'no carrier'})
        hs = self.store.put_all(resp, {})
        self.assertEqual(2, len(hs))
        self.assertIsNone(hs[1])
        self.assertTrue(re.match(r'g_', hs[0]))


if __name__ == '__main__':
    unittest.main()
