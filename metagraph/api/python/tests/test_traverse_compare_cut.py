"""graphlet_compare's page gives way in its base before its rows: the two evidence blocks,
then the export hint, are cut (named in fields_cut) as far as a page needs to carry its
first row whole. The search service found it live on staging: two time-budgeted walks of
one 200 bp seed on refseq33m (179,545 and 146,780 body bytes) answered result_too_large
at the default 2,048 bytes, the evidence blocks filling the page before any row. The pair
here is the real SRA one with the same shape: its uncut base alone (~1.75 KB) leaves no
room for a row at 2,048 or at 2,048 less the service's reserve."""

import gzip
import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib  # noqa: E402,F401

import metagraph.traverse.mcp_tools as M  # noqa: E402
from metagraph.traverse import Graphlet, GraphletStore  # noqa: E402
from metagraph.traverse.budget import LocalLimits  # noqa: E402
from metagraph.traverse.mcp_tools import GraphletTools, ToolLimits  # noqa: E402

SRA = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'traverse', 'real', 'sra')
ESSENTIAL = ('a', 'b', 'mode', 'comparable', 'equal', 'reason', 'depth_used', 'counts',
             'notes', 'support', 'scopes')
SERVICE_RESERVE = 256


def _load(name):
    with gzip.open(os.path.join(SRA, name + '.graphlet.json.gz')) as f:
        resp = json.loads(f.read())
    return Graphlet.from_response(resp['results'][0], resp)


def _size(obj):
    return len(json.dumps(obj, separators=(',', ':')))


class _Uncut:
    """The behaviour before the cut: no optional field of the base gives way."""

    def __enter__(self):
        self.saved = M._COMPARE_OPTIONAL
        M._COMPARE_OPTIONAL = ()

    def __exit__(self, *exc):
        M._COMPARE_OPTIONAL = self.saved


class TestCompareBaseCut(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.store = GraphletStore(os.path.join(self.tmp.name, 'sp'))
        self.tools = GraphletTools(self.store, {})
        self.a = self.store.put_graphlet(_load('sra_rand50_00__annotate_merge'))
        self.b = self.store.put_graphlet(_load('sra_rand50_00__cut_lists'))

    def tearDown(self):
        self.tmp.cleanup()

    def test_a_default_page_cuts_the_evidence_and_carries_rows(self):
        for mb in (None, 2048, 2048 - SERVICE_RESERVE):
            with _Uncut():
                self.assertEqual('result_too_large',
                                 self.tools.graphlet_compare(self.a, self.b, max_bytes=mb)['error'])
            out = self.tools.graphlet_compare(self.a, self.b, max_bytes=mb)
            self.assertNotIn('error', out, out)
            self.assertLessEqual(_size(out), mb or M.DEFAULT_MAX_BYTES)
            self.assertEqual(['evidence'], out['fields_cut'])
            self.assertNotIn('evidence', out)
            for k in ESSENTIAL:
                self.assertIn(k, out)
            self.assertEqual({'only_in_a': 224, 'only_in_b': 19, 'differ': 0}, out['counts'])
            self.assertTrue(out['rows'])
            self.assertNotIn('row_truncated', out)
            self.assertIn('next_cursor', out)

    def test_paging_yields_every_row_once_in_order(self):
        whole = self.tools.graphlet_compare(self.a, self.b, max_bytes=4 << 20)
        self.assertNotIn('fields_cut', whole)
        self.assertNotIn('next_cursor', whole)
        rows, cursor, pages = [], None, 0
        while True:
            out = self.tools.graphlet_compare(self.a, self.b, cursor=cursor)
            self.assertNotIn('error', out, out)
            self.assertLessEqual(_size(out), M.DEFAULT_MAX_BYTES)
            # each page states its own cut: never a field left out unnamed
            self.assertEqual('evidence' in out, 'fields_cut' not in out)
            rows.extend(out['rows'])
            pages += 1
            cursor = out.get('next_cursor')
            if not cursor:
                break
        self.assertGreater(pages, 1)
        self.assertEqual(whole['rows'], rows)

    def test_a_base_that_cannot_fit_is_refused_as_before(self):
        # below the essential base nothing is cut: cutting cannot make it fit
        with _Uncut():
            before = self.tools.graphlet_compare(self.a, self.b, max_bytes=300)
        after = self.tools.graphlet_compare(self.a, self.b, max_bytes=300)
        self.assertEqual('result_too_large', after['error'])
        self.assertEqual(before, after)

    def test_a_base_that_fits_is_answered_byte_for_byte_as_before(self):
        for mb in (4 << 20, 64 << 10, 8 << 10):
            with _Uncut():
                before = self.tools.graphlet_compare(self.a, self.b, max_bytes=mb)
            after = self.tools.graphlet_compare(self.a, self.b, max_bytes=mb)
            self.assertNotIn('fields_cut', after)
            self.assertEqual(json.dumps(before), json.dumps(after))
        # a compare of a handle with itself: a small base, unchanged at the default
        with _Uncut():
            before = self.tools.graphlet_compare(self.a, self.a)
        self.assertEqual(json.dumps(before), json.dumps(self.tools.graphlet_compare(self.a, self.a)))


class TestCompareStopCut(unittest.TestCase):
    """A comparison its local budget stops has no page: the same cut is made beside the
    local block that states the stop."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        d = self.tmp.name
        self.store = GraphletStore(os.path.join(d, 'sp'))
        self.tools = GraphletTools(self.store, {}, secret=b'k' * 32,
                                   export_dir=os.path.join(d, 'x'),
                                   local_limits=ToolLimits(view=LocalLimits(1),
                                                           heavy=LocalLimits(1)))
        self.a = self.store.put_graphlet(_load('sra_rand50_00__annotate_merge'))
        self.b = self.store.put_graphlet(_load('sra_rand50_00__cut_lists'))

    def tearDown(self):
        self.tmp.cleanup()

    def test_the_stop_answer_cuts_as_the_page_does(self):
        full = self.tools.graphlet_compare(self.a, self.b)
        self.assertEqual('unknown', full['comparable'])
        self.assertIn('local', full)
        self.assertNotIn('fields_cut', full)
        # where the answer fits (its local block compacted), it is answered unchanged
        for mb in (None, 1500, 1264):
            with _Uncut():
                before = self.tools.graphlet_compare(self.a, self.b, max_bytes=mb)
            self.assertNotIn('error', before)
            self.assertEqual(json.dumps(before),
                             json.dumps(self.tools.graphlet_compare(self.a, self.b, max_bytes=mb)))
        # where it was refused, the evidence (then the export hint) gives way
        for mb, cuts in ((1263, ['evidence']), (1100, ['evidence']), (500, ['evidence', 'export'])):
            with _Uncut():
                self.assertEqual('result_too_large',
                                 self.tools.graphlet_compare(self.a, self.b, max_bytes=mb)['error'])
            cut = self.tools.graphlet_compare(self.a, self.b, max_bytes=mb)
            self.assertEqual(cuts, cut['fields_cut'])
            self.assertLessEqual(_size(cut), mb)
            self.assertEqual(('unknown', None, []),
                             (cut['comparable'], cut['counts'], cut['rows']))
            self.assertEqual('work', cut['local']['stop']['resource'])


if __name__ == '__main__':
    unittest.main()
