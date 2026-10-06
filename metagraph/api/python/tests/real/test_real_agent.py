"""The agent tool layer (metagraph.traverse.mcp_tools, store, client) on REAL retrievals,
driven the way an agent drives it (DESIGN-traverse-graphlet.md §6, §14).

  * Scripted SESSIONS per index against the live servers (TestAgentSession*):
    traverse_fetch probe -> graphlet_summary -> graphlet_walks paged to exhaustion ->
    graphlet_walk -> graphlet_support -> graphlet_labels (name / ref selectors, an
    ambiguous bare name) -> graphlet_splits -> graphlet_claims paged -> graphlet_sequence
    slices -> graphlet_export fasta/gfa/json/mgt (a path outside the export dir refused)
    -> traverse_continue (execute and not) -> graphlet_compare -> graphlet_subtrie and
    the tools on the view -> graphlet_list / free / save / load, and the structured
    errors that need a backend (HTTP 400, an unreachable port, replay).
  * A SWEEP of the list tools over the cached cells (TestToolsOverCachedCells): every
    page within its ceiling, rows == total, no duplicate rows, the evidence block equal
    to the one derived from the server's own per-seed JSON summary, and the rows equal
    to the server's detail-full JSON where it has them (paths, splits, label_summary).
  * CURSORS and offline ERRORS on cached retrievals (no server needed).

Independent references: the server's per-seed JSON summary (outcome, limitations,
complete_to_bp, completeness_scope, labels_per_node) and detail-full JSON (paths,
segments, splits, label_summary) from the retrieval cache; the harness's spellers
(realdata.py, which read the G records only); the /resolve oracle. Every tool result
is checked on the way back by the Agent wrapper (size <= ceiling, labels as
{name, ref}, structured errors, the evidence block on local answers).

Skips cleanly without the servers or the cache:

    cd api/python/tests/real && python3 -m unittest -v test_real_agent
"""

import base64
import hashlib
import hmac
import json
import os
import re
import shutil
import socket
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import realdata as R  # noqa: E402  (puts api/python on sys.path)

from metagraph.traverse import (  # noqa: E402
    GraphletStore, GraphletView, TraverseClient, derive, dump, parse,
)
from metagraph.traverse.mcp_tools import (  # noqa: E402
    CAPABILITIES_MAX_BYTES, DEFAULT_MAX_BYTES, MIN_MAX_BYTES, SEQUENCE_MAX_BYTES, GraphletTools,
)

SPOOL_ROOT = os.path.join(R.SCRATCH, 'agent_spool')
EVIDENCE_KEYS = {'complete_to_bp', 'exact', 'support', 'reconverge', 'scope', 'outcome',
                 'limitations'}
LOCAL_ANSWERS = {'graphlet_summary', 'graphlet_walks', 'graphlet_walk', 'graphlet_support',
                 'graphlet_labels', 'graphlet_splits', 'graphlet_claims',
                 'graphlet_sequence', 'graphlet_subtrie'}
_REF = re.compile(r'(?:c:\d+|h:\d+:\d+)\Z')
CLAIM_KINDS = {'end', 'alive', 'diverged', 'merged', 'stretch', 'route_only', 'boundary'}
END_CLASSES = {'lost', 'dead_end', 'radius', 'capped', 'blocked', 'pruned', 'open'}

# the scripted sessions: per index the retrievals an agent makes (cells of the cache, so
# that the server's JSON is the independent reference; the live body is checked equal)
SESSIONS = {
    'mini_refseq': {'main': 'mini_ndm1__limit2', 'keep': 'mini_ndm1__limit2_keep',
                    'lower': 'mini_ndm1__cut_lists', 'extra': 'mini_ndm1__trace'},
    'uhgg': {'main': 'uhgg_conserved__limit2', 'keep': 'uhgg_conserved__limit2_keep',
             'lower': 'uhgg_16s__cut_lists', 'extra': 'uhgg_16s__cut_events'},
    'sra': {'main': 'sra_rand50_00__limit2', 'keep': 'sra_rand50_00__limit2_keep',
            'lower': 'sra_rand50_02__cut_lists', 'extra': 'sra_16s_PZ326290__cut_events'},
}
# a cell of another seed per index (graphlet_compare: comparable False)
OTHER_SEED = {'mini_refseq': 'mini_ndm1_rc__limit2', 'uhgg': 'uhgg_16s__limit2',
              'sra': 'sra_16s_PZ326290__limit2'}
# the sweep: bodies up to this size (the 64 MB-full "heavy" cells are excluded anyway);
# the detail-full JSON is compared where it is at most SWEEP_MAX_FULL
SWEEP_MAX_BODY = 1 << 20
SWEEP_MAX_FULL = 24 << 20
SWEEP_SMALL_PAGES = 4            # default-size pages checked per list in the sweep
# the sweep's reference pages: large enough for few pages (each page recomputes the
# whole list), small enough that a page's own cost stays low (_page re-measures the page
# after every row it adds, so one page of k rows costs O(k x page bytes))
SWEEP_BIG_PAGE = 1 << 16
SWEEP_WALK_PAGE = 1 << 18

STATS = {'tool_calls': 0, 'cells_swept': set(), 'session_cells': set(), 'pages': 0}


# ------------------------------------------------------------------ helpers

def new_spool(tag):
    os.makedirs(SPOOL_ROOT, exist_ok=True)
    return tempfile.mkdtemp(prefix=tag + '-', dir=SPOOL_ROOT)


def jsize(obj):
    """The bytes of a result as an MCP server sends it: compact JSON, UTF-8."""
    return len(json.dumps(obj, separators=(',', ':'), ensure_ascii=False).encode('utf-8'))


def iter_dicts(obj):
    """Every dict in a result, except the selectors a view echoes back ({ref} is an
    argument there, not a label description)."""
    if isinstance(obj, dict):
        yield obj
        for k, v in obj.items():
            if k != 'selectors':
                yield from iter_dicts(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from iter_dicts(v)


def check_labels(tc, out, where):
    """Labels leave the server as {name, ref}: never an id, never a bare name alone."""
    for d in iter_dicts(out):
        ref = d.get('ref')
        if isinstance(ref, str) and _REF.match(ref):
            tc.assertIn('name', d, (where, d))
            tc.assertNotIn('id', d, (where, d))


def closed_port():
    """A TCP port nothing listens on (bound, released, never listened)."""
    s = socket.socket()
    s.bind(('127.0.0.1', 0))
    port = s.getsockname()[1]
    s.close()
    return port


def informational(lim):
    """The scope statement of an arm on which no history was united (observed 0)."""
    return lim.get('kind') == 'scope' and lim.get('observed') == 0


# the label-class limitations (spec §7.0): lower bounds and overstatements
LABEL_CLASS = ('label_lists', 'inexact_counts', 'seed_labels', 'switch_sources',
               'greedy_losses', 'trace_record_boundaries')


def expected_evidence(result, envelope, side=None, delivery=None):
    """The evidence block derived from the SERVER's per-seed JSON summary alone (§6 as
    implemented: the outcome dimensions, the kinds of the limitations that apply per arm
    and per seed, exact iff the label evidence of every covered arm is complete and none
    of its lists was cut)."""
    arms = result['arms']
    sides = [s for s in ('left', 'right') if s in arms] if side is None else [side]
    lims = {}
    for s in sides:
        kinds = []
        for lim in arms[s].get('limitations') or []:
            if not informational(lim) and lim['kind'] not in kinds:
                kinds.append(lim['kind'])
        if kinds:
            lims[s] = kinds
    seed = []
    for lim in result.get('limitations') or []:
        if lim['kind'] not in seed:
            seed.append(lim['kind'])
    if result.get('resource_stop'):
        seed.append('resource_stop')
    if seed:
        lims['seed'] = seed
    outcome = dict(result['outcome'])
    if delivery is not None:
        outcome['delivery'] = delivery
    strategy = envelope['strategy']

    def arm_exact(s):
        # PER ARM (GPT review, finding 12): no label-class limitation of this arm or of
        # the seed, no cut list on it; one whose outcome is not complete with no such
        # limitation anywhere is exact nowhere
        if ((arms[s].get('labels_per_node') or {}).get('nodes_truncated') or 0) > 0:
            return False
        if any(l['kind'] in LABEL_CLASS for l in result.get('limitations') or []):
            return False
        stated = False
        for a, arm in arms.items():
            if any(l['kind'] in LABEL_CLASS for l in arm.get('limitations') or []):
                stated = True
                if a == s:
                    return False
        return stated or result['outcome']['label_evidence'] == 'complete'
    return {
        'complete_to_bp': {s: arms[s]['complete_to_bp'] for s in sides},
        'exact': all(arm_exact(s) for s in sides),
        'support': strategy.get('support', 'kmer'),
        'reconverge': strategy['branching']['on_reconverge'],
        'scope': {s: arms[s]['completeness_scope'] for s in sides},
        'outcome': outcome,
        'limitations': lims,
    }


def assert_evidence(tc, got, want, where=''):
    tc.assertIsInstance(got, dict, where)
    for k in ('complete_to_bp', 'exact', 'support', 'reconverge', 'scope', 'outcome'):
        tc.assertEqual(want[k], got[k], '%s: evidence.%s' % (where, k))
    tc.assertEqual(set(want['limitations']), set(got['limitations']),
                   '%s: evidence.limitations scopes' % (where,))
    for k, kinds in got['limitations'].items():
        tc.assertEqual(len(kinds), len(set(kinds)), '%s: duplicate kinds %r' % (where, kinds))
        tc.assertEqual(sorted(want['limitations'][k]), sorted(kinds),
                       '%s: evidence.limitations[%s]' % (where, k))


class Agent:
    """GraphletTools as an agent calls it. Every result is checked on the way back: its
    JSON size against the tool's ceiling (max(max_bytes, MIN_MAX_BYTES) when max_bytes is
    passed, else 2048; 16 KB for graphlet_sequence and for the request of
    traverse_continue(execute=False)), labels as {name, ref} without ids, an error as
    {error: str, message: str}, and the evidence block on every local answer."""

    def __init__(self, tc, tools):
        self.tc = tc
        self.tools = tools

    def __getattr__(self, name):
        fn = getattr(self.tools, name)
        if not name.startswith(('traverse_', 'graphlet_')):
            return fn

        def call(*args, **kw):
            out = fn(*args, **kw)
            STATS['tool_calls'] += 1
            tc = self.tc
            limit = kw.get('max_bytes')
            if isinstance(limit, int) and not isinstance(limit, bool):
                limit = max(limit, MIN_MAX_BYTES)
            else:
                dry_run = name == 'traverse_continue' and (
                    kw.get('execute', args[4] if len(args) > 4 else True) is False)
                limit = (self.tools.sequence_max_bytes
                         if name == 'graphlet_sequence' or dry_run else self.tools.max_bytes)
                if name == 'traverse_capabilities':
                    # the server's description of itself has its own default ceiling
                    limit = max(limit, CAPABILITIES_MAX_BYTES)
            tc.assertIsInstance(out, dict, name)
            n = jsize(out)
            tc.assertLessEqual(n, limit, '%s(%r, %r) returned %d bytes > %d'
                               % (name, args, kw, n, limit))
            check_labels(tc, out, name)
            if 'error' in out:
                tc.assertIsInstance(out['error'], str, out)
                if out['error'] != 'seed_failed':
                    tc.assertIsInstance(out.get('message'), str, out)
            elif name in LOCAL_ANSWERS and kw.get('detail', 'summary') == 'summary':
                ev = out.get('evidence')
                tc.assertIsInstance(ev, dict, '%s answered without evidence' % name)
                tc.assertTrue(EVIDENCE_KEYS <= set(ev), (name, sorted(ev)))
            return out
        return call


class Rows(list):
    """The rows of a paged answer; .cut maps the index of a row that came alone on a
    page with row_truncated to the cut fields that page named."""
    cut = None


def exhaust(tc, call, key=None, max_pages=50000):
    """Page |call| (a function of cursor=) to exhaustion: every page a page (no error),
    total constant, progress on every page, rows == total, no duplicate keys. A page
    with row_truncated holds exactly one row and names its cut fields."""
    out = call(cursor=None)
    tc.assertNotIn('error', out, out)
    total = out['total']
    rows, pages = Rows(), 0
    rows.cut = {}
    while True:
        pages += 1
        STATS['pages'] += 1
        tc.assertNotIn('error', out, out)
        tc.assertEqual(total, out['total'])
        if out.get('row_truncated'):
            tc.assertEqual(1, len(out['rows']), 'a truncated row comes alone')
            tc.assertTrue(out['cut_fields'], 'a truncated row names its cut fields')
            for f in out['cut_fields']:
                tc.assertIn(f, out['rows'][0])
            rows.cut[len(rows)] = list(out['cut_fields'])
        rows.extend(out['rows'])
        if 'next_cursor' not in out:
            break
        tc.assertTrue(out['rows'], 'a page with a next_cursor and no rows')
        tc.assertLess(pages, max_pages)
        out = call(cursor=out['next_cursor'])
    tc.assertEqual(total, len(rows), 'rows != total')
    if key is not None:
        keys = [key(r) for r in rows]
        tc.assertEqual(len(keys), len(set(keys)), 'duplicate rows across pages')
    return rows


def paged(tc, call, key=None, n_ok=False):
    """|call| (a function of cursor= and **kw) paged at its default size to exhaustion
    AND in one large page: the default pages hold exactly the reference rows, in order,
    each equal or cut where its page says so (so no row is repeated, skipped or moved
    across pages); |key| must be unique over the reference rows. -> (reference rows,
    paged rows)."""
    big = {'max_bytes': 1 << 22}
    if n_ok:
        big['n'] = 10 ** 7
    ref = exhaust(tc, lambda cursor: call(cursor=cursor, **big), key=key)
    tc.assertEqual({}, ref.cut)
    rows = exhaust(tc, lambda cursor: call(cursor=cursor))
    tc.assertEqual(len(ref), len(rows))
    for j, (r, f) in enumerate(zip(rows, ref)):
        assert_row_or_cut(tc, r, f, rows.cut.get(j), j)
    return ref, rows


def assert_row_or_cut(tc, row, full, cut, where=''):
    """|row| equals the uncut |full| row, except the named |cut| fields, which are
    prefixes of the full values (a string loses its tail, a list its last items)."""
    tc.assertEqual(set(full), set(row), where)
    for k, v in full.items():
        if cut and k in cut:
            tc.assertEqual(type(v), type(row[k]), (where, k))
            tc.assertEqual(v[:len(row[k])], row[k], (where, k))
            tc.assertLess(len(row[k]), len(v), (where, k))
        else:
            tc.assertEqual(v, row[k], (where, k))


def flank(molecule, seed, side):
    return molecule[len(seed):] if side == 'right' else molecule[:len(molecule) - len(seed)]


def full_paths(cell, i, side):
    return cell.full_result(i)['arms'][side]['paths']


def gfa_spell(s_seqs, names, k1):
    out = s_seqs[names[0]]
    for n in names[1:]:
        out += s_seqs[n][k1:]
    return out


# ------------------------------------------------------------------ live sessions

class _Session:
    """One scripted agent session on INDEX. setUpClass makes the session's retrievals
    live through traverse_fetch (the cached cell of the same seed and strategy is the
    independent reference; the live body is checked equal to it)."""

    INDEX = None

    @classmethod
    def setUpClass(cls):
        R.require_server(cls.INDEX)
        for cell in list(SESSIONS[cls.INDEX].values()) + [OTHER_SEED[cls.INDEX]]:
            R.require_cache(cls.INDEX, cell)
        if not R.live_matches_cache(cls.INDEX):
            raise unittest.SkipTest('%s: the live server is not the index the cache was '
                                    'filled from' % cls.INDEX)
        cls.spool = new_spool('session-' + cls.INDEX)
        cls.store = GraphletStore(cls.spool)
        cls.tools = GraphletTools(cls.store, {cls.INDEX: R.client(cls.INDEX)})
        cls.export_dir = os.path.realpath(os.path.join(cls.spool, 'exports'))
        cls.cells = {role: R.load_cell(cls.INDEX, cell)
                     for role, cell in SESSIONS[cls.INDEX].items()}
        cls.fetched, cls.handles = {}, {}
        for role, c in cls.cells.items():
            out = cls.tools.traverse_fetch(cls.INDEX, seed={'sequence': c.seed_sequence(0)},
                                           strategy=cls.strategy_of(c))
            if 'error' in out:
                raise AssertionError('%s/%s: traverse_fetch failed: %r'
                                     % (cls.INDEX, c.cell, out))
            cls.fetched[role] = out
            cls.handles[role] = out['handle']
            STATS['session_cells'].add('%s/%s' % (cls.INDEX, c.cell))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(getattr(cls, 'spool', ''), ignore_errors=True)

    @classmethod
    def strategy_of(cls, cell):
        return R.strategy(cell.strategy_name, cls.INDEX, cell.seed_record)

    def setUp(self):
        self.agent = Agent(self, self.tools)

    # -------------------------------------------------------------- shortcuts

    def cell(self, role):
        return self.cells[role]

    def h(self, role):
        return self.handles[role]

    def want(self, role, side=None, delivery=None):
        c = self.cell(role)
        return expected_evidence(c.graphlet_result(0), c.graphlet_response, side, delivery)

    def sides(self, role):
        return [s for s in ('left', 'right') if s in self.cell(role).graphlet_result(0)['arms']]

    def fresh(self, role):
        """The cached graphlet, parsed fresh (the harness's spellers read its G records)."""
        return self.cell(role).graphlet(0)

    def verifiable(self):
        return bool(R.identity(self.INDEX).get('index_fp'))

    def walk_with_continuation(self, role):
        g = self.store.graphlet(self.h(role))
        for side in self.sides(role):
            a = g.arms[side]
            for p in derive.paths(a):
                if a.segments[p.leaf].leaf.continuation is not None:
                    return side, p.id
        return None, None

    # -------------------------------------------------------------- the script

    def test_01_probe_then_fetch(self):
        """The probe (tight bounds, keep=False) stores nothing; each session retrieval
        is the cached body byte for byte, within 2 KB, with its server-derived evidence."""
        c = self.cell('main')
        before = (sorted(e['handle'] for e in self.store.list()),
                  sorted(os.listdir(os.path.join(self.spool, 'bodies'))))
        probe = self.agent.traverse_fetch(
            self.INDEX, seed={'sequence': c.seed_sequence(0)},
            strategy={'bounds': {'max_extension_bp': 60, 'max_live_paths': 20}}, keep=False)
        self.assertNotIn('error', probe, probe)
        self.assertNotIn('handle', probe)
        self.assertEqual(before, (sorted(e['handle'] for e in self.store.list()),
                                  sorted(os.listdir(os.path.join(self.spool, 'bodies')))))
        for side, to in probe['evidence']['complete_to_bp'].items():
            self.assertLessEqual(to, 60, side)
        self.assertEqual(probe['evidence']['outcome'], probe['summary']['outcome'])
        for role, out in self.fetched.items():
            cell = self.cell(role)
            self.assertLessEqual(jsize(out), DEFAULT_MAX_BYTES, role)
            self.assertEqual((self.INDEX, 'inline'), (out['index'], out['delivery']))
            self.assertEqual(cell.graphlet_result(0)['graphlet_bytes'], out['graphlet_bytes'])
            e = self.store.get(out['handle'])
            with open(os.path.join(self.spool, 'bodies', e.digest + '.mgt'),
                      encoding='utf-8', newline='') as f:
                self.assertEqual(cell.text(0), f.read(), '%s: live body != cached body' % role)
            self.assertEqual(hashlib.sha256(cell.text(0).encode()).hexdigest(), e.digest)
            assert_evidence(self, out['evidence'], self.want(role), role)
            s = out['summary']
            self.assertEqual(cell.graphlet_result(0)['outcome'], s['outcome'], role)
            self.assertEqual(self.verifiable(), s['index']['verifiable'], role)
            for side, arm in cell.graphlet_result(0)['arms'].items():
                if side in s['arms']:
                    self.assertEqual(arm['complete_to_bp'], s['arms'][side]['complete_to_bp'])

    def test_02_summary_and_limitations(self):
        for role in self.cells:
            h, result = self.h(role), self.cell(role).graphlet_result(0)
            out = self.agent.graphlet_summary(h)
            self.assertNotIn('error', out, out)
            assert_evidence(self, out['evidence'], self.want(role), role)
            for side in self.sides(role):
                one = self.agent.graphlet_summary(h, arm=side)
                assert_evidence(self, one['evidence'], self.want(role, side), (role, side))
                self.assertEqual([side], list(one['summary']['arms']))
                # every stated limitation in full, paged: the server's own list
                rows = exhaust(self, lambda cursor: self.agent.graphlet_summary(
                    h, arm=side, detail='limitations', cursor=cursor))
                want = [(l['kind'], l['knob'], l['limit'], l['observed'])
                        for l in (result['arms'][side].get('limitations') or [])
                        + (result.get('limitations') or [])]
                got = [(r['kind'], r['knob'], r['limit'], r['observed']) for r in rows]
                self.assertEqual(sorted(map(repr, want)), sorted(map(repr, got)), (role, side))
                for r in rows:
                    self.assertIn('effect', r)
                    self.assertEqual(r.get('informational', False),
                                     r['kind'] == 'scope' and r['observed'] == 0)

    def test_03_walks_paged_to_exhaustion(self):
        for role in self.cells:
            h, cell = self.h(role), self.cell(role)
            for side in self.sides(role):
                paths = full_paths(cell, 0, side)
                ev_want = self.want(role, side)
                ref = self.agent.graphlet_walks(h, side, n=10 ** 6, max_bytes=1 << 24)
                ref_rows = {r['walk']: r for r in ref['rows']}
                self.assertEqual(len(paths), ref['total'], (role, side))
                for rank in ('support', 'length', 'loss', 'id'):
                    def call(cursor, max_bytes=None, n=10, rank=rank):
                        out = self.agent.graphlet_walks(h, side, rank=rank, n=n, cursor=cursor,
                                                        max_bytes=max_bytes)
                        if 'error' not in out:
                            assert_evidence(self, out['evidence'], ev_want, (role, side, rank))
                        return out
                    rows = exhaust(self, call, key=lambda r: r['walk'])
                    ids = [r['walk'] for r in rows]
                    self.assertEqual(set(range(len(paths))), set(ids), (role, side, rank))
                    # stable: the same order again, with other page sizes and n=1
                    self.assertEqual(ids, [r['walk'] for r in exhaust(self, call)])
                    self.assertEqual(ids, [r['walk'] for r in exhaust(
                        self, lambda cursor: call(cursor, max_bytes=700))])
                    self.assertEqual(ids, [r['walk'] for r in exhaust(
                        self, lambda cursor: call(cursor, n=1))])
                    if rank == 'support':
                        keys = [(-r['n_full'], -r['n_alive'], -r['length_bp'], r['walk'])
                                for r in rows]
                        self.assertEqual(sorted(keys), keys, (role, side))
                    elif rank == 'length':
                        keys = [(-r['length_bp'], r['walk']) for r in rows]
                        self.assertEqual(sorted(keys), keys, (role, side))
                    elif rank == 'id':
                        self.assertEqual(sorted(ids), ids)
                    for j, r in enumerate(rows):
                        # rows equal the one-page reference unless cut, and a cut is named
                        assert_row_or_cut(self, r, ref_rows[r['walk']], rows.cut.get(j),
                                          (role, side, rank, r['walk']))
                        self.assertEqual(paths[r['walk']]['length_bp'], r['length_bp'])
                        self.assertLessEqual(len(r['labels_full']), 8)
                        self.assertEqual(r['n_full'] > 8, r.get('labels_full_truncated', False))
                        self.assertEqual(max(0, r['length_bp'] - ev_want['complete_to_bp'][side]),
                                         r['beyond_certified_bp'])
                        if r['complete']:
                            self.assertEqual(0, r['beyond_certified_bp'])

    def test_04_walks_spelled_and_truncated_rows(self):
        """spell=tail / full against the harness's speller; a row that does not fit its
        page comes alone, cut, and the cut fields are named (never silently shortened)."""
        for role in ('main', 'keep'):
            h, g = self.h(role), self.fresh(role)
            seed = g.seed.sequence
            for side in self.sides(role):
                mols = {pid: R.spell_path_with_seed(g, side, pid)
                        for pid in range(len(R.leaf_segments(g, side)))}
                tail = exhaust(self, lambda cursor: self.agent.graphlet_walks(
                    h, side, spell='tail', tail_bp=40, cursor=cursor), key=lambda r: r['walk'])
                for j, r in enumerate(tail):
                    fl = flank(mols[r['walk']], seed, side)
                    want = fl[-40:] if side == 'right' else fl[:40]
                    if 'sequence' in (tail.cut.get(j) or ()):
                        self.assertTrue(want.startswith(r['sequence']))
                    else:
                        self.assertEqual(want, r['sequence'], (role, side, r['walk']))
                rows = exhaust(self, lambda cursor: self.agent.graphlet_walks(
                    h, side, spell='full', cursor=cursor), key=lambda r: r['walk'])
                ref = {r['walk']: r for r in self.agent.graphlet_walks(
                    h, side, spell='full', n=10 ** 6, max_bytes=1 << 24)['rows']}
                for j, r in enumerate(rows):
                    fl = flank(mols[r['walk']], seed, side)
                    self.assertEqual(fl, ref[r['walk']]['sequence'], (role, side, r['walk']))
                    assert_row_or_cut(self, r, ref[r['walk']], rows.cut.get(j),
                                      (role, side, r['walk']))
                    if 'sequence' not in (rows.cut.get(j) or ()):
                        self.assertEqual(fl, r['sequence'], (role, side, r['walk']))
                STATS.setdefault('truncated_rows', 0)
                STATS['truncated_rows'] += len(rows.cut)

    def test_05_walk_chains(self):
        for role in self.cells:
            h, cell = self.h(role), self.cell(role)
            g = self.store.graphlet(h)
            for side in self.sides(role):
                paths = full_paths(cell, 0, side)
                segs = cell.full_result(0)['arms'][side]['segments']
                n = len(paths)
                for pid in sorted({0, n // 2, n - 1}):
                    rows, _ = paged(self, lambda cursor, **kw: self.agent.graphlet_walk(
                        h, side, pid, cursor=cursor, **kw), key=lambda r: r['segment'])
                    self.assertEqual(paths[pid]['segments'], [r['segment'] for r in rows])
                    pos = 0
                    for r in rows:
                        s = segs[r['segment']]
                        self.assertEqual((s['from_bp'], s['from_bp'] + s['length_bp']),
                                         (r['from_bp'], r['to_bp']))
                        self.assertEqual(pos, r['from_bp'])
                        pos = r['to_bp']
                        self.assertLessEqual(r['labels_in'], r['labels_in_total'])
                        self.assertLessEqual(r['labels_out'], r['labels_out_total'])
                        self.assertEqual(r['exact'], r['labels_in'] == r['labels_in_total']
                                         and r['labels_out'] == r['labels_out_total'])
                        self.assertEqual(min(r['labels_in'], 8), len(r['in']))
                        if len(s['parents']) > 1:
                            self.assertEqual(s['parents'], r['merge_of'])
                    self.assertEqual(paths[pid]['length_bp'], pos)
                    page = self.agent.graphlet_walk(h, side, pid)
                    assert_evidence(self, page['evidence'], self.want(role, side), (role, side))
                    leaf = g.arms[side].segments[R.leaf_segments(g, side)[pid]]
                    self.assertEqual(leaf.leaf.continuation is not None, 'continuation' in page)
                bad = self.agent.graphlet_walk(h, side, n)
                self.assertEqual('bad_argument', bad['error'])

    def test_06_support(self):
        for role in self.cells:
            h, cell = self.h(role), self.cell(role)
            for side in self.sides(role):
                paths = full_paths(cell, 0, side)
                for pid in sorted({0, len(paths) - 1}):
                    rows, _ = paged(self, lambda cursor, **kw: self.agent.graphlet_support(
                        h, side, pid, cursor=cursor, **kw), key=lambda r: r['from_bp'])
                    if paths[pid]['length_bp'] == 0:
                        self.assertEqual([], rows)
                        continue
                    pos = 0
                    for r in rows:
                        self.assertEqual(pos, r['from_bp'], (role, side, pid))
                        self.assertLess(r['from_bp'], r['to_bp'])
                        self.assertLessEqual(r['n_labels'], r['n_labels_total'])
                        pos = r['to_bp']
                    self.assertEqual(paths[pid]['length_bp'], pos, (role, side, pid))
                    assert_evidence(self, self.agent.graphlet_support(h, side, pid)['evidence'],
                                    self.want(role, side), (role, side, pid))
                    stepped = exhaust(self, lambda cursor: self.agent.graphlet_support(
                        h, side, pid, step=50, cursor=cursor))
                    self.assertTrue(all(r['to_bp'] - r['from_bp'] >= 50 or 'why' in r
                                        for r in stepped))
                    self.assertLessEqual(len(stepped), len(rows))

    def test_07_labels(self):
        for role in self.cells:
            h, cell = self.h(role), self.cell(role)
            full = cell.full_result(0)
            names = full['label_dict']
            for arm in [None] + self.sides(role):
                rows, _ = paged(self, lambda cursor, **kw: self.agent.graphlet_labels(
                    h, arm=arm, cursor=cursor, **kw), key=lambda r: (r['ref'], r['arm']),
                    n_ok=True)
                sides = self.sides(role) if arm is None else [arm]
                self.assertEqual(len(names) * len(sides), len(rows), (role, arm))
                assert_evidence(self, self.agent.graphlet_labels(h, arm=arm)['evidence'],
                                self.want(role, arm), (role, arm))
                keys = [(-r['direct_bp'], r['ref'], r['arm']) for r in rows]
                self.assertEqual(sorted(keys), keys)
                by_ref = {}
                for i, l in enumerate(names):
                    ref = 'h:%d:%d' % (l['column'], l['seq_id']) if l['kind'] == 'header' \
                        else 'c:%d' % l['column']
                    by_ref[ref] = (l['name'], full['label_summary'][i])
                for r in rows:
                    name, summ = by_ref[r['ref']]
                    self.assertEqual(name, r['name'])
                    s = summ[r['arm']]
                    self.assertEqual((s['direct_bp'], s['reach_bp'], s['reentries'], len(s['runs'])),
                                     (r['direct_bp'], r['reach_bp'], r['reentries'], r['runs']),
                                     (role, r['ref'], r['arm']))
                reach = exhaust(self, lambda cursor: self.agent.graphlet_labels(
                    h, arm=arm, rank='reach_bp', cursor=cursor))
                keys = [(-r['reach_bp'], r['ref'], r['arm']) for r in reach]
                self.assertEqual(sorted(keys), keys)
            # one label by name, by ref, bare name, bare ref: one answer; the cursor is
            # bound to the NORMALIZED selector, so it carries over between spellings
            g = self.store.graphlet(h)
            lab = g.labels[0]
            got = []
            for sel in ({'name': lab.name}, {'ref': lab.ref}, lab.name, lab.ref):
                rows = exhaust(self, lambda cursor: self.agent.graphlet_labels(
                    h, name=sel, n=1, cursor=cursor))
                got.append(rows)
            self.assertTrue(all(x == got[0] for x in got))
            first = self.agent.graphlet_labels(h, name={'name': lab.name}, n=1)
            if 'next_cursor' in first:
                again = self.agent.graphlet_labels(h, name=lab.ref, n=1,
                                                   cursor=first['next_cursor'])
                self.assertNotIn('error', again)
            self.assertEqual('unknown_label', self.agent.graphlet_labels(
                h, name={'ref': 'c:987654321'})['error'])
            self.assertEqual('unknown_label', self.agent.graphlet_labels(
                h, name='no such label at all')['error'])
            # at_bp: the labels alive at that outward base
            for side in self.sides(role):
                at = cell.graphlet_result(0)['arms'][side]['counts']['max_bp'] // 2
                alive = exhaust(self, lambda cursor: self.agent.graphlet_labels(
                    h, arm=side, at_bp=at, cursor=cursor))
                every = {(r['ref'], r['arm']) for r in exhaust(
                    self, lambda cursor: self.agent.graphlet_labels(h, arm=side, cursor=cursor))}
                self.assertTrue({(r['ref'], r['arm']) for r in alive} <= every)

    def test_08_ambiguous_bare_name(self):
        """No real index of this suite has two labels of one name (checked over every
        cached cell), so the ambiguity is made from the REAL main retrieval: label 1
        renamed to label 0's name, label 3 named like label 2's ref. A bare string
        then needs a ref; a tagged ref always resolves."""
        g = self.store.graphlet(self.h('main'))
        if len(g.labels) < 4:
            self.skipTest('fewer than 4 labels')
        g2 = parse(dump(g, envelope=True))
        dup = g2.labels[0].name
        g2.labels[1].name = dup
        g2.labels[3].name = g2.labels[2].ref
        h2 = self.store.put_graphlet(g2, self.store.get(self.h('main')).request,
                                     source=self.INDEX)
        side = self.sides('main')[0]
        for out in (self.agent.graphlet_labels(h2, name=dup),
                    self.agent.graphlet_labels(h2, name={'name': dup}),
                    self.agent.graphlet_labels(h2, name=g2.labels[2].ref),
                    self.agent.graphlet_walks(h2, side, label=dup),
                    self.agent.graphlet_claims(h2, side, labels=[dup]),
                    self.agent.graphlet_subtrie(h2, dup)):
            self.assertEqual('ambiguous_label', out['error'], out)
            self.assertIn('ref', out['hint'])
        for lab in g2.labels[:4]:
            out = self.agent.graphlet_labels(h2, name={'ref': lab.ref})
            self.assertNotIn('error', out, out)
            self.assertEqual((lab.ref, lab.name), (out['ref'], out['name']))
        self.tools.graphlet_free(h2)

    def test_09_splits(self):
        for role in self.cells:
            h, cell = self.h(role), self.cell(role)
            for side in self.sides(role):
                full = cell.full_result(0)['arms'][side]['splits']
                rows, _ = paged(self, lambda cursor, **kw: self.agent.graphlet_splits(
                    h, side, min_labels_before=0, cursor=cursor, **kw),
                    key=lambda r: r['segment'], n_ok=True)
                self.assertEqual(cell.graphlet_result(0)['arms'][side]['counts']['splits'],
                                 len(rows))
                assert_evidence(self, self.agent.graphlet_splits(h, side)['evidence'],
                                self.want(role, side), (role, side))
                want = sorted((s['segment'], s['at_bp'], s['kind'], s['labels_before'],
                               [(b['char'], b['labels_distinct']) for b in s['branches']])
                              for s in full)
                got = sorted((r['segment'], r['at_bp'], r['kind'], r['labels_before'],
                              [(b['char'], b['labels']) for b in r['branches']]) for r in rows)
                self.assertEqual(want, got, (role, side))
                two = exhaust(self, lambda cursor: self.agent.graphlet_splits(
                    h, side, cursor=cursor))
                self.assertEqual(sum(1 for r in rows if r['labels_before'] >= 2), len(two))

    def test_10_claims_paged(self):
        for role in self.cells:
            h, cell = self.h(role), self.cell(role)
            for side in self.sides(role):
                arm_json = cell.graphlet_result(0)['arms'][side]
                every, _ = paged(self, lambda cursor, **kw: self.agent.graphlet_claims(
                    h, side, route_consistent=False, cursor=cursor, **kw))
                first = self.agent.graphlet_claims(h, side)
                consistent, _ = paged(self, lambda cursor, **kw: self.agent.graphlet_claims(
                    h, side, cursor=cursor, **kw))
                # merge-entered claims (displayed only from evidence_from) and route_only
                # ones (none at a cut) are counted apart (GPT review, finding 11)
                filtered = sum((first.get('filtered') or {}).values())
                self.assertEqual(len(every), len(consistent) + filtered, (role, side))
                self.assertEqual(bool(filtered), 'hint' in first)
                cut = ((arm_json.get('labels_per_node') or {}).get('nodes_truncated') or 0) > 0
                self.assertEqual(cut and cell.graphlet_result(0)['label_mode'] == 'annotate',
                                 'caveat' in first, (role, side))
                assert_evidence(self, first['evidence'], self.want(role, side), (role, side))
                for r in every:
                    self.assertIn(r['kind'], CLAIM_KINDS)
                    self.assertIn(r['end_class'], END_CLASSES)
                    self.assertLessEqual(r['from_bp'], r['to_bp'])
                    if r['evidence_from'] is not None:
                        self.assertLessEqual(r['from_bp'], r['evidence_from'])
                        self.assertLessEqual(r['evidence_from'], r['to_bp'])
                    self.assertNotIn('arm', r)
                if every:
                    ref = every[0]['label']['ref']
                    one = exhaust(self, lambda cursor: self.agent.graphlet_claims(
                        h, side, labels=[{'ref': ref}], route_consistent=False, cursor=cursor))
                    self.assertEqual(sum(1 for r in every if r['label']['ref'] == ref), len(one))
                    self.assertTrue(all(r['label']['ref'] == ref for r in one))
                    longer = exhaust(self, lambda cursor: self.agent.graphlet_claims(
                        h, side, min_bp=100, route_consistent=False, cursor=cursor))
                    self.assertTrue(all(r['to_bp'] - (r['evidence_from'] if r['evidence_from']
                                                      is not None else r['to_bp']) >= 100
                                        for r in longer))

    def test_11_sequence_slices(self):
        h, g = self.h('main'), self.fresh('main')
        seed = g.seed.sequence
        for side in self.sides('main'):
            n = len(R.leaf_segments(g, side))
            for pid in sorted({0, n // 2, n - 1}):
                mol = R.spell_path_with_seed(g, side, pid)
                fl = flank(mol, seed, side)
                for with_seed, want in ((True, mol), (False, fl)):
                    out = self.agent.graphlet_sequence(h, side, walk=pid, with_seed=with_seed)
                    self.assertNotIn('error', out, out)
                    self.assertEqual((len(want), SEQUENCE_MAX_BYTES), (out['length'],
                                                                        out['ceiling_bytes']))
                    got, page = out['sequence'], out
                    while page.get('truncated'):
                        page = self.agent.graphlet_sequence(h, side, walk=pid,
                                                            with_seed=with_seed,
                                                            start=page['next_from'])
                        got += page['sequence']
                    self.assertEqual(want, got, (side, pid, with_seed))
                    cuts = [0, len(want) // 3, (2 * len(want)) // 3, len(want)]
                    pieces = [self.agent.graphlet_sequence(h, side, walk=pid, with_seed=with_seed,
                                                           start=a, end=b)
                              for a, b in zip(cuts, cuts[1:])]
                    self.assertEqual(want, ''.join(p['sequence'] for p in pieces))
                    for p, (a, b) in zip(pieces, zip(cuts, cuts[1:])):
                        self.assertEqual((a, b), (p['from'], p['to']))
                        assert_evidence(self, p['evidence'], self.want('main', side))
                    clipped = self.agent.graphlet_sequence(h, side, walk=pid, start=0,
                                                           end=len(fl) + 1000)
                    self.assertEqual(len(fl), clipped['to'])
                    empty = self.agent.graphlet_sequence(h, side, walk=pid, start=len(fl))
                    self.assertEqual('', empty['sequence'])
                    self.assertEqual('bad_argument', self.agent.graphlet_sequence(
                        h, side, walk=pid, start=len(fl) + 1)['error'])
                # walking order: outward base i at [i]
                walk_flank = R.walk_order_flank(g, side, R.leaf_segments(g, side)[pid])
                out = self.agent.graphlet_sequence(h, side, walk=pid, orientation='walk')
                self.assertEqual(walk_flank, out['sequence'])
                out = self.agent.graphlet_sequence(h, side, walk=pid, orientation='walk',
                                                   with_seed=True)
                self.assertEqual((seed if side == 'right' else seed[::-1]) + walk_flank,
                                 out['sequence'])
            a = g.arms[side]
            for sid in sorted({0, len(a.segments) - 1}):
                out = self.agent.graphlet_sequence(h, side, segment=sid)
                w = a.segments[sid].walk
                self.assertEqual(w if side == 'right' else w[::-1], out['sequence'])
                out = self.agent.graphlet_sequence(h, side, segment=sid, orientation='walk')
                self.assertEqual(w, out['sequence'])
        # the 16 KB ceiling holds for the whole result: a smaller ceiling pages
        small = GraphletTools(self.store, {}, sequence_max_bytes=1500)
        agent = Agent(self, small)
        side = self.sides('main')[-1]
        pid = max(range(len(R.leaf_segments(g, side))),
                  key=lambda p: len(R.spell_path_with_seed(g, side, p)))
        mol = R.spell_path_with_seed(g, side, pid)
        page = agent.graphlet_sequence(h, side, walk=pid, with_seed=True)
        got, pages = page['sequence'], 1
        while page.get('truncated'):
            self.assertEqual(page['to'], page['next_from'])
            page = agent.graphlet_sequence(h, side, walk=pid, with_seed=True,
                                           start=page['next_from'])
            got += page['sequence']
            pages += 1
        self.assertEqual(mol, got)
        self.assertGreater(pages, 1)

    def test_12_exports(self):
        h, cell, g = self.h('main'), self.cell('main'), self.fresh('main')
        seed = g.seed.sequence
        sides = self.sides('main')
        n_walks = sum(cell.graphlet_result(0)['arms'][s]['counts']['leaves'] for s in sides)
        n_segs = sum(cell.graphlet_result(0)['arms'][s]['counts']['segments'] for s in sides)
        # FASTA: one record per walk, every molecule the speller's
        out = self.agent.graphlet_export(h, 'fasta', path='walks.fa')
        self.assertNotIn('error', out, out)
        self.assertEqual(os.path.join(self.export_dir, 'walks.fa'), out['path'])
        self.assertEqual((n_walks, os.path.getsize(out['path'])), (out['records'], out['bytes']))
        with open(out['path']) as f:
            recs = f.read().split('>')[1:]
        for rec in recs:
            head, seq = rec.split('\n', 1)
            side, pid = head.split(' ')[0].split('_')
            self.assertEqual(R.spell_path_with_seed(g, side, int(pid)), seq.replace('\n', ''))
        one = self.agent.graphlet_export(h, 'fasta', arm=sides[0], walks=[0], path='one.fa')
        self.assertEqual(1, one['records'])
        # GFA: S per segment (+ seed), P per walk; every L overlaps by k-1 and every P
        # spells the walk's molecule
        out = self.agent.graphlet_export(h, 'gfa', path='trie.gfa')
        with open(out['path']) as f:
            lines = [l.split('\t') for l in f.read().splitlines()]
        S = {l[1]: l[2] for l in lines if l[0] == 'S'}
        L = [l for l in lines if l[0] == 'L']
        P = [l for l in lines if l[0] == 'P']
        self.assertEqual((n_segs + 1, n_walks), (len(S), len(P)))
        self.assertEqual(len(S) + len(L) + len(P), out['records'])
        k1 = g.k - 1
        for l in L:
            self.assertEqual(S[l[1]][-k1:], S[l[3]][:k1], l)
        for p in P:
            side, pid = p[1].split('_')
            names = [x[:-1] for x in p[2].split(',')]
            self.assertEqual(R.spell_path_with_seed(g, side, int(pid)), gfa_spell(S, names, k1))
        # JSON: the server's detail-full result after the §2.5 normalisation
        out = self.agent.graphlet_export(h, 'json', path='full.json')
        with open(out['path']) as f:
            mine = R.normalize_full_result(json.load(f))
        theirs = R.normalize_full_result(cell.full_result(0))
        self.assertEqual([], R.json_diff(theirs, mine))
        out = self.agent.graphlet_export(h, 'json', what='walk', arm=sides[0], walks=[0, 1],
                                         path='w.json')
        with open(out['path']) as f:
            w = json.load(f)
        self.assertEqual([0, 1], [x['walk'] for x in w['walks']])
        self.assertEqual(flank(R.spell_path_with_seed(g, sides[0], 1), seed, sides[0]),
                         w['walks'][1]['sequence'])
        self.assertEqual(full_paths(cell, 0, sides[0])[1]['segments'],
                         [r['segment'] for r in w['walks'][1]['segments']])
        # MGT: H, J, the unchanged body
        out = self.agent.graphlet_export(h, 'mgt', path='main.mgt')
        with open(out['path'], encoding='utf-8', newline='') as f:
            text = f.read()
        lines = text.split('\n')
        body = cell.text(0).split('\n')
        self.assertTrue(lines[1].startswith('J {'))
        # the body unchanged around the J line; Z counts the document's lines, J included
        self.assertEqual(body[:1] + body[1:-2], lines[:1] + lines[2:-2])
        self.assertEqual('Z %d' % (int(body[-2].split(' ')[1]) + 1), lines[-2])
        self.assertEqual(text, dump(parse(text)))
        # a default name; refused paths outside the export directory, nothing written
        out = self.agent.graphlet_export(h, 'gfa')
        self.assertEqual(os.path.join(self.export_dir, h + '.gfa'), out['path'])
        outside = os.path.join(self.spool, 'outside')
        os.makedirs(outside, exist_ok=True)
        link = os.path.join(self.export_dir, 'escape')
        if not os.path.islink(link):
            os.symlink(outside, link)
        for path in ('../outside/x.fa', os.path.join(outside, 'x.fa'), 'escape/x.fa',
                     'sub/../../outside/x.fa', '/etc/passwd'):
            for fmt in ('fasta', 'gfa', 'json', 'mgt'):
                self.assertEqual('path_not_allowed', self.agent.graphlet_export(
                    h, fmt, path=path)['error'], (fmt, path))
            self.assertEqual('path_not_allowed', self.agent.graphlet_save(h, path)['error'])
        self.assertEqual([], os.listdir(outside))
        for bad in ('', 'a\0b', 7):
            self.assertEqual('bad_argument', self.agent.graphlet_export(
                h, 'fasta', path=bad)['error'], bad)

    def test_13_continue(self):
        role = 'main'
        side, pid = self.walk_with_continuation(role)
        if side is None:
            role = 'lower'
            side, pid = self.walk_with_continuation(role)
        self.assertIsNotNone(side, 'no walk with a continuation in the session')
        h, g = self.h(role), self.fresh(role)
        mol = R.spell_path_with_seed(g, side, pid)
        chain = self.agent.graphlet_walk(h, side, pid)
        self.assertIn('continuation', chain)
        # the request, not run (an agent may inspect it first)
        req = self.tools.traverse_continue(h, side, pid, execute=False,
                                           overrides={'bounds': {'max_extension_bp': 200}})
        if 'error' not in req:
            self.assertEqual({'handle': h, 'arm': side, 'walk': pid,
                              'overlap_bp': len(req['request']['seeds'][0]['sequence'])},
                             req['parent'])
        # the backend runs it: a new, certified-on-its-own retrieval linked to its parent
        new = self.agent.traverse_continue(h, side, pid,
                                           overrides={'bounds': {'max_extension_bp': 200}})
        self.assertNotIn('error', new, new)
        parent = new['parent']
        self.assertEqual((h, side, pid), (parent['handle'], parent['arm'], parent['walk']))
        g2 = self.store.graphlet(new['handle'])
        cont = g2.seed.sequence
        self.assertEqual(parent['overlap_bp'], len(cont))
        self.assertEqual(chain['continuation']['length_bp'], len(cont))
        # the continuation is the molecule's outer end, natural orientation
        self.assertTrue(mol.endswith(cont) if side == 'right' else mol.startswith(cont))
        self.assertEqual([side], list(g2.arms))
        self.assertLessEqual(g2.arms[side].complete_to_bp, 200)
        self.assertTrue(R.oracle_present(self.INDEX, cont))
        if g2.mode == 'constrain':
            self.assertEqual(chain['continuation']['labels'], g2.seed.num_seed_labels)
        s = self.agent.graphlet_summary(new['handle'])
        self.assertEqual(parent, s['parent'])
        listed = [e for e in exhaust(self, lambda cursor: self.agent.graphlet_list(
            cursor=cursor)) if e['handle'] == new['handle']]
        self.assertEqual((h, parent), (listed[0]['derived_from'], listed[0]['parent']))
        # a bad override is the backend's 400, structured
        bad = self.agent.traverse_continue(h, side, pid,
                                           overrides={'bounds': {'max_extension_bp': 'far'}})
        self.assertEqual(('backend_error', 400), (bad['error'], bad.get('status')))
        self.assertEqual('bad_argument', self.agent.traverse_continue(
            h, side, pid, overrides=[1])['error'])
        self.tools.graphlet_free(new['handle'])

    def test_14_continue_request_is_deliverable(self):
        """traverse_continue(execute=False) returns the request an agent would run, under
        the 16 KB sequence ceiling it states (the request carries the continuation's
        bases and label names); max_bytes raises it."""
        found = 0
        for role in self.cells:
            g = self.store.graphlet(self.h(role))
            for side in self.sides(role):
                a = g.arms[side]
                for p in derive.paths(a):
                    if a.segments[p.leaf].leaf.continuation is None:
                        continue
                    out = self.agent.traverse_continue(self.h(role), side, p.id, execute=False)
                    self.assertNotIn('error', out, '%s %s walk %d: %r' % (role, side, p.id, out))
                    self.assertEqual(g.next_request(side, [p.id]), out['request'])
                    self.assertEqual(SEQUENCE_MAX_BYTES, out['ceiling_bytes'])
                    found += 1
                    break
        self.assertGreater(found, 0)

    def test_15_compare(self):
        idx = self.INDEX
        main, keep = self.h('main'), self.h('keep')
        # the same seed and strategy fetched again (a replay): equal where verifiable
        again = self.agent.traverse_fetch(replay=main)
        self.assertNotIn('error', again, again)
        twin = again['handle']
        self.assertEqual(self.store.get(main).digest, self.store.get(twin).digest)
        for mode in ('claims', 'walks', 'labels', 'prefix_subset'):
            out = self.agent.graphlet_compare(main, twin, mode=mode)
            if self.verifiable():
                self.assertEqual((True, True, 0), (out['comparable'], out['equal'], out['total']))
            else:
                self.assertEqual(('unverifiable', None), (out['comparable'], out['equal']))
                self.assertIn('index_fp', out['reason'])
            self.assertEqual({'a', 'b'}, set(out['evidence']))
            assert_evidence(self, out['evidence']['a'], self.want('main'))
        # merge vs keep of the same seed: qualified, never equal
        for mode in ('claims', 'walks'):
            first = self.agent.graphlet_compare(main, keep, mode=mode)
            self.assertIsNone(first['equal'])
            self.assertEqual('qualified' if self.verifiable() else 'unverifiable',
                             first['comparable'])
            rows, _ = paged(self, lambda cursor, **kw: self.agent.graphlet_compare(
                main, keep, mode=mode, cursor=cursor, **kw),
                key=lambda r: json.dumps(r, sort_keys=True))
            self.assertEqual(sum(first['counts'].values()), len(rows))
            for tag, k in (('a', 'only_in_a'), ('b', 'only_in_b'), ('differ', 'differ')):
                self.assertEqual(first['counts'][k], sum(1 for r in rows if r['side'] == tag))
            if 'next_cursor' in first:
                self.assertEqual('bad_cursor', self.agent.graphlet_compare(
                    main, twin, mode=mode, cursor=first['next_cursor'])['error'])
            ex = self.agent.graphlet_export(main, 'json', what='compare', other=keep,
                                            path='cmp-%s.json' % mode)
            if mode == 'claims':
                self.assertEqual(sum(first['counts'].values()), ex['records'])
        # another seed: not comparable
        oc = R.load_cell(idx, OTHER_SEED[idx])
        other = self.store.put(oc.graphlet_response, oc.request, source=idx)
        out = self.agent.graphlet_compare(main, other)
        self.assertEqual((False, None), (out['comparable'], out['equal']))
        self.assertIn('seed', out['reason'])
        self.tools.graphlet_free(other)
        self.tools.graphlet_free(twin)

    def test_16_subtrie_and_the_tools_on_the_view(self):
        h = self.h('main')
        g = self.store.graphlet(h)
        side = self.sides('main')[0]
        n = len(derive.paths(g.arms[side]))
        chosen = None
        for lab in g.labels:
            ws = exhaust(self, lambda cursor: self.agent.graphlet_walks(
                h, side, label={'ref': lab.ref}, route_consistent=False, cursor=cursor))
            if 0 < len(ws) < n:
                chosen, carrying = lab, {r['walk'] for r in ws}
                break
        if chosen is None:
            self.skipTest('no label carried by a strict subset of the walks')
        bodies = sorted(os.listdir(os.path.join(self.spool, 'bodies')))
        sub = self.agent.graphlet_subtrie(h, [{'ref': chosen.ref}])
        self.assertNotIn('error', sub, sub)
        v = sub['handle']
        self.assertEqual((h, h), (sub['of'], sub['view']['of']))
        self.assertEqual('for the selected labels', sub['evidence']['qualified'])
        # a view stores no new body: the backing body, the view in the entry
        self.assertEqual(bodies, sorted(os.listdir(os.path.join(self.spool, 'bodies'))))
        self.assertEqual(self.store.get(h).digest, self.store.get(v).digest)
        self.assertIsInstance(self.store.view(v), GraphletView)
        walks = exhaust(self, lambda cursor: self.agent.graphlet_walks(v, side, cursor=cursor),
                        key=lambda r: r['walk'])
        self.assertEqual(carrying, {r['walk'] for r in walks})
        page = self.agent.graphlet_walks(v, side)
        self.assertEqual(('for the selected labels', h),
                         (page['evidence']['qualified'], page['evidence']['view']['of']))
        want = self.want('main', side)
        assert_evidence(self, page['evidence'], want)
        labels = exhaust(self, lambda cursor: self.agent.graphlet_labels(v, cursor=cursor))
        self.assertEqual({chosen.ref}, {r['ref'] for r in labels})
        claims = exhaust(self, lambda cursor: self.agent.graphlet_claims(
            v, side, route_consistent=False, cursor=cursor))
        self.assertEqual({chosen.ref}, {r['label']['ref'] for r in claims} | {chosen.ref})
        summ = self.agent.graphlet_summary(v)
        self.assertEqual(len(carrying), summ['summary']['arms'][side]['walks'])
        outside = sorted(set(range(n)) - carrying)
        self.assertEqual('not_in_view', self.agent.graphlet_walk(v, side, outside[0])['error'])
        self.assertEqual('not_in_view', self.agent.graphlet_sequence(
            v, side, walk=outside[0])['error'])
        inside = sorted(carrying)[0]
        self.assertNotIn('error', self.agent.graphlet_sequence(v, side, walk=inside))
        self.assertNotIn('error', self.agent.graphlet_walk(v, side, inside))
        self.assertNotIn('error', self.agent.graphlet_support(v, side, inside))
        self.assertNotIn('error', self.agent.graphlet_splits(v, side))
        fa = self.agent.graphlet_export(v, 'fasta', path='view.fa')
        both = sum(summ['summary']['arms'][s]['walks'] for s in summ['summary']['arms'])
        self.assertEqual(both, fa['records'])
        for out in (self.agent.graphlet_export(v, 'gfa', path='view.gfa'),
                    self.agent.graphlet_export(v, 'json', path='view.json'),
                    self.agent.graphlet_compare(v, h),
                    self.agent.graphlet_compare(h, v),
                    self.agent.graphlet_subtrie(v, [{'ref': chosen.ref}])):
            self.assertEqual(('view_unsupported', h), (out['error'], out['of']))
        mgt = self.agent.graphlet_export(v, 'mgt', path='view.mgt')
        with open(mgt['path'], encoding='utf-8') as f:
            j = json.loads(f.read().split('\n')[1][2:])
        self.assertEqual([{'ref': chosen.ref}], j['view']['selectors'])
        loaded = self.agent.graphlet_load('view.mgt')
        self.assertNotIn('error', loaded, loaded)
        self.assertIsInstance(self.store.view(loaded['handle']), GraphletView)
        again = exhaust(self, lambda cursor: self.agent.graphlet_walks(
            loaded['handle'], side, cursor=cursor))
        self.assertEqual(carrying, {r['walk'] for r in again})
        self.assertEqual('not_in_view', self.agent.graphlet_labels(
            v, name={'ref': next(l.ref for l in g.labels if l.ref != chosen.ref)})['error'])
        for x in (v, loaded['handle']):
            self.tools.graphlet_free(x)

    def test_17_list_free_save_load(self):
        main = self.h('main')
        r1 = self.agent.traverse_fetch(replay=main)['handle']
        digest = self.store.get(main).digest
        self.assertEqual(digest, self.store.get(r1).digest)   # entries share one body
        bodies = sorted(os.listdir(os.path.join(self.spool, 'bodies')))
        listed = exhaust(self, lambda cursor: self.agent.graphlet_list(cursor=cursor),
                         key=lambda r: r['handle'])
        self.assertEqual({e['handle'] for e in self.store.list()}, {r['handle'] for r in listed})
        small = exhaust(self, lambda cursor: self.agent.graphlet_list(cursor=cursor,
                                                                      max_bytes=600))
        self.assertEqual([r['handle'] for r in listed], [r['handle'] for r in small])
        saved = self.agent.graphlet_save(r1, 'r1.mgt')
        self.assertEqual((os.path.join(self.export_dir, 'r1.mgt'),
                          os.path.getsize(saved['path'])), (saved['path'], saved['bytes']))
        loaded = self.agent.graphlet_load('r1.mgt')
        r2 = loaded['handle']
        self.assertEqual(digest, self.store.get(r2).digest)
        self.assertEqual(bodies, sorted(os.listdir(os.path.join(self.spool, 'bodies'))))
        self.assertEqual(self.agent.graphlet_summary(main)['summary'],
                         self.agent.graphlet_summary(r2)['summary'])
        self.assertEqual(self.agent.graphlet_walks(main, self.sides('main')[0])['rows'],
                         self.agent.graphlet_walks(r2, self.sides('main')[0])['rows'])
        self.assertEqual({'freed': r1}, self.agent.graphlet_free(r1))
        gone = self.agent.graphlet_summary(r1)
        self.assertEqual(('unknown_handle', False), (gone['error'], gone['replayable']))
        self.assertTrue(os.path.exists(os.path.join(self.spool, 'bodies', digest + '.mgt')))
        self.agent.graphlet_free(r2)
        self.assertEqual('unknown_handle', self.agent.graphlet_free(r2)['error'])
        self.assertIn(main, self.store)

    def test_18_backend_errors(self):
        idx, c = self.INDEX, self.cell('main')
        seq = c.seed_sequence(0)
        for kw in ({'seed': {'sequence': seq[:20]}},
                   {'seed': {'sequence': seq, 'labels': ['NO_SUCH_LABEL_%s' % idx]}},
                   {'seed': {'sequence': seq}, 'strategy': {'bounds': {'max_extension_bp': -1}}},
                   {'seed': {'sequence': seq}, 'strategy': {'bounds': {'no_such_knob': 1}}}):
            out = self.agent.traverse_fetch(idx, **kw)
            self.assertEqual(('backend_error', 400), (out['error'], out.get('status')), kw)
        self.assertEqual(('backend_error', 400), tuple(
            self.agent.traverse_resolve(idx, sequence=seq)[k] for k in ('error', 'status')))
        res = self.agent.traverse_resolve(idx, sequence=seq, discover={'max_labels': 50},
                                          limit=3)
        self.assertNotIn('error', res, res)
        self.assertLessEqual(len(res['labels']), 3)
        caps = self.agent.traverse_capabilities(idx)
        self.assertEqual(R.identity(idx)['index_meta_fp'], caps['capabilities']['index_meta_fp'])
        for bad in ({'seed': None}, {'seed': {'labels': []}}, {'seed': {'sequence': seq},
                                                               'strategy': []}):
            self.assertEqual('bad_request', self.agent.traverse_fetch(idx, **bad)['error'])
        self.assertEqual('unknown_index', self.agent.traverse_fetch(
            'nope', seed={'sequence': seq})['error'])
        # an unreachable backend port, also behind an existing handle (its source index)
        dead = GraphletTools(self.store, {idx: TraverseClient('127.0.0.1', closed_port(),
                                                              timeout=10)})
        agent = Agent(self, dead)
        side, pid = self.walk_with_continuation('main')
        if side is None:
            side, pid = self.walk_with_continuation('lower')
            h = self.h('lower')
        else:
            h = self.h('main')
        for out in (agent.traverse_fetch(idx, seed={'sequence': seq}),
                    agent.traverse_fetch(replay=self.h('main')),
                    agent.traverse_capabilities(idx),
                    agent.traverse_resolve(idx, sequence=seq, discover={'max_labels': 5}),
                    agent.traverse_continue(h, side, pid)):
            self.assertEqual('backend_unreachable', out['error'], out)
        # the local tools still answer from the store
        self.assertNotIn('error', agent.graphlet_summary(self.h('main')))

    def test_19_oversize(self):
        idx, c = self.INDEX, self.cell('main')
        body = c.graphlet_result(0)['graphlet_bytes']
        out = self.agent.traverse_fetch(idx, seed={'sequence': c.seed_sequence(0)},
                                        strategy=self.strategy_of(c),
                                        max_graphlet_mb=(body / 2) / (1 << 20))
        self.assertEqual('spooled', out['delivery'])
        h = out['handle']
        self.assertFalse(self.store.in_ram(h))
        self.assertEqual('spooled', out['summary']['outcome']['delivery'])
        assert_evidence(self, out['evidence'], self.want('main', delivery='spooled'))
        page = self.agent.graphlet_walks(h, self.sides('main')[0])
        assert_evidence(self, page['evidence'], self.want('main', self.sides('main')[0],
                                                          delivery='spooled'))
        e = self.store.get(h)
        with open(os.path.join(self.spool, 'bodies', e.digest + '.mgt'), encoding='utf-8',
                  newline='') as f:
            self.assertEqual(c.text(0), f.read())
        self.tools.graphlet_free(h)
        # the store's hard limit rejects the fetch explicitly; nothing is stored
        spool = new_spool('hardlimit')
        try:
            store = GraphletStore(spool, max_body_mb=(body / 2) / (1 << 20))
            agent = Agent(self, GraphletTools(store, {idx: R.client(idx)}))
            out = agent.traverse_fetch(idx, seed={'sequence': c.seed_sequence(0)},
                                       strategy=self.strategy_of(c))
            self.assertEqual('too_large', out['error'])
            self.assertEqual(([], []), (store.list(), os.listdir(os.path.join(spool, 'bodies'))))
        finally:
            shutil.rmtree(spool, ignore_errors=True)
        # an answer that cannot fit its max_bytes is explicit, with the knob to raise
        tiny = self.agent.graphlet_walks(self.h('main'), self.sides('main')[0], max_bytes=300)
        self.assertEqual(('result_too_large', 300), (tiny['error'], tiny['max_bytes']))
        self.assertGreater(tiny['bytes'], 300)

    def test_20_expired_handle_replays(self):
        idx, c = self.INDEX, self.cell('main')
        clock = [1000.0]
        spool = new_spool('expiry')
        try:
            store = GraphletStore(spool, clock=lambda: clock[0], ttl_disk_s=60)
            agent = Agent(self, GraphletTools(store, {idx: R.client(idx)}))
            h = agent.traverse_fetch(idx, seed={'sequence': c.seed_sequence(0)},
                                     strategy=self.strategy_of(c))['handle']
            first = agent.graphlet_walks(h, self.sides('main')[0], n=1)
            clock[0] += 120
            store.sweep()
            gone = agent.graphlet_summary(h)
            self.assertEqual(('unknown_handle', True), (gone['error'], gone['replayable']))
            self.assertIn('traverse_fetch(replay=', gone['hint'])
            self.assertEqual('unknown_handle', agent.graphlet_walks(
                h, self.sides('main')[0], cursor=first.get('next_cursor'))['error'])
            again = agent.traverse_fetch(replay=h)
            self.assertNotIn('error', again, again)
            with open(os.path.join(spool, 'bodies', store.get(again['handle']).digest + '.mgt'),
                      encoding='utf-8', newline='') as f:
                self.assertEqual(c.text(0), f.read())
        finally:
            shutil.rmtree(spool, ignore_errors=True)


class TestAgentSessionMiniRefseq(_Session, unittest.TestCase):
    INDEX = 'mini_refseq'


class TestAgentSessionUHGG(_Session, unittest.TestCase):
    INDEX = 'uhgg'


class TestAgentSessionSRA(_Session, unittest.TestCase):
    INDEX = 'sra'

    def test_14_continue_request_is_deliverable(self):
        # fixed BUG-3 (CONTINUE-DRYRUN-TOO-LARGE): on SRA the request carries the 1000 bp
        # continuation, the continuation's labels (SRA names are ~110-byte paths) and the
        # full normalized strategy: 2127 bytes, over the 2048 list ceiling it was held to,
        # with no max_bytes to raise and no export -- unobtainable. Repro of the old answer:
        # cell sra_rand50_00__limit2, traverse_continue(<handle>, 'right', 2,
        # execute=False) -> {'error': 'result_too_large', 'bytes': 2127, 'max_bytes': 2048}.
        super().test_14_continue_request_is_deliverable()


# ------------------------------------------------------------------ the sweep

def _sweep_rows(index=None):
    return [r for r in R.cells(index, max_body_bytes=SWEEP_MAX_BODY)]


@R.skip_unless_cache()
class TestToolsOverCachedCells(unittest.TestCase):
    """Every list tool on every arm of every cached result (bodies <= 1 MB): every row
    (large pages), the first default 2 KB pages equal to them, every page within its
    ceiling, rows == total, no duplicates, evidence equal to the server-derived block,
    totals equal to the server's counts, rows equal to the detail-full JSON
    (deterministic cells with a full JSON <= 24 MB)."""

    @classmethod
    def setUpClass(cls):
        cls.spool = new_spool('sweep')
        cls.store = GraphletStore(cls.spool, max_ram_mb=256)
        cls.tools = GraphletTools(cls.store, {})

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.spool, ignore_errors=True)

    def sweep(self, index):
        rows = _sweep_rows(index)
        if not rows:
            self.skipTest('no cached cells of %s' % index)
        agent = Agent(self, self.tools)
        for row in rows:
            c = R.load_cell(index, row)
            full_ok = c.deterministic and (row.get('full') or {}).get('json_bytes', 0) \
                <= SWEEP_MAX_FULL
            for i in range(c.n_results()):
                result = c.graphlet_result(i)
                if 'graphlet' not in result:
                    continue
                with self.subTest(cell=row['cell'], result=i):
                    self.check_result(agent, c, i, result, full_ok)
            STATS['cells_swept'].add('%s/%s' % (index, row['cell']))

    def cover(self, call, key=None, n_ok=False, where='', page=None):
        """Every row of a list tool: exhausted with large pages (each page costs O(total)
        in the tool, so 2 KB pages to exhaustion would be quadratic on the big cells),
        then the first default (2 KB) pages checked against them: the same rows in the
        same order, cut only where a page says so, a cursor exactly while rows remain."""
        big = {'max_bytes': page or SWEEP_BIG_PAGE}
        if n_ok:
            big['n'] = 10 ** 7
        rows = exhaust(self, lambda cursor: call(cursor=cursor, **big), key=key)
        self.assertEqual({}, rows.cut, where)
        out, got, pages = call(cursor=None), 0, 0
        while True:
            self.assertNotIn('error', out, (where, out))
            self.assertEqual(len(rows), out['total'], where)
            cut = out['cut_fields'] if out.get('row_truncated') else None
            for r in out['rows']:
                assert_row_or_cut(self, r, rows[got], cut, (where, got))
                got += 1
            pages += 1
            STATS['pages'] += 1
            if 'next_cursor' not in out:
                self.assertEqual(len(rows), got, where)
                break
            self.assertTrue(out['rows'], where)
            if pages >= SWEEP_SMALL_PAGES:
                break
            out = call(cursor=out['next_cursor'])
        return rows

    def check_result(self, agent, c, i, result, full_ok):
        h = self.store.put(c.graphlet_response, c.request, i, source=c.index)
        envelope = c.graphlet_response
        full = c.full_result(i) if full_ok else None
        where = '%s/%d' % (c.cell, i)
        try:
            out = agent.graphlet_summary(h)
            assert_evidence(self, out['evidence'], expected_evidence(result, envelope), where)
            lims = self.cover(lambda cursor, **kw: agent.graphlet_summary(
                h, detail='limitations', cursor=cursor, **kw), where=where)
            n_lims = sum(len(a.get('limitations') or []) for a in result['arms'].values()) \
                + len(result.get('limitations') or [])
            self.assertEqual(n_lims, len(lims), where)
            labels = self.cover(lambda cursor, **kw: agent.graphlet_labels(
                h, cursor=cursor, **kw), key=lambda r: (r['ref'], r['arm']), n_ok=True,
                where=where)
            n_labels = len(full['label_dict']) if full is not None \
                else result['seed']['num_labels']
            self.assertEqual(n_labels * len(result['arms']), len(labels), where)
            for side, arm in result['arms'].items():
                want = expected_evidence(result, envelope, side)
                counts = arm['counts']
                w = (where, side)
                walks = self.cover(lambda cursor, **kw: self.ev(agent.graphlet_walks(
                    h, side, cursor=cursor, **kw), want), key=lambda r: r['walk'], n_ok=True,
                    where=w)
                self.assertEqual(counts['leaves'], len(walks), w)
                self.assertEqual(set(range(counts['leaves'])), {r['walk'] for r in walks}, w)
                keys = [(-r['n_full'], -r['n_alive'], -r['length_bp'], r['walk']) for r in walks]
                self.assertEqual(sorted(keys), keys, w)
                splits = self.cover(lambda cursor, **kw: self.ev(agent.graphlet_splits(
                    h, side, min_labels_before=0, cursor=cursor, **kw), want),
                    key=lambda r: r['segment'], n_ok=True, where=w)
                self.assertEqual(counts['splits'], len(splits), w)
                every = self.cover(lambda cursor, **kw: self.ev(agent.graphlet_claims(
                    h, side, route_consistent=False, cursor=cursor, **kw), want), where=w)
                first = agent.graphlet_claims(h, side)
                consistent = self.cover(lambda cursor, **kw: agent.graphlet_claims(
                    h, side, cursor=cursor, **kw), where=w)
                self.assertEqual(len(every), len(consistent)
                                 + sum((first.get('filtered') or {}).values()), w)
                for pid in sorted({0, counts['leaves'] - 1}):
                    # one walk's chain / support is computed whole on every page (1.4 s per
                    # call on a 10 kb beam walk): one large reference page
                    chain = self.cover(lambda cursor, **kw: self.ev(agent.graphlet_walk(
                        h, side, pid, cursor=cursor, **kw), want), key=lambda r: r['segment'],
                        where=(w, pid), page=SWEEP_WALK_PAGE)
                    support = self.cover(lambda cursor, **kw: self.ev(agent.graphlet_support(
                        h, side, pid, cursor=cursor, **kw), want), where=(w, pid),
                        page=SWEEP_WALK_PAGE)
                    pos = 0
                    for r in chain:
                        self.assertEqual(pos, r['from_bp'], (w, pid))
                        pos = r['to_bp']
                    if pos == 0:
                        # a walk of no bases (the seed at a record end): no support run
                        self.assertEqual([], support, (w, pid))
                    else:
                        self.assertEqual((0, pos), (support[0]['from_bp'], support[-1]['to_bp']),
                                         (w, pid))
                        for a, b in zip(support, support[1:]):
                            self.assertEqual(a['to_bp'], b['from_bp'], (w, pid))
                    if full is not None:
                        p = full['arms'][side]['paths'][pid]
                        self.assertEqual(p['segments'], [r['segment'] for r in chain], (w, pid))
                        self.assertEqual(p['length_bp'], pos, (w, pid))
                if full is not None:
                    fa = full['arms'][side]
                    self.assertEqual([p['length_bp'] for p in fa['paths']],
                                     [r['length_bp'] for r in sorted(walks,
                                                                     key=lambda r: r['walk'])], w)
                    self.assertEqual(
                        sorted((s['segment'], s['at_bp'], s['labels_before'],
                                [(b['char'], b['labels_distinct']) for b in s['branches']])
                               for s in fa['splits']),
                        sorted((r['segment'], r['at_bp'], r['labels_before'],
                                [(b['char'], b['labels']) for b in r['branches']])
                               for r in splits), w)
            if full is not None:
                summ = {}
                for k, l in enumerate(full['label_dict']):
                    ref = 'h:%d:%d' % (l['column'], l['seq_id']) if l['kind'] == 'header' \
                        else 'c:%d' % l['column']
                    summ[ref] = full['label_summary'][k]
                for r in labels:
                    s = summ[r['ref']][r['arm']]
                    self.assertEqual((s['direct_bp'], s['reach_bp'], s['reentries'],
                                      len(s['runs'])),
                                     (r['direct_bp'], r['reach_bp'], r['reentries'], r['runs']),
                                     (where, r['ref'], r['arm']))
        finally:
            self.store.free(h)

    def ev(self, out, want):
        if 'error' not in out:
            assert_evidence(self, out['evidence'], want)
        return out

    def test_mini_refseq(self):
        self.sweep('mini_refseq')

    def test_uhgg(self):
        self.sweep('uhgg')

    def test_sra(self):
        self.sweep('sra')


# ------------------------------------------------------------------ cursors

def _tamper(cursor, **changes):
    """A cursor with its body edited and the MAC kept (what a client could do)."""
    b64, mac = cursor.rsplit('.', 1)
    body = json.loads(base64.urlsafe_b64decode(b64 + '=' * (-len(b64) % 4)))
    body.update(changes)
    raw = json.dumps(body, separators=(',', ':')).encode()
    return base64.urlsafe_b64encode(raw).decode().rstrip('=') + '.' + mac


CURSOR_CELLS = [('sra', 'sra_hub__strict'), ('uhgg', 'uhgg_16s__cut_lists'),
                ('mini_refseq', 'mini_ndm1__limit2')]


@R.skip_unless_cache()
class TestCursorsOnRealRetrievals(unittest.TestCase):
    """A cursor is opaque and bound to the handle, the tool and the normalized
    arguments: presented with other arguments, another handle or tool, edited, forged
    with another secret, or garbage, it is rejected as bad_cursor; n and max_bytes may
    change between pages; it survives a restart of the tool layer on the same spool."""

    def setUp(self):
        self.spool = new_spool('cursors')
        self.store = GraphletStore(self.spool)
        self.tools = GraphletTools(self.store, {})
        self.agent = Agent(self, self.tools)
        self.h = {}
        for idx, cell in CURSOR_CELLS:
            if R.cell_meta(idx, cell) is None:
                continue
            c = R.load_cell(idx, cell)
            self.h[cell] = (self.store.put(c.graphlet_response, c.request, source=idx), c)
        if not self.h:
            self.skipTest('no cursor cells cached')

    def tearDown(self):
        shutil.rmtree(self.spool, ignore_errors=True)

    def calls(self, h, side):
        """(tool, call(cursor=, **overrides), variants that must reject the cursor)."""
        A = self.agent
        return [
            ('walks', lambda cursor, **o: A.graphlet_walks(h, side, cursor=cursor, **o),
             [dict(rank='length'), dict(min_bp=1), dict(spell='tail'), dict(tail_bp=59),
              dict(route_consistent=False)]),
            ('claims', lambda cursor, **o: A.graphlet_claims(
                h, side, cursor=cursor, **dict({'route_consistent': False}, **o)),
             [dict(min_bp=1), dict(route_consistent=True)]),
            ('labels', lambda cursor, **o: A.graphlet_labels(h, cursor=cursor, **o),
             [dict(rank='reach_bp'), dict(min_direct_bp=1), dict(at_bp=0), dict(arm=side)]),
            ('splits', lambda cursor, **o: A.graphlet_splits(
                h, side, cursor=cursor, **dict({'min_labels_before': 0}, **o)),
             [dict(min_labels_before=1)]),
            ('support', lambda cursor, **o: A.graphlet_support(h, side, 0, cursor=cursor, **o),
             [dict(step=7)]),
            ('walk', lambda cursor, **o: A.graphlet_walk(h, side, 0, cursor=cursor, **o), []),
        ]

    def first_cursor(self, call):
        for mb in (2048, 1024, 700):
            out = call(None, max_bytes=mb)
            if 'next_cursor' in out:
                return out
        return None

    def test_cursors_are_bound_and_tamper_proof(self):
        tested = 0
        handles = [x[0] for x in self.h.values()]
        spool2 = new_spool('cursors2')
        self.addCleanup(shutil.rmtree, spool2, True)
        key2 = GraphletTools(GraphletStore(spool2), {})._key
        for cell, (h, c) in self.h.items():
            other = next((x for x in handles if x != h), None)
            for side in c.graphlet_result(0)['arms']:
                for tool, call, variants in self.calls(h, side):
                    first = self.first_cursor(call)
                    if first is None:
                        continue
                    cur = first['next_cursor']
                    tested += 1
                    where = (cell, side, tool)
                    # the same arguments: accepted, also with another max_bytes / n
                    self.assertNotIn('error', call(cur), where)
                    self.assertNotIn('error', call(cur, max_bytes=4096), where)
                    if tool in ('walks', 'labels', 'splits'):
                        self.assertNotIn('error', call(cur, n=3), where)
                    # other arguments: rejected
                    for v in variants:
                        self.assertEqual('bad_cursor', call(cur, **v)['error'], (where, v))
                    # tampering: the body edited under the old MAC, the MAC edited, parts
                    # missing, garbage, a non-string
                    b64, mac = cur.rsplit('.', 1)
                    raw = base64.urlsafe_b64decode(b64 + '=' * (-len(b64) % 4))
                    body = json.loads(raw)
                    for bad in (cur[:-1] + ('0' if cur[-1] != '0' else '1'),
                                _tamper(cur, o=body['o'] + 1), _tamper(cur, o=0),
                                _tamper(cur, t=body['t'] + 'x'), _tamper(cur, a='0' * 16),
                                _tamper(cur, h='g_000000000000'), b64, '.' + mac,
                                'not-a-cursor', '%%%.zz', 'A' * 200000 + '.' + 'b' * 16,
                                cur + 'x'):
                        self.assertNotEqual(cur, bad)
                        self.assertEqual('bad_cursor', call(bad)['error'], (where, bad[:40]))
                    for bad in (12345, ['x'], {'o': 1}):
                        self.assertEqual('bad_cursor', call(bad)['error'], (where, bad))
                    # the real body MACed with this server's key is accepted (control),
                    # with another spool's key rejected
                    same = b64 + '.' + hmac.new(self.tools._key, raw,
                                                hashlib.sha256).hexdigest()[:16]
                    self.assertNotIn('error', call(same), where)
                    forged = b64 + '.' + hmac.new(key2, raw, hashlib.sha256).hexdigest()[:16]
                    self.assertEqual('bad_cursor', call(forged)['error'], where)
                    # another handle, the same tool and arguments: rejected
                    if other is not None and tool == 'labels':
                        self.assertEqual('bad_cursor', self.agent.graphlet_labels(
                            other, cursor=cur)['error'])
        self.assertGreater(tested, 5)

    def test_cursors_cross_tools_are_rejected(self):
        (h, c), = list(self.h.values())[:1]
        side = next(iter(c.graphlet_result(0)['arms']))
        outs = {}
        for tool, call, _ in self.calls(h, side):
            first = self.first_cursor(call)
            if first is not None:
                outs[tool] = (call, first['next_cursor'])
        for t1, (call, _) in outs.items():
            for t2, (_, cur) in outs.items():
                if t1 != t2:
                    self.assertEqual('bad_cursor', call(cur)['error'], (t1, t2))

    def test_cursor_survives_a_restart_and_n_may_change(self):
        h, c = next(iter(self.h.values()))
        first = self.agent.graphlet_labels(h, n=2)
        self.assertIn('next_cursor', first)
        restarted = Agent(self, GraphletTools(GraphletStore(self.spool), {}))
        second = restarted.graphlet_labels(h, n=5, cursor=first['next_cursor'])
        self.assertNotIn('error', second)
        rows = exhaust(self, lambda cursor: self.agent.graphlet_labels(h, cursor=cursor))
        self.assertEqual(rows[2:7], second['rows'])

    def test_a_page_past_the_end_is_refused(self):
        h, c = next(iter(self.h.values()))
        side = next(iter(c.graphlet_result(0)['arms']))
        total = self.agent.graphlet_walks(h, side)['total']
        args = dict(arm=side, rank='support', min_bp=0, label=None, spell='none', tail_bp=60,
                    route_consistent=True, cursor=None)
        past = self.tools._cursor('graphlet_walks', h, args, total + 5)
        self.assertEqual('bad_cursor', self.agent.graphlet_walks(h, side, cursor=past)['error'])
        end = self.tools._cursor('graphlet_walks', h, args, total)
        out = self.agent.graphlet_walks(h, side, cursor=end)
        self.assertEqual(([], total), (out['rows'], out['total']))


# ------------------------------------------------------------------ offline errors

@R.skip_unless_cache('mini_refseq')
class TestStructuredErrorsOffline(unittest.TestCase):
    """The failures an agent meets without a backend, as structured results."""

    def setUp(self):
        R.require_cache('mini_refseq', 'mini_ndm1__limit2')
        self.spool = new_spool('errors')
        self.store = GraphletStore(self.spool)
        self.tools = GraphletTools(self.store, {'mini_refseq': TraverseClient(
            '127.0.0.1', closed_port(), timeout=5)})
        self.agent = Agent(self, self.tools)
        self.c = R.load_cell('mini_refseq', 'mini_ndm1__limit2')
        self.h = self.store.put(self.c.graphlet_response, self.c.request, source='mini_refseq')

    def tearDown(self):
        shutil.rmtree(self.spool, ignore_errors=True)

    def test_unknown_handles(self):
        A = self.agent
        for h in ('g_000000000000', 'g_' + 'f' * 12, '../../etc/passwd', 'g_', '',
                  'g_00000000000/../x'):
            for out in (A.graphlet_summary(h), A.graphlet_walks(h, 'left'),
                        A.graphlet_labels(h), A.graphlet_free(h), A.graphlet_save(h, 'x.mgt'),
                        A.graphlet_compare(self.h, h), A.traverse_fetch(replay=h)):
                self.assertEqual('unknown_handle', out['error'], (h, out))
                self.assertFalse(out['replayable'])
        for h in (None, 7, ['g_x'], {'h': 1}):
            self.assertEqual('bad_argument', A.graphlet_summary(h)['error'], h)
        self.assertEqual('bad_arm', A.graphlet_walks(self.h, 'up')['error'])

    def test_unknown_and_malformed_labels(self):
        A = self.agent
        for sel in ({'ref': 'c:999999'}, {'ref': 'h:0:999999'}, {'name': 'NZ_NOPE.1'}, 'NZ_NOPE.1',
                    {'id': 10 ** 6}):
            for out in (A.graphlet_labels(self.h, name=sel),
                        A.graphlet_walks(self.h, 'left', label=sel),
                        A.graphlet_claims(self.h, 'left', labels=[sel]),
                        A.graphlet_subtrie(self.h, [sel])):
                self.assertEqual('unknown_label', out['error'], (sel, out))
                self.assertIn('graphlet_labels', out['hint'])
        for sel in ({'ref': 'c:0', 'name': 'x'}, {'tag': 'x'}, 3.5, [[1]]):
            out = A.graphlet_labels(self.h, name=sel)
            self.assertEqual('bad_argument', out['error'], (sel, out))
        self.assertEqual('bad_argument', A.graphlet_subtrie(self.h, [])['error'])

    def test_missing_envelope(self):
        """A body-only graphlet (a spool body loaded by path) answers what the body
        holds and refuses what needs the envelope, explicitly."""
        A = self.agent
        body = os.path.join(self.spool, 'bodies', self.store.get(self.h).digest + '.mgt')
        loaded = A.graphlet_load(body)
        self.assertIsNone(loaded['summary'])
        b = loaded['handle']
        self.assertEqual('missing_envelope', A.graphlet_summary(b)['error'])
        self.assertEqual('missing_envelope', A.graphlet_export(b, 'json', path='b.json')['error'])
        g = self.store.graphlet(self.h)
        a = g.arms['right']
        pid = next(p.id for p in derive.paths(a)
                   if a.segments[p.leaf].leaf.continuation is not None)
        self.assertEqual('missing_envelope', A.traverse_continue(b, 'right', pid,
                                                                 execute=False)['error'])
        self.assertEqual('not_replayable', A.traverse_fetch(replay=b)['error'])
        walks = A.graphlet_walks(b, 'right')
        self.assertNotIn('error', walks)
        self.assertEqual(self.c.graphlet_result(0)['outcome'], walks['evidence']['outcome'])
        self.assertEqual(A.graphlet_walks(self.h, 'right')['rows'], walks['rows'])
        self.assertNotIn('error', A.graphlet_export(b, 'fasta', path='b.fa'))
        self.assertNotIn('error', A.graphlet_export(b, 'mgt', path='b.mgt'))

    def test_unreachable_backend(self):
        A = self.agent
        seq = self.c.seed_sequence(0)
        for out in (A.traverse_fetch('mini_refseq', seed={'sequence': seq}),
                    A.traverse_capabilities('mini_refseq'),
                    A.traverse_resolve('mini_refseq', sequence=seq, discover={'max_labels': 3}),
                    A.traverse_fetch(replay=self.h),
                    A.traverse_continue(self.h, 'right', next(
                        p.id for p in derive.paths(self.store.graphlet(self.h).arms['right'])
                        if self.store.graphlet(self.h).arms['right'].segments[p.leaf].leaf
                        .continuation is not None))):
            self.assertEqual('backend_unreachable', out['error'], out)
        self.assertEqual(1, len(self.store.list()))

    def test_result_too_large_names_the_knob(self):
        for mb in (200, 300, 500):
            for out in (self.agent.graphlet_walks(self.h, 'left', max_bytes=mb),
                        self.agent.graphlet_claims(self.h, 'left', max_bytes=mb),
                        self.agent.graphlet_summary(self.h, max_bytes=mb)):
                if 'error' in out:
                    self.assertEqual(('result_too_large', mb), (out['error'], out['max_bytes']))
                    self.assertIn('max_bytes', out['message'])

    def test_every_return_fits_even_a_tiny_max_bytes(self):
        # fixed BUG-4: "EVERY return is <= max_bytes" did not hold below ~180 bytes: the
        # result_too_large the wrapper substitutes was itself 180+ bytes and returned
        # unchecked (graphlet_walks(<mini_ndm1__limit2>, 'left', max_bytes=64) -> 182 B).
        # Now every return is <= max(max_bytes, MIN_MAX_BYTES); the Agent wrapper asserts it.
        for mb in (MIN_MAX_BYTES, 80, 128, 182):
            for out in (self.agent.graphlet_walks(self.h, 'left', max_bytes=mb),
                        self.agent.graphlet_summary(self.h, max_bytes=mb),
                        self.agent.graphlet_claims(self.h, 'left', max_bytes=mb)):
                if 'error' in out:
                    self.assertEqual('result_too_large', out['error'], out)
        for mb in (1, 10, MIN_MAX_BYTES - 1):
            out = self.agent.graphlet_walks(self.h, 'left', max_bytes=mb)
            self.assertEqual('bad_argument', out['error'], out)

    def test_continue_request_hint_is_actionable(self):
        # fixed BUG-3 (the other half): result_too_large told the agent to "raise
        # max_bytes", but traverse_continue had no max_bytes parameter: passing one raised
        # TypeError out of the tool. Now it takes one, and an argument a tool does not take
        # is a bad_argument result.
        g = self.store.graphlet(self.h)
        a = g.arms['right']
        pid = next(p.id for p in derive.paths(a)
                   if a.segments[p.leaf].leaf.continuation is not None)
        out = self.agent.traverse_continue(self.h, 'right', pid, execute=False,
                                           max_bytes=1 << 16)
        self.assertEqual(g.next_request('right', [pid]), out['request'])
        small = self.agent.traverse_continue(self.h, 'right', pid, execute=False,
                                             max_bytes=256)
        self.assertEqual('result_too_large', small['error'])
        self.assertIn('raise max_bytes', small['message'])
        self.assertEqual('bad_argument', self.agent.graphlet_free(self.h, max_bytes=4096)
                         ['error'])

    def test_loading_files(self):
        A = self.agent
        exports = os.path.join(self.spool, 'exports')
        os.makedirs(exports, exist_ok=True)
        for name, data in (('empty.mgt', b''), ('binary.mgt', b'\x00\xff\xfe' * 100),
                           ('latin1.mgt', self.c.text(0).encode('utf-8').replace(
                               b'NZ_', b'N\xd6_', 1)),
                           ('json.mgt', json.dumps(self.c.graphlet_response).encode()),
                           ('cut.mgt', self.c.text(0).encode()[:5000])):
            with open(os.path.join(exports, name), 'wb') as f:
                f.write(data)
            out = A.graphlet_load(name)
            self.assertEqual('format_error', out['error'], (name, out))
            self.assertNotIn('NZ_', out['message'])
        self.assertEqual('io_error', A.graphlet_load('missing.mgt')['error'])
        self.assertEqual('path_not_allowed', A.graphlet_load('/etc/hosts')['error'])
        self.assertEqual('path_not_allowed', A.graphlet_load('../../x.mgt')['error'])
        self.assertEqual(1, len(self.store.list()))


def tearDownModule():
    if STATS['tool_calls']:
        sys.stderr.write('\n[test_real_agent] %d tool calls, %d pages, %d cells swept, '
                         '%d session cells\n' % (STATS['tool_calls'], STATS['pages'],
                                                 len(STATS['cells_swept']),
                                                 len(STATS['session_cells'])))


if __name__ == '__main__':
    unittest.main()
