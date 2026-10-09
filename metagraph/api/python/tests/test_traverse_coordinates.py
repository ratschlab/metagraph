"""Record coordinates in the library (feature level 6; DESIGN §18 and §26), on real
mini_refseq responses committed under data/traverse/coords/
(scripts/traversal/traverse_coords_snapshots.py):

  * parsing: the block is validated eagerly against the body by from_response(), a saved
    file's J line and the store's first parse of an unparsed entry; the parsed index lives
    in the graphlet's cache, never in a slot, and survives every round trip;
  * every inconsistency is a GraphletFormatError 'coordinates: ...' (one test per rule);
  * claims (clipped to the cut), walks and label walks carry them -- CoordClaim, CoordWalk,
    CoordLabelWalk, only for graphlets with a block; compare() never clips them and says
    it does not compare them (coordinates beside the claims, never in their keys);
    to_json() equals the full result; FASTA headers state them 1-based closed (none in GFA);
  * requests: next_request() carries the setting and strips the cap when coordinates are
    dropped (the cap needs coordinates: true); build_request(), supports_coordinates() and
    the automatic rule of traverse() (the depth gate);
  * the MCP rows, the local-limits charges (W_COORD) and the traced peak, also of blocks of
    10^4 occurrences (widened()), and one charge for the index whether cached or not;
  * the partial delivery of a derived seed: a walked result with a derivation limitation is
    qualified, check_rules agrees, the summary names it; the majority-parent rule: from
    feature level 6 a merge's first parent carries the most labels (check_rules), and
    compare() qualifies a comparison across the rule change against the level-5 documents
    of data/traverse/r21; and an ambient budget under tools without local limits, a list
    stopped between rows included (GraphletView.walks' set).
"""

import copy
import gc
import gzip
import json
import os
import tempfile
import tracemalloc
import unittest

import traverse_testlib as T
import traverse_golden as TG

import metagraph.traverse as MT
from metagraph.traverse import coords, derive, export, ops
from metagraph.traverse._codec import UNLIMITED, GraphletFormatError
from metagraph.traverse.budget import LocalBudget, LocalBudgetExceeded, LocalLimits
from metagraph.traverse.client import (AUTO_COORDINATES_UNDER_MEMORY_BUDGET, TraverseClient,
                                       TraverseError, auto_coordinates,
                                       drop_coordinates_note)
from metagraph.traverse.mcp_tools import GraphletTools, ToolLimits
from metagraph.traverse.model import (Claim, CoordClaim, CoordLabelWalk, CoordWalk,
                                      LabelWalk, MissingEnvelope, Walk)
from metagraph.traverse.store import GraphletStore

COORDS = os.path.join(T.DATA, 'coords')
# the cells whose graphlet and full fetches stop at the same place (the partial-derivation
# cells state an elapsed time, and the memory stop is priced per detail)
DETERMINISTIC = ('rep0_cut1', 'rep0_unlimited_left', 'rep0_lower_bound', 'rep0_column',
                 'rep0_mixed', 'rep0_kmer_null')
BLOCKS = ('rep0_cut1', 'rep0_unlimited_left', 'rep0_lower_bound', 'rep0_column',
          'rep0_mixed', 'd3_memory_stop')


def resp(name, detail='graphlet'):
    with open(os.path.join(COORDS, '%s.%s.json.gz' % (name, detail)), 'rb') as f:
        return json.loads(gzip.decompress(f.read()))


def fresh(name, mutate=None):
    """from_response of a snapshot; |mutate|(result, response) edits copies first."""
    j = resp(name)
    r = j['results'][0]
    if mutate is not None:
        mutate(r, j)
    return MT.from_response(r, j)


def stripped(name):
    """The snapshot as it would be without coordinates: the block, the reason, the echo
    and the coordinates limitation (summary and K record) removed -- what an opt-out
    request returns (the stripped-equals-opt-out harness)."""
    j = resp(name)
    r = j['results'][0]
    r.pop('coordinates', None)
    r.pop('coordinates_reason', None)
    r['limitations'] = [l for l in r['limitations'] if l['kind'] != 'coordinates']
    out = j['strategy']['output']
    out.pop('coordinates', None)
    out.pop('max_coordinate_occurrences', None)
    lines = r['graphlet'].split('\n')
    keep = [x for x in lines if not x.startswith('K * coordinates ')]
    if len(keep) != len(lines):
        keep[-2] = 'Z %d' % (len(keep) - 1)
    r['graphlet'] = '\n'.join(keep)
    r['graphlet_bytes'] = len(r['graphlet'].encode())
    r['graphlet_lines'] = r['graphlet'].count('\n')
    return MT.from_response(r, j)


def widened(name, n, labels=None):
    """The snapshot's response with |n| occurrences in every list of |labels| (None: all
    seed labels) and the cap "unlimited": a block of 10^4 occurrences and more, as the
    wide fixture of DESIGN §26 gives, consistent with the body -- each seed-entered run
    continues every seed occurrence of its label, a switch-entered run has its own starts --
    so from_response() accepts it. The lists of the other labels stay."""
    j = resp(name)
    r = j['results'][0]
    c = r['coordinates']
    g = MT.from_response(r, j)
    sl = g.seed.length_bp
    starts = {}
    for e in c['seed']:
        if labels is None or e['label'] in labels:
            st = starts[e['label']] = [10 ** 7 * (e['label'] + 1) + i * 1009 for i in range(n)]
            e['occurrences'] = [[s, s + sl] for s in st]
            e.pop('occurrences_total', None)
    for side, lst in c['arms'].items():
        for e, run in zip(lst, g.arms[side].runs):
            if run.label not in starts:
                continue
            w = run.to_bp - run.from_bp
            e.pop('occurrences_total', None)
            if run.from_label is None and side == 'right':
                e['occurrences'] = [[s + sl + run.from_bp, s + sl + run.to_bp]
                                    for s in starts[run.label]]
            elif run.from_label is None:
                e['occurrences'] = [[s - run.from_bp - w, s - run.from_bp]
                                    for s in starts[run.label]]
            else:
                e['occurrences'] = [[5 * 10 ** 8 + i * 1013, 5 * 10 ** 8 + i * 1013 + w]
                                    for i in range(n)]
    c['max_occurrences'] = 'unlimited'
    j['strategy']['output']['max_coordinate_occurrences'] = 'unlimited'
    # no list cut: complete unless a run is a lower bound, and no K record
    c['complete'] = not c.get('runs_lower_bound')
    r['limitations'] = [l for l in r['limitations'] if l['kind'] != 'coordinates']
    lines = r['graphlet'].split('\n')
    keep = [x for x in lines if not x.startswith('K * coordinates ')]
    if len(keep) != len(lines):
        keep[-2] = 'Z %d' % (len(keep) - 1)
    r['graphlet'] = '\n'.join(keep)
    r['graphlet_bytes'] = len(r['graphlet'].encode())
    r['graphlet_lines'] = r['graphlet'].count('\n')
    return j


def _peak(fn):
    gc.collect()
    tracemalloc.start()
    base = tracemalloc.get_traced_memory()[0]
    try:
        out = fn()
        return tracemalloc.get_traced_memory()[1] - base, out
    finally:
        tracemalloc.stop()


# ======================================================================= parsing

class TestParse(unittest.TestCase):
    def test_every_cell(self):
        want = {'rep0_cut1': ('record', False, 16, 0),
                'rep0_unlimited_left': ('record', True, 0, 0),
                'rep0_lower_bound': ('record', False, 0, 1),
                'rep0_column': ('column', True, 0, 0),
                'rep0_mixed': ('mixed', False, 0, 1),
                'd3_memory_stop': ('record', True, 0, 0)}
        for name, (kind, complete, cut, lower) in want.items():
            g = fresh(name)
            c = g.coordinates
            self.assertEqual((kind, complete, cut, lower),
                             (c.kind, c.complete, c.lists_cut, c.runs_lower_bound), name)
            self.assertIsNone(g.coordinates_reason)
            self.assertEqual(set(g.arms), set(c.arms))
            for side, runs in c.arms.items():
                self.assertEqual(len(g.arms[side].runs), len(runs))
                for r, rc in zip(g.arms[side].runs, runs):
                    self.assertEqual((r.id, r.label, r.from_bp, r.to_bp),
                                     (rc.run, rc.label, rc.from_bp, rc.to_bp))
            self.assertEqual([], g.check_rules(), name)
        self.assertIs(UNLIMITED, fresh('rep0_unlimited_left').coordinates.max_occurrences)

    def test_reasons(self):
        self.assertEqual('support kmer', fresh('rep0_kmer_null').coordinates_reason)
        self.assertEqual('partial derivation', fresh('d3_trace_coordinates').coordinates_reason)
        self.assertEqual('not requested', fresh('d3_kmer').coordinates_reason)
        for name in ('rep0_kmer_null', 'd3_trace_coordinates', 'd3_kmer'):
            self.assertIsNone(fresh(name).coordinates)

    def test_the_index_is_cached_not_a_slot(self):
        g = fresh('rep0_cut1')
        self.assertIs(g.coordinates, g.cache[coords.CACHE_KEY])
        self.assertNotIn('coordinates', [f for f in MT.Graphlet.__slots__])
        # a graphlet without a block caches nothing: memory_bytes and the signature as
        # they always were
        o = fresh('d3_kmer')
        before = (ops.memory_bytes(o), ops.cache_signature(o))
        self.assertIsNone(o.coordinates)
        self.assertIsNone(o.coordinates_reason and None)
        self.assertEqual(before, (ops.memory_bytes(o), ops.cache_signature(o)))
        self.assertNotIn(coords.CACHE_KEY, o.cache)

    def test_round_trips(self):
        with tempfile.TemporaryDirectory() as d:
            for name in BLOCKS + ('rep0_kmer_null', 'd3_trace_coordinates'):
                g = fresh(name)
                text = g.dump(envelope=True)
                self.assertTrue(MT.is_canonical(text))
                g2 = MT.parse(text)
                self.assertEqual(g.coordinates, g2.coordinates, name)
                self.assertEqual(g.coordinates_reason, g2.coordinates_reason)
                if g.coordinates is not None:
                    # validated at the parse: in the cache before any query
                    self.assertIn(coords.CACHE_KEY, g2.cache)
                path = os.path.join(d, name + '.mgt')
                g.save(path)
                self.assertEqual(g.coordinates, MT.load(path).coordinates)
                j = resp(name)
                self.assertEqual(text, MT.standalone_text(j['results'][0]['graphlet'],
                                                          j['results'][0], j))

    def test_store(self):
        with tempfile.TemporaryDirectory() as d:
            store = GraphletStore(d, max_ram_mb=0)
            j = resp('rep0_lower_bound')
            h = store.put(j, {'seeds': []})
            g = store.graphlet(h)                     # parsed again (not resident)
            self.assertEqual(fresh('rep0_lower_bound').coordinates, g.coordinates)
            v = ops.GraphletView(g, [g.labels[0].ref])
            self.assertEqual(fresh('rep0_lower_bound').coordinates,
                             v.backing.coordinates)
            # an unparsed entry (a parse stopped on its budget) validates its block at its
            # first whole parse; an inconsistent one is refused there and stays unparsed
            r = copy.deepcopy(j['results'][0])
            r['coordinates']['k'] = 29
            bad = dict(j, results=[r])
            h2 = store.put_unparsed(r, bad, {'seeds': []})
            with self.assertRaises(GraphletFormatError):
                store.graphlet(h2)
            self.assertFalse(store.get(h2).parsed)
            h3 = store.put_unparsed(j['results'][0], j, {'seeds': []})
            self.assertEqual(fresh('rep0_lower_bound').coordinates,
                             store.graphlet(h3).coordinates)
            self.assertTrue(store.get(h3).parsed)

    def test_body_only(self):
        g = MT.parse(fresh('rep0_cut1').dump(envelope=False))
        self.assertIsNone(g.coordinates)
        self.assertEqual('no envelope', g.coordinates_reason)
        self.assertNotIsInstance(g.claims()[0], CoordClaim)
        with self.assertRaises(MissingEnvelope):
            g.to_fasta(coordinates=True)


# ============================================================ the rules (eager validation)

def _edit(fn):
    def mutate(r, j):
        fn(r['coordinates'], r, j)
    return mutate


class TestRejects(unittest.TestCase):
    """One inconsistency per rule: GraphletFormatError 'coordinates: ...', eagerly."""

    def assertRefused(self, name, fn, needle):
        with self.assertRaises(GraphletFormatError) as cm:
            fresh(name, _edit(fn))
        self.assertIn('coordinates: ', str(cm.exception))
        self.assertIn(needle, str(cm.exception))

    def first_run(self, c, side='right'):
        return c['arms'][side][0]

    def test_kind_k_and_cap(self):
        self.assertRefused('rep0_cut1', lambda c, r, j: c.update(kind='column'), 'kind')
        self.assertRefused('rep0_cut1', lambda c, r, j: c.update(kind='genome'), 'kind')
        self.assertRefused('rep0_cut1', lambda c, r, j: c.update(k=21), 'k 21')
        self.assertRefused('rep0_cut1', lambda c, r, j: c.update(k=True), 'k is not')
        self.assertRefused('rep0_cut1', lambda c, r, j: c.update(max_occurrences=0),
                           'max_occurrences')
        self.assertRefused('rep0_cut1', lambda c, r, j: c.update(max_occurrences=2),
                           'echoes')
        self.assertRefused('rep0_cut1', lambda c, r, j: c.update(complete='no'), 'complete')

    def test_seed_entries(self):
        self.assertRefused('rep0_cut1', lambda c, r, j: c['seed'].pop(), 'seed entries')
        self.assertRefused('rep0_cut1', lambda c, r, j: c['seed'].reverse(), 'ascending')
        self.assertRefused(
            'rep0_unlimited_left',
            lambda c, r, j: c['seed'][0]['occurrences'][0].__setitem__(1, 1), 'bases')

    def test_one_entry_per_run_in_order(self):
        self.assertRefused('rep0_cut1', lambda c, r, j: c['arms']['right'].pop(), 'entries')
        self.assertRefused('rep0_cut1', lambda c, r, j: c['arms'].pop('left'), 'arms')
        self.assertRefused('rep0_cut1', lambda c, r, j: self.first_run(c).update(run=1),
                           'R order')
        self.assertRefused('rep0_cut1', lambda c, r, j: self.first_run(c).update(label=5),
                           'label')
        self.assertRefused('rep0_cut1',
                           lambda c, r, j: self.first_run(c).update(to_bp=599), 'to_bp')

    def test_intervals(self):
        def length(c, r, j):
            iv = self.first_run(c)['occurrences'][0]
            iv[1] += 1
        self.assertRefused('rep0_unlimited_left', lambda c, r, j: c['arms']['left'][0][
            'occurrences'][0].__setitem__(1, -1), 'integer')
        self.assertRefused('rep0_cut1', length, 'bases, not')

        def descending(c, r, j):
            for e in c['arms']['left']:
                if len(e['occurrences']) > 1:
                    e['occurrences'].reverse()
                    return
            raise AssertionError('no list of two')
        self.assertRefused('rep0_unlimited_left', descending, 'ascend')
        self.assertRefused('rep0_cut1', lambda c, r, j: self.first_run(c).update(
            occurrences=[[1.5, 601.5]]), 'integer')

    def test_totals_exactly_on_the_cut_lists(self):
        def total_on_uncut(c, r, j):
            e = next(e for e in c['arms']['left'] if e['occurrences'])
            e['occurrences_total'] = len(e['occurrences']) + 1
        # "unlimited": no list can be cut
        self.assertRefused('rep0_unlimited_left', total_on_uncut, 'not cut')

        def total_below_the_list(c, r, j):
            e = next(e for e in c['arms']['left'] if 'occurrences_total' not in e
                     and e['occurrences'])
            e['occurrences_total'] = 1
        self.assertRefused('rep0_cut1', total_below_the_list, 'not cut')

        def cut_without_total(c, r, j):
            e = next(e for e in c['arms']['right'] if 'occurrences_total' in e)
            del e['occurrences_total']
        # the list is then whole, and the K record counts one list too many
        self.assertRefused('rep0_cut1', cut_without_total, 'limitation')

        def above_cap(c, r, j):
            c['max_occurrences'] = 1
            j['strategy']['output']['max_coordinate_occurrences'] = 1
        self.assertRefused('rep0_unlimited_left', above_cap, 'above the cap')

    def test_chains_ended_and_lower_bound(self):
        self.assertRefused('rep0_cut1', lambda c, r, j: self.first_run(c).update(
            chains_ended=0), 'chains_ended')
        self.assertRefused('rep0_cut1', lambda c, r, j: self.first_run(c).update(
            lower_bound=True), 'entered by the seed')
        self.assertRefused('rep0_cut1', lambda c, r, j: self.first_run(c).update(
            lower_bound=False), 'only as true')

        def unmark(c, r, j):
            for side in c['arms'].values():
                for e in side:
                    e.pop('lower_bound', None)
        self.assertRefused('rep0_lower_bound', unmark, 'runs_lower_bound')

        def recount(c, r, j):
            c['runs_lower_bound'] = 2
        self.assertRefused('rep0_lower_bound', recount, 'runs_lower_bound')

    def test_complete_and_the_limitation(self):
        self.assertRefused('rep0_cut1', lambda c, r, j: c.update(complete=True), 'complete')
        self.assertRefused('rep0_lower_bound', lambda c, r, j: c.update(complete=True),
                           'complete')

        def no_k(c, r, j):
            lines = r['graphlet'].split('\n')
            keep = [x for x in lines if not x.startswith('K * coordinates ')]
            keep[-2] = 'Z %d' % (len(keep) - 1)
            r['graphlet'] = '\n'.join(keep)
            r.pop('graphlet_bytes')
            r.pop('graphlet_lines')
        self.assertRefused('rep0_cut1', no_k, 'coordinates limitations')

        def other_count(c, r, j):
            r['graphlet'] = r['graphlet'].replace('lists_cut=i:16', 'lists_cut=i:15')
            r.pop('graphlet_bytes')
        self.assertRefused('rep0_cut1', other_count, 'does not state')

    def test_the_seed_continuity(self):
        def shift(c, r, j):
            # every occurrence of one seed-entered run moved by one record base (lengths
            # still right): it does not continue a seed occurrence
            for e in c['arms']['left']:
                if e['occurrences'] and e['from_bp'] == 0:
                    e['occurrences'] = [[s - 1, t - 1] for s, t in e['occurrences']]
                    return
        self.assertRefused('rep0_unlimited_left', shift, 'continues no seed occurrence')

    def test_presence_against_the_echo_and_the_retrieval(self):
        def unasked(c, r, j):
            j['strategy']['output'].pop('coordinates')
            j['strategy']['output'].pop('max_coordinate_occurrences')
        self.assertRefused('rep0_cut1', unasked, 'did not ask')

        def lost(r, j):
            r.pop('coordinates')
        with self.assertRaises(GraphletFormatError):
            fresh('rep0_cut1', lost)

        def both(c, r, j):
            r['coordinates_reason'] = 'support kmer'
        self.assertRefused('rep0_cut1', both, 'a block and')

        def reason_lost(r, j):
            r.pop('coordinates_reason')
        with self.assertRaises(GraphletFormatError):
            fresh('rep0_kmer_null', reason_lost)

        def wrong_reason(r, j):
            r['coordinates_reason'] = 'partial derivation'
        with self.assertRaises(GraphletFormatError):
            fresh('rep0_kmer_null', wrong_reason)

        def block_on_kmer(r, j):
            r['coordinates'] = resp('rep0_cut1')['results'][0]['coordinates']
            r.pop('coordinates_reason')
        with self.assertRaises(GraphletFormatError) as cm:
            fresh('rep0_kmer_null', block_on_kmer)
        self.assertIn('support trace only', str(cm.exception))

    def test_a_reason_the_library_does_not_know_is_kept(self):
        def later(r, j):
            r['coordinates_reason'] = 'a reason of a later level'
        g = fresh('rep0_kmer_null', later)
        self.assertEqual('a reason of a later level', g.coordinates_reason)

    def test_a_saved_file_is_checked_too(self):
        text = fresh('rep0_cut1').dump(envelope=True)
        head, j, rest = text.split('\n', 2)
        obj = json.loads(j[2:])
        obj['results'][0]['coordinates']['k'] = 30
        bad = '\n'.join([head, 'J ' + json.dumps(obj, sort_keys=True,
                                                  separators=(',', ':')), rest])
        with self.assertRaises(GraphletFormatError):
            MT.parse(bad)


class TestGolden(unittest.TestCase):
    """Every public operation and tool over the snapshots: the digests recorded with this
    library (data/traverse/golden/coordinates.json.gz; the opt-out gate over the committed
    real fixtures, committed.json.gz, is test_traverse_speedups's and is unchanged). An
    intended change of a coordinate output records them again (traverse_golden.py --out
    data/traverse/golden/coordinates.json.gz data/traverse/coords
    data/traverse/real/mini_refseq/*__trace_coords*.graphlet.json.gz)."""

    def test_digests(self):
        want = TG.read_json(os.path.join(TG.GOLDEN_DIR, 'coordinates.json.gz'))
        got = TG.all_digests(MT, TG.items_of(TG.coordinate_files()))
        self.assertEqual(sorted(want), sorted(got))
        bad = sorted(k for k in want if got[k] != want[k])
        self.assertEqual([], bad[:20], '%d of %d digests differ' % (len(bad), len(want)))


# ======================================================================= queries

class TestClaims(unittest.TestCase):
    def test_claims_carry_their_runs_coordinates(self):
        for name in BLOCKS:
            g = fresh(name)
            c = g.coordinates
            for side in g.arms:
                for cl in g.claims(side):
                    self.assertIsInstance(cl, CoordClaim)
                    self.assertEqual(c.run(side, cl.run), cl.coordinates)
                    d = cl.as_dict()
                    self.assertEqual([list(x) for x in cl.coordinates.intervals],
                                     d['coordinates']['intervals'])
                    self.assertEqual(coords.label_kind(g, cl.label.id),
                                     d['coordinates']['kind'])

    def test_a_cut_clips_both_arms(self):
        g = fresh('rep0_cut1')
        for side in g.arms:
            whole = {cl.run: cl for cl in g.claims(side)}
            for cl in g.claims(side, at_most_bp=100):
                rc, full = cl.coordinates, whole[cl.run].coordinates
                d = min(full.to_bp, 100) - full.from_bp
                self.assertEqual(full.from_bp + d, rc.to_bp)
                for (s, e), (fs, fe) in zip(rc.intervals, full.intervals):
                    self.assertEqual(d, e - s)
                    if side == 'right':
                        self.assertEqual(fs, s)       # walking outward = ascending
                    else:
                        self.assertEqual(fe, e)       # walking outward = descending
                self.assertEqual(full.total, rc.total)

    def test_claims_of_a_graphlet_without_coordinates_are_claims(self):
        for name in ('rep0_kmer_null', 'd3_kmer', 'd3_trace_coordinates'):
            for cl in fresh(name).claims(strict=False):
                self.assertIs(Claim, type(cl))
                self.assertNotIn('coordinates', cl.as_dict())

    def test_stripped_claims_equal_the_opt_out_claims(self):
        for name in ('rep0_cut1', 'rep0_lower_bound', 'rep0_mixed'):
            a = [c.as_dict() for c in fresh(name).claims()]
            for d in a:
                d.pop('coordinates')
            self.assertEqual([c.as_dict() for c in stripped(name).claims()], a)

    def test_walks_and_label_walks(self):
        for name in BLOCKS:
            g = fresh(name)
            c = g.coordinates
            for side in g.arms:
                a = g.arms[side]
                for w in g.walks(side):
                    self.assertIsInstance(w, CoordWalk)
                    ends = derive.end_labels(a, w.leaf)
                    self.assertEqual({g.labels[e.label].ref: c.run(side, e.run)
                                      for e in ends}, w.coordinates)
                    for cl in w.claims:
                        self.assertIsInstance(cl, CoordClaim)
                for lab in g.labels[:3]:
                    for lw in g.label_walks(lab, side):
                        self.assertIsInstance(lw, CoordLabelWalk)
                        self.assertEqual(c.run(side, lw.run), lw.coordinates)
        o = fresh('d3_kmer')
        self.assertIs(Walk, type(o.walks('right')[0]))
        self.assertIs(LabelWalk, type(o.label_walks(o.labels[0], 'right')[0]))

    def test_the_whole_seed_entered_walk_lies_after_the_seed(self):
        # a run entered by the seed at 0 continues a seed occurrence: right arm right
        # after it, left arm right before it (the interval rule)
        g = fresh('rep0_unlimited_left')
        c = g.coordinates
        n = g.seed.length_bp
        for rc in c.arms['left']:
            run = g.arms['left'].runs[rc.run]
            if run.from_label is None and rc.from_bp == 0:
                starts = {s for s, _ in c.seed[rc.label].intervals}
                for s, e in rc.intervals:
                    self.assertIn(e, starts)
        self.assertTrue(all(e - s == n for so in c.seed.values() for s, e in so.intervals))


class TestCompare(unittest.TestCase):
    def test_never_clipped_and_noted(self):
        for name in ('rep0_cut1', 'rep0_lower_bound', 'rep0_column'):
            a, b = fresh(name), fresh(name)
            sa, sb = stripped(name), stripped(name)
            for mode in ('claims', 'walks', 'labels', 'prefix_subset'):
                with_ = a.compare(b, mode=mode).as_dict()
                without = sa.compare(sb, mode=mode).as_dict()
                self.assertIn(ops.COORDINATES_NOTE, with_['notes'])
                with_['notes'].remove(ops.COORDINATES_NOTE)
                self.assertEqual(without, with_, (name, mode))
                self.assertEqual(sa.compare_cost(sb, mode=mode), a.compare_cost(b, mode=mode))
                # the same charge: compare() reads no coordinate
                b1, b2 = LocalBudget(), LocalBudget()
                a.compare(b, mode=mode, budget=b1)
                sa.compare(sb, mode=mode, budget=b2)
                self.assertEqual(b2.usage()['work_units'], b1.usage()['work_units'])
            # one side with them is enough for the note
            self.assertIn(ops.COORDINATES_NOTE, a.compare(sb).notes)
            self.assertNotIn(ops.COORDINATES_NOTE, sa.compare(sb).notes)

    def test_the_docstring_states_the_label_limit_and_the_note(self):
        doc = ops.compare.__doc__
        self.assertIn('max_labels_per_node', doc)
        self.assertIn('coordinates are not compared', doc)


# ======================================================================= exports

class TestToJsonEqualsFull(unittest.TestCase):
    def test_deterministic_cells(self):
        for name in DETERMINISTIC:
            g = fresh(name)
            full = resp(name, 'full')['results'][0]
            self.assertEqual(export.normalize_result(full),
                             export.normalize_result(g.to_json()), name)
            self.assertEqual(g.to_json(), g.to_json(budget=LocalBudget()))

    def test_a_copy(self):
        g = fresh('rep0_cut1')
        j = g.to_json()
        j['coordinates']['k'] = 0
        self.assertEqual(31, g.seed_summary['coordinates']['k'])


class TestFasta(unittest.TestCase):
    def headers(self, text):
        return [x for x in text.split('\n') if x.startswith('>')]

    def test_header_field(self):
        g = fresh('rep0_unlimited_left')
        a = g.arms['left']
        c = g.coordinates
        text = g.to_fasta()
        for p, h in zip(derive.paths(a), self.headers(text)):
            pairs = []
            for e in derive.end_labels(a, p.leaf):
                run = a.runs[e.run]
                if run.from_label is None and run.from_bp == 0 and run.to_bp == p.length_bp:
                    for s, t in c.run('left', e.run).intervals:
                        pairs.append((e.label, s, t + g.seed.length_bp))
            pairs.sort()
            if not pairs:
                self.assertNotIn('coords=', h)
                continue
            want = ','.join('%s:%d-%d' % (g.labels[l].name, s + 1, t)
                            for l, s, t in pairs[:8])
            self.assertIn(' coords=' + want, h)
            if len(pairs) > 8:
                self.assertIn(' coords_more=%d' % (len(pairs) - 8), h)
        # without the seed: the walk's own bases
        for h in self.headers(g.to_fasta(with_seed=False)):
            for item in h.split(' coords=')[1].split(' ')[0].split(',') \
                    if ' coords=' in h else ():
                s, e = map(int, item.rsplit(':', 1)[1].split('-'))
                length = int(h.split('length_bp=')[1].split(' ')[0])
                self.assertEqual(length, e - s + 1)

    def test_cut_lists_and_column_labels(self):
        g = fresh('rep0_cut1')
        hs = self.headers(g.to_fasta())
        self.assertTrue(any(' coords_cut=1' in h for h in hs))
        # column labels are never stated (global positions)
        self.assertFalse(any('coords=' in h for h in self.headers(
            fresh('rep0_column').to_fasta())))

    def test_byte_identical_without_coordinates(self):
        for name in ('rep0_cut1', 'rep0_lower_bound', 'rep0_mixed'):
            g = fresh(name)
            plain = stripped(name).to_fasta()
            self.assertEqual(plain, g.to_fasta(coordinates=False))
            self.assertEqual(plain, '\n'.join(
                x.split(' coords')[0] if x.startswith('>') else x
                for x in g.to_fasta().split('\n')))
        self.assertEqual(fresh('d3_kmer').to_fasta(), fresh('d3_kmer').to_fasta(
            coordinates=None))
        with self.assertRaises(ValueError):
            fresh('d3_kmer').to_fasta(coordinates=True)
        with self.assertRaises(ValueError):
            fresh('rep0_cut1').to_fasta(coordinates='yes')
        # 0 equals False and must not be dispatched as None, which would write the coords= it
        # is to drop
        for bad in (0, 1):
            with self.assertRaises(ValueError, msg=repr(bad)):
                fresh('rep0_cut1').to_fasta(coordinates=bad)

    def test_views_and_budgets(self):
        g = fresh('rep0_cut1')
        v = ops.GraphletView(g, [g.labels[0].ref])
        self.assertIn('coords=', v.to_fasta())
        self.assertNotIn('coords=', v.to_fasta(coordinates=False))
        self.assertEqual(g.to_fasta(), g.to_fasta(budget=LocalBudget()))

    def test_no_coordinates_in_gfa(self):
        # GFA has no field for them in v1
        self.assertEqual(stripped('rep0_cut1').to_gfa(), fresh('rep0_cut1').to_gfa())


# ======================================================================= requests

class TestNextRequest(unittest.TestCase):
    def test_carries_the_setting(self):
        g = fresh('rep0_cut1')
        leaf = derive.leaves(g.arms['right'])[0]
        req = g.next_request('right', [ops.path_id(g.arms['right'], leaf)])
        out = req['strategy']['output']
        self.assertEqual((True, 1), (out['coordinates'], out['max_coordinate_occurrences']))

    def test_dropping_coordinates_drops_the_cap(self):
        # drop_coordinates through a continuation: the server refuses the cap alone (400)
        g = fresh('rep0_cut1')
        pid = ops.path_id(g.arms['right'], derive.leaves(g.arms['right'])[0])
        req = g.next_request('right', [pid], output={'coordinates': False})
        out = req['strategy']['output']
        self.assertFalse(out['coordinates'])
        self.assertNotIn('max_coordinate_occurrences', out)
        req = g.next_requests('right', [pid], output={'coordinates': False})[0]
        self.assertNotIn('max_coordinate_occurrences', req['strategy']['output'])

    def test_deepen_and_traverse_continue(self):
        g = fresh('rep0_cut1')
        pid = ops.path_id(g.arms['right'], derive.leaves(g.arms['right'])[0])
        sent = []

        class Client(TraverseClient):
            def traverse_raw(self, request):
                sent.append(request)
                return resp('rep0_cut1')
        c = Client('h', 1)
        c.deepen(g, 'right', [pid], output={'coordinates': False})
        self.assertNotIn('max_coordinate_occurrences', sent[-1]['strategy']['output'])
        with tempfile.TemporaryDirectory() as d:
            store = GraphletStore(d)
            h = store.put(resp('rep0_cut1'), {'seeds': []}, source='mini')
            tools = GraphletTools(store, {'mini': c})
            out = tools.traverse_continue(h, 'right', pid, overrides={
                'output': {'coordinates': False}}, execute=False)
            self.assertNotIn('max_coordinate_occurrences', out['request']['strategy']['output'])


def _caps(supported=True, trace=True, block=True):
    caps = {'supports_trace': trace, 'has_coordinates': supported, 'feature_level': 6}
    if block:
        caps['coordinates'] = {'supported': supported, 'knob': 'output.coordinates',
                               'cap_knob': 'output.max_coordinate_occurrences',
                               'max_occurrences_default': 16}
    return caps


class _Client(TraverseClient):
    """A client whose server states |caps| and records what it is sent."""

    def __init__(self, caps, answer='rep0_cut1', status=200):
        super().__init__('h', 1)
        self.caps = caps
        self.sent = []
        self.caps_reads = 0
        self.answer = answer
        self.status = status

    def capabilities(self, graph=None, graph_path=None):
        self.caps_reads += 1
        if isinstance(self.caps, Exception):
            raise self.caps
        return self.caps

    def traverse_raw(self, request):
        self.sent.append(copy.deepcopy(request))
        if self.status != 200:
            raise TraverseError(400, "strategy.output: unknown field 'coordinates'")
        return resp(self.answer)


TRACE = {'support': 'trace', 'bounds': {'max_extension_bp': 300}}


class TestClient(unittest.TestCase):
    def test_build_request(self):
        c = TraverseClient('h', 1)
        self.assertNotIn('coordinates', c.build_request(['A'], TRACE)['strategy']['output'])
        out = c.build_request(['A'], TRACE, coordinates=True,
                              max_coordinate_occurrences='unlimited')['strategy']['output']
        self.assertEqual((True, 'unlimited'),
                         (out['coordinates'], out['max_coordinate_occurrences']))
        st = {'output': {'coordinates': True, 'max_coordinate_occurrences': 3}}
        out = c.build_request(['A'], st, coordinates=False)['strategy']['output']
        self.assertNotIn('coordinates', out)
        self.assertNotIn('max_coordinate_occurrences', out)
        # the cap is never sent alone (the server's 400: the cap needs coordinates: true)
        out = c.build_request(['A'], TRACE, max_coordinate_occurrences=3)['strategy']['output']
        self.assertNotIn('max_coordinate_occurrences', out)
        self.assertEqual({'output': {'coordinates': True, 'max_coordinate_occurrences': 3}},
                         st)                          # pure
        with self.assertRaises(ValueError):
            c.build_request(['A'], TRACE, coordinates='auto')
        # 0 and 1 equal False and True; dispatched by identity they would be taken for None,
        # 0 keeping coordinates and 1 dropping an explicit cap -- refused by type
        for bad in (0, 1, 0.0, 1.0):
            with self.assertRaises(ValueError, msg=repr(bad)):
                c.build_request(['A'], st, coordinates=bad)
            with self.assertRaises(ValueError, msg=repr(bad)):
                c.build_request(['A'], TRACE, coordinates=bad, max_coordinate_occurrences=5)
        # a cap the server would refuse (400, after it registered an attempt_id) is refused
        # before anything is sent
        for bad in (0, -1, True, 1.5, 'all', (1 << 64) - 1):
            with self.assertRaises(ValueError, msg=repr(bad)):
                c.build_request(['A'], TRACE, coordinates=True, max_coordinate_occurrences=bad)
        out = c.build_request(['A'], TRACE, coordinates=True,
                              max_coordinate_occurrences=(1 << 64) - 2)['strategy']['output']
        self.assertEqual((1 << 64) - 2, out['max_coordinate_occurrences'])

    def test_supports_coordinates(self):
        self.assertTrue(_Client(_caps()).supports_coordinates())
        self.assertFalse(_Client(_caps(supported=False)).supports_coordinates())
        self.assertFalse(_Client(_caps(trace=False)).supports_coordinates())
        self.assertFalse(_Client(_caps(block=False)).supports_coordinates())   # level 5
        self.assertFalse(_Client(TraverseError(404, 'no route')).supports_coordinates())
        c = _Client(_caps())
        c.supports_coordinates()
        c.supports_coordinates()
        self.assertEqual(1, c.caps_reads)             # kept per (graph, graph_path)
        c.supports_coordinates(refresh=True)
        self.assertEqual(2, c.caps_reads)

    def test_auto(self):
        c = _Client(_caps())
        r = c.traverse(['A'], TRACE)
        self.assertTrue(c.sent[-1]['strategy']['output']['coordinates'])
        self.assertEqual({'requested': True, 'reason': 'supported'}, r.coordinates_auto)
        self.assertIsNotNone(r.graphlets[0].coordinates)
        # not for kmer, annotate, a pinned setting or an index without them
        for st in ({'bounds': {}}, dict(TRACE, labels={'mode': 'annotate'}),
                   dict(TRACE, output={'coordinates': False})):
            c.traverse(['A'], st)
            # (a pinned false is sent as given: the server reads it as absent)
            self.assertFalse(c.sent[-1]['strategy']['output'].get('coordinates'))
        c.traverse(['A'], dict(TRACE, output={'coordinates': True}), coordinates=None)
        self.assertTrue(c.sent[-1]['strategy']['output']['coordinates'])
        u = _Client(_caps(supported=False))
        r = u.traverse(['A'], TRACE)
        self.assertNotIn('coordinates', u.sent[-1]['strategy']['output'])
        self.assertEqual('unsupported', r.coordinates_auto['reason'])
        self.assertEqual([], r.notes)

    def test_auto_under_a_memory_budget(self):
        # the branch where the depth gate passes (AUTO_COORDINATES_UNDER_MEMORY_BUDGET true):
        # asked under a memory budget too, and stated. At level 6 the gate fails for wide
        # columns (DESIGN §26.5), so the shipped default is false (the test below)
        self.assertFalse(AUTO_COORDINATES_UNDER_MEMORY_BUDGET)
        import metagraph.traverse.client as client_mod
        client_mod.AUTO_COORDINATES_UNDER_MEMORY_BUDGET = True
        self.addCleanup(setattr, client_mod, 'AUTO_COORDINATES_UNDER_MEMORY_BUDGET', False)
        c = _Client(_caps())
        st = dict(TRACE, bounds={'max_extension_bp': 300, 'max_memory_mb': 4})
        r = c.traverse(['A'], st)
        self.assertTrue(c.sent[-1]['strategy']['output']['coordinates'])
        self.assertEqual({'requested': True, 'reason': 'supported', 'memory_budget': True},
                         r.coordinates_auto)
        self.assertEqual([], r.notes)          # the snapshot did not stop on memory
        self.assertEqual((None, None), auto_coordinates(c, dict(st, output={
            'coordinates': True})))
        # never on an index that cannot report them (the null form only costs depth)
        u = _Client(_caps(supported=False))
        u.traverse(['A'], st)
        self.assertNotIn('coordinates', u.sent[-1]['strategy']['output'])

    def test_a_memory_stop_says_how_to_drop_them(self):
        # d3_memory_stop: a memory stop that offers drop_coordinates (coordinates asked for
        # automatically under a memory budget: the gate-passes branch)
        import metagraph.traverse.client as client_mod
        client_mod.AUTO_COORDINATES_UNDER_MEMORY_BUDGET = True
        self.addCleanup(setattr, client_mod, 'AUTO_COORDINATES_UNDER_MEMORY_BUDGET', False)
        c = _Client(_caps(), answer='d3_memory_stop')
        st = dict(TRACE, bounds={'max_extension_bp': 300, 'max_memory_mb': 4})
        r = c.traverse(['A'], st)
        stop = r.graphlets[0].resource_stop
        self.assertIn('drop_coordinates', stop.actions)
        self.assertEqual(1, len(r.notes))
        self.assertIn('coordinates=False', r.notes[0])
        self.assertIn('1 seed(s)', r.notes[0])
        # asked for explicitly, or without a memory budget: the caller chose, no note
        self.assertEqual([], c.traverse(['A'], st, coordinates=True).notes)
        self.assertEqual([], c.traverse(['A'], TRACE).notes)
        # a failed seed's stop counts too
        results = [{'error': 'x', 'resource_stop': {'actions': ['drop_coordinates']}},
                   {'resource_stop': {'actions': ['use_graphlet']}}, {}]
        stated = {'requested': True, 'reason': 'supported', 'memory_budget': True}
        self.assertIn('1 seed(s)', drop_coordinates_note(results, stated))
        self.assertIsNone(drop_coordinates_note(results, {'requested': True,
                                                          'reason': 'supported'}))

    def test_auto_under_a_memory_budget_off_when_the_gate_fails(self):
        # the rule's other branch (AUTO_COORDINATES_UNDER_MEMORY_BUDGET false): not asked,
        # and said so with how to ask
        import metagraph.traverse.client as client_mod
        saved = client_mod.AUTO_COORDINATES_UNDER_MEMORY_BUDGET
        client_mod.AUTO_COORDINATES_UNDER_MEMORY_BUDGET = False
        try:
            c = _Client(_caps())
            st = dict(TRACE, bounds={'max_extension_bp': 300, 'max_memory_mb': 4})
            r = c.traverse(['A'], st)
            self.assertNotIn('coordinates', c.sent[-1]['strategy']['output'])
            self.assertEqual('memory_budget', r.coordinates_auto['reason'])
            self.assertIn('coordinates=True', r.notes[-1])
            c.traverse(['A'], st, coordinates=True)
            self.assertTrue(c.sent[-1]['strategy']['output']['coordinates'])
        finally:
            client_mod.AUTO_COORDINATES_UNDER_MEMORY_BUDGET = saved

    def test_no_silent_retry(self):
        c = _Client(_caps(), status=400)
        with self.assertRaises(TraverseError):
            c.traverse(['A'], TRACE)
        self.assertEqual(1, len(c.sent))


# ======================================================================= the tools

class _ToolClient(_Client):
    def build_request(self, seeds, strategy=None, detail='graphlet', **kw):
        return super().build_request(seeds, strategy, detail=detail, **kw)


class TestTools(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.store = GraphletStore(self.tmp.name)
        self.client = _ToolClient(_caps())
        self.tools = GraphletTools(self.store, {'mini': self.client}, max_bytes=16384)

    def tearDown(self):
        self.tmp.cleanup()

    def put(self, name):
        return self.store.put(resp(name), {'seeds': []}, source='mini')

    def test_claims_rows(self):
        for name in ('rep0_cut1', 'rep0_lower_bound', 'rep0_mixed'):
            h = self.put(name)
            g = fresh(name)
            for side in g.arms:
                out = self.tools.graphlet_claims(h, side, route_consistent=False,
                                                 max_bytes=1 << 20)
                self.assertEqual('record' if name != 'rep0_mixed' else 'mixed',
                                 out['evidence']['coordinates']['kind'])
                claims = g.claims(side, strict=False)
                self.assertEqual(len(claims), len(out['rows']))
                for row, cl in zip(out['rows'], claims):
                    rc = cl.coordinates
                    self.assertEqual([list(x) for x in rc.intervals[:4]], row['coordinates'])
                    self.assertEqual(rc.total, row['coordinates_total'])
                    self.assertEqual(rc.total > min(4, len(rc.intervals)),
                                     row['coordinates_truncated'])
                    self.assertEqual(rc.lower_bound, row.get('coordinates_lower_bound', False))
                    if name == 'rep0_mixed':
                        self.assertEqual(coords.label_kind(g, cl.label.id),
                                         row['coordinates_kind'])
                    else:
                        self.assertNotIn('coordinates_kind', row)
        # rows of a retrieval without coordinates are as they were
        h = self.put('d3_kmer')
        out = self.tools.graphlet_claims(h, 'right')
        self.assertFalse(any('coordinates' in r for r in out['rows']))
        self.assertNotIn('coordinates', out['evidence'])

    def test_walks_and_labels_flag(self):
        h = self.put('rep0_cut1')
        out = self.tools.graphlet_walks(h, 'right', n=3)
        self.assertFalse(any('coordinates' in r for r in out['rows']))
        out = self.tools.graphlet_walks(h, 'right', n=3, coordinates=True)
        g = fresh('rep0_cut1')
        for row in out['rows']:
            w = g.walks('right', top=None, by='id')[row['walk']]
            for x in row['coordinates']:
                rc = w.coordinates[x['ref']]
                self.assertEqual([list(i) for i in rc.intervals[:4]], x['coordinates'])
        lab = g.labels[0].ref
        out = self.tools.graphlet_labels(h, name={'ref': lab}, coordinates=True)
        self.assertTrue(all('coordinates_total' in r for r in out['rows']))
        out = self.tools.graphlet_labels(h, coordinates=True)
        self.assertEqual('bad_argument', out['error'])
        h2 = self.put('d3_kmer')
        self.assertEqual('no_coordinates',
                         self.tools.graphlet_walks(h2, 'right', coordinates=True)['error'])

    def test_fetch_auto(self):
        out = self.tools.traverse_fetch(seed={'sequence': 'ACGT'}, strategy=TRACE)
        self.assertNotIn('error', out, out)
        self.assertTrue(self.client.sent[-1]['strategy']['output']['coordinates'])
        self.assertEqual('record', out['evidence']['coordinates']['kind'])
        st = dict(TRACE, bounds={'max_extension_bp': 300, 'max_memory_mb': 4})
        # the branch where the depth gate passes: under a memory budget too, and no note
        # without a stop (the shipped default is false at level 6, DESIGN §26.5)
        import metagraph.traverse.client as client_mod
        client_mod.AUTO_COORDINATES_UNDER_MEMORY_BUDGET = True
        self.addCleanup(setattr, client_mod, 'AUTO_COORDINATES_UNDER_MEMORY_BUDGET', False)
        out = self.tools.traverse_fetch(seed={'sequence': 'ACGT'}, strategy=st)
        self.assertTrue(self.client.sent[-1]['strategy']['output']['coordinates'])
        self.assertNotIn('coordinates_note', out)
        # a memory stop offering drop_coordinates: the note says how to walk further
        self.client.answer = 'd3_memory_stop'
        out = self.tools.traverse_fetch(seed={'sequence': 'ACGT'}, strategy=st)
        self.assertIn('coordinates=False', out['coordinates_note'])
        out = self.tools.traverse_fetch(seed={'sequence': 'ACGT'}, strategy=st,
                                        coordinates=True)
        self.assertNotIn('coordinates_note', out)
        self.client.answer = 'rep0_cut1'
        import metagraph.traverse.client as client_mod
        client_mod.AUTO_COORDINATES_UNDER_MEMORY_BUDGET = False
        try:
            out = self.tools.traverse_fetch(seed={'sequence': 'ACGT'}, strategy=st)
        finally:
            client_mod.AUTO_COORDINATES_UNDER_MEMORY_BUDGET = True
        self.assertNotIn('coordinates', self.client.sent[-1]['strategy']['output'])
        self.assertIn('coordinates=True', out['coordinates_note'])
        out = self.tools.traverse_fetch(seed={'sequence': 'ACGT'}, strategy=st,
                                        coordinates=True)
        self.assertTrue(self.client.sent[-1]['strategy']['output']['coordinates'])
        out = self.tools.traverse_fetch(seed={'sequence': 'ACGT'}, strategy=dict(
            TRACE, output={'coordinates': True, 'max_coordinate_occurrences': 3}),
            coordinates=False)
        self.assertNotIn('coordinates', self.client.sent[-1]['strategy']['output'])
        self.assertNotIn('max_coordinate_occurrences', self.client.sent[-1]['strategy'][
            'output'])
        self.assertEqual('bad_argument', self.tools.traverse_fetch(
            seed={'sequence': 'A'}, coordinates='yes')['error'])
        # strict, as every other flag (_bool): 0 and 1 are no silent no-op
        sent = len(self.client.sent)
        for bad in (0, 1):
            self.assertEqual('bad_argument', self.tools.traverse_fetch(
                seed={'sequence': 'A'}, strategy=TRACE, coordinates=bad)['error'])
        self.assertEqual(sent, len(self.client.sent))

    def test_budgeted_pages_equal_unbudgeted(self):
        tools = GraphletTools(self.store, {'mini': self.client}, local_limits=ToolLimits())
        h = self.put('rep0_lower_bound')
        for side in ('left', 'right'):
            a = self.tools.graphlet_claims(h, side, route_consistent=False, max_bytes=1 << 20)
            b = tools.graphlet_claims(h, side, route_consistent=False, max_bytes=1 << 20)
            self.assertEqual(a['rows'], b['rows'])
            a = self.tools.graphlet_walks(h, side, coordinates=True, max_bytes=1 << 20)
            b = tools.graphlet_walks(h, side, coordinates=True, max_bytes=1 << 20)
            self.assertEqual(a['rows'], b['rows'])


# ======================================================================= local limits

class TestStageL(unittest.TestCase):
    OPS = {
        'claims': lambda g, **k: g.claims(**k),
        'claims_cut': lambda g, **k: g.claims(at_most_bp=50, **k),
        'walks': lambda g, **k: g.walks('right', **k),
        'walks_id': lambda g, **k: g.walks('right', by='id', **k),
        'label_walks': lambda g, **k: g.label_walks(g.labels[0], **k),
        'to_fasta': lambda g, **k: g.to_fasta(**k),
        'to_json': lambda g, **k: g.to_json(**k),
        'summary': lambda g, **k: g.summary(**k),
    }

    def test_unlimited_budgets_answer_as_without(self):
        for name in BLOCKS:
            for op, fn in self.OPS.items():
                if op.startswith('walks') and 'right' not in fresh(name).arms:
                    continue
                a = fn(fresh(name))
                b = fn(fresh(name), budget=LocalBudget())
                self.assertEqual(a, b, (name, op))

    def test_coordinates_are_charged(self):
        g, s = fresh('rep0_lower_bound'), stripped('rep0_lower_bound')
        for op in ('claims', 'walks', 'label_walks', 'to_fasta'):
            b1, b2 = LocalBudget(), LocalBudget()
            self.OPS[op](g, budget=b1)
            self.OPS[op](s, budget=b2)
            self.assertGreater(b1.usage()['work_units'], b2.usage()['work_units'], op)

    def test_the_account_bounds_the_traced_peak(self):
        for name in ('rep0_unlimited_left', 'rep0_lower_bound', 'rep0_column'):
            for op, fn in self.OPS.items():
                if op.startswith('walks') and 'right' not in fresh(name).arms:
                    continue
                b = LocalBudget()
                fn(fresh(name), budget=b)
                g = fresh(name)
                peak, _ = _peak(lambda: fn(g))
                self.assertGreaterEqual(b.peak_bytes, peak, (name, op))
            j = resp(name)
            b = LocalBudget()
            MT.from_response(j['results'][0], j, budget=b)
            peak, _ = _peak(lambda: MT.from_response(j['results'][0], j))
            self.assertGreaterEqual(b.peak_bytes, peak, (name, 'from_response'))

    def test_the_account_bounds_the_traced_peak_of_a_wide_block(self):
        # on the mini_refseq cells a few occurrences hide behind the fixed overheads; over a
        # block of 15,000, a set of the seed starts per run built by the check of the seed
        # continuity, and a memo entry per copied list kept by to_json's deepcopy, would go
        # uncharged: 1.22x and 1.34x the account on the wide fixture. Two shapes: every label
        # at 300 occurrences (14,100 in all; to_json at 1.35x uncharged), and one seed label
        # at 5,000 with its one run (the wide fixture's shape, one big seed list checked
        # against its run: validate() at 1.14x the index price, from_response() 1.05x the
        # account, uncharged)
        cases = [('rep0_column', widened('rep0_column', 300), self.OPS),
                 ('rep0_unlimited_left/label 3',
                  widened('rep0_unlimited_left', 5000, labels={3}),
                  {k: self.OPS[k] for k in ('to_json', 'claims_cut')})]
        for name, j, ops_ in cases:
            r = j['results'][0]
            entries, occ = coords._block_counts(r['coordinates'])
            self.assertGreaterEqual(occ, 10000)
            b = LocalBudget()
            g = MT.from_response(r, j, budget=b)
            peak, _ = _peak(lambda: MT.from_response(r, j))
            self.assertGreaterEqual(b.peak_bytes, peak, (name, 'from_response'))
            # the validation alone within the price of the index it builds
            peak, _ = _peak(lambda: coords.validate(g))
            self.assertGreaterEqual(coords.index_price(entries, occ)[1], peak,
                                    (name, 'validate'))
            for op, fn in ops_.items():
                b = LocalBudget()
                fn(MT.from_response(r, j), budget=b)
                g = MT.from_response(r, j)
                peak, _ = _peak(lambda: fn(g))
                self.assertGreaterEqual(b.peak_bytes, peak, (name, op))

    def test_to_json_copies_the_block(self):
        # a structural copy (no deepcopy memo): equal, and sharing no container with the
        # graphlet's summary
        g = fresh('rep0_lower_bound')
        out = g.to_json()['coordinates']
        raw = g.seed_summary['coordinates']
        self.assertEqual(raw, out)
        seen = set()

        def ids(x):
            if isinstance(x, (dict, list)):
                seen.add(id(x))
                for v in (x.values() if isinstance(x, dict) else x):
                    ids(v)
        ids(raw)

        def shared(x):
            if isinstance(x, (dict, list)):
                return id(x) in seen or any(shared(v) for v in (
                    x.values() if isinstance(x, dict) else x))
            return False
        self.assertFalse(shared(out))

    def test_a_cold_index_is_charged_as_a_warm_one(self):
        # the index is charged once per call at its cold price, whether the cache holds it or
        # not: a cold build charging it unkeyed, and the other arm's use charging it again,
        # would count it twice
        for name in BLOCKS:
            for op in ('claims', 'claims_cut', 'walks', 'label_walks', 'to_fasta'):
                if op == 'walks' and 'right' not in fresh(name).arms:
                    continue
                got = []
                for cold in (False, True):
                    g = fresh(name)
                    if cold:
                        del g.cache[coords.CACHE_KEY]
                    b = LocalBudget()
                    self.OPS[op](g, budget=b)
                    got.append((b.usage()['work_units'], b.peak_bytes))
                    self.assertIn(coords.CACHE_KEY, g.cache)
                self.assertEqual(got[0], got[1], (name, op))

    def test_a_store_re_parse_attaches_the_index(self):
        # a parsed entry parsed again (its resident model expired or was evicted) is the
        # model from_response() made: the index cached, the same size and signature, and
        # a call on it charged as on the resident one
        with tempfile.TemporaryDirectory() as d:
            store = GraphletStore(d)
            h = store.put(resp('rep0_cut1'), {'seeds': []})
            g1 = store.graphlet(h)
            size, sig = g1.memory_bytes(), ops.cache_signature(g1)
            b1 = LocalBudget()
            g1.label_walks(g1.labels[0], budget=b1)
            store._drop_ram(h)
            g2 = store.graphlet(h)
            self.assertIsNot(g1, g2)
            self.assertIn(coords.CACHE_KEY, g2.cache)
            self.assertEqual((size, sig), (g2.memory_bytes(), ops.cache_signature(g2)))
            b2 = LocalBudget()
            g2.label_walks(g2.labels[0], budget=b2)
            self.assertEqual(b1.usage(), b2.usage())
            # charged like from_response(): the parse budget pays the index too
            pb = LocalBudget()
            store._drop_ram(h)
            store.graphlet(h, parse_budget=pb)
            j = resp('rep0_cut1')
            body = MT.from_response(j['results'][0], j).dump(envelope=False)
            b0 = LocalBudget()
            MT.parse(body, budget=b0)
            w, _ = coords.index_price(*coords._block_counts(j['results'][0]['coordinates']))
            self.assertEqual(b0.usage()['work_units'] + w, pb.usage()['work_units'])

    def test_a_validation_stop_keeps_the_body_unparsed(self):
        # the block is charged before it is built: a parse budget that admits the body but
        # not the index stops in it, and the store keeps the entry unparsed
        j = resp('rep0_unlimited_left')
        r = j['results'][0]
        b = LocalBudget()
        MT.from_response(r, j, budget=b)
        whole = b.usage()['work_units']
        entries, occ = coords._block_counts(r['coordinates'])
        w, _ = coords.index_price(entries, occ)
        with self.assertRaises(LocalBudgetExceeded):
            MT.from_response(r, j, budget=LocalBudget(work_units=whole - w))
        with tempfile.TemporaryDirectory() as d:
            store = GraphletStore(d, parse_limits=LocalLimits(work_units=whole - w))
            with self.assertRaises(LocalBudgetExceeded) as cm:
                store.put(j, {'seeds': []})
            self.assertFalse(store.get(cm.exception.handle).parsed)
            store.parse_limits = None
            self.assertIsNotNone(store.graphlet(cm.exception.handle).coordinates)

    def test_claims_stop_and_resume_with_coordinates(self):
        g = fresh('rep0_lower_bound')
        want = g.claims('right')
        b = LocalBudget()
        g.claims('right', budget=b)
        got, resume = [], None
        while True:
            try:
                got.extend(g.claims('right', budget=LocalBudget(
                    work_units=b.usage()['work_units'] * 3 // 4), resume=resume))
                break
            except LocalBudgetExceeded as e:
                self.assertTrue(e.partial.rows or got)
                got.extend(e.partial.rows)
                resume = e.partial.resume
        self.assertEqual(want, got)


# ======================================================================= partial derivation

class TestPartialDerivation(unittest.TestCase):
    def test_check_rules_accepts_the_qualified_walked_result(self):
        for name in ('d3_kmer', 'd3_trace_coordinates'):
            g = fresh(name)
            self.assertEqual(('partial', 'qualified'),
                             (g.outcome.walks, g.outcome.label_evidence))
            self.assertEqual([], g.check_rules(), name)
            lim = derive.partial_derivation(g)
            self.assertEqual('derivation', lim.kind)
            self.assertTrue(derive.qualifies(g, lim))
            # the evidence of such a walk is exact on no arm, and a comparison says why
            self.assertFalse(ops.evidence_block(g)['exact'])
            c = g.compare(fresh(name))
            # qualified by the derivation -- and unknown, weaker still: the walk stopped at
            # the seed, so no depth is certified on either side
            self.assertEqual('unknown', c.comparable)
            self.assertTrue(any('derivation (bounds.time_budget_ms)' in n for n in c.notes))

    def test_a_failed_derivation_qualifies_nothing(self):
        g = fresh('d3_kmer')
        g.outcome.walks = 'failed'
        g.outcome.label_evidence = 'complete'
        self.assertIsNone(derive.partial_derivation(g))
        self.assertFalse(any(f['field'] == 'outcome.label_evidence'
                             for f in derive._outcome_findings(g)))
        # and a walked one stating complete beside its derivation is still flagged
        g = fresh('d3_kmer')
        g.outcome.label_evidence = 'complete'
        self.assertTrue(any(f['field'] == 'outcome.label_evidence'
                            for f in g.check_rules()))

    def test_the_summary_names_it(self):
        s = fresh('d3_kmer').summary()
        self.assertEqual({'partial': True, 'kmers_read': 1, 'of': 370},
                         s['seed']['derivation'])
        self.assertNotIn('derivation', fresh('rep0_cut1').summary()['seed'])
        self.assertEqual({'reason': 'partial derivation'},
                         fresh('d3_trace_coordinates').summary()['coordinates'])


# ======================================================================= the majority-parent rule

class TestDisplayedParentAtMerges(unittest.TestCase):
    """From feature level 6 the server makes a merge's first parent -- the one the
    displayed walk follows -- the head carrying the most labels (ties in arrival order);
    the library displays the first parent as stored, check_rules states a level-6
    retrieval that breaks the rule (derive.carried_labels: the server's counts), and
    compare() qualifies a comparison whose sides may display a merge through different
    parents (derive.displayed_parent_rule, derive.majority_first)."""

    def merges(self, g):
        return [(side, s) for side, a in g.arms.items() for s in a.segments
                if len(s.parents) > 1]

    def test_the_counts_on_the_merge_fixtures(self):
        # (first, second) parent's labels at each merge: the majority first, ties kept
        want = {'merge': {3: [3, 2], 6: [3, 2]},          # at 62 the G allele, arrived 2nd
                'merge_ties': {3: [2, 2], 6: [2, 2]},     # ties: arrival order
                'annotate': {3: [3, 2], 6: [3, 2]}}       # true counts of cut lists
        for name, merges in want.items():
            g = T.graphlet(name)
            self.assertEqual(6, derive.feature_level(g))
            got = {s.id: [derive.carried_labels(g, g.arms[side], s, i)
                          for i in range(len(s.parents))] for side, s in self.merges(g)}
            self.assertEqual(merges, got, name)
            self.assertEqual([], g.check_rules(), name)
        # merge: both.fa came through both parents at 62; the minority's lineage was closed
        # there (an 'm' run) and is counted for its parent
        g = T.graphlet('merge')
        a = g.arms['right']
        self.assertEqual((5, 4), a.segments[6].parents)
        self.assertEqual([(3, 4, 62)], [(r.label, r.segment, r.to_bp) for r in a.runs
                                        if r.merged and r.segment == 4])

    def swapped(self, name, old, new, level):
        j = T.doc_json(name, 'graphlet')
        r = j['results'][0]
        self.assertIn(old, r['graphlet'])
        r['graphlet'] = r['graphlet'].replace(old, new)
        g = MT.from_response(r, j)
        g.envelope['capabilities']['feature_level'] = level
        return [f for f in g.check_rules() if f.get('field') == 'parents']

    def test_a_level_6_document_breaking_it_is_flagged(self):
        # annotate: the merge at 62 in arrival order (C allele, 2 labels, first) -- what a
        # level-5 server wrote, a violation at level 6 only
        old, new = 'G 5,4 62 8 ! 4', 'G 4,5 62 8 ! 4'
        (f,) = self.swapped('annotate', old, new, 6)
        self.assertEqual(('invariant', 6), (f['kind'], f['segment']))
        self.assertIn('parent 4 carries 2 labels, parent 5 carries 3', f['msg'])
        self.assertEqual([], self.swapped('annotate', old, new, 5))
        # constrain: the same order with the partition swapped along (its route stamps then
        # disagree too, which check_rules states separately)
        old, new = 'G 5,4 62 38 * * * 1-3|0', 'G 4,5 62 38 * * * 0|1-3'
        self.assertEqual(1, len(self.swapped('merge', old, new, 6)))
        self.assertEqual([], self.swapped('merge', old, new, 5))

    def test_cut_lists_in_constrain_mode_are_not_judged(self):
        # with recorded lists cut the partitions are lower bounds: no verdict
        g = T.graphlet('merge')
        a = g.arms['right']
        s = a.segments[6]
        s.partition[0], s.partition[1] = s.partition[1], s.partition[0]
        a.cache.clear()
        self.assertTrue([f for f in g.check_rules() if f.get('field') == 'parents'])
        a.labels_per_node.nodes_truncated = 1
        self.assertFalse([f for f in g.check_rules() if f.get('field') == 'parents'])

    def level5(self, name):
        """The level-5 document of the same locus (data/traverse/r21: the CLI fixture as a
        level-5 walker writes it, its 62 merge in arrival order)."""
        with open(os.path.join(T.DATA, 'r21', name + '_level5.graphlet.json'),
                  encoding='utf-8') as f:
            j = json.load(f)
        return MT.from_response(j['results'][0], j)

    @staticmethod
    def flagged(c):
        return [n for n in c.notes if n.startswith('the displayed parent at a merge')]

    def test_the_rules_and_the_level_5_documents(self):
        for name in ('merge', 'annotate'):
            old, new = self.level5(name), T.graphlet(name)
            self.assertEqual(('arrival', 'majority'), (derive.displayed_parent_rule(old),
                                                       derive.displayed_parent_rule(new)))
            # the level-5 document breaks the level-6 rule at its 62 merge; the level-6 one
            # keeps it
            self.assertIs(False, derive.majority_first(old, 'right'))
            self.assertIs(True, derive.majority_first(new, 'right'))
            self.assertIs(True, derive.majority_first(old, 'right', depth=62))
            body = MT.parse(new.dump(envelope=False))
            self.assertIsNone(derive.displayed_parent_rule(body))
            # capabilities that state no level: a server from before the levels were stated
            del old.envelope['capabilities']['feature_level']
            self.assertEqual('arrival', derive.displayed_parent_rule(old))
            old.envelope.pop('capabilities')
            self.assertIsNone(derive.displayed_parent_rule(old))
        g = T.graphlet('merge')
        g.arms['right'].labels_per_node.nodes_truncated = 1
        self.assertIsNone(derive.majority_first(g, 'right'))   # cut lists: not judged

    def test_compare_across_the_rule_change_is_qualified(self):
        # the same locus walked by a level-5 and a level-6 server: the walks through the 62
        # merge are spelled through different parents, so claims and walks differ where
        # the retrieval does not -- never a definite difference (a level-6 body without its
        # envelope against the level-5 retrieval must not be comparable True with definite
        # differences)
        for name in ('merge', 'annotate'):
            def pair(new_env, old_env):
                new, old = T.graphlet(name), self.level5(name)
                if not new_env:
                    new = MT.parse(new.dump(envelope=False))
                if not old_env:
                    old = MT.parse(old.dump(envelope=False))
                for g in (new, old):
                    g.index_fp = 'f' * 16       # one index on both sides (comparable)
                return new, old
            for new_env, old_env in ((True, True), (False, True), (True, False),
                                     (False, False)):
                new, old = pair(new_env, old_env)
                for a, b in ((new, old), (old, new)):
                    for mode in ('claims', 'walks', 'prefix_subset'):
                        c = a.compare(b, mode=mode)
                        what = (name, new_env, old_env, a is new, mode)
                        self.assertEqual('qualified', c.comparable, what)
                        self.assertTrue(self.flagged(c), what)
                        self.assertIsNone(c.equal, what)
                        cb = a.compare(b, mode=mode, budget=LocalBudget())
                        self.assertEqual(c.as_dict(), cb.as_dict(), what)
                    # the display is the only difference: the labels agree
                    c = a.compare(b, mode='labels')
                    self.assertFalse(self.flagged(c))
            # the counts it reads are charged, and compare_cost() prices them: its exact
            # estimate (a side without bases, nothing keyed) stays exact
            new, old = pair(True, True)
            w = ops._displayed_parent_price(new, old, ['right'], 1000)
            self.assertGreater(w, 0)
            for s in old.arms['right'].segments:
                s.walk = None
            for mode in ('claims', 'walks', 'prefix_subset'):
                est = ops.compare_cost(new, old, mode=mode)
                b = LocalBudget()
                new.compare(old, mode=mode, budget=b)
                self.assertTrue(est['exact'])
                self.assertEqual(b.used_work, est['work_units']['at_least'], (name, mode))
            # the comparison would otherwise have stated definite differences (annotate's
            # cut lists qualify it anyway)
            new, old = pair(True, True)
            old.envelope['capabilities']['feature_level'] = 6
            c = new.compare(old, mode='walks')
            self.assertFalse(self.flagged(c))
            self.assertTrue(c.only_in_a and c.only_in_b)
            if name == 'merge':
                self.assertEqual((True, False), (c.comparable, c.equal))

    def test_no_rule_difference_is_not_qualified(self):
        # (annotate's comparisons are qualified by its cut lists: only the note is checked)
        def same(c, name):
            self.assertFalse(self.flagged(c), name)
            if name == 'merge':
                self.assertEqual((True, True), (c.comparable, c.equal), name)
            else:
                self.assertFalse(c.only_in_a or c.only_in_b or c.differ, name)

        for name in ('merge', 'annotate'):
            new = T.graphlet(name)
            body = MT.parse(new.dump(envelope=False))
            for g in (new, body):
                g.index_fp = 'f' * 16
            # a level-6 retrieval and its body: both show the majority parent
            for a, b in ((new, body), (body, new), (body, body)):
                same(a.compare(b, mode='walks'), name)
            # a level-5 retrieval whose merges show their majority parent displays what
            # level 6 would: the levels differ, the displays do not
            five = T.graphlet(name)
            five.envelope['capabilities']['feature_level'] = 5
            five.index_fp = 'f' * 16
            same(new.compare(five, mode='walks'), name)
            # two retrievals of one older rule are not qualified, whatever they display
            o1, o2 = self.level5(name), self.level5(name)
            del o2.envelope['capabilities']['feature_level']      # a server before levels
            for g in (o1, o2):
                g.index_fp = 'f' * 16
            same(o1.compare(o2, mode='walks'), name)
        # a trie without merges is not flagged
        f1, f2 = T.graphlet('fork'), T.graphlet('fork')
        f2.envelope = None
        f2.seed_summary = None
        self.assertFalse(self.flagged(f1.compare(f2, mode='walks')))
        # a merge at or beyond the depth is not compared through
        new, old = T.graphlet('merge'), self.level5('merge')
        self.assertIsNone(ops._displayed_parent_rules(new, old, ['right'], 62))
        self.assertIsNotNone(ops._displayed_parent_rules(new, old, ['right'], 63))

    def test_the_displayed_walk_follows_the_first_parent(self):
        for name in ('merge', 'merge_ties', 'annotate'):
            g = T.graphlet(name)
            for side, s in self.merges(g):
                for w in g.walks(side):
                    if s.id in w.segments:
                        i = w.segments.index(s.id)
                        self.assertEqual(s.parents[0], w.segments[i - 1])
                        # its bases are the first parent's
                        p = g.arms[side].segments[s.parents[0]]
                        if w.sequence is not None and p.walk is not None and side == 'right':
                            self.assertEqual(p.walk, w.sequence[p.from_bp:p.end_bp])


# ======================================================================= ambient budgets

class TestAmbientBudgetUnderToolsWithoutLimits(unittest.TestCase):
    """Under a caller's ambient budget, tools without local limits must not store an entry
    and lose its handle when the summary after it stops, and a fetch whose parse stops must
    not raise AttributeError (no call to state the stop on)."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.store = GraphletStore(self.tmp.name)
        self.client = _ToolClient(_caps(supported=False), answer='rep0_cut1')
        self.tools = GraphletTools(self.store, {'mini': self.client}, max_bytes=16384)

    def tearDown(self):
        self.tmp.cleanup()

    def parse_units(self):
        j = resp('rep0_cut1')
        b = LocalBudget()
        MT.from_response(j['results'][0], j, budget=b)
        return b.usage()['work_units']

    def test_a_stopped_summary_keeps_the_handle(self):
        units = self.parse_units()
        with MT.local_budget(LocalLimits(work_units=units + 50)):
            out = self.tools.traverse_fetch(seed={'sequence': 'ACGT'}, strategy=TRACE)
        self.assertIn('handle', out, out)
        self.assertIn('summary_stop', out)
        self.assertNotIn('summary', out)
        self.assertEqual([out['handle']], [e['handle'] for e in self.store.list()])

    def test_a_stopped_parse_answers_without_raising(self):
        with MT.local_budget(LocalLimits(work_units=10)):
            out = self.tools.traverse_fetch(seed={'sequence': 'ACGT'}, strategy=TRACE)
        self.assertFalse(out['parsed'])
        self.assertEqual('work', out['stop']['resource'])
        self.assertEqual([out['handle']], [e['handle'] for e in self.store.list()])

    def test_graphlet_load(self):
        path = os.path.join(self.tmp.name, 'x.mgt')
        export_dir = self.tmp.name
        fresh('rep0_cut1').save(path)
        tools = GraphletTools(self.store, {}, export_dir=export_dir)
        with open(path, encoding='utf-8') as f:
            text = f.read()
        b = LocalBudget()
        g = MT.parse(text, budget=b)
        b2 = LocalBudget()
        g.dump(envelope=False, budget=b2)
        need = b.usage()['work_units'] + b2.usage()['work_units'] + 400
        with MT.local_budget(LocalLimits(work_units=need)):
            out = tools.graphlet_load(path)
        self.assertIn('handle', out, out)
        self.assertEqual([out['handle']], [e['handle'] for e in self.store.list()])
        if 'summary_stop' in out:
            self.assertIsNone(out['summary'])

    def test_graphlet_walks_stops_between_rows(self):
        # a lazy row stopped by the ambient budget after the first one must not write the
        # stop to the absent local call (AttributeError): the page states it, its rows whole,
        # within max_bytes, and its cursor goes on where it stopped
        h = self.store.put(resp('rep0_cut1'), {'seeds': []})
        for max_bytes in (60000, 3000):
            tools = GraphletTools(self.store, {}, max_bytes=max_bytes)
            whole = []
            page = tools.graphlet_walks(h, 'right', n=50, spell='full')
            while True:
                whole.extend(page['rows'])
                if 'next_cursor' not in page:
                    break
                page = tools.graphlet_walks(h, 'right', n=50, spell='full',
                                            cursor=page['next_cursor'])
            kinds = set()
            for units in range(200, 6000, 100):
                with MT.local_budget(LocalLimits(work_units=units)):
                    out = tools.graphlet_walks(h, 'right', n=50, spell='full')
                self.assertLessEqual(len(json.dumps(out, separators=(',', ':'))), max_bytes)
                if 'error' in out:
                    self.assertEqual('local_budget_exceeded', out['error'])
                    kinds.add('error')
                    continue
                self.assertEqual(whole[:len(out['rows'])], out['rows'])
                if 'stop' not in out:
                    kinds.add('page')
                    continue
                kinds.add('stopped')
                self.assertEqual('work', out['stop']['resource'])
                self.assertIs(False, out['complete'])
                got = list(out['rows'])
                page = out
                while 'next_cursor' in page:
                    page = tools.graphlet_walks(h, 'right', n=50, spell='full',
                                                cursor=page['next_cursor'])
                    got.extend(page['rows'])
                self.assertEqual(whole, got, units)
            self.assertEqual({'error', 'stopped', 'page'}, kinds, max_bytes)


class TestViewWalksBuildNoUnchargedSet(unittest.TestCase):
    def test_the_answers_and_charges(self):
        g = T.graphlet('merge')
        v = ops.GraphletView(g, [g.labels[0].ref])
        side = next(iter(v.segments))
        want = [w.path_id for w in g.walks(side) if w.path_id in set(v.path_ids[side])]
        self.assertEqual(want, [w.path_id for w in v.walks(side)])
        b1, b2 = LocalBudget(), LocalBudget()
        a = v.walks(side, budget=b1)
        v.walks(side, budget=b2)
        self.assertEqual(b1.usage(), b2.usage())
        self.assertEqual([w.path_id for w in a], want)
        with self.assertRaises(LocalBudgetExceeded):
            v.walks(side, budget=LocalBudget(work_units=0, memory_mb=0))


if __name__ == '__main__':
    unittest.main()
