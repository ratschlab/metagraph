"""Offline regression over real retrievals: the CI side of tests/real/offline.py.

Eighteen small retrievals of the live real-data suite (tests/real), committed under
data/traverse/real/ together with the /resolve answers their oracle checks need. Each is a
real server response from one of three indexes:

  sra           the primary regime (42 G nodes, RowDiff<BRWT>, SRA runs as column labels,
                no manifest: index identity unverifiable)
  uhgg          the basic regime (9.7 G nodes, 4,644 genome labels, unverifiable identity)
  mini_refseq   the refseq33m format (k-mer coordinates, header labels via .seqs, an index
                manifest: identity verifiable)

The live suite's own check code runs on them, unchanged, through private copies of
tests/real/{realdata,test_real_conformance,test_real_oracle}.py whose harness points at the
fixtures with every server off:

  TestRealOfflineConformance  roundtrip, counts, rules, to_json, invariants, summary,
                              identity, seed, save_load (+ full_invariants on the
                              time-budgeted cell, whose two fetches stopped apart)
  TestRealOfflineOracle       a_present, b_claims, b_lost, b_walks, b_stretch, b_maximal,
                              c_routes, d_annotate, e_direct, e_library, answered from the
                              recorded /resolve answers, and f_positions (record
                              coordinates) from the recorded source-record answers
  TestRealOfflineFixtures     the files and their sizes, the situations the selection is
                              meant to cover (read from the bodies), the hairpin and
                              reverse-complement facts of the primary regime, no socket

One test per (cell, check). A check that RAN when the fixtures were made (manifest
'expect') must still run: a skip of it fails, because a skip would hide a regression. A
/resolve payload without a recorded answer fails too (OfflineMiss): the oracle plan
changed, so the library's output did. If that change is intended, rebuild the fixtures
while the servers run (python3 tests/real/offline.py snapshot) and review the diff.

    python3 -m unittest discover -s api/python/tests -p 'test_traverse_real_offline.py'
"""

import collections
import gzip
import importlib.util
import json
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
API = os.path.dirname(HERE)
if API not in sys.path:
    sys.path.insert(0, API)


def _load_offline():
    path = os.path.join(HERE, 'real', 'offline.py')
    spec = importlib.util.spec_from_file_location('_mgt_traverse_real_offline', path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


offline = _load_offline()
H = offline.Harness()          # the committed fixtures; missing fixtures fail loudly here


def _make(check, row, strict):
    def test(self):
        try:
            getattr(self, 'check_' + check)(row)
        except unittest.SkipTest as e:
            if strict:
                self.fail('%s on %s/%s ran when the fixtures were made and now skips: %s'
                          % (check, row['index'], row['cell'], e))
            raise
    test.__doc__ = '%s: %s/%s%s' % (check, row['index'], row['cell'],
                                    '' if strict else ' (skipped when recorded)')
    return test


def _build_classes():
    attrs = {'conformance': {}, 'oracle': {}}
    for row in H.rows:
        exp = row.get('expect') or {}
        ran = exp.get('ran') or {}
        for kind, check in H.checks(row):
            strict = check in (ran.get(kind) or ())
            attrs[kind]['test_%s__%s' % (row['cell'], check)] = _make(check, row, strict)
    out = {}
    for name, kind, base in (
            ('TestRealOfflineConformance', 'conformance', offline.conformance_base(H)),
            ('TestRealOfflineOracle', 'oracle', offline.oracle_base(H))):
        cls = type(name, (base,), dict(attrs[kind], __module__=__name__))
        cls.__doc__ = 'The real suite\'s %s checks on %d fixture cells (%d tests)' % (
            kind, len(H.rows), len(attrs[kind]))
        out[name] = cls
    return out


globals().update(_build_classes())


# situations the selection must keep covering, per index (offline.situations reads them
# from the bodies): the primary regime with its reverse-complement specifics, merges and a
# switch on SRA, a contig end switching into the genome that continues past it and a
# multi-genome window on UHGG, header labels / trace / verifiable identity on mini_refseq,
# both label modes, a cut list, cut evidence, a time-budgeted partial result, a per-seed
# refusal in a batch
REQUIRED = {
    'sra': ['regime_primary', 'identity_unverifiable', 'labels_column', 'mode_constrain',
            'mode_annotate', 'reconverge_merge', 'merge', 'merge_closed_run', 'switch',
            'end_loss_budget', 'hairpin_skipped', 'end_dead_end_hairpin', 'hairpin_followed',
            'end_edge_reuse_rc', 'cut_list', 'cut_list_truncated_nodes',
            'label_evidence_lower_bound', 'cut_events', 'time_budget', 'walks_partial',
            'derivation_error'],
    'uhgg': ['regime_basic', 'identity_unverifiable', 'labels_column', 'mode_constrain',
             'mode_annotate', 'merge', 'switch', 'switch_into_extra', 'reconverge_keep',
             'multi_genome_seed'],
    'mini_refseq': ['regime_basic', 'identity_verifiable', 'labels_header', 'support_trace',
                    'mode_constrain', 'mode_annotate', 'switch', 'end_edge_reuse',
                    'reverse_complement_seed_pair', 'coordinates_record',
                    'coordinates_column'],
}


class TestRealOfflineFixtures(unittest.TestCase):
    """The fixture set itself: complete, small, covering what it claims, and offline."""

    @classmethod
    def setUpClass(cls):
        cls.found = offline.situations(H)

    def test_every_cell_has_its_files_and_each_file_is_small(self):
        self.assertGreaterEqual(len(H.rows), 12)
        for row in H.rows:
            for what in ('request', 'graphlet', 'full', 'oracle'):
                p = os.path.join(H.fixtures, row['index'], '%s.%s.json.gz' % (row['cell'], what))
                with self.subTest(cell=row['cell'], file=what):
                    self.assertTrue(os.path.exists(p), p)
        for dirpath, _, files in os.walk(H.fixtures):
            for fn in files:
                p = os.path.join(dirpath, fn)
                if fn.endswith('.gz'):
                    with self.subTest(file=os.path.relpath(p, H.fixtures)):
                        self.assertLess(os.path.getsize(p), offline.MAX_FILE_BYTES)

    def test_every_index_is_present_with_its_identity(self):
        caps = H.R.recorded_capabilities()
        self.assertEqual(sorted({r['index'] for r in H.rows}), sorted(H.R.INDEXES))
        self.assertEqual(caps['sra']['capabilities']['regime'], 'primary')
        self.assertEqual(caps['uhgg']['capabilities']['regime'], 'basic')
        self.assertIsNone(caps['sra']['capabilities']['index_fp'])
        self.assertIsNone(caps['uhgg']['capabilities']['index_fp'])
        self.assertRegex(caps['mini_refseq']['capabilities']['index_fp'], r'\A[0-9a-f]{64}\Z')
        servers = H.manifest()['servers']
        for index in H.R.INDEXES:
            for k in ('index_ns', 'index_fp', 'index_meta_fp', 'k', 'regime'):
                self.assertEqual(servers[index][k], caps[index]['capabilities'][k], (index, k))

    def test_the_selection_covers_its_situations(self):
        for index, wants in REQUIRED.items():
            for want in wants:
                with self.subTest(index=index, situation=want):
                    cells = [c for c in self.found.get(want, ()) if c.startswith(index + '/')]
                    self.assertTrue(cells, '%s: no fixture cell shows %s (found: %s)' % (
                        index, want, sorted(k for k, v in self.found.items()
                                            if any(c.startswith(index + '/') for c in v))))

    def test_every_check_area_runs_somewhere(self):
        """Each conformance check and each oracle area ran on at least one cell when the
        fixtures were made (the strict tests then keep it running)."""
        ran = collections.Counter()
        for row in H.rows:
            for kind, checks in (row['expect']['ran']).items():
                ran.update(checks)
        for check in tuple(H.C.CHECKS) + ('full_invariants',) + tuple(H.O.AREAS):
            with self.subTest(check=check):
                self.assertGreater(ran[check], 0)

    def test_recorded_skips_are_data_driven(self):
        """What may skip is only what the data does not contain (no merge-derived support,
        no label_lost inside a segment) or the non-deterministic cell's comparisons."""
        allowed = ('no run whose displayed support starts at a merge',
                   'no label_lost / switch end inside a segment', 'not deterministic')
        for row in H.rows:
            for check, why in row['expect']['skipped'].items():
                with self.subTest(cell=row['cell'], check=check):
                    self.assertTrue(any(a in why for a in allowed), why)
                    if 'not deterministic' in why:
                        self.assertFalse(row.get('deterministic'))
                        self.assertIn('full_invariants', row['expect']['ran']['conformance'])

    def test_a_hairpin_skipped_and_followed_on_the_primary_regime(self):
        """Spec §6.5 on real primary data, independently of the oracle: the hairpin step is
        a self-reverse-complementary (k+1)-mer; skipped, every label of the event ends
        dead_end/hairpin there; followed, the reverse-complement retrace one step later
        ends edge_reuse_rc."""
        R = H.R
        seen = {}
        for strategy, followed in (('left_only', False), ('left_only_follow', True)):
            row = R.cell_meta('sra', 'sra_hairpin__' + strategy)
            g = R.load_cell('sra', row).graphlet(0)
            self.assertEqual(g.regime, 'primary')
            a = g.arms['left']
            events = [(s, e) for s in a.segments for e in s.events if e.type == 'hairpin']
            self.assertEqual(len(events), 1, strategy)
            s, e = events[0]
            self.assertIs(e.followed, followed)
            mol = R.spell_with_seed(g, 'left', s.id)
            j = s.end_bp - e.at_bp                       # the node at outward at_bp - 1
            step = e.char + mol[j:j + g.k]               # the left arm prepends
            self.assertEqual(len(step), g.k + 1)
            self.assertEqual(R.revcomp(step), step, 'not a self-RC (k+1)-mer: %s' % step)
            labels = set(e.labels)
            self.assertEqual(len(labels), e.total)
            dh = [r for r in a.runs if r.end == 'Dh']
            v = [r for r in a.runs if r.end == 'V']
            if not followed:
                self.assertEqual({r.label for r in dh}, labels)
                self.assertTrue(all((r.reason, r.qualifier) == ('dead_end', 'hairpin')
                                    and r.segment == s.id and r.to_bp == e.at_bp for r in dh))
                self.assertEqual(v, [])
            else:
                self.assertEqual(dh, [])
                self.assertTrue(v)
                self.assertLessEqual({r.label for r in v}, labels)
                self.assertTrue(all(r.reason == 'edge_reuse_rc' and r.to_bp == e.at_bp + 1
                                    for r in v))
            seen[strategy] = (g.seed.sequence, e.at_bp, {g.labels[x].name for x in labels})
        self.assertEqual(seen['left_only'], seen['left_only_follow'],
                         'one seed, one decision point, one label set')

    def test_reverse_complement_seeds_on_a_basic_graph(self):
        """blaNDM-1 and its reverse complement on mini_refseq (basic: one strand): two
        seeds, each fully carried by header labels, both strands walkable."""
        R = H.R
        fwd = R.load_cell('mini_refseq', R.cell_meta('mini_refseq', 'mini_ndm1__trace'))
        rev = R.load_cell('mini_refseq', R.cell_meta('mini_refseq', 'mini_ndm1_rc__switch1'))
        self.assertEqual(R.revcomp(rev.seed_sequence()), fwd.seed_sequence())
        for c in (fwd, rev):
            g = c.graphlet(0)
            self.assertEqual(g.regime, 'basic')
            self.assertTrue(g.seed.num_seed_labels > 0)
            self.assertTrue(all(l.kind in ('h', 'header') for l in g.labels), c.cell)

    def test_the_contig_end_switches_into_the_genome_past_it(self):
        R = H.R
        row = R.cell_meta('uhgg', 'uhgg_contig_end__switch1_limit1_keep')
        cell = R.load_cell('uhgg', row)
        extra = set(cell.request['strategy']['labels']['extra'])
        g = cell.graphlet(0)
        seed_labels = {l.name for l in g.labels[:g.seed.num_seed_labels]}
        self.assertFalse(extra & seed_labels, 'the switch target does not carry the seed')
        into = collections.Counter()
        for side, a in g.arms.items():
            for s in a.segments:
                for e in s.events:
                    if e.type == 'switch':
                        into[g.labels[e.to_label].name in extra] += 1
            for r in a.runs:
                if r.from_label is not None:
                    self.assertIsNotNone(r.prev_run)
        self.assertGreater(into[True], 0)
        self.assertEqual(g.outcome.walks, 'complete')

    def test_nothing_reaches_a_server(self):
        R = H.R
        self.assertTrue(all(v is None for v in R.SERVERS.values()))
        self.assertFalse(R.server_up('sra'))
        with self.assertRaises(offline.OfflineMiss):
            R.http('sra', 'GET', 'traverse/capabilities')
        with self.assertRaises(offline.OfflineMiss):       # an unrecorded payload
            R.oracle_resolve('sra', 'ACGT' * 20)

    def test_the_oracle_answers_are_slimmed_copies_of_real_answers(self):
        for row in H.rows:
            with gzip.open(offline.bundle_path(H.fixtures, row['index'], row['cell'])) as f:
                b = json.load(f)
            with self.subTest(cell=row['cell']):
                self.assertEqual((b['index'], b['cell']), (row['index'], row['cell']))
                self.assertTrue(b['answers'])
                for key, e in b['answers'].items():
                    self.assertRegex(key, r'\A[0-9a-f]{64}\Z')
                    self.assertNotIn('candidates', e['response'])
                    if 'positions' in e['request']:
                        # a record-oracle answer (f_positions): per slice the digests of the
                        # bases there, per counted string its count (the strings kept as
                        # their length and digest)
                        n = len(e['request']['positions']['slices'])
                        self.assertEqual(n, len(e['response']['slices']))
                        self.assertTrue(all(isinstance(x, list) for x in e['response']['slices']))
                        self.assertTrue(all(isinstance(x, int) for x in e['response']['counts']))
                        continue
                    self.assertEqual(e['response']['num_kmers'],
                                     e['request']['sequence_bp'] - e['response']['k'] + 1)


if __name__ == '__main__':
    unittest.main()
