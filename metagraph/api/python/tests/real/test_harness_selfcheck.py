"""Self-check of the real-data harness (realdata.py): the cache is consistent, the two
independent spellers agree with each other and with the library, and the /resolve
oracle confirms what any walk must satisfy. Skips without a cache or servers.

    cd api/python/tests/real && python3 -m unittest -v test_harness_selfcheck
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import realdata as R  # noqa: E402


def _small_cells(index=None, **kw):
    return [c for c in R.cells(index, max_full_bytes=8 << 20, **kw)]


@R.skip_unless_cache()
class TestCacheConsistency(unittest.TestCase):

    def test_every_ok_cell_has_its_files_and_matching_requests(self):
        for row in R.cells(heavy=None):
            c = R.load_cell(row['index'], row)
            for what in ('request', 'graphlet', 'full'):
                self.assertTrue(os.path.exists(c.paths[what]), (row['cell'], what))
            req = c.request
            self.assertEqual(req['strategy']['output']['detail'], 'graphlet', row['cell'])
            self.assertIs(req['strategy']['output']['timing'], False, row['cell'])
            self.assertEqual(len(req['seeds']), row['n_results'], row['cell'])

    def test_full_and_graphlet_have_the_same_results_shape(self):
        for row in _small_cells():
            c = R.load_cell(row['index'], row)
            g, f = c.graphlet_response, c.full_response
            self.assertEqual(len(g['results']), len(f['results']), row['cell'])
            for i, (rg, rf) in enumerate(zip(g['results'], f['results'])):
                self.assertEqual('error' in rg, 'error' in rf, (row['cell'], i))
                self.assertEqual(g['strategy']['output']['detail'], 'graphlet')
                self.assertEqual(f['strategy']['output']['detail'], 'full')

    def test_matrix_rows_name_known_seeds_and_strategies(self):
        for row in R.cells(status=None, heavy=None):
            self.assertIn(row['strategy'], R.STRATEGIES, row['cell'])
            self.assertEqual(row['cell'], R.cell_id(row['seed'], row['strategy']))


@R.skip_unless_cache()
class TestSpellers(unittest.TestCase):

    def test_graphlet_speller_equals_full_json_speller_and_library(self):
        checked = 0
        for row in _small_cells(deterministic=True):
            if 'sequences_false' in row.get('tags', ()):
                continue
            c = R.load_cell(row['index'], row)
            for i in range(c.n_results()):
                g = c.graphlet(i)
                if g is None:
                    continue
                full = c.full_result(i)
                seed = c.seed_sequence(i)
                self.assertEqual(g.seed.sequence, seed.upper(), row['cell'])
                for side in g.arms:
                    leaves = R.leaf_segments(g, side)
                    self.assertEqual(len(leaves), len(full['arms'][side]['paths']),
                                     (row['cell'], side))
                    for pid in range(min(len(leaves), 50)):
                        mine = R.spell_path_with_seed(g, side, pid)
                        theirs = R.spell_full_with_seed(full, side, pid, seed.upper())
                        lib = g.spell(side, pid, with_seed=True)
                        self.assertEqual(mine, theirs, (row['cell'], i, side, pid))
                        self.assertEqual(mine, lib, (row['cell'], i, side, pid))
                        checked += 1
        self.assertGreater(checked, 0)


@R.skip_unless_server('mini_refseq')
class TestOracleOnMiniRefseq(unittest.TestCase):

    def test_seeds_are_present_and_carried(self):
        for s in R.catalog_seeds('mini_refseq'):
            if s['kind'] == 'batch':
                continue
            self.assertTrue(R.oracle_present('mini_refseq', s['sequence']), s['name'])
            runs = R.oracle_label_runs('mini_refseq', s['sequence'], s['labels'][:5])
            for name, r in runs.items():
                self.assertEqual(r, [(0, s['num_kmers'])], (s['name'], name))

    def test_a_shorter_than_k_sequence_raises(self):
        with self.assertRaises(ValueError):
            R.oracle_present('mini_refseq', 'ACGT')

    def test_exhaustive_leaf_labels_cover_their_molecules(self):
        """Under exhaustive (keep, per path) a label alive at a leaf, entered by the
        seed, is present at every node of the walk: /resolve must report it on every
        k-mer of seed + flank (natural orientation). Uses only the full JSON."""
        R.require_cache('mini_refseq', 'mini_ndm1__exhaustive')
        c = R.load_cell('mini_refseq', 'mini_ndm1__exhaustive')
        full = c.full_result(0)
        names = [l['name'] for l in full['label_dict']]
        seed = c.seed_sequence(0)
        k = 31
        for side in ('left', 'right'):
            for pid, p in enumerate(full['arms'][side]['paths']):
                mol = R.spell_full_with_seed(full, side, pid, seed)
                self.assertTrue(R.oracle_present('mini_refseq', mol), (side, pid))
                alive = [names[e['label']] for e in p['end_labels'] if e['route_bp'] == 0]
                if not alive:
                    continue
                runs = R.oracle_label_runs('mini_refseq', mol, alive)
                for name in alive:
                    self.assertEqual(runs[name], [(0, len(mol) - k + 1)], (side, pid, name))


if __name__ == '__main__':
    unittest.main()
