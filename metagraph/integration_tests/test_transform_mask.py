import filecmp
import os
import random
import re
import shlex
import shutil
import subprocess
import unittest

from base import PROTEIN_MODE, TestingBase, METAGRAPH, TEST_DATA_DIR

"""
`metagraph transform --mask-dummy` (docs/DESIGN-pattern-search.md §4): an existing succinct graph
gets the dummy-edge mask `metagraph build --mask-dummy` would have written, as the file
<graph>.edgemask beside it, and nothing else changes. The reference is always the build's own
mask: each graph is built with --mask-dummy, a copy of its .dbg without the mask is given one by
transform, and the two .edgemask files must be byte for byte the same. Also the counts the step
logs (against `stats --count-dummy`), the refusals, --force, and the server's --pattern-build-mask
and transform's -o rules on the command line.
"""

K = 11


def run(args):
    res = subprocess.run(shlex.split(METAGRAPH) + args, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    return res.returncode, (res.stdout + res.stderr).decode()


def dummy_counts(graph):
    """`metagraph stats --count-dummy`: (source dummies, sink dummies, real edges)."""
    code, out = run(['stats', '--count-dummy', graph])
    assert code == 0, out
    get = lambda key: int(re.search(rf'^{key}: (\d+)$', out, re.M).group(1))
    return get('dummy source edges'), get('dummy sink edges'), get('real edges')


@unittest.skipIf(PROTEIN_MODE, "the masks are compared on DNA graphs")
class TestTransformMaskDummy(TestingBase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        d = cls.tempdir.name
        # random records, with an N run and records starting alike (a branching dummy tree),
        # and the transcripts fixture: many records, many starts and ends
        rng = random.Random(11)
        records = [''.join(rng.choice('ACGT') for _ in range(rng.randint(30, 400)))
                   for _ in range(40)]
        records += ['GATTACA' + r[:50] for r in records[:5]]
        records[3] = records[3][:20] + 'NNNN' + records[3][24:]
        cls.fasta = os.path.join(d, 'records.fa')
        with open(cls.fasta, 'w') as f:
            for i, r in enumerate(records):
                f.write(f'>r{i}\n{r}\n')
        cls.inputs = [cls.fasta, TEST_DATA_DIR + '/transcripts_100.fa']

    def built_and_stripped(self, name, mode='basic', state='stat'):
        """(the graph built with --mask-dummy, a copy of its .dbg without the mask)"""
        d = os.path.join(self.tempdir.name, name)
        os.makedirs(d)
        built = os.path.join(d, 'built.dbg')
        self._build_graph(self.inputs, built, K, 'succinct', mode=mode,
                          extra_params=f'--mask-dummy --in-ram --state {state}')
        self.assertTrue(os.path.isfile(os.path.join(d, 'built.edgemask')))
        stripped = os.path.join(d, 'stripped.dbg')
        shutil.copyfile(built, stripped)
        return built, stripped

    def assertSameMask(self, built, stripped):
        self.assertTrue(filecmp.cmp(built[:-4] + '.edgemask', stripped[:-4] + '.edgemask',
                                    shallow=False))

    def test_every_mode_and_state(self):
        for mode in ('basic', 'canonical', 'primary'):
            for state in ('stat', 'fast', 'small'):
                with self.subTest(mode=mode, state=state):
                    built, stripped = self.built_and_stripped(f'{mode}_{state}', mode, state)
                    code, log = run(['transform', '--mask-dummy', '-p', '2', stripped])
                    self.assertEqual(0, code, log)
                    self.assertSameMask(built, stripped)
                    # the graph file is untouched and nothing else was written
                    self.assertTrue(filecmp.cmp(built, stripped, shallow=False))
                    self.assertEqual({'built.dbg', 'built.edgemask', 'stripped.dbg',
                                      'stripped.edgemask'},
                                     set(os.listdir(os.path.dirname(stripped))) - {
                                         'built.dbg.fasta.gz'})
                    # the counts it logs are those of stats --count-dummy
                    source, sink, real = dummy_counts(built)
                    m = re.search(r'Dummy edges marked in [\d.]+ s with \d+ threads: (\d+) edges, '
                                  r'(\d+) source dummies \(the main dummy edge included\), (\d+) '
                                  r'sink dummies, (\d+) k-mers', log)
                    self.assertIsNotNone(m, log)
                    self.assertEqual((source + sink + real, source, sink, real),
                                     tuple(int(x) for x in m.groups()))
                    self.assertGreater(source, 1)
                    self.assertGreater(sink, 0)
                    # the graph states the same k-mers as the built one
                    self.assertEqual(self._get_stats(built)['nodes (k)'],
                                     self._get_stats(stripped)['nodes (k)'])
                    self.assertEqual(str(real), self._get_stats(stripped)['nodes (k)'])

    def test_with_mmap_and_output_naming_the_graph(self):
        built, stripped = self.built_and_stripped('mmap')
        code, log = run(['transform', '--mask-dummy', '--mmap', '-o', stripped[:-4], stripped])
        self.assertEqual(0, code, log)
        self.assertSameMask(built, stripped)

    def test_existing_mask_and_force(self):
        built, stripped = self.built_and_stripped('force')
        mask = stripped[:-4] + '.edgemask'
        with open(mask, 'w') as f:
            f.write('not a mask')
        code, log = run(['transform', '--mask-dummy', stripped])
        self.assertEqual(1, code, log)
        self.assertIn('pass --force to build it again and replace it', log)
        with open(mask) as f:
            self.assertEqual('not a mask', f.read())
        # a broken mask makes the graph unloadable; --force rebuilds it without reading it
        code, log = run(['transform', '--mask-dummy', '--force', stripped])
        self.assertEqual(0, code, log)
        self.assertIn('Replaced', log)
        self.assertSameMask(built, stripped)
        self.assertEqual([], [f for f in os.listdir(os.path.dirname(stripped)) if '.tmp' in f])

    def test_bloom_filter_beside_the_graph(self):
        """DBGSuccinct::load reads a .bloom only together with a mask: stated when it applies."""
        built, stripped = self.built_and_stripped('bloom')
        code, log = run(['transform', '--initialize-bloom', '--bloom-fpp', '0.1',
                         '-o', stripped[:-4], stripped])
        self.assertEqual(0, code, log)
        self.assertTrue(os.path.isfile(stripped[:-4] + '.bloom'))
        code, log = run(['transform', '--mask-dummy', stripped])
        self.assertEqual(0, code, log)
        self.assertIn('.bloom exists beside the graph: DBGSuccinct::load reads a Bloom filter '
                      'only when the mask is present', log)
        self.assertSameMask(built, stripped)

    def test_refusals(self):
        built, stripped = self.built_and_stripped('refusals')
        other = os.path.join(self.tempdir.name, 'elsewhere')
        refused = [
            # another transformation asked with it would be dropped unseen
            (['--clear-dummy'], 'cannot be combined'),
            (['--state', 'fast'], 'cannot be combined'),
            (['--index-ranges', '5'], 'cannot be combined'),
            (['--to-adj-list'], 'cannot be combined'),
            (['--initialize-bloom'], 'cannot be combined'),
            (['--to-fasta'], 'cannot be combined'),
            (['--mode', 'primary'], 'cannot be combined'),
            # the loader reads the mask beside the graph only
            (['-o', other], 'writes the mask beside the graph'),
        ]
        for extra, message in refused:
            with self.subTest(extra=extra):
                code, log = run(['transform', '--mask-dummy'] + extra + [stripped])
                self.assertNotEqual(0, code, log)
                self.assertIn(message, log)
                self.assertFalse(os.path.exists(stripped[:-4] + '.edgemask'))
                self.assertFalse(os.path.exists(other + '.edgemask'))
        # --force only with --mask-dummy, --pattern-build-mask only where a graph is served
        code, log = run(['transform', '--force', '-o', other, stripped])
        self.assertNotEqual(0, code)
        self.assertIn('--force applies only to transform --mask-dummy', log)
        code, log = run(['stats', '--pattern-build-mask', stripped])
        self.assertNotEqual(0, code)
        self.assertIn('--pattern-build-mask applies to server_query with one graph', log)
        csv = os.path.join(self.tempdir.name, 'graphs.csv')
        with open(csv, 'w') as f:
            f.write(f'G,{stripped},{stripped}.column.annodbg\n')
        code, log = run(['server_query', '--pattern-build-mask', csv])
        self.assertNotEqual(0, code)
        self.assertIn('--pattern-build-mask applies to server_query with one graph', log)
        # only a succinct graph has a mask
        hashed = os.path.join(self.tempdir.name, 'hashed')
        self._build_graph(self.fasta, hashed, K, 'hash')
        code, log = run(['transform', '--mask-dummy', hashed + '.orhashdbg'])
        self.assertEqual(1, code, log)
        self.assertIn('is not a succinct graph', log)
        # nothing was written by any of them
        self.assertFalse(os.path.exists(stripped[:-4] + '.edgemask'))
        self.assertTrue(filecmp.cmp(built, stripped, shallow=False))


if __name__ == '__main__':
    unittest.main()
