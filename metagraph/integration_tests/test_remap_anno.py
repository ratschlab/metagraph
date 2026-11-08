import unittest
from parameterized import parameterized
import subprocess
from subprocess import PIPE
from tempfile import TemporaryDirectory
import os
from base import PROTEIN_MODE, TestingBase, METAGRAPH, MMAP_FLAG

"""Test moving annotation rows from a succinct graph to an SSHash graph"""

TEST_DATA_DIR = os.path.dirname(os.path.realpath(__file__)) + '/../tests/data'
INPUT = TEST_DATA_DIR + '/transcripts_100.fa'
NUM_THREADS = 4
K = 21

anno_file_extension = {'brwt': '.brwt.annodbg',
                       'row_diff_brwt': '.row_diff_brwt.annodbg'}


@unittest.skipIf(PROTEIN_MODE, 'SSHash graphs need DNA')
class TestRemapAnnotation(TestingBase):
    @classmethod
    def setUpClass(cls):
        cls.tempdir = TemporaryDirectory()
        cls.source = cls.tempdir.name + '/graph.dbg'
        cls.target = cls.tempdir.name + '/graph.sshashdbg'
        # Keep dummy nodes unmasked so that the source has rows without k-mers.
        cls._build_graph(INPUT, cls.source, K, 'succinct', 'basic')
        cls._build_sshash(INPUT, cls.target)

        cls.partial_input = cls.tempdir.name + '/partial.fa'
        with open(INPUT) as f:
            records = f.read().split('>')[1:]
        with open(cls.partial_input, 'w') as f:
            f.write(''.join('>' + r for r in records[:50]))
        cls.partial_target = cls.tempdir.name + '/partial.sshashdbg'
        cls._build_sshash(cls.partial_input, cls.partial_target)

    @staticmethod
    def _build_sshash(input, output):
        # SSHash needs each k-mer once, so build it from the unitigs of a succinct graph.
        prefix = output.removesuffix('.sshashdbg')
        TestingBase._build_graph(input, prefix + '_unitigs.dbg', K, 'succinct', 'basic')
        TestingBase._run_command(
            f'{METAGRAPH} transform -p {NUM_THREADS} --to-fasta -o {prefix}_unitigs'
            f' {prefix}_unitigs.dbg', 'Transform command', shell=True)
        TestingBase._run_command(
            f'{METAGRAPH} build -p {NUM_THREADS} --mode basic --graph sshash'
            f' -k {K} -o {output} {prefix}_unitigs.fasta.gz', 'Build command', shell=True)

    def _query(self, graph, anno, mode):
        res = subprocess.run(
            f'{METAGRAPH} query --query-mode {mode} -i {graph} -a {anno} {INPUT}' + MMAP_FLAG,
            shell=True, stdout=PIPE, stderr=PIPE)
        self.assertEqual(res.returncode, 0, res.stderr.decode())
        return sorted(res.stdout.decode().splitlines())

    def _remap(self, anno, out, target):
        return subprocess.run(
            f'{METAGRAPH} transform_anno -p {NUM_THREADS} --remap-graph {self.source}'
            f' -i {target} -o {out} {anno}',
            shell=True, stdout=PIPE, stderr=PIPE)

    @parameterized.expand([
        ('row_diff_brwt', False),
        ('row_diff_brwt', True),
        ('brwt', False),
    ])
    def test_remap_to_sshash(self, anno_repr, no_fork_opt):
        prefix = f'{self.tempdir.name}/{anno_repr}_{int(no_fork_opt)}'
        self._annotate_graph(INPUT, self.source, prefix, anno_repr,
                             no_fork_opt=no_fork_opt, num_threads=NUM_THREADS)
        anno = prefix + anno_file_extension[anno_repr]
        out = prefix + '_remapped'
        out_anno = out + anno_file_extension[anno_repr]

        res = self._remap(anno, out, self.target)
        self.assertEqual(res.returncode, 0, res.stderr.decode())

        stats = self._get_stats(f'-a {out_anno}')
        self.assertEqual(stats['returncode'], 0)
        self.assertEqual(stats['labels'], '100')
        target_stats = self._get_stats(self.target)
        self.assertEqual(stats['objects'], target_stats['nodes (k)'])

        for mode in ['labels', 'matches', 'signature']:
            expected = self._query(self.source, anno, mode)
            self.assertEqual(len(expected), 100)
            self.assertTrue(all(line.split('\t')[-1] for line in expected))
            self.assertEqual(expected, self._query(self.target, out_anno, mode))

        # A second run reuses the saved row map and writes the same annotation.
        with open(out_anno, 'rb') as f:
            first = f.read()
        self.assertTrue(os.path.exists(out + '.rowmap'))
        res = self._remap(anno, out, self.target)
        self.assertEqual(res.returncode, 0, res.stderr.decode())
        with open(out_anno, 'rb') as f:
            self.assertEqual(first, f.read())

        # The saved row map belongs to the first target and is rejected for another one.
        res = self._remap(anno, out, self.partial_target)
        self.assertNotEqual(res.returncode, 0)
        self.assertIn('different version', res.stderr.decode() + res.stdout.decode())

    def test_remap_rejects_different_kmers(self):
        prefix = f'{self.tempdir.name}/mismatch'
        self._annotate_graph(INPUT, self.source, prefix, 'row_diff_brwt',
                             num_threads=NUM_THREADS)
        res = self._remap(prefix + anno_file_extension['row_diff_brwt'],
                          prefix + '_remapped', self.partial_target)
        self.assertNotEqual(res.returncode, 0)
        self.assertIn('different k-mers', res.stderr.decode() + res.stdout.decode())


if __name__ == '__main__':
    unittest.main()
