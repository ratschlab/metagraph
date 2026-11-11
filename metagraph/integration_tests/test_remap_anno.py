import unittest
from parameterized import parameterized
import subprocess
from subprocess import PIPE
from tempfile import TemporaryDirectory
import glob
import gzip
import os
from base import PROTEIN_MODE, TestingBase, METAGRAPH, MMAP_FLAG

"""Test moving annotation rows from a succinct graph to an SSHash graph"""

TEST_DATA_DIR = os.path.dirname(os.path.realpath(__file__)) + '/../tests/data'
INPUT = TEST_DATA_DIR + '/transcripts_100.fa'
NUM_THREADS = 4
K = 21

anno_file_extension = {'brwt': '.brwt.annodbg',
                       'row_diff_brwt': '.row_diff_brwt.annodbg',
                       'row_diff_flat': '.row_diff_flat.annodbg',
                       'row_diff_disk': '.row_diff_disk.annodbg'}


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

        cls.primary_source, cls.primary_target, cls.flipped_target = cls._build_primary_graphs()

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

    @classmethod
    def _build_primary_graphs(cls):
        """A primary succinct graph, an SSHash graph built from its contigs, and an
        SSHash graph built from the reverse complements of those contigs."""
        prefix = cls.tempdir.name + '/primary'
        cls._build_graph(INPUT, prefix + '.dbg', K, 'succinct', 'primary')
        # The contigs of the primary graph spell every k-mer as the graph stores it.
        cls._run_command(
            f'{METAGRAPH} transform -p {NUM_THREADS} --to-fasta'
            f' -o {prefix}_contigs {prefix}.dbg', 'Transform command', shell=True)
        contigs = prefix + '_contigs.fasta.gz'
        cls._run_command(
            f'{METAGRAPH} build -p {NUM_THREADS} --mode primary --graph sshash'
            f' -k {K} -o {prefix}.sshashdbg {contigs}', 'Build command', shell=True)

        complement = str.maketrans('ACGT', 'TGCA')
        flipped = prefix + '_flipped.fa'
        with gzip.open(contigs, 'rt') as f, open(flipped, 'w') as out:
            for line in f:
                line = line.strip()
                out.write((line if line.startswith('>') else line.translate(complement)[::-1]) + '\n')
        cls._run_command(
            f'{METAGRAPH} build -p {NUM_THREADS} --mode primary --graph sshash'
            f' -k {K} -o {prefix}_flipped.sshashdbg {flipped}', 'Build command', shell=True)
        return prefix + '.dbg', prefix + '.sshashdbg', prefix + '_flipped.sshashdbg'

    def _query(self, graph, anno, mode):
        res = subprocess.run(
            f'{METAGRAPH} query --query-mode {mode} -i {graph} -a {anno} {INPUT}' + MMAP_FLAG,
            shell=True, stdout=PIPE, stderr=PIPE)
        self.assertEqual(res.returncode, 0, res.stderr.decode())
        return sorted(res.stdout.decode().splitlines())

    def _remap(self, anno, out, target, anno_type=None, in_memory=False, source=None):
        anno_type_flag = f' --anno-type {anno_type}' if anno_type else ''
        in_memory_flag = ' --disk-swap ""' if in_memory else ''
        return subprocess.run(
            f'{METAGRAPH} transform_anno -p {NUM_THREADS} --remap-graph {source or self.source}'
            f' -i {target}{anno_type_flag}{in_memory_flag} -o {out} {anno}',
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
        self.assertIn('does not match', res.stderr.decode() + res.stdout.decode())

    @parameterized.expand(['brwt', 'row_diff_brwt'])
    def test_remap_on_disk_and_in_memory(self, anno_repr):
        prefix = f'{self.tempdir.name}/disk_{anno_repr}'
        self._annotate_graph(INPUT, self.source, prefix, anno_repr,
                             num_threads=NUM_THREADS)
        anno = prefix + anno_file_extension[anno_repr]

        res = self._remap(anno, prefix + '_owned', self.target, in_memory=True)
        self.assertEqual(res.returncode, 0, res.stderr.decode())
        with open(prefix + '_owned' + anno_file_extension[anno_repr], 'rb') as f:
            expected = f.read()
        with open(prefix + '_owned.rowmap', 'rb') as f:
            expected_row_map = f.read()

        # By default the row map is built on disk. The second run reuses it from the mapped file.
        out = prefix + '_mapped'
        for _ in range(2):
            res = self._remap(anno, out, self.target)
            self.assertEqual(res.returncode, 0, res.stderr.decode())
            with open(out + anno_file_extension[anno_repr], 'rb') as f:
                self.assertEqual(expected, f.read())
            with open(out + '.rowmap', 'rb') as f:
                self.assertEqual(expected_row_map, f.read())
            self.assertEqual(glob.glob(out + '.rowmap.*'), [])
            self.assertEqual(glob.glob(f'{self.tempdir.name}/**/bucket.*', recursive=True), [])

    @parameterized.expand(['row_diff_flat', 'row_diff_disk'])
    def test_remap_and_convert(self, anno_type):
        prefix = f'{self.tempdir.name}/convert_{anno_type}'
        self._annotate_graph(INPUT, self.source, prefix, 'row_diff_brwt',
                             num_threads=NUM_THREADS)
        anno = prefix + anno_file_extension['row_diff_brwt']
        out = prefix + '_remapped'
        out_anno = out + anno_file_extension[anno_type]

        res = self._remap(anno, out, self.target, anno_type)
        self.assertEqual(res.returncode, 0, res.stderr.decode())
        self.assertFalse(os.path.exists(out + anno_file_extension['row_diff_brwt']))

        stats = self._get_stats(f'-a {out_anno}')
        self.assertEqual(stats['returncode'], 0)
        self.assertEqual(stats['objects'], self._get_stats(self.target)['nodes (k)'])
        for mode in ['labels', 'matches']:
            self.assertEqual(self._query(self.source, anno, mode),
                             self._query(self.target, out_anno, mode))

    @parameterized.expand([
        ('brwt', 'row_diff_flat'),
        ('row_diff_brwt', 'column'),
        ('row_diff_brwt', 'brwt'),
    ])
    def test_remap_rejects_conversion(self, anno_repr, anno_type):
        prefix = f'{self.tempdir.name}/reject_{anno_repr}_{anno_type}'
        self._annotate_graph(INPUT, self.source, prefix, anno_repr,
                             num_threads=NUM_THREADS)
        res = self._remap(prefix + anno_file_extension[anno_repr],
                          prefix + '_remapped', self.target, anno_type)
        self.assertNotEqual(res.returncode, 0)
        self.assertIn('--anno-type' if anno_repr == 'row_diff_brwt' else 'can convert only',
                      res.stderr.decode() + res.stdout.decode())

    @parameterized.expand(['brwt', 'row_diff_brwt'])
    def test_remap_primary_to_sshash(self, anno_repr):
        source, target = self.primary_source, self.primary_target
        prefix = f'{self.tempdir.name}/primary_{anno_repr}'
        self._annotate_graph(INPUT, source, prefix, anno_repr, num_threads=NUM_THREADS)
        anno = prefix + anno_file_extension[anno_repr]
        out = prefix + '_remapped'
        out_anno = out + anno_file_extension[anno_repr]

        res = self._remap(anno, out, target, source=source)
        self.assertEqual(res.returncode, 0, res.stderr.decode())

        # The SSHash graph also has a node, and so a row, for every reverse complement.
        stats = self._get_stats(f'-a {out_anno}')
        self.assertEqual(stats['objects'], self._get_stats(target)['nodes (k)'])
        for mode in ['labels', 'matches', 'signature']:
            expected = self._query(source, anno, mode)
            self.assertEqual(len(expected), 100)
            self.assertTrue(all(line.split('\t')[-1] for line in expected))
            self.assertEqual(expected, self._query(target, out_anno, mode))

        res = self._remap(anno, prefix + '_flipped', self.flipped_target, source=source)
        self.assertNotEqual(res.returncode, 0)
        self.assertIn('other orientation', res.stderr.decode() + res.stdout.decode())

    @parameterized.expand([
        ('basic', 'brwt', False),
        ('basic', 'row_diff_brwt', False),
        ('primary', 'row_diff_brwt', False),
        ('basic', 'row_diff_brwt', True),
        ('primary', 'row_diff_brwt', True),
    ])
    def test_row_map_from_extraction(self, mode, anno_repr, convert):
        """`transform --to-fasta --row-map` and `transform --graph sshash` write
        the row map that the remap computes."""
        source = self.source if mode == 'basic' else self.primary_source
        prefix = f'{self.tempdir.name}/extract_{mode}_{anno_repr}_{int(convert)}'
        target = prefix + '.sshashdbg'
        if convert:
            self._run_command(
                f'{METAGRAPH} transform -p {NUM_THREADS} --graph sshash'
                f' -o {prefix} {source}', 'Transform command', shell=True)
            self.assertEqual(glob.glob(f'{self.tempdir.name}/**/temp_convert_*', recursive=True), [])
        else:
            self._run_command(
                f'{METAGRAPH} transform -p {NUM_THREADS} --to-fasta --row-map'
                f' -o {prefix} {source}', 'Transform command', shell=True)
            self._run_command(
                f'{METAGRAPH} build -p {NUM_THREADS} --mode {mode} --graph sshash'
                f' -k {K} -o {target} {prefix}.fasta.gz', 'Build command', shell=True)
        self.assertEqual(glob.glob(f'{self.tempdir.name}/**/bucket.*', recursive=True), [])
        # In primary mode, the SSHash graph also has a node for every reverse complement.
        self.assertEqual(self._get_stats(target)['nodes (k)'],
                         self._get_stats(self.target if mode == 'basic'
                                         else self.primary_target)['nodes (k)'])

        self._annotate_graph(INPUT, source, prefix + '_anno', anno_repr,
                             num_threads=NUM_THREADS)
        anno = prefix + '_anno' + anno_file_extension[anno_repr]

        # Match the rows by k-mer, as without a row map.
        matched = prefix + '_matched'
        res = self._remap(anno, matched, target, source=source)
        self.assertEqual(res.returncode, 0, res.stderr.decode())
        with open(matched + '.rowmap', 'rb') as f:
            expected_row_map = f.read()
        with open(prefix + '.rowmap', 'rb') as f:
            self.assertEqual(expected_row_map, f.read())

        # The remap accepts the row map written with the contigs.
        res = self._remap(anno, prefix, target, source=source)
        self.assertEqual(res.returncode, 0, res.stderr.decode())
        with open(prefix + '.rowmap', 'rb') as f:
            self.assertEqual(expected_row_map, f.read())
        with open(matched + anno_file_extension[anno_repr], 'rb') as f:
            expected = f.read()
        with open(prefix + anno_file_extension[anno_repr], 'rb') as f:
            self.assertEqual(expected, f.read())
        for query_mode in ['labels', 'matches']:
            self.assertEqual(self._query(source, anno, query_mode),
                             self._query(target, prefix + anno_file_extension[anno_repr],
                                         query_mode))

    @parameterized.expand([(False, 1), (False, NUM_THREADS), (True, 1), (True, NUM_THREADS)])
    def test_convert_basic_succinct(self, clear_dummy, num_threads):
        """A basic succinct graph converts to the SSHash graph that a build from
        its contigs gives. Debug builds also compare the succ/pred bitmaps taken
        from the succinct graph with the ones SSHash computes."""
        prefix = f'{self.tempdir.name}/basic_{int(clear_dummy)}_{num_threads}'
        self._build_graph(INPUT, prefix + '.dbg', K, 'succinct', 'basic')
        source = prefix + '.dbg'
        if clear_dummy:
            self._run_command(
                f'{METAGRAPH} transform -p {NUM_THREADS} --clear-dummy'
                f' -o {prefix}_masked {source}', 'Transform command', shell=True)
            source = prefix + '_masked.dbg'
        self._run_command(
            f'{METAGRAPH} transform -p {num_threads} --graph sshash'
            f' -o {prefix}_converted {source}', 'Transform command', shell=True)
        self._run_command(
            f'{METAGRAPH} transform -p {num_threads} --to-fasta'
            f' -o {prefix}_contigs {source}', 'Transform command', shell=True)
        self._run_command(
            f'{METAGRAPH} build -p {NUM_THREADS} --mode basic --graph sshash'
            f' -k {K} -o {prefix}_built.sshashdbg {prefix}_contigs.fasta.gz',
            'Build command', shell=True)

        # the unitigs of an SSHash graph follow from its succ/pred bitmaps
        unitigs = []
        for graph in [f'{prefix}_converted', f'{prefix}_built']:
            self._run_command(
                f'{METAGRAPH} transform -p {NUM_THREADS} --to-fasta --unitigs'
                f' -o {graph}_unitigs {graph}.sshashdbg', 'Transform command', shell=True)
            with gzip.open(f'{graph}_unitigs.fasta.gz', 'rt') as f:
                unitigs.append(sorted(line for line in f if not line.startswith('>')))
        self.assertEqual(unitigs[0], unitigs[1])
        self.assertGreater(len(unitigs[0]), 0)

    @parameterized.expand(['--to-fasta --row-map', '--graph sshash'])
    def test_row_map_rejects_canonical(self, flags):
        prefix = self.tempdir.name + '/canonical'
        self._build_graph(INPUT, prefix + '.dbg', K, 'succinct', 'canonical')
        res = subprocess.run(
            f'{METAGRAPH} transform {flags} -o {prefix}_out {prefix}.dbg',
            shell=True, stdout=PIPE, stderr=PIPE)
        self.assertNotEqual(res.returncode, 0)
        self.assertIn('basic or primary', res.stderr.decode() + res.stdout.decode())

    @parameterized.expand([
        ('--graph sshash --to-fasta', '--to-fasta --row-map'),
        ('--graph sshash --initialize-bloom', 'cannot be combined'),
        ('--graph sshash --clear-dummy', 'cannot be combined'),
        ('--graph sshash --to-adj-list', 'cannot be combined'),
        ('--graph sshash --mode primary', 'keeps the mode'),
        ('--graph hash', 'only with'),
    ])
    def test_convert_rejects_options(self, flags, message):
        out = self.tempdir.name + '/rejected'
        res = subprocess.run(f'{METAGRAPH} transform {flags} -o {out} {self.source}',
                             shell=True, stdout=PIPE, stderr=PIPE)
        self.assertNotEqual(res.returncode, 0)
        self.assertIn(message, res.stderr.decode() + res.stdout.decode())
        self.assertEqual(glob.glob(out + '*'), [])

    def test_remap_rejects_mixed_modes(self):
        prefix = f'{self.tempdir.name}/mixed'
        self._annotate_graph(INPUT, self.source, prefix, 'brwt', num_threads=NUM_THREADS)
        res = self._remap(prefix + anno_file_extension['brwt'], prefix + '_remapped',
                          self.primary_target)
        self.assertNotEqual(res.returncode, 0)
        self.assertIn('basic mode', res.stderr.decode() + res.stdout.decode())

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
