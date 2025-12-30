import unittest
from parameterized import parameterized
import itertools
import subprocess
from subprocess import PIPE
from tempfile import TemporaryDirectory
import glob
import os
import gzip
from base import PROTEIN_MODE, TestingBase, METAGRAPH, TEST_DATA_DIR


"""Test graph construction"""

graph_file_extension = {'succinct': '.dbg',
                        'bitmap': '.bitmapdbg',
                        'hash': '.orhashdbg',
                        'hashfast': '.hashfastdbg',
                        'hashstr': '.hashstrdbg'}

build_params = {'succinct': ('succinct', '""'),
                'succinct_disk': ('succinct', '/tmp/'),  # build with disk swap
                'bitmap': ('bitmap', '""'),
                'hash': ('hash', '""'),
                'hashfast': ('hashfast', '""'),
                'hashstr': ('hashstr', '""')}

# Also test with swap in shm but only if it exists (Linux but not MacOS)
# (shm has a different filesystem, hence we're testing cross-device moves here)
if os.path.isdir("/dev/shm"):
    build_params['succinct_shm'] = ('succinct', '/dev/shm/')

BUILDS = [name for name, _ in build_params.items()]


class TestBuildWeighted(TestingBase):
    def setUp(self):
        self.tempdir = TemporaryDirectory()

    @parameterized.expand([repr for repr in BUILDS if not (repr == 'bitmap' and PROTEIN_MODE)])
    def test_simple_all_graphs(self, build):
        representation, tmp_dir = build_params[build]

        construct_command = '{exe} build --mask-dummy \
                --graph {repr} -k 20 --count-kmers --disk-swap {tmp_dir} -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=representation,
            tmp_dir=tmp_dir,
            outfile=self.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_graph = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[representation])
        self.assertEqual(stats_graph['returncode'], 0)
        self.assertEqual(stats_graph['k'], '20')
        self.assertEqual(stats_graph['nodes (k)'], '591997')
        self.assertEqual(stats_graph['mode'], 'basic')
        self.assertEqual(stats_graph['nnz weights'], '591997')
        self.assertEqual(stats_graph['avg weight'], '2.48587')

    @parameterized.expand([repr for repr in BUILDS if not (repr == 'bitmap' and PROTEIN_MODE)])
    def test_simple_all_graphs_contigs(self, build):
        representation, tmp_dir = build_params[build]

        construct_command = '{exe} build --mask-dummy \
                --graph {repr} -k 20 --count-kmers --disk-swap {tmp_dir} -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=representation,
            tmp_dir=tmp_dir,
            outfile=self.tempdir.name + '/graph_',
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        command = f'{METAGRAPH} transform --to-fasta \
                -o {self.tempdir.name}/graph_ {self.tempdir.name}/graph_{graph_file_extension[representation]}'

        res = subprocess.run([command], shell=True)
        self.assertEqual(res.returncode, 0)

        command = f'{METAGRAPH} build --mask-dummy \
                --graph {representation} -k 20 --count-kmers --disk-swap {tmp_dir} \
                -o {self.tempdir.name}/graph {self.tempdir.name}/graph_.fasta.gz'

        res = subprocess.run([command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_graph = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[representation])
        self.assertEqual(stats_graph['returncode'], 0)
        self.assertEqual(stats_graph['k'], '20')
        self.assertEqual(stats_graph['nodes (k)'], '591997')
        self.assertEqual(stats_graph['mode'], 'basic')
        self.assertEqual(stats_graph['nnz weights'], '591997')
        self.assertEqual(stats_graph['avg weight'], '2.48587')

    # TODO: add 'hashstr' once the canonical mode is implemented for it
    @parameterized.expand([repr for repr in BUILDS if repr != 'hashstr'])
    @unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets")
    def test_simple_all_graphs_canonical(self, build):
        representation, tmp_dir = build_params[build]

        construct_command = '{exe} build --mask-dummy \
                --graph {repr} --mode canonical --count-kmers --disk-swap {tmp_dir} -k 20 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=representation,
            tmp_dir=tmp_dir,
            outfile=self.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_graph = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[representation])
        self.assertEqual(stats_graph['returncode'], 0)
        self.assertEqual(stats_graph['k'], '20')
        self.assertEqual(stats_graph['nodes (k)'], '1159851')
        self.assertEqual(stats_graph['mode'], 'canonical')
        self.assertEqual(stats_graph['nnz weights'], '1159851')
        self.assertEqual(stats_graph['avg weight'], '2.53761')

    @parameterized.expand(BUILDS)
    def test_build_tiny_k(self, build):
        representation, tmp_dir = build_params[build]
        args = [METAGRAPH, 'build', '--mask-dummy', '--graph', representation,
                '--count-kmers',
                '--disk-swap', tmp_dir,
                '-k', '2',
                '-o', self.tempdir.name + '/graph',
                TEST_DATA_DIR + '/transcripts_1000.fa']
        construct_command = ' '.join(args)

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_graph = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[representation])
        self.assertEqual(stats_graph['returncode'], 0)
        self.assertEqual(stats_graph['k'], '2')
        self.assertEqual(stats_graph['nodes (k)'], '16')
        self.assertEqual(stats_graph['mode'], 'basic')
        self.assertEqual(stats_graph['nnz weights'], '16')
        self.assertEqual(stats_graph['avg weight'], '255')

    # TODO: add 'hashstr' once the canonical mode is implemented for it
    @parameterized.expand([repr for repr in BUILDS if repr != 'hashstr'])
    @unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets")
    def test_build_tiny_k_canonical(self, build):
        representation, tmp_dir = build_params[build]

        args = [METAGRAPH, 'build', '--mask-dummy', '--graph', representation, '--mode canonical',
                '--count-kmers',
                '--disk-swap', tmp_dir,
                '-k', '2',
                '-o', self.tempdir.name + '/graph',
                TEST_DATA_DIR + '/transcripts_1000.fa']
        construct_command = ' '.join(args)

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_graph = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[representation])
        self.assertEqual(stats_graph['returncode'], 0)
        self.assertEqual(stats_graph['k'], '2')
        self.assertEqual(stats_graph['nodes (k)'], '16')
        self.assertEqual(stats_graph['mode'], 'canonical')
        self.assertEqual(stats_graph['nnz weights'], '16')
        self.assertEqual(stats_graph['avg weight'], '255')

    @parameterized.expand(BUILDS)
    def test_build_from_kmc(self, build):
        representation, tmp_dir = build_params[build]

        construct_command = '{exe} build --mask-dummy \
                --graph {repr} --count-kmers --disk-swap {tmp_dir} -k 11 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=representation,
            tmp_dir=tmp_dir,
            outfile=self.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_1000_kmc_counters.kmc_suf'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_graph = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[representation])
        self.assertEqual(stats_graph['returncode'], 0)
        self.assertEqual(stats_graph['k'], '11')
        self.assertEqual(stats_graph['nodes (k)'], '469983')
        self.assertEqual(stats_graph['mode'], 'basic')
        self.assertEqual(stats_graph['nnz weights'], '469983')
        self.assertEqual(stats_graph['avg weight'], '3.15029')

    @parameterized.expand(BUILDS)
    def test_build_from_kmc_both(self, build):
        representation, tmp_dir = build_params[build]

        construct_command = '{exe} build --mask-dummy \
                --graph {repr} --count-kmers --disk-swap {tmp_dir} -k 11 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=representation,
            tmp_dir=tmp_dir,
            outfile=self.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_1000_kmc_counters_both_strands.kmc_suf'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_graph = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[representation])
        self.assertEqual(stats_graph['returncode'], 0)
        self.assertEqual(stats_graph['k'], '11')
        self.assertEqual(stats_graph['nodes (k)'], '802920')
        self.assertEqual(stats_graph['mode'], 'basic')
        self.assertEqual(stats_graph['nnz weights'], '802920')
        self.assertEqual(stats_graph['avg weight'], '3.68754')

    # TODO: add 'hashstr' once the canonical mode is implemented for it
    @parameterized.expand([repr for repr in BUILDS if repr != 'hashstr'])
    @unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets")
    def test_build_from_kmc_canonical(self, build):
        representation, tmp_dir = build_params[build]

        construct_command = '{exe} build --mask-dummy \
                --graph {repr} --count-kmers --disk-swap {tmp_dir} --mode canonical -k 11 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=representation,
            tmp_dir=tmp_dir,
            outfile=self.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_1000_kmc_counters.kmc_suf'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_graph = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[representation])
        self.assertEqual(stats_graph['returncode'], 0)
        self.assertEqual(stats_graph['k'], '11')
        self.assertEqual(stats_graph['nodes (k)'], '802920')
        self.assertEqual(stats_graph['mode'], 'canonical')
        self.assertEqual(stats_graph['nnz weights'], '802920')
        self.assertEqual(stats_graph['avg weight'], '3.68754')

    # TODO: add 'hashstr' once the canonical mode is implemented for it
    @parameterized.expand([repr for repr in BUILDS if repr != 'hashstr'])
    @unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets")
    def test_build_from_kmc_both_canonical(self, build):
        representation, tmp_dir = build_params[build]

        construct_command = '{exe} build --mask-dummy \
                --graph {repr} --count-kmers --disk-swap {tmp_dir} --mode canonical -k 11 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=representation,
            tmp_dir=tmp_dir,
            outfile=self.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_1000_kmc_counters_both_strands.kmc_suf'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_graph = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[representation])
        self.assertEqual(stats_graph['returncode'], 0)
        self.assertEqual(stats_graph['k'], '11')
        self.assertEqual(stats_graph['nodes (k)'], '802920')
        self.assertEqual(stats_graph['mode'], 'canonical')
        self.assertEqual(stats_graph['nnz weights'], '802920')
        self.assertEqual(stats_graph['avg weight'], '3.68754')

    @parameterized.expand(
        itertools.product(BUILDS,
                          [
                              (2, 3),
                              (3, 7),
                              (6, 63),
                              (8, 255),
                              (12, 3507.17),
                              (16, 5811.04),
                              (32, 5811.04),
                          ]
                          ))
    def test_kmer_count_width(self, build, width_result):
        representation, tmp_dir = build_params[build]
        count_width, avg_count_expected = width_result

        construct_command = '{exe} build --mask-dummy \
                --graph {repr} -k 4 --count-kmers --disk-swap {tmp_dir} --count-width {width} -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=representation,
            tmp_dir=tmp_dir,
            width=count_width,
            outfile=self.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_graph = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[representation])
        self.assertEqual(stats_graph['returncode'], 0)
        self.assertEqual(stats_graph['k'], '4')
        self.assertEqual(stats_graph['nodes (k)'], '256')
        self.assertEqual(stats_graph['mode'], 'basic')
        self.assertEqual(stats_graph['nnz weights'], '256')
        self.assertEqual(stats_graph['avg weight'], str(avg_count_expected))

    @parameterized.expand(itertools.chain(
        itertools.product(BUILDS,
                          [
                              (4, 2, 3),
                              (4, 6, 63),
                              (4, 8, 255),
                              (4, 12, 4095),
                              (4, 16, 65535),
                              (4, 32, 999998),

                              (29, 8, 255),
                              (29, 16, 65535),
                              (29, 32, 999986)
                          ]
                          ),
        itertools.product([repr for repr in BUILDS if repr != 'bitmap'],
                          [
                              (35, 8, 255),
                              (35, 16, 65535),
                              (35, 32, 999983),

                              (70, 8, 255),
                              (70, 16, 65535),
                              (70, 32, 999966),
                          ]
                          )
    ))
    @unittest.skipIf(PROTEIN_MODE, "Too large k-mer size for Protein alphabets")
    def test_kmer_count_width_large(self, build, k_width_result):
        representation, tmp_dir = build_params[build]
        k, count_width, avg_count_expected = k_width_result

        fasta_file = self.tempdir.name + '/CG_10_6.fasta.gz'
        with gzip.open(fasta_file, 'w') as f:
            f.write(b'>CG_10^6times\n')
            f.write(b'CG' * 10**6)

        construct_command = '{exe} build --mask-dummy \
                --graph {repr} -k {k} --count-kmers --disk-swap {tmp_dir} --count-width {width} -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=representation,
            tmp_dir=tmp_dir,
            k=k,
            width=count_width,
            outfile=self.tempdir.name + '/graph',
            input=fasta_file
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_graph = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[representation])
        self.assertEqual(stats_graph['returncode'], 0)
        self.assertEqual(stats_graph['k'], str(k))
        self.assertEqual(stats_graph['nodes (k)'], '2')
        self.assertEqual(stats_graph['mode'], 'basic')
        self.assertEqual(stats_graph['nnz weights'], '2')
        self.assertEqual(stats_graph['avg weight'], str(avg_count_expected))

    @parameterized.expand([repr for repr in BUILDS if not (repr == 'bitmap' and PROTEIN_MODE)])
    def test_header_abundance_counts(self, build):
        """Test --count-kmers with k-mer abundances from FASTA headers (Logan format) for all graph types"""
        representation, tmp_dir = build_params[build]
        fasta_path = os.path.join(TEST_DATA_DIR, 'logan_30.fa')
        outbase = os.path.join(self.tempdir.name, f'logan_graph_{representation}')
        cmd = f'{METAGRAPH} build --graph {representation} --count-kmers -k 31 -o {outbase} {fasta_path}'
        res = subprocess.run(cmd, shell=True)
        self.assertEqual(res.returncode, 0)
        weights_file = outbase + graph_file_extension[representation] + '.weights'
        self.assertTrue(os.path.exists(weights_file))
        stats_graph = self._get_stats(outbase + graph_file_extension[representation])
        self.assertEqual(stats_graph['returncode'], 0)
        self.assertEqual(stats_graph['k'], '31')
        self.assertEqual(stats_graph['nnz weights'], '728')
        self.assertEqual(stats_graph['avg weight'], '7.74863')
        self.assertEqual(stats_graph['mode'], 'basic')
        self.assertIn(stats_graph['nodes (k)'], ['728', '1079'])

if __name__ == '__main__':
    unittest.main()
