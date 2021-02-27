import unittest
from parameterized import parameterized
import itertools
import subprocess
from subprocess import PIPE
from tempfile import TemporaryDirectory
import glob
import os
import gzip


"""Test graph construction"""

METAGRAPH = './metagraph'
PROTEIN_MODE = os.readlink(METAGRAPH).endswith("_Protein")
TEST_DATA_DIR = os.path.dirname(os.path.realpath(__file__)) + '/../tests/data'

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

BUILDS = [name for name, _ in build_params.items()]


class TestBuildWeighted(unittest.TestCase):
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

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[representation],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 20', params_str[0])
        self.assertEqual('nodes (k): 591997', params_str[1])
        self.assertEqual('mode: basic', params_str[2])
        self.assertEqual('nnz weights: 591997', params_str[3])
        self.assertEqual('avg weight: 2.48587', params_str[4])

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

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[representation],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 20', params_str[0])
        self.assertEqual('nodes (k): 1159851', params_str[1])
        self.assertEqual('mode: canonical', params_str[2])
        self.assertEqual('nnz weights: 1159851', params_str[3])
        self.assertEqual('avg weight: 2.53761', params_str[4])

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

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[representation],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 2', params_str[0])
        self.assertEqual('nodes (k): 16', params_str[1])
        self.assertEqual('mode: basic', params_str[2])
        self.assertEqual('nnz weights: 16', params_str[3])
        self.assertEqual('avg weight: 255', params_str[4])

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

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[representation],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 2', params_str[0])
        self.assertEqual('nodes (k): 16', params_str[1])
        self.assertEqual('mode: canonical', params_str[2])
        self.assertEqual('nnz weights: 16', params_str[3])
        self.assertEqual('avg weight: 255', params_str[4])

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

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[representation],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 11', params_str[0])
        self.assertEqual('nodes (k): 469983', params_str[1])
        self.assertEqual('mode: basic', params_str[2])
        self.assertEqual('nnz weights: 469983', params_str[3])
        self.assertEqual('avg weight: 3.15029', params_str[4])

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

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[representation],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 11', params_str[0])
        self.assertEqual('nodes (k): 802920', params_str[1])
        self.assertEqual('mode: basic', params_str[2])
        self.assertEqual('nnz weights: 802920', params_str[3])
        self.assertEqual('avg weight: 3.68754', params_str[4])

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

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[representation],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 11', params_str[0])
        self.assertEqual('nodes (k): 802920', params_str[1])
        self.assertEqual('mode: canonical', params_str[2])
        self.assertEqual('nnz weights: 802920', params_str[3])
        self.assertEqual('avg weight: 3.68754', params_str[4])

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

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[representation],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 11', params_str[0])
        self.assertEqual('nodes (k): 802920', params_str[1])
        self.assertEqual('mode: canonical', params_str[2])
        self.assertEqual('nnz weights: 802920', params_str[3])
        self.assertEqual('avg weight: 3.68754', params_str[4])

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

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[representation],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 4', params_str[0])
        self.assertEqual('nodes (k): 256', params_str[1])
        self.assertEqual('mode: basic', params_str[2])
        self.assertEqual('nnz weights: 256', params_str[3])
        self.assertEqual('avg weight: {}'.format(avg_count_expected), params_str[4])

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

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[representation],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: {}'.format(k), params_str[0])
        self.assertEqual('nodes (k): 2', params_str[1])
        self.assertEqual('mode: basic', params_str[2])
        self.assertEqual('nnz weights: 2', params_str[3])
        self.assertEqual('avg weight: {}'.format(avg_count_expected), params_str[4])


if __name__ == '__main__':
    unittest.main()
