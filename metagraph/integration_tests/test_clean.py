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

GRAPH_TYPES = [graph_type for graph_type, _ in graph_file_extension.items()]


class TestCleanWeighted(unittest.TestCase):
    def setUp(self):
        self.tempdir = TemporaryDirectory()

    @parameterized.expand([repr for repr in GRAPH_TYPES if not (repr == 'bitmap' and PROTEIN_MODE)])
    def test_no_cleaning_contigs(self, representation):

        construct_command = '{exe} build --mask-dummy \
                --graph {repr} -k 20 --count-kmers -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=representation,
            outfile=self.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )

        res = subprocess.run([construct_command], shell=True, stdout=PIPE)
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

        clean_command = '{exe} clean \
                --to-fasta -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=representation,
            outfile=self.tempdir.name + '/contigs.fasta.gz',
            input=self.tempdir.name + '/graph' + graph_file_extension[representation]
        )

        res = subprocess.run([clean_command], shell=True)
        self.assertEqual(res.returncode, 0)

        reconstruct_command = '{exe} build --mask-dummy \
                --graph {repr} -k 20 --count-kmers -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=representation,
            outfile=self.tempdir.name + '/graph_clean',
            input=self.tempdir.name + '/contigs.fasta.gz'
        )

        res = subprocess.run([reconstruct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph_clean' + graph_file_extension[representation],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 20', params_str[0])
        self.assertEqual('nodes (k): 591997', params_str[1])
        self.assertEqual('mode: basic', params_str[2])
        self.assertEqual('nnz weights: 591997', params_str[3])
        self.assertEqual('avg weight: 2.48587', params_str[4])

    @parameterized.expand([repr for repr in GRAPH_TYPES if not (repr == 'bitmap' and PROTEIN_MODE)])
    def test_no_cleaning_contigs_2bit_counts(self, representation):

        construct_command = '{exe} build --mask-dummy \
                --graph {repr} -k 20 --count-kmers --count-width 2 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=representation,
            outfile=self.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )

        res = subprocess.run([construct_command], shell=True, stdout=PIPE)
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
        self.assertEqual('avg weight: 1.73589', params_str[4])

        clean_command = '{exe} clean \
                --to-fasta -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=representation,
            outfile=self.tempdir.name + '/contigs.fasta.gz',
            input=self.tempdir.name + '/graph' + graph_file_extension[representation]
        )

        res = subprocess.run([clean_command], shell=True)
        self.assertEqual(res.returncode, 0)

        reconstruct_command = '{exe} build --mask-dummy \
                --graph {repr} -k 20 --count-kmers -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=representation,
            outfile=self.tempdir.name + '/graph_clean',
            input=self.tempdir.name + '/contigs.fasta.gz'
        )

        res = subprocess.run([reconstruct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph_clean' + graph_file_extension[representation],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 20', params_str[0])
        self.assertEqual('nodes (k): 591997', params_str[1])
        self.assertEqual('mode: basic', params_str[2])
        self.assertEqual('nnz weights: 591997', params_str[3])
        self.assertEqual('avg weight: 1.73589', params_str[4])


@unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets")
class TestCleanWeightedCanonical(unittest.TestCase):
    def setUp(self):
        self.tempdir = TemporaryDirectory()

    # TODO: add 'hashstr' once the canonical mode is implemented for it
    @parameterized.expand(['succinct', 'bitmap', 'hash'])  # , 'hashstr']:
    def test_no_cleaning_contigs(self, representation):

        construct_command = '{exe} build --mask-dummy \
                --graph {repr} -k 31 --mode canonical --count-kmers -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=representation,
            outfile=self.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )

        res = subprocess.run([construct_command], shell=True, stdout=PIPE)
        self.assertEqual(res.returncode, 0)

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[representation],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 31', params_str[0])
        self.assertEqual('nodes (k): 1185814', params_str[1])
        self.assertEqual('mode: canonical', params_str[2])
        self.assertEqual('nnz weights: 1185814', params_str[3])
        self.assertEqual('avg weight: 2.4635', params_str[4])

        clean_command = '{exe} clean \
                --to-fasta -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=representation,
            outfile=self.tempdir.name + '/contigs.fasta.gz',
            input=self.tempdir.name + '/graph' + graph_file_extension[representation]
        )

        res = subprocess.run([clean_command], shell=True)
        self.assertEqual(res.returncode, 0)

        reconstruct_command = '{exe} build --mask-dummy \
                --graph {repr} -k 31 --mode canonical --count-kmers -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=representation,
            outfile=self.tempdir.name + '/graph_clean',
            input=self.tempdir.name + '/contigs.fasta.gz'
        )

        res = subprocess.run([reconstruct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph_clean' + graph_file_extension[representation],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 31', params_str[0])
        self.assertEqual('nodes (k): 1185814', params_str[1])
        self.assertEqual('mode: canonical', params_str[2])
        self.assertEqual('nnz weights: 1185814', params_str[3])
        self.assertEqual('avg weight: 2.4635', params_str[4])

    # TODO: add 'hashstr' once the canonical mode is implemented for it
    @parameterized.expand(['succinct', 'bitmap', 'hash'])  # , 'hashstr']:
    def test_no_cleaning_contigs_2bit_counts(self, representation):

        construct_command = '{exe} build --mask-dummy \
                --graph {repr} -k 31 --mode canonical --count-kmers --count-width 2 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=representation,
            outfile=self.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )

        res = subprocess.run([construct_command], shell=True, stdout=PIPE)
        self.assertEqual(res.returncode, 0)

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[representation],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 31', params_str[0])
        self.assertEqual('nodes (k): 1185814', params_str[1])
        self.assertEqual('mode: canonical', params_str[2])
        self.assertEqual('nnz weights: 1185814', params_str[3])
        self.assertEqual('avg weight: 1.72792', params_str[4])

        clean_command = '{exe} clean \
                --to-fasta -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=representation,
            outfile=self.tempdir.name + '/contigs.fasta.gz',
            input=self.tempdir.name + '/graph' + graph_file_extension[representation]
        )

        res = subprocess.run([clean_command], shell=True)
        self.assertEqual(res.returncode, 0)

        reconstruct_command = '{exe} build --mask-dummy \
                --graph {repr} -k 31 --mode canonical --count-kmers -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=representation,
            outfile=self.tempdir.name + '/graph_clean',
            input=self.tempdir.name + '/contigs.fasta.gz'
        )

        res = subprocess.run([reconstruct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph_clean' + graph_file_extension[representation],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 31', params_str[0])
        self.assertEqual('nodes (k): 1185814', params_str[1])
        self.assertEqual('mode: canonical', params_str[2])
        self.assertEqual('nnz weights: 1185814', params_str[3])
        self.assertEqual('avg weight: 1.72792', params_str[4])


if __name__ == '__main__':
    unittest.main()
