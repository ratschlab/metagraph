import unittest
from parameterized import parameterized
import itertools
import subprocess
from subprocess import PIPE
from tempfile import TemporaryDirectory
import glob
import os
import gzip
from base import TestingBase, METAGRAPH, TEST_DATA_DIR, NUM_THREADS


"""Test graph construction"""

PROTEIN_MODE = os.readlink(METAGRAPH).endswith("_Protein")

graph_file_extension = {'succinct': '.dbg',
                        'bitmap': '.bitmapdbg',
                        'hash': '.orhashdbg',
                        'hashfast': '.hashfastdbg',
                        'hashstr': '.hashstrdbg'}

GRAPH_TYPES = [graph_type for graph_type, _ in graph_file_extension.items()]


class TestCleanWeighted(TestingBase):
    def setUp(self):
        super().setUpClass()

    @parameterized.expand([repr for repr in GRAPH_TYPES if not (repr == 'bitmap' and PROTEIN_MODE)])
    def test_no_cleaning_contigs(self, representation):

        self._build_graph(input=TEST_DATA_DIR + '/transcripts_1000.fa',
                          output=self.tempdir.name + '/graph',
                          k=20, repr=representation,
                          extra_params="--mask-dummy --count-kmers")

        res = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[representation])
        out = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 20', out[0])
        self.assertEqual('nodes (k): 591997', out[1])
        self.assertEqual('mode: basic', out[2])
        self.assertEqual('nnz weights: 591997', out[3])
        self.assertEqual('avg weight: 2.48587', out[4])

        clean_fasta = self.tempdir.name + '/contigs.fasta.gz'
        self._clean(self.tempdir.name + '/graph' + graph_file_extension[representation],
                    output=clean_fasta,
                    extra_params='')  # no cleaning

        self._build_graph(input=clean_fasta,
                          output=self.tempdir.name + '/graph_clean',
                          k=20, repr=representation,
                          extra_params="--mask-dummy --count-kmers")

        res = self._get_stats(self.tempdir.name + '/graph_clean' + graph_file_extension[representation])
        out = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 20', out[0])
        self.assertEqual('nodes (k): 591997', out[1])
        self.assertEqual('mode: basic', out[2])
        self.assertEqual('nnz weights: 591997', out[3])
        self.assertEqual('avg weight: 2.48587', out[4])

    @parameterized.expand([repr for repr in GRAPH_TYPES if not (repr == 'bitmap' and PROTEIN_MODE)])
    def test_no_cleaning_contigs_2bit_counts(self, representation):

        self._build_graph(input=TEST_DATA_DIR + '/transcripts_1000.fa',
                          output=self.tempdir.name + '/graph',
                          k=20, repr=representation,
                          extra_params="--mask-dummy --count-kmers --count-width 2")

        res = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[representation])
        out = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 20', out[0])
        self.assertEqual('nodes (k): 591997', out[1])
        self.assertEqual('mode: basic', out[2])
        self.assertEqual('nnz weights: 591997', out[3])
        self.assertEqual('avg weight: 1.73589', out[4])

        clean_fasta = self.tempdir.name + '/contigs.fasta.gz'
        self._clean(self.tempdir.name + '/graph' + graph_file_extension[representation],
                    output=clean_fasta,
                    extra_params='')  # no cleaning

        self._build_graph(input=clean_fasta,
                          output=self.tempdir.name + '/graph_clean',
                          k=20, repr=representation,
                          extra_params="--mask-dummy --count-kmers")

        res = self._get_stats(self.tempdir.name + '/graph_clean' + graph_file_extension[representation])
        out = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 20', out[0])
        self.assertEqual('nodes (k): 591997', out[1])
        self.assertEqual('mode: basic', out[2])
        self.assertEqual('nnz weights: 591997', out[3])
        self.assertEqual('avg weight: 1.73589', out[4])

    @parameterized.expand([repr for repr in GRAPH_TYPES if not (repr == 'bitmap' and PROTEIN_MODE)])
    def test_clean_prune_tips_no_counts(self, representation):

        self._build_graph(input=TEST_DATA_DIR + '/transcripts_1000.fa',
                          output=self.tempdir.name + '/graph',
                          k=20, repr=representation,
                          extra_params="--mask-dummy")

        clean_fasta = self.tempdir.name + '/contigs.fasta.gz'
        self._clean(self.tempdir.name + '/graph' + graph_file_extension[representation],
                    output=clean_fasta,
                    extra_params='--prune-tips 60')

        self._build_graph(input=clean_fasta,
                          output=self.tempdir.name + '/graph_clean',
                          k=20, repr=representation,
                          extra_params="--mask-dummy")

        res = self._get_stats(self.tempdir.name + '/graph_clean' + graph_file_extension[representation])
        out = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 20', out[0])
        self.assertEqual('nodes (k): 589774', out[1])
        self.assertEqual('mode: basic', out[2])

    @parameterized.expand([repr for repr in GRAPH_TYPES if not (repr == 'bitmap' and PROTEIN_MODE)])
    def test_clean_prune_tips(self, representation):

        self._build_graph(input=TEST_DATA_DIR + '/transcripts_1000.fa',
                          output=self.tempdir.name + '/graph',
                          k=20, repr=representation,
                          extra_params="--mask-dummy --count-kmers")

        clean_fasta = self.tempdir.name + '/contigs.fasta.gz'
        self._clean(self.tempdir.name + '/graph' + graph_file_extension[representation],
                    output=clean_fasta,
                    extra_params='--prune-tips 60')

        self._build_graph(input=clean_fasta,
                          output=self.tempdir.name + '/graph_clean',
                          k=20, repr=representation,
                          extra_params="--mask-dummy --count-kmers")

        res = self._get_stats(self.tempdir.name + '/graph_clean' + graph_file_extension[representation])
        out = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 20', out[0])
        self.assertEqual('nodes (k): 589774', out[1])
        self.assertEqual('mode: basic', out[2])
        self.assertEqual('nnz weights: 589774', out[3])
        self.assertEqual('avg weight: 2.49001', out[4])

    @parameterized.expand([repr for repr in GRAPH_TYPES if not (repr == 'bitmap' and PROTEIN_MODE)])
    def test_cleaning_threshold_fixed(self, representation):

        self._build_graph(input=TEST_DATA_DIR + '/transcripts_1000.fa',
                          output=self.tempdir.name + '/graph',
                          k=20, repr=representation,
                          extra_params="--mask-dummy --count-kmers")

        clean_fasta = self.tempdir.name + '/contigs.fasta.gz'
        self._clean(self.tempdir.name + '/graph' + graph_file_extension[representation],
                    output=clean_fasta,
                    extra_params='--prune-unitigs 3')

        self._build_graph(input=clean_fasta,
                          output=self.tempdir.name + '/graph_clean',
                          k=20, repr=representation,
                          extra_params="--mask-dummy --count-kmers")

        res = self._get_stats(self.tempdir.name + '/graph_clean' + graph_file_extension[representation])
        out = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 20', out[0])
        self.assertEqual('nodes (k): 167395', out[1])
        self.assertEqual('mode: basic', out[2])
        self.assertEqual('nnz weights: 167395', out[3])
        self.assertEqual('avg weight: 5.52732', out[4])

    @parameterized.expand([repr for repr in GRAPH_TYPES if not (repr == 'bitmap' and PROTEIN_MODE)])
    def test_cleaning_prune_tips_threshold_fixed(self, representation):

        self._build_graph(input=TEST_DATA_DIR + '/transcripts_1000.fa',
                          output=self.tempdir.name + '/graph',
                          k=20, repr=representation,
                          extra_params="--mask-dummy --count-kmers")

        clean_fasta = self.tempdir.name + '/contigs.fasta.gz'
        self._clean(self.tempdir.name + '/graph' + graph_file_extension[representation],
                    output=clean_fasta,
                    extra_params='--prune-tips 60 --prune-unitigs 3')

        self._build_graph(input=clean_fasta,
                          output=self.tempdir.name + '/graph_clean',
                          k=20, repr=representation,
                          extra_params="--mask-dummy --count-kmers")

        res = self._get_stats(self.tempdir.name + '/graph_clean' + graph_file_extension[representation])
        out = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 20', out[0])
        self.assertEqual('nodes (k): 167224', out[1])
        self.assertEqual('mode: basic', out[2])
        self.assertEqual('nnz weights: 167224', out[3])
        self.assertEqual('avg weight: 5.52757', out[4])


@unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets")
class TestCleanWeightedCanonical(TestingBase):
    def setUp(self):
        super().setUpClass()

    # TODO: add 'hashstr' once the canonical mode is implemented for it
    @parameterized.expand(['succinct', 'bitmap', 'hash'])  # , 'hashstr']:
    def test_no_cleaning_contigs(self, representation):

        self._build_graph(input=TEST_DATA_DIR + '/transcripts_1000.fa',
                          output=self.tempdir.name + '/graph',
                          k=31, repr=representation, mode='canonical',
                          extra_params="--mask-dummy --count-kmers")

        res = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[representation])
        out = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 31', out[0])
        self.assertEqual('nodes (k): 1185814', out[1])
        self.assertEqual('mode: canonical', out[2])
        self.assertEqual('nnz weights: 1185814', out[3])
        self.assertEqual('avg weight: 2.4635', out[4])

        clean_fasta = self.tempdir.name + '/contigs.fasta.gz'
        self._clean(self.tempdir.name + '/graph' + graph_file_extension[representation],
                    output=clean_fasta,
                    extra_params='')  # no cleaning

        self._build_graph(input=clean_fasta,
                          output=self.tempdir.name + '/graph_clean',
                          k=31, repr=representation, mode='canonical',
                          extra_params="--mask-dummy --count-kmers")

        res = self._get_stats(self.tempdir.name + '/graph_clean' + graph_file_extension[representation])
        out = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 31', out[0])
        self.assertEqual('nodes (k): 1185814', out[1])
        self.assertEqual('mode: canonical', out[2])
        self.assertEqual('nnz weights: 1185814', out[3])
        self.assertEqual('avg weight: 2.4635', out[4])

    # TODO: add 'hashstr' once the canonical mode is implemented for it
    @parameterized.expand(['succinct', 'bitmap', 'hash'])  # , 'hashstr']:
    def test_no_cleaning_contigs_2bit_counts(self, representation):

        self._build_graph(input=TEST_DATA_DIR + '/transcripts_1000.fa',
                          output=self.tempdir.name + '/graph',
                          k=31, repr=representation, mode='canonical',
                          extra_params="--mask-dummy --count-kmers --count-width 2")

        res = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[representation])
        out = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 31', out[0])
        self.assertEqual('nodes (k): 1185814', out[1])
        self.assertEqual('mode: canonical', out[2])
        self.assertEqual('nnz weights: 1185814', out[3])
        self.assertEqual('avg weight: 1.72792', out[4])

        clean_fasta = self.tempdir.name + '/contigs.fasta.gz'
        self._clean(self.tempdir.name + '/graph' + graph_file_extension[representation],
                    output=clean_fasta,
                    extra_params='')  # no cleaning

        self._build_graph(input=clean_fasta,
                          output=self.tempdir.name + '/graph_clean',
                          k=31, repr=representation, mode='canonical',
                          extra_params="--mask-dummy --count-kmers")

        res = self._get_stats(self.tempdir.name + '/graph_clean' + graph_file_extension[representation])
        out = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 31', out[0])
        self.assertEqual('nodes (k): 1185814', out[1])
        self.assertEqual('mode: canonical', out[2])
        self.assertEqual('nnz weights: 1185814', out[3])
        self.assertEqual('avg weight: 1.72792', out[4])

    @parameterized.expand(['succinct', 'bitmap', 'hash'])  # , 'hashstr']:
    def test_clean_prune_tips_no_counts(self, representation):

        self._build_graph(input=TEST_DATA_DIR + '/transcripts_1000.fa',
                          output=self.tempdir.name + '/graph',
                          k=31, repr=representation, mode='canonical',
                          extra_params="--mask-dummy")

        clean_fasta = self.tempdir.name + '/contigs.fasta.gz'
        self._clean(self.tempdir.name + '/graph' + graph_file_extension[representation],
                    output=clean_fasta,
                    extra_params='--prune-tips 60')

        self._build_graph(input=clean_fasta,
                          output=self.tempdir.name + '/graph_clean',
                          k=31, repr=representation, mode='canonical',
                          extra_params="--mask-dummy")

        res = self._get_stats(self.tempdir.name + '/graph_clean' + graph_file_extension[representation])
        out = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 31', out[0])
        self.assertEqual('nodes (k): 1180802', out[1])
        self.assertEqual('mode: canonical', out[2])

    @parameterized.expand(['succinct', 'bitmap', 'hash'])  # , 'hashstr']:
    def test_clean_prune_tips(self, representation):

        self._build_graph(input=TEST_DATA_DIR + '/transcripts_1000.fa',
                          output=self.tempdir.name + '/graph',
                          k=31, repr=representation, mode='canonical',
                          extra_params="--mask-dummy --count-kmers")

        clean_fasta = self.tempdir.name + '/contigs.fasta.gz'
        self._clean(self.tempdir.name + '/graph' + graph_file_extension[representation],
                    output=clean_fasta,
                    extra_params='--prune-tips 60')

        self._build_graph(input=clean_fasta,
                          output=self.tempdir.name + '/graph_clean',
                          k=31, repr=representation, mode='canonical',
                          extra_params="--mask-dummy --count-kmers")

        res = self._get_stats(self.tempdir.name + '/graph_clean' + graph_file_extension[representation])
        out = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 31', out[0])
        self.assertEqual('nodes (k): 1180802', out[1])
        self.assertEqual('mode: canonical', out[2])
        self.assertEqual('nnz weights: 1180802', out[3])
        self.assertEqual('avg weight: 2.46882', out[4])

    @parameterized.expand(GRAPH_TYPES)
    def test_cleaning_threshold_fixed_both_strands(self, representation):

        self._build_graph(input=TEST_DATA_DIR + '/transcripts_1000.fa',
                          output=self.tempdir.name + '/graph',
                          k=31, repr=representation,
                          extra_params="--mask-dummy --count-kmers --fwd-and-reverse")

        # extract all unitigs, not only the primary ones
        clean_fasta = self.tempdir.name + '/contigs.fasta.gz'
        self._clean(self.tempdir.name + '/graph' + graph_file_extension[representation],
                    output=clean_fasta,
                    extra_params='--prune-unitigs 3')

        self._build_graph(input=clean_fasta,
                          output=self.tempdir.name + '/graph_clean',
                          k=31, repr=representation,
                          extra_params="--mask-dummy --count-kmers")

        res = self._get_stats(self.tempdir.name + '/graph_clean' + graph_file_extension[representation])
        out = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 31', out[0])
        self.assertEqual('nodes (k): 331452', out[1])
        self.assertEqual('mode: basic', out[2])
        self.assertEqual('nnz weights: 331452', out[3])
        self.assertEqual('avg weight: 5.52692', out[4])

    @parameterized.expand(['succinct', 'bitmap', 'hash'])  # , 'hashstr']:
    def test_cleaning_threshold_fixed(self, representation):

        self._build_graph(input=TEST_DATA_DIR + '/transcripts_1000.fa',
                          output=self.tempdir.name + '/graph',
                          k=31, repr=representation, mode='canonical',
                          extra_params="--mask-dummy --count-kmers")

        clean_fasta = self.tempdir.name + '/contigs.fasta.gz'
        self._clean(self.tempdir.name + '/graph' + graph_file_extension[representation],
                    output=clean_fasta,
                    extra_params='--prune-unitigs 3')

        self._build_graph(input=clean_fasta,
                          output=self.tempdir.name + '/graph_clean',
                          k=31, repr=representation, mode='canonical',
                          extra_params="--mask-dummy --count-kmers")

        res = self._get_stats(self.tempdir.name + '/graph_clean' + graph_file_extension[representation])
        out = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 31', out[0])
        self.assertEqual('nodes (k): 331452', out[1])
        self.assertEqual('mode: canonical', out[2])
        self.assertEqual('nnz weights: 331452', out[3])
        self.assertEqual('avg weight: 5.52692', out[4])

    @parameterized.expand(['succinct', 'bitmap', 'hash'])  # , 'hashstr']:
    def test_cleaning_prune_tips_threshold_fixed(self, representation):

        self._build_graph(input=TEST_DATA_DIR + '/transcripts_1000.fa',
                          output=self.tempdir.name + '/graph',
                          k=31, repr=representation, mode='canonical',
                          extra_params="--mask-dummy --count-kmers")

        clean_fasta = self.tempdir.name + '/contigs.fasta.gz'
        self._clean(self.tempdir.name + '/graph' + graph_file_extension[representation],
                    output=clean_fasta,
                    extra_params='--prune-tips 60 --prune-unitigs 3')

        self._build_graph(input=clean_fasta,
                          output=self.tempdir.name + '/graph_clean',
                          k=31, repr=representation, mode='canonical',
                          extra_params="--mask-dummy --count-kmers")

        res = self._get_stats(self.tempdir.name + '/graph_clean' + graph_file_extension[representation])
        out = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 31', out[0])
        self.assertEqual('nodes (k): 331266', out[1])
        self.assertEqual('mode: canonical', out[2])
        self.assertEqual('nnz weights: 331266', out[3])
        self.assertEqual('avg weight: 5.52728', out[4])


if __name__ == '__main__':
    unittest.main()
