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

        stats = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[representation])
        self.assertEqual('20', stats['k'])
        self.assertEqual('591997', stats['nodes (k)'])
        self.assertEqual('basic', stats['mode'])
        self.assertEqual('591997', stats['nnz weights'])
        self.assertEqual('2.48587', stats['avg weight'])

        clean_fasta = self.tempdir.name + '/contigs.fasta.gz'
        self._clean(self.tempdir.name + '/graph' + graph_file_extension[representation],
                    output=clean_fasta,
                    extra_params='')  # no cleaning

        self._build_graph(input=clean_fasta,
                          output=self.tempdir.name + '/graph_clean',
                          k=20, repr=representation,
                          extra_params="--mask-dummy --count-kmers")

        stats = self._get_stats(self.tempdir.name + '/graph_clean' + graph_file_extension[representation])
        self.assertEqual('20', stats['k'])
        self.assertEqual('591997', stats['nodes (k)'])
        self.assertEqual('basic', stats['mode'])
        self.assertEqual('591997', stats['nnz weights'])
        self.assertEqual('2.48587', stats['avg weight'])

    @parameterized.expand([repr for repr in GRAPH_TYPES if not (repr == 'bitmap' and PROTEIN_MODE)])
    def test_no_cleaning_contigs_2bit_counts(self, representation):

        self._build_graph(input=TEST_DATA_DIR + '/transcripts_1000.fa',
                          output=self.tempdir.name + '/graph',
                          k=20, repr=representation,
                          extra_params="--mask-dummy --count-kmers --count-width 2")

        stats = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[representation])
        self.assertEqual('20', stats['k'])
        self.assertEqual('591997', stats['nodes (k)'])
        self.assertEqual('basic', stats['mode'])
        self.assertEqual('591997', stats['nnz weights'])
        self.assertEqual('1.73589', stats['avg weight'])

        clean_fasta = self.tempdir.name + '/contigs.fasta.gz'
        self._clean(self.tempdir.name + '/graph' + graph_file_extension[representation],
                    output=clean_fasta,
                    extra_params='')  # no cleaning

        self._build_graph(input=clean_fasta,
                          output=self.tempdir.name + '/graph_clean',
                          k=20, repr=representation,
                          extra_params="--mask-dummy --count-kmers")

        stats = self._get_stats(self.tempdir.name + '/graph_clean' + graph_file_extension[representation])
        self.assertEqual('20', stats['k'])
        self.assertEqual('591997', stats['nodes (k)'])
        self.assertEqual('basic', stats['mode'])
        self.assertEqual('591997', stats['nnz weights'])
        self.assertEqual('1.73589', stats['avg weight'])

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

        stats = self._get_stats(self.tempdir.name + '/graph_clean' + graph_file_extension[representation])
        self.assertEqual('20', stats['k'])
        self.assertEqual('589774', stats['nodes (k)'])
        self.assertEqual('basic', stats['mode'])

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

        stats = self._get_stats(self.tempdir.name + '/graph_clean' + graph_file_extension[representation])
        self.assertEqual('20', stats['k'])
        self.assertEqual('589774', stats['nodes (k)'])
        self.assertEqual('basic', stats['mode'])
        self.assertEqual('589774', stats['nnz weights'])
        self.assertEqual('2.49001', stats['avg weight'])

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

        stats = self._get_stats(self.tempdir.name + '/graph_clean' + graph_file_extension[representation])
        self.assertEqual('20', stats['k'])
        self.assertEqual('167395', stats['nodes (k)'])
        self.assertEqual('basic', stats['mode'])
        self.assertEqual('167395', stats['nnz weights'])
        self.assertEqual('5.52732', stats['avg weight'])
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

        stats = self._get_stats(self.tempdir.name + '/graph_clean' + graph_file_extension[representation])
        self.assertEqual('20', stats['k'])
        self.assertEqual('167224', stats['nodes (k)'])
        self.assertEqual('basic', stats['mode'])
        self.assertEqual('167224', stats['nnz weights'])
        self.assertEqual('5.52757', stats['avg weight'])


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

        stats = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[representation])
        self.assertEqual('31', stats['k'])
        self.assertEqual('1185814', stats['nodes (k)'])
        self.assertEqual('canonical', stats['mode'])
        self.assertEqual('1185814', stats['nnz weights'])
        self.assertEqual('2.4635', stats['avg weight'])

        clean_fasta = self.tempdir.name + '/contigs.fasta.gz'
        self._clean(self.tempdir.name + '/graph' + graph_file_extension[representation],
                    output=clean_fasta,
                    extra_params='')  # no cleaning

        self._build_graph(input=clean_fasta,
                          output=self.tempdir.name + '/graph_clean',
                          k=31, repr=representation, mode='canonical',
                          extra_params="--mask-dummy --count-kmers")

        stats = self._get_stats(self.tempdir.name + '/graph_clean' + graph_file_extension[representation])
        self.assertEqual('31', stats['k'])
        self.assertEqual('1185814', stats['nodes (k)'])
        self.assertEqual('canonical', stats['mode'])
        self.assertEqual('1185814', stats['nnz weights'])
        self.assertEqual('2.4635', stats['avg weight'])

    # TODO: add 'hashstr' once the canonical mode is implemented for it
    @parameterized.expand(['succinct', 'bitmap', 'hash'])  # , 'hashstr']:
    def test_no_cleaning_contigs_2bit_counts(self, representation):

        self._build_graph(input=TEST_DATA_DIR + '/transcripts_1000.fa',
                          output=self.tempdir.name + '/graph',
                          k=31, repr=representation, mode='canonical',
                          extra_params="--mask-dummy --count-kmers --count-width 2")

        stats = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[representation])
        self.assertEqual('31', stats['k'])
        self.assertEqual('1185814', stats['nodes (k)'])
        self.assertEqual('canonical', stats['mode'])
        self.assertEqual('1185814', stats['nnz weights'])
        self.assertEqual('1.72792', stats['avg weight'])

        clean_fasta = self.tempdir.name + '/contigs.fasta.gz'
        self._clean(self.tempdir.name + '/graph' + graph_file_extension[representation],
                    output=clean_fasta,
                    extra_params='')  # no cleaning

        self._build_graph(input=clean_fasta,
                          output=self.tempdir.name + '/graph_clean',
                          k=31, repr=representation, mode='canonical',
                          extra_params="--mask-dummy --count-kmers")

        stats = self._get_stats(self.tempdir.name + '/graph_clean' + graph_file_extension[representation])
        self.assertEqual('31', stats['k'])
        self.assertEqual('1185814', stats['nodes (k)'])
        self.assertEqual('canonical', stats['mode'])
        self.assertEqual('1185814', stats['nnz weights'])
        self.assertEqual('1.72792', stats['avg weight'])

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

        stats = self._get_stats(self.tempdir.name + '/graph_clean' + graph_file_extension[representation])
        self.assertEqual('31', stats['k'])
        self.assertEqual('1180802', stats['nodes (k)'])
        self.assertEqual('canonical', stats['mode'])

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

        stats = self._get_stats(self.tempdir.name + '/graph_clean' + graph_file_extension[representation])
        self.assertEqual('31', stats['k'])
        self.assertEqual('1180802', stats['nodes (k)'])
        self.assertEqual('canonical', stats['mode'])
        self.assertEqual('1180802', stats['nnz weights'])
        self.assertEqual('2.46882', stats['avg weight'])

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

        stats = self._get_stats(self.tempdir.name + '/graph_clean' + graph_file_extension[representation])
        self.assertEqual('31', stats['k'])
        self.assertEqual('331452', stats['nodes (k)'])
        self.assertEqual('basic', stats['mode'])
        self.assertEqual('331452', stats['nnz weights'])
        self.assertEqual('5.52692', stats['avg weight'])

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

        stats = self._get_stats(self.tempdir.name + '/graph_clean' + graph_file_extension[representation])
        self.assertEqual('31', stats['k'])
        self.assertEqual('331452', stats['nodes (k)'])
        self.assertEqual('canonical', stats['mode'])
        self.assertEqual('331452', stats['nnz weights'])
        self.assertEqual('5.52692', stats['avg weight'])

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

        stats = self._get_stats(self.tempdir.name + '/graph_clean' + graph_file_extension[representation])
        self.assertEqual('31', stats['k'])
        self.assertEqual('331266', stats['nodes (k)'])
        self.assertEqual('canonical', stats['mode'])
        self.assertEqual('331266', stats['nnz weights'])
        self.assertEqual('5.52728', stats['avg weight'])


if __name__ == '__main__':
    unittest.main()
