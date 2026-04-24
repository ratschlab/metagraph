import unittest
from parameterized import parameterized
import subprocess
from subprocess import PIPE
from tempfile import TemporaryDirectory
import glob
import os
import re

from base import PROTEIN_MODE, DNA_MODE, TestingBase, METAGRAPH, TEST_DATA_DIR, NUM_THREADS, MMAP_FLAG


"""Test graph construction and alignment"""

graph_file_extension = {'succinct': '.dbg',
                        'bitmap': '.bitmapdbg',
                        'hash': '.orhashdbg',
                        'hashfast': '.hashfastdbg'}

GRAPH_TYPES = [graph_type for graph_type, _ in graph_file_extension.items()]


@unittest.skipUnless(DNA_MODE, "These alignment tests are only for the DNA4 alphabet")
class TestDNAAlign(TestingBase):
    def setUp(self):
        super().setUpClass()

    @parameterized.expand(GRAPH_TYPES)
    def test_simple_align_all_graphs(self, representation):

        self._build_graph(input=TEST_DATA_DIR + '/genome.MT.fa',
                          output=self.tempdir.name + '/genome.MT',
                          k=11, repr=representation,
                          extra_params="--mask-dummy")

        params = self._get_stats(self.tempdir.name + '/genome.MT' + graph_file_extension[representation])
        self.assertEqual('11', params['k'])
        self.assertEqual('16438', params['nodes (k)'])
        self.assertEqual('basic', params['mode'])

        stats_command = '{exe} align --align-only-forwards -i {graph} --align-min-exact-match 0.0 {reads}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/genome.MT' + graph_file_extension[representation],
            reads=TEST_DATA_DIR + '/genome_MT1.fq',
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().rstrip().split('\n')
        self.assertEqual(len(params_str), 7)
        self.assertEqual(params_str[0], 'MT-10/1\tAACAGAGAATAGTTTAAATTAGAATCTTAGCTTTGGGTGCTAATGGTGGAGTTAAAGACTTTTTCTCTGATTTGTCCTTGGAAAAAGGTTTTCATCTCCGGTTTACAAGACTGGTGTATTAGTTTATACTACAAGGACAGGCCCATTTGA\t+\tTAGAATCTTAG\t22\t11\t19S11=120S\t0')
        self.assertEqual(params_str[1], 'MT-8/1\tAAAACTAACCCCCTAATAAAATTAATTAACCACTCATTCATCGACCTCCCCACCCCATCCAACATCTCCGCATGATGAAACTTCGGCTCACTCCTTGGCGCCTGCCTGATCCTCCAAATCACCACAGGACTATTCCTAGCCATGCACTAC\t+\tAAAACTAACCCCCTAATAAAATTAATTAACCACTCATTCATCGACCTCCCCACCCCATCCAACATCTCCGCATGATGAAACTTCGGCTCACTCCTTGGCGCCTGCCTGATCCTCCAAATCACCACAGGACTATTCCTAGCCATGCACTAC\t310\t150\t150=\t0')
        self.assertEqual(params_str[2], 'MT-6/1\tATATGACTAGCTTACACAATAGCTTTTATAGTAAAGATACCTCTTTACGGACTCCACTTATGACTCCCTAAAGCCCATGTCGAAGCCCCCATCGCTGGGTCAATAGTACTTGCCGCAGTACTCTTAAAACTAGGCGGCTATGGTATAATA\t+\tATATGACTAGCTTACACAATAGCTTTTATAGTAAAGATACCTCTTTACGGACTCCACTTATGACTCCCTAAAGCCCATGTCGAAGCCCCCATCGCTGGGTCAATAGTACTTGCCGCAGTACTCTTAAAACTAGGCGGCTATGGTATAATA\t310\t150\t150=\t0')
        self.assertEqual(params_str[3], 'MT-4/1\tAGTATAGTAGTTCGCTTTGACTGGTGAAGTCTTAGCATGTACTGCTCGGAGGTTCGGTTCTGCTCCGAGGTCGCCCCAACCGAAATTTTTAATGCAGGTTTGGTAGTTTAGGACCTGTGGGTTTGTTAGGTACTGTTTGCATTAATAAAT\t*\t*\t0\t*\t*\t*')
        self.assertEqual(params_str[4], 'MT-2/1\tTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCACTTTCCACACAGACATCATAACAAAAAATTTCCACCAAACCCCCCCTCCCCCGCTTCTGGCCACAGCACTTAAACACATCTCTGCCAAAC\t+\tTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCACTTTCCACACAGACATCATAACAAAAAATTTCCACCAAACCCCCCCTCCCCCGCTTCTGGCCACAGCACTTAAACACATCTCTGCCAAAC\t310\t150\t150=\t0')
        last_split = params_str[5].split("\t")
        self.assertEqual(last_split[0], "MT-11/1")
        self.assertEqual(last_split[1], "AACAGAGAATTGTTTAAATTACAATCTTAGCTATGGGTGCTAAAGGTGGAGTTATAGACTTTTTCACTGATTTGTCGTTGGAAAAAGCTTTTCATCTCGGGTTTACAAGTCTGGTGTATTTGTTTATACTAGAAGGACAGGCGCATTTGA")
        self.assertEqual(last_split[4], "22")

    @parameterized.expand(GRAPH_TYPES)
    def test_simple_align_map_all_graphs(self, representation):

        self._build_graph(input=TEST_DATA_DIR + '/genome.MT.fa',
                          output=self.tempdir.name + '/genome.MT',
                          k=11, repr=representation,
                          extra_params="--mask-dummy")

        params = self._get_stats(self.tempdir.name + '/genome.MT' + graph_file_extension[representation])
        self.assertEqual('11', params['k'])
        self.assertEqual('16438', params['nodes (k)'])
        self.assertEqual('basic', params['mode'])

        stats_command = '{exe} align -i {graph} --map --count-kmers {reads}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/genome.MT' + graph_file_extension[representation],
            reads=TEST_DATA_DIR + '/genome_MT1.fq',
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().rstrip().split('\n')
        self.assertEqual(len(params_str), 7)
        self.assertEqual(params_str[0], 'MT-10/1\t1/140/1')
        self.assertEqual(params_str[1], 'MT-8/1\t140/140/140')
        self.assertEqual(params_str[2], 'MT-6/1\t140/140/140')
        self.assertEqual(params_str[3], 'MT-4/1\t0/140/0')
        self.assertEqual(params_str[4], 'MT-2/1\t140/140/140')
        self.assertEqual(params_str[5], 'MT-11/1\t1/140/1')
        self.assertEqual(params_str[6], 'MT-11/1\t1/140/1')

    @parameterized.expand(GRAPH_TYPES)
    def test_simple_align_map_all_graphs_subk(self, representation):

        self._build_graph(input=TEST_DATA_DIR + '/genome.MT.fa',
                          output=self.tempdir.name + '/genome.MT',
                          k=11, repr=representation,
                          extra_params="--mask-dummy")

        params = self._get_stats(self.tempdir.name + '/genome.MT' + graph_file_extension[representation])
        self.assertEqual('11', params['k'])
        self.assertEqual('16438', params['nodes (k)'])
        self.assertEqual('basic', params['mode'])

        stats_command = '{exe} align -i {graph} --map --count-kmers --align-length 10 {reads}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/genome.MT' + graph_file_extension[representation],
            reads=TEST_DATA_DIR + '/genome_MT1.fq',
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE, stderr=PIPE)
        if representation != 'succinct':
            self.assertNotEqual(res.returncode, 0)
            return

        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().rstrip().split('\n')
        self.assertEqual(len(params_str), 7)
        self.assertEqual(params_str[0], 'MT-10/1\t3/141/3')
        self.assertEqual(params_str[1], 'MT-8/1\t141/141/141')
        self.assertEqual(params_str[2], 'MT-6/1\t141/141/141')
        self.assertEqual(params_str[3], 'MT-4/1\t1/141/1')
        self.assertEqual(params_str[4], 'MT-2/1\t141/141/141')
        self.assertEqual(params_str[5], 'MT-11/1\t4/141/4')
        self.assertEqual(params_str[6], 'MT-11/1\t3/141/3')

    @parameterized.expand(GRAPH_TYPES)
    def test_simple_align_map_canonical_all_graphs(self, representation):

        self._build_graph(input=TEST_DATA_DIR + '/genome.MT.fa',
                          output=self.tempdir.name + '/genome.MT',
                          k=11, repr=representation, mode='canonical',
                          extra_params="--mask-dummy")

        params = self._get_stats(self.tempdir.name + '/genome.MT' + graph_file_extension[representation])
        self.assertEqual('11', params['k'])
        self.assertEqual('32782', params['nodes (k)'])
        self.assertEqual('canonical', params['mode'])

        stats_command = '{exe} align -i {graph} --map --count-kmers {reads}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/genome.MT' + graph_file_extension[representation],
            reads=TEST_DATA_DIR + '/genome_MT1.fq',
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().rstrip().split('\n')
        self.assertEqual(len(params_str), 7)
        self.assertEqual(params_str[0], 'MT-10/1\t140/140/140')
        self.assertEqual(params_str[1], 'MT-8/1\t140/140/140')
        self.assertEqual(params_str[2], 'MT-6/1\t140/140/140')
        self.assertEqual(params_str[3], 'MT-4/1\t129/140/129')
        self.assertEqual(params_str[4], 'MT-2/1\t140/140/139')
        self.assertEqual(params_str[5], 'MT-11/1\t2/140/2')
        self.assertEqual(params_str[6], 'MT-11/1\t140/140/140')

    @parameterized.expand(['succinct'])
    def test_simple_align_json_all_graphs(self, representation):

        self._build_graph(input=TEST_DATA_DIR + '/genome.MT.fa',
                          output=self.tempdir.name + '/genome.MT',
                          k=11, repr=representation,
                          extra_params="--mask-dummy")

        params = self._get_stats(self.tempdir.name + '/genome.MT' + graph_file_extension[representation])
        self.assertEqual('11', params['k'])
        self.assertEqual('16438', params['nodes (k)'])
        self.assertEqual('basic', params['mode'])

        stats_command = '{exe} align --align-only-forwards -i {graph} --align-min-exact-match 0.0 {reads}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/genome.MT' + graph_file_extension[representation],
            reads=TEST_DATA_DIR + '/genome_MT1.fq',
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().rstrip().split('\n')
        self.assertEqual(len(params_str), 7)

    @parameterized.expand(GRAPH_TYPES)
    def test_simple_align_fwd_rev_comp_all_graphs(self, representation):

        self._build_graph(input=TEST_DATA_DIR + '/genome.MT.fa',
                          output=self.tempdir.name + '/genome.MT',
                          k=11, repr=representation,
                          extra_params="--mask-dummy")

        params = self._get_stats(self.tempdir.name + '/genome.MT' + graph_file_extension[representation])
        self.assertEqual('11', params['k'])
        self.assertEqual('16438', params['nodes (k)'])
        self.assertEqual('basic', params['mode'])

        stats_command = '{exe} align -i {graph} --align-min-exact-match 0.0 {reads}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/genome.MT' + graph_file_extension[representation],
            reads=TEST_DATA_DIR + '/genome_MT1.fq',
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().rstrip().split('\n')
        self.assertEqual(len(params_str), 7)
        self.assertEqual(params_str[0], 'MT-10/1\tAACAGAGAATAGTTTAAATTAGAATCTTAGCTTTGGGTGCTAATGGTGGAGTTAAAGACTTTTTCTCTGATTTGTCCTTGGAAAAAGGTTTTCATCTCCGGTTTACAAGACTGGTGTATTAGTTTATACTACAAGGACAGGCCCATTTGA\t-\tTCAAATGGGCCTGTCCTTGTAGTATAAACTAATACACCAGTCTTGTAAACCGGAGATGAAAACCTTTTTCCAAGGACAAATCAGAGAAAAAGTCTTTAACTCCACCATTAGCACCCAAAGCTAAGATTCTAATTTAAACTATTCTCTGTT\t310\t150\t150=\t0')
        self.assertEqual(params_str[1], 'MT-8/1\tAAAACTAACCCCCTAATAAAATTAATTAACCACTCATTCATCGACCTCCCCACCCCATCCAACATCTCCGCATGATGAAACTTCGGCTCACTCCTTGGCGCCTGCCTGATCCTCCAAATCACCACAGGACTATTCCTAGCCATGCACTAC\t+\tAAAACTAACCCCCTAATAAAATTAATTAACCACTCATTCATCGACCTCCCCACCCCATCCAACATCTCCGCATGATGAAACTTCGGCTCACTCCTTGGCGCCTGCCTGATCCTCCAAATCACCACAGGACTATTCCTAGCCATGCACTAC\t310\t150\t150=\t0')
        self.assertEqual(params_str[2], 'MT-6/1\tATATGACTAGCTTACACAATAGCTTTTATAGTAAAGATACCTCTTTACGGACTCCACTTATGACTCCCTAAAGCCCATGTCGAAGCCCCCATCGCTGGGTCAATAGTACTTGCCGCAGTACTCTTAAAACTAGGCGGCTATGGTATAATA\t+\tATATGACTAGCTTACACAATAGCTTTTATAGTAAAGATACCTCTTTACGGACTCCACTTATGACTCCCTAAAGCCCATGTCGAAGCCCCCATCGCTGGGTCAATAGTACTTGCCGCAGTACTCTTAAAACTAGGCGGCTATGGTATAATA\t310\t150\t150=\t0')
        self.assertEqual(params_str[3], 'MT-4/1\tAGTATAGTAGTTCGCTTTGACTGGTGAAGTCTTAGCATGTACTGCTCGGAGGTTCGGTTCTGCTCCGAGGTCGCCCCAACCGAAATTTTTAATGCAGGTTTGGTAGTTTAGGACCTGTGGGTTTGTTAGGTACTGTTTGCATTAATAAAT\t-\tATTTATTAATGCAAACAGTACCTAACAAACCCACAGGTCCTAAACTACCAAACCTGCATTAAAAATTTCGGTTGGGGCGACCTCGGAGCAGAACCCAACCTCCGAGCAGTACATGCTAAGACTTCACCAGTCAAAGCGAACTACTATACT\t305\t149\t95=1X54=\t0')
        self.assertEqual(params_str[4], 'MT-2/1\tTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCACTTTCCACACAGACATCATAACAAAAAATTTCCACCAAACCCCCCCTCCCCCGCTTCTGGCCACAGCACTTAAACACATCTCTGCCAAAC\t+\tTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCACTTTCCACACAGACATCATAACAAAAAATTTCCACCAAACCCCCCCTCCCCCGCTTCTGGCCACAGCACTTAAACACATCTCTGCCAAAC\t310\t150\t150=\t0')
        last_split = params_str[5].split("\t")
        self.assertEqual(last_split[0], "MT-11/1")
        self.assertEqual(last_split[1], "AACAGAGAATTGTTTAAATTACAATCTTAGCTATGGGTGCTAAAGGTGGAGTTATAGACTTTTTCACTGATTTGTCGTTGGAAAAAGCTTTTCATCTCGGGTTTACAAGTCTGGTGTATTTGTTTATACTAGAAGGACAGGCGCATTTGA")
        self.assertEqual(last_split[4], "22")

    @parameterized.expand(GRAPH_TYPES)
    def test_simple_align_canonical_all_graphs(self, representation):

        self._build_graph(input=TEST_DATA_DIR + '/genome.MT.fa',
                          output=self.tempdir.name + '/genome.MT',
                          k=11, repr=representation, mode='canonical',
                          extra_params="--mask-dummy")

        params = self._get_stats(self.tempdir.name + '/genome.MT' + graph_file_extension[representation])
        self.assertEqual('11', params['k'])
        self.assertEqual('32782', params['nodes (k)'])
        self.assertEqual('canonical', params['mode'])

        stats_command = '{exe} align -i {graph} --align-min-exact-match 0.0 {reads}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/genome.MT' + graph_file_extension[representation],
            reads=TEST_DATA_DIR + '/genome_MT1.fq',
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().rstrip().split('\n')
        self.maxDiff = None
        self.assertEqual(len(params_str), 7)
        self.assertEqual(params_str[0], 'MT-10/1\tAACAGAGAATAGTTTAAATTAGAATCTTAGCTTTGGGTGCTAATGGTGGAGTTAAAGACTTTTTCTCTGATTTGTCCTTGGAAAAAGGTTTTCATCTCCGGTTTACAAGACTGGTGTATTAGTTTATACTACAAGGACAGGCCCATTTGA\t+\tAACAGAGAATAGTTTAAATTAGAATCTTAGCTTTGGGTGCTAATGGTGGAGTTAAAGACTTTTTCTCTGATTTGTCCTTGGAAAAAGGTTTTCATCTCCGGTTTACAAGACTGGTGTATTAGTTTATACTACAAGGACAGGCCCATTTGA\t310\t150\t150=\t0')
        self.assertEqual(params_str[1], 'MT-8/1\tAAAACTAACCCCCTAATAAAATTAATTAACCACTCATTCATCGACCTCCCCACCCCATCCAACATCTCCGCATGATGAAACTTCGGCTCACTCCTTGGCGCCTGCCTGATCCTCCAAATCACCACAGGACTATTCCTAGCCATGCACTAC\t+\tAAAACTAACCCCCTAATAAAATTAATTAACCACTCATTCATCGACCTCCCCACCCCATCCAACATCTCCGCATGATGAAACTTCGGCTCACTCCTTGGCGCCTGCCTGATCCTCCAAATCACCACAGGACTATTCCTAGCCATGCACTAC\t310\t150\t150=\t0')
        self.assertEqual(params_str[2], 'MT-6/1\tATATGACTAGCTTACACAATAGCTTTTATAGTAAAGATACCTCTTTACGGACTCCACTTATGACTCCCTAAAGCCCATGTCGAAGCCCCCATCGCTGGGTCAATAGTACTTGCCGCAGTACTCTTAAAACTAGGCGGCTATGGTATAATA\t+\tATATGACTAGCTTACACAATAGCTTTTATAGTAAAGATACCTCTTTACGGACTCCACTTATGACTCCCTAAAGCCCATGTCGAAGCCCCCATCGCTGGGTCAATAGTACTTGCCGCAGTACTCTTAAAACTAGGCGGCTATGGTATAATA\t310\t150\t150=\t0')
        self.assertEqual(params_str[3], 'MT-4/1\tAGTATAGTAGTTCGCTTTGACTGGTGAAGTCTTAGCATGTACTGCTCGGAGGTTCGGTTCTGCTCCGAGGTCGCCCCAACCGAAATTTTTAATGCAGGTTTGGTAGTTTAGGACCTGTGGGTTTGTTAGGTACTGTTTGCATTAATAAAT\t+\tAGTATAGTAGTTCGCTTTGACTGGTGAAGTCTTAGCATGTACTGCTCGGAGGTTGGGTTCTGCTCCGAGGTCGCCCCAACCGAAATTTTTAATGCAGGTTTGGTAGTTTAGGACCTGTGGGTTTGTTAGGTACTGTTTGCATTAATAAAT\t305\t149\t54=1X95=\t0')
        self.assertEqual(params_str[4], 'MT-2/1\tTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCACTTTCCACACAGACATCATAACAAAAAATTTCCACCAAACCCCCCCTCCCCCGCTTCTGGCCACAGCACTTAAACACATCTCTGCCAAAC\t+\tTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCACTTTCCACACAGACATCATAACAAAAAATTTCCACCAAACCCCCCCTCCCCCGCTTCTGGCCACAGCACTTAAACACATCTCTGCCAAAC\t310\t150\t150=\t0')
        last_split = params_str[5].split("\t")
        self.assertEqual(last_split[0], "MT-11/1")
        self.assertEqual(last_split[1], "AACAGAGAATTGTTTAAATTACAATCTTAGCTATGGGTGCTAAAGGTGGAGTTATAGACTTTTTCACTGATTTGTCGTTGGAAAAAGCTTTTCATCTCGGGTTTACAAGTCTGGTGTATTTGTTTATACTAGAAGGACAGGCGCATTTGA")
        self.assertEqual(last_split[4], "22")

    @parameterized.expand(['succinct'])
    def test_simple_align_canonical_subk_succinct(self, representation):

        self._build_graph(input=TEST_DATA_DIR + '/genome.MT.fa',
                          output=self.tempdir.name + '/genome.MT',
                          k=11, repr=representation, mode='canonical',
                          extra_params="--mask-dummy")

        params = self._get_stats(self.tempdir.name + '/genome.MT' + graph_file_extension[representation])
        self.assertEqual('11', params['k'])
        self.assertEqual('32782', params['nodes (k)'])
        self.assertEqual('canonical', params['mode'])

        stats_command = '{exe} align -i {graph} --align-min-exact-match 0.0 --align-min-seed-length 10 {reads}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/genome.MT' + graph_file_extension[representation],
            reads=TEST_DATA_DIR + '/genome_MT1.fq',
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().rstrip().split('\n')
        self.assertEqual(len(params_str), 7)
        self.assertEqual(params_str[0], 'MT-10/1\tAACAGAGAATAGTTTAAATTAGAATCTTAGCTTTGGGTGCTAATGGTGGAGTTAAAGACTTTTTCTCTGATTTGTCCTTGGAAAAAGGTTTTCATCTCCGGTTTACAAGACTGGTGTATTAGTTTATACTACAAGGACAGGCCCATTTGA\t+\tAACAGAGAATAGTTTAAATTAGAATCTTAGCTTTGGGTGCTAATGGTGGAGTTAAAGACTTTTTCTCTGATTTGTCCTTGGAAAAAGGTTTTCATCTCCGGTTTACAAGACTGGTGTATTAGTTTATACTACAAGGACAGGCCCATTTGA\t310\t150\t150=\t0')
        self.assertEqual(params_str[1], 'MT-8/1\tAAAACTAACCCCCTAATAAAATTAATTAACCACTCATTCATCGACCTCCCCACCCCATCCAACATCTCCGCATGATGAAACTTCGGCTCACTCCTTGGCGCCTGCCTGATCCTCCAAATCACCACAGGACTATTCCTAGCCATGCACTAC\t+\tAAAACTAACCCCCTAATAAAATTAATTAACCACTCATTCATCGACCTCCCCACCCCATCCAACATCTCCGCATGATGAAACTTCGGCTCACTCCTTGGCGCCTGCCTGATCCTCCAAATCACCACAGGACTATTCCTAGCCATGCACTAC\t310\t150\t150=\t0')
        self.assertEqual(params_str[2], 'MT-6/1\tATATGACTAGCTTACACAATAGCTTTTATAGTAAAGATACCTCTTTACGGACTCCACTTATGACTCCCTAAAGCCCATGTCGAAGCCCCCATCGCTGGGTCAATAGTACTTGCCGCAGTACTCTTAAAACTAGGCGGCTATGGTATAATA\t+\tATATGACTAGCTTACACAATAGCTTTTATAGTAAAGATACCTCTTTACGGACTCCACTTATGACTCCCTAAAGCCCATGTCGAAGCCCCCATCGCTGGGTCAATAGTACTTGCCGCAGTACTCTTAAAACTAGGCGGCTATGGTATAATA\t310\t150\t150=\t0')
        self.assertEqual(params_str[3], 'MT-4/1\tAGTATAGTAGTTCGCTTTGACTGGTGAAGTCTTAGCATGTACTGCTCGGAGGTTCGGTTCTGCTCCGAGGTCGCCCCAACCGAAATTTTTAATGCAGGTTTGGTAGTTTAGGACCTGTGGGTTTGTTAGGTACTGTTTGCATTAATAAAT\t+\tAGTATAGTAGTTCGCTTTGACTGGTGAAGTCTTAGCATGTACTGCTCGGAGGTTGGGTTCTGCTCCGAGGTCGCCCCAACCGAAATTTTTAATGCAGGTTTGGTAGTTTAGGACCTGTGGGTTTGTTAGGTACTGTTTGCATTAATAAAT\t305\t149\t54=1X95=\t0')
        self.assertEqual(params_str[4], 'MT-2/1\tTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCACTTTCCACACAGACATCATAACAAAAAATTTCCACCAAACCCCCCCTCCCCCGCTTCTGGCCACAGCACTTAAACACATCTCTGCCAAAC\t+\tTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCACTTTCCACACAGACATCATAACAAAAAATTTCCACCAAACCCCCCCTCCCCCGCTTCTGGCCACAGCACTTAAACACATCTCTGCCAAAC\t310\t150\t150=\t0')
        self.assertEqual(params_str[5], 'MT-11/1\tAACAGAGAATTGTTTAAATTACAATCTTAGCTATGGGTGCTAAAGGTGGAGTTATAGACTTTTTCACTGATTTGTCGTTGGAAAAAGCTTTTCATCTCGGGTTTACAAGTCTGGTGTATTTGTTTATACTAGAAGGACAGGCGCATTTGA\t+\tAACAGAGAATAGTTTAAATTAGAATCTTAGCTTTGGGTGCTAATGGTGGAGTTAAAGACTTTTTCTCTGATTTGTCCTTGGAAAAAGGTTTTCATCTCCGGTTTACAAGACTGGTGTATTAGTTTATACTACAAGGACAGGCCCATTTGA\t245\t137\t10=1X10=1X10=1X10=1X10=1X10=1X10=1X10=1X10=1X10=1X10=1X10=1X10=1X7=\t0')

    @parameterized.expand(GRAPH_TYPES)
    def test_simple_align_primary_all_graphs(self, representation):

        self._build_graph(input=TEST_DATA_DIR + '/genome.MT.fa',
                          output=self.tempdir.name + '/genome.MT.primary',
                          k=11, repr=representation, mode='primary',
                          extra_params="--mask-dummy")

        params = self._get_stats(self.tempdir.name + '/genome.MT.primary' + graph_file_extension[representation])
        self.assertEqual('11', params['k'])
        self.assertEqual('16391', params['nodes (k)'])
        self.assertEqual('primary', params['mode'])

        stats_command = '{exe} align -i {graph} --align-min-exact-match 0.0 {reads}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/genome.MT.primary' + graph_file_extension[representation],
            reads=TEST_DATA_DIR + '/genome_MT1.fq',
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().rstrip().split('\n')
        self.assertEqual(len(params_str), 7)
        self.assertEqual(params_str[0], 'MT-10/1\tAACAGAGAATAGTTTAAATTAGAATCTTAGCTTTGGGTGCTAATGGTGGAGTTAAAGACTTTTTCTCTGATTTGTCCTTGGAAAAAGGTTTTCATCTCCGGTTTACAAGACTGGTGTATTAGTTTATACTACAAGGACAGGCCCATTTGA\t+\tAACAGAGAATAGTTTAAATTAGAATCTTAGCTTTGGGTGCTAATGGTGGAGTTAAAGACTTTTTCTCTGATTTGTCCTTGGAAAAAGGTTTTCATCTCCGGTTTACAAGACTGGTGTATTAGTTTATACTACAAGGACAGGCCCATTTGA\t310\t150\t150=\t0')
        self.assertEqual(params_str[1], 'MT-8/1\tAAAACTAACCCCCTAATAAAATTAATTAACCACTCATTCATCGACCTCCCCACCCCATCCAACATCTCCGCATGATGAAACTTCGGCTCACTCCTTGGCGCCTGCCTGATCCTCCAAATCACCACAGGACTATTCCTAGCCATGCACTAC\t+\tAAAACTAACCCCCTAATAAAATTAATTAACCACTCATTCATCGACCTCCCCACCCCATCCAACATCTCCGCATGATGAAACTTCGGCTCACTCCTTGGCGCCTGCCTGATCCTCCAAATCACCACAGGACTATTCCTAGCCATGCACTAC\t310\t150\t150=\t0')
        self.assertEqual(params_str[2], 'MT-6/1\tATATGACTAGCTTACACAATAGCTTTTATAGTAAAGATACCTCTTTACGGACTCCACTTATGACTCCCTAAAGCCCATGTCGAAGCCCCCATCGCTGGGTCAATAGTACTTGCCGCAGTACTCTTAAAACTAGGCGGCTATGGTATAATA\t+\tATATGACTAGCTTACACAATAGCTTTTATAGTAAAGATACCTCTTTACGGACTCCACTTATGACTCCCTAAAGCCCATGTCGAAGCCCCCATCGCTGGGTCAATAGTACTTGCCGCAGTACTCTTAAAACTAGGCGGCTATGGTATAATA\t310\t150\t150=\t0')
        self.assertEqual(params_str[3], 'MT-4/1\tAGTATAGTAGTTCGCTTTGACTGGTGAAGTCTTAGCATGTACTGCTCGGAGGTTCGGTTCTGCTCCGAGGTCGCCCCAACCGAAATTTTTAATGCAGGTTTGGTAGTTTAGGACCTGTGGGTTTGTTAGGTACTGTTTGCATTAATAAAT\t+\tAGTATAGTAGTTCGCTTTGACTGGTGAAGTCTTAGCATGTACTGCTCGGAGGTTGGGTTCTGCTCCGAGGTCGCCCCAACCGAAATTTTTAATGCAGGTTTGGTAGTTTAGGACCTGTGGGTTTGTTAGGTACTGTTTGCATTAATAAAT\t305\t149\t54=1X95=\t0')
        self.assertEqual(params_str[4], 'MT-2/1\tTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCACTTTCCACACAGACATCATAACAAAAAATTTCCACCAAACCCCCCCTCCCCCGCTTCTGGCCACAGCACTTAAACACATCTCTGCCAAAC\t+\tTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCACTTTCCACACAGACATCATAACAAAAAATTTCCACCAAACCCCCCCTCCCCCGCTTCTGGCCACAGCACTTAAACACATCTCTGCCAAAC\t310\t150\t150=\t0')
        self.assertEqual(params_str[6].split("\t")[4], "310")
        last_split = params_str[5].split("\t")
        self.assertEqual(last_split[0], "MT-11/1")
        self.assertEqual(last_split[1], "AACAGAGAATTGTTTAAATTACAATCTTAGCTATGGGTGCTAAAGGTGGAGTTATAGACTTTTTCACTGATTTGTCGTTGGAAAAAGCTTTTCATCTCGGGTTTACAAGTCTGGTGTATTTGTTTATACTAGAAGGACAGGCGCATTTGA")
        self.assertEqual(last_split[4], "22")

    @parameterized.expand(['succinct'])
    def test_simple_align_primary_subk_succinct(self, representation):

        self._build_graph(input=TEST_DATA_DIR + '/genome.MT.fa',
                          output=self.tempdir.name + '/genome.MT.primary',
                          k=11, repr=representation, mode='primary',
                          extra_params="--mask-dummy")

        params = self._get_stats(self.tempdir.name + '/genome.MT.primary' + graph_file_extension[representation])
        self.assertEqual('11', params['k'])
        self.assertEqual('16391', params['nodes (k)'])
        self.assertEqual('primary', params['mode'])

        stats_command = '{exe} align -i {graph} --align-min-exact-match 0.0 --align-min-seed-length 10 {reads}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/genome.MT.primary' + graph_file_extension[representation],
            reads=TEST_DATA_DIR + '/genome_MT1.fq',
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().rstrip().split('\n')
        self.assertEqual(len(params_str), 7)
        self.assertEqual(params_str[0], 'MT-10/1\tAACAGAGAATAGTTTAAATTAGAATCTTAGCTTTGGGTGCTAATGGTGGAGTTAAAGACTTTTTCTCTGATTTGTCCTTGGAAAAAGGTTTTCATCTCCGGTTTACAAGACTGGTGTATTAGTTTATACTACAAGGACAGGCCCATTTGA\t+\tAACAGAGAATAGTTTAAATTAGAATCTTAGCTTTGGGTGCTAATGGTGGAGTTAAAGACTTTTTCTCTGATTTGTCCTTGGAAAAAGGTTTTCATCTCCGGTTTACAAGACTGGTGTATTAGTTTATACTACAAGGACAGGCCCATTTGA\t310\t150\t150=\t0')
        self.assertEqual(params_str[1], 'MT-8/1\tAAAACTAACCCCCTAATAAAATTAATTAACCACTCATTCATCGACCTCCCCACCCCATCCAACATCTCCGCATGATGAAACTTCGGCTCACTCCTTGGCGCCTGCCTGATCCTCCAAATCACCACAGGACTATTCCTAGCCATGCACTAC\t+\tAAAACTAACCCCCTAATAAAATTAATTAACCACTCATTCATCGACCTCCCCACCCCATCCAACATCTCCGCATGATGAAACTTCGGCTCACTCCTTGGCGCCTGCCTGATCCTCCAAATCACCACAGGACTATTCCTAGCCATGCACTAC\t310\t150\t150=\t0')
        self.assertEqual(params_str[2], 'MT-6/1\tATATGACTAGCTTACACAATAGCTTTTATAGTAAAGATACCTCTTTACGGACTCCACTTATGACTCCCTAAAGCCCATGTCGAAGCCCCCATCGCTGGGTCAATAGTACTTGCCGCAGTACTCTTAAAACTAGGCGGCTATGGTATAATA\t+\tATATGACTAGCTTACACAATAGCTTTTATAGTAAAGATACCTCTTTACGGACTCCACTTATGACTCCCTAAAGCCCATGTCGAAGCCCCCATCGCTGGGTCAATAGTACTTGCCGCAGTACTCTTAAAACTAGGCGGCTATGGTATAATA\t310\t150\t150=\t0')
        self.assertEqual(params_str[3], 'MT-4/1\tAGTATAGTAGTTCGCTTTGACTGGTGAAGTCTTAGCATGTACTGCTCGGAGGTTCGGTTCTGCTCCGAGGTCGCCCCAACCGAAATTTTTAATGCAGGTTTGGTAGTTTAGGACCTGTGGGTTTGTTAGGTACTGTTTGCATTAATAAAT\t+\tAGTATAGTAGTTCGCTTTGACTGGTGAAGTCTTAGCATGTACTGCTCGGAGGTTGGGTTCTGCTCCGAGGTCGCCCCAACCGAAATTTTTAATGCAGGTTTGGTAGTTTAGGACCTGTGGGTTTGTTAGGTACTGTTTGCATTAATAAAT\t305\t149\t54=1X95=\t0')
        self.assertEqual(params_str[4], 'MT-2/1\tTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCACTTTCCACACAGACATCATAACAAAAAATTTCCACCAAACCCCCCCTCCCCCGCTTCTGGCCACAGCACTTAAACACATCTCTGCCAAAC\t+\tTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCACTTTCCACACAGACATCATAACAAAAAATTTCCACCAAACCCCCCCTCCCCCGCTTCTGGCCACAGCACTTAAACACATCTCTGCCAAAC\t310\t150\t150=\t0')
        self.assertEqual(params_str[5], 'MT-11/1\tAACAGAGAATTGTTTAAATTACAATCTTAGCTATGGGTGCTAAAGGTGGAGTTATAGACTTTTTCACTGATTTGTCGTTGGAAAAAGCTTTTCATCTCGGGTTTACAAGTCTGGTGTATTTGTTTATACTAGAAGGACAGGCGCATTTGA\t+\tAACAGAGAATAGTTTAAATTAGAATCTTAGCTTTGGGTGCTAATGGTGGAGTTAAAGACTTTTTCTCTGATTTGTCCTTGGAAAAAGGTTTTCATCTCCGGTTTACAAGACTGGTGTATTAGTTTATACTACAAGGACAGGCCCATTTGA\t245\t137\t10=1X10=1X10=1X10=1X10=1X10=1X10=1X10=1X10=1X10=1X10=1X10=1X10=1X7=\t0')

    @parameterized.expand(['succinct'])
    def test_simple_align_fwd_rev_comp_json_all_graphs(self, representation):

        self._build_graph(input=TEST_DATA_DIR + '/genome.MT.fa',
                          output=self.tempdir.name + '/genome.MT',
                          k=11, repr=representation,
                          extra_params="--mask-dummy")

        params = self._get_stats(self.tempdir.name + '/genome.MT' + graph_file_extension[representation])
        self.assertEqual('11', params['k'])
        self.assertEqual('16438', params['nodes (k)'])
        self.assertEqual('basic', params['mode'])

        stats_command = '{exe} align --json -i {graph} --align-min-exact-match 0.0 {reads}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/genome.MT' + graph_file_extension[representation],
            reads=TEST_DATA_DIR + '/genome_MT1.fq',
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().rstrip().split('\n')
        self.assertEqual(len(params_str), 7)
        ref_align_str = [a.rstrip() for a in open(TEST_DATA_DIR + '/genome_MT1.align.json', 'r').readlines()]
        for [a, b] in zip(params_str, ref_align_str):
            self.assertEqual(a, b)

    @parameterized.expand(['succinct'])
    def test_simple_align_edit_distance_all_graphs(self, representation):

        self._build_graph(input=TEST_DATA_DIR + '/genome.MT.fa',
                          output=self.tempdir.name + '/genome.MT',
                          k=11, repr=representation,
                          extra_params="--mask-dummy")

        params = self._get_stats(self.tempdir.name + '/genome.MT' + graph_file_extension[representation])
        self.assertEqual('11', params['k'])
        self.assertEqual('16438', params['nodes (k)'])
        self.assertEqual('basic', params['mode'])

        stats_command = '{exe} align --json --align-edit-distance -i {graph} --align-min-exact-match 0.0 {reads}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/genome.MT' + graph_file_extension[representation],
            reads=TEST_DATA_DIR + '/genome_MT1.fq',
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().rstrip().split('\n')
        self.assertEqual(len(params_str), 7)
        ref_align_str = [a.rstrip() for a in open(TEST_DATA_DIR + '/genome_MT1.align.edit.json', 'r').readlines()]
        for [a, b] in zip(params_str, ref_align_str):
            self.assertEqual(a, b)


@unittest.skipUnless(PROTEIN_MODE, "These alignment tests are only for the Protein alphabet")
class TestProteinAlign(TestingBase):
    def setUp(self):
        super().setUpClass()

    # TODO: test alignment for protein sequences
    # @parameterized.expand(GRAPH_TYPES)
    # def test_simple_align_all_graphs(self, representation):
    #     pass


coord_anno_file_extension = {
    'column_coord': '.column_coord.annodbg',
    'row_diff_brwt_coord': '.row_diff_brwt_coord.annodbg',
}

COORD_ANNO_TYPES = list(coord_anno_file_extension.keys())


@unittest.skipUnless(DNA_MODE, "These alignment tests are only for the DNA4 alphabet")
class TestAlignCoordToHeader(TestingBase):
    """Test that 'metagraph align' resolves coordinates via .seqs index."""

    def setUp(self):
        self.tempdir = TemporaryDirectory()

    def tearDown(self):
        if getattr(self, 'tempdir', None) is not None:
            self.tempdir.cleanup()

    def _setup_graph(self, fa_files, anno_repr, k=5):
        """Build graph + coord annotation + .seqs index. Return (graph, anno) paths."""
        if isinstance(fa_files, str):
            fa_files = [fa_files]
        fa_input = ' '.join(fa_files)

        graph_base = self.tempdir.name + '/graph'
        graph = graph_base + '.dbg'
        anno_base = self.tempdir.name + '/annotation'
        anno = anno_base + coord_anno_file_extension[anno_repr]

        self._build_graph(fa_input, graph_base, k=k, repr='succinct', mode='basic')
        self._annotate_graph(fa_input, graph, anno_base, anno_repr, anno_type='filename')

        # The annotation column order may differ from the input file order
        # (e.g. due to parallel construction).  Query the actual column order
        # so the .seqs index matches.
        res = subprocess.run([f"{METAGRAPH} stats {anno} --print-col-names"],
                             shell=True, stdout=PIPE, stderr=PIPE)
        self.assertEqual(res.returncode, 0)
        columns = [l.strip() for l in res.stdout.decode().split('\n')[2:] if l.strip()]
        seqs_input = ' '.join(columns) if columns else fa_input

        cmd = (f"{METAGRAPH} annotate --anno-filename --index-header-coords -v "
               f"-i {graph} -o {anno_base} {seqs_input}" + MMAP_FLAG)
        res = subprocess.run([cmd], shell=True, stdout=PIPE, stderr=PIPE)
        self.assertEqual(res.returncode, 0,
                         f"Indexing headers failed: {res.stderr.decode()}")
        self.assertTrue(os.path.exists(anno_base + '.seqs'))
        return graph, anno

    def _run_align(self, graph, anno, query_fa, extra_flags='', only_forwards=True):
        """Run align and return parsed output lines as lists of tab-separated fields."""
        forwards = '--align-only-forwards ' if only_forwards else ''
        align_cmd = (f'{METAGRAPH} align {forwards}{extra_flags} '
                     f'-i {graph} -a {anno} {query_fa}' + MMAP_FLAG)
        res = subprocess.run([align_cmd], shell=True, stdout=PIPE, stderr=PIPE)
        self.assertEqual(res.returncode, 0, f"Align failed: {res.stderr.decode()}")
        output = res.stdout.decode().rstrip()
        if not output:
            return []
        return [line.split('\t') for line in output.split('\n')]

    def _write_fa(self, name, sequences):
        """Write a FASTA file. sequences is a list of (header, seq) tuples."""
        path = os.path.join(self.tempdir.name, name)
        with open(path, 'w') as f:
            for header, seq in sequences:
                f.write(f'>{header}\n{seq}\n')
        return path

    @parameterized.expand(COORD_ANNO_TYPES)
    def test_align_with_seqs_maps_coords_to_headers(self, anno_repr):
        """Core test: .seqs resolves global coords to per-sequence header:start-end."""
        #   seq1: GTATCGATCG       (10 bp)
        #   seq2: GCTAGCTAGCTAGCTA (16 bp)
        #   seq3: ATCGATCGAAAAACCCCCGGGGGTTTTT (28 bp)
        test_fa = self._write_fa('seqs.fa', [
            ('seq1', 'GTATCGATCG'),
            ('seq2', 'GCTAGCTAGCTAGCTA'),
            ('seq3', 'ATCGATCGAAAAACCCCCGGGGGTTTTT'),
        ])
        query_fa = self._write_fa('query.fa', [
            ('query1', 'TATCGATCG'),     # matches seq1
            ('query2', 'GCTAGCTAGCTAG'), # matches seq2
            ('query3', 'AAAAACCCCC'),    # matches seq3
        ])

        graph, anno = self._setup_graph(test_fa, anno_repr)
        rows = self._run_align(graph, anno, query_fa)
        self.assertEqual(len(rows), 3)

        # query1 -> seq1:2-10
        self.assertEqual(rows[0][0], 'query1')
        self.assertEqual(rows[0][1], 'TATCGATCG')
        self.assertEqual(rows[0][3], 'TATCGATCG')
        self.assertEqual(rows[0][6], '9=')
        self.assertEqual(rows[0][7], '0')
        self.assertEqual(rows[0][8], 'seq1:2-10')

        # query2 -> seq2:1-13
        self.assertEqual(rows[1][0], 'query2')
        self.assertEqual(rows[1][1], 'GCTAGCTAGCTAG')
        self.assertEqual(rows[1][3], 'GCTAGCTAGCTAG')
        self.assertEqual(rows[1][6], '13=')
        self.assertEqual(rows[1][7], '0')
        self.assertEqual(rows[1][8], 'seq2:1-13')

        # query3 -> seq3:9-18
        self.assertEqual(rows[2][0], 'query3')
        self.assertEqual(rows[2][1], 'AAAAACCCCC')
        self.assertEqual(rows[2][3], 'AAAAACCCCC')
        self.assertEqual(rows[2][6], '10=')
        self.assertEqual(rows[2][7], '0')
        self.assertEqual(rows[2][8], 'seq3:9-18')

    @parameterized.expand(COORD_ANNO_TYPES)
    def test_align_no_coord_mapping_flag(self, anno_repr):
        """--no-coord-mapping disables header resolution, keeps global file offsets."""
        test_fa = self._write_fa('seqs.fa', [
            ('seq1', 'GTATCGATCG'),
            ('seq2', 'GCTAGCTAGCTAGCTA'),
            ('seq3', 'ATCGATCGAAAAACCCCCGGGGGTTTTT'),
        ])
        query_fa = self._write_fa('query.fa', [
            ('query1', 'TATCGATCG'),
            ('query2', 'GCTAGCTAGCTAG'),
            ('query3', 'AAAAACCCCC'),
        ])

        graph, anno = self._setup_graph(test_fa, anno_repr)
        rows = self._run_align(graph, anno, query_fa, extra_flags='--no-coord-mapping')
        self.assertEqual(len(rows), 3)

        # With --no-coord-mapping, labels use global file path:offset format
        self.assertEqual(rows[0][0], 'query1')
        self.assertEqual(rows[0][6], '9=')
        self.assertEqual(rows[0][8], f'{test_fa}:2-10')

        self.assertEqual(rows[1][0], 'query2')
        self.assertEqual(rows[1][6], '13=')
        self.assertEqual(rows[1][8], f'{test_fa}:7-19')

        self.assertEqual(rows[2][0], 'query3')
        self.assertEqual(rows[2][6], '10=')
        self.assertEqual(rows[2][8], f'{test_fa}:27-36')

    @parameterized.expand(COORD_ANNO_TYPES)
    def test_align_shared_kmers_multiple_labels(self, anno_repr):
        """Query matching k-mers shared between sequences produces semicolon-separated labels."""
        test_fa = self._write_fa('seqs.fa', [
            ('seq1', 'GTATCGATCG'),
            ('seq2', 'GCTAGCTAGCTAGCTA'),
            ('seq3', 'ATCGATCGAAAAACCCCCGGGGGTTTTT'),
        ])
        # ATCGATCG appears in both seq1 (pos 3-10) and seq3 (pos 1-8)
        query_fa = self._write_fa('query.fa', [('shared_query', 'ATCGATCG')])

        graph, anno = self._setup_graph(test_fa, anno_repr)
        rows = self._run_align(graph, anno, query_fa)
        self.assertEqual(len(rows), 1)

        self.assertEqual(rows[0][0], 'shared_query')
        self.assertEqual(rows[0][1], 'ATCGATCG')
        self.assertEqual(rows[0][3], 'ATCGATCG')
        self.assertEqual(rows[0][6], '8=')
        self.assertEqual(rows[0][7], '0')
        # Both sequences reported, separated by ;
        self.assertEqual(rows[0][8], 'seq1:3-10;seq3:1-8')

    @parameterized.expand(COORD_ANNO_TYPES)
    def test_align_multiple_input_files(self, anno_repr):
        """Sequences from separate FASTA files get correct per-file header mapping."""
        fa1 = self._write_fa('file1.fa', [
            ('alpha', 'GTATCGATCGACGT'),
            ('beta', 'GCTAGCTAGCTAGCTA'),
        ])
        fa2 = self._write_fa('file2.fa', [
            ('gamma', 'ATCGATCGAAAAACCCCC'),
        ])
        query_fa = self._write_fa('query.fa', [
            ('q_alpha', 'TATCGATCG'),
            ('q_gamma', 'AAAAACCCCC'),
        ])

        graph, anno = self._setup_graph([fa1, fa2], anno_repr)
        rows = self._run_align(graph, anno, query_fa)
        self.assertEqual(len(rows), 2)

        # q_alpha -> alpha:2-10 (in file1.fa)
        self.assertEqual(rows[0][0], 'q_alpha')
        self.assertEqual(rows[0][3], 'TATCGATCG')
        self.assertEqual(rows[0][6], '9=')
        self.assertEqual(rows[0][8], 'alpha:2-10')

        # q_gamma -> gamma:9-18 (in file2.fa)
        self.assertEqual(rows[1][0], 'q_gamma')
        self.assertEqual(rows[1][3], 'AAAAACCCCC')
        self.assertEqual(rows[1][6], '10=')
        self.assertEqual(rows[1][8], 'gamma:9-18')

    @parameterized.expand(COORD_ANNO_TYPES)
    def test_align_with_mismatch(self, anno_repr):
        """Non-trivial CIGAR (substitution) still gets correct coord mapping."""
        test_fa = self._write_fa('seqs.fa', [
            ('seq1', 'GTATCGATCGACGTACGT'),
        ])
        # Position 12: G->A mismatch
        query_fa = self._write_fa('query.fa', [
            ('mismatch_q', 'GTATCGATCGACATACGT'),
        ])

        graph, anno = self._setup_graph(test_fa, anno_repr)
        rows = self._run_align(graph, anno, query_fa, extra_flags='--align-min-exact-match 0.0')
        self.assertEqual(len(rows), 1)

        self.assertEqual(rows[0][0], 'mismatch_q')
        self.assertEqual(rows[0][1], 'GTATCGATCGACATACGT')
        self.assertEqual(rows[0][3], 'GTATCGATCGACGTACGT')
        self.assertEqual(rows[0][6], '12=1X5=')
        self.assertEqual(rows[0][8], 'seq1:1-18')

    @parameterized.expand(COORD_ANNO_TYPES)
    def test_align_sequence_boundaries(self, anno_repr):
        """Matches at very start (coord 1) and very end of a sequence."""
        test_fa = self._write_fa('seqs.fa', [
            ('seq1', 'ACGTACGTACGT'),
            ('seq2', 'TTTTTAAAAATTTTT'),
        ])
        query_start = self._write_fa('q_start.fa', [('q_start', 'ACGTACGTA')])
        query_end = self._write_fa('q_end.fa', [('q_end', 'GTACGTACGT')])

        graph, anno = self._setup_graph(test_fa, anno_repr)

        # Match at start of seq1 -> coord starts at 1
        rows = self._run_align(graph, anno, query_start)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0][0], 'q_start')
        self.assertEqual(rows[0][6], '9=')
        self.assertEqual(rows[0][8], 'seq1:1-9')

        # Match at end of seq1 -> coord ends at len(seq1)=12
        rows = self._run_align(graph, anno, query_end)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0][0], 'q_end')
        self.assertEqual(rows[0][6], '10=')
        self.assertEqual(rows[0][8], 'seq1:3-12')

    @parameterized.expand(COORD_ANNO_TYPES)
    def test_align_cross_sequence_boundary(self, anno_repr):
        """Alignment path crosses from one indexed sequence into the next.

        When two sequences share the boundary k-1-mer, the graph connects
        them and the aligner follows a path whose coordinate range extends
        past the end of the first sequence.  The label is split across
        consecutive headers, joined by ';'.
        """
        # seq1 ends with ...ACGT, seq2 starts with ACGT... — they share the
        # k-1-mer 'ACGT' at the boundary, so the graph connects them.
        #   seq1: AAAAACGTACGT  (12 bp, k-mer coords 0-7)
        #   seq2: ACGTTTTTTTTT  (12 bp, k-mer coords 8-15)
        test_fa = self._write_fa('seqs.fa', [
            ('seq1', 'AAAAACGTACGT'),
            ('seq2', 'ACGTTTTTTTTT'),
        ])
        # The aligner finds a 16bp path crossing the boundary (8bp soft
        # clip at the start, since the first 8 query bp cannot extend into
        # a longer run than the cross-boundary one).
        query_fa = self._write_fa('query.fa', [
            ('q_concat', 'AAAAACGTACGTACGTTTTTTTTT'),
        ])

        graph, anno = self._setup_graph(test_fa, anno_repr)
        rows = self._run_align(graph, anno, query_fa)
        self.assertEqual(len(rows), 1)

        self.assertEqual(rows[0][0], 'q_concat')
        self.assertEqual(rows[0][6], '8S16=')
        # Range 8 nt in seq1 (local 5-12) + 8 nt in seq2 (local 1-8).
        self.assertEqual(rows[0][8], 'seq1:5-12;seq2:1-8')

    @parameterized.expand(COORD_ANNO_TYPES)
    def test_align_coord_to_header_matches_per_sequence_columns(self, anno_repr):
        """CoordToHeader output matches per-sequence-column annotation output.

        With CoordToHeader (one column + .seqs index) vs per-sequence columns
        (each sequence its own label via --anno-header), the aligner should
        produce the same labels with the same local coordinate ranges.
        """
        sequences = [
            ('seq1', 'GTATCGATCG'),
            ('seq2', 'GCTAGCTAGCTAGCTA'),
            ('seq3', 'ATCGATCGAAAAACCCCCGGGGGTTTTT'),
        ]
        test_fa = self._write_fa('seqs.fa', sequences)
        query_fa = self._write_fa('query.fa', [
            ('q1', 'TATCGATCG'),         # matches seq1 (and shared 'ATCGATCG' in seq3)
            ('q2', 'GCTAGCTAGCTAG'),     # matches seq2
            ('q3', 'AAAAACCCCC'),        # matches seq3
            ('q2_rc', 'CTAGCTAGCTAGC'),  # rc of a seq2 substring — hits on rev strand
        ])

        # Mode A: CoordToHeader — one column (filename), headers via .seqs.
        graph_a, anno_a = self._setup_graph(test_fa, anno_repr)
        rows_a = self._run_align(graph_a, anno_a, query_fa, only_forwards=False)

        # Mode B: per-sequence columns — each FASTA header is its own label.
        graph_b_base = self.tempdir.name + '/graph_b'
        graph_b = graph_b_base + '.dbg'
        anno_b_base = self.tempdir.name + '/anno_b'
        anno_b = anno_b_base + coord_anno_file_extension[anno_repr]
        self._build_graph(test_fa, graph_b_base, k=5, repr='succinct', mode='basic')
        self._annotate_graph(test_fa, graph_b, anno_b_base, anno_repr, anno_type='header')
        rows_b = self._run_align(graph_b, anno_b, query_fa, only_forwards=False)

        # The CIGAR (col 6) and label field (col 8) must match between modes.
        self.assertEqual(len(rows_a), len(rows_b))
        for row_a, row_b in zip(rows_a, rows_b):
            self.assertEqual(row_a[0], row_b[0])  # query name
            self.assertEqual(row_a[6], row_b[6],
                             f"CIGAR mismatch for {row_a[0]}: "
                             f"CoordToHeader={row_a[6]!r} vs per-sequence={row_b[6]!r}")
            # Normalize semicolon-separated labels (order may vary across modes).
            labels_a = sorted(row_a[8].split(';'))
            labels_b = sorted(row_b[8].split(';'))
            self.assertEqual(labels_a, labels_b,
                             f"label mismatch for {row_a[0]}: "
                             f"CoordToHeader={row_a[8]!r} vs per-sequence={row_b[8]!r}")


if __name__ == '__main__':
    unittest.main()
