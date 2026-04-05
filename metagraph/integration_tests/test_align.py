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
    """Test that 'metagraph align' resolves accession numbers via .seqs index (issue #608)."""

    def setUp(self):
        self.tempdir = TemporaryDirectory()

        self.test_fa = self.tempdir.name + '/test_sequences.fa'
        with open(self.test_fa, 'w') as f:
            f.write('>seq1\n')
            f.write('GTATCGATCG\n')
            f.write('>seq2\n')
            f.write('GCTAGCTAGCTAGCTA\n')
            f.write('>seq3\n')
            f.write('ATCGATCGAAAAACCCCCGGGGGTTTTT\n')

        self.query_fa = self.tempdir.name + '/query.fa'
        with open(self.query_fa, 'w') as f:
            f.write('>query1\n')
            f.write('TATCGATCG\n')

    def tearDown(self):
        if getattr(self, 'tempdir', None) is not None:
            self.tempdir.cleanup()

    def _index_headers(self, graph, anno_base):
        """Create the CoordToHeader (.seqs) index next to the annotation."""
        cmd = (f"{METAGRAPH} annotate --anno-filename --index-header-coords -v "
               f"-i {graph} -o {anno_base} {self.test_fa}" + MMAP_FLAG)
        res = subprocess.run([cmd], shell=True, stdout=PIPE, stderr=PIPE)
        self.assertEqual(res.returncode, 0,
                         f"Indexing headers failed: {res.stderr.decode()}")
        self.assertTrue(os.path.exists(anno_base + '.seqs'))

    @parameterized.expand(COORD_ANNO_TYPES)
    def test_align_without_seqs_shows_file_paths(self, anno_repr):
        """Without .seqs index, align output should contain file paths."""
        graph_base = self.tempdir.name + '/graph'
        graph = graph_base + '.dbg'
        anno_base = self.tempdir.name + '/annotation'
        anno = anno_base + coord_anno_file_extension[anno_repr]

        self._build_graph(self.test_fa, graph_base, k=5, repr='succinct', mode='basic')
        self._annotate_graph(self.test_fa, graph, anno_base, anno_repr, anno_type='filename')

        align_cmd = (f'{METAGRAPH} align --align-only-forwards '
                     f'-i {graph} -a {anno} {self.query_fa}' + MMAP_FLAG)
        res = subprocess.run([align_cmd], shell=True, stdout=PIPE, stderr=PIPE)
        self.assertEqual(res.returncode, 0, f"Align failed: {res.stderr.decode()}")

        output = res.stdout.decode()
        self.assertGreater(len(output.strip()), 0, "Alignment produced no output")
        self.assertIn(self.tempdir.name, output,
                      "Without .seqs, output should contain the file path")

    @parameterized.expand(COORD_ANNO_TYPES)
    def test_align_with_seqs_shows_headers(self, anno_repr):
        """With .seqs index, align output should contain header names, not file paths."""
        graph_base = self.tempdir.name + '/graph'
        graph = graph_base + '.dbg'
        anno_base = self.tempdir.name + '/annotation'
        anno = anno_base + coord_anno_file_extension[anno_repr]

        self._build_graph(self.test_fa, graph_base, k=5, repr='succinct', mode='basic')
        self._annotate_graph(self.test_fa, graph, anno_base, anno_repr, anno_type='filename')
        self._index_headers(graph, anno_base)

        align_cmd = (f'{METAGRAPH} align --align-only-forwards '
                     f'-i {graph} -a {anno} {self.query_fa}' + MMAP_FLAG)
        res = subprocess.run([align_cmd], shell=True, stdout=PIPE, stderr=PIPE)
        self.assertEqual(res.returncode, 0,
                         f"Align with .seqs failed: {res.stderr.decode()}")

        output = res.stdout.decode()
        self.assertGreater(len(output.strip()), 0, "Alignment produced no output")

        # File path should not appear when .seqs is loaded
        self.assertNotIn(self.tempdir.name, output,
                         "With .seqs, output should not contain file paths")

        # At least one header name should appear (query matches seq1 and seq3)
        has_header = any(h in output for h in ['seq1', 'seq3'])
        self.assertTrue(has_header,
                        f"With .seqs, output should contain header names. Got:\n{output}")

    @parameterized.expand(COORD_ANNO_TYPES)
    def test_align_no_coord_mapping_flag(self, anno_repr):
        """With --no-coord-mapping, align should show file paths even when .seqs exists."""
        graph_base = self.tempdir.name + '/graph'
        graph = graph_base + '.dbg'
        anno_base = self.tempdir.name + '/annotation'
        anno = anno_base + coord_anno_file_extension[anno_repr]

        self._build_graph(self.test_fa, graph_base, k=5, repr='succinct', mode='basic')
        self._annotate_graph(self.test_fa, graph, anno_base, anno_repr, anno_type='filename')
        self._index_headers(graph, anno_base)

        align_cmd = (f'{METAGRAPH} align --align-only-forwards --no-coord-mapping '
                     f'-i {graph} -a {anno} {self.query_fa}' + MMAP_FLAG)
        res = subprocess.run([align_cmd], shell=True, stdout=PIPE, stderr=PIPE)
        self.assertEqual(res.returncode, 0,
                         f"Align with --no-coord-mapping failed: {res.stderr.decode()}")

        output = res.stdout.decode()
        self.assertGreater(len(output.strip()), 0, "Alignment produced no output")

        # With --no-coord-mapping, should still show file paths
        self.assertIn(self.tempdir.name, output,
                      "With --no-coord-mapping, output should contain the file path")

    @parameterized.expand(COORD_ANNO_TYPES)
    def test_align_coords_match_query_coords(self, anno_repr):
        """Verify that align with .seqs produces headers consistent with query --query-mode coords."""
        graph_base = self.tempdir.name + '/graph'
        graph = graph_base + '.dbg'
        anno_base = self.tempdir.name + '/annotation'
        anno = anno_base + coord_anno_file_extension[anno_repr]

        self._build_graph(self.test_fa, graph_base, k=5, repr='succinct', mode='basic')
        self._annotate_graph(self.test_fa, graph, anno_base, anno_repr, anno_type='filename')
        self._index_headers(graph, anno_base)

        # Run query to get the expected header names (same mode as test_query coord tests)
        query_cmd = (f'{METAGRAPH} query --batch-size 0 --query-mode coords '
                     f'-i {graph} -a {anno} --min-kmers-fraction-label 0.0 '
                     f'{self.query_fa}' + MMAP_FLAG)
        res_query = subprocess.run([query_cmd], shell=True, stdout=PIPE, stderr=PIPE)
        self.assertEqual(res_query.returncode, 0)
        query_output = res_query.stdout.decode()

        # Run align
        align_cmd = (f'{METAGRAPH} align --align-only-forwards '
                     f'-i {graph} -a {anno} {self.query_fa}' + MMAP_FLAG)
        res_align = subprocess.run([align_cmd], shell=True, stdout=PIPE, stderr=PIPE)
        self.assertEqual(res_align.returncode, 0)
        align_output = res_align.stdout.decode()

        # Query lines look like: 0\tquery1\t<seq1>:...\t<seq3>:...
        query_headers = set(re.findall(r'<([^>]+)>', query_output))
        self.assertTrue(len(query_headers) > 0, f"Query should list headers:\n{query_output}")

        # Align prints labels as header:start-end (no angle brackets), possibly several paths
        # on one line. Match only headers that query also reported.
        align_headers = {h for h in query_headers if f'{h}:' in align_output}

        self.assertTrue(len(align_headers) > 0,
                        f"Align should reference at least one query header. Output:\n{align_output}")
        self.assertTrue(align_headers.issubset(query_headers),
                        f"Align headers {align_headers} should be a subset of query headers {query_headers}")


if __name__ == '__main__':
    unittest.main()
