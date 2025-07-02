import unittest
from parameterized import parameterized
import subprocess
from subprocess import PIPE
from tempfile import TemporaryDirectory
import glob
import os

from base import PROTEIN_MODE, DNA_MODE, TestingBase, METAGRAPH, TEST_DATA_DIR, NUM_THREADS


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
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        if representation != 'succinct':
            self.assertNotEqual(res.returncode, 0)
            return;

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


if __name__ == '__main__':
    unittest.main()
