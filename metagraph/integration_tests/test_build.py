import unittest
from parameterized import parameterized
import subprocess
from subprocess import PIPE
from tempfile import TemporaryDirectory
import glob
import os
from base import PROTEIN_MODE, DNA_MODE, TestingBase, METAGRAPH, TEST_DATA_DIR


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

succinct_states = {'stat', 'fast', 'small', 'dynamic'}

BUILDS = [name for name, _ in build_params.items()]


class TestBuild(TestingBase):
    def setUp(self):
        self.tempdir = TemporaryDirectory()

    @parameterized.expand([repr for repr in BUILDS if not (repr == 'bitmap' and PROTEIN_MODE)])
    def test_simple_all_graphs(self, build):
        representation, tmp_dir = build_params[build]

        construct_command = '{exe} build --mask-dummy --graph {repr} --disk-swap {tmp_dir} -k 20 -o {outfile} {input}'.format(
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
        self.assertEqual('20', stats_graph['k'])
        self.assertEqual('591997', stats_graph['nodes (k)'])
        self.assertEqual('basic', stats_graph['mode'])

    @parameterized.expand(succinct_states)
    def test_build_succinct_inplace(self, state):
        construct_command = f'{METAGRAPH} build -k 20 --graph succinct --state {state} \
                                            --inplace \
                                            -o {self.tempdir.name}/graph \
                                            {TEST_DATA_DIR}/transcripts_1000.fa'

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_graph = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension['succinct'])
        self.assertEqual(stats_graph['returncode'], 0)
        self.assertEqual('20', stats_graph['k'])
        self.assertEqual('597931', stats_graph['nodes (k)'])
        self.assertEqual('basic', stats_graph['mode'])
        self.assertEqual(state, stats_graph['state'])

    @parameterized.expand(['succinct'])
    def test_simple_bloom_graph(self, build):
        representation, tmp_dir = build_params[build]

        construct_command = '{exe} build --mask-dummy --graph {repr} --disk-swap {tmp_dir} -k 20 -o {outfile} {input}'.format(
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
        self.assertEqual('20', stats_graph['k'])
        self.assertEqual('591997', stats_graph['nodes (k)'])
        self.assertEqual('basic', stats_graph['mode'])

        convert_command = '{exe} transform -o {outfile} --initialize-bloom {bloom_param} {input}'.format(
            exe=METAGRAPH,
            outfile=self.tempdir.name + '/graph',
            bloom_param='--bloom-fpp 0.1',
            input=self.tempdir.name + '/graph.dbg'
        )
        res = subprocess.run([convert_command], shell=True)
        self.assertEqual(res.returncode, 0)

        convert_command = '{exe} transform -o {outfile} --initialize-bloom {bloom_param} {input}'.format(
            exe=METAGRAPH,
            outfile=self.tempdir.name + '/graph',
            bloom_param='',
            input=self.tempdir.name + '/graph.dbg'
        )
        res = subprocess.run([convert_command], shell=True)
        self.assertEqual(res.returncode, 0)

    # TODO: add 'hashstr' once the canonical mode is implemented for it
    @parameterized.expand([repr for repr in BUILDS if repr != 'hashstr'])
    @unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets")
    def test_simple_all_graphs_canonical(self, build):
        """
        Build simple canonical graphs
        """
        representation, tmp_dir = build_params[build]

        construct_command = '{exe} build --mask-dummy \
                --graph {repr} --mode canonical -k 20 -o {outfile} {input}'.format(
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
        self.assertEqual('20', stats_graph['k'])
        self.assertEqual('1159851', stats_graph['nodes (k)'])
        self.assertEqual('canonical', stats_graph['mode'])

    @parameterized.expand(BUILDS)
    def test_build_tiny_k(self, build):
        representation, tmp_dir = build_params[build]

        args = [METAGRAPH, 'build', '--mask-dummy', '--graph', representation,
                '-k', '2',
                '--disk-swap', tmp_dir,
                '-o', self.tempdir.name + '/graph',
                TEST_DATA_DIR + '/transcripts_1000.fa']
        construct_command = ' '.join(args)

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_graph = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[representation])
        self.assertEqual(stats_graph['returncode'], 0)
        self.assertEqual('2', stats_graph['k'])
        self.assertEqual('16', stats_graph['nodes (k)'])
        self.assertEqual('basic', stats_graph['mode'])

    # TODO: add 'hashstr' once the canonical mode is implemented for it
    @parameterized.expand([repr for repr in BUILDS if repr != 'hashstr'])
    @unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets")
    def test_build_tiny_k_canonical(self, build):
        representation, tmp_dir = build_params[build]

        args = [METAGRAPH, 'build', '--mask-dummy', '--graph', representation, '--mode canonical',
                '-k', '2',
                '--disk-swap', tmp_dir,
                '-o', self.tempdir.name + '/graph',
                TEST_DATA_DIR + '/transcripts_1000.fa']
        construct_command = ' '.join(args)

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_graph = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[representation])
        self.assertEqual(stats_graph['returncode'], 0)
        self.assertEqual('2', stats_graph['k'])
        self.assertEqual('16', stats_graph['nodes (k)'])
        self.assertEqual('canonical', stats_graph['mode'])

    @parameterized.expand(BUILDS)
    def test_build_tiny_k_parallel(self, build):
        representation, tmp_dir = build_params[build]

        construct_command = f'{METAGRAPH} build --mask-dummy --graph {representation} \
                                -k 2 -p 100 --disk-swap {tmp_dir} \
                                -o {self.tempdir.name}/graph \
                                {TEST_DATA_DIR}/transcripts_1000.fa'

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_graph = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[representation])
        self.assertEqual(stats_graph['returncode'], 0)
        self.assertEqual('2', stats_graph['k'])
        self.assertEqual('16', stats_graph['nodes (k)'])
        self.assertEqual('basic', stats_graph['mode'])

    # TODO: add 'hashstr' once the canonical mode is implemented for it
    @parameterized.expand([repr for repr in BUILDS if repr != 'hashstr'])
    @unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets")
    def test_build_tiny_k_parallel_canonical(self, build):
        representation, tmp_dir = build_params[build]

        construct_command = f'{METAGRAPH} build --mask-dummy --graph {representation} \
                                --mode canonical \
                                -k 2 -p 100 --disk-swap {tmp_dir} \
                                -o {self.tempdir.name}/graph \
                                {TEST_DATA_DIR}/transcripts_1000.fa'

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_graph = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[representation])
        self.assertEqual(stats_graph['returncode'], 0)
        self.assertEqual('2', stats_graph['k'])
        self.assertEqual('16', stats_graph['nodes (k)'])
        self.assertEqual('canonical', stats_graph['mode'])

    @parameterized.expand(BUILDS)
    def test_build_from_kmc(self, build):
        representation, tmp_dir = build_params[build]

        construct_command = '{exe} build --mask-dummy --graph {repr} --disk-swap {tmp_dir} -k 11 -o {outfile} {input}'.format(
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
        self.assertEqual('11', stats_graph['k'])
        self.assertEqual('469983', stats_graph['nodes (k)'])
        self.assertEqual('basic', stats_graph['mode'])

    @parameterized.expand(BUILDS)
    def test_build_from_kmc_both(self, build):
        representation, tmp_dir = build_params[build]

        construct_command = '{exe} build --mask-dummy --graph {repr} --disk-swap {tmp_dir} -k 11 -o {outfile} {input}'.format(
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
        self.assertEqual('11', stats_graph['k'])
        self.assertEqual('802920', stats_graph['nodes (k)'])
        self.assertEqual('basic', stats_graph['mode'])

    @parameterized.expand([repr for repr in BUILDS if repr != 'hashstr'])
    @unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets")
    def test_build_from_kmc_canonical(self, build):
        representation, tmp_dir = build_params[build]

        construct_command = '{exe} build --mask-dummy \
                --graph {repr} --disk-swap {tmp_dir} --mode canonical -k 11 -o {outfile} {input}'.format(
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
        self.assertEqual('11', stats_graph['k'])
        self.assertEqual('802920', stats_graph['nodes (k)'])
        self.assertEqual('canonical', stats_graph['mode'])

    @parameterized.expand([repr for repr in BUILDS if repr != 'hashstr'])
    @unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets")
    def test_build_from_kmc_both_canonical(self, build):
        representation, tmp_dir = build_params[build]

        construct_command = '{exe} build --mask-dummy \
                --graph {repr} --disk-swap {tmp_dir} --mode canonical -k 11 -o {outfile} {input}'.format(
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
        self.assertEqual('11', stats_graph['k'])
        self.assertEqual('802920', stats_graph['nodes (k)'])
        self.assertEqual('canonical', stats_graph['mode'])

    @parameterized.expand(['succinct', 'succinct_disk'])
    @unittest.skipUnless(DNA_MODE, "Need to adapt suffixes for other alphabets")
    def test_build_chunks_from_kmc(self, build):
        representation, tmp_dir = build_params[build]

        # Build chunks
        for suffix in ['$', 'A', 'C', 'G', 'T']:
            construct_command = '{exe} build --mask-dummy --disk-swap {tmp_dir} \
                                --graph {repr} -k 11 --suffix {suffix} -o {outfile} {input}'.format(
                exe=METAGRAPH,
                repr=representation,
                tmp_dir=tmp_dir,
                outfile=self.tempdir.name + '/graph',
                input=TEST_DATA_DIR + '/transcripts_1000_kmc_counters.kmc_suf',
                suffix=suffix
            )

            res = subprocess.run([construct_command], shell=True)
            self.assertEqual(res.returncode, 0)

        # Concatenate chunks
        construct_command = '{exe} concatenate --len-suffix 1 --graph {repr} -i {chunk_filebase} -o {outfile}'.format(
            exe=METAGRAPH,
            repr=representation,
            chunk_filebase=self.tempdir.name + '/graph',
            outfile=self.tempdir.name + '/graph_from_chunks',
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        # Check graph
        stats_graph = self._get_stats(self.tempdir.name + '/graph_from_chunks'
                                      + graph_file_extension[representation])
        self.assertEqual(stats_graph['returncode'], 0)
        self.assertEqual('11', stats_graph['k'])
        self.assertEqual('469983', stats_graph['nodes (k)'])
        self.assertEqual('basic', stats_graph['mode'])

    @parameterized.expand(['succinct', 'succinct_disk'])
    @unittest.skipUnless(DNA_MODE, "Need to adapt suffixes for other alphabets")
    def test_build_chunks_from_kmc_canonical(self, build):
        representation, tmp_dir = build_params[build]

        # Build chunks
        for suffix in ['$', 'A', 'C', 'G', 'T']:
            construct_command = '{exe} build --mask-dummy --disk-swap {tmp_dir} --graph {repr} --mode canonical -k 11 \
                    --suffix {suffix} -o {outfile} {input}'.format(
                exe=METAGRAPH,
                repr=representation,
                tmp_dir=tmp_dir,
                outfile=self.tempdir.name + '/graph',
                input=TEST_DATA_DIR + '/transcripts_1000_kmc_counters.kmc_suf',
                suffix=suffix
            )

            res = subprocess.run([construct_command], shell=True)
            self.assertEqual(res.returncode, 0)

        # Concatenate chunks
        construct_command = '{exe} concatenate --len-suffix 1 --graph {repr} -i {chunk_filebase} --mode canonical -o {outfile}'.format(
            exe=METAGRAPH,
            repr=representation,
            chunk_filebase=self.tempdir.name + '/graph',
            outfile=self.tempdir.name + '/graph_from_chunks',
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        # Check graph
        stats_graph = self._get_stats(self.tempdir.name + '/graph_from_chunks'
                                      + graph_file_extension[representation])
        self.assertEqual(stats_graph['returncode'], 0)
        self.assertEqual('11', stats_graph['k'])
        self.assertEqual('802920', stats_graph['nodes (k)'])
        self.assertEqual('canonical', stats_graph['mode'])


if __name__ == '__main__':
    unittest.main()
