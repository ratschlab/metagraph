import unittest
from parameterized import parameterized
import subprocess
from subprocess import PIPE
from tempfile import TemporaryDirectory
import glob
import os


"""Test graph construction"""

METAGRAPH = './metagraph'
TEST_DATA_DIR = os.path.dirname(os.path.realpath(__file__)) + '/../tests/data'

graph_type_to_data = {'succinct': ('succinct', '.dbg', 'vector'),
                      'succinct_disk': ('succinct', '.dbg', 'vector_disk'),
                      'bitmap': ('bitmap', '.bitmapdbg', 'vector'),
                      'hash': ('hash', '.orhashdbg', 'vector'),
                      'hashfast': ('hashfast', '.hashfastdbg', 'vector'),
                      'hashstr': ('hashstr', '.hashstrdbg', 'vector')}

GRAPH_TYPES = [graph_type for graph_type, _ in graph_type_to_data.items()]


class TestBuild(unittest.TestCase):
    def setUp(self):
        self.tempdir = TemporaryDirectory()

    def __get_stats(self, graph_filename):
        stats_command = METAGRAPH + ' stats ' + graph_filename
        res = subprocess.run(stats_command.split(), stdout=PIPE, stderr=PIPE)
        return res

    @parameterized.expand(GRAPH_TYPES)
    def test_simple_all_graphs(self, representation):

        construct_command = '{exe} build --graph {repr} --container {container} -k 20 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=graph_type_to_data[representation][0],
            container=graph_type_to_data[representation][2],
            outfile=self.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        res = self.__get_stats(self.tempdir.name + '/graph' + graph_type_to_data[representation][1])
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 20', params_str[0])
        self.assertEqual('nodes (k): 591997', params_str[1])
        self.assertEqual('canonical mode: no', params_str[2])

    @parameterized.expand(['succinct'])
    def test_simple_bloom_graph(self, representation):

        construct_command = '{exe} build --graph {repr} --container {container} -k 20 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=graph_type_to_data[representation][0],
            container=graph_type_to_data[representation][2],
            outfile=self.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        res = self.__get_stats(self.tempdir.name + '/graph' + graph_type_to_data[representation][1])
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 20', params_str[0])
        self.assertEqual('nodes (k): 591997', params_str[1])
        self.assertEqual('canonical mode: no', params_str[2])

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
    @parameterized.expand(['succinct', 'succinct_disk', 'bitmap', 'hash'])  # , 'hashstr']:
    def test_simple_all_graphs_canonical(self, representation):
        """
        Build simple canonical graphs
        """

        construct_command = '{exe} build \
                --graph {repr} --canonical -k 20 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=graph_type_to_data[representation][0],
            container=graph_type_to_data[representation][2],
            outfile=self.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        res = self.__get_stats(self.tempdir.name + '/graph' + graph_type_to_data[representation][1])
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 20', params_str[0])
        self.assertEqual('nodes (k): 1159851', params_str[1])
        self.assertEqual('canonical mode: yes', params_str[2])

    @parameterized.expand(GRAPH_TYPES)
    def test_build_tiny_k(self, representation):
        args = [METAGRAPH, 'build', '--graph', graph_type_to_data[representation][0],
                '-k', '2',
                '--container', graph_type_to_data[representation][2],
                '-o', self.tempdir.name + '/graph',
                TEST_DATA_DIR + '/transcripts_1000.fa']
        construct_command = ' '.join(args)

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        res = self.__get_stats(self.tempdir.name + '/graph' + graph_type_to_data[representation][1])
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 2', params_str[0])
        self.assertEqual('nodes (k): 16', params_str[1])
        self.assertEqual('canonical mode: no', params_str[2])

    # TODO: add 'hashstr' once the canonical mode is implemented for it
    @parameterized.expand(['succinct', 'succinct_disk', 'bitmap', 'hash'])  # , 'hashstr']:
    def test_build_tiny_k_canonical(self, representation):
        args = [METAGRAPH, 'build', '--graph', graph_type_to_data[representation][0], '--canonical',
                '-k', '2',
                '--container', graph_type_to_data[representation][2],
                '-o', self.tempdir.name + '/graph',
                TEST_DATA_DIR + '/transcripts_1000.fa']
        construct_command = ' '.join(args)

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        res = self.__get_stats(self.tempdir.name + '/graph' + graph_type_to_data[representation][1])
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 2', params_str[0])
        self.assertEqual('nodes (k): 16', params_str[1])
        self.assertEqual('canonical mode: yes', params_str[2])

    @parameterized.expand(GRAPH_TYPES)
    def test_build_from_kmc(self, representation):
        construct_command = '{exe} build --graph {repr} --container {container} -k 11 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=graph_type_to_data[representation][0],
            container=graph_type_to_data[representation][2],
            outfile=self.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_1000_kmc_counters.kmc_suf'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        res = self.__get_stats(self.tempdir.name + '/graph' + graph_type_to_data[representation][1])
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 11', params_str[0])
        self.assertEqual('nodes (k): 469983', params_str[1])
        self.assertEqual('canonical mode: no', params_str[2])

    @parameterized.expand(GRAPH_TYPES)
    def test_build_from_kmc_both(self, representation):
        construct_command = '{exe} build --graph {repr} --container {container} -k 11 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=graph_type_to_data[representation][0],
            container=graph_type_to_data[representation][2],
            outfile=self.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_1000_kmc_counters_both_strands.kmc_suf'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        res = self.__get_stats(self.tempdir.name + '/graph' + graph_type_to_data[representation][1])
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 11', params_str[0])
        self.assertEqual('nodes (k): 802920', params_str[1])
        self.assertEqual('canonical mode: no', params_str[2])

    @parameterized.expand(['succinct', 'succinct_disk', 'bitmap', 'hash'])  # , 'hashstr']:
    def test_build_from_kmc_canonical(self, representation):
        construct_command = '{exe} build \
                --graph {repr} --container {container} --canonical -k 11 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=graph_type_to_data[representation][0],
            container=graph_type_to_data[representation][2],
            outfile=self.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_1000_kmc_counters.kmc_suf'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        res = self.__get_stats(self.tempdir.name + '/graph' + graph_type_to_data[representation][1])
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 11', params_str[0])
        self.assertEqual('nodes (k): 802920', params_str[1])
        self.assertEqual('canonical mode: yes', params_str[2])

    @parameterized.expand(['succinct', 'succinct_disk', 'bitmap', 'hash'])  # , 'hashstr']:
    def test_build_from_kmc_both_canonical(self, representation):
        construct_command = '{exe} build \
                --graph {repr} --container {container} --canonical -k 11 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=graph_type_to_data[representation][0],
            container=graph_type_to_data[representation][2],
            outfile=self.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_1000_kmc_counters_both_strands.kmc_suf'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        res = self.__get_stats(self.tempdir.name + '/graph' + graph_type_to_data[representation][1])
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 11', params_str[0])
        self.assertEqual('nodes (k): 802920', params_str[1])
        self.assertEqual('canonical mode: yes', params_str[2])

    def test_build_chunks_from_kmc(self):
        representation = 'succinct'

        # Build chunks
        for suffix in ['$', 'A', 'C', 'G', 'T']:
            construct_command = '{exe} build \
                                --graph {repr} -k 11 --suffix {suffix} -o {outfile} {input}'.format(
                exe=METAGRAPH,
                repr=graph_type_to_data[representation][0],
                outfile=self.tempdir.name + '/graph',
                input=TEST_DATA_DIR + '/transcripts_1000_kmc_counters.kmc_suf',
                suffix=suffix
            )

            res = subprocess.run([construct_command], shell=True)
            self.assertEqual(res.returncode, 0)

        # Concatenate chunks
        construct_command = '{exe} concatenate --len-suffix 1 --graph {repr} -i {chunk_filebase} -o {outfile}'.format(
            exe=METAGRAPH,
            repr=graph_type_to_data[representation][0],
            chunk_filebase=self.tempdir.name + '/graph',
            outfile=self.tempdir.name + '/graph_from_chunks',
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        # Check graph
        res = self.__get_stats(self.tempdir.name + '/graph_from_chunks'
                               + graph_type_to_data[representation][1])
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 11', params_str[0])
        self.assertEqual('nodes (k): 469983', params_str[1])
        self.assertEqual('canonical mode: no', params_str[2])

    def test_build_chunks_from_kmc_canonical(self):
        representation = 'succinct'

        # Build chunks
        for suffix in ['$', 'A', 'C', 'G', 'T']:
            construct_command = '{exe} build --graph {repr} --canonical -k 11 \
                    --suffix {suffix} -o {outfile} {input}'.format(
                exe=METAGRAPH,
                repr=graph_type_to_data[representation][0],
                outfile=self.tempdir.name + '/graph',
                input=TEST_DATA_DIR + '/transcripts_1000_kmc_counters.kmc_suf',
                suffix=suffix
            )

            res = subprocess.run([construct_command], shell=True)
            self.assertEqual(res.returncode, 0)

        # Concatenate chunks
        construct_command = '{exe} concatenate --len-suffix 1 --graph {repr} -i {chunk_filebase} -o {outfile}'.format(
            exe=METAGRAPH,
            repr=graph_type_to_data[representation][0],
            chunk_filebase=self.tempdir.name + '/graph',
            outfile=self.tempdir.name + '/graph_from_chunks',
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        # Check graph
        res = self.__get_stats(self.tempdir.name + '/graph_from_chunks'
                               + graph_type_to_data[representation][1])
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 11', params_str[0])
        self.assertEqual('nodes (k): 802920', params_str[1])
        self.assertEqual('canonical mode: yes', params_str[2])


if __name__ == '__main__':
    unittest.main()
