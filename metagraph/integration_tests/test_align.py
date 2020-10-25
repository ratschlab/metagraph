import unittest
from parameterized import parameterized
import subprocess
from subprocess import PIPE
from tempfile import TemporaryDirectory
import glob
import os


"""Test graph construction and alignment"""

METAGRAPH = './metagraph'
TEST_DATA_DIR = os.path.dirname(os.path.realpath(__file__)) + '/../tests/data'

graph_file_extension = {'succinct': '.dbg',
                        'bitmap': '.bitmapdbg',
                        'hash': '.orhashdbg',
                        'hashfast': '.hashfastdbg',
                        'hashstr': '.hashstrdbg'}

GRAPH_TYPES = [graph_type for graph_type, _ in graph_file_extension.items()]


class TestAlign(unittest.TestCase):
    def setUp(self):
        self.tempdir = TemporaryDirectory()

    @parameterized.expand(GRAPH_TYPES)
    def test_simple_align_all_graphs(self, representation):

        self.maxDiff = None

        construct_command = '{exe} build --mask-dummy --graph {repr} -k 11 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=representation,
            outfile=self.tempdir.name + '/genome.MT',
            input=TEST_DATA_DIR + '/genome.MT.fa'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/genome.MT' + graph_file_extension[representation],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 11', params_str[0])
        self.assertEqual('nodes (k): 16438', params_str[1])
        self.assertEqual('canonical mode: no', params_str[2])

        stats_command = '{exe} align -i {graph} --align-vertical-bandwidth 1000000 --discovery-fraction 0.0 --align-local {reads}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/genome.MT' + graph_file_extension[representation],
            reads=TEST_DATA_DIR + '/genome_MT1.fq',
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')
        self.assertEqual(params_str[0], 'MT-10/1\tAACAGAGAATAGTTTAAATTAGAATCTTAGCTTTGGGTGCTAATGGTGGAGTTAAAGACTTTTTCTCTGATTTGTCCTTGGAAAAAGGTTTTCATCTCCGGTTTACAAGACTGGTGTATTAGTTTATACTACAAGGACAGGCCCATTTGA\t+\tTAGAATCTTAG\t22\t11\t19S11=120S\t0')
        self.assertEqual(params_str[1], 'MT-8/1\tAAAACTAACCCCCTAATAAAATTAATTAACCACTCATTCATCGACCTCCCCACCCCATCCAACATCTCCGCATGATGAAACTTCGGCTCACTCCTTGGCGCCTGCCTGATCCTCCAAATCACCACAGGACTATTCCTAGCCATGCACTAC\t+\tAAAACTAACCCCCTAATAAAATTAATTAACCACTCATTCATCGACCTCCCCACCCCATCCAACATCTCCGCATGATGAAACTTCGGCTCACTCCTTGGCGCCTGCCTGATCCTCCAAATCACCACAGGACTATTCCTAGCCATGCACTAC\t300\t150\t150=\t0')
        self.assertEqual(params_str[2], 'MT-6/1\tATATGACTAGCTTACACAATAGCTTTTATAGTAAAGATACCTCTTTACGGACTCCACTTATGACTCCCTAAAGCCCATGTCGAAGCCCCCATCGCTGGGTCAATAGTACTTGCCGCAGTACTCTTAAAACTAGGCGGCTATGGTATAATA\t+\tATATGACTAGCTTACACAATAGCTTTTATAGTAAAGATACCTCTTTACGGACTCCACTTATGACTCCCTAAAGCCCATGTCGAAGCCCCCATCGCTGGGTCAATAGTACTTGCCGCAGTACTCTTAAAACTAGGCGGCTATGGTATAATA\t300\t150\t150=\t0')
        self.assertEqual(params_str[3], 'MT-4/1\tAGTATAGTAGTTCGCTTTGACTGGTGAAGTCTTAGCATGTACTGCTCGGAGGTTCGGTTCTGCTCCGAGGTCGCCCCAACCGAAATTTTTAATGCAGGTTTGGTAGTTTAGGACCTGTGGGTTTGTTAGGTACTGTTTGCATTAATAAAT\t*\t*\t0\t*\t*\t*')
        self.assertEqual(params_str[4], 'MT-2/1\tTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCACTTTCCACACAGACATCATAACAAAAAATTTCCACCAAACCCCCCCTCCCCCGCTTCTGGCCACAGCACTTAAACACATCTCTGCCAAAC\t+\tTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCACTTTCCACACAGACATCATAACAAAAAATTTCCACCAAACCCCCCCTCCCCCGCTTCTGGCCACAGCACTTAAACACATCTCTGCCAAAC\t300\t150\t150=\t0')

    @parameterized.expand(GRAPH_TYPES)
    def test_simple_align_banded_all_graphs(self, representation):

        self.maxDiff = None

        construct_command = '{exe} build --mask-dummy --graph {repr} -k 11 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=representation,
            outfile=self.tempdir.name + '/genome.MT',
            input=TEST_DATA_DIR + '/genome.MT.fa'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/genome.MT' + graph_file_extension[representation],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 11', params_str[0])
        self.assertEqual('nodes (k): 16438', params_str[1])
        self.assertEqual('canonical mode: no', params_str[2])

        stats_command = '{exe} align -i {graph} --align-vertical-bandwidth 10 --discovery-fraction 0.0 --align-local {reads}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/genome.MT' + graph_file_extension[representation],
            reads=TEST_DATA_DIR + '/genome_MT1.fq',
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')
        self.assertEqual(params_str[0], 'MT-10/1\tAACAGAGAATAGTTTAAATTAGAATCTTAGCTTTGGGTGCTAATGGTGGAGTTAAAGACTTTTTCTCTGATTTGTCCTTGGAAAAAGGTTTTCATCTCCGGTTTACAAGACTGGTGTATTAGTTTATACTACAAGGACAGGCCCATTTGA\t+\tTAGAATCTTAG\t22\t11\t19S11=120S\t0')
        self.assertEqual(params_str[1], 'MT-8/1\tAAAACTAACCCCCTAATAAAATTAATTAACCACTCATTCATCGACCTCCCCACCCCATCCAACATCTCCGCATGATGAAACTTCGGCTCACTCCTTGGCGCCTGCCTGATCCTCCAAATCACCACAGGACTATTCCTAGCCATGCACTAC\t+\tAAAACTAACCCCCTAATAAAATTAATTAACCACTCATTCATCGACCTCCCCACCCCATCCAACATCTCCGCATGATGAAACTTCGGCTCACTCCTTGGCGCCTGCCTGATCCTCCAAATCACCACAGGACTATTCCTAGCCATGCACTAC\t300\t150\t150=\t0')
        self.assertEqual(params_str[2], 'MT-6/1\tATATGACTAGCTTACACAATAGCTTTTATAGTAAAGATACCTCTTTACGGACTCCACTTATGACTCCCTAAAGCCCATGTCGAAGCCCCCATCGCTGGGTCAATAGTACTTGCCGCAGTACTCTTAAAACTAGGCGGCTATGGTATAATA\t+\tATATGACTAGCTTACACAATAGCTTTTATAGTAAAGATACCTCTTTACGGACTCCACTTATGACTCCCTAAAGCCCATGTCGAAGCCCCCATCGCTGGGTCAATAGTACTTGCCGCAGTACTCTTAAAACTAGGCGGCTATGGTATAATA\t300\t150\t150=\t0')
        self.assertEqual(params_str[3], 'MT-4/1\tAGTATAGTAGTTCGCTTTGACTGGTGAAGTCTTAGCATGTACTGCTCGGAGGTTCGGTTCTGCTCCGAGGTCGCCCCAACCGAAATTTTTAATGCAGGTTTGGTAGTTTAGGACCTGTGGGTTTGTTAGGTACTGTTTGCATTAATAAAT\t*\t*\t0\t*\t*\t*')
        self.assertEqual(params_str[4], 'MT-2/1\tTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCACTTTCCACACAGACATCATAACAAAAAATTTCCACCAAACCCCCCCTCCCCCGCTTCTGGCCACAGCACTTAAACACATCTCTGCCAAAC\t+\tTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCACTTTCCACACAGACATCATAACAAAAAATTTCCACCAAACCCCCCCTCCCCCGCTTCTGGCCACAGCACTTAAACACATCTCTGCCAAAC\t300\t150\t150=\t0')

    @parameterized.expand(['succinct'])
    def test_simple_align_json_all_graphs(self, representation):

        construct_command = '{exe} build --mask-dummy --graph {repr} -k 11 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=representation,
            outfile=self.tempdir.name + '/genome.MT',
            input=TEST_DATA_DIR + '/genome.MT.fa'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/genome.MT' + graph_file_extension[representation],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 11', params_str[0])
        self.assertEqual('nodes (k): 16438', params_str[1])
        self.assertEqual('canonical mode: no', params_str[2])

        stats_command = '{exe} align -i {graph} --discovery-fraction 0.0 --align-local {reads}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/genome.MT' + graph_file_extension[representation],
            reads=TEST_DATA_DIR + '/genome_MT1.fq',
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')

    @parameterized.expand(GRAPH_TYPES)
    def test_simple_align_fwd_rev_comp_all_graphs(self, representation):

        self.maxDiff = None
        construct_command = '{exe} build --mask-dummy --graph {repr} -k 11 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=representation,
            outfile=self.tempdir.name + '/genome.MT',
            input=TEST_DATA_DIR + '/genome.MT.fa'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/genome.MT' + graph_file_extension[representation],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 11', params_str[0])
        self.assertEqual('nodes (k): 16438', params_str[1])
        self.assertEqual('canonical mode: no', params_str[2])

        stats_command = '{exe} align --align-both-strands -i {graph} --discovery-fraction 0.0 --align-local {reads}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/genome.MT' + graph_file_extension[representation],
            reads=TEST_DATA_DIR + '/genome_MT1.fq',
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')
        self.assertEqual(params_str[0], 'MT-10/1\tAACAGAGAATAGTTTAAATTAGAATCTTAGCTTTGGGTGCTAATGGTGGAGTTAAAGACTTTTTCTCTGATTTGTCCTTGGAAAAAGGTTTTCATCTCCGGTTTACAAGACTGGTGTATTAGTTTATACTACAAGGACAGGCCCATTTGA\t-\tTCAAATGGGCCTGTCCTTGTAGTATAAACTAATACACCAGTCTTGTAAACCGGAGATGAAAACCTTTTTCCAAGGACAAATCAGAGAAAAAGTCTTTAACTCCACCATTAGCACCCAAAGCTAAGATTCTAATTTAAACTATTCTCTGTT\t300\t150\t150=\t0')
        self.assertEqual(params_str[1], 'MT-8/1\tAAAACTAACCCCCTAATAAAATTAATTAACCACTCATTCATCGACCTCCCCACCCCATCCAACATCTCCGCATGATGAAACTTCGGCTCACTCCTTGGCGCCTGCCTGATCCTCCAAATCACCACAGGACTATTCCTAGCCATGCACTAC\t+\tAAAACTAACCCCCTAATAAAATTAATTAACCACTCATTCATCGACCTCCCCACCCCATCCAACATCTCCGCATGATGAAACTTCGGCTCACTCCTTGGCGCCTGCCTGATCCTCCAAATCACCACAGGACTATTCCTAGCCATGCACTAC\t300\t150\t150=\t0')
        self.assertEqual(params_str[2], 'MT-6/1\tATATGACTAGCTTACACAATAGCTTTTATAGTAAAGATACCTCTTTACGGACTCCACTTATGACTCCCTAAAGCCCATGTCGAAGCCCCCATCGCTGGGTCAATAGTACTTGCCGCAGTACTCTTAAAACTAGGCGGCTATGGTATAATA\t+\tATATGACTAGCTTACACAATAGCTTTTATAGTAAAGATACCTCTTTACGGACTCCACTTATGACTCCCTAAAGCCCATGTCGAAGCCCCCATCGCTGGGTCAATAGTACTTGCCGCAGTACTCTTAAAACTAGGCGGCTATGGTATAATA\t300\t150\t150=\t0')
        self.assertEqual(params_str[3], 'MT-4/1\tAGTATAGTAGTTCGCTTTGACTGGTGAAGTCTTAGCATGTACTGCTCGGAGGTTCGGTTCTGCTCCGAGGTCGCCCCAACCGAAATTTTTAATGCAGGTTTGGTAGTTTAGGACCTGTGGGTTTGTTAGGTACTGTTTGCATTAATAAAT\t-\tATTTATTAATGCAAACAGTACCTAACAAACCCACAGGTCCTAAACTACCAAACCTGCATTAAAAATTTCGGTTGGGGCGACCTCGGAGCAGAACCCAACCTCCGAGCAGTACATGCTAAGACTTCACCAGTCAAAGCGAACTACTATACT\t295\t149\t95=1X54=\t0')
        self.assertEqual(params_str[4], 'MT-2/1\tTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCACTTTCCACACAGACATCATAACAAAAAATTTCCACCAAACCCCCCCTCCCCCGCTTCTGGCCACAGCACTTAAACACATCTCTGCCAAAC\t+\tTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCACTTTCCACACAGACATCATAACAAAAAATTTCCACCAAACCCCCCCTCCCCCGCTTCTGGCCACAGCACTTAAACACATCTCTGCCAAAC\t300\t150\t150=\t0')

    @parameterized.expand(['succinct'])
    def test_simple_align_fwd_rev_comp_json_all_graphs(self, representation):

        construct_command = '{exe} build --mask-dummy --graph {repr} -k 11 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=representation,
            outfile=self.tempdir.name + '/genome.MT',
            input=TEST_DATA_DIR + '/genome.MT.fa'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/genome.MT' + graph_file_extension[representation],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 11', params_str[0])
        self.assertEqual('nodes (k): 16438', params_str[1])
        self.assertEqual('canonical mode: no', params_str[2])

        stats_command = '{exe} align -o {output} --json --align-both-strands -i {graph} --discovery-fraction 0.0 --align-local {reads}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/genome.MT' + graph_file_extension[representation],
            reads=TEST_DATA_DIR + '/genome_MT1.fq',
            output=self.tempdir.name + '/genome.MT' + graph_file_extension[representation] + '.align.json',
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = open(self.tempdir.name + '/genome.MT' + graph_file_extension[representation] + '.align.json', 'r').readlines()
        ref_align_str = open(TEST_DATA_DIR + '/genome_MT1.align.json', 'r').readlines()
        for [a, b] in zip(params_str, ref_align_str):
            self.assertEqual(a, b)

    @parameterized.expand(['succinct'])
    def test_simple_align_edit_distance_all_graphs(self, representation):

        construct_command = '{exe} build --mask-dummy --graph {repr} -k 11 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=representation,
            outfile=self.tempdir.name + '/genome.MT',
            input=TEST_DATA_DIR + '/genome.MT.fa'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/genome.MT' + graph_file_extension[representation],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 11', params_str[0])
        self.assertEqual('nodes (k): 16438', params_str[1])
        self.assertEqual('canonical mode: no', params_str[2])

        stats_command = '{exe} align -o {output} --json --align-both-strands --align-edit-distance -i {graph} --discovery-fraction 0.0 --align-local {reads}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/genome.MT' + graph_file_extension[representation],
            reads=TEST_DATA_DIR + '/genome_MT1.fq',
            output=self.tempdir.name + '/genome.MT' + graph_file_extension[representation] + '.align.json',
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = open(self.tempdir.name + '/genome.MT' + graph_file_extension[representation] + '.align.json', 'r').readlines()
        ref_align_str = open(TEST_DATA_DIR + '/genome_MT1.align.edit.json', 'r').readlines()
        for [a, b] in zip(params_str, ref_align_str):
            self.assertEqual(a, b)

    @parameterized.expand(['succinct'])
    def test_align_chain(self, representation):

        construct_command = '{exe} build --graph {repr} -k 31 -c -o {outfile} {input}'.format(
            exe=METAGRAPH,
            repr=representation,
            outfile=self.tempdir.name + '/refs',
            input=TEST_DATA_DIR + '/refs.fa'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/refs' + graph_file_extension[representation],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 31', params_str[0])
        self.assertEqual('nodes (k): 6346', params_str[1])
        self.assertEqual('canonical mode: yes', params_str[2])

        stats_command = '{exe} align -i {graph} --align-min-seed-length 15 --discovery-fraction 0.0 {reads}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/refs' + graph_file_extension[representation],
            reads=TEST_DATA_DIR + '/reads.fa',
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')
        self.assertEqual(len(res.stdout), 14752)

if __name__ == '__main__':
    unittest.main()
