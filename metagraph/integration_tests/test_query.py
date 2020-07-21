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

graph_file_extension = {'succinct': '.dbg',
                        'bitmap': '.bitmapdbg',
                        'hash': '.orhashdbg',
                        'hashfast': '.hashfastdbg',
                        'hashstr': '.hashstrdbg'}

anno_file_extension = {'column': '.column.annodbg',
                       'row': '.row.annodbg'}

GRAPH_TYPES = [graph_type for graph_type, _ in graph_file_extension.items()]

NUM_THREADS = 4


class TestQuery(unittest.TestCase):
    def setUp(self):
        self.tempdir = TemporaryDirectory()

    @parameterized.expand(GRAPH_TYPES)
    def test_query_all_graphs(self, graph_repr):

        construct_command = '{exe} build --mask-dummy -p {num_threads} \
                --graph {repr} -k 20 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            num_threads=NUM_THREADS,
            repr=graph_repr,
            outfile=self.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_100.fa'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 20', params_str[0])
        self.assertEqual('nodes (k): 46960', params_str[1])
        self.assertEqual('canonical mode: no', params_str[2])

        for anno_repr in ['row', 'column']:
            # build annotation
            annotate_command = '{exe} annotate --anno-header -i {graph} \
                    --anno-type {anno_repr} -o {outfile} {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                anno_repr=anno_repr,
                outfile=self.tempdir.name + '/annotation',
                input=TEST_DATA_DIR + '/transcripts_100.fa'
            )
            res = subprocess.run([annotate_command], shell=True)
            self.assertEqual(res.returncode, 0)

            # check annotation
            anno_stats_command = '{exe} stats -a {annotation}'.format(
                exe=METAGRAPH,
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
            )
            res = subprocess.run(anno_stats_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            params_str = res.stdout.decode().split('\n')[2:]
            self.assertEqual('labels:  100', params_str[0])
            # self.assertEqual('objects: 46960', params_str[1])
            self.assertEqual('representation: ' + anno_repr, params_str[3])

            # query graph
            query_command = '{exe} query -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_1000.fa'
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 137140)

            query_command = '{exe} query --count-labels -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_1000.fa'
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 136959)

            # query graph (fwd and reverse)
            query_command = '{exe} query --fwd-and-reverse -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_1000.fa'
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 261390)

            query_command = '{exe} query --fwd-and-reverse --count-labels -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_1000.fa'
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 260215)

            # align to graph
            query_command = '{exe} query --align -i {graph} -a {annotation} --discovery-fraction 0.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 12241)

            query_command = '{exe} query --align --count-labels -i {graph} -a {annotation} --discovery-fraction 0.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 12347)

            # query graph (multi-threaded)
            query_command = '{exe} query -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 1.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_1000.fa',
                num_theads=NUM_THREADS
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 137140)

            query_command = '{exe} query --count-labels -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 1.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_1000.fa',
                num_theads=NUM_THREADS
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 136959)

            # query graph (fwd and reverse, multi-threaded)
            query_command = '{exe} query --fwd-and-reverse -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 1.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_1000.fa',
                num_theads=NUM_THREADS
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 261390)

            query_command = '{exe} query --fwd-and-reverse --count-labels -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 1.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_1000.fa',
                num_theads=NUM_THREADS
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 260215)

            # align to graph (multi-threaded)
            query_command = '{exe} query --align -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 0.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa',
                num_theads=NUM_THREADS
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 12241)

            query_command = '{exe} query --align --count-labels -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 0.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa',
                num_theads=NUM_THREADS
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 12347)

            # align to graph (fwd and reverse multi-threaded)
            query_command = '{exe} query --fwd-and-reverse --align -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 0.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa',
                num_theads=NUM_THREADS
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 20522)

            query_command = '{exe} query --fwd-and-reverse --align --count-labels -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 0.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa',
                num_theads=NUM_THREADS
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 20636)

    @parameterized.expand(['succinct'])
    def test_query_graphs_bloom(self, graph_repr):

        construct_command = '{exe} build --mask-dummy -p {num_threads} \
                --graph {repr} -k 20 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            num_threads=NUM_THREADS,
            repr=graph_repr,
            outfile=self.tempdir.name + '/graph.bloom',
            input=TEST_DATA_DIR + '/transcripts_100.fa'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph.bloom' + graph_file_extension[graph_repr],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 20', params_str[0])
        self.assertEqual('nodes (k): 46960', params_str[1])
        self.assertEqual('canonical mode: no', params_str[2])

        convert_command = '{exe} transform -o {outfile} --initialize-bloom {bloom_param} {input}'.format(
            exe=METAGRAPH,
            outfile=self.tempdir.name + '/graph.bloom',
            bloom_param='--bloom-fpp 0.1',
            input=self.tempdir.name + '/graph.bloom' + graph_file_extension[graph_repr],
        )
        res = subprocess.run([convert_command], shell=True)
        self.assertEqual(res.returncode, 0)

        for anno_repr in ['row', 'column']:
            # build annotation
            annotate_command = '{exe} annotate --anno-header -i {graph} \
                    --anno-type {anno_repr} -o {outfile} {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph.bloom' + graph_file_extension[graph_repr],
                anno_repr=anno_repr,
                outfile=self.tempdir.name + '/annotation',
                input=TEST_DATA_DIR + '/transcripts_100.fa'
            )
            res = subprocess.run([annotate_command], shell=True)
            self.assertEqual(res.returncode, 0)

            # check annotation
            anno_stats_command = '{exe} stats -a {annotation}'.format(
                exe=METAGRAPH,
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
            )
            res = subprocess.run(anno_stats_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            params_str = res.stdout.decode().split('\n')[2:]
            self.assertEqual('labels:  100', params_str[0])
            self.assertEqual('objects: 46960', params_str[1])
            self.assertEqual('representation: ' + anno_repr, params_str[3])

            # query graph
            query_command = '{exe} query -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph.bloom' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_1000.fa'
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 137140)

            query_command = '{exe} query --count-labels -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph.bloom' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_1000.fa'
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 136959)

            # query graph (multi-threaded)
            query_command = '{exe} query -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 1.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph.bloom' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_1000.fa',
                num_theads=NUM_THREADS
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 137140)

            query_command = '{exe} query --count-labels -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 1.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph.bloom' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_1000.fa',
                num_theads=NUM_THREADS
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 136959)

    @parameterized.expand(GRAPH_TYPES)
    def test_query_all_graphs_batch(self, graph_repr):

        construct_command = '{exe} build --mask-dummy -p {num_threads} \
                --graph {repr} -k 20 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            num_threads=NUM_THREADS,
            repr=graph_repr,
            outfile=self.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_100.fa'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 20', params_str[0])
        self.assertEqual('nodes (k): 46960', params_str[1])
        self.assertEqual('canonical mode: no', params_str[2])

        for anno_repr in ['row', 'column']:
            # build annotation
            annotate_command = '{exe} annotate --anno-header -i {graph} \
                    --anno-type {anno_repr} -o {outfile} {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                anno_repr=anno_repr,
                outfile=self.tempdir.name + '/annotation',
                input=TEST_DATA_DIR + '/transcripts_100.fa'
            )
            res = subprocess.run([annotate_command], shell=True)
            self.assertEqual(res.returncode, 0)

            # check annotation
            anno_stats_command = '{exe} stats -a {annotation}'.format(
                exe=METAGRAPH,
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
            )
            res = subprocess.run(anno_stats_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            params_str = res.stdout.decode().split('\n')[2:]
            self.assertEqual('labels:  100', params_str[0])
            # self.assertEqual('objects: 46960', params_str[1])
            self.assertEqual('representation: ' + anno_repr, params_str[3])

            # query graph
            query_command = '{exe} query --fast -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_1000.fa'
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 137140)

            query_command = '{exe} query --fast --count-labels -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_1000.fa'
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 136959)

            # query graph (fwd and reverse)
            query_command = '{exe} query --fast --fwd-and-reverse -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_1000.fa'
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 261390)

            query_command = '{exe} query --fast --fwd-and-reverse --count-labels -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_1000.fa'
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 260215)

            # align to graph
            query_command = '{exe} query --align --fast -i {graph} -a {annotation} --discovery-fraction 0.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 12241)

            query_command = '{exe} query --align --fast --count-labels -i {graph} -a {annotation} --discovery-fraction 0.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 12347)

            # query graph (multi-threaded)
            query_command = '{exe} query --fast -i {graph} -a {annotation} -p {num_threads} --discovery-fraction 1.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_1000.fa',
                num_threads=NUM_THREADS
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 137140)

            query_command = '{exe} query --fast --count-labels -i {graph} -a {annotation} -p {num_threads} --discovery-fraction 1.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_1000.fa',
                num_threads=NUM_THREADS
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 136959)

            # query graph (fwd and reverse, multi-threaded)
            query_command = '{exe} query --fast --fwd-and-reverse -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 1.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_1000.fa',
                num_theads=NUM_THREADS
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 261390)

            query_command = '{exe} query --fast --fwd-and-reverse --count-labels -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 1.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_1000.fa',
                num_theads=NUM_THREADS
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 260215)

            # align to graph (multi-threaded)
            query_command = '{exe} query --align --fast -i {graph} -a {annotation} -p {num_threads} --discovery-fraction 0.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa',
                num_threads=NUM_THREADS
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 12241)

            query_command = '{exe} query --align --fast --count-labels -i {graph} -a {annotation} -p {num_threads} --discovery-fraction 0.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa',
                num_threads=NUM_THREADS
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 12347)

            # align to graph (fwd and reverse multi-threaded)
            query_command = '{exe} query --fast --fwd-and-reverse --align -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 0.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa',
                num_theads=NUM_THREADS
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 20522)

            query_command = '{exe} query --fast --fwd-and-reverse --align --count-labels -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 0.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa',
                num_theads=NUM_THREADS
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 20636)

    @parameterized.expand(GRAPH_TYPES)
    def test_query_all_graphs_tiny_batch(self, graph_repr):

        construct_command = '{exe} build --mask-dummy -p {num_threads} \
                --graph {repr} -k 20 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            num_threads=NUM_THREADS,
            repr=graph_repr,
            outfile=self.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_100.fa'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 20', params_str[0])
        self.assertEqual('nodes (k): 46960', params_str[1])
        self.assertEqual('canonical mode: no', params_str[2])

        for anno_repr in ['row', 'column']:
            # build annotation
            annotate_command = '{exe} annotate --anno-header -i {graph} \
                    --anno-type {anno_repr} -o {outfile} {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                anno_repr=anno_repr,
                outfile=self.tempdir.name + '/annotation',
                input=TEST_DATA_DIR + '/transcripts_100.fa'
            )
            res = subprocess.run([annotate_command], shell=True)
            self.assertEqual(res.returncode, 0)

            # check annotation
            anno_stats_command = '{exe} stats -a {annotation}'.format(
                exe=METAGRAPH,
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
            )
            res = subprocess.run(anno_stats_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            params_str = res.stdout.decode().split('\n')[2:]
            self.assertEqual('labels:  100', params_str[0])
            # self.assertEqual('objects: 46960', params_str[1])
            self.assertEqual('representation: ' + anno_repr, params_str[3])

            # query graph
            query_command = '{exe} query --fast --batch-size 100 -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_1000.fa'
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 137140)

            query_command = '{exe} query --fast --batch-size 100 --count-labels -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_1000.fa'
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 136959)


class TestQueryCanonical(unittest.TestCase):
    def setUp(self):
        self.tempdir = TemporaryDirectory()

    @parameterized.expand(['succinct', 'bitmap', 'hash'])  # 'hashstr'
    def test_query_all_graphs(self, graph_repr):

        construct_command = '{exe} build --mask-dummy --canonical -p {num_threads} \
                --graph {repr} -k 20 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            num_threads=NUM_THREADS,
            repr=graph_repr,
            outfile=self.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_100.fa'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 20', params_str[0])
        self.assertEqual('nodes (k): 91584', params_str[1])
        self.assertEqual('canonical mode: yes', params_str[2])

        for anno_repr in ['row', 'column']:
            # build annotation
            annotate_command = '{exe} annotate --anno-header -i {graph} \
                    --anno-type {anno_repr} -o {outfile} {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                anno_repr=anno_repr,
                outfile=self.tempdir.name + '/annotation',
                input=TEST_DATA_DIR + '/transcripts_100.fa'
            )
            res = subprocess.run([annotate_command], shell=True)
            self.assertEqual(res.returncode, 0)

            # check annotation
            anno_stats_command = '{exe} stats -a {annotation}'.format(
                exe=METAGRAPH,
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
            )
            res = subprocess.run(anno_stats_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            params_str = res.stdout.decode().split('\n')[2:]
            self.assertEqual('labels:  100', params_str[0])
            self.assertEqual('objects: 91584', params_str[1])
            self.assertEqual('representation: ' + anno_repr, params_str[3])

            # query graph
            query_command = '{exe} query -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_1000.fa'
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 137269)

            query_command = '{exe} query --count-labels -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_1000.fa'
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 137093)

            # align to graph
            query_command = '{exe} query --align -i {graph} -a {annotation} --discovery-fraction 0.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 12839)

            query_command = '{exe} query --align --count-labels -i {graph} -a {annotation} --discovery-fraction 0.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 12969)

    @parameterized.expand(['succinct'])  # 'hashstr'
    def test_query_graphs_bloom(self, graph_repr):

        construct_command = '{exe} build --mask-dummy --canonical -p {num_threads} \
                --graph {repr} -k 20 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            num_threads=NUM_THREADS,
            repr=graph_repr,
            outfile=self.tempdir.name + '/graph.bloom',
            input=TEST_DATA_DIR + '/transcripts_100.fa'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph.bloom' + graph_file_extension[graph_repr],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 20', params_str[0])
        self.assertEqual('nodes (k): 91584', params_str[1])
        self.assertEqual('canonical mode: yes', params_str[2])

        convert_command = '{exe} transform -o {outfile} --initialize-bloom {bloom_param} {input}'.format(
            exe=METAGRAPH,
            outfile=self.tempdir.name + '/graph.bloom',
            bloom_param='--bloom-fpp 0.1',
            input=self.tempdir.name + '/graph.bloom' + graph_file_extension[graph_repr],
        )
        res = subprocess.run([convert_command], shell=True)
        self.assertEqual(res.returncode, 0)

        for anno_repr in ['row', 'column']:
            # build annotation
            annotate_command = '{exe} annotate --anno-header -i {graph} \
                    --anno-type {anno_repr} -o {outfile} {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph.bloom' + graph_file_extension[graph_repr],
                anno_repr=anno_repr,
                outfile=self.tempdir.name + '/annotation',
                input=TEST_DATA_DIR + '/transcripts_100.fa'
            )
            res = subprocess.run([annotate_command], shell=True)
            self.assertEqual(res.returncode, 0)

            # check annotation
            anno_stats_command = '{exe} stats -a {annotation}'.format(
                exe=METAGRAPH,
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
            )
            res = subprocess.run(anno_stats_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            params_str = res.stdout.decode().split('\n')[2:]
            self.assertEqual('labels:  100', params_str[0])
            self.assertEqual('objects: 91584', params_str[1])
            self.assertEqual('representation: ' + anno_repr, params_str[3])

            # query graph
            query_command = '{exe} query -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph.bloom' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_1000.fa'
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 137269)

            query_command = '{exe} query --count-labels -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph.bloom' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_1000.fa'
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 137093)

    @parameterized.expand(['succinct', 'bitmap', 'hash'])  # 'hashstr'
    def test_query_all_graphs_batch(self, graph_repr):

        construct_command = '{exe} build --mask-dummy --canonical -p {num_threads} \
                --graph {repr} -k 20 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            num_threads=NUM_THREADS,
            repr=graph_repr,
            outfile=self.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_100.fa'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 20', params_str[0])
        self.assertEqual('nodes (k): 91584', params_str[1])
        self.assertEqual('canonical mode: yes', params_str[2])

        for anno_repr in ['row', 'column']:
            # build annotation
            annotate_command = '{exe} annotate --anno-header -i {graph} \
                    --anno-type {anno_repr} -o {outfile} {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                anno_repr=anno_repr,
                outfile=self.tempdir.name + '/annotation',
                input=TEST_DATA_DIR + '/transcripts_100.fa'
            )
            res = subprocess.run([annotate_command], shell=True)
            self.assertEqual(res.returncode, 0)

            # check annotation
            anno_stats_command = '{exe} stats -a {annotation}'.format(
                exe=METAGRAPH,
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
            )
            res = subprocess.run(anno_stats_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            params_str = res.stdout.decode().split('\n')[2:]
            self.assertEqual('labels:  100', params_str[0])
            self.assertEqual('objects: 91584', params_str[1])
            self.assertEqual('representation: ' + anno_repr, params_str[3])

            # query graph
            query_command = '{exe} query --fast -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_1000.fa'
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 137269)

            query_command = '{exe} query --fast --count-labels -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_1000.fa'
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 137093)

            # align to graph
            query_command = '{exe} query --align --fast -i {graph} -a {annotation} --discovery-fraction 0.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 12839)

            query_command = '{exe} query --align --fast --count-labels -i {graph} -a {annotation} --discovery-fraction 0.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 12969)

    @parameterized.expand(['succinct', 'bitmap', 'hash'])  # 'hashstr'
    def test_query_all_graphs_batch_tiny_batch(self, graph_repr):

        construct_command = '{exe} build --mask-dummy --canonical -p {num_threads} \
                --graph {repr} -k 20 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            num_threads=NUM_THREADS,
            repr=graph_repr,
            outfile=self.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_100.fa'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        params_str = res.stdout.decode().split('\n')[2:]
        self.assertEqual('k: 20', params_str[0])
        self.assertEqual('nodes (k): 91584', params_str[1])
        self.assertEqual('canonical mode: yes', params_str[2])

        for anno_repr in ['row', 'column']:
            # build annotation
            annotate_command = '{exe} annotate --anno-header -i {graph} \
                    --anno-type {anno_repr} -o {outfile} {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                anno_repr=anno_repr,
                outfile=self.tempdir.name + '/annotation',
                input=TEST_DATA_DIR + '/transcripts_100.fa'
            )
            res = subprocess.run([annotate_command], shell=True)
            self.assertEqual(res.returncode, 0)

            # check annotation
            anno_stats_command = '{exe} stats -a {annotation}'.format(
                exe=METAGRAPH,
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
            )
            res = subprocess.run(anno_stats_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            params_str = res.stdout.decode().split('\n')[2:]
            self.assertEqual('labels:  100', params_str[0])
            self.assertEqual('objects: 91584', params_str[1])
            self.assertEqual('representation: ' + anno_repr, params_str[3])

            # query graph
            query_command = '{exe} query --fast --batch-size 100 -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_1000.fa'
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 137269)

            query_command = '{exe} query --fast --batch-size 100 --count-labels -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                input=TEST_DATA_DIR + '/transcripts_1000.fa'
            )
            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(len(res.stdout), 137093)


if __name__ == '__main__':
    unittest.main()
