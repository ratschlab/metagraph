import unittest
import subprocess
from subprocess import PIPE
from tempfile import TemporaryDirectory
import glob
import os


"""Test graph annotation"""

METAGRAPH = './metagraph'
TEST_DATA_DIR = os.path.dirname(os.path.realpath(__file__)) + '/../tests/data'

graph_file_extension = {'succinct': '.dbg',
                        'bitmap': '.bitmapdbg',
                        'hash': '.orhashdbg',
                        'hashstr': '.hashstrdbg'}

anno_file_extension = {'column': '.column.annodbg',
                       'row': '.row.annodbg'}

NUM_THREADS = 4


class TestAnnotate(unittest.TestCase):
    def setUp(self):
        self.tempdir = TemporaryDirectory()

    def test_simple_all_graphs(self):
        """
        Simple annotation test
        """

        for graph_repr in ['succinct', 'bitmap', 'hash', 'hashstr']:

            construct_command = '{exe} build -p {num_threads} \
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
            res = subprocess.run(stats_command.split(), stdout=PIPE, stderr=PIPE)
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
                res = subprocess.run(anno_stats_command.split(), stdout=PIPE, stderr=PIPE)
                self.assertEqual(res.returncode, 0)
                params_str = res.stdout.decode().split('\n')[2:]
                self.assertEqual('labels:  100', params_str[0])
                self.assertEqual('objects: 46960', params_str[1])
                self.assertEqual('density: 0.0185072', params_str[2])
                self.assertEqual('representation: ' + anno_repr, params_str[3])

    def test_simple_all_graphs_canonical(self):
        """
        Simple annotation test
        """

        # TODO: add 'hashstr' once the canonical mode is implemented for it
        for graph_repr in ['succinct', 'bitmap', 'hash']:  # , 'hashstr']:

            construct_command = '{exe} build -p {num_threads} \
                    --graph {repr} --canonical -k 20 -o {outfile} {input}'.format(
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
            res = subprocess.run(stats_command.split(), stdout=PIPE, stderr=PIPE)
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
                res = subprocess.run(anno_stats_command.split(), stdout=PIPE, stderr=PIPE)
                self.assertEqual(res.returncode, 0)
                params_str = res.stdout.decode().split('\n')[2:]
                self.assertEqual('labels:  100', params_str[0])
                self.assertEqual('objects: 91584', params_str[1])
                self.assertEqual('density: 0.00948888', params_str[2])
                self.assertEqual('representation: ' + anno_repr, params_str[3])

    def test_simple_all_graphs_from_kmc(self):
        """
        Annotate non-canonical graph constructed from non-canonical KMC database
        """

        for graph_repr in ['succinct', 'bitmap', 'hash', 'hashstr']:

            construct_command = '{exe} build -p {num_threads} \
                    --graph {repr} -k 11 -o {outfile} {input}'.format(
                exe=METAGRAPH,
                num_threads=NUM_THREADS,
                repr=graph_repr,
                outfile=self.tempdir.name + '/graph',
                input=TEST_DATA_DIR + '/transcripts_1000_kmc_counters.kmc_suf'
            )

            res = subprocess.run([construct_command], shell=True)
            self.assertEqual(res.returncode, 0)

            stats_command = '{exe} stats {graph}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
            )
            res = subprocess.run(stats_command.split(), stdout=PIPE, stderr=PIPE)
            self.assertEqual(res.returncode, 0)
            params_str = res.stdout.decode().split('\n')[2:]
            self.assertEqual('k: 11', params_str[0])
            self.assertEqual('nodes (k): 469983', params_str[1])
            self.assertEqual('canonical mode: no', params_str[2])

            for anno_repr in ['row', 'column']:
                # build annotation
                annotate_command = '{exe} annotate --anno-label LabelName -i {graph} \
                        --anno-type {anno_repr} -o {outfile} {input}'.format(
                    exe=METAGRAPH,
                    graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                    anno_repr=anno_repr,
                    outfile=self.tempdir.name + '/annotation',
                    input=TEST_DATA_DIR + '/transcripts_1000_kmc_counters.kmc_suf'
                )
                res = subprocess.run([annotate_command], shell=True)
                self.assertEqual(res.returncode, 0)

                # check annotation
                anno_stats_command = '{exe} stats -a {annotation}'.format(
                    exe=METAGRAPH,
                    annotation=self.tempdir.name + '/annotation' + anno_file_extension[anno_repr],
                )
                res = subprocess.run(anno_stats_command.split(), stdout=PIPE, stderr=PIPE)
                self.assertEqual(res.returncode, 0)
                params_str = res.stdout.decode().split('\n')[2:]
                self.assertEqual('labels:  1', params_str[0])
                self.assertEqual('objects: 469983', params_str[1])
                self.assertEqual('density: 1', params_str[2])
                self.assertEqual('representation: ' + anno_repr, params_str[3])

    def test_simple_all_graphs_from_kmc_both(self):
        """
        Annotate non-canonical graph constructed from canonical KMC database
        """

        for graph_repr in ['succinct', 'bitmap', 'hash', 'hashstr']:

            construct_command = '{exe} build -p {num_threads} \
                    --graph {repr} -k 11 -o {outfile} {input}'.format(
                exe=METAGRAPH,
                num_threads=NUM_THREADS,
                repr=graph_repr,
                outfile=self.tempdir.name + '/graph',
                input=TEST_DATA_DIR + '/transcripts_1000_kmc_counters_both_strands.kmc_suf'
            )

            res = subprocess.run([construct_command], shell=True)
            self.assertEqual(res.returncode, 0)

            stats_command = '{exe} stats {graph}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
            )
            res = subprocess.run(stats_command.split(), stdout=PIPE, stderr=PIPE)
            self.assertEqual(res.returncode, 0)
            params_str = res.stdout.decode().split('\n')[2:]
            self.assertEqual('k: 11', params_str[0])
            self.assertEqual('nodes (k): 802920', params_str[1])
            self.assertEqual('canonical mode: no', params_str[2])

            for anno_repr in ['row', 'column']:
                # build annotation
                annotate_command = '{exe} annotate --anno-label LabelName -i {graph} \
                        --anno-type {anno_repr} -o {outfile} {input}'.format(
                    exe=METAGRAPH,
                    graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                    anno_repr=anno_repr,
                    outfile=self.tempdir.name + '/annotation_single',
                    input=TEST_DATA_DIR + '/transcripts_1000_kmc_counters.kmc_suf'
                )
                res = subprocess.run([annotate_command], shell=True)
                self.assertEqual(res.returncode, 0)

                # check annotation
                anno_stats_command = '{exe} stats -a {annotation}'.format(
                    exe=METAGRAPH,
                    annotation=self.tempdir.name + '/annotation_single' + anno_file_extension[anno_repr],
                )
                res = subprocess.run(anno_stats_command.split(), stdout=PIPE, stderr=PIPE)
                self.assertEqual(res.returncode, 0)
                params_str = res.stdout.decode().split('\n')[2:]
                self.assertEqual('labels:  1', params_str[0])
                self.assertEqual('objects: 802920', params_str[1])
                self.assertEqual('density: 0.5853423', params_str[2])
                self.assertEqual('representation: ' + anno_repr, params_str[3])

                # both strands
                annotate_command = '{exe} annotate --anno-label LabelName -i {graph} \
                        --anno-type {anno_repr} -o {outfile} {input}'.format(
                    exe=METAGRAPH,
                    graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                    anno_repr=anno_repr,
                    outfile=self.tempdir.name + '/annotation_both',
                    input=TEST_DATA_DIR + '/transcripts_1000_kmc_counters_both_strands.kmc_suf'
                )
                res = subprocess.run([annotate_command], shell=True)
                self.assertEqual(res.returncode, 0)

                # check annotation
                anno_stats_command = '{exe} stats -a {annotation}'.format(
                    exe=METAGRAPH,
                    annotation=self.tempdir.name + '/annotation_both' + anno_file_extension[anno_repr],
                )
                res = subprocess.run(anno_stats_command.split(), stdout=PIPE, stderr=PIPE)
                self.assertEqual(res.returncode, 0)
                params_str = res.stdout.decode().split('\n')[2:]
                self.assertEqual('labels:  1', params_str[0])
                self.assertEqual('objects: 802920', params_str[1])
                self.assertEqual('density: 1', params_str[2])
                self.assertEqual('representation: ' + anno_repr, params_str[3])

    def test_simple_all_graphs_from_kmc_both(self):
        """
        Annotate canonical graph with k-mers from KMC
        """

        # TODO: add 'hashstr' once the canonical mode is implemented for it
        for graph_repr in ['succinct', 'bitmap', 'hash']:  # , 'hashstr']:

            construct_command = '{exe} build -p {num_threads} \
                    --graph {repr} --canonical -k 11 -o {outfile} {input}'.format(
                exe=METAGRAPH,
                num_threads=NUM_THREADS,
                repr=graph_repr,
                outfile=self.tempdir.name + '/graph',
                input=TEST_DATA_DIR + '/transcripts_1000_kmc_counters.kmc_suf'
            )

            res = subprocess.run([construct_command], shell=True)
            self.assertEqual(res.returncode, 0)

            stats_command = '{exe} stats {graph}'.format(
                exe=METAGRAPH,
                graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
            )
            res = subprocess.run(stats_command.split(), stdout=PIPE, stderr=PIPE)
            self.assertEqual(res.returncode, 0)
            params_str = res.stdout.decode().split('\n')[2:]
            self.assertEqual('k: 11', params_str[0])
            self.assertEqual('nodes (k): 802920', params_str[1])
            self.assertEqual('canonical mode: yes', params_str[2])

            for anno_repr in ['row', 'column']:
                # build annotation
                annotate_command = '{exe} annotate --anno-label LabelName -i {graph} \
                        --anno-type {anno_repr} -o {outfile} {input}'.format(
                    exe=METAGRAPH,
                    graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                    anno_repr=anno_repr,
                    outfile=self.tempdir.name + '/annotation_single',
                    input=TEST_DATA_DIR + '/transcripts_1000_kmc_counters.kmc_suf'
                )
                res = subprocess.run([annotate_command], shell=True)
                self.assertEqual(res.returncode, 0)

                # check annotation
                anno_stats_command = '{exe} stats -a {annotation}'.format(
                    exe=METAGRAPH,
                    annotation=self.tempdir.name + '/annotation_single' + anno_file_extension[anno_repr],
                )
                res = subprocess.run(anno_stats_command.split(), stdout=PIPE, stderr=PIPE)
                self.assertEqual(res.returncode, 0)
                params_str = res.stdout.decode().split('\n')[2:]
                self.assertEqual('labels:  1', params_str[0])
                self.assertEqual('objects: 802920', params_str[1])
                self.assertEqual('density: 0.5', params_str[2])
                self.assertEqual('representation: ' + anno_repr, params_str[3])

                # both strands
                annotate_command = '{exe} annotate --anno-label LabelName -i {graph} \
                        --anno-type {anno_repr} -o {outfile} {input}'.format(
                    exe=METAGRAPH,
                    graph=self.tempdir.name + '/graph' + graph_file_extension[graph_repr],
                    anno_repr=anno_repr,
                    outfile=self.tempdir.name + '/annotation_both',
                    input=TEST_DATA_DIR + '/transcripts_1000_kmc_counters_both_strands.kmc_suf'
                )
                res = subprocess.run([annotate_command], shell=True)
                self.assertEqual(res.returncode, 0)

                # check annotation
                anno_stats_command = '{exe} stats -a {annotation}'.format(
                    exe=METAGRAPH,
                    annotation=self.tempdir.name + '/annotation_both' + anno_file_extension[anno_repr],
                )
                res = subprocess.run(anno_stats_command.split(), stdout=PIPE, stderr=PIPE)
                self.assertEqual(res.returncode, 0)
                params_str = res.stdout.decode().split('\n')[2:]
                self.assertEqual('labels:  1', params_str[0])
                self.assertEqual('objects: 802920', params_str[1])
                self.assertEqual('density: 0.5', params_str[2])
                self.assertEqual('representation: ' + anno_repr, params_str[3])


if __name__ == '__main__':
    unittest.main()
