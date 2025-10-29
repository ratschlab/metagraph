import unittest
from parameterized import parameterized
import subprocess
from subprocess import PIPE
from tempfile import TemporaryDirectory
import filecmp
import glob
import os
from base import PROTEIN_MODE, TestingBase, METAGRAPH, TEST_DATA_DIR, NUM_THREADS


"""Test graph annotation"""

graph_file_extension = {'succinct': '.dbg',
                        'bitmap': '.bitmapdbg',
                        'hash': '.orhashdbg',
                        # 'hashfast': '.hashfastdbg',
                        'hashstr': '.hashstrdbg'}

anno_file_extension = {'column': '.column.annodbg',
                       'row': '.row.annodbg'}

GRAPH_TYPES = [graph_type for graph_type, _ in graph_file_extension.items()]


class TestAnnotate(TestingBase):
    def setUp(self):
        self.tempdir = TemporaryDirectory()

    @parameterized.expand([repr for repr in GRAPH_TYPES if not (repr == 'bitmap' and PROTEIN_MODE)])
    def test_simple_all_graphs(self, graph_repr):

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

        stats_graph = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[graph_repr])
        self.assertEqual(stats_graph['returncode'], 0)
        self.assertEqual('20', stats_graph['k'])
        self.assertEqual('46960', stats_graph['nodes (k)'])
        self.assertEqual('basic', stats_graph['mode'])

        for anno_repr in ['row', 'column']:
            # build annotation
            annotate_command = f'{METAGRAPH} annotate --anno-header -p {NUM_THREADS} \
                                -i {self.tempdir.name}/graph{graph_file_extension[graph_repr]} \
                                --anno-type {anno_repr} \
                                -o {self.tempdir.name}/annotation \
                                {TEST_DATA_DIR}/transcripts_100.fa'

            res = subprocess.run([annotate_command], shell=True)
            self.assertEqual(res.returncode, 0)

            # check annotation
            stats_annotation = self._get_stats('-a ' + self.tempdir.name + '/annotation' + anno_file_extension[anno_repr])
            self.assertEqual(stats_annotation['returncode'], 0)
            self.assertEqual('100', stats_annotation['labels'])
            self.assertEqual(stats_graph['max index (k)'], stats_annotation['objects'])
            self.assertAlmostEqual(
                0.0185072 * (int(stats_graph['nodes (k)']) / int(stats_graph['max index (k)'])),
                float(stats_annotation['density']),
                places=6)
            self.assertEqual(anno_repr, stats_annotation['representation'])

    # TODO: add 'hashstr' once the canonical mode is implemented for it
    @parameterized.expand(['succinct', 'bitmap', 'hash'])  # , 'hashstr']:
    @unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets")
    def test_simple_all_graphs_canonical(self, graph_repr):

        construct_command = '{exe} build --mask-dummy -p {num_threads} \
                --graph {repr} --mode canonical -k 20 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            num_threads=NUM_THREADS,
            repr=graph_repr,
            outfile=self.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_100.fa'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_graph = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[graph_repr])
        self.assertEqual(stats_graph['returncode'], 0)
        self.assertEqual('20', stats_graph['k'])
        self.assertEqual('91584', stats_graph['nodes (k)'])
        self.assertEqual('canonical', stats_graph['mode'])

        for anno_repr in ['row', 'column']:
            # build annotation
            annotate_command = f'{METAGRAPH} annotate --anno-header -p {NUM_THREADS} \
                                -i {self.tempdir.name}/graph{graph_file_extension[graph_repr]} \
                                --anno-type {anno_repr} -o {self.tempdir.name}/annotation \
                                {TEST_DATA_DIR}/transcripts_100.fa'

            res = subprocess.run([annotate_command], shell=True)
            self.assertEqual(res.returncode, 0)

            # check annotation
            stats_annotation = self._get_stats('-a ' + self.tempdir.name + '/annotation' + anno_file_extension[anno_repr])
            self.assertEqual(stats_annotation['returncode'], 0)
            self.assertEqual('100', stats_annotation['labels'])
            self.assertEqual(stats_graph['max index (k)'], stats_annotation['objects'])
            self.assertAlmostEqual(
                0.00948888 * (int(stats_graph['nodes (k)']) / int(stats_graph['max index (k)'])),
                float(stats_annotation['density']),
                places=6)
            self.assertEqual(anno_repr, stats_annotation['representation'])

    @parameterized.expand(GRAPH_TYPES)
    def test_simple_all_graphs_from_kmc(self, graph_repr):
        """
        Annotate non-canonical graph constructed from non-canonical KMC database
        """

        construct_command = f'{METAGRAPH} build --mask-dummy -p {NUM_THREADS} \
                            --graph {graph_repr} -k 11 \
                            -o {self.tempdir.name}/graph \
                            {TEST_DATA_DIR}/transcripts_1000_kmc_counters.kmc_suf'

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_graph = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[graph_repr])
        self.assertEqual(stats_graph['returncode'], 0)
        self.assertEqual('11', stats_graph['k'])
        self.assertEqual('469983', stats_graph['nodes (k)'])
        self.assertEqual('basic', stats_graph['mode'])

        for anno_repr in ['row', 'column']:
            # build annotation
            annotate_command = f'{METAGRAPH} annotate --anno-label LabelName -p {NUM_THREADS} \
                                -i {self.tempdir.name}/graph{graph_file_extension[graph_repr]} \
                                --anno-type {anno_repr} -o {self.tempdir.name}/annotation \
                                {TEST_DATA_DIR}/transcripts_1000_kmc_counters.kmc_suf'

            res = subprocess.run([annotate_command], shell=True)
            self.assertEqual(res.returncode, 0)

            # check annotation
            stats_annotation = self._get_stats('-a ' + self.tempdir.name + '/annotation' + anno_file_extension[anno_repr])
            self.assertEqual(stats_annotation['returncode'], 0)
            self.assertEqual('1', stats_annotation['labels'])
            self.assertEqual(stats_graph['max index (k)'], stats_annotation['objects'])
            self.assertAlmostEqual(
                1 * (int(stats_graph['nodes (k)']) / int(stats_graph['max index (k)'])),
                float(stats_annotation['density']),
                places=6)
            self.assertEqual(anno_repr, stats_annotation['representation'])

    @parameterized.expand(GRAPH_TYPES)
    def test_simple_all_graphs_from_kmc_both(self, graph_repr):
        """
        Annotate non-canonical graph constructed from canonical KMC database
        """

        construct_command = f'{METAGRAPH} build --mask-dummy -p {NUM_THREADS} \
                            --graph {graph_repr} -k 11 \
                            -o {self.tempdir.name}/graph \
                            {TEST_DATA_DIR}/transcripts_1000_kmc_counters_both_strands.kmc_suf'

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_graph = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[graph_repr])
        self.assertEqual(stats_graph['returncode'], 0)
        self.assertEqual('11', stats_graph['k'])
        self.assertEqual('802920', stats_graph['nodes (k)'])
        self.assertEqual('basic', stats_graph['mode'])

        for anno_repr in ['row', 'column']:
            # build annotation
            annotate_command = f'{METAGRAPH} annotate --anno-label LabelName -p {NUM_THREADS} \
                                -i {self.tempdir.name}/graph{graph_file_extension[graph_repr]} \
                                --anno-type {anno_repr} -o {self.tempdir.name}/annotation_single \
                                {TEST_DATA_DIR}/transcripts_1000_kmc_counters.kmc_suf'

            res = subprocess.run([annotate_command], shell=True)
            self.assertEqual(res.returncode, 0)

            # check annotation
            stats_annotation = self._get_stats('-a ' + self.tempdir.name + '/annotation_single' + anno_file_extension[anno_repr])
            self.assertEqual(stats_annotation['returncode'], 0)
            self.assertEqual('1', stats_annotation['labels'])
            self.assertEqual(stats_graph['max index (k)'], stats_annotation['objects'])
            self.assertAlmostEqual(
                0.585342 * (int(stats_graph['nodes (k)']) / int(stats_graph['max index (k)'])),
                float(stats_annotation['density']),
                places=6)
            self.assertEqual(anno_repr, stats_annotation['representation'])

            # both strands
            annotate_command = f'{METAGRAPH} annotate --anno-label LabelName -p {NUM_THREADS} \
                                -i {self.tempdir.name}/graph{graph_file_extension[graph_repr]} \
                                --anno-type {anno_repr} -o {self.tempdir.name}/annotation_both \
                                {TEST_DATA_DIR}/transcripts_1000_kmc_counters_both_strands.kmc_suf'

            res = subprocess.run([annotate_command], shell=True)
            self.assertEqual(res.returncode, 0)

            # check annotation
            stats_annotation = self._get_stats('-a ' + self.tempdir.name + '/annotation_both' + anno_file_extension[anno_repr])
            self.assertEqual(stats_annotation['returncode'], 0)
            self.assertEqual('1', stats_annotation['labels'])
            self.assertEqual(stats_graph['max index (k)'], stats_annotation['objects'])
            self.assertAlmostEqual(
                1 * (int(stats_graph['nodes (k)']) / int(stats_graph['max index (k)'])),
                float(stats_annotation['density']),
                places=6)
            self.assertEqual(anno_repr, stats_annotation['representation'])

    # TODO: add 'hashstr' once the canonical mode is implemented for it
    @parameterized.expand(['succinct', 'bitmap', 'hash'])  # , 'hashstr']:
    @unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets")
    def test_simple_all_graphs_from_kmc_both_canonical(self, graph_repr):
        """
        Annotate canonical graph with k-mers from KMC
        """

        construct_command = f'{METAGRAPH} build --mask-dummy -p {NUM_THREADS} \
                            --graph {graph_repr} --mode canonical -k 11 \
                            -o {self.tempdir.name}/graph \
                            {TEST_DATA_DIR}/transcripts_1000_kmc_counters.kmc_suf'

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_graph = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[graph_repr])
        self.assertEqual(stats_graph['returncode'], 0)
        self.assertEqual('11', stats_graph['k'])
        self.assertEqual('802920', stats_graph['nodes (k)'])
        self.assertEqual('canonical', stats_graph['mode'])

        for anno_repr in ['row', 'column']:
            # build annotation
            annotate_command = f'{METAGRAPH} annotate --anno-label LabelName -p {NUM_THREADS} \
                                -i {self.tempdir.name}/graph{graph_file_extension[graph_repr]} \
                                --anno-type {anno_repr} -o {self.tempdir.name}/annotation_single \
                                {TEST_DATA_DIR}/transcripts_1000_kmc_counters.kmc_suf'

            res = subprocess.run([annotate_command], shell=True)
            self.assertEqual(res.returncode, 0)

            # check annotation
            stats_annotation = self._get_stats('-a ' + self.tempdir.name + '/annotation_single' + anno_file_extension[anno_repr])
            self.assertEqual(stats_annotation['returncode'], 0)
            self.assertEqual('1', stats_annotation['labels'])
            self.assertEqual(stats_graph['max index (k)'], stats_annotation['objects'])
            self.assertAlmostEqual(
                0.5 * (int(stats_graph['nodes (k)']) / int(stats_graph['max index (k)'])),
                float(stats_annotation['density']),
                places=6)
            self.assertEqual(anno_repr, stats_annotation['representation'])

            # both strands
            annotate_command = f'{METAGRAPH} annotate --anno-label LabelName -p {NUM_THREADS} \
                                -i {self.tempdir.name}/graph{graph_file_extension[graph_repr]} \
                                --anno-type {anno_repr} -o {self.tempdir.name}/annotation_both \
                                {TEST_DATA_DIR}/transcripts_1000_kmc_counters_both_strands.kmc_suf'

            res = subprocess.run([annotate_command], shell=True)
            self.assertEqual(res.returncode, 0)

            # check annotation
            stats_annotation = self._get_stats('-a ' + self.tempdir.name + '/annotation_both' + anno_file_extension[anno_repr])
            self.assertEqual(stats_annotation['returncode'], 0)
            self.assertEqual('1', stats_annotation['labels'])
            self.assertEqual(stats_graph['max index (k)'], stats_annotation['objects'])
            self.assertAlmostEqual(
                0.5 * (int(stats_graph['nodes (k)']) / int(stats_graph['max index (k)'])),
                float(stats_annotation['density']),
                places=6)
            self.assertEqual(anno_repr, stats_annotation['representation'])

    def test_annotate_with_disk_swap(self):
        graph_repr = 'succinct'
        anno_repr = 'column'

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

        stats_graph = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[graph_repr])
        self.assertEqual(stats_graph['returncode'], 0)
        self.assertEqual('20', stats_graph['k'])
        self.assertEqual('46960', stats_graph['nodes (k)'])
        self.assertEqual('basic', stats_graph['mode'])

        # build annotation
        annotate_command = f'{METAGRAPH} annotate --anno-header \
                            --disk-swap {self.tempdir.name} --mem-cap-gb 1e-6 \
                            -i {self.tempdir.name}/graph{graph_file_extension[graph_repr]} \
                            --anno-type {anno_repr} -o {self.tempdir.name}/annotation \
                            {TEST_DATA_DIR}/transcripts_100.fa'

        res = subprocess.run([annotate_command], shell=True)
        self.assertEqual(res.returncode, 0)

        # check annotation
        stats_annotation = self._get_stats('-a ' + f'{self.tempdir.name}/annotation{anno_file_extension[anno_repr]}')
        self.assertEqual(stats_annotation['returncode'], 0)
        self.assertEqual('100', stats_annotation['labels'])
        self.assertEqual(stats_graph['max index (k)'], stats_annotation['objects'])
        self.assertAlmostEqual(
            0.0185072 * (int(stats_graph['nodes (k)']) / int(stats_graph['max index (k)'])),
            float(stats_annotation['density']),
            places=6)
        self.assertEqual(anno_repr, stats_annotation['representation'])

    @parameterized.expand(GRAPH_TYPES)
    def test_annotate_coordinates(self, graph_repr):

        construct_command = f'{METAGRAPH} build --mask-dummy -p {NUM_THREADS} \
                              --graph {graph_repr} -k 11 \
                              -o {self.tempdir.name}/graph \
                              {TEST_DATA_DIR}/transcripts_100.fa'

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        # build annotation
        annotate_command = f'{METAGRAPH} annotate --anno-header -p {NUM_THREADS} --coordinates \
                            -i {self.tempdir.name}/graph{graph_file_extension[graph_repr]} \
                            -o {self.tempdir.name}/annotation \
                            {TEST_DATA_DIR}/transcripts_100.fa'

        res = subprocess.run([annotate_command], shell=True)
        self.assertEqual(res.returncode, 0)

    @parameterized.expand(GRAPH_TYPES)
    def test_annotate_coordinates_with_disk_swap(self, graph_repr):

        construct_command = f'{METAGRAPH} build --mask-dummy -p {NUM_THREADS} \
                              --graph {graph_repr} -k 11 \
                              -o {self.tempdir.name}/graph \
                              {TEST_DATA_DIR}/transcripts_100.fa'

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        # build annotation with disk swap
        annotate_command = f'{METAGRAPH} annotate --anno-header -p {NUM_THREADS} --coordinates \
                            --disk-swap {self.tempdir.name} --mem-cap-gb 1e-6 \
                            -i {self.tempdir.name}/graph{graph_file_extension[graph_repr]} \
                            -o {self.tempdir.name}/annotation \
                            {TEST_DATA_DIR}/transcripts_100.fa'

        res = subprocess.run([annotate_command], shell=True)
        self.assertEqual(res.returncode, 0)

        # build annotation in RAM
        annotate_command = f'{METAGRAPH} annotate --anno-header -p {NUM_THREADS} --coordinates \
                            -i {self.tempdir.name}/graph{graph_file_extension[graph_repr]} \
                            -o {self.tempdir.name}/annotation_ram \
                            {TEST_DATA_DIR}/transcripts_100.fa'

        res = subprocess.run([annotate_command], shell=True)
        self.assertEqual(res.returncode, 0)

        self.assertTrue(filecmp.cmp(f'{self.tempdir.name}/annotation.column.annodbg',
                                    f'{self.tempdir.name}/annotation_ram.column.annodbg'))
        self.assertTrue(filecmp.cmp(f'{self.tempdir.name}/annotation.column.annodbg.coords',
                                    f'{self.tempdir.name}/annotation_ram.column.annodbg.coords'))

if __name__ == '__main__':
    unittest.main()
