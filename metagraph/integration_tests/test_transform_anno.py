import unittest
from parameterized import parameterized, parameterized_class
import itertools
import subprocess
from subprocess import PIPE
from tempfile import TemporaryDirectory
import glob
import os
import gzip
from base import TestingBase, METAGRAPH, TEST_DATA_DIR, graph_file_extension, anno_file_extension
from helpers import get_test_class_name

"""Test operations on annotation columns"""

TEST_DATA_DIR = os.path.dirname(os.path.realpath(__file__)) + '/../tests/data'

NUM_THREADS = 4


class TestColumnOperations(TestingBase):
    @classmethod
    def setUpClass(cls):
        cls.tempdir = TemporaryDirectory()
        cls.graph_repr = 'succinct'
        cls.anno_repr = 'column'

        cls._build_graph(TEST_DATA_DIR + '/transcripts_100.fa',
                         cls.tempdir.name + '/graph',
                         20, cls.graph_repr, 'basic', '--mask-dummy')

        stats_graph = cls._get_stats(f'{cls.tempdir.name}/graph{graph_file_extension[cls.graph_repr]}')
        assert(stats_graph['returncode'] == 0)
        assert(stats_graph['k'] == '20')
        assert(stats_graph['nodes (k)'] == '46960')
        assert(stats_graph['mode'] == 'basic')

        cls.num_nodes = stats_graph['nodes (k)']
        cls.max_index = stats_graph['max index (k)']

        cls._annotate_graph(
            TEST_DATA_DIR + '/transcripts_100.fa',
            cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
            cls.tempdir.name + '/annotation',
            cls.anno_repr,
            extra_params='--count-kmers'
        )

    def setUp(self):
        # use relative paths to get rid of dir names
        self.old_cwd = os.getcwd()
        os.chdir(self.tempdir.name)
        self.annotation = f'annotation{anno_file_extension[self.anno_repr]}';

        # check annotation
        stats_annotation = self._get_stats(f'-a {self.annotation}')
        self.assertEqual(stats_annotation['returncode'], 0)
        self.assertEqual(stats_annotation['labels'], '100')
        self.assertEqual(stats_annotation['objects'], self.max_index)
        self.assertAlmostEqual(
            float(stats_annotation['density']),
            0.0185072 * int(self.num_nodes) / int(self.max_index),
            places=6)
        self.assertEqual(stats_annotation['representation'], self.anno_repr)

    def tearDown(self):
        os.chdir(self.old_cwd)

    def test_overlap(self):
        command = f'{METAGRAPH} transform_anno {self.annotation} \
            --compute-overlap {self.annotation} -o out -p {NUM_THREADS}'

        res = subprocess.run(command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(156421, len(res.stdout.decode()))

    def _check_aggregation_min(self, min_count, expected_density):
        command = f'{METAGRAPH} transform_anno {self.annotation} -p {NUM_THREADS} \
            --aggregate-columns --min-count {min_count} -o aggregated'

        res = subprocess.run(command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)

        stats_annotation = self._get_stats(f'-a aggregated{anno_file_extension[self.anno_repr]}')
        self.assertEqual(stats_annotation['returncode'], 0)
        self.assertEqual(stats_annotation['labels'], '1')
        self.assertEqual(stats_annotation['objects'], self.max_index)
        self.assertAlmostEqual(
            float(stats_annotation['density']),
            float(expected_density) * int(self.num_nodes) / int(self.max_index),
            places=5)
        self.assertEqual(stats_annotation['representation'], self.anno_repr)

    def test_aggregate_columns(self):
        self._check_aggregation_min(0, 1)
        self._check_aggregation_min(1, 1)
        self._check_aggregation_min(5, 0.0715077)
        self._check_aggregation_min(10, 0.00344974)
        self._check_aggregation_min(20, 0)

    def _check_aggregation_min_max_value(self, min_count, max_value, expected_density):
        command = f'{METAGRAPH} transform_anno {self.annotation} -p {NUM_THREADS} \
            --aggregate-columns --min-count {min_count} --max-value {max_value} -o aggregated'

        res = subprocess.run(command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)

        stats_annotation = self._get_stats(f'-a aggregated{anno_file_extension[self.anno_repr]}')
        self.assertEqual(stats_annotation['returncode'], 0)
        self.assertEqual(stats_annotation['labels'], '1')
        self.assertEqual(stats_annotation['objects'], self.max_index)
        self.assertAlmostEqual(
            float(stats_annotation['density']),
            float(expected_density) * int(self.num_nodes) / int(self.max_index),
            places=5)
        self.assertEqual(stats_annotation['representation'], self.anno_repr)

    def test_aggregate_columns_filtered(self):
        self._check_aggregation_min_max_value(0, 0, 0)
        self._check_aggregation_min_max_value(1, 0, 0)
        self._check_aggregation_min_max_value(2, 0, 0)
        self._check_aggregation_min_max_value(3, 0, 0)
        self._check_aggregation_min_max_value(5, 0, 0)

        self._check_aggregation_min_max_value(0, 1, 0.99704)
        self._check_aggregation_min_max_value(1, 1, 0.99704)
        self._check_aggregation_min_max_value(2, 1, 0.392994)
        self._check_aggregation_min_max_value(3, 1, 0.183305)
        self._check_aggregation_min_max_value(5, 1, 0.0715077)

        self._check_aggregation_min_max_value(0, 2, 0.998807)
        self._check_aggregation_min_max_value(1, 2, 0.998807)
        self._check_aggregation_min_max_value(2, 2, 0.394825)
        self._check_aggregation_min_max_value(3, 2, 0.183986)
        self._check_aggregation_min_max_value(5, 2, 0.0715077)

        self._check_aggregation_min_max_value(0, 5, 0.998999)
        self._check_aggregation_min_max_value(1, 5, 0.998999)
        self._check_aggregation_min_max_value(2, 5, 0.395315)
        self._check_aggregation_min_max_value(3, 5, 0.184817)
        self._check_aggregation_min_max_value(5, 5, 0.0715077)

        self._check_aggregation_min_max_value(0, 1000, 1)
        self._check_aggregation_min_max_value(1, 1000, 1)
        self._check_aggregation_min_max_value(2, 1000, 0.395336)
        self._check_aggregation_min_max_value(3, 1000, 0.184817)
        self._check_aggregation_min_max_value(5, 1000, 0.0715077)


if __name__ == '__main__':
    unittest.main()
